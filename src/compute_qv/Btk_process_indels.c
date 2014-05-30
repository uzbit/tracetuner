/**************************************************************************
 *
 * Copyright (c) 2004-2008, J. Craig Venter Institute. All rights reserved.
 *
 * This file is part of TraceTuner, the DNA sequencing quality value,
 * base calling and trace processing software.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

/*
 * Btk_process_indels.c
 */

/* static char fileIdentifier[] = "$Id: Btk_process_indels.c,v 1.15 2009/01/12 22:15:11 gdenisov Exp $"; */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <sys/param.h>
#include <stdint.h>

#include "Btk_qv.h"
#include "SFF_Toolkit.h"
#include "nr.h"
#include "util.h"
#include "Btk_lookup_table.h"
#include "Btk_qv_data.h"
#include "train.h"              /* needs util.h */
#include "Btk_process_peaks.h"  /* needs train.h */
#include "context_table.h"
#include "Btk_compute_qv.h"
#include "Btk_qv_funs.h"
#include "Btk_process_peaks.h"
#include "tracepoly.h"
#include "Btk_compute_tp.h"
#include "Btk_compute_tpars.h"  /* needs train.h */
#include "Btk_call_bases.h"
#include "Btk_match_data.h"
#include "Btk_qv_io.h"
#include "Btk_qv_funs.h"
#include "Btk_atod.h"
#include "Btk_process_raw_data.h"

#define ADD_SIGNAL_DOWNSTREAM                 0
#define CORR_FACTOR2                          40.
#define DEBUG                                 0
#define DP_WINDOW_SIZE                       30
#define DP_FRACTION                           0.6
#define EPS                                   0.0000001
#define H2_TO_AVE10_RATIO                     0.2
#define INDEL_RATIO                           0.2 
#define INDEL_RECALL_RATIO                    0.2
#define MAX_INDEL_SIZE                     1000
#define MAX_NUM_INDELS                        1
#define MAX_WIDTH_OF_PEAK                  1000
#define MIN_BASE_IND_START_INDEL            100
#define MIN_DIST_BETWEEN_INDELS              10
#define MIN_FRAC_MIXED_PEAKS_AHEAD            0.25
#define MIN_LENGTH_OF_INDEL                  10
#define MIN_LEN_INDEL                        30
#define MIN_LEN_STRETCH_PURE_PEAKS           15
#define MIN_LEN_SEARCH_STUTTER               12
#define MIN_LEN_STUTTER                       5
#define MIN_NUM_BASES_AHEAD                  20
#define MIN_PEAK_HEIGHT                     100
#define MIN_QV_JUMP                           7
#define MIN_WEIGHT                            0.8;
#define NUM_PARAMS                            4
#define NUM_PEAKS_CHECKED_AFTER_INDEL_START  20
#define NUM_QVS_UP                           21
#define NUM_QVS_DOWN                         21
#define POS_START_SEARCH_INDEL                5
#define SPACING_DEVIATION                     5
#define TYPE33_FRACTION                       0.2

typedef enum {
    DETECTION=1,             /* Dewtect indel */             
    SEPARATION,           /* Separate chromatograms in the indel region */              
} STEP;

typedef enum {
    LONG=1, 
    SHORT,
} CHROMAT_TYPE;

/*******************************************************************************
 * Function: Shape_idp
 * Purpose: compute the peak shape as a function of its model parameters
 *          and distance from the peak position
 *******************************************************************************
 */

double
Shape_idp(double A, double B, double x)
{
   double   x2 = x*x;
   if      (B <= EPS) { return A ; }
   else if (B > 1000) { return 0.; }

   return A * exp(-B*x2);
}

/*******************************************************************************
 * Function: is_dp_idp
 * Purpose:  check if a given peak is a dominant peak, that is,
 *           the observed peak signal at its position is greater than
 *           observed signal of any other color at this position
 *******************************************************************************
 */
int
is_dp_idp(Peak *pk, Data *data)
{
    int    j, height = pk->height,
           color  = pk->color_index,
           pos    = pk->pos;
    int    is_dp = 1;
    float  other_signal = -1., dominance_factor = DP_FRACTION;

    /* Determine if current peak is DP */
    for (j=0; j<NUM_COLORS; j++) {
        int new_pos = pos;

        if (j == color)
            continue;

        if ((new_pos <= data->pos_data_beg) || (new_pos >= data->pos_data_end))
            return (0);

        other_signal = (float)data->color_data[j].data[new_pos];

        if (other_signal > (float)height / dominance_factor) {
            is_dp = 0;
            break;
        }
    }

    return is_dp;
}

/*******************************************************************************
 * Function: colordata_find_peak_index_by_location_idp
 * Purpose: In the list of peaks of given color, find the peak such that a
 *          given position is located in its area. If the position does not
 *          fall into any peak's area, return peak_index equal to the negative
 *          index of the peak right to the given position. If there're no peaks
 *          to the right from the given position or the position is located
 *          left from the peak with index 0, return peak_index == NINF
 *******************************************************************************
 */
int
colordata_find_peak_index_by_location_idp(ColorData *color_data, int given_position,
    Peak *peak, int *peak_index)
{
    int lo, hi, mid;
    int beg = color_data->peak_list[0].beg;
    int end = color_data->peak_list[color_data->peak_list_len-1].end;

   *peak_index = -1;

    if (color_data->peak_list_len < 1) {
        return SUCCESS;
    }

    if ((given_position < beg) && (given_position > end))
    {
        return SUCCESS;
    }
    else {
        if (color_data->peak_list_len == 1) {
           *peak_index = 0;
           *peak = color_data->peak_list[*peak_index];
            return SUCCESS;
        }
    }

    lo = 0;
    hi = color_data->peak_list_len - 1;
    do {
        mid = (lo + hi) / 2;
        if (given_position < color_data->peak_list[mid].beg) {
            hi = mid - 1;
        }
        else if (given_position >= color_data->peak_list[mid].end) {
            lo = mid + 1;
        }
        else {
            *peak = color_data->peak_list[mid];
            *peak_index = mid;
            break;
        }
    } while (lo <= hi);

    return SUCCESS;
}

/*******************************************************************************
 * Function: get_highest_peaks
 * NOTE: take into account not only peaks found at the location of original 
 *       base but also peaks found in between the locations
 *******************************************************************************
 */
static void
get_highest_peaks(int base_ind, float *highest1_signal, float *highest2_signal, 
    int *pind1, int *pind2, int *color1, int *color2, int *pos1, int *pos2,
    Data *data, Options *options)
{
    int  i, peak_ind, loc, prev_loc, next_loc;
    float weight, min_weight, signal;
    float derivative;
    BtkMessage msg;
    Peak peak;
  
   *pind1  = -1;
   *pind2  = -1;   
   *color1 = -1;
   *color2 = -1;
   *pos1 = -1;
   *pos2 = -1;
   *highest1_signal = 0.;
   *highest2_signal = 0.;

    loc      = data->bases.coordinate[base_ind];
    prev_loc = (base_ind > 0) ?
               data->bases.coordinate[base_ind-1] : 0;
    next_loc = (base_ind < data->bases.length-1) ?
               data->bases.coordinate[base_ind+1] : data->length-1;

    /* Try peaks found at the original base locations */
    for (i=0; i<NUM_COLORS; i++)
    {
        ColorData *cd = &(data->color_data[i]);
        peak_ind = -1;

        colordata_find_peak_index_by_location(cd, loc, &peak, &peak_ind, &msg);

//      if ((peak_ind >= 0) && is_dye_blob(peak.pos, &peak, data, 0))
//          continue;

        if ((peak_ind > 0) && 
            (peak.ipos >  (loc + prev_loc)/2) &&
            (peak.ipos <= (loc + next_loc)/2))
            signal = peak.iheight;
        else
            signal = cd->data[loc];

        min_weight = MIN_WEIGHT;
        derivative = (1.-min_weight)/
            QVMAX((float)(loc-prev_loc)/2., (float)(next_loc-loc)/2.);

        if (*highest1_signal < signal)
        {
           *pos2     = *pos1;
           *pos1     =  
                ((peak_ind>0) &&
                 (peak.ipos >  (loc + prev_loc)/2) &&
                 (peak.ipos <= (loc + next_loc)/2)) ? peak.ipos : loc;
            weight   =  ((peak_ind>0) &&
                 (peak.ipos >  (loc + prev_loc)/2) &&
                 (peak.ipos <= (loc + next_loc)/2)) ? 
                (1.-derivative*(float)(abs(*pos1 - *pos2))) : 1.; 
           *highest2_signal = *highest1_signal * weight;
           *highest1_signal =  signal;                        
           *color2   = *color1;
           *color1   =  i;
           *pind2    = *pind1;
           *pind1    =  peak_ind;
        }
        else if ((*highest1_signal >= signal) && (*highest2_signal < signal))
        {
            *color2   = i;
            *pind2    = peak_ind;
            *pos2     = 
                 ((peak_ind>0) &&
                 (peak.ipos >  (loc + prev_loc)/2) &&
                 (peak.ipos <= (loc + next_loc)/2)) ? peak.pos : loc;
             weight   = ((peak_ind>0) &&
                 (peak.ipos >  (loc + prev_loc)/2) &&
                 (peak.ipos <= (loc + next_loc)/2)) ?
                (1.-derivative*(float)(abs(*pos1 - *pos2))) : 1.;
            *highest2_signal = (float)signal *weight;    
        }
#if 0
        fprintf(stderr, "    color=%d h=%d w=%f h1=%f h2=%f \n", i, peak.height, weight, *highest1_signal, *highest2_signal);
#endif
    }
}

/*******************************************************************************
 * Function: get_end_of_trace
 *******************************************************************************
 */
static int
get_end_of_trace(Data *data)
{
    int i, j, data_peak_ind_start, num_dps, num_type33, type;
    int base_index=-1;

    if (data->bases.length <= MIN_BASE_IND_START_INDEL)
        return data->length;

    data_peak_ind_start = 
        (data->bases.called_peak_list[MIN_BASE_IND_START_INDEL] != NULL) ?
         data->bases.called_peak_list[MIN_BASE_IND_START_INDEL]->data_peak_ind :
         data->bases.coordinate[MIN_BASE_IND_START_INDEL];

    for (i=data_peak_ind_start; i<data->peak_list_len; i++)
    {
        /* Count the fraction of dominant peaks of type 3 in a window of
         * 10 dominant peaks next to the current peak of type 3
         */
        if (data->peak_list[i]->is_called==0)
            continue;

        type = data->peak_list[i]->type;
        if (type !=33)
            continue;

        if (!is_dp_idp(data->peak_list[i], data))
            continue;

        if (data->peak_list[i]->base_index >=0)
            base_index = data->peak_list[i]->base_index;

        num_dps   = 1;
        num_type33 = 1;
        j = i;
#if 0
        fprintf(stderr, "trying peak at pos=%d, type=%d base_index=%d base=%c\n",
             data->peak_list[i]->pos, data->peak_list[i]->type,
             data->peak_list[i]->base_index, data->peak_list[i]->base);
#endif
        while ((num_dps < DP_WINDOW_SIZE) && (j < data->peak_list_len-1))
        {
            j++;
            if (is_dp_idp(data->peak_list[j], data))
            {
                num_dps++;
                type = data->peak_list[j]->type;
                if (type == 33)
                    num_type33++;
            }
        }
        if (num_type33 > DP_WINDOW_SIZE * TYPE33_FRACTION)
        {
#if 0
            int k;
            fprintf(stderr, "num_type3=%d end_pos=%d base_ind=%d\n",
                num_type3, data->peak_list[i]->pos, data->peak_list[i]->base_index);
            fprintf(stderr, "            Peaks found:\n");
            for (k=i; k< j; k++)
            {
                if (!is_dp_idp(data->peak_list[k], data))
                     continue;
                fprintf(stderr, "       Peak[%d] pos=%d type=%d base=%c is_called=%d base_ind=%d\n",
                    k, data->peak_list[k]->pos, data->peak_list[k]->type, data->peak_list[k]->base,
                    data->peak_list[k]->is_called, data->peak_list[k]->base_index);
            }
#endif
            data->bases.length = base_index+1;
            return data->peak_list[i]->pos;
        }
    }

    return data->length;
}


/*******************************************************************************
 * Function: get_autocorrelation
 * Purpose:  calculate and output correlation function
 *******************************************************************************
 */
static void
get_autocorrelation(int pos_indel_scan, int *chromatogram, 
    int num_datapoints, float *ans, int n)
{
    int    i;
    float *data1, *data2;

#if 0
    fprintf(stderr, "n=%d\n", n);
#endif
    data1 = CALLOC(float, n);
    data2 = CALLOC(float, n);

    for (i=0; i<n; i++)
    {
        if (i+pos_indel_scan < num_datapoints)
            data1[i] = (float)(chromatogram[i+pos_indel_scan]);
        else
            data1[i] = 0.;

        data2[i] = data1[i];
    }
    correl(data1-1, data2-1, n, ans-1);
    
    FREE(data1);
    FREE(data2);
}

/*******************************************************************************
 * Function: get_lag_scans 
 *******************************************************************************
 */
static int
get_lag_scans(int indel_scan, int shift_scans, float *fdata, int len, 
    float *max_peak_height, float *max_height_ratio)
{
    int   i, peak_ind=-1, max_peak_ind=-1, peak_pos=-1, max_peak_pos=-1;
    int   peak_beg = -1, peak_end = -1;
    float derivative_2, peak_height, height_ratio;
    int  *data;
 
   *max_peak_height  = -1;

    data = CALLOC(int, len);
    for (i=0; i<len; i++)
        data[i] = ROUND(fdata[i]);

    /* Find the index and the height of the highest peak */ 
    peak_height=-1.;
    height_ratio=-1.;
    peak_beg = -1;
    peak_end = -1;
    for (i = ((shift_scans > 0) ? QVMAX(7, shift_scans-6) : 0); 
         i < ((shift_scans > 0) ? QVMIN(shift_scans+6, len-1) : (len-1)); 
         i++) 
    {
        derivative_2 = fdata[i-1] - 2.*fdata[i] + fdata[i+1]; 

        if ((peak_beg < 0) && (derivative_2 < 0)) {
            peak_beg = i;
        }

        
        if ((i>1) && (peak_beg >= 0) && (peak_end < 0) &&
            is_maximum(i, data, len))
        {
            if (shift_scans > 0)
                peak_pos = i;
            else if ((shift_scans <= 0) && (peak_ind > 0))
                peak_pos = i;
        }

        /* Peak ends */
        if ((peak_beg >= 0) && (peak_end < 0) && (derivative_2 > 0))
        {
            peak_end = i;
            peak_ind++;
            if (peak_pos < 0)
                peak_pos = (peak_beg + peak_end)/2;
            peak_height = fdata[peak_pos];

            /* The highest peak */
            if ((shift_scans > 0) && (*max_peak_height < peak_height))
            {
                max_peak_pos    = peak_pos;
               *max_peak_height = peak_height;
            }
            else if ((shift_scans <= 0) && 
                     (peak_ind > 0) &&
                     (*max_peak_height < peak_height))
            {
                max_peak_ind    = peak_ind;
                max_peak_pos    = peak_pos;
               *max_peak_height = peak_height;
            }
            peak_beg = -1;
            peak_end = -1;
            peak_pos = -1;
        }
    }
       
#if 0
        fprintf(stderr, "i=%d der2=%f peak_ind=%d max_peak_ind=%d peak_hight=%f max_height=%f\n",
            i, derivative_2, peak_ind, max_peak_ind, peak_hight, max_height);
#endif
    /* Second, find the ratio of height of the highest peak and adjacent peak */
    peak_height=-1.;
    height_ratio=-1.;
    peak_beg = -1;
    peak_end = -1;
    peak_ind = -1;

    for (i = 1; i < len-1; i++)
    {
        derivative_2 = fdata[i-1] - 2.*fdata[i] + fdata[i+1];

        if ((peak_beg < 0) && (derivative_2 < 0)) {
            peak_beg = i;
        }

        if ((i>1)          && (peak_beg >= 0) && (peak_end < 0) &&
            (peak_ind > 0) && is_maximum(i, data, len))
            peak_pos = i;

        /* Peak ends */
        if ((peak_beg >= 0) && (peak_end < 0) && (derivative_2 > 0))
        {
            peak_end = i;
            peak_ind++;
            if (peak_pos < 0)
                peak_pos = (peak_beg + peak_end)/2;
            peak_height = fdata[peak_pos];

            /* Compare the highest peak with all other peaks of nonzero index */
            if (((peak_ind > 0) && (peak_ind == max_peak_ind+1)) ||
                ((peak_ind > 0) && (peak_ind == max_peak_ind-1)))
            {
                height_ratio = peak_height/(*max_peak_height);
                if (*max_height_ratio < height_ratio)
                    *max_height_ratio = height_ratio;
            }

#if 0
            fprintf(stderr, "  pind=%d, peak_h=%f max_ph=%f max_hr=%f \n",
                peak_ind, peak_height, *max_peak_height, *max_height_ratio);
#endif
            peak_beg = -1;
            peak_end = -1;
            peak_pos = -1;
        }
    }

#if 0
        fprintf(stderr, "i=%d der2=%f peak_ind=%d max_peak_ind=%d peak_hight=%f max_height=%f\n",
            i, derivative_2, peak_ind, max_peak_ind, peak_hight, max_height);
#endif

    FREE(data);
    return max_peak_pos;
}

/*******************************************************************************
 * Function: output_autocorrelation_functions
 *******************************************************************************
 */
static void
output_autocorrelation_functions(char *seq_name, float **ans, int n)
{
    int color, i;
    FILE *fp;
    
     /* Output normalized autocorrelations */
    fp = fopen("Autocorrelations.xgr", "w");
    fprintf(fp, "Title = Autocorrelations for %s\n", seq_name);
    fprintf(fp, "title_x = Lag  #\n");
    fprintf(fp, "title_y = Norm. autocorr. coeff\n");
    for (color=0; color<4; color++)
    {
        if (color == 0) fprintf(fp, "\ncolor = %d\n", 4);  /* green */
        if (color == 1) fprintf(fp, "\ncolor = %d\n", 3);  /* blue */
        if (color == 2) fprintf(fp, "\ncolor = %d\n", 7);  /* yellow */
        if (color == 3) fprintf(fp, "\ncolor = %d\n", 2);  /* red */
        if (color == 4) fprintf(fp, "\ncolor = %d\n", 9);  /* cyan */

        for (i=n/2; i<n; i++)
        {
            fprintf(fp, "%d %f\n", -(n-i), ans[color][i]);
        }
        for (i=0; i<n/2; i++)
        {
            fprintf(fp, "%d %f\n", i, ans[color][i]);
        }
        fprintf(fp, "next\n");
    }
    fclose(fp);
}

/*******************************************************************************
 * Function: correlate_data
 *******************************************************************************
 */
static int
count_identical_bases(int start, int end, int *c1, int *c2, float *h1, 
    float *h2, char *bases, float ratio)
{
    int i = start, n = 0, nmax = 0;
    float indel_ratio = (ratio > 0) ? ratio : INDEL_RATIO;
    if (start > end)
    {
        while(i > QVMAX(0, end))
        {
            if ((bases[i] == bases[i-1]) && 
                ((c2[i]==-1) || (h2[i]/h1[i] < indel_ratio)))
            {
                n++;
                if (nmax < n)
                    nmax = n;
            }
            else 
            {
                n = 0;
            }
            i--;
        }
    }
    else 
    {
        while (i < end) 
        {
            if ((bases[i] == bases[i+1]) &&
                ((c2[i]==-1) || (h2[i]/h1[i] < indel_ratio)))
            {
                n++;
                if (nmax < n)
                    nmax = n;
            }
            else 
            {
                n = 0;      
            }
            i++;
        }
    }
    return nmax;
}

/*******************************************************************************
 * Function: get_indel_size
 * Purpose:  determine the location and size of indels
 *******************************************************************************
 */
static int
get_autocorrelation_indel_size(int ind_loc, int shift_scans, int data_end_scan,
    char *seq_name, Data *data, Options options)
{
    int    i, j, color, best_color, n=1;
    int    lag_scans=0, best_lag_scans;
    float *ans[NUM_COLORS], height = 0., max_height = 0.; 
    float  height_ratio = 0., max_height_ratio = -1;

    if ((ind_loc <= 0) || (ind_loc >= data->length-1))
    {
        fprintf(stderr, "Position of indel start has not been determined\n");
        return -1;
    }
    else
    {
        /* Allocate memory */
#if 0
        data_end_scan = QVMIN(data_end_scan, indel1_start_scan + 600);
#endif
        while (n <  data_end_scan - ind_loc)
            n *= 2;

        for (i=0; i<NUM_COLORS; i++)
            ans[i] = CALLOC(float, n*2);

        /* ***********************************************************************
         * Calculating and normalize autocorrelation for each signal
         * ***********************************************************************
         */
        for (color=0; color<NUM_COLORS; color++)
        {
            float max = 0;

            get_autocorrelation(ind_loc, data->color_data[color].data,
                data_end_scan, ans[color], n);

            for (i=0; i<n/2; i++)
            {
               if (max < ans[color][i])
                   max = ans[color][i];
            }
            for (i=0; i<n; i++)
            {
               ans[color][i] /= max;
            }
        }

        best_color     = -1;
        best_lag_scans = 0;

        for (i=0; i<NUM_COLORS; i++)
        {
            /* Calculate the indel size */
            lag_scans = get_lag_scans(ind_loc, shift_scans, ans[i], n/2, &height,
                &height_ratio);

            if ((lag_scans > 0) && (max_height_ratio < height_ratio))
            {
                max_height = height;
                max_height_ratio = height_ratio;
                best_color = i;
                best_lag_scans = lag_scans;
            }

//           fprintf(stderr, "color=%d lag_scans=%d height_ratio=%f max_height_ratio=%f\n",
//              i, lag_scans, height_ratio, max_height_ratio);

            if (options.xgr)
                (void)output_autocorrelation_functions(seq_name, ans, n);
        }
    }
    for (j=0; j<4; j++)
        FREE(ans[j]);

    return best_lag_scans;
}


/*******************************************************************************
 * Function: correlated_colors
 *******************************************************************************
 */
static int
correlated_colors(int c1, int c2, int c1s, int c2s)
{
    if ((c1 == c1s) ||
        (c1 == c2s) ||
        (c2 == c1s) ||
       ((c2 == c2s) && (c2 >= 0)))
    {
        return 1;
    }
    return 0;
}

/*******************************************************************************
 * Function: correlated_3colors_and_heights
 * Purpose:  check if color of a peak at current position
 *           correlates with color of one of two peaks at a different position
 *******************************************************************************
 */
static int
correlated_3colors_and_heights(int c, int c1s, int c2s,
    float h1s, float h2s, float r2)
{
    if (  (c == c1s)                                     ||
        (((c == c1s) || (c == c2s)) && (h2s/h1s >= r2)))
    {
        return 1;
    }
    return 0;
}

/*******************************************************************************
 * Function: correlated_4colors_and_heights
 *******************************************************************************
 */
static int
correlated_4colors_and_heights(int c1, int c2, int c1s, int c2s,
    float h1, float h2, float h1s, float h2s, float r1, float r2)
{
    if (h2/h1 >= r1)
    {
        if (correlated_3colors_and_heights(c1, c1s, c2s, h1s, h2s, r2))
            return 1;

        if (correlated_3colors_and_heights(c2, c1s, c2s, h1s, h2s, r2))
            return 1;
    }
    else /* h2/h1 < r1 */
    {
        if (correlated_3colors_and_heights(c1, c1s, c2s, h1s, h2s, r2))
            return 1;
    }
    return 0;
}

/*******************************************************************************
 * Function: correlated_colors_and_heights
 *******************************************************************************
 */
static int
correlated_6colors_and_heights(
    int c1, int c2, int c1m, int c2m, int c1p, int c2p,
    float h1, float h2, float h1m, float h2m, float h1p, float h2p,
    float r1, float r2)
{
    if (h2/h1 >= r1)
    {
        if (correlated_3colors_and_heights(c1, c1m, c2m, h1m, h2m, r2) &&
            correlated_3colors_and_heights(c2, c1p, c2p, h1p, h2p, r2))
            return 1;

        if (correlated_3colors_and_heights(c2, c1m, c2m, h1m, h2m, r2) &&
            correlated_3colors_and_heights(c1, c1p, c2p, h1p, h2p, r2))
            return 1;
    }
    else if ((h2/h1 < r1) && (h2/h1 >= r2))
    {
        if (correlated_3colors_and_heights(c1, c1m, c2m, h1m, h2m, r2) &&
            correlated_3colors_and_heights(c2, c1p, c2p, h1p, h2p, r2))
            return 1;

        if (correlated_3colors_and_heights(c2, c1m, c2m, h1m, h2m, r2) &&
            correlated_3colors_and_heights(c1, c1p, c2p, h1p, h2p, r2))
            return 1;

        if (correlated_3colors_and_heights(c1, c1m, c2m, h1m, h2m, r2) &&
            correlated_3colors_and_heights(c1, c1p, c2p, h1p, h2p, r2))
            return 1;
    }
    else /* h2/h1 < r2 */
    {
        if (correlated_3colors_and_heights(c1, c1m, c2m, h1m, h2m, r2) &&
            correlated_3colors_and_heights(c1, c1p, c2p, h1p, h2p, r2))
            return 1;
    }
    return 0;
}

/*******************************************************************************
 * Function: correlated_4positions
 * Purpose:  check if there are no indel errors between two called positions
 *******************************************************************************
 */
static int
correlated_4positions(int s, int pos1, int pos2, int pos1s, int pos2s)
{
    int dev = SPACING_DEVIATION;
    int max_spacing = 12*s + dev;
    int min_spacing = 12*s - dev;
    int diff11 = abs(pos1 - pos1s), diff12 = abs(pos1 - pos2s),
        diff21 = abs(pos2 - pos1s), diff22 = abs(pos2 - pos2s);

    if ((diff11 >= min_spacing) && (diff11 <= max_spacing))
        return 1;

    if ((diff12 >= min_spacing) && (diff12 <= max_spacing))
        return 1;

    if ((diff21 >= min_spacing) && (diff21 <= max_spacing))
        return 1;

    if ((diff22 >= min_spacing) && (diff22 <= max_spacing))
        return 1;

    return 0;
}

/*******************************************************************************
 * Function: get_indel_shift
 * Purpose: determine the number of inserted/deleted bases
 *******************************************************************************
 */
static int
get_indel_shift(int i, int *c1, int *c2,  float *h1, float *h2,
    int *p1, int *p2, int len)
{
    int j;
    int s = 1;      /* shift */
    int num_peaks_checked_after_indel_start =
        NUM_PEAKS_CHECKED_AFTER_INDEL_START;
    int N = QVMIN(num_peaks_checked_after_indel_start -7, len-1-i);
    int best_indel_size = 0, min_num_unconfirmed = 10000;

    while ((s <= MAX_INDEL_SIZE)    &&
           (i + N + s < len))
    {
        int  num_confirmed = 0;
        int  num_correlated_in_signature = 0;

        num_peaks_checked_after_indel_start =
            QVMAX(num_peaks_checked_after_indel_start, s);
        N = QVMIN(num_peaks_checked_after_indel_start -1, len-1-i);
        for (j=i; j < i+N; j++)
        {
            float r1, r2;

            r1 = 0.2;
            r2 = 0.1;

            if (j+s >= len)
                continue;

            /* Forward confirmation */
            if (j - i < s)
            {
                if (!correlated_4positions(s, p1[j], p2[j], p1[j+s], p2[j+s]))
                {
                   N++;
                   num_confirmed++;
                   continue;
                }
                else if (correlated_4colors_and_heights(c1[j], c2[j], c1[j+s], c2[j+s],
                        h1[j], h2[j], h1[j+s], h2[j+s], r1, r2))
                {
                    num_confirmed++;
                    num_correlated_in_signature++;
                }
                else if ((N <= num_peaks_checked_after_indel_start) &&
                    (correlated_4colors_and_heights(
                         c1[j], c2[j], c1[j+s], c2[j+s],
                         h1[j], h2[j], h1[j+s], h2[j+s],
                         r1, r2) ||
                     correlated_4colors_and_heights(
                         c1[j], c2[j], c1[j+s], c2[j+s],
                         h1[j], h2[j], h1[j+s], h2[j+s],
                         r2, r2)))
                {
                    N++;
                    num_confirmed++;
                }

                else
                    break;
            }
            else
            /* Forward and backward confirmation */
            {
                if (!correlated_4positions(s, p1[j], p2[j], p1[j+s], p2[j+s]) ||
                    !correlated_4positions(s, p1[j], p2[j], p1[j-s], p2[j-s]))
                {
                   N++;
                   num_confirmed++;
                   continue;
                }

                if (correlated_6colors_and_heights(
                    c1[j], c2[j], c1[j-s], c2[j-s], c1[j+s], c2[j+s],
                    h1[j], h2[j], h1[j-s], h2[j-s], h1[j+s], h2[j+s], r1, r1))
                {
                    num_confirmed++;
                }
                else if (
                         (
                          N <= num_peaks_checked_after_indel_start)
                                                                           &&
                         correlated_colors(c1[j], c2[j], c1[j-s], c2[j-s]) &&
                         correlated_colors(c1[j], c2[j], c1[j+s], c2[j+s]))
                {
                    N++;
                    num_confirmed++;
                }
                else
                    break;
            }
        }
#if 0
        fprintf(stderr, "       i=%d loc=%d s=%d N=%d num_confirmed=%d\n",
            i, data->bases.coordinate[i]-6, s, N, num_confirmed);
#endif
        if (num_confirmed == N && num_correlated_in_signature == s)
        {
            return s;
        }

        if (min_num_unconfirmed > N - num_confirmed && 
            num_correlated_in_signature == s)
        {
            min_num_unconfirmed = N - num_confirmed;
            best_indel_size     = s;
        }

        s++;
    }
    
    if (min_num_unconfirmed <= 0)
        return best_indel_size;

    return -1;
}



/*******************************************************************************
 * Function: get_ave10_height           
 * Purpose: calculate an average height among 10 preceding highest peaks
 *******************************************************************************
 */
static float
get_ave10_height(int i, float *h)
{
    int j, n = QVMIN(10, i);
    float ave_h=0;

    if (i <= 0)
        return h[0];

    for (j=i; j>i-n; j--)
    {
        ave_h += h[j];
    }
    ave_h /= (float)n;

    return ave_h;
}

/*******************************************************************************
 * Function: get_stretch_pure_peaks             
 *******************************************************************************
 */
static int   
get_stretch_pure_peaks(int i, int *c1, int *c2, float *h1, float *h2, 
    Data *data, Options options)
{
    int   j;
    int   n = QVMIN(MIN_LEN_STRETCH_PURE_PEAKS, data->bases.length-i-1);
    float ave10_height = get_ave10_height(i-1, h1);
 
    for (j=i; j < i+n; j++)
    {
        if ((h2[j]/h1[j] > options.min_ratio) && 
            (h2[j]/ave10_height > H2_TO_AVE10_RATIO))
            return 0;      
    }
    return 1;
}

/*******************************************************************************
 * Function: compare_qvs 
 * Purpose:  Calculate the jump in average QV before and after an indel start
 *******************************************************************************
 */
static int      
compare_qvs(const void *qv1, const void *qv2)
{
    if (*(int *)qv1 < *(int *)qv2)
        return -1;
    else if (*(int *)qv1 > *(int *)qv2)
        return 1;

    return 0;
}


/*******************************************************************************
 * Function: get_qv_jump   
 * Purpose:  Calculate the jump in average QV before and after an indel start
 *******************************************************************************
 */
static float
get_qv_jump(int base_ind, uint8_t *quality_values, int num_bases)
{
    int    i, j=0, k=0, l=0, m=0;
    int   *qv_up, *qv_down;
    float  qv_jump, sum_qv_up = 0., sum_qv_down = 0.; 

    qv_up   = CALLOC(int, NUM_QVS_UP);
    qv_down = CALLOC(int, NUM_QVS_DOWN);

    for (i= base_ind-1; i>0 && j<NUM_QVS_UP; i--)
    {
        qv_up[j] = quality_values[i];
        j++;
    }
    qsort(qv_up, (size_t)j, sizeof(int), compare_qvs);
#if 0
    fprintf(stderr, "QVs up=\n");
    for(i=0; i<j; i++)
         fprintf(stderr, "%d ", qv_up[i]);
     fprintf(stderr, "\n");
#endif
    if (j >= 5)
    {
         for(i=(j-5)/2; i<j-(j-5)/2; i++)
         {
             sum_qv_up += qv_up[i];
             l++;        
         }
    }
    else
    {
        for(i=0; i<j; i++)
        {
            sum_qv_up += qv_up[i];
            l++;
        }
    }

    for (i= base_ind; i<num_bases && k< NUM_QVS_DOWN; i++)
    {
        if (i>=0)
        {
            qv_down[k] = quality_values[i];
            k++;
        }
    }
    qsort(qv_down, (size_t)k, sizeof(int), compare_qvs);
#if 0
    fprintf(stderr, "QVs down=\n");
    for(i=0; i<k; i++)
         fprintf(stderr, "%d ", qv_down[i]);
     fprintf(stderr, "\n");
#endif
    if (k >= 5)
    {
        for(i=(k-5)/2; i<k-(k-5)/2; i++)
        {
            sum_qv_down += qv_down[i];
            m++;
        }
    }
    else
    {
        for(i=0; i<k; i++)
        {
            sum_qv_down += qv_down[i];
            m++;
        }
    } 

    l = (l>0) ? l : 1;
    m = (m>0) ? m : 1;
    qv_jump = sum_qv_up/(float)l - sum_qv_down/(float)m;
//  if ((j==NUM_QVS_UP) && (k==NUM_QVS_DOWN))
//  {
//      qv_jump = qv_up[NUM_QVS_UP/2 + 1] - qv_down[NUM_QVS_DOWN/2 + 1];
#if 0
        fprintf(stderr, "qv_jump=%f qv_up=%f qv_down=%f\n", qv_jump,
            sum_qv_up/(float)l, sum_qv_down/(float)m);
//          (float)qv_up[NUM_QVS_UP/2 + 1], (float)qv_down[NUM_QVS_DOWN/2 + 1]);
#endif
//  }
    qv_jump = (qv_jump > 0) ? qv_jump : 0;

    FREE(qv_up);
    FREE(qv_down);
    return qv_jump;
}

/*******************************************************************************
 * Function: output_indels         
 *******************************************************************************
 */
static void
output_indels(char *seq_name, int *num_ident_bases, int *base_ind, int indcount, 
    int *indsize, int *indloc, char *bases, int num_bases, char *donor, 
    char *amplicon, uint8_t *quality_values, Options options) 
{
    int i, size;
    float qv_jump;

#if 0
    fprintf(stderr, "indcount=%d indsize= %d num_ident_bases= %d\n", 
        indcount, indsize[0], num_ident_bases[0]);
#endif
    qv_jump = get_qv_jump(*base_ind, quality_values, num_bases);

    if (indcount > 0)
    {
        for (i=0; i<indcount; i++)
        {
            if ((num_ident_bases[i] > MIN_LEN_STUTTER) && (indsize[i] == 12))
            {
                /* Stutter */
                fprintf(stdout,
                "%s: stutter_base_ind= %d lag_scans= %d loc= %d string= ?\n",
                    seq_name, base_ind[i], indsize[i], indloc[i]);
                if (options.Verbose != 1)
                    fprintf(stderr,"    donor= %s amplicon= %s qv_jump= %f\n",
                        donor, amplicon, qv_jump);
            }
            else /* Indel */
            {
                int j;
                char *string=NULL;
                size = (indsize[i]+3)/12;
                string = CALLOC(char, size+1);
                for (j=0; j<size; j++)
                {
                    string[j] = bases[base_ind[i] + j];
                } 
                string[size] = '\0';
                fprintf(stdout,
                "%s: indel_base_ind= %d lag_scans= %d loc= %d string= %s\n",
                    seq_name, base_ind[i], indsize[i], indloc[i], string);
                FREE(string);
                if (options.Verbose == 1)
                    fprintf(stdout,"\n");
                else
                    fprintf(stderr," donor= %s amplicon= %s qv_jump= %f\n", 
                         donor, amplicon, qv_jump);
            }
        }
    }
    else
    {
        if (options.Verbose == 1) { }
        else
            fprintf(stdout," donor= %s amplicon= %s qv_jump= %f\n", donor, amplicon,
               qv_jump);
    }
}

/*******************************************************************************
 * Function: detect_indels 
 * Purpose:  determine the location and size of indels 
 *******************************************************************************
 */
static void
detect_indels_from_chromatogram(int *indcount, int *num_ident_bases, 
    int **indloc, int **indbind, int **indsize, Data *data, 
    uint8_t *quality_values, Options options)
{
    int    i, stretch_of_pure_peaks_found = 0; 
    float *h1, *h2, *h3;                                   /* peak heights */
    float  ave10_height;
    int   *c1, *c2, *c3;                                   /* peak colors */
    int   *pi1, *pi2, *pi3;
    int   *pos1, *pos2, *pos3;
    int    shift, loc; 

    float  qv_jump;
    int    data_end_scan = data->bases.coordinate[data->bases.length-1];
    char  *seq_name;

    /* Allocate */
    c1      = CALLOC(int,   data->bases.length);
    c2      = CALLOC(int,   data->bases.length);
    c3      = CALLOC(int,   data->bases.length);
    pi1     = CALLOC(int,   data->bases.length);
    pi2     = CALLOC(int,   data->bases.length);
    pi3     = CALLOC(int,   data->bases.length);
    h1      = CALLOC(float, data->bases.length);
    h2      = CALLOC(float, data->bases.length);
    h3      = CALLOC(float, data->bases.length);
    pos1    = CALLOC(int,   data->bases.length);
    pos2    = CALLOC(int,   data->bases.length);
    pos3    = CALLOC(int,   data->bases.length);

   *indcount=0;

    /* Initialize */
    for (i=0; i< MAX_NUM_INDELS; i++)
      (*indloc)[i]  = (*indbind)[i] = (*indsize)[i] = -1;

    for (i=0; i< data->bases.length; i++)
    {
        c1[i] = c2[i] = pi1[i] = pi2[i] = -1;
        h1[i] = h2[i] = 0.;

        get_highest_peaks(i, &h1[i], &h2[i], &pi1[i], &pi2[i],
            &c1[i], &c2[i], &pos1[i], &pos2[i], data, &options);
#if 0
        fprintf(stderr, "i= %d c1=%d c2=%d \n", i, c1[i], c2[i]);
#endif
    }

    
    if ((seq_name = strrchr(options.file_name, '/')) != NULL) {
        seq_name++;
    }
    else {
        seq_name = options.file_name;
    }

    data->length = get_end_of_trace(data);

    /* Detect indel regions */
    for (i=POS_START_SEARCH_INDEL; i< data->bases.length; i++)
    {
        float ratio, ave_ratio;

        if ((pi1[i] < 0) || (pi2[i] < 0))
            continue;

        ave10_height = get_ave10_height(i-1, h1);
        ratio = h2[i]/h1[i];
        ave_ratio = h2[i]/ave10_height;

#if 0
        fprintf(stderr, "i=%d ratio= %f ave_ratio=%f\n", i, ratio, ave_ratio);
#endif
        if ((ratio     < INDEL_RATIO)        ||
            (ave_ratio < H2_TO_AVE10_RATIO))
        {
            if (stretch_of_pure_peaks_found == 0)
            {
                stretch_of_pure_peaks_found = get_stretch_pure_peaks(i, 
                    c1, c2, h1, h2, data, options);
                if (stretch_of_pure_peaks_found > 0) 
                {
                   *indcount   = 0;
                }
            }
            continue;
        }
        else
        {
            /* May be either SNP or indel */
            shift = get_indel_shift(i, c1, c2, h1, h2, pos1, pos2, 
                data->bases.length); 
            if (shift <= 0)
                continue;

            qv_jump = get_qv_jump(i, quality_values, data->bases.length);
            if (qv_jump <= MIN_QV_JUMP)
                continue;

            loc = data->bases.coordinate[i]-6;
            data_end_scan = QVMIN(data->bases.coordinate[data->bases.length-1],
                    loc + MIN_LEN_INDEL * 12);
            if (data_end_scan < loc+6)
                continue;
           
            if (*indcount >= MAX_NUM_INDELS)
                continue;
   
            num_ident_bases[*indcount] = count_identical_bases(i,
                QVMAX(i - MIN_LEN_SEARCH_STUTTER - 2, 0),
                c1, c2, h1, h2, data->bases.bases, options.min_ratio);

            /* Indel */
            if ((num_ident_bases[*indcount] < MIN_LEN_STUTTER) || 
                ((shift*12+2)/12 > 1))
            {
                /* Use autocorrelation curve to determine the indel size */
                (*indloc)[*indcount]  = loc;
                (*indsize)[*indcount] = get_autocorrelation_indel_size(loc, shift,
                        data_end_scan, options.file_name, data, options);
                if ((*indsize)[*indcount] <= 0)
                    (*indsize)[*indcount] = shift*12;
                (*indbind)[*indcount] = i;
                (*indcount)++;
            }
            else 
            {
                (*indloc)[*indcount]  = loc;                             
                (*indsize)[*indcount] = 12;
                (*indbind)[*indcount] = i;
                (*indcount)++;
            }
        }
    }

    if ((*indcount == 0) && ((*indbind)[0] > 0))
        (*indcount)++;

    FREE(c1);
    FREE(c2);
    FREE(c3);
    FREE(pi1);
    FREE(pi2);
    FREE(pi3);
    FREE(h1);
    FREE(h2);
    FREE(h3);
    FREE(pos1);
    FREE(pos2);
    FREE(pos3);

    return;
}

/*******************************************************************************
 * Function: get_intrinsic_signal             
 *******************************************************************************
 */
static void
get_intrinsic_signal(int h, float B, int *y, int *num_points, Options *options)
{
    int i;
    for (i=0; i<MAX_WIDTH_OF_PEAK; i++)
        y[i] = 0;

   *num_points = 0;
//  y[*num_points] = peak->height;
    y[*num_points] = h;
    while ((y[*num_points] >= 1) && ((*num_points) < MAX_WIDTH_OF_PEAK))
    {
        (*num_points)++;    
        y[*num_points] = INT_DBL(Shape_idp(h, B, *num_points));
    } 
}

/*******************************************************************************
 * Function: add_intrinsic_signal
 *******************************************************************************
 */
static int
add_intrinsic_signal(int pos, int *y, int num_points, 
    float frac, int *data, int len)
{
    int i;

    for (i=0; i<num_points; i++)
    {
        if (pos + i < len)
        {
            data[pos + i] += ROUND((float)y[i] * frac);
        }

        if ((i>0) && (pos - i >= 0))
        {
            data[pos - i] += ROUND((float)y[i] * frac);
        }
    }

    return SUCCESS;
}

/*******************************************************************************
 * Function: subtract_intrinsic_signal
 *******************************************************************************
 */
static int
subtract_intrinsic_signal(int pos, int *y, int num_points, float frac,
    int *data, int len)
{
    int i;

    for (i=0; i<num_points; i++)
    {
        if (i >= MAX_WIDTH_OF_PEAK)
            continue;

        if (pos + i < len)
        {
            data[pos + i] -= (int)(ROUND((float)y[i] * frac));
            if (data[pos + i] < 0)
                data[pos + i] = 0;
        }

        if ((i>0) && (pos - i >= 0))
        {
            data[pos - i] -= (int)ROUND((float)y[i] * frac);
            if (data[pos - i] < 0)
                data[pos - i] = 0;
        }
    }

    return SUCCESS;
}

/*******************************************************************************
 * Function: base2color
 *******************************************************************************
 */
static int
base2color(char base)
{
    if ((base == 'A') || (base == 'a')) { return 0; }
    if ((base == 'C') || (base == 'c')) { return 1; }
    if ((base == 'G') || (base == 'g')) { return 2; }
    if ((base == 'T') || (base == 't')) { return 3; }

    return -1;
}


static float
get_frac(int loc, int locm, int locp, ColorData *cd)
{
    int denom = (sqrt((float)cd->data[loc] * (float)cd->data[locm]) +
           sqrt((float)cd->data[loc] * (float)cd->data[locp]));
    if (denom > 0)
        return sqrt((float)cd->data[loc] * (float)cd->data[locm]) /
              (sqrt((float)cd->data[loc] * (float)cd->data[locm]) +
               sqrt((float)cd->data[loc] * (float)cd->data[locp]));
    else
        return 1;
}

/*******************************************************************************
 * Function: produce_new chromatogram
 *******************************************************************************
 */
static int 
produce_new_chromatogram(CHROMAT_TYPE type, int indsize_scans, int indbind,
    char *color2base, int **new_chromatogram, Data *data, 
    char *called_seq, Options options)
{
    int   i, *y, num_points=0;
    int   color1, color2, pind1, pind2, pos1, pos2;
    int   right_color, other_color, right_pind, other_pind;
    int   max_width_of_ipeak = MAX_WIDTH_OF_PEAK;
    float highest1_signal, highest2_signal; 
    float right_height, other_height, nearest_B= -1., nearest_h;
    Peak *peak;
    ColorData *cd;

    y = CALLOC(int, max_width_of_ipeak);

    /* Loop through all the called base locations */
    for (i=indbind; i<data->bases.length; i++)
    { 
        char called_base = called_seq[i];

        if (called_base == 'N')
        {
            continue;
        }

        color1 = -1;
        color2 = -1;
        get_highest_peaks(i, &highest1_signal, &highest2_signal,
            &pind1, &pind2, &color1, &color2, &pos1, &pos2, data,
            &options);

        if (highest2_signal/highest1_signal < INDEL_RATIO)
            color2 = -1;
 
//      fprintf(stderr, "i=%d color1=%d color2=%d\n", i, color1, color2); 
        right_color = base2color(called_base);
        other_color = (right_color == color1) ? color2 : color1;
        right_pind  = (right_color == color1) ? pind1  : pind2;
        other_pind  = (right_color == color1) ? pind2  : pind1;
        right_height= (right_color == color1) ? highest1_signal : highest2_signal;
        other_height= (right_color == color2) ? highest1_signal : highest2_signal;

        if ((right_color != color1) && (right_color != color2))
        {
            if (DEBUG)
            {
                fprintf(stderr, "i=%d color1=%d color2=%d\n", i,  color1, color2);
                fprintf(stderr, 
                    "i=%d base=%c chromat_type=%d right_color=%d other_color=%d color1=%d color2=%d\n",
                    i, called_base, type, right_color, other_color, color1, color2);
                fprintf(stderr, "Error in called base color!\n");
                return ERROR;
            }
            else
                continue;
        }

        /* Raise the right_peak */
        cd = &data->color_data[right_color];
        pos1 = data->bases.coordinate[i];
        if (right_pind >= 0)
        {
            peak = &cd->peak_list[right_pind]; 
            nearest_B = 1./peak->beta/peak->beta;
            pos1 = peak->pos;
        }
        nearest_h = (right_pind >= 0) ? peak->height : right_height;
        get_intrinsic_signal(nearest_h, nearest_B, y, &num_points, &options);
        if (add_intrinsic_signal(pos1, y, num_points,
            1.0, new_chromatogram[right_color], cd->length) != SUCCESS)
            return ERROR;

        if (other_color < 0)
            continue;

        /* Remove the other peak */
        cd = &data->color_data[other_color];
        pos2 = data->bases.coordinate[i];
        if (other_pind >= 0)
        {
            peak = &cd->peak_list[other_pind];
            nearest_B = 1./peak->beta/peak->beta;
            pos2 = peak->pos;
        }
        nearest_h = (other_pind >= 0) ? peak->height : other_height; 
        get_intrinsic_signal(nearest_h, nearest_B, y, &num_points, &options);
        if (subtract_intrinsic_signal(pos2, y, num_points, 
            1.0, new_chromatogram[other_color], cd->length) != SUCCESS)
            return ERROR;     
    }


    /* Process undetermined bases (those called 'N') */
    for (i=indbind; i<data->bases.length; i++)
    {
        int loc1, loc2, loc1m, loc2m, loc1p, loc2p;
        float frac;

        if (called_seq[i] != 'N')
            continue;

        color1 = -1;
        color2 = -1;

        get_highest_peaks(i, &highest1_signal, &highest2_signal,
            &pind1, &pind2, &color1, &color2, &pos1, &pos2, data,
            &options);

        if (highest2_signal/highest1_signal < INDEL_RATIO)
            color2 = -1;

        if (color1 < 0)
            continue;

        cd = &data->color_data[color1];
        loc1 = data->bases.coordinate[i];
        nearest_h = highest1_signal;
        if (pind1 >= 0)
        {
            peak = &cd->peak_list[pind1];
            nearest_h = peak->height;
            nearest_B = 1./peak->beta/peak->beta;
            loc1 = peak->pos;
        }
        loc1m = loc1 - indsize_scans;
        loc1p = loc1 + indsize_scans;
        if (type == LONG)
        {
            frac = get_frac(loc1, loc1m, loc1p, cd);
            if (frac < get_frac(loc1+3, loc1m+3, loc1p+3, cd))
                frac = get_frac(loc1+3, loc1m+3, loc1p+3, cd);
            if (frac < get_frac(loc1-3, loc1m-3, loc1p-3, cd))
                frac = get_frac(loc1-3, loc1m-3, loc1p-3, cd);
        }
        else
        {
            frac = 1. - get_frac(loc1, loc1m, loc1p, cd);
            if (frac < 1. - get_frac(loc1+3, loc1m+3, loc1p+3, cd))
                frac = 1. - get_frac(loc1+3, loc1m+3, loc1p+3, cd);
            if (frac < 1. - get_frac(loc1-3, loc1m-3, loc1p-3, cd))
                frac = 1. - get_frac(loc1-3, loc1m-3, loc1p-3, cd);
        }
        get_intrinsic_signal(nearest_h, nearest_B, y, &num_points, &options);
        (void)subtract_intrinsic_signal(loc1, y, num_points,
                1.0, new_chromatogram[color1], cd->length);
        (void)add_intrinsic_signal(loc1, y, num_points,
                2.*frac, new_chromatogram[color1], cd->length);

        if (color2 < 0) 
            continue;

        cd = &data->color_data[color2];
        loc2 = data->bases.coordinate[i];
        nearest_h = highest2_signal;
        if (pind2 >= 0)
        {
            peak = &cd->peak_list[pind2];
            nearest_h = peak->height;
            nearest_B = 1./peak->beta/peak->beta;
            loc2 = peak->pos;
        }
        loc2m = loc2 - indsize_scans;
        loc2p = loc2 + indsize_scans;
        
        if (type == LONG)
        {
            frac = get_frac(loc2, loc2m, loc2p, cd);
            if (frac < get_frac(loc2+3, loc2m+3, loc2p+3, cd))
                frac = get_frac(loc2+3, loc2m+3, loc2p+3, cd);
            if (frac < get_frac(loc2-3, loc2m-3, loc2p-3, cd))
                frac = get_frac(loc2-3, loc2m-3, loc2p-3, cd);
        }
        else
        {
            frac = 1. - get_frac(loc2, loc2m, loc2p, cd);
            if (frac < 1. - get_frac(loc2+3, loc2m+3, loc2p+3, cd))
                frac = 1. - get_frac(loc2+3, loc2m+3, loc2p+3, cd);
            if (frac < 1. - get_frac(loc2-3, loc2m-3, loc2p-3, cd))
                frac = 1. - get_frac(loc2-3, loc2m-3, loc2p-3, cd);
        }
        get_intrinsic_signal(nearest_h, nearest_B, y, &num_points, &options);
        (void)subtract_intrinsic_signal(loc2, y, num_points,
                1.0, new_chromatogram[color2], cd->length);
        (void)add_intrinsic_signal(loc2, y, num_points,
                2.*frac, new_chromatogram[color2], cd->length);
    }


    FREE(y);
    return SUCCESS;
}


/*******************************************************************************
 * Function: is_correlated
 * Purpose:  check if at two specified positions there are peaks 
 *           of the same color
 *******************************************************************************
 */
static int
is_correlated(int i, int m, int color, int *color1, int *color2,
    int pos, int *pos1, int *pos2, int len)
{
    int dev = SPACING_DEVIATION;

    if ((i+m < 0) || (i+m >= len))
        return 0;

    if (color < 0)
        return 0;

    if (color == color1[i+m])
    {
        if ((abs(m) == 1) &&
            (abs(pos - pos1[i+m]) <= 12*abs(m) + dev))
        {
            return 1;
        }
        if ((abs(m) > 1) &&
            (abs(pos - pos1[i+m]) >= 12*abs(m) - dev) &&
            (abs(pos - pos1[i+m]) <= 12*abs(m) + dev))
        {
            return 1;
        }
    }

    if (color == color2[i+m])
    {
        if ((abs(m) == 1) &&
            (abs(pos - pos2[i+m]) <= 12*abs(m) + dev))
        {
            return 1;
        }

        if ((abs(m) > 1) &&
            (abs(pos - pos2[i+m]) >= 12*abs(m) - dev) &&
            (abs(pos - pos2[i+m]) <= 12*abs(m) + dev))
        {
            return 1;
        }
    }
    return 0;
}

/*******************************************************************************
 * Function: make_first_pass    
 * Purpose: process the locations where there is only one choice
 *******************************************************************************
 */
static void
make_first_pass(int ibind, int m, int len, int *color1, int *color2, 
    int *pos1, int *pos2, int *called_color_short, int *called_color_long, 
    char *color2base)
{
    int i;

    /* Forward pass */
    for (i=ibind; i<len; i++)
    {

        if (called_color_long[i] >= 0)
            continue;

        if ((i>=ibind+m) && (called_color_short[i-m] >= 0))
            continue;

        if (color2[i] < 0)
        {
            /* Peak of only one color - there's no choice */
            called_color_long[i] = color1[i];
            if ((i>=ibind+m) && 
                is_correlated(i, -m, color1[i], color1, color2, pos1[i], pos1, pos2, len))
                called_color_short[i-m] = color1[i];
            continue;
        }
        else  /* color1 >= 0 and color2 >= 0 */
        {
            if ((i>=ibind+m) &&
                is_correlated(i, -m, color1[i], color1, color2, pos1[i], pos1, pos2, len) &&
               !is_correlated(i, -m, color2[i], color1, color2, pos2[i], pos1, pos2, len))
            {
                called_color_long[i]    = color1[i];
                called_color_short[i-m] = color1[i];
                continue;
            }
            else
            if ((i>=ibind+m) &&
                is_correlated(i, -m, color2[i], color1, color2, pos2[i], pos1, pos2, len) &&
               !is_correlated(i, -m, color1[i], color1, color2, pos1[i], pos1, pos2, len))
            {
                called_color_long[i]    = color2[i];
                called_color_short[i-m] = color2[i];
                continue;
            }
            else
                continue;
        }
    }

    /* Reverse pass */
    for (i=len-1; i>=ibind; i--)
    {
        if (called_color_short[i] >= 0)
            continue;

        if ((i<len-m) && (called_color_long[i+m] >= 0))
            continue;

        if (color2[i] < 0)
        {
            /* Peak of only one color - there's no choice */
            called_color_short[i] = color1[i];
            if ((i<len-m) &&
                is_correlated(i, m, color1[i], color1, color2, pos1[i], pos1, pos2, len))
                called_color_long[i+m] = color1[i];  // MAY BE A PROBLEM !!!
            continue;
        }

        if (i<len-m)
        {
            if (is_correlated(i, m, color1[i], color1, color2, pos1[i], pos1, pos2, len) &&
               !is_correlated(i, m, color2[i], color1, color2, pos2[i], pos1, pos2, len))
            {
                called_color_short[i]   = color1[i];
                called_color_long[i+m]  = color1[i];
                continue;
            }
            else
            if (is_correlated(i, m, color2[i], color1, color2, pos2[i], pos1, pos2, len) &&
               !is_correlated(i, m, color1[i], color1, color2, pos1[i], pos1, pos2, len))
            {
                called_color_short[i]   = color2[i];
                called_color_long[i+m]  = color2[i];
                continue;
            }
        }

    }  /* end of the 1st pass */
}

/*******************************************************************************
 * Function: make_second_pass
 * Purpose: process the locations where there are two choices
 *          with one of the peaks called
 *******************************************************************************
 */
static void
make_second_pass(int ibind, int m, int len, int *color1, int *color2, 
    int *pos1, int *pos2, int *called_color_short, int *called_color_long,
    float *highest1_signal, float *highest2_signal, char *color2base)
{
    int i;

    /* Forward pass */
    for (i=ibind; i<len; i++)
    {
        if (called_color_long[i] >= 0)
        {
            if ((called_color_long[i] != color1[i]) &&
                (called_color_long[i] != color2[i]))
            fprintf(stderr, "In make_second_pass: called_color_long=%d color1=%d color2=%d\n",
                   called_color_long[i], color1[i], color2[i]); 
            continue;
        }

        if ((i>=ibind+m) && (called_color_short[i-m] >= 0))
            continue;

        if ((i>=ibind+m) &&
            is_correlated(i, -m, color1[i], color1, color2, pos1[i], pos1, pos2, len) &&
            is_correlated(i, -m, color2[i], color1, color2, pos2[i], pos1, pos2, len))
        {
            if (called_color_long[i-m] >= 0 && called_color_long[i-m] == color1[i])
            {
                if (color2[i] >= 0)
                {
                    called_color_long[i]    = color2[i];
                    called_color_short[i-m] = color2[i];
                    continue;
                }
                else
                {
                    called_color_long[i]    = color1[i];
                    called_color_short[i-m] = color1[i];
                }
            }
            if (called_color_long[i-m] >= 0 && called_color_long[i-m] == color2[i])
            {
                called_color_long[i]    = color1[i];
                called_color_short[i-m] = color1[i];
                continue;
            }
        }
    }  

    /* Reverse pass */ 
    for (i=len-1; i>=ibind; i--)
    {
        if (called_color_short[i] >= 0)
            continue;

        if ((i<len-m) && (called_color_long[i+m] >= 0))
            continue;

        if ((i<len-m) &&
            is_correlated(i, m, color1[i], color1, color2, pos1[i], pos1, pos2, len) &&
            is_correlated(i, m, color2[i], color1, color2, pos2[i], pos1, pos2, len))
        {
            if (called_color_short[i+m] == color1[i])
            {
                if (color2[i] >= 0)
                {
                    called_color_short[i]  = color2[i];
                    called_color_long[i+m] = color2[i];
                    continue;
                }
                else
                {
                    called_color_short[i]  = color1[i];
                    called_color_long[i+m] = color1[i];
                    continue;
                }
            }
            else if (called_color_short[i+m] >= 0 && called_color_short[i+m] == color2[i])
            {
                called_color_short[i]  = color1[i];
                called_color_long[i+m] = color1[i];
                continue;
            }
        }
    } /* end of second pass */
}


static void
output_undetermined_bases_in_long_read(int m, int len, int *color1, 
    int *color2, int *pos1, int *pos2, int *called_color_long)
{
    int i, k=0;
    fprintf(stderr, "Long read:\n");
    for (i=0; i<len; i++)
    {
        if (called_color_long[i] >= 0)
            continue;
        k++;
        fprintf(stderr,
        "i=%3d c1-m=%2d c2-m=%2d called-m=%2d c1=%2d c2=%2d c1+m=%2d c2+m=%2d called+m=%2d  ",
               i,
               (i-m<0) ? -1 : color1[i-m], (i-m<0) ? -1 : color2[i-m],
               (i-m<0) ? -1 : called_color_long[i-m],
               color1[i], color2[i],
               (i+m>len-1) ? -1 : color1[i+m], (i+m>len-1) ? -1 : color2[i+m],
               (i+m>len-1) ? -1 : called_color_long[i+m]);
        if ((i>=m) &&
            !is_correlated(i, -m, color1[i], color1, color2, pos1[i], pos1, pos2, len) &&
            !is_correlated(i, -m, color2[i], color1, color2, pos2[i], pos1, pos2, len))
        {
            fprintf(stderr, "no backward corr\n");
        }
        else {
            fprintf(stderr, "\n");
        }
    }
    fprintf(stderr, "%d bases undetermined in the long  read\n", k);
}

static void
output_undetermined_bases_in_short_read(int m, int len, int *color1,
    int *color2, int *pos1, int *pos2, int *called_color_short)
{
    int i, k = 0;
    fprintf(stderr, "Short read:\n");
    for (i=0; i<len; i++)
    {
        if (called_color_short[i] >= 0)
            continue;

        k++;
        fprintf(stderr,
        "i=%3d c1-m=%2d c2-m=%2d called-m=%2d c1=%2d c2=%2d c1+m=%2d c2+m=%2d called+m=%2d  ",
               i,
               (i-m<0) ? -1 : color1[i-m], (i-m<0) ? -1 : color2[i-m],
               (i-m<0) ? -1 : called_color_short[i-m],
               color1[i], color2[i],
               (i+m>len-1) ? -1 : color1[i+m], (i+m>len-1) ? -1 : color2[i+m],
               (i+m>len-1) ? -1 : called_color_short[i+m]);
        if ((i<len-m) &&
            !is_correlated(i, m, color1[i], color1, color2, pos1[i], pos1, pos2, len) &&
            !is_correlated(i, m, color2[i], color1, color2, pos2[i], pos1, pos2, len))
        {
            fprintf(stderr, "no forward  corr\n");
        }
        else fprintf(stderr, "\n");
    }
    fprintf(stderr, "%d bases undetermined in the short read\n", k);
}

static void
compare_long_and_short_reads(int m, int len, int *color1,
    int *color2, int *pos1, int *pos2, 
    int *called_color_short, int *called_color_long)
{
    int i, k = 0;
    fprintf(stderr, "\nComparison of long and short reads:\n");
    for (i=0; i<len-m; i++)
    {
        if ((called_color_short[i] >= 0) &&
            (called_color_short[i] == called_color_long[i+m]))
            continue;
        k++;
        fprintf(stderr,
        "i=%3d called_color_short[i] = %d called_color_long[i+m] = %d\n",
         i, called_color_short[i], called_color_long[i+m]);
    }
    fprintf(stderr, "%d different bases in the long and short reads\n", k);
}


/*******************************************************************************
 * Function: call_long_and_short_sequences            
 * Purpose:  determine the sequence of called bases    
             m - indel size
 *******************************************************************************
 */
static int
call_long_and_short_sequences(int ibind, int isize, char *color2base,
    char *long_called_seq, char *short_called_seq,
    Data *data, Options options)
{
    int    i, len;
    int    indsize_bases = (isize+3)/12.;
    char  *base;
    int   *color1, *color2, *pind1, *pind2, *pos1, *pos2;
    int   *called_color_long, *called_color_short;
    float *highest1_signal, *highest2_signal;

    /* Determine the number of base positions in the indel region */
    len = data->bases.length;

    for (i=0; i<data->bases.length; i++)
    {
        if (i<ibind)
            long_called_seq[i] = short_called_seq[i] = data->bases.bases[i];
        else
            long_called_seq[i] = short_called_seq[i] = 'N';
    }
    long_called_seq[data->bases.length]  = '\0';
    short_called_seq[data->bases.length] = '\0';

    /* Allocate and initialize the arrays of bases and peak colors */
    base   = CALLOC(char, len+1);
    color1 = CALLOC(int,  len);     /* highest peak */
    color2 = CALLOC(int,  len);     /* 2nd highest peak */
    called_color_long  = CALLOC(int,  len);
    called_color_short = CALLOC(int,  len);
    pind1  = CALLOC(int,  len);     /* highest peak */
    pind2  = CALLOC(int,  len);     /* 2nd highest peak */
    highest1_signal = CALLOC(float,len);
    highest2_signal = CALLOC(float,len);
    pos1  = CALLOC(int,  len);
    pos2  = CALLOC(int,  len);


    /* Initialize all bases to 'N' */
    for (i=0; i<len; i++)
    {
        base[i]  = 'N';
        called_color_long[i]  = -1;
        called_color_short[i] = -1;
        get_highest_peaks(i, &(highest1_signal[i]), &(highest2_signal[i]),
            &(pind1[i]), &(pind2[i]), &(color1[i]), &(color2[i]),
            &(pos1[i]), &(pos2[i]), data, &options);

        if (highest2_signal[i]/highest1_signal[i] < INDEL_RATIO)
            color2[i] = -1;
#if 0
    fprintf(stderr, "i= %d colors= %d %d bases= %c %c \n", i, color1[i], color2[i],
     color1[i]>= 0 ? color2base[color1[i]] : 'N', color2[i]>= 0 ? color2base[color2[i]] : 'N');
#endif
    }
    base[len] = '\0';
#if 0
    fprintf(stderr, "ibind= %d indsize_bases= %d len= %d data->bases.length= %d\n", 
        ibind, indsize_bases, len, data->bases.length);
    fprintf(stderr, "Before the first pass long seq=\n%s\n", long_called_seq+ibind);
#endif

    /* First pass: only one color at location or only one correlation backwards */
    make_first_pass(ibind, indsize_bases, len, color1, color2, pos1, pos2, 
        called_color_short, called_color_long, color2base);

#if 0
    for (i=0; i<len; i++)
    {
        if (called_color_long[i] > 0 && called_color_long[i] < 4)
            long_called_seq[i] = color2base[called_color_long[i]];
        if (called_color_short[i] > 0 && called_color_short[i] < 4)
           short_called_seq[i] = color2base[called_color_short[i]]; 
    }    
    fprintf(stderr, "After the first pass long seq=\n%s\n", long_called_seq+ibind);
    fprintf(stderr, "After the first pass short seq=\n%s\n", short_called_seq+ibind);
#endif
#if 0
   {
       int j=0, k=0;
       for  (i=0; i<len; i++)
       {
           if (called_color_long[i] < 0)
               j++;
           if (called_color_short[i] < 0)
               k++;
       }
       fprintf(stderr, "Num of uncalled bases after the 1st pass: long = %d short = %d\n",
         j, k);
    }
#endif
    /* 2nd pass: two correlations backward; call long read */
     make_second_pass(ibind, indsize_bases, len, color1, color2, pos1, pos2, 
        called_color_short, called_color_long, highest1_signal, highest2_signal, 
        color2base);
#if 0
    for (i=0; i<len; i++)
    {
        if (called_color_long[i] > 0 && called_color_long[i] < 4)
            long_called_seq[i] = color2base[called_color_long[i]];
        if (called_color_short[i] > 0 && called_color_short[i] < 4)
           short_called_seq[i] = color2base[called_color_short[i]];
    }
    fprintf(stderr, "After second pass long seq=\n%s\n", long_called_seq+ibind);
    fprintf(stderr, "After second pass short seq=\n%s\n", short_called_seq+ibind);
#endif

#if 0

    {
       int j=0, k=0;
       for  (i=0; i<len; i++)
       {
           if (called_color_long[i] < 0)
               j++;
           if (called_color_short[i] < 0)
               k++;
       }
       fprintf(stderr, "Num of uncalled bases after the 2nd pass: long = %d short = %d\n",
         j, k);
    }

#endif

    /* Determine the color of deleted bases */
    for (i=ibind; i<ibind+indsize_bases; i++)
    {
        if (color2[i] < 0)
            called_color_long[i] = color1[i];
        else
        {
           if (called_color_short[i] == color1[i] && called_color_long[i] < 0)
               called_color_long[i]   = color2[i];
           if (called_color_short[i] == color2[i] && called_color_long[i] < 0)
               called_color_long[i]   = color1[i];
        }
#if 0
        fprintf(stderr, "called_color_long[%d] = %d\n", i, called_color_long[i]);
#endif
    }
#if 0
    fprintf(stderr, "Deleted bases= ");
    for (i=ibind; i<ibind+indsize_bases; i++)
    {
        fprintf(stderr, "%c", color2base[called_color_long[i]]);    
    }     
    fprintf(stderr, "\n");
#endif
    

#if 0
    {
       int j=0, k=0;
       for  (i=0; i<len; i++)
       {
           if (called_color_long[i] < 0)
               j++;
           if (called_color_short[i] < 0)
               k++;
       }
       fprintf(stderr, "Num of uncalled bases after the 3rd pass: long = %d short = %d\n",
         j, k);
    }
#endif

    if (0)
    {
        output_undetermined_bases_in_long_read(indsize_bases, len, color1, color2,
            pos1, pos2, called_color_long);
        output_undetermined_bases_in_short_read(indsize_bases, len,
            color1, color2, pos1, pos2, called_color_short);
        compare_long_and_short_reads(indsize_bases, len, color1, color2, pos1, pos2,
            called_color_short, called_color_long);
    }

    /* Determine the sequence of bases */
    for (i=ibind; i<len; i++)
    {
        if (i<ibind+indsize_bases)
        {
            if (called_color_long[i] >= 0)
                long_called_seq[i] = color2base[called_color_long[i]];
        }
        else
        {
            if ((called_color_long[i] >= 0) &&
                (called_color_long[i] == called_color_short[i-indsize_bases]))
            {
                long_called_seq[i] = color2base[called_color_long[i]];
                short_called_seq[i-indsize_bases] = 
                    color2base[called_color_short[i-indsize_bases]];
            }
        }
    }
#if 0
    fprintf(stderr, "Final long seq=\n%s\n", long_called_seq+ibind);
    fprintf(stderr, "Final short seq=\n%s\n", short_called_seq+ibind);
#endif
#if 0
    fprintf(stderr, "Base1=%s length=%d\n",  base, strlen(base));
#endif
    FREE(base);
    FREE(color1);
    FREE(color2);
    FREE(pind1);
    FREE(pind2);
    FREE(pos1);
    FREE(pos2);
    FREE(highest1_signal);
    FREE(highest2_signal);
    FREE(called_color_long);
    FREE(called_color_short);

    return SUCCESS;            
}

/*******************************************************************************
 * Function: get_base_index_by_position
 * Purpose: Find the index of the rightmost called base (= called peak)
 *          located to the left from the given position.
 *          If there is no such base, return -1
 *          If the given position equals the position of some base
 *          (= called peak), return the index of this base (called peak)
 *******************************************************************************
 */
int
get_base_index_by_position(int given_position, int *coordinate, int len)
{
    int i, base_index = -1;                          
    int loc, prev_loc, next_loc;
   
//  fprintf(stderr, "given_pos=%d len=%d first_coord=%d last_coord=%d\n",
//         given_position, len, coordinate[0], coordinate[len-1]); 

    if (len < 1) {
        base_index = -1;
        return base_index;
    }
    else if (len == 1) {
        base_index = 0;
        return base_index;
    }
    else if (len == 2) {
        if (given_position <= (coordinate[0] + coordinate[1])/2)
           base_index = 0;
        else 
           base_index = 1;
        return base_index;
    }

    if (given_position < (coordinate[0] + coordinate[1])/2)
    {
           base_index = 0;
           return base_index;
    }
    
    if (given_position > (coordinate[len-2] + coordinate[len-1])/2)
    {
           base_index = len-1;
           return base_index;
    }

    for (i=1; i< len-1; i++)
    {
        loc = coordinate[i];
        prev_loc = coordinate[i-1];
        next_loc = coordinate[i+1];
        if ((given_position > (loc + prev_loc)/2) &&
            (given_position <=(loc + next_loc)/2))
        {
            base_index = i;
            return base_index;
        }
    }

    return base_index;
}

/*******************************************************************************
 * Function: compute_tpars
 *           Without option -raw processing occurs as usual
 *           With -raw, raw data will be analyzed and then, depending on other
 *                output options, bases may or may not be called
 *           Main steps:
 *           1. Call mixed bases
 *           2. Detect indel location and size
 *           3. Output indel location and size
 *           4. Call long and short sequences 
 *           5. Split a mixed chromatogram
 *******************************************************************************
 */
int
Btk_process_indels(char *path, int *num_datapoints, int **chromatogram, 
    char *color2base, uint8_t *quality_values, Data *data, 
    ReadInfo *read_info, ContextTable *ctable, Options options, 
    BtkMessage *msg)
{
    int     i, indcount = 0;
    int    *indloc = NULL, *indsize_scans = NULL, *indbind = NULL, color;
    int    *num_ident_bases = NULL;
    double *params[NUM_PARAMS] = {NULL, NULL, NULL, NULL};
    char   *prefix_name = NULL, *scf_file_name = NULL;
    char   *seq_name = NULL; 
    char   *long_called_seq  = NULL, 
           *short_called_seq = NULL; 
    char    amplicon_name[MAXPATHLEN], donor_name[MAXPATHLEN];
    char    path1[MAXPATHLEN], path2[MAXPATHLEN];
    int     num_tokens = 0;
    char  **tokens;

    indloc          = CALLOC(int, MAX_NUM_INDELS);
    indbind         = CALLOC(int, MAX_NUM_INDELS);
    indsize_scans   = CALLOC(int, MAX_NUM_INDELS);
    num_ident_bases = CALLOC(int, MAX_NUM_INDELS);

    num_ident_bases[0] = 1;

    if ((seq_name = strrchr(options.file_name, '/')) != NULL) {
        seq_name++;
    }
    else {
        seq_name = options.file_name;
    }

    if (options.Verbose > 1)
    {
        /* Need this because subsequent processing will spoil the path */
        strcpy(path1, path);
        strcpy(path2, path);

        /* Count the number of tokens separated by '/' */
        if (strtok(path1, "/") != NULL)
        {
            num_tokens++;
            while (strtok(NULL, "/") != NULL)
            {
                num_tokens++;
            }

            if (num_tokens >= 4)
            {
                tokens = CALLOC(char *, num_tokens);
                for (i=0; i<num_tokens; i++)
                    tokens[i] = CALLOC(char, MAXPATHLEN);    

                strcpy(tokens[0], strtok(path2, "/"));

                for (i=1; i<num_tokens; i++)
                {
                    strcpy(tokens[i], strtok(NULL, "/"));
                }
                strcpy(amplicon_name, tokens[num_tokens-3]);
                strcpy(donor_name, tokens[num_tokens-4]);
                for (i=0; i<num_tokens; i++)
                    FREE(tokens[i]);
                FREE(tokens);
            }
            else
            {
                strcpy(amplicon_name, "?");
                strcpy(donor_name, "?");
            }
        }
        else
        {
            strcpy(amplicon_name, "?");
            strcpy(donor_name, "?"); 
        } 
    }

    /* 1. Call mixed bases */
//  call_mixed_bases(color2base, data, options);

    /* 2. Detect indel location and size */
    if ((options.indloc < 0) && (options.indsize <= 0))
    {
        (void)detect_indels_from_chromatogram(&indcount, num_ident_bases,
            &indloc, &indbind, &indsize_scans, data, quality_values,
            options);

//      fprintf(stderr, "Indloc = %d\n", indloc[0]);
    }
    else if ((options.indloc < 0) && (options.indsize > 0))
    {
        (void)detect_indels_from_chromatogram(&indcount, num_ident_bases,
            &indloc, &indbind, &indsize_scans, data, quality_values,
            options);
        
        indsize_scans[0] =
            get_autocorrelation_indel_size(indloc[0]+6, options.indsize,
            data->length, options.file_name, data, options);

    }
    else if ((options.indloc >= 0) && (options.indsize <= 0))
    {
        indloc[0] = options.indloc;
        indsize_scans[0] =
            get_autocorrelation_indel_size(indloc[0]+6, -1,
            data->length, options.file_name, data, options);

        indbind[0] = get_base_index_by_position(indloc[0]+6, 
            data->bases.coordinate, data->bases.length);
        if (indsize_scans[0] > 0) {
            indcount++;
        }
    }
    else
    {
        indcount++; 
        indloc[0]        = options.indloc;
        indsize_scans[0] = 
            get_autocorrelation_indel_size(indloc[0]+6, options.indsize,
            data->length, options.file_name, data, options);

        indbind[0] = get_base_index_by_position(indloc[0]+6, 
            data->bases.coordinate, data->bases.length);
    } /* detect indel location and size */
#if 0
    fprintf(stderr, "ibind= %d\n", indbind[0]);
#endif

//  fprintf(stderr, "options.indsize=%d indcount=%d indsize_scans = %d num_ident_bases = %d indbind=%d \n", 
//      options.indsize, indcount, indsize_scans[0], num_ident_bases[0], indbind[0]);

    /* 3. Output detected indel location and size */
    if (indcount == 0)
        output_indels(seq_name, num_ident_bases, indbind,
        indcount, indsize_scans, indloc, data->bases.bases, data->bases.length,
        donor_name, amplicon_name, quality_values, options);

    if ((indbind[0] < 0) || (indsize_scans[0] < 0))
    {
        FREE(indloc);
        FREE(indbind);
        FREE(indsize_scans);
        FREE(num_ident_bases);
        return SUCCESS;
    }
#if 0
    fprintf(stderr, "Original called mixed seq=\n");
    for (i=0; i<data->bases.length; i++)
        fprintf(stderr, "%c", data->bases.bases[i]);
    fprintf(stderr, "\n");
#endif

    /* 4. Call long and short sequences */
    long_called_seq  = CALLOC(char, data->bases.length+1);
    short_called_seq = CALLOC(char, data->bases.length+1); 
    call_long_and_short_sequences(indbind[0], indsize_scans[0], 
        color2base, long_called_seq, short_called_seq, data, 
        options);

#if 0
    fprintf(stderr, "Long seq=\n%s\n", long_called_seq+indbind[0]);
    fprintf(stderr, "Short_seq=\n%s\n", short_called_seq+indbind[0]);
#endif

    output_indels(seq_name, num_ident_bases, indbind, indcount, indsize_scans, 
        indloc, long_called_seq, strlen(long_called_seq),  
        donor_name, amplicon_name, quality_values, options);

#if 0
    fprintf(stderr, "Indbind = %d\n", indbind);
#endif
    
    if (options.indel_resolve)
    {
        int   prefix_len, indsize_bases = ROUND((float)indsize_scans[0]/12.);
        char *suffix = NULL;
        int  *long_chromatogram[NUM_COLORS] = {NULL, NULL, NULL, NULL};
        int  *short_chromatogram[NUM_COLORS] = {NULL, NULL, NULL, NULL};
        int  *long_locs=NULL, long_bases_len=0, long_data_len=0; 
  
        prefix_name   = CALLOC(char, MAXPATHLEN);
        scf_file_name = CALLOC(char, MAXPATHLEN);

        if ((suffix = strchr(options.file_name, '.')) != NULL)
            suffix++;
        prefix_len = strlen(options.file_name) - strlen(suffix);
        strncpy(prefix_name, options.file_name, prefix_len-1);
        prefix_name[prefix_len] = '\0';

#if 0
        fprintf(stderr, "File_name=%s Prefix=%s Suffix = %s\n", 
            options.file_name, prefix_name, suffix);
#endif
        if ((strncmp(suffix, "_long" , 5)  == 0) ||
            (strncmp(suffix, "_short", 6) == 0))
        {
            fprintf(stderr, "Chromatograms are separated already\n");
            fprintf(stderr, "No SCF file will be output\n");
            FREE(prefix_name);
            FREE(scf_file_name);
            FREE(indloc);
            FREE(indbind);
            FREE(indsize_scans);
            FREE(num_ident_bases);
            FREE(long_called_seq);
            FREE(short_called_seq);
            return SUCCESS;
        }

        for (color = 0; color < NUM_COLORS; color++)
        {
            ColorData *cd = &data->color_data[color];
            long_chromatogram[color] = CALLOC(int, data->color_data[0].length);
            short_chromatogram[color] = CALLOC(int, data->color_data[0].length);
            for (i=0; i<data->color_data[0].length; i++)
            {
                long_chromatogram[color][i] = cd->data[i];
                short_chromatogram[color][i] = cd->data[i];
            }
        }

        if (produce_new_chromatogram(LONG, indsize_scans[0], indbind[0], 
            color2base, long_chromatogram, data, long_called_seq, 
            options) != SUCCESS)
        {
            fprintf(stderr, "Error extracting long chromatogram\n");
            goto error;
        }

        for (color=0; color<NUM_COLORS; color++)
        {
            for (i=0; i<data->length; i++)
            {
                long_chromatogram[color][i] /= 2;
            }
        }

        /* Evaluate quality values for long chromatogram */ 
        for (i = 0; i < NUM_PARAMS; i++)
        {
            params[i] = CALLOC(double, data->bases.length);
        }

        for (i=0; i<data->bases.length; i++) 
        {
            quality_values[i] = 0;
        }

        /* Store "long" data */
        long_bases_len = data->bases.length;
        long_data_len  = data->length;
        long_locs  = CALLOC(int,  long_bases_len);
        for (i=0; i<long_bases_len; i++)
        {
            long_locs[i]  = data->bases.coordinate[i];
        }

        /* Process "short" chromatogram */
        if (produce_new_chromatogram(SHORT, indsize_scans[0], indbind[0], 
            color2base, short_chromatogram, data, short_called_seq,
            options) != SUCCESS)
        {
            fprintf(stderr, "Error extracting short chromatogram\n");
            goto error;
        }

        for (color=0; color<NUM_COLORS; color++)
        {
            for (i=0; i<data->length; i++)
            {
                short_chromatogram[color][i] /= 2;
            }
        }

        /* Adjust the long and short chromarogram to make sure
         * their sum exactly equals the original chromatogram
         */
        for (color=0; color<NUM_COLORS; color++)
        {
            for (i=0; i<data->length; i++)
            {
                int long_chromat  = long_chromatogram[color][i];
                int short_chromat = short_chromatogram[color][i];    
                int sum_chromats  = long_chromat + short_chromat;
                int residual = (chromatogram[color][i] - sum_chromats)/2;
                
                if ((long_chromat  + residual >= 0) &&
                    (short_chromat + residual >= 0))
                {
                    long_chromatogram[color][i] = long_chromat  + residual;
                    short_chromatogram[color][i]= short_chromat + residual;
                }
                else if ((long_chromat  + residual <  0) &&
                         (short_chromat + residual >= 0))
                {
                    long_chromatogram[color][i] = 0;
                    short_chromatogram[color][i]= short_chromat + residual +
                        (long_chromat  + residual);
                }
                else if ((long_chromat  + residual >= 0) &&
                         (short_chromat + residual <  0))
                {
                    short_chromatogram[color][i]= 0;
                    long_chromatogram[color][i] = short_chromat + residual +
                        (long_chromat  + residual);
                }
#if 0        
                residual = (chromatogram[color][i] - short_chromatogram[color][i]    
                                                   - long_chromatogram[color][i] );
                if (residual != 0)
                    fprintf(stderr, "Residual[[color=%d][scan=%d] = %d\n",
                      color, i, residual); 
#endif
#if 0
                fprintf(stderr, "color= %d pos= %d chromat= %d final_long= %d final_short=%d\n",
                    color, i, chromatogram[color][i], long_chromatogram[color][i],
                    short_chromatogram[color][i]);
#endif
            }
        }

        /* Output the results */
        sprintf(scf_file_name, "%s_long", prefix_name);
        if (output_scf_file(scf_file_name, options.scf_dir,
            long_called_seq, long_locs,
            quality_values, long_bases_len, long_data_len,
            long_chromatogram[0], long_chromatogram[1],
            long_chromatogram[2], long_chromatogram[3],
            color2base, options.chemistry) == ERROR)
        {
            fprintf(stderr, "Error producing long SCF file\n");
            goto error;
        }

        prefix_name[prefix_len-1] = '\0';
        sprintf(scf_file_name, "%s_short", prefix_name);
        if (output_scf_file(scf_file_name, options.scf_dir,
            short_called_seq, data->bases.coordinate,
            quality_values, data->bases.length - indsize_bases, 
            data->length - indsize_scans[0],
            short_chromatogram[0], short_chromatogram[1],      
            short_chromatogram[2], short_chromatogram[3],      
            color2base, options.chemistry) == ERROR)
        {
            fprintf(stderr, "Error producing short SCF file\n");
            goto error;
        }

        if (options.xgr)
        {
            output_chromatogram("1_Deconvolved_data_long.xgr",
            "Long deconvolved chromatogram",
            long_chromatogram[0], long_chromatogram[1],
            long_chromatogram[2], long_chromatogram[3],        
            long_data_len, data);

            output_chromatogram("1_Deconvolved_data_short.xgr",
            "Short deconvolved chromatogram",
            data->color_data[0].data, data->color_data[1].data,
            data->color_data[2].data, data->color_data[3].data,
            data->length, data);
       
            /* Output sum of separated chromatograms */ 
            for (color=0; color<NUM_COLORS; color++)
            {
                for (i=0; i<data->length; i++)
                {
                    long_chromatogram[color][i] += 
                        short_chromatogram[color][i];    
                }
            }

            output_chromatogram("1_Deconvolved_data_sum.xgr",
            "Sum of deconvolved chromatograms",
            long_chromatogram[0], long_chromatogram[1], 
            long_chromatogram[2], long_chromatogram[3], long_data_len, data);
        }

        for (i=0; i<4; i++)
        {
            FREE(long_chromatogram[i]);
            FREE(short_chromatogram[i]);
        }
        for (i = 0; i < NUM_PARAMS; i++) {
            FREE(params[i]);
        }
        FREE(long_locs);
        FREE(prefix_name);
        FREE(scf_file_name);
        FREE(indloc);
        FREE(indbind);
        FREE(indsize_scans);
        FREE(num_ident_bases);
        FREE(long_called_seq);
        FREE(short_called_seq);
    }
    else
    { 
        FREE(indloc);
        FREE(indbind);
        FREE(indsize_scans);
        FREE(num_ident_bases);
        FREE(long_called_seq);
        FREE(short_called_seq);
    }
    
    return SUCCESS;

error:
    for (i = 0; i < NUM_PARAMS; i++) 
    {
        FREE(params[i]);
    }
    FREE(prefix_name);
    FREE(scf_file_name);    
    FREE(long_called_seq);
    FREE(short_called_seq);
    return ERROR;
}
