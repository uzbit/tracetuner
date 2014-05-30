/**************************************************************************
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
 * Copyright (c) 1999-2003 Paracel, Inc.  All rights reserved.
 */

/*
 *  Btk_call_bases.c $Revision: 1.29 $
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

#include "Btk_qv.h"
#include "SFF_Toolkit.h"
#include "util.h"
#include "nr.h"
#include "Btk_qv_data.h"
#include "Btk_qv_funs.h"
#include "train.h"              /* needs util.h */
#include "Btk_process_peaks.h"  /* needs train.h */
#include "context_table.h"
#include "tracepoly.h"  
#include "Btk_lookup_table.h"
#include "Btk_call_bases.h"  
#include "Btk_compute_tpars.h"  /* needs train.h */
#include "Btk_compute_tp.h"
#include "Btk_get_mixed_bases.h"
#include "Btk_process_raw_data.h"

#define BASE_IND1                     75
#define BASE_IND2                    150
#define BC5_CALL_PEAKS                 0
#define BC5_MERGE_PEAKS                0
#define BC5_SPLIT_PEAKS                0
#define BC5_UNCALL_PEAKS               0
#define BC5_CALLED_NON_DIP             0
#define BC5_UNCALLED_DIP               1
#define BIF                            1.1   /* used in BC1-4 and BC5.3 */
#define CHECK_AFTER_LOOP_BC4           0
#define CHECK_BASE_INDEX               0
#define CHECK_CALLED_PEAKS             0
#define CHECK_CD_PEAK_IND              0
#define CHECK_CD_PEAK_LISTS            0
#define CHECK_DATA_BASE_INDEX          0
#define CHECK_NUM_CALLED_PEAKS         0
#define CHECK_PEAK_POSITION            0
#define CHECK_REORDERING               0
#define COMPLETE_MOB_SHIFTS            1  
#define DEBUG_CURTIS                   0 
#define DELETE_BASES                   1
#define DIP_FRACTION_FINAL             0.9 
#define DIP_FRACTION_INITIAL           0.7 
#define FRAC_SHORT_CALLED_PEAK         0.8
#define INSERT_BASES                   1
#define MAX_DELTA_SPACING              4
#define MAX_INSERTION_FACTOR           1.1
#define MAX_PEAK_WIDTH               100
#define MAX_POS_DYE_BLOBS            340
#define MAX_SHIFT_MAX                  2.0
#define MAX_SHIFT_INC                  0.5
#define MAX_SIZE_OF_SEARCH_REGION    100
#define MIN_DISTANCE_BETWEEN_BASES     2
#define MIN_DYE_BLOB_RESOLUTION       (0.3*ONE_PLUS)
#define MIN_RELATIVE_CHANNEL_WEIGHT    0.12
#define MIN_PEAK_SPACING               3
#define MIN_PEAK_SPLIT_RESOLUTION     (0.2*ONE_PLUS)
#define MIN_SIZE_OF_GOOD_REGION        6
#define MOB_POLYFIT_DEG                2
#define N_SPAC_MODEL_COEFF             3
#define NUM_PEAKS_PER_WINDOW          50
#define PEAK_SPLIT_FACTOR             (1.125*ONE_PLUS)
#define LOW_BOUNDARY_FRACTION         (0.9*ONE_MINUS)
#define LOW_BOUNDARY_FRACTION_SNP     (0.76*ONE_MINUS)
#define OUTPUT_SPACING_CURVE           0
#define POS_RAISE_INSERTION_RATE     400
#define PRINT_NORM_FOR_PLOTTING        0 /* flag for debug output of normalization */ 
#define RATIO_SHIFT_MAX_REL_SPC_ERR    1.5
#define REL_SPC_ERR_THRESH             0.25
#define REORDER_BASES                  0
#define SHIFT_MAX_THRES                2.0
#define SHIFT_ERR_GROW_FACTOR          1.2
#define SHOW_CASE                      0
#define SHOW_CHANGED_BASES             0
#define SHOW_PEAK_LIST                 0
#define SHOW_REVIEW                    0
#define SIGNAL_STRENGTH_FRACTION       0.65
#define SPACING_FROM_GOOD_REGION       1
#define SQRT_ENVELOPE_ASYMPTOTE        0.5  
#define TRUNCATED_HEIGHT            1600
#define USE_IS_DYE_BLOB_2001           1
#define USE_INSTINSIC_SIGNAL           0
#define WORST_BASE_NOT_DIP             0
#define WORST_BASE_NOT_TYPE11          0

/* crude estimate of peak spacing */
static const double	gSpacEst0   = 12.0;

extern double
mobShiftQuality(Peak * const peak[], int n_peaks);

extern int
collect_some_stats( Data *data, Options *options, BtkMessage *message,
                    Results *results );

static void
check_cd_peak_ind(Data *data)
{
    int i,j;
    fprintf(stderr, "CHECK_CD_PEAK_IND:\n");
    for (i=0; i<NUM_COLORS; i++)
    {
        ColorData *cd = &data->color_data[i];
        for (j=0; j<cd->peak_list_len; j++) 
        {
            if (j != cd->peak_list[j].cd_peak_ind)
                fprintf(stderr,
                "Error: peak(color=%d, ind=%d, base_ind=%d).cd_peak_ind =%d\n",
                i,j, cd->peak_list[j].base_index, cd->peak_list[j].cd_peak_ind);
        }
    }
    fprintf(stderr, "...done\n");
}

static void
check_cd_peak_lists(Data *data)
{
    int i,j;
    fprintf(stderr, "CHECK_CD_PEAK_LISTS:\n");
    for (i=0; i<NUM_COLORS; i++)
    {
        ColorData *cd = &data->color_data[i];
        for (j=1; j<cd->peak_list_len; j++)
        {
            double      ipos = cd->peak_list[j  ].ipos;
            double prev_ipos = cd->peak_list[j-1].ipos;
            Peak  *pl = cd->peak_list;

            if (ipos < prev_ipos)
            {
                fprintf(stderr,
                "Error: peak(color=%d, ind=%d, is_called=%d bind=%d).ipos = %.1f ", 
                i,j, cd->peak_list[j  ].is_called, 
                cd->peak_list[j  ].base_index, ipos);
                fprintf(stderr,
                "< peak(color=%d, ind=%d, is_called=%d bind=%d).ipos == %.1f\n", 
                i,j-1, cd->peak_list[j-1].is_called, 
                cd->peak_list[j-1].base_index, prev_ipos);
            }

            if ((pl[j].is_called == 0 && pl[j].base_index != -1) ||
                (pl[j].is_called != 0 && pl[j].base_index == -1))
            {
                fprintf(stderr, "Error: peak[%d](color=%d, is_called=%d).base_index = %d\n",
                j,i,pl[j].is_called,pl[j].base_index); 
            } 

            if (pl[j].end - pl[j].beg > MAX_PEAK_WIDTH)
                fprintf(stderr, "Error: too wide peak (color=%d ind=%d) beg=%d end=%d width = %d\n",
                    i, j, pl[j].beg, pl[j].end, pl[j].end - pl[j].beg);
        }
    }
    fprintf(stderr, "...done\n");
}

void
check_base_index(Data *data)
{
    int i = 0;
    Peak **cpl = data->bases.called_peak_list;

    fprintf(stderr, "CHECK_BASE_INDEX:\n");
    for (i=0; i<data->bases.length; i++) {
        int color, cd_peak_ind;

        if (data->bases.called_peak_list[i] == NULL)
            continue;

        color       = cpl[i]->color_index;
        cd_peak_ind = cpl[i]->cd_peak_ind;
        
        if (cpl[i]->is_called == 0)
            fprintf(stderr,
            "Error: called_peak[%d] points to uncalled color=%d ind=%d)\n",
            i, color, cd_peak_ind);

        if (cpl[i]->is_called > 0 && i != cpl[i]->base_index)
            fprintf(stderr,
            "Error: called_peak[%d]->base_index = %d (is_called= %d)\n",
            i, cpl[i]->base_index, cpl[i]->is_called);

        if (i != data->color_data[color].peak_list[cd_peak_ind].base_index)
        {
            fprintf(stderr,
            "Error: base_index==%4d != color_data[%d].peak[%3d].base_index=%d ",
            i, color, cd_peak_ind,
            data->color_data[color].peak_list[cd_peak_ind].base_index);
            fprintf(stderr, "(is_called = %d)\n", 
                data->color_data[color].peak_list[cd_peak_ind].is_called);
        }
    }
    fprintf(stderr, "...done\n");
}

void
check_data_base_index(Data *data)
{
    int i, base_index = 0;
    Peak **dpl = data->peak_list;
    fprintf(stderr, "CHECK_DATA_BASE_INDEX:\n");
    for (i=0; i<data->peak_list_len; i++) {
        if ((dpl[i]->is_called <= 0)     ||
            (dpl[i]->base_index < 0))
            continue;

        if (base_index != dpl[i]->base_index)
        {
            fprintf(stderr,
            "Error: base_index==%4d != data_base_index=%4d ",
                base_index,  dpl[i]->base_index);
            fprintf(stderr,
            "for data_peak %d base=%c is_called=%d ipos=%.1f\n",
                i, dpl[i]->base, dpl[i]->is_called, dpl[i]->ipos);
        }
        base_index++;
    }
    fprintf(stderr, "...done\n");
}

static int  
check_reordering(Data *data)
{
    int  i, j;
    Peak **cpl = data->bases.called_peak_list;
    char *bases = data->bases.bases;

    fprintf(stderr, "CHECK_REORDERING:\n");
    for (i=1; i<data->bases.length; i++) 
    {
         
        if (cpl[i  ] == NULL)
            continue;

        if (cpl[i]->base != bases[i])
            fprintf(stderr,
              "Error: called %d peak[%d]->base=%c != bases[%d]=%c, is_called=%d\n",
                cpl[i]->is_called, i, cpl[i]->base, i, bases[i], cpl[i]->is_called);
       
        if (cpl[i-1] == NULL)
            continue;
 
        if (cpl[i]->ipos < cpl[i-1]->ipos)
        {
            fprintf(stderr, 
            "Error: called %d peak_%c[%d]->ipos=%.1f < called %d peak_%c[%d]->ipos=%.1f\n",
            cpl[i]->is_called, cpl[i]->base, i, cpl[i]->ipos, 
            cpl[i-1]->is_called, cpl[i-1]->base, i-1, cpl[i-1]->ipos);
        }
    }

    for (j=0; j<NUM_COLORS; j++) {
        int base_index = -1;
        int peak_ind = -1;
        ColorData *cd = &data->color_data[j];

        for (i=0; i<data->color_data[j].peak_list_len; i++) {
            if (cd->peak_list[i].is_called <= 0) {
                continue;
            }
            if (base_index>=0 &&
                base_index >= cd->peak_list[i].base_index)
            {
                fprintf(stderr, "Error: peak[%d][%d].base_index=%d ", 
                    j, peak_ind, base_index);
                fprintf(stderr, ">= peak[%d][%d].base_index=%d\n",
                    j, i, data->color_data[j].peak_list[i].base_index);
            }
            base_index = cd->peak_list[i].base_index;
            peak_ind = i;
        }
    }
    fprintf(stderr, "...done\n");
    return SUCCESS;
}

static void 
check_called_peaks(Data *data)
{
    int i,j;
    char color2base[4]={'A','C','G','T'};
    char *bases = data->bases.bases;
    Peak **cpl = data->bases.called_peak_list;
    int  *coord = data->bases.coordinate;

    fprintf(stderr, "CHECK_CALLED_PEAKS:\n");
    for (i=0; i<data->bases.length; i++) 
    {
        int lbound = (i == 0) ? 0 : (coord[i] + coord[i-1])/2;
        int rbound = (i == data->bases.length-1) ?
            data->length-1 : (coord[i]+coord[i+1])/2;

        if (cpl[i] == NULL) continue;

        if (cpl[i]->ipos < lbound)
             fprintf(stderr,
                "Error: position %.1f of called %d peak %d < lbound = %d\n",
                cpl[i]->ipos, cpl[i]->is_called, i, lbound);

        if (cpl[i]->ipos >= rbound)
             fprintf(stderr,
                "Error: position %.f of called %d peak %d >= rbound = %d\n",
                cpl[i]->ipos, cpl[i]->is_called, i, rbound);

        if (i == 0) continue;
 
        for (j=1; j<data->bases.length; j++) 
        {
            if (j<=i) continue;

            if (cpl[j] != NULL && cpl[i] == cpl[j])
            { 
                fprintf(stderr,
                "Error: called_peaks %d and %d point to the same real peak\n",
                i, j);
            }
        }
        if (bases[i] != color2base[cpl[i]->color_index])
            fprintf(stderr,
                "Error: called base[%d]=%c != called_peak_base= %c\n",
                 i, bases[i], color2base[cpl[i]->color_index]);
    }
    fprintf(stderr, "...done\n");
}

static void
check_num_called_peaks(Data *data)
{
    int i, j, base_index;
    Peak **cpl = data->bases.called_peak_list;
    int num_called_peaks = data->bases.length;

    fprintf(stderr, "CHECK_NUM_CALLED_PEAKS:\n");
    for (j=0; j<NUM_COLORS; j++) 
    {
        for (i=0; i<data->color_data[j].peak_list_len; i++) 
        {
            if (data->color_data[j].peak_list[i].is_called > 0) 
            {
                base_index=data->color_data[j].peak_list[i].base_index;
                if (base_index < 0 || cpl[base_index] == NULL ||
                    cpl[base_index] != &data->color_data[j].peak_list[i])
                {
                    fprintf(stderr, 
                    "Base index %d of called peak %d (color=%d, is_called=%d) is incorrect\n",
                     base_index, i, j, data->color_data[j].peak_list[i].is_called);
                }
                if (i>0 && data->color_data[j].peak_list[i].base_index==
                    data->color_data[j].peak_list[i-1].base_index) {
                    fprintf(stderr,
                       "Peaks %d and %d of color %d have the same base index %d\n",
                      i-1, i, j, base_index);
                }
                if (i<data->color_data[j].peak_list_len-1 &&
                    data->color_data[j].peak_list[i  ].base_index==
                    data->color_data[j].peak_list[i+1].base_index) {
                    fprintf(stderr,
                       "Peaks %d and %d of color %d have the same base index %d\n",
                       i+1, i, j, base_index);
                }
            }
        }
    }

    for (i=0; i<data->bases.length; i++) {
        if (cpl[i] == NULL)
            num_called_peaks--;

        if (cpl[i]!= NULL &&
            (cpl[i]->is_called > 0)) {
            if (i>0 && data->bases.called_peak_list[i] != NULL &&
               cpl[i-1]!=NULL && cpl[i]->ipos <= cpl[i-1]->ipos) {
            }
            if (i<data->bases.length-1 &&
                cpl[i]!= NULL && cpl[i+1] != NULL &&
               cpl[i]->ipos >= cpl[i+1]->ipos) {
            }
        }
        else {
           fprintf(stderr, "Base %d remained uncalled", i);
           if (data->bases.called_peak_list[i]!= NULL)
               fprintf(stderr, ": is_called=%d\n",
               data->bases.called_peak_list[i]->is_called);
           else
               fprintf(stderr, "\n");
        }
    }
    if (num_called_peaks < data->bases.length)
        fprintf(stderr, "Warning: num_called_peaks=%d < num_bases=%d\n",
            num_called_peaks, data->bases.length);
    fprintf(stderr, "...done\n");
}

static double
get_spacing_from_good_region(int base_index, Data *data);

static void
output_spacing_curve(Data *data, Options *options)
{
    int i, l=0, base_index = 1, base_spacing;
    double sum_devs=0., curr_spacing1, curr_spacing2, curr_spacing;
    FILE *excelout, *xgraphout;

    excelout = fopen("Spacing_curve.excel", "w");
    fprintf(excelout,
    "base_ind   bpos   prev_bpos  scan  curr_spa  ");
    fprintf(excelout,
    "base_spa  curr_spa1 curr_spa2    delta\n");
    for (i=0; i<data->color_data[0].length; i++) {
        if ((base_index < data->bases.length-1) &&
            (i == data->bases.called_peak_list[base_index]->ipos))
        {
            base_index++;
            base_spacing =
                data->bases.called_peak_list[base_index  ]->ipos-
                data->bases.called_peak_list[base_index-1]->ipos;
        }

        if ((i == 0) || (i % 20) != 0)
            continue;

        curr_spacing1 = spacing_curve(i);
        curr_spacing2 = get_spacing_from_good_region(base_index, data);
        curr_spacing = (curr_spacing1 + curr_spacing2)/2.;

        l++;
        sum_devs += fabs(curr_spacing1-(double)base_spacing);

        fprintf(excelout,
        "%5d     %.1f %.1f  %5d    %5.3f  %5d       %5.3f    %5.3f    %6.3f\n",
            base_index,
            data->bases.called_peak_list[base_index-1]->ipos,
            data->bases.called_peak_list[base_index  ]->ipos,
            i,  curr_spacing,
            base_spacing, curr_spacing1, curr_spacing2,
            ((double)base_spacing-curr_spacing)/curr_spacing);
    }
    fclose(excelout);
    fprintf(stderr, "Average dev=%5.3f\n", (l>0) ? sum_devs/(double)l : sum_devs);
    xgraphout = fopen("Spacing_curve.xgraph", "w");
    fprintf(xgraphout, "\nTitleText: Spacing curves for %s\n",
        options->file_name);
    fprintf(xgraphout, "\n\"Curr_spacing \n");
    for (i=0; i<data->color_data[0].length; i++) {
        if ((base_index < data->bases.length-1) &&
            (i == data->bases.called_peak_list[base_index]->ipos))
        {
            base_index++;
            base_spacing =
                data->bases.called_peak_list[base_index  ]->ipos-
                data->bases.called_peak_list[base_index-1]->ipos;
        }

        if ((i == 0) || (i % 20) != 0)
            continue;

        curr_spacing1 = spacing_curve(i);
        curr_spacing2 = get_spacing_from_good_region(base_index, data);
        curr_spacing = (curr_spacing1 + curr_spacing2)/2.;

        fprintf(xgraphout, "%f   %f\n", (double)i, curr_spacing);
    }
    fprintf(xgraphout, "\n\"Base_spacing \n");
    base_index=1;
    for (i=0; i<data->color_data[0].length; i++) {
        if ((base_index < data->bases.length-1) &&
            (i == data->bases.called_peak_list[base_index]->ipos))
        {
            base_index++;
            base_spacing =
                data->bases.called_peak_list[base_index  ]->ipos-
                data->bases.called_peak_list[base_index-1]->ipos;
        }

        if ((i == 0) || (i % 100) != 0)
            continue;

        fprintf(xgraphout, "%f   %f\n", (double)i, (double)base_spacing);
    }

    fclose(xgraphout);
}

/*******************************************************************************
 * Function: colordata_find_peak_index_by_location
 * Purpose: In the list of peaks of given color, find the peak such that a 
 *          given position is located in its area. If the position does not 
 *          fall into any peak's area, return peak_index equal to the negative 
 *          index of the peak right to the given position. If there're no peaks 
 *          to the right from the given position or the position is located 
 *          left from the peak with index 0, return peak_index == NINF
 *******************************************************************************
 */
int
colordata_find_peak_index_by_location(ColorData *color_data, int given_position,
    Peak *peak, int *peak_index, BtkMessage *message)
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
 * Function: data_find_peak_index_by_position
 * Purpose:  in the list of pointers to all peaks, find the index of the peak
 *           (1) with position closest to the given position and
 *           (2) with peak_index between lo and hi.
 *           If no such peak is found, return -1
 *           If there are two ore more peaks equally distanced from the given
 *               position, return the smallest of their indexes
 *           Strategy:
 *           1) find a couple of peaks with indexes ind1 and ind2=ind1+1 and
 *              positions pos1 and pos2 such that pos1 <= given_position <= pos2
 *           2) determine which of the pos1 and pos2 is closer to the 
 *              given_position
 *           3) if pos2 is closer, return ind2; otherwise, scan peak list
 *              left of ind1 to find minimal index of peak with position pos1
 *              and return this index
 *******************************************************************************
 */
int
data_find_peak_index_by_position(Data *data, int given_position,
    int lo_ind, int hi_ind, int *peak_index, BtkMessage *message)
{
    int mid;

   *peak_index = -1;  
    if (data->peak_list_len < 1 ||
        given_position < data->peak_list[lo_ind]->ipos ||
        given_position > data->peak_list[hi_ind]->ipos) 
    {
        (void)sprintf(message->text, "Given position is out of bounds\n"); 
        return ERROR;        
    }

    if (given_position == data->peak_list[lo_ind]->ipos) {
       *peak_index = lo_ind;
        return SUCCESS;
    }
    if (given_position == data->peak_list[hi_ind]->ipos) {
        while (hi_ind > lo_ind && 
               data->peak_list[hi_ind  ]->ipos ==
               data->peak_list[hi_ind-1]->ipos)
        {
            hi_ind--;   
        }
       *peak_index = hi_ind;
        return SUCCESS;
    }
    do {
        mid = (lo_ind + hi_ind) / 2;
        if (given_position ==  data->peak_list[mid]->ipos) {
            while (mid > lo_ind && 
                data->peak_list[mid  ]->ipos ==
                data->peak_list[mid-1]->ipos)
            {
                mid--;
            }
           *peak_index = mid;
            break;
            return SUCCESS;
        }
        else if (given_position < data->peak_list[mid]->ipos) {
            hi_ind = mid;
        }
        else if (given_position > data->peak_list[mid]->ipos) {
            lo_ind = mid;
        }
    } while (lo_ind < hi_ind-1);
    
   *peak_index = ( given_position - data->peak_list[lo_ind]->ipos) <=
                 (-given_position + data->peak_list[hi_ind]->ipos) ? 
                 lo_ind : hi_ind;
    if (*peak_index==lo_ind) {
        while (lo_ind > 0 && data->peak_list[lo_ind  ]->ipos ==
                         data->peak_list[lo_ind-1]->ipos)
            {
                lo_ind--;
            }
           *peak_index = lo_ind; 
    }
    return SUCCESS;
}

/*******************************************************************************
 * Function: data_review_peak_list
 * Purpose: count and print out the number of
 *          - uncalled true peaks;
 *          - called noise peaks;
 *          - total true peaks, etc.
 *******************************************************************************
 */
#if SHOW_REVIEW
static void
data_review_peak_list(Data *data, BtkMessage *message)
{
    int  i, true_called=0, noise_called=0;
    int  true_uncalled=0, noise_uncalled=0, true=0, noise=0;
    for (i=0; i<data->peak_list_len; i++) {
        if (is_true_peak(*data->peak_list[i])) {
            true++;
            if (data->peak_list[i]->is_called > 0) {
                true_called++;
            }
            else {
                true_uncalled++;
            }
        }
        else {
            noise++;
            if (data->peak_list[i]->is_called > 0) {
                noise_called++;
            }
            else {
                noise_uncalled++;
            }
        }
    }
    fprintf(stderr, "Total peaks = %d\nTrue peaks=%d\nTrue called peaks=%d\n",
        data->peak_list_len, true, true_called);
    fprintf(stderr, "True uncalled peaks = %d\nNoise peaks=%d\nNoise called peaks=%d\n",
        true_uncalled, noise, noise_called);
    fprintf(stderr, "Noise uncalled peaks=%d\n",noise_uncalled);
}
#endif


/*******************************************************************************
 * Function: reorder_called_bases_and_peaks
 * Purpose: reorder bases and called_peak_list arrays so that
 *          the called peak position would not decrease with base number
 *******************************************************************************
 */
int
bc_reorder_called_bases_and_peaks(Data *data, BtkMessage *message)
{
    int     i, ordered, start_pos=0;
    char    tb;
    Peak   *tp;
    Peak  **cpl = data->bases.called_peak_list;
    char   *bases = data->bases.bases;

    fprintf(stderr, "REORDER_BASES\n"); 
    ordered=0;
    while (!ordered) 
    {
        ordered = 1;
        for (i=start_pos+1; i<data->bases.length; i++) 
        {
            if (cpl[i-1]->ipos <= cpl[i]->ipos)
            { 
                start_pos = i-1;
                continue;
            }

            /* Swap called_peaks and bases i and i-1 */
            ordered = 0;
            
            tb         = bases[i-1];
            bases[i-1] = bases[i  ];
            bases[i  ] = tb;

            tp       = cpl[i-1];
            cpl[i-1] = cpl[i  ];
            cpl[i  ] = tp;
            cpl[i-1]->base_index = i-1;
            cpl[i  ]->base_index = i  ;

            break;

        } /* end loop in i */
    } /* while (!ordered) */

    if (CHECK_REORDERING) check_reordering(data);

    return SUCCESS;
}

/*******************************************************************************
 * Function: bc_data_create_single_ordered_peak_list
 * Purpose: create an array of pointers to the peaks of any color,
 *          ordered by peak position
 *******************************************************************************
 */
int
bc_data_create_single_ordered_peak_list(Data *data, int *shift, 
    BtkMessage *message)
{
    int    i, color;
    int    best_color, best_base_index;
    double best_ipos, best_iheight, best_area;
    int    curr_index[NUM_COLORS] = {0, 0, 0, 0};
    /* current index in the list of peaks of a given color
     * that is put into a single peak list
     */
    FREE(data->peak_list);
    data->peak_list_len = 0;

#if 1
    for (i=0; i<NUM_COLORS; i++) {
        shift[i]=0;
    }
#endif
    i = data->color_data[0].peak_list_len +
        data->color_data[1].peak_list_len +
        data->color_data[2].peak_list_len +
        data->color_data[3].peak_list_len;

    data->peak_list = CALLOC(Peak *, i);

    /* Order peaks with respect to their position and put them into the list */
    i = 0;
    while ((curr_index[0] < data->color_data[0].peak_list_len) ||
           (curr_index[1] < data->color_data[1].peak_list_len) ||
           (curr_index[2] < data->color_data[2].peak_list_len) ||
           (curr_index[3] < data->color_data[3].peak_list_len))
    {
        /* Find the leftmost peak  */
        best_ipos       = INF;
        best_base_index = INF;
        best_iheight    = -1.;
        best_area       = -1;
        best_color      = -1;

        for (color=0; color < NUM_COLORS; color++)
        {
            double new_ipos = shift[color] +
                data->color_data[color].peak_list[curr_index[color]].ipos;
            int new_base_index =
                data->color_data[color].peak_list[curr_index[color]].base_index;
            double new_iheight =
                data->color_data[color].peak_list[curr_index[color]].iheight;
            double new_area =
                data->color_data[color].peak_list[curr_index[color]].area;
            
            if (curr_index[color] >= data->color_data[color].peak_list_len) {
                continue;
            }

            /* Peak with left position goes first */
            if (DBL_GT_DBL(best_ipos, new_ipos) ||

            /* Peak with lower base index goes first */
               (DBL_EQ_DBL(new_ipos, best_ipos) &&
                new_base_index < best_base_index)  ||

            /* Higher peak goes first */
               (DBL_EQ_DBL(new_ipos, best_ipos) &&
                new_base_index == best_base_index  &&
                DBL_GT_DBL(new_iheight, best_iheight)) ||                

            /* Peak with bigger area goes first */
               (DBL_EQ_DBL(new_ipos, best_ipos) &&
                new_base_index == best_base_index  &&
                DBL_EQ_DBL(new_iheight, best_iheight) &&
                DBL_GT_DBL(new_area, best_area)))
            {
                best_color      = color;
                best_ipos       = new_ipos;
                best_base_index = new_base_index;
                best_iheight    = new_iheight;
                best_area       = new_area;
            }
        } /* loo through all colors */

        if (best_color >= 0) 
        {
            Peak **dpl = data->peak_list;
            ColorData *cd = &data->color_data[best_color];
            dpl[i] = &cd->peak_list[curr_index[best_color]];
            dpl[i]->data_peak_ind  = i;
            dpl[i]->data_peak_ind2 = i;
            dpl[i]->cd_peak_ind    = curr_index[best_color];
            if (SHOW_PEAK_LIST)
            {
                fprintf(stderr, "data_peak[%d] color=%d cd_peak_ind=%d ipos=%.1f ",
                    i, dpl[i]->color_index, curr_index[best_color], dpl[i]->ipos);
                fprintf(stderr, "base_index=%d iheight=%.1f area=%.1f\n",
                    dpl[i]->base_index, dpl[i]->iheight, dpl[i]->area);
            }
            curr_index[best_color]++;
            i++;
        }
        else {
            fprintf(stderr, "Best color not found for i=%d:\n", i);
            fprintf(stderr, "   peak list length: %d %d %d %d\n",
                data->color_data[0].peak_list_len,
                data->color_data[1].peak_list_len,
                data->color_data[2].peak_list_len,
                data->color_data[3].peak_list_len);
            fprintf(stderr, "   curr_index[color] = %d %d %d %d\n",
                curr_index[0], curr_index[1], curr_index[2], curr_index[3]);
            fprintf(stderr, "   shift = %d %d %d %d \n\n",
                shift[0], shift[1], shift[2], shift[3]);
            break;
            break;
        }
    }
    data->peak_list_len = i;

    if (CHECK_REORDERING) check_reordering(data);

    return SUCCESS;
}


/*******************************************************************************
 * Function: call_peak       
 *******************************************************************************
 */
static int
call_peak(Data *data, int base_index, int data_peak_ind, int Case,
    BtkMessage *message)
{
    int j;

    if (Case <= 0) {
        fprintf(stderr, "Base called in wrong case %d\n", Case);
        return ERROR;
    }

    data->peak_list[data_peak_ind]->is_called  = Case;
    data->peak_list[data_peak_ind]->base_index = base_index;

    data->bases.length++;
    if (base_index > data->bases.length-1) {
        fprintf(stderr,
            "Inserting at base_index=%d > data->bases.length-1=%d\n",
        base_index, data->bases.length-1);
        return ERROR;
    }

    if (data->bases.length >= data->bases.max_length) {
        data->bases.max_length *= 2;
        data->bases.called_peak_list =
            REALLOC(data->bases.called_peak_list, Peak *,
            data->bases.max_length);
        MEM_ERROR(data->bases.called_peak_list);
        data->bases.bases = REALLOC(data->bases.bases, char,
            data->bases.max_length);
        MEM_ERROR(data->bases.bases);
        data->bases.coordinate = REALLOC(data->bases.coordinate, int,
            data->bases.max_length);
        MEM_ERROR(data->bases.coordinate);
    }

    /* Move bases, coordinates and called peaks one
     * spot to free space for the new element
     */
    if (base_index < data->bases.length-1) {
        (void)memmove(&data->bases.bases[base_index+1],
            &data->bases.bases[base_index],
            (data->bases.length-base_index-1) * sizeof(char));
        (void)memmove(&data->bases.called_peak_list[base_index+1],
            &data->bases.called_peak_list[base_index],
            (data->bases.length-base_index-1) * sizeof(Peak*));
        (void)memmove(&data->bases.coordinate[base_index+1],
            &data->bases.coordinate[base_index],
            (data->bases.length-base_index-1) * sizeof(int));
    }

    /* Increment base indexes of subsequent called peaks */
    for (j=base_index+1; j<data->bases.length; j++) {
        if (data->bases.called_peak_list[j] != NULL) {
            data->bases.called_peak_list[j]->base_index=j;
        }
        else {
            fprintf(stderr,
                "called_peak_list[%d] == NULL \n", j);
        }
    }

    /* Insert the new element */
    data->bases.bases[base_index] = data->peak_list[data_peak_ind]->base;
    data->bases.coordinate[base_index] = data->peak_list[data_peak_ind]->ipos;
    data->bases.called_peak_list[base_index] =
        data->peak_list[data_peak_ind];

    return SUCCESS;

error:
    return ERROR;
}


/*******************************************************************************
 * Function: adjust_called_base_position
 * Purpose:  make sure position of the called base does not coincide
 *           with position of previous or next called base
 *******************************************************************************
 */
static void
adjust_called_base_position(int base_index, Data *data)
{
    int ipos  = data->bases.called_peak_list[base_index]->ipos;

    if ((base_index > 0) &&
        (data->bases.called_peak_list[base_index-1] != NULL))
    {
        int iposm = (base_index > 0) ?
            data->bases.called_peak_list[base_index-1]->ipos : 0;
        if (ipos == iposm) data->bases.called_peak_list[base_index]->ipos++;
    }

    if ((base_index < data->bases.length-1) &&
        (data->bases.called_peak_list[base_index+1] != NULL)) 
    {
        int iposp = (base_index < data->bases.length-1) ?
            data->bases.called_peak_list[base_index+1]->ipos :
            data->bases.called_peak_list[data->bases.length-1]->ipos+1;

        if (ipos == iposp) data->bases.called_peak_list[base_index]->ipos--;
    }
}

/*******************************************************************************
 * Function: insert_and_resolve_peaks
 *
 * Purpose: 1) insert the peaks created upon splitting into peak lists
 *          2) reallocate memory
 *          3) resolve peaks
 *******************************************************************************
 */
int
insert_and_resolve_peaks(Data *data, int jc, int base_ind1, int base_ind2,
    int cd_peak_ind, int data_peak_ind, Peak peak1, Peak peak2, int Case, 
    Options *options, BtkMessage *message)
{
    int    j, base_index, il, ir;
    double max_resolution;

    peak1.base_index    = base_ind1;
    peak2.base_index    = base_ind2;

    if (peak1.pos < peak1.beg || peak1.pos > peak1.end) {
        (void)sprintf(message->text,
        "Position of peak1 (#%d) out of bounds (Case BC%d)\n",
            cd_peak_ind, Case);
        fprintf(stderr," Beg=%d Pos=%d End=%d\n", peak1.beg, peak1.pos, peak1.end);
        return ERROR;
    }

    if (peak2.pos < peak2.beg || peak2.pos > peak2.end) {
        (void)sprintf(message->text,
        "Position of peak2 (#%d) out of bounds (Case BC%d)\n",
            cd_peak_ind+1, Case);
        return ERROR;
    }
 
    /* Insert peak1 into peak list */
    peak1.cd_peak_ind = cd_peak_ind;
    data->color_data[jc].peak_list[cd_peak_ind] = peak1;
    data->bases.called_peak_list[base_ind1] =
        &data->color_data[jc].peak_list[cd_peak_ind];
    data->bases.called_peak_list[base_ind1]->base_index = base_ind1;
    data->bases.bases[base_ind1] = peak1.base;
        /* need this for truncated peaks */ 

    /* Reallocate memory for new peak, if needed */
    data->color_data[jc].peak_list_len++;
    if (data->color_data[jc].peak_list_len >
        data->color_data[jc].peak_list_max_len) {
        data->color_data[jc].peak_list_max_len *=2;
        data->color_data[jc].peak_list =
            REALLOC(data->color_data[jc].peak_list, Peak,
            data->color_data[jc].peak_list_max_len);
        MEM_ERROR(data->color_data[jc].peak_list);
    }

    /* Shift any peaks and bases to the right one spot to make room */
    cd_peak_ind++;
    (void)memmove(&data->color_data[jc].peak_list[cd_peak_ind+1],
        &data->color_data[jc].peak_list[cd_peak_ind],
        (data->color_data[jc].peak_list_len-1 - cd_peak_ind)
         * sizeof(Peak));

    /* Insert peak2 into peak list */
    peak2.cd_peak_ind = cd_peak_ind;
    peak2.base_index = base_ind2;
    data->color_data[jc].peak_list[cd_peak_ind]=peak2;

    /* Update cd_peak_ind of shifted peaks. 
     * NOTE: for now, the second peak remains uncalled  
     */
    for (j=cd_peak_ind+1; j<data->color_data[jc].peak_list_len; j++) 
        data->color_data[jc].peak_list[j].cd_peak_ind = j;

    /* Updata base_index and called peak list */
    if (data_peak_ind < 0) {

        /* We are in loop 1-4, so no new base will be inserted */
        data->bases.bases[base_ind2] = peak2.base;
        for (j=0; j<data->color_data[jc].peak_list_len; j++) {
            if (data->color_data[jc].peak_list[j].is_called > 0) {
                base_index = data->color_data[jc].peak_list[j].base_index;
                data->bases.called_peak_list[base_index] =
                    &data->color_data[jc].peak_list[j];
            }
        }
    }
    else {

        /* We are in loop 5; new base will be inserted 
         * For now, leave the 2nd peak uncalled and only
         * reset called_peak_list to peaks that were already called; 
         * later will call function call_peak like in case BC5.3
         */
        data->color_data[jc].peak_list[cd_peak_ind].is_called = 0;
        data->color_data[jc].peak_list[cd_peak_ind].base_index = -1;
        for (j=0; j<data->color_data[jc].peak_list_len; j++) {
            if (data->color_data[jc].peak_list[j].is_called > 0) {
                base_index = data->color_data[jc].peak_list[j].base_index;
                data->bases.called_peak_list[base_index] =
                    &data->color_data[jc].peak_list[j];
            }
        }
    }

    /* Update data peak list */
    if (data_peak_ind >= 0) {

        /* If data_peak_ind is defined (>=0), -> we are in 5th loop 
         * need to reallocate arrays, because we are adding new called peak
         */
        data->peak_list_len++;
        if (data->peak_list_len >= data->peak_list_max_len) {
            data->peak_list_max_len *= 2;
            data->peak_list = REALLOC(data->peak_list, Peak *,
                data->peak_list_max_len);
            MEM_ERROR(data->peak_list);
        }
       
        /* Members of data peak_list with indexes < data_peak_ind
         * will stay on their places and will
         * point to the same peaks as before, except those
         * of color jc, which should be "refreshed"
         * after reallocation of memory
         */
        for (j=0; j<cd_peak_ind; j++) {
            int k= data->color_data[jc].peak_list[j].data_peak_ind;
            data->peak_list[k] = &data->color_data[jc].peak_list[j];
        }

        /* Members with indexes >= data_peak_ind+1
         * should be shifted one spot to the right. 
         * Before the shift, those of these members that point to peaks 
         * of color jc should be refreshed
         */
        for (j=cd_peak_ind+1; j<data->color_data[jc].peak_list_len; j++) {
            int k= data->color_data[jc].peak_list[j].data_peak_ind;
            data->peak_list[k] = &data->color_data[jc].peak_list[j];
        }
        (void)memmove(&data->peak_list[data_peak_ind+2],
                      &data->peak_list[data_peak_ind+1],
        (data->peak_list_len-2 - data_peak_ind) * sizeof(Peak*));
        for (j=data_peak_ind+2; j<data->peak_list_len; j++) 
            data->peak_list[j]->data_peak_ind = 
                data->peak_list[j]->data_peak_ind2 = j;

        /* Finally, member data_peak_ind+1 should be set to 
         * point to the inserted peak
         */
        data->peak_list[data_peak_ind+1] =
            &data->color_data[jc].peak_list[cd_peak_ind];
        data->peak_list[data_peak_ind+1]->data_peak_ind =
             data->peak_list[data_peak_ind+1]->data_peak_ind2= data_peak_ind+1;

        /* Finally, call the inserted uncalled peak */
        call_peak(data, base_ind2, data_peak_ind+1, Case, message);
    }

    /* Resolve peaks created upon splitting */
    il = cd_peak_ind;
    while (il > 0 &&
        data->color_data[jc].peak_list[il].type > 20) {
        il--;
    }
    ir = cd_peak_ind;
    while (ir < data->color_data[jc].peak_list_len-1 &&
        data->color_data[jc].peak_list[ir].type%10 > 1) {
        ir++;
    }

    if (resolve_multiple_peaks(data, jc, il, ir-il+1,
        &max_resolution, options, message) != SUCCESS)
    {
        (void)fprintf(stderr,
            "Error calling resolve_multiple_peaks in BC%d\n", Case);
        return ERROR;
    }

    /* Make sure the positions of called peaks after resolving do not 
     * coincide with positions of adjacent called peaks 
     */
    for (j=il; j<=ir; j++) 
    {
        int base_index =  data->color_data[jc].peak_list[j].base_index;

        if (data->color_data[jc].peak_list[j].is_called == 0)
            continue;

        adjust_called_base_position(base_index, data);
    }

    return SUCCESS;

    error:
    return ERROR;
}


/*******************************************************************************
 * Function: uncall_peak        
 ******************************************************************************
 */
int
uncall_peak(int i, Data *data, BtkMessage *message)
{
    int j;

#if 0
    fprintf(stderr, "Deleting base %c at pos=%d\n",
        data->bases.bases[i], data->bases.called_peak_list[i]->ipos);
#endif

    if (data->bases.called_peak_list[i] != NULL) {
        data->bases.called_peak_list[i]->base_index = -1;
        data->bases.called_peak_list[i]->is_called  = 0;
    }
    if (i < data->bases.length-1) {

        /* Shift any called peaks and bases to the left one spot */
        (void)memmove(&data->bases.bases[i], &data->bases.bases[i+1],
            (data->bases.length-1-i) * sizeof(char));
        (void)memmove(&data->bases.called_peak_list[i],
            &data->bases.called_peak_list[i+1],
            (data->bases.length-1-i) * sizeof(Peak *));
        (void)memmove(&data->bases.coordinate[i],
            &data->bases.coordinate[i+1],
            (data->bases.length-1-i) * sizeof(int));
    }
    data->bases.length--;

    for (j=i; j<data->bases.length; j++)
        if (data->bases.called_peak_list[j] != NULL)
            data->bases.called_peak_list[j]->base_index = j;

    return SUCCESS;
}


/*******************************************************************************
 * Function: peak_with_low_boundaries
 *
 * MDC - 08/01/03 - This function returns a boolean indication of whether
 * the input peak has boundaries (beginning and ending heights) below some
 * magic value.  My change is only to add this explanation (no code changes).
 *******************************************************************************
 */
static int
peak_with_low_boundaries(Peak *peak, ColorData *cd, int het)
{
    double hbeg   = cd->data[peak->beg],
           hend   = cd->data[peak->end],
           max_height = cd->data[peak->max];

#if 0
    fprintf(stderr, "base=%c hbeg=%f hend=%f max_height*FRACTION=%f\n",
        peak->base, hbeg, hend, max_height*LOW_BOUNDARY_FRACTION_SNP);
#endif

    if (!het &&
        (hbeg < max_height * LOW_BOUNDARY_FRACTION) &&
        (hend < max_height * LOW_BOUNDARY_FRACTION))
        return 1;

    if (het &&
        (hbeg < max_height * LOW_BOUNDARY_FRACTION_SNP) &&
        (hend < max_height * LOW_BOUNDARY_FRACTION_SNP))
        return 1;

    return 0;
}

/*******************************************************************************
 * Function: is_poorly_resolved_peak
 *******************************************************************************
 */
static int
is_poorly_resolved_peak(Peak *peak, ColorData *cd, int het)
{
    if (!het &&
        ((peak->type == 11) || (peak->type == 22) || 
         (peak->type == 21) || (peak->type == 12))) 
        return 0;

    if ( peak_with_low_boundaries(peak, cd, het) &&
        (peak->resolution < MIN_DYE_BLOB_RESOLUTION ))
    return 0;

    return 1;
}

/*******************************************************************************
 * Function: get_dye_blob_multiplicity
 *******************************************************************************
 */
static int
get_dye_blob_multiplicity(Peak *peak, ColorData *cd, int het,
    int *num_poorly_resolved)
{
    int i = peak->cd_peak_ind, m = 1;

   *num_poorly_resolved = 0;

    if (peak->type == 11)
        return 1;

    /* Find the index of the first peak */
    if (peak->type > 20) {
        while ((i>0) && (cd->peak_list[i].type>20)) {
            i--;
        }
    }

    /* Now count all and the poorly resolved peaks */
    i++;   /* start from second peak */
    while ((i<cd->peak_list_len-1) && (cd->peak_list[i].type>20)) {
       m++;
       if (is_poorly_resolved_peak(&cd->peak_list[i], cd, het))
           (*num_poorly_resolved)++;
       i++;
    }

    return m;
}

/*******************************************************************************
 * Function: get_modified_height
 *******************************************************************************
 */
static float
get_modified_height(int pos, Peak *peak, ColorData *cd)
{
    float hbeg = cd->data[peak->beg];
    float hend = cd->data[peak->end];
    float modif_height;

    if ((pos < peak->beg) || (pos > peak->end))
        return 0.;

    modif_height = cd->data[pos] - 
        (hbeg*(peak->end - pos) - hend*(pos - peak->beg))
        /(peak->end - peak->beg);    

    if (modif_height < 0)
        modif_height = 0;

    return modif_height;
}

#if USE_IS_DYE_BLOB_2001
/*******************************************************************************
 * Function: is_dye_blob
 *******************************************************************************
 */
int
is_dye_blob(int pos, Peak *peak, Data *data, int het)
{
    int peak_ind =  peak->cd_peak_ind;
    ColorData *cd = &data->color_data[peak->color_index];

    if ((peak->beg < 0) || (peak->beg >= cd->length))
        return 1;
    if ((peak->end < 0) || (peak->end >= cd->length))
        return 1;

    /* Two or more truncated peaks in a row is a dye blob */
    if (peak->is_truncated &&
        (((peak->type    > 30) && (peak_ind > 0) &&
         (cd->peak_list[peak_ind-1].is_truncated)) ||
         ((peak->type%10 == 3) && (peak_ind < cd->peak_list_len-1) &&
          cd->peak_list[peak_ind+1].is_truncated)))
        return 1;

    if (!het &&
        ((peak->ipos > MAX_POS_DYE_BLOBS*DEFAULT_PEAK_SPACING) ||
          !is_poorly_resolved_peak(peak, cd, het) ||
         (peak_with_low_boundaries(peak, cd, het) &&
          (peak->resolution < MIN_DYE_BLOB_RESOLUTION))))
        return 0;
    else if (het &&
        ( !is_poorly_resolved_peak(peak, cd, het) ||
         (peak_with_low_boundaries(peak, cd, het) &&
          (peak->resolution < MIN_DYE_BLOB_RESOLUTION))))
        return 0;
    else {
        int p, m = get_dye_blob_multiplicity(peak, cd, het, &p);

        if ((m > 3) && (p > (m/2)))
//      if ((m > 3) && (p > 1)))
            return 1;

        if (het &&
            (((m >= 3) && (p > (m/2))) ||
             ((m >= 6) && (p >=(m/2))))
           )
            return 1;
    }

    return 0;
}

#else
/*******************************************************************************
 * Function: is_dye_blob
 *******************************************************************************
 */
int
is_dye_blob(int pos, Peak *peak, Data *data, int het)
{
    int i;
    int color        = peak->color_index;
    int modif_height; 
    int other_signal = 0, other_peak_ind=-1;
    ColorData *cd    = &data->color_data[color];

    if ((peak->beg < 0) || (peak->beg >= cd->length))
        return 0;

    if ((peak->end < 0) || (peak->end >= cd->length))
        return 0;

    if (pos > MAX_POS_DYE_BLOBS*DEFAULT_PEAK_SPACING)
        return 0;

    if ((peak->type == 11) && !peak->is_truncated)
        return 0;

    /* If all of the following apply, peak is dye blob: 
     * - current peak is truncated,
     * - its closest neibours of the same color on both sides are truncated
     * - there is a nonzero signal of other color undeneath of current peak
     */
    if (peak->is_truncated)
        modif_height = (float)cd->data[pos] * DYE_BLOB_FRACTION; 
    else
        modif_height = get_modified_height(pos, peak, cd);

    /* First, explore the possibility that there exists other peak 
     * at position pos
     */
    {
        Peak other_peak;
        BtkMessage message;
        other_signal = 0;
        for (i=0; i<NUM_COLORS; i++) 
        {
            int other_modif_height = 0;
            other_peak.pos    = 0;
            other_peak.height = 0;
            other_peak_ind = -1;

            if (i == color)
                continue;

            colordata_find_peak_index_by_location(&(data->color_data[i]), 
                pos, &other_peak, &other_peak_ind, &message);
                
            if ((other_peak_ind >= 0)  &&                           
                (other_peak.pos > pos - DEFAULT_PEAK_SPACING/2) &&
                (other_peak.pos < pos + DEFAULT_PEAK_SPACING/2)) 
            {
                if (other_peak.is_truncated)
                {
                    other_modif_height = data->color_data[i].data[pos]
                         * DYE_BLOB_FRACTION;
                }
                else 
                {
                    other_modif_height = 
                        get_modified_height(pos, &other_peak, 
                            &data->color_data[i]);
                }
            }
            if (other_signal < other_modif_height)
                other_signal = other_modif_height;
        }
        if (modif_height < other_signal)
            return 1;
    } 
       
    return 0;
}
#endif

/*******************************************************************************
 * Function: is_dip
 * Purpose:  check if a given peak is a dominant intrinsic peak
 *******************************************************************************
 */
int
is_dip(Peak **peak_list, int peak_list_len, int i, int *shift, Options *options,
    Data *data)
{
    int j, il, ir, base_index;
    double dip_fraction=DIP_FRACTION_INITIAL, other_signal, wiheight;
    int    color = peak_list[i]->color_index;

    wiheight = peak_list[i]->iheight;
        if (is_dye_blob(peak_list[i]->ipos, peak_list[i], data, 0))
                wiheight *= DYE_BLOB_FRACTION;
    /* Determine the range of competing peaks */
    il = i;
    while (il > 0) {
        int color1 = peak_list[il]->color_index;
        il--;
        if (peak_list[il]->ipos + peak_list[il]->orig_width/2. +
             2.*peak_list[il]->beta + shift[color1]   <
            peak_list[i ]->ipos + shift[color])
            break;
    }
    ir = i;
    while (ir < peak_list_len-1) {
        int color1 = peak_list[ir]->color_index;
        ir++;
        if (peak_list[ir]->ipos - peak_list[ir]->orig_width/2. -
            2.*peak_list[ir]->beta + shift[color1] >
            peak_list[i ]->ipos + shift[color])
        break;
    }

    /* Determine if the current peak is a DIP */
    base_index = 0;
    for (j = il; j <= ir; j++) {
        int color1 = peak_list[j]->color_index;

        if (peak_list[j]->is_called > 0)
            base_index = peak_list[j]->base_index;

        if (j == i)
            continue;

        other_signal = Shape(peak_list[j]->C0, peak_list[j]->orig_width/2.,
            peak_list[j]->beta,
            (double)(peak_list[i]->ipos + shift[color]
                   - peak_list[j]->ipos - shift[color1]),
            options);

        if (is_dye_blob(peak_list[j]->ipos, peak_list[j], data, 0))
            other_signal *= DYE_BLOB_FRACTION;


        /* NOTE:
         * the conditional below is not exact if a context table is
         * used, because this conditional uses wiheights computed
         * using cweights from the previous loop. These cweights,
         * however, might change due to insertion of new bases
         */
        if (base_index > BASE_IND2)
            dip_fraction = DIP_FRACTION_FINAL;
        else if ((base_index > BASE_IND1) && (base_index <=BASE_IND2)) {
            dip_fraction = DIP_FRACTION_INITIAL +
                (DIP_FRACTION_FINAL - DIP_FRACTION_INITIAL) /
                (double)(BASE_IND2 - BASE_IND1)   *
                (double)(base_index- BASE_IND1);
        }

        /* if (wiheight + EPS < dip_fraction * other_signal)  */
        if( DBL_LT_DBL( wiheight, dip_fraction * other_signal) )
            return 0;
    }

    return 1;
}

/*******************************************************************************
 * Function: get_truncated_width
 *******************************************************************************
 */
static int
get_truncated_width(Peak *peak, Data *data) 
{
    if (!peak->is_truncated)
		return 0;

	else {
        int i, truncated_width=-1, color = peak->color_index;
		ColorData *cd = &(data->color_data[color]);

        for (i = peak->beg; i <= peak->end; i++) {
            if (cd->data[i] == TRUNCATED_HEIGHT) 
				truncated_width++;
		}
        return truncated_width;
	}
}

/******************************************************************************
 * Function: get_representative_peak_signal
 ******************************************************************************
 */
static int
get_representative_peak_signal(int position, Peak *peak, Data *data, 
    Options *options)
{
    int        signal, color = peak->color_index;
    ColorData *cd = &data->color_data[color];

    if (position >= 0) {
        if (!peak->is_truncated) {
            signal = INT_DBL(Shape(peak->C0, peak->orig_width/2., peak->beta,
                     (double)(peak->ipos - position), options));
        }
        else  {
            signal = cd->data[position];
        }
    }
    else {
        signal = INT_DBL(peak->iheight);
    }

    return signal;
}

/*******************************************************************************
 * Function: is_better_peak
 * Purpose:  determine if the first peak is to be preferred to the second one
 *           when recalling N to the best guess.  Return 1 if yes, 0 if no.
 *
 * MDC - 08/01/03
 *
 * This function was VERY dependent on the color ordering, as it had many
 * places where it didn't check if the peaks had the same values; i.e. it
 * equated "equal to" with "less than" or "greater than", meaning the
 * function was not commutative (an issue when the color arrays are provided
 * in different orders).  I have changed it to break ties as follows, with
 * the first matching rule winning:
 * a.  Best peak with the former criteria, sans equality.
 * b.  Peak that is higher to the left, keep following equality until one
 *     wins or the end of the array is hit
 * c.  Peak that is higher to the right, keep following equality until one
 *     wins or the end of the array is hit
 *           <if we get here, we have 2 identical traces)
 * d.  First color in this list that appears:  ACGT.  Note that this order
 *     allows a simple character magnitude comparison.
 *
 *******************************************************************************
 */
static int
is_better_peak(Peak *peak1, Peak *peak2, int position, Data *data, 
    Options *options)
{
    int signal1, is_dye_blob1=0, lb1, width1;
    int signal2, is_dye_blob2=0, lb2, width2;
    int *d1;
    int *d2;
    int i;
    ColorData *cd1 = &data->color_data[peak1->color_index];
    ColorData *cd2 = &data->color_data[peak2->color_index];
	
    /* Set the peaks' signal value at the position.
     * If peak is not truncated, use its intrinsic signal;
     * otherwise, use apparent sinal
     */
    signal1 = get_representative_peak_signal(position, peak1, data, options);
    signal2 = get_representative_peak_signal(position, peak2, data, options);

    is_dye_blob1 = is_dye_blob(position, peak1, data, 0);
    is_dye_blob2 = is_dye_blob(position, peak2, data, 0);

    /* Compare the two peaks */
    if (!is_dye_blob1 && !is_dye_blob2) {
        if ((signal1 < TRUNCATED_HEIGHT) || (signal2 < TRUNCATED_HEIGHT)) {
            if (INT_GT_DBL(signal1, signal2 / SIGNAL_STRENGTH_FRACTION)) {
                return 1;
	    }
            if (INT_LT_DBL(signal1, signal2 * SIGNAL_STRENGTH_FRACTION)) {
                return 0;
	    }
	    lb1 = peak_with_low_boundaries(peak1, cd1,  0);
	    lb2 = peak_with_low_boundaries(peak2, cd2,  0);
	    if (lb1 && !lb2) {
		return 1;
	    }
	    if (!lb1 && lb2) {
		return 0;
	    }
	    if (signal1 > signal2) {
		return 1;
	    }
	    if (signal1 < signal2) {
		return 0;
	    }
	    /* fall-through to tie-breakers at end */
        }
        else { /* Both the peaks are truncated */
            width1 = get_truncated_width(peak1, data);
            width2 = get_truncated_width(peak2, data);
            if (width1 < width2) {
		return 1;
	    }
            if (width1 > width2) {
		return 0;
	    }
	    /* fall-through to tie-breakers at end */
        }
    }

    else if (!is_dye_blob1 && is_dye_blob2) {
        if (peak2->is_truncated) {
            if (INT_GT_DBL(peak1->iheight, signal2 * DYE_BLOB_FRACTION)) {
		return 1;
	    }
            if (INT_LT_DBL(peak1->iheight, signal2 * DYE_BLOB_FRACTION)) {
		return 0;
	    }
	    /* fall-through to tie-breakers at end */
	}
        else {
            if (INT_GT_DBL(signal1, signal2 * DYE_BLOB_FRACTION)) {
		return 1;
	    }
            if (INT_LT_DBL(signal1, signal2 * DYE_BLOB_FRACTION)) {
		return 0;
	    }
	    /* fall-through to tie-breakers at end */
	}
    }

    else if (is_dye_blob1 && !is_dye_blob2) {
        if (peak1->is_truncated) {
            if (DBL_GT_INT(signal1 * DYE_BLOB_FRACTION, peak2->iheight)) {
		return 1;
	    }
            if (DBL_LT_INT(signal1 * DYE_BLOB_FRACTION, peak2->iheight)) {
		return 0;
	    }
	}
        else {
            if (DBL_GT_INT(signal1 * DYE_BLOB_FRACTION, signal2)) {
		return 1;
	    }
            if (DBL_LT_INT(signal1 * DYE_BLOB_FRACTION, signal2)) {
		return 0;
	    }
	    /* fall-through to tie-breakers at end */
	}
    }
	
    else {	/* Both the peaks are dye blobs */
	if (signal1 > signal2) {
	    return 1;
	}
	if (signal1 < signal2) {
	    return 0;
	}
	width1 = get_truncated_width(peak1, data);
	width2 = get_truncated_width(peak2, data);
	if (width1 < width2) {
	    return 1;
	}
	if (width1 > width2) {
	    return 0;
	}
	/* fall-through to tie-breakers at end */
    }

    /*
     * Break the tie.
     */
    d1 = data->color_data[peak1->color_index].data;
    d2 = data->color_data[peak2->color_index].data;

    /* b. */
    for (i = position - 1; i >= 0; i--) {
	if (d1[i] > d2[i]) {
	    return 1;
	}
	if (d1[i] < d2[i]) {
	    return 0;
	}
    }

    /* c. */
    for (i = position + 1; i < data->color_data[peak1->color_index].length; i++)
    {
	if (d1[i] > d2[i]) {
	    return 1;
	}
	if (d1[i] < d2[i]) {
	    return 0;
	}
    }

    /* d. */
    return (peak1->base < peak2->base);
}

/*****************************************************************************
 * Function: create_DIP_list
 * Purpose:  Create an array containing pointers to dominant intrinsic peaks
 *           (i.e., peaks that are not smaller than
 *               DIP_FRACTION * max of all other peaks at center of the DIPs)
 *           Allocate memory for the output array
 * Inputs:
 *           peak_list  array of pointers to all peaks sorted by position
 *           n_peak     length of peak_list pointer array
 *
 * Outputs:
 *           dip        array of pointers to DIPs
 *           num_DIPs      length of DIP pointer array
 * Return:   Error flag
 *****************************************************************************/
int
create_DIP_list(Data *data, int beg_pos, int end_pos, int *shift,
    Peak ***dip, int *num_DIPs, Options *options, BtkMessage *message)
{
    int         n_peak = data->peak_list_len;
    int         i, *idip = CALLOC(int, n_peak);
    Peak        pli;

    MEM_ERROR(idip);
   *num_DIPs = 0;

    /* Create a list of indexes of DIPs */
    for ( i=0; i < n_peak; i++ ) {
        pli = *data->peak_list[i];

            if ( pli.is_called > 0 ) {
                idip[(*num_DIPs)++] = i; /* peak_list[i] is a DIP */
                continue;
            }

        if (is_dip(data->peak_list, n_peak, i, shift, options, data))
            idip[(*num_DIPs)++] = i; /* peak_list[i] is a DIP */
    }

    /* Build list of DIPs */
   *dip = CALLOC(Peak *, *num_DIPs);
    MEM_ERROR(dip);
    for ( i=0; i < *num_DIPs; i++ )
        (*dip)[i] = data->peak_list[idip[i]];

    FREE(idip);
    return SUCCESS;

 error:
    FREE(idip);
    FREE(dip);
    return ERROR;
}


#if PRINT_NORM_FOR_PLOTTING
/*******************************************************************************
 * Function: print_for_plotting       
 *
 *******************************************************************************
 */
static void
print_for_plotting( ReadInfo* read_info, 
             double xx[4][MaxNumBases], double yy[4][MaxNumBases],
             int num_xy[4] )
{
    int i, c_indx;
    static char* fnames[] = 
    { "out/a.dat", "out/c.dat", "out/g.dat", "out/t.dat" };

    for( c_indx=0; c_indx<4; c_indx++ ) {
        FILE* fp = fopen( fnames[c_indx], "w" );
        if( !fp ) { 
            fprintf( stderr, "can't open file %s\n", fnames[c_indx] );
            break;
        }
        for( i=0; i<num_xy[c_indx]; i++ ) {
            double 
                x1 = xx[c_indx][i],
                y1 = yy[c_indx][i],
                factor = readNormFactor( read_info, c_indx, x1 ),
                y2 = y1*factor,
                y3 = traceEvalFunc( &(read_info->read_polys[c_indx]), x1 ),
                y4 = readEvalAveFunc( read_info, x1 );
            y2 = QVMIN( 1600, y2 );
#if 0
            fprintf( fp, 
         "%3.0f %3.0f %3.0f %3.0f %3.0f %5.1f %5.1f %5.1f %5.1f %5.1f\n",
                     x1, x1, x1, x1, x1,
                     y1, y2, y3, y4, 100*factor );
#else
            y2=y2; y4=y4;/* stop compiler whining */
            fprintf( fp, 
         "%3.0f %3.0f %8.5f %8.5f\n",
                     x1, x1,
                     y1, y3 );            
#endif
        }
        fclose( fp );
    }
}
#endif /* PRINT_NORM_FOR_PLOTTING */
/*******************************************************************************
 * Function: get_coefficients   
 *
 *******************************************************************************
 */
static void
get_coefficients( Data* data, ReadInfo* read_info )
{
    double xx[4][MaxNumBases], yy[4][MaxNumBases];
    int l, num_xy[4];

    for( l=0; l<4; l++ ) { num_xy[l] = 0; }

    for( l=0; l<data->bases.length; l++ ) {
        int c_indx, p_indx;
        Peak* peak = data->bases.called_peak_list[l];
        if( peak==NULL ) { continue; }
        c_indx = peak->color_index;
        p_indx = num_xy[c_indx];
        xx[c_indx][p_indx] = peak->ipos;
        yy[c_indx][p_indx] = peak->iheight;        
        num_xy[c_indx]++;
    }

    readGetCoefficients( read_info, xx, yy, num_xy, 0 );

#if PRINT_NORM_FOR_PLOTTING
        print_for_plotting( read_info, xx, yy, num_xy );
#endif

}
/******************************************************************************* 
 * Function: get_weighted_peak_heights
 * Purpose:  for each intrinsic peak, compute the product intrinsic peak height, 
 *           normalization weight and context weight
 *           
 *******************************************************************************
 */
static int
get_weighted_peak_heights(Data *data, char *color2base, ContextTable *ctable,
    ReadInfo *read_info, Options *options, BtkMessage *message)
{
    int i, color, base_index=-1;
    double  nweight, cweight;
    char   *context = NULL;

    if (ctable != NULL)
        context = CALLOC(char, ctable->dimension);

    for (i=0; i<data->peak_list_len; i++) {
        if (data->peak_list[i]->is_called > 0) {
            base_index = data->peak_list[i]->base_index; 

            /* Get context (will be used for uncalled peaks) */
            if ((ctable != NULL) &&
                (base_index >= ctable->dimension-1) &&
                (base_index < data->bases.length-1))
            {
                memcpy(context, &data->bases.bases[base_index - 
                    (ctable->dimension-1) + 1], ctable->dimension-1);
            }
        }
             
        /* Normalization and context weights */    
        nweight = 1.;
        cweight = 1.;
        if (options->renorm && !options->het) {   
            nweight = readNormFactor( read_info, data->peak_list[i]->color_index,
                        data->peak_list[i]->ipos ); 
        }
        if ((ctable != NULL) && 
            (base_index >= ctable->dimension-1) &&
            (base_index < data->bases.length-1)) 
        {
            if (data->peak_list[i]->is_called > 0) {
                cweight = weight_from_context(
                    &data->bases.bases[base_index], ctable);
            }
            else {
                color = data->peak_list[i]->color_index;
                context[ctable->dimension-1] = color2base[color];
                cweight = weight_from_context(&context[ctable->dimension-1],
                    ctable);
            }
        }
        data->peak_list[i]->wiheight = 
            data->peak_list[i]->iheight * nweight * cweight;
        if (is_dye_blob(data->peak_list[i]->ipos, data->peak_list[i], data, 0))
            data->peak_list[i]->wiheight *= DYE_BLOB_FRACTION;
    }

    if (ctable != NULL)
        FREE(context);
    
    return SUCCESS;
}

/*******************************************************************************
 * Function: show_case
 *******************************************************************************
 */
void 
show_case(Data *data, int i, int Case, int jc, int peak_ind)
{
    if (data->bases.called_peak_list[i] != NULL) {
        fprintf(stderr, 
        "i=%3d base=%c coord=%d ipos=%.1f Case=BC%d jc=%d peak_ind=%d peak_res=%6.3f\n",
        i, data->bases.bases[i], data->bases.coordinate[i],
        data->bases.called_peak_list[i]->ipos, Case,
        jc, peak_ind, data->bases.called_peak_list[i]->resolution);
    }
    else {
        fprintf(stderr, 
        "i=%3d base=%c coord=%d              Case=BC%d jc=%d peak_ind=%d \n",
        i, data->bases.bases[i], data->bases.coordinate[i],
        Case, jc, peak_ind);
    }
}

/*******************************************************************************
 * Function: process_narrow_peak
 *******************************************************************************
 */
int
process_narrow_peak(Data *data, int jc, int peak_ind, int l, int i,
    int *peak_recalled)
{
    int position, lpos, rpos;

    /* Assumption: l < i */
    position = data->color_data[jc].peak_list[peak_ind].ipos;
    lpos = data->bases.coordinate[l];
    rpos = data->bases.coordinate[i];
   *peak_recalled = 0;

    /* Call that one of the two peaks which is closer to the original
     * base coordinate
     */
    if (abs(position-lpos) < abs(position-rpos)) {

        /* Select the left base */
        if (data->bases.called_peak_list[i] != NULL) {
            data->bases.called_peak_list[i]->is_called=0;
            data->bases.called_peak_list[i]->base_index=-1;
        }
        data->bases.called_peak_list[i] = NULL;
        data->bases.bases[i] = 'N';

        data->bases.called_peak_list[l]->is_called = 122;
        data->bases.called_peak_list[l]->base_index = l;
    }
    else if (abs(position-lpos) > abs(position-rpos)) {
        /* Select the right one */
        data->bases.called_peak_list[l]->is_called = 0;
        data->bases.called_peak_list[l]->base_index = -1;
        data->bases.called_peak_list[l] = NULL;
        data->bases.bases[l] = 'N';
        data->color_data[jc].peak_list[peak_ind].base_index =i;
        data->bases.called_peak_list[i]=
            &data->color_data[jc].peak_list[peak_ind];
        data->bases.called_peak_list[i]->is_called = 122;
       *peak_recalled = 1;
    }
    else /* equally distanced */ {
        if (data->color_data[jc].data[lpos] >
            data->color_data[jc].data[rpos]) 
        {
            if (data->bases.called_peak_list[i] != NULL) {
                 data->bases.called_peak_list[i]->is_called=0;
                 data->bases.called_peak_list[i]->base_index=-1;
            }
            data->bases.called_peak_list[i]=NULL;
            data->bases.bases[i] = 'N';
            data->bases.called_peak_list[l]->is_called = 122;            
        }
        else {
            data->bases.called_peak_list[l]->is_called = 0;
            data->bases.called_peak_list[l]->base_index = -1; 
            data->bases.called_peak_list[l]=NULL;
            data->bases.bases[l] = 'N';
            data->color_data[jc].peak_list[peak_ind].is_called = 122;
            data->color_data[jc].peak_list[peak_ind].base_index =i;
            data->bases.called_peak_list[i] =
             &data->color_data[jc].peak_list[peak_ind];
           *peak_recalled = 1;
        }
    }
    return SUCCESS; 
}


/*******************************************************************************
 * Function: get_base_pos
 * Purpose:  given base index, return its intrinsic position.
 *           If base for this index was not called, return ABI base location
*******************************************************************************
 */
static int
get_base_pos(int base_index, Data *data)
{
   int bpos = (data->bases.called_peak_list[base_index]==NULL) ?
        data->bases.coordinate[base_index] :
        data->bases.called_peak_list[base_index]->ipos;
   return bpos;
}
/*******************************************************************************
 * Function: getBIF51
*******************************************************************************
 */
static double
getBIF51(int base_index)
{
    double currBIF51, MAX_BIF=1., MIN_BIF=1.;
    int BIF_IND1=500, BIF_IND2=650;

    if      (base_index < BIF_IND1)
        return MAX_BIF;
    else if (base_index > BIF_IND2)
        return MIN_BIF;
    else
        currBIF51= MAX_BIF -
                (MAX_BIF - MIN_BIF) /
                (double)(BIF_IND2  - BIF_IND1)   *
                (double)(base_index- BIF_IND1);
    return currBIF51;
}

/*******************************************************************************
 * Function: getBDF
*******************************************************************************
 */
static double
getBDF(int base_index)
{
    double currBDF, MAX_BDF=1.18, MIN_BDF=1.28;
    int BDF_IND1=500, BDF_IND2=650;

    if      (base_index < BDF_IND1)
        return MAX_BDF;
    else if (base_index > BDF_IND2)
        return MIN_BDF;
    else
        currBDF= MAX_BDF -
                (MAX_BDF - MIN_BDF) /
                (double)(BDF_IND2  - BDF_IND1)   *
                (double)(base_index- BDF_IND1);
    return currBDF;
}

/*******************************************************************************
 * Function: get_base_insertion_factor
*******************************************************************************
 */
static double
get_base_insertion_factor(int base_index, Options *options)
{
    double insertion_factor, max_insertion_factor=1.;

    if (options->shift)
        max_insertion_factor = MAX_INSERTION_FACTOR;

    if      (base_index < BASE_IND1)
        return max_insertion_factor;
    else if (base_index > BASE_IND2)
        return 1;
    else
        insertion_factor =  max_insertion_factor -
                (max_insertion_factor - 1.) /
                (double)(BASE_IND2 - BASE_IND1)   *
                (double)(base_index- BASE_IND1);
    return insertion_factor;
}

/*******************************************************************************
 * Function: can_delete_base
 * Purpose:  determine if current base can be deleted
 *******************************************************************************
 */
static int
can_delete_base(Data *data, int prev_base_index, int next_base_index,
    double ave_spacing, double base_deletion_factor, Options *options)
{
    int prev_pos, next_pos;    /* Two spacings instead of 1 */
    int level = 0;
    double num_spacings;

    /* Level 1
     *********
     */
    num_spacings = 2.;
    next_pos = get_base_pos(next_base_index, data);
    prev_pos = get_base_pos(prev_base_index, data);
    if ( INT_LT_DBL( (next_pos - prev_pos),
        (num_spacings*ave_spacing/ base_deletion_factor) ) )
        level += 1000;

    /* Level 2
     *********
     */
    num_spacings = 3.;
    if (prev_base_index > 0) {
        next_pos = get_base_pos(next_base_index  , data);
        prev_pos = get_base_pos(prev_base_index-1, data);
        if ( INT_LT_DBL( (next_pos - prev_pos),
            (num_spacings*ave_spacing/ base_deletion_factor) ) )
            level += 100;
        }

    if (next_base_index < data->bases.length-1) {
        next_pos = get_base_pos(next_base_index+1, data);
        prev_pos = get_base_pos(prev_base_index  , data);
        if ( INT_LT_DBL( (next_pos - prev_pos),
            (num_spacings * ave_spacing/ base_deletion_factor) ) )
            level += 100;
        }

    /* Level 3
     *********
     */
    num_spacings =  4.;
    if (prev_base_index > 1) {
            next_pos = get_base_pos(next_base_index  , data);
        prev_pos = get_base_pos(prev_base_index-2, data);
        if ( INT_LT_DBL( (next_pos - prev_pos),
            (num_spacings * ave_spacing/ base_deletion_factor) ) )
            level += 10;
        }

    if ((prev_base_index > 0) &&
                (next_base_index < data->bases.length-1)) {
        next_pos = get_base_pos(next_base_index+1, data);
        prev_pos = get_base_pos(prev_base_index-1, data);
        if ( INT_LT_DBL( (next_pos - prev_pos),
            (num_spacings * ave_spacing/ base_deletion_factor) ) )
            level += 10;
        }

    if (next_base_index < data->bases.length-2) {
        next_pos = get_base_pos(next_base_index+2, data);
        prev_pos = get_base_pos(prev_base_index  , data);
        if ( INT_LT_DBL( (next_pos - prev_pos),
            (num_spacings * ave_spacing/ base_deletion_factor) ) )
            level += 10;
        }

    /* Level 4
     ************************************
     */
    num_spacings = 5.;
    if (prev_base_index > 2) {
            next_pos = get_base_pos(next_base_index  , data);
        prev_pos = get_base_pos(prev_base_index-3, data);
        if ( INT_LT_DBL( (next_pos - prev_pos),
            (num_spacings * ave_spacing/ base_deletion_factor) ) )
            level += 1;
        }

    if ((prev_base_index > 1) &&
                (next_base_index < data->bases.length-1)) {
        next_pos = get_base_pos(next_base_index+1, data);
        prev_pos = get_base_pos(prev_base_index-2, data);
        if ( INT_LT_DBL( (next_pos - prev_pos),
            (num_spacings * ave_spacing/ base_deletion_factor) ) )
            level += 1;
        }

    if ((prev_base_index > 0) &&
                (next_base_index < data->bases.length-2)) {
        next_pos = get_base_pos(next_base_index+2, data);
        prev_pos = get_base_pos(prev_base_index-1, data);
        if ( INT_LT_DBL( (next_pos - prev_pos),
            (num_spacings * ave_spacing/ base_deletion_factor) ) )
            level += 1;
        }

    if (next_base_index < data->bases.length-3) {
        next_pos = get_base_pos(next_base_index+3, data);
        prev_pos = get_base_pos(prev_base_index  , data);
        if ( INT_LT_DBL( (next_pos - prev_pos),
            (num_spacings * ave_spacing/ base_deletion_factor) ) )
            level += 1;
        }
#if 0
    if ((level>=1000) && (level%1000>=100) && (level%100>=10) && (level%10>=1))
        return 4;
    if ((level <1000) && (level%1000>=100) && (level%100>=10) && (level%10>=1))
        return 3;
    if ((level <1000) && (level%1000 <100) && (level%100>=10) && (level%10>=1))
        return 2;
    if ((level <1000) && (level%1000 <100) && (level%100 <10) && (level%10>=1))
        return 1;
#endif
    if (level >= 1111) return 4;
    if (level >=  211) return 3;
    if (level >=   31) return 2;
    if (level >=    4) return 1;

    return 0;
}

/*******************************************************************************
 * Function: can_insert_base      
 * Purpose:  determine if there's enough room for insertion of a new called base 
 * Returns:  0, 1, 2, 3 or 4
 *           4 - inequality satisfied at all levels
 *           3 -    -"-        -"-    at all but first level, etc.
 *******************************************************************************
 */
int
can_insert_base(Data *data, int prev_base_index, int next_base_index,     
    int num_add_peaks, double ave_spacing, double base_insertion_factor, 
    Options *options)
{
    int prev_pos, next_pos, num_spacings;    /* Two spacings instead of 1 */
    int level = 0;

    if (base_insertion_factor <= 0)
        base_insertion_factor = 
            get_base_insertion_factor(prev_base_index, options);

    /* Level 1 
     *********
     */
    num_spacings = num_add_peaks+1;
    next_pos = get_base_pos(next_base_index, data);
    prev_pos = get_base_pos(prev_base_index, data);    
    if ( INT_GT_DBL( (next_pos - prev_pos), 
        (num_spacings*ave_spacing/base_insertion_factor) ) )       
        level += 1000;

    /* Level 2
     *********
     */
    num_spacings = num_add_peaks+2;
    if (prev_base_index > 0) {
        next_pos = get_base_pos(next_base_index  , data);
        prev_pos = get_base_pos(prev_base_index-1, data);        
        if ( INT_GT_DBL( (next_pos - prev_pos), 
            (num_spacings*ave_spacing/base_insertion_factor) ) )           
            level += 100;
	}

    if (next_base_index < data->bases.length-1) {
        next_pos = get_base_pos(next_base_index+1, data);
        prev_pos = get_base_pos(prev_base_index  , data);        
        if ( INT_GT_DBL( (next_pos - prev_pos), 
            (num_spacings * ave_spacing/base_insertion_factor) ) ) 
            level += 100;
	}    

    /* Level 3
     *********
     */
    num_spacings =  num_add_peaks+3;

    if (prev_base_index > 1) {
	    next_pos = get_base_pos(next_base_index  , data);
        prev_pos = get_base_pos(prev_base_index-2, data);
        if ( INT_GT_DBL( (next_pos - prev_pos), 
            (num_spacings * ave_spacing/base_insertion_factor) ) )
            level += 10;
	}

    if ((prev_base_index > 0) &&
		(next_base_index < data->bases.length-1)) {
        next_pos = get_base_pos(next_base_index+1, data);
        prev_pos = get_base_pos(prev_base_index-1, data);
        if ( INT_GT_DBL( (next_pos - prev_pos), 
            (num_spacings * ave_spacing/base_insertion_factor) ) )
            level += 10;
	}

    if (next_base_index < data->bases.length-2) {
        next_pos = get_base_pos(next_base_index+2, data);
        prev_pos = get_base_pos(prev_base_index  , data);
        if ( INT_GT_DBL( (next_pos - prev_pos), 
            (num_spacings * ave_spacing/base_insertion_factor) ) )
            level += 10;
	}    

    /* Level 4
     ************************************     
     */
    num_spacings = num_add_peaks+4;
    if (prev_base_index > 2) {
	    next_pos = get_base_pos(next_base_index  , data);
        prev_pos = get_base_pos(prev_base_index-3, data);
        if ( INT_GT_DBL( (next_pos - prev_pos), 
            (num_spacings * ave_spacing/base_insertion_factor) ) )
            level += 1;
	}

    if ((prev_base_index > 1) &&
		(next_base_index < data->bases.length-1)) {
        next_pos = get_base_pos(next_base_index+1, data);
        prev_pos = get_base_pos(prev_base_index-2, data);
        if ( INT_GT_DBL( (next_pos - prev_pos), 
            (num_spacings * ave_spacing/base_insertion_factor) ) )
            level += 1;
	}

    if ((prev_base_index > 0) &&
		(next_base_index < data->bases.length-2)) {
        next_pos = get_base_pos(next_base_index+2, data);
        prev_pos = get_base_pos(prev_base_index-1, data);
        if ( INT_GT_DBL( (next_pos - prev_pos), 
            (num_spacings * ave_spacing/base_insertion_factor) ) )
            level += 1;
	}

    if (next_base_index < data->bases.length-3) {
        next_pos = get_base_pos(next_base_index+3, data);
        prev_pos = get_base_pos(prev_base_index  , data);
        if ( INT_GT_DBL( (next_pos - prev_pos), 
            (num_spacings * ave_spacing/base_insertion_factor) ) )
            level += 1;
	}
#if 0
    if ((level>=1000) && (level%1000>=100) && (level%100>=10) && (level%10>=1))
        return 4;
    if ((level <1000) && (level%1000>=100) && (level%100>=10) && (level%10>=1))
        return 3;
    if ((level <1000) && (level%1000 <100) && (level%100>=10) && (level%10>=1))
        return 2;
    if ((level <1000) && (level%1000 <100) && (level%100 <10) && (level%10>=1))
        return 1;
#endif
    if (level >= 1111) return 4;
    if (level >=  211) return 3;
    if (level >=   31) return 2;
    if (level >=    4) return 1;
    return 0;
}

/*******************************************************************************
 * Function: peak_outside_data_range
 *******************************************************************************
 */
static int
peak_outside_data_range(Peak *peak, Data *data)
{
    int       color = peak->color_index;
    ColorData *cd = &(data->color_data[color]);

    if ((peak->ipos < 0) || (peak->ipos > cd->length-1))
        return 1;

    return 0;
}

/*******************************************************************************
 * Function: check_spacing
 * Purpose:  perform a rough check of left and right spacings of current
 *           called peak by comparing them with expected spacing.
 *           If both of them are within 20% of expected_spacing, return 0
 *           If their sum is within 20% of 2 * expected_spacing, return 0
 *           If one of spacings is within 20% of expected_spacing,
 *              and the other is less than expected by > 20%, 
 *              return +1 (= potential insertion)
 *           If one of spacings is within 20% of expected_spacing,
 *              and the other is bigger than expected by > 20%, 
 *              return -1 (=potential deletion)
 *******************************************************************************
 */
static int
check_spacing(int data_peak_ind, Data *data, double expect_spacing)
{
   int base_index, ins=1, del=-1;
   double left_spacing, right_spacing, delta_left, delta_right;
   double DEL = .2;
 
   if (data->peak_list[data_peak_ind]->is_called <= 0) {
       fprintf(stderr, "Warning: checking spacing of uncalled peak %d\n",
           data_peak_ind);
       return 0;
   }
  
   if (expect_spacing <= 0) {
       fprintf(stderr, 
          "Warning: expected spacing = %f <= 0 for peak %d\n",
           expect_spacing, data_peak_ind);
       return ERROR;
   }

   base_index = data->peak_list[data_peak_ind]->base_index;

   /* Set left and right spacing */
   if (base_index==0) left_spacing = expect_spacing;
   else
       left_spacing  = 
          (double)(abs(data->bases.called_peak_list[base_index  ]->ipos-
                      data->bases.called_peak_list[base_index-1]->ipos));   
   if (base_index==data->bases.length-1) right_spacing = expect_spacing;   
   else
       right_spacing  =
          (double)(abs(data->bases.called_peak_list[base_index+1]->ipos-
                      data->bases.called_peak_list[base_index  ]->ipos));
  
   delta_left = 2*(left_spacing - expect_spacing) / 
                  (left_spacing + expect_spacing); 
   delta_right= 2*(right_spacing- expect_spacing) / 
                  (right_spacing+ expect_spacing);  
 
   if ((fabs(delta_left) < DEL) && (fabs(delta_right) < DEL)) return  0;
   if ((fabs(delta_left) < DEL) && (     delta_right >= DEL)) return del;
   if ((fabs(delta_left) < DEL) && (     delta_right <=-DEL)) return ins;
   if ((     delta_left >= DEL) && (fabs(delta_right) < DEL)) return del;
   if ((     delta_left <=-DEL) && (fabs(delta_right) < DEL)) return ins;
   if (      delta_left         +        delta_right <=-DEL ) return ins; 
         
   return del;
}

/*******************************************************************************
 * Function: get_worst_inserted_base
 *******************************************************************************
 */
static int
get_worst_inserted_base(int base_index, Data *data, double curr_spacing,
    Options *options)
{
    int i, j, worst_base=-1, spacing, spacing1, spacing2, worst_spacing=INF;
    int spacing_problems = 1;
    int shift[NUM_COLORS]={0, 0, 0, 0};

    if ((base_index <2) || (base_index > data->bases.length-3)) 
        return -1;

    j = base_index;
    i = data->bases.called_peak_list[base_index]->data_peak_ind;
    while (spacing_problems && 
           (j >= 2) &&
           (j < data->bases.length-2) &&
           (j-base_index < 3))
    {
        Peak *pkm2 = data->bases.called_peak_list[j-2];
        Peak *pkm1 = data->bases.called_peak_list[j-1];
        Peak *pk   = data->bases.called_peak_list[j  ];
        Peak *pkp1 = data->bases.called_peak_list[j+1];
        Peak *pkp2 = data->bases.called_peak_list[j+2];

        if (WORST_BASE_NOT_TYPE11 &
            (pk->type == 11) && 
            (pk->iheight > pkm1->iheight*FRAC_SHORT_CALLED_PEAK) && 
            (pk->iheight > pkp1->iheight*FRAC_SHORT_CALLED_PEAK))
        {
            i = data->bases.called_peak_list[j]->data_peak_ind;
            spacing_problems = check_spacing(i, data, curr_spacing);
            j++;
            continue;
        }

        if (WORST_BASE_NOT_DIP &&
            is_dip(data->peak_list, data->peak_list_len, i,
             shift, options, data) &&
            (pk->iheight > pkm1->iheight*FRAC_SHORT_CALLED_PEAK) &&
            (pk->iheight > pkp1->iheight*FRAC_SHORT_CALLED_PEAK))
        {
            i = data->bases.called_peak_list[j]->data_peak_ind;
            spacing_problems = check_spacing(i, data, curr_spacing);
            j++;
            continue;
        }

        if (can_insert_base(data, j-1, j+1, 1, curr_spacing, 
            getBIF51(j), options) >= 4)
        {
            i = data->bases.called_peak_list[j]->data_peak_ind;
            spacing_problems = check_spacing(i, data, curr_spacing);
            j++;
            continue;
        }

        spacing = pkp1->ipos - pkm1->ipos;
        spacing1= pkp2->ipos - pkm1->ipos;
        spacing2= pkp1->ipos - pkm2->ipos;

        if ((spacing + spacing1 + spacing2 < worst_spacing) &&
            (can_delete_base(data, j-1, j+1, curr_spacing, getBDF(j) , 
             options) >= 4))
        {
            worst_spacing = spacing + spacing1 + spacing2;
            worst_base    = j;
        } 
        j++; 
        i = data->bases.called_peak_list[j]->data_peak_ind;
        spacing_problems = check_spacing(i, data, curr_spacing); 
    }
   
    return worst_base;
}

/*******************************************************************************
 * Function: delete_uncalled_peak
 *******************************************************************************
 */
static int
delete_uncalled_peak(int cd_peak_ind, int data_peak_ind, ColorData *cd, 
    Data *data, BtkMessage *message)
{
    int i;

    if (cd->peak_list[cd_peak_ind].data_peak_ind != data_peak_ind) {
        fprintf(stderr, 
        "Misdefined data_peak_ind of cd peak %d in delete_uncalled_peak\n",
        cd_peak_ind);
        return ERROR;
    }

    if (data->peak_list[data_peak_ind]->cd_peak_ind != cd_peak_ind) {
        fprintf(stderr, 
        "Misdefined cd_peak_ind=%d of data peak %d in delete_uncalled_peak\n",
        data->peak_list[data_peak_ind]->cd_peak_ind, data_peak_ind);
        return ERROR;
    }

    if (data->peak_list[data_peak_ind]->data_peak_ind != data_peak_ind) {
        fprintf(stderr,
        "Misdefined data_peak_ind of data peak %d in delete_uncalled_peak\n",
        data_peak_ind);
        return ERROR;
    }

    /* Shift subsequent cd peaks to the left one spot to make room */
    (void)memmove(&cd->peak_list[cd_peak_ind  ],
                  &cd->peak_list[cd_peak_ind+1],
        (cd->peak_list_len-1 - cd_peak_ind) * sizeof(Peak));
    cd->peak_list_len--;

    /* Reset cd_peak_ind and data_peak_ind of shifted peaks
     * reset data peaks to shifted peaks 
     */
    for (i=cd_peak_ind; i<cd->peak_list_len; i++) {
        int k;
        cd->peak_list[i].cd_peak_ind = i;
        k = cd->peak_list[i].data_peak_ind;
        data->peak_list[k] = &cd->peak_list[i];
        if (cd->peak_list[i].is_called > 0) {
            k = cd->peak_list[i].base_index;
            data->bases.called_peak_list[k] =
                &cd->peak_list[i];
        }
    }

    /* Shift subsequent data peaks and update data_peak_ind */
    (void)memmove(&data->peak_list[data_peak_ind],
                  &data->peak_list[data_peak_ind+1],
        (data->peak_list_len-1 - data_peak_ind) * sizeof(Peak *));

    for (i=data_peak_ind; i<data->peak_list_len-1; i++) {
        data->peak_list[i]->data_peak_ind = 
            data->peak_list[i]->data_peak_ind2= i;
    }
    data->peak_list_len--;

    return SUCCESS;    
}

/*******************************************************************************
 * Function: merge_two_called_peaks
 *******************************************************************************
 */
static int
merge_two_called_peaks(int base_ind1, int base_ind2, Data *data,
    BtkMessage *message)
{
    int color = data->bases.called_peak_list[base_ind1]->color_index;
    int cd_peak_ind = data->bases.called_peak_list[base_ind1]->cd_peak_ind;
    ColorData *cd = &data->color_data[color];
    Peak peak1, peak2, peak;

    /* Assume: base_ind1 < base_ind2; */
    peak1 = *data->bases.called_peak_list[base_ind1];    
    peak2 = *data->bases.called_peak_list[base_ind2];
   
    if (peak1.color_index != peak2.color_index) {
        fprintf(stderr, "Attempt to merge two peaks of different colors\n");
        return ERROR;
    }

    if (peak1.cd_peak_ind != peak2.cd_peak_ind-1) {
        fprintf(stderr, "Attempt to merge two non-adjacent peaks\n");
        return ERROR;
    }
 
    peak         = peak1;
    peak.end     = peak2.end; 
    peak.iend    = peak2.iend;
    peak.area    = peak1.area + peak2.area;
    peak.relative_area = peak1.relative_area + peak2.relative_area;
    peak.pos     = get_peak_position(cd->data, peak.beg, peak.end, peak.area,
                   message);
//  peak.pos     = peak1.pos;
    peak.height  =         cd->data[peak.pos];
    peak.iheight = (double)cd->data[peak.pos];
    peak.ipos    = peak.pos;
    peak.spacing = peak1.spacing + (peak.pos - peak1.pos);
    peak.type    = peak1.type - peak1.type%10 + peak2.type%10; 
    peak.is_truncated = QVMAX(peak1.is_truncated, peak2.is_truncated);

    /* Insert new peak */
    cd->peak_list[cd_peak_ind] = peak;

    /* Now remove peak2, which is no longer needed.
     * Do it in 2 steps:
     *    - uncall the peak, remove the base
     *    - remove uncalled peak from colordata and data peak lists
     */
    if (uncall_peak(base_ind2, data, message) != SUCCESS)
        return ERROR;

    if (delete_uncalled_peak(peak2.cd_peak_ind, peak2.data_peak_ind,
        cd, data, message) != SUCCESS)
        return ERROR;     

    return SUCCESS;
}

/*******************************************************************************
 * Function: get_spacing_from_good_region
 *
 *******************************************************************************
 */
static double
get_spacing_from_good_region(int base_index, Data *data)
{
    int      i, j, new_spacing, max_spacing, max_pos, min_spacing, min_pos,
             num_spacings=0, sum_spacings=0, delta_spacing;
    double   ave_spacing;
    Peak   **Pl = data->bases.called_peak_list;

    i = base_index;

    if (base_index>=data->bases.length) {
        int ind = data->bases.length-1;
        return spacing_curve(data->bases.called_peak_list[ind]->ipos);
    }

    if (base_index <= MAX_SIZE_OF_SEARCH_REGION)
        return spacing_curve(data->bases.called_peak_list[base_index]->ipos);

    delta_spacing= 1;
    max_spacing  = 0;
    min_spacing  = INF;
    max_pos      = i;
    min_pos      = i;
    num_spacings = 0;
    sum_spacings = 0;

    j = i;
    while ((j             >  1) && 
           (i-j           <  MIN_SIZE_OF_GOOD_REGION) &&
           (base_index-j  <= MAX_SIZE_OF_SEARCH_REGION) &&
           (delta_spacing <= MAX_DELTA_SPACING)) 
    {
        j--;

        if (base_index - j >= MAX_SIZE_OF_SEARCH_REGION) {
            delta_spacing++;
            i = base_index;
            j = i;
            max_spacing  = 0;
            min_spacing  = INF;
            max_pos      = i;
            min_pos      = i;
            num_spacings = 0;
            sum_spacings = 0;
            continue;
        }
       
        new_spacing = Pl[j]->ipos - Pl[j-1]->ipos;
        if (new_spacing < 0          ) new_spacing =-new_spacing;
        if (max_spacing < new_spacing) {
            max_spacing = new_spacing;
            max_pos = j-1;
        }
        if (min_spacing > new_spacing) {
            min_spacing = new_spacing;
            min_pos = j-1;
        }

        if ((double)(max_spacing - min_spacing) <= delta_spacing)
        {
            num_spacings++;
            sum_spacings += new_spacing;

            if (i-j >= MIN_SIZE_OF_GOOD_REGION) {
                ave_spacing = (double)sum_spacings / (double)num_spacings;
                return ave_spacing;
            }
            else
            {
                continue;
            }
        }
        else {
            max_spacing  = 0;
            min_spacing  = INF;
            num_spacings = 0;
            sum_spacings = 0;
            i = QVMAX(max_pos, min_pos);
            j = i;
            max_pos = i;
            min_pos = i;
        }
    }

    if ((sum_spacings > 0.) && (j > 0))
        return sum_spacings/(double)j;

    new_spacing = QVMAX(MIN_PEAK_SPACING, (int)spacing_curve(Pl[base_index]->ipos));

    return new_spacing;
}

/*******************************************************************************
 * Function: split_observed_peak
 * Purpose: split observed peak into 2 virtual peaks.
 *          Do it a primitive way: set the position of the interface
 *          between the two peaks to the position of the original peak
 * Parameters:
 *     peak - original peak
 *     peak1, peak2 - pointers to virtual peaks
 *******************************************************************************
 */
static int
split_observed_peak(ColorData *cd, Peak peak, int prev_pos, int next_pos,
    Peak *peak1, Peak *peak2, BtkMessage *message)
{
    int bound_pos;

    /* Set position of the boundary between peak1 and peak2 */
    if (prev_pos > next_pos)
    {
        int temp_pos = next_pos;
        next_pos = prev_pos;
        prev_pos = temp_pos;
    }
    else if (prev_pos == next_pos)
    {
        if (peak.beg < prev_pos)
           prev_pos--;
        else if (next_pos < peak.end)
           next_pos++;
    }
    bound_pos = (prev_pos + next_pos)/2;
    if ((bound_pos <= peak.beg) || (bound_pos >= peak.end))
        bound_pos = (peak.beg + peak.end)/2;

    if (prev_pos > next_pos) {
        sprintf(message->text,
            "Swapped prev_pos and next_pos when splitting peak\n");
        return ERROR;
    }

    /* Reset peak1 */
   *peak1 = peak;
    peak1->type = peak.type - peak.type%10 +3;
    peak1->end = bound_pos;
    if (peak1->beg >= peak1->end) {
        peak1->beg  = peak1->end;
        peak1->end++;
        bound_pos = peak1->end;
    }
    peak1->max = get_peak_max(cd->data,  peak1->beg, peak1->end, message);
    peak1->area= get_peak_area(cd->data, peak1->beg, peak1->end, message);
    if (peak1->area < 0 ) {
        sprintf(message->text, "Error computing peak1.area while splitting\n");
        return ERROR;
    }
    peak1->pos = get_peak_position(cd->data, peak1->beg, peak1->end,
        peak1->area, message);
    peak1->ipos = peak1->pos;
    if ((peak1->pos < peak1->beg) || (peak1->pos > peak1->end)) {
        sprintf(message->text,
            "Error computing peak1.pos while splitting\n");
        return ERROR;
    }
    peak1->height = cd->data[peak1->pos];
    if (peak1->height==0) 
        peak1->height=1;
    peak1->iheight = peak1->height;

    /* Reset peak2 */
   *peak2 = peak;
    peak2->cd_peak_ind = peak.cd_peak_ind + 1;
    peak2->type = 30 + peak.type%10;
    peak2->beg = bound_pos;
    if (peak2->end <= peak2->beg) {
        peak2->end  = peak2->beg+1;
    }
    peak2->max = get_peak_max(cd->data,  peak2->beg, peak2->end, message);
    peak2->area= get_peak_area(cd->data, peak2->beg, peak2->end, message);
    if (peak2->area < 0 ) {
        sprintf(message->text,
            "Error computing peak2.area while splitting\n");
        return ERROR;
    }
    peak2->pos = get_peak_position(cd->data, peak2->beg, peak2->end, 
        peak2->area, message);
    peak2->ipos = peak2->pos;
    if ((peak2->pos < peak2->beg) || (peak2->pos > peak2->end)) {
        sprintf(message->text,
            "Error computing peak2.pos while splitting\n");
        return ERROR;
    }
    peak2->height = cd->data[peak2->pos];
    if (peak2->height==0)
        peak2->height=1;
    peak2->iheight = peak2->height;

    return SUCCESS;
}

/*******************************************************************************
 * Function: build_new_peak
 * Purpose: Build a new peak record in case where the base coordinate
 *          as provided by sample file is not within any existing peak's area
 *          Determine the closest left and right neibours of the new peak in
 *          the peak list as well as the index of the new peak.
 *          Set the new peak position to the base coordinate, peak.beg to the
 *          end of the previous peak in the peak list (left neibour) and
 *          peak.end to the beg next peak (right neibour)
 *******************************************************************************
 */
static int
build_new_peak(ColorData *cd, int beg, int pos, int end, 
    Peak *peak, int *peak_index, BtkMessage *message)
{
    int prev_type, next_type;

    if ((cd->dye_number < 1) || (cd->dye_number > 4)) {
        (void)printf("Dye Index=%d out of bounds: pos=%d\n",
                        cd->dye_number, pos);
        return ERROR;
    }
    peak->is_truncated = 0;
    peak->pos  = pos;
    peak->ipos = peak->pos;
    peak->max = peak->pos;
    peak->beg  = beg;
    peak->end  = end;
    peak->height = cd->data[pos];
    peak->iheight= peak->height;
    peak->area = get_peak_area(cd->data, peak->beg, peak->end, message);
    prev_type = (*peak_index > 0) ? cd->peak_list[*peak_index-1].type : 11; 
    next_type = (*peak_index < cd->peak_list_len-1) ? 
                cd->peak_list[*peak_index+1].type : 11;
    peak->type = (prev_type%10)*10 + (next_type-next_type%10)/10;
    peak->cd_peak_ind = *peak_index;

    /* Adjust boundaries of the previous peak  of the same color */
    if (*peak_index > 0 &&
        cd->peak_list[*peak_index-1].end > peak->beg)
    {
        cd->peak_list[*peak_index-1].end = peak->beg;
    }

    /* Adjust boundaries of the next peak of the same color */
    if (*peak_index <cd->peak_list_len-1 &&
        cd->peak_list[*peak_index].beg < peak->end)
    {
        cd->peak_list[*peak_index].beg = peak->end;
    }

    return SUCCESS;
}

/*******************************************************************************
 * Function: adjust_peak_position
 * Purpose: make sure the intrinsic position of the peak to be called is
 *          within the right range:
 *          1) between the previous and next called base
 *          2) between the two boundaries of the peak and
 *          3) between the positions of the previous and next peak 
 *             of the same color 
 *******************************************************************************
 */
void
adjust_peak_position(int base_ind, Peak *peak, int peak_index, int peak_beg, 
    int peak_end, int prev_called_pos, int next_called_pos, Data *data)
{
    ColorData *cd = &data->color_data[peak->color_index];
    double prev_ipos = (peak_index>0) ? cd->peak_list[peak_index-1].ipos : 1.;
    double next_ipos = (peak_index < cd->peak_list_len-1) ?
                        cd->peak_list[peak_index+1].ipos : data->length-1;

    double l_bound = QVMAX3(peak_beg, prev_ipos, prev_called_pos);
    double r_bound = QVMIN3(peak_end, next_ipos, next_called_pos);

    if (CHECK_REORDERING && l_bound >= r_bound)
        fprintf(stderr, 
            "Warning: can not adjust position of peak near base_ind=%d\n", 
            base_ind);
    peak->ipos = (peak->ipos > l_bound) ? peak->ipos : l_bound +1.;
    peak->ipos = (peak->ipos < r_bound) ? peak->ipos : r_bound -1.;
}

void
preset_base_calls(Data *data, char *color2base, Options options,
    BtkMessage *message)
{
    int i;
    char *bases = data->bases.bases;
    int len = data->bases.length;

    for (i=0; i<len; i++)
    {
        bases[i] = 'N';
    }
}

static int
get_best_peak_color_and_index_by_base_and_location(int *color_index,
    int *peak_index, int base, char *color2base, int location, 
    int prev_location, int next_location, Data *data, BtkMessage *message, 
    Options options)
{
    int    j, k, r, l, jc;
    int    new_peak_index, new_location = location;
    double height = 0., new_height;
    Peak   peak,        new_peak;

   *peak_index = -1;

    /* Find the color index of the base */
    for (j=0; j<NUM_COLORS; j++) {
        if (data->color_data[j].base == base)
            break;      /* j is the index */
    }

    /* Determine the index of the called peak.
     * Two cases:
     * A) base is A, C, G or T,   so that j <= 3 and
     * B) base is N or something, so that j == 3
     */   

    /* Case A : base is A, C, G or T */   
    if (base=='A' || base=='C' || base=='G' || base=='T') 
    {
        /*
         * Find the peak closest to the base coordinate in colordata
         * corresponding to the base
         */
        if (location <= data->color_data[j].length) {
            if (colordata_find_peak_index_by_location(&(data->color_data[j]),
                 location, &peak, peak_index, message) != SUCCESS)
            {
                return ERROR;
            }

            if (*peak_index >= 0)  /* we are done */
            {
                height = Shape(data->color_data[j].peak_list[*peak_index].C0,
                data->color_data[j].peak_list[*peak_index].orig_width/2.,
                data->color_data[j].peak_list[*peak_index].beta,
               (double)(data->color_data[j].peak_list[*peak_index].ipos -
                location), &options);
            }
            else {
                height = data->color_data[j].data[location];
            }
        }
       *color_index = j;
    }

    /* Case B: base is N */
    else
    {
       *color_index = jc = 0;
        colordata_find_peak_index_by_location(&(data->color_data[0]),
              location, &peak, peak_index, message);

        if (*peak_index < 0) {
            height = data->color_data[0].data[location];
        }
        else
        {
            height = Shape(data->color_data[0].peak_list[*peak_index].C0,
                data->color_data[0].peak_list[*peak_index].orig_width/2.,
                data->color_data[0].peak_list[*peak_index].beta,
               (double)(data->color_data[0].peak_list[*peak_index].ipos -
                location), &options);

            if (options.recallndb && is_dye_blob(location, &peak, data, 0))
            {
                if (peak.is_truncated)
                    height =
                        (float)data->color_data[0].data[location]
                        *DYE_BLOB_FRACTION;
                else
                {
                    height = get_modified_height(location, &peak,
                        &data->color_data[0]);
                }
            }
        }

        for (k=1; k<NUM_COLORS; k++) 
        {
            if (location > data->color_data[k].length) {
                continue;
            }

            r = colordata_find_peak_index_by_location(&(data->color_data[k]),
                    location, &new_peak, &new_peak_index,message);
            if(r==ERROR) return ERROR;

            if (new_peak_index < 0) {
                new_height = data->color_data[k].data[location];
            }
            else
            {
                new_height =
                        Shape(data->color_data[k].peak_list[new_peak_index].C0,
                    data->color_data[k].peak_list[new_peak_index].orig_width/2.,
                    data->color_data[k].peak_list[new_peak_index].beta,
                   (double)(data->color_data[k].peak_list[new_peak_index].ipos-
                    location), &options);

                if (options.recallndb &&
                    is_dye_blob(location, &new_peak, data, 0))
                {
                    if (new_peak.is_truncated)
                        new_height = (float)data->color_data[k].data[location]
                            *DYE_BLOB_FRACTION;
                    else
                    {
                        new_height = get_modified_height(location, &new_peak,
                            &data->color_data[k]);
                    }
                }
            }

            if (new_height > height) 
            {
               *color_index = k;
                jc = k;
               *peak_index=new_peak_index;
                height = new_height;
                peak = new_peak;
                new_location = location;
                continue;
            }
            else  /* new_height==height; */
            {  
                int llocation = location-1; 
                int rlocation = location+1;

                l = 1;
                while ( llocation > 0 || rlocation < data->color_data[k].length)
                {
                    if (data->color_data[jc].data[llocation] <
                        data->color_data[k ].data[llocation])
                    {
                        /* switch to the new color */
                        break;
                        jc = k;
                        new_location = llocation; 
                        continue;
                    }
                    if (data->color_data[jc].data[rlocation] <
                        data->color_data[k ].data[rlocation])
                    {
                        /* switch to the new color */
                        break;
                        jc = k;
                        new_location = rlocation;
                        continue;
                    }
                    l++;
                    llocation = location - l >= 0 ? location - l : 0;
                    rlocation = location + l < data->color_data[k].length ?
                                location + l : data->color_data[k].length-1;
                }
               *color_index = jc;
            }
        } /* end loop in k */
        if (colordata_find_peak_index_by_location(&(data->color_data[jc]),
            location, &peak, peak_index, message) != SUCCESS)
        {
            return ERROR;
        }
    } /* base is N */

    if (*peak_index >= 0 &&
         (peak.pos < prev_location || peak.pos > next_location))
           *peak_index = -1;

    return SUCCESS;
}

static void
check_peak_positions(Data *data)
{
    int jc, j;
    /* Check if there are peaks with location out of bounds */
    for (jc=0; jc<NUM_COLORS; jc++) 
    {
        for (j=0; j<data->color_data[jc   ].peak_list_len; j++) {
            if ((data->color_data[jc].peak_list[j].ipos <
                 data->color_data[jc].peak_list[j].beg)
                ||
                (data->color_data[jc].peak_list[j].ipos >
                 data->color_data[jc].peak_list[j].end))
            fprintf(stderr, "Before base calling: Color=%d Peak[%d] beg=%d pos=%f end=%d\n",
                jc, j, data->color_data[jc].peak_list[j].beg,
                data->color_data[jc].peak_list[j].ipos,
                data->color_data[jc].peak_list[j].end);
        }
    }
}

Peak
initialize_peak()
{
    Peak peak;
    peak.end = peak.beg = peak.area = peak.pos = peak.ipos = peak.max
             = peak.cd_peak_ind = peak.ipos_orig = peak.data_peak_ind = -1;
    peak.is_truncated = peak.type = peak.ave_sq_width2 = peak.ave_width2
             = peak.ave_width1 = peak.wiheight = 0;
    return peak;
}

int
call_peak_by_location_color_and_index(int location, int base_ind, int color, 
    char *color2base, int *peak_index, int prev_peak_index, int *Case,
    Data *data, BtkMessage *message, Options options)
{
    int i = base_ind, jc = color, j, l, r, k;
    static int peaks_added[4]={0,0,0,0};
    ColorData *cd = &data->color_data[jc];
    char    base = color2base[jc];
    Peak peak = initialize_peak(), peak1 = initialize_peak(),
                peak2 = initialize_peak();
    int *coord = data->bases.coordinate;
    int lbound = (i == 0) ? 0 : (coord[i] + coord[i-1])/2;
    int rbound = (i == data->bases.length-1) ?  
        data->color_data[jc].length-1 : (coord[i]+coord[i+1])/2;
    Peak **cpl = data->bases.called_peak_list;

    /* Initialization */
    peak.beg = 0;
    peak.end = INF;
    peak.area = 0.;
    peak.pos = -1;
    peak.ipos = -1;
    peak.max = -1;
    peak.is_truncated = 0;
    peak.type = -1;
    peak.cd_peak_ind = -1;

    if (*peak_index >= 0)
        peak = cd->peak_list[*peak_index];

    /* Case M1:
     **********
     * If the peak within which area the current base is located
     * is found and its index differs from the index of the most recently
     * called peak of the same color, mark the current peak called;
     * set the peak position equal to the called base coordinate
     */
    if ((*peak_index>=0) && (*peak_index != prev_peak_index) && 
        cd->peak_list[*peak_index].is_called == 0)
    {
        ColorData *cd = &data->color_data[jc];
        Peak *cdpl    = cd->peak_list;

       *Case = 1;
        peak.base   = base;
        peak.base_index = i;
        peak.is_called  = *Case;
        peak.area = get_peak_area(data->color_data[jc].data,
                peak.beg, peak.end, message);

        peak.height = cd->data[peak.pos];
        if (peak.ipos <  lbound) peak.ipos = lbound;     
        if (peak.ipos >= rbound) peak.ipos = rbound + 1;   

        if ((peak.ipos < peak.beg) || (peak.ipos > peak.end)) {
           (void)sprintf(message->text,
                "Position of peak (#%d) out of bounds (case M1)\n",
                         *peak_index);
            return ERROR;
        }
        peak.base = data->bases.bases[i] = color2base[jc];

        /* Insert peak into peak lists */
        cdpl[*peak_index] = peak;
        cpl[i]= &cdpl[*peak_index];
        cpl[i]->base_index = i;

#if SHOW_CASE
        fprintf(stderr,
            "i=%3d base=%c coord=%d ipos =%.1f Case= M%d rec_base=%c, peak_width=%f",
            i, base, data->bases.coordinate[i], peak.ipos,*Case,
            color2base[jc], peak.area/peak.height);
        fprintf(stderr," peak_index=%d\n", *peak_index);
#endif
    }  /* end of case M1 */

    /* Case M2:
     **********
     * If the peak within which area the base is located, is found, but the
     * most recent previously called base with the same letter was located
     * within area of this same peak, then split this observed peak into 2
     * different "real" peaks and create a new record in the list of peaks
     */
    else if ((*peak_index>=0) &&
             (*peak_index == prev_peak_index &&
              data->color_data[color].peak_list[*peak_index].is_called > 0))
    {
        double location2;

       *Case = 2;
        /* Modify the record for the previous "real" peak of the same color */
        j = jc;
        peak.is_called = *Case;
        l = data->color_data[j].peak_list[prev_peak_index].base_index;
        location2 = data->bases.called_peak_list[l]->ipos;

        r = split_observed_peak(cd, peak,
            (location2 < location) ? location2: location, 
            (location2 < location) ? location : location2,
            &peak1, &peak2, message);

        if (peak1.end <= peak1.beg) {
            peak1.end = peak1.beg+1;
            peak2.beg = peak1.end;
        }
        peak1.base = cd->peak_list[prev_peak_index].base;
        peak1.base_index = (l < i) ? l : i;
        peak1.is_truncated = peak.is_truncated;

        if (peak1.height==0)
        {
            peak1.height=1;
        }

        {
            int llbound = (l == 0) ? 0 : (coord[l] + coord[l-1])/2;
            int rrbound = l<data->bases.length-1 ?
                (coord[l]+coord[l+1])/2 : cd->length;
        
            if (l < i)
            {
                if (peak1.ipos < llbound) peak1.ipos = llbound;
                if (peak1.ipos >=rrbound) peak1.ipos = rrbound-1;
            }
            else
            {
                if (peak1.ipos <  lbound) peak1.ipos =  lbound;
                if (peak1.ipos >= rbound) peak1.ipos =  rbound-1;
            }
        }

        /* Insert peak1 into peak lists */
        data->color_data[j].peak_list[prev_peak_index] = peak1;
        data->bases.called_peak_list[l<i ? l : i] =
            &data->color_data[j].peak_list[prev_peak_index];
        data->bases.bases[l<i ? l : i] = base;

        /* Create new record for the current "real" peak
         * (This should work however many "real" peaks
         * merge into one "observed" peak, assuming that the
         * coordinates of the called bases are ordered;
         * each time we will add one more "real" peak)
         */
      (*peak_index)++;
        peak2.is_truncated = peak.is_truncated;
        peak2.base_index = i;
        peak2.base = base;
        peak2.cd_peak_ind = *peak_index;
        if (peak2.height==0)
        {
            peak2.height=1;
            peak2.iheight = peak2.height;
        }

        if (peak2.ipos < lbound) peak2.ipos = lbound;
        if (peak2.ipos >=rbound) peak2.ipos = rbound-1;
        {
            int llbound = (l == 0) ? 0 : (coord[l] + coord[l-1])/2;
            int rrbound = l<data->bases.length-1 ?
                (coord[l]+coord[l+1])/2 : cd->length;

            if (l > i)
            {
                if (peak2.ipos < llbound) peak2.ipos = llbound;
                if (peak2.ipos >=rrbound) peak2.ipos = rrbound-1;
            }
            else
            {
                if (peak2.ipos <  lbound) peak2.ipos =  lbound;
                if (peak2.ipos >= rbound) peak2.ipos =  rbound-1;
            }
        }

        cd->peak_list_len++;

        /* Shift any peaks to the right one spot to make room */
        (void)memmove(&data->color_data[j].peak_list[*peak_index+1],
                      &data->color_data[j].peak_list[*peak_index],
                      (data->color_data[j].peak_list_len - *peak_index - 1)
                          * sizeof(Peak));

        for (k=*peak_index+1; k<data->color_data[j].peak_list_len; k++)
        {
            data->color_data[j].peak_list[k].cd_peak_ind = k;
            if (data->color_data[j].peak_list[k].is_called > 0)
            {
                int base_ind=data->color_data[j].peak_list[k].base_index;  
                data->bases.called_peak_list[base_ind] =
                    &data->color_data[j].peak_list[k]; 
                data->bases.called_peak_list[base_ind]->base_index = base_ind;
            }
            else
                data->color_data[j].peak_list[k].base_index = -1;
        }

        /* Update base array */
        if (base != 'A' && base != 'C' && base != 'G' && base != 'T')
        {
            data->bases.bases[l<i?i:l] = color2base[j];
        }
        peak2.base = data->bases.bases[l<i?i:l];

        /* Insert peak2 into peak lists */
        cd->peak_list[*peak_index]=peak2;
        cpl[l<i?i:l] = &cd->peak_list[*peak_index];
        data->bases.bases[l<i?i:l] = base;
        peaks_added[j]++;

#if SHOW_CASE
        fprintf(stderr,
        "i=%3d base=%c coord=%d ipos1=%.1f Case= M%d rec_base=%c, peak_width=%f",
        l<i?l:i, base, data->bases.coordinate[l<i?l:i],peak1.ipos,*Case,
        color2base[jc], peak1.area/peak1.height);
        fprintf(stderr," peak_index=%d\n", prev_peak_index);

        fprintf(stderr,
        "i=%3d base=%c coord=%d ipos2=%.1f Case= M%d rec_base=%c, peak_width=%f",
        l<i?i:l,base,data->bases.coordinate[l<i?i:l],peak2.ipos,*Case,
        color2base[jc], peak2.area/peak2.height);
        fprintf(stderr," peak_index=%d\n", *peak_index);
#endif

    }  /* end of case M2 */

    /* Case M3:
     *********
     * If the peak within which area the base is located is not found,
     * create a new peak at the location of the called base,
     * then insert it into the peak_list and increment indices of all
     * subsequent peaks by 1
     */
    else
    {    /* peak_index < 0 */
       *Case=3;
        peak.color_index = jc;
       *peak_index = prev_peak_index + 1;
     
        if (location <  lbound) location  = lbound;
        if (location >= rbound) location  = rbound-1;

        build_new_peak(cd, lbound, location, rbound, &peak, peak_index, message);
        peaks_added[jc]++;

        peak.base   = color2base[jc];
        peak.base_index = i;
        peak.is_called = *Case;
        peak.resolution = 1.;
        peak.width2=peak.end-peak.beg;
        peak.height = cd->data[peak.pos];
        if (peak.height==0) peak.height=1;
        peak.iheight = peak.height;
        peak.beta = i>0?data->bases.called_peak_list[i-1]->beta:1.;
        peak.orig_width = i>0?data->bases.called_peak_list[i-1]->beta:1.;
        peak.ave_w02beta = i>0?data->bases.called_peak_list[i-1]->ave_w02beta:1.;
        peak.C0   = peak.iheight;

        if ((*peak_index==0) || (cd->peak_list[*peak_index-1].area == 0.0))
        {
            peak.relative_area = 1.;
        }
        else {
            peak.relative_area =
                cd->peak_list[*peak_index-1].relative_area /
                cd->peak_list[*peak_index-1].area * peak.area;    /* a rough estimate */
        }

        /* Shift any peaks to the right one spot to make room */
        cd->peak_list_len++;
        (void)memmove(&cd->peak_list[*peak_index+1],
            &cd->peak_list[*peak_index], (cd->peak_list_len - *peak_index - 1)
                     * sizeof(Peak));
        for (j=*peak_index+1; j<cd->peak_list_len; j++)
        {
            int bi = cd->peak_list[j].base_index;
            if (bi >= 0)
            {
                data->bases.called_peak_list[bi] = &cd->peak_list[j];
            }
            cd->peak_list[j].cd_peak_ind = j;
        }

        /* Insert the new peak into the list */
        cd->peak_list[*peak_index]=peak;
        data->bases.called_peak_list[i] =
            &cd->peak_list[*peak_index];
        data->bases.called_peak_list[i]->base_index = i;
        data->bases.called_peak_list[i]->cd_peak_ind = *peak_index;
        data->bases.bases[i] = color2base[jc];

#if SHOW_CASE
        fprintf(stderr,
        "i=%3d base=%c coord=%d ipos =%.1f Case= M%d rec_base=%c peak_index=%d ",
        i,base,data->bases.coordinate[i],peak.ipos,*Case,color2base[jc],
        *peak_index);
        fprintf(stderr," lbound=%d rbound=%d\n", lbound, rbound);
#endif
    } /* end of case M3 */
    return SUCCESS;
}

static int get_previous_peak_index(int i, int jc, int peak_ind,
    int location, Data *data, BtkMessage *msg)
{
    int j = i-1;
    ColorData *cd = &data->color_data[jc];

    /* case M1 or M2 */
    if (peak_ind >= 0)
    {
        if (i == 0) return -1;
        while (j > 0 && data->bases.called_peak_list[j]->color_index != jc)
            j--;
        return data->bases.called_peak_list[j]->cd_peak_ind;
    }

    /* Case M3 */
    if (location < cd->peak_list[0].ipos) 
        return -1;
    if (location >= cd->peak_list[cd->peak_list_len-1].ipos)
        return cd->peak_list_len-1;
    for (j=1; j<cd->peak_list_len; j++)
    {
        if (cd->peak_list[j  ].ipos >= location &&
            cd->peak_list[j-1].ipos < location)
            return j-1;
    }
    return -1;
}

/*******************************************************************************
 * Function: mark_called_peaks
 * Purpose: In the list of all peaks, mark those which are called
 *          This function is similar to the data_mark_called_peaks in
 *          TraceTuner 1.0, but differs from it in that intrinsic peak
 *          height is used when recalling 'N's.
 *          The algorith assumes that it will be able to find a difference
 *          the relevant peaks. If it's not so, then the trace corresponding
 *          to the base 'T' is selected
 *
 *          This code includes a simple base swapping avoider, to make sure
 *          the order of recalled bases is the same as that of original bases.
 *          Specifically, we require that
 *
 *          called_peak_list[current_peak]->ipos >=
 *              called_peak_list[previous_peak]ipos
 *
 *          and
 *
 *          called_peak_list[current_peak]->ipos <=  coordinate[next_peak]
 *
 *******************************************************************************
 */
int
mark_called_peaks(Data *data, char *color2base, Options options, 
    BtkMessage *message)
{
    char  *debug;
    char   base;
    double max_res;
    int    i, jc;
    int    peak_index;
    int    previous_peak_index[4] = {-1,-1,-1,-1};
    int    prev_color_index = -1;
    int    Case, cd_index, cd_length;
    int    shift[4] = {0, 0, 0, 0};
    int   *coord = data->bases.coordinate;
    Peak **cpl = data->bases.called_peak_list;

    debug = getenv("DEBUG");

    if (CHECK_PEAK_POSITION)
        check_peak_positions(data);

    /* Infer called peak from orig. base and its coordinate */
    for(i=0; i<QVMIN(data->bases.length, MAX_NUM_BASES); i++) 
    {
        int location = coord[i];
        int prev_location = (i== 0) ? 0 :
              (cpl[i-1] == NULL ?  coord[i-1]+1 : cpl[i-1]->ipos+1);
        int next_location = (i == data->bases.length-1) ?
               data->color_data[0].length-1 :
              (data->bases.called_peak_list[i+1] == NULL ?
               data->bases.coordinate[i+1]-1 :
               data->bases.called_peak_list[i+1]->ipos-1);

        base = data->bases.bases[i];
        if (islower((int)base)) base=toupper((int)base); /* always use upper */

        Case=0;
        cd_index=-1;
        cd_length = -1;

        if (debug != NULL) {
            if ((i>0) && (data->bases.coordinate[i]==data->bases.coordinate[i-1]))
                printf("Bases %d and %d are called at the same coordinate=%d\n",
                i-1,i,data->bases.coordinate[i]);
        }

        if (get_best_peak_color_and_index_by_base_and_location(&jc, &peak_index,
            base, color2base, location, prev_location, next_location,
            data, message, options) != SUCCESS)
            return ERROR;

        if (peak_index < 0)
            previous_peak_index[jc] = get_previous_peak_index(i, jc, peak_index,
                        location, data, message);

        if (call_peak_by_location_color_and_index(location, i, jc, color2base,
            &peak_index, previous_peak_index[jc], &Case,
            data, message, options) != SUCCESS)
            return ERROR;

        previous_peak_index[jc] = peak_index;
        prev_location = location;
        if (prev_color_index >= 0 &&
           (previous_peak_index[prev_color_index] >= 0) &&
            (prev_color_index != jc) && 
             (Case == 2 || Case == 3))
        {
            ColorData *cd = &data->color_data[prev_color_index];
            int pind = previous_peak_index[prev_color_index];
            int first_peak_ind = pind;
            int last_peak_ind  = pind;
            while ((first_peak_ind>0) && (cd->peak_list[first_peak_ind].type > 20))
                first_peak_ind--;
            while ((last_peak_ind<cd->peak_list_len-1) &&
                (cd->peak_list[last_peak_ind].type%10 > 1))
                last_peak_ind++;

            resolve_multiple_peaks(data, prev_color_index, first_peak_ind,
                   last_peak_ind-first_peak_ind+1, &max_res, &options, message);
        }    
        prev_color_index = jc;

        if (base != 'A' && base != 'C' && base != 'G' && base != 'T')
        {
            data->bases.bases[i] = color2base[jc];
            data->color_data[jc].peak_list[peak_index].base =
                color2base[jc];
        }
    }                                      /* end loop through all bases */

    /* Adjust peak ipos and iheight */
    for(i=0; i<QVMIN(data->bases.length, MAX_NUM_BASES); i++)
    {
        Peak **pk = data->bases.called_peak_list;

        if (pk[i]->iheight < 1.)
        {
            pk[i]->iheight = 1.;
        }

        // If after resolve_multiple_peaks ipos is out of range, i
        // keep the old base coordinate
        // Otherwise, set the coordinate to ipos
        if ((i > 0 && 
             pk[i]->ipos < pk[i-1]->ipos) 
            ||
            (i<QVMIN(data->bases.length, MAX_NUM_BASES)-1 &&
            pk[i]->ipos > pk[i+1]->ipos))
        {
            pk[i]->ipos = data->bases.coordinate[i];
            continue;
        }

        data->bases.coordinate[i] = 
            ROUND(data->bases.called_peak_list[i]->ipos);
    }

    /* This call is done in order to set data_peak_ind for each peak,
     * which will be used when computing the trace parameters */
    if (bc_data_create_single_ordered_peak_list(data, shift, message) != SUCCESS)
    {
        sprintf(message->text, "Error creating single peak list\n");
        return ERROR;
    }

    if (CHECK_CD_PEAK_IND)     check_cd_peak_ind(data);
    if (CHECK_CALLED_PEAKS)    check_called_peaks(data);
    if (CHECK_BASE_INDEX)      check_base_index(data);
    if (CHECK_DATA_BASE_INDEX) check_data_base_index(data);
    if (CHECK_REORDERING)      check_reordering(data);
    if (CHECK_CD_PEAK_LISTS)  check_cd_peak_lists(data);

    return SUCCESS;
}

/*******************************************************************************
 * Function: Btk_call_bases             
 *
 * Purpose: 
 *     1) In the list of all peaks, mark those which are called
 *     2) Recall those 'N's which are located in some peak's area
 *     3) If no appropriate peak is found for the called base,
 *        remove the base from the list of called bases
 *     4) Change called locations: set them to the positions of the called
 *        peaks (if the peak is found)
 *     5) If additional "good" but uncalled peaks are found,
 *        call them and extend the list of called bases
 * 
 *******************************************************************************
 */
int
Btk_call_bases(Data *data, char *color2base, ReadInfo *read_info, 
    ContextTable *ctable, Options *options, BtkMessage *message, 
    Results *results )
{
    char    *debug, base, *context=NULL;
    int      r, i, j, jc = -1, pind, l, k, lbound, rbound;
    int      location, base_index, type; 
    int      jc2, height;
    int      peak_ind, next_peak_ind;
    int      shift[NUM_COLORS] = {0, 0, 0, 0};
    double   height2, nweight, cweight;
    Peak     peak, new_peak;
    char    *orig_bases;
    int     *Case;
    int      case_len = 2*data->bases.length;
    clock_t  start_clock = clock(), curr_clock;
#if CHANGE_BASE_ORDER
    color2base_orig[NUM_COLORS]="GATC";  /* base order in sample file */
    char    recalled_base;
#endif

    if (SHOW_CASE)
        fprintf(stderr, "Loop BC1\n");

    debug = getenv("DEBUG");

    orig_bases = CALLOC(char, data->bases.length);
    Case       = CALLOC(int, case_len);               

    for (i=0; i<data->bases.length; i++) {
        data->bases.called_peak_list[i] = NULL;
        orig_bases[i] = data->bases.bases[i];
    }

    if (ctable != NULL) {
        context = CALLOC(char, ctable->dimension);
    }

    for (i=0; i<data->bases.length; i++) 
    {

        /* BC1: 1st loop
         ******************************************
         * process only the called locations which correspond
         * to bases A, C, G, T and fall into area of a peak
         * of corresponding color
         */
        int *coord = data->bases.coordinate;
        Peak **cpl = data->bases.called_peak_list;
        ColorData *cd;

        base = data->bases.bases[i];  
        Case[i]=-1;
      
        /* Always use upper case */
        if (islower((int)base)) 
        {
            base=toupper((int)base); 
            data->bases.bases[i]=base;
        }

        data->bases.called_peak_list[i] = NULL;

        if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
            continue;
        }

        /* Find the color index of the base */
        for (jc = 0; jc < NUM_COLORS; jc++) {
            if(data->color_data[jc].base == base) {
                break;
            }
        }
        cd = &data->color_data[jc];
        lbound = (i == 0) ? 0 : (coord[i] + coord[i-1])/2;   
        rbound = (i == data->bases.length-1) ?  
                  data->color_data[jc].length-1 :
                  (coord[i]+coord[i+1])/2; 
 
        /* Does the original called location
         * fall into area of peak of the "right" color ?
         */
        peak_ind =-1;
        location = coord[i];  

        if (location <= data->color_data[jc].length-1) {
            if (colordata_find_peak_index_by_location(&(data->color_data[jc]),
                location, &peak, &peak_ind, message)==ERROR) {
                fprintf(stderr,
                "Error calling find_peak_index_by_location\n");
                FREE(orig_bases);
                FREE(Case);
                return ERROR;
            }
        }
        peak.color_index = jc; /* need only when peak_ind=-1 */

        /* BC1.0: no peak of a given color found at the original 
         *        location
         *******************************************************
         * Leave the bae as is; it will be considered in loop 3
         */
        if (peak_ind < 0)
        { 
            Case[i]= 10;
            if (SHOW_CASE) fprintf(stderr,
                "i=%3d        coord=%d             Case=BC%d rec_base=N \n",
                i, data->bases.coordinate[i], Case[i]);
            data->bases.bases[i] = 'N';
            continue;
        }

        /* BC1.1: original location falls into area of
         *        uncalled peak of the right color
         * ********************************************
         */
        else if ((peak_ind>=0) && (peak.is_called <= 0)) 
        {

            /* (Skipped for now)
             * BC1.1.1: too short peak 
             * ************************************************
             * Rename base to N; will reprocess it in the 2nd loop
             */
	    if (0 && peak.iheight < MIN_CALLED_PEAK_HEIGHT)    
            {
                data->bases.bases[i] = 'N';
		Case[i]= 111;
                if (SHOW_CASE) 
                    fprintf(stderr,
                    "i=%3d base=%c coord=%d            Case=BC%d rec_base=N\n",
                    i, data->bases.bases[i], coord[i], Case[i]);
				continue;
	    }

            /* (Skipped for now)
             * BC1.1.2: dye blob (= poorly resolved peak near
             *          the beginning of read)
             * ***********************************************
             * Rename base to N; will reprocess it in the 2nd loop
             */
            if (0 && is_dye_blob(location, &peak, data, 0)) 
            {
                data->bases.bases[i] = 'N';
                Case[i]= 112;
                if (SHOW_CASE) 
                    fprintf(stderr,
                    "i=%3d base=N coord=%d             Case=BC%d type=%d beg=%d end=%d\n",
                    i, data->bases.coordinate[i], Case[i], peak.type, peak.beg,
                peak.end);
                continue;
            }

            jc2    = -1;
            height = 0;
            for (j=0; j<NUM_COLORS; j++) {
                if (j == jc)
                    continue;

                if (data->color_data[j].data[location] > height) {
                    height = data->color_data[j].data[location];
                    jc2 = j;
                }
            }

            /* (Skipped for now)
             * BC1.1.5: there's a non-dye-blob peak with higher
             *          intrinsic signal at ABI base location
             * ***************************************************
             * Rename base to N; will reprocess it in the 2nd loop
             */ 
	    if (0 && jc2 >=0 ) 
            {
		next_peak_ind = -1;
		if (colordata_find_peak_index_by_location(
			&(data->color_data[jc2]), location, 
			&new_peak, &next_peak_ind, message) == ERROR) 
                {
		    fprintf(stderr,        
			"Error calling find_peak_index_by_location\n");
                        FREE(orig_bases);
                        FREE(Case);
			return ERROR;
	        }
	        if (((next_peak_ind >= 0) && 
                     (peak.type == 11) && (new_peak.type == 11) && 
                     (data->color_data[jc].data[location] <
                      data->color_data[jc2].data[location]))
                        ||
                    /* Peaks are not of type 11 */
                     ((next_peak_ind >= 0) &&
                      (!is_dye_blob(location, &new_peak, data, 0)) &&
                      (peak_with_low_boundaries(&peak, 
                          &data->color_data[peak.color_index], 0) <=
                       peak_with_low_boundaries(&new_peak, 
                          &data->color_data[new_peak.color_index], 0)) &&
                       (Shape(    peak.C0,     peak.orig_width/2.,     peak.beta, 
                           (double)(    peak.ipos - location), options) < 
                        Shape(new_peak.C0, new_peak.orig_width/2., new_peak.beta, 
                           (double)(new_peak.ipos - location), options)))) 
                 {

                    data->bases.bases[i] = 'N';
                    Case[i]= 115;
                    if (SHOW_CASE) fprintf(stderr,
                        "i=%3d base=N coord=%d             Case=BC%d \n",
                        i, data->bases.coordinate[i], Case[i]);
                    continue;
                }
            }
 
            { 
                /* BC1.1.6: good peak! (default case)
                ************************************
                * mark the peak called; set the called location to the peak 
                * position. If the peak position is out of the boundaries
                * of vicinity of the original location, adjust it
                */
                peak.base   = base;
                peak.base_index = i;
                Case[i]=116;
                peak.is_called = Case[i];
                if (peak.ipos <  lbound)
                    peak.ipos =  lbound;
                if (peak.ipos >= rbound)
                    peak.ipos =  rbound;

               /* Modify the peak list record */
                cd->peak_list[peak_ind]=peak;
                cpl[i]= &cd->peak_list[peak_ind];
                cpl[i]->cd_peak_ind = peak_ind;
                cpl[i]->base_index = i;

                if (SHOW_CASE) fprintf(stderr,
                "i=%3d base=%c coord=%d ipos=%.1f Case=BC%d \n",
                i, data->bases.bases[i], data->bases.coordinate[i], 
                data->bases.called_peak_list[i]->ipos, Case[i]);
            }
            continue;
        }   /* end of Case BC1.1 */

        /* BC1.2: original location falls into area of
         *        already called peak of the right color
         * ********************************************
         */

        else if (peak_ind>=0 && cd->peak_list[peak_ind].is_called > 0)
        {
            ColorData *cd = &data->color_data[jc];
            l = cd->peak_list[peak_ind].base_index; 
            peak = cd->peak_list[peak_ind];
            
            /* BC1.2.1: the two base locs fall into are of the peak
             *          which is good enough for two bases     
             *******************************************************
             * Split the peak into 2 peaks and call each
             */  
            if (1) /* Force this for now */
            {
                Peak peak1, peak2;
         
                Case[i]= 121;
                peak.is_called = Case[i];

                /* In this case, l < i goes automatically */

                r = split_observed_peak(cd, peak, coord[l], coord[i], 
                    &peak1, &peak2, message);
                peak1.base_index = peak.base_index;
                peak2.base_index = i;

                if (insert_and_resolve_peaks(data, jc, l, i, peak_ind,
                    -1, peak1, peak2, Case[i], options, message) != SUCCESS)
                    return ERROR;
               
                /* Adjust peak positions after resolving */
                {
                    /* Assumption: l < i */
                    int llbound = (l ==  0) ? 0 : (coord[l] + coord[l-1])/2;
                    int rrbound = (coord[l]+coord[l+1])/2; 
                    Peak *prev_pk = &cd->peak_list[peak_ind];
                    Peak *curr_pk = &cd->peak_list[peak_ind+1];

                    if (prev_pk->ipos < llbound) prev_pk->ipos = llbound;
                    if (prev_pk->ipos >=rrbound) prev_pk->ipos = rrbound-1;    
                    if (curr_pk->ipos <  lbound) curr_pk->ipos =  lbound;
                    if (curr_pk->ipos >= rbound) curr_pk->ipos =  rbound-1;
                } 
                cpl[i]->iheight = cd->data[(int)cpl[i]->ipos];
                cpl[l]->iheight = cd->data[(int)cpl[l]->ipos];

                if (SHOW_CASE) fprintf(stderr,
                    "i=%3d base=%c coord=%d ipos=%.1f Case=BC%d \n",
                    l, data->bases.bases[l], coord[l], cpl[l]->ipos, Case[i]);

                if (SHOW_CASE) fprintf(stderr,
                    "i=%3d base=%c coord=%d ipos=%.1f Case=BC%d \n",
                    i, data->bases.bases[i], coord[i], cpl[i]->ipos, Case[i]);

                continue;
            }  /* end of Case BC1.2.1 */

            else 
            {
       
                /* BC1.2.2: peak is too narrow for two bases      
                 *************************************************
                 * assign one of the bases to a peak and convert
                 * the other base to 'N'
                 */
                int peak_recalled = 0;
                Case[i]=122;

                if (l < i)
                {
                    if (process_narrow_peak(data, jc, peak_ind, l, i, 
                        &peak_recalled) != SUCCESS)
                    return ERROR;
                }
                else
                {
                    if (process_narrow_peak(data, jc, peak_ind, i, l,
                        &peak_recalled) != SUCCESS)
                    return ERROR;
                }

                if (SHOW_CASE && !peak_recalled) fprintf(stderr,
                    "i=%3d base=%c coord=%d             Case=BC%d \n",
                    i, data->bases.bases[i], data->bases.coordinate[i],
                    Case[i]);

                if (SHOW_CASE && peak_recalled) fprintf(stderr,
                    "i=%3d base=%c coord=%d ipos=%.1f Case=BC%d \n",
                    i, data->bases.bases[i], data->bases.coordinate[i], 
                    data->bases.called_peak_list[i]->ipos, 
                    Case[i]);

            } /* end Case BC1.2.2 (two bases in narrow peak) */
        } /* end Case BC1.2 (peak_ind >= 0 && peak.is_called) */
        if (CHECK_CALLED_PEAKS)    check_called_peaks(data);
    } /* end loop BC1 */

    if (CHECK_REORDERING) check_reordering(data);
    if (CHECK_CD_PEAK_IND) check_cd_peak_ind(data);
    if (CHECK_CD_PEAK_LISTS)  check_cd_peak_lists(data);

    if (SHOW_CASE)   
        fprintf(stderr, "Loop BC2\n");

    if (options->renorm)
        get_coefficients( data, read_info );

    /* BC2: 2nd loop
     **************************************************
     * recall the 'N'-s that fall into some peak's area
     */
    for (i = 0; i < data->bases.length; i++) 
    {
        int *coord = data->bases.coordinate;
        base = data->bases.bases[i];
        Case[i] = -2;

        if (base == 'A' || base == 'C' || base == 'G' || base == 'T') {
            continue;
        }

        /* Determine the range between previous and next base
         * to make sure the current base is not outside of it 
         */
        lbound = (i == 0) ? 0 : (coord[i] + coord[i-1])/2;
        rbound = (i == data->bases.length-1) ?
                data->color_data[jc].length-1 :
               (coord[i]+coord[i+1])/2;

        /* Does the base location fall into a good peak's area? */
        location = data->bases.coordinate[i];

        /* Correct the base location if bases are swapped */
        {
            int lloc = 0;
            int rloc = data->length-1;
            if ((i>0) && (data->bases.called_peak_list[i-1] != NULL))
                lloc = data->bases.called_peak_list[i-1]->ipos;
            if (i<data->bases.length-1)
                rloc = data->bases.coordinate[i+1];
            if ((location <= lloc) || (location >= rloc))
            {
                location = QVMAX(lloc + 1, (lloc+rloc)/2);
            }
        }

        jc  = -1;
        peak_ind = -1;
        type = -1;
        if ((ctable != NULL) && (i >= ctable->dimension-1)) { 
            memcpy(context, &data->bases.bases[i -
                (ctable->dimension-1) + 1], ctable->dimension-1);
        }

        /* Find a peak of the best color 
         * (with strongest intrinsic signal)
         * at the original called location
         */
        for (j = 0; j < NUM_COLORS; j++) 
        {
            new_peak.ipos = -1;

            if (colordata_find_peak_index_by_location(&(data->color_data[j]),
                location, &new_peak, &pind, message)==ERROR) {
                fprintf(stderr,
                "Error calling find_peak_index_by_location\n");
                FREE(orig_bases);
                FREE(Case);
                return ERROR;
            }

            /* Consider only the case where new_peak is found 
             * at the ABI base location 
             */
            if (pind < 0) 
                continue;

            /* Ignore very short new_peaks */
            if (new_peak.iheight < MIN_CALLED_PEAK_HEIGHT) {
                continue;
            }

            /* At this step, we are now looking for peaks located within 
             * search region 
             */
            if ((new_peak.ipos < lbound) || (new_peak.ipos >= rbound)) 
            {
                pind = -1;
                continue;
            }

            if ((new_peak.is_called <= 0)
                 ||
                ((new_peak.is_called > 0) && 
                  (new_peak.area/new_peak.iheight > PEAK_SPLIT_FACTOR*
                   new_peak.ave_width2)
                )
               )
            {
                /* Select the new_peak with strongest signal at the ABI base location
                 */
                nweight = 1.;
                cweight = 1.;

                if (options->renorm && !options->het) {   
                    nweight = readNormFactor( read_info, 
                    data->color_data[j].peak_list[pind].color_index, 
                    location );
                }
                if ((ctable != NULL) && (i >= ctable->dimension-1)) {
                    context[ctable->dimension-1] = new_peak.base;
                    cweight = weight_from_context(&context[ctable->dimension-1],
                        ctable);
                } 

                if ((type < 0) || 
                    is_better_peak(&new_peak, &peak, location, data, options))
                {
                    peak_ind = pind;
                    jc       = j;
                    type     = new_peak.type;
                    peak     = new_peak;
                }
            } /* end if pind > 0 etc. */
        }   /* end loop in j */

        /* BC2.1: Peak of the best color found */
        if (peak_ind >= 0) 
        {
            ColorData *cd = &data->color_data[jc];
            Peak **cpl = data->bases.called_peak_list;

            base = color2base[jc];
            cd->peak_list[peak_ind] = peak;
            data->bases.bases[i] = base;
            peak.base = base;

            /* BC2.1.1: Peak of the best color is not called
             ***********************************************
             * Call it and include into called peak list
             */
            if (peak.is_called==0) 
            { 
                Case[i]=211;
                cd->peak_list[peak_ind].base = base;
                cd->peak_list[peak_ind].is_called = Case[i];
                cd->peak_list[peak_ind].base_index =i;
                
                cpl[i] = &cd->peak_list[peak_ind];
                cpl[i]->base_index = i;
                data->bases.bases[i] = base;

                /* Adjust peak positions after resolving */
                if (cpl[i]->ipos <  lbound) cpl[i]->ipos =  lbound;
                if (cpl[i]->ipos >= rbound) cpl[i]->ipos =  rbound-1;

                if (SHOW_CASE) fprintf(stderr,
                    "i=%3d base=%c coord=%d ipos=%.1f Case=BC%d lbound=%d rbound=%d\n",
                    i, data->bases.bases[i], data->bases.coordinate[i],
                    cpl[i]->ipos, Case[i], lbound, rbound);

            }  /* end Case BC2.1.1 (uncalled peak found) */

            /* BC2.1.2: the best peak is already called
             ******************************************
             * Split it and resolve the multiple peak
             */
            else 
            {

                /* BC2.1.2: the best peak is already called
                 ******************************************
                 * Split it and resolve the multiple peak             
                 */
                Peak peak1, peak2;

                /* Base index currently assigned to the called peak */
                l = data->color_data[jc].peak_list[peak_ind].base_index;
                 
                Case[i]=212;
                peak.is_called = Case[i];

                if (l < i)
                {
                    r = split_observed_peak(&data->color_data[jc],
                        peak, data->bases.coordinate[l],
                        data->bases.coordinate[i], &peak1, &peak2, message);

                    if (insert_and_resolve_peaks(data, jc, l, i, peak_ind,
                        -1, peak1, peak2, Case[i], options, message) != SUCCESS)
                        return ERROR;
                    
                    /* Adjust peak positions after resolving */
                    {
                        /* Assumption: l < i */
                        int llbound = (l ==  0) ? 0 : (coord[l] + coord[l-1])/2;
                        int rrbound = (coord[l]+coord[l+1])/2;
                        Peak *prev_pk = &cd->peak_list[peak_ind];
                        Peak *curr_pk = &cd->peak_list[peak_ind+1];

                        if (prev_pk->ipos < llbound) prev_pk->ipos = llbound;
                        if (prev_pk->ipos >=rrbound) prev_pk->ipos = rrbound-1;
                        if (curr_pk->ipos <  lbound) curr_pk->ipos =  lbound;
                        if (curr_pk->ipos >= rbound) curr_pk->ipos =  rbound-1;
                    }
                    cpl[i]->iheight = cd->data[(int)cpl[i]->ipos];
                    cpl[l]->iheight = cd->data[(int)cpl[l]->ipos];

                    if (SHOW_CASE) fprintf(stderr,
                        "i=%3d base=%c coord=%d ipos=%.1f Case=BC%d \n",
                        l, data->bases.bases[l], data->bases.coordinate[l],
                        data->bases.called_peak_list[l]->ipos, Case[i]);

                    if (SHOW_CASE) fprintf(stderr,
                        "i=%3d base=%c coord=%d ipos=%.1f Case=BC%d \n",
                        i, data->bases.bases[i], data->bases.coordinate[i],
                        data->bases.called_peak_list[i]->ipos, Case[i]);
                }
                else
                {
                    r = split_observed_peak(&data->color_data[jc],
                        peak, data->bases.coordinate[i], data->bases.coordinate[l],
                        &peak1, &peak2, message);

                    if (insert_and_resolve_peaks(data, jc, i, l, peak_ind,
                        -1, peak1, peak2, Case[i], options, message) != SUCCESS)
                        return ERROR;
                  
                    /* Adjust peak positions after resolving */
                    {
                        /* Assumption: l > i */
                        int llbound = (coord[l] + coord[l-1])/2;
                        int rrbound = l<data->bases.length-1 ? 
                                      (coord[l]+coord[l+1])/2 :
                                      cd->length;
                        Peak *prev_pk = &cd->peak_list[peak_ind];
                        Peak *curr_pk = &cd->peak_list[peak_ind+1];
    
                        if (prev_pk->ipos < llbound) prev_pk->ipos = llbound;
                        if (prev_pk->ipos >=rrbound) prev_pk->ipos = rrbound-1;
                        if (curr_pk->ipos <  lbound) curr_pk->ipos =  lbound;
                        if (curr_pk->ipos >= rbound) curr_pk->ipos =  rbound-1;
                    }
                    cpl[i]->iheight = cd->data[(int)cpl[i]->ipos];
                    cpl[l]->iheight = cd->data[(int)cpl[l]->ipos];
 
                    if (SHOW_CASE) fprintf(stderr,
                        "i=%3d base=%c coord=%d ipos=%.1f Case=BC%d \n",
                        i, data->bases.bases[i], data->bases.coordinate[i],
                        data->bases.called_peak_list[i]->ipos, Case[i]);

                    if (SHOW_CASE) fprintf(stderr,
                        "i=%3d base=%c coord=%d ipos=%.1f Case=BC%d \n",
                        l, data->bases.bases[l], data->bases.coordinate[l],
                        data->bases.called_peak_list[l]->ipos, Case[i]);
                
                }    
            } /* end BC2.1.2 (best peak already called) */
        } /* end BC2.1 (base location falls into peak's area) */
        
        else 
        {
 
            /* BC2.2: Base location falls into none peak's area
             *********************************************************
             * Leave it as is; it will be procrssed in loop 4
             */
             Case[i]=22;
             if (SHOW_CASE) fprintf(stderr,
                    "i=%3d base=%c coord=%d             Case=BC%d \n",
                    i, data->bases.bases[i], data->bases.coordinate[i], 
                    Case[i]);
             continue;
        }
        if (CHECK_CALLED_PEAKS)  check_called_peaks(data);
    } /* end loop BC2 */

    if (CHECK_REORDERING) check_reordering(data);
    if (CHECK_CD_PEAK_IND) check_cd_peak_ind(data);
    if (CHECK_CD_PEAK_LISTS)  check_cd_peak_lists(data);

    if (SHOW_CASE)   
        fprintf(stderr, "Loop BC3\n");

    /* BC3: 3rd loop
     ****************
     * process remaining A, C, G ,T: look for a good peak on
     * either side of the base coordinate, withing the search region
     */
    for (i = 0; i < data->bases.length; i++) 
    {
        ColorData *cd;
        int *coord = data->bases.coordinate;
        base = data->bases.bases[i];
        Case[i] = 3;

        /* Search region:
         */
        lbound = (i == 0) ? 0 : (coord[i] + coord[i-1])/2;
        rbound = (i == data->bases.length-1) ?
                  data->color_data[jc].length-1 :
                  (coord[i]+coord[i+1])/2;

        if ((data->bases.called_peak_list[i] != NULL &&
            (data->bases.called_peak_list[i]->is_called > 0)) ||
            (base != 'A' && base != 'C' && base != 'G' && base != 'T'))
        {
            continue;
        }
   
        /* Determine the color index of the base */
        for (jc=0; jc<NUM_COLORS; jc++) {
            if(data->color_data[jc].base == base) {
                break;
            }
        }
        cd = &data->color_data[jc];

        /* Search for a peak of the right color */
        location = lbound+1;
        peak.iheight = 0;
        peak_ind = -1;
        while (location < rbound)
        {
            pind = -1;
            if (colordata_find_peak_index_by_location(cd,
                location, &new_peak, &pind, message) == ERROR) {
                fprintf(stderr,
                "Error calling find_peak_index_by_location\n");
                FREE(orig_bases);
                FREE(Case);
                return ERROR;
            }

            if (pind > 0 && new_peak.iheight > peak.iheight)
            {
                peak = new_peak;
                peak_ind = pind;
                if (peak.ipos < lbound) 
                {
                    peak.ipos =  lbound;
                    peak.iheight = cd->data[(int)peak.ipos];
                }
                if (peak.ipos >= rbound)
                {
                    peak.ipos =  rbound-1;
                    peak.iheight = cd->data[(int)peak.ipos];
                }
            }
            location += 2;
        }

        /* BC3.1: best peak found and not called
         ***********************************
         * Assign the base to it
         */
        if (peak_ind >= 0 
//          &&  /* peak height > a given percentage 
//              of the average peak height */ 
            &&
            peak.is_called == 0) 
        {
            Case[i]= 31;
            cd->peak_list[peak_ind] = peak;
            cd->peak_list[peak_ind].is_called = Case[i];
            cd->peak_list[peak_ind].base_index = i;
            data->bases.called_peak_list[i] =
                &data->color_data[jc].peak_list[peak_ind];
            data->bases.bases[i] = peak.base;
            data->bases.called_peak_list[i]->base_index = i;

            if (SHOW_CASE) fprintf(stderr,
                    "i=%3d base=%c coord=%d ipos=%.1f Case=BC%d \n",
                    i, data->bases.bases[i], data->bases.coordinate[i], 
                    data->bases.called_peak_list[i]->ipos, Case[i]);
        }   /* end Case BC3.1 */
       
        /* BC3.2: best peak found and is already called
         **********************************************
         * Assign the base to it
         */
        else if (peak_ind >= 0 &&
//         ( /* peak height > a given percentage
//              of the average peak height */) &&
            peak.is_called > 0)
        {
            /* split already called peak */
        }


        /* BC3.3: No good peak is found on either side
             ***********************************************
             * call the base 'N'
             */
        else 
        {
            data->bases.bases[i] = 'N';
            Case[i]=32;
            if (SHOW_CASE) fprintf(stderr,
                "i=%3d base=%c coord=%d             Case=BC%d \n",
                i, data->bases.bases[i], data->bases.coordinate[i], Case[i]);
        } /* end Case BC3.2 */

        if (SHOW_CASE && (data->bases.called_peak_list[i] != NULL)) {
            fprintf(stderr, 
            "i=%3d base=%c coord=%d             Case=BC%d rec_base=%c peak_ind=%d peak_res=%6.3f\n",
            i, data->bases.bases[i], data->bases.coordinate[i],
            Case[i], color2base[jc], peak_ind, data->bases.called_peak_list[i]->resolution);
        }
        if (CHECK_CALLED_PEAKS) check_called_peaks(data);
    } /* end loop BC3 (in i through remaining A, C, G, T) */

    if (CHECK_REORDERING) check_reordering(data);
    if (CHECK_BASE_INDEX) check_base_index(data);
    if (CHECK_CD_PEAK_IND) check_cd_peak_ind(data);
    if (CHECK_NUM_CALLED_PEAKS) check_num_called_peaks(data);

    if (SHOW_CASE) 
        fprintf(stderr, "Loop BC4\n");


    /* BC4: 4th loop
     ***************************
     * recall the remaining 'N':
     * Find the highest good peak around
     */
    for (i = 0; i < data->bases.length; i++) 
    {
        int *coord = data->bases.coordinate;
        ColorData *cd;
        base = data->bases.bases[i];

        if (base == 'A' || base == 'C' || base == 'G' || base == 'T') {
            continue;
        }

        jc = -1;
        peak_ind = -1;
        if ((ctable != NULL) && (i >= ctable->dimension-1)) {
            memcpy(context, &data->bases.bases[i -
                (ctable->dimension-1) + 1], ctable->dimension-1);
        }

        /* Get the left / right boundary of the search region.
         * This will be the location of either previous / next original base
         * or previous /next called base
         */
        lbound = (i == 0) ? 0 : (coord[i] + coord[i-1])/2;
        rbound = (i == data->bases.length-1) ?
                  data->length-1 : (coord[i]+coord[i+1])/2;

        /* Find the best color and peak */
        for (j = 0; j < NUM_COLORS; j++) 
        {
            cd = &data->color_data[j];
            base = color2base[j];

            /* Scan the search region, find the best peak */
            location = lbound+1;
            while (location < rbound)
            {
                int height=0;

                if (cd->data[location] <= MIN_CALLED_PEAK_HEIGHT)
                {
                    location += 4;
                    continue;
                }
                pind = -1;
                if (colordata_find_peak_index_by_location(cd,
                    location, &new_peak, &pind, message) == ERROR) {
                    fprintf(stderr,
                    "Error calling find_peak_index_by_location\n");
                    FREE(orig_bases);
                    FREE(Case);
                    return ERROR;
                }

                /* Find the strongest trace signal at the peak position 
                 * to make sure it is not much greater than
                 * the peak's apparent height 
                 */
                if (pind > 0) {
                    height = 0;
                    for (jc2=0; jc2 < NUM_COLORS; jc2++) {
                        if (height < data->color_data[jc2].data[(int)new_peak.ipos])
                            height = data->color_data[jc2].data[(int)new_peak.ipos]; 
                    }
                }

                if ((pind > 0) && 
                    ( INT_GT_DBL(cd->data[(int)new_peak.ipos],
                                  height * SIGNAL_STRENGTH_FRACTION) ) &&
                    (new_peak.is_called==0)                            &&
                    !is_dye_blob(location, &new_peak, data, 0) &&
                    ((i==0) || 
                     (abs(new_peak.ipos-data->bases.called_peak_list[i-1]->ipos)
                     >=2))
                   )
                {
                    nweight = 1.; 
                    cweight = 1.;
                    if (options->renorm && !options->het) {
                        nweight = readNormFactor( read_info, 
                            new_peak.color_index, new_peak.ipos);
                    }
                    if ((ctable != NULL) && (i >= ctable->dimension-1)) {
                        context[ctable->dimension-1] = new_peak.base;
                        cweight = weight_from_context(
                            &context[ctable->dimension-1], ctable);
                    }
                    if ((pind > 0) && (new_peak.is_called <= 0) &&
                       ( (peak_ind < 0) 
                                       ||
                        ((pind != peak_ind) &&
                         is_better_peak(&new_peak, &peak, -1, data, 
                             options))))
                    {
                        peak = new_peak;
                        peak_ind = pind;
                        jc = j;
                        if (peak.ipos < lbound)
                            peak.ipos = lbound;
                        if (peak.ipos >= rbound)
                            peak.ipos =  rbound+1;
                    }
                }
                location += 2;
            } /* end search for a peak in the window */
        } /* end loop in j through all colors */

        if ((peak_ind >= 0) && 
            (data->color_data[jc].peak_list[peak_ind].iheight >= 
             MIN_CALLED_PEAK_HEIGHT)) 
        {
            Peak *pk;
            cd = &data->color_data[jc];

            /* BC4.1: best peak and color found 
             **********************************
             * Insert it into called peak list
             */
            Case[i] = 41;
            base = color2base[jc];
            pk   = &cd->peak_list[peak_ind];
            if (pk->ipos  < lbound) 
                pk->ipos  = lbound;
            if (pk->ipos >= rbound) 
                pk->ipos =  rbound-1;
            data->bases.bases[i] = base;
            pk->is_called = Case[i];
            pk->base_index = i;
            pk->cd_peak_ind = peak_ind;
            data->bases.called_peak_list[i] = pk;
 
            if (SHOW_CASE) fprintf(stderr,
                "i=%3d base=%c coord=%d ipos=%4.2f Case=BC%d \n",
                i, data->bases.bases[i], data->bases.coordinate[i], 
                data->bases.called_peak_list[i]->ipos, Case[i]);

        } /* end BC 4.1 */

        else 
        {
            if ((data->bases.called_peak_list[i] != NULL) &&
                (data->bases.called_peak_list[i]->is_called > 0)) {
                fprintf(stderr,
                    "ERROR: rejected base is called\n");
            }

            /* BC4.2: no good peak found around
             ****************************************************
             * Determine which trace has the highest intrinsic or 
             * apparent signal at the original base location and call
             * a peak in this trace.
             * NOTE: use apparent, not intrinsic signal!
             *       The intrinsic signal may not exist if
             *       the peak is too small to be detected
             */

            location = data->bases.coordinate[i];

	    /* Determine the best color */
	    height2 = (double)data->color_data[0].data[location];
	    jc      = 0;
            for (j=1; j<NUM_COLORS; j++)
            {
                if (location > data->color_data[j].length) {
                    continue;
                }

                /* Compare signals at the location */
                if (INT_GT_DBL(data->color_data[j].data[location], height2)) {
                    height2 = (double)data->color_data[j].data[location];
                    jc = j;
                    continue;
                }

                /* Compare signals at adjacent scans */
                {
                    int llocation = location-1, rlocation = location+1;
                    k = 1;
                    while ( llocation > 0 || rlocation < data->color_data[j].length)
                    {
                        if (data->color_data[jc].data[llocation] <
                            data->color_data[j ].data[llocation]
                            ||
                            data->color_data[jc].data[rlocation] <
                            data->color_data[j ].data[rlocation])
                        {
                            /* switch to the new color, height2 already right */
                            break;
                            jc = j;
                            continue;
                        }
                        k++;
                        llocation = location - k >= 0 ? location - k : 0;
                        rlocation = location + k < data->color_data[j].length ?
                                    location + k : data->color_data[j].length-1;
                    }
                }
            }

            /* Is there (any) peak in the best trace ? */
            if (jc >= 0) 
            {
                if (colordata_find_peak_index_by_location(&(data->color_data[jc]),
                    location, &peak, &peak_ind, message) != SUCCESS) {
                    fprintf(stderr,
                    "Error calling find_peak_index_by_location\n");
                    FREE(orig_bases);
                    FREE(Case);
                    return ERROR;
                }

                /* If the peak is found, make sure it is located between 
                 * the previous and next base
                 */
                if (peak_ind >=0) {
                    /* Determine the range between previous and next base
                     * to make sure the current base is not outside of it
                     */
                    lbound = (i == 0) ? 0 : (coord[i] + coord[i-1])/2;
                    rbound = (i == data->bases.length-1) ?
                        data->color_data[jc].length-1 :
                        (coord[i]+coord[i+1])/2;

                    if (peak.ipos < lbound)
                        peak.ipos = lbound;
                    if (peak.ipos >= rbound) 
                        peak.ipos =  rbound-1;
                }
            }

            if ((jc >= 0) && (peak_ind >= 0) && (peak.is_called < 0)) 
            {
                /* BC4.2.1: best trace is found and contains uncalled peak 
                 *          (maybe, not as good as we would like, but ...)
                 *********************************************************
                 * Call the peak
                 */
                Case[i]=421;
                peak.base   = color2base[jc];
                peak.base_index = i;
                peak.is_called = Case[i];
                peak.color_index = jc;
                peak.cd_peak_ind = peak_ind;
                if (peak.ipos < lbound)
                    peak.ipos = lbound;
                if (peak.ipos >= rbound)
                    peak.ipos =  rbound-1;

                data->color_data[jc].peak_list[peak_ind]=peak;
                data->bases.called_peak_list[i]=
                    &data->color_data[jc].peak_list[peak_ind];
                data->bases.bases[i] = peak.base;

                if (SHOW_CASE) fprintf(stderr,
                    "i=%3d base=%c coord=%d ipos=%4.2f Case=BC%d \n",
                    i, data->bases.bases[i], data->bases.coordinate[i], 
                    data->bases.called_peak_list[i]->ipos, Case[i]);
            }

            else if ((jc >=0) && 
                     (peak_ind >= 0) && 
                     (peak.is_called > 0) &&
                     (data->color_data[jc].data[
                      (data->bases.coordinate[i] > data->color_data[jc].length-1) ?
                       data->color_data[jc].length-1 : data->bases.coordinate[i]
                                               ] >= MIN_CALLED_PEAK_HEIGHT)
                    ) 
            {
                int llbound, rrbound;

                /* BC4.2.2: best trace is found and contains called peak
                 *************************************************************
                 * Split the peak 
                 */
                Peak peak1, peak2;
                Case[i]=422;
                cd = &data->color_data[jc];
                Peak **cpl = data->bases.called_peak_list;

                l = cd->peak_list[peak_ind].base_index;
                llbound = (l ==  0) ? 0 : (coord[l] + coord[l-1])/2;
                rrbound = l<data->bases.length-1 ?
                          (coord[l]+coord[l+1])/2 : cd->length;

                peak.is_called = Case[i];
                base = color2base[jc];
                peak.area = get_peak_area(cd->data, peak.beg,
                    peak.end, message);
                peak.pos = get_peak_position(cd->data,
                    peak.beg, peak.end, peak.area, message); 
                if (peak.pos < location)
                    (void)split_observed_peak(cd, peak,
                    peak.pos, location, &peak1, &peak2,
                    message);
                else
                    (void)split_observed_peak(&data->color_data[jc], peak,
                    location, peak.pos, &peak1, &peak2, message);

                if (l < i)
                {
                    if (insert_and_resolve_peaks(data, jc, l, i, peak_ind,
                    -1, peak1, peak2, Case[i], options, message) != SUCCESS)
                    return ERROR;

                    cpl[l] = &cd->peak_list[peak_ind];
                    cpl[i] = &cd->peak_list[peak_ind+1];
                    cpl[l]->base_index = l;
                    cpl[i]->base_index = i;
           
                    /* Adjust peak positions after resolving */
                    if (cpl[l]->ipos < llbound) cpl[l]->ipos = llbound;
                    if (cpl[l]->ipos >=rrbound) cpl[l]->ipos = rrbound-1;
                    if (cpl[i]->ipos <  lbound) cpl[i]->ipos =  lbound;
                    if (cpl[i]->ipos >= rbound) cpl[i]->ipos =  rbound-1;
                }
                else
                {
                    if (insert_and_resolve_peaks(data, jc, i, l, peak_ind,
                    -1, peak1, peak2, Case[i], options, message) != SUCCESS)
                    return ERROR;

                    cpl[i] = &cd->peak_list[peak_ind];
                    cpl[l] = &cd->peak_list[peak_ind+1];
                    cpl[l]->base_index = l;
                    cpl[i]->base_index = i;
 
                    /* Adjust peak positions after resolving */
                    if (cpl[l]->ipos < llbound) cpl[l]->ipos = llbound;
                    if (cpl[l]->ipos >=rrbound) cpl[l]->ipos = rrbound-1;
                    if (cpl[i]->ipos <  lbound) cpl[i]->ipos =  lbound;
                    if (cpl[i]->ipos >= rbound) cpl[i]->ipos =  rbound-1;
                }
                cpl[i]->iheight = cd->data[(int)cpl[i]->ipos];
                cpl[l]->iheight = cd->data[(int)cpl[l]->ipos];

                if (SHOW_CASE) 
                {
                    fprintf(stderr, "  Splitting called peak %d\n", l);
                    if (l < i)
                    {
                        fprintf(stderr,
                        "i=%3d base=%c coord=%d ipos=%4.2f Case=BC%d ",
                        l, data->bases.bases[l], data->bases.coordinate[l],
                        data->bases.called_peak_list[l]->ipos, Case[i]);
                        fprintf(stderr, "llbound= %d rrbound=%d \n",
                        llbound, rrbound);

                        fprintf(stderr,
                        "i=%3d base=%c coord=%d ipos=%4.2f Case=BC%d ",
                        i, data->bases.bases[i], data->bases.coordinate[i], 
                        data->bases.called_peak_list[i]->ipos, Case[i]);
                        fprintf(stderr, " lbound= %d  rbound=%d \n",
                        lbound, rbound);
                    }
                    else
                    {
                        fprintf(stderr,
                        "i=%3d base=%c coord=%d ipos=%4.2f Case=BC%d ",
                        i, data->bases.bases[i], data->bases.coordinate[i],
                        data->bases.called_peak_list[i]->ipos, Case[i]);
                        fprintf(stderr, " lbound= %d  rbound=%d \n",
                        lbound, rbound);
 
                        fprintf(stderr,
                        "i=%3d base=%c coord=%d ipos=%4.2f Case=BC%d ",
                        l, data->bases.bases[l], data->bases.coordinate[l],
                        data->bases.called_peak_list[l]->ipos, Case[i]);
                        fprintf(stderr, "llbound= %d rrbound=%d \n",
                        llbound, rrbound);
                    } 
                }
            }

            else 
            {
                /* BC4.2.4: there's no signal at the original base location
                 **********************************************************
                 */
                if (DELETE_BASES)
                {
                    /*
                     * BC4.2.4.1. Reject the base
                     */
                    Case[i]=424;
                    if (SHOW_CASE) fprintf(stderr,
                        "i=%3d base=%c coord=%d Case=BC%d: reject base \n",
                        i, data->bases.bases[i], data->bases.coordinate[i], Case[i]);

                    if (uncall_peak(i, data, message) != SUCCESS)
                        return ERROR; 

                    i--; /* will consider the current i again */
                } /* end BC4.2.4 (reject base) */
                else
                {
                    /* BC4.2.4.2. Recall to the best guess */
                    int prev_location = i>0 ? 
                        data->bases.called_peak_list[i-1]->ipos : 0.;
                    int next_location = i<data->bases.length-1 ?
                        data->bases.coordinate[i+1] : 
                        data->color_data[0].length;
                    int previous_peak_index;

                    if (get_best_peak_color_and_index_by_base_and_location(&jc, 
                        &peak_ind, base, color2base, location, prev_location, 
                        next_location, data, message, *options) != SUCCESS)
                        return ERROR;
                       
                    previous_peak_index = get_previous_peak_index(i, jc, peak_ind,
                        location, data, message); 
                 
                    if (call_peak_by_location_color_and_index(location, i, jc, 
                        color2base, &peak_ind, previous_peak_index, &Case[i],
                        data, message, *options) != SUCCESS)
                        return ERROR;
                }
            }
        } /* end BC4.2 (no good peak found around) */
        if (CHECK_BASE_INDEX)      check_base_index(data);
        if (CHECK_CALLED_PEAKS) check_called_peaks(data);
        if (CHECK_CD_PEAK_LISTS)  check_cd_peak_lists(data);
    } /* end BC4 (loop in i through all remaining Ns) */
    if (CHECK_NUM_CALLED_PEAKS) check_num_called_peaks(data);
    if (CHECK_CD_PEAK_LISTS)  check_cd_peak_lists(data);
    if (REORDER_BASES) (void)bc_reorder_called_bases_and_peaks(data, message);
    if (CHECK_CD_PEAK_LISTS)  check_cd_peak_lists(data);

    if (bc_data_create_single_ordered_peak_list(data, shift, message)
        != SUCCESS)
    {
        sprintf(message->text, "Error creating single peak list\n");
        FREE(orig_bases);
        FREE(Case);
        return ERROR;
    }
    if (CHECK_NUM_CALLED_PEAKS) check_num_called_peaks(data);
    if (CHECK_CD_PEAK_IND)  check_cd_peak_ind(data);
    if (CHECK_CALLED_PEAKS) check_called_peaks(data);
    if (CHECK_BASE_INDEX)   check_base_index(data);
    if (CHECK_DATA_BASE_INDEX) check_data_base_index(data);
    if (CHECK_REORDERING)   check_reordering(data);
    if (CHECK_CD_PEAK_LISTS)  check_cd_peak_lists(data);
 
    if (INSERT_BASES)
    {
        double  ave_spacing_abi=1., ave_spacing=1.;

        if (CHECK_NUM_CALLED_PEAKS) check_num_called_peaks(data);

        if (options->renorm) {
            get_weighted_peak_heights(data, color2base, ctable, NULL, options, 
			message); 
        }
        else {
            get_weighted_peak_heights(data, color2base, ctable, read_info, 
			options, message);
        }

        if (options->time) {        
            curr_clock = clock();
            fprintf(stderr, 
            "   First 4 base calling loops completed in %f sec. \n",
                (double)(curr_clock - start_clock)/(double)CLOCKS_PER_SEC);
            start_clock = curr_clock;
        }      

        get_peak_spacing(data, options, message);

        if (CHECK_CD_PEAK_IND) check_cd_peak_ind(data);

        if ( collect_some_stats( data, options, message, results )
             != SUCCESS ) {
            return ERROR;
        }

        if (options->time) {
            curr_clock = clock();
            fprintf(stderr, "   Mobility shifts corrected in %f sec. \n",
                (double)(curr_clock - start_clock)/(double)CLOCKS_PER_SEC);
            start_clock = curr_clock;
        }

        if (SHOW_CASE)
            fprintf(stderr, "Loop BC5\n");

        base_index = 0;

        for (i = 0; i < data->peak_list_len; i++) 
        {
            int    detect_spacing_problems = 0;
            double curr_spacing, curr_spacing1, curr_spacing2;

            if (data->peak_list[i]->is_called > 0)
                base_index = data->peak_list[i]->base_index + 1;            

            if (base_index > data->bases.length) 
                continue;
         
            curr_spacing1 = spacing_curve(data->peak_list[i]->ipos);
            curr_spacing2 = get_spacing_from_good_region(base_index, data);
         
            if (base_index < 500)       
                curr_spacing = curr_spacing2;
            else
                curr_spacing = (curr_spacing1 + curr_spacing2)/2.;

            if (SHOW_CASE && (data->peak_list[i]->is_called <= 0))
                fprintf(stderr,
                "      base=%c ipos=%.1f data_peak_ind=%d\n",
                    data->peak_list[i]->base, data->peak_list[i]->ipos, i);

            /* BC5: 5th loop
             ***************
             * Review all peaks and call more good / uncall more bad
             * peaks if possible
             */

            if (peak_outside_data_range(data->peak_list[i], data))
                continue;

             if (data->peak_list[i]->is_called > 0) {
                 detect_spacing_problems = check_spacing(i, data, curr_spacing); 
             }

            /* BC5.1: Called DIP 
             *************************
             * just update base_index
             */
            if ((data->peak_list[i]->is_called > 0) && 
                is_dip(data->peak_list, data->peak_list_len, i,  
                shift, options, data))
            {

                adjust_called_base_position(data->peak_list[i]->base_index, 
                    data);

                if (SHOW_CASE) {
                    if (i >= case_len)
                    {
                        case_len *= 2;
                        Case = REALLOC(Case, int, case_len);
                    }
                    Case[i]=51;
                    fprintf(stderr,
                    "i=%3d base=%c ipos=%.1f dpi=%4d is_dip=%d is_dye_blob=%d ",
                    data->peak_list[i]->base_index, data->peak_list[i]->base,
                    data->peak_list[i]->ipos, i, is_dip(data->peak_list, 
                    data->peak_list_len, i, shift, options, data), 
                    is_dye_blob(data->peak_list[i]->ipos,
                    data->peak_list[i], data, 0));
                    fprintf(stderr,
                    "type=%d Case=BC%d \n",  data->peak_list[i]->type, Case[i]);
                }

                /* BC5.1.1:
                 **********
                 * Last attempt to insert deleted base 
                 */
                if ((detect_spacing_problems < 0) && 
                    (((base_index > 1) &&
                     (can_insert_base(data, base_index-2, base_index-1, 1,
                        curr_spacing, getBIF51(base_index), options) == 4 )) 
                    ||
                    ((base_index > 2) &&
                     (can_insert_base(data, base_index-3, base_index-1, 2,
                        curr_spacing, getBIF51(base_index), options) == 4 )))
                   ) 
                {
                    /* If there is still a big gap to the left from current
                     * called peak, then build and call a new peak 
                     * at the middle of gap by choosing trace
                     * with highest apparent signal at the new peak position
                     */
                    int new_pos;
                    int jc=-1, data_peak_ind;
                    int height=0, best_height=0, best_pind=-1;
                    Peak peak1, peak2;

                    if ((base_index > 1) &&
                        (can_insert_base(data, base_index-2, base_index-1, 1,
                        curr_spacing, getBIF51(base_index), options) == 4 ))
                    {
                        new_pos = 
                        (data->peak_list[i]->ipos + 
                         data->bases.called_peak_list[base_index-2]->ipos)/2;
                    }
                    else {
                        int curr_pos= data->peak_list[i]->ipos,
                            prev1_pos = 
                            data->bases.called_peak_list[base_index-2]->ipos,
                            prev2_pos = 
                            data->bases.called_peak_list[base_index-3]->ipos;
             
                        if ((curr_pos-prev1_pos)>(prev1_pos-prev2_pos))
                            new_pos = (curr_pos+prev1_pos)/2;
                        else
                            new_pos = (prev1_pos+prev2_pos)/2;
                    }

                    for (j=0; j<NUM_COLORS; j++) 
                    {
                        pind =-1;

                        height = data->color_data[j].data[new_pos];
                        if (height <= 0) 
                            continue;

                        if (colordata_find_peak_index_by_location(&(
                            data->color_data[j]), new_pos, &new_peak,
                            &pind, message)==ERROR) 
                        {
                            fprintf(stderr,
                            "Error calling find_peak_index_by_location\n");
                            FREE(orig_bases);
                            FREE(Case);
                            return ERROR;
                        }

                        if (pind >= 0) {
                            int iheight =
                                (int)Shape(new_peak.C0, new_peak.orig_width/2.,
                                new_peak.beta, (double)(new_peak.ipos - new_pos),
                                options);
                            if ((iheight < height) && (iheight > 0))
                                height = iheight; 
                        }
            
                        if (height > best_height) {
                            jc = j;
                            best_pind   = pind;
                            best_height = height;
                            peak = new_peak;
                        }
                    }
                    pind = best_pind;
                    height = best_height;

                    /* BC5.1.1.1:
                     ***********
                     * Called peak detected at expected new location.
                     * Split and call it and insert both peaks
                     */
                    if ( BC5_SPLIT_PEAKS &&
                        (pind >= 0) && 
                        (peak.is_called > 0) &&
                        ((peak.base_index == base_index-2) ||
                         (peak.base_index == base_index-1))
                       ) 
                    {
                        int pos1, pos2;
                        int pbind = peak.base_index>0?peak.base_index-1:0;
                        int nbind = peak.base_index<data->bases.length-1?
                                    peak.base_index+1:data->bases.length-1;
 
                        pos1 = peak.beg + (peak.pos - peak.beg)/3;
                        if (pos1 <= data->bases.called_peak_list[pbind]->ipos)
                            pos1  = data->bases.called_peak_list[pbind]->ipos + 1; 
                        pos2 = peak.end - (peak.end - peak.pos)/3;
                        if (pos2 >= data->bases.called_peak_list[nbind]->ipos)
                            pos2  = data->bases.called_peak_list[nbind]->ipos-1; 
                        (void)split_observed_peak(&data->color_data[jc], peak,
                            pos1, pos2, &peak1, &peak2, message);
                        data_peak_ind = peak.data_peak_ind;

                        if (insert_and_resolve_peaks(data, jc, peak.base_index,
                            peak.base_index+1, pind, data_peak_ind, peak1, peak2, 
                            5111, options, message) != SUCCESS)
                            return ERROR;

                        if (SHOW_CASE && (base_index<data->bases.length)) 
                        {
                            if (i >= case_len)
                            {
                                case_len *= 2;
                                Case = REALLOC(Case, int, case_len);
                            }
                            Case[i]=5111;
                            fprintf(stderr,
                   "i=%3d base=%c ipos=%.1f dpi=%4d is_dip=%d is_dye_blob=%d type=%d Case=BC%d \n",
                            peak.base_index,  data->peak_list[data_peak_ind]->base,
                            data->peak_list[data_peak_ind]->ipos, 
                            data_peak_ind, is_dip(data->peak_list,
                            data->peak_list_len, data_peak_ind, shift, options, data),
                            is_dye_blob(data->peak_list[data_peak_ind]->ipos,
                            data->peak_list[data_peak_ind], data, 0),
                            data->peak_list[data_peak_ind]->type, Case[i]);

                            data_peak_ind =
                            data->bases.called_peak_list[peak.base_index+1]->data_peak_ind;
                            fprintf(stderr,
                   "i=%3d base=%c ipos=%.1f is_dip=%d is_dye_blob=%d type=%d Case=BC%d \n",
                            peak.base_index+1,  data->peak_list[data_peak_ind]->base,
                            data->peak_list[data_peak_ind]->ipos, is_dip(data->peak_list,
                            data->peak_list_len, data_peak_ind, shift, options, data),
                            is_dye_blob(data->peak_list[data_peak_ind]->ipos,
                            data->peak_list[data_peak_ind], data, 0),
                            data->peak_list[data_peak_ind]->type, Case[i]);
                        }

                        /* Because we inserted new called peak, base index
                         * should be incremented in case next peak is uncalled 
                         */
                        base_index++;

                        continue;
                    }

                    /* BC5.1.1.3:
                     ***********
                     * Non-called peak detected at expected new location.
                     * Call it 
                     */
                    if ( BC5_CALL_PEAKS &&
                        (pind >= 0) &&
                        (peak.is_called==0) &&
                        (peak.ipos < 
                            data->bases.called_peak_list[base_index-1]->ipos) &&
                        (peak.ipos > 
                            data->bases.called_peak_list[base_index-2]->ipos)
                       )
                    {
#if 0
                        fprintf(stderr, 
                        "      BC5.1.1.3: Calling Peak %d (base %c) at pos=%d ",
                           peak.data_peak_ind, peak.base, peak.ipos);
                        fprintf(stderr,
                        "between bases %d and %d\n", base_index-2, base_index-1);
#endif
                        if (i >= case_len)
                        {
                            case_len *= 2;
                            Case = REALLOC(Case, int, case_len);
                        }

                        Case[i]=5113;
                        call_peak(data, base_index-1, peak.data_peak_ind, Case[i],
                            message); 

                        if (SHOW_CASE && (base_index < data->bases.length)) {
                            int next_dpi = 
                                data->bases.called_peak_list[base_index]->data_peak_ind;
                            fprintf(stderr,
                            "i=%3d base=%c ipos=%.1f dpi=%4d is_dip=%d is_dye_blob=%d type=%d Case=BC%d \n",
                            base_index-1,  data->peak_list[peak.data_peak_ind]->base,
                            data->peak_list[peak.data_peak_ind]->ipos, 
                            peak.data_peak_ind, is_dip(data->peak_list,
                            data->peak_list_len, peak.data_peak_ind, shift, options, data),
                            is_dye_blob(data->peak_list[data_peak_ind]->ipos,
                            data->peak_list[peak.data_peak_ind], data, 0),
                            data->peak_list[peak.data_peak_ind]->type, Case[i]);
                            fprintf(stderr,
                            "i=%3d base=%c ipos=%.1f dpi=%4d is_dip=%d is_dye_blob=%d type=%d \n",
                            base_index,  data->peak_list[next_dpi]->base,
                            data->peak_list[next_dpi]->ipos, 
                            next_dpi, is_dip(data->peak_list,
                            data->peak_list_len, next_dpi, shift, options, data),
                            is_dye_blob(data->peak_list[next_dpi]->ipos,
                            data->peak_list[next_dpi], data, 0),
                            data->peak_list[next_dpi]->type);
                        }

                        /* Because we inserted new called peak, base index
                         * should be incremented in case next peak is uncalled 
                         */
                        base_index++;
                        continue;
                    }
                }

                /* BC5.1.2:
                 **********
                 * Attempt to delete inserted base
                 * Consider 3 called bases to the right from current, 
                 * select the worst one based on the spacing requirements.
                 * If the requirements are not met, delete the base
                 */
                else if ((BC5_MERGE_PEAKS || BC5_UNCALL_PEAKS) &&
                         (options->lut_type == MegaBACE) &&
                         (detect_spacing_problems > 0)
                        ) 
                {
                    int worst_base = -1;
                    int count_removed_bases = 0;

                    if (i >= case_len)
                    {
                        case_len *= 2;
                        Case = REALLOC(Case, int, case_len);
                    }
                    Case[i]=512;             
                    while ((count_removed_bases < 5) &&
                           (worst_base = get_worst_inserted_base(base_index-1,
                                data, curr_spacing, options)) >= 0)
                    {
                        /* Delete the worst base */
                        Peak *cp = data->bases.called_peak_list[worst_base  ];
                        Peak *cpm= data->bases.called_peak_list[worst_base-1];
                        Peak *cpp= data->bases.called_peak_list[worst_base+1];
               
                        count_removed_bases++; 

                        /* Should we merge two peaks or just uncall a peak? */
                        if ( BC5_MERGE_PEAKS &&
                             (cpm->type%10==3) && (cp->type > 30)  && 
                             (cpm->base        ==  cp->base)       &&
                             (cpm->cd_peak_ind ==  cp->cd_peak_ind)   ) 
                        {
                            /* Merge with previous called peak */ 
#if SHOW_CASE
                            fprintf(stderr, 
                        "   BC5.1.2.1: Merging peaks corr. to bases %d and %d; i=%d; ", 
                            worst_base-1, worst_base, i);

                            fprintf(stderr,
                            "merged peaks %d and %d\n", 
                            data->bases.called_peak_list[worst_base-1]->data_peak_ind,
                            data->bases.called_peak_list[worst_base  ]->data_peak_ind);
#endif
                            merge_two_called_peaks(worst_base-1, worst_base, data,
                                message); 
                            continue;
                        }
                        else if ( BC5_MERGE_PEAKS &&
                                 (cp->type%10==3) && (cpp->type > 30) &&
                                 (cp->base        ==  cpp->base)      &&
                                 (cp->cd_peak_ind ==  cpp->cd_peak_ind)) 
                        {
                            /* Merge with next called peak */ 
#if SHOW_CASE
                            fprintf(stderr, 
                       "    BC5.1.2.2: Merging peaks corr. to bases %d and %d; i=%d; \n", 
                            worst_base, worst_base+1, i);

                            fprintf(stderr,
                            "merged peaks %d and %d\n", 
                            data->bases.called_peak_list[worst_base  ]->data_peak_ind,
                            data->bases.called_peak_list[worst_base+1]->data_peak_ind);
#endif
                            merge_two_called_peaks(worst_base, worst_base+1, data,
                                message);
                            continue;
                        }
                        else if (BC5_UNCALL_PEAKS && (cp->type != 11))
                        {
                            /* Just uncall peak */
#if SHOW_CASE
                            fprintf(stderr,
                     "      BC5.1.2.3: Uncalling peak corr. to base %c/%d\n", 
                            data->bases.called_peak_list[worst_base]->base,
                            worst_base);
#endif
                            (void)uncall_peak(worst_base, data, message);
                             continue;
                        }
                    }
                }
            }
 
            if (base_index < 0) {
                continue;
            }

            /* BC5.2: Uncalled non-DIP or dye blob:
             **************************************
             * do nothinng
             */
            if ((data->peak_list[i]->is_called==0) &&
                (!is_dip(data->peak_list, data->peak_list_len, i, 
                 shift, options, data) || is_dye_blob(data->peak_list[i]->ipos,
                 data->peak_list[i], data, 0))) 
            {
                if (SHOW_CASE) 
                {
                    if (i >= case_len)
                    {
                        case_len *= 2;
                        Case = REALLOC(Case, int, case_len);
                    }
                    Case[i]=52;
                    fprintf(stderr,
              "     base=%c ipos=%.1f dpi=%4d is_dip=%d is_dye_blob=%d type=%d Case=BC%d \n",
                    data->peak_list[i]->base, data->peak_list[i]->ipos, 
                    i, is_dip(data->peak_list,
                    data->peak_list_len, i, shift, options, data), 
                    is_dye_blob(data->peak_list[i]->ipos,
                    data->peak_list[i], data, 0),
                    data->peak_list[i]->type, Case[i]);
                }
            }
 
            /* Ignore peaks that were shifted outside of data range */ 
            if (data->peak_list[i]->ipos > 
                data->color_data[data->peak_list[i]->color_index].length)
                continue;

            /* Uncalled DIP or called non-DIP: 
             * get prepared 
             ********************************
             */
            if (((data->peak_list[i]->is_called==0) &&
                  is_dip(data->peak_list, data->peak_list_len, i,
                      shift, options, data)                  &&
                  !is_dye_blob(data->peak_list[i]->ipos, 
                   data->peak_list[i], data, 0))
                                                         ||
                ((data->peak_list[i]->is_called > 0) &&
                 (!is_dip(data->peak_list, data->peak_list_len, i,
                  shift, options, data) || is_dye_blob(data->peak_list[i]->ipos, 
                   data->peak_list[i], data, 0))
                )
               )
            {
                 if (base_index > 1) {
                         ave_spacing = curr_spacing;
                 }
                 else {
                     ave_spacing = INF;      /* Will not insert after 1st base */
                 }
            }

            /* BC5.3: Uncalled DIP
             *********************
             * if there is room, call the peak and insert new base
             */
            if ( BC5_UNCALLED_DIP &&
                (base_index > 1)                                       &&
                (base_index < data->bases.length-1)                    &&
                (check_spacing(
                    data->bases.called_peak_list[base_index-1]->data_peak_ind, 
                    data, curr_spacing) < 0)                           &&
                (can_insert_base(data, base_index-1, base_index  , 1,
                 ave_spacing, BIF, options) >= 3) 
               )
            {
                if ((data->peak_list[i]->is_called==0)                 && 
                    is_dip(data->peak_list, data->peak_list_len, i,  
                        shift, options, data)                                 &&
                    !is_dye_blob(data->peak_list[i]->ipos, 
                   data->peak_list[i], data, 0) &&

                    /* Make sure the inserted base location does 
                     * not coincide with location of another, already
                     * called base 
                     */
                     (abs(data->peak_list[i]->ipos -
                      data->bases.called_peak_list[base_index]->ipos)
                      >= 2)  
                                                                        &&
                     (abs(data->peak_list[i]->ipos -
                      data->bases.called_peak_list[base_index-1]->ipos)
                      >= 2)
                    )
                {
                    if (i >= case_len)
                    {
                        case_len *= 2;
                        Case = REALLOC(Case, int, case_len);
                    }
                    Case[i]=53;
                    if (SHOW_CASE && (base_index<data->bases.length))
                    {
                        fprintf(stderr, 
                   "i=%3d base=%c ipos=%.1f dpi=%4d is_dip=%d is_dye_blob=%d type=%d Case=BC%d \n",
                        base_index,  data->peak_list[i]->base,      
                        data->peak_list[i]->ipos, i, is_dip(data->peak_list,
                        data->peak_list_len, i, shift, options, data), 
                        is_dye_blob(data->peak_list[i]->ipos, 
                        data->peak_list[i], data, 0),
                        data->peak_list[i]->type, Case[i]);
        
                        fprintf(stderr, 
                "       enough_room=%d enough_room_abi=%d data_peak_ind=%d",
                          ((base_index==0) || (base_index==data->bases.length) ||
                        can_insert_base(data, base_index-1, base_index, 1,
                          ave_spacing, -1., options)) ? 1:0,
                           ((base_index==0) || (base_index==data->bases.length) ||
                        can_insert_base(data, base_index-1, base_index, 1,
                          ave_spacing_abi, -1., options)) ? 1:0, i);
                        fprintf(stderr,
                          " peak_res=%6.3f\n", 
                            data->bases.called_peak_list[base_index]->resolution);
                    } /* end SHOW_CASE */

                    if (call_peak(data, base_index, i, Case[i], message) 
                        != SUCCESS) 
                    {
                        fprintf(stderr, "Error insering new base\n");
                        FREE(orig_bases);
                        FREE(Case);
                        return ERROR;
                    }

                    base_index++; /* since the current index is now called */
                }

                continue;
 
            } /* uncalled dip */


            /* BC5.4: Called non-DIP
             ***********************
             * uncall peak, delete called base
             */
            if (!BC5_CALLED_NON_DIP &&
                (data->peak_list[i]->is_called > 0)      &&
                !is_dip(data->peak_list, data->peak_list_len,
                i, shift, options, data)                        &&
                (data->peak_list[i]->type != 11))
            {
                 if (SHOW_CASE && (base_index-1 < data->bases.length-1)) {
                    fprintf(stderr,
                    "i=%3d base=%c ipos=%.1f dpi=%4d is_dip=%d is_dye_blob=%d ",
                    base_index-1, data->peak_list[i]->base,
                    data->bases.called_peak_list[base_index-1]->ipos, i,
                    is_dip(data->peak_list, data->peak_list_len,
                    i, shift, options, data),
                    is_dye_blob(data->peak_list[i]->ipos,
                    data->peak_list[i], data, 0));
                    fprintf(stderr,
                    "        Case=%d\n", data->peak_list[i]->is_called);
                }
            }

            if (BC5_CALLED_NON_DIP &&
                (data->peak_list[i]->is_called > 0)      &&
                !is_dip(data->peak_list, data->peak_list_len,
                i, shift, options, data)                        &&
                (data->peak_list[i]->type != 11))
            {
                int enough_room = 0, curr_base_index=base_index-1;
                /* NOTE: the base index of the current called peak
                 *       equals base_index-1 !!!!!
                 */
                if ((curr_base_index-1>=0) &&
                    (curr_base_index+1<=data->bases.length-1))
                    enough_room =
                        (can_delete_base(data, curr_base_index-1, curr_base_index+1, 
                            curr_spacing, getBDF(curr_base_index), 
                            options)>=3) ? 0 : 1;

                if (SHOW_CASE && (curr_base_index < data->bases.length-1)) {
                    int p;
                    if (i >= case_len)
                    {
                        case_len *= 2;
                        Case = REALLOC(Case, int, case_len);
                    }
                    Case[i]=54;
                    fprintf(stderr,
                    "i=%3d base=%c coord=%d ipos=%.1f dpi=%4d is_dip=%d is_dye_blob=%d\n",
                    curr_base_index, data->peak_list[i]->base,
                    data->bases.coordinate[curr_base_index],
                    data->bases.called_peak_list[curr_base_index]->ipos, i,
                    is_dip(data->peak_list, data->peak_list_len,
                    i, shift, options, data),
                    is_dye_blob(data->peak_list[i]->ipos,
                    data->peak_list[i], data, 0));

                    fprintf(stderr,
                    "       enough_room=%d enough_room_abi=%d Case=BC%d type=%d\n",
                    enough_room, can_insert_base(data, curr_base_index-1,
                    curr_base_index+1, 1, curr_spacing, -1., options), Case[i],
                    data->peak_list[i]->type);

                    fprintf(stderr,
                    "       is_truncated=%d multiplicity=%d, num_poor_res=%d\n",
                    data->peak_list[i]->is_truncated,
                    get_dye_blob_multiplicity(data->peak_list[i], 
                       &(data->color_data[data->peak_list[i]->color_index]), 0, 
                       &p), p);
                }

                if ((curr_base_index-1>=0) &&
                    (curr_base_index+1<=data->bases.length-1) &&
                     !enough_room)
                {
                    data->peak_list[i]->is_called = 0;
                    if (uncall_peak(curr_base_index, data, message) != SUCCESS) {
                        fprintf(stderr, "Error deleting called base\n");
                        FREE(orig_bases);
                        FREE(Case);
                        return ERROR;
                    }
                }
                else if ((curr_base_index == 0) ||
                         (curr_base_index >= data->bases.length-1) ||
                          !can_delete_base(data, curr_base_index-1,
                              curr_base_index+1, curr_spacing, 
                              getBDF(curr_base_index), options))
                {
                    /* Try to substitute base */
                    int l, jcl, peak_indl;
                    peak_indl = -1;
                    jcl = -1;
                    location = data->peak_list[i]->ipos;
                    for (l=0; j<NUM_COLORS; l++) {
                        if (l == peak.color_index)
                            continue;
                        if (colordata_find_peak_index_by_location(&(
                            data->color_data[l]),
                            location, &new_peak, &pind, message)==ERROR) {
                            fprintf(stderr,
                            "Error calling find_peak_index_by_location\n");
                            FREE(orig_bases);
                            FREE(Case);
                            return ERROR;
                         }
                         if ((pind > 0) &&
                             is_dip(data->peak_list, data->peak_list_len,
                                 new_peak.data_peak_ind,
                                 shift, options, data) &&
                             !is_dye_blob(new_peak.ipos, &new_peak, data, 0)) {
                             peak_indl = pind;
                             jcl = l;
                         }
                    }
                    if (peak_indl > 0) {
                        data->bases.bases[curr_base_index] =
                            data->color_data[jcl].peak_list[peak_indl].base;
                        data->peak_list[i]->is_called = 0;
                        data->peak_list[new_peak.data_peak_ind]->is_called = 1;
                        data->bases.called_peak_list[curr_base_index] =
                            data->peak_list[new_peak.data_peak_ind];
                    }
                }
                continue;
            }

            /* BC5.5: Two called bases at the same position
             ************************************************
             * remove the base corresponding to lower peak 
             */
            if ((data->peak_list[i]->is_called > 0)         &&
                (base_index > 1)                            &&
                (base_index < data->bases.length)           &&
                (can_delete_base(data, base_index-2, 
                    base_index, curr_spacing, 
                    getBDF(base_index), options)>=4)        &&

                 /* Check previous called base position */
                (((abs(data->peak_list[i]->ipos -
                   data->bases.called_peak_list[base_index-2]->ipos) <=
                   MIN_DISTANCE_BETWEEN_BASES) 
                                                            &&
                  is_better_peak(data->bases.called_peak_list[base_index-2],
                      data->peak_list[i], -1, data, options))
                                                                 ||
                  /* Check next called base position */
                  ((abs(data->peak_list[i]->ipos - 
                     data->bases.called_peak_list[base_index]->ipos) <=
                     MIN_DISTANCE_BETWEEN_BASES)            &&
                   is_better_peak(data->bases.called_peak_list[base_index],
                      data->peak_list[i], -1, data, options)) 
                )
               ) 
            {
                if (SHOW_CASE) {
                    if (i >= case_len)
                    {
                        case_len *= 2;
                        Case = REALLOC(Case, int, case_len);
                    }
                    Case[i]=55;
                    fprintf(stderr,
                    "i=%3d base=%c coord=%d ipos=%.1f removed Case=BC%d \n",
                    base_index-1, data->peak_list[i]->base,
                    data->bases.coordinate[base_index],
                    data->bases.called_peak_list[base_index]->ipos,
                    Case[i]);
                }
                data->peak_list[i]->is_called = 0;
                if (uncall_peak(base_index-1, data, message) != SUCCESS) 
                {
                    fprintf(stderr, "Error deleting called base\n");
                    FREE(orig_bases);
                    FREE(Case);
                    return ERROR;
                }
            }
        } /* loop in i through all peaks */ 

        if (REORDER_BASES) (void)bc_reorder_called_bases_and_peaks(data, message);

        for (r=0; r<NUM_COLORS; r++) shift[r] = 0;

        bc_data_create_single_ordered_peak_list(data, shift, message);

        if (OUTPUT_SPACING_CURVE) output_spacing_curve(data, options);
        if (CHECK_CD_PEAK_IND)  check_cd_peak_ind(data);
        if (CHECK_CALLED_PEAKS) check_called_peaks(data);
        if (CHECK_BASE_INDEX)   check_base_index(data);
        if (CHECK_DATA_BASE_INDEX) check_data_base_index(data);
        if (CHECK_REORDERING)   check_reordering(data);
    } /* if INSERT_BASES */

    if (SHOW_CHANGED_BASES)
    {
        for (i=0; i<data->bases.length; i++)
        {
            if (orig_bases[i] != data->bases.bases[i])
            {
                fprintf(stderr, "Replaced base %d/%c with %d/%c in Case=BC%d\n",
                    i, orig_bases[i], i, data->bases.bases[i], Case[i]); 
            }
        }
    }

    if (options->time) {
        curr_clock = clock();
        fprintf(stderr, "   New bases interted in in %f sec. \n",
            (double)(curr_clock - start_clock)/(double)CLOCKS_PER_SEC);
        start_clock = curr_clock;
    }

#if SHOW_REVIEW
    if (data_review_peak_list(data, message) != SUCCESS) {
        sprintf(message->text, "Error reviewing peak list\n");
        FREE(orig_bases);
        FREE(Case);
        return ERROR;
    }
#endif
    FREE(orig_bases);
    FREE(Case);
    return SUCCESS;
}
