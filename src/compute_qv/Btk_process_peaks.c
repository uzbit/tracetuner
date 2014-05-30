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
 *
 * 2.108 2003/11/06 18:18:39
 */

#include <stdio.h> 
#include <stdlib.h> 
#include <ctype.h> 
#include <math.h> 
#include <string.h> 
#include <assert.h> 
#include <time.h>
#include <float.h>

#include "Btk_qv.h" 
#include "nr.h"
#include "util.h" 
#include "Btk_qv_data.h" 
#include "Btk_qv_funs.h"
#include "train.h"              /* needs util.h */
#include "Btk_process_peaks.h"  /* needs train.h */
#include "context_table.h"
#include "tracepoly.h"
#include "Btk_lookup_table.h"
#include "Btk_call_bases.h" 

/* ONE_PLUS & ONE_MINUS to avoid difference between platforms !!! */ 

#define AVE_HEIGHT_FACTOR              (0.75*ONE_PLUS)
#define DEFAULT_ABI_SPACING            12.
#define DEFAULT_RATIO                   0.95
#define DEFAULT_WIDTH2                  8.0
#define DP_FRACTION                     0.6
#define HIGH_BOUNDARY_FACTOR           (1.0*ONE_PLUS)
#define INIT_BETA                      (3.6*ONE_PLUS)
#define INIT_ORIG_WIDTH                (3.0*ONE_PLUS)
#define LOW_BOUNDARY_FACTOR            (0.2*ONE_PLUS)
#define MAX_CONSENSUS_LEN          200000  
#define MAX_ITER                        5 
#define MAX_RELATIVE_WIDTH             (1.9*ONE_PLUS) /* threshold width for splitting */
#define MEDIAN_STACK_SIZE              11
#define MIN_DIST_BETWEEN_BOUNDARIES     3
#define MIN_PEAK_IHEIGHT                1.0
#define MIN_PEAK_RESOLUTION            (0.1*ONE_MINUS)
#define MIN_PEAK_WIDTH_FACTOR          (0.3)
#define MULTIPLICITY_THRESHOLD         (0.2*ONE_PLUS)
#define NUM_ABI_AVE_SPACING            20
#define OUTPUT_DETECTED_PEAK_LIST       0
#define PEAK_BOUND_SIGNAL_FRACTION      0.05
#define POS_START_SPLIT                 0
#define REPRESENTATIVE_HEIGHT_FRACTION (0.25*ONE_MINUS)
#define SHOW_MARKED                     0 
#define SHOW_MULTIPLE_RESOLUTION        0
#define SHOW_PEAK_MULTIPLICITY          0
#define SHOW_PEAK_CHARACTERISTICS       0
#define SHOW_SPLITTING                  0 
#define SHOW_REDEFINE_CASE              0
#define SPLIT_BIAS                      0.000001 
#define SPLIT_MAX                      10  
#define TRUNCATED_HEIGHT                1600
#define TRUNCATED_PEAK_RATIO            0.9999
#define USE_REPRESENTATIVE_WIDTH        0
#define WIDTH2_THRESHOLD               (0.2*ONE_MINUS)
#define TWO                             2.0

#define SQR(A)          ((A)*(A))

static const double EPSILON    =       (1.0e-6*ONE_PLUS);
static const double M_1_SQRTPI =        0.56418958;    /* 1.0 / sqrt(pi) */


/*******************************************************************************
 * The following several routines are an implementation of a simple stack
 * which knows how to find it's median value.
 *******************************************************************************
 */
typedef struct {
    double stack[MEDIAN_STACK_SIZE], temp[MEDIAN_STACK_SIZE];
    int size;
} MedianStack;

/*******************************************************************************
 * Functions: print_xxx
 * Purpose:  various debug and i/o functions
 *******************************************************************************
 */ 

static void
print_one_peak_header( FILE *fp )
{
    fprintf( fp ,
".  ic   ip ty    beg  ibeg    pos  ipos    end  iend    hite   ihite   reso\n"
            );
}

static void
print_one_peak( FILE *fp, const Peak *pk )
{
    fprintf( fp,"%c %3d %4d %2d  %5d %5d  %5d %4.2f %5d %5d   %5d %7.1f %6.3f\n",
             pk->base, pk->cd_peak_ind, pk->data_peak_ind, pk->type, 
             pk->beg, pk->ibeg, pk->pos, pk->ipos, pk->end, pk->iend,
             pk->height, pk->iheight, pk->resolution );
}

static void
print_multiple_peaks( FILE *fp, const Peak *peak, int n )
{
    int j;
    for( j=0; j<n; j++ ) {
        print_one_peak( fp, &(peak[j]) );
    }
}

static int
sane_peak_pos( Peak *pk )
{
    return ( pk->beg<pk->ipos && pk->ipos<pk->end );
}

static int
is_broken_peak_sub_list( FILE *fp, const Data *data, int color, 
                         int indx, int nn, BtkMessage *message )
{
    int j;
    const ColorData *cd = &data->color_data[color];

    for (j=indx; j<indx+nn; j++) {
        if( !sane_peak_pos(&(cd->peak_list[j])) ) { return 1; }
    }

    for (j=indx+1; j<indx+nn; j++) {
        if( cd->peak_list[j-1].end > cd->peak_list[j].beg ) { return 1; }
    }
    return 0;
}

static void
print_broken_peak_sub_list( FILE *fp, const Data *data, int color, 
                            int indx, int nn, BtkMessage *message)
{
    int j;
    const ColorData *cd = &data->color_data[color];
    
    fprintf( fp, " " );
    print_one_peak_header( fp ); 

    for (j=indx; j<indx+nn; j++) {
        Peak *pk = &(cd->peak_list[j]);
        if( !sane_peak_pos(pk) ) {
                fprintf( fp, " " );
                print_one_peak( fp, pk );
            }
        }

    for (j=indx+1; j<indx+nn; j++) {
        Peak *pk0 = &(cd->peak_list[j-1]), *pk1 = &(cd->peak_list[j]);
        if( pk0->end > pk1->beg ) {
            fprintf( fp, "0" ); 
            print_one_peak( fp, pk0 );
            fprintf( fp, "1" ); 
            print_one_peak( fp, pk1 );
        }
    }
}

static int
test_broken_peak_lists( FILE *fp, const Data *data, BtkMessage *message )
{
    int color, ret_broken=0;
    for( color=0; color<NUM_COLORS; color++ ) {
        int broken = is_broken_peak_sub_list( fp, data, color, 0, 
                          data->color_data[color].peak_list_len, message );
        if( broken ) {
            ret_broken = 1;
            print_broken_peak_sub_list( fp, data, color, 0, 
                          data->color_data[color].peak_list_len, message );
        }
    }
    return ret_broken;
}

/*******************************************************************************
 * Function: get_peak_area 
 * Purpose: Calculate the area under a peak.  The algorithm is simple: 
 *          the area is defined to be the total of the averages of the heights 
 *          between each pair of points that make up the peak. 
 *******************************************************************************
 */ 
double 
get_peak_area(int *data, int peak_beg, int peak_end, 
                BtkMessage *message) 
{ 
    int i, area;
    double peak_area; 
 
    if ((peak_beg <0) || (peak_end < 0) || (peak_end < peak_beg)) { 
        sprintf(message->text,  
        "Incorrect peak boundaries (%d, %d) passed to get_peak_area\n", 
        peak_beg, peak_end); 
        return ERROR; 
    } 
         
    if (peak_beg == peak_end) { 
        return (0.0); 
    } 
 
    if ((peak_end - peak_beg) == 1) { 
        return ((double)(data[peak_end] + data[peak_beg]) / TWO); 
    } 
 
    area = data[peak_beg] + data[peak_end]; 
    for (i = peak_beg + 1; i < peak_end; i++) { 
        area += 2 * data[i]; 
    } 
    peak_area = (double)area / TWO; 
    return (peak_area); 
} 
 
/*******************************************************************************
 * Function: get_peak_area_double 
 * Purpose: Same as get_peak_area, but with double, rather than integer, 
 *          data array. 
 *          Calculate the area under a peak.  The algorithm is simple: 
 *          the area is defined to be the total of the averages of the heights 
 *          between each pair of points that make up the peak. 
 ****************************************************************************** 
 */ 
double 
get_peak_area_double(double *data, int peak_beg, int peak_end, 
    BtkMessage *message) 
{ 
    int i;
    double area, peak_area; 
 
    if ((peak_beg <0) || (peak_end < 0)) { 
        sprintf(message->text, 
        "Incorrect peak boundaries (%d, %d) passed to get_peak_area\n", 
        peak_beg, peak_end); 
        return ERROR; 
    } 
        
    if (peak_beg == peak_end) { 
        return (0.0); 
    } 
 
    if ((peak_end - peak_beg) == 1) { 
        return ((double)(data[peak_end] + data[peak_beg]) / TWO ); 
    } 
     
    area = data[peak_beg] + data[peak_end]; 
    for (i = peak_beg + 1; i < peak_end; i++) { 
        area += 2 * data[i]; 
    } 
    peak_area = (double)area / TWO; 
    return (peak_area); 
} 
 
/****************************************************************************** 
 * Function: get_peak_position 
 * Purpose: Find location of observed peak as the position which bisects  
 *          its area. If the position is found, return it. Otherwise, return -1 
 ****************************************************************************** 
 */ 
int  
get_peak_position(int *data, int peak_beg, int peak_end, 
			    double peak_area, BtkMessage *message) 
{ 
    int j; 
    double area, old_area; 
 
    /* Check the easy (degenerate) case first */ 
    if (peak_end - peak_beg <= 1) { 
	return peak_end;        
    } 
 
    area = 0.0; 
    old_area = 0.0; 
    if ((peak_end - peak_beg) == 2) { 
       return peak_beg+1;       
    } 
    else { 
        for (j = peak_beg + 1; j <= peak_end; j++) { 
            area += (double)(data[j] + data[j-1]) / TWO; 
            /*if (area == peak_area / TWO) { */
            if( DBL_EQ_DBL( area, peak_area/TWO ) ) {
                return (j);          
            } 
            else if(  (area > ONE_PLUS * peak_area / 2.0) && 
                      (old_area <= ONE_PLUS * peak_area / 2.0)  ) {
                if ((area - peak_area / 2.0) < 
                    ONE_PLUS * (peak_area / 2.0 - old_area)) { 
                    return (j);
                } 
                else { 
                    return (j - 1); 
                } 
            } 
            old_area = area; 
        } 
    } 
 
    (void)sprintf(message->text, 
        "couldn't determine position of peak @ %d-%d", peak_beg, peak_end); 
    return ERROR; 
} 

/******************************************************************************
 * Function: get_peak_max
 * Purpose:  Find the position inside peak where signal has the highest value.
 ******************************************************************************
 */
int
get_peak_max(int *data, int peak_beg, int peak_end, BtkMessage *message)
{
    int j, max=-1, signal=-1;

    for (j=peak_beg; j<=peak_end; j++) {
        if (signal < data[j]) {
            signal = data[j];
            max = j;
        }
    } 
    return max;
}

#if 1 

/*******************************************************************************
 * Function: is_maximum 
 * Purpose: determine if the current element of data array is its internal 
 *          local maximum; return 1 if it is or 0 otherwise 
 *******************************************************************************
 */ 
int 
is_maximum(int n, int *data, int length) 
{ 
    if ((n > 1) && (n <length-2)) { 
        if      ((data[n]   > data[n-1]) && 
                 (data[n-1] > data[n-2]) && 
                 (data[n  ] > data[n+1]) && 
                 (data[n+1] > data[n+2])) 
        { 
            return 1; 
        } 
        else if ((data[n] > data[n-1]) && 
                 (data[n] ==data[n+1]) && 
                 (data[n] > data[n+2])) 
        { 
            return 1; 
        } 
        else if ((data[n] > data[n+1]) && 
                 (data[n] ==data[n-1]) && 
                 (data[n] > data[n-2])) 
        { 
            return 1; 
        } 
        else if ((data[n] == data[n-1]) && 
                 (data[n] == data[n+1]) && 
                 (data[n] >  data[n-2]) && 
                 (data[n] >  data[n+2])) 
        { 
            return 1; 
        } 
    } 
    return 0; 
} 

#endif
 
/*******************************************************************************
 * Function: get_default_peak 
 * Purpose:  create a peak with default characteristics 
 *******************************************************************************
 */ 
static int 
get_default_peak(ColorData *cd, Peak *peak) 
{ 
    peak->ibeg = 0; 
    peak->iend = 0; 
    peak->is_truncated = 0; 
    peak->height = 0; 
    peak->iheight = 0; 
    peak->base = cd->base; 
    peak->base_index = -1; 
    peak->color_index = cd->dye_number-1; 
    peak->data_peak_ind = 0; 
    peak->cd_peak_ind = 0; 
    peak->beg = 0; 
    peak->end = 0; 
    peak->is_called = 0;                     /* default */ 
    peak->spacing = 0; 
    peak->type = 0; 
    peak->width1 = 0.; 
    peak->width2 = 0.; 
    peak->C0 = 0.;
    peak->beta = 4.; 
    peak->orig_width = 0.; 
    peak->ave_w02beta = 0.01;   /* added by SLT */
    peak->resolution = -1; 
 
    return SUCCESS; 
} 

/*******************************************************************************
 * Function: is_minimum
 * Purpose: determine if the current element of data array is its internal
 *          local minimum; return 1 if it is or 0 otherwise
 *******************************************************************************
 */
static int
is_minimum(int n, int *data, int length)
{
    if ((n < 1) || (n > length-2)) 
        return 0;

    if ((data[n] < data[n-1]) && (data[n] < data[n+1]))
        return 1;

    if ((n > 1) && (n < length-2) && 
        (data[n  ] <= data[n-1]) && 
	(data[n-1] <  data[n-2]) &&
        (data[n  ] <= data[n+1]) && 
	(data[n+1] <  data[n+2]))
    {
        return 1;
    }

    return 0;
}

/*******************************************************************************
 * Function: is_big_peak
 * Purpose:  check if peak is "big" (this idea was suggested by David Ho)
 *******************************************************************************
 */
static int
is_big_peak(Peak peak, double ave_area, double prev_peak_area, int max_value)
{
    /* Too high peak */
    if      (peak.is_truncated)
        return 1;

    /* Peak of too big area */
    else if ((peak.area > ave_area       / AREA_FACTOR2) ||
             (peak.area > prev_peak_area / AREA_FACTOR1)   )
    {
        return 1;
    }

    return 0;
}

/*******************************************************************************
 * Function: is_dp
 * Purpose:  check if a given peak is a dominant peak, that is,
 *           the observed peak signal at its position is greater than
 *           observed signal of any other color at this position
 *******************************************************************************
 */
int
is_dp(Peak *pk, int *shift, Data *data)
{
    int    j, height = pk->height,
           color  = pk->color_index,
           pos    = pk->pos;
    int    is_dp = 1;
    float  other_signal = -1., dominance_factor = DP_FRACTION;

    /* Determine if current peak is DP */
    for (j=0; j<NUM_COLORS; j++) {
        int new_pos = pos + shift[color] - shift[j];

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
 * Function: colordata_detect_peaks
 * Purpose: Create a list of peaks of a given color 
 *******************************************************************************
 */ 
static int 
colordata_detect_peaks(Data *data, int color, Options *options, 
    BtkMessage *message) 
{ 
    ColorData *cd; 
    int i, j, k, l, derivative_2, endscan; 
    Peak peak = initialize_peak(), temp_peak = initialize_peak(); 
    double ave_peak_area=0., start_avg_peak_area, prev_peak_area; 
    double avg_peak_width=0., start_avg_peak_width, sum_peak_area, 
           sum_peak_width, average10_width, relative_width; 

    /* Initialization */
    peak.ave_sq_width2 = -1;
    peak.ave_width2 = -1;
    peak.ave_width1 = -1;
    peak.ipos_orig = -1;
    peak.wiheight = 0.;
    peak.data_peak_ind2 = -1;
    peak.max = -1;

    cd = &(data->color_data[color]); 
    cd->peak_list_len = 0; 
 
    get_default_peak(cd, &peak); 
 
    /*  Pre scan some data to get an average peak area to start; 
     * Start prescanning from position of base number 120 or, 
     * if total  number of bases is < 200, from the middle base. 
     * Continue prescanning through the position of the last base 
     * or untill 20 "good" peaks are found, whichever comes first 
     */ 
    j = 0; 
    temp_peak.ibeg = -1; 
    temp_peak.iend = -1; 
    temp_peak.is_truncated = 0;
    start_avg_peak_width = 0.0; 
    start_avg_peak_area = 0.0; 
    if (data->bases.length > 200) 
        k = 120; 
    else 
        k = (data->bases.length-1)/2; 
    l = data->bases.length-1; 
 
    /* Make sure cd->data[i+2] has sense */ 
    endscan = data->bases.coordinate[l] < cd->length-2 ? 
              data->bases.coordinate[l] : cd->length-2; 
 
    for (i = data->bases.coordinate[k]; 
        (i < endscan) && (j < 20); i++) { 
        derivative_2 = cd->data[i-1] - 2*cd->data[i] + cd->data[i+1]; 
 
        /* Temp peak begins */ 
        if ((derivative_2 < 0) && (temp_peak.ibeg < 0)) { 
            temp_peak.ibeg = i; 
        } 

        if ((temp_peak.ibeg >= 0) && (i > temp_peak.ibeg + 1) &&
            (cd->data[i  ] >= TRUNCATED_HEIGHT) &&
            (cd->data[i  ] == cd->data[i-1]) && 
            (cd->data[i-1] == cd->data[i-2]))
        {
            temp_peak.is_truncated = 1;
        }
 
        /* Temp peak ends */ 
        if ((temp_peak.ibeg >= 0) && 
        ((derivative_2 > 0 && cd->data[i]-2*cd->data[i+1]+cd->data[i+2] >= 0) ||
         is_minimum(i, cd->data, cd->length) ||
         (cd->data[i] <=0))) 
        { 
            temp_peak.iend = i; 
            if ( 
                /* Ignore oscillations of width 1 and narrow width!!! */ 
                 (temp_peak.iend-temp_peak.ibeg > 2) &&
                /* Idnore truncated peaks */
                 !temp_peak.is_truncated) 
            { 
                start_avg_peak_width += (temp_peak.iend-temp_peak.ibeg); 
                start_avg_peak_area += get_peak_area(cd->data,temp_peak.ibeg, 
                                                      temp_peak.iend,message); 
                j += 1; 
            } 
            temp_peak.ibeg = -1; 
            temp_peak.iend = -1; 
            temp_peak.is_truncated = 0;
        } 
    } 
 
    /* Initialize starting params */ 
    if (j>0) { 
        start_avg_peak_width /= j; 
        start_avg_peak_area /= j; 
    } 
    else { 
        start_avg_peak_width = 1.; 
        start_avg_peak_area  = 1.; 
    } 
    prev_peak_area = start_avg_peak_area; 
 
    /************************************************************************
     * Use same approach as Phred: first, find a couple of subsequent 
     * inflection points at one of which the 2nd derivative changes sign from 
     * + to - and at the other backwards; then, localize the peak position 
     * between the two points 
     ************************************************************************
     */ 
    peak.ibeg = -1;
    for (i = data->pos_data_beg + 1; i < data->pos_data_end - 2; i++) {
        derivative_2 = cd->data[i-1] - 2*cd->data[i] + cd->data[i+1]; 
 
        /* Peak begins */ 
        if ((cd->peak_list_len < MAX_NUM_OF_PEAK-1) &&
            (derivative_2 < 0) && (peak.ibeg < 0)) { 
            peak.ibeg = peak.max = i; 
        } 
 
        if ((peak.ibeg >= 0) && (i > peak.ibeg + 1) &&
            (cd->data[i  ] >= TRUNCATED_HEIGHT) &&
            (cd->data[i  ] == cd->data[i-1]) &&
            (cd->data[i-1] == cd->data[i-2]))
        {
            peak.is_truncated = 1;
            if (options->Verbose >= 5)
               fprintf(stderr, "Truncated peak detected at scan %d\n", i); 
        }
 
        if ((peak.ibeg >= 0) && (cd->data[peak.max] < cd->data[i])) {
            peak.max = i; 
        } 

        /* Peak ends */ 
        if ((peak.ibeg >= 0) &&
               (((derivative_2 > 0) && 
               INT_GT_DBL(i-peak.ibeg, MIN_PEAK_WIDTH_FACTOR*avg_peak_width)) ||
               is_minimum(i, cd->data, cd->length) ||
               (cd->data[i] <=0))
            )
        { 
            peak.iend = i; 
            /* Determine the peak's position between peak.ibeg and peak.iend */ 
            peak.area = get_peak_area(cd->data, peak.ibeg, peak.iend, message); 
            peak.pos = get_peak_position(cd->data, peak.ibeg, peak.iend, 
                peak.area, message); 
            peak.ipos = peak.pos; 
            if (peak.pos == ERROR) { 
                sprintf(message->text, "Error in get_peak_position\n"); 
                return ERROR; 
            } 
 
            /* Calculate average peak area and width. Exclude "too big" peaks" 
             * when calculating ave_peak_area and ave_peak_width
             */ 
            k = 0; 
            if (cd->peak_list_len < 11) { 
                ave_peak_area  = start_avg_peak_area; 
                avg_peak_width = start_avg_peak_width; 
            } 
            else 
            { 
                sum_peak_area  = 0.0; 
                sum_peak_width = 0.0; 
                for (j = cd->peak_list_len - 1; k < 10 && j > 0; j--) 
                { 
                    /* check for big peaks */ 
                    if (!is_big_peak(cd->peak_list[j], 
                        ave_peak_area, prev_peak_area, TRUNCATED_HEIGHT)) 
                    { 
                        ++k; 
                        sum_peak_area  += cd->peak_list[j].area; 
                        sum_peak_width += (cd->peak_list[j].iend - 
                                          cd->peak_list[j].ibeg); 
                    } 
                } 
                if (k < 10) { 
                   ave_peak_area  = start_avg_peak_area; 
                   avg_peak_width = start_avg_peak_width; 
                } 
                else { 
                   ave_peak_area  = sum_peak_area  / 10.0; 
                   avg_peak_width = sum_peak_width / 10.0; 
                } 
            } 

            /* Record the current peak.
             *
             * Do it only if 
             * 1) peak is wider than 2 scans
             * 2) its area is >= 10% of the average area of 10 preceeding peaks 
             * 3) its area is > 5% of the area of the
             *    immediately preceeding peak (of the same color).
             *
             * When checking these criteria, exclude "too big" peaks"
             * when calculating ave_peak_area and ave_peak_width 
             */
 
            if ( 
                /* Ignore oscillations of width 1 !!! */ 
                ((peak.iend-peak.ibeg > 2) || 
                 ((peak.iend-peak.ibeg==2) && (cd->data[peak.iend-1] > 0))) && 
 
                /* Ignore peaks of too small area, see Phred paper II, p.178 */ 
                  (peak.area >= AREA_FACTOR1 * prev_peak_area) && 
                  (peak.area >= AREA_FACTOR2 * ave_peak_area) 
               ) 
            { 
                peak.height = cd->data[peak.pos], 
                peak.iheight = (double)peak.height; 
                peak.beg = peak.ibeg; 
                peak.end = peak.iend; 
                peak.cd_peak_ind = cd->peak_list_len; 
                if (ave_peak_area  > 0.) { 
                    peak.relative_area = peak.area/ave_peak_area ; 
                } 
                else { 
                    peak.relative_area = 1.; 
                } 
                if (cd->peak_list_len >= cd->peak_list_max_len) { 
                    cd->peak_list_max_len *= 2; 
                    cd->peak_list = REALLOC(cd->peak_list, Peak, 
                    cd->peak_list_max_len); 
                    MEM_ERROR(cd->peak_list); 
                } 

                cd->peak_list[cd->peak_list_len] = peak; 
                cd->peak_list_len++; 

                /* Compute relative width for the current peak.
                 * Don't consider truncated peaks */ 
                average10_width = 0; 
                relative_width = 0; 
                {
                    int num_peaks=0;
                    for (j = QVMAX(0, cd->peak_list_len - 10); 
                         j < cd->peak_list_len; j++) 
                    {
                        if (cd->peak_list[j].is_truncated)
                            continue;

                        num_peaks++; 
                        average10_width += cd->peak_list[j].iend - 
                            cd->peak_list[j].ibeg; 
                    } 
                    num_peaks = QVMAX(num_peaks, 1);
                    average10_width /= num_peaks;
                }
                if (average10_width > 0) 
                    relative_width = (double) (peak.iend - peak.ibeg) / 
                                                    average10_width; 
                else 
                    relative_width = 1.0; 
 
                if (!is_big_peak(peak, ave_peak_area, prev_peak_area, 
                    TRUNCATED_HEIGHT)) 
                { 
                    prev_peak_area = peak.area; 
                } 
            } 
 
            /* Prepare for processing next peak */ 
            peak.ibeg = -1; 
            peak.iend = -1; 
            peak.is_truncated = 0; 
 
        } /* end processing current peak */ 
    } /* end loop in i */ 
 
#if SHOW_NUM_SPLIT_PEAKS 
    fprintf(stderr, "Number of splitted peaks=%d\n", num_splitted); 
#endif 
    return SUCCESS; 
 
error: 
    return ERROR; 
} 
 
/*******************************************************************************
 * Function: output_peak_list 
 * Purpose: write out the contents of the populated peak list to 
 *          a file for debugging purposes 
 *******************************************************************************
 */ 
static void 
output_peak_list(Data *data, BtkMessage *message) 
{ 
    int c, i; 

    FILE *fp = fopen( "Peaks", "w" );
    for (c = 0; c < NUM_COLORS; c++) { 
	for (i = 0; i < data->color_data[c].peak_list_len; i++) { 
	    (void)fprintf(fp, "peak color=%d ind=%d beg=%d pos=%d end=%d\n", 
                c, i,
	        data->color_data[c].peak_list[i].ibeg,
                data->color_data[c].peak_list[i].pos,
                data->color_data[c].peak_list[i].iend); 
	} 
    } 
    fclose(fp);
} 
 
/*******************************************************************************
 * Function: data_detect_peaks 
 * Purpose: for each of 4 colordata arrays, create a list of "observed" peaks 
 *******************************************************************************
 */ 
int 
data_detect_peaks(Data *data, Options *options, BtkMessage *message) 
{ 
    int r; 
    int color; 
 
    for (color = 0; color < NUM_COLORS; color++) { 
        if ((r = colordata_detect_peaks(data, color, options, message)) == ERROR) 
    { 
            return ERROR; 
        }
    } 
 
    if (OUTPUT_DETECTED_PEAK_LIST) {
        output_peak_list(data, message); 
    } 
 
    return SUCCESS; 
} 

 
/*******************************************************************************
 * Function: is_true_peak                    
 * Purpose: determine if the ratio of peak area (in the new definition) to its 
 *          height is close to the peak's width2, where width2 is the average   
 *          ratio of peak area to peak height among 10 "true" immediately  
 *          preceeding peaks 
 *******************************************************************************
 */ 
  
int 
is_true_peak(Peak peak) 
{ 
    if (!peak.is_truncated && 
        (peak.resolution >=0)  &&  
        (peak.resolution <= MIN_PEAK_RESOLUTION)) 
    { 
        return 1; 
    } 
    else { 
        return 0; 
    } 
} 
 
/*******************************************************************************
 * Function: get_left_half_width1 
 * Purpose: compute the peak width at a half of its apparent height   
 *******************************************************************************
 */ 
static double 
get_left_half_width1(int *data, int peak_beg, int peak_end, int peak_pos, 
    int peak_height, BtkMessage *message) 
{ 
    int    l; 
    double left, height2; 
 
    if (peak_end - peak_beg < 4) { 
        return (double)(peak_end - peak_beg)/2.; 
    } 
/*  if (peak_pos - peak_beg <= 3) { 
 *      return (peak_pos-1); 
 *  } 
 */ 
 
    l = peak_beg; 
    left = (double)peak_beg;   /* default value */ 
 
    height2 = (double)data[peak_pos]/2.; 
    
    if ( INT_LT_DBL( data[peak_beg], height2) ) { 
        while ( INT_LT_DBL( data[l], height2) ) { 
            //if ( INT_LT_DBL( data[peak_beg],height2) ) { 
            //    while ( INT_LT_DBL( data[l], height2) ) { 
            l++; 
        } 
        if ( INT_EQ_DBL(data[l],height2) ) { 
            left = (double)l; 
        } 
        else { 
            left = (double)(l-1)+ 
                   (double)(height2 - data[l-1])/ 
                   (double)(data[l] - data[l-1]);  
        } 
    } 
   
    if ((double)peak_pos - left < 0) { 
       fprintf(stderr, "Error in get_left_half_width1: \n"); 
       fprintf(stderr, "Peak_beg=%d peak_pos=%d peak_end=%d left_pos=%f\n", 
          peak_beg, peak_pos, peak_end, left  ); 
       return ERROR; 
    } 
  
    return ((double)peak_pos - left); 
} 
 
/*******************************************************************************
 * Function: get_right_half_width1 
 * Purpose: compute the peak width at a half of its height 
 *******************************************************************************
 */ 
 
static double 
get_right_half_width1(int *data, int peak_beg, int peak_end, int peak_pos,  
    int peak_height, BtkMessage *message) 
{ 
    int    r; 
    double right, height2; 
  
    if (peak_end - peak_beg < 4) { 
        return (double)(peak_end - peak_beg)/2.; 
    } 
/*  if (peak_end - peak_pos <= 3) { 
 *      return (peak_pos+1); 
 *  } 
 */ 
    r = peak_end; 
    right= (double)r; 
  
    height2 = (double)data[peak_pos]/2.; 
  
    if ( INT_LT_DBL( data[peak_end], height2 ) ) { 
        while ( INT_LT_DBL( data[r], height2 ) ) { 
            r--; 
        } 
        if ( INT_EQ_DBL( data[r], height2) ) { 
            right= (double)r; 
        } 
        else { 
            right= (double)r+ 
                   (double)(data[r]-height2       )/ 
                   (double)(data[r]- data[r+1]); 
        } 
    } 
    if (right - (double)peak_pos < 0) { 
       fprintf(stderr, "Error in get_right_half_width1: \n"); 
       fprintf(stderr,"Peak_beg=%d peak_pos=%d peak_end=%d r=%d right_pos=%f\n",
          peak_beg, peak_pos, peak_end, r, right );  
       return ERROR; 
    } 
 
    return (right - peak_pos); 
} 
 
/*******************************************************************************
 * Function: get_peak_width1 
 * Purpose: compute the peak width at a half of its height 
 *******************************************************************************
 */ 
  
static double 
get_peak_width1(ColorData *cd, int peak_index, BtkMessage *message) 
{ 
    double  width1=0.0, left, right; 
 
    if ((cd->peak_list[peak_index].type != 12) && 
        (cd->peak_list[peak_index].type != 21))  
    { 
        left = get_left_half_width1(cd->data, cd->peak_list[peak_index].beg, 
                cd->peak_list[peak_index].end, cd->peak_list[peak_index].pos, 
                cd->peak_list[peak_index].height, message); 
        if (left <0) { 
            sprintf(message->text, "Left half width1 <0\n"); 
            return ERROR; 
        } 
        right= get_right_half_width1(cd->data, cd->peak_list[peak_index].beg, 
                cd->peak_list[peak_index].end, cd->peak_list[peak_index].pos, 
                cd->peak_list[peak_index].height, message); 
        if (right <0) { 
            sprintf(message->text, "Right half width1 <0\n"); 
            return ERROR; 
        } 
        width1 =  left + right; 
    } 
    else if (cd->peak_list[peak_index].type == 12) {  
        left = get_left_half_width1(cd->data, cd->peak_list[peak_index].beg, 
                cd->peak_list[peak_index].end, cd->peak_list[peak_index].pos, 
                cd->peak_list[peak_index].height, message); 
        if (left <0) { 
            sprintf(message->text, "Left half width1 <0\n"); 
            return ERROR; 
        } 
        width1 = 2.*left; 
    } 
    else if (cd->peak_list[peak_index].type == 21) { 
        right= get_right_half_width1(cd->data, cd->peak_list[peak_index].beg, 
                cd->peak_list[peak_index].end, cd->peak_list[peak_index].pos, 
                cd->peak_list[peak_index].height, message); 
        if (right <0) { 
            sprintf(message->text, "Right half width1 <0\n"); 
            return ERROR; 
        } 
        width1 = 2.*right; 
    }  
    
    return width1; 
} 
 
 
/*******************************************************************************
 * Function: get_left_half_width2 
 * Purpose: compute the peaks left half width2  
 *******************************************************************************
 */ 
static double 
get_left_half_width2(int *data, int peak_beg, int peak_end, int peak_pos, 
    int peak_height, BtkMessage *message) 
{ 
    if (peak_end - peak_beg < 4) { 
        return get_peak_area(data, peak_beg, peak_end, message) 
            /2./(double)data[peak_pos]; 
    } 
 
    if ((peak_pos - peak_beg == 1) || (peak_pos - peak_beg == 0)) { 
        return (peak_pos - peak_beg); 
    } 
 
    if (data[peak_pos] <= 0) { 
       fprintf(stderr, "Error: data[peak_pos] = %d for peak_pos = %d\n", 
          data[peak_pos], peak_pos); 
       return ERROR; 
    } 
    return get_peak_area(data, peak_beg, peak_pos, message)/data[peak_pos]; 
} 
 
/*******************************************************************************
 * Function: get_right_half_width2 
 * Purpose: compute the peaks right half width2 
 *******************************************************************************
 */ 
static double 
get_right_half_width2(int *data, int peak_beg, int peak_end, int peak_pos, 
    int peak_height, BtkMessage *message) 
{ 
    if (peak_end - peak_beg < 4) { 
        return get_peak_area(data, peak_beg, peak_end, message)/2. 
            /data[peak_pos]; 
    } 
    if ((peak_end - peak_pos == 1) || (peak_end - peak_pos == 0)) { 
        return peak_end - peak_pos; 
    } 
    if (data[peak_pos] <= 0) { 
       fprintf(stderr, "Error: data[peak_pos] = %d for peak_pos = %d\n", 
          data[peak_pos], peak_pos); 
       return ERROR; 
    } 
    return get_peak_area(data, peak_pos, peak_end, message)/ 
        (double)data[peak_pos]; 
} 
 
/*******************************************************************************
 * Function: get_peak_width2 
 * Purpose: compute the peak width at a half of its height 
 *******************************************************************************
 */ 
 
static double 
get_peak_width2(ColorData *cd, int peak_index, BtkMessage *message) 
{ 
    double  width2=0.0, left_half_width2, right_half_width2; 
 
    if (cd->peak_list[peak_index].height <=0) { 
        fprintf(stderr, "Peak with height=%d passed to ", 
            cd->peak_list[peak_index].height); 
        fprintf(stderr, "get_peak_width2 function; color=%d, peak_ind=%d\n", 
            cd->dye_number-1, peak_index); 
        return ERROR; 
    } 
    if ((cd->peak_list[peak_index].type != 12) && 
        (cd->peak_list[peak_index].type != 21)) 
    { 
        if ((left_half_width2 = get_left_half_width2(cd->data,  
             cd->peak_list[peak_index].beg, 
             cd->peak_list[peak_index].end, cd->peak_list[peak_index].pos, 
             cd->peak_list[peak_index].height, message)) < 0) { 
            fprintf(stderr, "Error computing left_half_width2 for "); 
            fprintf(stderr,  
              "color=%d peak_index=%d beg=%d pos=%d end=%d area=%f height=%d\n",
                cd->dye_number-1, peak_index, cd->peak_list[peak_index].beg, 
                cd->peak_list[peak_index].pos, cd->peak_list[peak_index].end, 
                cd->peak_list[peak_index].area,  
                cd->data[cd->peak_list[peak_index].pos]); 
            return ERROR; 
        } 
        if ((right_half_width2 = get_right_half_width2(cd->data,  
             cd->peak_list[peak_index].beg, 
             cd->peak_list[peak_index].end, cd->peak_list[peak_index].pos, 
             cd->peak_list[peak_index].height, message)) <= 0) { 
            fprintf(stderr, "Error computing right_half_width2 for "); 
            fprintf(stderr, "color=%d, peak_index=%d\n", 
                cd->dye_number-1, peak_index); 
            return ERROR; 
        } 
        width2 = left_half_width2 + left_half_width2; 
    } 
    else if (cd->peak_list[peak_index].type == 12) { 
        if ((left_half_width2 = get_left_half_width2(cd->data, 
             cd->peak_list[peak_index].beg, 
             cd->peak_list[peak_index].end, cd->peak_list[peak_index].pos, 
             cd->peak_list[peak_index].height, message)) <= 0) { 
            fprintf(stderr, "Error computing left_half_width2 for "); 
            fprintf(stderr, "color=%d, peak_index=%d\n", 
                cd->dye_number-1, peak_index); 
            return ERROR; 
        } 
        width2 = 2.* left_half_width2; 
    } 
    else if (cd->peak_list[peak_index].type == 21) { 
        if ((right_half_width2 = get_right_half_width2(cd->data, 
             cd->peak_list[peak_index].beg, 
             cd->peak_list[peak_index].end, cd->peak_list[peak_index].pos, 
             cd->peak_list[peak_index].height, message)) <= 0) { 
            fprintf(stderr, "Error computing right_half_width2 for "); 
            fprintf(stderr, "color=%d, peak_index=%d\n", 
                cd->dye_number-1, peak_index); 
            return ERROR; 
        } 
        width2 = 2.* right_half_width2; 
    } 
    return width2; 
} 
 
/*******************************************************************************
 * Function: get_average_width2 
 * Purpose: compute the average width2 among num_upstream+num_downstream 
 *          non-truncated peaks of types 11, 12 and 21 upstream and downstream 
 *          from the current peak 
 *******************************************************************************
 */ 
double 
get_average_width(int w_ind, int peak_index, int num_upstream, 
    int num_downstream, double *stderr_width, 
    Data *data, BtkMessage *message) 
{ 
    int    i, j, cd_peak_ind; 
    int shift[NUM_COLORS] = {0, 0, 0, 0};
    double ave_width, ave_sq_width, width, sum = 0., sum2 = 0.; 
    double max_width, min_width; 
    ColorData *cd; 

    if (w_ind != 1 && w_ind != 2)
    {
        fprintf(stderr, "Wrong value of w_ind= %d passed to get_average_width\n",
            w_ind);
        exit(-1);
    }

    if (peak_index < 0) { 
        (void)sprintf(message->text, 
                "Negative peak_index passed to get_average_width%d\n", w_ind); 
        fprintf(stderr, "Negative peak_index passed to get_average_width%d\n",
            w_ind); 
        return ERROR; 
    }

    cd = &(data->color_data[data->peak_list[peak_index]->color_index]); 
    if (peak_index == 0) { 
        if ((w_ind == 1) &&
            (ave_width = get_peak_width1(cd, peak_index, message)) > 0) { 
           *stderr_width = 0.;
            return ave_width; 
        } 
        else if ((w_ind == 2) &&
            (ave_width = get_peak_width2(cd, peak_index, message)) > 0) {
           *stderr_width = 0.;
            return ave_width;
        }
        else { 
            fprintf(stderr, "Error calling get_peak_width%d \n", w_ind); 
            return ERROR; 
        } 
    } 
    else 
    { 
        i = peak_index; 
        j = 0; 
        max_width = NINF; 
        min_width = INF; 
        while (i >0 && j < num_downstream) { 
            cd = &(data->color_data[data->peak_list[i]->color_index]); 

            if (!is_dp(data->peak_list[i], shift, data))
            {
                i--;
                continue;
            }

            if (!data->peak_list[i]->is_truncated && 
                (data->peak_list[i]->type == 11)) 
            { 
                cd_peak_ind = data->peak_list[i]->cd_peak_ind; 
                width = (w_ind == 1) ?
                    get_peak_width1(cd, cd_peak_ind, message) :
                    get_peak_width2(cd, cd_peak_ind, message); 
                sum  += width; 
                sum2 += width * width;
                if (width > max_width) max_width = width; 
                if (width < min_width) min_width = width; 
                j++; 
            } 
            i--; 
        } 
        i = peak_index; 
        while (i<data->peak_list_len-1 && j < num_upstream+num_downstream) { 
            cd = &(data->color_data[data->peak_list[i]->color_index]); 

            if (!is_dp(data->peak_list[i], shift, data))
            {
                i++;
                continue;
            }

            if (!data->peak_list[i]->is_truncated && 
                ((data->peak_list[i]->type == 11) 
                                             || 
                ((data->peak_list[i]->type == 12) && 
                 (cd->data[data->peak_list[i]->end] < data->peak_list[i]->height
                 * AVE_HEIGHT_FACTOR)) 
                                             || 
                ((data->peak_list[i]->type == 21) && 
                 (cd->data[data->peak_list[i]->beg] < data->peak_list[i]->height
                 * AVE_HEIGHT_FACTOR))
            )) 
            { 
                cd_peak_ind = data->peak_list[i]->cd_peak_ind; 
                width = (w_ind == 1) ?
                        get_peak_width1(cd, cd_peak_ind, message) :
                        get_peak_width2(cd, cd_peak_ind, message); 
                sum  += width; 
                sum2 += width * width;
                if (width > max_width) max_width = width; 
                if (width < min_width) min_width = width; 
                j++; 
            } 
            i++; 
        } 
        if ((sum > 0) && (j>0)) { 
            if ((j > 2) && (sum-max_width-min_width > 0)) {
                ave_width = (sum-max_width-min_width)/(double)(j-2); 
                ave_sq_width = 
                    (sum2-max_width*max_width-min_width*min_width)/(double)(j-2);
               *stderr_width = sqrt(ave_sq_width - ave_width*ave_width);
            }
            else 
            {
                ave_width    = sum /(double)j; 
                ave_sq_width = sum2/(double)j;
               *stderr_width = sqrt(ave_sq_width - ave_width*ave_width);
            }
        } 
        else { 
            cd_peak_ind = data->peak_list[peak_index]->cd_peak_ind; 
            if ((w_ind == 1) &&
                (ave_width = get_peak_width1(cd, cd_peak_ind, message)) > 0) {
                ave_sq_width = ave_width*ave_width;
               *stderr_width = 0.;
                return ave_width;
            }
            else if ((w_ind == 2) && 
                (ave_width = get_peak_width2(cd, cd_peak_ind, message)) > 0) { 
                ave_sq_width = ave_width*ave_width;
               *stderr_width = 0.;
                return ave_width; 
            } 
            else { 
                fprintf(stderr, "Error calling get_peak_width2 \n"); 
                return ERROR; 
            } 
        } 
        return ave_width; 
    } 
} 

/*******************************************************************************
 * Function: get_abi_base_index_by_position
 * Purpose:  given a position in chromatogram, determine the base index
 *           of the closest ABI base call located to the left from the position
 *******************************************************************************
 */
static int
find_abi_base_index_by_position(int given_position, Data *data, 
    BtkMessage *message)
{
    int lo_ind, mid_ind, hi_ind;

    if ((given_position < 0) || (given_position > data->color_data[0].length-1))
    {
        (void)fprintf(stderr,        
        "Position passed to find_abi_base_index %d is out of bounds [0,%d]\n",
                      given_position, data->color_data[0].length-1 );
        return (-1);    
    }

    if (data->bases.length <= 1)
        return 0;

    if ((data->bases.length>1) && 
        (given_position < data->bases.coordinate[1])) {
        return 0;
    }

    if (given_position >= data->bases.coordinate[data->bases.length-1])
        return (data->bases.length-1);

    lo_ind  = 1;
    hi_ind  = data->bases.length-1;

    do {
        mid_ind = (lo_ind + hi_ind) / 2;
        if      (given_position == data->bases.coordinate[mid_ind]) 
            return mid_ind;
        else if (given_position <  data->bases.coordinate[mid_ind])
            hi_ind = mid_ind;
        else
            lo_ind = mid_ind;  
    } while (lo_ind < hi_ind-1);

    return lo_ind;
}


/*******************************************************************************
 * Function: get_average_abi_spacing
 * Purpose:  for a given position in chromatogram, determine the average
 *           spacing between a specified number of ABI base calls
 *******************************************************************************
 */
double
get_average_abi_spacing(int position, Data *data, BtkMessage *message)
{
    int    i, j=0, k;

    i = find_abi_base_index_by_position(position, data, message);

    if (i < 0) {
        fprintf(stderr, 
        "Error in computing abi_base_index_by_position: i=%d\n", i);
        return (-1);
    }

    k = i;
    if (k==0)
        return DEFAULT_ABI_SPACING;

    while ((k>0) && (j<NUM_ABI_AVE_SPACING)) {
        k--;
        if (data->bases.coordinate[k] < data->bases.coordinate[k+1])
            j++;
    }
    if (j==0)
        return DEFAULT_ABI_SPACING;

    return (double)(data->bases.coordinate[i] -
                    data->bases.coordinate[k])/(double)j;
}

/*******************************************************************************
 * Function: is_first_poorly_resolved_peak
 *******************************************************************************
 */
static int
is_first_poorly_resolved_peak_sub( const Peak pk[], int peak_list_len,
                                   const int data_array[], BtkMessage *message)
{
    int        beg = pk->beg,
               end = pk->end,
               type = pk->type;
    double     height = pk->height,
               bheight= data_array[beg], 
               eheight= data_array[end];    

    if (type    == 13) return 1;  
    if (type%10 ==  1) return 0;

    if( pk->cd_peak_ind <= 0 ) /* beginning of read */
        return 1;
    if( pk->cd_peak_ind >= peak_list_len-1 )
        return 0;

    if (type == 12) {
        double nextheight = pk[+1].height;
        if ((eheight >= HIGH_BOUNDARY_FACTOR *     height) ||
            (eheight >= HIGH_BOUNDARY_FACTOR * nextheight))
        return 1;
    }
    
    if ((type == 23) && 
        (bheight < HIGH_BOUNDARY_FACTOR * height) &&
        (bheight < HIGH_BOUNDARY_FACTOR * (double)pk[-1].height))
        return 1;

    if (type == 22) {
        double prevheight = pk[-1].height; 
        double nextheight = pk[+1].height; 
        if ((bheight <  HIGH_BOUNDARY_FACTOR *     height) &&
            (bheight <  HIGH_BOUNDARY_FACTOR * prevheight) &&
            (eheight >= HIGH_BOUNDARY_FACTOR *     height) &&
            (eheight >= HIGH_BOUNDARY_FACTOR * nextheight))
        return 1;
    }

    return 0;
}

static int
is_first_poorly_resolved_peak(const Data *data, int jc, int peak_cd_ind,
    BtkMessage *message)
{
    const ColorData *cd = &(data->color_data[jc]);
    const Peak* pk = &(cd->peak_list[peak_cd_ind]);
    int* data_array = cd->data;
    return is_first_poorly_resolved_peak_sub( pk, cd->peak_list_len, 
                                              data_array, message );
}

/*******************************************************************************
 * Function: is_last_poorly_resolved_peak
 *******************************************************************************
 */
static int
is_last_poorly_resolved_peak_sub( const Peak pk[], int peak_list_len, 
                                  const int data_array[], BtkMessage *message )

{
    int        end = pk->end,
               type = pk->type;
    double     height = pk->height,
               eheight= data_array[end];

    if (type % 10 == 1) return 1;
    if (type % 10 == 3) return 0;
    if (pk->cd_peak_ind >= peak_list_len-1 ) /* end of read */
        return 1;

    if (type % 10 == 2) {
        double nextheight = pk[+1].height; 

        if ((eheight < HIGH_BOUNDARY_FACTOR *     height) &&
            (eheight < HIGH_BOUNDARY_FACTOR * nextheight))
        return 1;
    }

    return 0;
}

static int
is_last_poorly_resolved_peak( const Data *data, int jc, int peak_cd_ind,
                              BtkMessage *message )
{
    const ColorData *cd = &(data->color_data[jc]);
    const Peak* pk = &(cd->peak_list[peak_cd_ind]);
    int* data_array = cd->data;
    return is_last_poorly_resolved_peak_sub( pk, cd->peak_list_len,
                                             data_array, message );
}

/*******************************************************************************
 * Function: get_expected_multiplicity
 * Purpose:  determine expected number of individual peaks that a given
 *           multiple peak should be resolved to
 *******************************************************************************
 */
static double 
get_expected_multiplicity(Data *data, int jc, int first_peak_cd_ind, 
    int num_peaks_observed, BtkMessage *message, double* spacing )
{
    int    last_peak_cd_ind;
    double ave_num_peaks_expected, ave_spacing;
    Peak *peak = data->color_data[jc].peak_list;
    *spacing = ave_spacing = get_average_abi_spacing(
        peak[first_peak_cd_ind].pos, data, message);
    last_peak_cd_ind = first_peak_cd_ind+num_peaks_observed-1;
    ave_num_peaks_expected = 1 + 
        (double)(peak[last_peak_cd_ind].pos - peak[first_peak_cd_ind].pos) /
        ave_spacing;
    return ave_num_peaks_expected;
}
 
/*******************************************************************************
 * Function: colordata_expand_peaks 
 * Purpose: modify members of a peak list of given color using a new definition 
 *          of a peak, which results in an extension of the peak boundaries 
 *******************************************************************************
 */ 
static int 
colordata_expand_peaks(ColorData *cd, BtkMessage *message) 
{ 
    int        i, j; 
    int        height_min, midpoint;             
    int Case = 0; 

    /* Loop through peak list; make preliminary estimates  
     * for "outer" peak boundaries  
     */ 

    /* 1st "iteration": set lower and upper limits for peak boundaries  */    

    cd->peak_list[0].beg = 0;
    if (cd->peak_list_len > 1) {
        cd->peak_list[0].end = cd->peak_list[1].ibeg;
    }
    else {
        cd->peak_list[0].end = cd->length-1;
    }

    if (cd->peak_list_len > 1) {
        for (i=1; i<cd->peak_list_len-1; i++) { 
            cd->peak_list[i].beg = cd->peak_list[i-1].iend; 
            cd->peak_list[i].end = cd->peak_list[i+1].ibeg; 
        }         /* end loop in i  */  
        cd->peak_list[cd->peak_list_len-1].beg = 
            cd->peak_list[cd->peak_list_len-2].iend;
        cd->peak_list[cd->peak_list_len-1].end = cd->length-1;
    }
    else
        return SUCCESS;
   

    /* 2nd "iteration": expand peak boundaries */ 

    for (i=0; i<cd->peak_list_len; i++) {    
        /* left boundary */  
        if (cd->peak_list[i].beg < cd->peak_list[i].ibeg)
            for (j =  cd->peak_list[i].ibeg; j > cd->peak_list[i].beg; j--) 
            { 
                if ((cd->data[j] <= 0) || is_minimum(j, cd->data, cd->length)) 
                { 
	            cd->peak_list[i].beg = j; 
                    break; 
                } 
            }     
        /* right boundary */ 
        if (cd->peak_list[i].iend < cd->peak_list[i].end)
            for (j =  cd->peak_list[i].iend; j < cd->peak_list[i].end; j++) 
            { 
                if ((cd->data[j] <= 0) || is_minimum(j, cd->data, cd->length)) 
                { 
                    cd->peak_list[i].end   = j; 
                    break; 
                } 
            }
    }         /* end loop in i  */ 
 
    /* Loop through all peaks but the last once again;  
     * revise the peak boundaries and determine each peak's type        
     */ 
    if (cd->peak_list[0].beg > 0) { 
       cd->peak_list[0].type = 1; 
    } 
    else { 
       cd->peak_list[0].type = 3; 
    } 
 
    for (i=1; i <= cd->peak_list_len-1; i++) { 
#if 0
        if (cd->peak_list[i].max < 0)
            fprintf(stderr, 
                "In expand_peaks: color=%d peak[%d].max==%d\n",
                cd->peak_list[i].color_index, i, cd->peak_list[i].max);
#endif
    
        /* Case R1:  
         *********** 
         * Boundaries of adjacent peaks are distinct (type 1);  
         * accept them as they are              
         */ 
        if ((cd->peak_list[i].beg - cd->peak_list[i-1].end >
               MIN_DIST_BETWEEN_BOUNDARIES)
              ||
            ((cd->peak_list[i].beg - cd->peak_list[i-1].end > 0) &&
            (
             (
              ((double)cd->data[cd->peak_list[i  ].beg] <
               (double)cd->peak_list[i  ].height        * LOW_BOUNDARY_FACTOR) &&
              ((double)cd->data[cd->peak_list[i-1].end] <
               (double)cd->peak_list[i-1].height        * LOW_BOUNDARY_FACTOR)
             )
              ||
             /* Peak may be very wide and short, so its area will be big 
              * enough for the peak to be included into peak list, but
              * height is not big enough to meet the requirements above.
              * However, we still want identify its boundaries as type 1.
              * For this case, we add two alternative requirements below
              */
             ((double)cd->data[cd->peak_list[i  ].beg] < MIN_PEAK_IHEIGHT)
              ||
             ((double)cd->data[cd->peak_list[i-1].end] < MIN_PEAK_IHEIGHT)
            ))
           ) 
        { 
            cd->peak_list[i-1].type = 1 + 10*cd->peak_list[i-1].type; 
            cd->peak_list[i  ].type = 1; 
            cd->peak_list[i-1].area =  
                get_peak_area(cd->data, cd->peak_list[i-1].beg,  
                              cd->peak_list[i-1].end, message); 
            cd->peak_list[i-1].pos =  
                get_peak_position(cd->data, cd->peak_list[i-1].beg, 
                    cd->peak_list[i-1].end, cd->peak_list[i-1].area, message); 
            cd->peak_list[i-1].max =
                get_peak_max(cd->data, cd->peak_list[i-1].beg,
                                  cd->peak_list[i-1].end, message);
            if (cd->peak_list[i-1].pos == ERROR) { 
                return ERROR; 
            } 
            if ((i > 1) && (cd->peak_list[i-1].pos < cd->peak_list[i-2].pos)) { 
                (void)sprintf(message->text, 
                    "Wrong order of peak position #%d in case R1\n", i-1); 
                return ERROR; 
            } 
            cd->peak_list[i-1].height = cd->data[cd->peak_list[i-1].pos]; 
            cd->peak_list[i-1].iheight = (double)cd->peak_list[i-1].height; 
 
            cd->peak_list[i-1].width1 = get_peak_width1(cd, i-1, message); 
            if (cd->peak_list[i-1].width1 < 0) { 
/*              sprintf(message->text, 
 *                  "Width1 <0 for peak[%d], trace %d, case R1\n",  
 *                  i-1, cd->dye_number-1); 
 */ 
                return ERROR; 
            } 
 
            if ((i==1) || (cd->peak_list[i-1].type == 11)) { 
                if ((cd->peak_list[i-1].width2 =  
                    get_peak_width2(cd, i-1, message)) <= 0.) { 
                    return ERROR; 
                } 
            } 
            else { 
                cd->peak_list[i-1].width2 = cd->peak_list[i-2].width2; 
            }     
            if (cd->peak_list[i-1].width2 < 0) { 
                return ERROR; 
            } 
            if ((cd->peak_list[i-1].pos < cd->peak_list[i-1].beg) || 
                (cd->peak_list[i-1].pos > cd->peak_list[i-1].end))   
            { 
                (void)sprintf(message->text, 
                    "Position of peak #%d out of bounds in case R1\n", i-1); 
                return ERROR; 
            } 
            if (SHOW_REDEFINE_CASE) { 
                Case = 1; 
                fprintf(stderr, 
                "In expand_peaks: color=%d ind=%d type=%d pos=%d Case=R%d\n", 
                    cd->peak_list[i-1].color_index, i-1, 
                    cd->peak_list[i-1].type, cd->peak_list[i-1].pos,
                    Case); 
            }
            continue; 
        } 
     
        /* Case R2:  
         ********* 
         */  
        else { 
            height_min = QVMIN(cd->peak_list[i].height,  
                cd->peak_list[i-1].height); 
 
            /* Case R2.1 
             *********** 
             * signal is small at the junction (type 1) -> accept the boundaries
             * as if it was zero 
             */ 
            if ((cd->peak_list[i].beg - cd->peak_list[i-1].end <= 0) && 
                ((double)(cd->data[cd->peak_list[i  ].beg]) <     
                 (double)(height_min) * LOW_BOUNDARY_FACTOR) && 
                ((double)(cd->data[cd->peak_list[i-1].end]) <     
                 (double)(height_min) * LOW_BOUNDARY_FACTOR)) 
            { 
                cd->peak_list[i-1].end = cd->peak_list[i].beg; 
                cd->peak_list[i-1].type = 1 + 10*cd->peak_list[i-1].type; 
                cd->peak_list[i  ].type = 1; 
                cd->peak_list[i-1].area =  
                    get_peak_area(cd->data, cd->peak_list[i-1].beg,  
                                            cd->peak_list[i-1].end, message); 
                cd->peak_list[i-1].pos = get_peak_position(cd->data,  
                    cd->peak_list[i-1].beg, cd->peak_list[i-1].end, 
                    cd->peak_list[i-1].area, message); 
                cd->peak_list[i-1].max =
                    get_peak_max(cd->data, cd->peak_list[i-1].beg,
                                  cd->peak_list[i-1].end, message);
                if (cd->peak_list[i-1].pos == ERROR) { 
                    return ERROR; 
                } 
                if ((i>1) && (cd->peak_list[i-1].pos < cd->peak_list[i-2].pos)){
                    (void)sprintf(message->text, 
                    "Wrong order of peak position #%d in case R2.1\n", i-1); 
                   return ERROR; 
                } 
                cd->peak_list[i-1].height = cd->data[cd->peak_list[i-1].pos]; 
                cd->peak_list[i-1].iheight = (double)cd->peak_list[i-1].height; 
 
                cd->peak_list[i-1].width1 = get_peak_width1(cd, i-1, message); 
                if (cd->peak_list[i-1].width1 < 0) { 
                    sprintf(message->text,"Width1 <0 in case R2.1\n");    
                    return ERROR; 
                } 
 
                if ((i==1) || (cd->peak_list[i-1].type == 11)) { 
                    if ((cd->peak_list[i-1].width2 =  
                        get_peak_width2(cd, i-1, message)) <= 0) { 
                        return ERROR; 
                    } 
                } 
                else { 
                    cd->peak_list[i-1].width2 = cd->peak_list[i-2].width2; 
                } 
 
                if ((cd->peak_list[i-1].pos < cd->peak_list[i-1].beg) || 
                    (cd->peak_list[i-1].pos > cd->peak_list[i-1].end))   
                { 
                    (void)sprintf(message->text, 
                        "Position of peak #%d out of bounds in case R2.1\n",  
                        i-1); 
                    return ERROR; 
                } 
                if (SHOW_REDEFINE_CASE) { 
                    Case = 21; 
                    fprintf(stderr, 
                    "In   expand_peaks: color=%d ind=%d type=%d pos=%d Case=R%d\n", 
                    cd->peak_list[i-1].color_index, i-1, 
                    cd->peak_list[i-1].type, cd->peak_list[i-1].pos, Case); 
                }
                continue; 
            } 
             
            /* Case R2.2 
             *********** 
             * peak boundaries at local minimum (type 2) - accept them          
             */  
            else if ((cd->peak_list[i].beg - cd->peak_list[i-1].end >=0) && 
                     ( 
                      ((double)(cd->data[cd->peak_list[i  ].beg]) >=
                          (double)(height_min) * LOW_BOUNDARY_FACTOR) || 
                      ((double)(cd->data[cd->peak_list[i-1].end]) >=
                          (double)(height_min) * LOW_BOUNDARY_FACTOR) 
                     )                                                   &&
                     /* Take care of the case where ibeg==beg or iend==end,
                      * so no extension was needed. This is boundary of type
                      * 3, not 2. It may appaer, for example, as the result of 
                      * splitting truncated peak immediately after detection
                      */
                     (cd->peak_list[i  ].beg < cd->peak_list[i  ].ibeg)  &&
                     (cd->peak_list[i-1].end > cd->peak_list[i-1].iend)
                    ) 
            { 
                cd->peak_list[i-1].type = 2 + 10*cd->peak_list[i-1].type; 
                cd->peak_list[i  ].type = 2; 
                cd->peak_list[i-1].area =  
                    get_peak_area(cd->data, cd->peak_list[i-1].beg,  
                                            cd->peak_list[i-1].end, message); 
                cd->peak_list[i-1].pos =  
                    get_peak_position(cd->data, cd->peak_list[i-1].beg, 
                    cd->peak_list[i-1].end, cd->peak_list[i-1].area, message); 
                cd->peak_list[i-1].max =
                    get_peak_max(cd->data, cd->peak_list[i-1].beg,
                    cd->peak_list[i-1].end, message);
                if (cd->peak_list[i-1].pos == ERROR) { 
                    return ERROR; 
                } 
                if ((i>1) && (cd->peak_list[i-1].pos < cd->peak_list[i-2].pos)){
                    (void)sprintf(message->text, 
                        "Wrong order of peak position #%d in case R2.2\n",  
                        i-1); 
                   return ERROR; 
                } 
                cd->peak_list[i-1].height = cd->data[cd->peak_list[i-1].pos]; 
                cd->peak_list[i-1].iheight = (double)cd->peak_list[i-1].height; 
                cd->peak_list[i-1].width1 = get_peak_width1(cd, i-1, message); 
                if (cd->peak_list[i-1].width1 < 0) { 
                    sprintf(message->text,"Width1 <0 in case R2.1\n"); 
                    return ERROR; 
                } 
                if (i==1) { 
                   if ((cd->peak_list[i-1].width2 =  
                       get_peak_width2(cd, i-1, message)) <= 0.) { 
                       return ERROR; 
                   } 
                } 
                else { 
                    cd->peak_list[i-1].width2 = cd->peak_list[i-2].width2; 
                } 
                if ((cd->peak_list[i-1].pos < cd->peak_list[i-1].beg) || 
                    (cd->peak_list[i-1].pos > cd->peak_list[i-1].end))   
                { 
                    (void)sprintf(message->text, 
                        "Position of peak #%d out of bounds in case R2.2\n",  
                        i-1); 
                    return ERROR; 
                } 
                     
                if (SHOW_REDEFINE_CASE) {
                    Case = 22; 
                    fprintf(stderr, 
                    "In   expand_peaks: color=%d ind=%d type=%d pos=%d Case=R%d\n", 
                    cd->peak_list[i-1].color_index,  i-1, 
                    cd->peak_list[i-1].type, cd->peak_list[i-1].pos, Case); 
                }
                continue; 
            } 
 
            /* Case R2.3 
             *********** 
             * boundaries overlap (type 3) - set them to the midpoint  
             * between iend of the previous and ibeg of the next peak 
             */ 
            else { 
                cd->peak_list[i-1].type = 3 + 10*cd->peak_list[i-1].type; 
                cd->peak_list[i  ].type = 3; 
                midpoint = (cd->peak_list[i-1].iend + cd->peak_list[i].ibeg)/2; 
                cd->peak_list[i-1].end = midpoint; 
                cd->peak_list[i  ].beg = midpoint; 
                cd->peak_list[i-1].area = 
                    get_peak_area(cd->data, cd->peak_list[i-1].beg,  
                                            cd->peak_list[i-1].end, message);  
                cd->peak_list[i-1].pos =  
                    get_peak_position(cd->data, cd->peak_list[i-1].beg, 
                    cd->peak_list[i-1].end, cd->peak_list[i-1].area, message); 
                cd->peak_list[i-1].max =
                    get_peak_max(cd->data, cd->peak_list[i-1].beg,
                    cd->peak_list[i-1].end, message);
                cd->peak_list[i-1].height = cd->data[cd->peak_list[i-1].pos]; 
                cd->peak_list[i-1].iheight = (double)cd->peak_list[i-1].height; 
                if (i==1) { 
                    if ((cd->peak_list[i-1].width2 =  
                       get_peak_width2(cd, i-1, message)) <= 0.) { 
                       return ERROR; 
                   } 
                } 
                else if (i>1) { 
                    cd->peak_list[i-1].width2 = cd->peak_list[i-2].width2; 
                } 
                if ((i>1) && (cd->peak_list[i-1].pos < cd->peak_list[i-2].pos)){
                    (void)sprintf(message->text, 
                    "Wrong order of peak1 position #%d in case R2.3\n",  
                        i-1); 
                    return ERROR; 
                } 
                if (cd->peak_list[i].pos < cd->peak_list[i-1].pos) { 
                    (void)sprintf(message->text, 
                    "Wrong order of peak2 (#%d) position in case R2.3\n", i); 
                   return ERROR; 
                } 
            }       /* split peaks */ 
        }           /* outer boundaries either coincide or overlap */ 
 
        if ((cd->peak_list[i-1].pos < cd->peak_list[i-1].beg) || 
            (cd->peak_list[i-1].pos > cd->peak_list[i-1].end))   
        { 
            (void)sprintf(message->text, 
                "Position of peak #%d out of bounds in case R2.3\n", i-1); 
            return ERROR; 
        } 
        if (cd->peak_list[i-1].height <= 0) { 
            fprintf(stderr, "Error in   expand_peaks: "); 
            fprintf(stderr, "Peak(color=%d, index=%d).height=%d\n", 
                cd->dye_number-1, i-1, cd->peak_list[i-1].height);
            return ERROR; 
        } 
        if (SHOW_REDEFINE_CASE) { 
            Case = 23; 
            fprintf(stderr, 
            "In   expand_peaks: color=%d ind=%d type=%d pos=%d Case=R%d\n", 
            cd->peak_list[i-1].color_index,
            i-1, cd->peak_list[i-1].type, cd->peak_list[i-1].pos, Case); 
        }
    }               /* loop through all but the last peak */ 
 
    /* Case R3 
     *********** 
     * last peak  
     */ 
    i = cd->peak_list_len-1; 
    /* Force the last peak to be of type *1; need this for 
     * resolving multiple peaks  
     */ 
    cd->peak_list[i].type = 1 + 10*cd->peak_list[i].type; 
    cd->peak_list[i].area = get_peak_area(cd->data, cd->peak_list[i].beg,  
        cd->peak_list[i].end, message); 
    cd->peak_list[i].pos = get_peak_position(cd->data, cd->peak_list[i].beg,  
        cd->peak_list[i].end, cd->peak_list[i].area, message); 
    cd->peak_list[i].max =
        get_peak_max(cd->data, cd->peak_list[i].beg,
        cd->peak_list[i].end, message);
    if (cd->peak_list[i].pos == ERROR) { 
        return ERROR; 
    } 
    if ((cd->peak_list_len > 1) &&  
        (cd->peak_list[i].pos < cd->peak_list[i-1].pos))  
    { 
        (void)sprintf(message->text, 
         "Wrong order of peak #%d position in case R3\n",  
            i); 
         return ERROR; 
    } 
    cd->peak_list[i].height = cd->data[cd->peak_list[i].pos]; 
    cd->peak_list[i].iheight = (double)cd->peak_list[i].height; 
 
    cd->peak_list[i].width1 = get_peak_width1(cd, i, message); 
    if (cd->peak_list[i].width1 <= 0) { 
        sprintf(message->text,"Peak[%d].width1 = %f in case R3; peak.beg=%d peak.end=%d\n",  
        i, cd->peak_list[i].width1, cd->peak_list[i].beg, cd->peak_list[i].end);    
        return ERROR; 
    } 
 
    if (i==0) { 
        if ((cd->peak_list[i].width2 =  
            get_peak_width2(cd, i, message)) <= 0.) { 
            return ERROR; 
        } 
    } 
    else { 
        cd->peak_list[i].width2 = cd->peak_list[i-1].width2; 
    } 
    if ((cd->peak_list[i].pos < cd->peak_list[i].beg) || 
        (cd->peak_list[i].pos > cd->peak_list[i].end))   
    { 
        (void)sprintf(message->text, 
            "Position of peak #%d out of bounds in case R3\n", i); 
        return ERROR; 
    } 
    
     /* In the above processing, type 1 boundaries are sometimes 
      * mis-classified as type 3.  Loop over the data once again to fix this
      */
    for ( i=0; i<cd->peak_list_len-1; i++ ) {
        Peak 
            *pk0 = &cd->peak_list[i],
            *pk1 = &cd->peak_list[i+1];
        if( pk0->end == pk1->beg ) { 
            double 
                mean_hite = 0.5*(pk0->height+pk1->height),
                mid_hite = cd->data[pk0->end];
            if( mid_hite < mean_hite *  LOW_BOUNDARY_FACTOR ) {
                int
                    type0 = pk0->type%10,
                    type1 = pk1->type/10;
                if( type0!=type1 ) { fprintf(stderr,"pk0->type!=pk1->type\n"); }
                if( type0!=1 ) {
                    pk0->type = 10*(pk0->type/10) + 1;
                    pk1->type = pk1->type%10 + 10;
                } 
            }
        }
    }
    
    return SUCCESS; 
}

#if 0  /* debug routine */
#define PP( field,format) fprintf( fp, #field "=%" #format "\n" , pk->##field );
static void 
print_peak( const Peak *pk, FILE *fp )
{
#if 0
    PP( area, g );
    PP( base_index, d );
    PP( color_index, d );
    PP( data_peak_ind, d );
    PP( data_peak_ind2, d );
    PP( ibeg, d );
    PP( iend, d );
#endif
    PP( base, c );
    PP( cd_peak_ind, d );
    PP( beg, d );
    PP( pos, d );
    PP( end, d );
    PP( height, d );
#if 0
    PP( iheight, g );
    PP( wiheight, g );
    PP( is_called, d );
    PP( is_truncated, d );
    PP( ipos, d );
    PP( ipos_orig, d );
    PP( max, d );
    PP( relative_area, g );
    PP( spacing, d );
    PP( type, d );
    PP( width1, g );
    PP( width2, g );
    PP( ave_width1, g );
    PP( ave_width2, g );
    PP( orig_width, g );
    PP( beta, g );
    PP( ave_w02beta, g );
    PP( C0, g );
    PP( resolution, g );
#endif
}
#endif

/*******************************************************************************
 * Function: fix_peak_number
 * Purpose:  Given that you currently have a group of peaks:
 *           (num_old peaks beginning at indx)
 *           and you now wish to have a new group of peaks
 *           (num_new peaks beginning at indx),
 *           this routine takes care of the bookeeping:
 *           (1) Allocating memory if needed.
 *           (2) Moving peaks in the array which are to the right of 
 *               this group of peaks.
 *           (3) Adjusting cd_peak_ind for peaks to the right.
 *******************************************************************************
 */
static int
fix_peak_number( Data *data, int color, int indx, int num_old, int num_new,
                 BtkMessage *message )
{
    ColorData *cd = &data->color_data[color];
    Peak *peak;
    int 
        num_move = cd->peak_list_len-(indx+num_old), 
        num_add = num_new-num_old, j;
    if( num_new > num_old ) {
        int resize=0;
        while( cd->peak_list_len+num_add > cd->peak_list_max_len ) {
            cd->peak_list_max_len *= 2;
            resize = 1;
        }
        if( resize ) {
            cd->peak_list = REALLOC( cd->peak_list, Peak,cd->peak_list_max_len);
            MEM_ERROR(cd->peak_list);
        }
    }
    peak = &cd->peak_list[indx];
    cd->peak_list_len += num_add;
    memmove( peak+num_new, peak+num_old, num_move*sizeof(Peak) );
    for( j=indx+num_new; j<cd->peak_list_len; j++ ) {
        cd->peak_list[j].cd_peak_ind += num_add;
    }
    return SUCCESS;
 error:
    return ERROR;
}

/****************************************************************************** 
 * Function: populate_group_of_peaks
 * Purpose:  Given a group of peaks (n peaks beginning at indx)
 *           This routine resets all their fields so as to have
 *           n peaks evenly spaced between the limits described by
 *           first->pos and last->ipos
 ******************************************************************************
 */
static void
populate_group_of_peaks( Data *data, int color, int indx, int n, 
                         const Peak *first, const Peak *last, 
                         BtkMessage *message, int debug )
{
    ColorData *cd = &data->color_data[color];
    Peak *peak = &cd->peak_list[indx];
    static char base[] = { 'A', 'C', 'G', 'T' };
    int i, m=n-1;
    double  del = n==1 ? 0.0 : (double)(last->pos-first->pos)/(double)(m);
    
    if( debug ) {
        fprintf( stderr, "first and last peak:\n" );
        print_one_peak( stderr, first );
        print_one_peak( stderr, last );
    }

    for( i=0; i<n; i++ ) {
        int type10, type1;
        double pos = first->pos + i*del;
        Peak *pk = &peak[i];
        
        if( debug ) {
            fprintf( stderr, "top of loop, i=%d\n", i );
            print_one_peak( stderr, pk );
        }
        
        if( n==1 ) {
            pos = 0.5 * (first->pos+last->pos);
            pk->beg = first->beg;
            pk->end = last->end;
        } else {
            pk->beg = ROUND( pos-0.5*del );
            pk->end = ROUND( pos+0.5*del );
            pk->beg = QVMAX( pk->beg, 0 );
            pk->end = QVMIN( pk->end, cd->length-1 );
            pk->beg = QVMAX( pk->beg, first->beg );
            pk->end = QVMIN( pk->end, last->end );
            if( i==0 ) { pk->beg = first->beg; }
            if( i==m ) { pk->end = last->end; }
        }

        /* -1 ==> don't care */
        pk->area= get_peak_area(cd->data, pk->beg, pk->end, message);
        pk->pos = get_peak_position(cd->data, pk->beg, pk->end, pk->area,
                  message);
        pk->max = get_peak_max(cd->data, pk->beg, pk->end, message);
        pk->ipos_orig      = pk->pos;
        pk->ipos           = pk->pos;
        pk->base           = base[color];
        pk->base_index     = -1;
        pk->color_index    = color;
        pk->cd_peak_ind    = indx+i;
        pk->data_peak_ind  = -1;
        pk->data_peak_ind2 = -1;
        pk->height         = cd->data[pk->pos];
        pk->iheight        = pk->height;
        pk->wiheight       = pk->height;
        pk->is_called      = 0;
        pk->is_truncated   = ( pk->height >= TRUNCATED_HEIGHT );
        pk->relative_area  = 1.0;
        pk->spacing        = ROUND(del);
        /* Make peaks in this group of peaks to be of types: x3 33 ... 33 3y 
         * where x & y are determined by the 1st and last peaks
         */
        type10 = (i==0 ? 10*(first->type/10) : 30);
        type1  = (i==m ? last->type%10       :  3);
        pk->type = type10 + type1;

        pk->width1 = pk->width2     = pk->ave_width1 = pk->ave_width2 
                   = pk->orig_width = pk->beta       = pk->ave_w02beta
                   = pk->C0         = pk->resolution = -1.0;

        if( debug ) {
            fprintf( stderr, "bot of loop, i=%d\n", i );
            print_one_peak( stderr, pk );
        }
    }
}


static void
print_peak_sub_list( FILE *fp, const Data *data, int color, 
                     int indx, int nn, BtkMessage *message )
{
    int j;
    const ColorData *cd = &data->color_data[color];
    fprintf( fp, "f l "); print_one_peak_header(fp);
    for (j=indx; j<indx+nn; j++) {
        int 
            f = is_first_poorly_resolved_peak(data, color, j, message),
            l = is_last_poorly_resolved_peak(data, color, j, message );
        Peak* pk = &(cd->peak_list[j]);
        fprintf( fp, "%1d %1d ", f, l );
        print_one_peak(fp, pk );
    }
}

#if 0 /* debug routine to print out info on all peaks */
static void
print_peak_list( FILE *fp, const Data *data, int color, BtkMessage *message )
{
    print_peak_sub_list( fp, data, color, 0, 
                         data->color_data[color].peak_list_len, message );
}
static void
print_all_peak_lists( FILE *fp, const Data *data, BtkMessage *message )
{
    int i;
    for( i=0; i<4; i++ ) {
        print_peak_list( fp, data, i, message );
    }
}
void
my_print_peak_list( const char* fname, const Data* data, int color,
                    BtkMessage *message )
{
    FILE *fp = fopen( fname, "w" );
    assert(fp);
    print_peak_list( fp, data, color, message );
}
#endif

static void 
show_peak_multiplicity( ColorData *cd, int color,
                    int j, 
                    int num_peaks_observed,
                    double ave_num_peaks_expected,
                    int num_peaks_expected, double spacing )
{
    int i;
    double eps = ave_num_peaks_expected - num_peaks_expected;
    int n = num_peaks_observed;
    Peak *peak = cd->peak_list;
    fprintf(stderr,
            "%c i=%d,%d pos %d - %d = %2d sp=%4.1f "
            "pks: obs=%d %c exp=%d = %5.4f%+7.4f\n",            
            peak[j].base, peak[j].cd_peak_ind, peak[j].data_peak_ind,
            peak[j].pos, peak[j+n-1].pos,
            peak[j+n-1].pos-peak[j].pos,
            spacing,
            num_peaks_observed,
            (num_peaks_observed==num_peaks_expected ? '=' :
            (num_peaks_observed<num_peaks_expected ? '<': '>')),
            num_peaks_expected, ave_num_peaks_expected,
            -eps );
    
    fprintf(stderr, "         Observed peak types=");
    for (i=0; i<n; i++) {
        fprintf(stderr, "%d ", peak[j+i].type);
    }
    fprintf(stderr, "\n");
}

static int
compute_num_peaks_expected( int num_peaks_observed, 
                            double ave_num_peaks_expected )
{
    int num_peaks_expected;
    double
        eps = fabs(ave_num_peaks_expected - ROUND(ave_num_peaks_expected));
    if( eps< MULTIPLICITY_THRESHOLD ) {
        num_peaks_expected = ROUND(ave_num_peaks_expected);
    } else {
        if( ave_num_peaks_expected > num_peaks_observed ) {
            num_peaks_expected = (int)floor(ave_num_peaks_expected);
        } else {
            num_peaks_expected = (int)ceil(ave_num_peaks_expected);
        }
    }
    return num_peaks_expected;
}
double
array_inverse_CDF( const int h[], int size, double p )
{
    double sum, v0, v1, vtest, d0, d1, retval;
    int i, i0, i1;
    assert( p>=0 );
    assert( p<=1.0 );

    sum = 0;
    for( i=0; i<size; i++ ) { sum += h[i]; }
    
    vtest = p*sum;

    if( p<0.5 ) {
        v0 = v1 = 0;
        i1 = -1;
        do {
            i1++;
            v0 = v1;
            v1 += h[i1];
        } while( v1 < vtest );
        i0 = i1 - 1;
        d0 = i0+1.0;
        d1 = i1+1.0;
    } else {
        v0 = v1 = sum;
        i0 = size;
        do {
            i0--;
            v1 = v0;
            v0 -= h[i0];
        } while( v0 >= vtest );
        i1 = i0 + 1;
        d0 = i0;
        d1 = i1;
    }

    retval = ( d0*(v1-vtest) + d1*(vtest-v0) ) / (v1-v0);    
    return retval - 0.5;
}

void
flesh_out_peak( int data[], int color, int cd_peak_ind,
                Peak *pk, BtkMessage *message )
{
    static char base[] = { 'A', 'C', 'G', 'T' };
    pk->area = get_peak_area(data, pk->beg, pk->end, message);
    pk->pos  = get_peak_position(data, pk->beg, pk->end, pk->area,
                                message);
    pk->max  = get_peak_max(data, pk->beg, pk->end, message);
    pk->ipos_orig      = pk->pos;
    pk->ipos           = pk->pos;
    pk->base           = base[color];
    pk->base_index     = -1;
    pk->color_index    = color;
    pk->cd_peak_ind    = cd_peak_ind;
    pk->data_peak_ind  = -1;
    pk->data_peak_ind2 = -1;
    pk->height         = data[pk->pos];
    pk->iheight        = pk->height;
    pk->wiheight       = pk->height;
    pk->is_called      = 0;
    pk->is_truncated   = ( pk->height >= TRUNCATED_HEIGHT );
    pk->relative_area  = 1.0;
    if( cd_peak_ind==0 ) {
        pk->spacing = (int)DEFAULT_ABI_SPACING;
    } else {
        pk->spacing = pk->pos - pk[-1].pos;
    }
    pk->width1 = pk->width2     = pk->ave_width1 = pk->ave_width2 
               = pk->orig_width = pk->beta       = pk->ave_w02beta
               = pk->C0         = pk->resolution = -1.0;        
}

static void
add_one_peak( Data *data, int color, int indx, int num_old,
              BtkMessage *message )
{
    ColorData *cd = &data->color_data[color];
    int j, j_max_beg_end=0, t0, t1;

    Peak peak0, peak1, *p0, *p1, *p2;

    /* find widest peak */
    {
        int max_beg_end=-1;
        for( j=0; j<num_old; j++ ) {
            Peak *p0 = &cd->peak_list[indx+j];
            int beg_end = p0->end - p0->beg;
            if( beg_end > max_beg_end ) {
                max_beg_end = beg_end;
                j_max_beg_end = j;
            }
        }
    }

    peak0=cd->peak_list[indx+j_max_beg_end];
    peak1=cd->peak_list[indx+j_max_beg_end+1];

    fix_peak_number( data, color, indx+j_max_beg_end+1, 0, 1, message );

    p0 = &cd->peak_list[indx+j_max_beg_end+0];
    p1 = &cd->peak_list[indx+j_max_beg_end+1];
    p2 = &cd->peak_list[indx+j_max_beg_end+2];

    {
        int *local_data = &cd->data[p0->beg];
        int size = p0->end - p0->beg + 1;
        double boundary = p0->beg + array_inverse_CDF( local_data, size, 0.5 );
        p1->end = p0->end;
        p0->end = ROUND(boundary);
        p1->beg = p0->end;
    }
    t0 = p0->type/10;
    t1 = p0->type%10;
    p0->type = 10*t0 + 3;
    p1->type = 30 + t1;

    flesh_out_peak( cd->data, color, indx+j_max_beg_end+0, p0, message );
    flesh_out_peak( cd->data, color, indx+j_max_beg_end+1, p1, message );
}

static void
rem_one_peak( Data *data, int color, int indx, int num_old,
              BtkMessage *message )
{
    ColorData *cd = &data->color_data[color];
    int j, j_min_beg_end=0, beg, end, new_type;
    Peak *p0, *p1;
    /* find narrowest peak pair */
    {
        int min_beg_end=10000; /* infinity */
        for( j=0; j<num_old-1; j++ ) {
            int beg_end;
            p0 = &cd->peak_list[indx+j];
            p1 = &cd->peak_list[indx+j+1];
            beg_end = p1->end - p0->beg;
            if( beg_end < min_beg_end ) {
                min_beg_end = beg_end;
                j_min_beg_end = j;
            }
        }
    }

    p0 = &cd->peak_list[indx+j_min_beg_end+0];
    p1 = &cd->peak_list[indx+j_min_beg_end+1];

    new_type = 10*(p0->type/10) + p1->type%10;
    beg = p0->beg;
    end = p1->end;

    fix_peak_number( data, color, indx+j_min_beg_end+1, 1, 0, message );

    p0->type = new_type;
    p0->beg = beg;
    p0->end = end;

    flesh_out_peak( cd->data, color, indx+j_min_beg_end+0, p0, message );
}

/*******************************************************************************
 * Function: set_peak_widths      
 *******************************************************************************
 */

static void 
set_peak_widths(Data *data, BtkMessage *message)
{
    int i, j;

    for (i=0; i<NUM_COLORS; i++)
    {
        for (j=0; j<data->color_data[i].peak_list_len; j++)
        {
            Peak *peak = &data->color_data[i].peak_list[j];
            ColorData *cd = &data->color_data[i];
 
            if (peak->width1 < EPSILON)
                peak->width1 = get_peak_width1(cd, peak->cd_peak_ind, message);
            
            if (peak->width2 < EPSILON)
                peak->width2 = get_peak_width2(cd, peak->cd_peak_ind, message);
        }
    }
}

/*******************************************************************************
 * Function: data_expand_peaks              
 * Purpose: modify peaks in each of 4 peak lists using a new peak definition, 
 *          which extends peaks; expand   peak area, peak position and peak  
 *          height                 
 *******************************************************************************
 */ 
int 
data_expand_peaks(Data *data, Options *options, BtkMessage *message) 
{ 
    int color;

    for (color = 0; color < NUM_COLORS; color++) { 
        ColorData *cd = &data->color_data[color]; 
        if (cd->peak_list_len == 0) { 
            continue; 
        } 
        if (colordata_expand_peaks(cd, message) == ERROR) { 
            return ERROR; 
        }  
        set_peak_widths(data, message);
    }

#if 0
    /*for( color=0; color<NUM_COLORS; color++ )*/
    fprintf( stderr, "before\n" );
    color = 0;
    print_peak_list( stderr, data, color, message );
#endif
    
    if (!options->respace ) 
        { return SUCCESS; }

    for (color = 0; color < NUM_COLORS; color++) {
        int j;
        ColorData *cd = &data->color_data[color];
        for (j=0; j<cd->peak_list_len; j++) {
            int n, num_peaks_expected, num_peaks_observed;
            double ave_num_peaks_expected, spacing;
            int debug=0;
#if 0
            Peak *ref = &(cd->peak_list[j]);
            /*debug = ( color==3 && ref->pos>7430 && ref->pos<7470 );*/
            debug = ( color==0 && ref->pos>5900 && ref->pos<6000 );
#endif
            /* Detect the first poorly resolved peak in the row */
            if (!is_first_poorly_resolved_peak(data, color, j, message)) 
                continue;
            
            /* Detect the last poorly resolved peak in the row */
            n = 2;
            while (!is_last_poorly_resolved_peak(data, color, j+n-1, 
                            message) && j+n < cd->peak_list_len)
            {
                n++;
            }

            if( debug ) {
                fprintf( stderr, "before:\n" );
                print_peak_sub_list(stderr, data, color, j, n, message);
            }

            num_peaks_observed = n;
            ave_num_peaks_expected = get_expected_multiplicity(data,color,j,
                                        n, message, &spacing );

            num_peaks_expected = compute_num_peaks_expected( 
                                   num_peaks_observed, ave_num_peaks_expected );

            if( debug ) {
                show_peak_multiplicity( cd, color, j, num_peaks_observed,
                                        ave_num_peaks_expected,
                                        num_peaks_expected, spacing );
            }
            
            if(0) /* (num_peaks_expected != num_peaks_observed ) { */
            {
                int jj;
                fprintf( stderr, "color=%d obs=%d exp=%d ", 
                         color, 
                         num_peaks_observed, num_peaks_expected );
                for( jj=j; jj<j+num_peaks_observed; jj++ ) {
                    Peak *p0=&(cd->peak_list[jj]);
                    fprintf( stderr, " (%d %d %d %d t=%d) ", 
                             p0->beg, p0->pos, p0->end, p0->end-p0->beg,
                             p0->type);
                }
                fprintf( stderr, "\n" );
            }

            if (num_peaks_expected != num_peaks_observed ) {
                int new_num_peaks;
                Peak first, last;
                
                /* debug stuff: print out info on poorly resolved peaks */
#if SHOW_PEAK_MULTIPLICITY
                show_peak_multiplicity( cd, color, j, num_peaks_observed,
                                        ave_num_peaks_expected,
                                        num_peaks_expected, spacing );
#endif
                first = cd->peak_list[j];
                last  = cd->peak_list[j+n-1];
                
                if( 1 ) {       /*  Change by one peak at most */

                    if( num_peaks_expected > num_peaks_observed ) {
                        add_one_peak( data, color, j, num_peaks_observed,
                                    message );
                        new_num_peaks = num_peaks_observed+1;
                    } else {
                        if( num_peaks_expected < num_peaks_observed )
                        {
                            rem_one_peak( data, color, j, num_peaks_observed,
                                          message );
                            new_num_peaks = num_peaks_observed-1;
                        } else {
                            new_num_peaks = num_peaks_observed;
                        }
                    }

                } else { /* distribute peaks on a regular grid */
                    if( fix_peak_number( data, color, j, n, num_peaks_expected,
                                         message ) == ERROR ) return ERROR;
                    
                    populate_group_of_peaks( data, color, j, num_peaks_expected,
                                             &first, &last, message, debug );
                    new_num_peaks = num_peaks_expected;
                }
                if( debug ) {
                    fprintf( stderr, "after:\n" );
                    print_peak_sub_list(stderr, data, color, j, 
                                        num_peaks_expected, message );
                }

                if(0)
                {
                    int jj;
                    fprintf( stderr, " ---> " );
                    for( jj=j; jj<j+new_num_peaks; jj++ ) {
                        Peak *p0=&(cd->peak_list[jj]);
                        fprintf( stderr, " (%d %d %d %d) ", 
                                p0->beg, p0->pos, p0->end, p0->end-p0->beg );
                    }
                    fprintf( stderr, "\n" );
                }

            }                
        }
    }
#if 0
    /*for( color=0; color<NUM_COLORS; color++ )*/
    fprintf( stderr, "after\n" );
    color = 3;
    print_peak_list( stderr, data, color, message );
#endif
    return SUCCESS; 
} 

/*******************************************************************************
 * Function: single_peak_resolution
 * Purpose: Calculate the resolution parameter for the single peak fit.
 *******************************************************************************
 */ 
static double
single_peak_resolution( Data* data, const Peak* peak, Options* options )
{
    int j;
    double x, mdata, model, resolution;
    ColorData* cd = &data->color_data[peak->color_index];

    /* Compute the peak resolution */ 
    resolution = 0.; 
    for (j=peak->beg; j<peak->end; j++) { 
        x = (double)j + 0.5; 
        mdata = (double)(cd->data[j]+cd->data[j+1])/2.; 
        model= Shape( peak->C0, peak->orig_width/2., peak->beta, 
                      (x - (double)peak->pos), options); 
        resolution += fabs(mdata - model);
    }
    return resolution / peak->area;
}

static float
update_ipos(float B, int pos,  int height, int lpos, int lheight,
                     int rpos, int rheight)
{
    if      (lheight== rheight) return (float)pos;
    else if (height == lheight) return (float)(pos+lpos)/2.;
    else if (height == rheight) return (float)(pos+rpos)/2.;
    else if (height  < rheight) return (float)rpos;
    else if (height  < lheight) return (float)lpos;
    /* lheigh < height and rheight < height */

    /* Interpolate the three points with quadratic parabola:
     * y = y0 - A*(x-x0)^2
     */
    float y1 = lheight;
    float x1 = lpos;
    float y2 = height;
    float x2 = pos;
    float y3 = rheight;
    float x3 = rpos;
    return 0.5*(2.*(y2-y1)*(x3+x1)-(y3-y1)*(x2+x1)) 
              /(2.*(y2-y1)        -(y3-y1)        );
}


/*******************************************************************************
 * Function: fit_single_peak 
 * Purpose: estimate model parameters by fitting a  peak of type 11; 
 *          compute the peak resolution 
 *******************************************************************************
 */
static int 
fit_single_peak(Peak *pk, double stderr_w1, double stderr_w2,
    Data *data, Options *options, BtkMessage *message) 
{ 
    double w1, w2, width_ratio;    
    int color = pk->color_index;
    ColorData *cd = &data->color_data[color];
    double aw1, aw2, beta1, beta2;

    w1 = pk->width1;
    w2 = pk->width2;
    aw1 = pk->ave_width1;
    aw2 = pk->ave_width2;

    if (w1 < EPSILON) 
        pk->width1 = get_peak_width1(cd, pk->cd_peak_ind, message);
  
    w1    = ((w1 > aw1-stderr_w1) && (w1 < aw1+stderr_w1)) ? w1 : aw1; 
    w2    = ((w2 > aw2-stderr_w2) && (w2 < aw2+stderr_w2)) ? w2 : aw2;
    beta1 = w1 / TWO_SQRT_LN2;
    beta2 = w2 * M_1_SQRTPI;
    width_ratio = w1/w2;

    if( options->gauss ) {
        /* Gauss model: C = (C0/(4*PI*beta**2))*exp(-x**2/(4*beta**2)) */
        int i;
        double sum_exp = 0., sum_exp2 = 0., x2, B;
        pk->beta = (beta1 + beta2)/2.;        
        pk->orig_width = 0.0;
        B = 1./4./pk->beta/pk->beta;

        /* Set ipos */
        {
            int pos = pk->pos, height = cd->data[pos]; 
            int lpos = QVMAX(pos-1, 0), lheight = cd->data[lpos];
            int rpos = QVMIN(pos+1, data->length-1), rheight = cd->data[rpos];
            pk->ipos = 
                update_ipos(B, pos, height, lpos, lheight, rpos, rheight);
        }
        /* Set iheight */
        for (i=pk->beg; i<=pk->end; i++)
        {
            x2 = (double)((i-pk->ipos)*(i-pk->ipos));  
            sum_exp  += (double)(cd->data[i])
                       *Shape(1., 0.0, pk->beta, (double)(i-pk->ipos), options);
            sum_exp2 += Shape(1., 0.0, pk->beta, (double)(i-pk->ipos), options)
                       *Shape(1., 0.0, pk->beta, (double)(i-pk->ipos), options);
        }
        pk->C0 = pk->iheight = (sum_exp2) > 0 ? sum_exp / sum_exp2 : pk->height;
//      fprintf(stderr, "single ipos= %f beg=%d end=%d\n", pk->ipos, pk->beg, pk->end);
    } else {
        double w02beta;
        /* width_ratio = calc_ok_ratio( data, pk, width_ratio ); */
        w02beta = Phi(width_ratio);
        pk->beta   = w1 / 2./W1(w02beta); 
        pk->orig_width = 2.*w02beta * pk->beta; 
        pk->C0 = pk->iheight / Erf(w02beta);
    }

    pk->resolution = single_peak_resolution( data, pk, options );
#if 0
    fprintf(stderr, 
        "Fitting peak color=%d ipos=%d C0=%f beta=%f resolution=%f\n",
       pk->color_index, pk->ipos, pk->C0, pk->beta, pk->resolution); 
#endif
    return SUCCESS;    
}

/*******************************************************************************
 * Function: compute_deviation
 *  
 * Purpose: compute the deviation between the data and the model
 *          for the multiple peak fit
 *  
 *******************************************************************************
 */ 
static void
compute_deviation( const Peak peak_list[], const int data_array[], int j, int n,
                   const Options* options, double dev[] )
{
    int beg, end, m, m1;
    /* Compute the deviation of the data from the model  
     * Deviation is the data minus (all) the intrinsic peak  
     * signals but the current one 
     */ 
    beg = peak_list[j].beg; 
    end = peak_list[j].end; 
    if (j > 0) { 
        beg = peak_list[j-1].beg; 
    } 
    if (j < n-1) { 
        end = peak_list[j+1].end; 
    } 
    
    for (m1=beg; m1<=end; m1++) { 
        const Peak *pk;
        m = m1 - peak_list[0].beg; 
        dev[m] = (double)data_array[m1];  // initializing dev array 
        if (j > 0) {
            pk = &peak_list[j-1];
            dev[m] -= Shape(pk->C0, pk->orig_width/TWO, pk->beta, 
                            m1-pk->ipos, options );
        } 
        if (j > 1) { 
            pk = &peak_list[j-2];
            dev[m] -= Shape(pk->C0, pk->orig_width/TWO, pk->beta, 
                            m1-pk->ipos, options );
        } 
        if (j < n-1) { 
            pk = &peak_list[j+1];
            dev[m] -= Shape(pk->C0, pk->orig_width/TWO, pk->beta, 
                            m1-pk->ipos, options );
        } 
        if (j < n-2) {
            pk = &peak_list[j+2];
            dev[m] -= Shape(pk->C0, pk->orig_width/TWO, pk->beta, 
                            m1-pk->ipos, options );
        } 
        if (dev[m] < 0) { 
            dev[m] = 0.; 
        } 
    }
}

/****************************************************************************** 
 * Function: compute_peak_resolution
 * Purpose: compute the sum of the absolute value of the difference between
 *          the model and the data.
 *******************************************************************************
 */
static double
compute_peak_resolution( const Peak* pk, int begin, 
                         const double dev[], const Options* options )
{
    int m, m1;
    double dev_area;

    dev_area = 0.; 
    for (m1 = pk->beg; m1 < pk->end; m1++)
    { 
        m = m1 - begin; 
        dev_area += fabs((double)(dev[m]+dev[m+1])/TWO
                  - Shape(pk->C0, pk->orig_width/TWO, pk->beta,
                         (double)(2*m1 + 1 - 2*pk->ipos)/TWO, options));
    }
    if (pk->area > 0.) {  dev_area /= pk->area;  } 
    return dev_area;
}

/*******************************************************************************
 * Function: peak_weight
 * Purpose: For the functions below, which attempt to align the peaks with a 
 * regular lattice, I want to weight the peaks somehow.
 * The exact method is somewhat arbitrary.
 *******************************************************************************
 */
static double
peak_weight( const Peak* pk )
{
    double temp = pk->iheight;
    /* height^2 seems like a reasonable weight, but other choices make sense */
    return temp; /* NOTICE!! temp*temp; */
}

/*******************************************************************************
 *Function: calc_peak_cost
 * Calculate a measure of how close the peaks are to a regular grid.
 * I'm using a weighted least squares cost function.
 *******************************************************************************
 */
static double
calc_peak_cost( const Peak pk[], int nn, double spacing, double shift )
{
    int i;
    double ref_pos=pk[0].ipos, cost=0.0;
    for( i=0; i<nn; i++ ) {
        double 
            indx_i = i,  /* i --> indx[i] in general */
            x = pk[i].ipos,
            y = ( ref_pos + spacing * indx_i ),
            wt = peak_weight(&pk[i]),
            temp = x-y-shift;
        cost += wt*temp*temp;
    }
    return cost;
}

/*******************************************************************************
 * Function: solve_special_2x2
 * Purpose: Solve the following eqs. for (x,y) given (a,b,c,e,f)
 *    a*x + b*y = e
 *    b*x + c*y = f 
 *******************************************************************************
 */
static int 
solve_special_2x2( double a, double b, double c, double e, double f,
                   double* x, double* y )
{
    double ac=a*c, bb=b*b, det=ac-bb;
    if( fabs(det) <= EPSILON*bb ) {
        *x = *y = 0.0;
        return 0;       /* not OK */
    }
    *x = (c*e-b*f)/det;
    *y = (a*f-b*e)/det;
    return 1;           /* OK */
}

/*******************************************************************************
 * Function: calc_peak_shift_and_cost
 * Purpose: Calculate the shift needed to minimize the cost of aligning the 
 * peaks to a regular grid.
 * Peak spacing is assumed constant.
 * shift = shift w.r.t. pk[0]
 *******************************************************************************
 */
static int
calc_peak_shift_and_cost( const Peak pk[], int nn, double spacing, 
                          double* shift, double* cost )
{
    int i;
    double numer=0.0, denom=0.0;
    double ref_pos = pk[0].ipos;
    for( i=0; i<nn; i++ ) {
        double wt = peak_weight(&pk[i]),
               indx_i = i;  /* i --> indx[i] in general */
        numer += wt * ( pk[i].ipos - ( ref_pos + spacing * indx_i ) );
        denom += wt;
    }
    if( denom <= EPSILON*numer ) {
        fprintf( stderr, "numerical problem solving for shift\n" );
        *shift = *cost = 0.0;
        return 0;
    }
    *shift = numer/denom;

    *cost = calc_peak_cost( pk, nn, spacing, *shift );
    return 1;
}

/*******************************************************************************
 * Function: calc_peak_spacing_shift_cost
 * Purpose: Calculate the shift and spacing needed to minimize the cost of 
 * aligning the peaks to a regular grid.
 * shift = shift w.r.t. pk[0]
 *******************************************************************************
 */
static int
calc_peak_spacing_shift_cost( const Peak pk[], int nn, 
                              double* spacing, double* shift, double* cost )
{
    int i, ok;
    double a, b, c, e, f, ref_pos = pk[0].ipos;

    a = b = c = e = f = 0.0;
    for( i=0; i<nn; i++ ) {
        double wt = peak_weight(&pk[i]),
               indx_i = i;  /* i --> indx[i] in general */
        a += wt;
        b += wt*indx_i;      
        c += wt*indx_i*indx_i;
        e += wt*( pk[i].ipos - ref_pos );
        f += wt*indx_i*( pk[i].ipos - ref_pos );
    }

    ok = solve_special_2x2( a, b, c, e, f, shift, spacing );
    if( !ok ) {
#if 0
        fprintf( stderr, 
                 "numerical problem solving 2x2 for spacing and shift\n" );
        fprintf( stderr, "ref_pos=%f a=%f b=%f c=%f e=%f f=%f\n",
                 ref_pos, a, b, c, e, f );
        print_one_peak_header( stderr );
        print_multiple_peaks( stderr, pk, nn );
#endif
        *cost = 0.0;
        return 0;
    }

    *cost = calc_peak_cost( pk, nn, *spacing, *shift );
    return 1;
}

/*******************************************************************************
 * Function: set_peak_lattice_sub
 * Purpose: Given a peak spacing 'spacing', and a shift (w.r.t. pk[0]),
 * align the peaks onto a regular lattice.
 *******************************************************************************
 */
static void
set_peak_lattice_sub( Peak pk[], int nn, double spacing, double shift )
{
    int i;
    double ref_pos = pk[0].ipos;
    for( i=0; i<nn; i++ ) {
        double temp = ref_pos + shift + spacing * i;
        pk[i].ipos = ROUND(temp);
        pk[i].ipos = QVMAX(pk[i].ipos,pk[i].beg);
        pk[i].ipos = QVMIN(pk[i].ipos,pk[i].end);
    }
#if 0   /* NOTICE: it might make sense to include the following: */
    for( i=0; i<nn-1; i++ ) {
        double temp = 0.5*(pk[i].pos+pk[i+1].pos);
        pk[i].end = pk[i].iend = pk[i+1].beg = pk[i+1].ibeg = ROUND(temp);
    }
#endif
}

/*******************************************************************************
 * Function: set_peak_lattice
 * Purpose: Calculate the relevant lattice parameters and then
 * align the peaks onto a regular lattice.
 *******************************************************************************
 */
static void
set_peak_lattice( Peak pk[], int b_indx, int e_indx, double spacing_in )
{
    double spacing=spacing_in, shift, cost;
    Peak* pr_pk = pk+b_indx;        /* poorly resolved peaks */
    int ok, num_peaks = (e_indx-b_indx+1);
    const int strategy=1;

    if( num_peaks==1 ) { return; }
    if( strategy==1 ) { /* calculate spacing */
        ok = calc_peak_spacing_shift_cost( pr_pk, num_peaks, 
                                           &spacing, &shift, &cost );
    } else {            /* just use spacing given as input */
        ok = calc_peak_shift_and_cost( pr_pk, num_peaks, 
                                       spacing, &shift, &cost );
    }
    if( ok ) {
        set_peak_lattice_sub( pr_pk, num_peaks, spacing, shift );
    }
    /* fprintf( stderr, "(%d->%d [%g,%g])", b_indx, e_indx, shift, cost ); */
}


#define END_STATE -1
#define BEGIN_STATE 0
#define MIDDLE_STATE 1

/*******************************************************************************
 * Function: fix_peak_positions
 * Purpose: Loop over the peaks within a group so as to:
 *          (1) Break this group into smaller "poorly resolved" groups.
 *          (2) Align each of these smaller groups onto a regular lattice.
 *******************************************************************************
 */
static void
fix_peak_positions( Peak pk[], int peak_list_len, const int data_array[], 
                    int nn, double ave_spacing, BtkMessage* message )
{
    int indx, state = END_STATE, b_indx=0;

    for( indx=0; indx<nn; indx++ ) {
        int
            first = is_first_poorly_resolved_peak_sub( pk+indx, peak_list_len,
                                                       data_array, message ),
            last  = is_last_poorly_resolved_peak_sub( pk+indx, peak_list_len,
                                                      data_array, message );

        switch (state) {

        case END_STATE:
            if( last ) {                /* Process current peak as singlet.*/
                state = END_STATE;
            } else if ( first ) {       /* Start new group.*/
                b_indx = indx;
                state = BEGIN_STATE;
            } else {                    /* Process current peak as singlet.*/
                state = END_STATE;
            }
            break;

        case BEGIN_STATE:
            if( last ) {                /* Process current group.*/
                set_peak_lattice( pk, b_indx, indx, ave_spacing );
                state = END_STATE;
            } else if ( first ) {       /* Process last peak as singlet.*/
                b_indx = indx;          /* Start new group. */
                state = BEGIN_STATE;
            } else {                    /* Continue current group. */
                state = MIDDLE_STATE;
            }
            break;

        case MIDDLE_STATE:
            if( last ) {                /* Process current group */
                set_peak_lattice( pk, b_indx, indx, ave_spacing );
                state = END_STATE;
            } else if ( first ) {       /* Process old group */
                                        /* and start new group. */
                set_peak_lattice( pk, b_indx, indx-1, ave_spacing );
                b_indx = indx;
                state = BEGIN_STATE;
            } else {                    /* Continue current group */
                state = MIDDLE_STATE;
            }
            break;

        default:
            fprintf( stderr, "bogus state= %d\n", state );
            break;
        }
    }
}

#undef END_STATE 
#undef BEGIN_STATE
#undef MIDDLE_STATE

#if 0
static void
print_peak_info( Peak pk[], int peak_list_len, int data_array[], int nn,
                 double ave_spacing, BtkMessage* message )
{
    int i;

    fprintf( stderr, "as=%g ", ave_spacing );
    for( i=0; i<nn; i++ ) {
        fprintf( stderr, "(%2d,%5d)", pk[i].type, pk[i].ipos );
    }

    fprintf( stderr, " : " );

    for( i=0; i<nn; i++ ) {
        int first, last, flag;
        first = is_first_poorly_resolved_peak_sub( pk+i, peak_list_len,
                                                   data_array, message );
        last  = 2*is_last_poorly_resolved_peak_sub( pk+i, peak_list_len,
                                                    data_array, message );
        flag = first | last;
        fprintf( stderr, " %d", flag );
    }
    /*fprintf( stderr, "\n" );*/
    fprintf( stderr, " : " );
}
#endif

/*******************************************************************************
 * Function: get_intrinsic_position_and_height
 * Purpose:  Utility routine for resolve_sub().
 *           Given an array of deviations, estimate:
 *            (1) intrinsic peak position, (2) intrinsic peak height
 *           by finding the maximum value of the array.
 *******************************************************************************
 */
static void
get_intrinsic_position_and_height( double array[], int beg, int end, 
    double beta, double *ipos, double *iheight, int debug, Options *options )
{
    int ctr, i;         // ctr = counter of points with positive signal 
    int width = end - beg;
    double sum_exp = 0., sum_exp2 = 0., x;
    double B = 1./4./beta/beta;

    if( debug ) fprintf( stderr, "pmax:" );
   *ipos = (double)(width-1)*0.5;
   *iheight = 0.0;

    for(ctr=0,i=0; i<width; i++ ) 
    { 
        if( array[i]>0.0 ) ctr++; 
    }

    if( ctr && width>1 ) {
        // Determine peak position
       *iheight = -1.0;
        for( i=0; i<width; i++ ) {
            if( debug ) fprintf( stderr, " %2d,%8.2f", i, array[i]);
            if( DBL_GT_DBL( array[i], *iheight ) ) {
                if(debug) fprintf( stderr, "_" );
               *ipos = i;
               *iheight = array[i];
            }
        }
        i = *ipos;
//       fprintf(stderr, "Initial ipos= %f\n", *ipos);
       *ipos = update_ipos(B, i,           array[i], 
                              i-1>0?i-1:0, i-1>0?array[i-1]:array[0],
                              i+1<width-1?i+1:width-1, 
                              i+1<width-1?array[i+1]:array[width-1]); 
        if (*ipos <= 0) *ipos = 1.;
        if (*ipos >= width-1) *ipos = width - 1.;

        // Determine intrinsic height
        for( i=0; i<width; i++ ) {
            x = *ipos - (double)i;
            sum_exp  += (double)(array[i])*Shape(1., 0.0, beta, x, options);
            sum_exp2 += Shape(1., 0.0, beta, x, options)
                       *Shape(1., 0.0, beta, x, options);
//          if (beg==65)
//              fprintf(stderr, "beg=%d width=%d ipos=%f beta=%f x= %f shape=%f sum_exp=%f sum_exp2=%f\n", beg, width, *ipos, beta, x, Shape(1., 0.0, beta, x, options), sum_exp, sum_exp2);
        }
//      fprintf(stderr, "   array= ");
//      for( i=0; i<width; i++ ) 
//          fprintf(stderr, "%f ", array[i]);
//      fprintf(stderr, "\n");
//      fprintf(stderr, "iheight_old = %f iheight_new = %f\n", *iheight, sum_exp / sum_exp2);
       *iheight = sum_exp / sum_exp2;
    }
    
    // Determine peak iheight
    
    if( debug ) fprintf( stderr, " ipos=%f iheight=%f\n", *ipos, *iheight );
   *ipos += beg;
//  fprintf(stderr, "multiple ipos = %f beg=%d end=%d \n", *ipos, beg, end);
}

/*******************************************************************************
 * Function: compute_max_resolution
 * Purpose:  Utility routine for resolve_sub().
 *           The return value is used as a convergence criterion.
 *******************************************************************************
 */
static double
compute_max_resolution( Peak *peak_list, int n )
{
    int j;
    double
        max_resolution = 0.0,
        resolution = 0.0, 
        area = 0.0; 
    for (j = 0; j < n; j++) { 
        if ((peak_list[j].type < 30    ) && 
            (peak_list[j].type % 10 < 3)) 
        { 
            if (max_resolution < peak_list[j].resolution) { 
                max_resolution = peak_list[j].resolution; 
            } 
            resolution = 0.; 
            area = 0.; 
        } 
        else { 
            resolution += peak_list[j].resolution * peak_list[j].area; 
            area       += peak_list[j].area; 
            if ((peak_list[j].type == 32) || 
                (peak_list[j].type == 31)) 
            { 
                resolution /= area; 
                if (max_resolution < resolution) { 
                    max_resolution = resolution; 
                } 
            } 
        } 
    } 
    return max_resolution;
}


/*******************************************************************************
 * Function: resolve_multiple_peaks & resolve_sub   
 *  
 * Purpose: fit n overlapping peaks in a row, starting from the peak with  
 *          index i, with mathematical model. For each peak, compute resolution,
 *          intrinsic peak boundaries, area, position and height.
 *          Peak boundaries are conditionally changed
 *            (depending upon options->respace)  
 *  
 *******************************************************************************
 */ 
static int
resolve_sub( Peak* peak_list, int n, int peak_list_len, int color, 
             int* data_array, int length, 
             double rep_width2, 
             double rep_w02beta, double rep_spacing, 
             double* max_resolution, Options *options, BtkMessage *message )
{
    int        j, k=0;  
    double     *dev=NULL;
    int  debug=0, fix_option = 1; /* NOTICE */
    double C0_factor, beta_factor, orig_width_factor;
#if 0
    Peak *pk = &(peak_list[0]);
    debug = (color==0 && pk->pos>8590 && pk->pos<8650 );
#endif
    /* print_peak_info( peak_list, peak_list_len, data_array, n, ave_spacing,
                        message ); */

    /* Allocate memory */
    {
        int size = peak_list[n-1].end-peak_list[0].beg + 1000;
        if (size <= 0) { 
            fprintf(stderr,"Tried to malloc %d bytes in resolve_sub\n", size);
            goto error;
        } else {
            /* fprintf( stderr, "CALLOC of %d bytes\n", size ); */
            dev= CALLOC(double, size);
            MEM_ERROR(dev); 
        }
    }

    if( options->gauss ) {
        C0_factor = 1.0;
        beta_factor = M_1_SQRTPI;
        orig_width_factor = 0.0;
    } else {
        beta_factor = 1.0 / TWO / W1(rep_w02beta);
        C0_factor = 1.0 / Erf(rep_w02beta);
        orig_width_factor = TWO*rep_w02beta;
    }

    {
        int ok;
        double *diag, *off_diag, *hites, *ihites;
        double 
            ref_hite = 100.0, /* somewhat arbitrary */
            ref_beta = rep_width2 * beta_factor,
            ref_orig_width = ref_beta * orig_width_factor,
            shape0;

        diag     = (double*)malloc( n*sizeof(double) );
        off_diag = (double*)malloc( (n-1)*sizeof(double) );
        hites    = (double*)malloc( n*sizeof(double) );
        ihites   = (double*)malloc( n*sizeof(double) );
        if ( !diag || !off_diag || !hites || !ihites ) {
            ok = 0;
        } else {
            int k;
            shape0 = Shape( ref_hite*C0_factor, ref_orig_width/TWO,
                            ref_beta, 0.0, options );
            for ( k=0; k<n; k++ ) {
                diag[k]  = 1.0;
                hites[k] = peak_list[k].height;
                if( k<n-1 ) {
                    off_diag[k] = Shape( ref_hite*C0_factor, ref_orig_width/TWO,
                                         ref_beta, 
                                         peak_list[k+1].ipos-peak_list[k].ipos,
                                         options ) / shape0;
                    /*fprintf( stderr, "%6.2f", off_diag[k] );*/
                }
            }
            /*fprintf(  stderr, "\n" );*/
            ok = solve_sym_tridiag( diag, off_diag, hites, ihites, n );
        }

        /* Initialization */
        for (j=0; j<n; j++) { 
            peak_list[j].ave_width2  = rep_width2;
            peak_list[j].beta        = ref_beta;
            peak_list[j].orig_width  = ref_orig_width;
            peak_list[j].iheight = ok ? ihites[j] : peak_list[j].height;
            peak_list[j].C0          = peak_list[j].iheight * C0_factor;
            peak_list[j].ave_w02beta = rep_w02beta; 
//          peak_list[j].ibeg        = peak_list[j].beg; 
            peak_list[j].ipos        = peak_list[j].pos;
            /*fprintf(stderr,"%4d %f\n",
              peak_list[j].height,peak_list[j].iheight );*/
        }
        FREE( diag );
        FREE( off_diag );
        FREE( hites );
        FREE( ihites );
    }

    if( debug ) {
        fprintf( stderr, "fo=%d (begin of resolve)\n", fix_option );
        print_multiple_peaks( stderr, peak_list, n );
    }
    if ( options->respace && (fix_option==3) ) {
        fix_peak_positions( peak_list, peak_list_len, data_array, n, 
                            rep_spacing, message );
        if( debug ) {
            fprintf( stderr, "fo=%d (before loop, after lattice)\n",fix_option);
            print_multiple_peaks( stderr, peak_list, n );
        }
    }

    /* Iterating ...  
     *************** 
     */  
    for (k=0; k < MAX_ITER; k++) 
    {  

        if( debug ) fprintf( stderr, "iteration #%d:\n", k );

        if ( options->respace && (fix_option==1) ) {
            if( debug ) {
                fprintf( stderr, "fo=%d k=%d (top of loop, before lattice)\n",
                         fix_option, k );
                print_multiple_peaks( stderr, peak_list, n );
            }
//          fix_peak_positions( peak_list, peak_list_len, data_array, n, 
//                              rep_spacing, message );
            if( debug ) {
                fprintf( stderr, "fo=%d k=%d (top of loop, after lattice)\n",
                         fix_option, k );
                print_multiple_peaks( stderr, peak_list, n );
            }
        }

        /* Loop through all peaks */ 
        for (j=0; j<n; j++) 
        {
            int
                offset = peak_list[j].beg - peak_list[0].beg,
                size = peak_list[j].end - peak_list[j].beg + 1;
            double pos, hite;
            int debug2 = 0; /* ( debug && j>=2 && j<=4 ); */

            if( offset < 0 || size < 0 ) {
                if (options->Verbose > 1)
                    fprintf( stderr, "resolve_sub: offset=%d size=%d\n", 
                         offset, size );
                goto error;
            }
            /* Compute the deviation of the data from the model.
             * The deviation is the data signal minus (all) the intrinsic peak  
             * signals but the current (j-th) one 
             */     
            compute_deviation( peak_list, data_array, j, n, options, dev );

            get_intrinsic_position_and_height(dev+offset, peak_list[j].beg, 
                peak_list[j].end, peak_list[j].beta, &pos, &hite, debug2,
                options );

//          hite = (hite>MIN_PEAK_IHEIGHT*ONE_PLUS) ? hite : 
//              QVMAX((double)data_array[peak_list[j].ipos]/2., MIN_PEAK_IHEIGHT);
            peak_list[j].ipos    = ROUND(pos);
            peak_list[j].iheight = hite;             
            peak_list[j].C0      = hite * C0_factor;

            peak_list[j].resolution = 
                compute_peak_resolution( &peak_list[j], peak_list[0].beg,
                                         dev, options );

            if (debug) { print_one_peak( stderr, &(peak_list[j]) ); }
        } /* end loop in j */ 
 
        /* Compute max_resolution */ 
        *max_resolution = compute_max_resolution( peak_list, n );

        if ( options->respace && ((fix_option==2)||(fix_option==3)) ) {
            if( debug ) {
                fprintf( stderr,"fo=%d k=%d (bottom of loop, before lattice)\n",
                         fix_option, k );
                print_multiple_peaks( stderr, peak_list, n );
            }
            fix_peak_positions( peak_list, peak_list_len, data_array, n, 
                                rep_spacing, message );
            if( debug ) {
                fprintf( stderr,"fo=%d k=%d (bottom of loop, after lattice)\n",
                         fix_option, k );
                print_multiple_peaks( stderr, peak_list, n );
            }
        }
		/* Stop the iterations if they "converged" */
        if (*max_resolution < MIN_PEAK_RESOLUTION) { 
            break; 
        } 
 
    } /* end loop in k */  

/*  printf("\nTitleText: Peaks %d-%d, max_res=%f, Trace %d, Base %c\n\n\n", 
 *     i, i+n-1, *max_resolution, color, color_to_base[color] ); 
 */ 
 
#if SHOW_MULTIPLE_PEAK_FIT 
    {
        static char color_to_base[] = { 'A', 'C', 'G', 'T' };
        printf(
          "\nTitleText: Peaks %d-%d, max_res=%.4f, Trace %d, Base %c\n",
          peak_list[0].cd_peak_ind, peak_list[n-1].cd_peak_ind,
          *max_resolution, color, color_to_base[color] ); 
        printf("\n\"Data \n");
        for (m=peak_list[0].beg; 
             m<peak_list[n-1].end; m++) 
        { 
            printf("%d   %d\n", m, data_array[m]); 
            dev[m-peak_list[0].beg] = 0; 
        } 
        for (j=0; j<n; j++) {
            int beg, end;
            beg = peak_list[j].beg; 
            end = peak_list[j].end; 
            if (j>0) { 
                beg=peak_list[j-1].beg; 
            } 
            if (j<n-1) { 
                end = peak_list[j+1].end; 
            } 
            printf("\n\"Model Peak %d\n", peak_list[j].cd_peak_ind ); 
            for (m=beg; m<=end; m++) 
            { 
                double x, x1, model;
                x = (double)m; 
                model= Shape(peak_list[j].C0,
                             peak_list[j].orig_width/TWO, 
                             peak_list[j].beta, 
                             (x - (double)peak_list[j].ipos), options ); 
                x1 = x;     
                printf("%f   %f\n", x1, model); 
                dev[m-peak_list[0].beg] += model; 
            } 
        } 
        printf("\n\"Model Overall \n"); 
        for (m=peak_list[0].beg; m<peak_list[n-1].end; m++) { 
            printf("%d   %g\n", m, dev[m-peak_list[0].beg]); 
        } 
    }
#endif 
 
#if SHOW_MULTIPLE_RESOLUTION 
    fprintf(stderr, "_max_res=%.3f ITER=%d\n", 
    *max_resolution, k);
#endif
    if( debug ) {
        fprintf( stderr, "f0=%d (end of resolve)\n", fix_option );
        print_multiple_peaks( stderr, peak_list, n );
    }

    FREE(dev); 
 
    return SUCCESS; 
 
 error: 
    FREE(dev); 
 
    return ERROR;
}

/*******************************************************************************
 * Function: resolve_multiple_peaks & resolve_sub   
 *  
 * Purpose: fit n overlapping peaks in a row, starting from the peak with  
 *          index i, with mathematical model. For each peak, compute resolution,
 *          intrinsic peak boundaries, area, position and height.  
 *          Do not change the apparent peak boundaries.  
 *  
********************************************************************************
 */ 
int 
resolve_multiple_peaks(Data *data, int color, int i, int n,  
    double *max_resolution, Options *options, BtkMessage *message) 
{ 
    ColorData *cd;
    Peak* short_peak_list;
    int*  data_array;
    int   length=0;
    double rep_w02beta, rep_spacing, rep_width1, rep_width2, rep_width_ratio;
#if SHOW_MULTIPLE_RESOLUTION 
    fprintf(stderr, "C=%d Peaks %3d-%3d:", color, i, i+n-1); 
#endif 
 
    cd = &data->color_data[color]; 
    short_peak_list = cd->peak_list+i;
    data_array = cd->data;
    length = cd->length;

    rep_spacing = get_average_abi_spacing( short_peak_list[0].pos,
                                           data, message);

    /* Initiation  
     ************ 
     */ 
    rep_width1 = cd->peak_list[i].ave_width1;
    rep_width2 = cd->peak_list[i].ave_width2;
    rep_width_ratio = rep_width1/rep_width2;                      
    rep_w02beta = options->gauss ? 0.0 : Phi(rep_width_ratio);

    return resolve_sub( short_peak_list, n, cd->peak_list_len,
                        color, data_array, length,
                        rep_width2,
                        rep_w02beta, rep_spacing, max_resolution, 
                        options, message );
}

/*******************************************************************************
 * Function: data_resolve_single_peaks 
 * Purpose: fit each single observed (or apparent) peak with a model,  
 *          determine the model's parameters  
 * Notes:   1. Fit all peaks, including the truncated ones, since 
 *             we need the peak resolution parameter for all peaks 
 *          2. Upon the resolution, peak.ipos will be used instead of peak.pos 
 *******************************************************************************
 */ 
static int 
data_resolve_single_peaks( Data *data, Options *options, BtkMessage *message) 
{ 
    int        i;
    double     stderr_w1, stderr_w2;

    for (i=0; i < data->peak_list_len; i++) {
        Peak* peak = data->peak_list[i];
        peak->ave_width1 = get_average_width(1, peak->data_peak_ind,
            NUM_PEAK_UPSTREAM, NUM_PEAK_DOWNSTREAM,
            &stderr_w1, data, message);
        peak->ave_width2 = get_average_width(2, peak->data_peak_ind,
            NUM_PEAK_UPSTREAM, NUM_PEAK_DOWNSTREAM,
           &stderr_w2, data, message);
        if ((peak->type   ==11        ) || 
            (peak->type%10==1  && i==0) || 
            (peak->type<20     && i==data->peak_list_len-1)) 
        {
            if (fit_single_peak(peak, stderr_w1, stderr_w2,
                data, options, message) != SUCCESS)
                return ERROR;
        }
    }             /* loop in i */ 
    return SUCCESS; 
} 

/*******************************************************************************
 * Function: data_resolve_peaks 
 *******************************************************************************
 */ 
int 
data_resolve_peaks( Data *data, Options *options, BtkMessage *message) 
{ 
    int        i, j, k, n, color; 
    int        shift[NUM_COLORS] = {0, 0, 0, 0};
    double     max_resolution; 
    ColorData *cd; 
    Peak      *pk, *pk_n;

    /* Create a list of all peaks */
    if (bc_data_create_single_ordered_peak_list(data, shift, message) != SUCCESS) 
    {
        sprintf(message->text, "Error creating single peak list\n");
        return ERROR;
    }

    set_peak_widths(data, message);
#if 0
    /* Check setting of width1 and 2 */
    for (j = 0; j < data->peak_list_len; j++)
    {
        if ((data->peak_list[j]->width1 < EPSILON) ||
            (data->peak_list[j]->width1 < EPSILON))
            fprintf(stderr, "Peak%d type=%d width1=%f width2=%f\n",
                 j, data->peak_list[j]->type, 
                 data->peak_list[j]->width1,  data->peak_list[j]->width2);
    }
#endif
 
    if( data_resolve_single_peaks( data, options, message ) != SUCCESS ) {
        return ERROR;
    }
    
    /* print_peak_list( stderr, data, 0, message ); */

    /* Scan the list of all peaks, find multiple peaks and resolve them */ 
    for (j = 0; j < data->peak_list_len; j++) 
    { 
 
        /* Set ave_width_ratio and ave_width2*/ 
        color = data->peak_list[j]->color_index; 
        cd = &data->color_data[color]; 
        i = data->peak_list[j]->cd_peak_ind; 
        pk = &cd->peak_list[i];

        /* Skip already resolved peaks */ 
        if (data->peak_list[j]->resolution >0) 
            continue; 
         
        /* Skip peaks of type != 12 and != 13 */ 
        if ((                             pk->type  ==11)  || 
            ((i==0                  ) && (pk->type%10==1)) || 
            ((i==cd->peak_list_len-1) && (pk->type   <20)))  
        { 
            continue; 
        } 
 
        /* Determine the cd_peak_ind, i+k-1, of the last peak  
         * in the row of overlapping peaks  
         */ 
        k = 1; 
        for (n=1; j+n < data->peak_list_len; n++) { 
            if (data->peak_list[j+n]->color_index != color) { 
                continue; 
            } 
            k++; 
            if ((data->peak_list[j+n]->type%10 == 1)              || 
                (i+k-1 == data->color_data[color].peak_list_len-1)) 
            { 
                break; 
            } 
        } 
 
        /* Resolve peaks of a given color with indices i, i+1, ...,  i+k-1 */ 
        for (n=0; n<k; n++) {
            pk_n = &cd->peak_list[i+n];
            pk_n->ave_width1 = pk->ave_width1;
            pk_n->ave_width2 = pk->ave_width2;                
        }
        if (resolve_multiple_peaks(data, color, i, k, &max_resolution,  
            options, message) != SUCCESS)  
            return ERROR; 
    } 
#if 0
    for (j = 0; j < data->peak_list_len; j++) {
        if (data->peak_list[j]->resolution < 0.)
        fprintf(stderr, "Peak %d color=%d cd_peak_ind=%d unresolved\n",
                j, data->peak_list[j]->color_index, 
                data->peak_list[j]->cd_peak_ind);   
    }
#endif
    return SUCCESS; 
} 

/*******************************************************************************
 * Function: Btk_process_peaks 
 * Purpose: create peak list, resolve peaks, compute apparent and intrinsic 
 *          peak characteristics 
 *******************************************************************************
 */ 
int 
Btk_process_peaks(Data *data, Options *options, BtkMessage *message)  
{ 
    clock_t start_clock = clock(), curr_clock;

    if (data_detect_peaks(data, options, message) != SUCCESS) {
#if 0
        fprintf(stderr, "Error calling data_detect_peaks\n");
#endif
        goto error;
    }

    /* print_peak_list( stderr, data, 3, message ); */

    /* Check the number of peaks in traces */
    {
        int color;
        for (color = 0; color < NUM_COLORS; color++) {
            ColorData *cd = &data->color_data[color];
            if (cd->peak_list_len <= 10) {
                fprintf(stderr, "Note: number of peaks in trace %d is %d\n",
                    color, cd->peak_list_len);
            }
        }
        if (data->color_data[0].peak_list_len +
            data->color_data[1].peak_list_len +
            data->color_data[2].peak_list_len +
            data->color_data[3].peak_list_len < 1) {
            fprintf(stderr, "Total number of peaks < 1. Exit\n");
            return SUCCESS;   
        }    
    }

    if (options->time) {
        curr_clock = clock();
        fprintf(stderr, "   Peaks detected in %f sec. \n",
            (float)(curr_clock - start_clock)/(float)CLOCKS_PER_SEC);
        start_clock = curr_clock;
    }         

    if (data_expand_peaks(data, options, message) != SUCCESS) {
#if 0
        fprintf(stderr, "Error calling data_expand_peaks\n");
#endif
        goto error;
    }

    if (options->time) {  
        curr_clock = clock();
        fprintf(stderr, "   Peaks expanded in %f sec. \n",
            (float)(curr_clock - start_clock)/(float)CLOCKS_PER_SEC);
        start_clock = curr_clock;
    }  

    if (data_resolve_peaks(data, options, message) != SUCCESS) {
#if 0
        fprintf(stderr, "Error calling data_resolve_peaks\n");
#endif
        goto error;
    } 
   
    if( 0 ) {   /* debug stuff */
        if( test_broken_peak_lists( stderr, data, message ) ) {
            fprintf( stderr, "broken peaks after resolve\n" );
            return (-1);
        }
    }

    if (options->time) {  
        curr_clock = clock();
        fprintf(stderr, "   Peaks resolved in %f sec. \n",
            (float)(curr_clock - start_clock)/(float)CLOCKS_PER_SEC);
        start_clock = curr_clock;
    }   

    /* print_peak_sub_list( stderr, data, 3, 210, 6, message ); */

    return SUCCESS;

    error:
    data_release(data);

    return ERROR;
}

