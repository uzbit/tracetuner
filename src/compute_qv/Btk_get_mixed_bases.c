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
 * $Id: Btk_get_mixed_bases.c,v 1.14 2009/01/12 22:15:11 gdenisov Exp $                      
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <stdint.h>

#include "Btk_qv.h"
#include "SFF_Toolkit.h"
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
#include "Btk_get_mixed_bases.h"
#include "Btk_compute_tpars.h"  /* needs train.h */
#include "Btk_call_bases.h"
#include "Btk_match_data.h"
#include "Btk_qv_io.h"
#include "Btk_qv_funs.h"
#include "Btk_atod.h" 
#include "Btk_process_raw_data.h"

#define BASE_MERGE_FACTOR 0.85
#define INSERT_BASES     1
#define MERGE_BASE_CALLS 0
#define MIN_SIGNAL_VALUE 20
#define NUM_PEAK_AVERAGE 10
#define PERFORM_SECOND_LOOP 1
        /* Disabled 2nd loop temporarily, because it swaps bases;
         * for example see file ~/cvs/qv/rel/x86-sunos-5.7/Jun_unordered_bases_May10.ab1,
         * near position 336
         */
#define SHOW_ALL_BASES     0
#define SHOW_CALLED_BASES  0
#define SHOW_CASE          0
#define SHOW_INPUT_OPTIONS 0
#define SHOW_MARKED        0   
#define SHOW_MIXED_BASES   0
#define SHOW_PEAK_SEARCH   0
#define STORE_IS_RESOLVED  1
#define STORE_CASE         0
#define USE_BEST_BASE_POS  0
#define WINDOW_7 7

extern unsigned int max_colordata_value;

/*******************************************************************************
 * Function: get_mixed_base_position
 * Purpose:  for mixed base, return the position of the peak with the "worst" 
 *           spacing. By this, we mean the one of the two mixed base's peaks
 *           which has the highest value of the peak spacing ratio
 *******************************************************************************
 */
int
get_mixed_base_position(Data *data, int data_peak_ind, int data_peak_ind2,
    Options *options, BtkMessage *message)
{
    if (data_peak_ind == data_peak_ind2) {
        return data->peak_list[data_peak_ind]->ipos;
    }
    else if ((data_peak_ind==0) || (data_peak_ind2==0) ||
        (data_peak_ind==data->peak_list_len-1) || 
        (data_peak_ind2==data->peak_list_len-1)) {
        return (data->peak_list[data_peak_ind ]->ipos +
                data->peak_list[data_peak_ind2]->ipos)/2;
    }
    else {
#if USE_BEST_BASE_POS
        /* Find indexes of preceeding and following called peak */
        int pos = data->peak_list[data_peak_ind ]->ipos;
        int pos2= data->peak_list[data_peak_ind2]->ipos;
        int base_index = QVMAX(data->peak_list[data_peak_ind ]->base_index,
                               data->peak_list[data_peak_ind2]->base_index);
        if ((base_index==0) || (base_index==data->bases.length-1)) 
            return (pos + pos2)/2;
        else {
            int prev_pos = data->bases.called_peak_list[base_index-1]->ipos;
            int next_pos = data->bases.called_peak_list[base_index+1]->ipos;
            return (abs((prev_pos+next_pos)/2-pos) < 
                    abs((prev_pos+next_pos)/2-pos2)) ? pos : pos2;
        }
#endif
#if !USE_BEST_BASE_POS
        return (data->peak_list[data_peak_ind ]->ipos +
                data->peak_list[data_peak_ind2]->ipos)/2;
#endif
    }
}

/*******************************************************************************
 * Function: mixed_base 
 * Purpose:  convert two regular input base characters into a mixed base IUB 
 *           character
 *******************************************************************************
 */
char
mixed_base(char b1, char b2)      
{
    if      (b1 == b2) return b1;
    else if (((b1 == 'A') && (b2 == 'C')) || ((b1 == 'C') && (b2 == 'A'))) 
        return 'M';
    else if (((b1 == 'A') && (b2 == 'G')) || ((b1 == 'G') && (b2 == 'A'))) 
        return 'R';
    else if (((b1 == 'A') && (b2 == 'T')) || ((b1 == 'T') && (b2 == 'A'))) 
        return 'W'; 
    else if (((b1 == 'C') && (b2 == 'G')) || ((b1 == 'G') && (b2 == 'C'))) 
        return 'S'; 
    else if (((b1 == 'C') && (b2 == 'T')) || ((b1 == 'T') && (b2 == 'C'))) 
        return 'Y'; 
    else if (((b1 == 'G') && (b2 == 'T')) || ((b1 == 'T') && (b2 == 'G'))) 
        return 'K'; 
    
        return 'N';
}

/*******************************************************************************
 * Function: is_mixed_base
 * Purpose:  determine if base is mixed
 *******************************************************************************
 */
int 
is_mixed_base(char b)
{
    if   ((b == 'S') || (b == 'K') || (b == 'M') || 
          (b == 'W') || (b == 'R') || (b == 'Y'))
        return 1;

    return 0;
}

/*******************************************************************************
 * Function: get_base
 *******************************************************************************
 */
static char
get_base(Data *data, int data_peak_ind, int data_peak_ind2, char *color2base)
{
    char b1 = color2base[data->peak_list[data_peak_ind]->color_index],
         b2 = color2base[data->peak_list[data_peak_ind2]->color_index];

    char b = mixed_base(b1, b2);

    return b;
}

/*******************************************************************************
 * Function: store_alternative_base_call
 ******************************************************************************* 
 */
static int 
store_alternative_base_call(Data *data, char *color2base, int base_index, 
    int qv, int data_peak_ind, int data_peak_ind2, int *curr_num2, 
    int *max_num2, AltBase **altbases, Options *options, BtkMessage *message)
{
    int temp_base_index;

     /* If current alternative base call coinsides with previous one 
      * for the same base index, but has higher quality, delete the previous 
      * record before recording the current alternative base call
      */
    if (( *curr_num2 >= 1)                                    &&
        ((*altbases)[*curr_num2-1].base_index == base_index)  &&
        ((*altbases)[*curr_num2-1].base       ==
            get_base(data, data_peak_ind, data_peak_ind2,
                color2base))                                  &&
        ((*altbases)[*curr_num2-1].qv       < qv))
        (*curr_num2)--;

    if (( *curr_num2 > 0)                                     &&
        ((*altbases)[*curr_num2-1].base_index == base_index)  &&
        ((*altbases)[*curr_num2-1].base       ==
             get_base(data, data_peak_ind, data_peak_ind2, 
                 color2base)))
        return SUCCESS;

    temp_base_index = data->peak_list[data_peak_ind]->base_index;
    data->peak_list[data_peak_ind]->base_index = base_index;

   (*altbases)[*curr_num2].base = get_base(data, data_peak_ind,
       data_peak_ind2, color2base);
   (*altbases)[*curr_num2].pos = get_mixed_base_position(data, data_peak_ind, 
       data_peak_ind2, options, message);
   data->peak_list[data_peak_ind]->base_index = temp_base_index;
   (*altbases)[*curr_num2].qv = qv;
   (*altbases)[*curr_num2].base_index = base_index;
   (*altbases)[*curr_num2].data_peak_ind = data_peak_ind;
   (*altbases)[*curr_num2].data_peak_ind2= data_peak_ind2;
   (*curr_num2)++;

    /* Reallocate memory if the actual size of arrays of alternative
     * base calls approaches limit
     */
    if (*curr_num2 >= *max_num2) {
        *max_num2 *= 2;
        *altbases = REALLOC(*altbases, AltBase, *max_num2);
    }
    return SUCCESS;
}

/*******************************************************************************
 * Function: get_average_called_peak_spacing 
 * Purpose:  return the average height of called peak; averaging should
 *           be done among NUM_PEAK_AVERAGE called peaks immediately
 *           preceding the current peak
 *******************************************************************************
 */
static double
get_average_called_peak_spacing(Data *data, int i)
{

    int    j=0, k=i;
    double sum=0.;

    if ((i < 0) || (i>data->bases.length-1)) {
        fprintf(stderr,
            "Index i=%d passed to get_average_called_peak_spacing ", i);
        fprintf(stderr, "is ouside of bounds [0, %d]\n", data->bases.length-1);
        return ERROR;
    }

    if ((data->bases.called_peak_list[k] != NULL) && (k==0))
        return DEFAULT_PEAK_SPACING;

    while ((k>0) && (j<NUM_PEAK_AVERAGE)) {
        if ((data->bases.called_peak_list[k  ] == NULL) ||
            (data->bases.called_peak_list[k-1] == NULL))
            continue;

        sum += (data->bases.called_peak_list[k  ]->ipos -
                data->bases.called_peak_list[k-1]->ipos);
        k--;
        j++;
    }

    if (j==0) return DEFAULT_PEAK_SPACING;

    return sum/(double)j;

}

/*******************************************************************************
 * Function: get_average_called_peak_height
 * Purpose:  return the average height of called peak; averaging should
 *           be done among NUM_PEAK_AVERAGE called peaks immediately 
 *           preceding the current peak
 *******************************************************************************
 */
double
get_average_called_peak_height(Data *data, int i)
{

    int    j=0, k=i;
    int    left_ind = (i - 5 >= 0) ? (i - 5) : 0;
    int    right_ind =(i + 5 <= data->bases.length -1) ? (i + 5) :
                                data->bases.length -1;
    int    num = right_ind - left_ind + 1;
    double min = 10000., second_min = 10000.;
    double max= 0.,     second_max= 0.;
    double sum = 0.;

    if ((i < 0) || (i>data->bases.length-1)) {
        fprintf(stderr, 
            "Index i=%d passed to get_average_called_peak_height ", i);
        fprintf(stderr, "is ouside of bounds [0, %d]\n", data->bases.length-1);
        return ERROR;
    }

    if (data->bases.called_peak_list[k]==NULL)
        return 1.;

    for (j = left_ind; j <= right_ind; j++)
    {
        double iheight = data->bases.called_peak_list[j]->iheight;
        if (max < iheight) 
        { 
            second_max = max;
            max = iheight; 
        }
        else if (max >= iheight && second_max < iheight)
        {
            second_max = max;
        }

        if (min  > iheight)
        {
            second_min  = min ;
            min  = iheight;
        }
        else if (min  <= iheight && second_min  > iheight)
        {
            second_min  = min ;
        }

        if (data->bases.called_peak_list[k] == NULL) 
            continue;

        sum += iheight;    
    }

    return (sum - max - second_max - min - second_min)/
        (double)(num - 4);
}

/*******************************************************************************
 * Function:get_constraints
 *******************************************************************************
 */
static void
get_constraints(int *MIN_QV, double *MAX_PHR3, double *MAX_PHR7, 
    double *MAX_PSR7, double *MAX_PRES, Options *options)
{
   *MIN_QV   = 10;
   *MAX_PHR3 = 4.5;
   *MAX_PHR7 = 29.; 
   *MAX_PSR7 = 100.;
   *MAX_PRES = 0.55;
}

/*******************************************************************************
 * Function: get_bounds_of_the_search_resion
 *******************************************************************************
 */
static int 
get_right_bound_of_search_region(Data *data, int i)
{
    int pos = data->bases.called_peak_list[i  ]->ipos, right_bound,
        next_pos = (i < data->bases.length - 1) ?
              data->bases.called_peak_list[i+1]->ipos : data->length-1,
        right_height = (i < data->bases.length - 1) ?
            INT_DBL(data->bases.called_peak_list[i+1]->iheight) : INF;
  
    if (0)
        right_bound = (i < data->bases.length - 1) ?
           INT_DBL((double)pos + (double)(next_pos - pos)/
            (right_height + data->bases.called_peak_list[i]->iheight) *
            right_height +1) : pos;
    else 
        right_bound = INT_DBL((pos + next_pos) / 2.) + 2;

    return right_bound;
}

/*******************************************************************************
 * Function: find alternative_peaks
 *           for each original base call, alternative
 *           peaks are searched for between the left bound and right bound
 ******************************************************************************* 
 */
static int
find_alternative_peaks(Data *data, int i, int *found_peak_ind,
    int *found_peak_is_called, int *num_peaks_found, 
    int *next_called_peak_found, int *prev_right_bound, double max_pres, 
    double min_peak_height, Options *options, BtkMessage *message) 
{
    int        j, jc, k, pind, pos[3]; 
    static int left_bound=0, right_bound=0;
    double     iheight = data->bases.called_peak_list[i]->iheight,
               ave_spacing;
    Peak       peak;
    Peak     **cpl = data->bases.called_peak_list;

    /* Setting the left bound of the search region */
    left_bound = (i > 0) ? 
        QVMAX(get_right_bound_of_search_region(data, i-1)-4, 0) : 0;
    if ((i>0) && (left_bound <= cpl[i-1]->ipos))
         left_bound  = cpl[i-1]->ipos+1;
    if (left_bound < *prev_right_bound) 
        left_bound = *prev_right_bound;

    /* Setting the right bound of the search region */
    right_bound = get_right_bound_of_search_region(data, i);
    if ((i<data->bases.length-1) && (right_bound >= cpl[i+1]->ipos)) 
         right_bound  = cpl[i+1]->ipos;  
   *prev_right_bound = right_bound;

    if (SHOW_ALL_BASES)
        fprintf(stderr, "    left_bound=%d right_bound=%d\n", 
            left_bound, right_bound);

    /* Set the three locations of search */
    pos[0] = left_bound;
    pos[1] = data->bases.called_peak_list[i]->ipos;
    pos[2] = right_bound;

    /* Originally called peak is the first found one */
    found_peak_ind[0] = cpl[i]->data_peak_ind;
    found_peak_is_called[0] = 1;
    jc = data->bases.called_peak_list[i]->color_index;
   *num_peaks_found = 0;
   *next_called_peak_found = 0;     
    ave_spacing = get_average_called_peak_spacing(data, i); 

    /* Find all possible alternative peaks by search through
     * four colors and three locations */
    for (j=0; j<NUM_COLORS; j++) {
        for (k=0; k<3; k++) {
            if (j == jc)
                continue;

            if (data->color_data[j].data[pos[k]] <= MIN_SIGNAL_VALUE)
                continue;

            pind = -1;
            if (colordata_find_peak_index_by_location(&(data->color_data[j]),
                pos[k], &peak, &pind, message) != SUCCESS) {
                return ERROR;
            }
           
            if (SHOW_PEAK_SEARCH) {
                if (pind < 0) 
                    fprintf(stderr, 
                    "    i=%d j=%d k=%d search_loc=%d pind=%d\n",
                    i, j, k, pos[k], pind);
                else {
                    fprintf(stderr,
               "    i=%d j=%d k=%d search_loc=%d pind=%d base=%c ipos=%f res=%f\n",
                    i, j, k, pos[k], pind, peak.base, peak.ipos, peak.resolution);
                    fprintf(stderr,
               "         is_called=%d is_blob=%d data_peak_ind=%d found_peak_ind=%d\n",
                    peak.is_called, is_dye_blob(peak.ipos, &peak, data, 1), 
                       peak.data_peak_ind, found_peak_ind[*num_peaks_found]);
                    fprintf(stderr,
                    "         left_bound=%d right_bound=%d ratio=%f\n",
                    left_bound, right_bound, 
                    (double)peak.iheight/(double)iheight);
                }
            }

            if (pind < 0)
                continue;

            if (is_dye_blob(peak.ipos, &peak, data, 1))
                continue;

            if (peak.resolution > max_pres)
                continue;

            if (peak.iheight < min_peak_height)
                continue;

            if ((num_peaks_found > 0) &&
                (peak.data_peak_ind == found_peak_ind[*num_peaks_found]))
                continue;
 
            if ( (peak.is_called == 0) &&
                ((peak.ipos < left_bound) || (peak.ipos >= right_bound)))
                continue;

             /* Seach second called peak only at the middle position
              * (k==1), that is, only if it "overlaps" with
              * original called peak
              */
            if ( (peak.is_called > 0) &&
                ((peak.base_index >= i+2) || (k != 1) ||
                 (peak.base_index < i)    || (peak.ipos < pos[1]) ||
                ( INT_GT_DBL( (peak.ipos-pos[k]), ave_spacing/2.) )))
                continue;

            if (peak.iheight / iheight < options->min_ratio)
                continue;

            /* Don't consider next called peak as potential
             * candidate for merge with current called peak if 
             * spacing considerations suggest they should not merge
             */
            if ((peak.is_called > 0) &&
                (i>0) &&
                (i<data->bases.length-1) &&
                 ( can_insert_base(data, i-1, i+1, 1,
                   spacing_curve(data->bases.called_peak_list[i]->ipos) *
                   BASE_MERGE_FACTOR, -1., options) ||
                  ((i<data->bases.length-2) &&
                   can_insert_base(data, i, i+2, 1,
                   spacing_curve(data->bases.called_peak_list[i+1]->ipos) *
                   BASE_MERGE_FACTOR, -1., options))))
                continue;

            if (peak.is_called > 0)
               *next_called_peak_found = 1;
    
          /* Store the found peak index */       
          (*num_peaks_found)++;
            found_peak_ind[*num_peaks_found] =
                peak.data_peak_ind;
            if (peak.is_called)
                found_peak_is_called[*num_peaks_found] = 1;
            else
                found_peak_is_called[*num_peaks_found] = 0;

            if (SHOW_ALL_BASES) {
                fprintf(stderr,
                    "    Found new peak %c is_called=%d data_peak_ind=%d",
                    data->color_data[j].peak_list[pind].base,
                    data->color_data[j].peak_list[pind].is_called,
                    data->color_data[j].peak_list[pind].data_peak_ind);
                fprintf(stderr," ipos=%f iheight=%.1f is_dye_blob=%d\n",
                    data->color_data[j].peak_list[pind].ipos,
                    peak.iheight, is_dye_blob(peak.ipos, &peak, data, 1));
            }
        } /* loop through 3 locations */
    } /* loop through 4 colors; now all new peaks should be found */
  
    if (SHOW_ALL_BASES)
        fprintf(stderr, "    num_peaks_found=%d\n", *num_peaks_found);
 
    return SUCCESS; 
}

/*******************************************************************************
 * Function: get_mixed_qv
 *******************************************************************************
 */
static int
get_mixed_qv(int qv1, int qv2) 
{
    int qv;
    double prob1, prob2, prob;
    
    prob1 = 1. - pow(10., -(double)qv1/10.);
    prob2 = 1. - pow(10., -(double)qv2/10.);
    prob  = prob1 * prob2;
    qv = -10. * log10(1. - prob);
    
    return qv;
}

/*******************************************************************************
 * Function: get_pure_base_quality
 * Purpose:  return a quality value of a regular or mixed base
 * Note:     if data_peak_ind2 >= 0, then this peak will be ignored
 *           when computing trace parameters
 *******************************************************************************
 */
int
get_pure_base_quality(Data *data, int data_peak_ind, int data_peak_ind2,
    double frac, double min_uncalled, BtkLookupTable *table, ReadInfo *read_info, 
    Options *options, BtkMessage *message)
{
    double phr3, phr7, psr7, pres;

    phr3 = get_peak_height_ratio(data, data_peak_ind, data_peak_ind2, 3, frac,
        min_uncalled, options, message);
    phr7 = get_peak_height_ratio(data, data_peak_ind, data_peak_ind2, 7, frac,
        min_uncalled, options, message);
    psr7 = get_peak_spacing_ratio(data, data_peak_ind, data_peak_ind2, 7,
        options, message);
    pres = get_peak_resolution(data, data_peak_ind, options, message);

    return get_quality_value(phr3, phr7, psr7, pres, table);
}

/*******************************************************************************
 * Function: get_quality_of_alternative_call
 *           i    - base index of called base
 *           j, k - data indexes of called peak(s) 
 *******************************************************************************
 */
static int
get_quality_of_alternative_call(int i, int j, int k, int *found_peak_ind, 
    int *found_peak_is_called, Data *data, BtkLookupTable *table, 
    ReadInfo *read_info, Options *options, BtkMessage *message) 
{
    int    qv = 0, qv1 =0, qv2 = 0;
    Peak **dpl = data->peak_list;
    

    /* Pure base call */                  
    if (found_peak_ind[j] == found_peak_ind[k]) {
        dpl[found_peak_ind[j]]->is_called = 1;
        dpl[found_peak_ind[j]]->base_index= i;

        qv = get_pure_base_quality(data, found_peak_ind[j], -1, 1.,
            0., table, read_info, options, message); 
    }

    /* Mixed base call */
    else 
    {
        double h1 = data->peak_list[found_peak_ind[j]]->wiheight;
        double h2 = data->peak_list[found_peak_ind[k]]->wiheight;
        double h  = get_average_called_peak_height(data, i);
        double delta = 0.;

        dpl[found_peak_ind[j]]->is_called = 1;
        dpl[found_peak_ind[j]]->base_index= i;
        dpl[found_peak_ind[k]]->is_called = 0;
        dpl[found_peak_ind[k]]->base_index= -1;

        /* Determine QV of pure central base in the 1-st descending trace */
        if (options->het) {
            delta = h1 > h2 ? 0. : (h2 - h1);
            if (h1 > h/2. && h2 > h/2.) { delta = h2 - h/2; }
            qv1 = get_pure_base_quality(data, found_peak_ind[j], found_peak_ind[k],
            0.5, delta, table, read_info, options, message);
        }
        else if (options->mix) {
            double ratio = h1/h2;
            qv1 = get_pure_base_quality(data, found_peak_ind[j], found_peak_ind[k], 
            1./(1.+ratio), delta, table, read_info, options, message);
        }

        dpl[found_peak_ind[k]]->is_called = 1;
        dpl[found_peak_ind[k]]->base_index= i;
        dpl[found_peak_ind[j]]->is_called = 0;
        dpl[found_peak_ind[j]]->base_index= -1;

        if (options->het) {
            delta = h1 > h2 ? (h1 - h2) : 0.;
            if (h1 > h/2. && h2 > h/2.) { delta = h1 - h/2; }
            qv2 = get_pure_base_quality(data, found_peak_ind[k], found_peak_ind[j],
            0.5, delta, table, read_info, options, message);
        }
        else if (options->mix) {
            double ratio = h1/h2;
            qv2 = get_pure_base_quality(data, found_peak_ind[k], found_peak_ind[j],
            ratio/(1.+ratio), delta, table, read_info, options, message);
        }
        qv = get_mixed_qv(qv1, qv2);
    }
    /* Unset alternative base call */
    dpl[found_peak_ind[j]]->is_called = 0;
    dpl[found_peak_ind[k]]->is_called = 0;
    dpl[found_peak_ind[j]]->base_index = -1;
    dpl[found_peak_ind[k]]->base_index = -1;
    if ((k>0) && found_peak_is_called[k]) {
        dpl[found_peak_ind[k]]->is_called = 1;
        dpl[found_peak_ind[k]]->base_index= i+1;
    }

    return qv;
}

/*******************************************************************************
 * Function: set_best_base_call
 *******************************************************************************
 */
static int
set_best_base_call(int i, int dpi_best, int dpi2_best, int *data_peak_ind1, 
    int *data_peak_ind2, char *color2base, Data *data, Options *options) 
{
    char b1, b2, base;
    Peak **dpl = data->peak_list;

    data->bases.called_peak_list[i] = dpl[dpi_best];
    dpl[dpi_best ]->is_called  = 1;
    dpl[dpi2_best]->is_called  = 1;
    dpl[dpi_best ]->base_index = i;
    if (dpi_best != dpi2_best)
        dpl[dpi2_best]->base_index =-1;
    dpl[dpi_best ]->data_peak_ind = dpi_best;
    dpl[dpi_best ]->data_peak_ind2= dpi2_best;
    dpl[dpi2_best]->data_peak_ind = dpi2_best;
    dpl[dpi2_best]->data_peak_ind2= dpi_best;

    b1=color2base[dpl[dpi_best ]->color_index];
    b2=color2base[dpl[dpi2_best]->color_index];
    base = b1;
    if (!options->poly) {
        base = mixed_base(b1, b2);
        data->bases.bases[i] = base;
        dpl[dpi_best ]->base = b1;
        dpl[dpi2_best]->base = b2;
    }
    else { 
        dpl[dpi_best ]->base = data->bases.bases[i];
        dpl[dpi2_best]->base = data->bases.bases[i];
    }

    /* Warning: we accept the following convention:
     * if we call mixed base, then the mixed base character
     * will be assigned to data->bases.bases[i], but not to any
     * of peaks forming the mixed base call
     */

    data->bases.called_peak_list[i] = dpl[dpi_best ];
    data->bases.called_peak_list[i]->data_peak_ind = dpi_best;
    data->bases.called_peak_list[i]->data_peak_ind2= dpi2_best;
    data_peak_ind1[i] = dpi_best;
    data_peak_ind2[i] = dpi2_best;
    if (dpi_best == dpi2_best)
        data_peak_ind2[i] = -1;

    return SUCCESS;
}

/*******************************************************************************
 * Function: call_mixed_bases  
 * Purpose:  Re-evaluate all base calls at their locations, trying to call
 *           mixed bases where necessary.
 *           At each called location, consider a few alternative base calls, 
 *           estimate quality value for each of them and make the call with 
 *           the highest quality value. Return the highest quality value base 
 *           call and all alternatives, together with their locations and 
 *           quality values
 * 
 *           DO NOT reevaluate the base call, if:
 *           - peak resolution > MAX_PRES
 *           - quality value   < MIN_QV
 *******************************************************************************
 */
static int
call_mixed_bases(Data *data, char *color2base, uint8_t *quality_values,
    int *max_num2, AltBase **altbases, int **data_peak_ind1, 
    int **data_peak_ind2, ReadInfo *read_info, BtkLookupTable *table, 
    ContextTable *ctable, Options *options, BtkMessage *message)
{
    char   base, base2,*context = NULL;
    int    i, j, k, n, pos, jbest, kbest, curr_num2, prev_right_bound=0;
    int    found_peak_ind[12], found_peak_is_called[12], num_peaks_found;
    uint8_t qv, orig_qv=0, qv_best=0; 
    int    data_peak_ind_best=-1, 
           data_peak_ind2_best=-1, next_called_peak_found, best_base_used;
    double ratio, ave_wiheight, wiheight, nweight=1., cweight=1.; 
    double ratio_best=0., min_ratio = options->min_ratio;
    int    MIN_QV;
    double MAX_PHR3, MAX_PHR7, MAX_PSR7, MAX_PRES;
    Peak **dpl = data->peak_list;

    curr_num2 = 0;    /* initialize the number of alternative bases */

    get_constraints(&MIN_QV, &MAX_PHR3, &MAX_PHR7, &MAX_PSR7, &MAX_PRES, options);

    if (ctable != NULL)
        context = CALLOC(char, ctable->dimension);

    for (j=0; j<NUM_COLORS; j++) {
        found_peak_ind[j] = 0;
        found_peak_is_called[j] = 0;
    }
 
    for (i=0; i<data->bases.length; i++) {

        if ((ctable != NULL) && (i < ctable->dimension))
            continue;

        base2 = 'N';
        pos = data->bases.called_peak_list[i]->ipos;
        wiheight = data->bases.called_peak_list[i]->wiheight;
        data_peak_ind_best = data_peak_ind2_best = -1;

        /* Normalization and Context weight */
        nweight = 1.;
        cweight = 1.;
        if (options->renorm) {
//          nweight = readNormFactor( read_info, dpl[i]->color_index,
//                      dpl[i]->ipos );
        }
        if ((ctable != NULL) &&
            (i >= ctable->dimension-1) &&
            (i < data->bases.length-1))
        {
//          cweight = weight_from_context(&data->bases.bases[i], ctable);
        }

        if (SHOW_ALL_BASES) {
            fprintf(stderr,
            "i=%d \n  Orig. base=%c ipos=%f data_peak_ind=%d\n", i,
            data->bases.bases[i], data->bases.called_peak_list[i]->ipos, 
            data->bases.called_peak_list[i]->data_peak_ind);
        }

        /* Find all alternative peaks */
        find_alternative_peaks(data, i, found_peak_ind, found_peak_is_called, 
            &num_peaks_found, &next_called_peak_found, &prev_right_bound, 
            MAX_PRES, wiheight * options->min_ratio, options, message);

        /* If no high uncalled peaks found, proceed to the next base call */
        base = data->bases.bases[i];
        
        if (num_peaks_found == 0) 
        {
            data_peak_ind_best = data_peak_ind2_best = 
                data->bases.called_peak_list[i]->data_peak_ind;

            /* Collect info for .poly file */
            (*data_peak_ind1)[i] = found_peak_ind[0];
            (*data_peak_ind2)[i] = -1;
 
            if (SHOW_CALLED_BASES) {
                qv_best = get_pure_base_quality(data, 
                    data->bases.called_peak_list[i]->data_peak_ind, -1,
                    1., 0., table, read_info, options, message);
                fprintf(stderr, "  Final base1=%c ipos=%f qv=%2d data_peak_ind=%d",
                    data->peak_list[data_peak_ind_best]->base, 
                    data->bases.called_peak_list[i]->ipos,
                    qv_best, data_peak_ind_best);
                fprintf(stderr, " data_peak_ind2=%d\n", data_peak_ind2_best); 
            }
            quality_values[i] = get_pure_base_quality(data,
                data->bases.called_peak_list[i]->data_peak_ind, -1,
                1., 0., table, read_info, options, message);
            continue;
        }
        else 
        {
            if (SHOW_ALL_BASES)
                fprintf(stderr, "    %d new peaks found\n", num_peaks_found);

            /* Uncall the current originally called peak;
             * get prepared for evaluating alternative calls
             */
            data->bases.called_peak_list[i]->is_called=0;
            data->bases.called_peak_list[i]->base_index=-1;
            qv_best = 0;
            best_base_used =0;
            jbest = 0;
            kbest = 0;
            data_peak_ind_best  = found_peak_ind[0];
            data_peak_ind2_best = found_peak_ind[0];
            
            /* Find the best regular or mixed base. Use the symmetry of 
             * the function get_quality with respect to its arguments, 
             * which results from the symmetry of the function get_window
             */
            for (n = 0; n <= QVMIN(next_called_peak_found, MERGE_BASE_CALLS); 
                 n++) 
            {
                /* Try to consider the next, originally called peak.
                 * In the n-loop, n == 0 or 1:
                 * n == 0 means that next originally called peak is still 
                 *        considered called as regular base
                 * n == 1 means we will try to couple the next originally called 
                 *        peak it with current originally called peak
                 *        to form a mixed base
                 */
                for (j=0; j<=num_peaks_found; j++) {
                    for (k=j; k<=num_peaks_found; k++) { 

                        if ((n==0) && (j!=0) && (k!=j))
                            continue;

                        if ((n==0) && (k!=0) && found_peak_is_called[k])
                            continue;

                        if ((n==1) && ((j==k) || (j!=0) || 
                            !found_peak_is_called[k]))
                            continue;

                        qv = get_quality_of_alternative_call(i, j, k, 
                             found_peak_ind, found_peak_is_called, 
                             data, table, read_info, options, message);
                        
                        if ((j==0) && (k==0))
                            orig_qv=qv;
 
                        if (SHOW_ALL_BASES) {
                            fprintf(stderr, 
                            "\n    n=%d j=%d k=%d qv(%c)=qv(%c,%c)=%2d qv_best=%2d ",
                                n, j, k, mixed_base(
                                dpl[found_peak_ind[j]]->base,
                                dpl[found_peak_ind[k]]->base),
                                dpl[found_peak_ind[j]]->base,
                                dpl[found_peak_ind[k]]->base,
                                qv, qv_best);
                        }

                        ratio = dpl[found_peak_ind[k]]->iheight/
                                dpl[found_peak_ind[j]]->iheight;
                        
                        if ( qv >  qv_best ||
                            (qv == qv_best && j == k && jbest != kbest)) 
                        {
                            if ((qv_best > 0) &&
                                /* Don't store the orig. base call for now */
                                ((data->bases.called_peak_list[i]->data_peak_ind 
                                                       !=data_peak_ind_best) ||
                                 (data->bases.called_peak_list[i]->data_peak_ind2
                                                       !=data_peak_ind2_best))) 
                            {
                                (void)store_alternative_base_call(data, 
                                    color2base, i, qv_best, data_peak_ind_best, 
                                    data_peak_ind2_best, &curr_num2, max_num2,
                                    altbases, options, message);  
                                if (SHOW_CALLED_BASES) {
                                   fprintf(stderr,
                              "\n    store1 alternative base %c/%d data_peak_ind=%d",
                                   mixed_base(
                                   dpl[data_peak_ind_best ]->base,
                                   dpl[data_peak_ind2_best]->base), i,
                                   data_peak_ind_best);
                                   fprintf(stderr, " data_peak_ind2=%d\n",
                                   data_peak_ind2_best);
                                }
                            } 
                            jbest = j;
                            kbest = k;
                            qv_best    = qv;
                            ratio_best = ratio;
                            data_peak_ind_best  = found_peak_ind[j];
                            data_peak_ind2_best = found_peak_ind[k];
                        }
                        else {
                            /* Store the base call if 
                             * it is different from the orig. base call 
                             */
                            if ((data->bases.called_peak_list[i]->data_peak_ind 
                                   != found_peak_ind[j]) ||
                                (data->bases.called_peak_list[i]->data_peak_ind2
                                   != found_peak_ind[k])) 
                            {
                                (void)store_alternative_base_call(data, 
                                    color2base, i, qv, found_peak_ind[j], 
                                    found_peak_ind[k], &curr_num2, max_num2, 
                                    altbases, options, message);
                            }
                            if (SHOW_CALLED_BASES) {
                                fprintf(stderr,
                            "\n   store2 alternative base %c/%d data_peak_ind=%d",
                                    mixed_base( dpl[found_peak_ind[j]]->base,  
                                    dpl[found_peak_ind[k]]->base), i,
                                    found_peak_ind[j]);
                                fprintf(stderr, " data_peak_ind2=%d\n",
                                    found_peak_ind[k]);
                            }
                        }
                    }           /* loop in k */
                }               /* loop in j */
            }                   /* loop in n */

            ave_wiheight = get_average_called_peak_height(data, i);

            /* Call the best base if it meets constraints */
            if ((qv_best                          >= MIN_QV) && 
                (QVMIN(ratio_best, 1./ratio_best) >= min_ratio)) 
            {
                if ((data->bases.called_peak_list[i]->data_peak_ind 
                         != data_peak_ind_best) ||
                    (data->bases.called_peak_list[i]->data_peak_ind2
                         != data_peak_ind2_best)) 
                {
                   (void)store_alternative_base_call(data, color2base,
                       i, orig_qv, 
                       data->bases.called_peak_list[i]->data_peak_ind, 
                       data->bases.called_peak_list[i]->data_peak_ind2,
                       &curr_num2, max_num2, altbases, options, message);
                   if (SHOW_CALLED_BASES) {
                       fprintf(stderr, 
                       "    store3 alternative base %c/%d data_peak_ind=%d",
                       mixed_base(data->bases.called_peak_list[i]->base, 
                                  data->bases.called_peak_list[i]->base),
                       i, data->bases.called_peak_list[i]->data_peak_ind);
                       fprintf(stderr, " data_peak_ind2=%d\n",
                          data->bases.called_peak_list[i]->data_peak_ind2);
                   }
                }
#if 0
                fprintf(stderr, "  qv_best=%d MIN_QV=%d\n", qv_best, 
                    MIN_QV);
#endif          
                best_base_used= 1;
                quality_values[i] = qv_best;

                set_best_base_call(i, data_peak_ind_best, data_peak_ind2_best, 
                   *data_peak_ind1, *data_peak_ind2, color2base, data,
                    options);
                if (SHOW_CALLED_BASES) {
                    fprintf(stderr, "  Final base2=%c ipos=%f qv=%2d \n",
                        data->bases.called_peak_list[i]->base, 
                        data->bases.called_peak_list[i]->ipos, qv_best);
                }
            }

            /* Otherwise, use the original base call */
            else {
                if (SHOW_MIXED_BASES || SHOW_ALL_BASES) {
                    fprintf(stderr, 
                     "  Base %c did not meet constraints\n",
                    mixed_base(
                    dpl[data_peak_ind_best ]->base,
                    dpl[data_peak_ind2_best]->base));
                    fprintf(stderr, 
                    "    qv_best=%d MIN_QV=%d \n", qv_best, MIN_QV);
                    fprintf(stderr,
                    "    ratio_best=%f min_ratio=%f\n", ratio_best, min_ratio);
                    fprintf(stderr,
                    "    iheigh=%f iheight2=%f ave_wiheight=%f\n",
                    dpl[data_peak_ind_best ]->iheight,
                    dpl[data_peak_ind2_best]->iheight,
                    ave_wiheight);
                }
                data->bases.called_peak_list[i  ]->is_called = 1;
                data->bases.called_peak_list[i  ]->base_index= i;
                next_called_peak_found=0;
                if (i < data->bases.length-1) {
                    data->bases.called_peak_list[i+1]->is_called = 1;
                    data->bases.called_peak_list[i+1]->base_index= i+1;
                }
                data->bases.called_peak_list[i]->base = 
                    data->bases.bases[i];

                quality_values[i] = get_pure_base_quality(data,
                    data->bases.called_peak_list[i]->data_peak_ind, -1,
                    1., 0., table, read_info, options, message);
 
                if ((data->bases.called_peak_list[i]->data_peak_ind 
                         != data_peak_ind_best) ||
                    (data->bases.called_peak_list[i]->data_peak_ind2
                         != data_peak_ind2_best)) 
                {
                   (void)store_alternative_base_call(data, color2base,
                       i, qv_best, data_peak_ind_best, data_peak_ind2_best, 
                       &curr_num2, max_num2, altbases, options, message);
                   if (SHOW_CALLED_BASES) {
                       fprintf(stderr,
                       "    store4 alternative base %c/%d data_peak_ind=%d", 
                       mixed_base(
                       dpl[data_peak_ind_best ]->base,
                       dpl[data_peak_ind2_best]->base), i,
                       data_peak_ind_best);
                       fprintf(stderr, " data_peak_ind2=%d\n",
                       data_peak_ind2_best);
                   }
                }
                if (SHOW_CALLED_BASES) 
                fprintf(stderr, "  Final base3=%c ipos=%f qv=%2d \n",
                    data->bases.called_peak_list[i]->base, 
                    data->bases.called_peak_list[i]->ipos, orig_qv);
            }

            /* Handle the next originally called peak */
            if (next_called_peak_found && best_base_used) {
                if ((jbest!=kbest) &&
                    (found_peak_is_called[jbest]==1) &&
                    (found_peak_is_called[kbest]==1))
                {    
                    uncall_peak(i+1, 
                        data, message);
                }
                else {
                    if (i < data->bases.length-1) {
                        data->bases.called_peak_list[i+1]->is_called = 1;
                        data->bases.called_peak_list[i+1]->base_index = i+1;
                    }
                }
            }
        } /* num_peaks_found > 0 */
    } /* 1st loop (pass) in i */
   *max_num2 = curr_num2;
   *altbases = REALLOC(*altbases, AltBase, curr_num2);
    
#if PERFORM_SECOND_LOOP
    /* 2nd loop
     *
     * Make the 2nd pass through all base calls. Don't recall bases, 
     * but only reevaluate the quality of already called bases and all 
     * their alternatives. The new quality values should, generally, be higher 
     * than the original ones, because there will be fewer uncalled peaks to 
     * the right from the current called peak. If it turns out that, upon 
     * reevaluation, QV of the current base call is lower than QV of any of its 
     * alternative calls, then warning message will be output
     *
     * NOTE: at the second pass, don't merge called peaks
     *
     */
    for (i=0; i<data->bases.length; i++) {

        /* Store information about current base call */
        int orig_data_peak_ind, orig_data_peak_ind2;
        orig_data_peak_ind = data->bases.called_peak_list[i]->data_peak_ind;
        orig_data_peak_ind2= data->bases.called_peak_list[i]->data_peak_ind2;
     
        /* Re-evaluate quality value of current base call */
        found_peak_ind[0]       = orig_data_peak_ind;
        found_peak_ind[1]       = orig_data_peak_ind2;
        found_peak_is_called[1] = 0;
        quality_values[i] = get_quality_of_alternative_call(i, 0, 1,
            found_peak_ind, found_peak_is_called, data, table, read_info, 
            options, message);
        qv_best = quality_values[i];

        /* Unset current base call temporarily */
        data->peak_list[orig_data_peak_ind ]->is_called = 0;
        data->peak_list[orig_data_peak_ind2]->is_called = 0;
        data->bases.called_peak_list[i]->base_index = -1;
                 
        /* Reevaluate quality value of alternative base calls */
        for (j=0; j<curr_num2; j++) {

            if ((*altbases)[j].base_index == i) {

                /* Set alternative call */
                data->peak_list[(*altbases)[j].data_peak_ind ]->is_called = 1;
                data->peak_list[(*altbases)[j].data_peak_ind2]->is_called = 1;
                data->bases.called_peak_list[i] = 
                    data->peak_list[(*altbases)[j].data_peak_ind ];
                data->bases.called_peak_list[i]->base_index = i;

                /* Reevaluate its quality */
                found_peak_ind[0] = (*altbases)[j].data_peak_ind ;
                found_peak_ind[1] = (*altbases)[j].data_peak_ind2;
                (*altbases)[j].qv = get_quality_of_alternative_call(i, 0, 1,
                    found_peak_ind, found_peak_is_called, data, table, read_info,
                    options, message);
                if ((*altbases)[j].qv > qv_best) 
                {
                    int temp;
                    char temp_base;

                    if (options->Verbose > 1) {
                        fprintf(stderr, 
                        "Warning: new mixed base call detected in 2nd pass:\n");
                        fprintf(stderr,
                        "Base index=%d qv called=%d qv alternative=%d\n",
                        i, qv_best, (*altbases)[j].qv);
                    }

                    /* Make the alternative base call current 
                     * and original current base call alternative
                     */
                    temp = (*altbases)[j].data_peak_ind;
                    (*altbases)[j].data_peak_ind = orig_data_peak_ind;
                    orig_data_peak_ind = temp;

                    temp = (*altbases)[j].data_peak_ind2;
                    (*altbases)[j].data_peak_ind2 = orig_data_peak_ind2;
                    orig_data_peak_ind2 = temp;

                    (*data_peak_ind1)[i] = orig_data_peak_ind;
                    (*data_peak_ind2)[i] = orig_data_peak_ind2;

                    temp = (*altbases)[j].qv;
                    (*altbases)[j].qv = quality_values[i];
                    quality_values[i] = temp;

                    temp_base = (*altbases)[j].base;
                    (*altbases)[j].base = data->bases.bases[i];
                    data->bases.bases[i] = temp_base;
                } 

                /* Unset alternative call */
                data->peak_list[(*altbases)[j].data_peak_ind ]->is_called = 0;
                data->peak_list[(*altbases)[j].data_peak_ind2]->is_called = 0;
                data->bases.called_peak_list[i]->base_index = -1;
            }   
        } 

        /* Set back the original base call */
        data->peak_list[orig_data_peak_ind ]->is_called = 1;
        data->peak_list[orig_data_peak_ind2]->is_called = 1;
        data->bases.called_peak_list[i] =
            data->peak_list[orig_data_peak_ind ];
        data->bases.called_peak_list[i]->base_index = i;
    } /* 2nd loop (pass) in i */
#endif

    return SUCCESS;
}

/*******************************************************************************
 * Function: total_number_of_peaks
 *******************************************************************************
 */
static int
total_number_of_peaks(Data *data)
{
    return data->color_data[0].peak_list_len +
           data->color_data[1].peak_list_len +
           data->color_data[2].peak_list_len +
           data->color_data[3].peak_list_len;
}

/*******************************************************************************
 * Function: Btk_get_mixed_bases
 *******************************************************************************
 */
int
Btk_get_mixed_bases(int *num_bases, char **bases, int **peak_locs, 
    int num_datapoints, int **chromatogram, char *color2base, 
    uint8_t **quality_values, ReadInfo *read_info, BtkLookupTable *table, 
    ContextTable *ctable, Options options, BtkMessage *message, 
    Results *results)
{
    int      i, num2, renorm;
    int     *data_peak_ind1=NULL, *data_peak_ind2=NULL; /* for .poly file */
    AltBase *altbases;
    Data     data; 
    clock_t start_clock = clock(), curr_clock;

    if (SHOW_INPUT_OPTIONS)
        show_input_options(&options);

    if (data_create(&data, num_datapoints, *num_bases, color2base, message)
	!= SUCCESS)
    {
        sprintf(message->text, "Error calling data_create\n");
        return ERROR;
    }
 
    if (data_populate(num_bases, bases, options.edited_bases,
        peak_locs, num_datapoints, chromatogram, color2base, &data, 
        &options, message) != SUCCESS)
    {
        sprintf(message->text, "Error calling data_populate\n");
        fprintf(stderr, "Error calling  data_populate\n");
        goto error;
    }

    if (options.time) {  
        curr_clock = clock();
        fprintf(stderr, "Data structure populated in %f sec. \n",
            (float)(curr_clock - start_clock)/(float)CLOCKS_PER_SEC);
        start_clock = curr_clock;
    }

    if (Btk_process_peaks(&data, &options, message) != SUCCESS) {
        if (total_number_of_peaks(&data) >= 10) {
            sprintf(message->text, "Error calling Btk_process_peaks\n");
            fprintf(stderr, "Error calling Btk_process_peaks\n");
        }
        goto error;
    }
 
    if (options.time) {
        curr_clock = clock();
        fprintf(stderr, "Peaks processed in %f sec. \n",
            (float)(curr_clock - start_clock)/(float)CLOCKS_PER_SEC);
        start_clock = curr_clock; 
    }  

#if INSERT_BASES
    if (Btk_call_bases(&data, color2base, read_info, ctable,
        &options, message, results ) != SUCCESS) {
        fprintf(stderr, "Error calling bases\n");
        goto error;
    }
#else
    if (mark_called_peaks(&data, color2base, options, message) != SUCCESS)
    {
        fprintf(stderr, "Error calling pure_bases\n");
        goto error;
    }
#endif

    /* Output .poly file */
    if (options.poly) {
        if (Btk_output_poly_file(&data, &options, message) != SUCCESS)
        {
            fprintf(stderr, "Error creating  poly file\n");
            goto error;
        }
    }

    // Multiply peak iheights by context weight
    {
        int base_index = -1;
        char context[4] = {' ', ' ', ' ', '\0'};
        double cweight;
        for (i=0; i<data.peak_list_len; i++)
        {
            if (data.peak_list[i]->base_index >= 2) {
                base_index = data.peak_list[i]->base_index;
                context[0] = data.bases.bases[base_index - 2];
                context[1] = data.bases.bases[base_index - 1];
            }
            context[2] = data.peak_list[i]->base;

            if (base_index >= 2) {
                cweight = get_context_weight(context);
                data.peak_list[i]->wiheight = 
                    data.peak_list[i]->iheight * cweight;
            }
        }
    }

    if (options.time) {  
            curr_clock = clock();
            fprintf(stderr, "Bases called in %f sec. \n",
                (float)(curr_clock - start_clock)/(float)CLOCKS_PER_SEC);
            start_clock = curr_clock;
    }                          
        
    num2           = data.bases.length * 6;
    altbases       = CALLOC(AltBase, num2);
    data_peak_ind1 = CALLOC(int,  num2);
    data_peak_ind2 = CALLOC(int,  num2);
   *quality_values = REALLOC(*quality_values, uint8_t, data.bases.length);
    MEM_ERROR(*quality_values);

    if (options.het || options.mix)
    {
        if (options.Verbose > 1) {
            fprintf(stderr, "Calling mixed bases; min_ratio==%.2f\n", 
                options.min_ratio);
        }
      
        /* Don't renormalize in call_mixed_bases */ 
        renorm = options.renorm;
        options.renorm = 0; 
        if (call_mixed_bases(&data, color2base, *quality_values, 
            &num2, &altbases, &data_peak_ind1, &data_peak_ind2,
            read_info, table, ctable, &options, message) != SUCCESS) 
        {
            fprintf(stderr, "Error calling mixed bases\n");
            goto error;
        }
         
        options.renorm = renorm;
    }
    else {
        num2 = 0;
    }

    if (*num_bases != data.bases.length) {
        *num_bases  = data.bases.length;

       *bases = REALLOC(*bases, char, data.bases.length);
        MEM_ERROR(*bases);
       *peak_locs = REALLOC(*peak_locs, int, data.bases.length);
        MEM_ERROR(*peak_locs);
    }

    for (i = 0; i < data.bases.length; i++) {
        (*bases)[i] = data.bases.bases[i];
        if (data.bases.called_peak_list[i]->data_peak_ind ==
            data.bases.called_peak_list[i]->data_peak_ind2)
            (*peak_locs)[i] = data.bases.called_peak_list[i]->ipos;
        else
            (*peak_locs)[i] = get_mixed_base_position(&data, 
            data.bases.called_peak_list[i]->data_peak_ind,
            data.bases.called_peak_list[i]->data_peak_ind2, 
            &options, message);
    }


    if (options.time) {   
        curr_clock = clock();
        fprintf(stderr, "Trace parameters computed in %f sec. \n",
            (float)(curr_clock - start_clock)/(float)CLOCKS_PER_SEC);
        start_clock = curr_clock;
    }     

    if (SHOW_MIXED_BASES || SHOW_ALL_BASES) 
    {
        double ratio, wratio, res, w1, w2;
        int data_peak_ind, data_peak_ind2, jc;
        char base, b1, b2;
        for (i=0; i<data.bases.length; i++) {
            if (!is_mixed_base(data.bases.called_peak_list[i]->base))
                continue;

            data_peak_ind = data.bases.called_peak_list[i]->data_peak_ind;
            data_peak_ind2= data.bases.called_peak_list[i]->data_peak_ind2;

            ratio = QVMIN(
                  data.peak_list[data_peak_ind]->iheight/
                  data.peak_list[data_peak_ind2]->iheight,
                  data.peak_list[data_peak_ind2]->iheight/
                  data.peak_list[data_peak_ind]->iheight);
            res = QVMAX(
                  data.peak_list[data_peak_ind]->resolution,
                  data.peak_list[data_peak_ind2]->resolution);
            base = data.bases.bases[i];

            jc = data.peak_list[data_peak_ind]->color_index;
            data.bases.bases[i] = color2base[jc];
            w1 = 1.;
            if ((ctable != NULL) && (i < ctable->dimension)) {
                w1 = weight_from_context(&data.bases.bases[i], ctable);
            }
            jc = data.peak_list[data_peak_ind2]->color_index;
            data.bases.bases[i] = color2base[jc];
            w2 = 1.;
            if ((ctable != NULL) && (i < ctable->dimension)) {
                w2 = weight_from_context(&data.bases.bases[i], ctable);
            }
            data.bases.bases[i] = base;

            wratio = QVMIN(
                  data.peak_list[data_peak_ind ]->iheight * w1/
                  data.peak_list[data_peak_ind2]->iheight / w2,
                  data.peak_list[data_peak_ind2]->iheight * w2/
                  data.peak_list[data_peak_ind ]->iheight / w1);
            b1 = data.peak_list[data_peak_ind ]->base;
            b2 = data.peak_list[data_peak_ind2]->base;

            fprintf(stderr,
                "Detected mixed base[ %d ]= %c =(%c,%c) phr= %f ",
                i+1, mixed_base(b1, b2), b1, b2, ratio);
        }
    }

    if (options.tip_dir[0] != '\0') {
        if (Btk_output_tip_file(&data, color2base, options)
            != SUCCESS)
        {
            fprintf(stderr, "Error calling Btk_output_tip_file\n");
            goto error;
        }
    }

    if (options.tab_dir[0] != '\0') {
        /* Output alternative base calls */
        if (Btk_output_tab_file(num2, altbases, options) != SUCCESS)
        {
            fprintf(stderr, "Error writing .tab file\n");
            goto error;
        }
    }
    data_release(&data);
    FREE(altbases);
    FREE(data_peak_ind1);
    FREE(data_peak_ind2);

    return SUCCESS;

error:
    data_release(&data);
    FREE(altbases);
    FREE(data_peak_ind1);
    FREE(data_peak_ind2);

    return ERROR;
}
