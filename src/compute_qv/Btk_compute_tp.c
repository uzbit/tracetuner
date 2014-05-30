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
 * 1.73 2003/11/06 18:18:26
 */

/*
 *  Btk_compute_tp.c  $Revision: 1.18 $
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdint.h>

#include "Btk_qv.h"
#include "SFF_Toolkit.h"
#include "util.h"
#include "Btk_qv_data.h"
#include "Btk_lookup_table.h"
#include "context_table.h"
#include "train.h"
#include "Btk_compute_qv.h"
#include "Btk_qv_funs.h"
#include "Btk_process_peaks.h"  /* needs train.h */
#include "tracepoly.h"
#include "Btk_compute_tpars.h"
#include "Btk_compute_tp.h"
#include "Btk_call_bases.h"
#include "Btk_get_mixed_bases.h"

#define CALLED_HEIGHT_THRESHOLD 50.
#define CHECK_BASE_INDEX         0
#define CHECK_PEAK_INDEX         0
#define CHECK_WINDOW             0
#define DEBUG                    0
#define MAX_SPACING_FACTOR       1.2
#define MIN_SPACING_FACTOR       0.8
#define OUTPUT_TRACE_PARAMS      0
#define SHOW_SPACING             0
#define SHOW_SCALING_FACTORS     0
#define SHOW_WINDOW              0
#define TRUNCATED_HEIGHT      1600
#define USE_MAX_RESOLUTION       1
#define USE_OLD_PHR              0
#define USE_OLD_PSR              0
#define USE_OLD_PRES             0


/****************************************************************************
 * Function: get_window
 * Purpose: determine indexes of peaks corresponding to the boundaries of
 *          a window of given size, centered at the current peak.
 *          The window consists of a given number of called peaks.  
 *          When determining the window, peaks with index data_peak_ind2
 *          should be ignored, even if it is called
 ****************************************************************************
 */
static void   
get_window(Data *data, int data_peak_ind,  int data_peak_ind2, int window, 
    int *il, int *ir, int *lwin, int *rwin, Options *options)
{
    int  i=data_peak_ind, half_window=window/2;

    if (CHECK_WINDOW && !data->peak_list[data_peak_ind]->is_called) {
        fprintf(stderr, 
        "Passed uncalled peak %d to get_window function\n",
        data_peak_ind);
        exit_message(options, -1);
    }

    /* Find the index of the leftmost peak in the window */
   *lwin = 0;                        /* counter of called peaks */
   *il = (i>0) ? (i-1) : i;
    while ((*il >= 0) && (*lwin < half_window)) {

        if ( data->peak_list[*il]->is_called &&
            *il != data_peak_ind2 &&

            /* If a window contains mixed peaks, count them only once */
            (data->peak_list[*il]->base_index >= 0)) 
        {
            (*lwin)++;
        }
        (*il)--;
    }
    /* Correct il after the while loop */
    (*il)++;

    /* Find the index of the rightmost peak in the window */
   *rwin = 0;                        /* counter of called peaks */
   *ir = (i < data->peak_list_len-1) ? (i+1) : i;
    while ((*ir <= data->peak_list_len-1) && (*rwin < half_window)) {

        if ( data->peak_list[*ir]->is_called &&
            *ir != data_peak_ind2 &&

            /* If a window contains mixed peaks, count them only once */
            (data->peak_list[*ir]->base_index >= 0)) 
        {
            (*rwin)++;
        }
        (*ir)++;
    }
    /* Correct ir after the while loop */
    (*ir)--;
}

/****************************************************************************
 * Function: get_max_uncalled_height
 * Purpose:  Determine the heihght of the largest uncalled peak in the
 *           specified window. 
 * Arguments: 
 ****************************************************************************
 */
double
get_max_uncalled_height(Data *data, int data_peak_ind2, double frac, 
    int il, int ir, double *max_uncalled_ipos, Options *options)
{
   int    j;
   double max_uncalled_iheight = 0., iheight, ipos,
          min_called_peak_height=MIN_CALLED_PEAK_HEIGHT;

   for (j = il; j <= ir; j++) {

        if (j != data->peak_list[j]->data_peak_ind) {
            fprintf(stderr, "Error in data_peak_ind\n");
            return ERROR;
        }

        if (data->peak_list[j]->is_called || 
            j == data_peak_ind2 ||
            data->peak_list[j]->base_index >= 0 ||
            data->peak_list[j]->data_peak_ind != j) 
            continue;

        iheight = data->peak_list[j]->iheight;
        ipos    = data->peak_list[j]->ipos;
        if (!options->recalln && !options->recallndb && !options->ladder) {
            if (is_dye_blob(data->peak_list[j]->ipos,
                data->peak_list[j], data, 
                QVMAX(options->het, options->mix)))
            iheight *= DYE_BLOB_FRACTION;
            if (iheight < min_called_peak_height)
                iheight = min_called_peak_height;
        }

        if (j != data_peak_ind2 && max_uncalled_iheight < iheight*frac)
        {
            max_uncalled_iheight = iheight*frac;
           *max_uncalled_ipos    = ipos;
        }

    }  /* loop through the window */

    return max_uncalled_iheight;
}

/****************************************************************************
 * Function: get_min_max_called_height
 * Purpose:  Determine the height of the shortest and highest called peak 
 *           in the window of specified size. If mixed base is called, use 
 *           a sum of the two peak heights
 *           data_peak_ind - index of the central called peak,
 *                           which height should not be scaled
 *           data_peak_ind2- index of peak which should be ignored 
 *                           (taken off from chromatogram)
 *           frac          - scaling factor
 ****************************************************************************
 */
static int   
get_min_max_called_height(int data_peak_ind, int data_peak_ind2,
    double *min_called, double *max_called, double *min_called_ipos,
    double frac, Data *data, int il, int ir, Options *options)
{
   int    j;
   double iheight; 
   Peak **pl = data->peak_list;

   *min_called = pl[data_peak_ind]->iheight;
   *max_called = *min_called;
   *min_called_ipos = pl[data_peak_ind]->ipos;
   
   for (j = il; j <= ir; j++) {
        int j2 = data->peak_list[j]->data_peak_ind2;

        if (j != data->peak_list[j]->data_peak_ind) {
            fprintf(stderr, "Error in data_peak_ind\n");
            return ERROR;
        }

        if ((data->peak_list[j]->is_called == 0) ||
            (j == data_peak_ind2) ||
            (data->peak_list[j]->base_index < 0)) 
            continue;

        /* Get the height of the current called peak */      
        iheight = data->peak_list[j]->iheight;
         
        /* Scale the peak height */
        if (j != data_peak_ind)
            iheight *= frac;
        else
            iheight = QVMAX(iheight, 1);

        if (!options->recalln && !options->recallndb && !options->ladder) {
            if (is_dye_blob(data->peak_list[j]->ipos,
                    data->peak_list[j ], data, 
                    QVMAX(options->het, options->mix)) ||
                is_dye_blob(data->peak_list[j2]->ipos,
                    data->peak_list[j2], data,
                    QVMAX(options->het, options->mix)))
            {
                iheight *= DYE_BLOB_FRACTION;
                if (iheight < MIN_CALLED_PEAK_HEIGHT)
                    iheight = MIN_CALLED_PEAK_HEIGHT;
            }
        }

        /* Update min_called and max_called iheight */
        if (iheight < *min_called) {
           *min_called = iheight;
           *min_called_ipos = data->peak_list[j]->ipos;
        }
        if (iheight > *max_called) {
           *max_called = iheight;
        }

    }  /* end loop through the window */

    return SUCCESS; 
}

static int
get_min_max_spacing(int data_peak_ind, int data_peak_ind2, double *min_spacing, 
    double *max_spacing, double *min_spacing_ipos, double *max_spacing_ipos, 
    Data *data, int il, int ir, Options *options, BtkMessage *message)
{
    int    j, prev_called_ind;
    double called_pos, prev_called_pos;
    double spacing;
    Peak **pl  = data->peak_list;
    Peak **cpl = data->bases.called_peak_list;
    int base_ind = pl[data_peak_ind]->base_index;

    if (DEBUG && base_ind < 0) {
        fprintf(stderr, "Error in get_min_max_spacing: central peak %d is_called= %d base_ind<0\n",
        data_peak_ind, pl[data_peak_ind]->is_called);
        return -1;
    }

   *min_spacing = base_ind > 0 ? cpl[base_ind]->ipos - cpl[base_ind-1]->ipos : 
                                 cpl[base_ind+1]->ipos - cpl[base_ind]->ipos;
   *max_spacing = *min_spacing;
   *min_spacing_ipos = *max_spacing_ipos = pl[data_peak_ind]->ipos;
    prev_called_ind = il;
    prev_called_pos = data->peak_list[il]->ipos;

    /* Correction for the case where leftmost called peak is mixed */
    if ((options->het || options->mix) &&
        data->peak_list[il]->data_peak_ind!=data->peak_list[il]->data_peak_ind2)
        prev_called_pos = get_mixed_base_position(data, il,
            data->peak_list[il]->data_peak_ind2, options, message);

    for (j = il+1; j <= ir; j++)
    {
        if ((data->peak_list[j]->is_called == 0) ||
            (j == data_peak_ind2) ||
            (data->peak_list[j]->base_index < 0))
            continue;
         
        if (data->peak_list[j]->data_peak_ind2 == prev_called_ind )
            continue;
            

        called_pos = data->peak_list[j]->ipos;
        spacing = called_pos - prev_called_pos;

        /* Make sure the spacing between the two called peaks
         * corresponding to the mixed base is not counted
         */
        if (*min_spacing > spacing) {
            *min_spacing = spacing;
            *min_spacing_ipos = j<=data_peak_ind ? pl[j]->ipos : pl[j-1]->ipos;
        }
        if (*max_spacing < spacing) {
            *max_spacing = spacing;
            *max_spacing_ipos = j<=data_peak_ind ? pl[j]->ipos : pl[j-1]->ipos;
        }

        prev_called_pos = called_pos;
        prev_called_ind = j;
    }
    return SUCCESS;
}

/****************************************************************************
 * Function: get_peak_height_ratio in the window of specified number
 *           of called peaks. Both the peaks passed to this function
 *           should be called.
 *           data_peak_ind - index of the current called peak;
 *                           its height should not be scaled
 *           data_peak_ind2- index of the peak which should be ignored
 *                           (that is, taken off from chromatogram)
 *           frac          - multiplier of peak height for non-current peaks
 ****************************************************************************
 */
double 
get_peak_height_ratio(Data *data, int data_peak_ind, int data_peak_ind2, 
    int window, double frac, double min_uncalled, Options *options, 
    BtkMessage *message)
{
    int il, ir; 
    int lwin, rwin;  // # of called peaks in a window to the left and 
                     // right from the central peak
    double max_uncalled_height=min_uncalled, min_called_height, 
           max_called_height, phr;
    Peak **pl = data->peak_list;
    double ipos = pl[data_peak_ind]->ipos;
    double max_uncalled_ipos = -1, min_called_ipos = -1, bad_ipos;

    get_window(data, data_peak_ind, data_peak_ind2, window, &il, &ir, &lwin, 
        &rwin, options);

    if (CHECK_WINDOW) {
        fprintf(stderr, 
        "il=%d il_called=%d ilpos=%f ir=%d ir_called=%d irpos=%f window=%d\n",
        il, pl[il]->is_called, pl[il]->ipos,
        ir, pl[ir]->is_called, pl[ir]->ipos,
        abs(pl[il]->base_index - pl[ir]->base_index)+1);
    } 

    max_uncalled_height = QVMAX(min_uncalled, get_max_uncalled_height(data, 
        data_peak_ind2, frac, il, ir, &max_uncalled_ipos, options));

    if (get_min_max_called_height(data_peak_ind, data_peak_ind2,
        &min_called_height, &max_called_height, &min_called_ipos,
        frac, data, il, ir, options) != SUCCESS)
        return ERROR;

    bad_ipos = abs(ipos - min_called_ipos) < abs(ipos - max_uncalled_ipos) ?
               min_called_ipos : max_uncalled_ipos;

    phr = (min_called_height > 0) ? 
           max_uncalled_height / min_called_height : 100.0;
  
    if (SHOW_SPACING) {
        int base_ind = pl[data_peak_ind]->base_index;
        if (window == 3)
        fprintf(stderr,
        "    i=%d min_called=%f max_uncalled=%f phr3=%5.3f\n",
             base_ind, min_called_height, max_uncalled_height, phr);
        if (window == 7)
        fprintf(stderr,
        "    i=%d min_called=%f max_uncalled=%f phr7=%5.3f\n",
             base_ind, min_called_height, max_uncalled_height, phr);
    }
    return phr;      
}

/*******************************************************************************
 * Function: get_peak_spacing_ratio
 * Purpose: Calculate the ratio of highest to lowest spacing between called
 *          peaks in the window of 7 called peaks centered at the current one
 *******************************************************************************
 */
double
get_peak_spacing_ratio(Data *data, int data_peak_ind, int data_peak_ind2,
    int window, Options *options, BtkMessage *message)
{
    int    j, il, ir, lwin, rwin;
    double min_spacing, max_spacing;
    double min_spa_ipos, max_spa_ipos;
    double psr, bad_ipos, ipos;
    Peak **pl = data->peak_list;

    /* Get indexes of the boundary peaks of the window.
     * If the boundary peak corresponds to a mixed base, then il and ir
     * are indexes of the peak that has base_index > 0
     */
    get_window(data, data_peak_ind, data_peak_ind2, window, &il, &ir, &lwin, 
        &rwin, options);

    if (SHOW_WINDOW)
        fprintf(stderr, 
        "        window: il=%d ir=%d\n", il, ir);
       
    if (DEBUG) {
        for (j = il+1; j <= ir; j++) {
            if ((j > 0) && 
                (data->peak_list[j]->ipos < data->peak_list[j-1]->ipos))
                fprintf(stderr, 
                "Disordered peaks: peak[%d].ipos=%f < peak[%d].ipos=%f data_peak_ind 1,2= %d %d\n",
                j, data->peak_list[j]->ipos, j-1, data->peak_list[j-1]->ipos,
                data_peak_ind, data_peak_ind2);
        }
    } 

    if (get_min_max_spacing(data_peak_ind, data_peak_ind2, &min_spacing, 
        &max_spacing, &min_spa_ipos, &max_spa_ipos, data, il, ir, options, 
        message) != SUCCESS)
        return ERROR;

    ipos = pl[data_peak_ind]->ipos;
    bad_ipos = abs(ipos - max_spa_ipos) < abs(ipos - min_spa_ipos) ?
               max_spa_ipos: min_spa_ipos; 

    /* Calculate the peak spacing ratio */
    psr = max_spacing / (min_spacing>0?min_spacing:1.);

    if (SHOW_SPACING) {
        int base_ind = pl[data_peak_ind]->base_index;
        fprintf(stderr, 
        "    i=%d min_spacing=%f max_spacing=%f psr=%5.3f\n",
             base_ind, min_spacing, max_spacing, psr);
    } 
    return psr;
}

/*******************************************************************************
 * Function: get_peak_resolution
 * Purpose:  for regular base, just return its peak resolution;
 *           for mixed base, return the bigger of the two peak resolutions
 *******************************************************************************
 */
double
get_peak_resolution(Data *data, int data_peak_ind, Options *options,
    BtkMessage *message)
{
    double res = data->peak_list[data_peak_ind]->resolution;

    return (data->peak_list[data_peak_ind]->height > CALLED_HEIGHT_THRESHOLD) ?
           res : 1.;
}

/*******************************************************************************
 * Function: get_peak_resolution_parameter
 * Purpose:  for regular base, just return its peak resolution; 
 *           for mixed base, return the bigger of the two peak resolutions
 *******************************************************************************
 */
static double
get_peak_resolution_parameter(Data *data, int data_peak_ind, Options *options, 
    BtkMessage *message)
{
    int i;
    int base_index = data->peak_list[data_peak_ind]->base_index;
    int win_beg = base_index - WINDOW_7/2 ;
    int win_end = base_index + WINDOW_7/2 ;
    double res = data->peak_list[data_peak_ind]->resolution;
    double max_res = NINF;

    if (base_index < 0)
        return 1.;

    res = NINF;
    if (win_beg < 0)
        win_beg = 0;
    if (win_end > data->bases.length-1)
        win_end = data->bases.length-1;
      
    for (i=win_beg; i <= win_end; i++)
    {
        res = (data->bases.called_peak_list[i]->height >
               CALLED_HEIGHT_THRESHOLD) ?  
               data->bases.called_peak_list[i]->resolution : 1.;

        if (max_res < res)
            max_res = res;
    }

    return res;    
}


/*******************************************************************************
 * Function: get_trace_parameters_of_pure_bases             
 * Purpose: for each called base, compute 4 trace parameters
 * Note: this function will be called only if none of -het and -mix
 *       options is specified 
 *******************************************************************************
 */
int
get_trace_parameters_of_pure_bases(int base_index, Data *data, char *color2base, 
    Options *opts, BtkMessage *message)
{
    int i;

#if SHOW_REVIEW
    data_review_peak_list(data, message);
#endif
    if (base_index < 0 || base_index > data->bases.length-1)
    {
        fprintf(stderr, "Base index out of range\n");
        return ERROR;
    }
    
    i = data->bases.called_peak_list[base_index]->data_peak_ind; 
           
    data->trace_parameters.phr3[base_index] =     
        get_peak_height_ratio(data, i, -1, WINDOW_3, 1., 0., opts, message);
    data->trace_parameters.phr7[base_index] =
        get_peak_height_ratio(data, i, -1, WINDOW_7, 1., 0., opts, message);
    data->trace_parameters.psr7[base_index] =
        get_peak_spacing_ratio(data, i, -1, WINDOW_7, opts, message);
    data->trace_parameters.pres[base_index] =
        get_peak_resolution_parameter(data, i, opts, message);
 
    return SUCCESS;
}

double
get_context_weight(char *context)
{
    if (strcmp(context, "GCC") == 0) return 0.959166283895999;
    if (strcmp(context, "AGT") == 0) return 1.42537877662666;
    if (strcmp(context, "TGA") == 0) return 0.917418796329302;
    if (strcmp(context, "TGT") == 0) return 1.25654181086771;
    if (strcmp(context, "CGA") == 0) return 0.841592444429528;
    if (strcmp(context, "ATC") == 0) return 0.993893558355747;
    if (strcmp(context, "AAC") == 0) return 1.05654744889751;
    if (strcmp(context, "AGC") == 0) return 1.18098395806666;
    if (strcmp(context, "TAC") == 0) return 1.08362343932169;
    if (strcmp(context, "ACA") == 0) return 1.19862003804775;
    if (strcmp(context, "TCG") == 0) return 1.5240544867462;
    if (strcmp(context, "CCG") == 0) return 1.06012104834425;
    if (strcmp(context, "CTG") == 0) return 0.908899999494062;
    if (strcmp(context, "GCA") == 0) return 1.17590500872575;
    if (strcmp(context, "GTG") == 0) return 1.00091503789467;
    if (strcmp(context, "AAG") == 0) return 1.18096908535966;
    if (strcmp(context, "CAC") == 0) return 0.992981200364162;
    if (strcmp(context, "GTT") == 0) return 1.06576084606065;
    if (strcmp(context, "AGA") == 0) return 0.826648235645059;
    if (strcmp(context, "ACC") == 0) return 0.944609342994373;
    if (strcmp(context, "CCA") == 0) return 1.26458658869927;
    if (strcmp(context, "TGG") == 0) return 1.12601644195801;
    if (strcmp(context, "CGC") == 0) return 0.984679305072416;
    if (strcmp(context, "CTC") == 0) return 0.978809857198512;
    if (strcmp(context, "TTG") == 0) return 1.07508351963171;
    if (strcmp(context, "TAA") == 0) return 1.12094141287651;
    if (strcmp(context, "CAG") == 0) return 1.27053380942843;
    if (strcmp(context, "ACG") == 0) return 1.23756760079144;
    if (strcmp(context, "ATG") == 0) return 0.978948205319427;
    if (strcmp(context, "AAA") == 0) return 0.993778161758255;
    if (strcmp(context, "GTA") == 0) return 1.20762324447585;
    if (strcmp(context, "TAG") == 0) return 1.4880979348256;
    if (strcmp(context, "CTT") == 0) return 0.985210564542276;
    if (strcmp(context, "GGA") == 0) return 0.968137129189855;
    if (strcmp(context, "GTC") == 0) return 1.03912433651762;
    if (strcmp(context, "TGC") == 0) return 1.11010715213778;
    if (strcmp(context, "TCA") == 0) return 1.27934352035566;
    if (strcmp(context, "ATT") == 0) return 1.04294054731527;
    if (strcmp(context, "TAT") == 0) return 0.854983207664805;
    if (strcmp(context, "AAT") == 0) return 1.08490902268228;
    if (strcmp(context, "ACT") == 0) return 1.02438215597705;
    if (strcmp(context, "CAA") == 0) return 1.08755893009568;
    if (strcmp(context, "GAC") == 0) return 1.0055742087977;
    if (strcmp(context, "GGT") == 0) return 1.54160745527815;
    if (strcmp(context, "TCC") == 0) return 1.00679673688028;
    if (strcmp(context, "TTT") == 0) return 1.07154651507132;
    if (strcmp(context, "AGG") == 0) return 1.01235739074533;
    if (strcmp(context, "CGT") == 0) return 1.36323936614225;
    if (strcmp(context, "ATA") == 0) return 1.32048288218383;
    if (strcmp(context, "CGG") == 0) return 0.98245359837622;
    if (strcmp(context, "CAT") == 0) return 0.996356010149638;
    if (strcmp(context, "CCC") == 0) return 1.00003821006626;
    if (strcmp(context, "GGG") == 0) return 1.16862795443708;
    if (strcmp(context, "TTA") == 0) return 1.18016843163403;
    if (strcmp(context, "GAG") == 0) return 1.17601780194004;
    if (strcmp(context, "CTA") == 0) return 1.1069205981111;
    if (strcmp(context, "GAT") == 0) return 1.01828760021087;
    if (strcmp(context, "TCT") == 0) return 1.08595471272615;
    if (strcmp(context, "TTC") == 0) return 1.02873478854932;
    if (strcmp(context, "GCG") == 0) return 1.07917396588071;
    if (strcmp(context, "GGC") == 0) return 1.33801056267355;
    if (strcmp(context, "GCT") == 0) return 1.01651951359486;
    if (strcmp(context, "GAA") == 0) return 0.926244912611492;
    if (strcmp(context, "CCT") == 0) return 1.15021622512442;
    return 1.;
}

/*******************************************************************************
 * Function: populate_params_array   
 * Purpose: populate the array of trace parameters which the function
 *          Btk_compute_tp will output
 *******************************************************************************
 */
int
populate_params_array(int i, Data* data, double** params)
{
    if ((data->trace_parameters.phr3[i] < 0) || 
        (data->trace_parameters.phr7[i] < 0)) {
        printf("i=%d phr3, phr7= %f %f\n",i, 
        data->trace_parameters.phr3[i],
        data->trace_parameters.phr7[i]);
    }
    params[0][i] = data->trace_parameters.phr3[i];
    params[1][i] = data->trace_parameters.phr7[i];
    params[2][i] = data->trace_parameters.psr7[i];
    params[3][i] = data->trace_parameters.pres[i];
    
    return SUCCESS;
}

/*******************************************************************************
 * Function: Btk_compute_tp
 * Purpose: Calculate trace parameters 
 * Note: this function will be called only if none of -het and -mix
 *       options is specified
 * 
 *******************************************************************************
 */
int
Btk_compute_tp(Data *data, char *color2base, int NUM_PARAMS, 
    double *params[NUM_PARAMS], ReadInfo *read_info, ContextTable *ctable, 
    Options *options, BtkMessage *message)
{
    int i;

    if (CHECK_PEAK_INDEX) {
        int i, j;
        for (j=0; i<NUM_COLORS; j++) {
            for (i=0; i<data->color_data[j].peak_list_len; i++) {
                if (i != data->color_data[j].peak_list[i].cd_peak_ind)
                    fprintf(stderr,
                    "peak_ind=%d != cd_peak_ind=%d for color j=%d\n",
                    i, data->color_data[j].peak_list[i].cd_peak_ind, j);
            }
        }
    }

    for (i=0; i<data->bases.length; i++)
    {    
        if (get_trace_parameters_of_pure_bases(i, data, color2base, options,
            message) != SUCCESS) {
            sprintf(message->text, "Error calling get_trace_parameters_of_pure_bases\n");
	    goto error;
        }

        if (populate_params_array(i, data, params) != SUCCESS) {
            sprintf(message->text, "Error calling populate_params_array\n");
            goto error;
        }
    }

    if (CHECK_BASE_INDEX) check_base_index(data);

    return SUCCESS;

error:
    return ERROR;
}
