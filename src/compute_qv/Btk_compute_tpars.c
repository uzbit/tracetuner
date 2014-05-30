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
 */

/*
 *  Btk_compute_tpars.c  $Revision: 1.21 $
 */

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
#include "nr.h"
#include "util.h"
#include "Btk_lookup_table.h"
#include "Btk_qv_data.h"
#include "SFF_Toolkit.h"
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
#include "Btk_process_indels.h"
#include "Btk_get_mixed_bases.h"

#define MIN_NUM_PEAKS                5
#define SHOW_INPUT_OPTIONS           0
#define SORT_ORIG_LOCATIONS          1
#define STORE_IS_RESOLVED            0
#define STORE_CASE                   0

extern unsigned int max_colordata_value;

/*******************************************************************************
 * Function: colordata_release
 * Purpose: free the memory allocated by colordata_create()
 *******************************************************************************
 */
void
colordata_release(ColorData *color_data)
{
    FREE(color_data->data);
    FREE(color_data->peak_list);
}
 
/*******************************************************************************
 * Function: bases release
 * Purpose: free the memory allocated by bases_create()
 *******************************************************************************
 */
void
bases_release(TT_Bases *bases)
{
    FREE(bases->coordinate);
    FREE(bases->bases);
    FREE(bases->called_peak_list);
    bases->length = 0;
}
 
/*******************************************************************************
 * Function: trace_parameters_release
 * Purpose: free the memory allocated by trace_parameters_create()
 *******************************************************************************
 */
void
trace_parameters_release(TraceParameters *tp)
{
    FREE(tp->phr3);
    FREE(tp->phr7);
    FREE(tp->psr7);
    FREE(tp->pres);
}
/*******************************************************************************
 * Function: data_release
 * Purpose: free the memory allocated by data_create()
 *******************************************************************************
 */
void
data_release(Data *data)
{
    int i;
 
    for (i = 0; i < NUM_COLORS; i++) {
        colordata_release(&data->color_data[i]);
    }
    bases_release(&data->bases);
    trace_parameters_release(&data->trace_parameters);
 
    FREE(data->peak_list);
}
 
/*******************************************************************************
 * Function: colordata_create
 * Purpose: allocate memory for colordata structure
 *******************************************************************************
 */
int
colordata_create(ColorData *color_data, int length, int color,
    char *color2base, BtkMessage *message)
{
    color_data->length = length;
    color_data->peak_list_len = 0;
    color_data->peak_list_max_len = MAX_NUM_OF_PEAK;

    color_data->peak_list = CALLOC(Peak, color_data->peak_list_max_len);
    MEM_ERROR(color_data->peak_list);

    color_data->data = CALLOC(int, color_data->length);
    MEM_ERROR(color_data->data);
 
    color_data->dye_number = color + 1;
    color_data->base = color2base[color];
 
    return SUCCESS;
 
error:
    FREE(color_data->data);
    FREE(color_data->peak_list);
    return ERROR;
}
 
/*******************************************************************************
 * Function: bases_create
 * Purpose: allocate memory for bases structure
 *******************************************************************************
 */
int
bases_create(TT_Bases *bases, int length, BtkMessage *message)
{

    bases->length      = length;
    bases->orig_length = length;
    bases->max_length  = length * 2; 

    bases->coordinate = CALLOC(int, bases->max_length);
    MEM_ERROR(bases->coordinate);
 
    bases->bases = CALLOC(char, bases->max_length);
    MEM_ERROR(bases->bases);
 
    bases->called_peak_list = CALLOC(Peak *, bases->max_length);
    MEM_ERROR(bases->called_peak_list);
 
    return SUCCESS;
 
error:
    FREE(bases->coordinate);
    FREE(bases->bases);
    FREE(bases->called_peak_list);
    return ERROR;
}
 
/*******************************************************************************
 * Function: trace_parameters_create
 * Purpose: allocate memory for trace_parameters structure
 *******************************************************************************
 */
int
trace_parameters_create(TraceParameters *tp, int length, BtkMessage *message)
{
    tp->length = 0;      
    tp->phr3   = NULL;
    tp->phr7   = NULL;
    tp->psr7   = NULL;
    tp->pres   = NULL;

    return SUCCESS;
}
 
/*******************************************************************************
 * Function: data_create
 * Purpose: allocate memory for data structure
 *******************************************************************************
 */
int
data_create(Data *data, int length_cd, int length_bs, char *color2base, 
    BtkMessage *message)
{
    int i, r;
    (void)memset(data, 0, sizeof(data));

    data->length = 0; 
    for (i = 0; i < NUM_COLORS; i++) {
        if ((r = colordata_create(&data->color_data[i], length_cd, i, 
            color2base, message)) != SUCCESS)
        {
            goto error;
        }
        if (data->color_data[i].length > data->length) 
            data->length = data->color_data[i].length;
    }
 
    if ((r = bases_create(&data->bases, length_bs, message)) != SUCCESS) {
        goto error;
    }
 
    if ((r = trace_parameters_create(&data->trace_parameters, length_bs,
                                        message))
        != SUCCESS)
    {
        goto error;
    }

    data->peak_list_len     = length_cd * 4;
    data->peak_list_max_len = length_cd * 8;
 
    data->peak_list = CALLOC(Peak *, data->peak_list_max_len);
    MEM_ERROR(data->peak_list);
    
 
    return SUCCESS;
 
error:
    data_release(data);
    return ERROR;
}

/*******************************************************************************
 * Function: partition_locations
 * Purpose:  partition an array locs[p..r], that is, reorder the locs 
 *           and determine such an index q (p <= q <= r) that any location 
 *           with index <= q is smaller or equal to the location of any base 
 *           with index >q. Return the value q
 *******************************************************************************
 */
unsigned long
partition(int *locs, int p, int r)
{
    int  x, itemp;
    int  i, q;

    x = locs[p];
    i = p - 1;
    q = r + 1;
    while (1) {
        /* scan from right until (locs[q] <= x) is found */
        do {
            q = q - 1;
        } while (locs[q] > x);
        /* scan from left until (locs[i] >= x) is found */
        do {
            i = i + 1;
        } while (locs[i] < x);
        /* scan and exchange elements until scan indices cross,
           i.e., (i >= q) */
        if (i < q) {
          /* exchange elements of the array */
            itemp = locs[i];
            locs[i] = locs[q];
            locs[q] = itemp;
        }
        else
        {
            return q;
        }
    }
}

/*******************************************************************************
 * Function: quicksort_bases
 * Purpose:  Sort bases with respect to their locations, using
 *           a quicksort algorithm. Perform this
 *           for bases and locations with indeces between p and r
 *******************************************************************************
 */
void
quicksort_locations(int *locs, int p, int r)
{
    int q;

    if (p < r) {
        q = partition(locs, p, r);
        quicksort_locations(locs, p,   q);
        quicksort_locations(locs, q+1, r);
    }
}
 
/*******************************************************************************
 * Function: bases_populate
 * Purpose: populate the arrays allocated by bases_create()
 *******************************************************************************
 */
int
bases_populate(int *num_bases, char **bases, int edited_bases,
    int **locs, Data *data, Options *options, BtkMessage *message)
{
    int i, m=0;
 
    if (data->bases.length != *num_bases) {
        (void)sprintf(message->text,
                "Internal error: data->bases.length=%d; num_bases=%d",
                data->bases.length, *num_bases);
        return ERROR;
    }

    /* Sort bases with respect to their locations,
     * if they are edited (and thus can be "swapped")
     */
    if (((SORT_ORIG_LOCATIONS > 0) || (edited_bases==1)) && 
         !options->recalln && !options->recallndb && !options->ladder) 
    {
        quicksort_locations(*locs, 0, *num_bases-1);

        /* Determine the number of bases having zero locations */
        if ((*locs)[0] == 0) {
            while ((*locs)[m] == 0) {
                m++;
            }
        }
       *num_bases = *num_bases - m;
        data->bases.length = *num_bases;
    }
    else {
        m = 0;
    }

    /* Insert the bases with nonzero peak locations */
    for (i = 0; i < data->bases.length; i++) {
        data->bases.bases[i] = (*bases)[i+m];
        data->bases.coordinate[i] = (*locs)[i+m];
        data->bases.called_peak_list[i] = NULL;
    }

    /* Now bases and locs contain a garbage! */
    if (edited_bases) {
        data->bases.bases = REALLOC(data->bases.bases, char, *num_bases);
        data->bases.coordinate = REALLOC(data->bases.coordinate,
            int, *num_bases);
        data->bases.called_peak_list = REALLOC(data->bases.called_peak_list,
            Peak *, *num_bases);
    }

 
    return SUCCESS;
}
 
/*******************************************************************************
 * Function: colordata_populate
 * Purpose: populate the arrays allocated by colordata_create()
 *******************************************************************************
 */
int
colordata_populate(int num_datapoints, int **chromatogram, char *color2base,
                      Data *data, BtkMessage *message)
{
    int i, j, max_value;

    max_colordata_value = 0;

    for (j = 0; j < NUM_COLORS; j++) {
        data->color_data[j].base = color2base[j];
        /*
         * Fill in the data points.  Along the way, find the maximum data
         * value for this color.
         */
        max_value = 0;
        for (i = 0; i < data->color_data[j].length; i++) {
            if (chromatogram[j][i]>=0) {
                data->color_data[j].data[i] = chromatogram[j][i];
            }
            else {
//              data->color_data[j].data[i] = 0;
#if 1
                data->color_data[j].data[i] = - chromatogram[j][i];
#endif

            }
            if (data->color_data[j].data[i] > max_value) {
                max_value = data->color_data[j].data[i];
            }
        }

        data->color_data[j].max_value = max_value;

        if (max_value > (int)max_colordata_value)
            max_colordata_value = max_value;
    }

    return SUCCESS;
}

/*******************************************************************************
 * Function: output_data
 * Purpose: write out the contents of the populated Data structure to
 *          a file for debugging purposes
 *******************************************************************************
 */
static void
output_data(Data *data, BtkMessage *message)
{
    FILE *fp;
    int c, i;
 
    if ((fp = fopen("Data", "w")) == NULL) {
        (void)perror("Data");
        return;
    }
 
    (void)fprintf(fp, "base  coord   pos\n");
    for (i = 0; i < data->bases.length; i++) {
        (void)fprintf(fp, "   %c  %5d\n", data->bases.bases[i],
                        data->bases.coordinate[i]);
    }
 
    for (c = 0; c < NUM_COLORS; c++) {
        (void)fprintf(fp, "\ncolor %d => base '%c'; max value = %d; values:\n",
                        c, data->color_data[c].base,
                        data->color_data[c].max_value);
        for (i = 0; i < data->color_data[c].length; i++) {
            (void)fprintf(fp, "%6d\n", data->color_data[c].data[i]);
        }
    }
 
    (void)fclose(fp);
}
 
/*******************************************************************************
 * Function: data_populate
 * Purpose: populate the arrays allocated by data_create()
 *******************************************************************************
 */
int
data_populate(int *num_bases, char **bases, int edited_bases,
    int **locs, int num_datapoints, int **chromatogram, 
    char *color2base, Data *data, Options *options, BtkMessage *message)
{
    int r;
 
    if ((r = bases_populate(num_bases, bases, edited_bases,
        locs, data, options, message)) == ERROR)
    {
        return r;
    }
 
    if ((r = colordata_populate(num_datapoints, chromatogram, color2base, data,
        message)) == ERROR)
    {
        return r;
    }
 
    if (getenv("DEBUG") != NULL) {
        output_data(data, message);
    }

    data->pos_data_beg = 0;
    data->pos_data_end = num_datapoints;
 
    return SUCCESS;
}

/*******************************************************************************
 * Function: is_resolved
 * Purpose: Set peak.is_resolved parameter for each peak in the list
 *          of peaks of given color (for the definition of resolved peak,
 *          see Phred paper, part II, page 188)
 *******************************************************************************
 */
int
is_resolved(ColorData *color_data, int cd_peak_ind)
{
    int j;
    int min_height;
    int i=(cd_peak_ind>0)?cd_peak_ind:1; 
    int left_resolved=0, right_resolved=0;
 
    min_height = QVMIN(color_data->peak_list[i-1].height,
                       color_data->peak_list[i  ].height);
    for (j = color_data->peak_list[i-1].pos+1;
         j < color_data->peak_list[i  ].pos  ;
         j++)
    {
        if (color_data->data[j] < min_height) {
            left_resolved=1;
            break;
        }
    }
    if ((cd_peak_ind==0) || 
        (cd_peak_ind==color_data->peak_list_len-1))
        return left_resolved;

    min_height = QVMIN(color_data->peak_list[i+1].height,
                       color_data->peak_list[i  ].height);
    for (j = color_data->peak_list[i  ].pos+1;
         j < color_data->peak_list[i+1].pos  ;
         j++)
    {
        if (color_data->data[j] < min_height) {
            right_resolved = 1;
            break;
        }
    }
    if (left_resolved && right_resolved)
        return 1;

    return 0;
}

/*******************************************************************************
 * Function: populate_iheight_array
 *
 *******************************************************************************
 */
static int
populate_iheight_array(Data* data, double* iheight, double* iheight2,
     double* ave_iheight, BtkMessage *msg)
{
    int i, base_index, j, jc;
    int position = -1;
    Peak peak;
    int peak_index;

    for (i=0; i<data->peak_list_len; i++) 
    {
        if (!data->peak_list[i]->is_called) {
            continue;
        }
        base_index = data->peak_list[i]->base_index;
    
        /* Determine iheight */ 
        if (STORE_IS_RESOLVED && !STORE_CASE) {
            if (base_index >= 0) {
                int j = data->peak_list[i]->color_index;
                iheight[base_index] = is_resolved(&data->color_data[j],
                   data->peak_list[i]->cd_peak_ind);
            }
        }
        else if (!STORE_IS_RESOLVED && STORE_CASE) {
            if (base_index >= 0)
                iheight[base_index] = (double)data->peak_list[i]->is_called;
        }
        else {
            if (base_index >= 0) 
                iheight[base_index] = data->peak_list[i]->iheight;
        }

        /* Determine iheight2 */
        position = data->peak_list[i]->ipos;
        jc = data->peak_list[i]->color_index;
        iheight2[base_index] = 0.;
        for (j=0; j<NUM_COLORS; j++) 
        {
            if (j == jc)
                continue;

            if (colordata_find_peak_index_by_location(&(data->color_data[j]),
                position, &peak, &peak_index, msg) != SUCCESS) {
                return ERROR;
            }

            if (peak_index < 0)
                continue;

            if (is_dye_blob(peak.ipos, &peak, data, 0))
                continue;

            if ((base_index<data->bases.length-1) &&
                (peak.pos >= (position +
                    data->bases.called_peak_list[base_index+1]->pos)/2))
                continue;

            if ((base_index>0) && (peak.pos < 
               (position + data->bases.called_peak_list[base_index-1]->pos)/2))
                continue;

            if (iheight2[base_index] < peak.iheight)
                iheight2[base_index] = peak.iheight;
        }
        if (base_index>data->bases.length-1) base_index=data->bases.length-1;
        ave_iheight[base_index] = get_average_called_peak_height(data, 
            base_index);
    }
    return SUCCESS;
}

/*******************************************************************************
 * Function: show_input_options
 *******************************************************************************
 */
void
show_input_options(Options *options)
{
    fprintf(stderr, "\nInput options:\nFile_name=%s lut_type=%d path=%s\n",
        options->file_name, options->lut_type, options->path);
    fprintf(stderr, "edited_bases=%d gauss=%d nocall=%d respace=%d\n",
        options->edited_bases, options->gauss, options->nocall,
        options->respace);
    fprintf(stderr,
    "recalln=%d renorm=%d shift=%d het=%d tab_dir=%s tal_dir=%s tip_dir=%s Verbose=%d\n",
        options->recalln, options->renorm, options->shift, options->het,
        options->tab_dir, options->tal_dir, options->tip_dir, options->Verbose);
    fprintf(stderr, "scaling factors = %f %f %f %f\n\n",
        options->sf[0], options->sf[1], options->sf[2], options->sf[3]);
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


static int
get_trace_parameters_of_mixed_base(int i, int data_peak_ind, int data_peak_ind2, 
    double **params, double r, Options *op, Data *data, BtkMessage *msg)
{
    double psr1, psr2;
    
    params[0][i] = QVMAX(
        get_peak_height_ratio(data, data_peak_ind, data_peak_ind2, 3, 
                   1./(1.+r), 0., op, msg),
        get_peak_height_ratio(data, data_peak_ind2, data_peak_ind, 3, 
                   r/(1.+r), 0., op, msg));
    params[1][i] = QVMAX(
        get_peak_height_ratio(data, data_peak_ind, data_peak_ind2, 7,
                   1./(1.+r), 0., op, msg),
        get_peak_height_ratio(data, data_peak_ind2, data_peak_ind, 7,
                   r/(1.+r), 0., op, msg));
        psr1 = get_peak_spacing_ratio(data, data_peak_ind, data_peak_ind2, 7,
                   op, msg);
        psr2 = get_peak_spacing_ratio(data, data_peak_ind2, data_peak_ind, 7,
                   op, msg);
    params[2][i] = QVMAX(psr1, psr2);
#if 0
    fprintf(stderr, "mixed_base= %c data_peak_ind 1,2 = %d %d psr[%d] = %f %f\n",
        mixed_base(data->peak_list[data_peak_ind]->base,
        data->peak_list[data_peak_ind2]->base), data_peak_ind, data_peak_ind2, i, psr1, psr2);
#endif
    params[3][i] = QVMAX(
        get_peak_resolution(data, data_peak_ind , op, msg),
        get_peak_resolution(data, data_peak_ind2, op, msg));
#if 0
    fprintf(stderr, "params = %f %f %f %f \n",
        params[0][i], params[1][i], params[2][i], params[3][i]);
#endif

    return SUCCESS;
}

static int
get_data_peak_ind_of_second_base(Data *data, int data_peak_ind, char base2,
    char *color2base, int *color2)
{
    Peak peak;
    int pind;
    BtkMessage msg;
   *color2 = (color2base[0] == base2) ? 0 :
             (color2base[1] == base2) ? 1 :
             (color2base[2] == base2) ? 2 : 3;
    colordata_find_peak_index_by_location(&(data->color_data[*color2]),
        data->peak_list[data_peak_ind]->pos, &peak, &pind, &msg);
    if (pind >= 0)
        return data->color_data[*color2].peak_list[pind].data_peak_ind;
    else
        return -1;
}

/*******************************************************************************
 * Function: Btk_compute_tpars_Sanger
 *           Without option -raw processing occurs as usual
 *           With -raw, raw data will be analyzed and then, depending on other
 *                output options, bases may or may not be called
 *******************************************************************************
 */
int
Btk_compute_tpars_Sanger(int *num_bases, char **bases, int **peak_locs, 
    int *num_datapoints, int **chromatogram, char *color2base, int NUM_PARAMS,
    double **params0, double **params1, double **params2, double **params3,  
    double **iheight, double **iheight2, double **ave_iheight, 
    ReadInfo *read_info, BtkLookupTable *table, ContextTable *ctable,
    Options options, BtkMessage *message, Results *results )
{
    int      i;
    Data     data; 
    double  *params[4] = {NULL, NULL, NULL, NULL}; 
    clock_t  start_clock = clock(), curr_clock;
    int     *quality_values = NULL;
    char    *orig_bases = "";

    if (SHOW_INPUT_OPTIONS)
        show_input_options(&options);

    if (*num_bases > MAX_NUM_BASES)
    {
        fprintf(stderr,
                "Number of called bases exceeds allowed maximum\n");
        return ERROR;
    }

    /* Store original bases */
    if (options.het || options.mix)
    {
        orig_bases = CALLOC(char, *num_bases);
        for (i=0; i<*num_bases; i++)
        {
            orig_bases[i] = (*bases)[i];
        }
    }

    if (data_create(&data, *num_datapoints, *num_bases, 
        color2base, message) != SUCCESS)
    {
        sprintf(message->text, "Error calling data_create\n");
        return ERROR;
    }

    if (data_populate(num_bases, bases, options.edited_bases,
        peak_locs, *num_datapoints, chromatogram, color2base, &data,
        &options, message) != SUCCESS)
    {
        sprintf(message->text, "Error calling data_populate\n");
        fprintf(stderr, "Error calling  data_populate - 1\n");
        goto error;
    }

    if (options.raw_data) {
        if (Btk_process_raw_data(num_datapoints, chromatogram, "ACGT",
                 &data, options, message) == ERROR)
        {
            fprintf(stderr, "Error processing raw data\n");
            goto error;
        }

        if (data_populate(num_bases, bases, options.edited_bases,
            peak_locs, *num_datapoints, chromatogram, color2base, &data,
            &options, message) != SUCCESS)
        {
            sprintf(message->text, "Error calling data_populate\n");
            fprintf(stderr, "Error calling  data_populate - 2\n");
            goto error;
        }
    }
    else {
        if (options.xgr) {
            output_chromatogram("0_Orig_analyzed_data.xgr",
            "Analyzed data", 
            chromatogram[0], chromatogram[1], chromatogram[2], chromatogram[3], 
           *num_datapoints, &data);
        }
    }

    if (options.process_bases) 
    {
        if (options.time) {
            curr_clock = clock();
            fprintf(stderr, "Data structure populated in %f sec. \n",
                (float)(curr_clock - start_clock)/(float)CLOCKS_PER_SEC);
            start_clock = curr_clock;
        }
    
        if (Btk_process_peaks(&data, &options, message) != SUCCESS) {
            if (total_number_of_peaks(&data) >= 10) {
                sprintf(message->text, "Error calling process_peaks\n");
                fprintf(stderr, "Error calling process_peaks\n");
            }
	    fprintf(stderr, "Error processing peaks\n");
            goto error;
        }

        /* Return if data contain very few peaks */
        for (i=0; i<NUM_COLORS; i++)
        {
            if ((&data.color_data[i] == NULL) ||
                ( data.color_data[i].peak_list_len < MIN_NUM_PEAKS))
            {
                data_release(&data);
                fprintf(stderr, "This trace is classified as junk. Exit.\n");
                return ERROR;
            }
        }

        if (options.time) {
            curr_clock = clock();
            fprintf(stderr, "Peaks processed in %f sec. \n",
                (float)(curr_clock - start_clock)/(float)CLOCKS_PER_SEC);
            start_clock = curr_clock;
        }

        if (!options.recalln && !options.recallndb && !options.ladder &&
            !options.het     && !options.mix) 
            // the last two options may be passed only from train.c
        {
            if (Btk_call_bases(&data, color2base, read_info, ctable,
                &options, message, results ) != SUCCESS) {
                fprintf(stderr, "Error calling bases\n");
                goto error;
            }
            if (options.time) {
                curr_clock = clock();
                fprintf(stderr, "Bases called in %f sec. \n",
                    (float)(curr_clock - start_clock)/(float)CLOCKS_PER_SEC);
                start_clock = curr_clock;
            }
        }
        else if (options.recalln || options.recallndb || options.ladder ||
                 options.het     || options.mix)
        {
            if (options.recalln && (options.Verbose > 1)) {
                fprintf(stderr, "Using original base calls with Ns recalled ");
                fprintf(stderr, "to the best guess\n");
            }
            else if (options.recallndb && (options.Verbose > 1)) {
                fprintf(stderr,
                    "Using original base calls with Ns and dye blobs recalled ");
                fprintf(stderr, "to the best guess\n");
            }
             else if (options.ladder && (options.Verbose > 1)) {
                fprintf(stderr,
                    "Recalling all the bases to the best guess ");
                fprintf(stderr, "at their original locations");
            }

            if (options.ladder)
            {
               preset_base_calls(&data, color2base, options, message);
            }

            if (mark_called_peaks(&data, color2base, options, message)
                != SUCCESS)
            {
                goto error;
            }
            for (i=0; i<data.bases.length; i++)
            {
                if ((data.bases.called_peak_list[i]->ipos < 1) ||
                    (data.bases.called_peak_list[i]->pos  < 1))
                {
                    data.bases.called_peak_list[i]->ipos = 1;
                    data.bases.called_peak_list[i]->pos  = 1;
                }
            }
            if (options.time) {
                curr_clock = clock();
                fprintf(stderr, "Bases called in %f sec. \n",
                    (float)(curr_clock - start_clock)/(float)CLOCKS_PER_SEC);
                start_clock = curr_clock;
            }
        }

        // Retain the original mixed bases (this will only be used by train)
        if (options.het || options.mix)
        {
            for (i=0; i< *num_bases; i++)
            {
                if (is_mixed_base(orig_bases[i]))
                {
                    data.bases.bases[i] = orig_bases[i];
                }
            }
        }

        data.trace_parameters.length = data.bases.length;
        data.trace_parameters.phr3 = CALLOC(
            double, data.trace_parameters.length);
        MEM_ERROR(data.trace_parameters.phr3);
        data.trace_parameters.phr7 = CALLOC(
            double, data.trace_parameters.length);
        MEM_ERROR(data.trace_parameters.phr7);
        data.trace_parameters.psr7 = CALLOC(
            double, data.trace_parameters.length);
        MEM_ERROR(data.trace_parameters.psr7);
        data.trace_parameters.pres = CALLOC(
            double, data.trace_parameters.length);
        MEM_ERROR(data.trace_parameters.pres);

        for (i = 0; i < NUM_PARAMS; i++)
        {
            params[i] = CALLOC(double, data.bases.length);
            MEM_ERROR(params[i]);
        }
   
        if (*num_bases != data.bases.length) {
            *num_bases  = data.bases.length;
            data.trace_parameters.length = data.bases.length;

//         *bases = REALLOC(*bases, char, data.bases.length);
//          MEM_ERROR(*bases);
//         *peak_locs = REALLOC(*peak_locs, int, data.bases.length);
//          MEM_ERROR(*peak_locs);
        }
 
       *params0 = CALLOC(double, data.bases.length);
       *params1 = CALLOC(double, data.bases.length);
       *params2 = CALLOC(double, data.bases.length);
       *params3 = CALLOC(double, data.bases.length);
       *iheight = CALLOC(double, data.bases.length);
       *iheight2= CALLOC(double, data.bases.length);
       *ave_iheight = CALLOC(double, data.bases.length);
       *bases   = REALLOC(*bases, char, data.bases.length);
       *peak_locs = REALLOC(*peak_locs, int, data.bases.length);

        for (i = 0; i < data.bases.length; i++) {
            (*bases)[i] = data.bases.bases[i];
            (*peak_locs)[i] = data.bases.called_peak_list[i]->ipos;
        }

        if (!options.indel_detect && !options.indel_resolve) 
        {
            if (options.het || options.mix)
            {
                Peak **dpl = data.peak_list;
                Peak **cpl = data.bases.called_peak_list;

                // Update trace parameters for mixed bases
                for (i = 0; i < data.bases.length; i++) 
                {
                    if (is_mixed_base(data.bases.bases[i]))
                    {
                        char mbase = data.bases.bases[i];
                        int data_peak_ind = cpl[i]->data_peak_ind;
                        int color2;
                        int pos = cpl[i]->ipos;
                        char base = cpl[i]->base;
                        char base2= 
                           (mixed_base(base, 'A') == mbase) ? 'A' :
                          ((mixed_base(base, 'C') == mbase) ? 'C' :
                          ((mixed_base(base, 'G') == mbase) ? 'G' :
                          ((mixed_base(base, 'T') == mbase) ? 'T' :
                           'N')));         
                        int data_peak_ind2 = 
                            get_data_peak_ind_of_second_base(&data, data_peak_ind, 
                                base2, color2base, &color2);
                        data.peak_list[data_peak_ind ]->data_peak_ind2 = data_peak_ind2;
                        data.peak_list[data_peak_ind2]->data_peak_ind2 = data_peak_ind;
                        data.peak_list[data_peak_ind2]->base_index = i;
                        data.peak_list[data_peak_ind2]->is_called = 1;
                        double h1 = dpl[data_peak_ind ]->iheight;
                        double h2 = (data_peak_ind2 >= 0) ? dpl[data_peak_ind2]->iheight : 
                                    (color2 >= 0) ? data.color_data[color2].data[pos] : 1.;
                        double ratio = h1/h2;
#if 0 
                        fprintf(stderr, "First peak ind= %d Second peak ind= %d\n", data_peak_ind, data_peak_ind2);
#endif
 
                        if (data_peak_ind < 0 || data_peak_ind2 < 0)
                            fprintf(stderr, 
                                "Error: in Btk_compute_tpars: data_peak_ind=%d data_peak_ind2=%d\n",
                                data_peak_ind, data_peak_ind2);

                        get_trace_parameters_of_mixed_base(i, data_peak_ind,
                            data_peak_ind2, params, ratio, &options, &data, message);
#if 0
                        fprintf(stderr, "i= %d params_mixed = %f %f %f %f \n", i, params[0][i], params[1][i], params[2][i], params[3][i]);
#endif
                    } 
                }
#if 0
           fprintf(stderr, "params_mixed_0 = %f %f %f %f \n", params[0][i], params[1][i], params[2][i], params[3][i]);
#endif
                for (i = 0; i < data.bases.length; i++)
                {
                    if (is_mixed_base(data.bases.bases[i]))
                        continue;

                    if (get_trace_parameters_of_pure_bases(i, &data, color2base, &options,
                        message) != SUCCESS) {
                            sprintf(message->text, "Error calling get_trace_parameters_of_pure_bases\n");
                            goto error;
                    }

                    if (populate_params_array(i, &data, params) != SUCCESS) {
                        sprintf(message->text, "Error calling populate_params_array\n");
                        goto error;
                    }
                }
#if 0
           fprintf(stderr, "params_mixed_1 = %f %f %f %f \n", params[0][i], params[1][i], params[2][i], params[3][i]);
#endif
            }
            else
            {
                for (i = 0; i < data.bases.length; i++)
                {
                    if (get_trace_parameters_of_pure_bases(i, &data, color2base, &options,
                            message) != SUCCESS) {
                        sprintf(message->text, "Error calling get_trace_parameters_of_pure_bases\n");
                        goto error;
                    }
                    if (populate_params_array(i, &data, params) != SUCCESS) {
                        sprintf(message->text, "Error calling populate_params_array\n");
                        goto error;
                    }
                }
            }
        }
        else
        {
            if (Btk_compute_tp(&data, color2base, 4, params, read_info, ctable,
                &options, message) != SUCCESS) {
                fprintf(stderr, "Error calling Btk_compute_tp\n");
                goto error;
            }

            quality_values = CALLOC(int, data.bases.length);
            if (Btk_compute_tp(&data, color2base, NUM_PARAMS, params, read_info, 
                ctable, &options, message) != SUCCESS) {
                fprintf(stderr, "Error calling Btk_compute_tp\n");
                goto error;
            }

            for (i=0; i<data.bases.length; i++)
            {
                quality_values[i] = get_quality_value(params[0][i],
                                params[1][i], params[2][i],
                                params[3][i], table);
            }
            if (Btk_process_indels(options.file_name, num_datapoints, 
                chromatogram, color2base, quality_values, &data, read_info, 
                ctable, options, message) != SUCCESS)
                goto error;
 
            for (i = 0; i < data.bases.length; i++) {
                (*bases)[i] = data.bases.bases[i];
                (*peak_locs)[i] = data.bases.called_peak_list[i]->ipos;
            }

            FREE(quality_values);
        }
#if 0
        fprintf(stderr, "params_mixed_2 = %f %f %f %f \n", params[0][357], params[1][357], params[2][357], params[3][357]);
#endif
        /* Output .poly file */
        if (options.poly) {
            if (Btk_output_poly_file(&data, &options, message) != SUCCESS)
            {
                fprintf(stderr, "Error creating  poly file\n");
                goto error;
            }
        }

        if (options.time) {   
            curr_clock = clock();
            fprintf(stderr, "Trace parameters computed in %f sec. \n",
                (float)(curr_clock - start_clock)/(float)CLOCKS_PER_SEC);
            start_clock = curr_clock;
        }     

        if (populate_iheight_array(&data, *iheight, *iheight2, *ave_iheight,
            message) != SUCCESS) 
        {
            fprintf( stderr, "Error calling populate_iheight_array\n");
	    goto error;
        }
    
        for (i=0; i<*num_bases; i++) {
            (*params0)[i] = params[0][i];
            (*params1)[i] = params[1][i];
            (*params2)[i] = params[2][i];
            (*params3)[i] = params[3][i];
        }


#if 0
        fprintf(stderr, "params_mixed_3 = %f %f %f %f \n", (*params0)[357], (*params1)[357], (*params2)[357], (*params3)[357]);
#endif
        /* Count the number of peaks with shoulders */
        {
            int j, num_QV20_with_shoulders = 0;
            int num_QV20 = 0;
            for (i=0; i<NUM_COLORS; i++)
            {
                for (j=1; j<data.color_data[i].peak_list_len; j++)
                {
                    if (data.color_data[i].peak_list[j-1].is_called)
                    {
                        int base_ind =
                            data.color_data[i].peak_list[j-1].base_index;
                        int qv = get_quality_value((*params0)[base_ind],
                                                   (*params1)[base_ind],
                                                   (*params2)[base_ind],
                                                   (*params3)[base_ind],
                                                    table);
                        if (qv >= 20) num_QV20++;
                        if ((qv >= 20) &&
                            !data.color_data[i].peak_list[j].is_called &&
                            (data.color_data[i].peak_list[j].type >= 31))
                        {
                            num_QV20_with_shoulders++;
                        }
                    }
                }
            }
//          for (i=0; i<data.bases->length; i++)
//          {
//              int qv = get_quality_value((*params0)[i], (*params1)[i], (*params2)[i],
//                           (*params3)[i], table);
//              if (qv >= 20) num_QV20++;
//          }
           results->frac_QV20_with_shoulders = (float)num_QV20_with_shoulders*100.
                                              /(float)(num_QV20>0 ? num_QV20 : 1);
        }

#if 0
        fprintf(stderr, "params_mixed_last = %f %f %f %f \n", (*params0)[357], (*params1)[357], (*params2)[357], (*params3)[357]);
#endif

        if (options.tip_dir[0] != '\0') {
            if (Btk_output_tip_file(&data, color2base, options)
                != SUCCESS)
            {
                fprintf(stderr, "Error calling Btk_output_tip_file\n");
                goto error;
            }
        }
        for (i = 0; i < NUM_PARAMS; i++) {
            FREE(params[i]);
        }
    }
    data_release(&data);
    if (options.het || options.mix)
        FREE(orig_bases);

    return SUCCESS;

error:
    data_release(&data);
    for (i = 0; i < NUM_PARAMS; i++) {
        FREE(params[i]);
    }
    if (options.het || options.mix)
        FREE(orig_bases);

    return ERROR;
}

static int
base2color(char b)
{
    if      (b == 'A') return 0;
    else if (b == 'C') return 1;
    else if (b == 'G') return 2;
    else if (b == 'T') return 3;
    return -1;
}

static float
get_read_noise(float *call0_flows, int num0, float *call1_flows, int num1)
{
    int   i;
    float m0 =.0, m1 = 0.; 
    float s0 =.0, s1 = 0.;

    for (i=0; i<num0; i++) 
    {
        m0 += call0_flows[i];
        s0 += call0_flows[i]*call0_flows[i];
    }
 
    for (i=0; i<num1; i++)
    {
        m1 += call1_flows[i];
        s1 += call1_flows[i]*call1_flows[i];
    }

    m0 /= (float)num0;
    s0 /= (float)num0;
    m1 /= (float)num1;
    s1 /= (float)num1;

    s0 = sqrt(s0 - m0*m0);
    s1 = sqrt(s1 - m1*m1);

//  printf("m0= %f m1= %f s0= %f s1= %f\n", m0, m1, s0, s1);

    return (m1 - m0)/(s1 + s0);
}

/*******************************************************************************
 * Function: Btk_compute_tpars_454
 *
 * Trace parameters for 454 data are:
 *
 * 1 - noise in a given flow
 * 2 - max noise in a radius of 10 flows around the current one
 * 3 - read noise = (m1 _ m0)/(s0 + s1)
 * 4 - current homopolymer count 
 * 5 - previous homopolymer count for the same base as current one
 * 6 - position on read
 *
 *******************************************************************************
 */

int
Btk_compute_tpars_454(sffHeader *h, sffRead *r, double *params[6],
    double *flow_per_base, uint8_t *orig_qvs, Options options,
    BtkMessage *message )
{
    int    i, j, k, n, bind, num_hprs = 1;
    int    win_beg, win_end;
    float *noise_per_flow;
    float *max_noise;
    float *call0_flows = CALLOC(float, h->number_of_flows_per_read);
    float *call1_flows = CALLOC(float, h->number_of_flows_per_read);
    int    num_call0_flows = 0, num_call1_flows = 0;
    int    prev_hp_count[4] = {0, 0, 0, 0};
    float  read_noise = 0.;   

    // Count hprs
    for (i=1; i<r->number_of_bases; i++)
        if (r->bases[i] != r->bases[i-1])
            num_hprs++;

    // Compute noise per flow
    noise_per_flow  = CALLOC(float, num_hprs);
    bind = 0;
    j = 0;
    for (i=0; i<h->number_of_flows_per_read; i++)
    {
        if (h->flow_chars[i] == r->bases[bind])
        {
            n = 1;
            while (r->bases[bind+n] == r->bases[bind])
                n++;

            noise_per_flow[j] = fabs((float)r->flowgram_values[i]/100.0 - (float)n);
            j++;
            bind += n;
            
            if (n == 1)
            {
                call1_flows[num_call1_flows] = (float)r->flowgram_values[i]/100.0;
                num_call1_flows++;
            }
        }
        else
        {
            call0_flows[num_call0_flows] = (float)r->flowgram_values[i]/100.0;
            num_call0_flows++;
        }
    }

    // Compute max noise in a window of 21 flows
    max_noise = CALLOC(float, num_hprs);
    for (i=0; i<num_hprs; i++)
    {
        win_beg = i;
        while (i - win_beg <= 10 && win_beg > 0)
            win_beg--;

        win_end = i;
        while (win_end - i >= 10 && win_end < num_hprs-1)
            win_end++;

        max_noise[i] = noise_per_flow[i];
        for (j=win_beg; j<= win_end; j++)
        {
            if (max_noise[i] < noise_per_flow[j])
                max_noise[i] = noise_per_flow[j];
        }
    }

    // Compute the noise per read
    read_noise = get_read_noise(call0_flows, num_call0_flows, call1_flows, num_call1_flows);

    // Compute predictors of quality
    j = 0;       // index of a homopolymer run
    bind = 0;    // base index
    for (i=0; i<h->number_of_flows_per_read; i++)
    {
        if (toupper(h->flow_chars[i]) == r->bases[bind])
        {
            params[0][bind] = noise_per_flow[j];
            params[1][bind] = max_noise[j];
            params[2][bind] = read_noise;  // noise in a read
            params[4][bind] = base2color(r->bases[bind]) >= 0 ?
                prev_hp_count[base2color(r->bases[bind])] : 10;
            params[5][bind] = bind; // position on read
            flow_per_base[bind] = (double)r->flowgram_values[i]/100.0;
            orig_qvs[bind]      = (double)r->quality_values[bind];
            n = 1;
            while (r->bases[bind] == r->bases[bind+n])
            {
                params[0][bind+n] = params[0][bind];
                params[1][bind+n] = params[1][bind];
                params[2][bind+n] = params[2][bind];
                params[4][bind+n] = params[4][bind];
                params[5][bind+n] = bind+n;
                orig_qvs[bind+n]  = (double)r->quality_values[bind+n];
                n++;
            }

            flow_per_base[bind] /= (float)n;
            for (k=0; k<n; k++) {
                params[3][bind+k] = (double)n;    // current homopolymer count
                if (k > 0)
                    flow_per_base[bind+k] = flow_per_base[bind]; 
            }
            prev_hp_count[base2color(r->bases[bind])] = n;

            bind += n;
            j++;
        }
    }
    FREE(noise_per_flow);
    FREE(max_noise);
    FREE(call0_flows);
    FREE(call1_flows);
    return SUCCESS;
}

