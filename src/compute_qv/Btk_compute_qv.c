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
 * $Id: Btk_compute_qv.c,v 1.17 2009/01/12 22:15:11 gdenisov Exp $
 */

/*
 *  Btk_compute_qv.c  $Revision: 1.17 $
 */

#include <math.h>   
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <limits.h>
#include <stdint.h>

#include "Btk_qv.h"
#include "util.h"
#include "Btk_lookup_table.h"
#include "Btk_default_table.h" 
#include "Btk_qv_data.h"
#include "context_table.h"     
#include "train.h"
#include "Btk_compute_qv.h"
#include "tracepoly.h"
#include "Btk_get_mixed_bases.h"
#include "SFF_Toolkit.h"
#include "Btk_compute_tpars.h"  /* needs train.h */
#include "Btk_atod.h"

/* Staff needed for license */
#define PRODNAMELEN (200)
#define MAX_NAME_LEN 1000
#define NUM_PARAMS  (4)

/*******************************************************************************
 * Function: quality_value
 *******************************************************************************
 */
uint8_t
get_quality_value(double phr3, double phr7, double psr7, double pres,
    BtkLookupTable *table)
{
    int i;

    if (table->num_lut_entries < 1) {
        return 0;
    }

    for (i = 0; i < table->num_lut_entries; i++) {
        if ((phr3 <= table->tpar[(int)table->entries[i].phr3i].phr3t) &&
            (phr7 <= table->tpar[(int)table->entries[i].phr7i].phr7t) &&
            (psr7 <= table->tpar[(int)table->entries[i].psr7i].psr7t) &&
            (pres <= table->tpar[(int)table->entries[i].presi].prest))
        {
            return table->entries[i].qval;
        }
    }

    /* If the right entry is not found, return the last entry in the table */
    return table->entries[table->num_lut_entries-1].qval;
}

/*******************************************************************************
 * Function: assign_quality
 *******************************************************************************
 */
int
assign_quality(int num_called_bases, uint8_t *quality_values, 
    double **params, Options *options, BtkLookupTable *table)
{
    int i;

    /* Calculate the quality values and output the results */
    for (i = 0; i < num_called_bases; i++) {
	quality_values[i] = get_quality_value(
        params[0][i], params[1][i], params[2][i], params[3][i], 
        table);
    }
    return SUCCESS;
}

/*******************************************************************************
 * Function: Btk_compute_qv
 *******************************************************************************
 */
int
Btk_compute_qv(int *num_called_bases, char **called_bases, 
               int **called_peak_locs, int *num_datapoints, int **chromatogram, 
               char *color2base, BtkLookupTable *table, 
#if USE_CONTEXT_TABLE
               ContextTable *ctable,
#endif
               uint8_t **quality_values, Options options, BtkMessage *message, 
               Results *results)
{
    int           i, r;
    double       *params[NUM_PARAMS]; 
    double       *iheight  = NULL;
    double       *iheight2 = NULL;
    double       *ave_iheight = NULL;
    ReadInfo      read_info;
#if !USE_CONTEXT_TABLE
    ContextTable *ctable = NULL;
#endif

    for (i = 0; i < NUM_PARAMS; i++) {
        params[i] = NULL;
    }

    /* Reset parameters depending on the lookup table used 
     * Default is ABI3100 table
     */
    options.sf[0] = 0.8; 
    options.sf[1] = 3.5;
    options.sf[2] = 0.1;
    options.sf[3] = 1.4;

    if ((!options.het && !options.mix) || options.recalln || options.recallndb 
        || options.ladder) 
    { 
        if (Btk_compute_tpars_Sanger(num_called_bases, called_bases, called_peak_locs, 
            num_datapoints, chromatogram, color2base, NUM_PARAMS,
            &params[0], &params[1], &params[2], &params[3], &iheight, &iheight2,
            &ave_iheight, &read_info, table, ctable, options, message, results) 
            != SUCCESS) 
        {
                goto error;
        }

        if (options.process_bases) 
        {
           *quality_values = REALLOC(*quality_values, uint8_t, *num_called_bases);
            MEM_ERROR(*quality_values);

            if ((r = assign_quality(*num_called_bases, *quality_values,
                params, &options, table)) != SUCCESS) {
                goto error;
            }

        } 
    }
  
    if (options.het || options.mix) 
    {
        /* Notes: 
         * 1) function Btk_get_mixed_bases will repeat most of the calls made by 
         *    Btk_compute_tpars
         * 2) If options.het > 0 AND options.poly > 0, then heterozygotes 
         *    will be processed as usual, but instead of mixed base they
         *    will be assigned a pure base character corresponding to the 
         *    highest peak (at Tim's request). These pure bases will be
         *    output to PHD file
         */
        if (Btk_get_mixed_bases(num_called_bases, called_bases, 
            called_peak_locs, *num_datapoints, chromatogram, 
            color2base, quality_values, &read_info, table, ctable, 
            options, message, results) != SUCCESS)
        {
            goto error;
        }
    }

    for (i = 0; i < NUM_PARAMS; i++) {
        FREE(params[i]);
    }
    FREE(iheight);
    FREE(iheight2);
    FREE(ave_iheight);
    return SUCCESS;

error:
    for (i = 0; i < NUM_PARAMS; i++) {
        FREE(params[i]);
    }
    FREE(iheight);
    FREE(iheight2);
    FREE(ave_iheight);

    return ERROR;
}
