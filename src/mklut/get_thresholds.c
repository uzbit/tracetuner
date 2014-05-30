
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

/**     Copyright (c) 1999 Paracel Inc.  All rights reserved.
 **
 **/

/*
 *   $Id: get_thresholds.c,v 1.6 2009/01/16 15:14:23 gdenisov Exp $
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>

#include "select.h"
#include "params.h"
#include "lut.h"
#include "Btk_atod.h"

#define MIN_BIN_FRACTION   0.25
#define BLOCKSIZE          1024
#define DISPLAY_THRESHOLD_DETAILS 0

extern int Verbose;     /* How much status info to print, if any */
extern int Compress;    /* Whether to compress thresholds */

/***************************************************************************
 * Function: get_thresholds
 *
 * Purpose: determine a set of thresholds for each parameter such that an
 *          even number of bases fall between each threshold.
 *
 * Called by: main in file lut.c
 * Calls:     randomized_select
 *
 ***************************************************************************/
void
get_thresholds(BASE *base, unsigned long base_count, PARAMETER parameter[],
    int PARAMETER_COUNT)
{
    int i, j, jj, k, kk, m, compressed, remain_weight;
    int *weight_value;
    int params_per_threshold, thresholds_per_value, unique_count, unique_sum;
    double *parameter_value;

    parameter_value = (double *) malloc(sizeof(double) * base_count);

    if (parameter_value == NULL)
    {
        puts("couldn't malloc parameter_value");
        exit(1);
    }

    for (i = 0; i < PARAMETER_COUNT; i++)
    {
        weight_value = (int *) malloc(sizeof(*weight_value) * BLOCKSIZE);
        if (weight_value == NULL)
        {
            puts("couldn't malloc weight_value");
            exit(1);
        }

        for (j = 0; j < (int)base_count; j++)
            parameter_value[j] = base[j].parameter[i];

        quicksort(parameter_value, 0, base_count - 1);
#if 0
        for (j=1; j< (int)base_count; j++) {
            if (parameter_value[j] < parameter_value[j-1])
                fprintf(stderr, "Unsorted values found: par[%d]=%f < par[%d]=%f\n",
                j, parameter_value[j], j-1, parameter_value[j-1]);
            if (j < 20 || j> (int)base_count-20 )
                fprintf(stderr, "param_value[%d]=%f\n", j, parameter_value[j]);
        }      
#endif

        for (j = 1; j < (int)base_count; j++) 
        {
            if (parameter_value[j] < parameter_value[j - 1])
            {
                fprintf(stderr, 
                    "\nparameter %d: %f is smaller than parameter %d: %f\n"
                    "\nFix quicksort and try again.\n\n",
                    j, parameter_value[j], j - 1, parameter_value[j - 1]);
                exit(1);
            }

        }
        
        /* Now find the UNIQUE parameter values, and collapse the array. */
        weight_value[0] = 1;

        for (j = 1, unique_count = 0; j < (int)base_count; j++) 
        {
            if (parameter_value[j] != parameter_value[j - 1])
            {
                if (++unique_count % BLOCKSIZE == 0)
                    weight_value = (int*)realloc(weight_value, 
                          sizeof(*weight_value) * (unique_count + BLOCKSIZE));

                parameter_value[unique_count] = parameter_value[j];
                weight_value[unique_count] = 1;
            }
            else
                weight_value[unique_count]++;
        }

        
#if 0
        for (j=0; j<=unique_count; j++) {
            if (i==3 && parameter_value[j] <= parameter_value[j-1])
                fprintf(stderr, "Unsorted values found: par[%d]=%f < par[%d]=%f\n",
                j, parameter_value[j], j-1, parameter_value[j-1]);
            if (i==3 && (j < 10 || j>unique_count-10) )
                fprintf(stderr, "param_value[%d]=%f\n", j, parameter_value[j]);
        }
#endif

        if (Verbose > 1)
            fprintf(stderr, "parameter #%d: %lu values, %d unique\n",
                                         i,  base_count, unique_count + 1);

        parameter[i].max = parameter_value[unique_count];

        for (m = base_count, k = -1, j = 0; 
             (j < parameter[i].threshold_count) && (k < unique_count-1); 
             j++)
        {
            /* Calculate number per bin with current parameters left 
             * and bins left 
             */
            params_per_threshold = (int)rint(((double) m) / 
                                         (parameter[i].threshold_count - j));
            unique_sum = weight_value[++k];

            if (unique_sum < params_per_threshold)
            {
                while (unique_sum < params_per_threshold)
                    unique_sum += weight_value[++k];

                if ((params_per_threshold - unique_sum + weight_value[k] >
                    unique_sum - params_per_threshold)
                ||
                    params_per_threshold - unique_sum + weight_value[k - 1] <
                    (rint) (((double) params_per_threshold) * MIN_BIN_FRACTION))
                {
                    parameter[i].threshold[j] = (parameter_value[k] 
                        + parameter_value[k + 1]) / 2;
                    m -= unique_sum;
                }
                else
                {
                    parameter[i].threshold[j] = (parameter_value[k - 1] 
                        + parameter_value[k]) / 2;
                    m -= unique_sum - weight_value[k--];
                }
            }
            else
            {
                thresholds_per_value = (int)floor(((double) weight_value[k]) / 
                                                  params_per_threshold);
                remain_weight = weight_value[k] 
                    - thresholds_per_value * params_per_threshold;

                for (jj = j; jj < j + thresholds_per_value; jj++)
                {
                    parameter[i].threshold[j] = (parameter_value[k] 
                        + parameter_value[k + 1]) / 2;
                    m -= params_per_threshold;
                }

                if (remain_weight > 0)
                {
                    parameter[i].threshold[j] = (parameter_value[k] 
                        + parameter_value[k + 1]) / 2;
                    m -= remain_weight;
                }
            }
        }
        parameter[i].threshold_count = j;

        parameter[i].threshold[parameter[i].threshold_count - 1] 
            = parameter[i].max;

        if (Compress)
        {
            /*
             * Compress the thresholds if any are duplicated.  Compress instead
             * of redistributing across the remaining values to limit the
             * skewness of the resulting bins.  For example, if we're picking
             * three thresholds from the data (0 0 0 0 0 0 1 2 3), we'd get 0,
             * 0, and 3.  Since the first two are 0, we want two buckets
             * instead, one with a value of 0 (6 bases), and the other with a
             * value of 3 (3 bases).  These buckets are somewhat skewed, but
             * they're the best we can do with this data.  If, instead, we
             * chose to redistribute, we'd get three buckets all right, one
             * with a value of 0 (6 bases), another with a value of 2 (2
             * bases), and the last one with a value of 3 (1 base), which is
             * even more skewed.
             */

            for (kk = 1, compressed = j = 0; 
                 j < parameter[i].threshold_count - 1; 
                 j++)
            {
                while ((j + kk < parameter[i].threshold_count-1) &&
                    (parameter[i].threshold[j] == 
                     parameter[i].threshold[j + kk])) {
                    kk++;
                }

                if (kk > 1)   /* got some dups; shift down and continue */
                {
                    assert(parameter[i].threshold_count - (j + kk) >= 0);
                    memmove(&parameter[i].threshold[j + 1], 
                        &parameter[i].threshold[j + kk],
                        (parameter[i].threshold_count - (j + kk)) *
                        sizeof(parameter[i].threshold[0]));
                    parameter[i].threshold_count -= kk - 1;
                    compressed += kk - 1;
                }
            }

            if (Verbose > 1 && compressed > 0)
                fprintf(stderr, 
                    "compressed down to %d thresholds for parm #%d\n",
                    parameter[i].threshold_count, i);
        }
        free(weight_value);
    }

    free(parameter_value);
}

/***************************************************************************
 * Function: get_thresholds2
 * 
 * Purpose: determine a set of thresholds for each parameter such that an 
 *          even number of bases fall between each threshold.
 *
 * Called by: main in file lut.c
 *
 ***************************************************************************
 */
void
get_thresholds2(PARAMETER parameter[], int PARAMETER_COUNT)
{
    int i, j=0, k=0;
    int expected_bin_size, sum, sum1, num_thresholds;     

    for (i = 0; i < PARAMETER_COUNT; i++)
    {

        fprintf(stderr, 
            "   parameter %d, num_val=%d ", 
            i, parameter[i].num_val);
        sum = 0;
        for (j = 0; j < parameter[i].num_val; j++) {
            sum += parameter[i].weight[j];
            if (j < parameter[i].num_val-1 &&
                    parameter[i].value[j] >= parameter[i].value[j+1]) {
                fprintf(stderr, "In get_thresholds2: parameters disordered!!!\n");
            }
        }

        /* Set the last threshold */
        num_thresholds = parameter[i].threshold_count;
        parameter[i].threshold[num_thresholds-1] = 
            parameter[i].value[parameter[i].num_val-1];
        parameter[i].max = parameter[i].threshold[num_thresholds-1];
        fprintf(stderr,
            " last threshold=%f\n", parameter[i].threshold[num_thresholds-1]);
        /* Compute the rest of the thresholds */
        j = 0;
        k = 0;
        while (k < parameter[i].threshold_count-1) {
            expected_bin_size = sum / num_thresholds;

            /* Case 1: one unique value per bin 
             **********************************
             */
            if (j < parameter[i].num_val-1 &&
                parameter[i].weight[j] >= expected_bin_size) 
            {
//              fprintf(stderr, "Case 1\n");
                parameter[i].threshold[k] = 
                    (parameter[i].value[j] + parameter[i].value[j+1])/2.; 
                k++;
                j++;
                num_thresholds--;
                sum -= parameter[i].weight[j];
#if DISPLAY_THRESHOLD_DETAILS
                fprintf(stderr, "Threshold=%d, param=%d, tot_weight=%d exp_size=%d...\n",
                    k, j, parameter[i].weight[j], expected_bin_size);
#endif
            }
            
            /* Case 2: one or more unique values per bin 
             *******************************************
             */
            else if (j < parameter[i].num_val-1 &&
                parameter[i].weight[j] < expected_bin_size) 
            {
//              fprintf(stderr, "Case 2\n");
                sum1 = parameter[i].weight[j];
                while (sum1 < expected_bin_size &&
                    j < parameter[i].num_val-1) 
                {
                    j++;
                    sum1 += parameter[i].weight[j];
                }
//              fprintf(stderr, "   j=%d sum1=%d\n", j, sum1);
                sum -= sum1;
                if (j < parameter[i].num_val-1 && sum1 >=expected_bin_size) {
                    parameter[i].threshold[k] =
                        (parameter[i].value[j] + parameter[i].value[j+1])/2.;
#if DISPLAY_THRESHOLD_DETAILS
                    fprintf(stderr, 
                        "Threshold=%d, param=%d, tot_weight=%d exp_size=%d...\n",
                        k, j, sum1, expected_bin_size);
#endif
                    k++;
                    num_thresholds--;
                }
                else if (j == parameter[i].num_val-1) {
                    parameter[i].threshold_count = k;
                    parameter[i].threshold[k-1] = parameter[i].value[j];
                    parameter[i].threshold = 
                        (double*)realloc(parameter[i].threshold, k*sizeof(double)); 
#if DISPLAY_THRESHOLD_DETAILS
                    fprintf(stderr, 
                        "Threshold=%d, param=%d, tot_weight=%d exp_size=%d...\n",
                        k, j, sum1, expected_bin_size);
#endif
                }        
            }
            /* No parameters left */
            else if (j == parameter[i].num_val-1) {
                parameter[i].threshold_count = k;
                parameter[i].threshold[k] = parameter[i].value[j];
                parameter[i].threshold =
                    (double*)realloc(parameter[i].threshold, k*sizeof(double));
            }
        } /* end loop in k */
    }     /* end loop in i */
}
