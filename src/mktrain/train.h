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
 * 2.14 2003/11/06 19:57:21
 */


#ifndef TRAIN_H
#define TRAIN_H 

#define MAX_CHUNKS      30
#define SCANS_PER_CHUNK 500
#define MIN_PEAKS_PER_CHUNK 30
#define OUTPUT_VECTOR_ALIGNMENT 0
#define VECTOR_MATCH_MULTIPLIER 3

/* Structure used to accumulate the statistics of errors
 * produced for a single read
 */
typedef struct {
    int align_length;           /* Length of the best-align region 
                                 * for a given read 
                                 */
    int count_ins;              /* Num of ins errs in the best-align region */
    int count_del;              /* Num of del errs in the best-align region */
    int count_sub;              /* Num of sub errs in the best-align region */
    int count_err;              /* Num of all errs in the best-align region */ 
    int count_hetero;           /* Number of heterozygotes in ref sequence  */  
    int count_tp_hetero;        /* Number of true  positive heterozygotes   */
    int count_fp_hetero;        /* Number of false positive heterozygotes   */
    int count_fn_hetero;        /* Number of false negative heterozygotes   */
    int count_mm_hetero;        /* Number of mismatched     heterozygotes   */
    int count_tp_qv_hetero;     /* Number of true  pos heterozygotes of QV>=20 */  
    int count_fp_qv_hetero;     /* Number of false pos heterozygotes of QV>=20 */
    int count_fn_qv_hetero;     /* Number of false neg heterozygotes of QV>=20 */
    int count_chunks;           /* Used size of arrays below */
    int min_qv_cor;             /* Lowest QV assigned to correctly called heterozygote */
    double max_phr3_cor;        /* Highest phr3 of correctly called heterozygote */
    double max_phr7_cor;        /* Highest phr7 of correctly called heterozygote */     
    double max_psr7_cor;        /* Highest psr7 of correctly called heterozygote */     
    double max_pres_cor;        /* Highest pres of correctly called heterozygote */     
    Stats stats[MAX_CHUNKS];    /* mean & std. dev. of normalized spacing */
    Histo histo[MAX_CHUNKS];    /* histogram of normalized spacing */
    double frac_QV20_with_shoulders; /* Fraction of QV>=20 peaks with shoulders */  
} Results;

#endif
