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
 *   $Id: lut.h,v 1.5 2009/01/16 15:15:23 gdenisov Exp $     
 */

#define BASE_COUNT_SCALE   (100000)
#define QVMAX(a,b)  (((a)>(b))?(a):(b))

typedef enum {
    CURRENT,
    PREVIOUS
} TIME;

typedef struct {
    char   is_match; /* 1 if base call is correct, 0 otherwise */
    char   schar;
    double parameter[6]; /* parameter values for a particular base */
    int    qv;
    char   new_frag_beg;
} BASE;

typedef struct {
    int     num_val;     /* length of the array of unique values */
    double *value;       /* array of unique parameter values */
    int    *weight;      /* array of weights (numbers of duplicates of a given
                          * unique value);
                          */
    double *threshold;   /* all possible threshold values for a given parameter */
    int threshold_count; /* number of thresholds (e.g., 50) for a parameter */
    int dimension;       /* precomputed dimension() */
    double max;
} PARAMETER;

typedef struct {
    unsigned long correct_base_call_count;
    unsigned long incorrect_base_call_count;
    unsigned long total_base_call_count;
    double        error_rate;
    int           quality_value;
    double        parameter[6];   /* threshold values for the cut */
    int           index[6];       /* threshold indices for the cut */
    int           sum_of_indices; /* sum of threshold indices from "index" array */
} HIGHEST_QV_CUT;

typedef struct {
    unsigned long correct;   /* number of correct base calls */
    unsigned long incorrect; /* number of incorrect base calls */
} BIN;

typedef struct {
    unsigned long correct;   /* number of correct base calls */
    unsigned long incorrect; /* number of incorrect base calls */
} CUT;

typedef struct {
    BIN *bin;
    CUT *previous_cube; /* the cut scores for n-1 dimensions */
    CUT *current_cube;
    CUT boundary_cut;
    int parameter_count;
    PARAMETER *parameter;
    int dimension[6];
} INFO;
