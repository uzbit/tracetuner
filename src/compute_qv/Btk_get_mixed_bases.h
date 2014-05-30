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
 * 2.3 2003/11/06 18:18:38
 */

#if 0   /* The following files need to be included: */
#include "util.h"    // needed by train.h
#include "train.h"    // typedef Results
#endif

extern int get_mixed_base_position(Data *, int, int, Options *, BtkMessage *);
extern int is_mixed_base(char);
extern char mixed_base(char, char);
extern double get_average_called_peak_height(Data *, int);
extern int 
Btk_get_mixed_bases(
    int *,		/* pointer to input len. of the array of called bases */
    char **,		/* pointer to input array of called bases */
    int **,		/* pointer to input array of the called peak locations*/
    int,		/* input length of the chromatograms array */
    int **,		/* input arrays which store chromatographic data for
			 * each of dyes, assumed to be 4 colors 
                         */
    char *,		/* array of bases corresponding to the colors */
    uint8_t **,         /* pointer to output array of quality values  */
    ReadInfo *,         /* pointer to data used in trace renormalization 
                         * procedure 
                         */
    BtkLookupTable *,   /* pointer to a lookup table */
    ContextTable * ,    /* pointer to a context table */ 
    Options options,    /* structure including file_name, nocall, etc. */
    BtkMessage *,	/* error code and descriptive text */
    Results *           /* statistical results used by train (not ttuner) */
);
