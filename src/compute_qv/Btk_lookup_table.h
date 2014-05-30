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
 * 1.8 2003/11/06 18:18:39
 */

typedef struct _btk_quality_lookup_entry {
    char  qval ;    /* quality value */
    char  phr3i;    /* index of peak height ratio 3 threshold */
    char  phr7i;    /* index of peak height ratio 7 threshold */
    char  psr7i;    /* index of peak spacing  ratio threshold */
    char  presi;    /* index of peak resolution     threshold */
} BtkLookupEntry;

typedef struct _trace_param_threshold_entry {
    double phr3t;    /* peak height ratio 3 threshold */
    double phr7t;    /* peak height ratio 7 threshold */
    double psr7t;    /* peak spacing  ratio threshold */
    double prest;    /* peak resolution     threshold */
} TraceParamEntry;

typedef struct _btk_quality_lookup_table {
    int              num_tpar_entries;
    TraceParamEntry *tpar;
	int              num_lut_entries;
	BtkLookupEntry  *entries;
} BtkLookupTable;


extern BtkLookupTable *Btk_read_lookup_table(char * /*path*/);
extern void Btk_destroy_lookup_table(BtkLookupTable * /*table*/);
