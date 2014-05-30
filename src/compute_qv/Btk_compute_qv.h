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
 * 2.16 2003/11/06 18:18:26
 */

/*
 *  Btk_compute_qv.h  $Revision: 1.8 $
 */

extern int assign_quality(int, uint8_t *, double **, Options *, 
    BtkLookupTable *);
extern int Btk_compute_qv(
    int *,	       /* pointer to input length of the array of called bases*/
    char **,	       /* pointer to input array of called bases */
    int **,	       /* pointer to input array of the called peak locations */
    int *,	       /* input length of the chromatograms array */
    int **,	       /* input arrays which store chromatographic data for
			* each of dyes, assumed to be 4 colors 
                        */
    char *,	       /* array of bases corresponding to the colors */
    BtkLookupTable *,  /* pointer to a populated structure which represents
		        * a lookup table 
                        */
#if USE_CONTEXT_TABLE
    ContextTable *,    /* pointer to a populated structure which represents
                        * a context table
                        */
#endif
    uint8_t **,        /* pointer to output array of quality values */
    Options,           /* options structure: includes file_name, nocall, etc. */
    BtkMessage *,      /* error code and descriptive text */
    Results *
);
extern uint8_t get_quality_value(double, double, double, double, BtkLookupTable *);
