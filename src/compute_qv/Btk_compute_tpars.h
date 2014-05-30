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
 * 2.23 2003/11/06 18:18:27
 */

/*
 *  Btk_compute_tpars.h  $Revision: 1.10 $
 */

#if 0   /* The following files need to be included: */
#include "util.h"    // needed by train.h
#include "train.h"    // typedef Results
#endif

extern int is_resolved(ColorData *, int );
extern void colordata_release(ColorData *color_data);
extern void bases_release(TT_Bases *bases);
extern void trace_parameters_release(TraceParameters *tp);
extern void data_release(Data *data);
extern int  colordata_create(ColorData *, int, int, char *, BtkMessage *);
extern int bases_create(TT_Bases *, int, BtkMessage *);
extern int trace_parameters_create(TraceParameters *, int, BtkMessage *);
extern int data_create(Data *, int, int, char *, BtkMessage *);
extern int bases_populate(int *, char **, int, int **, Data *, Options *,
    BtkMessage *);
extern int colordata_populate(int, int **, char *, Data *, BtkMessage *);
extern int data_populate(int *, char **, int,
    int **, int, int **, char *, Data *, Options *, BtkMessage *);
extern void data_nelease(Data *);
extern int bc_reorder_called_bases_and_peaks(Data *, BtkMessage *);
extern void show_input_options(Options *);

extern int 
Btk_compute_tpars_Sanger(
    int *,		/* pointer to input len. of the array of called bases */
    char **,		/* pointer to input array of called bases */
    int **,		/* pointer to input array of the called peak locations*/
    int *,		/* input length of the chromatograms array */
    int **,		/* input arrays which store chromatographic data for
			 * each of dyes, assumed to be 4 colors 
                         */
    char *,		/* array of bases corresponding to the colors */
    int ,               /* number of trace parameters (= 4 for Sanger data) */
    double **par0,      /* pointer to output array of trace parameters */
    double **par1,      /* pointer to output array of trace parameters */
    double **par2,      /* pointer to output array of trace parameters */
    double **par3,      /* pointer to output array of trace parameters */
    double **iheight,   /* pointer to output array of intrinsic peak height */
    double **iheight2,  /* pointer to output array of secondary peak height */
    double **ave_iheight,
    ReadInfo *,         /* pointer to data used in trace renormalization 
                         * procedure 
                         */
    BtkLookupTable *,   /* pointer to a lookup table */
    ContextTable * ,    /* pointer to a context table */ 
    Options options,    /* structure including file_name, nocall, etc. */
    BtkMessage *,	/* error code and descriptive text */
    Results *           /* statistical results used by train (not ttuner) */
);

extern int
Btk_compute_tpars_454(
    sffHeader *h, 
    sffRead *r, 
    double *params[6],
    double *flow_per_base, 
    uint8_t *orig_qvs, 
    Options options,
    BtkMessage *message
);
