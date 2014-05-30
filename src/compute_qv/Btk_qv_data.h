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
 * 2.44 2003/11/06 18:18:41
 */

#ifndef __BTK_QV_DATA__
#define __BTK_QV_DATA__

#define MAX_ERROR_TOLERABLE 2
#define NUM_COLORS     (4)     /* # of dyes */
#define WINDOW_3 3
#define WINDOW_7 7
#define QVMIN(a,b)  (((a)<(b))?(a):(b))
#define QVMAX(a,b)  (((a)>(b))?(a):(b))
#define QVMIN3(a,b,c)  (QVMIN((a),QVMIN((b),(c))))
#define QVMAX3(a,b,c)  (QVMAX((a),QVMAX((b),(c))))     
#define INF    200000000       /* infinity */
#define NINF  -200000000       /* negative infinity */
#define MAX_NUM_OF_PEAK 6000   /* default value for the number of peaks */
                               /* corresponding to a particular color data */
#define MAX_NUM_BASES   4000

#define WIDTH_FACTOR1 1.5
#define WIDTH_FACTOR2 0.75 
/*#define MIN_PEAK_RESOLUTION 0.3*/
/*#define RESOLUTION_FACTOR 0.00001*/
#define MERGE_PEAKS 0
#define MAX_NAME_LENGTH 256

extern double Erf(double);
extern double F(double);
extern double Phi(double);

/* data structure for peak */
typedef struct {
    double area;		/* peak area */
    char   base;		/* base (if any) which corresponds to the peak */
    int    base_index;		/* index of the called base in bases array */
    int    color_index;         /* 0, 1, 2, ..., NUM_COLORS-1 */
    int    cd_peak_ind;         /* index of peak in the colordata peak list */
    int    data_peak_ind;       /* index of peak in the data peak list */
    int    data_peak_ind2;      /* index of the second peak in the data peak list 
                                 * in the case of mixed base
                                 */
    int    ibeg;	        /* position of the left  inflection point of a peak */
    int    iend;		/* position of the right inflection point of a peak */
    int    beg;                 /* left outer boundary of a peak */
    int    end;                 /* right outer boundary of a peak */
    int    height;		/* peak apparent hight (=signal at the peak's position) */
    double iheight;             /* peak's intrinsic height (as computed from the model) */
    double wiheight;            /* weighed intrinsic peak height (=iheight multiplied by
                                 * the normalization and context weight */ 
    int    is_called;		/* 1 for yes, 0 for no */
    int    is_truncated;        /* is_truncated */
    int    pos;	        	/* apparent peak position on a chromatogram */
    float  ipos;                /* intrinsic peak position on a chromatogram */
    int    ipos_orig;		/* original intrinsic peak position */
				/* (before mobility shift correction) */
    int    max;                 /* position of peak's maximum or -1 */
    double relative_area;	/* ratio of peak area to the average area of
                                 * 10 preceeding peaks 
                                 */
    int    spacing;		/* distance between the current and previous peak's pos */
    int    type;                /* depends on the types of its boundaries; assumes
                                 * 9 possible values: 11,12,13,21,22,23,31,32 and 33 
                                 */
    double width1;              /* peak width at half-height */
    double width2;              /* ratio of peak area to its height */   
    double ave_width1;          /* average peak width1/width2 */
    double ave_width2;          /* average peak width2 */
    double ave_sq_width2;       /* average square of peak width2 */
    double orig_width;          /* original width of a peak as computed from
                                 * mathematical model */
    double beta;                /* square root of the product of time of peaks's 
                                 * evolution aand the DNA fragment diffusion coefficient,
                                 * as computed from mathematical model */
    double ave_w02beta;         /* averaged ratio of orig_width to beta */
    double C0;                  /* original height of a peak */
    double resolution;          /* integral of the difference between the original 
                                 * and intrinsic peak within the peak's bounds
                                 */
} Peak;

typedef struct {
    char  *bases;		/* array of called bases */
    Peak **called_peak_list;	/* array of pointers to called peaks;
                                 * Note: called_peak_list[i].pos, which we compute
				 * from color_data, may differ from coordinate[i]  
                                 */
    int   *coordinate;		/* list of called peak locations as obtained
                                 * directly from the sample file (or database)
                                 */
    int    length;		/* length of the array "bases" as read in from
                                 * the sample file 
                                 */
    int    max_length;          /* maximal (allocated) length of the array "bases" */
    int    orig_length;
} TT_Bases;

typedef struct {
    char  base;			/* the base, 'A', 'C', 'G' or 'T', is linked
				 * to the color_number depending on the sample
				 * file 
                                 */
    int  *data;			/* data points for a particular color_number */
    int   dye_number;		/* 1, 2, ..., NUM_COLORS  */
    int   length;		/* length of the data array */       
    int   max_value;		/* maximum value of signal for the trace */
    Peak *peak_list;		/* list of all peak corresponding to a
				 * particular color_number 
                                 */
    int   peak_list_len;	/* actual length of the array of all peaks;
				 * this is determined from the color_data */
    int   peak_list_max_len;	/* default peak list length, MAX_NUM_OF_PEAK, which
				 * is used when initializing the structure 
                                 */
} ColorData;

typedef struct {
    int     length;		/* current length trace parameters array */
    double *phr3;		/* peak height ratio in a window of 3 called peaks */
    double *phr7;		/* peak height ratio in a window of 7 called peaks */
    double *pres;		/* peak resolution */
    double *psr7;		/* peak distance (or spacing) ratio */
} TraceParameters;

typedef struct {
    TT_Bases      bases;
    ColorData  color_data[NUM_COLORS];	/* chromatograms */
    Peak     **peak_list;		/* array of pointers to all peaks */
    int        peak_list_len;           /* total number of peaks (of any color) */
    int        peak_list_max_len;       /* max allowed total number of peaks */
    int        length;                  /* max of the 4 colordata lengths */
    char       color2base[NUM_COLORS];
    int        pos_data_beg;            /* used when processing raw data */
    int        pos_data_end;            /* used when processing raw data */
    TraceParameters trace_parameters;
    char       chemistry[MAX_NAME_LENGTH];
} Data;

typedef struct {
    char  *chemistry;                       
    int    dev;
    int    edited_bases;      /* use edited, rather than called, bases */
    char   file_name[MAX_NAME_LENGTH]; /* sample file name        */
    int    gauss;             /* use gaussian peak shape model */
    int    het;               /* call heterozygotes */
    int    indel_detect;      /* detect indels, recall bases */
    int    indel_resolve;     /* det. indels, recall bases, separate chroamts */
    int    indsize;
    int    indloc;
    int    inp_phd;           /* read orig. bases & peak locs from .phd file */
    char   inp_phd_dir[MAX_NAME_LENGTH];  /* read input phd file(s) from this dir */
    int    ladder;
    int    lut_type;          /* lookup table type */
    float  min_ratio;         /* threshold peak height ratio */  
    int    mix;               /* whether to call mixed bases */
    int    multicomp;         /* Whether to multicomponent the (raw) data */
    int    nocall;            /* skip recalling bases */
    int    respace;           /* skip "fixing" of the multiple peak positions */
    double sf[4];             /* Scaling factors; used for debugging */
    char   scf_dir[MAX_NAME_LENGTH];   /* output directory name   */
    int    scf_version;       /* output scf file version (= 2 or 3) */
    char   path[MAX_NAME_LENGTH];      /* path to the sample file */
    int    poly;
    char   poly_dir[MAX_NAME_LENGTH];   /* output directory name   */
    int    process_bases;
    int    raw_data;          /* use raw colordata */
    int    recalln;           /* use TT 1.0's base recalling procedure */
    int    recallndb;         /* only recall Ns and dye blobs */
    int    renorm;            /* skip renormalization of traces */ 
    int    shift;             /* skip mobility shift correction */    
    char   tab_dir[MAX_NAME_LENGTH];   /* output directory name   */
    char   tal_dir[MAX_NAME_LENGTH];
    char   hpr_dir[MAX_NAME_LENGTH];
    int    time;              /* show the time of processing */
    char   tip_dir[MAX_NAME_LENGTH];   /* output directory name   */
    int    Verbose;           /* whether and how much status info to print */
    int    xgr;               /* output results in xgraph-readable format */
} Options;

/* Alternative base call */
typedef struct {
    char  base;             /* base character */
    int   pos;              /* alternative base position */
    int   qv;               /* alternative base quality value */
    int   base_index;       /* index of original base in the read */
    int   data_peak_ind;    /* indexes of peaks forming alternative */
    int   data_peak_ind2;   /* base call in data peak list  */
} AltBase;

typedef enum {
    ABI3700pop5=1,      /* ABI 3700 Pop-5 table */
    ABI3700pop6,        /* ABI 3700 Pop-6 table */
    ABI3100,            /* ABI 3100 Pop-6 table */
    ABI3730pop7,        /* ABI 3730 Pop-7 BDv3. table */
    MegaBACE            /* MegaBACE table */
} TABLETYPE;

extern void
exit_message(Options *, int);
#endif
