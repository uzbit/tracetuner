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
 * 1.12 2003/11/06 18:18:40
 */

#if 0   /* The following files need to be included: */
//#include "Stats.h"    // needed by train.h
//#include "Histo.h"    // needed by train.h
#endif
 
#define NUM_PEAK_UPSTREAM 3 
#define NUM_PEAK_DOWNSTREAM 3   
#define AREA_FACTOR1 (0.01*ONE_PLUS)
#define AREA_FACTOR2 (0.02*ONE_PLUS)

extern int is_dp(Peak *, int *, Data *); 
extern int Btk_process_peaks(Data *, Options *, BtkMessage *); 
extern int data_create_single_ordered_peak_list(Data *, BtkMessage *); 
extern int data_recall_bases(Data *, char *, char *, BtkMessage *); 
extern double get_average_width1(Data *, int, int, int, BtkMessage *); 
extern double get_average_width2(Data *, int, int, int, double *, BtkMessage *); 
extern double get_average_ratio(Data *, int, int, int, BtkMessage *); 
extern double get_peak_area(int *, int, int, BtkMessage *); 
extern double get_average_abi_spacing(int , Data *, BtkMessage *);
extern int get_peak_position(int *, int, int, double, BtkMessage *); 
extern int is_true_peak(Peak); 
extern int get_peak_max(int *, int, int, BtkMessage *);
extern int resolve_multiple_peaks(Data *, int, int, int, double *,  
    Options *, BtkMessage *); 
extern void data_release(Data *); 
extern int is_maximum(int n, int *data, int length);
extern int data_detect_peaks(Data *, Options *, BtkMessage *);
extern int data_expand_peaks(Data *, Options *, BtkMessage *);
extern int data_resolve_peaks(Data *, Options *, BtkMessage *);
extern double get_average_width(int, int, int, int, double *,
    Data *, BtkMessage *);
