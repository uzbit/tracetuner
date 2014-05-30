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
 * 1.11 2003/11/06 18:18:26
 */


/*
 *  Btk_compute_tp.h  $Revision: 1.7 $
 */

extern int Btk_compute_tp(Data *, char *, int, double *[], ReadInfo *,
    ContextTable *, Options *, BtkMessage *);
extern double get_peak_height_ratio(Data *, int, int, int, double,
    double,  Options *, BtkMessage *);
extern double get_peak_spacing_ratio(Data *, int, int, int, 
    Options *, BtkMessage *);
extern double get_peak_resolution(Data *, int, 
    Options *, BtkMessage *);
extern int get_trace_parameters_of_pure_bases(int, Data *, char *, 
    Options *, BtkMessage *);
extern int populate_params_array(int, Data* , double** );
extern double get_context_weight(char *context);
