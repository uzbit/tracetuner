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
 * 1.31 2003/11/06 18:18:25
 */

/*
 *  Btk_call_bases.h $Revision: 1.9 $
 */

#define QV_PI 3.14159265358979323846
#define DYE_BLOB_FRACTION 0.05
#define MIN_CALLED_PEAK_HEIGHT 5.

extern int bc_data_create_single_ordered_peak_list(Data *, int *, BtkMessage *);
extern int colordata_find_peak_index_by_location(ColorData *, int,
    Peak *, int *, BtkMessage *);
extern int uncall_peak(int, Data *, BtkMessage *);
extern int can_insert_base(Data *, int, int, int,
    double, double, Options *);
extern int is_dye_blob(int pos, Peak *, Data *, int); 
extern int mk_dip_list(Peak **, int, Peak ***, int *, 
    Options *, BtkMessage *, Data *);
extern int Btk_call_bases(Data *, char *, ReadInfo *, 
    ContextTable *, Options *, BtkMessage *, Results * );
extern int mark_called_peaks(Data *, char *, Options, BtkMessage *);
extern int create_DIP_list(Data *, int, int, int *, Peak ***, int *, 
    Options *, BtkMessage *);
extern void preset_base_calls(Data *, char *, Options, BtkMessage *);
extern void check_base_index(Data *data);
extern Peak initialize_peak();
