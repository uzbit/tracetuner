/**************************************************************************
 * This file is part of TraceTuner, the DNA sequencing quality value,
 * base calling and trace processing software.
 *
 * Copyright (c) 2003-2008 Gennady Denisov.  All rights reserved.
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

#define DEFAULT_PEAK_SPACING 12
#define POLYFIT_DEGREE 5
#define NUM_MULTICOMP_ITER 16

extern int get_peak_spacing(Data *, Options *, BtkMessage *);
extern int Btk_process_raw_data(int *, int **, char *, Data *, 
    Options, BtkMessage *);
extern void output_analyzed_data(char *, char *, int *, int *, int *, int *,
    int, int, int *, int, Data *);
extern void output_chromatogram(char *, char *, int *, int *, int *, int *,
    int, Data *);
extern int multicomponent(int **, int, Options *, BtkMessage *);
extern int prebaseline(int, int **, Options *, BtkMessage *);
extern double spacing_curve(int);
