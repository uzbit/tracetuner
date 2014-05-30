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
 * 2.2 2003/11/06 18:18:23
 */

/*
 *  Btk_align_peaks.h $Revision: 1.5 $
 */

#define NUM_COLORS     (4)     /* # of dyes */
#define QVMIN(a,b)  (a<b)?a:b
#define QVMAX(a,b)  (a>b)?a:b
#define INF    200000000       /* infinity */
#define NINF  -200000000       /* negative infinity */
#define MAX_NAME_LENGTH 200

extern int 
Btk_align_peaks(int **, int, int **, long *, int, 
    int, char *, char *, BtkMessage *);

