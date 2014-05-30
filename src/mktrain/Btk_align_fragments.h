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
 * 2.4 2003/11/05 22:59:12
 */

/* Data structure for base (used only in the current file)
 */
typedef struct {
    char   base;
    int    pos;             /* base position */
    int    ind;             /* base index in a fragment */
} Base;

/* A GAlign structure contains information returned from an alignment:
 * Note: the indexs (aind, rind) start from 1 instead of 0
 */
typedef struct {
    double score;       // alignment score
    int   *aind;        // alignment vector containing a-index
    int   *rind;        // alignment vector containing r-index
    int   *apos;        // alignment vector containing a-position
    int   *rpos;        // alignment vector containing r-positions
    char  *achar;       // alignment vector containing a-characters
    char  *rchar;       // alignment vector containing r-characters
    char  *mchar;       // '|' for match and  ' ' for mismatch
    int    len;         // length of alignment vectors
    int    max_len;
    int    num_gaps;
} GAlign;

extern long 
g3(char);

extern long 
p3(char, char, int);

extern int
GAlign_init(GAlign*, BtkMessage*);

extern int
GAlign_release(GAlign*, BtkMessage*);

extern int
align_fragments(Base*, int, int*, Base*, int, long (*g)(char), 
		long (*p)(char, char, int), GAlign*, FILE*, BtkMessage*);
