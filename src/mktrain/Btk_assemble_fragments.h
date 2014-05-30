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
 * 2.4 2003/11/05 22:59:13
 */

#include "Btk_match_data.h"

#define MAX_ASSEMBLED_VECTOR_LEN 400 // normally the assemebled vector is
				     // below 100 bases long.  400 is just
				     // a big number to ensure enough space
				     // is allocated

/**
 * A WeightedBase structure wraps a base character with its match
 * weight and insertion-deletion weight.
 */
typedef struct {
  char base;
  int  weight_match;
  int  weight_insdel;
} WeightedBase;


/**
 * A WeightedBases structure contains a contigous WeightedBase sequence.
 */
typedef struct {
  WeightedBase* bases;
  int           length;		// size of the WeghtedBase sequence
  int           max_length;	// maximum size allowed
} WeightedBases;

extern int
weightedbases_create(WeightedBases*, Contig*, BtkMessage*);

extern int
weightedbases_release(WeightedBases*, BtkMessage*);

extern int
weightedbases_print(WeightedBases*, FILE*, BtkMessage*);

extern int
extract_bases_from_WB(WeightedBases*, Base*, int*, BtkMessage*);

extern int
extract_bases_from_contig(Contig*, Base*, BtkMessage*);

extern int
compute_weightedbases(WeightedBases*, GAlign*, int, FILE*, BtkMessage*);

extern int
assemble_fragments(WeightedBases*, Contig*, int, FILE*, BtkMessage*);
