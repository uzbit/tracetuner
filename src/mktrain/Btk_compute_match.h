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
 * $Id: Btk_compute_match.h,v 1.5 2008/11/27 13:23:33 gdenisov Exp $                  
 */

#ifndef BTK_COMPUTE_MATCH_H
#define BTK_COMPUTE_MATCH_H

#define KTUP		6	/* K-tuples used in indexing. */
#define KEEP_BEST	10	/* Keep at most this many possible regions. */
#define MAXDELTA	100	/* Combine regions that are less than this
				 *	far apart. 
				 */
#define VECTORDELTA     15      /* If the alignment ends more than this far
				 * away from the restriction site, ignore it.
				 */

/*
 * The following three are vector fudge factors:
 *
 * MAX_PRE_CUT_POS and MIN_POST_CUT_POS help us rule out illegitimate
 *	alignments that don't match at the positions we think they should.
 * MIN_POST_CUT_LEN helps get rid of very short post-cut alignments that
 *	may mess things up.
 * For legitimate vector alignments, the pre-cut alignment may be very short,
 * 	so there is no MIN_PRE_CUT_LEN.
 */
#define MAX_PRE_CUT_POS  50
    /*
     * If the pre-cut part of the vector aligns with the sample after this
     * position, ignore that alignment.
     * E.g. if the beginning of the vector aligns with position 350 (by random
     * chance), ignore.
     */
#define MIN_POST_CUT_POS 400
    /* The analogous number for post-cut alignment. */

#define MIN_POST_CUT_LEN	30
    /*
     * If the post-cut part of the vector aligns with the sample,
     * but the alignment is less than this long, ignore the alignment.
     * E.g. if the post-cut part of the vector aligns with 15 bases,
     * ignore.
     */

typedef struct {
    int diag; /* The diagonal = lib_pos - query_pos + query_size - 1 */
    int beg;  /* The beginning position of the match, in query coordinates. */
    int end;  /* The ending position of the match, in query coordinates. */
    int cbeg; /* The beginning position of the match,
	       * in consensus (library) coordinates. */
    int cend; /* The ending position of the match,
	       * in consensus (library) coordinates. */
    int score;
} Region;

extern int
Btk_compute_match( Align_params* align_pars, Align_params* align_parsIUB,
    Contig* lib_seq, Contig* lib_rev_comp,  Contig* query_seq,
    int *num_good_alignments, Range *align_range, double repeatFraction,
    Align* best_alignment, Vector* vector, Align* start, Align* finish,
    Range* clearRange, int read_direction, BtkMessage* message);


#endif
