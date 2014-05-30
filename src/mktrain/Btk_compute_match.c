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
 * $Id: Btk_compute_match.c,v 1.9 2009/01/13 20:44:52 gdenisov Exp $                  
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "Btk_qv.h"
#include "util.h"
#include "Btk_qv_data.h"
#include "Btk_match_data.h"
#include "Btk_sw.h"
#include "Btk_compute_match.h"

#define QUERY_LEN_MULTIPLIER 0.1 

/* This function saves a region in the 'best-regions' list if
 *     there are open spaces in the best-regions list
 *         OR
 *     the score of the current region is better than some existing score.
 *         In this case, the existing region is removed, and the new
 *         region is saved.
 * The function's synopsis is:
 *
 * result = save_region( best-regions, region )
 *
 * where
 *      best-regions  is the list of regions
 *      region        is the address of the region to be saved
 *
 *      result        is 0 on success, !0 if an error occurs.
 */
static int
save_region( Region* best, Region* reg )
{
    int i, min;

    min=0;	
    for(i=1; i<KEEP_BEST; i++) {
        if ( best[i].score < best[min].score ) {
            min=i;
        }
    }

    if ( reg->score > best[min].score ) {
        best[min].score = reg->score;
        best[min].beg = reg->beg;
        best[min].end = reg->end;
        best[min].cbeg = reg->cbeg;
        best[min].cend = reg->cend;
        best[min].diag = reg->diag;
    }
    return SUCCESS;
}

/* This function finds regions where the data sequence and the
 * query sequence are most similar (without gaps). Its synopsis is:
 *
 * result = find_regions( ap, best, lib_contig, query_contig, clearRange, 
 *                        message)
 *
 * where
 *      ap            is the address of an Align_params data structure,
 *                    created by set_alignment_parameters
 *      best          is the list of best regions (of size KEEP_BEST)
 *      lib_contig    is the address of the library contig, whose
 *                    lookup table should already be populated
 *      query_contig  is the address of the query contig
 *      clearRange    is the range (of the query) over which to do the
 *                    alignment
 *      message       is the address of a BtkMessage, where information
 *                    about an error will be put, if any
 *
 *      result        is 0 on success, !0 if an error occurs.
 *
 * Algorithm:
 *     This function only considers sequences without gaps, so similarity
 *     is defined in terms of matches and mismatches.
 *
 *     It is based on the FastA algorithm.  See 'Rapid Sequence Comparison ..'
 *     by William R. Pearson.  Meth Enzym v183, 1990 p 72 for more info.
 *
 *     The basic idea is this:
 *         For each contigous k-tuple in the query
 *              Calculate an index (see contig_make_fasta_lookup_table for
 *                      details of how to calculate such an index)
 *              Use the library's lookup table to find each occurrence
 *                      of the k-tuple in the library
 *              For each occurance of the k-tuple in the library
 *                  Note that there's a match in the diagonal defined by
 *                      lib_position - query_position
 *                  If there is no region along that diagonal, define
 *                      a new region with score ktup*match
 *                  If there already is a region along that diagonal,
 *                      figure out whether it's better to extend the existing
 *                      region, or save the existing region and start a new
 *                      region.
 */
static int
find_regions( Align_params *ap, Region *best, Contig *lib_contig,
              Contig* query_contig, Range clearRange, BtkMessage *message)
{
    Region *diagonal;
    int r, i, j, query_pos, lib_pos, curr_diag;
    int index;
    int lib_size, query_size, size;
    int dist;
    int tmp_score;

    lib_size = lib_contig->length-KTUP+1;
    query_size = query_contig->length-KTUP+1;

    size = query_size + lib_size;

    diagonal= (Region*) malloc( sizeof(Region)*size);
    MEM_ERROR(diagonal);

    for( i=0; i<KEEP_BEST; i++) {
	best[i].score = -1;
    }

    for( i=0; i<size; i++) {
	diagonal[i].beg=-1;
	diagonal[i].cbeg=-1;
	diagonal[i].diag=i;
    }

    index=0;

    /* For each contigous k-tuple in the query */
    for( i=0, query_pos=clearRange.begin;
         query_pos<clearRange.end+1;
         query_pos++, i++ )
    {
	/* Calculate an index (see contig_make_fasta_lookup_table for
         * details of how to calculate such an index)
         */
	index *= ALPHABET_SIZE;
	index = index % lib_contig->lut.length;
	index += base2int( query_contig->sequence[query_pos] );
	if ( i < KTUP-1 ) { continue; } /* Not enough bases so far. */

	/* Use the library's lookup table to find each occurrence
	 * of the k-tuple in the library.
	 * For each occurance of the k-tuple in the library
	 */
	for( j=0; j<lib_contig->lut.sizes[index]; j++)
	{
	    lib_pos = lib_contig->lut.data[index][j];
	    curr_diag = lib_pos - query_pos + query_size - 1;

	    /* Note that there's a match in the diagonal defined by
	     * lib_position - query_position
	     */
	    if ( diagonal[curr_diag].beg == -1 )
	    {
		/* If there is no region along that diagonal, define
		 * a new region with score ktup*match
		 */
		diagonal[curr_diag].beg=query_pos-KTUP+1;
		diagonal[curr_diag].end=query_pos;
		diagonal[curr_diag].cbeg=lib_pos-KTUP+1;
		diagonal[curr_diag].cend=lib_pos;
		diagonal[curr_diag].score=KTUP*ap->matrix[0][0];
	    }
	    else	 /* There is already a region. */
	    {
		/* If there already is a region along that diagonal,
		 * figure out whether it's better to extend the existing
		 * region, or save the existing region and start a new
		 * region.
		 */
		dist = query_pos-diagonal[curr_diag].end;
		if ( dist <= KTUP ) {
		    tmp_score = diagonal[curr_diag].score
				+ dist * ap->matrix[0][0];
		} else {
		    tmp_score = diagonal[curr_diag].score
				+ KTUP * ap->matrix[0][0]
				+ (dist - KTUP) * ap->matrix[0][1];
		}
		if ( tmp_score < diagonal[curr_diag].score ) {
		    /* Better to save the existing region and start a new
		     * region.
		     */
		    save_region( best, &diagonal[curr_diag] );

		    diagonal[curr_diag].beg=query_pos-KTUP+1;
		    diagonal[curr_diag].end=query_pos;
		    diagonal[curr_diag].cbeg=lib_pos-KTUP+1;
		    diagonal[curr_diag].cend=lib_pos;
		    diagonal[curr_diag].score=KTUP*ap->matrix[0][0];
		} else {
		    /* Better to extend the existing region. */
		    diagonal[curr_diag].cend=lib_pos;
		    diagonal[curr_diag].end=query_pos;
		    diagonal[curr_diag].score=tmp_score;
		}
	    }
			
			
	}
    }

    /* Save any regions we are currently working on. */
    for( i=0; i<size; i++) {
	if ( diagonal[i].beg != -1 ) {
	    save_region( best, &diagonal[i] );
	}
    }

    r=SUCCESS;
    goto cleanup;

    error:
	r=ERROR;

    cleanup:
	if ( diagonal!= NULL ) { free( diagonal); }
	
	return r;
}

#if 0

/* This function prints the best regions to the specified file stream.
 * Its synopsis is:
 *
 * fprint_regions( out, best )
 *
 * where
 *      out           is the specified FILE *
 *      best          is the list of best regions (of size KEEP_BEST)
 *
 */
static void
fprint_regions( FILE *fout, Region *best )
{
    int i;

    fprintf(fout, "     i  score   diag  begin    end cbegin   cend\n"); 
    for (i=0; i<KEEP_BEST; i++)
    {
	if ( best[i].score != -1 )
	{
	    fprintf(fout, "%6d %6d %6d %6d %6d %6d %6d\n", i, 
		best[i].score,
		best[i].diag,
		best[i].beg,
		best[i].end,
		best[i].cbeg,
		best[i].cend
	    );
	}
    }
    fprintf(fout, "\n");
}

static void
fprint_sorted_regions( FILE *fout, Region* ptrs[] )
{
    int i;

    fprintf(fout, "     i  score   diag  begin    end cbegin   cend\n"); 
    for (i=0; i<KEEP_BEST; i++)
    {
        Region* r = ptrs[i];
	if ( r->score != -1 )
	{
	    fprintf(fout, "%6d %6d %6d %6d %6d %6d %6d\n", i, 
		r->score,
		r->diag,
		r->beg,
		r->end,
		r->cbeg,
		r->cend
	    );
	}
    }
    fprintf(fout, "\n");
}

#endif


/* This function joins regions that overlap.  Its synopsis is:
 *
 * result = join_regions( best, newbest, length, liblength, message )
 *
 * where
 *      best          is the original list of best regions (of size KEEP_BEST)
 *      newbest       is the list of best regions after joining
 *      length        is the query length
 *      liblength     is the library length
 *      message       is the address of a BtkMessage, where information
 *                    about an error will be put, if any
 *
 *      result        is 0 on success, !0 if an error occurs.
 *
 * Algorithm:
 *     Calculates the beginning and ending position of the full fragment
 *         with the consensus.  For example, if the fragment is 500 bases
 *         long, and positions 25-100 of the fragment match with positions
 *         6025-6100 of the consensus, then the full fragment will match
 *         positions 6000-6500 of the consensus.
 *     Sorts by consensus position.
 *     Joins regions whose consensus start position differs by less than
 *         MAXDELTA.
 *
 * Assumes:
 *     The best regions that are defined occupy positions 0..x of the
 *     array.  (I.e., not 0, 4, 6, but rather a contiguous range. )
 */
static int
join_regions( Region *best, Region *newbest, int length, int liblength,
              BtkMessage *message )
{
    Region* ptrs[KEEP_BEST];
    Region* tmp;
    int i, j, max;
#if 0
    fprintf( stderr, "BEFORE: best:\n" );
    fprint_regions( stderr, best );
#endif
    for(max=0; max<KEEP_BEST; max++)
    {
        Region* r = &best[max];
        int bot, top;
	if ( r->score == -1 ) {
	    break;
	}
	ptrs[max] = r;
        bot = QVMIN( r->beg, r->cbeg );
        /* following "liblength-1" was "liblength" OLD: off by 1 (SLT) */
        top = QVMIN( length - 1 - r->end, liblength - 1 - r->cend );
        r->beg  -= bot;
        r->cbeg -= bot;
        r->end  += top;
        r->cend += top;
    }

    for(i=0; i<KEEP_BEST; i++) {
	newbest[i].score = -1;
    }

    if ( max == 0 ) {
	return SUCCESS;
    }

    /* Sort in order of increasing offset. */
    for(i=0; i<max; i++)
    {
	/* Invariant: list is sorted up to and including i-1. */
	for(j=i; j>0; j--)
	{
            if ( ptrs[j]->diag < ptrs[j-1]->diag)
	    {
		tmp=ptrs[j-1];
		ptrs[j-1]=ptrs[j];
		ptrs[j]=tmp;
	    }
	}
    }

    newbest[0] = *ptrs[0];

#if 1
    for(i=1, j=0; i<max; i++)
    {
        if( (ptrs[i]->diag - newbest[j].diag ) < MAXDELTA )
	{
	    /* Combine these two. */
            /* We don't use .cend in any subsequent processing (SLT)
             *  newbest[j].cend = ptrs[i]->cend; // original
             *  newbest[j].cend = -1000000;      // debug
             */
	    newbest[j].score += ptrs[i]->score;
	} else {
	    j++;
            newbest[j] = *ptrs[i];
	}
    }
#else
    for( i=0; i<max; i++ ) {
        newbest[i] = best[i];  
    }
#endif

#if 0
     fprintf( stderr, "AFTER: best:\n" );
     fprint_regions( stderr, best );
     fprint_sorted_regions( stderr, ptrs );
     fprintf( stderr, "newbest:\n" );
     fprint_regions( stderr, newbest );
#endif
	
    return SUCCESS;
}


/* This function rescores regions, possibly expanding or contracting
 * to maximize score.  Its synopsis is:
 *
 * result = rescore_regions( ap, best, newbest, library, query,
 * 			     message)
 *
 * where
 *      ap            is the address of an Align_params data structure,
 *                    created by set_alignment_parameters
 *      best          is the original list of best regions (of size KEEP_BEST)
 *      newbest       is the list of best regions after rescoring
 *      library       is the address of the library Contig
 *      query         is the address of the query Contig
 *      message       is the address of a BtkMessage, where information
 *                    about an error will be put, if any
 *
 *      result        is 0 on success, !0 if an error occurs.
 *
 * Algorithm:
 *     We use a simple dynamic programming algorithm:
 *        best_score(i) = MAX( best_score(i-1) + score(i), 0 )
 *        best_score(overall) = MAX(over i) { best_score(i) }
 */
static int
rescore_regions( Align_params *ap, Region *best, Region *newbest,
                 Contig *lib_contig, Contig* query_contig, BtkMessage *message)
{
    int lib_offset, query_start, i,j,k;
    int best_score, best_start, best_end;
	/* best_score is the score from best_start to best_end. */
    int curr_score, curr_start;
	/* curr_score is the score from curr_start to j. */

    for( k=0; k<KEEP_BEST; k++) {
	newbest[k].score = -1;
    }

    for( i=0, k=0; i< KEEP_BEST; i++) {
	if ( best[i].score > -1 )
	{
	    query_start = best[i].beg,
	    lib_offset = best[i].cbeg - best[i].beg;

	    best_score=0;
	    best_start=query_start;
	    best_end=query_start-1;
	    curr_score=0;
	    curr_start=query_start;
	    for( j=query_start; j<= best[i].end; j++)
	    {
	        /* best_score(i-1) + score(i) */
	        if ( query_contig->sequence[j] ==
		     lib_contig->sequence[j+lib_offset] ) {
		    curr_score += ap->matrix[0][0];
		} else {
		    curr_score += ap->matrix[0][1];
		}

	        /* best_score(i) = MAX( best_score(i-1) + score(i), 0 ) */
		if ( curr_score < 0 ) {
		    /* We'd get a better score by dumping everything
		     * and starting here.
		     */
		    curr_score=0;
		    curr_start=j+1;
		}

		/* best_score(overall) = MAX(over i) { best_score(i) } */
		if ( curr_score > best_score ) {
		    best_score=curr_score;
		    best_start=curr_start;
		    best_end=j;
		}
	    }

	    /* Now, extend the region if that helps. */

	    /* First, extend the end. */
	    curr_score = best_score;
	    for( j=best_end+1;
		 (j < query_contig->length) &&
		 ( j+lib_offset < lib_contig->length) &&
		 (curr_score > 0 );
		 j++)
	    {
		if ( query_contig->sequence[j] ==
		     lib_contig->sequence[j+lib_offset]) {
		    curr_score += ap->matrix[0][0];
		} else {
		    curr_score += ap->matrix[0][1];
		}
		if ( curr_score < 0 ) {
		    /* We're not going to get a better score which contains
		     * the original region.
		     */
		    break;
		}
		if ( curr_score > best_score ) {
		    best_score=curr_score;
		    best_end=j;
		}
	    }

	    /* Now, extend the beginning. */
	    curr_score = best_score;
	    for(j=best_start-1;
		     (j >= 0) && ( j+lib_offset >= 0 ) && (curr_score > 0 );
	        j--)
	    {
		if ( query_contig->sequence[j] ==
		     lib_contig->sequence[j+lib_offset]) {
		    curr_score += ap->matrix[0][0];
		} else {
		    curr_score += ap->matrix[0][1];
		}
		if ( curr_score < 0 ) {
		    /* We're not going to get a better score which contains
		     * the original region.
		     */
		    break;
		}
		if ( curr_score > best_score ) {
		    best_score=curr_score;
		    best_start=j;
		}
	    }

	    /* Finally, assign the best scores. */
	    newbest[k].beg = best_start;
	    newbest[k].end = best_end;
	    newbest[k].score = best_score;
	    newbest[k].diag = best[i].diag;
	    newbest[k].cbeg = best_start+lib_offset;
	    newbest[k].cend = best_end+lib_offset;
	    k++;
	}
    }

    return SUCCESS;
}


/* This function performs the first several steps of the FastA algorithm
 * on two sequences.  Its synopsis is:
 *
 * result = fasta_preprocess( ap, contig1, contig2, best, clearRange, message)
 *
 * where
 *      ap            is the address of an Align_params data structure,
 *                    created by set_alignment_parameters
 *      contig1       is the library contig, with its lookup table set
 *      contig2       is the query contig
 *      best          is the returned list of best regions (of size KEEP_BEST)
 *      clearRange    is the range of the query on which to perform the
 *                    preprocessing
 *      message       is the address of a BtkMessage, where information
 *                    about an error will be put, if any
 *
 *      result        is 0 on success, !0 if an error occurs.
 *
 * Algorithm:
 *     Calculate regions with highest density of matches (without gaps).
 *     Rescore those regions.
 *     Join close-by regions together.
 */
static int
fasta_preprocess( Align_params *ap, Contig* contig1, Contig* contig2,
                  Region *best, Range clearRange, BtkMessage *message)
{
    int r;
    Region tmpbest[KEEP_BEST];


    /* Find regions with highest density of perfect matches. */
    r=find_regions( ap, best, contig1, contig2, clearRange, message);
    if ( r==ERROR ) { return ERROR; }

/*
    fprintf(stderr, "Best regions:\n");
    fprint_regions( stderr, best, contig2->length );
*/

    /* Recalculate scores for those highest regions. */
    r=rescore_regions( ap, best, tmpbest, contig1, contig2, message);
    if ( r==ERROR ) { return ERROR; }

/*
    fprintf(stderr, "Best regions, rescored:\n");
    fprint_regions( stderr, tmpbest, contig2->length );
*/

    /* Join regions (with gaps). */
    r=join_regions( tmpbest, best, contig2->length, contig1->length,
                    message);
    if ( r==ERROR ) { return ERROR; }

/*
    fprintf(stderr, "Best regions after joining:\n");
    fprint_regions( stderr, best, contig2->length );
*/

    /* The actual FastA algorithm would do : */

    /* Create alignment using dynamic programming, for those
     * residues within band of size 32 around regions. */

    return SUCCESS;

}


/* This function returns the index of the highest scoring region.
 * Its synopsis is:
 *
 * highest = find_best_region(best, size)
 *
 * where
 *      best          is the list of best regions 
 *      size          is the length of the 'best' array
 *
 *      highest       is the index of best with the highest score
 *
 */
static int
find_best_region(Region *best, int size, int *read_direction)
{
    int i, max;
    
    max=0;
    for(i=0; i<size; i++) {
        if ( best[i].score > best[max].score )
        {
            if (i<size/2) *read_direction =  1;
            else          *read_direction = -1;
	    max=i;
	}
    }
    return max;
}


/* This function aligns the ends of the fragment with the vector
 * and extracts the clear range.  Its synopsis is:
 *
 * result = strip_vector( align_pars, query, vector, clear_range, start,
			  finish, message )
 *
 * where
 *      align_pars    is the address of an Align_params data structure,
 *                    created by set_alignment_parameters
 *      query         is the query contig
 *      vector        is the address of the vector, or NULL if there is none
 *      clear_range   is the address of a Range structure, which will be
 *		      set to the clear range - the unaligned range between
 *                    where the pre-Cut and post-Cut regions align
 *      start         is the address of an Align structure, which will be
 *                    filled if the pre-Cut part of the vector aligns
 *                    If the vector does not align, the score field will
 *                    be set to 0.
 *      finish        is the address of an Align structure, which will be
 *                    filled if the post-Cut part of the vector aligns
 *                    If the vector does not align, the score field will
 *                    be set to 0.
 *      message       is the address of a BtkMessage, where information
 *                    about an error will be put, if any
 *
 *      result        is 0 on success, !0 if an error occurs.
 *
 */
static int
strip_vector ( Align_params* align_pars, Contig* query, Vector* vector,
	       Range* clear_range, Align* start, Align* finish,
	       BtkMessage* message )
{
    int r, i;
    Range tmp;
    SubContig shortContig;

    clear_range->begin=0;
    /* clear_range->end=query->length; OLD: off by 1 (SLT) */
    clear_range->end=query->length -1;

    if (vector==NULL) {
	start->score = -1;
	finish->score = -1;
	return SUCCESS;
    }

    r=align_create(start, query->length, 0, 0, message);
    if ( r==ERROR ) { return ERROR; }
    r=align_create(finish, query->length, 0, 0, message);
    if ( r==ERROR ) { return ERROR; }

    /* Align with start. */
    start->num_gaps = QVMAX( query->length, vector->preCut.length );

    tmp.begin = 0;
    /*tmp.end = QVMIN( query->length, MAX_PRE_CUT_POS ); OLD: off by 1 (SLT) */
    tmp.end = QVMIN( query->length-1, MAX_PRE_CUT_POS );

    r=Btk_sw_alignment(align_pars, query, &(vector->preCut), NULL,
	start, tmp, 0, 1, message);
    if ( r==ERROR ) { return ERROR; }

    /* If the alignment is good, fix its coordinates, 
     * otherwise, set its score to 0.
     */
    if ( ( start->trace_len > 0 ) &&
         ( (vector->preCut.length - start->trace_dpos[ start->trace_len - 1 ] )
             < VECTORDELTA ) )
    {
        clear_range->begin = start->trace_qpos[ start->trace_len - 1 ] + 1;
	if ( vector->is_reverse ) {
	    for(i=0; i<start->trace_len; i++) {
	        start->trace_dpos[i] = vector->primerStart 
				     - start->trace_dpos[i];
	    }
	} else {
	    for(i=0; i<start->trace_len; i++) {
	        start->trace_dpos[i] += vector->primerStart;
	    }
	}
#if 0
        fprintf( stderr, "Vector beginning:\n");
        r=align_fprint( stderr, start, 70, message ); 
        if ( r==ERROR ) { return ERROR; }
#endif

    } else {
	start->score = 0;    /* Do not use. */
    }


    if ((vector->postCut).sequence == NULL) {
        /* Only short vector is specified by the user. */
        finish->score = -1;
        return SUCCESS;
    }


    /* Align with finish.  Because the postCut part may be very large,
     * we should be clever to speed things up.  
     */
    tmp.begin=0;
    tmp.end=100;
    r=contig_create_sub(&shortContig, &(vector->postCut), 0, 100, message);
    if ( r==ERROR ) { return ERROR; }

    finish->num_gaps = QVMAX( query->length, shortContig.length );

    r=Btk_sw_alignment(align_pars, query, &(shortContig), NULL,
	finish, tmp, 0, 1, message);
    if ( r==ERROR ) { return ERROR; }

    if ( finish->score <= 0 ) {
	/* No good matches. */
	return SUCCESS;
    }

    tmp.begin= QVMAX( (finish->trace_qpos[0]-10), MIN_POST_CUT_POS );
    /* tmp.end=query->length;    OLD: off by 1 (SLT) */
    tmp.end=query->length -1;
    finish->num_gaps = 0;    /* Force percentage of gaps, not whole number. */

    r=Btk_sw_alignment(align_pars, query, &(vector->postCut), NULL,
	finish, tmp, 0, 1, message);
    if ( r==ERROR ) { return ERROR; }

    /* If the alignment is good, fix its coordinates, 
     * otherwise, set its score to 0.
     */
    if ( ( finish->trace_len > 0 ) &&
         ( finish->trace_len > MIN_POST_CUT_LEN ) &&
         ( finish->trace_dpos[0]  < VECTORDELTA ) ) {
        clear_range->end = finish->trace_qpos[ 0 ] - 1;
	if ( vector->is_reverse ) {
	    for(i=0; i<finish->trace_len; i++) {
	        finish->trace_dpos[i] = vector->primerStart 
		  		        - vector->preCut.length
		  		        - finish->trace_dpos[i];
	    }
	} else {
	    for(i=0; i<finish->trace_len; i++) {
	        finish->trace_dpos[i] += vector->primerStart
		  		         + vector->preCut.length;
	    }
	}
#if 1
        fprintf( stderr, "Vector end:\n");
        r=align_fprint( stderr, finish, 70, message ); 
        if ( r==ERROR ) { return ERROR; }
#endif

    } else {
	finish->score = 0;    /* Do not use. */
    }

    return SUCCESS;
}

/* This function finds possible good alignments using the FastA
 * (heuristic) algorithm.  It then uses the Smith-Waterman (exact)
 * algorithm on the best several to get exact scores.  It returns
 * the number of high scoring alignments, and if there is only one,
 * the alignment itself.  Its synopsis is:
 *
 * result = Btk_compute_match( align_pars, align_parsIUB, lib_seq,
 *      lib_rev_comp,  query_seq, num_good, align_range, repeatFrac, best,
 *      vector, start, finish, clearRange, message)
 *
 * where
 *      align_pars    is the address of an Align_params data structure,
 *                    created by set_alignment_parameters, for use in 
 *                    aligning consensus with sample
 *      align_parsIUB is the address of an Align_params data structure,
 *                    created by set_alignment_parameters_IUB, for use in 
 *                    aligning vector with sample
 *      lib_seq       is the address of the library contig.
 *                    If its lookup table is
 *                    set, fine, otherwise this function sets it.
 *      lib_rev_comp  is the address of the reverse complement (library) contig.
 *                    If its lookup table is set, fine, otherwise this
 *                    function sets it.
 *      query_seq     is the query contig
 *      num_good      is the address of an int, which will be set to
 *                    the number of good alignments
 *      align_range   is the address of a Range structure, which will be
 *		      set to the sub-range of the query that aligns best
 *                    with the consensus
 *      repeatFrac    is a threshold.  Scores above this fraction of the
 *                    high score are assumed to be repeats.
 *      best          is the address of the return Align structure, which
 *                    will be set if there is a unique best score
 *      vector        is the address of the vector, or NULL if there is none
 *      start         is the address of an Align structure, for alignment
 *                    with vector before the consensus
 *                    start->score will be 0 if there is no good alignment
 *      finish        is the address of an Align structure, for alignment
 *                    with vector after the consensus
 *                    finish->score will be 0 if there is no good alignment
 *      clear_range   is the address of a Range structure, which will be
 *		      set to the sub-range of the query between the
 *		      start and end parts of the vector (if found)
 *      message       is the address of a BtkMessage, where information
 *                    about an error will be put, if any
 *
 *      result        is 0 on success, !0 if an error occurs.
 *
 * Algorithm:
 *     Use FastA to find regions with good alignments.
 *     Use Smith-Waterman on the best regions to find the exact scores.
 *     Discard repeats.
 *     If there is a unique best region, i.e., not repeats, return it.
 */
int
Btk_compute_match( Align_params* align_pars, Align_params* align_parsIUB,
    Contig* lib_seq, Contig* lib_rev_comp,  Contig* query_seq,
    int *num_good_alignments, Range *align_range, double repeatFraction,
    Align *best_alignment, Vector* vector, Align* start, Align* finish,
    Range* clearRange, int phd, BtkMessage* message)
{
    Region best[2*KEEP_BEST];
    int r, i, final;
    Align best_aligns[2*KEEP_BEST];
    int maxscore=-1;
    int read_direction;

    best_alignment->score=-1;

    align_range->begin=0;
    /* align_range->end = query_seq->length; OLD: off by 1 (SLT) */
    align_range->end = query_seq->length - 1;
#if 0
    fprintf( stderr, "lib_len= %d %d   q_len=%d\n",
             lib_seq->length, lib_rev_comp->length, query_seq->length );
#endif
    /* If the FastA lookup tables haven't been created yet, do so. */
    if( lib_seq != NULL && lib_seq->lut.length == 0 ) {
        r=contig_make_fasta_lookup_table( lib_seq, KTUP, message);
        if ( r==ERROR ) { 
            fprintf(stderr, "Error: Cannot make lookup_table1\n");
            return ERROR; 
        }
    }
    if( lib_rev_comp != NULL && lib_rev_comp->lut.length == 0 ) {
        r=contig_make_fasta_lookup_table( lib_rev_comp, KTUP, message);
        if ( r==ERROR ) { 
            fprintf(stderr, "Error: Cannot make lookup_table2\n");
            return ERROR; 
        }
    }

    for( i=0; i<2*KEEP_BEST; i++) {
        best_aligns[i].score = -1;
    }

    /* Strip vector. */
    r=strip_vector( align_parsIUB, query_seq, vector, clearRange, start,
		    finish, message );
    if ( r==ERROR) { 
        return ERROR; 
    }

    /* Fasta Forward. */
    if (lib_seq != NULL && lib_seq->sequence != NULL)
    {
        r=fasta_preprocess( align_pars, lib_seq, query_seq, best, *clearRange,
			message ); 
        if ( r==ERROR ) { return ERROR; }
    }

    /* FastA Backward. */
    if (lib_rev_comp != NULL && lib_rev_comp->sequence != NULL)
    {
        r=fasta_preprocess( align_pars, lib_rev_comp, query_seq, &best[KEEP_BEST],
                        *clearRange, message ); 
        if ( r==ERROR ) { 
            fprintf(stderr, "Error: Cannot process fasta file\n");
            return ERROR; 
        }
    }

    final = find_best_region(best, 2*KEEP_BEST, &read_direction);

    /* Do Smith-Waterman alignment for those regions whose score is
     * at least half of the highest score.
     */
    if ( best[final].score > 0 ) {
        for( i=0; i<2*KEEP_BEST; i++) {
            if ( best[i].score > best[final].score/2 ) 
            {
                r=align_create(&best_aligns[i], query_seq->length,
                    //best[i].cbeg==0 ? -best[i].beg : best[i].cbeg,
                      best[i].cbeg,
                      i/KEEP_BEST, message);
                if ( r==ERROR ) {
                    fprintf(stderr, "Error: Cannot create alignment\n");
                    return ERROR;
                }

                r=Btk_sw_alignment(align_pars, query_seq, lib_seq,
                                   lib_rev_comp, &best_aligns[i], *clearRange,
                                   read_direction, phd, message);
                if ( r==ERROR ) {
                    fprintf(stderr, "Error: Cannot produce SW alignment\n");
                    return ERROR;
                }

                if ( best_aligns[i].score > maxscore ) {
		    maxscore = best_aligns[i].score;
                }
            }
        }
    }

    *num_good_alignments=0;
    /* Print out best alignments, based on Smith-Waterman score. */
    if ( maxscore != -1 )
    {
        for( i=0; i<2*KEEP_BEST; i++) 
        {
	    /* Only print if they are within a factor of 2 of the absolute
             * best, AND the trace_len is at least 10% of the sequence length.
             */
            if ((best_aligns[i].score > (int) (repeatFraction*maxscore) ) && 
                (best_aligns[i].trace_len > (int) (
                    QUERY_LEN_MULTIPLIER*query_seq->length)  ) )
            {
		(*num_good_alignments)++;
		final=i;
/*
                fprintf( stderr, "i=%d, Score = %d, Dir=%c\n", i, 
                         best_aligns[i].score, (i/KEEP_BEST)?'B':'F' );
                r=align_fprint( stderr, &best_aligns[i], 70, message);
                if ( r==ERROR ) {
                    return ERROR;
                }
*/
            }
        }
    }

    if( *num_good_alignments == 1 )
    {
	/* There is a unique best match, so this is valid. */
        align_range->begin = best_aligns[final].trace_qpos[0];
        align_range->end = 
	    best_aligns[final].trace_qpos[ best_aligns[final].trace_len - 1];
        *best_alignment = best_aligns[final];
	best_aligns[final].score = -1;
    }
    else
        goto error;

    /* Clean up. */
    for( i=0; i<2*KEEP_BEST; i++) {
        if ( best_aligns[i].score > -1 ) {
            r=align_release( &best_aligns[i], message );
            if ( r==ERROR ) {
                fprintf(stderr, "Error: Cannot release alignment\n");
                return ERROR;
            }
        }
    }

    return SUCCESS;
error:
    /* Clean up. */
    for( i=0; i<2*KEEP_BEST; i++) {
        if ( best_aligns[i].score > -1 ) {
            r=align_release( &best_aligns[i], message );
            if ( r==ERROR ) {
                fprintf(stderr, "Error: Cannot release alignment\n");
                return ERROR;
            }
        }
    }
    return ERROR;
}

