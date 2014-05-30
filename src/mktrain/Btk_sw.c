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
 * $Id: Btk_sw.c,v 1.9 2009/01/13 20:46:47 gdenisov Exp $                  
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "Btk_qv.h"
#include "util.h"
#include "Btk_qv_data.h"
#include "Btk_match_data.h"
#include "Btk_sw.h"

#define MAX_SIZE_OF_BASES  5000
#define  FRACTION_GAPS     0.2	/* Allow no more than 20% gaps. */
#define DESIRED_WIDTH 100

typedef struct{
    int width; /* half-band-width -- LZ */
    int maxq,maxd;
    char * path;
    int nmax1;
} passed_vars;


/* This function returns the maximum of four numbers and an indication
 * of which argument achieved the maximum value.  Its synopsis is:
 *
 * max = align_max_4( x0, x1, x2, x3, best )
 *
 * where
 *      x0..x3              are integers
 *      best                is the address of the outputted argument,
 *                          i.e., if x0 is the max, *best = 0, etc.
 *      max                 is the maximum of x0..x3
 *
 * Ex.:
 *     max = align_max_4( 10, 30, 20, 10, &c)
 * returns
 *     max = 30
 *     c = '\1'
 */
static int
align_max_4(int x0, int x1, int x2, int x3, char* best_dir, int dir, int phd)
{
    int max;

    if (dir > 0) /* shift gaps to the end of sequence */
    {
        max=x2;
        (*best_dir) = (char) 2;

        if ( x1 > max )
        {
            max=x1;
            (*best_dir) = (char) 1;
        }

        if ( x0 > max )
        {
            max=x0;
            (*best_dir) = (char) 0;
        }
    }
    else 
    {
        if (phd)
        {
            max=x0;
            (*best_dir) = (char) 0;

            if ( x1 > max )
            {
                max=x1;
                (*best_dir) = (char) 1;
            }

            if ( x2 > max )
            {
                max=x2;
                (*best_dir) = (char) 2;
            }
        }
        else
        {
            max=x1;
            (*best_dir) = (char) 1;

            if ( x2 > max )
            {
                max=x2;
                (*best_dir) = (char) 2;
            }

            if ( x0 > max )
            {
                max=x0;
                (*best_dir) = (char) 0;
            }
        }
    }

    if ( x3 > max )
    {
        max=x3;
        (*best_dir) = (char) 3;
    }

    return max;
}

/******************************************************************************
 * Function: are_similar_bases
 * Purpose:   return 1 or 0 depending on whether bases are "similar" or not.
 ******************************************************************************
 */
static int
are_similar_bases(char qbase, char dbase)
{
    if (qbase == dbase)
        return 1;

    else if (dbase == 'M') {
        if ((qbase == 'A') || (qbase == 'C') || (qbase == 'M'))
            return 1;
    }
    else if (dbase == 'R') {
        if ((qbase == 'A') || (qbase == 'G') || (qbase == 'R'))
            return 1;
    }
    else if (dbase == 'K') {
        if ((qbase == 'G') || (qbase == 'T') || (qbase == 'K'))
            return 1;
    }
    else if (dbase == 'S') {
        if ((qbase == 'G') || (qbase == 'C') || (qbase == 'S'))
            return 1;
    }
    else if (dbase == 'W') {
        if ((qbase == 'A') || (qbase == 'T') || (qbase == 'W'))
            return 1;
    }    
    else if (dbase == 'Y') {
        if ((qbase == 'T') || (qbase == 'C') || (qbase == 'Y'))
            return 1;
    }
    return 0;
}

/* This function performs the Smith-Waterman dynamic programming
 * algorithm on the two sequences. Its synopsis is:
 *
 * result = align_pair_1( passed, align, params, query, data, message )
 *
 * where
 *      passed        is the address of a passed_vars structure, which
 *                    contains information used by both this function
 *                    and align_pair_2
 *      align         is the address of the Align structure.  The only
 *                    change to it will be in the score field.
 *      params        is the address of an Align_params structure,
 *                    which contains a scoring matrix and ins/del penalties
 *      query         is the address of a Contig containing the query
 *                    sequence
 *      data          is the address of a Contig containing the data
 *                    sequence
 *      message       is the address of a BtkMessage, where information
 *                    about an error will be put, if any
 *
 *      result        is 0 on success, !0 if an error occurs.
 *
 * The algorithm:
 * Creates the dynamic programming matrix using a band of size passed->width:
 *     ---------> m (data)
 *  |  xxxx
 *  |   xxxx
 *  |    xxxx
 *  \/    xxxx
 *  n      xxxx
 * (query)
 *
 * (In the diagram shown above, the passed->width = 2. -- LZ)
 *
 *     Imagine an extra row and an extra column for empty matches,
 *     (i.e., gaps on the far left end):
 *         ---g               or              actg
 *         actg                               ---g
 *     these would be row -1 and column -1.
 * 
 * Fills in by columns.
 *
 * Stores the path from the current position back to the previous best.
 * The indexing is: path[m, n] =>  path[ (nmax+1)*(m+1) + (n+1) ]
 *     where the matrix 'path' has been compressed into an array, and
 *     the '+1's are due to the two imaginary rows.
 *
 * Note: the function aligns the whole sequences, and returns the
 *     best local sub-alignment.  To align sub-sequences, create
 *     SubContigs and then call the function.
 */
static int
align_pair_1(passed_vars *pv, Align* align, Align_params* align_pars,
    Contig* query, Contig* data, int dir, int phd, BtkMessage* message)
{
    int np;
    int mmax, nmax;

    /* High scores, for current column and previous column. */
    int *sh_current=NULL, *sh_prev=NULL, *tmp=NULL;
    int nstart,nend;
    int m=0,n=0;
    int i;
 
    int cur_cell, up_cell, left_cell, upleft_cell;
 
    int iq,id;

    int r;

    align->score=0;
 
    pv->nmax1=0;
    pv->maxq=0;
    pv->maxd=0;

    mmax = QVMIN( (query->length+pv->width), data->length );
    nmax = QVMIN( (mmax+pv->width) , query->length+1 );
 
    pv->nmax1 = nmax+1;
 
    if(nmax < 0 || mmax < 0)  {
        message->code=1;
        return ERROR;
    }

    sh_prev = (int*) malloc(sizeof(int)*(nmax+1));
    MEM_ERROR(sh_prev);
    sh_current = (int*) malloc(sizeof(int)*(nmax+1));
    MEM_ERROR(sh_current);
    pv->path = (char*) malloc(sizeof(char)*(nmax+1)*(mmax+1));
    MEM_ERROR(pv->path);

    for(i=0; i<(nmax+1)*(mmax+1); i++)
	pv->path[i]=25;
 
    /* Enforce zero boundary condition */
    for(i=0;i<(nmax+1);i++){
        sh_prev[i] = 0;
    }
 
 
    for(m=0;m<mmax;m++){
        sh_current[0] = 0;
        nstart = QVMAX( 0, (m - pv->width + 1) );
        nend = QVMIN( (m+pv->width-1) , query->length );

        for(n=nstart;n<nend;n++){
 
            np = n-nstart;

	    /* Calculate index of current cell and its neighbors.
	     * The only trickiness is if we've just started, and
	     * the size of sh_prev is less than 2*width.
	     */
            cur_cell = np + 1;
            up_cell = cur_cell-1;
            if(nstart == 0) {
                left_cell = np + 1;
            } else {
                left_cell = np + 2;
            }
            upleft_cell = left_cell - 1;

            iq = (int) query->sequence[n];
            id = (int) data->sequence[m];

    	    sh_current[cur_cell] = align_max_4 (
		sh_prev[upleft_cell]+align_pars->matrix[iq][id],
		sh_current[up_cell]+align_pars->gap_init,
		sh_prev[left_cell]+align_pars->gap_ext,
		0,
		&pv->path[ (m+1)*(nmax+1) + (np+1) ],
                dir, phd
	    );

#ifdef ALIGN_GLOBAL
            if((n==query->length-1) || (m == mmax-1)) {
#endif
               if(align->score < sh_current[cur_cell]){
                    align->score = sh_current[cur_cell];
                    pv->maxq = n;
                    pv->maxd = m;
                }
#ifdef ALIGN_GLOBAL
            }
#endif
      } /* end inner loop */
	sh_current[nend-nstart+1] = NINF; /* So that left_cell will be defined
		   	                   * in the next iteration.
                                           */
	tmp=sh_prev;
	sh_prev=sh_current;
        sh_current=tmp;
    } /* end outer loop */


    /* End of SW. */
 
    r=SUCCESS;
    goto cleanup;

    error: 
        r = ERROR;
        FREE(pv->path);

    cleanup:
        FREE(sh_current);
        FREE(sh_prev);

    return r;
}

/* This function uses the results of the Smith-Waterman dynamic programming
 * algorithm to generate the alignment.  Its synopsis is:
 *
 * result = align_pair_2( passed, align, query, data, message )
 *
 * where
 *      passed        is the address of a passed_vars structure, which
 *                    contains information used by both this function
 *                    and align_pair_2
 *      align         is the address of the Align structure
 *      query         is the address of a Contig containing the query
 *                    sequence
 *      data          is the address of a Contig containing the data
 *                    sequence
 *      message       is the address of a BtkMessage, where information
 *                    about an error will be put, if any
 *
 *      result        is 0 on success, !0 if an error occurs.
 *
 * The algorithm:
 * Starts at n = passed->maxq, m = passed->maxd.
 *
 * Uses the stored path.
 * The indexing is: path[m, n] =>  path[ (nmax+1)*(m+1) + (n+1) ]
 *     where the matrix 'path' has been compressed into an array, and
 *     the '+1's are due to the two imaginary rows. (See align_pair_1.)
 *
 * Traces its way back to the beginning of the best local alignment.
 *
 * Sets the qpos and dpos arrays.  (Positions are with respect to the
 *     data and query Contigs, so if those are actually SubContigs or
 *     reverse complements, the calling function will have to adjust those.)
 * Sets the trace_dir array:
 *     Match     => '|'
 *     Mismatch  => ' '
 *     Deletion  => '1'
 *     Insertion => '2'
 * Sets the qchar, mchar, and dchar arrays.
 */
static int
align_pair_2(passed_vars *pv, Align* align, Contig* query,
             Contig* data, BtkMessage* message)
{
    int np;
    char c1,c2,c3;
    int k;

    int state;

    int countSame=0;
    int countMismatch=0;
    int countDelGap=0;
    int countInsGap=0;

    int mxnp;

    int m=0,n=0;
    int i;


    if(align->trace_max_len < 2*query->length+pv->width){
	/* trace_max_len too small */
        align->trace_max_len *= 2;
        align->trace_qpos = (int*) realloc(align->trace_qpos,
                                           sizeof(int)*align->trace_max_len);
        MEM_ERROR(align->trace_qpos);

        align->trace_dpos = (int*) realloc(align->trace_dpos,
                                           sizeof(int)*align->trace_max_len);
        MEM_ERROR(align->trace_dpos);

        align->trace_dir = (char*) realloc(align->trace_dir,
                                           align->trace_max_len);
        MEM_ERROR(align->trace_dir);

        align->trace_qchar = (char*) realloc(align->trace_qchar,
                                             align->trace_max_len);
        MEM_ERROR(align->trace_qchar);

        align->trace_mchar = (char*) realloc(align->trace_mchar,
                                             align->trace_max_len);
        MEM_ERROR(align->trace_mchar);

        align->trace_dchar = (char*) realloc(align->trace_dchar,
                                             align->trace_max_len);
        MEM_ERROR(align->trace_dchar);
    }

    n = pv->maxq;
    m = pv->maxd;

    k = 0;

    /* Loop through and figure out size of alignment region. */

    while(m >= 0 && n >= 0){
        np = n - (QVMAX( 0, (m - pv->width + 1) ));
        mxnp = (m+1)*pv->nmax1 + (np+1);
        state = pv->path[mxnp];

        switch (state) {
        case 0 :	/* Match or Mismatch. */
            n--;
            m--;
            break;
        case 1 :	/* Deletion. */
            n--;
            break;
        case 2 :	/* Insertion. */
            m--;
            break;
        case 3 :	/* End (for local alignment.) */
            m=-1;
            n=-1;
            k--;	/* Do not include the current position in the
                         * count 'k'. */
            break;
        default:
            sprintf(message->text,
                "Alignment backtrace accessing illegal array reference\n");
            goto error;
        }

        /* check if out of bounds */
        if(n == (m-pv->width) || n == (m+pv->width))
        {
            sprintf(message->text, "Alignment backtrace hit boundary\n");
            goto error;
        }

        k++;
    } /* end of while loop to trace back */

    align->trace_len = k;

    n = pv->maxq;
    m = pv->maxd;
    k--;

    /* Loop through again and set all array values. */


    while( k >= 0 ){
        np = n - (QVMAX( 0, (m - pv->width + 1) ));
        mxnp = (m+1)*pv->nmax1 + (np+1);
        state = pv->path[mxnp];

        align->trace_qpos[k] = n;
        align->trace_dpos[k] = m;

        switch (state) {
        case 0 :	/* Match or Mismatch. */
            if(are_similar_bases(query->sequence[n],
                                  data->sequence[m])){
                align->trace_dir[k] = '|';
                countSame++;
            } else {
                align->trace_dir[k] = ' ';
                countMismatch++;
            }
            n--;
            m--;
            break;
        case 1 :	/* Deletion. */
            align->trace_dir[k] = 1+ '0';
            n--;
            countDelGap++;
            break;
        case 2 :	/* Insertion. */
            align->trace_dir[k] = 2+'0';
            m--;
            countInsGap++;
            break;
        case 3 :	/* End (for local alignment.) */
            m=-1;
            n=-1;
            break;
        default:
            sprintf(message->text,
                "Alignment backtrace accessing illegal array reference\n");
            goto error;
        }


        /* check if out of bounds */
        if(n == (m-pv->width) || n == (m+pv->width))
        {
            sprintf(message->text, "Alignment backtrace hit boundary\n");
            goto error;
        }

        k--;
    } /* end of while loop to set array values */

    for(i=0;i<align->trace_len;i++){
        c1 = query->sequence[align->trace_qpos[i]];
        c2 = data->sequence[align->trace_dpos[i]];
        c3 = align->trace_dir[i];
        if(c3 == '1'){
            c2 = '-';
            c3 = ' ';
        }
        if(c3 == '2'){
            c1 = '-';
            c3 = ' ';
        }
        align->trace_qchar[i] = c1;
        align->trace_mchar[i] = c3;
        align->trace_dchar[i] = c2;
    }

    return SUCCESS;

    error: 
    fprintf( stderr, "%s\n", message->text );
        return ERROR;
}

/* This function performs the Smith-Waterman dynamic programming
 * algorithm on two sequences and sets the appropriate arrays in
 * an Align data structure. Its synopsis is:
 *
 * result = align_pair( align, params, query, data, message )
 *
 * where
 *      align         is the address of the Align structure
 *      params        is the address of an Align_params structure,
 *                    which contains a scoring matrix and ins/del penalties
 *      query         is the address of a Contig containing the query
 *                    sequence
 *      data          is the address of a Contig containing the data
 *                    sequence
 *      message       is the address of a BtkMessage, where information
 *                    about an error will be put, if any
 *
 *      result        is 0 on success, !0 if an error occurs.
 *
 * Note: the function aligns the whole sequences, and returns the
 *     best local sub-alignment.  To align sub-sequences, create
 *     SubContigs and then call the function.
 * 
 * If align->num_gaps > 0, use that number of gaps,
 *    align->num_gaps == 0, use FRACTION_GAPS*min_length,
 *    align->num_gaps < 0, no preprocssing was used before Btk_sw_alignment()
 *    was called, use max(data->length, query->length) to cover the full
 *    matrix.
 * 
 * Therefore if Btk_sw_alignment() needs to be called without pre-processing
 * (such as FastA preprocessing), one can assign a negative value to
 * align->num_gaps before calling Btk_sw_alignment() to flag the sw algorithm
 * implemented here to compute the entire matrix.  In this case the pv.width
 * is calculated as the max of data->length and query->length.
 */
static int
align_pair(Align* align, Align_params* align_pars, Contig* query, 
    Contig* data, int dir, int phd, BtkMessage* message)
{   
    passed_vars pv;
    int min_length;
    int r;

    min_length = QVMIN( query->length, data->length );
    /* Negative align->num_gaps value means that no pre-processing was used
     * to find best regions before Btk_sw_alignment() was called.  So in
     * order to cover the entire matrix, the pv.width needs to be the max
     * of query length and data length. -- LZ
     */
    if ( align->num_gaps < 0 ) {
    	pv.width = QVMAX(query->length, data->length);
    } else {
    	if ( align->num_gaps == 0 ) {
            pv.width = (int)(min_length*FRACTION_GAPS);
    	} else {
	    pv.width = 2*align->num_gaps;
    	}

    	/* added by SLT (May 28 2001) */
    	if ( pv.width < DESIRED_WIDTH ) {
            pv.width = QVMIN( min_length, DESIRED_WIDTH );
    	}

        /* The following is somewhat of a kludge to deal with the fact
         * that calling with negative offsets apparently cause purify to choke.
     	 * Because we can't use negative offsets to specify exactly where
     	 * to start the search, we just open up the search window.
     	 * I realize this negates all the code above, but eventually
     	 * we may want to fix this the right way.
     	 */
    	pv.width = min_length;   /* kludge added by SLT (May 29 2001) */
    }

    r = align_pair_1(&pv, align, align_pars, query, data, dir, phd, message);
    if (r == ERROR) {
	return ERROR;
    }

    align_pair_2(&pv, align, query, data, message);

    FREE(pv.path);

    return SUCCESS;
}


/* This function performs the Smith-Waterman dynamic programming
 * algorithm on two sequences and sets the appropriate arrays in
 * an Align data structure. Its synopsis is:
 *
 * result = Btk_sw_alignment( params, query, library, rev_comp,
 *                            align, range, message )
 *
 * where
 *      params        is the address of an Align_params structure,
 *                    which contains a scoring matrix and ins/del penalties
 *      query         is the address of a Contig containing the query
 *                    sequence
 *      library       is the address of a Contig containing the library
 *                    sequence
 *      rev_comp      is the address of a Contig containing the reverse
 *                    complement of the library sequence
 *      align         is the address of the Align structure.
 *                    Its contig_offset and base_is_reverse fields 
 *                    should be set when the function is called.
 *                    At the end, all of the relevant fields in align
 *                    will be set.
 *		      NOTE: a negative value of align->num_gaps means
 *                          the entire matrix for aligning sub-sequences
 *                          will be computed.  This is useful when one
 *                          needs to call this Btk_sw_alignment() function
 *                          directly (skipping any pre-processing procedures
 *                          that probe for best regions to run sw), eg.
 *                          when aligning two short sequences.
 *      range         is the range of query that should be aligned.
 *      message       is the address of a BtkMessage, where information
 *                    about an error will be put, if any
 *
 *      result        is 0 on success, !0 if an error occurs.
 *
 * Note: the function aligns sub-sequences, and returns the
 *     best local sub-alignment.  It sets trace_qpos and
 *     trace_dpos with respect to the original sequences, not the
 *     sub-sequences.  (Previous versions of training code set 
 *     trace_dpos with respect either the forward or reverse sequence,
 *     this version always sets it with respect to forward.)
 */
int 
Btk_sw_alignment(Align_params *ap, Contig *query, Contig *library,
    Contig *rev_comp_library, Align *align,
    Range align_range, int dir, int phd, BtkMessage* message)
{
    int    r=0,i;
    int    align_range_len, align_range_begin;
    SubContig lib_ranges, query_ranges;

    /* Get the "align" part of the fragment sequence 
     * or its reverse compliment
     */
    /* I'm not sure (SLT), but the following test is probably only needed 
     * because of bugs in 
     * "strip_vector()", "join_regions()", "Btk_compute_match()"
     * which did things like:
     *    align_range.end = query->length
     * instead of:
     *    align_range.end = query->length-1
     * These bugs are now fixed, but for paranoia's sake, I'll keep the
     * following debug message in for now.
     */
#if 1 /* NEW: a (hopefully unneeded) debug test */
    {
        int len1 = align_range.end - align_range.begin + 1,
            len2 = query->length   - align_range.begin;
        align_range_len = len1;
        if( len2 < len1 ) {
            fprintf( stderr, 
                     "We really did need this test in Btk_sw_alignment\n" );
            fprintf(  stderr, "end=%d beg=%d len=%d\n",
                      align_range.end, align_range.begin,query->length );
            align_range_len = len2;
        }
    }
#else /* OLD: off by 1 (SLT) */
    align_range_len = QVMIN( align_range.end - align_range.begin + 1,
    		     query->length   - align_range.begin      );
#endif
    align_range_begin = align_range.begin;

    /* Create sub-contig for query. */
    r=contig_create_sub( &query_ranges, query, align_range_begin,
                         align_range_len, message );
    if(r==ERROR) {
        return ERROR;
    }


    /* Create sub-contig for library or its reverse complement, 
     * whichever is appropriate.
     */
    if(align->base_is_reverse) {
        r=contig_create_sub( &lib_ranges, rev_comp_library, 
                             align->contig_offset,
                             rev_comp_library->length - align->contig_offset,
                             message );
    } else {
        r=contig_create_sub( &lib_ranges, library, 
                             align->contig_offset,
                             library->length - align->contig_offset,
                             message );
    }
    if(r==ERROR) {
        return ERROR;
    }

    /* Align sub-range of library with sub-range of fragment sequence. */
    r = align_pair( align, ap, &query_ranges, &lib_ranges, dir, phd, message );

    if(r==ERROR) {
        goto error;
    }

    /* Change trace_dpos and trace_qpos to specify global position,
     * i.e., position w.r.t. 0, not w.r.t. alignment range.
     */
    for(i=0; i<align->trace_len; i++) {
        align->trace_qpos[i] += align_range_begin;
        align->trace_dpos[i] += align->contig_offset;
        if(align->base_is_reverse) {
            align->trace_dpos[i] = library->length - align->trace_dpos[i] - 1;
        }
    }

    return SUCCESS;

    error:
        return ERROR;
}
