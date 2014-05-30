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
 * 2.8 2003/11/05 22:59:12
 */


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include "util.h"
#include "Btk_qv.h"
#include "Btk_align_fragments.h"
#include "Btk_assemble_fragments.h"

/* This function creates a weightedbase array from a contig.
 * The wights (weight_match and weight_insdel) are assigned to value 1.
 *
 * Its synopsis is:
 * result = weightedbases_create(weightedbases, contig, message)
 * where
 *      weightedbases    is the address of the weighted base array to be 
 *                       created
 *      contig           is the address of the source contig
 *      message          is the address of the Btk message object, where an
 *                       error information will be stored, if any. 
 *	result		is 0 on success, !0 if an error occurs.
 */
int
weightedbases_create(WeightedBases* wbases, Contig* c, BtkMessage* message)
{
    int i;

    if (wbases == NULL) {
        sprintf(message->text, 
		"first argument of weightedbases_create() is NULL");
        return ERROR;
    }

    if (c == NULL) {
        sprintf(message->text, 
		"second argument of weightedbases_create() is NULL");
        return ERROR;
    }

    /* the length of the wbases is equal to the length of the contig.
     * Note if contig is longer than the allowed maximum length of the
     * WeightedBases, the rest of the contig will be truncated.
     */
    wbases->length = c->length;
    if (wbases->length > wbases->max_length) {
        wbases->length = wbases->max_length;
    }

    /* The base of the weighted base will be the base of the contig.
     * The weights are initilaized to value of 1.
     */
    for (i=0; i<wbases->length; i++) {
        (wbases->bases[i]).base = c->sequence[i];
	(wbases->bases[i]).weight_match = 1;
	(wbases->bases[i]).weight_insdel = 1;
    }

    return SUCCESS;
}


/* This function releases the resource held by the specified weighted base
 * array.  Its synopsis is:
 * 
 * result = weightedbases_release(weightedbases, message)
 * where
 *      weightedbases    is the address of the weighted base array of which
 *                       the resource is to be released
 *      message          is the address of the btk message object, where an
 *                       error information will be stored, if any.
 *	result		is 0 on success, !0 if an error occurs.
 */
int
weightedbases_release(WeightedBases* wbases, BtkMessage* message) {
    if (wbases == NULL)
        return SUCCESS;

    free(wbases->bases);
    wbases->bases = NULL;
    wbases->length = 0;
    wbases->max_length = 0;

    return SUCCESS;
}


/**
 * This function prints out the base characters, weights of the specified
 * weightedbase array to the specified file.
 *
 * Its synopsis is:
 *
 * result = weightedbases_print(weightedbases, fout)
 * where
 *      weightedbases    is the address of the weighted base array to
 *                       be printed
 *      fout             is the specified FILE *
 *      message          is the address of the btk message object, where an
 *                       error information will be stored, if any.
 *	result		is 0 on success, !0 if an error occurs.
 */
int
weightedbases_print(WeightedBases* wbases, FILE* fout, BtkMessage* message) {
    int i, j, line_num, index;
    int bases_per_line = 10;  // print 10 bases per line

    if (wbases == NULL) {
        sprintf(message->text, 
		"first argument of weightedbases_print() is NULL");
        return ERROR;
    }

    if (fout == NULL) {
        sprintf(message->text, 
		"second argument of weightedbases_print() is NULL");
        return ERROR;
    }

    line_num = (int) ceil((double) (wbases->length) / bases_per_line);

    for (j = 0; j < line_num; j++) {

        fprintf(fout, " ind:\t");
	for (i = 0; i < bases_per_line; i++) {
	    index = j * bases_per_line + i;
	    if (index < wbases->length) {
	        fprintf(fout, "%5d  ", index + 1);
	    }
	}

	fprintf(fout, "\nbase:\t");
	for (i = 0; i < bases_per_line; i++) {
	    index = j * bases_per_line + i;
	    if (index < wbases->length) {
	        fprintf(fout, "    %c  ", (wbases->bases[index]).base);
	    }
	}

	fprintf(fout, "\n w_m:\t");
	for (i = 0; i < bases_per_line; i++) {
	    index = j * bases_per_line + i;
	    if (index < wbases->length) {
	        fprintf(fout, "%5d  ", (wbases->bases[index]).weight_match);
	    }
	}

	fprintf(fout, "\nw_id:\t");
	for (i = 0; i < bases_per_line; i++) {
	    index = j * bases_per_line + i;
	    if (index < wbases->length) {
	        fprintf(fout, "%5d  ", (wbases->bases[index]).weight_insdel);
	    }
	}
	fprintf(fout, "\n");
	if (j < (line_num - 1)) {
	    fprintf(fout, "\n");
	}
    }
    
    /* Output current approximation for short vector */
    fprintf(fout, "\nCurrent short vector:\n");
    for (i=0; i<wbases->length; i++)
        fprintf(fout, "%c", wbases->bases[i].base);
    fprintf(fout, "\n\n\n");
     
    return SUCCESS;
}


/**
 * This function extracts a Base array from the specified WeightedBases (a
 * weighted base array).
 *
 * Its synopsis is:
 * result = extract_bases_from_WB(weightedbases, bases, message)
 * where
 *         weightedbases    is the address of the specified WeightedBases
 *                          from which the Base array is extracted
 *         bases            is the address of the Base array
 *         message          is the address of the btk message object, where an
 *                          error information will be stored, if any.
 *	result		is 0 on success, !0 if an error occurs.
 */
int
extract_bases_from_WB(WeightedBases* wbases, Base* bases, 
		      int* weights, BtkMessage* message)
{
    int i;

    if (wbases == NULL) {
        sprintf(message->text, 
		"first argument of extract_bases_from_WB() is NULL");
        return ERROR;
    }

    if (bases == NULL) {
        sprintf(message->text, 
		"second argument of extract_bases_from_WB() is NULL");
        return ERROR;
    }

    if (weights == NULL) {
        sprintf(message->text, 
		"third argument of extract_bases_from_WB() is NULL");
        return ERROR;
    }

    for (i = 0; i<wbases->length; i++) {
        bases[i].base = (wbases->bases[i]).base;
	bases[i].ind = i;
	weights[i] = (wbases->bases[i]).weight_match
	  		+ (wbases->bases[i]).weight_insdel;
	fprintf(stderr, "    WB: base = %c, index= %d, wights= %d\n", 
		bases[i].base, bases[i].ind, weights[i]);
    }
    fprintf(stderr, "\n");

    return SUCCESS;
}

/**
 * This function extracts a Base array from the specified Contig.
 *
 * Its synopsis is:
 * result = extract_bases_from_contig(contig, bases, message)
 * where
 *      contig          is the address of the specified Contig
 *                      from which the Base array is extracted
 *      bases           is the address of the Base array
 *      message         is the address of the btk message object, where
 *                      an error information will be stored, if any.
 *	result		is 0 on success, !0 if an error occurs.
 */
int
extract_bases_from_contig(Contig* c, Base* bases, BtkMessage* message)
{
    int i;

    if (c == NULL) {
        sprintf(message->text, 
		"first argument of extract_bases_from_contig() is NULL");
        return ERROR;
    }

    if (bases == NULL) {
        sprintf(message->text, 
		"second argument of extract_bases_from_contig() is NULL");
        return ERROR;
    }

    for (i = 0; i<c->length; i++) {
        bases[i].base = c->sequence[i];
	bases[i].ind = i;
	fprintf(stderr, "contig: base = %c, index= %d\n", 
		bases[i].base, bases[i].ind);
    }
    fprintf(stderr, "\n");

    return SUCCESS;
}


/**
 * This function checks if the specified base is among 'ACGT'.
 *
 * Its synopsis is:
 * result = isACGT(base)
 * where
 * 	base		is the specified base character
 *	result		is 1 if the specified base is A or C or G or T
 *			or a or c or g or t; 0 otherwise.
 */
int
isACGT(char base) {
    char upper;
    upper = (char) toupper((int) base);
    if (upper == 'A' || upper == 'T' || upper == 'C' || upper == 'G') {
        return 1;
    } else {
        return 0;
    }
}


/**
 * This function checks if the specified position of the specified party
 * of the specified alignment is a terminating gap.
 *
 * Its synopsis is:
 * result = isTerminatingGap(party, position, align)
 * where
 *	party		is the character denoting one of the two parties
 *			in the alignment. either 'a' or 'r'.
 *	position	is the position in the alignment.
 *	align		is the address of the alignment.
 *	result		is 1 if the given position is a terminating gap,
 *			0 otherwise.
 */
int
isTerminatingGap(char party, int pos, GAlign* align) {
    int  i;
    char base;
    int  numOfGaps = 0;
    int  numOfBases = 0;

    // Go through the leading positions, if they are all gaps, return 1.
    for (i = 1; i < pos; i++) {
      	numOfBases++;
        if (party == 'a')
	    base = align -> achar[i];
        else
	    base = align -> rchar[i];
	if (base == '-')
	    numOfGaps++;
    }
    if (numOfGaps == numOfBases)
        return 1;

    /* The gap is not one of those leading gaps, if any.
     * Go through the trailing positions (the ones after the specified pos),
     * if they are all gaps, return 1; otherwise return 0.
     */
    numOfGaps = 0;
    numOfBases = 0;
    for (i = pos + 1; i <= (align->len); i++) {
      	numOfBases++;
        if (party == 'a')
	    base = align -> achar[i];
        else
	    base = align -> rchar[i];
	if (base == '-')
	    numOfGaps++;
    }
    if (numOfGaps == numOfBases)
        return 1;
    else
      	return 0;
}

/**
 * This function computes the result WeightedBases from the
 * specified WeightedBases and the specified alignment.
 * 
 * Its synopsis is:
 * result = compute_weightedbases(weightedbases, align, 
 *				  gapsAreDifferent, fout, message)
 * where
 *	weightedbases	is the address of the specified WeightedBases,
 *			also the address of the result WeightedBases.
 *	align		is the address of the specified alignment
 *	gapsAreDifferent is the flag indicating whether or not the
 *			terminating gaps are treated differently from
 *			the internal gaps in the alignment.
 *	fout		is the specified FILE*
 *	message		is the address of the BtkMessage that holds an
 *			error information string, if any.
 *	result		is 0 on success, !0 if an error occurs.
 *	
 * The algorithm:
 * Each weighted base has a base (base), a match weight (weight_match)
 * and a insertion-deletion weight (weight_insdel).  Global alignment
 * algorithm was used to align the bases in the weighted base array and
 * another contig to make the specified alignment.  Then from the weighted
 * base array and the alignment, we are able to calculate the new weighted
 * base array (including the result bases and result weights) using the
 * following algorithm:
 *
 * Consider at an arbitrarate position (pos) in the alignment:
 *	wb_old		as the weighted bases in the alignment
 *	wb_new		as the result weighted bases
 *	b1 		as the character at pos in wb_old
 *	b2 		as the character at pos in the contig
 *      w_m_old 	as match weight of base b1 in wb_old, if b1 != '-',
 *			which means at pos, there is a gap in the wb_old
 *	w_id_old    	as insertion-deletion weight of base b1 in wb_old,
 *			if b1 != '-'
 *	b_new		as the result base character from this alignment
 *			between b1 and b2
 *	w_m_new		as the match weight of b1 in wb_new if b1 is not
 *			deleted, or the match weight of b2 in wb_new if
 *			b1 is '-' (b2 is added to the weighted base array).
 *	w_id_new	as the insertion-deletion weight of b1 in wb_new
 *			if b1 is not deleted, or the match weight of b2 
 *			in wb_new if b1 is '-' (b2 is added to the weighted
 *			base array).
 *
 * 1) If b1 == b2 and b1 is one of ACGT and b2 is one of ACGT:
 *	b_new = b1; w_m_new = w_m_old + 1; w_id_new = w_id_new + 1;
 *	Note: it is impossible to have b1 == b2 and b1 or b2 to be either
 *	      '-' or 'N'
 *
 * 2) If b1 != b2:
 *    a) if b1 is one of ACGT and b2 is one of ACGT:
 *    	   w_m_new = w_m_old - 1; w_id_new = w_id_new - 1;
 *	   if w_m_new <= 0, b_new = 'N'; else b_new = b1
 *    b) if b1 is 'N' and b2 is one of ACGT:
 *	   b_new = b2; w_m_new = w_m_old + 1; w_id_new = w_id_new + 1;
 *    c) if b1 is '-':
 *	   i) if b1 is one of the terminating gaps and terminating gaps are
 *            treated differently
 *	   	b_new = b2; w_m_new = 1; w_id_new = 1
 *	   ii) otherwise:   
 *	   	b_new = b2; w_m_new = 1; w_id_new = 0
 *    d) if b2 is '-':
 *         i) if b2 is one of the terminating gaps and terminating gaps are
 *            treated differently
 *	        b_new = b1; w_m_new = w_m_old; w_id_new = w_id_old;
 *	   ii) otherwise:
 *	        w_m_new = w_m_old; w_id_new = w_id_old - 1;
 *	        if w_id_old < 0, b1 is deleted from the weighted base array;
 *	        otherwise, b_new = b1.
 *
 *
 * For an example:
 *   The weighted base array before alignment was:
 *	(base)		 C A C T G A C T
 *	(weight_match)   2 3 2 4 1 2 1 2
 * 	(weight_insdel)	 1 0 2 1 2 0 3 1
 *   The config used to form the specified alignment was: CCGGTAG
 *   The global alignment between the weighted bases and the contig is: *
 *
 *		1 2 3 4 5 6 7 8	9 <-- pos	
 *		C A C T G - A C	T <-- wb_old
 *		|   |   |   |
 *		C - C G G T A G	- <-- contig
 *
 *   The resulting weighted base array (wb_new) will be:
 *	(base)		 C C T G T A N T
 *	(weight_match)   3 3 3 2 1 3 0 2
 * 	(weight_insdel)	 2 3 2 3 0 1 4 1
 *
 */
int
compute_weightedbases(WeightedBases* wbases, GAlign* align,
		      int gapsAreDifferent, FILE* fout, BtkMessage* message)
{
    /* "temp" is used to temporarily hold the computed weighted bases,
     *  after the computation finished, its value is assigned back to
     *  "wbases" and its resource is released.
     */
    WeightedBases temp;
    int index_wbases, index_align, index_temp, i;
    int w_m, w_id;
    char base_a, base_r;

    if (wbases == NULL) {
        sprintf(message->text, 
		"first argument of compute_weightedbases() is NULL");
        return ERROR;
    }

    if (align == NULL) {
        sprintf(message->text, 
		"second argument of compute_weightedbases() is NULL");
        return ERROR;
    }

    if (align == NULL) {
        sprintf(message->text, 
		"third argument of compute_weightedbases() is NULL");
        return ERROR;
    }

    temp.max_length= MAX_ASSEMBLED_VECTOR_LEN;
    temp.length = align->len;
    temp.bases = CALLOC(WeightedBase, temp.length);
    MEM_ERROR(temp.bases);

    /* Iterate through the alignment vector containing a-index,
     * find its index in wbases, assign the new wbases, compute
     * new wights 
     */
    index_wbases = 0;
    index_temp = 0;

    for (index_align = (align->len); index_align > 0; index_align--) {
//        fprintf(stderr, 
//		"\nindex_align = %d  index_wbases = %d  index_temp = %d\n",
//		index_align, index_wbases, index_temp);

	base_a = (char) toupper((int) (align->achar[index_align]));
	base_r = (char) toupper((int) (align->rchar[index_align]));

//	fprintf(stderr, "base_a = %c  base_r = %c\n",
//		base_a, base_r);

//	fprintf(stderr, "wbases: base = %c, w_m = %d, w_id = %d\n",
//                (wbases->bases[index_wbases]).base,
//                (wbases->bases[index_wbases]).weight_match,
//                (wbases->bases[index_wbases]).weight_insdel);

	if (base_a == base_r ) {
	    if (isACGT(base_a) && isACGT(base_r)) {
	      // if it is a match in the alignment
	      (temp.bases[index_temp]).base = base_a;
	      (temp.bases[index_temp]).weight_match =
		  (wbases->bases[index_wbases]).weight_match + 1;
	      (temp.bases[index_temp]).weight_insdel =
		  (wbases->bases[index_wbases]).weight_insdel + 1;
	      index_temp++;
//	      fprintf(stderr, "temp: base = %c, w_m = %d, w_id = %d\n",
//		      (temp.bases[index_temp - 1]).base,
//		      (temp.bases[index_temp - 1]).weight_match,
//		      (temp.bases[index_temp - 1]).weight_insdel);
	      index_wbases ++;
	    }
	    // It is not possible that both are '-' or both are 'N'
	} else {
	    // base_a != base_r
	    if (isACGT(base_a) && isACGT(base_r)) {
	        w_m = (wbases->bases[index_wbases]).weight_match - 1;
		w_id = (wbases->bases[index_wbases]).weight_insdel + 1;

		if (w_m <= 0) {
		    (temp.bases[index_temp]).base = 'N';
		} else {
		    (temp.bases[index_temp]).base = 
		        (wbases->bases[index_wbases]).base;
		}

		(temp.bases[index_temp]).weight_match = w_m;
		(temp.bases[index_temp]).weight_insdel = w_id;
		index_temp++;
//		fprintf(stderr, "temp: base = %c, w_m = %d, w_id = %d\n",
//			(temp.bases[index_temp - 1]).base,
//			(temp.bases[index_temp - 1]).weight_match,
//			(temp.bases[index_temp - 1]).weight_insdel);
		index_wbases ++;
      	    } else if (base_a == 'N' && isACGT(base_r)) {
	        (temp.bases[index_temp]).base = base_r;
		(temp.bases[index_temp]).weight_match =
		    (wbases->bases[index_wbases]).weight_match + 1;
		(temp.bases[index_temp]).weight_insdel =
		    (wbases->bases[index_wbases]).weight_insdel + 1;
		index_temp++;
//		fprintf(stderr, "temp: base = %c, w_m = %d, w_id = %d\n",
//			(temp.bases[index_temp - 1]).base,
//			(temp.bases[index_temp - 1]).weight_match,
//			(temp.bases[index_temp - 1]).weight_insdel);
		index_wbases ++;
	    } else if (base_r == '-') {
	        if (gapsAreDifferent
			&& isTerminatingGap('r', index_align, align)) {
		    (temp.bases[index_temp]).base = base_a;
	            (temp.bases[index_temp]).weight_match
		        = (wbases->bases[index_wbases]).weight_match;
	            (temp.bases[index_temp]).weight_insdel
		        = (wbases->bases[index_wbases]).weight_insdel;
		    index_temp++;
//		    fprintf(stderr,
//			    "temp: base = %c, w_m = %d, w_id = %d\n",
//			    (temp.bases[index_temp - 1]).base,
//			    (temp.bases[index_temp - 1]).weight_match,
//			    (temp.bases[index_temp - 1]).weight_insdel);
		} else {
	            w_m = (wbases->bases[index_wbases]).weight_match;
		    w_id = (wbases->bases[index_wbases]).weight_insdel - 1;

		    if (w_id >= 0) {
		        (temp.bases[index_temp]).base = base_a;
		        (temp.bases[index_temp]).weight_match = w_m;
		        (temp.bases[index_temp]).weight_insdel = w_id;
		        index_temp++;
//		        fprintf(stderr,
//				"temp: base = %c, w_m = %d, w_id = %d\n",
//			    	(temp.bases[index_temp - 1]).base,
//			    	(temp.bases[index_temp - 1]).weight_match,
//			    	(temp.bases[index_temp - 1]).weight_insdel);
			}
			/* Otherwise (if w_id<0), this base is deleted */
		}
		index_wbases ++;
	    } else if (base_a == '-') {
	        if (gapsAreDifferent
			&& isTerminatingGap('a', index_align, align)) {
		    (temp.bases[index_temp]).base = base_r;
		    (temp.bases[index_temp]).weight_match = 1;
		    (temp.bases[index_temp]).weight_insdel = 1;
		} else {
	            (temp.bases[index_temp]).base = base_r;
		    (temp.bases[index_temp]).weight_match = 1;
		    (temp.bases[index_temp]).weight_insdel = 0;
		}
		index_temp++;
//		fprintf(stderr, "temp: base = %c, w_m = %d, w_id = %d\n",
//			(temp.bases[index_temp - 1]).base,
//			(temp.bases[index_temp - 1]).weight_match,
//			(temp.bases[index_temp - 1]).weight_insdel);
	    }
	    // other cases are not possible
	}
	
    } // for

    temp.length =  index_temp;
  
    // assign the value of "temp" back to "wbases"
    wbases->length = temp.length;
    for (i=0; i<temp.length; i++) {
        (wbases->bases[i]).base = (temp.bases[i]).base;
	(wbases->bases[i]).weight_match = (temp.bases[i]).weight_match;
	(wbases->bases[i]).weight_insdel = (temp.bases[i]).weight_insdel;
    }

    FREE(temp.bases);

    return SUCCESS;

  error:
    return ERROR;
}


/**
 * This function assembles the specified WeightedBases and the specified
 * Contig into a WeightedBases (a weighted base array).
 *
 * Its synopsis is:
 * result = assemble_fragments(weightedbases, contig, 
 *			       gapsAreDifferent, fout, message)
 * where
 *	weightedbases	is the address of the specified WeightedBases,
 *			also the address of the result WeightedBases.
 *	contig		is the address of the specified contig
 *	gapsAreDifferent is the flag indicating whether or not the
 *			terminating gaps are treated differently from
 *			the internal gaps in the alignment.
 *	fout		is the specified FILE*
 *	message		is the address of the BtkMessage that holds an
 *			error information string, if any.
 *	result		is 0 on success, !0 if an error occurs.
 *	
 */
int
assemble_fragments(WeightedBases* vectorw, Contig* fragc,
	           int gapsAreDifferent, FILE* fout, BtkMessage* message) 
{
    Base *vectorb, *fragb;
    int  *weights;
    GAlign align;
    int r;

    if (vectorw == NULL) {
        sprintf(message->text, 
		"first argument of assemble_fragments() is NULL");
        return ERROR;
    }

    if (fragc == NULL) {
        sprintf(message->text, 
		"second argument of assemble_fragments() is NULL");
        return ERROR;
    }

    if (fout == NULL) {
        sprintf(message->text, 
		"third argument of assemble_fragments() is NULL");
        return ERROR;
    }

    if (vectorw->length == -1) {
        // no bases in "vectorw"
        r = weightedbases_create(vectorw, fragc, message);

	if (r == ERROR) {
	    fprintf(stderr, "%s\n", message->text);
	}

	return r;
    } 

    fragb = NULL;
    vectorb = NULL;
    weights = NULL;
    GAlign_init(&align, message);

    fragb = CALLOC(Base, fragc->length);
    MEM_ERROR(fragb);
    r = extract_bases_from_contig(fragc, fragb, message);
    if (r == ERROR) {
          fprintf(stderr, "%s\n", message->text);
	  goto error;
    }

    vectorb = CALLOC(Base, vectorw->length);
    MEM_ERROR(vectorb);
    weights = CALLOC(int, vectorw->length);
    MEM_ERROR(weights);
    r = extract_bases_from_WB(vectorw, vectorb, weights, message);
    if (r == ERROR) {
          fprintf(stderr, "%s\n", message->text);
	  goto error;
    }

    // Globally align the bases from the "vectorw" and the "fragc"
    r = align_fragments(vectorb, vectorw->length, weights,
			fragb, fragc->length,
			g3, p3, &align, fout, message);
    if (r == ERROR) {
          fprintf(stderr, "%s\n", message->text);
	  goto error;
    }

    // Compute the new weightedbases from the global alignment
    r = compute_weightedbases(vectorw, &align, gapsAreDifferent, fout, message);
    if (r == ERROR) {
          fprintf(stderr, "%s\n", message->text);
	  goto error;
    }
    
    GAlign_release(&align, message);
    FREE(vectorb);
    FREE(weights);
    FREE(fragb);
    return SUCCESS;

  error:
    GAlign_release(&align, message);
    FREE(vectorb);
    FREE(weights);
    FREE(fragb);
    return ERROR;
}
