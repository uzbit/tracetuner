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
 * $Id: Btk_match_data.h,v 1.7 2009/01/01 16:41:02 gdenisov Exp $                   
 */

#ifndef BTK_MATCH_DATA_H
#define BTK_MATCH_DATA_H

/**********************************************************************/
/*                Contig definition and functions.                    */
/**********************************************************************/

/* A Contig is a data structure consisting of a contiguous sequence.
 * It may also contain a lookup table corresponding to that sequence,
 * or it may not.
 */

#define ALPHABET_SIZE (5)
    
    typedef struct {
       int    length;
       int   *sizes;	/* Sizes of each data array. */
       long **data;
    } lookup_table;
    
    typedef struct 
    {
       char*	    sequence;
       uint8_t     *qv;
       int    	    max_length;
       int    	    length;
       lookup_table lut;
    } Contig;
    
    typedef Contig SubContig;

/**********************************************************************/
/*             Align_params definition and functions.                 */
/**********************************************************************/

/* An Align_params structure contains the substitution matrix, insert
 * penalty and delete penalty necessary to do an alignment.
 */

    typedef struct {
       int**  matrix;
       int    matrix_row_len;
       int    gap_init;
       int    gap_ext;
    } Align_params;
    
/**********************************************************************/
/*                 Align definition and functions.                    */
/**********************************************************************/

/* An Align structure contains information returned from an alignment:
 *     score             score of the alignment
 *     base_is_reverse   0 if the data aligns with the query
 *                       1 if its reverse complement aligns
 *     trace_qpos        position in the query of the character in trace_qchar
 *     trace_dpos        position in the query of the character in trace_dchar
 *     trace_dir         the directions taken in the SW algorithm:
 *                       match/mismatch, insert, delete
 *     trace_qchar       the characters (including gaps) of the query
 *     trace_mchar       '|'s at positions where query and data match
 *     trace_dchar       the characters (including gaps) of the data
 *     trace_len         the length of the alignment (counts gaps)
 *     trace_max_len     the amount of space in the various arrays
 *     contig_offset     the position in the data that corresponds to 
 *                       query position 0
 *     num_gaps          the maximum allowable number of gaps.
 */
    typedef struct 
    { 
       int    score;
       int    base_is_reverse;
       int*   trace_qpos;
       int*   trace_dpos;
       char*  trace_dir;
       char*  trace_qchar;
       char*  trace_mchar;
       char*  trace_dchar;
       int    trace_len;
       int    trace_max_len;
       int    contig_offset;
       int    num_gaps;
    } Align;
    
/**********************************************************************/
/*                        Range definition.                           */
/**********************************************************************/
    typedef struct {
       int begin;
       int end;
    } Range;

/**********************************************************************/
/*                        Vector definition.                          */
/**********************************************************************/

/* A Vector structure contains vector information, namely
 *
 * preCut	the sequence from primer to restriction site
 * postCut	the sequence after restriction site
 * is_reverse	0 if the primer and restriction site match the vector
 * 		1 if the reverse complement
 * primerStart	the index in the vector where the primer starts
 *		NOTE: always w.r.t. forward vector
 */

    typedef struct {
	Contig preCut;
	Contig postCut;
	int is_reverse;
	int primerStart;
    } Vector;

    extern int local_read_fasta(char *, Contig *, BtkMessage *);
    extern int readVector( char *, char *, char *, Vector*, Align_params *,
        BtkMessage *);
    extern int readShortVector(char *, Vector *, BtkMessage *);
    extern int contig_init(Contig *, BtkMessage*);
    extern int contig_create(Contig *, char*, int, uint8_t *, BtkMessage*);
    extern int contig_copy(Contig *, Contig *, BtkMessage*);
    extern int contig_release(Contig *, BtkMessage*);
    extern int contig_make_fasta_lookup_table( Contig*, int, BtkMessage*);
    extern int contig_complement(Contig*, BtkMessage*);
    extern int contig_reverse(Contig*, BtkMessage*);
    extern int contig_get_reverse_comp(Contig*, Contig*, BtkMessage*);
    extern int contig_create_sub(SubContig*, Contig*, int, int, BtkMessage*);
    extern int contig_fprint(FILE*, Contig*, BtkMessage*);
    extern int base2int( char );
    extern int set_alignment_parameters(Align_params *, int, int, int, int, 
        BtkMessage*);
    extern int set_alignment_parameters_IUB(Align_params *, int ,
        int, int, int, BtkMessage*);
    extern int alignment_parameters_release(Align_params *, BtkMessage *);
    extern int align_init( Align*, BtkMessage*);
    extern int align_create( Align*, int, int, int, BtkMessage*);
    extern int align_fprint( FILE*, Align*, int, BtkMessage*);
    extern int align_release( Align*, BtkMessage* );
#endif 
