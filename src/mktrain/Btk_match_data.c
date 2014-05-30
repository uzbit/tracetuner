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
 * $Id: Btk_match_data.c,v 1.11 2009/01/06 18:18:12 gdenisov Exp $                 
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "Btk_qv.h"
#include "util.h"
#include "Btk_match_data.h"

#define LOOKUP_SIZE       (10)
#define MAX_SIZE_OF_BASES (5000)
#define FASTA_LEN 1000

/**********************************************************************/
/*                     Vector functions:                              */
/**********************************************************************/

/**
 * This function reads a single FastA formatted entry from a file, and
 * constructs the specified contig with the bases read from the file..
 * FastA format: lines beginning with '>' signify comments
 *               other lines signify data.
 *
 * Its synopsis is:
 * result = local_read_fasta(file_name, contig, message)
 *      file_name       is the address of the specified FASTA file name
 *      contig          is the address of the specified contig
 *      message         is the address of a BtkMessage where information
 *                      about an error will be put, if any
 *      result          is 0 on success, !0 if an error occurs
 */
int
local_read_fasta(char* file_name, Contig* contig, BtkMessage *message)
{
    FILE            *infile;
    char            buffer[BUFLEN];
    int             curr_len, in_data, line_len;
    char           *sequence;
    uint8_t        *quality_values = NULL;
    int             i, r;
    int             length;

    if ((infile=fopen(file_name,"r"))== NULL) {
        sprintf(message->text, "Unable to open FastA file '%s'\n", file_name);
        return ERROR;
    }

    length =  FASTA_LEN;
    sequence = (char*) malloc( sizeof(char)*FASTA_LEN );
    MEM_ERROR( sequence );

    curr_len=0;
    sequence[curr_len] = '\0';
    in_data = -1;

    while ((fgets(buffer, BUFLEN, infile)) != NULL) {
        if ( buffer[0] == '>') {
            /* Header line. */
            if ( in_data==1 ) {
                /* We've read one sequence, and we hit the next
                 * header line.
                 */
                break;
            }
            in_data=0;
        } else {
            if ( in_data == -1 ) {
                fprintf(stderr, "Warning: sequence did not contain a ");
                fprintf(stderr, "valid '>' header line! - file = %s\n",
                        file_name);
            }
            in_data=1;

            /* Take care of carriage returns. */
            line_len  = strlen(buffer);
            if (buffer[line_len-1] == '\n')  {
                buffer[line_len-1] = '\0';
                line_len--;
            }

            curr_len += line_len;

            /* Allocate more memory, if necessary. */
            if (curr_len >= length) {
                length += FASTA_LEN;
                sequence = REALLOC(sequence, char, length);
                MEM_ERROR( sequence );
            }

            strncpy(&sequence[curr_len-line_len],buffer,line_len+1);
        }
    }

    length = strlen(sequence);
    for (i=0; i<length; i++)
        sequence[i] = (char)toupper((int)sequence[i]);    

    contig_create( contig, sequence, length, quality_values, message );

    r=SUCCESS;
    goto cleanup;

error:
    r=ERROR;

cleanup:
    fclose(infile);
    FREE( sequence );

    return r;
}

/* This function reads in vector information and creates a Vector data
 * structure.  Its synopsis is:
 *
 * result = readVector( VectorName, PrimerName, SiteName, vector, ap, message)
 *
 * where
 *      VectorName      is the name of a file containing the vector
 *      PrimerName      is the name of a file containing the primer
 *      SiteName        is the name of a file containing the restriction site
 *      vector          is the address of the Vector to be created
 *      ap              is the address of an Align_params data structure,
 *                      correctly set for IUB matching
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 */
int
readVector( char* VectorName, char* PrimerName, char* SiteName, Vector* vector,
            Align_params *ap, BtkMessage* message)
{
    Contig vectorContig;
    Contig consensusrc;
    Contig site;
    Contig primer;
    int    r;
    int    cutPoint;
    int    siteStart;
    char   *s, *s2;
    uint8_t *quality_values = NULL;

    /* Read and complement vector. */
    r=local_read_fasta( VectorName, &vectorContig, message );
    if ( r==ERROR ) { return ERROR; }
    r=contig_get_reverse_comp( &consensusrc, &vectorContig, message );
    if ( r==ERROR ) { return ERROR; }

    /* Read primer and restriction site. */
    r=local_read_fasta( PrimerName, &primer, message );
    if ( r==ERROR ) { return ERROR; }
    r=local_read_fasta( SiteName, &site, message );
    if ( r==ERROR ) { return ERROR; }

    /* Site will be of the form XXX^YYY.  Find the cut-point (i.e., ^),
     * store its location, and remove it from the string.
     */
    s=strchr( site.sequence, '^');
    if ( s == NULL) {
        sprintf(message->text, "No ^ found in restriction site.\n");
        return ERROR;
    } else {
        cutPoint = (int) (s-site.sequence)/sizeof(char);
        memmove(s, s+1, (strlen(s)-1)*sizeof(char) );
        site.length--;
        site.sequence[ site.length ] = '\0';
    }
    if ( strchr( site.sequence, '^') != NULL ) {
        sprintf(message->text,
                "At least two ^'s found in restriction site.\n");
        return ERROR;
    }

    /* Find primer and site in vector or its complement.
     * Here we assume it only occurs in one or the other.
     */
    s = strstr( vectorContig.sequence, primer.sequence );
    if ( s != NULL ) {
        vector->primerStart = (int) (s - vectorContig.sequence)/sizeof(char);
        vector->is_reverse=0;
    } else {
        s = strstr( consensusrc.sequence, primer.sequence );
        if ( s == NULL ) {
            sprintf(message->text, "Primer not found in vector.\n");
            return ERROR;
        }
        vector->primerStart = consensusrc.length
                                - (int) (s - consensusrc.sequence)/sizeof(char)
                                - 1;
        vector->is_reverse=1;
    }

    s2 = strstr( s, site.sequence );
    if ( s == NULL ) {
        sprintf(message->text, "Restriction site not found in vector.\n");
        return ERROR;
    }
    siteStart = (int) (s2 - s)/sizeof(char);

    /* Split vector into pre-,post-Cut parts. */
    if (contig_create( &(vector->preCut), s, siteStart+cutPoint, 
        quality_values, message) != SUCCESS)
        return ERROR;
    
    if (contig_create( &(vector->postCut), &s2[cutPoint], strlen(&s2[cutPoint]),
        quality_values,  message) != SUCCESS)
        return ERROR; 

    /* Clean up. */
    r=contig_release( &vectorContig, message );
    if ( r==ERROR ) { return ERROR; }
    r=contig_release( &consensusrc, message );
    if ( r==ERROR ) { return ERROR; }
    r=contig_release( &site, message );
    if ( r==ERROR ) { return ERROR; }
    r=contig_release( &primer, message );
    if ( r==ERROR ) { return ERROR; }

    return SUCCESS;
}

/******************************************************************************
 * This function reads in short vector information and creates a Vector
 * data structure.  Its synopsis is:
 *
 * result = readShortVector(vectorName, vector, message)
 * where
 *      vectorName      is the name of a file containing the short vector
 *      vector          is the address of the Vector to be created
 *      message          is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 ******************************************************************************
 */
int
readShortVector(char *vectorName, Vector *vector, BtkMessage *message)
{
    Contig shortVectorContig;
    int r;
    uint8_t *qv;

    contig_init(&shortVectorContig, message);

    /* Read the short vector */
    r = local_read_fasta(vectorName, &shortVectorContig, message);
    if (r == ERROR) {
        contig_release(&shortVectorContig, message);
        return ERROR;
    }
    qv = CALLOC(uint8_t, shortVectorContig.length);

    r = contig_create(&(vector->preCut), shortVectorContig.sequence,
        shortVectorContig.length, qv, message);
    if (r == ERROR) {
        contig_release(&shortVectorContig, message);
        return ERROR;
    }
    FREE(qv);

    /* When a user provides the short vector information, postCut
       will not be available. */
    (vector->postCut).sequence = NULL;
    vector->primerStart = 0;
    vector->is_reverse = 0;

    /* Clean up. */
    contig_release(&shortVectorContig, message);
    return SUCCESS;
}


/**********************************************************************/
/*                     Contig functions:                              */
/**********************************************************************/

/**
 * This function initializes a Contig.  Its synopsis is:
 * 
 * result = contig_init(contig, message)
 * where
 *      contig          is the address of the Contig to be initialized
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 */
int
contig_init(Contig *contig, BtkMessage* message)
{
    contig->sequence = NULL;
    contig->length = 0;
    contig->max_length= 0;
    (contig->lut).length = 0;
    return SUCCESS;
} 

/*
 * This function creates a Contig from a string.  Its synopsis is:
 *
 * result = contig_create(contig, bases, length, message)
 *
 * where
 *      contig          is the address of the Contig to be created
 *      bases           is a string of bases that will become the sequence
 *                      of the Contig
 *      length          is the length of the sequence 'bases'
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 * Note:
 *     The Contig's sequence field will be made uppercase and will end
 *         with '\0'.
 *     The Contig's length does not include the '\0'.
 *     The Contig's lookup table will be initialized to NULL; its length
 *         will be set to 0.
 *     This function should not be called on a Contig that already has
 *         data stored in it!
 */
int
contig_create(Contig * contig, char* bases, int length, 
    uint8_t *quality_values, BtkMessage* message)
{
    int i;

    contig->sequence = (char*) malloc(sizeof(char)*(length+1));
    MEM_ERROR(contig->sequence);

    contig->qv = (uint8_t*) malloc(sizeof(uint8_t)*length);
    MEM_ERROR(contig->qv);

    contig->max_length = length;
    contig->length = length;

    for (i=0; i<length; i++) {
        contig->sequence[i] = (char) toupper( (int) bases[i] );
        if (quality_values != NULL)
            contig->qv[i] = quality_values[i];
        else
            contig->qv[i] = 0;
    }
    contig->sequence[length] = '\0';

    contig->lut.sizes=NULL;
    contig->lut.data=NULL;
    contig->lut.length=0;

    return SUCCESS;

    error: 
        return ERROR;
}


/*
 * This function creates a Contig by copying an existing Contig.
 * Its synopsis is:
 *
 * result = contig_create(to, from, message)
 *
 * where
 *      to              is the address of the Contig to be created
 *      from            is the address of the existing Contig
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 * Note:
 *     The Contig's sequence field will be made uppercase.
 *     The Contig's lookup table will be initialized to NULL; its length
 *         will be set to 0.
 *     This function should not be called on a 'to' Contig that already has
 *         data stored in it!
 */
int
contig_copy(Contig * to, Contig* from, BtkMessage* message)
{
    return contig_create( to, from->sequence, from->length, from->qv, 
        message);
}


/*
 * This function frees the memory held by a Contig.  Its synopsis is:
 *
 * result = contig_release(contig, message)
 *
 * where
 *      contig          is the address of the Contig to be released
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 * Note:
 *     This function will check explicitly whether lut has any data in it.
 *     If so, it will free that memory.  If not, it will not even try.
 */
int 
contig_release(Contig * contig, BtkMessage* message)
{
    int i;

    FREE(contig->sequence);
    FREE(contig->qv);
    contig->max_length=0;
    contig->length=0;

    if ( contig->lut.length > 0 )
    {
	if (contig->lut.data != NULL) 
	{
            for( i=0; i < contig->lut.length; i++) 
	    {
             	FREE( contig->lut.data[i] );
            }
	}
        FREE( contig->lut.data );
 
        FREE( contig->lut.sizes );
     }

    return SUCCESS;
}


/*
 * This function returns a number corresponding to a base, namely:
 *     A = 0, C = 1, G = 2, T = 3, other = 4.
 * Its synopsis is:
 *
 * num = base2int( base )
 *
 * where
 *      base            is a base, in character form
 *
 *      num             is the resulting integer
 *
 */
int
base2int( char base )
{
    switch ( base ) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default : return 4;
    }
}


/*
 * This function populates the lut field with a lookup table, for
 * use in the FastA algorithm.  Its synopsis is:
 *
 * result = contig_make_fasta_lookup_table(contig, ktup, message)
 *
 * where
 *      contig          is the address of the Contig whose lookup table
 *                      needs to be created
 *      ktup            is an integer specifying what size k-tuples should
 *                      be used.  For example, if ktup=6 (recommended for DNA)
 *                      6-tuples will be used.
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 * The algorithm:
 * Consider each consecutive K-tuple, which we will reference
 *     by ending positions.
 *     For exmple, if ktup = 6, consider 6-tuples:
 *         0123456789...
 *         ACAGCTGCTCGAACGATGATAGCACAAACGCCCCCTGTG
 *         555555
 *          666666
 *           777777
 *            888888
 *             etc.
 *
 * Convert the k-tuple into an integer index,
 *    where index(x) is
 *        (( x0*ALPHABET_SIZE + x1 )*ALPHABET_SIZE + x2)... + xk.
 *    Assume, for example, that ktup = 6 and ALPHABET_SIZE = 5 
 *    then, ktuple #5,
 *        => ACAGCT
 *        => ( 0, 1, 0, 2, 1, 3 )  * using base2int *
 *        => ((((( 0 * 5 ) + 1 ) * 5 + 0 ) * 5 + 2 ) * 5 + 1 ) * 5 + 3
 *        => 683
 *
 * Put the position of the current K-tuple into the lookup table at position
 *     index(x).
 *
 *     In this example, put the value 5 into the lookup table at position 683 
 */
int
contig_make_fasta_lookup_table( Contig *contig, int ktup, BtkMessage *message)
{
    int *lookup_table_actual=NULL;
    int index, i, r;

    if( contig->lut.length != 0 )
    {
        sprintf(message->text,
                "There is already a lookup table for this Contig.\n");
        return ERROR;
    }

    for (i=0, contig->lut.length=1; i<ktup; i++)
    {
        contig->lut.length *= ALPHABET_SIZE;
    }

    contig->lut.data = (long**) malloc(sizeof(long*)*(contig->lut.length));
    MEM_ERROR(contig->lut.data);

    contig->lut.sizes = (int*) malloc( sizeof(int)*(contig->lut.length) );
    MEM_ERROR(contig->lut.sizes);

    lookup_table_actual = (int*) malloc( sizeof(int)*(contig->lut.length) );
    MEM_ERROR(lookup_table_actual);

    if( contig->lut.length < ktup )
    {
        sprintf(message->text, "Sequence is too small.\n");
        return ERROR;
    }

    /* Initialize the lookup table. */
    for( i=0; i< contig->lut.length; i++)
    {
        lookup_table_actual[i]=0;
        contig->lut.sizes[i]=0;
        contig->lut.data[i]=NULL;
    }

    index=0;
    for( i=0; i< contig->length; i++ )
    {
        index *= ALPHABET_SIZE;
        index = index % (contig->lut.length);
        index += base2int( contig->sequence[i] );
        if ( i < ktup-1 ) continue; /* Not enough bases so far. */
        contig->lut.sizes[ index ] ++;
        if ( contig->lut.sizes[ index ] > lookup_table_actual[ index ])
        {
            lookup_table_actual[index] += LOOKUP_SIZE;
            contig->lut.data[index] =
                (long*) realloc(contig->lut.data[index],
                                sizeof(long) * (lookup_table_actual[index]) );
        }
        contig->lut.data[index][contig->lut.sizes[index]-1]=i;
    }


    r=SUCCESS;
    goto cleanup;

    error:
        r=ERROR;

    cleanup:
        FREE( lookup_table_actual );
	
    return r;
}


/*
 * This function changes the sequence of a Contig to its complement.
 * Its synopsis is:
 *
 * result = contig_complement(contig, message)
 *
 * where
 *      contig          is the address of the Contig to complement
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 * Note:
 *     This will ruin the lookup table, if it has already been populated.
 */
int 
contig_complement(Contig* contig, BtkMessage* message)
{
    int i;
    char c;

    for(i=0;i< contig->length;i++) {
        if      (contig->sequence[i]=='A') c='T';
        else if (contig->sequence[i]=='T') c='A';
        else if (contig->sequence[i]=='C') c='G';
        else if (contig->sequence[i]=='G') c='C';
        else if (contig->sequence[i]=='O') c='J';
        else if (contig->sequence[i]=='J') c='O';

      /* IUB symbols */
        else if (contig->sequence[i]=='M') c='K';
        else if (contig->sequence[i]=='R') c='Y';
        else if (contig->sequence[i]=='W') c='W';

        else if (contig->sequence[i]=='S') c='S';
        else if (contig->sequence[i]=='Y') c='R';
        else if (contig->sequence[i]=='K') c='M';
        else if (contig->sequence[i]=='V') c='B';
        else if (contig->sequence[i]=='H') c='D';
        else if (contig->sequence[i]=='D') c='H';
        else if (contig->sequence[i]=='B') c='V';
  
        else c=contig->sequence[i];

	contig->sequence[i]=c;
    }

    return SUCCESS;
}


/*
 * This function reverses the sequence of a Contig.  Its synopsis is:
 *
 * result = contig_reverse(contig, message)
 *
 * where
 *      contig          is the address of the Contig to reverse
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 * Note:
 *     This will ruin the lookup table, if it has already been populated.
 */
int 
contig_reverse(Contig* contig, BtkMessage* message)
{
    int i;
    char c;

    for( i=0; i < contig->length/2; i++ )
    {
        c=contig->sequence[i];
        contig->sequence[i] = contig->sequence[contig->length-i-1];
        contig->sequence[contig->length-i-1] = c;
    }
    return SUCCESS;
}


/*
 * This function creates a new Contig that is the reverse complement
 * of a given Contig.  Its synopsis is:
 *
 * result = contig_get_reverse_comp(to, from, message)
 *
 * where
 *      to              is the address of the new Contig
 *      from            is the address of the original Contig
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 * Note:
 *     All provisos from contig_copy, contig_reverse and contig_complement
 *     apply here as well.
 */
int 
contig_get_reverse_comp(Contig* rev_comp, Contig* contig, BtkMessage* message)
{
    int r;

    r=contig_copy( rev_comp, contig, message );
    if ( r==ERROR ) {
        return ERROR;
    }

    r=contig_complement( rev_comp, message );
    if ( r==ERROR ) {
        return ERROR;
    }

    r=contig_reverse( rev_comp, message );
    if ( r==ERROR ) {
        return ERROR;
    }

    return SUCCESS;
}

/*
 * This function creates a SubContig, a structure that points to,
 * but does not copy, the original contig.  Its synopsis is:
 *
 * result = contig_create_sub(subcontig, contig, offset, length, message)
 *
 * where
 *      subcontig       is the address of the new SubContig
 *      contig          is the address of the original Contig
 *      offset          is the offset in the original Contig where the
 *                      SubContig should start
 *      length          is the length of the SubContig (if this length
 *                      exceeds the end of the original Contig, an error
 *                      is returned)
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 * Note:
 *     A SubContig only points to the original structure, it does not
 *         copy it, so a SubContig should not be free'd or released.
 *
 *     The lookup table for the original Contig is not valid for the
 *         SubContig, so the new lookup table is initialized to NULL,
 *         length = 0.
 */
int 
contig_create_sub(SubContig* subcontig, Contig* contig, int offset,
    int length, BtkMessage* message)
{
    subcontig->sequence = contig->sequence + offset;
    if ( length > contig->length - offset ) {
        sprintf(message->text, "Sub-contig length is too large.\n");
        return ERROR;
    } else {
        subcontig->length = length;
    }
    subcontig->qv = contig->qv + offset;

    subcontig->max_length = -1;
    subcontig->lut.sizes=NULL;
    subcontig->lut.data=NULL;
    subcontig->lut.length=0;

    return SUCCESS;
}


/*
 * This function prints a Contig to the specified output stream.
 * Its synopsis is:
 *
 * result = contig_fprint(out, contig, message)
 *
 * where
 *      out             is the address of a FILE
 *      contig          is the address of the Contig
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 * Note:
 *     Prints lines of length 70.  One might consider parameterizing that.
 */
int 
contig_fprint(FILE* fout, Contig* contig, BtkMessage* message)
{
    int i;

    for(i=0; i < contig->length; i++)
    {
       if ( (i+1)%70 == 0 ) fprintf(fout, "%c\n", contig->sequence[i]);
       else                 fprintf(fout, "%c",   contig->sequence[i]);
    }
    if ( i%70 != 0 ) fprintf(fout, "\n");

    return SUCCESS;
}

/**********************************************************************/
/*                     Align_params functions:                        */
/**********************************************************************/


/*
 * This function sets the substitution matrix for an Align_params
 * structure.  Its synopsis is:
 *
 * result = align_set_matrix(align_par, matrix, row_len, message)
 *
 * where
 *      align_par       is the address of the Align_params structure
 *      matrix          is the matrix (i.e., **int)
 *      row_len         is the size of the matrix
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 * Note:
 *      Assigns only, does not copy.
 */
static int 
align_set_matrix(Align_params* align_par,int** matrix, int matrix_row_len,
                        BtkMessage* message)
{
    align_par->matrix = matrix;
    align_par->matrix_row_len = matrix_row_len;

    return SUCCESS;
}


/*
 * This function sets the insert and delete scores for an Align_params
 * structure.  Its synopsis is:
 *
 * result = align_set_indel(align_par, ins, del, message)
 *
 * where
 *      align_par       is the address of the Align_params structure
 *      ins             is the insert score (e.g. -16)
 *      del             is the delete score (e.g. -16)
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 */
static int 
align_set_indel(Align_params* align_params, int init, int ext, 
                       BtkMessage* message)
{
    align_params->gap_init = init;
    align_params->gap_ext  = ext;

    return SUCCESS;
}


/*
 * This function sets the parameters for an Align_params structure.
 * Its synopsis is:
 *
 * result = set_alignment_parameters(align_par, match, mismatch, ins, del,
                                     message)
 *
 * where
 *      align_par       is the address of the Align_params structure
 *      match           is the match score (e.g. 10)
 *      mismatch        is the mismatch score (e.g. -20)
 *      ins             is the insert score (e.g. -16)
 *      del             is the delete score (e.g. -16)
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 */
int
set_alignment_parameters( Align_params *ap, int match, int mismatch,
                          int ginit, int gext, BtkMessage* message)
{
    int    i,j;
    int    r=0;
    int    matrix_row_len = 256;

/* Create a matrix of size 256x256 */
    ap->matrix_row_len = matrix_row_len;
    ap->matrix = CALLOC(int *, ap->matrix_row_len);
    for(i=0;i<ap->matrix_row_len;i++) {
        ap->matrix[i] = CALLOC(int, ap->matrix_row_len);
    }

    for(i=0; i<256; i++) {
        for(j=0;j<256;j++) {
            if(i==j) {
                ap->matrix[i][j] = match;
            }
            else {
                ap->matrix[i][j] = mismatch;
            }
        }
    }

    r = align_set_matrix(ap, ap->matrix, ap->matrix_row_len, message);
    if(r==ERROR) {
        return ERROR;
    }

    r = align_set_indel(ap, ginit, gext, message);
    if(r==ERROR) {
        return ERROR;
    }

    return SUCCESS;
}


/*
 * This function sets the parameters for an Align_params structure according
 * to IUB conventions. Its synopsis is:
 *
 * result = set_alignment_parameters_IUB(align_par, match, mismatch, ins, del,
                                     message)
 *
 * where
 *      align_par       is the address of the Align_params structure
 *      match           is the match score (e.g. 10)
 *      mismatch        is the mismatch score (e.g. -20)
 *      ins             is the insert score (e.g. -16)
 *      del             is the delete score (e.g. -16)
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 * IUB code	Meaning	
 *    A           A
 *    C           C
 *    G        	  G
 *   T/U          T
 *    M         A or C
 *    R         A or G
 *    W         A or T
 *    S         C or G
 *    Y         C or T
 *    K         G or T
 *    V         A or C or G
 *    H         A or C or T
 *    D         A or G or T
 *    B         C or G or T
 *    N         A or C or G or T
 *
 * We will assume at least one of the two sequences has actual ACTG,
 * i.e. N matches with A, but N does not match with N.
 */
int
set_alignment_parameters_IUB( Align_params *ap, int match, int mismatch,
                              int ginit, int gext, BtkMessage* message)
{
    int    i,j;
    int    r=0;
    int**  matrix;
    int    matrix_row_len;

/* Create a matrix of size 256x256 */
    matrix_row_len = 256;
    matrix = (int**) malloc(sizeof(int*)*matrix_row_len);
    for(i=0;i<matrix_row_len;i++) {
        matrix[i] = (int*) malloc(sizeof(int)*matrix_row_len);
    }

    for(i=0; i<256; i++) {
        for(j=0;j<256;j++) {
            if(i==j) {
                matrix[i][j] = match;
            }
            else {
                matrix[i][j] = mismatch;
            }
        }
    }

    /* Singles */
    matrix['T']['U'] = matrix['U']['T'] = match;

    /* Doubles */
    matrix['M']['A'] = matrix['A']['M'] = match;
    matrix['M']['C'] = matrix['C']['M'] = match;

    matrix['R']['A'] = matrix['A']['R'] = match;
    matrix['R']['G'] = matrix['G']['R'] = match;

    matrix['W']['A'] = matrix['A']['W'] = match;
    matrix['W']['T'] = matrix['T']['W'] = match;

    matrix['S']['C'] = matrix['C']['S'] = match;
    matrix['S']['G'] = matrix['G']['S'] = match;

    matrix['Y']['C'] = matrix['C']['Y'] = match;
    matrix['Y']['T'] = matrix['T']['Y'] = match;

    matrix['K']['G'] = matrix['G']['K'] = match;
    matrix['K']['T'] = matrix['T']['K'] = match;

    /* Triples */
    matrix['V']['A'] = matrix['A']['V'] = match;
    matrix['V']['C'] = matrix['C']['V'] = match;
    matrix['V']['G'] = matrix['G']['V'] = match;

    matrix['H']['A'] = matrix['A']['H'] = match;
    matrix['H']['C'] = matrix['C']['H'] = match;
    matrix['H']['T'] = matrix['T']['H'] = match;

    matrix['D']['A'] = matrix['A']['D'] = match;
    matrix['D']['G'] = matrix['G']['D'] = match;
    matrix['D']['T'] = matrix['T']['D'] = match;

    matrix['B']['C'] = matrix['C']['B'] = match;
    matrix['B']['G'] = matrix['G']['B'] = match;
    matrix['B']['T'] = matrix['T']['B'] = match;

    /* Quadruples. */
    matrix['N']['A'] = matrix['A']['N'] = match;
    matrix['N']['C'] = matrix['C']['N'] = match;
    matrix['N']['G'] = matrix['G']['N'] = match;
    matrix['N']['T'] = matrix['T']['N'] = match;

    r = align_set_matrix(ap, matrix,matrix_row_len,message);
    if(r==ERROR) {
        return ERROR;
    }

    r = align_set_indel(ap, ginit, gext, message);
    if(r==ERROR) {
        return ERROR;
    }

    return SUCCESS;
}



/*
 * This function releases the memory held by an Align_params structure.
 * Its synopsis is:
 *
 * result = alignment_parameters_release(ap, message)
 *
 * where
 *      align_par       is the address of the Align_params structure
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 */
int
alignment_parameters_release(Align_params* ap, BtkMessage* message)
{
    int i;

    for(i=0;i<ap->matrix_row_len;i++) {
        FREE(ap->matrix[i]); 
    }
    FREE(ap->matrix); 
    return SUCCESS;
}

/**********************************************************************/
/*                     Align functions:                               */
/**********************************************************************/

/**
 * This function releases the memory resource held by the specified Align
 * struncture.  Its synopsis is:
 *
 * result = align_release(align, message)
 *
 * where
 *      align           is the address of the Align structure
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 */
int
align_init(Align *align, BtkMessage *message)
{
    if (align == NULL) {
        return ERROR;
    }

    align->trace_qpos = NULL;
    align->trace_dpos = NULL;
    align->trace_dir = NULL;
    align->trace_qchar = NULL;
    align->trace_mchar = NULL;
    align->trace_dchar = NULL;

    align->trace_len = 0;
    align->trace_max_len = 0;
    align->contig_offset = 0;

    return SUCCESS;
}



/*
 * This function creates an Align structure. Its synopsis is:
 *
 * result = align_create(align, num_bases, offset, is_complement, message)
 *
 * where
 *      align           is the address of the Align structure
 *      num_bases       is the length of the query sequence
 *                      ? Evidently, this is not used.
 *      offset          is the position in the data sequence
 *                      that corresponds to query position 0
 *      is_complement   is 1 if the query matches the reverse complement,
 *                      or 0 if the query matches the original sequence
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 */
int
align_create(Align * align, int num_called_bases, 
             int offset, int is_complement, BtkMessage* message)
{
    align->trace_qpos = (int*) malloc(sizeof(int)*MAX_SIZE_OF_BASES);
    MEM_ERROR(align->trace_qpos);
 
    align->trace_dpos = (int*) malloc(sizeof(int)*MAX_SIZE_OF_BASES);
    MEM_ERROR(align->trace_dpos);
 
    align->trace_dir = (char*) malloc(MAX_SIZE_OF_BASES);
    MEM_ERROR(align->trace_dir);
 
    align->trace_qchar = (char*) malloc(MAX_SIZE_OF_BASES);
    MEM_ERROR(align->trace_qchar);
 
    align->trace_mchar = (char*) malloc(MAX_SIZE_OF_BASES);
    MEM_ERROR(align->trace_mchar);
 
    align->trace_dchar = (char*) malloc(MAX_SIZE_OF_BASES);
    MEM_ERROR(align->trace_dchar);
 
    align->trace_len = 0;
    align->trace_max_len = MAX_SIZE_OF_BASES;

    /* Insert offset info */ 
    align->contig_offset = offset;
    align->base_is_reverse = is_complement;
    align->num_gaps = 0;   
 
    return SUCCESS;

    error:
        return ERROR;
}


/*
 * This function prints an Align structure to the specified file stream.
 * Its synopsis is:
 *
 * result = align_fprint(out, align, display_size, message)
 *
 * where
 *      out             is the address of a FILE
 *      align           is the address of the Align structure
 *      display_size    is the amount of alignment to show per line,
 *                      e.g. display_size=70 will show 70 characters of
 *                      alignment, and several other of header, per line
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 */
int
align_fprint(FILE *fout, Align* align, int display_size, BtkMessage* message)
{
    char  *seq1=NULL, *seq2=NULL, *matchchar=NULL, *tickmarks=NULL;
    int   i,r;

    if ( align->trace_len < 1 )
        return SUCCESS;

    seq1 = (char*) malloc( (display_size+1) * sizeof(char) );
    MEM_ERROR( seq1 );

    seq2 = (char*) malloc( (display_size+1) * sizeof(char) );
    MEM_ERROR( seq2 );

    matchchar = (char*) malloc( (display_size+1) * sizeof(char) );
    MEM_ERROR( matchchar );

    tickmarks = (char*) malloc( (display_size+1) * sizeof(char) );
    MEM_ERROR( tickmarks );

    for(i=0; i<display_size; i++)
    {
        if ((i+1)%10 == 0) {
	    tickmarks[i]=':';
	} else if ((i+1)%5 == 0) {
	    tickmarks[i]='.';
	} else {
	    tickmarks[i]=' ';
	}
    }
    tickmarks[display_size] = '\0';

    for (i=0; i<align->trace_len; i++) {
	seq1[i%display_size]=align->trace_dchar[i];
	seq1[i%display_size+1]='\0';
	seq2[i%display_size]=align->trace_qchar[i];
	seq2[i%display_size+1]='\0';
	matchchar[i%display_size]=align->trace_mchar[i];
	matchchar[i%display_size+1]='\0';
        if ((i+1)%display_size == 0){
	    fprintf(fout, "         %s\n", tickmarks);
	    fprintf(fout, "%8d %s\n",
                    align->trace_dpos[i-display_size+1]+1, seq1);
	    fprintf(fout, "         %s\n", matchchar);
	    fprintf(fout, "%8d %s\n",
                    align->trace_qpos[i-display_size+1]+1, seq2);
	    fprintf(fout, "\n");
	}
    }

    /* If we have a partial line left, print it. */
    if (i%display_size != 0){
        tickmarks[i%display_size] = '\0';
	i=(i/display_size)*display_size;

	/* This will be < i, because i%display_size != 0. */
	fprintf(fout, "         %s\n", tickmarks);
	fprintf(fout, "%8d %s\n", align->trace_dpos[i]+1, seq1);
	fprintf(fout, "         %s\n", matchchar);
	fprintf(fout, "%8d %s\n", align->trace_qpos[i]+1, seq2);
	fprintf(fout, "\n");
    }

    r = SUCCESS;
    goto cleanup;

    error:
        r = ERROR;

    cleanup:
        FREE( seq1 );
        FREE( seq2 );
        FREE( matchchar );
	FREE( tickmarks );

    return r;
}


/*
 * This function releases the memory held by an Align structure.
 * Its synopsis is:
 *
 * result = align_release( align, message)
 *
 * where
 *      align           is the address of the Align structure
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 */
int 
align_release(Align * align, BtkMessage* message)
{
    FREE(align->trace_qpos);
    FREE(align->trace_dpos);
    FREE(align->trace_dir);
    FREE(align->trace_qchar);
    FREE(align->trace_mchar);
    FREE(align->trace_dchar);

    align->trace_max_len=0;
    align->trace_len=0;
    align->contig_offset=0;

    return SUCCESS;
}

