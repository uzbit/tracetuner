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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include "util.h"
#include "Btk_qv.h"
#include "Btk_qv_data.h"
#include "Btk_align_fragments.h"

#define PENALTYCOEF 100000

/**
 * This function releases the memory resource held by the specified GAlign 
 * struncture.  Its synopsis is:
 *
 * result = GAlign_release(align, message)
 * 
 * where
 *	align		is the address of the GAlign structure
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 */
int
GAlign_init(GAlign *align, BtkMessage *message)
{
    if (align == NULL) {
	return ERROR;
    }

    align->aind = NULL;
    align->rind = NULL;
    align->apos = NULL;
    align->rpos = NULL;
    align->achar = NULL;
    align->rchar = NULL;
    align->mchar = NULL;

    align->len = 0;
    align->max_len = 0;
    align->num_gaps = 0;

    return SUCCESS;
}

/**
 * This function releases the memory resource held by the specified GAlign 
 * struncture.  Its synopsis is:
 *
 * result = GAlign_release(align, message)
 * 
 * where
 *	align		is the address of the GAlign structure
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 */
int
GAlign_release(GAlign *align, BtkMessage *message)
{
    if (align == NULL) {
	return SUCCESS;
    }
    FREE(align->aind);
    FREE(align->rind);
    FREE(align->apos);
    FREE(align->rpos);
    FREE(align->achar);
    FREE(align->rchar);
    FREE(align->mchar);

    align->len = 0;
    align->max_len = 0;
    align->num_gaps = 0;

    return SUCCESS;
}

/*******************************************************************************
 * Function: p3
 * Purpose:  compute a match premium / mismatch penalty for characters
*******************************************************************************
 */
long    
p3(char abase, char rbase, int aweight)
{
    if (abase == rbase) {
        return (10 * PENALTYCOEF + (long) aweight);
    } else {
        return (-20 * PENALTYCOEF);
    }
}
 
/*******************************************************************************
 * Function: g3
 * Purpose:  compute gap penalty for a character
*******************************************************************************
 */
long   
g3(char base)
{
    return (-40 * PENALTYCOEF);
}
 

/*******************************************************************************
 * Function: populate_fragment_similarity_matrix
 * Purpose:  compute elements of a similarity matrix
 *******************************************************************************
 */
int
populate_fragment_similarity_matrix(Base *afrag, Base *rfrag, int m, int *aw, 
    int n, long **a, long (*g)(char), long (*p)(char, char, int),
    BtkMessage *message)
{
    int i, j;

    /* Populate 0-th row and 0-th column */
    a[0][0] = 0;
    for (i=1; i<=m; i++) {
        a[i][0] = a[i-1][0] + (*g)(afrag[i-1].base);
    }
    for (j=1; j<=n; j++) {
        a[0][j] = a[0][j-1] + (*g)(rfrag[j-1].base);
    }

    /* Populate the rest of the elements */
    for (i=1; i<=m; i++)
    {
        for (j=1; j<=n; j++) {
            a[i][j] = a[i-1][j-1]   + (*p)(afrag[i-1].base, 
					   rfrag[j-1].base,
					   aw[i-1]);
            if (a[i][j] < a[i-1][j] + (*g)(afrag[i-1].base))
                a[i][j] = a[i-1][j] + (*g)(afrag[i-1].base);
            if (a[i][j] < a[i][j-1] + (*g)(rfrag[j-1].base))
                a[i][j] = a[i][j-1] + (*g)(rfrag[j-1].base);
        }
    }


#if 0
    /* Output the similarity matrix */
    for (i=0; i<=m; i++)
    {
        printf("\n i==%d\n", i);
        for (j=0; j<=n; j++) {
            printf("%9d  ", a[i][j]);
            if (j%10 == 9) printf("\n");
        }
    }
#endif

    return SUCCESS;
}

/*******************************************************************************
* Function: itoc                               
* Purpose:  return integer in the range from 1 to 9 as a character symbol
*******************************************************************************
*/
char
itoc(int i)
{
    if      (i <= 1) return '1';
    else if (i == 2) return '2';
    else if (i == 3) return '3';
    else if (i == 4) return '4';
    else if (i == 5) return '5';
    else if (i == 6) return '6';
    else if (i == 7) return '7';
    else if (i == 8) return '8';
    
    return '9';
}


/*******************************************************************************
* Function: compute_optimal_fragment_alignment
* Purpose:  compute the alignment vectors corresponding to the the alignment
*           with highest score
*******************************************************************************
*/
int
compute_optimal_fragment_alignment(Base *afrag, Base *rfrag, int i, int* aw,
				   int j, long **a, long (*g)(char),
				   long (*p)(char, char, int), GAlign *align,
    				   BtkMessage *message)
{
    int len = 0;
    align->score = a[i][j];
    align->num_gaps = 0;
  
    align->aind[len] =  afrag[i-1].ind;  
    align->rind[len] =  rfrag[j-1].ind;   

    while (i > 0 || j > 0) {
#if 0
        fprintf(stderr,
            "In compute_optimal_fragment_alignment: i=%d, j=%d\n",
            i, j);
#endif
        if (i>0 && a[i][j] == a[i-1][j]+(*g)(afrag[i-1].base))
        {
            len++;
            align->achar[len] =  afrag[i-1].base;
            align->rchar[len] =  '-';
            align->mchar[len] =  ' ';
            align->aind[len] =  afrag[i-1].ind;
            align->rind[len] =  len>0 ? align->rind[len-1] : align->rind[0];
            i--;
            align->num_gaps++;
        }
        else if (i>0 && j >0 &&
                 a[i][j] == a[i-1][j-1]+(*p)(afrag[i-1].base, 
					     rfrag[j-1].base,
					     aw[i-1]))
        {
            len++;
            align->achar[len] =  afrag[i-1].base;
            align->rchar[len] =  rfrag[j-1].base;
            if (afrag[i-1].base==rfrag[j-1].base) {
                align->mchar[len] =  '|';
            }
            else {
                align->mchar[len] =  ' ';
            }
            align->aind[len]  =  afrag[i-1].ind;
            align->rind[len]  =  rfrag[j-1].ind;
            i--;
            j--;
        }
        else if (j >0 && a[i][j] == a[i][j-1]+(*g)(rfrag[j-1].base))
            /* has to be j >0 && a[i][j] == a[i][j-1] + (*g)(rfrag[j-1].base) */
        {
            len++;
            align->achar[len] =  '-';
            align->rchar[len] =  rfrag[j-1].base;
            align->mchar[len] =  ' ';
            align->aind[len]  =  len>0 ? align->aind[len-1] : align->aind[0];  
            align->rind[len]  =  rfrag[j-1].ind;
            j--;
            align->num_gaps++;
        }

        if (i==0 && j==0) {
            align->len = len;
            return SUCCESS;
        }
    }

    return SUCCESS;
}


/*******************************************************************************
 * Function: align_fragments 
 * Purpose:  align two DNA fragments   
 *******************************************************************************
 */
int
align_fragments(Base *afrag, int afrag_len, int *aweights, Base *rfrag, 
		int rfrag_len, long (*g)(char), long (*p)(char, char, int),
		GAlign *align, FILE* fout, BtkMessage *message)
{
    int      i, m, n, aind, rind;
    long **a;                     // similarity matrix
    int      sum_len, printed, to_print, rest_len;

    /* Allocate memory */
    m = afrag_len;
    n = rfrag_len;
    fprintf(stderr, "\nAfrag_len = %d, Rfrag_len = %d\n\n", afrag_len, rfrag_len);
    sum_len = m + n;

    a = CALLOC(long *, m+1);
    MEM_ERROR(a);
    /* m+1-th element corresponds to a gap */

    for (i=0; i<=m; i++) {
       a[i] = CALLOC(long, n+1);
    }
    align->achar = CALLOC(char, sum_len);
    MEM_ERROR(align->achar);
    align->mchar = CALLOC(char, sum_len);
    MEM_ERROR(align->mchar);
    align->rchar = CALLOC(char, sum_len);
    MEM_ERROR(align->rchar);
    align->aind  = CALLOC(int , sum_len);
    MEM_ERROR(align->aind);
    align->rind  = CALLOC(int , sum_len);
    MEM_ERROR(align->rind);

    fprintf(stderr, "    Populating fragment similarity matrix ...\n");
    if (populate_fragment_similarity_matrix(afrag, rfrag,
        m, aweights, n, a, g, p, message) != SUCCESS) {
        goto error;
    }
    fprintf(stderr, "    done...\n\n");

    fprintf(stderr, "    Computing_optimal_fragment_alignment ...\n");
    if (compute_optimal_fragment_alignment(afrag, rfrag,
        m, aweights, n, a, g, p, align, message) != SUCCESS) {
        goto error;
    }
    fprintf(stderr, "    done...\n\n");

    /* Output the alignment results */
    rest_len = align->len;  /* num of characters left to print */
    fprintf(stderr, "align->len=%d\n", align->len);
    fprintf(fout, "align->len=%d\n", align->len);
    to_print = QVMIN(70, align->len); /* num of characters to print in the current line */
    printed = 1;            /* num of characters already printed */
    aind = align->aind[printed];
    rind =  align->rind[printed];
    while (rest_len > 0) {
        fprintf(stderr, "%4d  ", aind);
        fprintf(fout, "%4d  ", aind);
        for (i=printed; i<printed+to_print; i++) {
            fprintf(stderr, "%c", align->achar[i]);
            fprintf(fout, "%c", align->achar[i]);
            aind = align->aind[i];
        }
        fprintf(stderr, "\n");
        fprintf(fout, "\n");
        fprintf(stderr, "      ");
        fprintf(fout, "      ");
        for (i=printed; i<printed+to_print; i++) {
            fprintf(stderr, "%c", align->mchar[i]);
            fprintf(fout, "%c", align->mchar[i]);
        }
        fprintf(stderr, "\n");
        fprintf(fout, "\n");
        fprintf(stderr, "%4d  ", rind);
        fprintf(fout, "%4d  ", rind);
        for (i=printed; i<printed+to_print; i++) {
            fprintf(stderr, "%c", align->rchar[i]);
            fprintf(fout, "%c", align->rchar[i]);
            rind = align->rind[i];
        }
        fprintf(stderr, "\n\n");
        fprintf(fout, "\n\n");
        rest_len -= to_print;
        printed += to_print;
        if (rest_len > 0) {
            to_print = rest_len > to_print ? to_print : rest_len;
        }
    }

    fprintf(stderr, "\nAlignment score=%9.3f\n",  align->score / PENALTYCOEF);
    fprintf(fout, "\nAlignment score=%9.3f\n",  align->score / PENALTYCOEF);
    fprintf(stderr, "Alignment length=%d\n",  align->len);
    fprintf(fout, "Alignment length=%d\n",  align->len);
    fprintf(stderr, "Num_gaps=%d\n", align->num_gaps);
    fprintf(fout, "Num_gaps=%d\n", align->num_gaps);

    /* Compute 5th dye mobility shift */

    for (i=0; i<=m; i++) {
       FREE(a[i]);
    }
    FREE(a);

    return SUCCESS;

error:
    FREE(align->achar);
    FREE(align->mchar);
    FREE(align->rchar);
    FREE(align->apos );
    FREE(align->rpos );
    FREE(align->aind );
    FREE(align->rind );
    for (i=0; i<=m; i++) {
       FREE(a[i]);
    }
    FREE(a);

    return ERROR;
}
