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
 * 2.5 2003/11/06 18:18:23
 */

/*
 *  Btk_align_peaks.c $Revision: 1.5 $
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "Btk_qv.h"
#include "util.h"
#include "Btk_align_peaks.h"
#include "Btk_process_raw_data.h"
#include "Btk_qv_data.h"

#define AREA_FACTOR1 0.0100
#define AREA_FACTOR2 0.0200
#define WIDTH_FACTOR 2.000
// No More (SLT 5/30/2001): 
/* these should end with ...0001 to avoid difference between platforms !!! */

#define HALF_WINDOW_SIZE 5
#define MAX_CONSENSUS_LEN 200000 
//#define PEAK_HEIGHT_FACTOR 0.0001
#define SHOW_MULTIPLE_PEAK_FIT 0
#define SHOW_MULTIPLE_RESOLUTION 0
#define SHOW_PEAK_STATS 0
#define WEIGHT 0.

/* Primitive data structure for peak (used only in the current file)
 */
typedef struct {
    double area;            /* peak area */
    int    beg;             /* position of the left boundary of a peak */
    int    end;             /* position of the right boundary of a peak */
    int    height;          /* peak's hight (=signal at the peak's position) */
    int    pos;             /* peak's position */
    int    ind;             /* peak's index in the peak list */
    int    is_truncated;
    double relative_area;
    double relative_height;
    double relative_spacing;
} PeakA;

/* Data structure for base (used only in the current file)
 */
typedef struct {
    char   base;
    int    pos;             /* base position */
    int    ind;             /* base index in a fragment */
} Base;

/* An Align structure contains information returned from an alignment:
 */
    typedef struct
    {
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
    } Align;

/*******************************************************************************
 * Function: p           
 * Purpose:  compute a match premium / mismatch penalty
*******************************************************************************
 */
double
p1(PeakA apeak, PeakA rpeak)
{
    double height, spacing;

    if (apeak.relative_height >= rpeak.relative_height)
        height = apeak.relative_height - rpeak.relative_height;
    else 
        height = rpeak.relative_height - apeak.relative_height;

    if (apeak.relative_spacing>= rpeak.relative_spacing)
        spacing= apeak.relative_spacing- rpeak.relative_spacing;
    else
        spacing= rpeak.relative_spacing- apeak.relative_spacing;

    return (-height - WEIGHT * spacing);

#if 0
    return (QVMIN(apeak.relative_height,  rpeak.relative_height) -
            QVMAX(apeak.relative_height,  rpeak.relative_height) 
                                                                   +
            QVMIN(apeak.relative_spacing, rpeak.relative_spacing)-
            QVMAX(apeak.relative_spacing, rpeak.relative_spacing));
#endif
}

/*******************************************************************************
 * Function: p2
 * Purpose:  compute a match premium / mismatch penalty
*******************************************************************************
 */
double
p2(PeakA apeak, PeakA bpeak)
{
    double spacing;

    if (apeak.relative_spacing>= bpeak.relative_spacing)
        spacing = bpeak.relative_spacing > 0 ? 
            apeak.relative_spacing/bpeak.relative_spacing : 1;
    else
        spacing = apeak.relative_spacing > 0 ?
            bpeak.relative_spacing/apeak.relative_spacing : 1;

    return -spacing;
}
 
/*******************************************************************************
 * Function: p3
 * Purpose:  compute a match premium / mismatch penalty for characters
*******************************************************************************
 */
int    
p3(char abase, char rbase)
{
    if (abase == rbase)
        return 1;
    else 
        return (-1);
}
 
/*******************************************************************************
 * Function: g1           
 * Purpose:  compute gap penalty for a peak using relative_height
*******************************************************************************
 */
double
g1(PeakA peak)
{
    return (-peak.relative_height);
}    

/*******************************************************************************
 * Function: g2
 * Purpose:  compute gap penalty for a peak using relative_spacing
*******************************************************************************
 */
double
g2(PeakA peak)
{
    if (peak.relative_spacing > 1) 
        return (-peak.relative_spacing);
    else 
        return (-1./peak.relative_spacing);
}

/*******************************************************************************
 * Function: g3
 * Purpose:  compute gap penalty for a character
*******************************************************************************
 */
int   
g3(char base)
{
    return (-2);
}
 
/*******************************************************************************
 * Function: populate_peak_list_similarity_matrix     
 * Purpose:  compute elements of a similarity matrix 
 *******************************************************************************
 */
int
populate_peak_list_similarity_matrix(PeakA *alist, PeakA *rlist, int m,
    int n, double **a, double (*g)(PeakA), double (*p)(PeakA, PeakA),
    BtkMessage *message)
{
    int i, j;
    
    /* Populate 0-th row and 0-th column */
    a[0][0] = 0.;
    for (i=1; i<=m; i++) {
        a[i][0] = a[i-1][0] + (*g)(alist[i-1]);
#if 0
            if (g(alist[i-1]) > 0)
                printf("g(%d)==%5.3f\n",i-1, (*g)(alist[i-1]));
#endif
    }
    for (j=1; j<=n; j++) {
        a[0][j] = a[0][j-1] + (*g)(rlist[j-1]);
#if 0
            if (g(alist[j-1]) > 0)
                printf("g(%d)==%5.3f\n",j-1, g(rlist[j-1]));
#endif
    }
     
    /* Populate the rest of the elements */ 
#if 0
    fprintf(stderr, "In populate ...: m=%d, n=%d\n", m, n);
#endif
    
    for (i=1; i<=m; i++)
    {
        for (j=1; j<=n; j++) {
            a[i][j] = a[i-1][j-1]   + (*p)(alist[i-1], rlist[j-1]);
            if (a[i][j] < a[i-1][j] + (*g)(alist[i-1]))
                a[i][j] = a[i-1][j] + (*g)(alist[i-1]);
            if (a[i][j] < a[i][j-1] + (*g)(rlist[j-1]))
                a[i][j] = a[i][j-1] + (*g)(rlist[j-1]);
//          a[i][j] = QVMAX(a[i][j], a[i-1][j] + g(alist[i-1]));
//          a[i][j] = QVMAX(a[i][j], a[i][j-1] + g(rlist[j-1]));
#if 0
            if (p(alist[i-1], rlist[j-1]) > 0)
                printf("p(%d,%d)==%15.8f  heights=%5.3f %5.3f spacings=%5.3f %5.3f\n", 
                    i-1,j-1, p(alist[i-1], rlist[j-1]),
                    alist[i-1].relative_height, rlist[j-1].relative_height,
                    alist[i-1].relative_spacing, rlist[j-1].relative_spacing);
#endif
        }
    }

    
#if 0
    /* Output the similarity matrix */
    for (i=0; i<=m; i++)
    {
        printf("\n i==%d\n", i);
        for (j=0; j<=n; j++) {
            printf("%9.3f  ", a[i][j]);
            if (j%10 == 9) printf("\n");
        }
    }
#endif
    return SUCCESS;
}

/*******************************************************************************
 * Function: populate_fragment_similarity_matrix
 * Purpose:  compute elements of a similarity matrix
 *******************************************************************************
 */
int
populate_fragment_similarity_matrix(Base *afrag, Base *rfrag, int m,
    int n, int **a, int (*g)(char), int (*p)(char, char),
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
            a[i][j] = a[i-1][j-1]   + (*p)(afrag[i-1].base, rfrag[j-1].base);
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
* Function: compute_optimal_peak_list_alignment
* Purpose:  compute the alignment vectors corresponding to the the alignment
*           with highest score 
*******************************************************************************
*/
int
compute_optimal_peak_list_alignment(PeakA *alist, PeakA *rlist, int i, int j,
    double **a, double (*g)(PeakA), double (*p)(PeakA, PeakA), Align *align, 
    double max_height, double max_spacing, BtkMessage *message)
{
    int len = 0;
    double eps=0.00000001;
    align->score = a[i][j];

    align->num_gaps = 0;
    while (i > 0 || j > 0) {
        if (i>0 && fabs(a[i][j] -a[i-1][j]-(*g)(alist[i-1])) < eps) 
        {
            len++;
            if (max_height > 0) {
                align->achar[len] = 
                    itoc(INT_DBL(alist[i-1].relative_height/max_height * 10.));
            }
            if (max_spacing > 0) {
                align->achar[len] = 
                   itoc(INT_DBL(alist[i-1].relative_spacing/max_spacing * 10.));
            }
            align->rchar[len] =  '0';
            align->mchar[len] =  ' ';           
            align->apos[len] =  alist[i-1].pos;
            align->rpos[len] =  -1;
            align->aind[len] =  alist[i-1].ind;
            align->rind[len] =  -1;          
            i--;
            align->num_gaps++;
        }
        else if (i>0 && j >0 && 
                 fabs(a[i][j] - a[i-1][j-1]-(*p)(alist[i-1], rlist[j-1])) < eps) 
        {
            len++;
            if (max_height > 0) {
                align->achar[len] = 
                    itoc(INT_DBL(alist[i-1].relative_height/max_height * 10.));
                align->rchar[len] = 
                    itoc(INT_DBL(rlist[j-1].relative_height/max_height * 10.));
            }
            if (max_spacing > 0) {
                align->achar[len] = 
                   itoc(INT_DBL(alist[i-1].relative_spacing/max_spacing * 10.));
                align->rchar[len] = 
                   itoc(INT_DBL(rlist[j-1].relative_spacing/max_spacing * 10.));
            }
            align->mchar[len] =  '|';
            align->apos[len] =  alist[i-1].pos;
            align->rpos[len] =  rlist[j-1].pos;
            align->aind[len] =  alist[i-1].ind;
            align->rind[len] =  rlist[j-1].ind;
            i--;
            j--;
        }
        else if (j >0 && fabs(a[i][j] - a[i][j-1]-(*g)(rlist[j-1])) < eps)
            /* has to be j >0 && a[i][j] == a[i][j-1] + (*g)(rlist[j-1]) */
        {
            len++;
            align->achar[len] =  '0';
            if (max_height > 0) {
                align->rchar[len] = 
                    itoc(INT_DBL(rlist[j-1].relative_height/max_height * 10.));
            }
            if (max_spacing > 0) {
                align->rchar[len] = 
                   itoc(INT_DBL(rlist[j-1].relative_spacing/max_spacing * 10.));
            }
            align->mchar[len] =  ' ';
            align->apos[len] =  -1;
            align->rpos[len] =  rlist[j-1].pos;
            align->aind[len] =  -1;
            align->rind[len] =  rlist[j-1].ind;
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
* Function: compute_optimal_fragment_alignment
* Purpose:  compute the alignment vectors corresponding to the the alignment
*           with highest score
*******************************************************************************
*/
int
compute_optimal_fragment_alignment(Base *afrag, Base *rfrag, int i, int j,
    int **a, int (*g)(char), int (*p)(char, char), Align *align,
    BtkMessage *message)
{
    int len = 0;
    align->score = a[i][j];
    align->num_gaps = 0;
  
    align->aind[len] =  afrag[i-1].ind;  
    align->rind[len] =  afrag[j-1].ind;   

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
                 a[i][j] == a[i-1][j-1]+(*p)(afrag[i-1].base, rfrag[j-1].base))
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
 * Function: align_peak_lists                
 * Purpose:  align the two lists of peaks;
 *******************************************************************************
 */
int
align_peak_lists(PeakA *alist, int alist_len, PeakA *rlist, int rlist_len,
    double (*g)(PeakA), double (*p)(PeakA, PeakA), Align *align, 
    double max_height, double max_spacing, BtkMessage *message)
{
    int      i, m, n, aind, rind;
    double **a;                     // similarity matrix
    int      sum_len, printed, to_print, rest_len;

    /* Allocate memory */ 
    m = alist_len;
    n = rlist_len;
    sum_len = m + n;

    a = CALLOC(double *, m+1);
    /* m+1-th element corresponds to a gap */

    for (i=0; i<=m; i++) {
       a[i] = CALLOC(double, n+1); 
    }
    align->achar = CALLOC(char, sum_len);
    align->mchar = CALLOC(char, sum_len);
    align->rchar = CALLOC(char, sum_len);
    align->apos  = CALLOC(int , sum_len);
    align->rpos  = CALLOC(int , sum_len);
    align->aind  = CALLOC(int , sum_len);
    align->rind  = CALLOC(int , sum_len);

    if (populate_peak_list_similarity_matrix(alist, rlist, 
        m, n, a, g, p, message) != SUCCESS) {
        goto error;  
    }

    if (compute_optimal_peak_list_alignment(alist, rlist, 
        m, n, a, g, p, align, max_height, max_spacing,
        message) != SUCCESS) {
        goto error;  
    }

    /* Output the alignment results */
    rest_len = align->len;  /* num of characters left to print */
    to_print = 50;          /* num of characters to print in the current line */
    printed = 1;            /* num of characters already printed */
    aind = align->aind[printed];
    rind =  align->rind[printed];
    while (rest_len > 0) {
        fprintf(stderr, "%4d  ", aind);
        for (i=printed; i<printed+to_print; i++) {
            fprintf(stderr, "%c", align->achar[i]);
            if (align->aind[i] >=0) 
                aind = align->aind[i];
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "      ");
        for (i=printed; i<printed+to_print; i++) {
            fprintf(stderr, "%c", align->mchar[i]);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "%4d  ", rind);
        for (i=printed; i<printed+to_print; i++) {
            fprintf(stderr, "%c", align->rchar[i]);
            if (align->rind[i] >=0)
                rind = align->rind[i];
        }
        fprintf(stderr, "\n\n");
        rest_len -= to_print;
        printed += to_print;
        if (rest_len > 0) {
            to_print = rest_len > to_print ? to_print : rest_len;
        }
    }

    fprintf(stderr, "\nAlignment score=%9.3f\n",  align->score); 
    fprintf(stderr, "Alignment length=%d\n",  align->len);
    fprintf(stderr, "Num_gaps=%d\n", align->num_gaps);
   
    /* Compute 5th dye mobility shift */

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

/*******************************************************************************
 * Function: align_fragments 
 * Purpose:  align two DNA fragments   
 *******************************************************************************
 */
int
align_fragments(Base *afrag, int afrag_len, Base *rfrag, int rfrag_len,
    int (*g)(char), int (*p)(char, char), Align *align, BtkMessage *message)
{
    int   i, m, n, aind, rind;
    int **a;                     // similarity matrix
    int   sum_len, printed, to_print, rest_len;

    /* Allocate memory */
    m = afrag_len;
    n = rfrag_len;
    fprintf(stderr, "\nAfrag_len = %d, Rfrag_len = %d\n\n", afrag_len, rfrag_len);
    sum_len = m + n;

    a = CALLOC(int *, m+1);
    /* m+1-th element corresponds to a gap */

    for (i=0; i<=m; i++) {
       a[i] = CALLOC(int, n+1);
    }
    align->achar = CALLOC(char, sum_len);
    align->mchar = CALLOC(char, sum_len);
    align->rchar = CALLOC(char, sum_len);
    align->aind  = CALLOC(int , sum_len);
    align->rind  = CALLOC(int , sum_len);

    fprintf(stderr, "    Populating fragment similarity matrix ...\n");
    if (populate_fragment_similarity_matrix(afrag, rfrag,
        m, n, a, g, p, message) != SUCCESS) {
        goto error;
    }
    fprintf(stderr, "    done...\n\n");

    fprintf(stderr, "    Computing_optimal_fragment_alignment ...\n");
    if (compute_optimal_fragment_alignment(afrag, rfrag,
        m, n, a, g, p, align, message) != SUCCESS) {
        goto error;
    }
    fprintf(stderr, "    done...\n\n");

    /* Output the alignment results */
    rest_len = align->len;  /* num of characters left to print */
    to_print = 70;          /* num of characters to print in the current line */
    printed = 1;            /* num of characters already printed */
    aind = align->aind[printed];
    rind =  align->rind[printed];
    while (rest_len > 0) {
        fprintf(stderr, "%4d  ", aind);
        for (i=printed; i<printed+to_print; i++) {
            fprintf(stderr, "%c", align->achar[i]);
            aind = align->aind[i];
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "      ");
        for (i=printed; i<printed+to_print; i++) {
            fprintf(stderr, "%c", align->mchar[i]);
        }
        fprintf(stderr, "\n");
        fprintf(stderr, "%4d  ", rind);
        for (i=printed; i<printed+to_print; i++) {
            fprintf(stderr, "%c", align->rchar[i]);
            rind = align->rind[i];
        }
        fprintf(stderr, "\n\n");
        rest_len -= to_print;
        printed += to_print;
        if (rest_len > 0) {
            to_print = rest_len > to_print ? to_print : rest_len;
        }
    }

    fprintf(stderr, "\nAlignment score=%9.3f\n",  align->score);
    fprintf(stderr, "Alignment length=%d\n",  align->len);
    fprintf(stderr, "Num_gaps=%d\n", align->num_gaps);

    /* Compute 5th dye mobility shift */

    FREE(align->achar);
    FREE(align->mchar);
    FREE(align->rchar);
    FREE(align->aind );
    FREE(align->rind );
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

/*******************************************************************************
 * Function: get_peak_area
 * Purpose: Calculate the area under a peak.  The algorithm is simple:
 *          the area is defined to be the total of the averages of the heights
 *          between each pair of points that make up the peak.
 *******************************************************************************
 */
static double
get_peak_area(int *data, int peak_beg, int peak_end,
                BtkMessage *message)
{
    int i, area;
    double peak_area;

    if (peak_beg <0 || peak_end < 0 || peak_end < peak_beg) {
        sprintf(message->text,
        "Incorrect peak boundaries (%d, %d) passed to get_peak_area\n",
        peak_beg, peak_end);
        return ERROR;
    }
       
    if (peak_beg == peak_end) {
        return (0.0);
    }

    if ((peak_end - peak_beg) == 1) {
        return ((double)(data[peak_end] + data[peak_beg]) / 2.0);
    }

    area = data[peak_beg] + data[peak_end];
    for (i = peak_beg + 1; i < peak_end; i++) {
        area += 2 * data[i];
    }
    peak_area = (double)area / 2.0;
    return (peak_area);
}

/*******************************************************************************
 * Function: get_peak_position
 * Purpose: Find location of observed peak as the position which bisects
 *          its area. If the position is found, return it. Otherwise, return -1.
 *******************************************************************************
 */
static int
get_peak_position(int *data, int peak_beg, int peak_end,
                            double peak_area, BtkMessage *message)
{
    int j;
    double area, old_area;

    /* Check the easy (degenerate) case first */
    if (peak_end - peak_beg <= 1) {
        return peak_end;
    }

    area = 0.0;
    old_area = 0.0;
    if ((peak_end - peak_beg) == 2) {
       return peak_beg+1;
    }
    else {
        for (j = peak_beg + 1; j <= peak_end; j++) {
            area += (double)(data[j] + data[j-1]) / 2.0;
            if (area == peak_area / 2.0) {
                return (j);
            }
            else if ((area > peak_area / 2.0) && (old_area <= peak_area / 2.0)) {
                if ((area - peak_area / 2.0) < (peak_area / 2.0 - old_area)) {
                    return (j);
                }
                else {
                    return (j - 1);
                }
            }
            old_area = area;
        }
    }

    (void)sprintf(message->text,
        "couldn't determine position of peak @ %d-%d", peak_beg, peak_end);
    return ERROR;
}

/*******************************************************************************
 * Function: is_big_peak
 * Purpose:  check if peak is "big" in David Ho's definition
 *******************************************************************************
 */
static int
is_big_peak(PeakA peak, double ave_width, double ave_area, int max_value)
{
#if 1
    if      (INT_GT_DBL( (peak.end-peak.beg), (ave_width * WIDTH_FACTOR) ) &&
                      peak.is_truncated)
    {
        return 1;
    }
    else if (INT_GT_DBL( (peak.end-peak.beg), (ave_width * WIDTH_FACTOR) ) &&
                      peak.area          > ave_area  / AREA_FACTOR1 &&
                      peak.area          > ave_width * max_value)
    {
        return 1;
    }
#else
    if      ((double)(peak.end-peak.beg) > ave_width * WIDTH_FACTOR &&
                      peak.is_truncated)
    {
        return 1;
    }
    else if ((double)(peak.end-peak.beg) > ave_width * WIDTH_FACTOR &&
                      peak.area          > ave_area  / AREA_FACTOR1 &&
                      peak.area          > ave_width * max_value)
    {
        return 1;
    }
#endif
    return 0;
}

/*******************************************************************************
 * Function: create_peak_list
 * Purpose:  create a "primitive" list of peaks for a given trace using
 *           the inflection points; include into the list all peaks
 *           of width > 1
 *******************************************************************************
 */
static int
create_peak_list(int *peak_list_len, int *peak_list_max_len, 
    PeakA *peak_list, int *data, int length, double *max_relative_height, 
    double *max_relative_spacing, BtkMessage *message)
{
    int i, j, k=0, il, ir, lind, rind; 
    int derivative_2, startscan, endscan, max_value=0;
    PeakA peak, temp_peak;
    double ave_peak_area, start_avg_peak_area, prev_peak_area, ave_height;
    double avg_peak_width, start_avg_peak_width, sum_peak_area, sum_peak_width;

  (*peak_list_len) = 0;

    /* Determine the maximum value of signal */
    for (i=0; i<length; i++) {
        if (max_value < data[i]) 
            max_value = data[i];
    }

    /*  Pre scan some data to get an average peak area to start;
     * Start prescanning from position of base number 120 or,
     * if total  number of bases is < 200, from the middle base.
     * Continue prescanning through the position of the last base
     * or untill 20 "good" peaks are found, whichever comes first
     */
    j = 0;
    temp_peak.beg = 0;
    temp_peak.end = 0;
    start_avg_peak_width = 0.0;
    start_avg_peak_area = 0.0;

    startscan = 1000;
    endscan   = 2000;
    if (endscan > length-2) {
        startscan = length/4;
        endscan   = length/2;
    }

    for (i = startscan; (i < endscan) && (j < 20); i++) {
        derivative_2 = data[i-1] - 2*data[i] + data[i+1];

        /* Temp peak begins */
        if ((derivative_2 < 0) && !temp_peak.beg) {
            temp_peak.beg = i;
        }

        /* Temp peak ends */
        if (temp_peak.beg &&
        ((derivative_2 > 0 && data[i]-2*data[i+1]+data[i+2] >= 0)
             || data[i] <=0))
        {
            temp_peak.end = i;
            if (
                /* Ignore oscillations of width 1 and narrow width!!! */
                 temp_peak.end-temp_peak.beg > 2 )
            {
                start_avg_peak_width += (temp_peak.end-temp_peak.beg);
                start_avg_peak_area += get_peak_area(data,temp_peak.beg,
                                                      temp_peak.end,message);
                j += 1;
            }
            temp_peak.beg = 0;
            temp_peak.end = 0;
        }
    }

    /* Initialize starting params */
    start_avg_peak_width /= j;
    start_avg_peak_area  /= j;
    prev_peak_area = start_avg_peak_area;

    /*
     * Use same approach as Phred: first, find a couple of subsequent
     * inflection points at one of which the 2nd derivative changes sign from
     * + to - and at the other backwards; then, localize the peak position
     * between the two points
     */
    for (i = 1; i < length - 2; i++) {
        derivative_2 = data[i-1] - 2*data[i] + data[i+1];

        /* Peak begins */
        if ((derivative_2 < 0) && !peak.beg) {
            peak.beg = i;
        }

        if (peak.beg && (data[i] >= max_value)) {
            peak.is_truncated = 1;
        }

        /* Peak ends */
        if (peak.beg &&
         ((derivative_2 > 0 && data[i]-2*data[i+1]+data[i+2] >= 0)
             || data[i] <=0))
        {
            peak.end = i;

            /* Determine the peak's position between peak.beg and peak.end */
            peak.area = get_peak_area(data, peak.beg, peak.end, message);
            peak.pos = get_peak_position(data, peak.beg, peak.end,
                peak.area, message);
            if (peak.pos == ERROR) {
                sprintf(message->text, "Error in get_peak_position\n");
                return ERROR;
            }

            /*
             * Record the peak only if its area is at least 10% of the average
             * area of 10 preceeding peaks and 5% of the area of the
             * immediately preceeding peak (of the same color).
             * When checking these criteria, do not count "too big" peaks".
             * To this end, calculate the ave_peak_area and ave_peak_width among
             * previous 10 peaks, not including the "big" ones
             */
            k = 0;
            if (*peak_list_len < 11) {
                ave_peak_area  = start_avg_peak_area;
                avg_peak_width = start_avg_peak_width;
            }
            else {
                sum_peak_area  = 0.0;
                sum_peak_width = 0.0;
                for (j = *peak_list_len - 1; k < 10 && j > 0; j--)
                {
                    /* check for big peaks */
                    if (!is_big_peak(peak_list[j], avg_peak_width,
                        ave_peak_area, max_value))
                    {
                        ++k;
                        sum_peak_area  += peak_list[j].area;
                        sum_peak_width += (peak_list[j].end -
                                          peak_list[j].beg);
                    }
                }
                if (k < 10) {
                   ave_peak_area  = start_avg_peak_area;
                   avg_peak_width = start_avg_peak_width;
                }
                else {
                   ave_peak_area  = sum_peak_area  / 10.0;
                   avg_peak_width = sum_peak_width / 10.0;
                }
            }

            if (
                /* Ignore oscillations of width 1 !!! */
                (peak.end-peak.beg > 2 ||
                 (peak.end-peak.beg==2 && data[peak.end-1] > 0)) &&

            /* Ignore peaks of too small area, see Phred paper II, p.178 */
                  peak.area >= AREA_FACTOR1 * prev_peak_area &&
                  peak.area >= AREA_FACTOR2 * ave_peak_area
               )
            {
                peak.height = data[peak.pos],
                peak.ind    = *peak_list_len; 
                if (ave_peak_area  > 0.) {
                    peak.relative_area = peak.area/ave_peak_area ;
                }
                else {
                    peak.relative_area = 1.;
                }
                if (*peak_list_len >= *peak_list_max_len) {
                  (*peak_list_max_len) *= 2;
                    peak_list = REALLOC(peak_list, PeakA, *peak_list_max_len);
                    MEM_ERROR(peak_list);
                }
                peak_list[*peak_list_len] = peak;
              (*peak_list_len)++;

                /* Update the value of previous_peak_area parameter for next
                 * round only if the current peak is not "too big" peak;
                 * that is, only not very big peks will be used in the
                 * criiteria for checking peak into peak list
                 */
                if (!is_big_peak(peak, avg_peak_width, ave_peak_area,
                    max_value))
                    prev_peak_area = peak.area;
            }

            /* Prepare for processing next peak */
            peak.beg = 0;
            peak.end = 0;
            peak.is_truncated = 0;

        } /* end processing current peak */
    } /* end loop in i */

    /* Determine relative height and spacing for each peak in the list */
    for (i=0; i<*peak_list_len; i++) {

        /* Find the left boundary of the window */
        il = i;
        j  = 0;
        while (il>0 && j<HALF_WINDOW_SIZE) {
            il--;
            j++;
        }
        
        /* Find the right boundary of the window */
        ir = i;
        j  = 0;
        while (ir<*peak_list_len-1 && j<HALF_WINDOW_SIZE) {
            ir--;
            j++;
        }

        /* Compute the average peak height in the window */
        ave_height = 0;
        k = 0;
        for (j = il; j<=ir; j++) {
            if (!peak_list[j].is_truncated) {
                k++;
                ave_height += peak_list[j].height;
            } 
        }
        ave_height /= (double)k;

        if (i==il) {
            peak_list[i].relative_spacing = 1.;

            /* Find the high closest right neighbour of the current peak */
            rind  = i+1;
            j = i+1;
            while (j < ir) {
                if (peak_list[j].height > ave_height) {
                    rind  = j;
                    break;
                }
                j++;
            }
            peak_list[i].relative_height =
                (double)peak_list[i].height / (double)peak_list[rind ].height;
        }
        else if (i==*peak_list_len-1) {
            peak_list[i].relative_spacing = 1.;

            /* Find the high closest left neighbour of the current peak */
            lind = i-1;
            j = i-1;
            while (j > il) {
                if (peak_list[j].height > ave_height) {
                    lind  = j;
                    break;
                }
                j--;
            }
            peak_list[i].relative_height =
                (double)peak_list[i].height / (double)peak_list[lind].height;
        }
        else {

             /* Find the high closest right neighbour of the current peak */
            rind  = i+1;
            j = i+1;
            while (j < ir) {
                if (peak_list[j].height > ave_height) {
                    rind  = j;
                    break;
                }
                j++;
            }

            /* Find the high closest left neighbour of the current peak */
            lind = i-1;
            j = i-1;
            while (j > il) {
                if (peak_list[j].height > ave_height) {
                    lind  = j;
                    break;
                }
                j--;
            }
            peak_list[i].relative_spacing =
                (double)(peak_list[i    ].pos - peak_list[lind].pos) /
                (double)(peak_list[rind ].pos - peak_list[i   ].pos);
            peak_list[i].relative_height =
                (double)peak_list[i].height * 2. /
               ((double)peak_list[lind].height+(double)peak_list[rind ].height);
        }
        if (peak_list[i].relative_height  > (*max_relative_height))
           *max_relative_height = peak_list[i].relative_height;
        if (peak_list[i].relative_spacing > (*max_relative_spacing))
           *max_relative_spacing = peak_list[i].relative_spacing;
    }
    return SUCCESS;

error:
    return ERROR;
}

/*******************************************************************************
 * Function: get_vector_fragment
 *
 * Purpose:  read the pGEM reverse vector sequence
 *******************************************************************************
 */
int
get_vector_fragment(char **vector, int *vect_len)
{
    int  pGEMrev_len;
    char pGEMrev[] = 
                    "ATGATTACGCCAAGCTATTTAGGTGACACTATAGAATACTCAAGCTTGCATGCCTGCAG\
GTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTCGCCCTATAGTGAGTCGTATTACAATTCACTGGCCGTCGTTT\
TACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGT\
AATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGACGCGCCCTGTAGCGGC\
GCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGCTACACTTGCCAGCGCCCTAGCGCCCGCTCCTTTCGC\
TTTCTTCCCTTCCTTTCTCGCCACGTTCGCCGGCTTTCCCCGTCAAGCTCTAAATCGGGGGCTCCCTTTAGGGTTCCGAT\
TTAGTGCTTTACGGCACCTCGACCCCAAAAAACTTGATTAGGGTGATGGTTCACGTAGTGGGCCATCGCCCTGATAGACG\
GTTTTTCGCCCTTTGACGTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCAAACTGGAACAACACTCAACCCTAT\
CTCGGTCTATTCTTTTGATTTATAAGGGATTTTGCCGATTTCGGCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAAT\
TTAACGCGAATTTTAACAAAATATTAACGCTTACAATTTCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTT\
CACACCGCATATGGTGCACTCTCAGTACAATCTGCTCTGATGCCGCATAGTTAAGCCAGCCCCGACACCCGCCAACACCC\
GCTGACGCGCCCTGACGGGCTTGTCTGCTCCCGGCATCCGCTTACAGACAAGCTGTGACCGTCTCCGGGAGCTGCATGTG\
TCAGAGGTTTTCACCGTCATCACCGAAACGCGCGAGACGAAAGGGCCTCGTGATACGCCTATTTTTATAGGTTAATGTCA\
TGATAATAATGGTTTCTTAGACGTCAGGTG";

    pGEMrev_len = strlen(pGEMrev);
   *vector = CALLOC(char , pGEMrev_len+1);
    strcpy(*vector, pGEMrev); 
   *vect_len = pGEMrev_len;
 
    return SUCCESS;
}
/*******************************************************************************
 * Function: Btk_align_peaks      
 * Purpose: - create list of all peaks from analyzed data;
 *          - reanalyze raw data;
 *          - create list of all peaks from reanalyzed data;
 *          - align the two lists of peaks, find correspondence between the
 *            peak positions;
 *          - interpolate the map between peak positions and compute
 *            new base locations;
 *******************************************************************************
 */
int
Btk_align_peaks(
    int   **chromatogram,  
    int     num_points,
    int   **rchromatogram, 
    long    rnum_points[NUM_COLORS],
    int     js,             /* index of "selected" trace */
    int     num_called_bases, 
    char   *called_bases, 
    char   *color2base, 
    BtkMessage *message) 
{
    int   i, j, alist_len, alist_max_len, rlist_len, rlist_max_len;
    int   vector_len, blist_len;
    PeakA *alist, *rlist, *blist;
    char *vector;
    Base *vectfrag, *callfrag;
    Align align1, align2, align3;
    double max_relative_height=0., max_relative_spacing=0.;

    fprintf(stderr, "Color2base=%s js=%d Selected base = %c\n\n", 
        color2base, js, color2base[js]);

    /* Count the number of "selected" bases in the vector sequence */
    get_vector_fragment(&vector, &vector_len);
    blist_len = 0;
    for (i=0; i<vector_len; i++) {
        if (vector[i] == color2base[js])
            blist_len++;
    }

    fprintf(stderr, "rnum_points = %ld %ld %ld %ld \n",
        rnum_points[0], rnum_points[1], rnum_points[2], rnum_points[3]); 
    if (process_raw_data(rchromatogram, rnum_points) != SUCCESS) {
        goto error;
    }

    /* Allocate memory for peak lists */
    alist_max_len = rlist_max_len = MAX_NUM_OF_PEAK;
    alist = CALLOC(PeakA, alist_max_len);
    rlist = CALLOC(PeakA, rlist_max_len);
    blist = CALLOC(PeakA, blist_len);

    /* Create peak list for the analyzed "selected" trace */
    if (create_peak_list(&alist_len, &alist_max_len, alist, 
        chromatogram[js], num_points, &max_relative_height,
        &max_relative_spacing, message) != SUCCESS) {
        goto error;
    }

    /* Create peak list for the reanalyzed 5th dye trace */
    if (create_peak_list(&rlist_len, &rlist_max_len, rlist, 
        rchromatogram[js], (int)rnum_points[js], 
        &max_relative_height, &max_relative_spacing, message) != SUCCESS) {
        goto error;
    }
    fprintf(stderr, "\n\nmax_relative_height=%f, max_relative_spacing=%f\n\n",
        max_relative_height, max_relative_spacing);

    fprintf(stderr, "Aligning analyzed and reanalyzed peak lists ...\n\n");
    fprintf(stderr, "Alist_len = %d, Rlist_len = %d\n\n", alist_len, rlist_len);
    if (align_peak_lists(alist, alist_len, rlist, rlist_len,
        g1, p1, &align1, max_relative_height/4., -1., message) != SUCCESS) {
        goto error;
    }
    fprintf(stderr, "done ...\n\n");

    /* Initialize array 'blist' which corresponds to "selected"
     * bases in the vector 
     */ 
    j=0;
    for (i=0; i<vector_len; i++) {
        if (vector[i]==color2base[js]) {
            blist[j].pos = 12 * i;
            blist[j].ind = j;
            j++;
        }
    }
    for (j=0; j<blist_len; j++) {
        if (j==0 || j==blist_len-1) {
            blist[j].relative_spacing = 1.;
        }
        else {
            blist[j].relative_spacing =
                (double)(blist[j  ].pos-blist[j-1].pos)/
                (double)(blist[j+1].pos-blist[j  ].pos);
        }
    }

    /* Align array of "selected" bases against array of analyzed peaks */
    fprintf(stderr,
        "Aligning array of '%c's against array of analyzed peaks ...\n\n",
        color2base[js]);
    fprintf(stderr, "Alist_len = %d, Blist_len = %d\n\n", alist_len, blist_len);
    if (align_peak_lists(alist, alist_len, blist, blist_len,
            g2, p2, &align2, -1., max_relative_spacing/2., message) != SUCCESS) {
            goto error;
    }
    fprintf(stderr, "done ...\n\n");

    /* Align called fragment against correct sequence */
    fprintf(stderr,
        "Aligning called fragment against correct sequence ...\n\n");
#if 1
    fprintf(stderr,"pGEMrev=\n");
    for (i=0; i<vector_len; i++) {
        fprintf(stderr,"%c", vector[i]);
    }
    fprintf(stderr,"\n");
#endif
    vectfrag = CALLOC(Base, vector_len);
    callfrag = CALLOC(Base, num_called_bases);
    for (i=0; i<vector_len; i++) {
        vectfrag[i].base = vector[i];
        vectfrag[i].ind  = i;
    }
    for (i=0; i<num_called_bases; i++) {
        callfrag[i].base = called_bases[i];
        callfrag[i].ind  = i;
    }
    if (align_fragments(vectfrag, vector_len, callfrag, num_called_bases,
        g3, p3, &align3, message) != SUCCESS) {
            goto error;
    }
    fprintf(stderr,"done ...\n\n");

    FREE(alist);
    FREE(rlist);
    FREE(blist);
    FREE(vector);
    FREE(vectfrag);
    FREE(callfrag);
    return SUCCESS;   

error:
    return ERROR;   
}
