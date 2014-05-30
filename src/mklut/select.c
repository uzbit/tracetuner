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

/***************************************************************************
 * 
 * $Id: select.c,v 1.4 2008/11/27 13:09:30 gdenisov Exp $        
 *
 * source: p.187 of the book "Algorithms" by Thomas H. Cormen,
 * Charles E. Leiserson, and Ronald L. Rivest
 *
 ***************************************************************************/

 
#include "select.h"

/*******************************************************************************
 * Function: random_unsigned_long
 * Purpose:  return a random int between a and min(a+32767, b)
 *******************************************************************************
 */
unsigned long
random_unsigned_long(unsigned long a, unsigned long b)
{
    unsigned long r = (unsigned long) rand();
/* returns an int between 0 and 32767 */
    if (a+r > b)
        return a;
    else
        return (a+r);
}

/*******************************************************************************
 * Function: partition
 * Purpose: partition an array A[p..r] so that everything smaller than x is
 * on the left side of the array and everything larger than
 * x is on the right side of the array
 *******************************************************************************
 */
unsigned long
partition(double *A, unsigned long p, unsigned long r)
{
    double x, temp;
    unsigned long i, j;

    x = A[p];
    i = p - 1;
    j = r + 1;
    while (1) {
        /* scan from right until (A[j] <= x) is found */
        do {
            j = j - 1;
        } while (A[j] > x);
        /* scan from left until (A[i] >= x) is found */
        do {
            i = i + 1;
        } while (A[i] < x);
        /* scan and exchange elements until scan indices cross,
           i.e., (i >= j) */
        if (i < j) {
          /* exchange elements of the array */
            temp = A[i];
            A[i] = A[j];
            A[j] = temp;
        }
        else
            return j;
    }
}

/*******************************************************************************
 * Function: randomized_partition
 *******************************************************************************
 */
unsigned long
randomized_partition(double *A, unsigned long p, unsigned long r)
{
    unsigned long i;
    double temp;

    i = random_unsigned_long(p, r);
    temp = A[p];
    A[p] = A[i];
    A[i] = temp;
    return (partition(A,p,r));
}

/*******************************************************************************
 * Function: randomized_select   
 *******************************************************************************
 */
double
randomized_select(double *A, unsigned long p, unsigned long r,
    unsigned long i)
{
    unsigned long q, k;

    if (p == r)
        return (A[p]);
    q = randomized_partition(A, p, r);
    k = q - p + 1;
    if (i <= k)
        return (randomized_select(A, p, q, i));
    else
        return (randomized_select(A, q+1, r, i-k));
}

/*******************************************************************************
 * Function: quicksort
 *******************************************************************************
 */
void
quicksort(double *A, unsigned long p, unsigned long r)
{
    unsigned long q;

    if (p < r) {
        q = randomized_partition(A, p, r);
        quicksort(A, p, q);
        quicksort(A, q+1, r);
    }
}

/*******************************************************************************
 * Function: partition2
 * Purpose: partition an array A[p..r] so that everything smaller than x is
 *          on the left side of the array and everything larger than
 *          x is on the right side of the array
 *******************************************************************************
 */
unsigned long
partition2(double *A, int *w, unsigned long p, unsigned long r)
{
    double x, temp, itemp;
    unsigned long i, j;

    x = A[p];
    i = p - 1;
    j = r + 1;
    while (1) {
        /* scan from right until (A[j] <= x) is found */
        do {
            j = j - 1;
        } while (A[j] > x);
        /* scan from left until (A[i] >= x) is found */
        do {
            i = i + 1;
        } while (A[i] < x);
        /* scan and exchange elements until scan indices cross,
           i.e., (i >= j) */
        if (i < j) {
          /* exchange elements of the array */
            temp = A[i];
            A[i] = A[j];
            A[j] = temp;
            itemp = w[i];
            w[i]  = w[j];
            w[j]  = (int)itemp;
        }
        else
            return j;
    }
}

/*******************************************************************************
 * Function: randomized_partition2
 *******************************************************************************
 */ 
unsigned long 
randomized_partition2(double *A, int *w, unsigned long p, unsigned long r)
{
    unsigned long i;
    double temp;
    int    itemp;
 
    i = random_unsigned_long(p, r);
    temp = A[p];
    A[p] = A[i];
    A[i] = temp;
    itemp = w[p];
    w[p]  = w[i];
    w[i]  = itemp;
    return (partition2(A, w, p, r));
} 
 
/*******************************************************************************
 * Function: quicksort2
 *******************************************************************************
 */
void
quicksort2(double *A, int *w, unsigned long p, unsigned long r)
{
    unsigned long q;

    if (p < r) {
	q = randomized_partition2(A, w, p, r);
	quicksort2(A, w, p, q);
	quicksort2(A, w, q+1, r);
    }
}
