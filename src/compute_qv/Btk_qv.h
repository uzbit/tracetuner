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
 * 1.33 2003/11/06 18:18:41
 */

#define TT_VERSION "TT_4.0" /* Official program version */

#define BTKMESSAGE_LENGTH       (256)
#define	NUM_COLORS		(4)
#define MONITOR			0    /* controls some progress output */
#define USE_OLD_TPARS           0

#define ABI 1
#define SCF 2
#define ZTR 3
#define SFF 4

#define ONE_PLUS (1.0+16.0*DBL_EPSILON) 
#define ONE_MINUS (1.0-16.0*DBL_EPSILON) 

#define F_ONE_PLUS (1.0+16.0*FLT_EPSILON) 
#define F_ONE_MINUS (1.0-16.0*FLT_EPSILON) 
#define QV_PI 3.14159265358979323846
#define TWO_SQRT_LN2 1.6651092223

/* For platform independence, use the following macros instead of casting a 
 * double (or float) to an int
 */
#define INT_DBL(X) ((int)((1.0+16.0*DBL_EPSILON)*(X)))
#define INT_FLT(X) ((int)((1.0+16.0*FLT_EPSILON)*(X)))
#define ROUNDPOS(X)  ( INT_DBL((X)+0.5) )
#define ROUND(X)     ( ((X)>0.0) ? ROUNDPOS(X) : -ROUNDPOS(-(X)) )

#if 1
/* Given:
 *   int i;
 *   double d;
 * For platform independence, use the following macros for comparison 
 * instead of casting the int to a double
 */

/* "fudge" up and down macros */
#define F_UP(INT) ( (INT)>0 ? (+16.0*DBL_EPSILON) : (-16.0*DBL_EPSILON) )
#define F_DN(INT) ( (INT)>0 ? (-16.0*DBL_EPSILON) : (+16.0*DBL_EPSILON) )

#define INT_LT_DBL(I,D)   ((1.0+F_UP(I))*(double)(I)< (D)) /* (double)i< d */
#define INT_LE_DBL(I,D)   ((1.0+F_DN(I))*(double)(I)<=(D)) /* (double)i<=d */
#define INT_GE_DBL(I,D)   ((1.0+F_UP(I))*(double)(I)>=(D)) /* (double)i>=d */
#define INT_GT_DBL(I,D)   ((1.0+F_DN(I))*(double)(I)> (D)) /* (double)i> d */


#define DBL_LT_INT(D,I)   ((D)< (1.0+F_DN(I))*(double)(I)) /* d< (double)i */
#define DBL_LE_INT(D,I)   ((D)<=(1.0+F_UP(I))*(double)(I)) /* d<=(double)i */
#define DBL_GE_INT(D,I)   ((D)>=(1.0+F_DN(I))*(double)(I)) /* d>=(double)i */
#define DBL_GT_INT(D,I)   ((D)> (1.0+F_UP(I))*(double)(I)) /* d> (double)i */

#define INT_EQ_DBL(I,D)   ( fabs((double)(I)-(D)) < 16.0*DBL_EPSILON )
#define DBL_EQ_INT(D,I)   INT_EQ_DBL(I,D)

#define G_UP(DBL) ( (DBL)>0.0 ? (+16.0*DBL_EPSILON) : (-16.0*DBL_EPSILON) )
#define G_DN(DBL) ( (DBL)>0.0 ? (-16.0*DBL_EPSILON) : (+16.0*DBL_EPSILON) )

#define DBL_LT_DBL(A,B)   ((1.0+G_UP(A))*(A) < (B))  
#define DBL_GT_DBL(A,B)   ((1.0+G_DN(A))*(A) > (B))

#else

/* Make things like they used to be: */
#define INT_LT_DBL(I,D)   ((double)(I)< (D))
#define INT_LE_DBL(I,D)   ((double)(I)<=(D))
#define INT_GE_DBL(I,D)   ((double)(I)>=(D))
#define INT_GT_DBL(I,D)   ((double)(I)> (D))

#define DBL_LT_INT(D,I)   ((D)< (double)(I))
#define DBL_LE_INT(D,I)   ((D)<=(double)(I))
#define DBL_GE_INT(D,I)   ((D)>=(double)(I))
#define DBL_GT_INT(D,I)   ((D)> (double)(I))

#define INT_EQ_DBL(I,D)   ((double)(I)==(D))
#define DBL_EQ_INT(D,I)   INT_EQ_DBL(I,D)

#define DBL_LT_DBL(A,B)   ( (A) < (B))  
#define DBL_GT_DBL(A,B)   ( (A) > (B))

#endif

#define DBL_EQ_DBL(A,B)   (fabs((A)-(B))<16.0*DBL_EPSILON)

/* Controls whether or not ContextTable appears in signature of
 * Btk_compute_qv()
 */
#define USE_CONTEXT_TABLE       0

#define CALLOC(type,num) (type *) calloc((unsigned)(num),(sizeof (type)))

#define REALLOC(ptr,type,num) \
        (type *) realloc((type *)(ptr), (unsigned)((num)*sizeof(type)))

#define FREE(x) if((x)!=NULL) {free((char *)(x));(x)=NULL;}

typedef struct {
   int  code;
   char text[BTKMESSAGE_LENGTH];
} BtkMessage;
