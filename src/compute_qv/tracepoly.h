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
 * 2.8 2003/11/06 18:18:50
 */

#ifndef _TRACEPOLY_H_
#define _TRACEPOLY_H_

/* The degree of the poly = NumCoef-1
 * I.e., NumCoef=3  ===> quadratic (x^2)
 * 
 * Don't set this to anything but 1 or 2 unless you know what you are doing */

#define NumCoef 2       /* for even poly ===> quadratic (x^2) */
#define MaxNumBases 4000
#define TELL_ME_IF_I_EXCEED_THE_LIMITS 1 /* debug */

/*******************************************************************************
 *   typedefs for the prototypes below 
 ******************************************************************************/

typedef struct {
    int num_coef_;
    double coef_[NumCoef];
    double x_min_, x_max_, y_max_, x_cutoff_, y_cutoff_;
} TracePoly;

typedef TracePoly ReadPolynomials[4];

typedef struct {
    ReadPolynomials read_polys;
    double x_min_, x_max_;
} ReadInfo;

/*******************************************************************************
 *   prototypes
 ******************************************************************************/

/* Need the following for mklut directory (context.c) */
double traceEvalFunc(TracePoly *, double);
void   readPrint(ReadInfo *, FILE *);
double readEvalAveFunc(ReadInfo *, int);

int readGetCoefficients(ReadInfo *, double [][MaxNumBases],
                        double [][MaxNumBases], int [], int);
double readNormFactor(ReadInfo *, int, int);

#endif /* _TRACEPOLY_H_ */




