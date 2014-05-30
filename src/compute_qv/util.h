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
 * 1.20 2003/11/06 18:24:29
 */

#ifndef BTK_UTIL_H_
#define BTK_UTIL_H_

#define ERROR -1
#define SUCCESS 0
#define ENDOFLIST -2
#define BUFLEN 1000

#define MEM_ERROR(a) \
    if ((a) == NULL) { \
	(void)sprintf(message->text, \
			"insufficient memory at file=%s,line=%d\n", \
			__FILE__, __LINE__); \
	goto error; \
    }

#define MIN2(a,b)	((a) < (b) ? (a) : (b))
#define MAX2(a,b)	((a) > (b) ? (a) : (b))

void cg_stat(int *data, int *weight, int n, float *mean, float *std);
void polyfit(float x[], float y[], float yerr[], int n_data,
	     float a[], int na);

void fpoly(float x, float p[], int np);
int qv_round(double x);
int qv_isnan(double x);


/*******************************************************************************
 * Class Stats
 * This struct is intended to be used as an "object"
 * i.e. the data elements are intended to be "private"
 * and access is intended to be only through the routines provided here.
 ******************************************************************************/

#define USE_SIMPLE_STATS_ 0

#if USE_SIMPLE_STATS_
/* Simple (accumulate sum and sum of squares)
 *        (potential numerical problems for variance calculation)
 */
typedef struct
{
    double min_, max_, sum_, sum_sqr_;
    double count_;
} Stats;
#else   /* Robust Stats */
/* Robust (accumulate mean and variance directly)
 *        (This uses a recursion relation which is more numerically
 *         stable than for simple stats)
 */
typedef struct
{
    double min_, max_, mean_, var_;
    double count_;
} Stats;
#endif

/*******************************************************************************
 * Prototypes for "methods"
 ******************************************************************************/

void
statsReset( Stats* s );

void
statsUpdate( Stats* s, double x );

void
statsUpdateNTimes( Stats* s, double x, double n );

void
statsUpdateFromStats( Stats* s, const Stats *ss );

double
statsMin(    const Stats* s );

double
statsMax(    const Stats* s );

double
statsMean(   const Stats* s );

double
statsVar(    const Stats* s );

double
statsStdDev( const Stats* s );

double
statsCount(  const Stats* s );

void
statsPrint(  const Stats* s, FILE* fp );


/*******************************************************************************
 * Class Histo
 * This struct is intended to be used as an "object"
 * i.e. the data elements are intended to be "private"
 * and access is intended to be only through the routines provided here.
 ******************************************************************************/

#define HISTO_SIZE 1000

typedef struct {
    double x_min_, x_max_, bin_to_x_, x_to_bin_;
    double array_[HISTO_SIZE];
    double count_;
} Histo;


/*******************************************************************************
 * Prototypes for "methods"
 ******************************************************************************/

void
histoReset( Histo* h, double xmin, double xmax );

void
histoUpdate( Histo* h, double x );

void
histoUpdateNTimes( Histo* h, double x, double n );

void
histoUpdateFromHisto( Histo* h, const Histo* hh );

double
histoCount( const Histo* h );

double
histoInverseCDF( const Histo* h, double p );

double 
histoMedian( const Histo* h );

double
histoMinusSigma( const Histo* h );

double 
histoPlusSigma( const Histo* h );

double 
histoPercentileStdDev( const Histo* h );

void
histoPrintMedian( const Histo* h, FILE* fp );

void
histoPrint( const Histo* h, FILE* fp );

#endif /* include guard */
