/**************************************************************************
 * This file is part of TraceTuner, the DNA sequencing quality value,
 * base calling and trace processing software.
 *
 * Copyright (c) 2003-2008 Gennady Denisov.  All rights reserved.
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

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "util.h"

/* The following allows me to use this file stand-alone with its test drivers
 * which is useful for development, debugging, testing, etc.
 */
#if 0
#define DRIVER_STATS   /* or -DDRIVER in Makefile */
#endif
#if 0
#define DRIVER_HISTO   /* or -DDRIVER in Makefile */
#endif

#if !( defined(DRIVER_STATS) || defined(DRIVER_HISTO) ) 

#include "Btk_qv.h"
#include "Btk_qv_data.h"
#include "Btk_qv_funs.h"
#include "nr.h"
#include "train.h"
#include "tracepoly.h"
#include "context_table.h"
#include "Btk_call_bases.h"
#include "Btk_process_raw_data.h"

#define MAX_PEAK_HEIGHT 1600

/* This is also defined in call_bases.c
 * It is not an "extern" because C and C++ have slightly different
 * ways of treating "extern", and I want to compile with either compiler
 */
static const float	gSpacEst0 = 12.0;


/*******************************************************************************
 * Function: cg_stat
 * Purpose: Computed robust or weighted average and standard deviation of
 *          integer data
 *******************************************************************************
 */
void cg_stat(int *data, int *weight, int n, float *mean, float *std)
	/* compute weighted mean and standard deviation */
        /* if weight == NULL, compute robust mean and standard deviation */
	/* Curtis Gehman, 2000.11.13 */
{
  const double eps = 0.10;
	    /* probability that a legitimate datum will be discarded */
	    /* as an outlier is < eps, assuming normal distribution */

  int i, j, m;
  float sum = 0, sumsqr = 0;
  double x, p;
  float	tot_wgt = 0;

	    /* preliminary estimates */
  if ( weight == NULL ) {
    for ( i=0; i < n; ++i ) {
      sum += data[i];
      sumsqr += data[i] * data[i];
    }
    *mean = sum / n;
    *std = sumsqr / n - *mean * *mean;
    *std = sqrt(*std);
  } else {
    for ( i=0; i < n; ++i ) {
      tot_wgt += weight[i];
      sum += weight[i] * data[i];
      sumsqr += weight[i] * data[i] * data[i];
    }
    *mean = sum / tot_wgt;
    *std = sumsqr / tot_wgt - *mean * *mean;
    *std = sqrt(*std);
  }

	    /* check for outliers */
  if ( weight == NULL ) {
    x = 1.41421; /* square root of 2 */
    x *= *std;
    p = 1 - eps;
    p = 1 - pow(p, 1.0/n);  /* probability needed for individual instance */
    for ( j = INT_FLT(*std); (1.0 - Erf(j/x)) > p; ++j ) ;
    for ( i=m=0; i < n; ++i )
      if ( fabs(data[i] - *mean) > j ) {
	m++;  /* keep track of number of data discarded */
	sum -= data[i];
	sumsqr -= data[i]*data[i];
      }
    if ( m > n/5 )
      fprintf(stderr,
	      "cg_stat: Warning - discarding %3d of %4d data.\n", m, n);
    m = n - m;
	/* refine statistics */
    *mean = (float) sum / m;
    *std = (float) sumsqr / m - *mean * *mean;
    *std = sqrt(*std);
  }

  return;
}


void polyfit(float x[], float y[], float sig[], int n_data, float a[], int na)

   /************************************************************************
    * Description: Find a polynomial fit to data in x[] and y[].
    *
    * Inputs:   x, y    arrays of data
    *		sig	array of uncertainties of y's
    *           n_data  length of both data arrays
    *           na      order of polynomial + 1 (e.g., 4 for a cubic fit)
    * Outputs:  a       polynomial coefficients; a[n] is coefficient of x^n
    * Return:   void
    * Comments: If no uncertainties are available pass a NULL pointer for sig.
    *		Uses NR singular value decomposition fitting routine.
    */

{
    int     i;
    // float   **u, *w;
    float   **v;
    int     *ia;
    float   chisq;
    int		free_sig = 0;
    // float   **matrix(long nrl, long nrh, long ncl, long nch);  // NR
    // void    free_matrix(float **m, long nrl, long nrh, long ncl, long nch);

    if ( MONITOR > 2 ) fprintf(stderr, "  polyfit() ...\n");

    // Allocate memory required by NR (NR's unusual index convention)
    // u = matrix(1, n_data, 1, na);
    v = matrix(1, na, 1, na);
    // w = (float *) malloc(na*sizeof(float));

    if ( sig == NULL ) {	// Create dummy uncertainties
	sig = (float *) malloc(n_data*sizeof(float));
	free_sig = 1;
	for ( i=0; i < n_data; ++i ) sig[i] = 1.0;
    }

    // svdfit(x-1, y-1, sig-1, n_data, a-1, na, u, v, w-1, &chisq, fpoly);
    if ( MONITOR > 1 ) {
    // fprintf(stderr, "    singular values = (");
    // for ( i=0; i < na; ++i )
        // fprintf(stderr, "%12.3e%s", w[i], i==na-1 ? ")\n" : ",");
    }

    ia = (int *) malloc(na*sizeof(int));
    for ( i=0; i < na; ++i ) ia[i] = 1;
    lfit(x-1, y-1, sig-1, n_data, a-1, ia-1, na, v, &chisq, fpoly);

    free(ia);
    // free_matrix(u, 1, n_data, 1, na);
    free_matrix(v, 1, na, 1, na);
    // free(w);
    if ( free_sig ) free(sig);
}

#endif

void fpoly(float x, float powers[], int np)

   /************************************************************************
    * Description: Compute powers of x up to np.
    *
    * Inputs:	x
    *		np	highest power of x to compute
    * Outputs:	p	NR array of powers of x (base offset of 1), length np+1
    * Return:
    * Comments:
    */

{
    int j;

    powers[1] = 1.0;
    for ( j=2; j <= np; ++j ) 
        powers[j] = powers[j-1]*x;
}


/**********************************************************************
 * Round off to the nearest integer.
 *
 * Inputs:	x	Floating point value to qv_round off.
 *
 * Outputs:	
 * Return:	Integer nearest x.
 * Comments:	
 */
int qv_round(double x)
{
  if ( x > 0 ) return INT_DBL(x + 0.5);
  else if ( x < 0 ) return INT_DBL(x - 0.5);
  else return INT_DBL(x);
}


/*******************************************************************************
 * Implementation of class Stats (as close to a class as I can get in "C")
 *
 * This code implements a class for 
 * accumulating statistics with one of two possible implementations:
 * The implementations differ in the way in which they trade off numerical 
 * accuracy vs. speed.
 * simple implementation: (accumulate sum and sum of squares)
 *           (potential numerical problems for variance calculation)
 * robust implementation: (accumulate mean and variance directly)
 *           (This uses a recursion relation which is more numerically
 *           stable than above, but it is somewhat more expensive)
 *
 * Usage:
 * First declare an object of type Stats:
 *     Stats stats;
 * Make sure it is reset before using it:
 *     statsReset(&stats);
 * Inside a loop in which "rv" represents a random variable for
 * which we wish to collect statistics:
 *     statsUpate(&stats,rv);
 * After collecting statistics, we may use one of several "methods"
 * to get the (sample) min, max, mean, standard deviation, count,
 * or to print out the various statistics.
 *
 * The interface for Stats is as follows 
 * (the 1st argument of type Stats* is suppressed):
 *    statsReset()  --> sets all the statistical "accumulators" to zero
 *    statsUpdate(rv) --> increments the various "accumulators"
 *    statsUpdateNTimes(rv,n) --> equivalent to calling statsUpdate(rv) n times
 *                                (note that n need not be an integer)
 *    statsUpdateFromStats(Stats) --> increments accumulators using Stats
 *    statsMin()    --> returns the min value of the r.v.
 *    statsMax()    --> returns the max value of the r.v.
 *    statsMean()   --> returns the mean value of the r.v.
 *    statsVar()    --> returns the variance of the r.v.
 *    statsStdDev() --> returns the standard deviation of the r.v.
 *    statsCount(fp)--> returns the "count" to "fp"
 *                      (note that "count" need not be an integer
 *                       see statsUpdateNTimes() above)
 *    statsPrint(fp)--> prints min, max, mean, std. dev., count
 *
 ******************************************************************************/


#define MYMIN(A,B)      ( (A)<(B) ? (A) : (B) )
#define MYMAX(A,B)      ( (A)>(B) ? (A) : (B) )
#define MYSQR(X)        ( (X)*(X) )


/*******************************************************************************
 * Implementation of Class Stats
 ******************************************************************************/

#if USE_SIMPLE_STATS_ /********************************************************/

/* Faster but less robust implementation */

void
statsReset( Stats* s )
{
    s->min_     =  DBL_MAX;
    s->max_     = -DBL_MAX;
    s->sum_     = 0.0;
    s->sum_sqr_ = 0.0;
    s->count_   = 0;
}

void
statsUpdate( Stats* s, double x )
{
    s->min_      = MYMIN( s->min_, x );
    s->max_      = MYMAX( s->max_, x );
    s->sum_     += x;
    s->sum_sqr_ += x*x;
    s->count_   += 1;
}

void
statsUpdateNTimes( Stats* s, double x, double n )
{
    s->min_      = MYMIN( s->min_, x );
    s->max_      = MYMAX( s->max_, x );
    s->sum_     += n*x;
    s->sum_sqr_ += n*x*x;
    s->count_   += n;
}

void
statsUpdateFromStats( Stats* s, const Stats *ss )
{ 
    s->min_      = MYMIN( s->min_, ss->min_ );
    s->max_      = MYMAX( s->max_, ss->max_ );
    s->sum_     += ss->sum_;
    s->sum_sqr_ += ss->sum_sqr_;
    s->count_   += ss->count_;
}

double
statsMean( const Stats* s )
{
    if( s->count_ > 0 ) {  return (s->sum_)/(double)(s->count_);  } 
    else                {  return 0.0;  }
}

double
statsVar( const Stats* s )
{ 
    if( s->count_ > 1 ) {
        double mean = statsMean(s), dcount=s->count_;
        return ( (s->sum_sqr_ - dcount*MYSQR(mean)) / (dcount-1.0) );
    } else {
        return 0.0;
    }
}


#else  /* USE_SIMPLE_STATS_ ***************************************************/

/* More robust but slower implementation */

void
statsReset( Stats* s )
{
    s->min_   =  DBL_MAX;
    s->max_   = -DBL_MAX;
    s->mean_  = 0.0;
    s->var_   = 0.0;
    s->count_ = 0;
}

void
statsUpdate( Stats* s, double x )
{
    double old_mean = s->mean_, n;
    s->min_    = MYMIN( s->min_, x );
    s->max_    = MYMAX( s->max_, x );
    s->count_ += 1;

    n = s->count_;
    s->mean_ = s->mean_*((n-1.0)/n) + x/n;	/* recursion for mean */
    {   /* recursion for variance */
        double temp1 = (x-old_mean)/n, temp2 = x-s->mean_;        
        s->var_ = (n-1.0)*s->var_ + (n-1.0)*MYSQR(temp1) + MYSQR(temp2);
        s->var_ /= n;
    }
}
 
void
statsUpdateNTimes( Stats* sn, double m_x, double m_count )
{
    double n, m, un, um, up, vn, vm, vp, tn, tm;
    /* Mnemonics of variable names: 
     * n = count for sn
     * m = count for sm
     * p = n+m
     * u = Greek letter mu (for mean)
     * v = variance
     * t = temporary
     */

    /* Take care of degenerate cases */
    if( m_count    == 0.0 ) { return; }
    if( sn->count_ == 0.0 ) { sn->mean_=m_x; sn->count_=m_count; return; }

    sn->min_ = MYMIN( sn->min_, m_x );
    sn->max_ = MYMAX( sn->max_, m_x );

    n  = sn->count_;
    un = sn->mean_;
    vn = sn->var_;

    m  = m_count;
    um = m_x;
    vm = 0.0;

    up = ( n*un + m*um ) / (n+m);	/* recursion for mean */
    tn = un - up;
    tm = um - up;

    /* recursion for variance */
    vp = ( n*vn + n*MYSQR(tn) + m*vm + m*MYSQR(tm) ) / (n+m);    

    sn->count_ += m_count;
    sn->mean_   = up;
    sn->var_    = vp;
}

void
statsUpdateFromStats( Stats* sn, const Stats* sm )
{
    double n, m, un, um, up, vn, vm, vp, tn, tm;
    /* Mnemonics of variable names: 
     * n = count for sn
     * m = count for sm
     * p = n+m
     * u = Greek letter mu (for mean)
     * v = variance
     * t = temporary
     */

    /* Take care of degenerate cases */
    if( sm->count_ == 0.0 ) { return; }
    if( sn->count_ == 0.0 ) { *sn = *sm; return; }

    sn->min_ = MYMIN( sn->min_, sm->min_ );
    sn->max_ = MYMAX( sn->max_, sm->max_ );
    
    n  = sn->count_;
    un = sn->mean_;
    vn = sn->var_;

    m  = sm->count_;
    um = sm->mean_;
    vm = sm->var_;

    up = ( n*un + m*um ) / (n+m);	/* recursion for mean */
    tn = un - up;
    tm = um - up;

    /* recursion for variance */
    vp = ( n*vn + n*MYSQR(tn) + m*vm + m*MYSQR(tm) ) / (n+m);
    
    sn->count_ += sm->count_;
    sn->mean_   = up;
    sn->var_    = vp;
}

double
statsMean(   const Stats* s )  {  return (s->mean_);  }

double
statsVar(    const Stats* s )  
{
    if( s->count_ > 1.0 ) {
        return s->var_*s->count_/(s->count_-1.0); 
    } else {
        return 0.0;
    }
}

#endif /* USE_SIMPLE_STATS_ ***************************************************/

double
statsMin(    const Stats* s )  {  return s->min_;  }

double
statsMax(    const Stats* s )  {  return s->max_;  }

double
statsStdDev( const Stats* s )  {  return sqrt(statsVar(s));  }

double
statsCount(  const Stats* s )  {  return s->count_;  }

void
statsPrint(  const Stats* s, FILE* fp )
{
    fprintf( fp,
      /*"(min,max)=(%25.16e,%25.16e) (mean,std)=(%25.16e,%25.16e) ctr=%d\n",*/
             "(min,max)=(%11.4g,%11.4g) (mean,std)=(%11.4g,%11.4g) ctr=%f\n",
             statsMin(s), statsMax(s),  
             statsMean(s), statsStdDev(s), statsCount(s) );
}

/*******************************************************************************
 * The following code serves as both:
 * (1) An example of how to use the class Stats, and
 * (2) A debugging test to see if it is implemented correctly
 * This is a stand alone driver which will compile and run as-is.
 * (No other .c files need be linked)
 ******************************************************************************/

#ifdef DRIVER_STATS

//#include "Stats.h"

/* From Numerical Recipes (some code deleted to avoid static variables) */
double gasdev( unsigned short seed16v[3] )
{
    double fac,rsq,v1,v2;
    
    do {
        v1=2.0*erand48(seed16v)-1.0;
        v2=2.0*erand48(seed16v)-1.0;
        rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    return v2*fac;
}

void 
test_stats1(void)
{
    Stats         x_stats;
    int exponent;
    
    for( exponent=0; exponent<1; exponent++ )
    {
        double mean = pow(10.0,exponent);
        int trials;
        fprintf( stdout, "exponent = %d\n", exponent );

        for( trials=100; trials<=100; trials++ )
        {
            unsigned short seed16v[3];
            int t;
            
            fprintf( stdout, "trials = %d\n", trials );
            /* state for erand48() random number generator */
            seed16v[0]=0x0123; seed16v[1]=0x4567; seed16v[2]=0x89ab; 
            statsReset( &x_stats );
            
            for( t=0; t<trials; t++ )
            {   // NOTICE
                double x = mean + gasdev(seed16v);                
                statsUpdate( &x_stats, x );
            }
            statsPrint( &x_stats, stdout );
        }
    }
}

void
test_stats2(void)
{
    Stats stats0,  stats1,  stats2,  stats3, 
          stats0a, stats1a, stats2a, stats3a;
    int trials = 3;
    unsigned short seed16v[3];
    int t;
    
    //for( trials = 0; trials<=6; trials++ )
    trials = 4000;
    {
        seed16v[0]=0x0123; seed16v[1]=0x4567; seed16v[2]=0x89ab; 
        
        statsReset( &stats0 );
        statsReset( &stats0a );
        statsReset( &stats1 );
        statsReset( &stats1a );
        statsReset( &stats2 );
        statsReset( &stats2a );
        statsReset( &stats3 );
        statsReset( &stats3a );
        
        for( t=0; t<trials; t++ )
        {
            double x = gasdev(seed16v);
            //if ( t%2 == 0 ) { x=-1.0; }
            //else            { x=+1.0; }
            //printf( " %7.4f", x );
#if 0
            statsUpdate( &stats3a, x );
            if ( t%2 == 0 ) {  statsUpdate( &stats1a, x );  } 
            else            {  statsUpdate( &stats2a, x );  }
#else
            double n = 0.3; //0.3;
            statsUpdateNTimes( &stats3a, x, n );
            if ( t%2 == 0 ) {  statsUpdateNTimes( &stats1a, x, n );  } 
            else            {  statsUpdateNTimes( &stats2a, x, n );  }
#endif
        }
        //printf( "\n" );

        statsUpdateFromStats( &stats0, &stats3a );
#if 1
        statsUpdateFromStats( &stats1, &stats0a );
        statsUpdateFromStats( &stats1, &stats1a );
        //stats1 = stats1a;
        statsUpdateFromStats( &stats1, &stats0a );
        statsUpdateFromStats( &stats1, &stats2a );
        statsUpdateFromStats( &stats1, &stats0a );
        
        statsUpdateFromStats( &stats2, &stats0a );
        statsUpdateFromStats( &stats2, &stats2a );
        //stats2 = stats2a;
        statsUpdateFromStats( &stats2, &stats0a );
        statsUpdateFromStats( &stats2, &stats1a );
        statsUpdateFromStats( &stats2, &stats0a );
        
        statsUpdateFromStats( &stats3, &stats0a );
        statsUpdateFromStats( &stats3, &stats3a );
        statsUpdateFromStats( &stats3, &stats0a );
#endif
        //statsPrint( &stats1a, stdout );
        //statsPrint( &stats2a, stdout );
        statsPrint( &stats3a, stdout );
        statsPrint( &stats0, stdout );
#if 1
        statsPrint( &stats1, stdout );
        statsPrint( &stats2, stdout );
        statsPrint( &stats3, stdout );
#endif
    }
}

int
main( int argc, char* argv[] )
{
    /* printf( "(min,max)= %24.12e %24.12e\n", -DBL_MAX, DBL_MAX ); */
    test_stats2();
    return 0;
}

#endif          /* DRIVER_STATS */


/*******************************************************************************
 * Implementation of class Histo (as close to a class as I can get in "C")
 *
 * Usage:
 * First declare an object of type Histo:
 *     Histo histo;
 * Make sure it is reset before using it:
 *     histoReset(&histo);
 * Inside a loop in which "rv" represents a random variable for
 * which we wish to collect statistics:
 *     histoUpate(&histo,rv);
 * After collecting statistics, we may use one of several "methods"
 * to get the median, x-th percentile, count, etc.,
 * or to print out the various statistics.
 *
 * The interface for Histo is as follows 
 * (the 1st argument of type Histo* is suppressed):
 *    histoReset(xmin,xmax) --> reset histogram bins and counter to zero.
 *                              set min and max of the bins.
 *    histoUpdate(rv)   --> update internal state with new random variable
 *                           (i.e., increment appropriate bin)
 *    histoUpdateNTimes(rv,n) --> equivalent to calling histoUpdate(rv) n times
 *                                (note that n need not be an integer)
 *    histoUpdateFromHisto(Histo) --> update internal state from another Histo
 *    histoCount()  --> returns the "count" 
 *                      (note that "count" need not be an integer
 *                       see histoUpdateNTimes() above)
 *    histoInverseCDF(x)--> histoInverseCDF(0.5) returns median value
 *                      --> histoInverseCDF(0.1587) returns -sigma (for normal)
 *                      --> histoInverseCDF(0.8413) returns +sigma (for normal)
 *    histoMedian()     --> returns histoInverseCDF(0.5)
 *    histoMinusSigma() --> returns histoInverseCDF(0.1587)
 *    histoPlusSigma()  --> returns histoInverseCDF(0.8413)
 *    histoPercentileStdDev()->returns 0.5*(histoPlusSigma()-histoMinusSigma())
 *    histoPrintDbg(fp) --> prints out internal state (mainly for debugging)
 *    histoPrint(fp)    --> prints out some of the above statsistics to "fp"
 *
 * Assumptions:
 *   rv = random variable
 *   xmin <= rv <= xmax
 * If rv is outside this range it is "clamped" to put in the histogram.
 * I.e., 
 *       if(rv < xmin) { rv=xmin; }
 *       if(rv > xmax) { rv=xmax; }
 * The resolution (size of the bins) is given by (xmax-xmin)/HISTO_SIZE
 *
 ******************************************************************************/

#define ROUNDPOS(X)     ( INT_DBL((X)+0.5) )
#define ROUND(X)        ( ((X)>0.0) ? ROUNDPOS(X) : -ROUNDPOS(-(X)) )
#define MIN(A,B)        ( ((A)<(B)) ? (A) : (B) )
#define MAX(A,B)        ( ((A)>(B)) ? (A) : (B) )


void
histoReset( Histo* h, double xmin, double xmax )
{
    int i;
    if (xmin>=xmax) {
        fprintf( stderr, "screwy histogram: min=%g max=%g\n", xmin, xmax );
        exit(-1);
    }
    h->x_min_ = xmin;
    h->x_max_ = xmax;
    h->bin_to_x_   = (xmax-xmin)/HISTO_SIZE;
    h->x_to_bin_   = HISTO_SIZE/(xmax-xmin);
    //printf( "bin_to_x=%f  x_to_bin=%f\n",  h->bin_to_x_, h->x_to_bin_  );
    for( i=0; i<HISTO_SIZE; i++ ) { h->array_[i] = 0; }
    h->count_ = 0;
}


void
histoUpdateNTimes( Histo* h, double xx, double n )
{
    double x = MAX( h->x_min_, MIN( xx, h->x_max_ ) );
    double temp = (x-h->x_min_)*h->x_to_bin_ - 0.5;
    int bin = ROUND(temp);
    bin = MAX( 0, bin );
    bin = MIN( HISTO_SIZE-1, bin );
    h->array_[bin] += n;
    h->count_ += n;
}

void
histoUpdate( Histo* h, double xx )
{
    histoUpdateNTimes( h, xx, 1.0 );
}


void
histoUpdateFromHisto( Histo* h, const Histo *hh )
{
    int i;
    for ( i=0; i<HISTO_SIZE; i++ ) {
        h->array_[i] += hh->array_[i];
    }
    h->count_ += hh->count_;
}


void
histoPrintDbg( const Histo* h, FILE* fp )
{
    int i;
    for( i=0; i<HISTO_SIZE; i++ ) {
        fprintf( fp, "  %6.2f %2f", 
                 (double)i*h->bin_to_x_+h->x_min_, h->array_[i] );
    }    
    fprintf( fp, "\n" );
}


double
histoInverseCDF( const Histo* h, double p )
{
    double v0, v1, vtest, d0, d1, retval;
    int i0, i1;

    assert( p>=0 );
    assert( p<=1.0 );
    vtest = p*h->count_;

    if( p<0.5 ) {
        /*if( 0 ) { */
        v0 = v1 = 0;
        i1 = -1;
        do {
            i1++;
            v0 = v1;
            v1 += h->array_[i1];
        } while( v1 < vtest );
        i0 = i1 - 1;
        d0 = i0+1.0;
        d1 = i1+1.0;
    } else {
        v0 = v1 = h->count_;
        i0 = HISTO_SIZE;
        /*fprintf( stderr, "v0=%f i0=%d vtest=%f\n", v0, i0, vtest ); */
        do {
            i0--;
            v1 = v0;
            v0 -= h->array_[i0];
            /*fprintf( stderr, 
              "v1=%f v0=%f i0=%d  array=%d\n", v1, v0, i0, h->array_[i0] ); */
        } while( v0 >= vtest );
        i1 = i0 + 1;
        d0 = i0;
        d1 = i1;
    }

    retval = ( d0*(v1-vtest) + d1*(vtest-v0) ) / ((v1-v0)*h->x_to_bin_);    
    /*fprintf( stderr, 
    "p = %f d = %f %f  v = %f %f %f  r= %f\n", 
    p, d0, d1, v0, vtest, v1, retval ); */
    return retval + h->x_min_;
}


double histoMedian( const Histo* h )
{
    return histoInverseCDF( h, 0.5 );
}


double histoMinusSigma( const Histo* h )
{
    return histoInverseCDF( h, 0.1587 );
}


double histoPlusSigma( const Histo* h )
{
    return histoInverseCDF( h, 0.8413 );
}

double histoPercentileStdDev( const Histo* h )
{
    return 0.5 * ( histoPlusSigma(h) - histoMinusSigma(h) );
}


double
histoCount( const Histo* h )
{
    return h->count_;
}

void
histoPrint( const Histo* h, FILE* fp )
{
    double med, msig, psig, sig;
    med  = histoMedian(h);
    msig = histoMinusSigma(h);
    psig = histoPlusSigma(h);
    sig  = histoPercentileStdDev(h);
    fprintf( fp, "(15.87%%,50%%,84.13%%)=(%4.2f,%4.2f,%4.2f) 'sig'=%4.2f\n",
             msig, med, psig, sig );
}

/*******************************************************************************
 * The following code serves as both:
 * (1) An example of how to use the class Histo, and
 * (2) A debugging test to see if it is implemented correctly
 * This is a stand alone driver which will compile and run as-is.
 * (No other .c files need be linked)
 ******************************************************************************/


#ifdef DRIVER_HISTO

void
test1(void)
{
    Histo h;
    const int size = 5;
    int i;
    double x;

    histoReset( &h, 0.0, 4.0 );
    for( x=0.0; x<1.01; x+=0.25 ) {
        histoUpdate( &h, MIN(0.99999,x) );
    }
    histoUpdate( &h, 0 );
    /* histoPrintDbg( &h, stderr ); */

    fprintf( stderr, "count= %f\n", histoCount(&h) );

    for( i=0; i<size; i++ ) 
    {
        double p =  (double)i/(double)(size-1);
        /* double p = 1.0; */
        double x = histoInverseCDF( &h, p );
        /* fprintf( stderr, "p=%g  p-1=%g\n", p, p-1.0 ); */
        fprintf( stderr, "   %4.2f %6.3f\n", p, x );
    }
}

/* From Numerical Recipes (some code deleted to avoid static variables) */
double gasdev( unsigned short seed16v[3] )
{
    double fac,rsq,v1,v2;
    
    do {
        v1=2.0*erand48(seed16v)-1.0;
        v2=2.0*erand48(seed16v)-1.0;
        rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    return v2*fac;
}


void
test2(void)
{
    Histo h, h1, h2, h3;
    unsigned short seed16v[3];
    int trials=1000, t;
    double mean = 0.0;

    seed16v[0]=0x0123; seed16v[1]=0x4567; seed16v[2]=0x89ab; 
    histoReset(&h,  -2.0, 2.0 );
    histoReset(&h1, -2.0, 2.0);
    histoReset(&h2, -2.0, 2.0);
    histoReset(&h3, -2.0, 2.0);

    for( t=0; t<trials; t++ )
    {
        double x = mean + gasdev(seed16v);
#if 0
        histoUpdate( &h, x );
        if ( t%2 == 0 ) { histoUpdate( &h1, x ); }
        else            { histoUpdate( &h2, x ); }
#else
        double n = 3.14159;
        histoUpdateNTimes( &h, x, n );
        if ( t%2 == 0 ) { histoUpdateNTimes( &h1, x, n ); }
        else            { histoUpdateNTimes( &h2, x, n ); }
#endif
    }
    
    histoUpdateFromHisto( &h3, &h1 );
    histoUpdateFromHisto( &h3, &h2 );

    histoPrint( &h, stdout );
    histoPrint( &h1, stdout );
    histoPrint( &h2, stdout );
    histoPrint( &h3, stdout );
    //histoPrintDbg( &h, stdout );
}

int
main(void)
{
    fprintf( stderr, "version 1\n" );
    test2();
    return 0;
}

#endif          /* DRIVER */


#if !( defined(DRIVER_STATS) || defined(DRIVER_HISTO) )

/*****************************************************************************
 * Function: mobShiftQuality(Peak peak[], int n_peaks)
 * Purpose:  Return a statistical measure of the quality of corrected ipos
 *           Measure the quality of mobility shift corrections by examining
 *           statistics of the spacings of the corrected and original positions
 *           of intrinsic peaks.
 *
 * Inputs:
 *	peak		array of pointers to peaks to consider
 *	n_peaks		length of array of peaks
 *
 * Outputs:		none
 * Return:		measure of quality of shifted peak positions
 * Comments:	        rv means "random variable" in comments below
 */
double
mobShiftQuality(Peak * const peak[], int n_peaks)
{
    int         i, chunk_num, ctr=0;
    const Peak *pk0, *pk1;
    Histo       spac_histo_old, spac_histo_new;
    double      sum_hite[MAX_CHUNKS],
                spacing_old[  MAX_CHUNKS],   spacing_new[  MAX_CHUNKS],
                norm_spac_old[MAX_CHUNKS],   norm_spac_new[MAX_CHUNKS];
    int         count[MAX_CHUNKS];
    double      sum_h=0., sum_height_tot=0;

    /* Initialization */
    spac_histo_new.count_ = 0;
    spac_histo_new.x_to_bin_ = 0;
    spac_histo_new.x_min_ = INF;
    spac_histo_new.x_max_ = -1;
    spac_histo_old.count_ = 0;
    spac_histo_old.x_to_bin_ = 0;
    spac_histo_old.x_min_ = INF;
    spac_histo_old.x_max_ = -1;

    /* return 0; */ /* kludgy way of turning this routine off  */

    /* 1st pass over data to get stats on peak spacing (old & new) and height 
     */
    for ( i=1, pk0=peak[0], chunk_num=-1; i<n_peaks; i++, pk0=pk1 ) {
        pk1 = peak[i];
        if ( pk1->ipos / SCANS_PER_CHUNK > chunk_num ) {
            if ( chunk_num>=MAX_CHUNKS-1 ) { break; }
            chunk_num++;
            if( chunk_num > 0 ) {
                /* store some stats from last chunk */
                count[   chunk_num-1] = ctr;
                sum_hite[chunk_num-1] = sum_h;
                spacing_old[chunk_num-1] = histoMedian( &spac_histo_old );
                spacing_new[chunk_num-1] = histoMedian( &spac_histo_new );
            }
            /* Reset all statisics accumulators */ 
            ctr = 0;
            sum_h = 0.0;
            histoReset( &spac_histo_old, 0.0, 5.0*gSpacEst0 ); /* min --> max */
            histoReset( &spac_histo_new, 0.0, 5.0*gSpacEst0 ); /* min --> max */
        }
        /* Accumulate statistics for this chunk */
        ctr++;
        sum_h += pk1->iheight;
        sum_height_tot += pk1->iheight;
        /* (1) rv = peak spacing (not normalized)
         * (2) stats are accumulated according to "weight" given by iheight
         */
        histoUpdateNTimes( &spac_histo_old, pk1->ipos-pk0->ipos,
                           pk1->iheight );
        histoUpdateNTimes( &spac_histo_new, pk1->ipos     -pk0->ipos,
                           pk1->iheight );
    }
    /* store some stats from last chunk */
    count[   chunk_num] = ctr;
    sum_hite[chunk_num] = sum_h;
    spacing_old[chunk_num] = histoMedian( &spac_histo_old );
    spacing_new[chunk_num] = histoMedian( &spac_histo_new );


    /* 2nd pass over data to get stats on normalized peak spacing (old & new)
     */
    for ( i=1, pk0=peak[0], chunk_num=-1; i<n_peaks; i++, pk0=pk1 ) {
        pk1 = peak[i];
        if ( pk1->ipos / SCANS_PER_CHUNK > chunk_num ) {
            if ( chunk_num+1>=MAX_CHUNKS ) { break; }
            chunk_num++;
            if( chunk_num > 0 ) {
                /* store peak spacing stats from last chunk */
                norm_spac_old[chunk_num-1] = 
                    histoPercentileStdDev( &spac_histo_old );
                norm_spac_new[chunk_num-1] = 
                    histoPercentileStdDev( &spac_histo_new );
            }
            /* Reset all statisics accumulators */ 
            histoReset( &spac_histo_old, 0.0, 5.0 ); /* min --> max */
            histoReset( &spac_histo_new, 0.0, 5.0 ); /* min --> max */
        }
        /* Accumulate peak spacing statistics for this chunk:
         * (1) rv = peak spacing normalized by representative spacing
         * (2) stats for rv are accumulated by "weight" given by iheight
         */
        histoUpdateNTimes( &spac_histo_old, 
                     (pk1->ipos - pk0->ipos)/spacing_old[chunk_num],
                           pk1->iheight );
        histoUpdateNTimes( &spac_histo_new, 
                     (pk1->ipos      - pk0->ipos     )/spacing_new[chunk_num],
                           pk1->iheight );
    }
    /* store peak spacing stats from last chunk */
    norm_spac_old[chunk_num] = histoPercentileStdDev( &spac_histo_old );
    norm_spac_new[chunk_num] = histoPercentileStdDev( &spac_histo_new );

    /* "Integrate" over read to get a single statistic comparing 
     *  old vs. new peak spacing
     */
    {
        double sum_tot=0.0;

        for( i=0; i<=chunk_num; i++ ) {
            double temp = (norm_spac_old[i]-norm_spac_new[i]) / 
                (norm_spac_old[i]+norm_spac_new[i]),
                weight = sum_hite[i] / sum_height_tot;
#if 0 /* debug */
            fprintf( stderr, "%2d %2d %6.3f %6.3f  %6.3f * %6.3f = %6.3f %2d\n",
                     i, i, norm_spac_old[i], norm_spac_new[i], temp, weight, 
                     temp*weight, count[i] );
#endif
            sum_tot += temp*weight;
        }
        /* fprintf( stderr, "sum=%f\n", sum_tot ); */
        return sum_tot;    
    }
}


/*******************************************************************************
 * Function: collect_some_stats
 * Purpose:  Accumulate statistics on normalized peak spacing
 *******************************************************************************
 */
int
collect_some_stats( Data *data, Options *options, BtkMessage *message,
                    Results *results )
{
    int i, chunk_num = 0;
    Peak *pk0, *pk1;
    Stats stats;
    Histo histo;
    Peak  **dip;
    double ave_spacing[MAX_CHUNKS], med_spacing[MAX_CHUNKS];
    int    count2[MAX_CHUNKS];
    int	n_dip = 0, ctr=0;
    int shift[NUM_COLORS] = {0, 0, 0, 0};

    /* Initialization */
    histo.count_ = 0;
    histo.x_to_bin_ = 0;
    histo.x_min_ = INF;
    histo.x_max_ = -1;
    stats.mean_ = -1;
    stats.count_ = 0;
    
    results->count_chunks = 0;
    if( !results ) { return SUCCESS; }

    /* Create a list of all peaks */
    if (bc_data_create_single_ordered_peak_list(data, shift, message) != SUCCESS) {
        sprintf(message->text, "Error creating single peak list\n");
        return ERROR;
    }

#if 0
    for( i=0; i<data->peak_list_len; i++ ) {
        Peak *pk = data->peak_list[i];
        pk->wiheight = pk->iheight; /* need this for create_DIP_list */
    }
#endif

    if ( create_DIP_list(data, 0, data->length, shift, &dip, &n_dip, options, 
        message)
         != SUCCESS ) 
        goto error;

    
    /* Make first pass over data to collect ave_spacing and med_spacing 
     */
    for ( i=1, pk0=dip[0], chunk_num=-1; i<n_dip; i++, pk0=pk1 ) {
        pk1 = dip[i];
        if ( pk1->ipos / SCANS_PER_CHUNK > chunk_num ) {
            if ( chunk_num>=MAX_CHUNKS-1 ) { break; }
            chunk_num++;
            if( chunk_num > 0 ) {
                /* store peak spacing stats from last chunk */
                count2[chunk_num-1]      = ctr;
                ave_spacing[chunk_num-1] = statsMean(   &stats );
                med_spacing[chunk_num-1] = histoMedian( &histo );
            }
            /* Reset all statisics accumulators */ 
            ctr = 0;
            statsReset( &stats );
            histoReset( &histo, 0.0, 5.0*gSpacEst0 ); /* min --> max */
        }
        /* Accumulate peak spacing statistics for this chunk */
        ctr++;
        /* (1) rv = peak spacing (not normalized)
         * (2) stats are accumulated according to "weight" given by iheight
         */
        statsUpdateNTimes( &stats, pk1->ipos - pk0->ipos, pk1->iheight );
        histoUpdateNTimes( &histo, pk1->ipos - pk0->ipos, pk1->iheight );
    }
    /* store some stats from last chunk */
    count2[chunk_num]      = ctr;
    ave_spacing[chunk_num] = statsMean(   &stats );
    med_spacing[chunk_num] = histoMedian( &histo );

#if 0
    fprintf( stderr, "some stats\n" );
    for ( i=0; i<=chunk_num; i++ ) {
        fprintf( stderr, "%2d %6.3f %6.3f %6.3f  %2d\n",
                i, ave_spacing[i], med_spacing[i],
                ave_spacing[i]-med_spacing[i], count2[i] );
    }
#endif


    /* Make 2nd pass over data to get stats on normalized peak spacing
     */
    for ( i=1, pk0=dip[0], chunk_num=-1; i<n_dip; i++, pk0=pk1 ) {
        pk1 = dip[i];
        if ( pk1->ipos / SCANS_PER_CHUNK > chunk_num ) {
            if ( chunk_num+1>=MAX_CHUNKS ) { break; }
            if ( count2[chunk_num+1] < MIN_PEAKS_PER_CHUNK ) { break; }
            chunk_num++;
        }
        /* Accumulate peak spacing statistics for this chunk:
         * (1) rv = peak spacing normalized by representative spacing
         * (2) stats for rv are accumulated by "weight" given by iheight
         */
        statsUpdateNTimes( &results->stats[chunk_num], 
                     (pk1->ipos - pk0->ipos)/ave_spacing[chunk_num],
                           pk1->iheight );
        histoUpdateNTimes( &results->histo[chunk_num], 
                     (pk1->ipos - pk0->ipos)/med_spacing[chunk_num],
                           pk1->iheight );
    }

    results->count_chunks = QVMAX( results->count_chunks, chunk_num+1 );

    free( dip );
    return SUCCESS;

 error:
    free( dip );
    return ERROR;
}

#endif
