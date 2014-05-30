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
#include <string.h>
#include <unistd.h>


#include "Btk_atod.h"
#include "tracepoly.h"


/*  1.7
 *
 *  context - Produce an output similar to the check program using on
 *               Tracetuner train files except this version should run on
 *               phred train files.
 */


/* The following used to be separate files, they are now incorporated 
 * into this file
*/
/* #include "Stats.h" */
/* #include "Hcube.h" */
/* #include "Histo.h" */
/* #include "median.h" */


#if 1
/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 * Here's where I begin incorporating what used to be the separate file: 
 * "Stats.h:"
 vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/
#ifndef STATS_H__
#define STATS_H__


/*******************************************************************************
 * SStats == Simple Statistics (accumulate sum and sum of squares)
 *           (potential numerical problems for variance calculation)
 * RStats == Robust Statistics (accumulate mean and variance directly)
 *           (This uses a recursion relation which is more numerically
 *           stable than for Sstats)
 * These structs are intended to be used as "objects"
 * i.e. the data elements are intended to be "private"
 * and access is intended to be through the routines provided here.
 *
 * The interfaces are essentially identical, so it should suffice to 
 * give an example only of the use of RStats:
 * 
 ******************************************************************************/


/*******************************************************************************
 * Declaration of Class SStats
 ******************************************************************************/

typedef struct
{
    double min_, max_, sum_, sum_sqr_;
    int count_;
} SStats;

void
sStatsReset( SStats* s );

void
sStatsUpdate( SStats* s, double x );

double
sStatsMin( SStats* s );

double
sStatsMax( SStats* s );

double
sStatsMean( SStats* s );

double
sStatsVar( SStats* s );

double
sStatsStdDev( SStats* s );

int
sStatsCount( SStats* s );

void
sStatsPrint( SStats* s, FILE* fp );


/*******************************************************************************
 * Declaration of Class RStats
 ******************************************************************************/

typedef struct
{
    double min_, max_, mean_, var_;
    int count_;
} RStats;

void
rStatsReset( RStats* s );

void
rStatsUpdate( RStats* s, double x );

double
rStatsMin( RStats* s );

double
rStatsMax( RStats* s );

double
rStatsMean( RStats* s );

double
rStatsVar( RStats* s );

double
rStatsStdDev( RStats* s );

int
rStatsCount( RStats* s );

void
rStatsPrint( RStats* s, FILE* fp );

#endif /* STATS_H__ */

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 * Here's where I end incorporating what used to be the separate file: 
 * "Stats.h:"
 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 * Here's where I begin incorporating what used to be the separate file: 
 * "Stats.c:"
 vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/

#if 0   /* these includes are needed only if this is a separate file */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Stats.h"
#endif

/*******************************************************************************
 * Conventions to emulate C++ classes in C:
 * E.g., to emulate a C++ class called "MyType" I use the conventions:
 *   - Two files:    MyType.h and MyType.c
 *   - In the file MyType.h is the declaration:  typedef struct { .... } MyType;
 *   - Methods are emulated by using functions whose names start with "myType".
 *     The 1st argument to the "method" is of type: MyType*
 *   - Data members ("instance variables") are not accessed directly
 *     (i.e., they are assumed to be "private").
 *     They are accessed through "methods".
 *     E.g., the "instance variable" (field of struct) called "number"
 *     is accessed using a function called "myTypeNumber()"
 *     The C compiler cannot enforce this convention as can the C++ compiler
 *     which is one reason I use the convention of putting a trailing
 *     underscore on the instance variable names.
 *     (It makes it easier to grep for any code which breaks the convention.)
 *   - Constructors and Destructors are emulated using functions with names
 *     like: myTypeAllocate() and myTypeFree(), or
 *           myTypeInit() if no resource allocation is done
 *           (and thus no freeing need be done)
 *       
 * This file (and the corresponding .h file) implements two related "classes"
 * The classes are essentially identical in function and interface.
 * They differ in the way in which they trade off numerical accuracy
 * vs. speed.
 * SStats == Simple Statistics (accumulate sum and sum of squares)
 *           (potential numerical problems for variance calculation)
 * RStats == Robust Statistics (accumulate mean and variance directly)
 *           (This uses a recursion relation which is more numerically
 *           stable than for Sstats, but it is more expensive)
 *
 * Because the interfaces are essentially identical, it will suffice to
 * demonstrate the use of RStats.
 * First declare an object of type RStats:
 *     RStats stats;
 * Make sure it is reset before using it:
 *     rStatsReset(&stats);
 * Inside a loop in which "x" represents a random variable (r.v.) for
 * which we wish to collect statistics:
 *     rStatsUpate(&stats,x);
 * After collecting statistics, we may use one of several "methods"
 * to get the (sample) min, max, mean, standard deviation, count,
 * or to print out the various statistics.
 *
 * The interface for RStats is as follows:
 *    rStatsReset()  --> sets all the statistical "accumulators" to zero
 *    rStatsUpdate() --> increments the various "accumulators"
 *    rStatsMin()    --> returns the min value of the r.v.
 *    rStatsMax()    --> returns the max value of the r.v.
 *    rStatsMean()   --> returns the mean value of the r.v.
 *    rStatsVar()    --> returns the variance of the r.v.
 *    rStatsStdDev() --> returns the standard deviation of the r.v.
 *    rStatsCount()  --> returns the "count" 
 *                       (how many times rStatsUpdate() was called since
 *                       an rStatsReset() was done)
 *    rStatsPrint()  --> prints min, max, mean, std. dev, count
 *
 * The interface for SStats is as follows:
 *    sStatsReset()  --> sets all the statistical "accumulators" to zero
 *    sStatsUpdate() --> increments the various "accumulators"
 *    sStatsMin()    --> returns the min value of the r.v.
 *    sStatsMax()    --> returns the max value of the r.v.
 *    sStatsMean()   --> returns the mean value of the r.v.
 *    sStatsVar()    --> returns the variance of the r.v.
 *    sStatsStdDev() --> returns the standard deviation of the r.v.
 *    sStatsCount()  --> returns the "count" 
 *                       (how many times sStatsUpdate() was called since
 *                       an sStatsReset() was done)
 *    sStatsPrint()  --> prints min, max, mean, std. dev, count
 *
 *
 ******************************************************************************/


#define MYMIN(A,B)      ( (A)<(B) ? (A) : (B) )
#define MYMAX(A,B)      ( (A)>(B) ? (A) : (B) )
#define MYSQR(X)        ( (X)*(X) )


/*******************************************************************************
 * Implementation of Class SStats
 ******************************************************************************/

void
sStatsReset( SStats* s )
{
    s->min_     =  DBL_MAX;
    s->max_     = -DBL_MAX;
    s->sum_     = 0.0;
    s->sum_sqr_ = 0.0;
    s->count_   = 0;
}

void
sStatsUpdate( SStats* s, double x )
{
    s->min_     = MYMIN( s->min_, x );
    s->max_     = MYMAX( s->max_, x );
    s->sum_     += x;
    s->sum_sqr_ += x*x;
    s->count_   += 1;
}

double
sStatsMin( SStats* s )
{
    return s->min_;
}

double
sStatsMax( SStats* s )
{
    return s->max_;
}

double
sStatsMean( SStats* s )
{
    if( s->count_ > 0 ) {
        return (s->sum_)/(double)(s->count_);
    } else {
        return 0.0;
    }
}

double
sStatsVar( SStats* s )
{ 
    if( s->count_ > 1 ) {
        double mean = sStatsMean(s), dcount=s->count_;
        return ( (s->sum_sqr_ - dcount*MYSQR(mean)) / (dcount-1.0) );
    }
    else {
        return 0.0;
    }
}

double
sStatsStdDev( SStats* s )
{
    return sqrt(sStatsVar(s));
}

int
sStatsCount( SStats* s )
{
    return s->count_;
}

void
sStatsPrint( SStats* s, FILE* fp )
{
    fprintf( fp,
         "(min,max)=(%25.16e,%25.16e) (mean,std)=(%25.16e,%25.16e) ctr=%d\n",
            sStatsMin(s), sStatsMax(s),  
             sStatsMean(s), sStatsStdDev(s), sStatsCount(s) );
}



/*******************************************************************************
 * Implementation of Class RStats
 ******************************************************************************/

void
rStatsReset( RStats* s )
{
    s->min_   =  DBL_MAX;
    s->max_   = -DBL_MAX;
    s->mean_  = 0.0;
    s->var_   = 0.0;
    s->count_ = 0;
}

void
rStatsUpdate( RStats* s, double x )
{
    double old_mean = s->mean_, n;
    s->min_    = MYMIN( s->min_, x );
    s->max_    = MYMAX( s->max_, x );
    s->count_ += 1;

    n = s->count_;
    s->mean_ = s->mean_*((n-1.0)/n) + x/n;	/* recursion for mean */
    if( s->count_>1 )                           /* recursion for variance */
    {
        double temp1 = (x-old_mean)/n, temp2 = x-s->mean_;
        s->var_ = s->var_*((n-2.0)/(n-1.0)) + MYSQR(temp1) 
                + MYSQR(temp2)/(n-1.0);
        /* var_ = var_*((n-2.0)/(n-1.0)) + sqr((x-old_mean)/n) 
                + sqr(x-mean_)/(n-1.0); */
    }
}

double
rStatsMin( RStats* s )
{
    return s->min_;
}

double
rStatsMax( RStats* s )
{
    return s->max_;
}

double
rStatsMean( RStats* s )
{
    return (s->mean_);
}

double
rStatsVar( RStats* s )
{
    return (s->var_);
}

double
rStatsStdDev( RStats* s )
{
    return sqrt(rStatsVar(s));
}

int
rStatsCount( RStats* s )
{
    return s->count_;
}

void
rStatsPrint( RStats* s, FILE* fp )
{
    fprintf( fp, 
          "(min,max)=(%25.16e,%25.16e) (mean,std)=(%25.16e,%25.16e) ctr=%d\n",
            rStatsMin(s), rStatsMax(s),  
            rStatsMean(s), rStatsStdDev(s), rStatsCount(s) );
}



/*******************************************************************************
 * The following code serves as both:
 * (1) An example of how to use the classes SStats and RStats, and
 * (2) A debugging test to see if they are implemented correctly
 * This is a stand alone driver which will compile and run as-is.
 * (No other .c files need be linked)
 ******************************************************************************/

/* #define DRIVER */  /* or -DDRIVER in Makefile */
#ifdef DRIVER

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "Stats.h"

/* From Numerical Recipes (some code commented out to avoid static variables) */
double gasdev( unsigned short seed16v[3] )
{
    /*double ran1(long *idum);*/
    /*static int iset=0;*/
    /*static double gset;*/
    double fac,rsq,v1,v2;
    
    /*if  (iset == 0) {*/
    do {
        v1=2.0*erand48(seed16v)-1.0;
        v2=2.0*erand48(seed16v)-1.0;
        rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    /*gset=v1*fac;*/
    /*iset=1;*/
    return v2*fac;
    /*} else {*/
    /*    iset=0;*/
    /*    return gset;*/
    /*}*/
}


void 
test_stats1(void)
{
    SStats         x_stats1;
    RStats         x_stats2;
    int exponent;
    
    for( exponent=0; exponent<9; exponent++ )
    {
        double mean = pow(10.0,exponent);
        int trials;
        fprintf( stdout, "exponent = %d\n", exponent );

        for( trials=0; trials<4; trials++ )
        {
            unsigned short seed16v[3];
            int t;
            
            fprintf( stdout, "trials = %d\n", trials );
            /* state for erand48() random number generator */
            seed16v[0]=0x0123; seed16v[1]=0x4567; seed16v[2]=0x89ab; 
            sStatsReset( &x_stats1 );
            rStatsReset( &x_stats2 );
            
            for( t=0; t<trials; t++ )
            {
                double x = mean + gasdev(seed16v);
                
                sStatsUpdate( &x_stats1, x );
                rStatsUpdate( &x_stats2, x );
            }
            sStatsPrint( &x_stats1, stdout );
            rStatsPrint( &x_stats2, stdout );
        }
    }
}
  

int
main( int argc, char* argv[] )
{
    /* printf( "(min,max)= %24.12e %24.12e\n", -DBL_MAX, DBL_MAX ); */
    test_stats1();
    return 0;
}

#endif          /* DRIVER */

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 * Here's where I end incorporating what used to be the separate file: 
 * "Stats.c:"
 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 * Here's where I begin incorporating what used to be the separate file: 
 * "Hcube.h:"
 vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/
#ifndef HCUBE_H__
#define HCUBE_H__

/* Either one of the following is needed to get typedef size_t */
/* #include <stdlib.h> */
/* #include <stddef.h> */



/*******************************************************************************
 * Class Hcube declaration
 ******************************************************************************/

typedef struct {
    void *data_;
    int num_elem_, num_dim_;
    /* need <stdlib.h> or <stddef.h> for typedef size_t */
    size_t sizeof_data_elem_;   
} Hcube;


void
hcubeAllocate( Hcube* hc, size_t sizeof_data_elem, 
               int num_elements, int num_dim );

void
hcubeFree( Hcube* hc );

int 
hcubeNumDim( const Hcube* hc );

int 
hcubeNumElem( const Hcube* hc );

int 
hcubeNumData( const Hcube* hc );

void*
hcubeData( Hcube* hc );

void*
hcubeRef( Hcube* hc, const int indices[] );


#define HCUBE_ALLOCATE( HCUBE_PTR, TYPE, NUM_ELEM, NUM_DIM ) \
        hcubeAllocate( HCUBE_PTR, sizeof(TYPE), NUM_ELEM, NUM_DIM )



/*******************************************************************************
 * Class HcubeIter declaration
 ******************************************************************************/

typedef struct {
    /* Nothing magic about the number 32.  But consider the following:
     * If num_elem_==2 (probably the smallest sized hypercube that one
     * would consider), then the size of the data for the hypercube
     * would be:
     * (2^32)*sizeof(element) >= 2^32 >= address space for a 32-bit machine.
     * So, ..., 32 seems like a realistic upper bound on the size of
     * the array "indices_".
     * Before running into a problem imposed by the number 32 below,
     * you will run into the more stringent constraints imposed by 
     * the machine address space
     * (and also the maximum int allowed by the machine).
     */
    int indices_[32];
    int num_elem_, num_dim_;
} HcubeIter;


void
hcubeIterReset( HcubeIter* hci );

void
hcubeIterInit( HcubeIter* hci, Hcube* hc );

const int * 
hcubeIterIndices( HcubeIter* hci );

int
hcubeIterIncrement( HcubeIter* hci );


#endif          /* HCUBE_H__ */


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 * Here's where I end incorporating what used to be the separate file: 
 * "Hcube.h:"
 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 * Here's where I begin incorporating what used to be the separate file: 
 * "Hcube.c:"
 vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/
#if 0   /* these includes are needed only if this is a separate file */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "Hcube.h"
#endif

/*******************************************************************************
 * Conventions to emulate C++ classes in C:
 * E.g., to emulate a C++ class called "MyType" I use the conventions:
 *   - Two files:    MyType.h and MyType.c
 *   - In the file MyType.h is the declaration:  typedef struct { .... } MyType;
 *   - Methods are emulated by using functions whose names start with "myType".
 *     The 1st argument to the "method" is of type: MyType*
 *   - Data members ("instance variables") are not accessed directly
 *     (i.e., they are assumed to be "private").
 *     They are accessed through "methods".
 *     E.g., the "instance variable" (field of struct) called "number"
 *     is accessed using a function called "myTypeNumber()"
 *     The C compiler cannot enforce this convention as can the C++ compiler
 *     which is one reason I use the convention of putting a trailing
 *     underscore on the instance variable names.
 *     (It makes it easier to grep for any code which breaks the convention.)
 *   - Constructors and Destructors are emulated using functions with names
 *     like: myTypeAllocate() and myTypeFree(), or
 *           myTypeInit() if no resource allocation is done
 *           (and thus no freeing need be done)
 *       
 * This file (and the corresponding .h file) implements two related "classes"
 *
 * (1) Class Hcube implements a "hypercube" of data
 *     The data elements may be builtin types (e.g., int), or structs.
 *     The user declares an object of type "Hcube" as follows:
 *         Hcube h;
 *     The Hcube object 'h' must have it's internal data allocated 
 *     and freed as follows:
 *         HCUBE_ALLOCATE( &h, <Type>, <num_elem>, <num_dim> );
 *         hcubeFree( &h );
 *     Conceptually, the object 'h' is an array with num_dim dimensions.
 *     The size in each dimension is the same (=num_elem).
 *     The total number of elements is thus num_elem^num_dim
 *     E.g., a hypercube "h" with Type=double, num_elem=4 and num_dim=3 is 
 *     (conceptually) the same as the C declaration:
 *         double h[4][4][4];
 *     To get or set an element of the C-style array, 
 *     one could write the following:
 *         temp1 = h[ii][jj][kk]; 
 *         h[ii][jj][kk]= temp2;
 *
 *     Note that (in C) the data is stored such that:
 *         &(h[ii][jj][kk])+1 == &(h[ii][jj][kk+1])
 *     
 *     To get (or set) an element of a similar hypercube "hc" one would
 *     do the following: 
 *     (A) Declare an array of ints, e.g.,
 *             int indices[3];
 *     (B) Set the values in "indices" appropriately,
 *         There are two possible conventions.
 *         The (compile time) choice is determined by the
 *         macro "INDICES_CONVENTION" below.
 *         For the example given above:
 *         INDICES_CONVENTION==1 implies:
 *             indices[2]=ii; indices[1]=[jj]; indices[0]=kk;
 *         INDICES_CONVENTION==0 implies:
 *             indices[0]=ii; indices[1]=[jj]; indices[2]=kk;
 *     (C) Call hcubeRef:
 *         double* ref = (double*)hcubeRef(&hc,indices);
 *     (D) Dereference "ref" to either get or set the value:
 *             temp1 = *ref;
 *             *ref = temp2;
 *
 *     The "methods" are:
 *         hcubeAllocate()   --> you need to call this before using an Hcube
 *         HCUBE_ALLOCATE()  --> macro interface to hcubeAllocate
 *         hcubeFree()       --> youe need to call this when done with an Hcube
 *         hcubeNumDim()     --> returns num_dim
 *         hcubeNumElem()    --> returns num_elem
 *         hcubeNumData()    --> returns num_elem^num_dim
 *         hcubeData()       --> returns a pointer to the 1st data element
 *         hcubeRef()        --> returns a pointer to an arbitrary data element
 *
 * (2) Class HcubeIter implements an "iterator" to an Hcube object.
 *     (This is a common C++ construct, as used e.g., in the 
 *     C++ Standard Template Library (STL))
 *     This class is basically just a wrapper around an array of indices
 *     which is used as described above for the Hcube class.
 *
 *     One would declare an object of type HcubeIter as follows:
 *          HcubeIter hci;
 *     Assuming the existance of an object "hc" of type Hcube, one would
 *     initialize "hci" as follows:
 *          HcubeIterInit( &hci, &hc );
 *     To increment the iterator:
 *          hcubeIterIncrement( &hci );
 *     To use the iterator call "hcubeIterIndices"
 *     This returns an int* which may be used in a call to "hcubeRef()", e.g.,
 *          double* ref = (double*)hcubeRef( &hc, hcubeIterIndices(&hci) );
 *
 *     Also, "method"  hcubeIterIncrement returns a boolean specifying whether 
 *     or not the data is exhausted, so that a useful idiom is:
 *
 *     do {
 *        <--- whatever --->
 *      } while( hcubeIterIncrement(&hci) );
 *
 *     The "method" hcubeIterReset( HcubeIter* hci ) may be called
 *     to set the indices array back to zero.
 *
 *     The "methods" are:
 *       hcubeIterInit      --> user needs to call this before using
 *       hcubeIterIndices   --> returns an int* which can be used in
 *                                a call to hcubeRef()
 *       hcubeIterReset     --> sets "indices" array to zero
 *       hcubeIterIncrement --> increments indices in what could be considered:
 *                              lexicographic order (analogy with C-style array)
 *                              (i.e., we march through memory monotonically)
 *
 ******************************************************************************/



/*******************************************************************************
 * Class Hcube implementation
 ******************************************************************************/



#define INDICES_CONVENTION      0       /* 0 makes more sense to me */


void
hcubeAllocate( Hcube*   hc, 
               size_t   sizeof_data_elem, 
               int      num_elements, 
               int      num_dim )
{
    assert( num_dim>0 );
    assert( num_dim<=32 );
    hc->sizeof_data_elem_ = sizeof_data_elem;
    hc->num_elem_         = num_elements;
    hc->num_dim_          = num_dim;
    hc->data_  = calloc( hcubeNumData(hc), sizeof_data_elem );
}


void
hcubeFree( Hcube* hc )
{
    if( hc->data_ ) { free( hc->data_ ); }
    hc->data_ = NULL;
}


int 
hcubeNumDim( const Hcube* hc )
{
    return hc->num_dim_;
}


int 
hcubeNumElem( const Hcube* hc )
{
    return hc->num_elem_;
}


int 
hcubeNumData( const Hcube* hc )
{
    int n, prod=0;
    for( n=1,prod=1; n<=hc->num_dim_; n++ ) {
        prod *= hc->num_elem_;
    }
    return prod;        /* (hc->num_elem_)^(hc->num_dim_) */
}


void* 
hcubeData( Hcube* hc )
{
    return hc->data_;
}


void*
hcubeRef( Hcube*       hc, 
          const int    indices[] )
{
    char* ptr = (char*)hc->data_;
    int n, ctr;
/* Example to help out below, think of a C-style array:
 *     int a[ii][jj][kk];
 */
#if INDICES_CONVENTION
    /* indices[2]=ii; indices[1]=jj; indices[0]=kk; */
    for( n=hc->num_dim_-2,ctr=indices[hc->num_dim_-1]; n>=0; n-- )
#else
    /* indices[0]=ii; indices[1]=jj; indices[2]=kk; */
    for( n=1,ctr=indices[0]; n<=hc->num_dim_-1; n++ )
#endif
    {
        ctr = (ctr*hc->num_elem_)+indices[n];
    }
    return ptr + ctr*hc->sizeof_data_elem_;
}



/*******************************************************************************
 * Class HcubeIter implementation
 ******************************************************************************/

void
hcubeIterReset( HcubeIter* hci )
{
    int i;
    for( i=0; i<hci->num_dim_; i++ ) { hci->indices_[i]=0; }
}


void
hcubeIterInit( HcubeIter* hci, Hcube* hc )
{
    hci->num_elem_ = hc->num_elem_;
    hci->num_dim_  = hc->num_dim_;
    hcubeIterReset( hci );
}


const int * 
hcubeIterIndices( HcubeIter* hci )
{
    return hci->indices_;
}


int
hcubeIterIncrement( HcubeIter* hci )
{
    int n, more_data=0;

/* Example to help out below, think of a C-style array:
 *     int a[ii][jj][kk];
 */
#if INDICES_CONVENTION
    /* indices[2]=ii; indices[1]=jj; indices[0]=kk; */
    for( n=0; n<hci->num_dim_; n++ )
#else
    /* indices[0]=ii; indices[1]=jj; indices[2]=kk; */
    for( n=hci->num_dim_-1; n>=0; n-- )
#endif
    {
        if( hci->indices_[n] < hci->num_elem_-1 ) {
            hci->indices_[n]++;
            more_data = 1;
            break;
        } else {
            hci->indices_[n] = 0;
        }
    }
    /* return boolean indicating whether or not we've hit the end */
    return more_data;
}



/*******************************************************************************
 * The following code serves as both:
 * (1) An example of how to use the classes Hcube and HcubeIter, and
 * (2) A debugging test to see if they are implemented correctly
 * This is a stand alone driver which will compile and run as-is.
 * (No other .c files need be linked)
 ******************************************************************************/

#if 0
#define DRIVER  /* or -DDRIVER in Makefile */
#endif

#ifdef DRIVER

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "Hcube.h"

typedef struct {
    double x,y;
    char pad[27];       /* just for the heck of it to test things out */
} Point;


void
make_consecutive_points( Hcube* hc )
{
    Point* data = (Point*)hcubeData(hc);
    Point* temp;
    int size = hcubeNumData( hc );
    int i;

    fprintf( stdout, "data='%8p'\n", data ); 
    for( i=0; i<size; i++ ) {
        temp = data+i;
        temp->x = 10*i;
        temp->y = 20*i;
        fprintf( stdout, "i=%4d,  temp='%8p' %8.2f   %8.2f \n", 
                 i, temp, temp->x, temp->y ); 
    }
}


void
print_point( Hcube* hc, const int indices[] )
{
    Point* point = (Point*)hcubeRef( hc, indices );
    int n;

    for( n=hcubeNumDim(hc)-1; n>=0; n-- ) {
        fprintf( stdout, "%3d  ", indices[n] );
    }
    fprintf( stdout, "point='%8p' %8.2f   %8.2f\n",
             point, point->x, point->y );
}


void
print_all_points( Hcube* hc )
{
    HcubeIter hci;

    hcubeIterInit( &hci, hc );

    do {
        print_point( hc, hcubeIterIndices(&hci) );
    } while( hcubeIterIncrement(&hci) );

}


void
test_hcube( int max_num_dim )
{
    Hcube hcubes[32];
    int n, num_elem;
    int min_dim=1;
    num_elem = 2;

    for( n=min_dim; n<=max_num_dim; n++ ) {
        HCUBE_ALLOCATE( &(hcubes[n]), Point, num_elem, n );
    }

    for( n=min_dim; n<=max_num_dim; n++ ) {
        make_consecutive_points( &(hcubes[n]) );        
    }


    for( n=min_dim; n<=max_num_dim; n++ ) {
        print_all_points( &(hcubes[n]) );
    }

    for( n=min_dim; n<=max_num_dim; n++ ) {
        hcubeFree( &(hcubes[n]) );
    }
}


int
main( void )
{
    fprintf( stderr, "version 2\n" );
    test_hcube(3);
    return 0;
}



#endif         /* DRIVER */

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 * Here's where I end incorporating what used to be the separate file: 
 * "Hcube.c:"
 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 * Here's where I begin incorporating what used to be the separate file: 
 * "Histo.h:"
 vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/

#ifndef HISTO_H__
#define HISTO_H__

#if 1
#define HISTO_SIZE 600
#define HISTO_X_TO_BIN 100.0
#else
#define HISTO_SIZE 6
#define HISTO_X_TO_BIN 5.0
#endif


typedef struct {
    int array_[HISTO_SIZE];
    int count_;
} Histo;


void
histoReset( Histo* h );

void
histoUpdate( Histo* h, double x );

int
histoCount( Histo* h );

double
histoInverseCDF( Histo* h, double p );

double 
histoMedian( Histo* h );

double
histoMinusSigma( Histo* h );

double 
histoPlusSigma( Histo* h );

#endif  /* HISTO_H__ */

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 * Here's where I end incorporating what used to be the separate file: 
 * "Histo.h:"
 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 * Here's where I begin incorporating what used to be the separate file: 
 * "Histo.c:"
 vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/
#if 0   /* these includes are needed only if this is a separate file */
#include <assert.h>
#include <stdio.h>

#include "Histo.h"
#endif

/*******************************************************************************
 * Histogram "class".
 *
 * Assumptions:
 *   rv = random variable
 *   0 <= rv <= HISTO_SIZE/HISTO_X_TO_BIN
 * If rv is outside this range it is "clamped" to put in the histogram.
 * I.e., 
 *       if(rv < 0) { rv=0; }
 *       if(rv > HISTO_SIZE/HISTO_X_TO_BIN) { rv=HISTO_SIZE/HISTO_X_TO_BIN; }
 * The resolution (size of the bins) is given by HISTO_BIN_TO_X
 *
 * Methods: (1st argument of type Histo* is ignored below)
 *     histoReset()      --> reset histogram bins and counter to zero
 *     histoUpdate(rv)   --> update internal state with new random variable
 *                           (i.e., increment appropriate bin)
 *     histoCount()      --> return total count so far
 *     histoInverseCDF() --> histoInverseCDF(0.5) returns median value
 *                       --> histoInverseCDF(0.1587) returns -sigma (for normal)
 *                       --> histoInverseCDF(0.8413) returns +sigma (for normal)
 *     histoMedian()     --> returns histoInverseCDF(0.5)
 *     histoMinusSigma() --> returns histoInverseCDF(0.1587)
 *     histoPlusSigma()  --> returns histoInverseCDF(0.8413)
 *     histoPrint()      --> prints out internal state (mainly for debugging)
 ******************************************************************************/

#define HISTO_BIN_TO_X (1.0/HISTO_X_TO_BIN)

#define ROUNDPOS(X)     ( (int)((X)+0.5) )
#define ROUND(X)        ( ((X)>0.0) ? ROUNDPOS(X) : -ROUNDPOS(-(X)) )
#define MIN(A,B)        ( ((A)<(B)) ? (A) : (B) )
#define MAX(A,B)        ( ((A)>(B)) ? (A) : (B) )


void
histoReset( Histo* h )
{
    int i;
    for( i=0; i<HISTO_SIZE; i++ ) { h->array_[i] = 0; }
    h->count_ = 0;
}


void
histoUpdate( Histo* h, double x )
{
    double temp = x*HISTO_X_TO_BIN - 0.5;
    int bin = ROUND(temp);
    bin = MAX( 0, bin );
    bin = MIN( HISTO_SIZE-1, bin );
    h->array_[bin]++;
    h->count_++;
}


void
histoPrint( Histo* h, FILE* fp )
{
    int i;
    for( i=0; i<HISTO_SIZE; i++ ) {
        fprintf( fp, "   %f %d", (double)i*HISTO_BIN_TO_X, h->array_[i] );
    }    
    fprintf( stderr, "\n" );
}


double
histoInverseCDF( Histo* h, double p )
{
    double v0, v1, vtest, d0, d1, retval;
    int i0, i1;

    assert( p>=0 );
    assert( p<=1.0 );
    vtest = p*h->count_;

    if( h->count_ == 0 ) { return 0.0; }        /* somewhat arbitrary */

    if( p<0.5 ) {
        /*if( 0 ) { */  /* for debug */
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
        /* fprintf( stderr, "v0=%f i0=%d vtest=%f\n", v0, i0, vtest ); */
        do {
            i0--;
            v1 = v0;
            v0 -= h->array_[i0];
            /*fprintf( stderr, 
             * "v1=%f v0=%f i0=%d  array=%d\n", v1, v0, i0, h->array_[i0] );
             */
        } while( v0 >= vtest );
        i1 = i0 + 1;
        d0 = i0;
        d1 = i1;
    }

    retval = ( d0*(v1-vtest) + d1*(vtest-v0) ) / ((v1-v0)*HISTO_X_TO_BIN);    
    return retval;
}


double histoMedian( Histo* h )
{
    return histoInverseCDF( h, 0.5 );
}


double histoMinusSigma( Histo* h )
{
    return histoInverseCDF( h, 0.1587 );
}


double histoPlusSigma( Histo* h )
{
    return histoInverseCDF( h, 0.8413 );
}


int
histoCount( Histo* h )
{
    return h->count_;
}


#if 0
#define DRIVER   /* or -DDRIVER in Makefile */
#endif

#ifdef DRIVER

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void
test1(void)
{
    Histo h;
    const int size = 5;
    int i;
    double x;

    histoReset( &h );
    for( x=0.0; x<1.01; x+=0.25 ) {
        histoUpdate( &h, MIN(0.99999,x) );
    }
    histoUpdate( &h, 0 );
    /* histoPrint( &h, stderr ); */

    fprintf( stderr, "count= %d\n", histoCount(&h) );

    for( i=0; i<size; i++ ) 
    {
        double p =  (double)i/(double)(size-1);
        /* double p = 1.0; */
        double x = histoInverseCDF( &h, p );
        /* fprintf( stderr, "p=%g  p-1=%g\n", p, p-1.0 ); */
        fprintf( stderr, "   %4.2f %6.3f\n", p, x );
    }

}

/* From Numerical Recipes (some code commented out to avoid static variables) */
double gasdev( unsigned short seed16v[3] )
{
    /*double ran1(long *idum);*/
    /*static int iset=0;*/
    /*static double gset;*/
    double fac,rsq,v1,v2;
    
    /*if  (iset == 0) {*/
    do {
        v1=2.0*erand48(seed16v)-1.0;
        v2=2.0*erand48(seed16v)-1.0;
        rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    /*gset=v1*fac;*/
    /*iset=1;*/
    return v2*fac;
    /*} else {*/
    /*    iset=0;*/
    /*    return gset;*/
    /*}*/
}


void
test2(void)
{
    Histo h;
    double x;
    unsigned short seed16v[3];
    int trials=1000000, t;
    double mean = 2.5;

    seed16v[0]=0x0123; seed16v[1]=0x4567; seed16v[2]=0x89ab; 
    histoReset(&h);

    for( t=0; t<trials; t++ )
    {
        double x = mean + gasdev(seed16v);
        histoUpdate( &h, x );
    }
    
    {
        double med, msig, psig;
        med = histoMedian(&h);
        msig = histoMinusSigma(&h);
        psig = histoPlusSigma(&h);
        fprintf( stderr, " msig,median,psig = %f %f %f \n", 
                 msig, med, psig );
    }  
}

int
main(void)
{
    fprintf( stderr, "version 1\n" );
    test2();
    return 0;
}

#endif          /* DRIVER */

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 * Here's where I end incorporating what used to be the separate file: 
 * "Histo.c:"
 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 * Here's where I begin incorporating what used to be the separate file: 
 * "median.h:"
 vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/
#ifndef MEDIAN_H__
#define MEDIAN_H__

void 
heap_sort( int n, double ra[]);

double
array_median( int n, double array[], double temp[] );

#endif /*  MEDIAN_H__ */

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 * Here's where I end incorporating what used to be the separate file: 
 * "median.h:"
 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 * Here's where I begin incorporating what used to be the separate file: 
 * "median.c:"
 vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/
#if 0   /* these includes are needed only if this is a separate file */
#include <assert.h>

#include "median.h"
#endif

/*******************************************************************************
 * Heapsort algorithm from Numerical Recipes 2nd edition.
 * Modified so as to accept 0-based array of doubles.
 ******************************************************************************/

void 
heap_sort( int n, double ra_temp[] ) /* ra_temp assumed to be 0-based */
{
    int i,ir,j,l;
    double rra, *ra;

    ra = ra_temp - 1;   /* convert to 1-based array for N.R. */

    if (n < 2) return;
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
        if (l > 1) {
            rra=ra[--l];
        } else {
            rra=ra[ir];
            ra[ir]=ra[1];
            if (--ir == 1) {
                ra[1]=rra;
                break;
            }
        }
        i=l;
        j=l+l;
        while (j <= ir) {
            if (j < ir && ra[j] < ra[j+1]) j++;
            if (rra < ra[j]) {
                ra[i]=ra[j];
                i=j;
                j <<= 1;
            } else j=ir+1;
        }
        ra[i]=rra;
    }
}


double
array_median( int n, double array[], double temp[] )
{
    int i;
    assert( n>0 );
    for( i=0; i<n; i++ ) {
        temp[i] = array[i];
    }
    heap_sort( n, temp );
    if( n%2 == 0 ) {
        int i1,i2;
        i2 = n/2;
        i1 = i2-1;
        return 0.5 * ( temp[i1] + temp[i2] );
    } else {
        int ii;
        ii = n/2;
        return temp[ii];
    }
}

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 * Here's where I end incorporating what used to be the separate file: 
 * "median.c:"
 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 * This is the beginning of the original "checkheight.c" file
 vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/

#endif

#define MAXLINE		(1000)


int ACGT_to_int[256];
static const char Int_to_ACGT[4] = { 'A', 'C', 'G', 'T' };
static const int Not_ACGT = -1000;


static int
new_getbase( FILE* input_fp, int *spos, char* schar, int* cpos, char* cchar, 
         double *phr3, double *phr7, double *psr7, double *pres, 
         double* iheight, int *is_match, int *set_filename, char* filename )
{
    char linebuf[MAXLINE], *s;
    double m;
    int debug = 0;

    *set_filename = 0;
    *is_match = 0;  /* This handles the case where there are no parameters */
    if(debug) printf( "entering getbase\n" );
    for (;;) {
	if (fgets(linebuf, MAXLINE, input_fp) == NULL) {
	    return 0;
	}

	/* Ignore all white space lines */
	if (strspn(linebuf, " \t\r\n") == strlen(linebuf)) {
	    continue;
	}

        /* Extract filename from comment if present */
        if(linebuf[0] == '#') {

            s = strtok(linebuf, " \t\n"); /* grab '#' */
            if( s==NULL ) { continue; } /* This can't happen */

            s = strtok(NULL, " \t\n");
            if( s==NULL ) { continue; }

            if(    !strcmp(s,"File:") 
                || !strcmp(s,"file:") 
                || !strcmp(s,"FILE:") ) {
                s = strtok(NULL, " \t\n");
                if( s!=NULL ) {
                    strcpy( filename, s );
                    *set_filename = 1;
                    return 1;   /* NOTICE */
                }
                continue;
            }
        }

        /* Ignore all other comments */
	if ((linebuf[0] == '#')
            || ((linebuf[0] == '/') && (linebuf[1] == '*'))
            || (linebuf[0] == ';')) {
	    continue;
	}

	/* Get rid of 1) sample position, 2) sample base,
	 * 3) consensus position and 4) consensus base.
	 */
	s=strtok(linebuf, " \t\n");
	if ( s==NULL ) { continue; }
	*spos=atoi(s);

	s=strtok(NULL, " \t\n");
	if ( s==NULL ) { continue; }
	*schar=s[0];

	s=strtok(NULL, " \t\n");
	if ( s==NULL ) { continue; }
	*cpos=atoi(s);

	s=strtok(NULL, " \t\n");
	if ( s==NULL ) { continue; }
	*cchar=s[0];

	s+=strlen(s)+1;

	/* If sample == '-', there are no training parameters. */
	if ( (*schar) == '-' ) {            
            return 1;         
        }

	if ((Btk_atod(&s, &m)  != 1)
	|| (Btk_atod(&s, phr3) != 1)
	|| (Btk_atod(&s, phr7) != 1)
	|| (Btk_atod(&s, psr7) != 1)
	|| (Btk_atod(&s, pres) != 1)
	|| (Btk_atod(&s, iheight) != 1))
	{
	    return 0;
	}
	*is_match = (int)m;
	return 1;
    }
}


static int
get_next_read( FILE* input_fp,
               int spos_a[], char schar_a[],
               int cpos_a[], char cchar_a[], 
               double iheight_a[], int is_match_a[],
               char filename[], char filename_next[], int* num_bases )
{
    /* arguments to new_getbase() */
    double pres, phr3, phr7, psr7, iheight;
    int    spos,  cpos, is_match;
    char   schar, cchar;
    int    set_filename;

    *num_bases = 0;

    strcpy( filename, filename_next );
    while( new_getbase( input_fp, &spos, &schar, &cpos, &cchar, 
                        &phr3, &phr7, &psr7,
                        &pres, &iheight, &is_match,
                        &set_filename, filename_next ) )
    {
        if( set_filename ) {
            if( *num_bases ) {
                return 1;
            } else {
                *num_bases = 0;
                strcpy( filename, filename_next );
                continue;
            }
        }

        spos_a[*num_bases]     = spos;
        schar_a[*num_bases]    = schar;
        cpos_a[*num_bases]     = cpos;
        cchar_a[*num_bases]    = cchar;
        iheight_a[*num_bases]  = iheight;
        is_match_a[*num_bases] = is_match;

        (*num_bases)++; 
    }
    return ( (*num_bases>0) ? 1 : 0 );
}


/* Initialize the global variable "ACGT_to_int" */
static void
set_ACGT_to_int(void)
{
    int i;
    /* Set most of the array to some big number.
     * This is not strictly needed, 
     * it just enables a simple sanity check.
     * I.e., if anybody tries to use the big number, it will
     * immediately cause an "array out of bounds" (with purify)
     */
    for( i=0; i<256; i++ ) { ACGT_to_int[i]=Not_ACGT; }
    ACGT_to_int['A']=0; 
    ACGT_to_int['C']=1;   
    ACGT_to_int['G']=2;    
    ACGT_to_int['T']=3;
}


/*******************************************************************************
 * The next 2 routines deal with a FIFO data structure
 * (FIFO = First In First Out)
 ******************************************************************************/

static void
push_char_fifo( char fifo[], int size, char base_code )
{
    int i;
    for( i=size-1; i>0; i-- )    {  fifo[i] = fifo[i-1];  }
    fifo[0] = base_code;
}


static void
push_int_fifo( int fifo[], int size, int base_code )
{
    int i;
    for( i=size-1; i>0; i-- )    {  fifo[i] = fifo[i-1];  }
    fifo[0] = base_code;
}


static void
push_double_fifo( double fifo[], int size, double iheight )
{
    int i;
    for( i=size-1; i>0; i-- )    {  fifo[i] = fifo[i-1];  }
    fifo[0] = iheight;
}


#if 0
static double
average_iheight( double fifo[], int size )
{
    double sum = 0.0;
    int i;
    for( i=0; i<size; i++ ) {
        sum += fifo[i];
    }
    return sum / (double)size;
}
#endif


/*******************************************************************************
 * reset_stats() & reset_histos()
 * resets all the fields of the hypercube.
 * so they are ready to accumulate stats (mean & std. dev.) & histogram
 ******************************************************************************/
static void
reset_stats( Hcube* hc )
{
    HcubeIter hci;
    RStats* stats;

    hcubeIterInit( &hci, hc );
    do {
        stats = (RStats*)hcubeRef( hc, hcubeIterIndices(&hci) );
        rStatsReset( stats );
    } while( hcubeIterIncrement(&hci) );
}

static void
reset_histos( Hcube* hc )
{
    HcubeIter hci;
    Histo* histo;

    hcubeIterInit( &hci, hc );
    do {
        histo = (Histo*)hcubeRef( hc, hcubeIterIndices(&hci) );
        histoReset( histo );
    } while( hcubeIterIncrement(&hci) );
}


/*******************************************************************************
 * update_stats() & update_histo() updates the fields of the 
 * appropriate element of the hypercube 
 *    RStats ---> accumulate stats for mean & std. dev.
 *    Histo  ---> accumulate histogram
 ******************************************************************************/
static void
update_stats( Hcube* hc, int indices[], double value )
{
    RStats* stats = (RStats*)hcubeRef( hc, indices );
    rStatsUpdate( stats, value );
}

static void
update_histo( Hcube* hc, int indices[], double value )
{
    Histo* histo = (Histo*)hcubeRef( hc, indices );
    histoUpdate( histo, value );
}


/*******************************************************************************
 * print_stats() prints the mean and std. dev. of height statistic
 * for each of the combinations of sequences of ACGT in a hypercube
 ******************************************************************************/
static void
print_stats_and_histos( Hcube* stats_hc, Hcube* histo_hc, FILE* fp )
{
    int d, dim = hcubeNumDim(stats_hc), total=0;
    HcubeIter stats_hci;
    RStats rs_min, rs_14p, rs_med, rs_86p, rs_max, 
        rs_mean, rs_std, rs_std_prime, rs_count;
    Histo h_min, h_14p, h_med, h_86p, h_max, 
          h_mean, h_std, h_std_prime, h_count;

    assert( dim == hcubeNumDim(histo_hc) );
    assert( hcubeNumElem(stats_hc) == hcubeNumElem(histo_hc) );

    rStatsReset( &rs_min );
    rStatsReset( &rs_14p );
    rStatsReset( &rs_med );
    rStatsReset( &rs_86p );
    rStatsReset( &rs_max );
    rStatsReset( &rs_mean );
    rStatsReset( &rs_std );
    rStatsReset( &rs_std_prime );
    rStatsReset( &rs_count );

    histoReset( &h_min );
    histoReset( &h_14p );
    histoReset( &h_med );
    histoReset( &h_86p );
    histoReset( &h_max );
    histoReset( &h_mean );
    histoReset( &h_std );
    histoReset( &h_std_prime );
    histoReset( &h_count );


    /* print header */
    fprintf( fp, "\n\nContext: %1d base%c",
             dim, ( dim==1 ? ' ' : 's' ) );
    fprintf( fp, "\n" );
    //for( d=0; d<dim; d++ ) { fprintf( fp, "  " ); }
    //fprintf( fp  "                                      --- log base e ---\n" );

    if( dim==1 ) {
        fprintf( fp, "  " );
    } else {
        fprintf( fp, "5'" );
        for( d=1; d<dim-1; d++ ) { fprintf( fp, "  " ); }
        fprintf( fp, "3'" );
    }
    fprintf( fp, 
             "  min   14%%   med   86%%      max    mean   std  std'  count\n" );
    //%%   med   86%%\n" );
              
    hcubeIterInit( &stats_hci, stats_hc );
    do {
        const int* indices = hcubeIterIndices( &stats_hci );
        RStats* st = (RStats*)hcubeRef( stats_hc, indices );
        Histo*  hi = (Histo*) hcubeRef( histo_hc, indices );

        for( d=dim-1; d>=0; d-- ) {
        /*for( d=0; d<dim; d++ ) { */
            fprintf( fp, "%c ", Int_to_ACGT[indices[d]] );
        }
        {
            double
                percentile14 = histoMinusSigma(hi),
                median       = histoMedian(hi),
                percentile86 = histoPlusSigma(hi),
                mean         = rStatsMean(st),
                stddev       = rStatsStdDev(st),
                std_prime    = 0.5*(percentile86-percentile14),
                dmin         = rStatsMin(st),
                dmax         = rStatsMax(st),
                log_p14, log_med, log_p86;
            int
                hi_count     = histoCount(hi),
                st_count     = histoCount(hi);
            assert( hi_count == st_count );
            if( st_count == 0 ) {
                dmin = percentile14 = median = percentile86 = dmax 
                    = mean = stddev = std_prime = -1.0;
                log_p14 = log_med = log_p86 = 0.0;
            } else {    /* was log10() */
                log_p14 = log(percentile14);
                log_med = log(median);
                log_p86 = log(percentile86);
            }
            
            rStatsUpdate( &rs_min,   dmin );
            rStatsUpdate( &rs_14p,   percentile14 );
            rStatsUpdate( &rs_med,   median );
            rStatsUpdate( &rs_86p,   percentile86 );
            rStatsUpdate( &rs_max,   dmax );
            rStatsUpdate( &rs_mean,  mean  );
            rStatsUpdate( &rs_std,   stddev );
            rStatsUpdate( &rs_std_prime,   std_prime );
            rStatsUpdate( &rs_count, st_count );

            histoUpdate( &h_min,   dmin );
            histoUpdate( &h_14p,   percentile14 );
            histoUpdate( &h_med,   median );
            histoUpdate( &h_86p,   percentile86 );
            histoUpdate( &h_max,   dmax );
            histoUpdate( &h_mean,  mean  );
            histoUpdate( &h_std,   stddev );
            histoUpdate( &h_std_prime,   std_prime );
            histoUpdate( &h_count, st_count );

            fprintf( fp, 
//"%5.2f %5.2f %5.2f %5.2f %8.2f   %5.2f%6.2f %7d %5.2f %5.2f %5.2f  %5.2f\n",
"%5.2f %5.2f %5.2f %5.2f %8.2f   %5.2f%6.2f %5.2f %7d\n",
                     dmin, percentile14, median, percentile86, dmax,
                     mean, stddev, std_prime, st_count
                     //log_p14, log_med, log_p86,
                      );
            total += st_count;
        }
    } while( hcubeIterIncrement(&stats_hci) );

    for( d=0; d<dim; d++ ) { fprintf( fp, "--" ); }
    fprintf( fp, "------------------------------------------------------------\n" );

    /* print min statistics */
    if( dim==1 ) {
        fprintf( fp, "mn" );
    } else {
        fprintf( fp, "min:" );
    }
    for( d=2; d<dim; d++ ) { fprintf( fp, "  " ); }
    fprintf( fp, 
             "%5.2f %5.2f %5.2f %5.2f %8.2f   %5.2f%6.2f %5.2f %7d\n",
             rStatsMin(&rs_min),
             rStatsMin(&rs_14p),
             rStatsMin(&rs_med),
             rStatsMin(&rs_86p),
             rStatsMin(&rs_max),
             rStatsMin(&rs_mean),
             rStatsMin(&rs_std),
             rStatsMin(&rs_std_prime),
             (int)rStatsMin(&rs_count) );

    /* print 14% statistics */
    if( dim==1 ) {
        fprintf( fp, "14" );
    } else {
        fprintf( fp, "14%%:" );
    }
    for( d=2; d<dim; d++ ) { fprintf( fp, "  " ); }
    fprintf( fp, 
             "%5.2f %5.2f %5.2f %5.2f %8.2f   %5.2f%6.2f %5.2f\n",
             histoMinusSigma(&h_min),
             histoMinusSigma(&h_14p),
             histoMinusSigma(&h_med),
             histoMinusSigma(&h_86p),
             histoMinusSigma(&h_max),
             histoMinusSigma(&h_mean),
             histoMinusSigma(&h_std),
             histoMinusSigma(&h_std_prime)
        );
    /* (int)histoMinusSigma(&h_count) ); */

    /* print median statistics */
    if( dim==1 ) {
        fprintf( fp, "md" );
    } else {
        fprintf( fp, "med:" );
    }
    for( d=2; d<dim; d++ ) { fprintf( fp, "  " ); }
    fprintf( fp, 
             "%5.2f %5.2f %5.2f %5.2f %8.2f   %5.2f%6.2f %5.2f\n",
             histoMedian(&h_min),
             histoMedian(&h_14p),
             histoMedian(&h_med),
             histoMedian(&h_86p),
             histoMedian(&h_max),
             histoMedian(&h_mean),
             histoMedian(&h_std),
             histoMedian(&h_std_prime)
        );
    /* (int)histoMedian(&h_count) ); */

    /* print 86% statistics */
    if( dim==1 ) {
        fprintf( fp, "86" );
    } else {
        fprintf( fp, "86%%:" );
    }
    for( d=2; d<dim; d++ ) { fprintf( fp, "  " ); }
    fprintf( fp, 
             "%5.2f %5.2f %5.2f %5.2f %8.2f   %5.2f%6.2f %5.2f\n",
             histoPlusSigma(&h_min),
             histoPlusSigma(&h_14p),
             histoPlusSigma(&h_med),
             histoPlusSigma(&h_86p),
             histoPlusSigma(&h_max),
             histoPlusSigma(&h_mean),
             histoPlusSigma(&h_std),
             histoPlusSigma(&h_std_prime)
 );
    /* (int)histoPlusSigma(&h_count) ); */

    /* print max statistics */
    if( dim==1 ) {
        fprintf( fp, "mx" );
    } else {
        fprintf( fp, "max:" );
    }
    for( d=2; d<dim; d++ ) { fprintf( fp, "  " ); }
    fprintf( fp, 
             "%5.2f %5.2f %5.2f %5.2f %8.2f   %5.2f%6.2f %5.2f %7d\n",
             rStatsMax(&rs_min),
             rStatsMax(&rs_14p),
             rStatsMax(&rs_med),
             rStatsMax(&rs_86p),
             rStatsMax(&rs_max),
             rStatsMax(&rs_mean),
             rStatsMax(&rs_std),
             rStatsMax(&rs_std_prime),
             (int)rStatsMax(&rs_count) );

    /* print mean statistics */
    if( dim==1 ) {
        fprintf( fp, "av" );
    } else {
        fprintf( fp, "ave:" );
    }
    for( d=2; d<dim; d++ ) { fprintf( fp, "  " ); }
    fprintf( fp, 
             "%5.2f %5.2f %5.2f %5.2f %8.2f   %5.2f%6.2f %5.2f %7d\n",
             rStatsMean(&rs_min),
             rStatsMean(&rs_14p),
             rStatsMean(&rs_med),
             rStatsMean(&rs_86p),
             rStatsMean(&rs_max),
             rStatsMean(&rs_mean),
             rStatsMean(&rs_std),
             rStatsMean(&rs_std_prime),
             (int)rStatsMean(&rs_count) );

    /* print standard deviation statistics */
    if( dim==1 ) {
        fprintf( fp, "sd" );
    } else {
        fprintf( fp, "std:" );
    }
    for( d=2; d<dim; d++ ) { fprintf( fp, "  " ); }
    fprintf( fp, 
             "%5.2f %5.2f %5.2f %5.2f %8.2f   %5.2f%6.2f %5.2f %7d\n",
             rStatsStdDev(&rs_min),
             rStatsStdDev(&rs_14p),
             rStatsStdDev(&rs_med),
             rStatsStdDev(&rs_86p),
             rStatsStdDev(&rs_max),
             rStatsStdDev(&rs_mean),
             rStatsStdDev(&rs_std),
             rStatsStdDev(&rs_std_prime),
             (int)rStatsStdDev(&rs_count) );
    
    /* print "standard deviation" statistics */
    if( dim==1 ) {
        fprintf( fp, "s'" );
    } else {
        fprintf( fp, "sd':" );
    }
    for( d=2; d<dim; d++ ) { fprintf( fp, "  " ); }
    fprintf( fp, 
             "%5.2f %5.2f %5.2f %5.2f            %5.2f%6.2f %5.2f\n",
             0.5*(histoPlusSigma(&h_min)  - histoMinusSigma(&h_min)),
             0.5*(histoPlusSigma(&h_14p)  - histoMinusSigma(&h_14p)),
             0.5*(histoPlusSigma(&h_med)  - histoMinusSigma(&h_med)),
             0.5*(histoPlusSigma(&h_86p)  - histoMinusSigma(&h_86p)),
             //0.5*(histoPlusSigma(&h_max)  - histoMinusSigma(&h_max)),
             0.5*(histoPlusSigma(&h_mean) - histoMinusSigma(&h_mean)),
             0.5*(histoPlusSigma(&h_std)  - histoMinusSigma(&h_std)),
             0.5*(histoPlusSigma(&h_std_prime) - histoMinusSigma(&h_std_prime)) );
     
    /* total */
    fprintf( fp, "total:" );
    for( d=0; d<dim; d++ ) { fprintf( fp, "  " ); }
    fprintf( fp, "                                              %8d\n", total );
    
}


/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 * begin new stuff
 *vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/

#if 0
static double
rel_diff( double a, double b )
{
    double    sum   = a+b,     diff  = a-b;
    return ( diff==0.0 ? 0.0 : fabs(2.0*diff/sum) );

}
#endif


static void
tracePrintPeaks( ReadInfo* read_info, 
                 double xx[4][MaxNumBases], double yy[4][MaxNumBases],
                 int num_xy[4] )
{
    int i, b_code;
    static char* fnames[] = 
    { "out/a.dat", "out/c.dat", "out/g.dat", "out/t.dat" };

    for( b_code=0; b_code<4; b_code++ ) {
        FILE* fp = fopen( fnames[b_code], "w" );
        assert(fp);
        for( i=0; i<num_xy[b_code]; i++ ) {
            double 
                x1 = xx[b_code][i],
                y1 = yy[b_code][i],
                factor = readNormFactor( read_info, b_code, (int)x1 ),
                y2 = y1*factor,
                y3 = traceEvalFunc( &(read_info->read_polys[b_code]), x1 ),
                y4 = readEvalAveFunc( read_info, (int)x1 );
#if 0
            fprintf( fp, 
         "%3.0f %3.0f %3.0f %3.0f %3.0f %5.1f %5.1f %5.1f %5.1f %5.1f\n",
                     x1, x1, x1, x1, x1,
                     y1, y2, y3, y4, 100*factor );
#else
            y2=y2; y4=y4;
            fprintf( fp, 
         "%3.0f %3.0f %8.5f %8.5f\n",
                     x1, x1,
                     y1, y3 );            
#endif
        }
        fclose( fp );
    }

    for( b_code=0; b_code<4; b_code++ ) {
        FILE* fp = fopen( fnames[b_code], "a" );
        assert(fp);
        for( i=num_xy[b_code]-1; i>=0; i-- ) {
            double 
                x1 = 2.0*read_info->x_max_ - xx[b_code][i],
                y1 = yy[b_code][i],
                factor = readNormFactor( read_info, b_code, (int)x1 ),
                y2 = y1*factor,
                y3 = traceEvalFunc( &(read_info->read_polys[b_code]), x1 ),
                y4 = readEvalAveFunc( read_info, (int)x1 );
#if 0
            fprintf( fp, 
         "%3.0f %3.0f %3.0f %3.0f %3.0f %5.1f %5.1f %5.1f %5.1f %5.1f\n",
                     x1, x1, x1, x1, x1,
                     y1, y2, y3, y4, 100*factor );
#else
            y2=y2; y4=y4;
            fprintf( fp, 
         "%3.0f %3.0f %8.5f %8.5f\n",
                     x1, x1,
                     y1, y3 );            
#endif
        }
        fclose( fp );
    }

}

void
readPrintColorFunc( ReadPolynomials rp, FILE* fp, int base_code )
{
    int pos;
    TracePoly *tp = &(rp[base_code]); 

    for( pos=0; pos<=1000; pos+=20 ) {
        double val;
        if( pos <= tp->x_max_ ) {
            val = traceEvalFunc( tp, pos );
        } else { val = 0.0; }
        fprintf( fp, " %4.2f", val );
    }
    fprintf( fp, "\n" );
}


void
tracePrintXY( TracePoly* tp, FILE* fp,  double xx[], double yy[], int n  )
{  
    int i;
    //printf( "print_xy: nc=%d\n", nc );
    for( i=0; i<n; i++ ) {
        //if( i%10==0 ) { fprintf( fp, "\n" ); }
        fprintf( fp, "%3.0f %3.0f %5.1f %5.1f\n", 
                 xx[i], xx[i], yy[i], traceEvalFunc( tp, xx[i] ) );
    }
    //fprintf( fp, "\n" );
}


void
tracePrintNormBeginPoly( TracePoly* tp, FILE* fp )
{
    double height_begin = traceEvalFunc( tp, 0.0 );
    int pos;

    for( pos=0; pos<=1000; pos+=20 ) {
        double height;
        if( pos <= tp->x_max_ ) { height = traceEvalFunc( tp, pos ); }
        else { height = 0.0; }
        fprintf( fp, " %4.2f", height/height_begin );
    }
    fprintf( fp, "\n" );
}

#if 0
void
tracePrintNormColorFunc( ReadPolynomials rp, FILE* fp, int base_code )
{
    int pos;
    TracePoly
        *tp_ref = &(rp[4]),     /* NOTICE */
        *tp     = &(rp[base_code]);  /* NOTICE 0 vs 4 for ref. */

    for( pos=0; pos<=1000; pos+=20 ) {
        double height_ref, height, ratio;
        if( pos <= tp->x_max_ ) {
            ratio = traceEvalFunc( tp_ref, pos );
        } else { ratio = 0.0; }
        fprintf( fp, " %4.2f", ratio );
    }
    fprintf( fp, "\n" );
}
#endif

#if 0
static void
print_spos( int spos[], char schar[], int num_bases )
{
    int i;
    for( i=0; i<num_bases; i++ ) {
        printf( "%4d %c\n", spos[i], schar[i] );
    }
}
#endif


static int
move_data( int spos[], char schar[],
           double iheight[], int is_match[],
           int num_bases,
           double xx[4][MaxNumBases], double yy[4][MaxNumBases], 
           int num_xy[4] )
{
    int i, base_code, sum;
    for( base_code=0; base_code<4; base_code++ ) { num_xy[base_code] = 0; }

    for( i=0; i<num_bases; i++ ) {
        int indx;
        base_code = ACGT_to_int[(int)schar[i]];
        if( base_code == Not_ACGT ) { continue; }
        if( (is_match && !is_match[i]) ) { continue; }
        indx = num_xy[base_code];
        xx[base_code][indx] = spos[i];
        yy[base_code][indx] = iheight[i];
        num_xy[base_code]++;
    }
    sum = 0;
    for( base_code=0; base_code<4; base_code++ ) { sum += num_xy[base_code]; }
    return sum;
}


static void
height_stats( FILE* input_fp, int window_width, int dim_max, int verbose, 
                  int skip_length, double end_fraction, 
                  double threshold_debug, int xperimental )
{
    /* arguments to get_next_read() */
    double iheight[MaxNumBases];
    int    spos[MaxNumBases],  cpos[MaxNumBases], is_match[MaxNumBases];
    char   schar[MaxNumBases], cchar[MaxNumBases];
    char   filename[MAXLINE], filename_next[MAXLINE];
    int num_bases;

    int read_num=0, base_indx;

    /* other important state info (controls processing flow) */
    int end_state, dim;
    int num_matches_begin, num_matches_mid;
    int consec_matches;
    double ref_height;

    /* statistics gathering */
    int base_code_fifo[window_width];
    double iheight_fifo[window_width], temp[window_width];
    char base_char_fifo[window_width];

    /* debug info (getting rid of this stuff doesn't hurt) */
    int total_num_bases, 
        tot_matches, 
        tot_begin_matches, 
        tot_mid1_matches, tot_mid2_matches, 
        tot_end_matches,
        broken, violation;

    double xx[4][MaxNumBases], yy[4][MaxNumBases];
    int num_xy[4];

    Hcube   /* Ignore 0th element of arrays (treat as 1-based) */
        array_of_stats_cubes[dim_max+1], 
        array_of_histo_cubes[dim_max+1];
    {
        int i;
        for( i=0; i<window_width; i++ ) {
            base_char_fifo[i] = 'Z';
        }
    }
    for( dim=1; dim<=dim_max; dim++ ) {
        HCUBE_ALLOCATE( &(array_of_stats_cubes[dim]), RStats, 4, dim );
        reset_stats( &(array_of_stats_cubes[dim]) );
        HCUBE_ALLOCATE( &(array_of_histo_cubes[dim]), Histo, 4, dim );
        reset_histos( &(array_of_histo_cubes[dim]) );
    }

    /* debug initialization */
    tot_matches = tot_begin_matches = tot_mid1_matches
                = tot_mid2_matches  = tot_end_matches = 0;
    total_num_bases = 0;
    broken = 0;

    /* The following initialization (to arbitrary values)
     * is NOT needed. I'm doing only to stop compiler warnings.
     */
    end_state = violation = -1;
    num_matches_begin = num_matches_mid = consec_matches = -1;
    ref_height = -1.0;

    strcpy( filename_next, "" );        /* need to initialize */
    while( get_next_read( input_fp, spos, schar, cpos, cchar,
                          iheight, is_match, filename, filename_next,
                          &num_bases ) )

    {
        int ok_read, num_matched;
        ReadInfo read_info;
        int bad_read_num = -1; //103; //104; //328; //854;
#if TELL_ME_IF_I_EXCEED_THE_LIMITS
        /* For debug purposes, tracepoly.c wants these globals */
        extern char* Filename;
        extern int ReadNum, BadReadNum, ReportNum;
        Filename = filename;
        ReadNum = read_num+1;
        BadReadNum = bad_read_num;
        ReportNum = 0;
#endif
        read_num++;
        if (read_num % 10 == 0) { fprintf(stderr, "\rread#=%d ", read_num); }
        
        //if( total_num_bases > 710000 ) break;
        //if( read_num>913 ) break;
        //fprintf( stderr, "%d %s\n", read_num, filename );
        num_matched = move_data( spos, schar, iheight, is_match,
                                 num_bases, xx, yy, num_xy );

        ok_read = readGetCoefficients( &read_info, xx, yy, num_xy, read_num );
        if( !ok_read ) { continue; }

#if 1   /* This is useful debug stuff for printing arrays out to be plotted */
        if( xperimental && read_num==bad_read_num ) {
            fprintf( stderr, "read#=%d file=%s\n", read_num, filename );
            //print_spos( spos, schar, num_bases );
            tracePrintPeaks( &read_info, xx, yy, num_xy );
            //for( base_indx=0; base_indx<5; base_indx++ ) 
            //{ readPrintColorPoly( read_polys, stdout, base_indx ); }  
            readPrint( &read_info, stderr );
            break;
        }
#endif

        /* start of new read */
        end_state = 0;
        num_matches_begin = num_matches_mid = 0;
        consec_matches = 0;
        violation = 0;
        
        for( base_indx=0; base_indx<num_bases; base_indx++ ) {
            int base_code;
            double norm_iheight;

#if 0
            if (++total_num_bases % 10000 == 0) { 
                fprintf(stderr, "\r%d bases read so far...", 
                        total_num_bases);
            }
#endif
            /* For training, I wish to consider only unsaturated matches */
            /* NOTE: the following test should be consistent with the
             * one in the routine "readGetCoefficients"
             */
             if( !is_match[base_indx] || iheight[base_indx]>=1599 ) {
                consec_matches=0; 
                goto end_for;
            }  /* ------> */

            consec_matches++;
            tot_matches++;

            if( end_state ) {    /* Don't process rest of read */
                tot_end_matches++;
                goto end_for;                           /* ------> */
            }

            if( num_matches_begin <= skip_length ) { /* begin of read */
                num_matches_begin++;
                tot_begin_matches++;
                goto end_for;                           /* ------> */
            }

            /* Now we are processing matches in the "middle" of read 
             * We must wait "window_width-1" matches for FIFO to fill
             */
            num_matches_mid++;
            
            base_code = ACGT_to_int[(int)schar[base_indx]];

            norm_iheight = iheight[base_indx]
                * readNormFactor( &read_info, base_code, 
                                  spos[base_indx] );

            push_char_fifo( base_char_fifo, window_width, schar[base_indx] );
            push_int_fifo(  base_code_fifo, window_width, base_code );
            push_double_fifo( iheight_fifo, window_width, norm_iheight );
           
            if( num_matches_mid >= window_width ) {
                double ave, value;
                ave = array_median( window_width, iheight_fifo, temp);
                /* ave = average_iheight( base_height_fifo, size ); */
                // new_height *= weight_from_reverse_context( base_char_fifo );
                value = ave/norm_iheight;
                
                if( value > threshold_debug ) {
                    if( !violation ) {
                        fprintf( stderr, "read_num=%d %s\n", 
                                 read_num, filename );
                        printf( "read_num=%d %s\n", read_num, filename );
                    }
                    violation = 1;
                    fprintf( stderr, "%c %d  %d   %d   %f = %f / %f\n", 
                             schar[base_indx], spos[base_indx], cpos[base_indx], 
                             consec_matches, value, ave, norm_iheight );
                    printf( "%c %d  %d   %d   %f = %f / %f\n", 
                            schar[base_indx], spos[base_indx], cpos[base_indx], 
                            consec_matches, value, ave, norm_iheight );
                }
                if( num_matches_mid == window_width ) { 
                    ref_height = ave; 
                }
                if( ave < end_fraction*ref_height ) {
                    end_state = 1;          /* ignore rest of read */
                    tot_end_matches++;
                    goto end_for;                      /* ------> */
                } else {
                    tot_mid2_matches++;     /* FIFO is full */
                }
                for( dim=1; dim<=dim_max; dim++ ) {
                    if( consec_matches >= dim ) {
                        //printf( "%d\n", total_num_bases );    /* NOTICE */
                        if(verbose) { printf( "." ); }
                        update_stats( &(array_of_stats_cubes[dim]), 
                                      base_code_fifo, value );
                        update_histo( &(array_of_histo_cubes[dim]), 
                                      base_code_fifo, value );
                    }
                }
            } else {
                tot_mid1_matches++; /* Waiting for FIFO to fill */
            }
            
        end_for:
            {  /* This is just debug stuff. It is NOT needed otherwise */
                if( tot_matches != tot_begin_matches + tot_mid1_matches
                    + tot_mid2_matches + tot_end_matches ) {
                    broken = 1; 
                    break;
                }
            }
            
        } /* for ( base_indx=0; ... ) */
    } /* while( get_next_read() ) */

    if( dim_max>0 )    /* print out the overall statistics */
    {
        FILE* fp = stdout; /* fopen( "out/test.out", "w" ); */
        assert( fp );
        fprintf( fp, "\n\ntotal matches:\n" );
        fprintf( fp,"  begin    mid1    mid2     end     tot     sum\n");
        fprintf( fp, "%7d %7d %7d %7d %7d %7d\n",
                 tot_begin_matches, tot_mid1_matches, tot_mid2_matches, 
                 tot_end_matches, tot_matches,
                 tot_begin_matches+tot_mid1_matches
                 +tot_mid2_matches+tot_end_matches );
        if( broken ) {
            fprintf( fp,     "Something' wrong\n" );
            fprintf( stderr, "Something' wrong\n" );
            exit(-1);
        }
        for( dim=1; dim<=dim_max; dim++ ) {
            print_stats_and_histos( &(array_of_stats_cubes[dim]),  
                                    &(array_of_histo_cubes[dim]), fp );
        }
    }

    for( dim=1; dim<=dim_max; dim++ ) {
        hcubeFree( &(array_of_stats_cubes[dim]) );
        hcubeFree( &(array_of_histo_cubes[dim]) );
    }
        //fprintf( stdout, "read_num=%d\n", read_num );
        //print_all_coef( stdout, coef );
}


/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 * end new stuff
 *^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/


/*******************************************************************************
 * show_usage()
 ******************************************************************************/
void
show_usage(char *argv[],
           int width, int dim_max, int skip,
           double fraction_end, double threshold_debug, 
           int verbose, int xperimental )
{
    fprintf(stderr, "usage:\n %s ", argv[0] );
    fprintf(stderr, 
  "[-w <width>]\
   [-d <dim_max>]\
   [-s <skip>]\
   [-f <fract>]\
   [-t <threshold_debug>]\
   [-v] [-x] < <input file> > <output file>\n");
    fprintf(stderr, "or\n%s -h\n", argv[0] );

    fprintf( stderr,
"meaning of options and default values:\n\
-h  help (this usage message)\n\
-w  filter width for smoothing peaks heights (default=%d)\n\
-d  dimension_max for table (default=%d)\n\
-s  number of matched bases to skip at begin of read (default=%d)\n\
-f  stop processing read if ave_peak_height falls below this fract of\n\
    peak height at beginning of read (default=%g)\n\
-t  print out debug info if ratio: ave_peak_height / current_peak_height\n\
    exceeds this value (default=%g)\n\
-v  verbose (for debug only, not very useful right now) (default=%s)\n\
-x  experimental (for debug only -- no tables will be written out)\n\
    (default=%s)\n",
             width, dim_max, skip, fraction_end, threshold_debug,
             (verbose ?     "1 {true}" : "0 {false}" ), 
             (xperimental ? "1 {true}" : "0 {false}") );
}


/*******************************************************************************
 * main()
 ******************************************************************************/
int
main(int argc, char *argv[])
{
    int i, width, dim_max, verbose, skip, xperimental;
    double fraction_end, threshold_debug;

    //fprintf( stderr, "version 3\n" ); /*  debug */

    /* Set defaults. */
    width   = 10;
    dim_max = 8;        /* NOTICE */
    skip    = 0; //30;
    fraction_end    = 0.0; //0.5;
    threshold_debug = 1000000.0; /* DBL_MAX; */
    xperimental = 0;
    verbose = 0;

#if 0
    if (argc == 1) {
	show_usage(argv);
	exit(2);
    }
#endif

    while( (i = getopt(argc, argv, "w:d:s:f:t:vhx") ) != EOF ) {
        switch(i) {
        case 'w':
            width = atoi(optarg);
            if ( width <= 0 ) {
                fprintf(stderr, "Error: illegal width (-w %d)\n", 
                        width );
                exit(-1);
            }
            break;
        case 's':
            skip = atoi(optarg);
            if ( skip < 0 ) {
                fprintf(stderr, "Error: illegal skip (-s %d)\n", 
                        skip );
                exit(-1);
            }
            break;
        case 'f':
            fraction_end = atof(optarg);
            if ( fraction_end < 0.0 || fraction_end >1.0 ) {
                fprintf(stderr, "Error: illegal fraction_end (-f %f)\n", 
                        fraction_end );
                exit(-1);
            }
            break;
        case 't':
            threshold_debug = atof(optarg);
            if ( threshold_debug < 0.0 ) {
                fprintf(stderr, "Error: illegal threshold (-f %f)\n", 
                        threshold_debug );
                exit(-1);
            }
            break;
        case 'd':
            dim_max = atoi(optarg);
            /* 16 is somewhat arbitrary */
            if( dim_max <=0 || dim_max >= 16 ) { 
                fprintf(stderr, "Error: illegal dim_max (-d %d)\n", 
                        dim_max );
                exit(-1);
            }
            break;
	case 'v':
	    verbose = 1;
	    break;
	case 'x':
	    xperimental = 1;
	    break;
	case 'h':
	default:
	    show_usage(argv, width, dim_max, skip,fraction_end, threshold_debug,
                       verbose, xperimental );
	    exit(2);
	}
    }
    
    if( xperimental ) {
        dim_max = 0;
    } else {
        printf( "\
    width      = %d\n\
    dim_max    = %d\n\
    skip       = %d\n\
    fraction   = %f\n\
    threshold  = %g\n\
    verbose    = %c\n",
                width, dim_max, skip, fraction_end, 
                threshold_debug, (verbose ? 'T' : 'F' ) );
    }
    fflush( stdout );

    set_ACGT_to_int();/* Initialize global variable ACGT_to_int */

    height_stats( stdin, width, dim_max, verbose, skip,
                  fraction_end, threshold_debug, xperimental );

    return(0);
}
