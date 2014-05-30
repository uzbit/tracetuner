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
 * 2.6 2003/11/06 18:18:46
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>   
#include <float.h>

#include "Btk_qv.h"
#include "Btk_atod.h"
#include "context_table.h"

#define CONTEXT_DEBUG 0
#define MAXLINE       (1000)
#define CHUNK         (1000)

#if 0
#define CONTEXT_DRIVER 1  /* or -DCONTEXT_DRIVER in Makefile */
#include "Btk_atod.c"
#endif

/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
 * Here's where I begin incorporating what used to be the separate file:
 * "Hcube.h:"
 vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/

#ifndef HCUBE_H__
#define HCUBE_H__

/*******************************************************************************
 * Class Hcube declaration 
 *******************************************************************************/

typedef struct {
    void *data_;
    int my_data_, num_elem_, num_dim_;
    size_t sizeof_data_elem_; /*need <stdlib.h>or<stddef.h> for typedef size_t*/
} Hcube;


void
hcubeAllocate( Hcube* hc, size_t sizeof_data_elem, int num_elements, 
               int num_dim );

void
hcubeInit( Hcube* hc, size_t sizeof_data_elem, int num_elements, 
           int num_dim, void* data );

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

#define HCUBE_INIT( HCUBE_PTR, TYPE, NUM_ELEM, NUM_DIM, PTR ) \
        hcubeInit( HCUBE_PTR, sizeof(TYPE), NUM_ELEM, NUM_DIM, PTR )


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

void
hcubeIterInitElemDim( HcubeIter* hci, int elem, int dim );

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
/* #include <assert.h> */
/* #include <stdio.h> */
/* #include <stdlib.h> */

/* #include "Hcube.h" */


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
 *         hcubeAllocate()   --> You need to call one of these
 *         hcubeInit()       --> two routines before using an Hcube.
 *         HCUBE_ALLOCATE()  --> macro interface to hcubeAllocate
 *         HCUBE_INIT()      --> macro interface to hcubeAllocate
 *         hcubeFree()       --> you need to call this when done with an Hcube
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
 *     Methods: (1st argument of type Hcube* is ignored below)
 *            call one of the next two routines before using
 *       hcubeIterInit(Hcube* hc)  --> use Hcube to initialize
 *       hcubeIterInitElemDim(int elem, int dim) --> "manual" initialization
 *       hcubeIterIndices()   --> returns an int* which can be used in
 *                                a call to hcubeRef()
 *       hcubeIterReset()     --> sets "indices" array to zero
 *       hcubeIterIncrement() --> increments indices in what could be considered:
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
    hc->my_data_          = 1;
    hc->sizeof_data_elem_ = sizeof_data_elem;
    hc->num_elem_         = num_elements;
    hc->num_dim_          = num_dim;
    hc->data_  = calloc( hcubeNumData(hc), sizeof_data_elem );
}


void
hcubeInit( Hcube*   hc, 
           size_t   sizeof_data_elem, 
           int      num_elements, 
           int      num_dim,
           void*    data )
{
    assert( num_dim>0 );
    assert( num_dim<=32 );
    hc->my_data_          = 0;
    hc->sizeof_data_elem_ = sizeof_data_elem;
    hc->num_elem_         = num_elements;
    hc->num_dim_          = num_dim;
    hc->data_             = data;
}


void
hcubeFree( Hcube* hc )
{
    /* No free is needed if hcubeInit() was called */
    if( hc->my_data_ ) {
        /* No free is needed if we've already done so */
        if( hc->data_ ) { 
            free( hc->data_ );
        }
    }
    hc->data_ = NULL;   /* Record the fact that we've freed memory */
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


void
hcubeIterInitElemDim( HcubeIter* hci, int elem, int dim )
{
    hci->num_elem_ = elem;
    hci->num_dim_  = dim;
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



/*******************************************************************************
 * The following code serves as both:
 * (1) An example of how to use the function "weight_from_reverse_context()",
 * and
 * (2) A debugging test to see if it is implemented correctly
 * This is a stand alone driver which will compile and run as-is.
 * The only other .c file needed is "default_context.c"
 ******************************************************************************/

#ifdef CONTEXT_DRIVER   /* or -DCONTEXT_DRIVER in Makefile */
/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/

static char Int_to_ACGT[4] = { 'A', 'C', 'G', 'T' };

static void
print_weight_using_ints( Hcube* hc, const int indices[] )
{
    double* weight = (double*)hcubeRef( hc, indices );
    int n;

    for( n=hcubeNumDim(hc)-1; n>=0; n-- ) {
        fprintf( stdout, "%1c ", Int_to_ACGT[indices[n]] );
    }
    fprintf( stdout, " %4.2f\n", *weight );
}


static void
print_all_weights_using_ints( Hcube* hc )
{
    HcubeIter hci;

    hcubeIterInit( &hci, hc );

    do {
        print_weight_using_ints( hc, hcubeIterIndices(&hci) );
    } while( hcubeIterIncrement(&hci) );
}


/* This routine loops over the entire table and prints to stdout.
 * It uses simple indexing (no interpolation for mixed bases).
 */
static void
test_reverse_using_int_indices(char *context_table, ContextTable *ctable )
{
    Hcube hcube;
    //ContextTable *ctable;

    //ctable = read_context_table(context_table);
    HCUBE_INIT( &hcube, double, 4, ctable->dimension, ctable->weights);

    print_all_weights_using_ints( &hcube );

    hcubeFree( &hcube );
}


static void
print_reverse_weight_using_acgt( const char base_code[], ContextTable *ctable )
{
    int dim = ctable->dimension;
    double weight;
    int n;

    for( n=dim-1; n>=0; n-- ) {
        fprintf( stdout, "%c ", base_code[n] );
    }
    fprintf( stdout, "\n" );

    weight = weight_from_reverse_context( base_code, ctable );
    fprintf( stdout, " weight = %4.2f\n", weight );
}


static void
print_weight_using_acgt( const char base_code[], ContextTable *ctable )
{
    int dim = ctable->dimension;
    double weight;
    int n;

    for( n=dim-1; n>=0; n-- ) {
        fprintf( stdout, "%c ", base_code[-n] );
    }
    fprintf( stdout, "\n" );

    weight = weight_from_context( base_code, ctable );
    fprintf( stdout, " weight = %4.2f\n", weight );
}


#if 0
static char Int_to_NAC[] = { 'A', 'C', 'G', 'T',
                               'R', 'Y', 'K', 'M', 'S', 'W',
                               'B', 'D', 'H', 'V', 'N' };
static char Int_to_NAC[] = { 'G', 'Y', 'B' };
#else
static char Int_to_NAC[] = { 'A', 'C', 'G', 'T' };

#endif

/* This routine loops over a subset of the table and prints to stdout.
 * (The subset is defined by the above static array)
 * It does interpolation for mixed bases.
 * (Using call to "weight_from_reverse_context"
 */
#define NUM_ELEMENTS(A) (sizeof(A)/sizeof(A[0]))
static void
test_reverse_using_char_indices(ContextTable *ctable)
{
    HcubeIter hci;
    int d, dim = ctable->dimension;
    char char_context[32];
    
    /* We are getting set up so that every array index loops
     *    over the values in "Int_to_NAC".
     * The total number of times through the "do {} while()"
     *    loop below is  NUM_ELEMENTS(Int_to_NAC)^dim
     */
    hcubeIterInitElemDim( &hci, NUM_ELEMENTS(Int_to_NAC), dim );

    printf( "test reverse\n" );
    do {
        for( d=0; d<dim; d++ ) {
            char_context[d] = Int_to_NAC[hcubeIterIndices(&hci)[d]];
        }
        print_reverse_weight_using_acgt( char_context, ctable );
    } while( hcubeIterIncrement(&hci) );
}

static void
test_using_char_indices(ContextTable *ctable) 
{
    HcubeIter hci;
    int d, dim = ctable->dimension;
    char temp[32];
    char* char_context = temp+31;
    
    /* We are getting set up so that every array index loops
     *    over the values in "Int_to_NAC".
     * The total number of times through the "do {} while()"
     *    loop below is  NUM_ELEMENTS(Int_to_NAC)^dim
     */
    hcubeIterInitElemDim( &hci, NUM_ELEMENTS(Int_to_NAC), dim );

    printf( "test\n" );
    do {
        for( d=0; d<dim; d++ ) {
            char_context[-d] = Int_to_NAC[hcubeIterIndices(&hci)[d]];
        }
        print_weight_using_acgt( char_context, ctable );
    } while( hcubeIterIncrement(&hci) );
}


#if 0
static void
test_loop_over_nac(void)
{
    static const int dim = 5;
    static char context[] = { 'A', 'R', 'D', 'N', 'Y' };

    ContextIter ci;
    contextIterInit( &ci, context, dim );

    do {
        contextIterPrintCharIndices( &ci );
        printf( "\n" );
    } while( contextIterIncrement(&ci) );
}
#endif

int
main( int argc, char* argv[] )
{
    char* fname;
    ContextTable *ctable;
    assert( argc==2 );
    fname = argv[1];
    fprintf( stderr, "fname=%s\n", fname );
    ctable = read_context_table(fname);
#if 0
    set_NucleicAcidCode_to_ACGT();
    set_ACGT_to_int();
#endif
    fprintf( stderr, "version 1\n" );
    test_reverse_using_char_indices(ctable);
#if 0
    test_using_char_indices(ctable);
    test_reverse_using_int_indices(ctable);
    test_loop_over_nac();
#endif
    destroy_context_table( ctable );
    return 0;
}

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
#endif          /* DRIVER */


/*******************************************************************************
 * "NucleicAcidCode_to_ACGT" is a data structure encoding the mixed base info:
 * i.e., R=AG,  Y=CT,  K=GT,  M=AC,  B=CGT,  ...,  N=ACGT
 * "ACGT_to_int" is a mapping from A,C,G,T to 0,1,2,3
 *    (This info can be pulled out of "NucleicAcidCode_to_ACGT", but it's
 *     more cryptic)
 ******************************************************************************/

typedef struct {
    int num;
    char map[4];
}  NucleicAcidCode;

NucleicAcidCode NucleicAcidCode_to_ACGT[256];

static int  ACGT_to_int[256];


/*******************************************************************************
 * Initialize the global variable "ACGT_to_int"
 ******************************************************************************/

static void 
set_ACGT_to_int(void)
{
    int i;
    /* Set most of the array to some bogus number.
     * This is not strictly needed, 
     * it just enables a simple sanity check.
     * I.e., if anybody tries to use the bogus number, it will
     * immediately cause an "array out of bounds" (with purify)
     */
    for( i=0; i<256; i++ ) {
        ACGT_to_int[i] = -1;
    }
    ACGT_to_int['A'] = 0;
    ACGT_to_int['C'] = 1;
    ACGT_to_int['G'] = 2;
    ACGT_to_int['T'] = 3;
}



/*******************************************************************************
 * Initialize the global variable "NucleicAcidCode_to_ACGT"
 ******************************************************************************/

static void 
set_NucleicAcidCode_to_ACGT(void)
{
    int i,j;
    NucleicAcidCode* nac;
    /* Set most of the array to some bogus number.
     * This is not strictly needed, 
     * it just enables a simple sanity check.
     * I.e., if anybody tries to use the bogus number, it will
     * immediately cause an "array out of bounds" (with purify)
     */
    for( i=0; i<256; i++ ) {
        nac = &NucleicAcidCode_to_ACGT[i];
        nac->num = -1;  /* bogus num */
        for( j=0; j<4; j++ ) { nac->map[j] = 'Z'; } /* bogus char */
    }
    
    nac = &NucleicAcidCode_to_ACGT['A'];
    nac->num = 1;
    nac->map[0] = 'A';

    nac = &NucleicAcidCode_to_ACGT['C'];
    nac->num = 1;
    nac->map[0] = 'C';

    nac = &NucleicAcidCode_to_ACGT['G'];
    nac->num = 1;
    nac->map[0] = 'G';

    nac = &NucleicAcidCode_to_ACGT['T'];
    nac->num = 1;
    nac->map[0] = 'T';

    nac = &NucleicAcidCode_to_ACGT['R'];
    nac->num = 2;
    nac->map[0] = 'A';
    nac->map[1] = 'G';

    nac = &NucleicAcidCode_to_ACGT['Y'];
    nac->num = 2;
    nac->map[0] = 'C';
    nac->map[1] = 'T';

    nac = &NucleicAcidCode_to_ACGT['K'];
    nac->num = 2;
    nac->map[0] = 'G';
    nac->map[1] = 'T';

    nac = &NucleicAcidCode_to_ACGT['M'];
    nac->num = 2;
    nac->map[0] = 'A';
    nac->map[1] = 'C';

    nac = &NucleicAcidCode_to_ACGT['S'];
    nac->num = 2;
    nac->map[0] = 'C';
    nac->map[1] = 'G';

    nac = &NucleicAcidCode_to_ACGT['W'];
    nac->num = 2;
    nac->map[0] = 'A';
    nac->map[1] = 'T';

    nac = &NucleicAcidCode_to_ACGT['B'];
    nac->num = 3;
    nac->map[0] = 'C';
    nac->map[1] = 'G';
    nac->map[2] = 'T';

    nac = &NucleicAcidCode_to_ACGT['D'];
    nac->num = 3;
    nac->map[0] = 'A';
    nac->map[1] = 'G';
    nac->map[2] = 'T';

    nac = &NucleicAcidCode_to_ACGT['H'];
    nac->num = 3;
    nac->map[0] = 'A';
    nac->map[1] = 'C';
    nac->map[2] = 'T';

    nac = &NucleicAcidCode_to_ACGT['V'];
    nac->num = 3;
    nac->map[0] = 'A';
    nac->map[1] = 'C';
    nac->map[2] = 'G';

    nac = &NucleicAcidCode_to_ACGT['N'];
    nac->num = 4;
    nac->map[0] = 'A';
    nac->map[1] = 'C';
    nac->map[2] = 'G';
    nac->map[3] = 'T';

}


/*******************************************************************************
 * Declaration of ContextIter
 ******************************************************************************/

typedef struct {
    /* Nothing magic about the number 32.  But 4^32 is a very big number!! */
    int  int_indices_[32],  int_limits_[32], public_int_indices_[32];
    char char_indices_[32], char_init_[32], base_code_[32];
    int  num_dim_;
} ContextIter;

#define INDICES_CONVENTION      0       /* 0 makes more sense to me */

static void
contextIterReset( ContextIter* ci );

static void
contextIterInit( ContextIter* ci, const char base_code[], int dim );

static const int * 
contextIterIntIndices( ContextIter* ci );

static int
contextIterIncrement( ContextIter* ci );


#if 0
/*******************************************************************************
 * Class ContextIter implementation
 ******************************************************************************/

static void 
contextIterPrint( ContextIter* ci, FILE* fp )
{
    int i;
    int dim = ci->num_dim_;
    
    fprintf( fp, "base code         " );
    for(i=dim-1;i>=0;i--){ fprintf( fp, " %c", ci->base_code_[i] ); } 
    printf("\n");

    fprintf( fp, "char_init         " );
    for(i=dim-1;i>=0;i--){ fprintf( fp, " %c", ci->char_init_[i] ); } 
    printf("\n");

    fprintf( fp, "char_indices      " );
    for(i=dim-1;i>=0;i--){ fprintf( fp, " %c", ci->char_indices_[i] ); } 
    printf("\n");

    fprintf( fp, "int_limits        " );
    for(i=dim-1;i>=0;i--){ fprintf( fp, " %d", ci->int_limits_[i] ); } 
    printf("\n");

    fprintf( fp, "int_indices       " );
    for(i=dim-1;i>=0;i--){ fprintf( fp, " %d", ci->int_indices_[i] ); } 
    printf("\n");

    fprintf( fp, "public_int_indices" );
    for(i=dim-1;i>=0;i--){ fprintf( fp, " %d", ci->public_int_indices_[i] ); } 
    printf("\n");
}
#endif

static void
contextIterReset( ContextIter* ci )
{
    int i;
    char base;
    for( i=0; i<ci->num_dim_; i++ ) { 
        ci->int_indices_[i] = 0;
        base = ci->char_init_[i];
        ci->char_indices_[i] = base;
        ci->public_int_indices_[i] = ACGT_to_int[(int)base];
    }
}


static void
contextIterInit( ContextIter* ci, const char base_code[], int dim )
{
    int i;
    ci->num_dim_  = dim;
    for( i=0; i<ci->num_dim_; i++ ) { 
        ci->int_limits_[i] = NucleicAcidCode_to_ACGT[(int)base_code[i]].num;
        ci->base_code_[i] = base_code[i];
        ci->char_init_[i] =  NucleicAcidCode_to_ACGT[(int)base_code[i]].map[0];
    }
    contextIterReset( ci );
}


static const int * 
contextIterIntIndices( ContextIter* ci )
{
    return ci->public_int_indices_;
}

static int
contextIterIncrement( ContextIter* ci )
{
    int n, more_data=0;

/* Example to help out below, think of a C-style array:
 *     int a[ii][jj][kk];
 */
#if INDICES_CONVENTION
    /* indices[2]=ii; indices[1]=jj; indices[0]=kk; */
    for( n=0; n<ci->num_dim_; n++ )
#else
    /* indices[0]=ii; indices[1]=jj; indices[2]=kk; */
    for( n=ci->num_dim_-1; n>=0; n-- )
#endif
    { 
        /* int limit = NucleicAcidCode_to_ACGT[(int)ci->base_code_[n]].num; */
        if( ci->int_indices_[n] < ci->int_limits_[n]-1 ) 
        {
            ci->int_indices_[n]++;
            more_data = 1;
            break;
        } else {
            ci->int_indices_[n] = 0;
        }
    }
    for( n=0; n<ci->num_dim_; n++ ) {
        char acgt
     = NucleicAcidCode_to_ACGT[(int)ci->base_code_[n]].map[ci->int_indices_[n]];
        ci->char_indices_[n] = acgt;
        ci->public_int_indices_[n] = ACGT_to_int[(int)acgt];
        /*     = NucleicAcidCode_to_ACGT[(int)acgt].map[0]; */
    }
    
    /* return boolean indicating whether or not we've hit the end */
    return more_data;
}

/*******************************************************************************
 * Weight_from_reverse_context()
 * Inputs:
 *     base_code[0]   = current base
 *     base_code[1]   = previous base
 *          .
 *        (etc.)
 *          .
 *     base_code[n-1] = "oldest" base
 * where "base" may be any one of: A,C,G,T, R,Y,K,M,S,W, B,D,H,V, N
 ******************************************************************************/
double 
weight_from_reverse_context( const char base_code[], ContextTable *ctable )
{
    static int initialized = 0;
    static Hcube hcube;
    int dim = ctable->dimension;
    const int max_dim = 32; /* should be OK; 4^32 is a big number!! */
    typedef double EntryType;
    double sum;

    if( !initialized ) {
        /* initialize global variables */
        set_ACGT_to_int();
        set_NucleicAcidCode_to_ACGT();

        /* Initialize static variables */
        HCUBE_INIT( &hcube, EntryType, 4 /* 4 bases: ACGT */, dim, 
            ctable->weights );
        assert( hcubeNumDim( &hcube ) <= max_dim );
        initialized = 1;
    }
    {
        int count = 0;
        ContextIter ci;
        contextIterInit( &ci, base_code, dim );

        /*printf( "init:\n" );*/
        /*contextIterPrint( &ci, stdout );*/

        sum = 0.0;
        do {
            const int* i_indices = contextIterIntIndices(&ci);
            double val = *(EntryType*)hcubeRef( &hcube, i_indices );
            if( val < 0 ) {
                val = 1.0;
                fprintf( stderr, 
            "Missing entry in context table (assume 1.0 for now)\n"
            "Better Solution: More training data or smaller context length\n" );
            }
                sum += val;
                count++;
#if CONTEXT_DEBUG
                contextIterPrintCharIndices(&ci);
                printf( "  %4.2f  %6.2f\n", val, sum ); 
#endif
        } while( contextIterIncrement(&ci) );
        sum /= count;
    }
    return sum;
}


/*******************************************************************************
 * Weight_from_context()
 * Inputs:
 *     base_code[0]    = current base
 *     base_code[-1]   = previous base
 *          .
 *        (etc.)
 *          .
 *     base_code[-n+1] = "oldest" base
 * where "base" may be any one of: A,C,G,T, R,Y,K,M,S,W, B,D,H,V, N
 ******************************************************************************/
double 
weight_from_context( const char base_code[], ContextTable *ctable )
{
    int dim, d;
    char context[32];   /* 4^^32 is a big number */

    dim = ctable->dimension;
    for( d=0; d<dim; d++ ) {
        context[d] = base_code[-d];
    }
    return weight_from_reverse_context( context, ctable );
}

/* ****************************************************************************
 * This function reads in and parses a context table.  Its synopsis is:
 *
 * ctable = read_context_table(path)
 *
 * where
 *      path    is the name of the file containing the context table
 *
 *      ctable  is a CntextTable allocated and populated by this function
 *              on success, NULL otherwise.
 *
 * If this function returns a table, then the caller should free up the
 * resources (when through with it) by calling destroy_context_table().
 * 
 * Each line in the ctable should consist of 16 float numbers
 ******************************************************************************/
ContextTable *
read_context_table(char *path)
{
    FILE         *fp;
    int           i=0, j, size;
    char         *s, linebuf[MAXLINE];  
    ContextTable *ctable = NULL;

    if ((fp = fopen(path, "r")) == NULL) {
        return(ctable);
    }

    ctable = CALLOC(ContextTable, 1);
    size = CHUNK;
    ctable->weights = (double *)malloc(sizeof(double) * size);
    while (fgets(linebuf, MAXLINE, fp) != NULL) {
        /* Ignore all white space lines and comments */
        if (strspn(linebuf, " \t\r\n") == strlen(linebuf)) {
            continue;
        }

        if (((linebuf[0] == '/') && (linebuf[1] == '*'))
          || (linebuf[0] == ';') || (linebuf[0] == '#'))
        {
            continue;
        }

        s = linebuf;
        j = 0;
        while (Btk_atod(&s, &ctable->weights[i]) == 1) {
            i++;
            j++;
        }
        if (j == 0) {
            goto error_return;
        }

        if (i >= size/2) {
            size += CHUNK;
            ctable->weights = REALLOC(ctable->weights, double, size);
        }
    }
    if (ferror(fp)) {
        goto error_return;
    }

    /* Close file and give back unused space */
    (void)fclose(fp);
    ctable->weights = REALLOC(ctable->weights, double, i);

    /* Determine the dimension of the ctable */
    if (i%4 != 0) {
        fprintf(stderr, "Incorrect size of the context table\n");
        goto error_return; 
    }
    j = 1;
    size = 4;
    while (size < i) {
        j++;
        size *= 4;
        if (size > i) {
            fprintf(stderr, "Incorrect size of the context table\n");
            goto error_return;
        }         
    }
    ctable->dimension = j;

    return(ctable);

error_return:
    (void)fclose(fp);
        FREE(ctable->weights);
        FREE(ctable);
    return(NULL);
}

/*****************************************************************************
 * Function: destroy_context_table
 * Purpose:  reclaim the storage taken up by a context table previously
 *           returned by read_context_table()
 * Synopsis:
 *
 *     destroy_context_table(ctable)
 *
 * where
 *     ctable  is a ContextTable previously returned by
 *             read_context_table().
 */
void
destroy_context_table(ContextTable *ctable)
{
    if (ctable == NULL) {
        return;
    }
    FREE(ctable->weights);
}
