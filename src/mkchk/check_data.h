/* 
 * Copyright (c) 1999-2003 Paracel, Inc.  All rights reserved.
 *
 * $Id: check_data.h,v 1.7 2009/01/12 15:28:40 gdenisov Exp $                 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAXLINE         (1000)
#define CHUNK           (100)
#define MAXQVALUE       (100)
#define MAXNUMBINS      (5000)

extern BASE *get_bases(unsigned long , unsigned long *, int, char);
extern int   getbase(int *, char* , int *, int* , char* , 
    double *, double *, double *, double *, 
    double *, double *, double *);
extern int   getbase_phd(int *, char* , int *, int* , char*,
                         int *);
extern double calc_qv(int , int );
extern double calc_qvP(double , int );
extern int is_mixed_base(char );
