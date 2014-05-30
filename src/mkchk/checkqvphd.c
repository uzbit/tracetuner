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
 * $Id: checkqvphd.c,v 1.7 2009/01/12 15:28:40 gdenisov Exp $                     
 *
 * Produce an output similar to the check program using on
 * Tracetuner train files except this version should run on
 * phred train files.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "Btk_atod.h"
#include "Btk_lookup_table.h"
#include "Btk_qv_data.h"
#include "params.h"
#include "lut.h"
#include "check_data.h"
#include "Btk_qv.h"  

static int BinSize;
static int pyrosequencing;

/*******************************************************************************
 * Function: crossValidate
 *******************************************************************************
 */
static int
crossValidate(void)
{
    unsigned long    nbases, ngood_bases, errbases, nmsbases; 
    int    qv, i;
    int    wrong[MAXQVALUE+1], total[MAXQVALUE+1];
    double weightederrori, observed;
    int    cum_total=0;
    int    max_qv = 0;

    BASE          *bases = NULL;
    unsigned long  initial_base_room = BASE_COUNT_SCALE;
    unsigned long  base_count;

    (void)memset(wrong, 0, sizeof(wrong));
    (void)memset(total, 0, sizeof(total));

    nbases = ngood_bases = errbases = nmsbases = 0;

    fprintf(stderr, "Reading bases from stdout ... \n");
    bases = get_bases(initial_base_room, &base_count, 4, 1);
    nbases = base_count;
    (void)fprintf(stderr, "\r%lu bases total.        \n", nbases);

    fprintf(stderr, "nbases= %lu ngoodbases=%lu\n", nbases, ngood_bases);

    for (i=0; i<nbases; i++)
    {
        if (bases[i].schar != '-')
            ngood_bases++;

        qv = bases[i].qv;
        if (max_qv < qv) { max_qv = qv; }

        total[qv/BinSize]++;
        if (bases[i].is_match == (char)0)
        {
            if (pyrosequencing && i > 0 && i < nbases-1)
            {
                if (bases[i].schar == bases[i-1].schar && bases[i].schar != '-')
                {
                    wrong[bases[i-1].qv/BinSize]++;
                    errbases++;
                }
                else if (bases[i].schar == bases[i+1].schar && bases[i].schar != '-')
                {
                    wrong[bases[i+1].qv/BinSize]++;
                    errbases++;
                }
                else
                {
                    wrong[bases[i].qv/BinSize]++;
                    errbases++;
                }
            }
            else
            {
                wrong[bases[i].qv/BinSize]++;
                errbases++;
            }
        }
    }
    (void)fprintf(stderr, "\r%lu bases total.        \n", nbases);
    (void)printf("predict  observe  #correct  #incorrect  #total");
    (void)printf(" frac_qv_or_better\n");
    weightederrori = 0;

    for (i = 0; i <= max_qv; i++) 
    {
        cum_total += total[i];
	    if (total[i] > 0) 
        {
	        observed = calc_qv(wrong[i], total[i]);
	        (void)printf("%7d  %5d  %7d  %8d   %7d      %6.3f\n",
			i, (int)rint(observed), total[i] - wrong[i], wrong[i],
			total[i],
                        (float)(nbases-cum_total)/(float)nbases);
	        weightederrori += total[i] * fabs(i - rint(observed));
        }
        else {
            (void)printf( "%7d                                   0      %6.3f\n",
            i, (float)(nbases-cum_total)/(float)nbases);
                                }

    }

    (void)printf("\n%lu bases total\n", nbases); 
    (void)printf("\nNumber of matched or substituted bases=%lu\n",
        nmsbases);
    (void)printf("ngood_bases=%lu", ngood_bases);
    (void)printf("Fraction of errors = %f\n",
        (double)errbases/(double)nbases);
    if( ngood_bases ) {
        weightederrori /= ngood_bases;
        (void)printf("Accuracy:  %6.4f\n", weightederrori);
    } else {
        printf( "Accuracy: undefined\n" );
    }
    FREE(bases);
    return(0);
}

/*******************************************************************************
 * Function: show_usage    
 *******************************************************************************
 */
void
show_usage(char *argv[])
{
    (void)fprintf(stderr, "\nVersion: %s\n", TT_VERSION);
    (void)fprintf(stderr, "usage: %s [ -h ] \
   <    <alignment_file>\n", argv[0]);
}

/*******************************************************************************
 * Function: main       
 *******************************************************************************
 */
int
main(int argc, char *argv[])
{
    int i;

    /* Set defaults. */
    BinSize = 1;
    pyrosequencing = 0;

    while( (i = getopt(argc, argv, "b:ph") ) != EOF ) {

        switch(i) {
	        case 'b':
	        BinSize = atoi(optarg);
	        if ( BinSize <= 0 ) {
		        fprintf(stderr, "Error: illegal bin size %d\n", BinSize);
		        exit(-1);
	        }
	        break;
            case 'p':
                pyrosequencing++;
                break;
            case 'h':
                show_usage(argv);
                exit(2);
        	default:
            break;
	    }
    }

    setbuf(stderr, NULL);	/* unbuffered */

    crossValidate();

    return(0);
}

