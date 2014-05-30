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
 * $Id: checkqv.c,v 1.11 2009/01/12 15:28:40 gdenisov Exp $                   
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>

#include "Btk_lookup_table.h"
#include "Btk_atod.h"
#include "Btk_qv.h"
#include "util.h"
#include "train.h"
#include "Btk_qv_data.h"
#include "Btk_compute_qv.h"
#include "params.h"
#include "lut.h"
#include "check_data.h"

static int BinSize;

/*******************************************************************************
 * Function: crossValidate     
 *******************************************************************************
 */
static int
crossValidate( char *tablefile )
{
    BtkLookupTable *table;
    unsigned long  i, nbases, ngood_bases, errbases; 
    int    qv;
    int    wrong[MAXQVALUE+1], total[MAXQVALUE+1];
    double weightederrori, observed;
    int    cum_total=0;
    int    max_qv = 0;

    BASE          *base = NULL;
    unsigned long  initial_base_room = BASE_COUNT_SCALE;
    unsigned long  base_count;

    fprintf(stderr, "Cross-validating ...\n");

    if ((table = Btk_read_lookup_table(tablefile)) == NULL) {
	(void)fprintf(stderr, "%s: couldn't read lookup table\n", tablefile);
	exit(1);
    }
    fprintf(stderr, "%s: %d entries\n", tablefile, table->num_lut_entries);

    (void)memset(wrong, 0, sizeof(wrong));
    (void)memset(total, 0, sizeof(total));

    nbases = ngood_bases = errbases = 0;

    fprintf(stderr, "Reading bases from stdout ... \n");
    base = get_bases(initial_base_room, &base_count, 4, 0);
    nbases = base_count;
    (void)fprintf(stderr, "\r%lu bases total.        \n", nbases);

    for (i=0; i<nbases; i++)
    {
        if (base[i].schar != '-')            
            ngood_bases++;

	    qv = get_quality_value(base[i].parameter[0], 
                               base[i].parameter[1],
                               base[i].parameter[2],
                               base[i].parameter[3], table);
        total[qv/BinSize]++;
        if (max_qv < qv) { max_qv = qv; }

        if (base[i].is_match == (char)0)
        {
            wrong[qv/BinSize]++;
            errbases++;
        }
        if (i > 0 && i%100000 == 0)
        {
            fprintf(stderr, "\r%lu bases processed", i);
        }
    }
    (void)fprintf(stderr, "\r%lu bases total.        \n", nbases);

    (void)printf("pred_QV  obs_QV    #correct  #incorrect  #total");
    (void)printf(" frac_given_qv_or_better\n");
    weightederrori = 0;
    for (i = 0; i <= max_qv; i++) {
        if ( i*BinSize > MAXQVALUE ) { break; }
  	    if (total[i] > 0) {
            cum_total += total[i];
	        observed = calc_qv(wrong[i], total[i]);
	        (void)printf("%5lu  %5d      %7d  %7d    %7d     %6.3f\n",
                         i*BinSize, (int)rint(observed), 
                         total[i] - wrong[i], wrong[i],
                         total[i], 
                         (float)(nbases-cum_total)/(float)nbases);
	        weightederrori += total[i] * fabs(i*BinSize - rint(observed));
	    }
        else {
            (void)printf( "%5lu                                       0     %6.3f\n",
                i*BinSize, (float)(nbases-cum_total)/(float)nbases);
        }

    }
    (void)printf("\n%lu bases total\n", nbases);
#if 0
        printf("ngood_bases=%d", ngood_bases);
#endif
    (void)printf("Fraction of errors = %f\n",
        (double)errbases/(double)nbases);
    if( ngood_bases ) {
        weightederrori /= ngood_bases;
        (void)printf("Accuracy:  %6.4f\n", weightederrori);
    } else {
        printf( "Accuracy: undefined\n" );
    }
    Btk_destroy_lookup_table(table);
    FREE(base);
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
    (void)fprintf(stderr, "usage: %s    [ -h ] \
  <lookuptable>   <    <alignment_file>\n", 
    argv[0]);
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

    if (argc == 1) {
	show_usage(argv);
	exit(2);
    }

    while( (i = getopt(argc, argv, "b:") ) != EOF ) {
        switch(i) {
	        case 'b':
	            BinSize = atoi(optarg);
	            if ( BinSize <= 0 ) {
		        fprintf(stderr, "Error: illegal bin size %d\n", BinSize);
		        exit(-1);
	            }
	            break;
	        default:
	            show_usage(argv);
	            exit(2);
	    }
    }

    setbuf(stderr, NULL);	/* unbuffered */

	if ( optind >= argc ) {
	    show_usage(argv);
	    exit(2);
	} else {
            crossValidate(argv[optind]);
	}
	    
    return(0);
}
