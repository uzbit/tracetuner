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
 * 1.9 2003/11/05 22:58:56
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "Btk_lookup_table.h"
#include "Btk_atod.h"
#include "Btk_qv.h"
#include "Btk_qv_data.h"
#include "util.h"
#include "train.h"
#include "Btk_compute_qv.h"
#include "params.h"
#include "lut.h"
#include "check_data.h"

typedef enum{
    CROSS_VALIDATE=1,    /* do the standard check table. */
    ERROR_STATS,         /* raw error statistics. */
    CUTOFF_STATS,        /* error statistics with QV cutoff. */
} REPORT_TYPE;

static REPORT_TYPE ReportType;
static int CutOff;
static int BinSize;

static int
errorStatswithQVs( char *tablefile , int minqv )
{
    BtkLookupTable *table;
    double pres, phr3, phr7, psr7, iheight;
    int nbases, match, qv, i;
    int wrong[MAXNUMBINS+1], total[MAXNUMBINS+1], total_cutoff[MAXNUMBINS+1];
    int maxbin;
    int spos, cpos;
    char schar, cchar;

    double mean_phr3[MAXNUMBINS+1], mean_phr7[MAXNUMBINS+1], mean_psr7[MAXNUMBINS+1], 
      mean_pres[MAXNUMBINS+1];
    int Num_phr3[MAXNUMBINS+1], Num_phr7[MAXNUMBINS+1], Num_psr7[MAXNUMBINS+1], 
      Num_pres[MAXNUMBINS+1];






    if ((table = Btk_read_lookup_table(tablefile)) == NULL) {
	(void)fprintf(stderr, "%s: couldn't read lookup table\n", tablefile);
	exit(1);
    }
    fprintf(stderr, "%s: %d entries\n", tablefile, table->num_lut_entries);

    (void)memset(wrong, 0, sizeof(wrong));
    (void)memset(total, 0, sizeof(total));
    (void)memset(total_cutoff, 0, sizeof(total_cutoff));

    maxbin=0;
    nbases = 0;

    for (i=0; i<MAXNUMBINS+1; i++){
      mean_phr3[i]=0.;
      mean_phr7[i]=0.;
      mean_pres[i]=0.;
      mean_psr7[i]=0.;
      Num_phr3[i]=0;
      Num_phr7[i]=0;
      Num_pres[i]=0;
      Num_psr7[i]=0;
    }

    while (getbase(&spos, &schar, &cpos, &cchar, &phr3,
                   &phr7, &psr7, &pres, &match, &iheight)) {
	if (++nbases % 10000 == 0) {
	    (void)fprintf(stderr, "\r%d bases read so far...", nbases);
	}
	total[(spos-1)/BinSize]++;
        if ( (spos-1)/BinSize + 1 > maxbin )
            maxbin = (spos-1)/BinSize + 1;

	  
	qv=get_quality_value(phr3, phr7, psr7, pres, table, USE_SORTED_LUT);
	/*	printf("qv=%d\n", qv);*/
	/*printf("Strange!\n");*/

	  if ( qv < minqv ) { continue; }



	mean_phr3[(spos-1)/BinSize]=mean_phr3[(spos-1)/BinSize]+phr3;
	Num_phr3[(spos-1)/BinSize]++;

	mean_phr7[(spos-1)/BinSize]=mean_phr7[(spos-1)/BinSize]+phr7;
	Num_phr7[(spos-1)/BinSize]++;

	mean_psr7[(spos-1)/BinSize]=mean_psr7[(spos-1)/BinSize]+psr7;
	Num_psr7[(spos-1)/BinSize]++;
	
	mean_pres[(spos-1)/BinSize]=mean_pres[(spos-1)/BinSize]+pres;
	Num_pres[(spos-1)/BinSize]++;


    }
    (void)fprintf(stderr, "\r%d bases total.        \n", nbases);


    (void)printf("coord   total   phr3     phr7     psr7      pres\n");
    for (i = 0; i < maxbin; i++) {

	mean_phr3[i]=mean_phr3[i]/Num_phr3[i];
	mean_phr7[i]=mean_phr7[i]/Num_phr7[i];
	mean_psr7[i]=mean_psr7[i]/Num_psr7[i];
	mean_pres[i]=mean_pres[i]/Num_pres[i];

	(void)printf("%4d  %9d  %6.3f   %6.3f   %6.3f   %6.3f\n",
		i*BinSize+BinSize/2, 
    		total[i], 
             	mean_phr3[i], mean_phr7[i], mean_psr7[i], mean_pres[i]);
    }


    Btk_destroy_lookup_table(table);
    return(0);
}

void
show_usage(char *argv[])
{
    (void)fprintf(stderr, "\nVersion: %s\n", TT_VERSION);
    (void)fprintf(stderr, "usage: %s \n\
  [ -e ] \n\
  [ -c <cutoff> ] \n\
  [ -b <binsize> ] \n\
  <lookuptable> < <basefile>\n", 
    argv[0]);
}

int
main(int argc, char *argv[])
{
    int i, cutoffspecified;

    /* Set defaults. */
    ReportType = CROSS_VALIDATE;
    CutOff = -1;
    cutoffspecified=0;
    BinSize = 1; /* Bases 0-99 in one bin, 100-199 in a second, etc. */

    if (argc == 1) {
	show_usage(argv);
	exit(2);
    }

    while( (i = getopt(argc, argv, "b:c:e") ) != EOF ) {

        switch(i) {
        case 'b':
            BinSize = atoi(optarg);
            if ( BinSize <= 0 ) {
                fprintf(stderr, "Error: illegal bin size %d\n", BinSize);
                exit(-1);
            }
            break;
        case 'c':
            ReportType = CUTOFF_STATS;
            CutOff = atoi(optarg);
            cutoffspecified++;
            break;
	case 'e':
	    ReportType = ERROR_STATS;
	    break;
	default:
	    show_usage(argv);
	    exit(2);
	}
	
  }

    fprintf(stderr, "BinSize = %d\n", BinSize);
    setbuf(stderr, NULL);	/* unbuffered */

     errorStatswithQVs(argv[optind], CutOff);
    /*    crossValidate(argv[optind]);*/

/* 
    switch(ReportType) {
    case CROSS_VALIDATE:

      	if ( optind >= argc ) {
	    show_usage(argv);
	    exit(2);
	} else {
      
	crossValidate(argv[optind]);
           
	   	}
	break;
    case ERROR_STATS:
        errorStats();
	break;
    case CUTOFF_STATS:
	if ( (optind >= argc) || (!cutoffspecified) ) {
	    show_usage(argv);
	    exit(2);
	} else {
	printf("Mark 3 \n");
            errorStatswithQVs(argv[optind], CutOff);
	}
	break;
    }
*/
    return(0);
}

