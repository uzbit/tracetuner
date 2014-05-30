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
 * $Id: checkbc.c,v 1.11 2009/01/05 21:29:31 gdenisov Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "Btk_atod.h"
#include "Btk_qv.h"
#include "Btk_qv_data.h"
#include "util.h"
#include "train.h"
#include "params.h"
#include "lut.h"
#include "check_data.h"
 
static int BinSize = 1;

/*******************************************************************************
 * Function: errorStatswithQVs
 *******************************************************************************
 */
static int
errorStatswithQVs( )
{
    char   schar, cchar, prev_schar = ' ';
    int    nbases, match, nmsbases,
           Nqvp[MAXNUMBINS+1], i; 
    int    wrong[MAXNUMBINS+1], total[MAXNUMBINS+1]; 
    int    maxbin;
    int    spos, cpos;
    int    Ins=0, Del=0, Sub=0; 
    int    ninsertion[MAXNUMBINS+1], ndeletion[MAXNUMBINS+1], 
           nsubstitution[MAXNUMBINS+1];
    double phr3, phr7, psr7, pres, iheight, iheight1, iheight2;

    (void)memset(wrong, 0, sizeof(wrong));
    (void)memset(total, 0, sizeof(total));

    maxbin = 0;
    nbases = 0;
    nmsbases = 0;
    while (getbase(&cpos, &cchar,  &match, &spos, &schar, &phr3,
                   &phr7, &psr7, &pres, &iheight, &iheight1, &iheight2)) {
	if (++nbases % 10000 == 0) {
	    (void)fprintf(stderr, "\r%d bases read so far...", nbases);
	}

        if (cchar == 'X')
            continue;
        if ((cchar == 'N') && (schar != '-'))
            continue;

	total[(spos-1)/BinSize]++;
        if (schar != cchar ) {
            wrong[(spos-1)/BinSize]++;
        }

        if ( (spos-1)/BinSize + 1 > maxbin )
            maxbin = (spos-1)/BinSize + 1;

        if ( schar == '-' ) { 
           Del++;
	   ndeletion[(spos-1)/BinSize]++;
           prev_schar = schar;
	   continue; 
        }
        else
        {
	    if ( cchar == '-' ) 
            {
                Ins++;
	        ninsertion[(spos-1)/BinSize]++;
    	    } 
            else if (schar != cchar) 
            {
                Sub++;
    	        nsubstitution[(spos-1)/BinSize]++;
	    }
            if (prev_schar == '-')
            {
                Nqvp[(spos-2)/BinSize]++;
            }
	 
            if (cchar != '-')
                nmsbases++;
        }
 
        prev_schar = schar;
    }
    (void)fprintf(stderr, "\r%d bases total.        \n", nbases);

    (void)printf("coord  Err_tot  Err_del  Err_ins  Err_sub "); 
    (void)printf(" #Bases  #Del   #Ins  #Sub   #Err  \n");
    for (i = 0; i < maxbin; i++) {
	(void)printf("%4d    %6.3f  %6.3f   %6.3f   %6.3f %5d  %5d  %5d  %5d  %5d  \n",
		i*BinSize+BinSize/2, 
        100.*(double)(ndeletion[i]+ninsertion[i]+nsubstitution[i])/(total[i]>0?(double)total[i]:1.),
        100.*(double)ndeletion[i]/(total[i]>0?(double)total[i]:1.),
        100.*(double)ninsertion[i]/(total[i]>0?(double)total[i]:1.),
        100.*(double)nsubstitution[i]/(total[i]>0?(double)total[i]:1.),
        total[i], ndeletion[i], ninsertion[i], nsubstitution[i],
        ndeletion[i]+ninsertion[i]+nsubstitution[i]);
    }

    (void)printf("\nNumber of bases=%d\n", nbases);
    (void)printf("\nNumber of matched or substituted bases=%d\n",
        nmsbases);
    (void)printf("Number of errors: Ins=%d Del=%d Sub=%d Tot=%d\n",
        Ins, Del, Sub, Ins+Del+Sub);

    return(0);
}

void
show_usage(char *argv[])
{
    (void)fprintf(stderr, "\nVersion: %s\n",TT_VERSION);
    (void)fprintf(stderr, 
        "usage: %s    [ -h ]     <    <alignment_file>\n", 
    argv[0]);
}

int
main(int argc, char *argv[])
{
    int i;

    while( (i = getopt(argc, argv, "b:h") ) != EOF ) {

        switch(i) {
        case 'h':
            show_usage(argv);
            exit(2);
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

    fprintf(stderr, "Reading bases from stdout ...\n");
    setbuf(stderr, NULL);	/* unbuffered */

     errorStatswithQVs(argv[optind]);
    /*    crossValidate(argv[optind]);*/
    return(0);
}

