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

/**
 ** chkorder.c
 **
 ** Test program to make sure base calls and quality values don't depend
 ** the order color traces are presented to the code.  This code was taken,
 ** in large part, from example.c.
 **/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdint.h>

#include "Btk_qv.h"
#include "Btk_qv_data.h"
#include "Btk_lookup_table.h"
#include "Btk_default_table.h"
#include "util.h"
#include "train.h"
#include "Btk_compute_qv.h"
#include "Btk_match_data.h"
#include "Btk_qv_io.h"


/*
 * Call some bases and qvs, and (optionally) check against a "standard".
 *
 * r = callbases(filename, nbases, bases, locations, nvals, vals, colororder,
 *			table, options,
 *			firstnbases, firstbases, firstlocations, firstqv,
 *			callednbases, calledbases, calledlocations, qv,
 *			msg)
 *
 * filename		is the root filename of the sample file
 * nbases		is the length of the bases and locations arrays
 * bases		is the address of the first element of a char array
 *			containing the original base calls
 * locations		is the address of the first element of an int array
 *			containing the positions in the traces of the bases
 *			called in bases
 * nvals		is the length of the arrays contained in vals
 * vals			is the address of an array of int arrays containing
 *			the trace values (4 in all).  The traces MUST be in
 *			ACGT order
 * colororder		is a string with 'A', 'C', 'G', and 'T' in the order
 *			that the trace arrays (from vals) should be put into
 *			to call bases and qvs
 * table		is the address of the BtkLookupTable to use
 * options		is an Options structure
 * firstnbases		is the length of the firstbases, firstlocations, and
 *			firstqv arrays, if any (0 => don't compare)
 * firstbases		is the address of an array of quality values to
 *			to compare to, if any (NULL => don't compare)
 * firstlocations	is the address of an array of quality values to
 *			to compare to, if any (NULL => don't compare)
 * firstqv		is the address of an array of quality values to
 *			to compare to, if any (NULL => don't compare)
 * callednbases		is the address of an int var where this function will
 *			put the length of the calledbases, calledlocations,
 *			and qv arrays
 * calledbases		is the address of a char* var where this function will
 *			put the final base calls
 * calledlocations	is the address of an int* var where this function will
 *			put the locations where calledbases were found
 * qv			is the address of an int* var where this function will
 *			put the quality values corresponding to calledbases
 * msg			is the address of a BtkMessage structure into which
 *			a message will be put on error
 *
 * r			is zero on success, non-zero on failure
 */
static int
callbases(
    char *filename,
    int nbases,
    char *bases,
    int *locations,
    int nvals,
    int *vals[4],		/* must be in ACGT order */
    char *colororder,
    BtkLookupTable *table,
    Options options,
    int firstnbases,
    char *firstbases,
    int *firstlocations,
    uint8_t *firstqv,
    int *callednbases,
    char **calledbases,
    int **calledlocations,
    uint8_t **qv,
    BtkMessage *msg)
{
    int *tmpvals[4], r, basesdiff, i;
    Results results;

    /* set up the trace arrays */
    switch(colororder[0]) {
    case 'A':
	 tmpvals[0] = vals[0];
	 break;
    case 'C':
	 tmpvals[0] = vals[1];
	 break;
    case 'G':
	 tmpvals[0] = vals[2];
	 break;
    case 'T':
	 tmpvals[0] = vals[3];
	 break;
    default:
	 (void)sprintf(msg->text, "unknown color '%c'", colororder[0]);
	 return -1;
    }
    switch(colororder[1]) {
    case 'A':
	 tmpvals[1] = vals[0];
	 break;
    case 'C':
	 tmpvals[1] = vals[1];
	 break;
    case 'G':
	 tmpvals[1] = vals[2];
	 break;
    case 'T':
	 tmpvals[1] = vals[3];
	 break;
    default:
	 (void)sprintf(msg->text, "unknown color '%c'", colororder[1]);
	 return -1;
    }
    switch(colororder[2]) {
    case 'A':
	 tmpvals[2] = vals[0];
	 break;
    case 'C':
	 tmpvals[2] = vals[1];
	 break;
    case 'G':
	 tmpvals[2] = vals[2];
	 break;
    case 'T':
	 tmpvals[2] = vals[3];
	 break;
    default:
	 (void)sprintf(msg->text, "unknown color '%c'", colororder[2]);
	 return -1;
    }
    switch(colororder[3]) {
    case 'A':
	 tmpvals[3] = vals[0];
	 break;
    case 'C':
	 tmpvals[3] = vals[1];
	 break;
    case 'G':
	 tmpvals[3] = vals[2];
	 break;
    case 'T':
	 tmpvals[3] = vals[3];
	 break;
    default:
	 (void)sprintf(msg->text, "unknown color '%c'", colororder[3]);
	 return -1;
    }

    *callednbases = nbases;
    *calledbases  = malloc(nbases + 1);
    (void)memcpy(*calledbases, bases, nbases + 1);
    *calledlocations = (int *)malloc(nbases * sizeof(locations[0]));
    (void)memcpy(*calledlocations, locations, nbases * sizeof(locations[0]));
    *qv = (uint8_t *)malloc(nbases * sizeof((*qv)[0]));

    (void)fprintf(stderr, "%s.%s: ", filename, colororder);
    (void)fflush(stderr);
    if ((r = Btk_compute_qv(callednbases, calledbases, calledlocations, &nvals,
			    tmpvals, colororder, table, qv, options, msg,
             &results))
	!= 0)
    {
	(void)fprintf(stderr, "%s\n", msg->text);
	return r;
    }

    /* Check against the "standard" */
    if (firstnbases > 0) {
	if (*callednbases != firstnbases) {
	    (void)fprintf(stderr, "%d bases vs. %d\n", *callednbases,
			    firstnbases);
	}
	else {
	    basesdiff = 0;
	    for (i = 0; (i < *callednbases) && (i < firstnbases); i++) {
		if ((*calledbases)[i] != firstbases[i]) {
		    basesdiff++;
		}
	    }
	    if (basesdiff > 0) {
		(void)fprintf(stderr, "%d bases different\n", basesdiff);
	    }
	    else {
		(void)fprintf(stderr, "ok\n");
	    }
	}
    }
    else {
	(void)fprintf(stderr, "%d bases\n", *callednbases);
    }

    return 0;
}


int
main(int argc, char *argv[])
{
    int i, n;
    char *lookup_table, *smp, *smptail;
    BtkLookupTable *table;
    char *chemistry = "", status_code[100];
    int nvals, *vals[4], filetype = -1;
    uint8_t *quality_values = NULL;
    BtkMessage  msg;
    Options options;

    /* original values read from the sample file */
    char *origbases;
    int orignbases, *origlocations;

    /* values computed the first time, for comparison purposes */
    char *firstbases;
    int firstnbases, *firstlocations;

    /* values computed each time */
    char *bases;
    int nbases, *locations; 
    uint8_t *qv, *firstqv;

    /* all permutations of colors */
    char *orders[] = {
	"ACGT",
	"ACTG",
	"AGCT",
	"AGTC",
	"ATCG",
	"ATGC",
	"CAGT",
	"CATG",
	"CGAT",
	"CGTA",
	"CTAG",
	"CTGA",
	"GACT",
	"GATC",
	"GCAT",
	"GCTA",
	"GTAC",
	"GTCA",
	"TACG",
	"TAGC",
	"TCAG",
	"TCGA",
	"TGAC",
	"TGCA"
    };
    int numorders = sizeof(orders) / sizeof(orders[0]);


    if (argc < 2) {
        (void)fprintf(stderr, "usage: %s <samplefiles...>\n", argv[0]);
        exit(1);
    }

    /*
     *  If a non-standard lookup table is specified as an environment
     *  variable use it, otherwise NULL gets the default.
     */
    lookup_table = getenv("LOOKUP_TABLE");

    if (lookup_table != NULL) {
        if (strcmp(lookup_table, "3730") == 0) {
            table = Btk_get_3730pop7_table();
        }
        else if (strcmp(lookup_table, "3700pop5") == 0) {
            table = Btk_get_3700pop5_table();
        }
        else if (strcmp(lookup_table, "3700pop6") == 0) {
            table = Btk_get_3700pop6_table();
        }
        else {
            if ((table = Btk_read_lookup_table(lookup_table)) == NULL) {
                fprintf(stderr, "Couldn't read lookup table '%s'.\n", 
                        lookup_table);
                exit(1);
            }
        }
    }
    else {
        table = NULL;
    }

    for (n = 1; n < argc; n++) {
        smp = argv[n];
        if ((smptail = strrchr(smp, '/')) != NULL) {
            smptail++;
	}
        else {
            smptail = smp;
	}

        /* Setting default options 
         */
        options.nocall         = 0;
	options.respace        = 0;
        options.recalln        = 0;
        options.edited_bases   = 0;
        options.gauss          = 1;
        options.tal_dir[0]     = '\0';
        options.tab_dir[0]     = '\0';
        options.het            = 0;
        options.time           = 0;
        options.min_ratio      = (float)0.15;
        options.Verbose        = 3;
        options.inp_phd        = 0;
        options.inp_phd_dir[0] = '\0';
        options.file_name[0]   = '\0';
        options.lut_type       = ABI3730pop7;

        if (Btk_read_sample_file(smp, &orignbases, &origbases, 0,
		    &origlocations, &quality_values, &nvals, &vals[0], &vals[1],
		    &vals[2], &vals[3], NULL, &chemistry, status_code, &filetype, 
                    options, &msg)
	    != 0)
        {
            (void)fprintf(stderr, "%s: couldn't read sample file\n", smptail);
            continue;
        }

        (void)strcpy(options.file_name, smptail);

        if (table == NULL) {
            if (strstr(chemistry, "3730")) {
                table = Btk_get_3730pop7_table();
            }
            else if (strstr(chemistry, "POP5")) {
                table = Btk_get_3700pop5_table();
            }
            else if (strstr(chemistry, "POP6") && 
                     strstr(chemistry, "3700")) {
                table = Btk_get_3700pop6_table();
            }
            else if (strstr(chemistry, "3100")) {
                table = Btk_get_3100pop6_table();
            }
            else {
		(void)fprintf(stderr,
				"Can't select the lookup table automatically.\n"
				"Using built-in ABI 3730 Pop-7 table.\n");
		table = Btk_get_3730pop7_table();
            }
        }

	/*
	 * Do it the first time, but remember the answers.
	 */
	if (callbases(smptail, orignbases, origbases, origlocations,
			nvals, vals, orders[0], table, options,
			0, NULL, NULL, NULL,
			&firstnbases, &firstbases, &firstlocations, &firstqv,
			&msg)
	    != 0)
	{
            (void)fprintf(stderr, "%s\n", msg.text);
            goto cleanup_a_file;
        }

	for (i = 1; i < numorders; i++) {
	    if (callbases(smptail, orignbases, origbases, origlocations,
			    nvals, vals, orders[i], table, options,
			    firstnbases, firstbases, firstlocations, firstqv,
			    &nbases, &bases, &locations, &qv, &msg)
		!= 0)
	    {
		(void)fprintf(stderr, "%s\n", msg.text);
		goto cleanup_a_file;
	    }
	}

    }

cleanup_a_file:
    FREE(qv);
    FREE(quality_values);
    FREE(quality_values);
    FREE(bases);
    FREE(locations);
    for (i=0; i<4; i++) {
	FREE(vals[i]);
    }
    Btk_destroy_lookup_table(table);

    return(0);
}
