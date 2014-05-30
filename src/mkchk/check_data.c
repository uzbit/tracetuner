/* 
 * Copyright (c) 1999-2003 Paracel, Inc.  All rights reserved.
 *
 * $Id: check_data.c,v 1.13 2009/01/12 15:28:40 gdenisov Exp $                   
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "params.h"
#include "lut.h"
#include "Btk_lookup_table.h"
#include "Btk_atod.h"
#include "Btk_qv.h"
#include "util.h"

#define MAXLINE		(1000)
#define CHUNK		(100)
#define MAXQVALUE	(60)
#define MAXNUMBINS      (5000)

/***************************************************************************
 * get_bases
 * purpose: read bases from stdin to populate the pre-malloced array
 * <base> of length <base_count>. stdin should have the following format:
 * returns: the number of vali
 * id is_match parameter1 parameter2 parameter3 parameter4 is_valid
 *
 * called by: main
 * calls: none
 *
 * ASSUMES: 4 parameters
 *          input lines consist of
 *             sample_pos sample_char consensus_pos cons_char training_params
 *
 ***************************************************************************/
BASE *
get_bases(unsigned long initial_base_room, unsigned long *base_count,
    int parameter_count, char phred)
{
    int           i, j;
    unsigned long base_room, bases_used = 0;
    BASE         *bases;
    int           linenum = 0, prev_good_base = 0;
    char          linebuf[1024], *s, new_frag_beg = 0;
    double        is_match, qv;
    int           spos, cpos;   /* Sample and consensus positions. */
    char          schar, cchar, prev_schar='X'; /* Sample and consensus characters. */
    BtkMessage    msg, *message = &msg;;

    base_room = initial_base_room;
    bases = (BASE *) malloc(base_room * sizeof(bases[0]));

    while (fgets(linebuf, sizeof(linebuf), stdin) != NULL)
    {
        linenum++;

/* Ignore all white space lines and comments */
        if (strspn(linebuf, " \t\r\n") == strlen(linebuf))
            continue;

        if (linebuf[0] == '#' && linebuf[2] == 'F')
            new_frag_beg = 1;

        if (linebuf[0] == '#' || linebuf[0] == ';' ||
           (linebuf[0] == '/' && linebuf[1] == '*'))
            continue;

/* Get rid of 1) sample position, 2) sample base,
 * 3) consensus position and 4) consensus base. */

        s = strtok(linebuf, " \t\n");
        if (s == NULL) continue;
        cpos = atoi(s);

        s = strtok(NULL, " \t\n");
        if (s == NULL) continue;
        cchar = s[0];

        s=strtok(NULL, " \t\n");
        if ( s==NULL ) { continue; }
        is_match=atoi(s);

        s = strtok(NULL, " \t\n");
        if (s == NULL) continue;
        spos = atoi(s);

        s = strtok(NULL, " \t\n");
        if (s == NULL) continue;
        schar = s[0];

        s += strlen(s) + 1;

        if ((cchar == 'X') || (cchar == 'N') || (schar == 'N')) {
            /* Consensus originally had '-', as a result of a polymorphism */
            continue;
        }
        if (schar == '-')  /* deleted base */
        {
            /* Consensus may have 'N' at this position!
             * We define fictitious values of * pres and phr3,
             * and use real phr7 and psr7 of previous base
             */
            if ((bases_used < 1) || new_frag_beg) {
                continue;
            }
            else {
                is_match = 0; 
                bases[bases_used].schar = schar;
                bases[bases_used].new_frag_beg = new_frag_beg;
                bases[bases_used].is_match = is_match;
                bases_used++;
                prev_schar = schar;
                if ((bases_used % 100000) == 0) {
                    fprintf(stderr, "\r%lu bases have been read", bases_used);
                }
                if (bases_used == base_room)
                {
                    base_room += BASE_COUNT_SCALE;
                    bases = (BASE*)realloc(bases, base_room * sizeof(bases[0]));
                }
                new_frag_beg = 0;
                continue;
            }
        }
        else  /* undeleted base */
        {
            if (!phred)
            {
                if (parameter_count == 4 &&
                   ((Btk_atod(&s, &bases[bases_used].parameter[0]) != 1)
                ||  (Btk_atod(&s, &bases[bases_used].parameter[1]) != 1)
                ||  (Btk_atod(&s, &bases[bases_used].parameter[2]) != 1)
                ||  (Btk_atod(&s, &bases[bases_used].parameter[3]) != 1)))
                {
                    fprintf(stderr,
                        "line %d:\n%s\nmissing/garbled parameter; skipping\n",
                        linenum, linebuf);
                    continue;
                }
                else if (parameter_count == 6 &&
                   ((Btk_atod(&s, &bases[bases_used].parameter[0]) != 1)
                ||  (Btk_atod(&s, &bases[bases_used].parameter[1]) != 1)
                ||  (Btk_atod(&s, &bases[bases_used].parameter[2]) != 1)
                ||  (Btk_atod(&s, &bases[bases_used].parameter[3]) != 1)
                ||  (Btk_atod(&s, &bases[bases_used].parameter[4]) != 1)
                ||  (Btk_atod(&s, &bases[bases_used].parameter[5]) != 1)))
                {
                    fprintf(stderr,
                        "line %d:\n%s\nmissing/garbled parameter; skipping\n",
                        linenum, linebuf);
                    continue;
                }

                if ((prev_schar == '-') && (bases_used > 0))
                {
                    i = bases_used-1;
                    while ((i>=0) && (bases[i].is_match == (char)0))
                    {
                        for (j=0; j<parameter_count; j++)
                            bases[i].parameter[j] =
                                QVMAX(bases[prev_good_base].parameter[j],
                                      bases[bases_used    ].parameter[j]);
                        i--;
                    }
                }
            }
            else
            {
                if (Btk_atod(&s, &qv) != 1) 
                {
                    fprintf(stderr,
                        "line %d:\n%s\nmissing/garbled parameter; skipping\n",
                        linenum, linebuf);
                    continue;    
                }
                if (phred)
                {
                    bases[bases_used].qv = (int)(qv + 0.0001);
                }

                if ((prev_schar == '-') && (bases_used > 0))
                {
                    int i = bases_used-1;
                    while ((i>=0) && (bases[i].is_match == (char)0))
                    {
                        bases[i].qv = QVMAX(bases[prev_good_base].qv,
                                            bases[bases_used    ].qv);
                        i--;
                    }
                }
            }
        }

        prev_schar     = schar;
        prev_good_base = bases_used;
        bases[bases_used].is_match = (char)is_match;
        bases[bases_used].schar = schar;
        bases[bases_used].new_frag_beg = new_frag_beg;
        bases_used++;
        new_frag_beg = 0;

        if ((bases_used % 100000) == 0) {
            fprintf(stderr, "\r%lu bases have been read", bases_used);
        }
        if (bases_used == base_room)
        {
            base_room += BASE_COUNT_SCALE;
            bases = (BASE*)realloc(bases, base_room * sizeof(bases[0]));
            MEM_ERROR(bases);
        }
    }
/* Give back any unused space */

    if (bases_used != base_room)
        bases = (BASE*)realloc(bases, bases_used * sizeof(bases[0]));

   *base_count = bases_used;

    return bases;
error:
    exit(-1);     
}

/*******************************************************************************
 * Function: getbase
 * This function reads a single valid base's parameters from the training file
 * (stdin).  Its synopsis is:
 *
 * success = getbase(spos, schar, cpos, cchar, phr3, phr7, psr7, pres, is_match)
 *
 * where
 *      spos            is where the sample position should go
 *      schar           is where the sample character (ACTG-) should go
 *      cpos            is where the consensus position should go
 *      cchar           is where the consensus character (ACTG-) should go
 *	phr3		are where the double parameters should go (by ref)
 *	phr7
 *	psr7
 *	pres
 *	is_match	is where the int parameter should go (by ref)
 *
 *	success		is 1 if a valid base is read, 0 otherwise
 *
 * ASSUMES: input lines consist of
 *             sample_pos sample_char consensus_pos cons_char training_params
 *
 * NOTE: if schar == '-', the training parameters are not set.
 *******************************************************************************
 */
int
getbase(int *cpos, char* cchar, int *is_match, int* spos, char* schar, 
        double *phr3, double *phr7, double *psr7, double *pres, 
        double *iheight, double *iheight1, double *iheight2)
{
    char linebuf[MAXLINE], *s;
    static int line = 0;

    for (;;) 
    {
	if (fgets(linebuf, MAXLINE, stdin) == NULL) {
	    return 0;
	}
        
        line++;

	/* Ignore all white space lines and comments */
	if (strspn(linebuf, " \t\r\n") == strlen(linebuf)) {
	    continue;
	}
	if ((linebuf[0] == '#')
	|| ((linebuf[0] == '/') && (linebuf[1] == '*'))
	|| (linebuf[0] == ';'))
	{
	    continue;
	}

	/* Get rid of 1) sample position, 2) sample base,
	 * 3) consensus position and 4) consensus base.
	 */
	s=strtok(linebuf, " \t\n");
	if ( s==NULL ) { continue; }
       *cpos=atoi(s);

	s=strtok(NULL, " \t\n");
	if ( s==NULL ) { continue; }
       *cchar=s[0];

        s=strtok(NULL, " \t\n");
        if ( s==NULL ) { continue; }
       *is_match = atoi(s);  

	s=strtok(NULL, " \t\n");
	if ( s==NULL ) { continue; }
       *spos=atoi(s);

	s=strtok(NULL, " \t\n");
	if ( s==NULL ) { continue; }
       *schar=s[0];

	s+=strlen(s)+1;
        
        /* Ignore base positions where consensus base is 'N' or 'X' */
        if (((*cchar) == 'N') || ((*cchar) == 'X')) 
        {
            continue;
        }

        /* If sample == '-', there are no training parameters. */
        if ( (*schar) == '-' ) { 
           *is_match = 0;
            return 1; 
        }

        if ((Btk_atod(&s, phr3)     != 1) ||
            (Btk_atod(&s, phr7)     != 1) ||
            (Btk_atod(&s, psr7)     != 1) ||
            (Btk_atod(&s, pres)     != 1))
        { 
            fprintf(stderr,
                "line %d: missing/garbled parameter; skipping\n", line);
            continue;
        }

        return 1;
    }
    return 1;
}

/*******************************************************************************
 * Function: calc_qv
 * This function calculates a quality value, given statistics about a group
 * of base calls.  Its synopsis is:
 *
 * qv = calc_qv(wrong, total)
 *
 * where
 *	wrong	is the number of incorrect base calls
 *	total	is the total number of right and wrong base calls
 *
 *	qv	is the computed quality value, WITHOUT rounding
 *******************************************************************************
 */
double
calc_qv(int wrong, int total)
{
    double d;
    double qv;


    if (wrong == 0) {
	d = 1.0 / (total + 1.0);
    }
    else {
	d = (double)wrong / total;
    }
    qv = -10.0 * log10(d);

    return(qv);
}

/*******************************************************************************
 * Function: calc_qvP
 * Purpose:  calculated average predicted quality value in a bin
 *******************************************************************************
 */
double
calc_qvP(double prob, int Nprob)
{
    long double d;
    double qv;


    if (prob == 0.000000) {
        d = 1.0 / (Nprob + 1.0);
    }
    else {
        d = (double)prob / Nprob;
    }
    qv = -10.0 * log10(d);

    return(qv);
}

/*******************************************************************************
 * Function: getbase_phred
 * This function reads a single valid base's parameters from the training file
 * (stdin).  Its synopsis is:
 *
 * success = getbase(spos, schar, cpos, cchar, is_match, quality_values)
 *
 * where
 *      spos            is where the sample position should go
 *      schar           is where the sample character (ACTG-) should go
 *      cpos            is where the consensus position should go
 *      cchar           is where the consensus character (ACTG-) should go
 *	phr3		are where the double parameters should go (by ref)
 *	phr7
 *	psr
 *	pre
 *	is_match	is where the int parameter should go (by ref)
 *
 *	success		is 1 if a valid base is read, 0 otherwise
 *
 * ASSUMES: input lines consist of
 *             sample_pos sample_char consensus_pos cons_char training_params
 *
 * NOTE: if schar == '-', the training parameters are not set.
 *******************************************************************************
 */
int
getbase_phd(int* cpos, char* cchar, int *is_match, int *spos, char* schar,
    int *quality_values)
{
    char linebuf[MAXLINE], *s;
    double qv;

    for (;;) {
	if (fgets(linebuf, MAXLINE, stdin) == NULL) {
	    return 0;
	}

	/* Ignore all white space lines and comments */
	if (strspn(linebuf, " \t\r\n") == strlen(linebuf)) {
	    continue;
	}
#if 0
        if ((linebuf[0] == '#') && (linebuf[2]=='F'))
            fprintf(stderr, "%s", linebuf);
#endif
	if ((linebuf[0] == '#')
	|| ((linebuf[0] == '/') && (linebuf[1] == '*'))
	|| (linebuf[0] == ';'))
	{
	    continue;
	}

	/* Get rid of 1) sample position, 2) sample base,
	 * 3) consensus position and 4) consensus base, 5)quality value.
	 */
        s=strtok(linebuf, " \t\n");
        if ( s==NULL ) { continue; }
       *cpos=atoi(s);

        s=strtok(NULL, " \t\n");
        if ( s==NULL ) { continue; }
       *cchar=s[0];

        s=strtok(NULL, " \t\n");
        if ( s==NULL ) { continue; }
       *is_match=atoi(s);
 
	    s=strtok(NULL, " \t\n");
    	if ( s==NULL ) { continue; }
       *spos=atoi(s);

    	s=strtok(NULL, " \t\n");
    	if ( s==NULL ) { continue; }
       *schar=s[0];

    	s+=strlen(s)+1;

        if (Btk_atod(&s, &qv) != 1) {
            continue;
	}
       *quality_values = (int)qv;
	return 1;
    }
}

