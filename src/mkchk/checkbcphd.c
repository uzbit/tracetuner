/* 
 * Copyright (c) 1999-2003 Paracel, Inc.  All rights reserved.
 *
 * $Id: checkbcphd.c,v 1.11 2008/12/15 14:07:24 gdenisov Exp $                   
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
#include "params.h"
#include "lut.h"
#include "check_data.h"

static int use_consensus_pos = 0;
static int max_var = 0;
static int qv_weights = 0;
static char color2base[6]={'A', 'C', 'G', 'T', 'N', '-'};
static int pyrosequencing = 0;

static int BinSize;

static int 
base2color(char b)
{
    if      (b == 'A' || b == 'a') return 0;
    else if (b == 'C' || b == 'c') return 1;
    else if (b == 'G' || b == 'g') return 2;
    else if (b == 'T' || b == 't') return 3;
    else if (b == 'N' || b == 'n') return 4;
    else                           return 5;  // gap '-'
}

/*******************************************************************************
 * Function: errorStatswithQVs
 *******************************************************************************
 */
static int
errorStatswithQVs()
{
    char   schar, cchar, prev_schar='-', prev_cchar=' '; 
    int    quality_value=0, Ins=0, Del=0, Sub=0;
    int    i, nbases, match, qv=0, nmsbases, tot_err,
           Nqvp[MAXNUMBINS+1], qvp[MAXNUMBINS+1], qvo[MAXNUMBINS+1];
    int    wrong[MAXNUMBINS+1], total[MAXNUMBINS+1]; 
    int    maxbin, line=1;
    int    spos, cpos, prev_qv=0, pos=0;
    double ninsertion[MAXNUMBINS+1], 
           ndeletion[MAXNUMBINS+1], 
           nsubstitution[MAXNUMBINS+1];
    int    nsubstitution_max_var[6][6][MAXNUMBINS+1];
    int    ninsertion_max_var[6][MAXNUMBINS+1];
    int    ndeletion_max_var[6][MAXNUMBINS+1];
    double Pp[MAXNUMBINS+1]; 
    int    weight;
    int    homopolymer=0, prev_homopolymer=0;

    /*  
    if ((table = Btk_read_lookup_table(tablefile)) == NULL) {
	(void)fprintf(stderr, "%s: couldn't read lookup table\n", tablefile);
	exit(1);
    }
    
    (void)fprintf(stderr, "%s: %d entries\n", tablefile, table->num_entries);
    */
    (void)memset(wrong, 0, sizeof(wrong));
    (void)memset(total, 0, sizeof(total));

    for (i=0; i<=MAXNUMBINS; i++) {
        Pp[i]=0.;
        total[i]=0;
    }

    maxbin = 0;
    nbases = 0;
    nmsbases = 0;   /* matched or substituted bases */
    fprintf(stderr, "Reading bases from stdout ... \n");

    while (getbase_phd(&cpos, &cchar, &match, &spos, &schar,
                         &quality_value)) 
    {
        line++;
        pos = spos;
        if (schar == prev_schar) homopolymer = 1;
        if (use_consensus_pos)   pos = cpos;

	    if (++nbases % 10000 == 0) {
    	    (void)fprintf(stderr, "\r%d bases read so far...", nbases);
    	}

        weight = qv_weights ? quality_value : 1;

        if ((cchar == 'X') || (cchar == 'N' && schar != '-') || (schar == 'N'))
            continue;

        total[(pos-1)/BinSize] += weight;       
        if (schar != cchar ) {
            wrong[(pos-1)/BinSize] += weight;       
        }

        if ((pos-1)/BinSize + 1 > maxbin)
            maxbin = (pos-1)/BinSize + 1;

        if ( schar == '-' ) {  // Deletion
            Del++;
            ndeletion[(pos-1)/BinSize] += 1.;
            if (max_var)
            {
                int ic = base2color(cchar);
                (ndeletion_max_var[ic][(pos-1)/BinSize]) += weight;
            }
        }
        else
        {
            if ( cchar == '-' ) // Insertion
            {
                Ins++;
                ninsertion[(pos-1)/BinSize] += 1.;
                if (max_var)
                {
                    int is = base2color(schar);
                    ninsertion_max_var[is][(pos-1)/BinSize] += weight;       
                }
            }
            else  // Substitution
            {
                nmsbases++;
                if (schar != cchar) 
                {
                    Sub++;
                    nsubstitution[(pos-1)/BinSize] += 1.;
                    if (max_var)
                    {
                        int ic = base2color(cchar); 
                        int is = base2color(schar);
                        (nsubstitution_max_var[ic][is][(pos-1)/BinSize]) += 
                            weight; 
                    }
                }
                if (prev_schar == '-')
                {
                    qv = QVMIN(quality_value, prev_qv);
                    Pp[(pos-2)/BinSize] = Pp[(pos-2)/BinSize]+exp(-qv/10.0*log(10));
                    Nqvp[(pos-2)/BinSize]++;
                }
            }
        }

        if (pyrosequencing)
        {
            /* Substitutions should be treated differently in case of 454 data.
             * For example, situation 
             * ...AAC... - consensus     or  ...GAC... - consensus
             * ...ACC... - read              ...AAC... - read
             * must be treated as 0.5 insertion and 0.5 deletion,
             * rather than 1 substitution
             */
            if (((homopolymer && !prev_homopolymer) ||
                 (!homopolymer &&  prev_homopolymer))
                && 
                prev_schar != prev_cchar &&
                prev_schar != '-'        &&
                prev_cchar != '-')
            {
                nsubstitution[(pos-1)/BinSize] -= 1.0;
                ninsertion[(pos-1)/BinSize]    += 0.5;
                ndeletion[(pos-1)/BinSize]     += 0.5;
            }
        }

	qv = quality_value;

      	Pp[(pos-1)/BinSize] = Pp[(pos-1)/BinSize]+exp(-qv/10.0*log(10)); 
	Nqvp[(pos-1)/BinSize]++;

        prev_qv    = quality_value;
        prev_schar = schar;
        prev_cchar = cchar;
        prev_homopolymer = homopolymer;
        homopolymer = 0;
    }
    (void)fprintf(stderr, "\r%d bases total.        \n", nbases);

    if (!max_var)
    {
        (void)printf("coord  Err_tot  Err_del  Err_ins  Err_sub ");
        (void)printf(" #Bases  #Del   #Ins  #Sub   #Err  \n");
    }
    else
    {
        char c = '\%';
        (void)printf("coord  #Total ");
        (void)printf("  #Del  #Ins  #Sub  #Var_tot %cVar_del %cVar_ins %cVar_sub %cVar_tot #Sub_tot ID_v S_v\n",
        c, c, c, c);
    }
    for (i = 0; i < maxbin; i++) 
    {
        qvp[i]=(int)(calc_qvP(Pp[i], Nqvp[i])+0.5);
        qvo[i]=(int)(calc_qv(wrong[i], total[i])+0.5);
        if (!max_var)
        {
            (void)printf("%4d   %6.3f   %6.3f   %6.3f   %6.3f %5d     %3d    %3d   %3d   %4d\n",
                i*BinSize+BinSize/2+1,
                total[i]> 0 ? 100.*(ndeletion[i]+ninsertion[i]+nsubstitution[i])/
                (double)total[i] : 0.,
                total[i]> 0 ? 100.*ndeletion[i]/(double)total[i] : 0.,
                total[i]> 0 ? 100.*ninsertion[i]/(double)total[i]: 0.,
                total[i]> 0 ? 100.*nsubstitution[i]/(double)total[i] : 0.,
                total[i], (int)ndeletion[i], (int)ninsertion[i], (int)nsubstitution[i],
                (int)(ndeletion[i]+ninsertion[i]+nsubstitution[i]));
        }
        else
        {
            int j, k; 
            int max_var_substitution=0, max_var_insertion=0,  max_var_deletion=0;
            int sc=-1, ss=-1;
            int is=-1, dc=-1;
            
            for (j=0; j<4; j++) {
                if (max_var_insertion < ninsertion_max_var[j][i])
                {
                    max_var_insertion = ninsertion_max_var[j][i];
                    is = j;
                }
                if (max_var_deletion  <  ndeletion_max_var[j][i])
                {
                    max_var_deletion  =  ndeletion_max_var[j][i];
                    dc = j;
                }
                for (k=0; k<4; k++) {
                    if (max_var_substitution < nsubstitution_max_var[j][k][i])
                    {
                        max_var_substitution = nsubstitution_max_var[j][k][i];
                        sc = j;
                        ss = k;
                    }
                }
            }
            tot_err = max_var_deletion+max_var_insertion+max_var_substitution;
            (void)printf("%4d  %5d  %5d  %5d  %5d %5d   %7.3f  %7.3f   %7.3f  %7.3f  %5d     %c/%c  %c/%c\n",
                i*BinSize+BinSize/2+1,
                total[i], 
                max_var_deletion, max_var_insertion, max_var_substitution,
                max_var_deletion+max_var_insertion+max_var_substitution,
                total[i]> 0 ? 100.*(double)max_var_deletion/(double)total[i] : 0.,
                total[i]> 0 ? 100.*(double)max_var_insertion/(double)total[i] : 0.,
                total[i]> 0 ? 100.*(double)max_var_substitution/(double)total[i] : 0., 
                total[i]> 0 ? 100.*(double)tot_err/total[i] : 0.,
                (int)nsubstitution[i], 
                max_var_substitution > 0 ? color2base[sc] : ' ',
                max_var_substitution > 0 ? color2base[ss] : ' ',
                max_var_deletion > max_var_insertion ? '-' : (max_var_deletion == max_var_insertion ? ' ' : color2base[is]),
                max_var_deletion > max_var_insertion ? color2base[dc] : (max_var_deletion == max_var_insertion ? ' ' : '-'));
        }
    }

    (void)printf("\nNumber of bases=%d\n", nbases);
    (void)printf("\nNumber of matched or substituted bases=%d\n", 
        nmsbases); 
    (void)printf("Number of errors: Ins=%d Del=%d Sub=%d Tot=%d\n",
        Ins, Del, Sub, Ins+Del+Sub);
    return(0);
}

/*******************************************************************************
 * Function: show_usage
 *******************************************************************************
 */
void
show_usage(char *argv[])
{
    (void)fprintf(stderr, "\nVersion: %s\n",TT_VERSION);
    (void)fprintf(stderr, "usage: %s [ -h ] "   
 "    <    <alignment_file>\n", argv[0]);
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
    BinSize = 1; /* Bases 0-99 in one bin, 100-199 in a second, etc. */

    while( (i = getopt(argc, argv, "b:hcpvw") ) != EOF ) 
    {
        switch(i) 
        {
            case 'b':
                BinSize = atoi(optarg);
                if ( BinSize <= 0 ) {
                    fprintf(stderr, "Error: illegal bin size %d\n", BinSize);
                    exit(-1);
                }
                break;
            case 'c':
                use_consensus_pos++;
                break;
            case 'p':
                pyrosequencing++;
                break;
            case 'v':
                max_var++;         
                break;
            case 'w':
                qv_weights++;
                break;
	        default:
	            show_usage(argv);
	            exit(2);
	    }
    }
    setbuf(stderr, NULL);	/* unbuffered */
    errorStatswithQVs();

    return(0);
}

