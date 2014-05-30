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
 * 2.3 2003/11/05 22:59:16
 *
 * This program takes as input 
 * - two phd files or
 * - one or more phd files and a directory of phd files or
 * - two directories of phd files,
 * compares the phd file and outputs the results
 */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>
#include <sys/param.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#ifdef OSF1
extern __const char *__const sys_errlist[];
#endif

#include "util.h"
#include "Btk_qv_data.h"
#include "Btk_qv.h"
#include "Btk_match_data.h"    
#include "Btk_qv_io.h"
#include "Btk_lookup_table.h"
#include "context_table.h"
#include "tracepoly.h"
#include "train.h" 
#include "Btk_compute_tpars.h" 
#include "Btk_compute_match.h"
#include "Btk_default_table.h"
#include "Btk_compute_qv.h"

#define BUFLEN 1000
#define DEBUG 0  
#define FASTA_LEN 1000
#define FP_TO_FN_RATIO 4
#define MAX_FILENAME_LEN 100
#define MAX_COMMAND_LEN 200
#define MAX_BASES_LEN 4000

typedef enum {
     INP_FILES=1,       /* input will come from files named on the cmd line */
     INP_DIR,           /* input will come from all files in a directory */
}  INPTYPE;

static char DirName[BUFLEN]; /* Path of dir or project file. */
static  INPTYPE InputType;
static char OutputName[BUFLEN];
static char Dir2Name[BUFLEN];

/* Scores and penalties. */
static int Match;
static int Substitution;
static int GapInit;
static int GapExt;

static double RepeatFraction = 0.85;
static int    OutputSpecified;

/***********************************************************************/
static void usage(char** argv)
{
    fprintf(stderr, "\nVersion: %s\n", TT_VERSION);
    fprintf(stderr, "\n"
    "usage: %s [ -h ]\n" 
    "     [ -d2 <dir> ]        [ -o <output_file> ]\n"
    "     [ -M <match> ]       [ -X <mismatch> ]\n"
    "     [ -E <gap_ext> ]     [ -I <gap_init > ]\n"
    "     <sample_file(s)>     | -d <input_dir>\n", argv[0]);
}

/***********************************************************************/
static void
help_message(int argc, char *argv[])
{
    /* NOTE: when comparing Phred and TT PHD files, 
     *       Phred file should be specified first and TT - second.
     *       For example,
     *           comparePhdFiles <Phred_phd> <TT_phd>
     *       or
     *           comparePhdFiles -d2 <TT_dir> -d <Phred_dir>
     *
     *       In this case, the output base positions and qvs
     *       will agree with what is seen in the viewer
     *       when command
     *           ttview  <file> <Phred_phd> <TT_phd>
     *       These output positions will correspond to 
     *       Phred data 
     */
    (void)fprintf(stderr, 
    "-h                    (Help) This message\n"
    "-o <output_file>      Specify the name of the output file. By default,\n"
    "-M <match>            Specify the match premium     (default is +10)\n"
    "-X <mismatch>         Specify the mismatch penalty  (default is -20)\n"
    "-E <gap_ext>          Specify the gap extension penalty  (default is -30)\n"
    "-I <gap_init>>        Specify the gap initiation penalty (default is -40)\n"
    "-d2 <dir2>            Compare the input PHD files with PHD files having\n" 
    "                      the same names, but located in directory <dir2> \n"
    "-d  <dir>             Read the input PHD files from specified directory\n");
}


/*
 * This function prints as formatted an error message as it can, then exits
 * if the Force flag isn't set.  Its synopsis is:
 *
 * error(label, msg, err)
 *
 * where
 *      label   is a char* to be printed at the beginning of the line.
 *              Usually the name of the file being processed at the time of
 *              the error, or the name of the program running (argv[0]).
 *      msg     is a char* to a message to print.  Usually describes what
 *              was attempted.
 *      err     is an int error code.  Usually errno.  If 0, then no error
 *              code is printed, nor is any system-provided textual
 *              description of the error.
 */
static void
error(char *label, char *msg, int err)
{
    if (err > 0) {
        fprintf(stderr, "%s: %s: %d: %s\n", label, msg, err,
                  strerror(err));
    } else
        fprintf(stderr, "%s: %s\n", label, msg);
}


/********************************************************************************
 * Function: compare_two_phd_files(file1, file2, fout, message)
 *
 ********************************************************************************
 */
static int
compare_two_phd_files(char *file1, char *file2, FILE *fout,
    BtkMessage *message)
{
    int          num_bases, num_bases2;
    char        *bases = NULL, *bases2 = NULL;
    int         *qv    = NULL, *qv2    = NULL; 
    int         *locs  = NULL, *locs2  = NULL;
    Align_params align_pars, align_pars_IUB;
    Contig       consensus, rev_comp;
    int          r=0, i, count_diff_qv=0;
    int          num_align;
    Contig       fragment;
    Range        align_range;
    Range        clear_range;
    Align        best_alignment;
    Align        vector_start;
    Align        vector_end;
    char        *phd_name1, *phd_name2;

#ifdef __WIN32
    if ((phd_name1 = strrchr(file1, '\\')) != NULL) {
#else
    if ((phd_name1 = strrchr(file1, '/')) != NULL) {
#endif
        phd_name1++;
    }
    else {
        phd_name1 = file1;
    }

#ifdef __WIN32
    if ((phd_name2 = strrchr(file2, '\\')) != NULL) {
#else
    if ((phd_name2 = strrchr(file2, '/')) != NULL) {
#endif
        phd_name2++;
    }
    else {
        phd_name2 = file2;
    }
  
    /* Set alignment parameters globally. */
    if (set_alignment_parameters( &align_pars, Match, Substitution, GapInit  ,
        GapExt, message ) != SUCCESS)
    {
       fprintf(stderr, message->text);
       exit (ERROR);
    }

    if (set_alignment_parameters_IUB( &align_pars_IUB, Match, Substitution,
        GapInit, GapExt, message ) != SUCCESS)
    {
       fprintf(stderr, message->text);
       exit (ERROR);
    }

    /* Allocated memory for compared sequences */ 
    num_bases = num_bases2 = MAX_BASES_LEN;

    bases = CALLOC(char, num_bases);
    MEM_ERROR(bases);
    bases2 = CALLOC(char, num_bases2);
    MEM_ERROR(bases2);

    locs = CALLOC(int, num_bases);
    MEM_ERROR(locs);
    locs2 = CALLOC(int, num_bases2);
    MEM_ERROR(locs2);

    qv  = CALLOC(int, num_bases);
    MEM_ERROR(qv); 
    qv2 = CALLOC(int, num_bases2);
    MEM_ERROR(qv2);

    if (Btk_read_phd_file(file1, bases, qv, locs, &num_bases, message) 
        != SUCCESS) 
    {
        fprintf(stderr, 
        "Cannot read phd file %s\n", file1);
        return r;
    }
    
    if (Btk_read_phd_file(file2, bases2, qv2, locs2, &num_bases2, message) 
        != SUCCESS) 
    {
        fprintf(stderr, 
        "Cannot read phd file %s\n", file2);
        return r;
    }
#if 0
    fprintf(stderr, "bases=\n%s\n", bases);
    fprintf(stderr, "bases2=\n%s\n", bases2);
#endif 
    /* Set 'consensus' to be the second sequence */
    contig_create(&consensus, bases2, strlen(bases2), qv2, message);
    consensus.sequence      = bases2;
    consensus.qv            = qv2;
    consensus.length        = num_bases2;
    consensus.max_length    = MAX_BASES_LEN;

    /* Create reverse complement of consensus. */
    if (contig_get_reverse_comp( &rev_comp, &consensus, message )
        != SUCCESS)
    {
        fprintf(stderr, message->text);
        exit (ERROR);
    }

    /* Initialize with default score */
    align_init(&best_alignment, message);
    align_init(&vector_start, message);
    align_init(&vector_end, message);

    best_alignment.score = -1;
    vector_start.score = -1;
    vector_end.score = -1;
   
    if (contig_create(&fragment, bases, num_bases, qv, message)
        != SUCCESS)
    {
        fprintf(stderr, "Cannot create contig\n");
        return ERROR;
    }

    if (Btk_compute_match(&align_pars, &align_pars_IUB, &consensus, 
        &rev_comp, &fragment, &num_align, &align_range, 
        RepeatFraction, &best_alignment, NULL, 
        &vector_start, &vector_end, &clear_range, 0, message)
        != SUCCESS)
    {
        fprintf(stderr, "Cannot compute match\n");
        return ERROR;
    } 

    if (num_align == 0) {
        fprintf(fout, "\n# No good alignments - ignore file %s\n",
               file1);
        fprintf(stderr, "\nNo good alignments - ignore file %s\n",
               file1);
        return SUCCESS;
    }
    else if ( (num_align != 1) || (best_alignment.score == -1) ) {
        fprintf(fout, 
               "\nFile: %s\nPossible repeat. Alignment failed.\n",
               phd_name1);
        fprintf(stderr, "\nFile: %s\n",
               phd_name1);
        return SUCCESS;
    }
    else {
        /* Count bases with significant differences in QVs */
        int read_pos[BUFLEN];
        for (i=0; i<best_alignment.trace_len; i++) {
            if (best_alignment.trace_qpos[i] <= 0)
                continue;
            if (best_alignment.trace_dpos[i] <= 0)
                continue;
#if 0
            fprintf(stderr, 
                "i=%d trace_qpos[i]=%d trace_dpos[i]=%d\n", 
                 i, best_alignment.trace_qpos[i], 
                    best_alignment.trace_dpos[i]);
#endif
            if (best_alignment.trace_qchar[i] ==
                best_alignment.trace_dchar[i]) 
            {
                int qpos = best_alignment.trace_qpos[i];
                int dpos = best_alignment.trace_dpos[i]; 
                if (((qv[qpos]  <= 15) || (qv2[dpos] <= 15)) 
                      &&
                    ((qv[ qpos] - qv2[dpos] >= 10) ||
                     (qv2[dpos] - qv[ qpos] >= 10))
                   )
                { 
                    if (qv[qpos]  - qv2[dpos] >= 10)
                        read_pos[count_diff_qv] = -qpos;
                    if (qv2[dpos] - qv[qpos]  >= 10)
                        read_pos[count_diff_qv] = qpos;
                    count_diff_qv++;
#if 0
                    fprintf(stderr, 
                       "count_diff_qv=%d qpos=%d\n", 
                        count_diff_qv, qpos);
#endif
                }
            }
        }
   
        /* Print the collected positions */
        fprintf(fout,   "\nFile: %s \n", phd_name1);
        fprintf(stderr, "\nFile: %s \n", phd_name1);
        fprintf(fout, 
            "Total number of bases with significantly different qvs: %d\n",
            count_diff_qv); 
        for (i=0; i<count_diff_qv; i++) {
            if (read_pos[i] >= 0)
                fprintf(fout, "%d(%d), ", read_pos[i]+1, qv[abs(read_pos[i])]);
            if (read_pos[i] < 0)
                fprintf(fout, "%d(%d), ", read_pos[i]-1, qv[abs(read_pos[i])]);
            if ((i > 0) && (i%10 == 0)) 
                fprintf(fout, "\n");
        }
        fprintf(fout, "\n"); 
    }

    FREE(bases);
    FREE(bases2);
    FREE(locs);
    FREE(locs2);
    FREE(qv);
    FREE(qv2);
 
    return SUCCESS;

error:
    return ERROR;
}

/*
 * This function processes all files contained in the specified directory.
 * Its synopsis is:
 *
 * result = process_dir(fout, DirName, Dir2Name, message)
 *
 * where
 *      fout            is the (already open) output stream
 *      DirName         is the name (path) of the directory from where
 *                      to read the input PHD files
 *      Dir2Name        is the name (path) of the directory containing
 *                      PHD files to compare with
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 */
static int process_dir(
    FILE       *fout,
    char       *DirName,
    char       *Dir2Name,
    BtkMessage *message)
{
    DIR *d;
    struct dirent *de;
    char path[MAXPATHLEN], path2[MAXPATHLEN];
    struct stat statbuf;
    int r;

    message->text[0] = 0;

    if ((d = opendir(DirName)) == NULL) {
    	error(DirName, "couldn't open dir", errno);
    	return ERROR;
    }

    while ((de = readdir(d)) != NULL) {
    	(void)sprintf(path, "%s/%s", DirName, de->d_name);
        (void)sprintf(path2, "%s/%s", Dir2Name, de->d_name);
    	if (stat(path, &statbuf) != 0) {
      	    error(path, "can't stat", errno);
      	    closedir(d);
      	    return ERROR;
    	}
    	if (statbuf.st_mode & S_IFDIR) {
      	    /* skip subdirectories */
            continue;
    	}

        r = compare_two_phd_files(path, path2, fout, message);

        if (r != SUCCESS ) {
            fprintf(stderr,"Error processing input: %s\n", message->text);
            fprintf(fout,"# Error processing input: %s\n", message->text);
        }
    }

    (void)closedir(d);
    return SUCCESS;
}

/*******************************************************************************
 * Function: initialize_results
 *******************************************************************************
 */
void
initialize_results(Results *results)
{
    int i;
    results->align_length = 0;
    results->count_ins = 0;
    results->count_del = 0;
    results->count_sub = 0;
    results->count_err = 0;
    results->count_hetero = 0;
    results->count_tp_hetero    = 
    results->count_fp_hetero    = results->count_fn_hetero    = 0;
    results->count_tp_qv_hetero = results->count_fn_qv_hetero = 
    results->count_fp_qv_hetero = 0;
    results->min_qv_cor  = 100;
    results->max_phr3_cor = results->max_phr7_cor =
    results->max_psr7_cor = results->max_pres_cor = 0.; 
    results->count_chunks = 0;
    for ( i=0; i<MAX_CHUNKS; i++ ) {
        statsReset( &results->stats[i] );
        /* assumption: normalized spacing is mostly in [0,2] */
        histoReset( &results->histo[i], 0.0, 2.0 );
    }
}

/*******************************************************************************
 * Function: get_file_name
 * Purpose:  extract the file name from full path
 *
 *******************************************************************************
 */
char *get_file_name(char *file_path)
{
    char *file_name = NULL;

#ifdef __WIN32
    if ((file_name = strrchr(file_path, '\\')) != NULL) {
#else
    if ((file_name = strrchr(file_path, '/')) != NULL) {
#endif
        file_name++;
    }
    else {
        file_name = file_path;
    }
   
    return file_name;
}

/*******************************************************************************
 * Function: main
 *******************************************************************************
 */
int main(int argc, char** argv)
{
    char         *args, *file_name, second_file_name[BUFLEN];
    FILE         *fout;
    int           r, i, j, optind;
    BtkMessage    message;
 
    if (argc == 1) {
    	usage(argv);
    	exit(2);
    }

    /* Set defaults. */
    OutputName[0]       ='\0';
    OutputSpecified     = 0;
    DirName[0]          ='\0';
    InputType           =  INP_FILES;
    Match               =  20;
    Substitution        = -5;
    GapInit             = -4;
    GapExt              = -3;
    Dir2Name[0]         ='\0';
    second_file_name[0] ='\0';

    /*
     * This argument-processing code probably looks a bit weird; if so, it's
     * because it was converted from getopt(), which isn't supported under
     * the Win32-compatible version of gcc we have.
     */
    for (optind = 1; optind < argc; optind++) {
        if (argv[optind][0] != '-') {
            break;
        }
        /* Make sure there is an option argument for those options which are
         * supposed to have an argument
         */
        if ((strcmp(argv[optind], "-d2") == 0) 
            &&
            (optind==argc-1 || argv[optind+1][0]=='-'))
        {
            usage(argv);
            exit(2);
        }

        /* Parse the current argument string (args), starting from the 2nd
         * character (j==1), which is not '-'.
         */
        for (j = 1, args = argv[optind]; args[j] != '\0'; j++) {
            i = args[j];
            switch(i) {

            case 'E':
                GapExt  =atoi(argv[++optind]);
                if ( GapExt   > 0 ) {
                    fprintf(stderr, "Error: GapExt   must not be positive.\n");
                    exit(2);
                }
                break;
           
            case 'd':
                if (strcmp(args, "-d2") == 0) {
                    (void)strncpy(Dir2Name, argv[++optind], sizeof(Dir2Name));
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-d") == 0) {
                    InputType =  INP_DIR;
                    (void)strncpy(DirName, argv[++optind], sizeof(DirName));
                    break;
                }
                else {
                    usage(argv);
                    exit(2);
                }
 
            case 'h':
                help_message(argc, argv);
                exit(2);

            case 'I':
                GapInit  =atoi(argv[++optind]);
                if ( GapInit   > 0 ) {
                    fprintf(stderr, "Error: GapInit   must not be positive\n");
                    exit(2);
                }
                break;
            
            case 'M':
                Match=atoi(argv[++optind]);
                if ( Match < 0 ) {
                    fprintf(stderr, "Error: Match must not be negative\n");
                    exit(2);
                }
                break;

            case 'o':
                OutputSpecified++;
                (void)strncpy(OutputName, argv[++optind], sizeof(OutputName));
                break;

             case 'X':
                Substitution=atoi(argv[++optind]);
                if ( Substitution > 0 ) {
                    fprintf(stderr, 
                    "Error: Substitution penalty must not be positive.\n");
                    exit(2);
                }
                break;

            case '?':
            default:
                usage(argv);
                exit(2);
            }
        }
    }

    /*
     * Set line buffering on the status output so that someone monitoring
     * progress can actually see incremental progress.
     */
    (void)setbuf(stderr, NULL);

    if (OutputSpecified) {
    	/* Set up output file. */
    	if ((fout=fopen(OutputName, "w"))== NULL) {
      	    fprintf(stderr, "Cannot open output file '%s'\n", OutputName);
      	    exit(ERROR);
    	}
    } else {
    	fprintf(stderr, "No output file specified.  Using stdout.\n");
    	fout=stdout;
    }

    switch(InputType) {
    case  INP_FILES:

        if (optind >= argc) {
            fprintf(stderr, 
            "No input phd file(s) or directory of files specified. Exit\n");
            exit(ERROR);
        }
        else if ((optind == argc-1) && (Dir2Name[0] == '\0')) {
            fprintf(stderr, 
            "Only one PHD file and no compare directory specified. Exit\n");
            exit(ERROR);
        }
        else if ((optind < argc-2) && (Dir2Name[0] == '\0')) {
            fprintf(stderr, 
            "More than two PHD files and no compare directory specified. Exit\n");
            exit(ERROR);
        }
        else if ((optind == argc-2) && (Dir2Name[0] == '\0')) {
            r = compare_two_phd_files(argv[argc-2], argv[argc-1], 
                fout, &message);
            if (r != SUCCESS ) {
                fprintf(stderr,"Error processing input: %s\n", message.text);
                fprintf(fout,"# Error processing input: %s\n", message.text);
            }
        }
    	else {
    	    for (i = optind; i < argc; i++) {
                file_name = get_file_name(argv[i]);
                sprintf(second_file_name, "%s/%s", Dir2Name, file_name);

   	        r = compare_two_phd_files(argv[i], second_file_name,
                    fout, &message);

                if (r != SUCCESS ) {
                    fprintf(stderr,"Error processing input: %s\n", message.text);
                    fprintf(fout,"# Error processing input: %s\n", message.text);
                }
            }
    	}
    	break;

    case  INP_DIR:

    	if (process_dir(fout, DirName, Dir2Name, &message) != SUCCESS) {
      	    fprintf(stderr,"Processing error: %s\n", message.text);
  	    exit (ERROR);
    	}
    	break;

    default:
    	fprintf(stderr, "%s: internal error: name type %d\n", argv[0],
		  InputType);
    	exit(ERROR);
    }

    exit(SUCCESS);
}

