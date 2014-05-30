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
 * $Id: trainphd.c,v 1.26 2009/01/13 20:44:52 gdenisov Exp $
 *
 * This file takes as input one or more consensus files and one or more
 * phd files, reads bases from the phd file(s), aligns the fragment(s)   
 * against the respective consensus(es), and outputs a trainphred file.  
 *
 * The trainphred file consists of one line for each base in the alignment,
 * and includes the information:
 *
 * Fragment_position Fragment_character Consensus_position Consensus_character
 *     quality value, is_match(0/1) 
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
#include <unistd.h>
#include <stdint.h>

#ifdef OSF1
extern __const char *__const sys_errlist[];
#endif

#include "Btk_qv.h"
#include "util.h"
#include "Btk_qv_data.h"
#include "Btk_lookup_table.h"
#include "context_table.h"
#include "tracepoly.h"
#include "train.h"              /* needs util.h */
#include "SFF_Toolkit.h"
#include "Btk_compute_tpars.h"  /* needs train.h */
#include "Btk_match_data.h"
#include "Btk_compute_match.h"
#include "Btk_qv_io.h"
#include "train_data.h"

#define DEBUG 0
#define BUFLEN 1000
#define FASTA_LEN 1000
#define MAX_FILENAME_LEN 100
#define MAX_COMMAND_LEN 200
#define NAME_PROJECTFILE 3

static char InputName[BUFLEN]; /* Path of dir or project file. */
static int InputType;

static char ConsensusName[BUFLEN]; /* Name of the Consensus file. */
static int  ConsensusSpecified;    /* Whether the user has specified a name. */

static char OutputName[BUFLEN];    /* Name of the Output file. */
static int  OutputSpecified;       /* Whether the user has specified a name. */

static char PrimerName[BUFLEN];    /* Name of the Primer file. */
static int  PrimerSpecified;       /* Whether the user has specified a name. */

static char VectorName[BUFLEN];    /* Name of the Vector file. */
static int  VectorSpecified;       /* Whether the user has specified a name. */

static char SiteName[BUFLEN];      /* Name of the Site file. */
static int  SiteSpecified;         /* Whether the user has specified a name. */

static double RepeatFraction;      /* Any sequence whose score is at least
                                    * RepeatFraction*HighScore
                                    * is considered a repeat.
                                    */
static double MaxFractionOfErrors; /* Fragments which alignment with consensus
                                    * contain more than this fraction of errors,
                                    * will be rejected
                                    */
static int hist; /* output "Histogram" of # of bases of a given QV */
static int read_direction;
static double min_portion_aligned; /* at least this part of read sdhould 
                                    * be aligned for the read bases to be used 
                                    */
static int min_read_length;
static int pyrosequencing;
static int output_alternative_bcalls;

/* Scores and penalties. */
static int Match;
static int MisMatch;
static int GapInit;
static int GapExt;

static char tab_file_name[BUFLEN];
static char tab_dir_name[BUFLEN];

/* *****************************************************************************
 */
static void show_usage(char** argv)
{
    fprintf(stderr, "\nVersion: %s\n",TT_VERSION);
    (void)fprintf(stderr, "\
usage: %s [ -h ]\n\
        [ -C <ref_seq_file> ]    [ -V <vector> ]\n\
        [ -P <primer> ]          [ -S <site> ]\n\
        [ -M <match_premium> ]   [ -X <subst_penalty> ] [ -G <gap_penalty> ]\n\
        [ -r <repeat_fraction> ] [ -f <max_fraction_of_errors> ] \n\
        [ -a <min_portion_aligned> ] [ -l <min_read_length> ] \n\
        [ -t <tab_dir> ]         [ -o <output_file> ]\n\
        <phd_file(s)> || -d <input_phd_dir> || -j <project_file>\n",
                   argv[0]);
}

static void
help_message(int argc, char *argv[])
{
    (void)fprintf(stderr, "\t-h                    (Help) This message\n");
    (void)fprintf(stderr, "\t-C <consensusfile>    Specify the name of the FASTA file which contains \n");
    (void)fprintf(stderr, "\t                      the consensus sequence \n");
    (void)fprintf(stderr, "\t-V <vector>           Specify the name of the FASTA file which contains the \n");
    (void)fprintf(stderr, "\t                      vector sequence \n");
    (void)fprintf(stderr, "\t-P <primer>           Specify the name of the FASTA file which contains the \n");
    (void)fprintf(stderr, "\t                      primer sequence \n");
    (void)fprintf(stderr, "\t-S <site>             Specify the name of the FASTA file which contains the  \n");
    (void)fprintf(stderr, "\t                      restriction site sequence \n");
    (void)fprintf(stderr, "\t-M <match>            Specify the match premium     (default is 10) \n");
    (void)fprintf(stderr, "\t-X <mismatch>         Specify the mismatch penalty  (default is 20) \n");
    (void)fprintf(stderr, "\t-G <gap_penalty >     Specify the gap initiation or extension penalty (default is 40) \n");
    (void)fprintf(stderr, "\t-r <repeat_fraction>  Specify the repeat_fraction   (default is 0.85) \n");
    (void)fprintf(stderr, "\t-f <max_frac_of_err>  Specify the allowable fraction of errors within the  \n");
    (void)fprintf(stderr, "\t                      best alignment region. Default is 0.1. If actual  \n");
    (void)fprintf(stderr, "\t                      fraction of errors exceeds this vale, the fragment  \n");
    (void)fprintf(stderr, "\t                      will be rejected (=not used in training process) \n");
    (void)fprintf(stderr, "\t-a <min_portion_aligned> Instructs trainphd to ignore a read if less than \n");
    (void)fprintf(stderr, "\t                      the specified portion of its bases is aligned with \n");
    (void)fprintf(stderr, "\t                      reference sequence (default is 0)\n");
    (void)fprintf(stderr, "\t-l <min_read_length>  Instructs trainphd to ignore a read of length less \n");
    (void)fprintf(stderr, "\t                      than specified number of basepairs (0 by default) \n");
    (void)fprintf(stderr, "\t-o <output_file>      Specify the name of the output file. By default, \n");
    (void)fprintf(stderr, "\t                      the output will be made to stdout   \n");
    (void)fprintf(stderr, "\t-t <tab_dir>          Instructs trainphd extract alternative base calls from TAB \n");
    (void)fprintf(stderr, "\t                      files stored in directory tab_dir and to output these \n");
    (void)fprintf(stderr, "\t                      calls, together with quality value for each alternative call \n");
    (void)fprintf(stderr, "\t-d <dir>              Read the input PHD files from specified directory \n");
    (void)fprintf(stderr, "\t-j <projectfile>      Specify the name of the projectfile which comprises \n");
    (void)fprintf(stderr, "\t                      two columns: the full path to the FASTA file which contains \n");
    (void)fprintf(stderr, "\t                      the consensus sequence, followed by the full path to the\n");
    (void)fprintf(stderr, "\t                      sample file \n");
}

/*******************************************************************************
 * File: release
 *******************************************************************************/
static int 
release3(char* called_bases, int* called_peak_locs, uint8_t* quality_values )
{
/* int i; */
 
   FREE(called_bases);
   FREE(called_peak_locs);
   FREE(quality_values);
return SUCCESS;
}

/*******************************************************************************
 * File: process_phd_file
 *
 *******************************************************************************/
static int 
process_phd_file(Align_params *alp, Align_params *alpIUB, char* phd_file_name, 
    char* tab_file_name, FILE *fout, char * consensus_name, 
    Contig* consensusc, Contig* consensusrc, Vector* vector, BtkMessage* message)
{
    int            r=0, i, count_sub, count_del, count_ins, reject_read;
    int            num_called_bases, num_alt_bases=0;
    char          * called_bases=NULL;
    int           *called_peak_locs=NULL;
    uint8_t       *quality_values=NULL;
    char          *bases2=NULL;
    int           *peak_locs2=NULL;
    uint8_t       *qvs2=NULL;
    int           *bind2=NULL;
    int            num_align;
    Contig         fragment;
    Range          align_range;
    Range          clear_range;
    Align          best_alignment;
    Align          vector_start;
    Align          vector_end;
    int            alignment_size;

    /* Initialize with default score (S.L.T. 8/20/2001) */
    best_alignment.score = -1;
    vector_start.score = -1;
    vector_end.score = -1;

    if (phd_file_name[0] > 0)
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "File: %s\n", phd_file_name);
        fprintf(stderr, "Reference: %s\n", consensus_name);
        fprintf(fout, "# File: %s\n", phd_file_name);
        fprintf(fout, "# Reference: %s\n", consensus_name);
	
	/* Allocate memory for data array */
	num_called_bases=MAX_BASES_LEN;
	called_bases = CALLOC(char, num_called_bases);
	MEM_ERROR(called_bases);
 
	called_peak_locs = CALLOC(int, num_called_bases);
	MEM_ERROR(called_peak_locs);
 
	quality_values = CALLOC(uint8_t, num_called_bases);
	MEM_ERROR(quality_values);

	if (Btk_read_phd_file(phd_file_name, called_bases, quality_values, 
            called_peak_locs, &num_called_bases, message) != SUCCESS) {
            fprintf(stderr, "Cannot read phd file %s\n", phd_file_name);
	    return r;       
        }

        // Special case for 454 reads: if read contains "n", reject it
        reject_read = 0;
        for (i=0; i<num_called_bases; i++)
        {
            if (called_bases[i] == 'N' || called_bases[i] == 'n')
                reject_read = 1;
        }

        if (contig_create(&fragment, called_bases, num_called_bases, 
            quality_values, message)
            != SUCCESS)
        {
            fprintf(stderr, "Cannot create contig\n");
            return ERROR;
        }

        if (Btk_compute_match(alp, alpIUB, consensusc, consensusrc, &fragment,
            &num_align, &align_range, RepeatFraction, &best_alignment,
            vector, &vector_start, &vector_end, &clear_range, 
            1, message) != SUCCESS)
        {
            fprintf(stderr, "Cannot compute match\n");
            return ERROR;
        }

        if (num_align == 0) {
            fprintf(fout, "# No good alignments - ignore file %s\n",
                   phd_file_name);
            fprintf(stderr, "No good alignments - ignore file %s\n",
                   phd_file_name);
        } 
        else if ( (num_align != 1) || (best_alignment.score == -1) ) {
            fprintf(fout, "# Possible repeat !!! - ignore file %s\n",
                   phd_file_name);
            fprintf(stderr, "Possible repeat !!! - ignore file %s\n",
                   phd_file_name);
        }
        else if (pyrosequencing && reject_read) // use only reads with Ns
        {
            fprintf(fout, "# 454 read contains Ns - ignore file %s\n",
                   phd_file_name);
            fprintf(stderr, "454 read contains Ns - ignore file %s\n",
                   phd_file_name);
        } 
        else {
            /* Count number of mismatches */
            count_sub = count_del = count_ins = 0; 
            alignment_size=0;
            for (i=0; i<best_alignment.trace_len; i++) {
                if (best_alignment.trace_qchar[i]=='-')
                    count_del++;
                else if (best_alignment.trace_dir[i] == '1') 
                    count_ins++;
                else if (best_alignment.trace_qchar[i] !=
                         best_alignment.trace_dchar[i])
                    count_sub++;
                alignment_size++;
            }

            /* Append the results to train file only if number of
             * mismatches is < MaxFractionOfErrors
             */
            if (((double)(count_sub+count_ins+count_del) <
                 (double)alignment_size * MaxFractionOfErrors)
                                                                &&
                 ((double)alignment_size >=
                  (double)num_called_bases * min_portion_aligned)
                                                                &&
                 (num_called_bases >= min_read_length)) 
            {
                if (tab_file_name[0] > 0)
                {
                    num_alt_bases = MAX_BASES_LEN;
                    bases2 = CALLOC(char, num_alt_bases);
                    MEM_ERROR(bases2);

                    peak_locs2 = CALLOC(int, num_alt_bases);
                    MEM_ERROR(peak_locs2);

                    qvs2 = CALLOC(uint8_t, num_alt_bases);
                    MEM_ERROR(qvs2);

                    bind2 = CALLOC(int, num_alt_bases);
                    MEM_ERROR(bind2);

                    if (Btk_read_tab_file(tab_file_name, bases2, qvs2,
                        peak_locs2, bind2, &num_alt_bases, message) != SUCCESS) {
                        fprintf(stderr, "Cannot read tab file %s\n", tab_file_name);
                        return r;
                    }
                }

                if (append_to_trainphred_data(fout, num_called_bases, 
                    called_peak_locs, clear_range, align_range, 
                    &vector_start, &best_alignment, &vector_end, 0,
                    QVMAX(best_alignment.trace_dpos[0], 
                    best_alignment.trace_dpos[best_alignment.trace_len-1]),
                    quality_values, num_alt_bases, bases2,
                    qvs2, bind2, message) != SUCCESS)
                {
                    fprintf(stderr, "Cannot append_to_train_data\n");
                    return ERROR;
                }
                if (tab_file_name[0] > 0)
                {
                    FREE(bases2);
                    FREE(peak_locs2);
                    FREE(qvs2);
                    FREE(bind2);
                }
            }
            else {
                fprintf(fout, "# BAD PROCESSING!!! - ignore file %s\n",
                        phd_file_name);
                fprintf(stderr, "BAD PROCESSING!!! - ignore file %s\n",
                        phd_file_name);
            }
        }
 
        error:
        ; /* make compiler happy */
        }
  
    /* added quality_values to list of vars. to be freed (S.L.T. 8/20/2001) */ 
        if (release3(called_bases, called_peak_locs, 
                     quality_values) != SUCCESS) {
            fprintf(stderr, "Cannot release\n");
            return ERROR;
        }
 
    if (contig_release( &fragment, message ) != SUCCESS) {
       fprintf(stderr, "Cannot release contig\n");
       return ERROR;
    }

    if ( (num_align == 1) && (best_alignment.score != -1) ) {
        if (align_release( &best_alignment, message ) != SUCCESS) {
            fprintf(stderr, "Cannot release align\n");
            return ERROR;
        }
    }

    /* Empirical observation made on 8/20/2001 that:
     * if input parm "vector" == NULL, then vector_start.score 
     * was uninitialized at this point which caused a core dump.
     */
    if (vector_start.score != -1) {
        if (align_release( &vector_start, message ) != SUCCESS) {
            fprintf(stderr, "Cannot release align\n");
            return ERROR; 
        }
    }

    if (vector_end.score != -1) {
        if (align_release( &vector_end, message ) != SUCCESS) {
            fprintf(stderr, "Cannot release align\n");
            return ERROR; 
        }
    }

    return SUCCESS;
}

/******************************************************************************
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
 *******************************************************************************/
static void
error(char *label, char *msg, int err)
{
     if (err > 0) {
        fprintf(stderr, "%s: %s: %d: %s\n", label, msg, err,
                  strerror(err));
    } else
        fprintf(stderr, "%s: %s\n", label, msg);
}

/*******************************************************************************
 * File: populate_histogram   
 *
 *******************************************************************************/
static int
populate_histogram(char* phd_file_name, long *histogram, BtkMessage* message)
{
    int            j;
    int            num_called_bases;
    char          *called_bases=NULL;
    int           *called_peak_locs=NULL;
    uint8_t       *quality_values=NULL;

    if (phd_file_name[0] > 0)
    {
        /* Allocate memory for data array */
        num_called_bases=MAX_BASES_LEN;
        called_bases = CALLOC(char, num_called_bases);
        MEM_ERROR(called_bases);

        called_peak_locs = CALLOC(int, num_called_bases);
        MEM_ERROR(called_peak_locs);

        quality_values = CALLOC(uint8_t, num_called_bases);
        MEM_ERROR(quality_values);

        fprintf(stderr, "%s\n", phd_file_name);

        if (Btk_read_phd_file(phd_file_name, called_bases, quality_values,
            called_peak_locs, &num_called_bases, message) != SUCCESS) 
        {
            fprintf(stderr, "Cannot read phd file %s\n", phd_file_name);
            return ERROR;
        }

        for (j=0; j<num_called_bases; j++)
        {
            if (quality_values[j] < 100)
                histogram[quality_values[j]]++;
            else
                histogram[99]++;
        }
    }
    return SUCCESS;
error:
    return ERROR;
}

/******************************************************************************
 *      File: populate_histogram_from_phd_dir                 
 ******************************************************************************/
static int
populate_histogram_from_phd_dir(char *dir, long *histogram, BtkMessage *message)
{
    DIR    *d;
    struct  dirent *de;
    char    path[MAXPATHLEN];
    struct  stat statbuf;
    int     r=0;

    if ((d = opendir(dir)) == NULL) {
        error(dir, "couldn't open dir", errno);
        return ERROR;
    }

    while ((de = readdir(d)) != NULL) 
    {
        (void)sprintf(path, "%s/%s", dir, de->d_name);
        if (stat(path, &statbuf) != 0) 
        {
            error(path, "can't stat", errno);
            goto error;
             r = ERROR;
        }
        if (statbuf.st_mode & S_IFDIR) 
        {
            /* skip subdirectories */
            continue;
        }
        populate_histogram(path, histogram, message);
    }

    (void)closedir(d);
    return SUCCESS;

error:
    (void)closedir(d);
    return r;
}

/*******************************************************************************
 * File: output_histogram
 *
 *******************************************************************************/
static void 
output_histogram(long *histogram)
{
    int i, cum=0, sum=0;
    FILE *fp = fopen( "Histogram", "w" );
    for (i = 0; i < 100; i++)
    {
        sum += histogram[i];
    }
 
    for (i = 0; i < 100; i++) 
    {
        cum += histogram[i];
        (void)fprintf(fp, "%2d %6ld %2.3f\n", i, histogram[i], (float)cum/(float)sum);
    }
    fclose(fp);
}



/******************************************************************************
 * This function processes all files contained in the specified directory.
 * Its synopsis is:
 *
 * result = process_dir(align_pars, align_pars_IUB, fout, cons, rev_comp, dir, vector, message)
 *
 * where
 *      align_pars      is the address of an Align_params data structure
 *      fout            is the (already open) output stream
 *      cons            is the address of the consensus Contig
 *      rev_comp        is the address of the reverse complement Contig
 *      dir             is the name (path) of a directory
 *      vector          is the address of the vector, or NULL if there is none
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 ******************************************************************************/
static int
process_dir( Align_params *ap, Align_params *apIUB, FILE* fout,
    char *consensus_name, Contig* consensus, Contig* rev_comp, char *dir, 
    char *tab_dir_name, Vector *vector, BtkMessage *message)
{
    DIR *d;
    struct dirent *de;
    char path[MAXPATHLEN];
    char path2[MAXPATHLEN];
    char string[MAXPATHLEN];
    char *name;
    struct stat statbuf;
    int r=0;

    if ((d = opendir(dir)) == NULL) {
        error(dir, "couldn't open dir", errno);
        return ERROR;
    }

    while ((de = readdir(d)) != NULL) {
        if (strstr(de->d_name, ".phd.1") == NULL)
            continue; 
        (void)sprintf(path, "%s/%s", dir, de->d_name);
        if (stat(path, &statbuf) != 0) {
            error(path, "can't stat", errno);
            goto error;
             r = ERROR;
        }
        if (statbuf.st_mode & S_IFDIR) {
            /* skip subdirectories */
            continue;
        }

        if (output_alternative_bcalls)
        {
#ifdef __WIN32
            if ((name = strrchr(path, '\\')) != NULL) {
#else
            if ((name = strrchr(path, '/')) != NULL) {
#endif
                name++;
            }
            (void)sprintf(string, "%s/%s", tab_dir_name, name); // name ends in "phd.1"
            if ((name = strstr(string, "phd.1")) != NULL)
            {
                printf("name= %s len= %d\n", name, (int)strlen(name));
                strncpy(path2, string, strlen(string)-strlen(name)); // remove the phd.1 suffix
                path2[strlen(string)-strlen(name)] = '\0';
            }
            strcat(path2, "tab");
        }
        fprintf(stderr, "path= %s path2= %s\n", path, path2);
        if (process_phd_file(ap, apIUB, path, path2,
            fout, ConsensusName, consensus, rev_comp, vector, message) 
            != SUCCESS )
        {
            goto error;
        }
    }

    (void)closedir(d);
    return SUCCESS;

error:
    (void)closedir(d);
    return r;
}

/******************************************************************************
 * This function processes all files contained in the specified project file.
 * Its synopsis is:
 *
 * result = process_projectfile(align_pars, align_pars_IUB, fout, projectFile, 
 *          vector, message)
 *
 * where
 *      align_pars      is the address of an Align_params data structure
 *      align_pars_IUB  is the address of an Align_params data structure for
 *                      aligning vector with sample
 *      fout            is the (already open) output stream
 *      projectFile     is the name (path) of a project file
 *      vector          is the address of the vector, or NULL if there is none
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *******************************************************************************/
static int
process_projectfile(Align_params *ap, Align_params *apIUB, FILE* fout,
    char *projectFile, Vector* vector, BtkMessage *message)
{
    char    buffer1[BUFLEN];
    char    buffer2[BUFLEN];
    char    previousConsensusName[BUFLEN];
    char   *current, *name;
    Contig  consensus;
    Contig  rev_comp;
    char    phdFileName[MAXPATHLEN];
    char    tabFileName[MAXPATHLEN];
    char    string[MAXPATHLEN];
    int     have_consensus;
    FILE   *infile;

    if ((infile=fopen(projectFile,"r"))== NULL) {
        sprintf(message->text, "Unable to open project file '%s'\n",
                projectFile);
        return ERROR;
    }

    have_consensus=0;
    buffer2[0]='\0';
    previousConsensusName[0] = '\0';
    current  = buffer1;
    while ((fgets(current, BUFLEN, infile)) != NULL) {
	/* Ignore all white space lines and comments */
	if (strspn(current, " \t\r\n") == strlen(current)) {
	    continue;
	}
	if ((current[0] == '#')                         ||
	   ((current[0] == '/') && (current[1] == '*')) ||
            (current[0] == ';'))
	{
	    continue;
	}

	sscanf(current, "%s %s", ConsensusName, phdFileName);

	fprintf( stderr, "Consensus = %s\n", ConsensusName );

	/* If consensus is different from previous, read it in. */
	if ( strcmp(ConsensusName, previousConsensusName) != 0 ) {
	    /* Release previous, if necessary. */
	    if ( have_consensus ) 
            {
		if (contig_release( &consensus, message ) != SUCCESS) {
                    sprintf(message->text, "Can not release consensus\n");
                    return ERROR;
                }
		if (contig_release( &rev_comp, message ) != SUCCESS) {
                    sprintf(message->text, "Can not release rev_consensus\n");
                    return ERROR;
                }
		have_consensus=0;
            }

	    /* Read consensus file. */
	    if (local_read_fasta( ConsensusName, &consensus, message ) != SUCCESS) {
                sprintf(message->text, "Can not read consensus\n");
                return ERROR;
	    }

	    /* Create reverse complement of consensus. */
            if (contig_get_reverse_comp( &rev_comp, &consensus, message ) 
                != SUCCESS) 
            {
                sprintf(message->text, "Can not get reverse consensus\n");
                return ERROR;
            }
	    have_consensus=1;
        }

        strcpy(previousConsensusName, ConsensusName);

        if (output_alternative_bcalls)
        {
#ifdef __WIN32
            if ((name = strrchr(phdFileName, '\\')) != NULL) {
#else
            if ((name = strrchr(phdFileName, '/')) != NULL) {
#endif
                name++;
            }
            (void)sprintf(string, "%s/%s", tab_dir_name, name); // name ends in "phd.1"
            if ((name = strstr(string, "phd.1")) != NULL)
            {
                strncpy(tabFileName, string, strlen(string)-strlen(name)); // remove the phd.1 suffix
                tabFileName[strlen(string)-strlen(name)] = '\0';
            }
            strcat(tabFileName, "tab");
        }

	fprintf( fout, "# Reference = %s\n", ConsensusName );
        if (process_phd_file(ap, apIUB, phdFileName, tabFileName, fout, 
            ConsensusName, &consensus, &rev_comp, vector, message) 
            != SUCCESS) {
            sprintf(message->text, "Can not process phd file\n");
            return ERROR; 
        }
    }

    fclose(infile);

    /* Clean up. */
    if ( have_consensus ) {
	if (contig_release( &consensus, message ) != SUCCESS) {
            return ERROR; 
        }
	if (contig_release( &rev_comp, message ) != SUCCESS) {
            return ERROR; 
        }
    }

    return SUCCESS;
}

/******************************************************************************
 * File: main
 ******************************************************************************/
int main(int argc, char** argv)
{
    FILE  *fout = stdout;
    Align_params align_pars, align_pars_IUB;
    int r=0, i;
    BtkMessage message;
    Contig consensus, rev_comp;
    Vector vector;
    Vector *vecP = NULL;
 
    if (argc == 1) {
        show_usage(argv);
        exit(2);
    }

    /* Set defaults. */
    ConsensusName[0]    = '\0';
    ConsensusSpecified  = 0;
    OutputName[0]       = '\0';
    OutputSpecified     = 0;
    PrimerName[0]       = '\0';
    PrimerSpecified     = 0;
    VectorName[0]       = '\0';
    VectorSpecified     = 0;
    SiteName[0]         = '\0';
    SiteSpecified       = 0;
    InputName[0]        = '\0';
    InputType           = NAME_FILES;
    tab_dir_name[0]     = '\0';
    RepeatFraction      = 0.85;
    MaxFractionOfErrors = 0.2;
    Match               = 10;
    MisMatch            = -20;
    GapInit             = -40;
    GapExt              = -40;
    read_direction      = 0;
    min_portion_aligned = 0.;
    min_read_length     = 0;
    output_alternative_bcalls = 0;
    static char prev_option = '\0';

    /* Parse the command line */
    while ( (i = getopt(argc, argv, "P:V:S:M:X:G:C:o:d:r:f:j:v:a:l:t:ghp")) != EOF ) {
        switch(i) {
        case 'P':
            PrimerSpecified++;
            (void)strncpy(PrimerName, optarg, sizeof(PrimerName));
	    break;
        case 'V':
            VectorSpecified++;
            (void)strncpy(VectorName, optarg, sizeof(VectorName));
	    break;
        case 'S':
            SiteSpecified++;
            (void)strncpy(SiteName, optarg, sizeof(SiteName));
            break;
        case 'M':
            Match=atoi(optarg);
            if ( Match < 0 ) {
                fprintf(stderr, "Error: Match must not be negative.\n");
                exit(2);
            }
            break;
        case 'X':
            MisMatch=atoi(optarg);
            if ( MisMatch > 0 ) {
                fprintf(stderr, "Error: MisMatch must not be positive.\n");
                exit(2);
            }
            break;
        case 'G':
            GapInit=atoi(optarg);
            GapExt=atoi(optarg);
            if ( GapInit > 0 ) {
                fprintf(stderr, "Error: GapInit must not be positive.\n");
                exit(2);
            }
            prev_option = 'G';
            break;
        case 'C':
            ConsensusSpecified++;
            (void)strncpy(ConsensusName, optarg, sizeof(ConsensusName));
            break;
        case 'f':
            MaxFractionOfErrors=atof(optarg);
            break;
        case 'o':
            OutputSpecified++;
            (void)strncpy(OutputName, optarg, sizeof(OutputName));
            break;
        case 'd':
            InputType = NAME_DIR;
            (void)strncpy(InputName, optarg, sizeof(InputName));
            break;
        case 'j':
            InputType = NAME_PROJECTFILE;
            (void)strncpy(InputName, optarg, sizeof(InputName));
            break;
        case 'p':
            pyrosequencing++;
        case 'r':
            RepeatFraction=atof(optarg);
            if ( (RepeatFraction < 0 ) || (RepeatFraction > 1.0) )
            {
                fprintf(stderr,
                    "Error: RepeatFraction must be in the range [0,1.0]\n");
                exit(2);
            }
            break;
        case 'v':
            read_direction=atoi(optarg);
            if ( read_direction != 0 ) {
                read_direction = read_direction / abs(read_direction);
            }
            break; 
         case 'l':
            min_read_length=atoi(optarg);
            break;
        case 'a':
            min_portion_aligned=atof(optarg);
            break;
        case 'g':
            hist++;
            break;
        case 'h':
            help_message(argc, argv);
            exit(2);
        case 't':
            (void)strncpy(tab_dir_name, optarg, sizeof(tab_dir_name));
            output_alternative_bcalls++;
            break;
        default:
            if (prev_option != 'G')
            {
                show_usage(argv);
                exit(2);
            }
        }
    }

    /*
     * Set line buffering on the status output so that someone monitoring
     * progress can actually see incremental progress.
     */
    (void)setbuf(stderr, NULL);

    if (!hist)
    {
        if (OutputSpecified) {
            /* Set up output file. */
            if ((fout=fopen(OutputName, "w"))== NULL) {
                fprintf(stderr, "Cannot open output file '%s'\n", OutputName);
                exit(ERROR);
            }
        } else {
            (void) fprintf(stderr, "No output file specified.  Using stdout.\n");
            fout=stdout;
        }

        /* Set parameters for alignment with consensus */  
        if (set_alignment_parameters( &align_pars, Match, MisMatch, GapInit,
            GapExt, &message ) != SUCCESS) 
        {
           fprintf(stderr, message.text);
           exit (ERROR);
        }
   
        /* Set parameters for alignment with vector */  
        if (set_alignment_parameters_IUB( &align_pars_IUB, 
            Match * VECTOR_MATCH_MULTIPLIER, MisMatch, 
            GapInit, GapExt, &message ) != SUCCESS)
        {
           fprintf(stderr, message.text);
           exit (ERROR);
        }

        /* Print parameter values to output file. */
            fprintf( fout, "# Match     = %d\n", Match );
            fprintf( fout, "# MisMatch  = %d\n", MisMatch );
            fprintf( fout, "# GapInit   = %d\n", GapInit);
            fprintf( fout, "# GapExt    = %d\n", GapExt);
            fprintf( fout, "# RepeatFraction = %f\n", RepeatFraction );
            fprintf( fout, "#\n");

        vecP=NULL;
        if ( VectorSpecified && PrimerSpecified && SiteSpecified ) {
        	fprintf(fout, "# Vector = %s\n", VectorName);
        	fprintf(fout, "# Primer = %s\n", PrimerName);
        	fprintf(fout, "# Site = %s\n",   SiteName);
        	if (readVector( VectorName, PrimerName, SiteName, &vector, 
    	        &align_pars_IUB, &message) ==ERROR ) {
      	        fprintf(stderr, message.text);
      	        exit (ERROR);
    	    }
    	    vecP=&vector;
        } else if ( VectorSpecified 
    	      && (! PrimerSpecified) && (!SiteSpecified) ) {
        	fprintf(fout, "# Short Vector = %s\n", VectorName);
        	if (readShortVector( VectorName, &vector, &message) == ERROR ) {
          	    fprintf(stderr, message.text);
     	        exit (ERROR);
    	    }
    	    vecP=&vector;
        } else if (PrimerSpecified || SiteSpecified) {
        	fprintf(stderr, "Error: incomplete specification of vector.\n");
        	fprintf(stderr, "Two ways of specifying vector:\n");
        	fprintf(stderr, 
    		" 1. Vector, Primer and Site.\n 2. Short Vector only.\n");
        	exit(2);
        }
    }

    switch(InputType) {
    case NAME_FILES:

	if (!ConsensusSpecified && !hist) {
	    (void) fprintf(stderr, "%s: no consensus file specified.\n",
			   argv[0]);
	    show_usage(argv);
	    exit(2);
	}
        else if (hist)
        {
            /* Output a histogram of # of bases having a given QV */
            int    j;
            long   histogram[100];
         
            for (j=0; j<100; j++)
                 histogram[j] = 0;

            for(i = optind; i < argc; i++)
            {
                populate_histogram(argv[i], histogram, &message);
            }
            output_histogram(histogram); 
            break;
        }
        else 
        {
            fprintf( fout, "#\n");
	    /* Read consensus file. */
	    if (local_read_fasta( ConsensusName, &consensus, &message ) != SUCCESS)
            {
	        fprintf(stderr, message.text);
	        exit (ERROR);
	    }
	    /* Create reverse complement of consensus. */
            if (contig_get_reverse_comp( &rev_comp, &consensus, &message ) 
                != SUCCESS) 
            {
	        fprintf(stderr, message.text);
	        exit (ERROR);
	    }

            /* Loop for each fragment: */
            for(i = optind; i < argc; i++)
            {
                char *name;
                char string[MAXPATHLEN];
#ifdef __WIN32
                if ((name = strrchr(argv[i], '\\')) != NULL) {
#else
                if ((name = strrchr(argv[i], '/')) != NULL) {
#endif
                    name++;
                }
                
                (void)sprintf(string, "%s/%s", tab_dir_name, name); 
                if (tab_dir_name[0] != '\0')
                {
                    if ((name = strstr(string, "phd.1")) != NULL)
                    {
                        // remove the phd.1 suffix
                        strncpy(tab_file_name, string, strlen(string)-strlen(name));
                        tab_file_name[strlen(string)-strlen(name)] = '\0';
                    }
                    strcat(tab_file_name, "tab");
                }
            
                if (process_phd_file(&align_pars, &align_pars_IUB, argv[i], tab_file_name, fout, 
                    ConsensusName, &consensus, &rev_comp, vecP, &message) != SUCCESS)  
                {
                    fprintf(stderr, "Error code=%d\n", r);
                    fprintf(stderr,"error:%s\n",message.text);
                    return ERROR;
                }
            }
            break;
        }

    case NAME_DIR:

	if (!ConsensusSpecified && !hist) {
	    (void) fprintf(stderr, "%s: no consensus file specified.\n",
			   argv[0]);
	    show_usage(argv);
	    exit(2);
	}
        else if (hist)
        {
            /* Output a histogram of # of bases having a given QV */
            int    j;
            long   histogram[100];

            for (j=0; j<100; j++)
                 histogram[j] = 0;

            populate_histogram_from_phd_dir(InputName, histogram, &message);
            output_histogram(histogram);
            break;
        }
        else
        {
            fprintf( fout, "#\n");
	    /* Read consensus file. */
	    if (local_read_fasta( ConsensusName, &consensus, &message ) != SUCCESS) {
	        fprintf(stderr, message.text);
	        exit (ERROR);
	    }
	    /* Create reverse complement of consensus. */
	    if (contig_get_reverse_comp( &rev_comp, &consensus, &message )
                != SUCCESS) 
            {
	        fprintf(stderr, message.text);
	        exit (ERROR);
	    }

            if (process_dir(&align_pars, &align_pars_IUB, fout, ConsensusName, 
                &consensus, &rev_comp, InputName, tab_dir_name, vecP, &message) != SUCCESS) 
            {
                (void)fprintf(stderr, "%s: %s\n", argv[0], message.text);
	        exit (ERROR);
            }
            break;
        }

    case NAME_PROJECTFILE:
        if (process_projectfile(&align_pars, &align_pars_IUB, fout, 
            InputName, vecP, &message) != SUCCESS) 
        {
            (void)fprintf(stderr, "%s: %s\n", argv[0], message.text);
	    exit (ERROR);
        }
        break;
	
    default:
        (void)fprintf(stderr, "%s: internal error: name type %d\n", argv[0],
                        InputType);
        exit(ERROR);
    }

    /* Clean up. */
    if (!hist)
    {
        if (alignment_parameters_release( &align_pars, &message ) != SUCCESS) {
           fprintf(stderr, message.text);
           exit (ERROR);
        }
 
        if (alignment_parameters_release( &align_pars_IUB, &message ) != SUCCESS) {
           fprintf(stderr, message.text);
           exit (ERROR);
        }

        if ( ( InputType == NAME_FILES ) || ( InputType == NAME_DIR ) ) {
    	    if (contig_release( &consensus, &message ) != SUCCESS) {
	        fprintf(stderr, message.text);
	        exit (ERROR);
	    }

            if (contig_release( &rev_comp, &message ) != SUCCESS) {
	        fprintf(stderr, message.text);
	        exit (ERROR);
            }
        }

        fclose(fout);
    }
    exit(SUCCESS);
}

