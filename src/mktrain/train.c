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
 * $Id: train.c,v 1.22 2009/01/12 22:22:31 gdenisov Exp $                   
 *
 * This file takes as input one or more consensus files and one or more
 * sample files, calls bases in the sample files, aligns the sample files
 * with the respective consensus, and outputs a training file.  The
 * training file consists of one line for each base in the alignment,
 * and includes the information:
 *     Sample_position Sample_character Consensus_position Consensus_character
 *         is_match(0/1) Training parameters
 */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <dirent.h>
#include <stdint.h>
#include <sys/param.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <float.h> 

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
#include "train.h"              /* needs util.h */
#include "SFF_Toolkit.h"
#include "Btk_compute_tpars.h" 
#include "Btk_get_mixed_bases.h"
#include "Btk_compute_match.h"
#include "Btk_default_table.h"
#include "Btk_atod.h"
#include "Btk_compute_qv.h"
#include "train_data.h"
#include "ABI_Toolkit.h"
#include "Btk_process_peaks.h"

#define BUFLEN            1000
#define CALIBRATION          1
#define DEBUG                0  
#define FASTA_LEN         1000
#define FP_TO_FN_RATIO       4
#define MAX_FILENAME_LEN   100
#define MAX_COMMAND_LEN    200
#define MIN_LEN_SWEET_SPOT 300

typedef enum {
     INP_FILES=1,       /* input will come from files named on the cmd line */
     INP_DIR,           /* input will come from all files in a directory */
     INP_PROJECTFILE,   /* input will come a file with one consensus and 
                         * one sample file per line 
                         */
}  INPTYPE;

static int  NUM_PARAMS = 4;
static char InputName[BUFLEN]; /* Path of dir or project file. */
static  INPTYPE InputType;

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
static double RepeatFraction;      /* Any sequence whose score >=RepeatFraction*
                                    * HighScore is considered a repeat 
                                    */
static double MaxFractionOfErrors; /* If the fraction of errors within 
                                    * best-alignment region exceeds this value,
                                    * then the whole fragment sequence will be
                                    * rejected, that is, will not be used 
                                    * in the training process
                                    */
static char status_code[BUFLEN];
static int nocall;                 /* Whether to use ABI base calls */
static int recalln;                /* Whether to only recall 'N's to the best guess*/
static int recallndb;              /* Whether to only recall 'N's and dye blobs 
                                      to the best guess */
static int ladder;                 
static int dev;
static int edited_bases;           /* Whether to sed edited bases */
static int shift;                  /* whether skip mobiliti shift corrections */
static int gauss;                  /* whether to use gaussian peak shape model*/
static int renorm;                 /* whether skip renormalization of traces */ 
static int respace;                /* whether to respace poorly resolved peaks
                                    */
static int het;                    /* whether to call heterozygotes */
static int mix;                    /* whether to call mixed bases */
static int lut;                    /* which lookup table to use */
static int inp_phd;                /* whether to read original bases and locations
                                    *from phd file rather than from sample file
                                    */
static char inp_phd_dir[BUFLEN];   /* Name of input phd directory */
static float min_ratio;            /* a threshold ratio of the current to 
                                    * called peak height under which the 
                                    * current peak will be considered noise
                                    */
static double p1, p2, p3, p4;      /*Some parameters, which are passed through*/
                                   /* the command line */

static int min_read_length;
static double min_portion_aligned; /* at least this part of read sdhould
                                    * be aligned for the read bases to be used
                                    */
/* Scores and penalties. */
static int Match;
static int MisMatch;
static int GapInit;  
static int GapExt;

/* 
 *  Global variables to summarize number of input files.  These are set in
 *  process_sample_file and output as comments at the end of the run.
 */
int Count_input_files = 0;  
int Count_file_errors = 0;  
int Count_processing_errors = 0;
int Count_possible_repeats  = 0;
int Count_no_alignments     = 0;

/* ------------------------------------------------------------------- */
static void usage(char** argv)
{
    (void)fprintf(stderr, "\nVersion: %s\n",TT_VERSION);
    (void)fprintf(stderr, "usage: %s [ -h ] \n", argv[0]);
    (void)fprintf(stderr, "\t[ -nocall ] [ -recalln ] [ -edited_bases ] \n");
    (void)fprintf(stderr, "\t[ -het    ] [ -mix ]     [ -min_ratio <ratio> ]\n");
    (void)fprintf(stderr, "\t[ -C <ref_seq_file>]     [ -V <vector> ]\n");
    (void)fprintf(stderr, "\t[ -P <primer> ]          [ -S <site> ]\n");
    (void)fprintf(stderr, "\t[ -M <match_premium> ]   [ -X <subst_penalty> ] [ -G <gap_penalty> ] \n");
    (void)fprintf(stderr, "\t[ -fr <repeat_frac> ]    [ -fe <max_frac_of_errors> ]\n");
    (void)fprintf(stderr, "\t[ -a <min_portion_aligned> ] [ -l <min_read_length> ]\n");
    (void)fprintf(stderr, "\t[ -ipd <dir> ]         [ -o <output_file> ]\n");
    (void)fprintf(stderr,
        "        <sample_file(s)> || -d <input_dir(s)> || -p <project_file>\n");
}

/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
static void usage_dev(char** argv)
{
    (void)fprintf(stderr, "\nVersion: %s\n",TT_VERSION);
    (void)fprintf(stderr, "usage: %s [ -h ] [ -Q ] \n", argv[0]);
    (void)fprintf(stderr, "\t[ -nocall | -recalln | -recallndb | -ladder]\n");
    (void)fprintf(stderr, "\t[ -shift ][ -convolved ] [ -edited_bases ]\n");
    (void)fprintf(stderr, "\t[ -renorm ]            [ -respace ]\n");
    (void)fprintf(stderr, "\t[ -het ]      [ -mix ] [ -min_ratio <ratio> ]\n");
    (void)fprintf(stderr, "\t[ -p1 <double_value> ] [ -p2 <double_value> ]\n");
    (void)fprintf(stderr, "\t[ -p3 <double_value> ] [ -p4 <double_value> ]\n");
    (void)fprintf(stderr, "\t[ -C <ref_seq_file>]  [ -V <vector> ]\n");
    (void)fprintf(stderr, "\t[ -P <primer> ]        [ -S <site> ]\n");
    (void)fprintf(stderr, "\t[ -M <match> ]         [ -X <mismatch> ] [ -G <gap_penalty> ]\n");
    (void)fprintf(stderr,
        "\t[ -fr <repeat_fraction> ] [ -fe <max_fraction_of_errors> ]\n");
    (void)fprintf(stderr, "\t[ -ipd <dir> ]         [ -o <output_file> ]\n");
    (void)fprintf(stderr,
        "      <sample_file(s)> || -d <input_dir(s)> || -p <project_file>\n");
}

/* ------------------------------------------------------------------- */
static void
help_message(int argc, char *argv[])
{
    (void)fprintf(stderr, "\t-h                    (Help) This message\n");
    (void)fprintf(stderr, "\t-nocall               Disable base recalling and just use the original \n");
    (void)fprintf(stderr, "\t                      called bases read from the input sample file \n");
    (void)fprintf(stderr, "\t-recalln              Disable adding bases to or deleting from the \n");
    (void)fprintf(stderr, "\t                      original called sequence. Only recall Ns \n");
    (void)fprintf(stderr, "\t-recallndb            Disable adding bases to or deleting from the \n");
    (void)fprintf(stderr, "\t                      original called sequence. Only recall Ns and dye blobs \n");
    (void)fprintf(stderr, "\t-ladder               Similar to -recalln, but all bases, not only Ns\n");
    (void)fprintf(stderr, "\t                      will be recalled from original locations\n");
    (void)fprintf(stderr, "\t-edited_bases         Start base recalling from the ABI's edited bases \n");
    (void)fprintf(stderr, "\t-het                  Call heterozygous bases \n");
    (void)fprintf(stderr, "\t-mix                  Call mixed        bases \n");
    (void)fprintf(stderr, "\t-min_ratio <ratio>    Override the default threshold ratio of heights of \n");
    (void)fprintf(stderr, "\t                      the lowest peak to the highest peak at a given \n");
    (void)fprintf(stderr, "\t                      position \n");
    (void)fprintf(stderr, "\t-fr <repeat_fraction> Specify the repeat_fraction   (default is 0.85) \n");
    (void)fprintf(stderr, "\t-fe <max_frac_of_err> Specify the allowable fraction of errors within the  \n");
    (void)fprintf(stderr, "\t                      best alignment region. Default is 0.1. If actual  \n");
    (void)fprintf(stderr, "\t                      fraction of errors exceeds this vale, the fragment  \n");
    (void)fprintf(stderr, "\t                      will be rejected (=not used in training process) \n");
    (void)fprintf(stderr, "\t-o <output_file>      Specify the name of the output file. By default, \n");
    (void)fprintf(stderr, "\t                      the output will be made to stdout   \n");
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
    (void)fprintf(stderr, "\t-d <dir>              Read the input sample files from specified directory \n");
    (void)fprintf(stderr, "\t-p <projectfile>      Specify the name of the projectfile which consists of \n");
    (void)fprintf(stderr, "\t                      two columns. Each line in this file contains the  \n");
    (void)fprintf(stderr, "\t                      full path to the FASTA file which contains the \n");
    (void)fprintf(stderr, "\t                      consensus sequence, followed by the full path to the\n");
    (void)fprintf(stderr, "\t                      sample file \n");
}

/**
 * This function releases the resources held by the specified objects.
 */
static int
release2(Contig *fragment, Align *best_alignment, Align *vector_start_align, 
    Align *vector_end_align, BtkMessage *message) {

    if (fragment != NULL) {
    	if (contig_release( fragment, message) == ERROR) return ERROR;
    }
    
    if (best_alignment != NULL) {
        if (align_release( best_alignment, message) == ERROR) return ERROR;
    }

    if (vector_start_align != NULL) {
        if (align_release( vector_start_align, message) == ERROR) return ERROR;
    }

    if (vector_end_align != NULL) {
      if (align_release( vector_end_align, message) == ERROR) return ERROR;
    }

    return SUCCESS;
}


/*******************************************************************************
 * Function: is_heterozygote
 *******************************************************************************
 */
static int
is_heterozygote(char b)
{
    if ((b == 'S') || (b == 'K') || (b == 'M') ||
        (b == 'W') || (b == 'R') || (b == 'Y'))
        return 1;

    return 0;
}

/*******************************************************************************
 * Function: is_acgt         
 *******************************************************************************
 */
static int
is_acgt(char b)
{
    if ((b == 'A') || (b == 'C') || (b == 'G') || (b == 'T'))
        return 1;

    return 0;
}

static void
count_stats(Align *best_alignment, int *count_sub, int *count_del, int *count_ins, 
    int *count_tp_hetero, int *count_fp_hetero, int *count_fn_hetero,
    int *count_tp_qv_hetero, int *count_fp_qv_hetero, int *count_fn_qv_hetero,
    int *count_hetero, int *count_mm_hetero, uint8_t *quality_values,
    int *interval, int *min_num_err, Results *results, Options *options)
{
    int i;
    int count_beg = 0, count_end = best_alignment->trace_len;
    int *cumul_error = CALLOC(int, best_alignment->trace_len);

   *count_sub          = *count_del          = *count_ins          = 
   *count_tp_hetero    = *count_fp_hetero    = *count_fn_hetero    =
   *count_tp_qv_hetero = *count_fp_qv_hetero = *count_fn_qv_hetero =
   *count_hetero       = *count_mm_hetero    = 0;

   *min_num_err = best_alignment->trace_len;
    if (*interval > best_alignment->trace_len) 
        *interval = best_alignment->trace_len;

    for (i=count_beg; i<count_end; i++)
    {
        int bpos   = best_alignment->trace_qpos[i],
            qv     = quality_values[bpos],
            qv_min = 20;

        if ((is_acgt(best_alignment->trace_dchar[i]) ||
             is_heterozygote(best_alignment->trace_dchar[i])) &&
            (best_alignment->trace_qchar[i]=='-'))  {
            (*count_del)++;
        }
        else if (best_alignment->trace_dir[i] == '1') {
            (*count_ins)++;
        }
        else if ((is_acgt(best_alignment->trace_qchar[i]) ||
                    is_heterozygote(best_alignment->trace_qchar[i]))
                    &&
                   (is_acgt(best_alignment->trace_dchar[i]) ||
                    is_heterozygote(best_alignment->trace_dchar[i]))
                    &&
                   (best_alignment->trace_qchar[i] !=
                    best_alignment->trace_dchar[i])) {
            (*count_sub)++;
        }

        cumul_error[i] = *count_ins + *count_del + *count_sub;

        if (is_heterozygote(best_alignment->trace_dchar[i]))
            count_hetero++;
        if (is_heterozygote(best_alignment->trace_qchar[i]) &&
            (best_alignment->trace_dchar[i] ==
             best_alignment->trace_qchar[i]))
            count_tp_hetero++;
        if (is_heterozygote(best_alignment->trace_qchar[i]) &&
            (best_alignment->trace_dchar[i] ==
             best_alignment->trace_qchar[i])                &&
            (qv >= qv_min))
            count_tp_qv_hetero++;
        if (is_heterozygote(best_alignment->trace_dchar[i]) &&
            (best_alignment->trace_dchar[i] !=
             best_alignment->trace_qchar[i]))
        {
            count_fn_hetero++;
        }
        if (is_heterozygote(best_alignment->trace_dchar[i]) &&
            (best_alignment->trace_dchar[i] !=
             best_alignment->trace_qchar[i])                &&
            (qv >= 20))
        {
            count_fn_qv_hetero++;
        }
        if (is_heterozygote(best_alignment->trace_qchar[i]) &&
            (best_alignment->trace_dchar[i] !=
             best_alignment->trace_qchar[i])) {
            count_fp_hetero++;
        }
        if (is_heterozygote(best_alignment->trace_qchar[i]) &&
            (best_alignment->trace_dchar[i] !=
             best_alignment->trace_qchar[i])                &&
            (qv >= qv_min)) {
            count_fp_qv_hetero++;
        }
        if (is_heterozygote(best_alignment->trace_dchar[i]) &&
            is_heterozygote(best_alignment->trace_qchar[i]) &&
            (best_alignment->trace_dchar[i] !=
             best_alignment->trace_qchar[i]))
        {
            count_mm_hetero++;
        }
        if (is_heterozygote(best_alignment->trace_qchar[i]) &&
            (best_alignment->trace_dchar[i] ==
             best_alignment->trace_qchar[i])) {
            if (results->min_qv_cor   > qv)
                results->min_qv_cor   = qv;
        }
    }
    for (i=count_beg; i<count_end - (*interval); i++) 
    {
        if (*min_num_err > cumul_error[*interval + i] - cumul_error[i])
            *min_num_err = cumul_error[*interval + i] - cumul_error[i]; 
    }
    FREE(cumul_error);
}

static void
extract_bases_from_alignment(int num_bases, char *bases, Align *best_alignment)
{
    int j, bind;
    for (j=0; j< best_alignment->trace_len; j++) 
    {
        if (best_alignment->trace_dir[j] == '2') {
            continue;
        }

        bind = best_alignment->trace_qpos[j];
        if (is_heterozygote(best_alignment->trace_dchar[j]))
        {
            bases[bind] = best_alignment->trace_dchar[j]; 
            best_alignment->trace_mchar[j] = '|';
        }
    }         
    fprintf(stderr, "\n");
}

static void
update_alignment(int num_bases, char *bases, Align *best_alignment)
{
    int j, bind;
    for (j=0; j< best_alignment->trace_len; j++)
    {
        if (best_alignment->trace_dir[j] == '2') {
            continue;
        }

        bind = best_alignment->trace_qpos[j];
        if (is_heterozygote(bases[bind]) || bases[bind]== 'N')
        {
            best_alignment->trace_qchar[j] = bases[bind];
            if (best_alignment->trace_dchar[j] == bases[bind])
                best_alignment->trace_mchar[j] = '|';
            else 
                best_alignment->trace_mchar[j] = ' ';
        }
    }    
}

static int 
process_sffRead(
    sffHeader    *header, 
    sffRead      *read, 
    Align_params *alp,
    Align_params *alpIUB,
    Contig       *consensus,
    Contig       *consensusrc,
    Vector       *vector,
    double       *params[NUM_PARAMS],
    FILE         *fout,
    BtkMessage   *message,
    Options      *options)
{
    int           num_align; 
    int           alignment_size;
    double       *flows = CALLOC(double, read->number_of_bases);
    uint8_t      *orig_qvs = CALLOC(uint8_t, read->number_of_bases);
    Contig        fragment;
    Range         align_range;
    Range         clear_range;
    Align         best_alignment;
    Align         vector_start_align;
    Align         vector_end_align;

    contig_init(&fragment, message);
    align_init(&best_alignment, message);
    align_init(&vector_start_align, message);
    align_init(&vector_end_align, message);

    if (Btk_compute_tpars_454(header, read, 
        params, flows, orig_qvs, *options, message) != SUCCESS)
    {
        release1_454(read, NUM_PARAMS, params);
        return ERROR;
    }

    if (options->Verbose > 0) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Read: %s\n", read->name);
    }
    if (options->Verbose >= 0) {
        fprintf(fout, "# Read: %s\n", read->name);
    }
   
    if (contig_create(&fragment, read->bases, read->number_of_bases, 
        read->quality_values, message) != SUCCESS)
    {
        release2(&fragment, &best_alignment, &vector_start_align, &vector_end_align,
            message);
        return ERROR;
    }

    if (Btk_compute_match(alp, alpIUB, consensus, consensusrc, &fragment,
            &num_align, &align_range, RepeatFraction, &best_alignment, vector,
            &vector_start_align, &vector_end_align, &clear_range, 0,
            message) == ERROR)
    {
        Count_processing_errors++;
        if (options->Verbose > 0) {
            fprintf(fout, "# Error calling Btk_compute_match: file %s\n",
                read->name);
            fprintf(stderr, "Error calling Btk_compute_match: file %s\n",
                read->name);
        }
        goto error;
    }

    alignment_size=best_alignment.trace_len;

    if (((double)alignment_size >= 
         (double)read->number_of_bases * min_portion_aligned)
                                                                &&
                (read->number_of_bases >= min_read_length))
    {
        if (append_to_train_data(fout, read->number_of_bases, NULL, NUM_PARAMS, 
            params, orig_qvs, flows, NULL, NULL, clear_range, align_range, 
            &vector_start_align, &best_alignment, &vector_end_align,
            0, QVMAX(best_alignment.trace_dpos[0],
            best_alignment.trace_dpos[best_alignment.trace_len-1]),
            1, options, message) == ERROR)
        {
            Count_processing_errors++;
            release2(&fragment, &best_alignment, 
            &vector_start_align, &vector_end_align, message);
            return ERROR;
        }
    }
    release2(&fragment, &best_alignment, &vector_start_align, &vector_end_align,
        message);
    FREE(flows);
    return SUCCESS;
error: 
    release2(&fragment, &best_alignment, &vector_start_align, &vector_end_align,
        message);
    FREE(flows);
    return ERROR; 
}

static int
process_sff_file(Align_params *alp,
    Align_params   *alpIUB,
    char           *sfffile,
    FILE           *fout,
    Contig         *consensus,
    Contig         *consensusrc,
    Vector         *vector,
    BtkMessage     *message,
    Options        *options)
{
    int i, j;
    char *sff_name;
    FILE *sff;
    sffHeader   *h   = CALLOC(sffHeader,   1);
    sffManifest *m   = CALLOC(sffManifest, 1);
    sffRead     *r   = CALLOC(sffRead,     1);
    double      *params[6]= {NULL, NULL, NULL, NULL, NULL, NULL};

    options->dev = dev;

    if (sfffile[0] == '\0') {
        FREE(h);
        FREE(m);
        FREE(r);
        return SUCCESS;
    }
    message->text[0] = '\0';

#ifdef __WIN32
    if ((sff_name = strrchr(sfffile, '\\')) != NULL) {
#else
    if ((sff_name = strrchr(sfffile, '/')) != NULL) {
#endif
        sff_name++;

    }
    else {
        sff_name = sfffile;
    }

    strcpy(options->file_name, sff_name);

    if ((sff = fopen(sfffile, "rb")) == NULL) {
        fprintf(stderr, "Error:  Unable to open SFF file:  %s\n", sfffile);
        goto error;
    }

    message->text[0] = '\0';  /* message may have been set with no error */

    if (readsff_header(sff, h, m) == kWrongFileType)
    {
        FREE(h);
        FREE(m);
        FREE(r);
        fclose(sff);
        return kWrongFileType;
    }

    for (i=0; i < h->number_of_reads; i++) 
    {
        readsff_read(sff, h, r);
        for (j=0; j<NUM_PARAMS; j++)
        {
            params[j] = CALLOC(double, r->number_of_bases);
        }
        if (process_sffRead(h, r, alp, alpIUB, consensus, consensusrc, vector, 
            params, fout, message, options) != SUCCESS)
        {
            fprintf(stderr, "Error processing read %s\n", r->name);
        }
        for (j=0; j<NUM_PARAMS; j++)
        {
            FREE(params[j]);
        }
    }

    //  Read the manifest?
    //
    if (m->manifest_length == 0) {
        //  We haven't read it yet.
        readsff_manifest(sff, h, m);
    }

    if (m->manifest != NULL) {
        //  Make sure that the reads have been rescored.

        if (strstr(m->manifest, "<qualityScoreVersion>1.1.03</qualityScoreVersion>") == NULL) {
            fprintf(stderr, "WARNING:  Fragments not rescored!\n");
        }
    }
    fclose(sff);

    FREE(h->data_block);
    FREE(h);
    FREE(m->manifest);
    FREE(m);
    FREE(r->data_block);
    FREE(r);
    return SUCCESS;

error:
    fclose(sff);
    FREE(h->data_block);
    FREE(h);
    FREE(m->manifest);
    FREE(m);
    FREE(r->data_block);
    FREE(r);
    return ERROR;
}

/*******************************************************************************
 * Function: process_sample_file
 *******************************************************************************
 */
static int 
process_sample_file(Align_params *alp, 
		      Align_params   *alpIUB,
		      char           *frag_name, 
		      FILE           *fout, 
		      Contig         *consensus, 
		      Contig         *consensusrc,
		      Vector         *vector, 
                      BtkLookupTable *table,
                      ContextTable   *ctable,
		      BtkMessage     *message,
                      Options         options, 
                      Results        *results)
{
    int             r,i, count_sub, count_del, count_ins, 
                    count_hetero,    count_mm_hetero,
                    count_tp_hetero, count_tp_qv_hetero, 
                    count_fn_hetero, count_fn_qv_hetero,
                    count_fp_hetero, count_fp_qv_hetero;
    int             num_bases; 
    uint8_t        *quality_values=NULL;
    char           *bases     = NULL, *chemistry = NULL, *call_method = NULL;
    int            *peak_locs = NULL;
    int             num_datapoints;
    int            *chromatogram[NUM_COLORS];
    double         *params[4]= {NULL, NULL, NULL, NULL};
    double         *iheight  = NULL;
    double         *iheight2 = NULL;
    double         *ave_iheight = NULL;
    int             num_align, filetype = -1;
    Contig          fragment;
    Range           align_range;
    Range           clear_range;
    Align           best_alignment;
    Align           vector_start_align;
    Align           vector_end_align;
    int             alignment_size, interval = MIN_LEN_SWEET_SPOT, min_num_err;
    BtkLookupTable *lookup_tbl=NULL;
    ReadInfo        read_info;
    char	   *consensus_seq = NULL;

    if (!ConsensusSpecified) {
        contig_init(consensus, message);
        contig_init(consensusrc, message);
    }
    contig_init(&fragment, message);
    align_init(&best_alignment, message);
    align_init(&vector_start_align, message);
    align_init(&vector_end_align, message);

    if (frag_name[0] == 0) {
    	Count_file_errors++;
	release2(&fragment, &best_alignment, 
	         &vector_start_align, &vector_end_align, message);
        return ERROR;
    }
   
    if (options.Verbose > 0) { 
        fprintf(stderr, "\n");
        fprintf(stderr, "File: %s\n", frag_name);
    }
    if (options.Verbose >= 0) {
        fprintf(fout, "# File: %s\n", frag_name);
    }
    Count_input_files++;

    /* Setting members of Options structure */
    options.nocall = 0;     
    options.recalln = recalln;
    options.recallndb = recallndb;
    options.ladder = ladder;       
    options.process_bases = 1;
    options.xgr = 0; 
    options.dev = dev;
    options.edited_bases = edited_bases;
    options.gauss = gauss;
    options.het = het;
    options.mix = mix;
    options.shift = shift;
    options.renorm = renorm;
    options.respace = respace;
    options.lut_type     = ABI3730pop7;
    options.min_ratio = min_ratio;
    options.time = 0;
    options.tip_dir[0]   = '\0';
    options.tab_dir[0]   = '\0';
    options.tal_dir[0]   = '\0';
    options.raw_data = 0;
    options.xgr = 0;
    options.sf[0] = p1;
    options.sf[1] = p2;
    options.sf[2] = p3;
    options.sf[3] = p4;
    options.inp_phd = inp_phd;
    options.indel_detect = 0;
    options.indel_resolve = 0;
    options.poly         = 0;
    options.poly_dir[0]  = '\0';
    strcpy(options.inp_phd_dir, inp_phd_dir);
    strcpy(options.file_name, frag_name);
    strcpy(options.path, "");


    if ((r = Btk_read_sample_file(frag_name, &num_bases, &bases,
        edited_bases, &peak_locs, &quality_values, &num_datapoints, 
        &chromatogram[0], &chromatogram[1], 
        &chromatogram[2], &chromatogram[3], 
      	&call_method, &chemistry, status_code, &filetype, options,  /* global */
        message)) != SUCCESS)
    {
        if (filetype != ABI && filetype != SCF)
        {
            Btk_release_file_data(bases, peak_locs, quality_values,
                chromatogram, &call_method, &chemistry);
            return kWrongFileType;
        }
        message->text[0] = '\0';  /* message may have been set with no error */
      	Count_file_errors++;
	Btk_release_file_data(bases, peak_locs, quality_values,
                chromatogram, &call_method, &chemistry);
      	return ERROR;
    }

    /* consensus not sepcified externally, and no consensus read from sample
       file. */
    if (!ConsensusSpecified && consensus_seq == NULL) {
        fprintf(stderr, "\nNo consensus sequence specified.\n");
        fprintf(stderr, "Please, specify consensus file \n");
        exit(-1);
    }

    /* consensus read from sample file. */
    if (!ConsensusSpecified) {
	r = contig_create(consensus, consensus_seq,
			  strlen(consensus_seq), quality_values, message);
	FREE(consensus_seq);
    	if (r == ERROR) {
    	    Count_file_errors++;
            Btk_release_file_data(bases, peak_locs, quality_values,
                chromatogram, &call_method, &chemistry);
	    release2(&fragment, &best_alignment,
		     &vector_start_align, &vector_end_align, message);
    	    return ERROR;
    	}

	/* Create reverse complement of consensus. */
    	r = contig_get_reverse_comp(consensusrc, consensus, message);

    	if (r == ERROR) {
    	    Count_file_errors++;
            Btk_release_file_data(bases, peak_locs, quality_values,
                chromatogram, &call_method, &chemistry);
	    release2(&fragment, &best_alignment,
		     &vector_start_align, &vector_end_align, message);
    	    return ERROR;
    	}
    }

    /*  Check that we got some called bases from the input file */
    if (num_bases <= 0) {
        Btk_release_file_data(bases, peak_locs, quality_values,
                chromatogram, &call_method, &chemistry);
    	release1_sanger(bases, peak_locs, chromatogram, quality_values, 
            NUM_PARAMS, NULL, NULL);
	release2(&fragment, &best_alignment, 
	         &vector_start_align, &vector_end_align, message);
    	sprintf(message->text,"No called base information in input file");
    	Count_file_errors++;
    	return ERROR;
    }
    if (num_datapoints <= 0) {
        Btk_release_file_data(bases, peak_locs, quality_values,
                chromatogram, &call_method, &chemistry);
    	release1_sanger(bases, peak_locs, chromatogram, quality_values, 
            NUM_PARAMS, NULL, NULL);
	release2(&fragment, &best_alignment, 
	         &vector_start_align, &vector_end_align, message);
    	sprintf(message->text,"No peak information in input file");
    	Count_file_errors++;
    	return ERROR;
    }

    if (options.inp_phd > 0)
        fprintf(stderr,
            "Reading original bases and locations from phd file in dir: %s\n",
            options.inp_phd_dir);

    /* Removing the suffix ".Z" or ".gz" from the name of compressed sample */
    if (frag_name[strlen(frag_name) - 3] == '.' &&
        frag_name[strlen(frag_name) - 2] == 'g' &&
        frag_name[strlen(frag_name) - 1] == 'z') {
        frag_name[strlen(frag_name) - 3] = '\0';
    }
    if (frag_name[strlen(frag_name) - 2] == '.' &&
        frag_name[strlen(frag_name) - 1] == 'Z') {
        frag_name[strlen(frag_name) - 2] = '\0';
    }

    strcpy(options.file_name, frag_name);
    strcpy(options.path, "");      
 
    /*
     *  If external lookup table not specified, determine from PDMF
     *  string which internal table to use (POP5, POP6, etc.).
     */
    if ((options.Verbose > 0) && (chemistry != NULL)) {
        (void)fprintf(stderr, "Mobility file = %s\n", chemistry);
    }
    lookup_tbl = table;
    if (table == NULL) {
        if ((lut == ABI3730pop7) || (!lut && (chemistry != NULL) &&
             strstr(chemistry, "POP5")))
        {
            lookup_tbl = Btk_get_3730pop7_table();
            if ((options.Verbose > 1) && (options.het || options.mix)) {
                (void)fprintf(stderr, "Using built-in ABI3730 Pop-7 table.\n");
            }
        }
        else if ((lut == ABI3700pop5) || (!lut && (chemistry != NULL) && 
             strstr(chemistry, "POP5"))) 
        {
            lookup_tbl = Btk_get_3700pop5_table();
            if ((options.Verbose > 1) && (options.het || options.mix)) {
                (void)fprintf(stderr, "Using built-in ABI3700 Pop-5 table.\n");
            }
        }
        else if ((lut == ABI3700pop5) ||
                (!lut &&  (chemistry != NULL) &&
                 strstr(chemistry, "3700") && 
                 strstr(chemistry, "POP6")))
        {
            lookup_tbl = Btk_get_3700pop6_table();
            if ((options.Verbose > 0) && (options.het || options.mix)) {
                (void)fprintf(stderr, "Using built-in ABI3700 Pop-6 table.\n");
            }
        }
        else if ((lut == ABI3100) ||
                (!lut &&  (chemistry != NULL) &&
                 strstr(chemistry, "3100"))) {
            lookup_tbl = Btk_get_3100pop6_table();
            if ((options.Verbose > 0) && (options.het || options.mix)) {
                fprintf(stderr,
                             "Using built-in 3100 Pop-6 table for this run.\n");
            }
        }
        else {
            lookup_tbl = Btk_get_3700pop5_table();
            if ((options.het || options.mix) && options.Verbose > 1) {
                (void)fprintf(stderr,
                "Can't select the lookup table automatically. \n");
                (void)fprintf(stderr,
                "Using built-in ABI 3700 Pop-5 table.\n");
            }
        }
    }
    table = lookup_tbl;

    if (quality_values == NULL) {
        quality_values = CALLOC(uint8_t, num_bases);
        MEM_ERROR(quality_values);
    }

    /* In default mode, determine new bases from electropherogram
     * In -mix or -het mode, will use original bases for now
     */
    if (!options.het && !options.mix)
    {
        if (Btk_compute_tpars_Sanger(&num_bases, &bases, &peak_locs, 
                &num_datapoints, chromatogram, "ACGT", NUM_PARAMS, 
                &params[0], &params[1], &params[2], &params[3], 
                &iheight, &iheight2, &ave_iheight, &read_info, table, ctable, 
                options, message, results ) != SUCCESS)
        {
            Count_processing_errors++;
            release1_sanger(bases, peak_locs, chromatogram, quality_values, 
                NUM_PARAMS, params, iheight);
            FREE(call_method);
            FREE(iheight2);
            FREE(ave_iheight);
            release2(&fragment, &best_alignment, 
                  &vector_start_align, &vector_end_align, message);
            return ERROR;
        }
    }
    
    r = contig_create(&fragment, bases, num_bases, quality_values, message);
    if (r == ERROR) {
    	Count_processing_errors++;
        release1_sanger(bases, peak_locs, chromatogram, quality_values, 
            NUM_PARAMS, params, iheight);
        FREE(call_method);
        FREE(iheight2);
        FREE(ave_iheight);
        release2(&fragment, &best_alignment, &vector_start_align, &vector_end_align, 
            message);
     return ERROR;
    }
 
    if (Btk_compute_match(alp, alpIUB, consensus, consensusrc, &fragment,
	    &num_align, &align_range, RepeatFraction, &best_alignment, vector, 
            &vector_start_align, &vector_end_align, &clear_range, 0,
            message) == ERROR) 
    {
    	Count_processing_errors++;
        if (options.Verbose > 0) {
            fprintf(fout, "# Error calling Btk_compute_match: file %s\n",
                frag_name);
            fprintf(stderr, "Error calling Btk_compute_match: file %s\n",
                frag_name);
        }
	goto error;
    }

    /* Recall all the bases except mixed ones 
     * Don't change the number of bases
     */
    if (options.het || options.mix)
    {
        extract_bases_from_alignment(num_bases, bases, &best_alignment);
        if (Btk_compute_tpars_Sanger(&num_bases, &bases, &peak_locs, 
                &num_datapoints, chromatogram, "ACGT", NUM_PARAMS,
                &params[0], &params[1], &params[2], &params[3],
                &iheight, &iheight2, &ave_iheight, &read_info, table, ctable, 
                options, message, results ) != SUCCESS)
        {
            Count_processing_errors++;
            release1_sanger(bases, peak_locs, chromatogram, quality_values, 
                NUM_PARAMS, params, iheight);
            FREE(call_method);
            FREE(iheight2);
            FREE(ave_iheight);
            release2(&fragment, &best_alignment,
                &vector_start_align, &vector_end_align, message);
            return ERROR;
        }
        update_alignment(num_bases, bases, &best_alignment);
    }

    if (num_align == 0) {
        if (options.Verbose > 0) {
    	    fprintf(fout, "# No good alignments - ignore file %s\n",
	        frag_name);
    	    fprintf(stderr, "No good alignments - ignore file %s\n",
	        frag_name);
        }
    	Count_no_alignments++;
    } else if ( (num_align != 1) || (best_alignment.score == -1) ) {
        if (options.Verbose > 0) { 
    	    fprintf(fout, "# Possible repeat !!! - ignore file %s\n",
	        frag_name);
    	    fprintf(stderr, "Possible repeat !!! - ignore file %s\n",
	        frag_name);
        }
    	Count_possible_repeats++;
    } 
    else 
    {
        
    	/* Count number of mismatches */
        count_stats(&best_alignment, &count_sub, &count_del, &count_ins, 
            &count_tp_hetero, &count_fp_hetero, &count_fn_hetero,
            &count_tp_qv_hetero, &count_fp_qv_hetero, &count_fn_qv_hetero,
            &count_hetero, &count_mm_hetero, quality_values, 
            &interval, &min_num_err, results, &options);

        alignment_size=best_alignment.trace_len;

        if (
#if CALIBRATION
    	/* Append the results to train file only if the fraction of
         * mismatches is < MaxFractionOfErrors       
         */
    	    ((double)min_num_err < (double)interval * MaxFractionOfErrors)
#else 
        // use a reference sequence instead of consensus sequence
            ((double)(count_sub+count_ins+count_del) <
                 (double)alignment_size * MaxFractionOfErrors)
#endif
                                                                &&
                 ((double)alignment_size >=
                  (double)num_bases * min_portion_aligned)
                                                                &&
                 (num_bases >= min_read_length))
        {
       	    if (append_to_train_data(fout, num_bases, peak_locs, NUM_PARAMS, 
                params, quality_values, iheight, iheight2, ave_iheight, 
                clear_range, align_range, 
                &vector_start_align, &best_alignment, &vector_end_align, 0, 
                QVMAX(best_alignment.trace_dpos[0],
                best_alignment.trace_dpos[best_alignment.trace_len-1]),
                options.Verbose, &options, message) == ERROR) 
            {
                Count_processing_errors++;
                release1_sanger(bases, peak_locs, chromatogram, quality_values, 
                    NUM_PARAMS, params, iheight);
                FREE(call_method);
                release2(&fragment, &best_alignment, &vector_start_align, 
                    &vector_end_align, message);
                return ERROR;
      	    }
            results->align_length += alignment_size;
            results->count_del += count_del;
            results->count_ins += count_ins;
            results->count_sub += count_sub;
            results->count_err += count_del + count_ins + count_sub;
            results->count_hetero       += count_hetero;
            results->count_tp_hetero    += count_tp_hetero;
            results->count_fn_hetero    += count_fn_hetero;
            results->count_fp_hetero    += count_fp_hetero;
            results->count_mm_hetero    += count_mm_hetero;
            results->count_fn_qv_hetero += count_fn_qv_hetero;
            results->count_fp_qv_hetero += count_fp_qv_hetero;
            results->count_tp_qv_hetero += count_tp_qv_hetero;   
    	} else {
            if (options.Verbose > 0) {
                fprintf(stderr, 
                    "Alignment_size=%d\n", alignment_size);
                fprintf(stderr, 
                    "tot_err=%d MaxFractionOfErrors=%f\n",
                    count_sub+count_ins+count_del, MaxFractionOfErrors);

                fprintf(fout, "# BAD PROCESSING!!! - ignore file %s\n",
	            frag_name);
      	        fprintf(stderr, "BAD PROCESSING!!! - ignore file %s\n",
	            frag_name);
            }
            Count_processing_errors++;
        }
        if (options.Verbose > 0)
    	    fprintf(stderr, "Align range length:%4d\n", alignment_size);
        if (options.Verbose > 0)
    	fprintf(stderr, "Number of matched bases: %d\n", 
	        alignment_size-count_sub-count_ins-count_del);
    }

    /* consensus was read from the sample file */
    if (!ConsensusSpecified) {
	if (consensus != NULL) {
	    contig_release(consensus, message);
	}
	if (consensusrc != NULL) {
	    contig_release(consensusrc, message);
	}
    }
    release1_sanger(bases, peak_locs, chromatogram, quality_values,
        NUM_PARAMS, params, iheight);
    FREE(call_method);
    FREE(iheight2);
    FREE(ave_iheight);
    release2(&fragment, &best_alignment, &vector_start_align, &vector_end_align, 
        message);
 
    return SUCCESS;
error:
    for (i = 0; i < NUM_PARAMS; i++) {
        FREE(params[i]);
    }
    /* consensus was read from the sample file */
    if (!ConsensusSpecified) {
	if (consensus != NULL) {
	    contig_release(consensus, message);
	}
	if (consensusrc != NULL) {
	    contig_release(consensusrc, message);
	}
    }

    release1_sanger(bases, peak_locs, chromatogram, quality_values,
        NUM_PARAMS, params, iheight);
    FREE(call_method);
    FREE(iheight2);
    FREE(ave_iheight);
    release2(&fragment, &best_alignment,
         &vector_start_align, &vector_end_align, message);
    return ERROR;
}

static int
process_file(Align_params *alp,
    Align_params   *alpIUB,
    char           *frag_name,
    FILE           *fout,
    Contig         *consensus,
    Contig         *consensusrc,
    Vector         *vector,
    BtkLookupTable *table,
    ContextTable   *ctable,
    BtkMessage     *message,
    Options        *options,
    Results        *results)
{
    int   result;

    NUM_PARAMS = 6;
    if (process_sff_file(alp, alpIUB, frag_name, fout,
        consensus, consensusrc, vector, message,
        options) != kWrongFileType)
    {
           fprintf(stderr, "This is SFF file\n");
    }
    else
    {
        NUM_PARAMS = 4; 
        if (process_sample_file(alp, alpIUB, frag_name, fout,
            consensus, consensusrc, vector, table, ctable, message,
            *options, results) != kWrongFileType)
        {
               fprintf(stderr, "This is ABI or SCF file\n");
        }
        else {
           fprintf(stderr, "Wrong file type: not an ABI, SCF or SFF\n");
        }
    }

    return result;
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

/*
 * This function processes all files contained in the specified directory.
 * Its synopsis is:
 *
 * result = process_dir(align_pars, align_pars_IUB, fout, cons, consensusrc, 
 *                      dir, vector, message, results)
 *
 * where
 *      align_pars      is the address of an Align_params data structure for
 *                      aligning consensus with sample
 *      align_pars_IUB  is the address of an Align_params data structure for
 *                      aligning vector with sample
 *      fout            is the (already open) output stream
 *      cons            is the address of the consensus Contig
 *      consensusrc     is the address of the reverse complement Contig
 *      dir             is the name (path) of a directory
 *      vector          is the address of the vector, or NULL if there is none
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      results         is the address of Results structure where information
 *                      about processing resuts is accumulated
 */
static int process_dir(
               Align_params   *ap,
	       Align_params   *apIUB,
	       FILE           *fout,
	       Contig         *consensus,
	       Contig         *consensusrc,
	       char           *dir,
	       Vector         *vector,
               BtkLookupTable *table,
               ContextTable   *ctable,
	       BtkMessage     *message,
               Options         options,
               Results        *results)
{
    DIR *d;
    struct dirent *de;
    char path[MAXPATHLEN];
    struct stat statbuf;
    int r;

    message->text[0] = 0;

    if ((d = opendir(dir)) == NULL) {
    	error(dir, "couldn't open dir", errno);
    	return ERROR;
    }

    while ((de = readdir(d)) != NULL) {
    	(void)sprintf(path, "%s/%s", dir, de->d_name);
    	if (stat(path, &statbuf) != 0) {
      	    error(path, "can't stat", errno);
      	    closedir(d);
      	    return ERROR;
    	}
    	if (statbuf.st_mode & S_IFDIR) {
      	    /* skip subdirectories */
            continue;
    	}

        if ((r = process_file(ap, apIUB, path, fout, consensus,
            consensusrc, vector, table, ctable, message, &options,
            results)) != SUCCESS) 
        {
      	    fprintf(stderr,"Error processing input file: %s\n",message->text);
      	    fprintf(fout,"# Error processing input file: %s\n",message->text);
    	}
    }

    (void)closedir(d);
    return SUCCESS;
}

/*
 * This function processes all files contained in the specified project file.
 * Its synopsis is:
 *
 * result = process_projectfile(align_pars, align_pars_IUB, fout, projectFile,
	     vector, message)
 *
 * where
 *      align_pars      is the address of an Align_params data structure for
 *                      aligning consensus with sample
 *      align_pars_IUB  is the address of an Align_params data structure for
 *                      aligning vector with sample
 *      fout            is the (already open) output stream
 *      projectFile     is the name (path) of a project file
 *      vector          is the address of the vector, or NULL if there is none
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 */
static int
process_projectfile(
            Align_params   *ap,
	    Align_params   *apIUB,
	    FILE           *fout,
	    char           *projectFile,
	    Vector         *vector,
            BtkLookupTable *table, 
            ContextTable   *ctable,
	    BtkMessage     *message,
            Options         options,
            Results        *results)
{
    char    buffer1[BUFLEN];
    char    buffer2[BUFLEN];
    char    *current, *previous, *tmp;
    Contig  consensus;
    Contig  consensusrc;
    char    fragmentName[MAXPATHLEN];
    int     r, have_consensus;
    FILE    *infile;

    if ((infile=fopen(projectFile,"r"))== NULL) {
    	sprintf(message->text, "Unable to open project file '%s'\n",
	        projectFile);
    	return ERROR;
    }

    have_consensus=0;
    buffer2[0]='\0';
    previous = buffer2;
    current  = buffer1;
    while ((fgets(current, BUFLEN, infile)) != NULL) {
    	/* Ignore all white space lines and comments */
    	if (strspn(current, " \t\r\n") == strlen(current)) {
      	    continue;
    	}
    	if ((current[0] == '#')
		|| ((current[0] == '/') && (current[1] == '*'))
		|| (current[0] == ';')) {
	    continue;
      	}

    	sscanf(current, "%s %s", ConsensusName, fragmentName);

    	/* If consensus is different from previous, read it in. */
    	if ( strcmp(current, previous) != 0 ) {
      	    /* Release previous, if necessary. */
            if ( have_consensus ) {
	    	r=contig_release( &consensus, message );
		if ( r==ERROR ) { return ERROR; }
		r=contig_release( &consensusrc, message );
		if ( r==ERROR ) { return ERROR; }
		have_consensus=0;
      	    }

            /* Read consensus file. */
            r=local_read_fasta( ConsensusName, &consensus, message );
      	    if ( r==ERROR ) { return ERROR; }

      	    /* Create reverse complement of consensus. */
      	    r=contig_get_reverse_comp( &consensusrc, &consensus, message );
      	    if ( r==ERROR ) { return ERROR; }
      	    have_consensus=1;
    	}

    	tmp=previous;
    	previous=current;
    	current=tmp;

        if (options.Verbose >0) {
      	    fprintf(stderr, "Reference = %s\n", ConsensusName );
    	    fprintf(fout, "# Reference = %s\n", ConsensusName );
        }

    	if (process_file(ap, apIUB, fragmentName, fout, &consensus,
	    &consensusrc, vector, table, ctable, message, &options,
            results) == ERROR) 
        {
      	    fprintf(stderr,"Error processing input file: %s\n",fragmentName);
            fprintf(fout,"# Error processing input file: %s\n",fragmentName);
    	}
    }

    fclose(infile);

    /* Clean up. */
    if ( have_consensus ) {
    	r=contig_release( &consensus, message );
    	if ( r==ERROR ) { return ERROR; }
    	r=contig_release( &consensusrc, message );
    	if ( r==ERROR ) { return ERROR; }
    }

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
    results->count_chunks = 0;
    for ( i=0; i<MAX_CHUNKS; i++ ) {
        statsReset( &results->stats[i] );
        /* assumption: normalized spacing is mostly in [0,2] */
        histoReset( &results->histo[i], 0.0, 2.0 );
    }
}

static void
print_results_stats( const Results *results, FILE *fp )
{
    int i;
    for( i=0; i<results->count_chunks; i++ ) {
        const Stats *s = &results->stats[i];
        const Histo *h = &results->histo[i];
        if ( statsCount(s) < MIN_PEAKS_PER_CHUNK ) break;
        if ( histoCount(h) < MIN_PEAKS_PER_CHUNK ) break;
        fprintf( fp, "%2d %2d %5.3f %5.3f  %g\n",
                 i, i, statsStdDev(s), histoPercentileStdDev(h),
                 statsCount(s) );
    }    
}

/*******************************************************************************
 * Function: main
 *******************************************************************************
 */
int main(int argc, char** argv)
{
    char           *args, *lookup_table, *context_table;
    int             i, j, optind;
    FILE           *fout;
    Align_params    align_pars, align_pars_IUB;
    BtkMessage      message;
    BtkLookupTable *table  = NULL;
    ContextTable   *ctable = NULL;
    Contig          consensus, consensusrc;
    Options         options;
    Results         results;
    Vector          vector;
    Vector         *vecP;
 
    if (argc == 1) {
    	usage(argv);
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
    InputType           = INP_FILES;
    RepeatFraction      = 0.85;
    MaxFractionOfErrors = 0.05;
    Match               = 10;
    MisMatch            = -20;
    GapInit             = -40;
    GapExt              = -40;
    nocall              = 0;
    recalln             = 0;
    recallndb           = 0;
    ladder              = 0;
    edited_bases        = 0;
    gauss               = 1;
    shift               = 0;
    respace             = 0;
    renorm              = 0;
    het                 = 0;
    mix                 = 0;
    min_ratio           = 0.1;
    p1                  = 0.3;
    p2                  = 0.5;
    lookup_table        = NULL;
    context_table       = NULL;
    options.Verbose     = 1;
    dev                 = 0;
    options.process_bases = 1;
    inp_phd             = 0;
    inp_phd_dir[0]      ='\0';
    min_portion_aligned = 0.;
    min_read_length     = 0;
    initialize_results(&results);

    /*
     * This argument-processing code probably looks a bit weird; if so, it's
     * because it was converted from getopt(), which isn't supported under
     * the Win32-compatible version of gcc we have.
     */
    for (optind = 1; optind < argc; optind++) 
    {
        if (argv[optind][0] != '-') {
            break;
        }

        /* Make sure there is an option argument for those options which are
         * supposed to have an argument
         */
        if (((strcmp(argv[optind], "-min_ratio") == 0) ||
             (strcmp(argv[optind], "-ct")        == 0) ||
             (strcmp(argv[optind], "-C" )        == 0) ||
             (strcmp(argv[optind], "-V" )        == 0) ||
             (strcmp(argv[optind], "-P" )        == 0) ||
             (strcmp(argv[optind], "-S" )        == 0) ||
             (strcmp(argv[optind], "-M" )        == 0) ||
             (strcmp(argv[optind], "-X" )        == 0) ||
             (strcmp(argv[optind], "-G" )        == 0) ||
             (strcmp(argv[optind], "-fe")        == 0) ||
             (strcmp(argv[optind], "-fr")        == 0) ||
             (strcmp(argv[optind], "-o" )        == 0) ||
             (strcmp(argv[optind], "-p" )        == 0) ||
             (strcmp(argv[optind], "-p1")        == 0) ||
             (strcmp(argv[optind], "-p2")        == 0) ||
             (strcmp(argv[optind], "-p3")        == 0) ||
             (strcmp(argv[optind], "-p4")        == 0) ||
             (strcmp(argv[optind], "-d" )        == 0) ||
             (strcmp(argv[optind], "-ipd" )      == 0)
            )
            &&
            (optind==argc-1 || argv[optind+1][0]=='-'))
        {
            usage(argv);
            exit(2);
        }

        /* Make sure that all the flags are from alowed list
         */
        if (((optind == argc-1) || (argv[optind + 1][0] == '-')) &&
            ((strcmp(argv[optind], "-h")            != 0) &&
             (strcmp(argv[optind], "-dev")          != 0) &&
             (strcmp(argv[optind], "-Q")            != 0) &&
             (strcmp(argv[optind], "-V")            != 0) &&
             (strcmp(argv[optind], "-QQ")           != 0) &&
             (strcmp(argv[optind], "-VV")           != 0) &&
             (strcmp(argv[optind], "-nocall")       != 0) &&
             (strcmp(argv[optind], "-recalln")      != 0) &&
             (strcmp(argv[optind], "-recallndb")    != 0) &&
             (strcmp(argv[optind], "-ladder")       != 0) &&
             (strcmp(argv[optind], "-edited_bases") != 0) &&
             (strcmp(argv[optind], "-raw")          != 0) &&
             (strcmp(argv[optind], "-mc")           != 0) &&
             (strcmp(argv[optind], "-xgr")          != 0) &&
             (strcmp(argv[optind], "-shift")        != 0) &&
             (strcmp(argv[optind], "-convolved")    != 0) &&
             (strcmp(argv[optind], "-het")          != 0) &&
             (strcmp(argv[optind], "-mix")          != 0)))
        {
            usage(argv);
            fprintf(stderr, "\nInvalid flag specified.\n");
            exit(2);
        }

        /* Parse the current argument string (args), starting from the 2nd
         * character (j==1), which is not '-'.
         */
        for (j = 1, args = argv[optind]; args[j] != '\0'; j++) {
            i = args[j];
            switch(i) {
           
            case 'c':
                if ((strcmp(args, "-ct") != 0) &&
                    (strcmp(args, "-convolved") != 0)) {
                    if (dev) usage_dev(argv);
                    else     usage(argv);
                    exit(2);
                }
                else if (strcmp(args, "-convolved") == 0) {
                    gauss = 0;
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-ct") != 0) {
                    context_table = argv[++optind];
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }

            case 'C':
                ConsensusSpecified++;
                (void)strncpy(ConsensusName, argv[++optind], 
                    sizeof(ConsensusName));
                break;

            case 'd':
                if (strcmp(args, "-dev") == 0) {
                    if (argc == 2)
                    {
                        usage_dev(argv);
                        exit(2);
                    }
                    else
                    {
                       dev++;
                       j = strlen(args) - 1;   /* break out of inner loop */
                       break;
                    }
                }
                else if (strcmp(args, "-d") == 0) {
                    InputType =  INP_DIR;
                    (void)strncpy(InputName, argv[++optind], sizeof(InputName));
                }
                break;

            case 'G':
                GapExt=-atoi(argv[++optind]);
                GapInit = GapExt;
                break;
 
            case 'e':
                if (strcmp(args, "-edited_bases") == 0) {
                    edited_bases++;
                    j = strlen(args)-1;  /* break out of inner loop */
                    break;
                }

            case 'f':
                if ((strcmp(args, "-fr") != 0) &&
                    (strcmp(args, "-fe") != 0)) {
                    if (dev) usage_dev(argv);
                    else     usage(argv);
                    exit(2);
                }
                else if (strcmp(args, "-fr") == 0) {
                    RepeatFraction=atof(argv[++optind]); 
                    j = strlen(args) - 1;  /* break out of inner loop */  
                    break;
                }
                else if (strcmp(args, "-fe") == 0) {
                    MaxFractionOfErrors=atof(argv[++optind]);
                    j = strlen(args) - 1;  /* break out of inner loop */  
                    break;
                }

            case 'F':  /* Version 1.0 had a -F flag which is now the default */
                break; /* Catch this case here so it's not flagged as an error*/

            case 'h':
                if (strcmp(args, "-het") == 0) {
                    het++;
                    MaxFractionOfErrors = 1.;
                    RepeatFraction = 0.9999;
                    j = strlen(args)-1;  /* break out of inner loop */
                    break;
                }
                help_message(argc, argv);
                exit(2);

            case 'i':
                if (strcmp(args, "-ipd") == 0) {
                    inp_phd = 1;
                    (void)strncpy(inp_phd_dir, argv[++optind],
                                  sizeof(inp_phd_dir));
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else {
                    if (dev) usage_dev(argv);
                    else     usage(argv);
                    exit(2);
                }

            case 'm':
                if ((strcmp(args, "-min_ratio") != 0) &&
                    (strcmp(args, "-mix")       != 0)) {
                    if (dev) usage_dev(argv);
                    else     usage(argv);
                    exit(2);
                } 
                else if (strcmp(args, "-min_ratio") == 0) {
                    min_ratio = atof(argv[++optind]);
                    j = strlen(args) - 1;  /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-mix") == 0) {
                    mix++;
                    Match = 20;
                    MaxFractionOfErrors = 1.;
                    RepeatFraction = 0.9999;
                    j = strlen(args)-1;  /* break out of inner loop */
                    break;
                }
           
            case 'M':
                Match=atoi(argv[++optind]);
                break; 

            case 'n':
                if (strcmp(args, "-nocall"  ) != 0) {   
                    if (dev) usage_dev(argv);
                    else     usage(argv);
                    exit(2);
                }
                else if (strcmp(args, "-nocall" ) == 0)
                    nocall++;
                j = strlen(args) - 1;   /* break out of inner loop */
                break;

            case 'o':
                OutputSpecified++;
                (void)strncpy(OutputName, argv[++optind], sizeof(OutputName));
                break;

            case 'p':
                if ((strcmp(args, "-p1") != 0) &&
                    (strcmp(args, "-p2") != 0) &&
                    (strcmp(args, "-p3") != 0) &&
                    (strcmp(args, "-p4") != 0) &&
                    (strcmp(args, "-p"   ) != 0)) {
                    if (dev) usage_dev(argv);
                    else     usage(argv);
                    exit(2);
                }
                else if (strcmp(args, "-p1") == 0) {
                    p1 = atof(argv[++optind]); 
                    j = strlen(args) - 1;   /* break out of inner loop */
                }
                else if (strcmp(args, "-p2") == 0) {
                    p2 = atof(argv[++optind]);                  
                    j = strlen(args) - 1;   /* break out of inner loop */
                }
                else if (strcmp(args, "-p3") ==0) {
                    p3 = atof(argv[++optind]);                  
                    j = strlen(args) - 1;   /* break out of inner loop */
                }
                else if (strcmp(args, "-p4") == 0) {
                    p4 = atof(argv[++optind]);                  
                    j = strlen(args) - 1;   /* break out of inner loop */
                }
                else if (strcmp(args, "-p") == 0) {
                   InputType =  INP_PROJECTFILE;
                   ConsensusSpecified++;
                   (void)strncpy(InputName, argv[++optind], sizeof(InputName));
                } 
                break;

            case 'P':
                PrimerSpecified++;
                (void)strncpy(PrimerName, argv[++optind], sizeof(PrimerName));
                break;

            case 'a':
                min_portion_aligned=atof(argv[++optind]);
                break;
            case 'l':
                if (strcmp(args, "-ladder") == 0)
                {
                    ladder++;
                }
                else if (strcmp(args, "-l") == 0) {
                    min_read_length=atoi(argv[++optind]);
                    break;
                }
                else {
                    if (dev) usage_dev(argv);
                    else     usage(argv);
                    exit(2);
                }
                j = strlen(args) - 1;  /* break out of inner loop */
                break;

            case 'r':
                if ((strcmp(args, "-recalln") != 0) &&
                    (strcmp(args, "-recallndb") != 0) &&
                    (strcmp(args, "-renorm" ) != 0) &&
                    (strcmp(args, "-respace") != 0)) {
                    if (dev) usage_dev(argv);
                    else     usage(argv);
                    exit(2);
                }
                else if (strcmp(args, "-recalln") == 0) {
                    recalln++;
                    j = strlen(args)-1;  /* break out of inner loop */
                }
                else if (strcmp(args, "-recallndb") == 0) {
                    recallndb++;
                    j = strlen(args)-1;  /* break out of inner loop */
                }
                else if (strcmp(args, "-renorm") == 0) {
                    renorm++;
                    j = strlen(args)-1;  /* break out of inner loop */
                }
                else if (strcmp(args, "-respace") == 0) {
                    respace++;
                    j = strlen(args)-1;  /* break out of inner loop */
                }
                break;

             case 'Q':
                options.Verbose--;
                break;

            case 's':
                if (strcmp(args, "-shift") != 0) {
		    if (dev) usage_dev(argv);
                    else     usage(argv);
                    exit(2); 
                } 
                else if (strcmp(args, "-shift") == 0) 
                { 
                    shift++; j = strlen(args)-1;  /* break out of inner loop */ 
                } 
                break;

	    case 'S': 
                SiteSpecified++; 
                (void)strncpy(SiteName, argv[++optind], sizeof(SiteName)); break;

	    case 't': 
                if (strcmp(args, "-t") != 0) { 
                    if (dev) usage_dev(argv);
                    else     usage(argv);
                    exit(2);
                }
                lookup_table = argv[++optind];
                break;

            case 'V':
                VectorSpecified++;
                (void)strncpy(VectorName, argv[++optind], sizeof(VectorName));
                break;
            
            case 'X':
                MisMatch=-atoi(argv[++optind]);
                break;

            case '?':
            default:
                usage(argv);
                exit(2);
            }
        }
    }

    if ((nocall && recalln) || (nocall && recallndb) || (recalln && recallndb) ||
        (nocall && ladder)  || (recalln && ladder)   || (recallndb && ladder))
    {
        fprintf(stderr,
        "\ntrain: Please, specify only one of the options -nocall, -recalln");
        fprintf(stderr,
        "or -recallndb\n");
         exit(-1);
    }

    if ((RepeatFraction < 0 ) ||
        (RepeatFraction > 1.0) ) {
        fprintf(stderr,
        "Error: RepeatFraction must be in the range [0,1.0]\n");
        exit(2);
    }

    if ((MaxFractionOfErrors < 0 ) || (MaxFractionOfErrors > 1.0) ) {
        fprintf(stderr,
        "Error: MaxFractionOfErrors must be in the range [0,1.0]\n" );
        exit(2);
    }
    
    /*
     * Set line buffering on the status output so that someone monitoring
     * progress can actually see incremental progress.
     */
    (void)setbuf(stderr, NULL);

    /* Read the lookup table */
    if (lookup_table != NULL) 
    {
        if ((table = Btk_read_lookup_table(lookup_table)) == NULL) {
            (void)fprintf(stderr,
                "Couldn't read the lookup table '%s'\n", lookup_table);
            exit(1);
        }
        else {
            if (options.Verbose > 1) {
                (void)fprintf(stderr,
                "Using lookup table '%s' for this run.\n",
                lookup_table);
            }
        }
    }

    /* Read the context table */
    if (context_table != NULL) {
        if ((ctable = read_context_table(context_table)) == NULL) {
            (void)fprintf(stderr,
                "Couldn't read the context table '%s'\n", context_table);
            exit(1);
        }
        else {
            if (options.Verbose > 1) {
                (void)fprintf(stderr,
                "Using context table '%s' of dimension %d for this run.\n",
                context_table, ctable->dimension);
            }
        }
    }

    if (OutputSpecified) {
    	/* Set up output file. */
    	if ((fout=fopen(OutputName, "w"))== NULL) {
      	    fprintf(stderr, "Cannot open output file '%s'\n", OutputName);
      	    exit(ERROR);
    	}
    	fprintf(stderr, "Software Version: %s\n", TT_VERSION);
    } else {
    	if (options.Verbose > 1) 
            fprintf(stderr, "No output file specified.  Using stdout.\n");
    	fout=stdout;
    }
  
    /* Set parameters for alignment with consensus */
    if (set_alignment_parameters( &align_pars, Match, MisMatch, GapInit,
	GapExt  , &message ) ==ERROR ) {
    	fprintf(stderr, message.text);
    	exit (ERROR);
    }

    /* Set parameters for alignment with vector */
    if (set_alignment_parameters_IUB( &align_pars_IUB, 
        Match * VECTOR_MATCH_MULTIPLIER, MisMatch, GapInit, GapExt, 
        &message ) ==ERROR ) {
    	fprintf(stderr, message.text);
    	exit (ERROR);
    }

    /* Print parameter values to output file and stderr. */
    fprintf( fout, "# Version   = %s\n", TT_VERSION);
    fprintf( fout, "# Match     = %d\n", Match );
    fprintf( fout, "# MisMatch  = %d\n", MisMatch );
    fprintf( fout, "# GapInit   = %d\n", GapInit);
    fprintf( fout, "# GapExt    = %d\n", GapExt);
    fprintf( fout, "# RepeatFraction = %f\n", RepeatFraction );
    fprintf( fout, "#\n");

    fprintf(stderr, " Version   = %s\n", TT_VERSION);
    fprintf(stderr, " Match     = %d\n", Match);
    fprintf(stderr, " MisMatch  = %d\n", MisMatch);
    fprintf(stderr, " GapInit   = %d\n", GapInit);
    fprintf(stderr, " GapExt    = %d\n", GapExt);
    fprintf(stderr, " RepeatFraction = %f\n", RepeatFraction );
    fprintf(stderr, "\n");

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

    switch(InputType) {
    case  INP_FILES:

    	if (ConsensusSpecified) {
      	    fprintf( fout, "# Reference = %s\n", ConsensusName );
    	    fprintf( fout, "#\n");
            
            fprintf(stderr, " Reference = %s\n", ConsensusName );
            fprintf(stderr, "\n");

       	    /* Read consensus file. */
    	    if (local_read_fasta( ConsensusName, &consensus, &message)
                == ERROR) {
      	        fprintf(stderr, message.text);
   	        exit (ERROR);
            }
   	    /* Create reverse complement of consensus. */
    	    if (contig_get_reverse_comp( &consensusrc, &consensus, &message)
                == ERROR) {
  	        fprintf(stderr, message.text);
 	        exit (ERROR);
            }
    	}

            if (optind >= argc)
            {
                fprintf(stderr, "No input was given to train. Please add to the command line a list\n");
                fprintf(stderr, "of sample file names, a sample file directory, or a project file.\n");
                exit(ERROR);
            }

    	/* Loop for each fragment: */
    	for (i = optind; i < argc; i++) {
   	    if (process_file(&align_pars, &align_pars_IUB, argv[i],
                fout, &consensus, &consensusrc, vecP, table, ctable, &message, 
                &options, &results) != SUCCESS) 
            {
                fprintf(stderr,"Error processing input file: %s\n",
		message.text);
	        fprintf(fout,"# Error processing input file: %s\n",
	    	message.text);
    	    }
        }
    	break;
        
        case  INP_DIR:

    	if (ConsensusSpecified) {
    	    fprintf( fout, "# Reference = %s\n", ConsensusName );
        	fprintf( fout, "#\n");

    	    /* Read consensus file. */
    	    if (local_read_fasta( ConsensusName, &consensus, &message) 
                == ERROR) {
  	        fprintf(stderr, "Error reading consensus file: %s\n",
                        message.text);
  	        exit (ERROR);
    	    }

    	    /* Create reverse complement of consensus. */
    	    if (contig_get_reverse_comp( &consensusrc, 
	    			    &consensus, &message) == ERROR) {
  	        fprintf(stderr, "Error getting reverse compliment: %s\n",
		    message.text);
  	        exit(ERROR);
            }
    	}

    	for (i = optind; i < argc; i++) {
            if (process_dir(&align_pars, &align_pars_IUB, fout, &consensus,
                &consensusrc, argv[i], vecP, table, ctable, &message,
                options, &results)   != SUCCESS) {
                fprintf(stderr,"Error processing input files: %s\n", message.text);
                exit (ERROR);
            }
        }
    	break;

    case  INP_PROJECTFILE:
    	if (process_projectfile(&align_pars, &align_pars_IUB, fout,
			InputName, vecP, table, ctable, &message,
                        options, &results) != SUCCESS) {
  	     fprintf(stderr, "%s: %s\n", argv[0], message.text);
 	     exit (ERROR);
    	}
    	break;
	
    default:
    	fprintf(stderr, "%s: internal error: name type %d\n", argv[0],
		  InputType);
    	exit(ERROR);
    }

    /* Clean up. */
 
    if (alignment_parameters_release( &align_pars, &message) == ERROR) {
    	fprintf(stderr, "Error releasing memory: %s\n",message.text);
    	exit(ERROR);
    }

    if (alignment_parameters_release( &align_pars_IUB, &message) == ERROR) {
    	fprintf(stderr, "Error releasing memory: %s\n",message.text);
    	exit(ERROR);
    }

    if ((InputType ==  INP_FILES ) || ( InputType ==  INP_DIR ) ) {
    	if(contig_release( &consensus, &message)) {
      	    fprintf(stderr, "Error releasing memory: %s\n",message.text);
  	    exit (ERROR);
    	}

    	if (contig_release( &consensusrc, &message) == ERROR) {
  	    fprintf(stderr, "Error releasing memory: %s\n",message.text);
  	    exit (ERROR);
    	}
    }

    /*
     *   Print out file count/processing summary.
     */
    fprintf(fout,"#\n");
    fprintf(fout,"# File count summary for this run: \n");
    fprintf(fout,"#  %6d files input\n", Count_input_files);
    fprintf(fout,"#  %6d files processed\n", Count_input_files -
	    (Count_file_errors + Count_processing_errors + 
	    Count_possible_repeats + Count_no_alignments));
    fprintf(fout,"#  %6d files rejected with file errors\n", 
	    Count_file_errors);
    fprintf(fout,"#  %6d files rejected with processing errors\n",
	    Count_processing_errors);
    fprintf(fout,"#  %6d files rejected with possible repeats\n",
	    Count_possible_repeats);
    fprintf(fout,"#  %6d files rejected with no alignments\n",
	    Count_no_alignments);
    fclose(fout);

    /*
     *   Print out Results summary
     */
    fprintf(stderr, "\nSUM OF ALIGN LENGTHS:%4d\n", 
         results.align_length);
    fprintf(stderr, "TOTAL MATCHED BASES:  %4d\n",
        results.align_length-
        results.count_sub-
        results.count_ins-
        results.count_del);
    fprintf(stderr, "NUMBERS OF ERRORS: \n");
    fprintf(stderr,
        "   SUB=%d,   INS=%d,   DEL=%d   TOT=%d\n",
         results.count_sub, results.count_ins, results.count_del,
        (results.count_sub +results.count_ins +results.count_del));
    printf("\n# NUMBERS OF ERRORS:   SUB=%d,   INS=%d,   DEL=%d   TOT=%d\n",
         results.count_sub, results.count_ins, results.count_del,
        (results.count_sub +results.count_ins +results.count_del));
    fprintf(stderr, "FRACTIONS OF ERRORS: \n");
    fprintf(stderr,
        "   SUB=%f,   INS=%f,   DEL=%f   TOT=%f\n",
        (double)results.count_sub/QVMAX((double)results.align_length, 1),
        (double)results.count_ins/QVMAX((double)results.align_length, 1),
        (double)results.count_del/QVMAX((double)results.align_length, 1),
        (double)(results.count_sub+results.count_ins+results.count_del)
        /QVMAX((double)results.align_length, 1));
    printf("\n# FRACTIONS OF ERRORS: SUB=%f,   INS=%f,   DEL=%f   TOT=%f\n",
        (double)results.count_sub/QVMAX((double)results.align_length, 1),
        (double)results.count_ins/QVMAX((double)results.align_length, 1),
        (double)results.count_del/QVMAX((double)results.align_length, 1),
        (double)(results.count_sub+results.count_ins+results.count_del)
        /QVMAX((double)results.align_length, 1));

    if (results.count_hetero > 0) 
    {
        int qv_min =20;
        fprintf(stderr, "\nSCALING FACTORS=%f %f %f %f QV_MIN=%d\n",
            options.sf[0], options.sf[1], options.sf[2], options.sf[3],
            qv_min);

        fprintf(stderr, 
          "\nCORRECT      HETEROZYGOTES: %3d\n",
            results.count_hetero);

        fprintf(stderr, 
            "\nTP HETEROZYGOTES: %3d\n",
            results.count_tp_hetero);

        fprintf(stderr, 
            "FN HETEROZYGOTES: %3d FN_RATIO=%f",
            results.count_fn_hetero,
            (double) results.count_fn_hetero/
            (double)(results.count_fn_hetero + results.count_tp_hetero +
                     results.count_mm_hetero));
        fprintf(stderr, "   AVE FN PER ONE PROCCESSED FILE: %4.2f\n",
            (double)results.count_fn_hetero/
            (double)(Count_input_files -
            (Count_file_errors + Count_processing_errors +
            Count_possible_repeats + Count_no_alignments)));

        fprintf(stderr,
            "FP HETEROZYGOTES: %3d FP_RATIO=%f",
            results.count_fp_hetero,
            (double) results.count_fp_hetero/
            (double)(results.count_fp_hetero + results.count_tp_hetero +
                     results.count_mm_hetero));
        fprintf(stderr, "   AVE FP PER ONE PROCCESSED FILE: %4.2f\n",
            (double)results.count_fp_hetero/
            (double)(Count_input_files -
            (Count_file_errors + Count_processing_errors +
            Count_possible_repeats + Count_no_alignments)));    
        fprintf(stderr, "MM HETEROZYGOTES: %3d\n",
            results.count_mm_hetero);
        fprintf(stderr, 
            "ERRORS IN HETEROZYGOTES: %3d ",
            results.count_fn_hetero+results.count_fp_hetero); 
        fprintf(stderr, "     AVE PER ONE PROCCESSED FILE: %4.2f\n",
            (double)(results.count_fp_hetero+results.count_fn_hetero)/
            (double)(Count_input_files -
            (Count_file_errors + Count_processing_errors +
            Count_possible_repeats + Count_no_alignments)));

        fprintf(stderr,
            "FRACTION OF FN HETEROZYGOTES: %f\n",
            (double)results.count_fn_hetero/
            (double)results.count_hetero);

        fprintf(stderr, "\n");
        fprintf(stderr, 
            "TP HETEROZYGOTES                      OF QV >= %2d: %3d\n",
            qv_min, results.count_tp_qv_hetero);
        fprintf(stderr, 
            "FN HETEROZYGOTES SUBSTITUTED BY BASES OF QV >= 20: %3d\n",
            results.count_fn_qv_hetero);
        fprintf(stderr, 
            "FP HETEROZYGOTES                      OF QV >= %2d: %3d\n",
            qv_min, results.count_fp_qv_hetero);
        fprintf(stderr, 
            "ERRORS IN HETEROZYGOTES               OF QV >= %2d: %3d\n",
            qv_min, results.count_fn_qv_hetero+
                    results.count_fp_qv_hetero);
        fprintf(stderr,
            "FRAC  OF HETEROZYGOTES SUBST BY BASES OF QV >= %d: %f\n\n",
            qv_min, (double)results.count_fn_qv_hetero/
            (double)results.count_hetero);
        fprintf(stderr,
            "N_TP-N_FP-N_FN = %d   FRAC1=%f\n\n",
            results.count_tp_hetero - results.count_fn_hetero -
                results.count_fp_hetero/FP_TO_FN_RATIO,
            (double)(results.count_tp_hetero - results.count_fn_hetero -
                     results.count_fp_hetero/FP_TO_FN_RATIO)/
            (double)(results.count_hetero));
        fprintf(stderr,
            "N_TP_QV20-N_FP_QV20-N_FN_QV20 = %d   FRAC2=%f\n\n",
            results.count_tp_qv_hetero - results.count_fn_qv_hetero -
            results.count_fp_qv_hetero,
            (double)(results.count_tp_qv_hetero - results.count_fn_qv_hetero -
                     results.count_fp_qv_hetero) /
            (double)(results.count_hetero));
    }

    {
        FILE *fp = fopen( "stats.dat", "w" );
        if (fp) {
            print_results_stats( &results, fp );
        } else {
            fprintf( stderr, "error opening stats.dat\n" );
        }
        fclose(fp);
    }
    exit(SUCCESS);
}
