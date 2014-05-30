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
 * $Id: main.c,v 1.23 2009/01/10 21:40:16 gdenisov Exp $
 */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#ifndef __WIN32
#include <unistd.h>
#endif
#include <math.h>
#ifdef __WIN32
#define MAXPATHLEN      (255)
#else
#include <sys/param.h>
#endif
#include <time.h>
#ifdef __DEVSTUDIO
#include <io.h>
#else
#include <dirent.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <float.h>

#include "Btk_qv.h"
#include "util.h"
#include "Btk_qv_data.h"
#include "Btk_lookup_table.h"
#include "context_table.h"
#include "train.h"
#include "Btk_compute_qv.h"
#include "Btk_match_data.h"
#include "Btk_qv_io.h"
#include "Btk_default_table.h"
#include "Btk_process_raw_data.h"
#include "SFF_Toolkit.h"
#include "ABI_Toolkit.h"

#define MAXBIN 1024     /* Max number of bins for quality report */
#define CHECK_LICENSE 0
#define MAX_NAME_LEN 1000
#define MAX_LENGTH_FLOWGRAM 10000
#define SUP(a) (((a)>0) ? (1) : (0))
typedef enum {
    TEXT_MODE = 0,FASTA_MODE, QUAL_MODE, FLOWGRAM_MODE, MANIFEST_MODE, 
    ACCNO_MODE, PHD_MODE 
} MODE;
int numManifest;

static int Verbose;	/* Whether and how much status info to print */

static char InputName[BUFLEN];	/* path of dir or file-of-files */
static int InputType;

static int ConsensusSpecified;		/* whether consensus is specified */
static char ConsensusName[BUFLEN];	/* path of consensus file */

static int AlnType;
static int HprType;

static int OutputFasta;		/* whether to write FASTA files */
static char FastaDirName[BUFLEN];	/* path of dir */
static int FastaType;

static int OutputFlow;

static int OutputFourMultiFastaFiles;
static char MultiFastaFilesDirName[BUFLEN];
static char multiseqsFileName[BUFLEN];
static char multiqualFileName[BUFLEN];
static char multilocsFileName[BUFLEN];
static char multistatFileName[BUFLEN];
static char status_code[BUFLEN];

static int dev = 0;
static int opts = 0;
static int file_type = -1;

static int OutputSCF;		/* whether to write Staden SCF files */
static char SCFDirName[BUFLEN];
static int SCFType;

static int OutputPhd;		/* whether to write .phd.1 files */
static char PhdDirName[BUFLEN];	/* path of dir */
static int PhdType;

static int OutputQual;		/* whether to write .qual files */
static char QualDirName[BUFLEN];	/* path of dir */
static int QualType;

static int OutputQualRpt;	/* whether to write quality report */
static char QualRptName[BUFLEN];/* filename for quality report */
static struct Qual_data_struct {
    int number_over_20[MAXBIN];   /* Number of bases in read with QV >= 20 */
    int trimmed_lengths[MAXBIN];  /* Read lengths after trimming           */
    int sum_number;               /* Totals to compute averages            */
} *Qual_data = NULL;

static char multiseqFileName[BUFLEN];
static char multiqualFileName[BUFLEN];
static int trim_window = 10;      /* width of trimming window */
static float trim_threshold = 20; /* when average of trim window goes above
                                   * this, trimming stops 
                                   */
static int left_trim_point, right_trim_point;

/* Any sequence whose score >= RepeatFraction*HighScore is considered a repeat*/
static double RepeatFraction;

/* Scores and penalties. */
static int Match;
static int MisMatch;
static int Insertion;
static int Deletion;

clock_t start_clock, curr_clock;

static void
usage(int argc, char *argv[])
{
    fprintf( stderr, 
    "\nVersion: %s\n"
    "usage: %s\n"
    "    [ -h ]     [ -Q ]      [ -V ]\n"
    "    [ -nocall] [ -recalln ][ -edited_bases ] \n" 
    "    [ -het   ] [ -mix     ][ -min_ratio <phr>     ]\n"  
    "    [ -trim_window  <size>][ -trim_threshold <qv> ][ -ipd <dir>]\n"
    "    [ -t <lookup_table>   ][ -C <consensus_file>  ][ -cv3   ]\n"
    "    [ -indel_detect ][ -indel_resolve ][ -indloc <loc> ][ -indsize <size> ]\n"
    "    [ -3730 ][ -3700pop5 ][ -3700pop6 ][ -3100 ][ -mbace]\n"   
    "    [ -p   | -pd   <dir> ][ -s | -sd <dir> ][ -q   | -qd   <dir> ]\n"  
    "    [ -d   | -dd   <dir> ][ -c | -cd <dir> ][ -tab | -tabd <dir> ]\n" 
    "    [ -hpr | -hprd <dir> ][ -qr     <file> ][ -tal | -tald <dir> ]\n"
    "    [ -sa         <file> ][ -qa     <file> ][ -f   ][ -o   <dir> ]\n"
    "    { <sample_file(s)>    | -id     <dir>  | -if  <fileoffiles> }\n"
             , TT_VERSION, argv[0] );
}

static void
usage_dev(int argc, char *argv[])
{
    fprintf( stderr,
    "\nVersion: %s\n"
    "usage: %s\n"
    "    [ -h ] [ -Q ] [ -V ] [ -opts ] [ -dev] \n"
    "    [ -nocall][ -recalln   ][ -recallndb ][ -edited_bases ][ -ladder] \n"  
    "    [ -het   ][ -mix       ][ -min_ratio <phr>     ]\n"
    "    [ -trim_window  <size> ][ -trim_threshold <qv> ]\n" 
    "    [ -t <lookup_table>    ][ -ct <context_table>  ]\n"
    "    [ -cv3   ] [ -time     ][ -C <consensus_file>  ]\n"
    "    [ -convolved ][ -shift ][ -renorm ][ -respace ]\n"
    "    [ -raw ] [ -xgr ][ -mc ] \n"
    "    [ -indel_detect ][ -indel_resolve ][ -indloc <loc> ][ -indsize <size> ]\n"
    "    [ -3730][ -3700pop5][ -3700pop6][ -3100][ -mbace]\n"
    "    [ -p  | -pd  <dir> ] [ -s | -sd <dir> ] [ -tip | -tipd <dir> ]\n"
    "    [ -q  | -qd  <dir> ] [ -c | -cd <dir> ] [ -tal | -tald <dir> ]\n"
    "    [ -d  | -dd  <dir> ] [ -qr     <file> ] [ -tab | -tabd <dir> ]\n"
    "    [ -ipd <dir> ]   [ -hpr | -hprd <dir> ] [ -f          <file> ]\n"
    "    [ -sa       <file> ] [ -qa     <file> ] [ -o           <dir> ]\n"
    "    { <sample_file(s)>   | -id     <dir>    | -if  <fileoffiles> }\n"
             , TT_VERSION, argv[0] );
}

static void
help_message(int argc, char *argv[])
{
    fprintf( stderr, 
"    -h                   (Help) This message\n"
"    -Q                   (Quiet) Turn off status messages\n"
"    -V                   (Verbose) Output more status messages\n"
"    -nocall              Disable base recalling and just use the original\n"
"                         called bases read from the input sample file\n"
"    -recalln             Disable adding bases to or deleting from the\n"
"                         original called sequence. Only recall Ns\n"
"    -het                 (For Sanger data only) Call hetezygotes. The ratio\n"
"                         of two alleles is assumed to be 1:1\n"
"    -mix                 (For Sanger data only) Call mixed bases. No assumption\n"  
"                         is made about the alleles ratio\n" 
"    -min_ratio <ratio>   (For Sanger data only) Override the default threshold\n"
"                         ratio of heights of\n"
"    -trim_window <size>  Set the trimming window size for averaging quality\n"
"                         to the specified value. The default is 10.\n"
"    -trim_threshold <qv> Set the average quality value used in trimming to\n"
"    -C <consensusfile>   Specify the name of the FASTA file which contains\n"
"                         the consensus sequence\n"
"    -edited_bases        (For Sanger data only) Start base recalling from \n"
"                         ABI's edited bases\n"
"    -t <table>           Use specified lookup table. This option overrides\n"
"                         the default (automatic choice of the lookup table)\n"
"                         as well as the options -3700pop5, -3700pop6, -3100,\n"
"                         and -mbace. To get a message showing \n"
"                         which table was used, specify -V option\n"
"    -indel_detect        (For Sanger data only) Detect heterozygous indels\n"
"                         and report their location, size and string to stderr\n"
"    -indel_resolve       (For Sanger data only) Detect heterozygous indels\n"
"                         and deconvolve the mixed electropherogram into two\n"
"                         pure electropherograms corresponding to the short and \n"
"                         long haplotypes\n"
"    -indloc <loc>        (For Sanger data only) Specify the location in electropherogram\n"
"                         where heterozygous indel is supposed to be detected\n"
"    -indsize <size>      (For Sanger data only) Specify the size of indel\n"
"                         which is supposed  \n"
"    -3730                (For Sanger data only) Use the built-in ABI 3730-pop7 lookup table\n"
"    -3700pop5            (For Sanger data only) Use the built-in ABI 3700-pop5 lookup table\n"
"    -3700pop6            (For Sanger data only) Use the built-in ABI 3700-pop6 lookup table\n"
"    -3100                (For Sanger data only) Use the built-in ABI 3100-pop6 lookup table\n"
"    -mbace               (For Sanger data only) Use the built-in MegaBACE lookup table\n"
"    -c                   (For Sanger data only) Output SCF file(s) in the \n"
"                         current directory\n"
"    -cd <dir>            (For Sanger data only) Output SCF file(s) in the \n"
"                         specified directory\n"
"    -cv3                 (For Sanger data only) Use version 3 for output \n"
"                         SCF file. Default is version 2.\n"
"    -o <dir>             (For Sanger data only) Output multi-fasta files of \n"
"                         bases (tt.seq), their locations (tt.pos), quality \n"
"                         values (tt.qual) and status reports (tt.status) to \n"
"                         directory <dir>\n"
"    -p                   Output .phd.1 file(s), in the current directory,\n"
"                         one file per sample\n"
"    -pd <dir>            Output .phd.1 file(s), in the specified directory,\n"
"                         one file per sample\n"
"    -q                   For Sanger data, output .qual file(s) in the current \n"
"                         directory; for 454 data, output multi-FASTA file of\n"
"                         quality values to stdout\n"
"    -qa <file>           Append .qual file(s) to <file>\n"
"    -qd <dir>            Output .qual file(s), in the specified directory\n"
"    -s                   For Sanger data, output .seq file(s) in FASTA format \n"
"                         in the current directory; for 454 data, output \n"
"                         multi-FASTA sequence file to stdout\n"
"    -sa <file>           Append .seq file(s) in multi-FASTA format to <file>\n"
"    -sd <dir>            Output .seq file(s) in FASTA format in the specified\n"
"                         directory\n"
"    -qr <file>           Output a quality report that gives data for a\n"
"                         histogram on the number of reads with quality\n"
"                         values >= 20, to the specified file\n"
"    -if <file>           Read the input sample filenames from the specified\n"
"                         file\n"
"    -id <dir>            Read the input sample files from specified directory\n"
"    -tab                 (For Sanger data only) Call heterozygotes or mixed \n"
"                         bases and output .tab file(s) in the  current directory\n"
"    -tabd <dir>          (For Sanger data only) Call mixed bases and output \n"
"                         .tab file(s), in the specified directory\n"
"    -tal                 (For Sanger data only) Output .tal file(s) in the \n"
"                         current directory\n"
"    -tald <dir>          (For Sanger data only) Output .tal file(s) in the \n"
"                         specified directory <dir> \n"
"    -hpr                 Output a homopolymer runs file in current directory\n"
"    -hprd <dir>          Output a homopolymer runs file(s),in the specified directory\n"
"    -d                   (For Sanger data only) Output .poly file(s) in the \n"
"                         current directory\n"
"    -dd  <dir>           (For Sanger data only) Output .poly file(s),in the \n"
"                         specified directory <dir>\n"
"    -ipd <dir>           (For Sanger data only) Input the original bases and \n"
"                         peak locations from PHD file in the specified directory.\n"
        );
}


static void
output_options(Options options)
{
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "\tchemistry = %s\n", options.chemistry);
    fprintf(stderr, "\tedited_bases = %d\n", options.edited_bases);
    fprintf(stderr, "\tfile_name = %s\n", options.file_name);                                        
    fprintf(stderr, "\tgauss = %d\n", options.gauss); 
    fprintf(stderr, "\thet = %d\n",options.het );  
    fprintf(stderr, "\tindel_detect = %d\n", options.indel_detect);
    fprintf(stderr, "\tindel_resolve = %d\n", options.indel_resolve);
    fprintf(stderr, "\tindsize = %d\n", options.indsize);
    fprintf(stderr, "\tindloc = %d\n", options.indloc);
    fprintf(stderr, "\tinp_phd = %d\n", options.inp_phd);
    fprintf(stderr, "\tinp_phd_dir = %s\n", options.inp_phd_dir);     
    fprintf(stderr, "\tlut_type = %d\n", options.lut_type);
    fprintf(stderr, "\tmin_ratio = %f\n", options.min_ratio);
    fprintf(stderr, "\tmix = %d\n", options.mix );    
    fprintf(stderr, "\tmulticomp = %d\n", options.multicomp); 
    fprintf(stderr, "\tnocall = %d\n", options.nocall); 
    fprintf(stderr, "\trespace = %d\n", options.respace);
    fprintf(stderr, "\tscf_dir = %s\n", options.scf_dir); 
    fprintf(stderr, "\tscf_version = %d\n", options.scf_version);
    fprintf(stderr, "\tpath = %s\n", options.path);   
    fprintf(stderr, "\tpoly_dir = %s\n", options.poly_dir);
    fprintf(stderr, "\tprocess_bases = %d\n", options.process_bases);
    fprintf(stderr, "\traw_data = %d\n", options.raw_data);
    fprintf(stderr, "\trecalln = %d\n", options.recalln);
    fprintf(stderr, "\trecallndb = %d\n", options.recallndb); 
    fprintf(stderr, "\trecallndb = %d\n", options.ladder);         
    fprintf(stderr, "\trenorm = %d\n", options.renorm); 
    fprintf(stderr, "\tshift = %d\n", options.shift);  
    fprintf(stderr, "\ttab_dir = %s\n", options.tab_dir);  
    fprintf(stderr, "\ttal_dir = %s\n", options.tal_dir);
    fprintf(stderr, "\thpr_dir = %s\n", options.hpr_dir);
    fprintf(stderr, "\ttime = %d\n", options.time); 
    fprintf(stderr, "\ttip_dir = %s\n", options.tip_dir);
    fprintf(stderr, "\tVerbose = %d\n", options.Verbose);
    fprintf(stderr, "\txgr = %d\n", options.xgr);   

    return;
}

/*
 *  accum_qual_report()
 *
 *  Calculate the number of quality values over 20, longer continuous read
 *  over 20, and last base position over 20 and tally the appropriate bins
 *  in the quality report structure.  Allocate the structure the first
 *  time through.
 */

void accum_qual_report(struct Qual_data_struct **qual_data, 
                       uint8_t *quality_values, int num_called_bases,
                       int trimmed_read_length) 
{
    static int mem_error = 0;
    int i, number_over_20 = 0;

    if (mem_error) return;  /* Only try alloc memory and report problem once */

    /* Allocate memory for the quality data structure. */

    if (*qual_data == NULL) {
        *qual_data = CALLOC(struct Qual_data_struct, 1);
        if (*qual_data == NULL) {
            mem_error = 1;
            fprintf(stderr,"Unable to allocate memory for quality report.\n");
            return;
        }
    }

    for (i = 0; i < num_called_bases; i++)
        if (quality_values[i] >= 20)
            number_over_20++;

    if (number_over_20 / 10 < MAXBIN)
    {
        (*qual_data)->number_over_20[number_over_20 / 10]++;
        (*qual_data)->trimmed_lengths[trimmed_read_length / 10]++;
    }

    (*qual_data)->sum_number += number_over_20;
}


/*
 *  output_qual_report()
 *
 *  After all the sample files have been processed, this routine will
 *  output the data accumulated in for a quality value histogram.
 */

void output_qual_report(struct Qual_data_struct *qual_data, 
			char *filename)
{
  FILE *f;
  int i, last_bucket=0, num_of_files = 0;

  if (qual_data == NULL)
    return;

  /* Figure out last bucket to report */
  for (i=0; i < MAXBIN-1; i++) {
    if (qual_data->number_over_20[i] > 0 || qual_data->trimmed_lengths[i] > 0)
      last_bucket = i+1;
    num_of_files += qual_data->number_over_20[i];
  }
  /* If there's no data, don't generate a report */
  if (num_of_files == 0) return;

  if ((f = fopen(filename, "w")) == NULL) {
    fprintf(stderr,"Unable to open file '%s': %s\n",
	    filename, strerror(errno));
    return;
  }

  fprintf(f,"TraceTuner Quality Value Histogram Data\n");
  fprintf(f,"Number of Files Reported = %d\n\n",num_of_files);
  fprintf(f,"            Num of Reads     Trimmed Read Lengths\n");
  fprintf(f,"Bucket      with QV >= 20    window=%d, threshold=%d\n",
		  trim_window, INT_FLT(trim_threshold));
  fprintf(f,"---------   --------------   -----------------------\n");

  for (i=0; i < last_bucket; i++) {

    fprintf(f,"%4d-%-4d         %8d                  %8d\n",
	    i*10, (i+1)*10-1,
	    qual_data->number_over_20[i], qual_data->trimmed_lengths[i]);
  }
  fprintf(f,"\nAverage           %8d\n",
	  qual_data->sum_number / num_of_files);
  fclose(f);
}


/*
 * This function prints as formatted an error message.  Its synopsis is:
 *
 * error(label, msg, err)
 *
 * where
 *	label	is a char* to be printed at the beginning of the line.
 *		Usually the name of the file being processed at the time of
 *		the error, or the name of the program running (argv[0]).
 *	msg	is a char* to a message to print.  Usually describes what
 *		was attempted.
 *	err	is an int error code.  Usually errno.  If 0, then no error
 *		code is printed, nor is any system-provided textual
 *		description of the error.
 */
static void
error(char *label, char *msg, int err)
{
    char *sysmsg;

    if ((sysmsg = strerror(err)) != NULL) {
	(void)fprintf(stderr, "%s: %s: %d: %s\n", label, msg, err, sysmsg);
    }
    else if (err == 0) {
	(void)fprintf(stderr, "%s: %s\n", label, msg);
    }
    else {
	(void)fprintf(stderr, "%s: %s: %d\n", label, msg, err);
    }
}


/*
 * This function parses the specified mobility file name and finds the
 * appropriate built-in lookup table.  Its synopsis is:
 *
 * result = parse_mobility_file_name(query, options, message)
 *
 * where
 *	query		is the address of the mobility file name
 *	options		is the address of the Options
 *	message		is the address of a BtkMessage where information about
 *			an error will be put, if any
 *	result		the found lookup table
 */
static BtkLookupTable* 
parse_mobility_file_name(char* query, Options* options, BtkMessage* message,
    BtkLookupTable* table) 
{
    BtkLookupTable* lookup_table = table;
    if (query == NULL)
        return Btk_get_3730pop7_table();
    if (strstr(query, "Mega") || strstr(query, "BACE")) {
        /* the sample file was generated from MegaBACE or LI-COR sequencing
           machine. */
        options->lut_type = MegaBACE;    
        if (table == NULL) {
            lookup_table = Btk_get_mbace_table();
            (void)fprintf(stderr,
                "Using a default built-in MegaBACE table.\n");
        }
        return lookup_table;
    } else if (strstr(query, "LI-COR")) {
        /* the sample file was generated from MegaBACE or LI-COR sequencing
           machine. */
        options->lut_type = ABI3730pop7;
        if (table == NULL) {
            lookup_table = Btk_get_3730pop7_table();
            (void)fprintf(stderr,
                "Can't select the lookup table automatically.\n");
            (void)fprintf(stderr,
                "Using a default built-in ABI3730 Pop-7 table.\n");
        }
        return lookup_table;
    } else if (strstr(query, "POP5")) {
        options->lut_type = ABI3700pop5;
        if (table == NULL) {
            lookup_table = Btk_get_3700pop5_table();
            if (options->Verbose > 1) {
                (void)fprintf(stderr,
                    "Using a built-in ABI3700 Pop-5 table.\n");
            }
        }
        return lookup_table;
    } else if (strstr(query, "3700") && strstr(query,"POP6")) {
        options->lut_type = ABI3700pop6;
        if (table == NULL) {
            lookup_table = Btk_get_3700pop6_table();
            if (options->Verbose > 1) {
                (void)fprintf(stderr,
                    "Using a built-in ABI3700 Pop-6 table.\n");
            }
        }
        return lookup_table;
    } else if (strstr(query, "3100") && strstr(query, "POP6")) {
        options->lut_type = ABI3100;
        if (table == NULL) {
            lookup_table = Btk_get_3100pop6_table();
            if (options->Verbose > 1)
                (void)fprintf(stderr,
                "Using a built-in ABI3100 Pop-6 table.\n");
        }
        return lookup_table;
    }
    else if (strstr(query, "3730") || strstr(query, "KB") ||
        strstr(query, "POP7"))
    {
        options->lut_type = ABI3730pop7;
        if (table == NULL) {
            lookup_table = Btk_get_3730pop7_table();
            if (options->Verbose > 1)
                (void)fprintf(stderr,
                "Using a built-in ABI3730 Pop-7 table.\n");
        }
        return lookup_table;
    } else if (strstr(query, "ET")) {
        if (strstr(query, "ET Terminators") || strstr(query, "ET Primers")) {
            /* The mobility file name indicates the sample file was
               generated from MolDyn_MegaBACE machine. */
            options->lut_type = ABI3730pop7;
            if (table == NULL) {
                lookup_table = Btk_get_3700pop5_table();
                (void)fprintf(stderr,
                    "Can't select the lookup table automatically.\n");
                (void)fprintf(stderr,
                    "Using a default built-in ABI3730 Pop-7 table.\n");
            }
            return lookup_table;
        }
    } else {
        options->lut_type = ABI3730pop7;
        if (table == NULL) {
            lookup_table = Btk_get_3700pop5_table();
            (void)fprintf(stderr,
                "Can't select the lookup table automatically.\n");
            (void)fprintf(stderr,
                "Using a default built-in ABI3730 Pop-7 table.\n");
        }
        return lookup_table;
    }
    return lookup_table;
}

static int
process_sample_file(BtkLookupTable *table, ContextTable *ctable,
    char *path, char *ConsensusName, char *ConsensusSeq,
    Options *options, BtkMessage *message)
{
    char *seq_name, *called_bases, *call_method;
    int   r, j, num_called_bases=0, *called_peak_locs, num_datapoints;
    int   trimmed_read_length, file_type = -1;
    int  *chromatogram[NUM_COLORS]; 
    uint8_t *quality_values;
    int   consFromSample = 0; // whether consensus sequence is from
                              // the sample file
    Results results;

    if (path[0] == '\0') {
        return SUCCESS;
    }
    message->text[0] = '\0';

#ifdef __WIN32
    if ((seq_name = strrchr(path, '\\')) != NULL) {
#else
    if ((seq_name = strrchr(path, '/')) != NULL) {
#endif
        seq_name++;
    }
    else {
        seq_name = path;
    }

    strcpy(options->file_name, seq_name);

    for (j = 0; j < NUM_COLORS; j++) {
        chromatogram[j] = NULL;
    }
    called_bases     = NULL;
    called_peak_locs = NULL;
    call_method      = NULL;
    quality_values   = NULL;

    if ((!ConsensusSpecified) || ConsensusSeq == NULL) {
        consFromSample = 1;
    }

    if (Btk_read_sample_file(path, &num_called_bases, &called_bases,
        options->edited_bases, &called_peak_locs, &quality_values,
        &num_datapoints, &chromatogram[0], &chromatogram[1],
        &chromatogram[2], &chromatogram[3], &call_method, &(options->chemistry),
        status_code, &file_type, *options, message) != SUCCESS)
    {
        if (file_type != ABI && file_type != SCF)
        {
            return kWrongFileType;
        }
        else if (status_code[0] == '\0')
        {
            strcpy(status_code, "ABIFILE_FAILURE");
        }
        else
        {
            /* status_code == "PHREDFILE_FAILURE */
            ;
        }
        output_four_multi_fasta_files(multiseqsFileName, multiqualFileName,
            multilocsFileName, multistatFileName, num_called_bases,
            called_bases, quality_values, called_peak_locs,
            results.frac_QV20_with_shoulders, status_code, *options);
        status_code[0] = '\0';
        goto error;
    }
    message->text[0] = '\0';  /* message may have been set with no error */

    if (opts)
        output_options(*options);
#if 0
    fprintf(stderr, "After Btk_read_sample_file: options->file_name=%s\n", options->file_name);
#endif

    /*
     *  Check that we read some valid bases and peak locations
     */
    if ((num_called_bases <= 0) && !options->raw_data) {
        sprintf(message->text,"Can't process - no base calls in file.");
        goto error;
    }
    if (num_datapoints <= 0) {
        sprintf(message->text,"Can't process - no peak locations in file.");
        goto error;
    }

    if (SUP(options->het)+SUP(options->mix) > 1)
    {
        fprintf(stderr, "\nPlease, specify no more than one of the options");
        fprintf(stderr, " -het and -mix\n");
        exit_message(options, -1);
    }
    else if (SUP(options->nocall)+SUP(options->recalln) > 1)
    {
        fprintf(stderr, "\nPlease, specify no more than one of the options");
        fprintf(stderr, " -nocall and -recalln\n");
        exit_message(options, -1);
    }
    else if (
        (SUP(options->nocall)+SUP(options->recalln)+SUP(options->recallndb)+
         SUP(options->ladder) > 1))
    {
        fprintf(stderr, "\nPlease, specify no more than one of the options");
        fprintf(stderr, " -nocall, -recalln, -recallndb and -ladder\n");
        exit_message(options, -1);
    }
    /* "Secret" option */
    else if (SUP(options->het)+SUP(options->recalln)-SUP(options->poly) > 1)
    {
        fprintf(stderr, "\nPlease, specify no more than one of the options");
        fprintf(stderr, " -het and -recalln\n");
        exit_message(options, -1);
    }
    else if ((SUP(options->het)+SUP(options->mix)==0) &&
             (options->tab_dir[0] != '\0'))
    {
        fprintf(stderr,
        "\nOption -tab or -tabd <dir> can only be used with ");
        fprintf(stderr, "either -het or -mix\n");
        exit_message(options, -1);
    }

    if (options->Verbose <= 1) {
        fprintf(stderr, "%s\n", options->file_name);
    }
    else
    {
        fprintf(stderr, "%s: %d points, %d bases originally\n",
            seq_name, num_datapoints, num_called_bases);
    }

    /* Removing the suffix ".Z" or ".gz" from the name of compressed sample */
    if (path[strlen(path) - 3] == '.' &&
        path[strlen(path) - 2] == 'g' &&
        path[strlen(path) - 1] == 'z') {
        path[strlen(path) - 3] = '\0';
    }
    if (path[strlen(path) - 2] == '.' &&
        path[strlen(path) - 1] == 'Z') {
        path[strlen(path) - 2] = '\0';
    }

#if 0
    fprintf(stderr, "\n===============================================\n");
    fprintf(stderr, "File name = %s\n", options->file_name);
    fprintf(stderr, "\n===============================================\n");
#endif


    /*
     *  If external lookup table not specified, determine from PDMF
     *  string which internal table to use (POP5, POP6, etc.).
     */
    if (table == 0)
    {
        if (Verbose > 2)
            (void)fprintf(stderr, "Mobility file name: %s\n", options->chemistry);

        table = parse_mobility_file_name(options->chemistry, options, message,
            table);

        if (table == 0)
        {
            table = Btk_get_3730pop7_table();
            options->lut_type = ABI3730pop7;
            (void)fprintf(stderr,
            "Can't select the lookup table automatically. \n");
            (void)fprintf(stderr,
            "Using a built-in ABI 3730 Pop-7 table.\n");
        }
    }

    /* Don't call compute_qv (that is, use original bases, locs and QVs) if:
     * -nocall option used AND there are QVs in the original sample file.
     */
#if 0
    fprintf(stderr, "Before Btk_compute_qv: options->file_name=%s\n", options->file_name);
#endif
    if (!options->nocall || (quality_values == NULL))
    {

        if (!options->inp_phd && (quality_values == NULL)) {
            quality_values = CALLOC(uint8_t, num_called_bases);
            MEM_ERROR(quality_values);
        }

        if ( Btk_compute_qv(&num_called_bases, &called_bases, &called_peak_locs,
            &num_datapoints, chromatogram, "ACGT", table,
#if USE_CONTEXT_TABLE
                        ctable,
#endif
                        &quality_values, *options, message, &results ) == ERROR)
        {
            sprintf(status_code, "%s", "TT_TRASH");
            if (OutputFourMultiFastaFiles)
                output_four_multi_fasta_files(multiseqsFileName, multiqualFileName,
                multilocsFileName, multistatFileName, num_called_bases,
                called_bases, quality_values, called_peak_locs,
                results.frac_QV20_with_shoulders, status_code, *options);
            if (Verbose > 1)
                fprintf(stderr, "0 bases finally\n");
            status_code[0] = '\0';
            goto error;
        }
    }
    else if (options->nocall && (Verbose != 0))
        fprintf(stderr, "Using original base calls and quality values\n");

    if ((Verbose > 1) && options->process_bases) {
        fprintf(stderr, "%d bases finally\n", num_called_bases);
    }
    message->text[0] = '\0';  /* message may have been set with no error */

    /* Count and output number of bases of QV>=20 */
    {
        int count, num_qv20=0;
        for (count=0; count<num_called_bases; count++)
            if (quality_values[count] >= 20)
                num_qv20++;
        if (Verbose > 1)
            fprintf(stderr, "%d bases of QV >= 20\n", num_qv20);
    }

    trimmed_read_length = find_trim_points(num_called_bases, quality_values,
        trim_window, trim_threshold, &left_trim_point, &right_trim_point);

    if ((options->tal_dir[0] != '\0') && !options->indel_resolve) {
        if (Btk_output_tal_file(path     ,
            AlnType == NAME_DIR ? options->tal_dir : NULL,
            ConsensusName, ConsensusSeq,
            called_bases, num_called_bases,
            Match, MisMatch, Insertion, Deletion, (float)RepeatFraction,
            Verbose) == ERROR)
        {
            goto error;
        }
    }

    if ((options->hpr_dir[0] != '\0') && !options->indel_resolve) {
        if ((r = Btk_output_hpr_file(path,
            HprType == NAME_DIR ? options->hpr_dir : NULL,
            called_bases, called_peak_locs, quality_values, num_called_bases,
            num_datapoints, Verbose))
            == ERROR)
        {
            goto error;
        }
    }

    if (OutputPhd && !options->indel_resolve) {
        if ((r = Btk_output_phd_file(path     ,
            PhdType == NAME_DIR ? PhdDirName : NULL,
            called_bases, called_peak_locs, quality_values, num_called_bases,
            num_datapoints, options->nocall, options->chemistry,
            left_trim_point, right_trim_point, trim_threshold, Verbose))
            == ERROR)
        {
            goto error;
        }
    }

    if (OutputQual && !options->indel_resolve) {
        if ((r = Btk_output_quality_values(QualType, path,
            QualDirName, multiqualFileName,
            quality_values, num_called_bases, left_trim_point, right_trim_point,
            Verbose)) == ERROR)
        {
            goto error;
        }
    }

    if (OutputFasta && !options->indel_resolve) {
        if ((r = Btk_output_fasta_file(FastaType, path     ,
            FastaDirName, multiseqFileName,
            called_bases, num_called_bases, left_trim_point, right_trim_point,
            Verbose)) == ERROR)
        {
            goto error;
        }
    }

    if (OutputSCF && !options->indel_resolve) {
        if (output_scf_file(path, SCFType == NAME_DIR ? SCFDirName : ".",
            called_bases, called_peak_locs, quality_values,
            num_called_bases, num_datapoints, chromatogram[0],
            chromatogram[1], chromatogram[2], chromatogram[3],
            "ACGT", options->chemistry) == ERROR)
        {
            goto error;
        }
    }

    if (OutputQualRpt && !options->indel_resolve) {
        accum_qual_report(&Qual_data, quality_values, num_called_bases,
                        trimmed_read_length);
    }

    if (OutputFourMultiFastaFiles)
    {
        sprintf(status_code, "%s", "TT_SUCCESS");

        output_four_multi_fasta_files(multiseqsFileName, multiqualFileName,
            multilocsFileName, multistatFileName, num_called_bases,
            called_bases, quality_values, called_peak_locs,
            results.frac_QV20_with_shoulders, status_code, *options);
        status_code[0] = '\0';
    }

    Btk_release_file_data(called_bases, called_peak_locs, quality_values,
        chromatogram, &call_method, &(options->chemistry));

    if (options->time > 0) {
        curr_clock = clock();
        fprintf(stderr,
        "File %s processed in %f sec. \n",
            path     ,
            (float)(curr_clock - start_clock)/(float)CLOCKS_PER_SEC);
        start_clock = curr_clock;
    }

    if (consFromSample) {
        FREE(ConsensusSeq);
    }
    return SUCCESS;

error:
    if (chromatogram[0] != NULL) {
        Btk_release_file_data(called_bases, called_peak_locs, quality_values,
            chromatogram, &call_method, &options->chemistry);
    }

    if (consFromSample) {
        FREE(ConsensusSeq);
    }

    return ERROR;
}

static void
process_sffRead(sffHeader *h, sffRead *r,  Options *options)
{
    int i, j;
    int left, right;
    float leftCut, rightCut, sig;

    if (OutputFasta)
    {
        int length = r->clip_quality_right-r->clip_quality_left+1;
        printf(">%s length=%d\n", r->name, length); 
  
#if 0 
        for (i=r->clip_quality_left,j=0; i <= r->clip_quality_right; i++) {
            fputc(( i >= r->clip_quality_left && i <= r->clip_quality_right
               ? toupper(r->bases[i-1]) : tolower(r->bases[i-1]) ), stdout);
            if (++j == 60) {
                fputc('\n', stdout);
                j = 0;
            }
        }
        if (length % 60 != 0)
            fputc('\n', stdout);
#else
        for (i=0; i<r->number_of_bases; i++)
        {
            fputc(toupper(r->bases[i]), stdout);
            if (++j == 60) {
                fputc('\n', stdout);
                j = 0;
            }
        }
        if (length % 60 != 0) {
            fputc('\n', stdout); 
        }
#endif
    }
    else if (OutputQual)
    {
        int length = 0;
        printf(">%s\n", r->name);
        
        for (i=r->clip_quality_left,j=0; i <= r->clip_quality_right; i++) {
            printf("%2d", r->quality_values[i-1]);
            if (++j == 20) {
                fputc('\n', stdout);
                j = 0;
                length = 0;
            } else if (i < r->clip_quality_right) {
                fputc(' ', stdout);
                length += 3;
            }
        }
        if ((length+3)%60 > 3) 
            fputc('\n', stdout);
    }
    else if (OutputFlow)
    {
        int *flowIndex = CALLOC(int, r->number_of_bases);
        printf(">%s\n", r->name);

//      printf("Number of bases= %d Flow index per base:\n", r->number_of_bases);
        for (i=0; i<r->number_of_bases; i++)
        {
            flowIndex[i] = (int)r->flow_index_per_base[i];
        }
        for (i=0; i<r->number_of_bases; i++)
        {
            if (i > 0) flowIndex[i] += flowIndex[i-1];
//          printf("%d ", flowIndex[i]);
        }
//      printf("\n");

        left = 0;
        leftCut = 0.0;
        right = h->number_of_flows_per_read - 1;
        rightCut = 0.0;

#if 0
        left = flowIndex[r->clip_quality_left-1] - 1;
        leftCut = 0.0;
        for (i=2; r->clip_quality_left - i >= 0 && flowIndex[r->clip_quality_left-i] - 1 == left; i++) {
            leftCut += 1.0;
        }
        right = flowIndex[r->clip_quality_right-1] - 1;
        rightCut = 0.0;

        for (i=0; r->clip_quality_right + i < r->number_of_bases && 
                  flowIndex[r->clip_quality_right+i] - 1 == right; i++) 
        {
            rightCut += 1.0;
        }
#endif

//      printf("Number of bases= %d Flow index per base:\n", r->number_of_bases);
        
        j = 0;
        for (i=0; i<h->number_of_flows_per_read; i++)
        {
            sig = (float)r->flowgram_values[i]/100.0 
                - (i == left ? leftCut : (i == right ? rightCut : 0.0));
            printf("%c,%.2f", h->flow_chars[i], sig); 
            if (++j == 10) {
                printf("\n");
                j = 0;
            } else if (i < h->number_of_flows_per_read) {
                printf(" ");
            }
        }  
        printf("\n");
    }
}

static int
process_sff_file(char *sfffile, Options *options, BtkMessage *message)
{
    int i;
    char *sff_name;
    FILE *sff;
    sffHeader   *h   = CALLOC(sffHeader,   1);
    sffManifest *m   = CALLOC(sffManifest, 1);
    sffRead     *r   = CALLOC(sffRead,     1);

    if (sfffile[0] == '\0') {
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

    /* Exclude invalid command line options for 454 data */
    if (options->het || options->mix)
    {
        fprintf(stderr, "\nCan not use option -het or -mix with 454 data\n");
        exit(2);
    }
    if (options->edited_bases)
    {
        fprintf(stderr, "\nCan not use option -edited_bases with 454 data\n");
        exit(2);
    }
    if (options->indel_detect || options->indel_resolve)
    {
        fprintf(stderr, "\nCan not use option -indel_detect or -indel_resolve with 454 data\n");
        exit(2);
    }
    if (OutputSCF)           
    {
        fprintf(stderr, "\nCan not use option -c or -cd <dir> with 454 data\n");
        exit(2);
    } 
    if (options->tab_dir[0] != '\0')
    {
        fprintf(stderr, "\nCan not use option -tab or -tabd <dir> with 454 data\n");
        exit(2);
    }
    if (options->tal_dir[0] != '\0')
    {
        fprintf(stderr, "\nCan not use option -tal or -tald <dir> with 454 data\n");
        exit(2);
    }
    if (options->tip_dir[0] != '\0')
    {
        fprintf(stderr, "\nCan not use option -tip or -tipd <dir> with 454 data\n");
        exit(2);
    }
    if (options->poly || options->poly_dir[0] != '\0')
    {
        fprintf(stderr, "\nCan not use option -d or -dd <dir> with 454 data\n");
        exit(2);
    }
    if (options->lut_type == ABI3730pop7 || options->lut_type == ABI3700pop5 ||
        options->lut_type == ABI3700pop6 || options->lut_type == ABI3100     ||
        options->lut_type == MegaBACE)
    {
        fprintf(stderr, "\nCan not use the specified Sanger lookup table with 454 data\n");
        exit(2);
    } 
 
    for (i=0; i < h->number_of_reads; i++) {
        readsff_read(sff, h, r);

        process_sffRead(h, r, options);
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

    if (options->time > 0) {
        curr_clock = clock();
        fprintf(stderr,
        "File %s processed in %f sec. \n",
            sfffile  ,
            (float)(curr_clock - start_clock)/(float)CLOCKS_PER_SEC);
        start_clock = curr_clock;
    }

    FREE(h->data_block);
    FREE(h);
    FREE(m->manifest);
    FREE(m);
    FREE(r->data_block);
    FREE(r);
    return SUCCESS;

error:
    FREE(h->data_block);
    FREE(h);
    FREE(m->manifest);
    FREE(m);
    FREE(r->data_block);
    FREE(r);
    return ERROR;
}

/*
 * This function processes a single sample file.  Its synopsis is:
 *
 * result = process_file(table, ctable, file_name, ..., message)
 *
 * where
 *	table		is the address of a BtkLookupTable returned by
 *			Btk_read_lookup_table()
 *	file_name	is the name (path) of the sample file
 *	message		is the address of a BtkMessage where information about
 *			an error will be put, if any
 *
 *	result		is 0 on success, !0 if an error occurs
 */
static int
process_file(BtkLookupTable *table, ContextTable *ctable,
    char *path, int argind, int num_arg, char *ConsensusName, char *ConsensusSeq,
    Options *options, BtkMessage *message)
{
    char *seq_name;

    file_type = -1;

    if (path[0] == '\0') {
	return SUCCESS;
    }
    message->text[0] = '\0';

#ifdef __WIN32
    if ((seq_name = strrchr(path, '\\')) != NULL) {
#else
    if ((seq_name = strrchr(path, '/')) != NULL) {
#endif
        seq_name++;
    }
    else {
        seq_name = path;
    }
 
    strcpy(options->file_name, seq_name);

    if (process_sff_file(path, options, message) !=
        kWrongFileType)
    {
        // this is SFF file
    }
    else if (process_sample_file(table, ctable, path, ConsensusName,
            ConsensusSeq, options, message) != kWrongFileType) 
    {
        // this is either ABI or SCF file
    }
    else { 
       fprintf(stderr, "%s: not an ABI, SCF or SFF file\n", path);
    }

    return SUCCESS;
}


/*
 * This function processes all sample files listed in the specified file.
 * Its synopsis is:
 *
 * result = process_fileoffiles(table, ctable, ..., message)
 *
 * where
 *	table		is the address of a BtkLookupTable returned by
 *			Btk_read_lookup_table()
 *      ctable          is the address of a CntaxtTable returned by
 *                      read_context_table()
 *	fileoffiles	is the name (path) of a file with one sample file
 *			per line
 *	message		is the address of a BtkMessage where information about
 *			an error will be put, if any
 *
 *	result		is 0 on success, !0 if an error occurs
 */
static int
process_fileoffiles(BtkLookupTable *table, ContextTable *ctable, 
    char *fileoffiles, char *ConsensusName, char *ConsensusSeq,
    Options *options, BtkMessage *message)
{
    FILE *fp;
    char line[BUFLEN], *s;
    int r;
    struct stat statbuf;


    if (Verbose > 1) {
	(void)fprintf(stderr, "%s: opening file-of-files\n", fileoffiles);
    }
    message->text[0] = '\0';
    if ((fp = fopen(fileoffiles, "r")) == NULL) {
	error(fileoffiles, "couldn't open", errno);
	return ERROR;
    }

    /* Read one line of the file at a time */
    while (fgets(line, sizeof(line), fp) != NULL) {
	/* Trim the newline, if any */
	if ((s = strchr(line, '\n')) != NULL) {
	    *s = '\0';
	}

	if (stat(line, &statbuf) != 0) {
	    error(line, "can't stat", errno);
	    continue;
	}
	if (statbuf.st_mode & S_IFDIR) {
	    /* skip subdirectories */
	    fprintf(stderr, "%s: skipping subdirectory\n", line);
	    continue;
	}

	if (process_file(table, ctable, line, -1, -1, ConsensusName, 
            ConsensusSeq, options, message) != SUCCESS) 
        {
	    fprintf(stderr, "%s: %s\n", line, message->text);
	}
	fprintf(stderr,"\n");
    }
    if (ferror(fp)) {
	error(fileoffiles, "couldn't read", errno);
	r = ERROR;
	goto error;
    }

    (void)fclose(fp);
    if (Verbose > 1) {
	(void)fprintf(stderr, "%s: closed file-of-files\n", fileoffiles);
    }
    return SUCCESS;

error:
    (void)fclose(fp);
    return r;
}

/*
 * This function processes all files contained in the specified directory.
 * Its synopsis is:
 *
 * result = process_dir(table, ctable, dir, ...)
 *
 * where
 *	table		is the address of a BtkLookupTable returned by
 *			Btk_read_lookup_table()
 *      ctable          is the address of a ContextTable returned by
 *                      read_context_table()
 *	dir		is the name (path) of a directory
 *
 *	Routine prints out its own error messages as needed.
 */
static void
process_dir(BtkLookupTable *table, ContextTable *ctable,
    char *dir, char *ConsensusName, char *ConsensusSeq,
    Options *options, BtkMessage *message)
{
#ifndef __DEVSTUDIO
    DIR *d;
    struct dirent *de;
    struct stat statbuf;
#else
    long handle;
    char filespec[MAXPATHLEN];
    struct _finddata_t fileinfo;
    struct _stat buffer;
#endif
    char path_and_name[MAXPATHLEN];

#ifdef __DEVSTUDIO
    sprintf(filespec, "%s\\*.*", dir);
    handle = _findfirst(filespec, &fileinfo);

    do {
        sprintf(path_and_name, "%s\\%s", dir, fileinfo.name);

        if (_stat(path_and_name, &buffer) != 0) {
            error(path_and_name, "can't stat", errno);
            continue;
        }

        if (buffer.st_mode & _S_IFDIR) {    /* skip subdirectories */
            if (Verbose > 2) {
                fprintf(stderr, "%s: skipping subdirectory\n", path_and_name);
            }
            continue;
        }

        if (process_file(table, ctable, path_and_name, -1, -1, ConsensusName, 
            ConsensusSeq, options, message) != SUCCESS)
        {
            fprintf(stderr, "%s: %s\n\n", path_and_name, message->text);
        }

        fprintf(stderr, "\n");
    } while ((_findnext(handle, &fileinfo)) == 0);

    _findclose(handle);
#else
    if (Verbose > 1) {
        fprintf(stderr, "%s: opening dir\n", dir);
    }

    if ((d = opendir(dir)) == NULL) {
        error(dir, "couldn't open dir", errno);
        return;
    }

    while ((de = readdir(d)) != NULL) {
#ifdef __WIN32
        sprintf(path_and_name, "%s\\%s", dir, de->d_name);
#else
        sprintf(path_and_name, "%s/%s", dir, de->d_name);
#endif

        if (stat(path_and_name, &statbuf) != 0) {
            error(path_and_name, "can't stat", errno);
            continue;
        }

        if (statbuf.st_mode & S_IFDIR) {    /* skip subdirectories */
            if (Verbose > 2) {
                fprintf(stderr, "%s: skipping subdirectory\n", path_and_name);
            }
            continue;
        }

        if (process_file(table, ctable, path_and_name, -1, -1, ConsensusName, 
            ConsensusSeq, options, message) != SUCCESS)
        {
            fprintf(stderr, "%s: %s\n\n", path_and_name, message->text);
        }

    }

    closedir(d);

    if (Verbose > 1) {
        fprintf(stderr, "%s: closed dir\n", dir);
    }
#endif

    return;
}

/*******************************************************************************
 * Function: validateDirectory
 *******************************************************************************
 */
void validateDirectory(char *targetDirName, Options *options)
{
    struct stat targetDirStat;
    FILE *fp;
    char fname[1024];

    if (stat(targetDirName, &targetDirStat))
        switch (errno)
        {
            case ENOENT:
                fprintf(stderr, "\nError: Target directory %s does not exist."
                        "\nPlease create it (or check your typing), "
                        "and try again.\n\n", targetDirName);
                exit_message(options, EXIT_FAILURE);
            case ENOTDIR:
                fprintf(stderr, 
                        "\nError: Path specified for target directory %s "
                                "is invalid.\n\n", targetDirName);
                exit_message(options, EXIT_FAILURE);
            case EACCES:
                fprintf(stderr, 
                      "\nError: Access to target directory %s not allowed.\n\n",
                        targetDirName);
                exit_message(options, EXIT_FAILURE);
            default:
                fprintf(stderr, "\nError: Could not open target directory %s, "
                                "errno = %d\n\n", targetDirName, errno);
                exit_message(options, EXIT_FAILURE);
        }

    if (!(targetDirStat.st_mode & S_IFDIR))
    {
        fprintf(stderr, "\nError: Target directory %s is a file.\n\n", 
                targetDirName);
        exit_message(options, EXIT_FAILURE);
    }
#ifdef __WIN32
    sprintf(fname, "%s\\A__hIgHlY___unLIkeLY_NamE", targetDirName);
#else
    sprintf(fname, "%s/A__hIgHlY___unLIkeLY_NamE", targetDirName);
#endif
    fp = fopen(fname, "w");
    if (fp == NULL)
    {
        fprintf(stderr, "\nError: Cannot write to target directory %s due "
                        "to insufficient permissions.\n\n", targetDirName);
        exit_message(options, EXIT_FAILURE);
    }
    else {
        fclose(fp);
        unlink(fname);
    }
}

int
main(int argc, char *argv[])
{
    char           *args, *lut_name = NULL, *context_table = NULL;
    int             i, j, optind, listtype=0;
    struct stat     statbuf;
    BtkMessage      message;
    BtkLookupTable *table = NULL;
    ContextTable   *ctable = NULL;
    char	   *ConsensusSeq = NULL;
    Options         options;

    /* Set time */ 
    start_clock = clock();

    if (argc == 1) {
	usage(argc, argv);
	exit(2);
    }

    /* Set defaults */
    dev             = 0;
    opts            = 1;
    Verbose         = 1;
    options.Verbose = 1;
    options.process_bases = 0;
    InputName[0]    = '\0';
    InputType       = NAME_FILES;
    options.inp_phd = 0;
    options.inp_phd_dir[0]='\0';
    ConsensusSpecified = 0;
    ConsensusName[0]='\0';
    AlnType         = NAME_FILES;
    HprType         = NAME_FILES;
    OutputFasta     = 0;
    FastaDirName[0] = '\0';
    FastaType       = NAME_NONE;
    OutputSCF       = 0;
    SCFDirName[0]   = '\0';
    SCFType         = NAME_FILES;
    OutputPhd       = 0;
    PhdDirName[0]   = '\0';
    PhdType         = NAME_FILES;
    OutputQual      = 0;
    OutputQualRpt   = 0;
    QualDirName[0]  = '\0';
    QualType        = NAME_NONE;
    lut_name        = NULL;
    context_table   = NULL;
    options.chemistry    = NULL;
    options.indel_detect = 0;
    options.indel_resolve= 0;
    options.indloc       = -1;
    options.indsize      = 0;
    options.lut_type     = 0;
    options.nocall       = 0;
    options.recalln      = 0;
    options.recallndb    = 0;
    options.ladder       = 0;
    options.edited_bases = 0;
    options.gauss        = 1;
    options.shift        = 0;
    options.renorm       = 0;
    options.respace      = 0;
    options.raw_data     = 0;
    options.multicomp    = 0;
    options.xgr          = 0;
    options.scf_dir[0]   = '\0';
    options.scf_version  = 0;
    options.time         = 0;
    options.tip_dir[0]   = '\0';
    options.tab_dir[0]   = '\0';
    options.tal_dir[0]   = '\0';
    options.hpr_dir[0]   = '\0';
    options.poly         = 0;
    options.poly_dir[0]  = '\0';
    options.het          = 0;
    options.mix          = 0;
    options.min_ratio    = (float)0.1;           /* default */
    RepeatFraction       = 0.9;
    Match                = 20;
    MisMatch             = -5;
    Insertion            = -3;
    Deletion             = -3;
    multiqualFileName[0] = '\0';
    multiseqFileName[0]  = '\0';
    MultiFastaFilesDirName[0] = '\0';
    status_code[0]            = '\0';
    OutputFourMultiFastaFiles = 0;


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
        if (((optind == argc - 1) || (argv[optind + 1][0] == '-')) &&
            ((strcmp(argv[optind], "-C" )             == 0) ||
             (strcmp(argv[optind], "-cd")             == 0) ||
             (strcmp(argv[optind], "-ct")             == 0) ||
             (strcmp(argv[optind], "-dd"  )           == 0) ||
             (strcmp(argv[optind], "-ipd")            == 0) ||
             (strcmp(argv[optind], "-id")             == 0) || 
             (strcmp(argv[optind], "-if")             == 0) ||
             (strcmp(argv[optind], "-indloc")         == 0) ||
             (strcmp(argv[optind], "-indsize")        == 0) ||
             (strcmp(argv[optind], "-min_ratio")      == 0)  ||
             (strcmp(argv[optind],  "-o")             == 0) ||
             (strcmp(argv[optind], "-pd")             == 0) ||
             (strcmp(argv[optind], "-qd")             == 0) ||
             (strcmp(argv[optind], "-qa")             == 0) ||
             (strcmp(argv[optind], "-qr")             == 0) ||
             (strcmp(argv[optind], "-sd")             == 0) ||
             (strcmp(argv[optind], "-sa")             == 0) ||
             (strcmp(argv[optind], "-t" )             == 0) ||
             (strcmp(argv[optind], "-tipd")           == 0) ||
             (strcmp(argv[optind], "-tabd")           == 0) ||
             (strcmp(argv[optind], "-tald")           == 0) ||
             (strcmp(argv[optind], "-hprd")           == 0) ||
             (strcmp(argv[optind], "-trim_window")    == 0) ||
             (strcmp(argv[optind], "-trim_threshold") == 0)))
        {
            usage(argc, argv);
            fprintf(stderr, "\nInvalid option specified.\n");
            exit(2);
        }

        /* Make sure that all the flags are from alowed list
         */
        if (((optind == argc-1) || (argv[optind + 1][0] == '-')) &&
            ((strcmp(argv[optind], "-h")            != 0) &&
             (strcmp(argv[optind], "-dev")          != 0) && 
             (strcmp(argv[optind], "-opts")         != 0) &&
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
             (strcmp(argv[optind], "-indel_detect") != 0) &&
             (strcmp(argv[optind], "-indel_resolve")!= 0) &&
             (strcmp(argv[optind], "-mc")           != 0) &&
             (strcmp(argv[optind], "-xgr")          != 0) &&
             (strcmp(argv[optind], "-shift")        != 0) &&
             (strcmp(argv[optind], "-convolved")    != 0) &&
             (strcmp(argv[optind], "-het")          != 0) &&
             (strcmp(argv[optind], "-mix")          != 0) &&
             (strcmp(argv[optind], "-cv3")          != 0) &&
             (strcmp(argv[optind], "-3700pop5")     != 0) &&
             (strcmp(argv[optind], "-3700pop6")     != 0) &&
             (strcmp(argv[optind], "-3730")         != 0) &&
             (strcmp(argv[optind], "-3100")         != 0) &&
             (strcmp(argv[optind], "-mbace")        != 0) &&
             (strcmp(argv[optind], "-p")            != 0) &&
             (strcmp(argv[optind], "-q")            != 0) &&
             (strcmp(argv[optind], "-c")            != 0) &&
             (strcmp(argv[optind], "-d")            != 0) &&
             (strcmp(argv[optind], "-s")            != 0) &&
             (strcmp(argv[optind], "-f")            != 0) &&
             (strcmp(argv[optind], "-tip")          != 0) &&
             (strcmp(argv[optind], "-tal")          != 0) &&
             (strcmp(argv[optind], "-hpr")          != 0) &&
             (strcmp(argv[optind], "-tab")          != 0)))
        {
            usage(argc, argv);
            fprintf(stderr, "\nInvalid flag specified.\n");
            exit(2);
        }

        /* Parse the current argument string (args), starting from the 2nd
         * character (j==1), which is not '-'.
         */
        for (j = 1, args = argv[optind]; args[j] != '\0'; j++) 
        {
            i = args[j];
            switch (i) {
            case '3':
                if (strcmp(args, "-3730") == 0) {
                    if (lut_name == NULL)
                        options.lut_type = ABI3730pop7;
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                if (strcmp(args, "-3700pop5") == 0) {
                    if (lut_name == NULL)
                        options.lut_type = ABI3700pop5;
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                if (strcmp(args, "-3700pop6") == 0) {
                    if (lut_name == NULL)
                        options.lut_type = ABI3700pop6;
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                if (strcmp(args, "-3100") == 0) {
                    if (lut_name == NULL)
                        options.lut_type = ABI3100;
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else {
                    usage(argc, argv);
                    fprintf(stderr, "\nInvalid option specified.\n");
                    exit(2);
                }
            case 'C':
                if (strcmp(args, "-C") == 0) {
                    ConsensusSpecified++;
                    (void)strncpy(ConsensusName, argv[++optind],
                                  sizeof(ConsensusName));
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else {
                    usage(argc, argv);
                    fprintf(stderr, "\nInvalid option specified.\n");
                    exit(2);
                }

            case 'd':
                if (strcmp(args, "-dev") == 0) {
                    if (argc == 2)
                    {
                        usage_dev(argc, argv);
                        exit(2);
                    }
                    else
                       dev++;
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-d") == 0) {
                    options.poly=1;
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-dd") == 0) {
                    options.poly=2;
                    (void)strncpy(options.poly_dir, argv[++optind],
                        sizeof(options.poly_dir));
                    validateDirectory(options.poly_dir, &options);
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }

                switch (listtype) {
                case 'i':
                    InputType = NAME_DIR;
                    (void)strncpy(InputName, argv[++optind],
                                  sizeof(InputName));
                    break;
                case 'p':
                    PhdType = NAME_DIR;
                    (void)strncpy(PhdDirName, argv[++optind],
                                  sizeof(PhdDirName));
                    validateDirectory(PhdDirName, &options);
                    break;
                case 'c':
                    SCFType = NAME_DIR;
                    (void)strncpy(SCFDirName, argv[++optind],
                                  sizeof(SCFDirName));
                    validateDirectory(SCFDirName, &options);
                    (void)strcpy(options.scf_dir, ".");
                    break;
                default:
                    usage(argc, argv);
                    fprintf(stderr, "\nInvalid option specified.\n");                 
                    exit(2);
                }
                break;

            case 'e':
                if (strcmp(args, "-edited_bases") == 0) {
                    options.edited_bases++;
                    j = strlen(args)-1;  /* break out of inner loop */
                    break;
                }

            case 'f':
                if (strcmp(args, "-f") == 0) {
                    OutputFlow++;
                    j = strlen(args)-1;  /* break out of inner loop */
                    break;
                }
                switch (listtype) {
                case 'i':
                    InputType = NAME_FILEOFFILES;
                    (void)strncpy(InputName, argv[++optind],
                                  sizeof(InputName));
                    break;
                default:
                    usage(argc, argv);
                    fprintf(stderr, "\nInvalid option specified.\n");
                    exit(2);
                }
                break;

            case 'F':  /* Version 1.0 had a -F flag which is now the default */
                break; /* Catch this case here so it's not flagged as an error*/

            case 'h':
                if (strcmp(args, "-het") == 0) {
                    options.het++;
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-hpr") == 0) {
                    HprType = NAME_DIR;
                    strcpy(options.hpr_dir, ".");
                    validateDirectory(options.hpr_dir, &options);
                    HprType = NAME_DIR;
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-hprd") == 0) {
                    HprType = NAME_DIR;
                    (void)strncpy(options.hpr_dir, argv[++optind],
                                  sizeof(options.hpr_dir));
                    validateDirectory(options.hpr_dir, &options);
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                help_message(argc, argv);
                exit(2);

            case 'i':
                if (strcmp(args, "-indel_detect") == 0){
                    options.indel_detect = 1;
                    options.recalln = 1;
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-indel_resolve") == 0){
                    options.indel_detect++;
                    options.indel_resolve++;  
                    options.recalln = 1;
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-indloc") == 0){
                    options.indloc = (int)atoi(argv[++optind]);
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-indsize") == 0){
                    options.indsize = (int)atoi(argv[++optind]);
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-ipd") == 0) {
                    options.inp_phd = 1;
                    (void)strncpy(options.inp_phd_dir, argv[++optind],
                              sizeof(options.inp_phd_dir));
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                listtype = i;
                break;

            case 'l':
                if (strcmp(args, "-ladder") == 0)
                    options.ladder++;
                else {
                    usage(argc, argv);
                    exit(2);
                }
                j = strlen(args) - 1;  /* break out of inner loop */
                break;

            case 'm':
                if (strcmp(args, "-mix") == 0) {
                    options.mix++;
                }
                else if (strcmp(args, "-min_ratio") == 0) {
                    options.min_ratio = (float)atof(argv[++optind]);   
                }
                else if (strcmp(args, "-mc") == 0) {
                    options.multicomp++;
                }
                else if (strcmp(args, "-mbace") == 0) {
                    if (lut_name == NULL)
                        options.lut_type = MegaBACE;
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else
                {
                    usage(argc, argv);
                    exit(2);
                }
                j = strlen(args) - 1;  /* break out of inner loop */
                break;

            case 'n':
                if (strcmp(args, "-nocall" ) == 0)
                    options.nocall++;
                else {
                    usage(argc, argv);
                    exit(2);
                }
                j = strlen(args) - 1;   /* break out of inner loop */
                break;

             case 'o':
                if (strcmp(args, "-opts" ) == 0)
                {
                    opts++;
                    j = strlen(args) - 1;   /* break out of inner loop */
                }
                else if (strcmp(args, "-o" ) == 0)
                {
                    OutputFourMultiFastaFiles++;
                    strncpy(MultiFastaFilesDirName, argv[++optind],
                         sizeof(MultiFastaFilesDirName));
                    validateDirectory(MultiFastaFilesDirName, &options);
                    j = strlen(args) - 1;   /* break out of inner loop */
                }
                else  {
                    usage(argc, argv);
                    exit(2);
                }
                break;

            case 'p':
                listtype = i;
                OutputPhd++;
                break;

            case 'q':
                if (strcmp(args, "-q") == 0) {
                    OutputQual++;
                    QualType |= NAME_FILES;
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-qa") == 0) {
                    OutputQual++;
                    QualType |= NAME_MULTI;
                    strncpy(multiqualFileName, argv[++optind],
                            sizeof(multiqualFileName));
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-qd") == 0) {
                    OutputQual++;
                    QualType |= NAME_DIR;
                    strncpy(QualDirName, argv[++optind], sizeof(QualDirName));
                    validateDirectory(QualDirName, &options);
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-qr") == 0) {
                    OutputQualRpt++;
                    strncpy(QualRptName, argv[++optind], sizeof(QualRptName));
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else 
                {
                    usage(argc, argv);
                    exit(2);
                }
                break;

            case 'r':
                if (strcmp(args, "-respace") == 0)
                    options.respace++;
                else if (strcmp(args, "-renorm") == 0)
                    options.renorm++;
                else if (strcmp(args, "-recalln") == 0)
                    options.recalln++;
                else if (strcmp(args, "-recallndb") == 0)
                    options.recallndb++;
                else if (strcmp(args, "-raw") == 0) 
                    options.raw_data++;         
                else {
                    usage(argc, argv);
                    exit(2);
                }
                j = strlen(args) - 1;  /* break out of inner loop */
                break;

            case 'Q':
                Verbose = 0;
                options.Verbose = 0;
                break;

            case 's':
                if (strcmp(args, "-shift") == 0) {
                    options.shift++;
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-s") == 0) {
                    OutputFasta++;
                    FastaType |= NAME_FILES;
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-sa") == 0) {
                    OutputFasta++;
                    FastaType |= NAME_MULTI;
                    strncpy(multiseqFileName, argv[++optind],
                            sizeof(multiseqFileName));
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-sd") == 0) {
                    OutputFasta++;
                    FastaType |= NAME_DIR;
                    strncpy(FastaDirName, argv[++optind], sizeof(FastaDirName));
                    validateDirectory(FastaDirName, &options);
                    j = strlen(args) - 1;   /* break out of inner loop */
                }
                break;

            case 't':
                if (strcmp(args, "-time") == 0) {
                    options.time=1;
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-tip") == 0) {
                    strcpy(options.tip_dir, ".");
                    validateDirectory(options.tip_dir, &options);
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-tipd") == 0) {
                    (void)strcpy(options.tip_dir, argv[++optind]);
//                      sizeof(options.tip_dir));
                    validateDirectory(options.tip_dir, &options);
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-tab") == 0) {
                    strcpy(options.tab_dir, ".");
                    validateDirectory(options.tab_dir, &options);
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-tabd") == 0) {
                    (void)strncpy(options.tab_dir, argv[++optind], 
                        sizeof(options.tab_dir));
                    validateDirectory(options.tab_dir, &options);
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-tal") == 0) {
                    AlnType = NAME_DIR;
                    strcpy(options.tal_dir, ".");
                    validateDirectory(options.tal_dir, &options);
                    AlnType = NAME_DIR;
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-tald") == 0) {
                    AlnType = NAME_DIR;
                    (void)strncpy(options.tal_dir, argv[++optind],
                                  sizeof(options.tal_dir));
                    validateDirectory(options.tal_dir, &options);
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else if (strcmp(args, "-trim_window") == 0) {
                    trim_window = atoi(argv[++optind]);
                    j = strlen(args) - 1;  /* break out of inner loop */
                    if (trim_window <= 0) {
                        usage(argc, argv);
                        exit(2);
                    }
                    break;
                }
                else if (strcmp(args, "-trim_threshold") == 0) {
                    trim_threshold = (float)atof(argv[++optind]);
                    j = strlen(args) - 1;  /* break out of inner loop */
                    if (trim_threshold <= 0) {
                        usage(argc, argv);
                        exit(2);
                    }
                    break;
                }
                else if (strcmp(args, "-t") == 0){
                    lut_name = argv[++optind];
                    if ((strstr(lut_name, "3700") != NULL) ||
                        (strstr(lut_name,  "Pop5") != NULL) ||
                        (strstr(lut_name,  "POP5") != NULL)) {
                        options.lut_type = ABI3700pop5;
                    } else if ((strstr(lut_name, "Pop6") != NULL) &&
                        (strstr(lut_name,  "POP6") != NULL)) {
                        options.lut_type = ABI3700pop6;
                    } else if (strstr(lut_name, "3100") != NULL) {
                        options.lut_type = ABI3100;
                    } else if (strstr(lut_name, "3730") != NULL) {
                        options.lut_type = ABI3730pop7;
                    } else if (strstr(lut_name, "Mega") != NULL ||
                               strstr(lut_name, "mega") != NULL ||
                               strstr(lut_name, "BACE") != NULL ||
                               strstr(lut_name, "bace") != NULL)
                    {
                        options.lut_type = MegaBACE;
                    }
                    else {
                        fprintf(stderr,
                        "Type of external table can not be determined; ");
                        fprintf(stderr,
                        "assuming it is 3730pop7\n");
                         options.lut_type = ABI3730pop7;
                    }

                }
                else {
                    usage(argc, argv);
                    exit(2);
                }
                break;

            case 'V':
                Verbose++;
                options.Verbose++;
                break;

            case 'c':
                if (strcmp(args, "-ct") == 0){
                    context_table = argv[++optind];
                    j = strlen(args) - 1;   /* break out of inner loop */
                }
                else if (strcmp(args, "-cv3") == 0) {
                    options.scf_version = 3;
                    j = strlen(args) - 1;   /* break out of inner loop */
                }
                else if (strcmp(args, "-convolved") == 0) {
                    options.gauss = 0;
                    j = strlen(args) - 1;   /* break out of inner loop */
                }
                else {
                    listtype = i;
                    OutputSCF++;
                    options.scf_version = 2;
                }
                break;

            case 'x':
                if (strcmp(args, "-xgr") == 0){
                    options.xgr++;
                    j = strlen(args) - 1;   /* break out of inner loop */
                    break;
                }
                else {
                    usage(argc, argv);
                    exit(2);
                }

            case '?':
            default:
                usage(argc, argv);
                exit(2);
            }
        }
    }

    if ((dev == 0) &&
        (options.shift || options.renorm || options.respace ||
         options.raw_data || options.xgr || options.multicomp ||
         options.recallndb || (options.gauss == 0) ||
         (options.tip_dir[0] != '\0')))
    {
        usage(argc, argv);
        exit(2);
    }

    if (OutputPhd || OutputQual || OutputFasta || 
        OutputQualRpt || OutputSCF || OutputFourMultiFastaFiles ||
        (options.tal_dir[0] != '\0') || (options.tip_dir[0] != '\0') ||
        (options.tab_dir[0] != '\0') || (options.hpr_dir[0] != '\0') || 
        options.mix || options.poly || options.indel_detect || 
        options.het || options.indel_resolve) 
    {
        options.process_bases = 1;
    }

    /* Make sure that data was provided */
    if ((optind == argc-1) && (argv[optind][0] == '-') &&
        (strcmp(argv[optind], "-h") != 0)
//      && (strcmp(argv[optind], "-dev") != 0)
                                              ) {
        usage(argc, argv);
        fprintf(stderr, "\nNo data was provided.\n");
        exit(2);
    }

#if 0
    fprintf(stderr, "OutputPhd=%d OutputQual=%d OutputFasta=%d OutputSCF=%d \n", 
        OutputPhd, OutputQual, OutputFasta, OutputSCF);
    fprintf(stderr, "OutputQualRpt=%d OutputAln=%d tip=%d tab=%d het=%d mix=%d poly=%d \n",
        OutputQualRpt, OutputAln, tip, tab, het, mix, poly);
#endif
    if (!OutputPhd     && !OutputQual  && !OutputFasta  && !OutputSCF && 
        !OutputQualRpt && (options.tal_dir[0] == '\0')  && !OutputFlow &&
         (options.tip_dir[0] == '\0')  && (options.tab_dir[0] == '\0') && 
         (options.hpr_dir[0] == '\0') && !options.poly && 
         !options.indel_detect && !options.indel_resolve &&
        !options.raw_data && !options.xgr && !OutputFourMultiFastaFiles) 
    {
        usage(argc, argv);
	(void)fprintf(stderr, "%s: no output type specified\n", argv[0]);
	exit(2);
    }

    if ((options.tab_dir[0] != '\0') && !options.het && !options.mix)
    {
        (void)fprintf(stderr, "Can not produce TAB file without -het or -mix option\n");
        exit(2);
    }

    if (multiqualFileName[0] != '\0') {
        unlink(multiqualFileName);
    }

    if (multiseqFileName[0] != '\0') {
        unlink(multiseqFileName);
    }

    /*
     * Set line buffering on the status output so that someone monitoring
     * progress can actually see incremental progress.
     */
    (void)setvbuf(stderr, NULL, _IOLBF, 0);

    /*
     * Table specified on command-line may refer explicitly to the '3700pop5',
     * '3700pop6', '3100' or 'mbace' built-in tables or 
     * to an external table.  
     *
     * If no table specified on command-line, try the environment.  
     *
     * If nothing there, the default table will later be selected 
     * for each file to match the PDMF string.
     * 
     * Order of preference:
     * 1) external table specified with "-t" option
     * 2) built-in table specified with one of options -3730, -3700pop5, etc.
     * 3) external table pointed to by environmental variable LOOKUP_TABLE
     * 4) built-in table determined automatically from the name of mobility file
     * 5) default built-in lookup table ABI3730-pop7
     */

    if (!options.nocall)
    {
        /* External lookup table was not passed 
         * but its type was specified through a command line option
         */
        if ((lut_name == NULL) && (options.lut_type != 0)) 
        {
            if (options.lut_type == ABI3700pop5) {
                table = Btk_get_3700pop5_table();
                if (Verbose > 1) {
                    fprintf(stderr,
                        "Using a built-in 3700 Pop-5 table for this run.\n");
                }
            }
            else if (options.lut_type == ABI3700pop6) {
                table = Btk_get_3700pop6_table();
                if (Verbose > 1) {
                    fprintf(stderr,
                             "Using a built-in 3700 Pop-6 table for this run.\n");
                }
            }
            else if (options.lut_type == ABI3730pop7) {
                table = Btk_get_3730pop7_table();
                if (Verbose > 1) {
                    fprintf(stderr,
                             "Using a built-in 3700 Pop-6 table for this run.\n");
                }
            }
            else if (options.lut_type == ABI3100) {
                table = Btk_get_3100pop6_table();
                if (Verbose > 1) {
                    fprintf(stderr,
                             "Using a built-in 3100 Pop-6 table for this run.\n");
                }
            }
            else if (options.lut_type == MegaBACE) {
                table = Btk_get_mbace_table();
                if (Verbose > 1) {
                    fprintf(stderr,
                    "Using a built-in MegaBACEtable for this run.\n");
                }
            }
        }
        else if ((lut_name == NULL) && (options.lut_type == 0))
        {
    	    lut_name = getenv("LOOKUP_TABLE");
            if (lut_name != NULL) 
            {
                if ((strstr(lut_name, "3700") != NULL) ||
                    (strstr(lut_name, "Pop5") != NULL) ||
                    (strstr(lut_name, "POP5") != NULL)) {
                    options.lut_type = ABI3700pop5;
                } else if ((strstr(lut_name, "Pop6") != NULL) &&
                    (strstr(lut_name,  "POP6") != NULL)) {
                    options.lut_type = ABI3700pop6;
                } else if (strstr(lut_name, "3100") != NULL) {
                    options.lut_type = ABI3100;
                } else if (strstr(lut_name, "3730") != NULL) {
                    options.lut_type = ABI3730pop7;
                } else if (strstr(lut_name, "BACE") != NULL ||
                           strstr(lut_name, "bace") != NULL ||
                           strstr(lut_name, "Mega") != NULL ||
                           strstr(lut_name, "mega") != NULL)
                {
                    options.lut_type = MegaBACE;
                }
                else {
                    fprintf(stderr,
                    "Type of the table stored in LOOKUP_TABLE ");
                    fprintf(stderr,
                    "can not be determined; assuming it is 3730\n");
                     options.lut_type = ABI3730pop7;
                }
            }
        }

        /* Read the named lookup table: either external or 
         * specified by LOOKUP_TABLE variable
         */
        if ((lut_name != NULL) && (table == NULL))
        {
            if ((table = Btk_read_lookup_table(lut_name)) == NULL) {
                fprintf(stderr, 
                   "Couldn't read lookup table '%s'.\n",
                      lut_name);
                exit_message(&options, 1);
            }
            else 
            {
                if (Verbose > 1) {
                    fprintf(stderr, 
                    "Using external lookuptable '%s' for this run.\n",
                        lut_name);
                }
            }
        }
        /* At this point, table may be == NULL if and only if
         * - valid table name was not passed with -t option,    
         * - table type was not specified through a command line option and
         * - valid table was not passed through a LOOKUP_TABLE variable 
         */
    }

    /* Read the context table */
    if (context_table != NULL) {
        if ((ctable = read_context_table(context_table)) == NULL) {
            fprintf(stderr, 
                "Couldn't read the context table '%s'\n", context_table);
            exit_message(&options, 1);
        }
        else {
            if (Verbose > 1) {
                fprintf(stderr, 
                "Using context table '%s' of dimension %d for this run.\n", 
                context_table, ctable->dimension);
            }
        }
    }

    if (ConsensusSpecified && ConsensusName != NULL) {
        if (read_sequence_from_fasta(ConsensusName, &ConsensusSeq, &message)
	    	== ERROR) {
	    FREE(ConsensusSeq);
	}
    }

    if (OutputFourMultiFastaFiles)
    {
        sprintf(multiseqsFileName, "%s/tt.seq",    MultiFastaFilesDirName);
        sprintf(multiqualFileName, "%s/tt.qual",   MultiFastaFilesDirName);
        sprintf(multilocsFileName, "%s/tt.pos",    MultiFastaFilesDirName);
        sprintf(multistatFileName, "%s/tt.status", MultiFastaFilesDirName);

        unlink(multiseqsFileName);
        unlink(multiqualFileName);
        unlink(multilocsFileName);
        unlink(multistatFileName);
    }

    if (optind == argc)
        fprintf(stderr, "No input data is specified\n");

    switch (InputType) {
    case NAME_FILES:
	for (i = optind; i < argc; i++) 
        {
            if (argv[i][0] == '-')
            {
                usage(argc, argv);
                fprintf(stderr, 
                    "\nCommand line option(s) or flag(s) specified after the data.");
                fprintf(stderr," They should precede the data\n\n");
                exit_message(&options, 2);
            }

	    if (stat(argv[i], &statbuf) != 0) {
		error(argv[i], "can't stat", errno);
		continue;
	    }
	    if (statbuf.st_mode & S_IFDIR) {
		error(argv[i], "is a directory!", 0);
		continue;
	    }
            strcpy(options.path, argv[i]);
	    if (process_file(table, ctable, argv[i], i, argc, ConsensusName, 
                ConsensusSeq, &options, &message) != SUCCESS) 
            {
		if (message.text[0] != '\0') {
		    fprintf(stderr, "%s: %s\n", argv[i], message.text);
		}
                fprintf(stderr,"\n");
		continue;
	    }
            fprintf(stderr,"\n");
	}
	break;

    case NAME_FILEOFFILES:
        if (file_type == SFF)
        {
            fprintf(stderr, "Can not use -if option with SFF files\n");
            exit (-1);
        }
	if (process_fileoffiles(table, ctable, InputName, ConsensusName, 
            ConsensusSeq, &options, &message) != SUCCESS) 
        {
	    if (message.text[0] != '\0') {
		fprintf(stderr, "%s: %s\n", argv[0], message.text);
	    }
	}
	break;

    case NAME_DIR:
        if (file_type == SFF)
        {
            fprintf(stderr, "Can not use -id option with SFF files\n");
            exit (-1);
        }
	process_dir(table, ctable, InputName, ConsensusName, 
        ConsensusSeq, &options, &message);
	break;

    default:
	fprintf(stderr, "%s: internal error: name type %d\n", argv[0],
			InputType);
	exit_message(&options, 1);
    }

    if (OutputQualRpt) {
      /*  Routine outputs its own error messages to stderr if necessary */
      output_qual_report(Qual_data, QualRptName);
      FREE(Qual_data);
    }

    /* Clean up. */
    if (!options.nocall)
        Btk_destroy_lookup_table(table);
    destroy_context_table(ctable);
    FREE(ConsensusSeq);

    return SUCCESS;
}
