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
 * Copyright (c) 2000 Paracel Inc.  All rights reserved.
 *
 * 2.23 2003/11/05 22:59:15
 *
 * This file takes as input one consensus file and one or more sample
 * files, aligns the sequences with the consensus, extracts the front
 * fragments of the sequence, which most likely contain the vector
 * sequence.  Then the file aligns all these fragments into a "consensus"
 * short vector sequence using a scoring algorithm.  (The details of this
 * algorithm can be found in Btk_assemble_fragments.c)
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

#ifdef OSF1
extern __const char *__const sys_errlist[];
#endif

#include "util.h"
#include "Btk_qv_data.h"
#include "Btk_qv.h"
#include "Btk_match_data.h"   
#include "Btk_qv_io.h"
#include "Btk_align_fragments.h"
#include "Btk_assemble_fragments.h"
#include "Btk_lookup_table.h"
#include "context_table.h"
#include "tracepoly.h"
#include "train.h"              /* needs util.h */
#include "Btk_compute_tpars.h"  /* needs train.h */
#include "Btk_compute_match.h"
#include "train_data.h"

#define DEBUG 0
#define BUFLEN 1000
#define FASTA_LEN 1000
#define NO_GOOD_ALIGN -20
#define BAD_PROCESSING -21
#define POSSIBLE_REPEAT -22

typedef enum {
    /* input will come from files named on the cmd line */
     INP_FILES=1,       
    /* input will come from all files in a directory */
     INP_DIR,           
    /* input will come a file with one consensus and one sample file 
     * per line
     */
     INP_PROJECTFILE,   
}  INPTYPE;

static char InputName[BUFLEN]; /* Path of dir. */
static  INPTYPE InputType;

// Name of the Consensus file.
static char ConsensusName[BUFLEN]; 
// Whether the user has specified a consensus.
static int  ConsensusSpecified;    

// Name of the initial vector file.
static char InitVectorName[BUFLEN]; 
// Whether the user has specified an intial vector.
static int  InitVectorSpecified;    

// Name of the Output file.
static char OutputName[BUFLEN];
// Whether the user has specified the output.
static int  OutputSpecified;

/* Whether the terminating gaps in the alignment are treated
   differently from the internal gaps. */
static int  GapsAreDifferent;

/* Any sequence whose score >= RepeatFraction*HighScore is 
 * considered a repeat 
 */
static double RepeatFraction;

static int nocall;            /* Whether to use ABI base calls */
static int recalln;           /* Whether to recall 'N's to the best guess */
static int edited_bases;      /* Whether to sed edited bases */
static Options options;       /* Will be passed to Btk_compute_tpars */
static char status_code[BUFLEN];

/* Scores and penalties. */
static int Match;
static int MisMatch;
static int GapInit;
static int GapExt;

/* Min and max lengths of the fragments used to assemble the short vector. */
static int MinFragLen;
static int MaxFragLen;

/* 
 *  Global variables to summarize number of input files.  These are set in
 *  process_fragment_file and output as comments at the end of the run.
 */
int Count_input_files = 0;  
int Count_file_errors = 0;  
int Count_processing_errors = 0;
int Count_possible_repeats  = 0;
int Count_no_alignments     = 0;
int Count_long_frag     = 0;
int Count_short_frag     = 0;
int Count_assembling_errors     = 0;


/*
 *  Global variables to store result contigs generated from aligning
 *  sequence with consensus -- the front portion of the sequence that
 *  might contain the vector.
 */
Contig* seqVector, finalVector;

/**
 * This function prints the usage message.
 */
static void show_usage(char** argv)
{
    fprintf(stderr, "Version: %s\n",TT_VERSION);
    fprintf(stderr, "\nusage: %s [ -o <output_file> ] \n", argv[0]);
    fprintf(stderr, "        [ -M <match> ]           (change default match premium)\n");
    fprintf(stderr, "        [ -X <mismatch> ]        (change default mismatch penalty)\n");
    fprintf(stderr, "        [ -E <gap_ext > ]        (change default gap extension penalty)\n");
    fprintf(stderr, "        [ -I <gap_init > ]       (change default gap initiation penalty)\n");
    fprintf(stderr, "        [ -l <minfraglen> ]      (change default min fragment length)\n");
    fprintf(stderr, "        [ -L <maxfraglen> ]      (change default max fragment length)\n");
    fprintf(stderr, "	[ -g <initvector_name>]  (set file name of the initial short vector sequence)\n");
    fprintf(stderr, "	[ -G ]			 (distinguish terminating gaps from internal gaps)\n");
    fprintf(stderr, "        -C <consensus_name> \n"); 
    fprintf(stderr, "        <sample_file(s)> | -d <input_dir> | -p <project_file>\n");
}

/**
 * This function uses default settings of ttuner to process the specified
 * sample file, then reads in bases into the specified Contig object.
 *
 * Its synopsis is:
 * result = read_bases(file_name, contig, message)
 * where
 *	file_name	is the address of the specified sample file name
 *	contig		is the address of the specified contig
 *	message		is the address of a BtkMessage where information
 *			about an error will be put, if any
 *	result		is 0 on success, !0 if an error occurs
 */
static int
read_bases(char* frag_name, Contig* contig, 
           Contig* consensusc, Contig* consensusrc,
           BtkMessage *message) 
{

    int             num_bases;
    char*           bases     = NULL, *chemistry = NULL, *call_method = NULL;
    int*            peak_locs = NULL, *quality_values;
    int             num_datapoints;
    int            *chromatogram[NUM_COLORS];
    double         *params[NUM_PARAMS];
    double         *iheight;
    double         *iheight2;
    double         *ave_iheight;
    int	            i, r;
    BtkLookupTable *table = NULL; 
    ContextTable   *ctable = NULL;
    ReadInfo        read_info;
  
    if (frag_name == NULL) {
        sprintf(message->text, "first argument of read_bases() is NULL");
	Count_file_errors++;
	return ERROR;
    }

    if ((r = Btk_read_sample_file(frag_name, &num_bases, &bases,
	  edited_bases, &peak_locs, &quality_values, &num_datapoints,
    	  &chromatogram[0], &chromatogram[1],
	  &chromatogram[2], &chromatogram[3],
    	  &call_method, &chemistry, status_code, options,              /* global */
          message)) != SUCCESS)
    {
      	message->text[0] = '\0';  /* message may have been set with no error */
      	Count_file_errors++;
      	return ERROR;
    }

    /*  Check that we got some called bases from the input file */
    if(num_bases <= 0) {
    	release1(bases, peak_locs, chromatogram, quality_values, NULL, NULL);
    	sprintf(message->text,"No called base information in input file");
    	Count_file_errors++;
    	return ERROR;
    }
    if(num_datapoints <= 0) {
    	release1(bases, peak_locs, chromatogram, quality_values, NULL, NULL);
    	sprintf(message->text,"No peak information in input file");
    	Count_file_errors++;
    	return ERROR;
    }

    for(i=0; i < NUM_PARAMS; i++) {
    	params[i] = CALLOC(double, num_bases);
    	if(params[i] == NULL) {
      	    release1(bases, peak_locs, chromatogram, quality_values, params, iheight);
            FREE(iheight2);
            FREE(ave_iheight);
      	    Count_processing_errors++;
      	    return ERROR;
    	}
    }
    iheight = CALLOC(double, num_bases);
    iheight2= CALLOC(double, num_bases);
    ave_iheight = CALLOC(double, num_bases);
    if(iheight==NULL) {
        release1(bases, peak_locs, chromatogram, quality_values, params, iheight);
        FREE(iheight2);
        FREE(ave_iheight);
        Count_processing_errors++;
        return ERROR;
    }

    if (nocall && recalln) {
      	fprintf(stderr, 
	    "\nPlease, specify only one of the options -nocall, -recalln\n");
      	exit(-1);
    }

    options.nocall = nocall;
    options.recalln = recalln;
    options.edited_bases = edited_bases;
    options.Verbose = 1;
    options.het = 0;
    strcpy(options.file_name, frag_name);

    if (Btk_compute_tpars(&num_bases, &bases, &peak_locs, &num_datapoints,
        chromatogram, "ACGT", &params[0], &params[1], &params[2], &params[3], 
        &iheight, &iheight2, &ave_iheight, &read_info, table, ctable, options, 
        message, NULL) != SUCCESS) {
        Count_processing_errors++;
        release1(bases, peak_locs, chromatogram, quality_values, params, iheight);
        FREE(iheight2);
        FREE(ave_iheight);
        return ERROR;
    }

    // make the contig from the computed bases
    r = contig_create(contig, bases, num_bases, quality_values, message);
    if(r==ERROR) {
    	Count_processing_errors++;
    	release1(bases, peak_locs, chromatogram, quality_values, params, iheight);
        FREE(iheight2);
        FREE(ave_iheight);
    	return ERROR;
    }

    release1(bases, peak_locs, chromatogram, quality_values, params, iheight);
    FREE(iheight2);
    FREE(ave_iheight);
    return SUCCESS;
  
}

/**
 * This function releases the resources held by the specified objects.
 */
static int
release2(Contig *frag_front, Contig *fragment, 
        Align *best_alignment, Align *vector_start, Align *vector_end,
        BtkMessage *message) {

    if (frag_front != NULL) {
    	if (contig_release( frag_front, message) == ERROR) return ERROR;
    }
    
    if (fragment != NULL) {
    	if (contig_release( fragment, message) == ERROR) return ERROR;
    }
    
    if (best_alignment != NULL) {
        if (align_release( best_alignment, message) == ERROR) return ERROR;
    }

    if (vector_start != NULL) {
        if (align_release( vector_start, message) == ERROR) return ERROR;
    }

    if (vector_end != NULL) {
      if (align_release( vector_end, message) == ERROR) return ERROR;
    }

    return SUCCESS;
}

/**
 * This function first processes the specified sample file.  Then it reads
 * the resulted bases, aligns it to the specified consensus, then takes the
 * front piece of the sample config outside of the best alignment, assembles
 * this piece with the specified WeightedBases.
 *
 * Its synopsis is:
 * result = process_fragment_file(align_params, align_paramsIUB,
 *				  file_name, fout,
 *				  consensus, consensusrc,
 *				  result_vector, message)
 * where
 *	align_params	is the address of an Align_params data structure
 *			for aligning consensus with sample
 *	align_paramsIUB	is the address of an Align_params data structure
 *			for aligning vector with sample  (Note: no vector
 *			is provided for assembling short vector.)
 *	file_name	is the address of the specified sample file name.
 *	fout		is the (already open) output stream
 *	consensus	is the address of the consensus contig
 *	consensusrc	is the address of the reverse complement contig
 *			of the consensus
 *	result_vector	is the address of the WeightedBases object holding
 *			the result of the assembled short vector sequence
 *	message		is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *      result          is 0 on success, !0 if an error occurs
 */
static int process_fragment_file(Align_params *alp, 
				 Align_params *alpIUB,
				 char* frag_name, 
				 FILE *fout, 
				 Contig* consensusc, 
				 Contig* consensusrc,
				 WeightedBases* result_vector,
				 BtkMessage* message) 
{
    int            r,i, count_sub, count_del, count_ins;
    int            num_align;
    Contig         fragment;
    Contig 	   frag_front;
    Range          align_range;
    Range          clear_range;
    Align          best_alignment;
    Align          vector_start;
    Align          vector_end;
    int            alignment_size;
 
    num_align = 0;
    contig_init(&fragment, message);
    contig_init(&frag_front, message);
    align_init(&best_alignment, message);
    align_init(&vector_start, message);
    align_init(&vector_end, message);

    if(frag_name[0] == 0) {
        Count_file_errors++;
        goto error;
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "File: %s\n", frag_name);
    fprintf(fout, "# File: %s\n", frag_name);

    Count_input_files++;

    // processes the sample file and reads the bases 
    if ((r = read_bases(frag_name, &fragment, 
                        consensusc, consensusrc,                       
                        message)) != SUCCESS) {
        message->text[0] = '\0';  /* message may have been set with no error */
        Count_processing_errors++;
        goto error;
    }

    /* aligns the bases with the consensus to obtain the front piece of the
     * sample bases which will be used for short vector assembling
     */
    r = Btk_compute_match(alp, alpIUB, consensusc, consensusrc, &fragment,
	 		  &num_align, &align_range, RepeatFraction, 
			  &best_alignment,
			  NULL, &vector_start, &vector_end, 
			  &clear_range, 0, message);
    if(r == ERROR) {
        fprintf(stderr,"Error processing input file: %s\n", frag_name);
        fprintf(fout,"# Error processing input file: %s\n", frag_name);
        goto error;
    }
    if (num_align == 0) {

        fprintf(fout, "# No good alignments - ignore file %s\n",
	        frag_name);
        fprintf(stderr, "No good alignments - ignore file %s\n",
	        frag_name);
        Count_no_alignments++;
        goto error;
    } else if ( (num_align != 1) || (best_alignment.score == -1) ) {

        fprintf(fout, "# Possible repeat !!! - ignore file %s\n",
	        frag_name);
        fprintf(stderr, "Possible repeat !!! - ignore file %s\n",
	        frag_name);
        Count_possible_repeats++;
        goto error;
    } else {

        /* Count number of mismatches */
        count_sub = 0;
        count_del = 0;
        count_ins = 0; 
        alignment_size=best_alignment.trace_len;
        for (i=best_alignment.trace_len/3; 
             i<best_alignment.trace_len/3*2; i++) {
            if (best_alignment.trace_dir[i] == '2') {
	        count_del++;
            } else if (best_alignment.trace_dir[i] == '1') {
	        count_ins++;
            } else if (best_alignment.trace_dir[i] == ' ') {
	        count_sub++;
            }
        }

        /* Append the results to output file only if number of
         * mismatches is < 10%
	 */
	if ((count_sub+count_ins+count_del) * 10 >= alignment_size/3) {
	    fprintf(fout, "# BAD PROCESSING!!! - ignore file %s\n",
	            frag_name);
            fprintf(stderr, "BAD PROCESSING!!! - ignore file %s\n",
	            frag_name);
            Count_processing_errors++;
            goto error;
        }
	fprintf(stderr, "Fraction of errors: sub=%f, ins=%f, del=%f tot=%f\n",
		(double)count_sub/(double)alignment_size,
		(double)count_ins/(double)alignment_size,
		(double)count_del/(double)alignment_size,
		(double)(count_sub+count_ins+count_del)
		    /(double)alignment_size);
	fprintf(stderr, "Align range length:%4d\n", alignment_size);
	fprintf(stderr, "Align range begin:%4d\n", align_range.begin);
	fprintf(stderr, "Align range end:%4d\n", align_range.end);
	fprintf(stderr, "Number of matched bases: %d\n", 
		alignment_size-count_sub-count_ins-count_del);

	fprintf(stderr, "Extracting the front piece...\n");
  
	/* extract the front piece (from fragment beginning to best alignment
         * beginning) from the fragment 
	 */
	r = contig_create(&frag_front, fragment.sequence,
			      align_range.begin, fragment.qv, message);
	if (r == ERROR) {
	    fprintf(stderr, "Error extracting frag_front: %s\n", 
		    message->text);
	    goto error;
	}

	if (frag_front.length > MaxFragLen) {
	    fprintf(stderr,"Frag_front too long!\n");
	    fprintf(fout,"# Frag_front too long!\n");
	    Count_long_frag++;
	    goto error;
	} else if (frag_front.length < MinFragLen) {
	    fprintf(stderr,"Frag_front too short!\n");
	    fprintf(fout,"# Frag_front too short!\n");
	    Count_short_frag++;
	    goto error;
	}

	// assemble the frag_front with the result vector
	frag_front.sequence[frag_front.length] = '\0';
	fprintf(stderr, "frag_front sequence = %s, len = %d\n",
		frag_front.sequence, frag_front.length);
	r = assemble_fragments(result_vector, &frag_front, 
			       GapsAreDifferent, fout, message);
	if (r == ERROR) {
	    fprintf(stderr,"Error assembling frag_front!\n");
	    fprintf(fout,"# Error assembling frag_front!\n");
	    Count_assembling_errors++;
	    goto error;
	}

	/* output result_vector to fout */
	r = weightedbases_print(result_vector, fout, message);
	if (r == ERROR) {
	    fprintf(stderr, "%s\n", message->text);
	    goto error;
	}
    
	/* extract the end piece (from best alignment end to fragment end)
	   from the fragment */
	/*
	if (frag_end != NULL) {
	    r = contig_create(frag_end, 
				  fragment.sequence + align_range.end + 1,
				  fragment.length - (align_range.end + 1),
				  message);
	    if (r == ERROR) {
	        fprintf(stderr, "Error extracting frag_end: %s\n", 
			message->text);
		goto error;
	    }
	}
	*/

    }

    release2(&frag_front, &fragment, &best_alignment, 
	     &vector_start, &vector_end, message);
    return SUCCESS;

 error:
    release2(&frag_front, &fragment, &best_alignment, 
	     &vector_start, &vector_end, message);
    return ERROR;
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
 * result = process_dir(align_pars, align_pars_IUB, fout, cons, rev_comp, dir,
		 	result_vector, message)
 *
 * where
 *      align_pars      is the address of an Align_params data structure for
 *                      aligning consensus with sample
 *      align_pars_IUB  is the address of an Align_params data structure for
 *                      aligning vector with sample
 *      fout            is the (already open) output stream
 *      cons            is the address of the consensus Contig
 *      rev_comp        is the address of the reverse complement Contig
 *      dir             is the name (path) of a directory
 *	result_vector	is the address of the WeightedBases object holding
 *			the result of the assembled short vector sequence
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *      result          is 0 on success, !0 if an error occurs
 */
static int 
process_dir(Align_params *ap,
		       Align_params *apIUB,
		       FILE* fout,
		       Contig* consensus,
		       Contig* rev_comp,
		       char *dir,
		       WeightedBases* result_vector,
		       BtkMessage *message)
{
    DIR *d;
    struct dirent *de;
    char path[MAXPATHLEN];
    struct stat statbuf;

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

	process_fragment_file(ap, apIUB, path, fout, consensus,
			      rev_comp, 
			      result_vector, message);
    }

    (void)closedir(d);
    return SUCCESS;
}


/*
 * This function processes all files contained in the specified project file.
 * Its synopsis is:
 *
 * result = process_projectfile(align_pars, align_pars_IUB, fout, projectFile,
             			result_vector, message)
 *
 * where
 *      align_pars      is the address of an Align_params data structure for
 *                      aligning consensus with sample
 *      align_pars_IUB  is the address of an Align_params data structure for
 *                      aligning vector with sample
 *      fout            is the (already open) output stream
 *      projectFile     is the name (path) of a project file
 *	result_vector	is the address of the WeightedBases object holding
 *			the result of the assembled short vector sequence
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *      result          is 0 on success, !0 if an error occurs
 */
static int
process_projectfile(
                    Align_params *ap,
                    Align_params *apIUB,
                    FILE* fout,
                    char *projectFile,
                    WeightedBases* result_vector,
                    BtkMessage *message)
{
    char    buffer1[BUFLEN];
    char    buffer2[BUFLEN];
    char    *current, *previous, *tmp;
    Contig  consensus;
    Contig  rev_comp;
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
	      r=contig_release( &rev_comp, message );
	      if ( r==ERROR ) { return ERROR; }
	      have_consensus=0;
	  }

	  /* Read consensus file. */
	  r=local_read_fasta( ConsensusName, &consensus, message );
	  if ( r==ERROR ) { return ERROR; }

	  /* Create reverse complement of consensus. */
	  r=contig_get_reverse_comp( &rev_comp, &consensus, message );
	  if ( r==ERROR ) { return ERROR; }
	  have_consensus=1;
      }

      tmp=previous;
      previous=current;
      current=tmp;

      fprintf(stderr, "Reference = %s\n", ConsensusName );
      fprintf(fout, "# Reference = %s\n", ConsensusName );

      process_fragment_file(ap, apIUB, fragmentName, fout, &consensus,
			    &rev_comp, result_vector, message);
    }

    fclose(infile);

    /* Clean up. */
    if ( have_consensus ) {
        r=contig_release( &consensus, message );
	if ( r==ERROR ) { return ERROR; }
	r=contig_release( &rev_comp, message );
	if ( r==ERROR ) { return ERROR; }
    }

    return SUCCESS;
}


/* ------------------------------------------------------------------- */
int 
main(int argc, char** argv)
{
    extern char    *optarg;
    extern int      optind;
    FILE           *fout;
    Align_params    align_pars, align_pars_IUB;
    int             r, i;
    BtkMessage      message;
    Contig          consensus, rev_comp;
    Contig          init_vector; // first guess of the short vector sequence.
    WeightedBases   result_vector;

    if (argc == 1) {
        show_usage(argv);
	exit(2);
    }

    /* Set defaults. */
    ConsensusName[0]    = '\0';
    ConsensusSpecified  = 0;
    InitVectorName[0]   = '\0';
    InitVectorSpecified = 0;
    OutputName[0]       = '\0';
    OutputSpecified     = 0;
    InputName[0]        = '\0';
    InputType           =  INP_FILES;
    GapsAreDifferent	= 0;
    RepeatFraction      = 0.85;
    Match               = 10;
    MisMatch            = -20;
    GapInit             = -40;
    GapExt              = -30;
    MinFragLen	        = 5;
    MaxFragLen	        = 50;
    nocall              = 0;
    recalln             = 0;
    edited_bases        = 0;


    while ( (i = getopt(argc, argv, "M:X:I:E:l:L:C:o:d:p:g:G")) != EOF ) {
        switch(i) {
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
	case 'I':
	    GapInit=atoi(optarg);
	    if ( GapInit > 0 ) {
	        fprintf(stderr, "Error: GapInit   must not be positive.\n");
		exit(2);
	    }
	    break;
	case 'E':
	    GapExt=atoi(optarg);
	    if ( GapExt > 0 ) {
	        fprintf(stderr, "Error: GapExt must not be positive.\n");
		exit(2);
	    }
	    break;
	case 'l':
	    MinFragLen=atoi(optarg);
	    if ( MinFragLen < 0 ) {
	        fprintf(stderr, "Error: MinFragLen must not be negative.\n");
		exit(2);
	    }
	    break;
	case 'L':
	    MaxFragLen=atoi(optarg);
	    if ( MaxFragLen <= MinFragLen ) {
	        fprintf(stderr, 
		      "Error: MaxFragLen must be greater than MinFragLen.\n");
	        exit(2);
	    }
	    if ( MaxFragLen > 200 ) {
	        fprintf(stderr, 
			"Error: MaxFragLen must not be greater than 200.\n");
		exit(2);
	    }
	    break;
	case 'C':
	    ConsensusSpecified++;
	    (void)strncpy(ConsensusName, optarg, sizeof(ConsensusName));
	    break;
	case 'o':
	    OutputSpecified++;
	    (void)strncpy(OutputName, optarg, sizeof(OutputName));
	    break;
	case 'd':
	    InputType =  INP_DIR;
	    (void)strncpy(InputName, optarg, sizeof(InputName));
	    break;
	case 'p':
	    InputType =  INP_PROJECTFILE;
	    (void)strncpy(InputName, optarg, sizeof(InputName));
	    break;
	case 'g':
	    InitVectorSpecified++;
	    (void)strncpy(InitVectorName, optarg, sizeof(InitVectorName));
	    break;
	case 'G':
	    GapsAreDifferent++;
	    break;
	default:
	    show_usage(argv);
	    exit(2);
	}
    }

    /*
     * Set line buffering on the status output so that someone monitoring
     * progress can actually see incremental progress.
     */
    (void)setbuf(stderr, NULL);

    if (OutputSpecified) {
        /* Set up output file. */
        if ((fout = fopen(OutputName, "w")) == NULL) {
	  fprintf(stderr, "Cannot open output file '%s'\n", OutputName);
	  exit(ERROR);
	}
	fprintf(stderr, "Software Version: %s\n", TT_VERSION);
    } else {
        fprintf(stderr, "No output file specified.  Using stdout.\n");
	fout=stdout;
    }

    /* Set parameters for alignment with consensus */   
    r = set_alignment_parameters(&align_pars, Match, MisMatch, GapInit,
			       GapExt, &message );
    if ( r == ERROR ) {
        fprintf(stderr, message.text);
	exit (ERROR);
    }

    /* Set parameters for alignment with vector */  
    r = set_alignment_parameters_IUB(&align_pars_IUB, 
        Match * VECTOR_MATCH_MULTIPLIER, MisMatch, 
        GapInit, GapExt, &message );
    if ( r == ERROR ) {
        fprintf(stderr, message.text);
	exit (ERROR);
    }

    /* Print parameter values to output file. */
    fprintf( fout, "# Version   = %s\n", TT_VERSION);
    fprintf( fout, "# Match     = %d\n", Match );
    fprintf( fout, "# MisMatch  = %d\n", MisMatch );
    fprintf( fout, "# GapInit   = %d\n", GapInit );
    fprintf( fout, "# GapExt    = %d\n", GapExt );
    fprintf( fout, "# RepeatFraction = %f\n", RepeatFraction );
    fprintf( fout, "# MinFragLen  = %d\n", MinFragLen );
    fprintf( fout, "# MaxFragLen  = %d\n", MaxFragLen );
    fprintf( fout, "#\n");

    result_vector.length = -1;
    result_vector.max_length = MAX_ASSEMBLED_VECTOR_LEN;
    result_vector.bases = CALLOC(WeightedBase, result_vector.max_length);
    if (result_vector.bases == NULL) {
        fprintf(stderr, "Insufficient memory at file=%s, line=%d\n",
		__FILE__, __LINE__);
	exit(ERROR);
    }

    if (InitVectorSpecified) {
        contig_init(&init_vector, &message);

        if (local_read_fasta(InitVectorName, &init_vector, &message)
		== ERROR) {
	    exit(ERROR);
	}
	
	r = weightedbases_create(&result_vector, &init_vector, &message);
	if ( r == ERROR) {
	    fprintf(stderr, 
		    "Error when creating the initial weighted bases: %s\n",
		    message.text);
	    exit(ERROR);
	}
	if (&init_vector != NULL) {
	    contig_release(&init_vector, &message);
	}
	fprintf( fout, "# Initial Vector = %s\n", InitVectorName );
    } else {
        fprintf( fout, "# No Initial Vector defined.\n");
    }

    switch(InputType) {
    case  INP_FILES:
        if (!ConsensusSpecified) {
	    fprintf(stderr, "%s: no consensus file specified.\n",argv[0]);
	    show_usage(argv);
	    exit(2);
	}
	fprintf( fout, "# Reference = %s\n", ConsensusName );
	fprintf( fout, "#\n");

	/* Read consensus file. */
	if(local_read_fasta( ConsensusName, &consensus, &message) == ERROR) {
	    fprintf(stderr, message.text);
	    exit (ERROR);
	}

	/* Create reverse complement of consensus. */
	r = contig_get_reverse_comp( &rev_comp, &consensus, &message );
	if ( r == ERROR ) {
	    fprintf(stderr, message.text);
	    exit (ERROR);
	}

	/* Loop for each fragment: */
	for(i = optind; i < argc; i++) {
	    r = process_fragment_file(&align_pars, &align_pars_IUB, argv[i],
				      fout, &consensus, &rev_comp,
				      &result_vector, &message);
	}
	break;

    case  INP_DIR:

        if (!ConsensusSpecified) {
	    fprintf(stderr, "%s: no consensus file specified.\n",argv[0]);
	    show_usage(argv);
	    exit(2);
	}
	fprintf( fout, "# Reference = %s\n", ConsensusName );
	fprintf( fout, "#\n");

	/* Read consensus file. */
	if(local_read_fasta( ConsensusName, &consensus, &message) == ERROR) {
	    fprintf(stderr, "Error reading consensus file: %s\n",message.text);
	    exit (ERROR);
	}

	/* Create reverse complement of consensus. */
	if(contig_get_reverse_comp( &rev_comp, &consensus, &message) 
	        == ERROR) {
	    fprintf(stderr, "Error getting reverse compliment: %s\n",
		    message.text);
	    exit(ERROR);
	}

	if(process_dir(&align_pars, &align_pars_IUB, fout, &consensus,
		       &rev_comp, InputName,
		       &result_vector, &message)   != SUCCESS) {
	    fprintf(stderr,"Error processing input files: %s\n", InputName);
	    exit (ERROR);
	}
	break;

    case  INP_PROJECTFILE:
        if (process_projectfile(&align_pars, &align_pars_IUB, fout,
			       InputName, 
			       &result_vector, &message) != SUCCESS) {
	    fprintf(stderr, "Error processing project file: %s\n", InputName);
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

    if (weightedbases_release( &result_vector, &message) == ERROR) {
        fprintf(stderr, "Error releasing memory: %s\n",message.text);
	exit(ERROR);
    }

    if ((InputType ==  INP_FILES ) || ( InputType ==  INP_DIR ) ) {
        if (contig_release( &consensus, &message)) {
	    fprintf(stderr, "Error releasing memory: %s\n",message.text);
	    exit (ERROR);
	}

	if (contig_release( &rev_comp, &message) == ERROR) {
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
	     Count_possible_repeats + Count_no_alignments +
	     Count_long_frag + Count_short_frag));
    fprintf(fout,"#  %6d files rejected with file errors\n", 
	    Count_file_errors);
    fprintf(fout,"#  %6d files rejected with processing errors\n",
	    Count_processing_errors);
    fprintf(fout,"#  %6d files rejected with possible repeats\n",
	    Count_possible_repeats);
    fprintf(fout,"#  %6d files rejected with no alignments\n",
	    Count_no_alignments);
    fprintf(fout,"#  %6d files rejected with too long front fragment\n",
	    Count_long_frag);
    fprintf(fout,"#  %6d files rejected with too short front fragment\n",
	    Count_short_frag);
    fprintf(fout,"#  %6d files rejected with assembling errors\n",
	    Count_assembling_errors);
    fclose(fout);
    
    exit(SUCCESS);
}



