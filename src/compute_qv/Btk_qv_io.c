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

/**	Copyright (c) 1999 Paracel Inc.  All rights reserved.
 **/

/*
 *  $Id: Btk_qv_io.c,v 1.21 2009/01/04 16:46:21 gdenisov Exp $
 */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#ifdef __WIN32
#else
#include <unistd.h>
#endif
#include <math.h>
#ifdef __WIN32
#define MAXPATHLEN      (255)
#else
#include <sys/param.h>
#endif
#include <sys/types.h>
#include <errno.h>
#include <float.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>

#include "ABI_Toolkit.h"
#include "SCF_Toolkit.h"
#include "SFF_Toolkit.h"
#include "FileHandler.h"
#include "Btk_qv.h"
#include "Btk_qv_data.h"
#include "util.h"
#include "Btk_match_data.h" 
#include "Btk_qv_io.h"
#include "Btk_qv_funs.h"
#include "Btk_compute_match.h"
#include "train.h"
#include "context_table.h"
#include "tracepoly.h"
#include "Btk_call_bases.h"
#include "Btk_scf_data.h"

#define DEBUG                        0
#define BAD_PROCESSING_MULTIPLIER    0.5
#define FASTA_LEN                 1000
#define MAX_FILE_NAME_LENGTH      1024
#define BTK_FASTA_WIDTH		   (50)
#define ADJUST_SECONDARY_PEAK        0
#define SHOW_SUBSTITUTIONS           0
#define USE_DEFAULT_CHEMISTRY        0

unsigned int max_colordata_value;

void
exit_message(Options *op, int errlevel)
{
    fprintf(stderr,
        "TraceTuner exiting abruptly with errorlevel %d while processing trace %s\n",
        errlevel, op->file_name);
    exit(errlevel);
}

/********************************************************************************
 * This function prints as formatted an error message as it can.  Its
 * synopsis is:
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
 ********************************************************************************
 */
static void
error(char *label, char *msg, int err)
{
    if (err > 0)
        fprintf(stderr, "%s: %s: %d: %s\n", label, msg, err, strerror(err));
    else if (err == 0)
        fprintf(stderr, "%s: %s\n", label, msg);
    else
        fprintf(stderr, "%s: %s: %d\n", label, msg, err);
}


/* ------------------------------------------------------------------- */
/* Reads a single FastA formatted entry from a file.
 * FastA format: lines beginning with '>' signify comments
 *               other lines signify data.
 */
int
read_sequence_from_fasta(char* file_name, char** sequence, BtkMessage *message)
{
    FILE            *infile;
    char            buffer[BUFLEN];
    int             curr_len, in_data, line_len;
    int             r;
    int             length;

    if ((infile=fopen(file_name,"r"))== NULL) {
        sprintf(message->text, "Unable to open FastA file '%s'\n", file_name);
        return ERROR;
    }

    length =  FASTA_LEN;
    *sequence = (char*) malloc( sizeof(char)*FASTA_LEN );
    MEM_ERROR( *sequence );

    curr_len=0;
    *sequence[curr_len] = '\0';
    in_data = -1;

    while ((fgets(buffer, BUFLEN, infile)) != NULL) {
        if ( buffer[0] == '>') {
        /* Header line. */
            if ( in_data==1 ) {
                /* We've read one sequence, and we hit the next header line. */
                break;
            }
            in_data=0;
        } else {
            if ( in_data == -1 ) {
                fprintf(stderr, "Warning: sequence did not contain a ");
                fprintf(stderr, "valid '>' header line! - file = %s\n",
                        file_name);
            }
            in_data=1;

            /* Take care of carriage returns. */
            line_len  = strlen(buffer);
            if (buffer[line_len-1] == '\n')  {
                buffer[line_len-1] = '\0';
                line_len--;
            }

            curr_len += line_len;

            /* Allocate more memory, if necessary. */
            if (curr_len >= length) {
                length += FASTA_LEN;
                *sequence = (char*) realloc((void*)*sequence, sizeof(char)*(length));
                MEM_ERROR( *sequence );
            }

            strncpy(&((*sequence)[curr_len-line_len]),buffer,line_len+1);
        }
    }

    r=SUCCESS;
    goto cleanup;

error:
    r=ERROR;

cleanup:
    fclose(infile);

    return r;
}

/* ------------------------------------------------------------------- */
/* Reads a single FastA formatted entry from a file.
 * FastA format: lines beginning with '>' signify comments
 *               other lines signify data.
 */
int
read_fasta(char* file_name, Contig* contig, BtkMessage *message)
{
    char* sequence;
    int r, length;
    uint8_t *quality_values = NULL;
    r = read_sequence_from_fasta(file_name, &sequence, message);
    if (r == SUCCESS) {
	length = strlen(sequence);
        quality_values = CALLOC(uint8_t, length);
	contig_create( contig, sequence, length, quality_values, message );
    }
    FREE( sequence );
    FREE( quality_values );
    return r;
}

/*******************************************************************************
 * Function: read_abi_nums
 * Purpose: read the number of called bases and the number of
 *          datapoints in each trace from sample file using the ABI Toolkit
 *******************************************************************************
 */
static int
read_abi_nums(int *num_called_bases, int use_edited_bases, int *num_datapoints,
Options *options)
{
    ABIError   r;
    int        j;
    long       num_bases;
    long       num_peaks;
    long       num_points;
    short      dye_number;

    if (options->inp_phd == 0) {
        /* Read the number of called bases */
        if (use_edited_bases) {
            if ((r = ABI_NumEditedBases(&num_bases)) != kNoError) {
                return r;
            }
        }
        else {
            if ((r = ABI_NumCalledBases(&num_bases)) != kNoError) {
                return r;
            }
        }
       *num_called_bases = num_bases;

        /* Read the number of called peak locations */
        if (use_edited_bases) {
            if ((r = ABI_NumEditedPeakLocations(&num_peaks)) != kNoError) {
                return r;
            }
        }
        else {
            if ((r = ABI_NumCalledPeakLocations(&num_peaks)) != kNoError) {
                return r;
            }
        }
        if (num_peaks != num_bases) {
            (void)fprintf(stderr, "Number of bases != number of peaks\n");
            return ERROR;
        }
    }

    /* Choose any color number */
    *num_datapoints = INF;
    for (j = 0; j < NUM_COLORS; j++) {
        dye_number = j + 1;

        /* Read the number of data points for this color */
        r = ABI_NumAnalyzedData(0, dye_number, &num_points);
        if (r != kNoError) {
            return r;
        }
        if (*num_datapoints > num_points) {
           *num_datapoints = num_points;
        }
    }

    return SUCCESS;
}

/********************************************************************************
 * Function: read_scf_nums
 ********************************************************************************
 */
static int
read_scf_nums(int *num_called_bases, int *num_datapoints, Options *options)
{
     long num_bases;
     long num_points;

     if (options->inp_phd == 0) {
         SCF_NumBases(&num_bases);

        /*  We skip the cross-checking of number of bases vs. number
         *  of peak locations that the corresponding ABI routine does,
         *  since in ABI chromats, called base arrays and called peak
         *  arrays are stored as separate data blobs, whereas in SCF files
         *  each called base is bound together with its peak location.
         */
       *num_called_bases = num_bases;
    }

    SCF_NumAnalyzedData(&num_points);
   *num_datapoints = num_points;

    return SUCCESS;
}

/****************************************************************************
 * Function: read_abi_bases_locs_and_quality_values
 * Purpose: read the array of ABI called bases and array of their locations
 *          from sample file using the ABI Toolkit
 ****************************************************************************
 */
static int
read_abi_bases_locs_and_quality_values(
    int   num_bases,
    char *called_bases,
    int   use_edited_bases,
    int  *called_locs,
    uint8_t  *quality_values,
    BtkMessage *message)
{
    short     *peak_locs;
    short      num_orig_qvs=0;
    char      *orig_qv;
    int        i;
    ABIError   r;

    /* Read the called or edited bases */
    if (use_edited_bases) {
        if ((r = ABI_EditedBases(called_bases)) != kNoError) {
            return r;
        }
    }
    else {
        if ((r = ABI_CalledBases(called_bases)) != kNoError) {
            return r;
        }
    }

    /* Create an auxiliary array of called peak locations */
    peak_locs = CALLOC(short, num_bases);
    MEM_ERROR(peak_locs);

    /* Read the called or edited peak locations
     * and FREE the auxiliary array
     */
    if (use_edited_bases) {
        if ((r = ABI_EditedPeakLocations(peak_locs)) != kNoError) {
            goto error;
        }
    }
    else {
        if ((r = ABI_CalledPeakLocations(peak_locs)) != kNoError) {
            goto error;
        }
    }
    for (i = 0; i < num_bases; i++) {
        called_locs[i] = peak_locs[i];
    }

    /* Read original quality values */
    num_orig_qvs = 0;
    r = kNoError;
    if (((r = ABI_NumQualityValues(&num_orig_qvs)) != kNoError) ||
        ((r == kNoError) && (num_orig_qvs == 0)))
    {
        num_orig_qvs = 0;
    }

    if ((int)num_orig_qvs == num_bases)
    {
        orig_qv = CALLOC(char, num_bases);
        if ((r = ABI_BasecallerQualityValues(orig_qv)) != kNoError) {
            goto error;
        }
        for (i = 0; i < num_bases; i++) {
            quality_values[i] = (uint8_t)orig_qv[i];
        }
#if 0
        {
           int i;
           fprintf(stderr, "orig_qvs=\n");
           for (i=0; i<num_bases; i++) {
              fprintf(stderr, "%2d ", orig_qv[i]);
              if ((i>0) && (i%30 == 0))
                  fprintf(stderr, "\n");
           }
           fprintf(stderr, "\n");
        }
#endif
        FREE(orig_qv);
    }

    FREE(peak_locs);
    return SUCCESS;

error:
    FREE(peak_locs);
    return r;
}


/********************************************************************************
 * Function: read_scf_bases_and_locs
 ********************************************************************************
 */
static int
read_scf_bases_and_locs(
    int   num_bases,
    char *called_bases,
    int  *called_locs,
    BtkMessage *message)
{
     short     *peak_locs;
     int        i;

     SCF_Bases(called_bases);

     peak_locs = CALLOC(short, num_bases);
     MEM_ERROR(peak_locs);

     SCF_PeakLocations(peak_locs);

     for (i = 0; i < num_bases; i++)
          called_locs[i] = peak_locs[i];

     FREE(peak_locs);

     return SUCCESS;

 error:
     FREE(peak_locs);

     return kMemoryFull;
}

/******************************************************************************
 * Function: read_consensus_from_sample_file
 ******************************************************************************
 */
int
read_consensus_from_sample_file(char **consensus_seq, int verbose)
{
    ABIError   r;
    long cons_len;

    /* Read the number of bases in consensus sequence*/
    if ((r = ABI_NumConsensusBases(&cons_len)) != kNoError) {
        return r;
    }

#if 0
        fprintf(stderr, "Cons_length = %ld\n", cons_len);
#endif

    *consensus_seq = CALLOC(char, cons_len + 1);
    if ((r = ABI_ConsensusBases(*consensus_seq)) != kNoError) {
        return r;
    }

    (*consensus_seq)[cons_len] = '\0';
    return r;
}

/****************************************************************************
 * Function: read_color_data
 * Purpose: read the data array for each trace and the color2base array
 *          which shows the base corresponding to each trace
 ****************************************************************************
 */
static int
read_abi_color_data(
    int   num_values,
    int **chromatogram,
    char *color2base,
    BtkMessage *message)
{
    int  *analyzed_data;
    short dye_number;
    char  base;
    int   i, j;
    ABIError r;

    analyzed_data = CALLOC(int, num_values);
    MEM_ERROR(analyzed_data);

    for (j = 0; j < NUM_COLORS; j++) {
        dye_number = j + 1;

        /* Which base corresponds to the selected dye number? */
        if ((r = ABI_DyeIndexToBase(dye_number, &base)) != kNoError) {
            FREE(analyzed_data);
            return r;
        }
        color2base[j] = base;

        /* Read the chromatograms and free the auxiliary array */
        r = ABI_AnalyzedData(0, dye_number, analyzed_data);
        if (r != kNoError) {
            FREE(analyzed_data);
            return r;
        }
        for (i = 0; i < num_values; i++) {
            chromatogram[j][i] = analyzed_data[i];
        }
    }

    FREE(analyzed_data);
    return SUCCESS;
error:
    FREE(analyzed_data);
    return ERROR;
}

/********************************************************************************
 * Function: read_scf_color_data
 ********************************************************************************
 */
static int
read_scf_color_data(int num_values,
                    int **chromatogram,
                    char *color2base,
                    BtkMessage *message)
{
     int  *analyzed_data;
     short dye_number;
     int   i;

     analyzed_data = CALLOC(int, num_values);
     MEM_ERROR(analyzed_data);

     color2base[0] = 'A';
     color2base[1] = 'C';
     color2base[2] = 'G';
     color2base[3] = 'T';

     for (dye_number = 0; dye_number < NUM_COLORS; dye_number++)
     {
          SCF_AnalyzedData(dye_number, analyzed_data);

          for (i = 0; i < num_values; i++)
               chromatogram[dye_number][i] = analyzed_data[i];
     }

     FREE(analyzed_data);

     return SUCCESS;

 error:
     return kMemoryFull;
}

/*
 * This function extracts base calls and chromatogram traces from a single
 * sample file.  Its synopsis is:
 *
 * result = Btk_read_sample_file(file_name, num_bases, called_bases,
 *	called_locs, num_values, avals, cvals,
 *	gvals, tvals, call_method, chemistry,
 *	options, message, verbose)
 *
 * where
 *	file_name	is the name (path) of the sample file
 *	num_bases	is the address where the number of bases will be put
 *	called_bases	is the address where a pointer to the array of called
 *			bases will be put
 *      edited_bases      integer which indicates whether to use the array of
 *                      called (0) or edited bases (1) from sample file
 *	called_locs	is the address where a pointer to the array of called
 *			base locations will be put
 *	num_values	is the address where the number of trace points will
 *			be put
 *	avals		is the address where a pointer to the array of A trace
 *			points will be put
 *	cvals		is the address where a pointer to the array of C trace
 *			points will be put
 *	gvals		is the address where a pointer to the array of G trace
 *			points will be put
 *	tvals		is the address where a pointer to the array of T trace
 *			points will be put
 *	call_method	is an optional address where a string that describes
 *			the method used to call the bases will be put
 *	chemistry	is an optional address where a string that describes
 *			the chemistry used (primer or terminator) will be put
 *	message		is the address of a BtkMessage where information about
 *			an error will be put, if any
 *	verbose		is whether to write status messages to stderr, and
 *			how verbosely
 *
 *	result		is 0 on success, !0 if an error occurs
 *
 * The caller is responsible for free()ing the arrays returned; the convenience
 * function Btk_release_file_data() can be used for this purpose.
 */
int
Btk_read_sample_file(
    char       *file_name,
    int        *num_bases,
    char      **called_bases,
    int         use_edited_bases,
    int       **called_locs,
    uint8_t   **quality_values,
    int        *num_values,
    int       **avals,
    int       **cvals,
    int       **gvals,
    int       **tvals,
    char      **call_method,
    char      **chemistry,
    char       *status_code,
    int        *fileType,
    Options     options,
    BtkMessage *message)
{
    char *seq_name;
    int   i, r=0;
    int  *chromatogram[NUM_COLORS] = {NULL, NULL, NULL, NULL};
    void *p;
    long  fileSize, n;
    char  color2base[5];
    unsigned char magic[2];
    FILE *fp;
    char  tempFileName[MAX_FILE_NAME_LENGTH] = "";
    char  command[2048];
    char  phd_file_name[1000];

    for (i = 0; i < NUM_COLORS; i++) {
        chromatogram[i] = NULL;
    }
   *called_bases = NULL;
   *called_locs  = NULL;

    /* Set seq_name to just the name of the file, no leading path. */

#ifdef __WIN32
    if ((seq_name = strrchr(file_name, '\\')) != NULL) {
#else
    if ((seq_name = strrchr(file_name, '/')) != NULL) {
#endif
        seq_name++;
    }
    else {
        seq_name = file_name;
    }

    if ((fp = fopen(file_name, "r")) == NULL)
    {
        sprintf(message->text,
                "could not open file for reading, errno = %d", errno);
        return ERROR;
    }
    fread(magic, 1, 2, fp);
    fclose(fp);

#ifdef __WIN32
    if (magic[0] == 0x1f && magic[1] == 0x8b)   /* gzipped file */
    {
        sprintf(tempFileName, "TTUNERQQ.TMP");

        sprintf(command, "gzip -d -f -c %s > TTUNERQQ.TMP", file_name);

        if (system(command) != 0)
        {
            sprintf(message->text,
                    "gunzip command failed on file %s, errno = %d", file_name,
                    errno);
            return ERROR;
        }
    }
    else if (magic[0] == 0x1f && magic[1] == 0x9d)   /* UNIX compressed file */
    {
        if (seq_name[strlen(seq_name) - 2] != '.' ||
            seq_name[strlen(seq_name) - 1] != 'Z')
        {
            sprintf(message->text,
                    "cannot uncompress file %s, name must end in .Z",
                    file_name);
            return ERROR;
        }

        sprintf(tempFileName, "TTUNERQQ.TMP");

        sprintf(command, "gzip -d -f -c %s > TTUNERQQ.TMP", file_name);

        if (system(command) != 0)
        {
            sprintf(message->text,
                    "gunzip command failed on file %s, errno = %d", file_name,
                    errno);
            return ERROR;
        }
    }
#else
    if (magic[0] == 0x1f && magic[1] == 0x8b)   /* gzipped file */
    {
        sprintf(tempFileName, "/tmp/%s", seq_name);
        if (seq_name[strlen(seq_name) - 3] == '.' &&
            seq_name[strlen(seq_name) - 2] == 'g' &&
            seq_name[strlen(seq_name) - 1] == 'z')
            tempFileName[strlen(tempFileName) - 3] = '\0';

        sprintf(command, "gunzip -f -c %s > %s", file_name, tempFileName);
        if (system(command) != 0)
        {
            sprintf(message->text,
                    "gunzip command failed on file %s, errno = %d", file_name,
                    errno);
            return ERROR;
        }
    }
    else if (magic[0] == 0x1f && magic[1] == 0x9d)   /* UNIX compressed file */
    {
        if (seq_name[strlen(seq_name) - 2] != '.' ||
            seq_name[strlen(seq_name) - 1] != 'Z')
        {
            sprintf(message->text,
                    "cannot uncompress file %s, name must end in .Z",
                    file_name);
            return ERROR;
        }

        sprintf(tempFileName, "/tmp/%s", seq_name);
                tempFileName[strlen(tempFileName) - 2] = '\0';

        sprintf(command, "uncompress -f -c %s > %s", file_name, tempFileName);
        if (system(command) != 0)
        {
            sprintf(message->text,
                    "gunzip command failed on file %s, errno = %d", file_name,
                    errno);
            return ERROR;
        }
    }
#endif

    if ((r = F_Open(tempFileName[0] != '\0' ? tempFileName : file_name,
                    &p, &fileSize, fileType)) != kNoError) {
        sprintf(message->text, "Error opening file: %s", 
                  ABI_ErrorString((ABIError)r));
        goto error;
    }

    if (tempFileName[0] != '\0')
        unlink(tempFileName);

    if (*fileType == ABI)
    {
         if ((r = read_abi_nums(num_bases, use_edited_bases, num_values,
             &options)) != kNoError) {
              strcpy(status_code, "ABIFILE_FAILURE");
              sprintf(message->text, "Error reading file: %s",
                      ABI_ErrorString((ABIError)r));
              goto error;
         }
    }
    else if (*fileType == SCF)
    {
         if ((r = read_scf_nums(num_bases, num_values, &options))
             != kNoError) {
             strcpy(status_code, "ABIFILE_FAILURE");
             sprintf(message->text, "Error reading file: %s",
                 ABI_ErrorString((ABIError)r));
             goto error;
         }
    }
    else {
        return kWrongFileType;
    }

    if ((options.Verbose > 1) && (*fileType == ABI && use_edited_bases) &&
        (options.inp_phd == 0)) {
        fprintf(stderr, "Starting from edited bases\n");
    }

    if ((options.inp_phd == 0) && (num_bases > 0)) {
        /* Allocate memory for data arrays */
       *called_bases = CALLOC(char, *num_bases);
        MEM_ERROR(*called_bases);
       *called_locs = CALLOC(int, *num_bases);
        MEM_ERROR(*called_locs);
       *quality_values = CALLOC(uint8_t, *num_bases);
        MEM_ERROR(*quality_values);
    }

    for (i = 0; i < NUM_COLORS; i++) {
         chromatogram[i] = CALLOC(int, *num_values * 3);
         MEM_ERROR(chromatogram);
    }

    if ((*fileType == ABI) && (options.inp_phd == 0))
    {
         if ((r = read_abi_bases_locs_and_quality_values(*num_bases,
                *called_bases, use_edited_bases, *called_locs,
                *quality_values, message))
             != SUCCESS)
         {
              strcpy(status_code, "ABIFILE_FAILURE");
              sprintf(message->text, "Error reading file: %s",
                      ABI_ErrorString((ABIError)r));
              goto error;
         }
    }
    else if ((*fileType == SCF) && (options.inp_phd == 0))
    {
         if ((r = read_scf_bases_and_locs(*num_bases, *called_bases,
            *called_locs, message))
             != SUCCESS)
         {
              strcpy(status_code, "ABIFILE_FAILURE");
              sprintf(message->text, "Error reading file: %s",
                      ABI_ErrorString((ABIError)r));
              goto error;
         }
    }


    if (*num_bases > MAX_NUM_BASES)
    {
       *num_bases = MAX_NUM_BASES;
    }

    if (*fileType == ABI)
    {
         if ((r = read_abi_color_data(*num_values, chromatogram, color2base,
                                      message))
             != SUCCESS)
         {
             strcpy(status_code, "ABIFILE_FAILURE");
             sprintf(message->text, "Error reading file: %s",
                 ABI_ErrorString((ABIError)r));
             goto error;
         }
    }
    else if (*fileType == SCF)
    {
         if ((r = read_scf_color_data(*num_values, chromatogram, color2base,
                                      message))
             != SUCCESS)
         {
              strcpy(status_code, "ABIFILE_FAILURE");
              sprintf(message->text, "Error reading file: %s",
                      ABI_ErrorString((ABIError)r));
              goto error;
         }
    }

    if (options.Verbose > 2) {
         (void)fprintf(stderr, "Base order: %c%c%c%c\n", color2base[0],
            color2base[1], color2base[2], color2base[3]);
    }

    if (call_method != NULL)
    {
        *call_method = CALLOC(char, BTKMESSAGE_LENGTH);
         MEM_ERROR(*call_method);

         /*
          *  If the analysis version isn't available, just leave the string
          *  null and it will get printed out as unknown later on.
          */
         if (*fileType == ABI)
         {
              if ((r = ABI_AnalysisVersion(BTKMESSAGE_LENGTH,
                                           *call_method, &n))
                  != kNoError)
              {
                   FREE(*call_method);
              } else {
                   (*call_method)[n] = '\0';
              }
         }
    }

    if (*chemistry != NULL) {
        *chemistry = CALLOC(char, BTKMESSAGE_LENGTH);
         MEM_ERROR(*chemistry);

         /*
          *  If the chemistry isn't available, just leave the string
          *  null and it will get printed out as unknown later on.
          */
         if (*fileType == ABI)
         {
              if ((r = ABI_MobilityFile(BTKMESSAGE_LENGTH, *chemistry, &n))
                  != kNoError)
              {
                   FREE(*chemistry);
              } else {
                   (*chemistry)[n] = '\0';
              }
         }
    }

    F_Close(p, *fileType);
    p = NULL;

    /* If -ipd <dir> option is used, read original bases and
     * locations from phd file, rather than from sample file
     */
    if (options.inp_phd) {
        strcpy(tempFileName, seq_name);
        if (seq_name[strlen(seq_name) - 3] == '.' &&
            seq_name[strlen(seq_name) - 2] == 'g' &&
            seq_name[strlen(seq_name) - 1] == 'z') {
            tempFileName[strlen(tempFileName) - 3] = '\0';
        }
        if (seq_name[strlen(seq_name) - 2] == '.' &&
            seq_name[strlen(seq_name) - 1] == 'Z') {
            tempFileName[strlen(tempFileName) - 2] = '\0';
        }
#ifdef __WIN32
        sprintf(phd_file_name, "%s\\%s.phd.1",
            options.inp_phd_dir, tempFileName);
#else
        sprintf(phd_file_name, "%s/%s.phd.1",
            options.inp_phd_dir, tempFileName);
#endif
       *num_bases = -1;
        if ((*num_bases = get_phd_num_bases(phd_file_name, message)) < 0)
        {
            fprintf(stderr,"%s", "PHREDFILE_FAILURE");
            sprintf(status_code, "%s", "PHREDFILE_FAILURE");
            sprintf(message->text, "Could not read phd file: %s\n",
                phd_file_name);
            FREE(*call_method);
            return ERROR;
        }

        if (*num_bases > MAX_NUM_BASES)
        {
            *num_bases = MAX_NUM_BASES;
        }

        if (*num_bases > 0) {
           *called_bases = CALLOC(char, *num_bases);
            MEM_ERROR(called_bases);

           *called_locs = CALLOC(int, *num_bases);
            MEM_ERROR(called_locs);

           *quality_values = CALLOC(uint8_t, *num_bases);
            MEM_ERROR(quality_values);

            if (Btk_read_phd_file(phd_file_name, *called_bases, *quality_values,
               *called_locs, num_bases, message) != SUCCESS)
            {
                sprintf(status_code, "%s", "PHREDFILE_FAILURE");
                FREE(called_bases);
                FREE(called_locs);
                FREE(quality_values);
                FREE(*call_method);
                sprintf(message->text, "Could not read phd file: %s\n",
                    phd_file_name);
                return ERROR;
            }
#if 0
            fprintf(stderr, "options.Verbose=%d options.inp_phd =%d\n",
               options.Verbose, options.inp_phd);
#endif
            if ((options.Verbose > 1) && (options.inp_phd > 0))
            fprintf(stderr,
               "Reading original bases and locations from phd file: %s\n",
               phd_file_name);
        }
    }

    /* Set the correct input vars to the correct colors */
    switch(color2base[0]) {
    case 'A':
         *avals = chromatogram[0];
         break;
    case 'C':
         *cvals = chromatogram[0];
         break;
    case 'G':
         *gvals = chromatogram[0];
         break;
    case 'T':
         *tvals = chromatogram[0];
         break;
    default:
         (void)sprintf(message->text, "unknown color '%c'", color2base[0]);
         r = ERROR;
         goto error;
    }
    switch(color2base[1]) {
    case 'A':
         *avals = chromatogram[1];
         break;
    case 'C':
         *cvals = chromatogram[1];
         break;
    case 'G':
         *gvals = chromatogram[1];
         break;
    case 'T':
         *tvals = chromatogram[1];
         break;
    default:
         (void)sprintf(message->text, "unknown color '%c'", color2base[1]);
         r = ERROR;
         goto error;
    }
    switch(color2base[2]) {
    case 'A':
         *avals = chromatogram[2];
         break;
    case 'C':
         *cvals = chromatogram[2];
         break;
    case 'G':
         *gvals = chromatogram[2];
         break;
    case 'T':
         *tvals = chromatogram[2];
         break;
    default:
         (void)sprintf(message->text, "unknown color '%c'", color2base[2]);
         r = ERROR;
         goto error;
    }
    switch(color2base[3]) {
    case 'A':
         *avals = chromatogram[3];
         break;
    case 'C':
         *cvals = chromatogram[3];
         break;
    case 'G':
         *gvals = chromatogram[3];
         break;
    case 'T':
         *tvals = chromatogram[3];
         break;
    default:
         (void)sprintf(message->text, "unknown color '%c'", color2base[3]);
         r = ERROR;
         goto error;
    }

    return SUCCESS;

    error:
    if (chromatogram[0] != NULL) {
         Btk_release_file_data(*called_bases, *called_locs, *quality_values,
             chromatogram, call_method, (chemistry != NULL) ? chemistry : NULL);
    }
    if (p != NULL) {
         F_Close(p, *fileType);
         p = NULL;
    }
    return r;
}

/*
 * This function frees up all the arrays allocated by Btk_read_sample_file().
 * Its synopsis is:
 *
 * Btk_release_file_data(called_bases, called_locs, chromatogram)
 *
 * where
 *	called_bases	is the pointer to called bases originally set by
 *			Btk_read_sample_file()
 *	called_locs	is the pointer to called base locations originally set
 *			by Btk_read_sample_file()
 *	chromatogram	is the array of pointers to color-specific trace points
 *			originally set by Btk_read_sample_file()
 *	call_method	is the string originally set by Btk_read_sample_file(),
 *			if any
 *	chemistry	is the string originally set by Btk_read_sample_file(),
 *			if any
 */
void
Btk_release_file_data(
    char  *called_bases,
    int   *called_locs,
    uint8_t   *quality_values,
    int  **chromatogram,
    char **call_method,
    char **chemistry)
{
    int i;


    FREE(called_bases);
    FREE(called_locs);
    FREE(quality_values);
    for (i = 0; i < NUM_COLORS; i++) {
	FREE(chromatogram[i]);
    }
    FREE(*call_method);
    FREE(*chemistry);
}

/*******************************************************************************
 * Function: find_trim_points
 *******************************************************************************
 */
int
find_trim_points(int num_called_bases, uint8_t *quality_values, int trim_window,
    float trim_threshold, int *left_trim_point, int *right_trim_point)
{
    int i, j, qualifying_range_found = 0;
    float sum, average;

    for (i = 0; i <= num_called_bases - trim_window; i++)
    {
        for (j = i, sum = 0; j < trim_window + i; j++)
            sum += quality_values[j];

        average = sum / trim_window;

        if (average >= trim_threshold)
        {
           *left_trim_point = i;
            qualifying_range_found = 1;
            break;
        }
    }

    if (!qualifying_range_found)
    {
       *left_trim_point = *right_trim_point = -1;
        return 0;
    }

    for (i = num_called_bases - 1; i >= trim_window - 1; i--)
    {
        for (j = i, sum = 0; j > i - trim_window; j--)
            sum += quality_values[j];

        average = sum / trim_window;

        if (average >= trim_threshold)
        {
           *right_trim_point = i;
            break;
        }
    }

    if (*right_trim_point - *left_trim_point + 1 < 20)
    {
       *left_trim_point = *right_trim_point = -1;
        return 0;
    }
    else
        return (*right_trim_point - *left_trim_point + 1);
}


/********************************************************************************
 * This function writes out the quality values (only) as a .qual file.  Its
 * synopsis is:
 *
 * success = Btk_output_quality_values(file_name, path, quality_values,
 *     num_values, left_trim_point, right_trim_point, verbose)
 *
 * where
 *	file_name	is the name of the sample file
 *	path		is the path name of the directory in which to write
 *			the .qual file, if any (NULL means current dir)
 *	quality_values	is an array of quality values
 *	num_values	is the number of elements in quality_values
 *	verbose		is whether to write status messages to stderr, and
 *			how verbosely
 *
 *	success		is SUCCESS or ERROR
 ********************************************************************************
 */
int
Btk_output_quality_values(
    int QualType,
    char *file_name,
    char *path,
    char *multiqualFileName,
    uint8_t *quality_values,
    int num_values,
    int left_trim_point,
    int right_trim_point,
    int verbose)
{
    int i;
    char *seq_name, qual_file_name[MAXPATHLEN];
    FILE *qv_out = NULL, *dir_out = NULL, *multi_out = NULL;

    /* Use the name of the sample file, sans path, as the sequence name */
#ifdef __WIN32
    if ((seq_name = strrchr(file_name, '\\')) != NULL) {
#else
    if ((seq_name = strrchr(file_name, '/')) != NULL) {
#endif
        seq_name++;
    }
    else {
        seq_name = file_name;
    }

    if (QualType & NAME_DIR) {
#ifdef __WIN32
        sprintf(qual_file_name, "%s\\%s.qual", path, seq_name);
#else
        sprintf(qual_file_name, "%s/%s.qual", path, seq_name);
#endif
        if ((dir_out = fopen(qual_file_name, "w")) == NULL) {
            error(qual_file_name, "couldn't open", errno);
            return ERROR;
        }
    }

    if (QualType & NAME_FILES) {
        sprintf(qual_file_name, "%s.qual", seq_name);

        if ((qv_out = fopen(qual_file_name, "w")) == NULL) {
            error(qual_file_name, "couldn't open", errno);
            return ERROR;
        }
    }

    if (QualType & NAME_MULTI) {
        if ((multi_out = fopen(multiqualFileName, "a")) == NULL) {
            error(multiqualFileName, "couldn't open", errno);
            fclose(qv_out);
            return ERROR;
        }
    }

    /* Output a header into the quality file.  Fields are taken following
     * phred's example: filename, number of quality values, number of bases
     * trimmed, number after trimming, and sample file source. */
    if (qv_out) {
        fprintf(qv_out, ">%s %d %d %d\n", seq_name, num_values, left_trim_point,
                right_trim_point - left_trim_point + 1);
    }

    if (dir_out) {
        fprintf(dir_out, ">%s %d %d %d\n", seq_name, num_values, left_trim_point,
                right_trim_point - left_trim_point + 1);
    }

    if (multi_out) {
        fprintf(multi_out, ">%s %d %d %d\n", seq_name, num_values, left_trim_point,
                right_trim_point - left_trim_point + 1);
    }

    for (i = 0; i < num_values; i++) {
        if (qv_out) {
            fprintf(qv_out, " %2d", quality_values[i]);

            if ((i % 17 == 16) || (i == num_values - 1)) {
                putc('\n', qv_out);
            }
        }

        if (dir_out) {
            fprintf(dir_out, " %2d", quality_values[i]);

            if ((i % 17 == 16) || (i == num_values - 1)) {
                putc('\n', dir_out);
            }
        }

        if (multi_out) {
            fprintf(multi_out, " %2d", quality_values[i]);

            if ((i % 17 == 16) || (i == num_values - 1)) {
                putc('\n', multi_out);
            }
        }
    }

    if (qv_out) {
        fclose(qv_out);
    }

    if (dir_out) {
        fclose(dir_out);
    }

    if (multi_out) {
        fclose(multi_out);
    }

    return SUCCESS;
}

/**
 * This function releases the resources held by the specified objects.
 */
static int
release(Contig *consensus, Contig *consensusrc, Contig *fragment, 
        Align *best_alignment, Align *vector_start, Align *vector_end,
	Align_params *align_pars, Align_params *align_pars_IUB,
	BtkMessage *message) {

    if (consensus != NULL) {
        if (contig_release( consensus, message) == ERROR) return ERROR;
    }

    if (consensusrc != NULL) {
        if (contig_release( consensusrc, message) == ERROR) return ERROR;
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

    if (align_pars != NULL) {
	if (alignment_parameters_release(align_pars, message) == ERROR)
            return ERROR;
    }

    if (align_pars_IUB != NULL) {
    	if (alignment_parameters_release(align_pars_IUB, message) == ERROR)
            return ERROR;
    }

    return SUCCESS;
}

	/********************************************************************************
	 * This function writes out a FASTQ file with the (potentially) recalled
	 * bases.  Its synopsis is:
	 *
	 * success = Btk_output_fastq_file(file_name, path, called_bases, quality_values,
	 *                  num_bases, verbose)
	 *
	 * where
	 *	file_name	is the name of the sample file
	 *	path		is the path name of the directory in which to write
	 *			the .seq file, if any (NULL means current dir)
	 *	called_bases	is an array of base calls
	 *	quality_values	is an array of quality values
	 *	num_bases	is the number of elements in called_bases
	 *	verbose		is whether to write status messages to stderr, and
	 *			how verbosely
	 *
	 *	success		is SUCCESS or ERROR
	 ********************************************************************************
	 */
	int
	Btk_output_fastq_file(
						  int FastqType,
						  char *file_name,
						  char *path,
						  char *multifastqFileName,
						  char *called_bases,
						  uint8_t *quality_values,
						  int num_bases,
						  int left_trim_point,
						  int right_trim_point,
						  int verbose)
	{
		int i;
		char qchar;
		char *seq_name, fastq_file_name[MAXPATHLEN];
		char sequence_string[num_bases + 1];
		char quality_string[num_bases + 1];
		FILE *fastq_out = NULL, *dir_out = NULL, *multi_out = NULL;

		/* Use the name of the sample file, sans path, as the sequence name. */
#ifdef __WIN32
		if ((seq_name = strrchr(file_name, '\\')) != NULL) {
#else
			if ((seq_name = strrchr(file_name, '/')) != NULL) {
#endif
				seq_name++;
			}
			else {
				seq_name = file_name;
			}

			/* Ignore any trailing N characters in called_bases */
			strncpy( sequence_string, called_bases, num_bases );
			sequence_string[num_bases] = 0;
			
			for(i=0;i<num_bases;i++)
			{
				/* Quality values are already integers */
				/* Write a Sanger FASTQ file, PHRED scores with ASCII offset 33 */
				qchar = 33 + quality_values[i];
				//assert( (33 <= qchar) && (qchar <= 126) );
				quality_string[i] = qchar;
				//fprintf(stderr, "%c, %d\n", quality_string[i], quality_values[i]);
			}
			quality_string[num_bases] = 0; /* Explicit null terminator */
			//fprintf(stderr, "Sequence string:\n%d, %d, %d, %s\n", (int)num_bases, (int)strlen(sequence_string), (int)strlen(quality_string), sequence_string); /* enable to print sequence string to stderr */
			assert( num_bases == strlen(quality_string) );
			assert( num_bases == strlen(sequence_string) );
			assert( strlen(quality_string) == strlen(sequence_string) );

			if (FastqType & NAME_DIR) {
#ifdef __WIN32
				sprintf(fastq_file_name, "%s\\%s.fastq", path, seq_name);
#else
				sprintf(fastq_file_name, "%s/%s.fastq", path, seq_name);
#endif
				if ((dir_out = fopen(fastq_file_name, "w")) == NULL) {
					error(fastq_file_name, "couldn't open", errno);
					return ERROR;
				}
			}

			if (FastqType & NAME_FILES) {
				sprintf(fastq_file_name, "%s.fastq", seq_name);
				if ((fastq_out = fopen(fastq_file_name, "w")) == NULL) {
					error(fastq_file_name, "couldn't open", errno);
					return ERROR;
				}
			}

			if (FastqType & NAME_MULTI) {
				if ((multi_out = fopen(multifastqFileName, "a")) == NULL) {
					error(multifastqFileName, "couldn't open", errno);
					fclose(fastq_out);
					return ERROR;
				}
			}

			if (fastq_out) {
				fprintf(fastq_out, "@%s %d %d %d\n%s\n+\n%s\n", seq_name, num_bases,
						left_trim_point, right_trim_point - left_trim_point + 1,
						sequence_string, quality_string);
			}

			if (dir_out) {
				fprintf(dir_out, "@%s %d %d %d\n%s\n+\n%s\n", seq_name, num_bases,
						left_trim_point, right_trim_point - left_trim_point + 1,
						sequence_string, quality_string);
			}

			if (multi_out) {
				fprintf(multi_out, "@%s %d %d %d\n%s\n+\n%s\n", seq_name, num_bases,
						left_trim_point, right_trim_point - left_trim_point + 1,
						sequence_string, quality_string);
			}

			if (fastq_out) {
				fclose(fastq_out);
			}

			if (dir_out) {
				fclose(dir_out);
			}

			if (multi_out) {
				fclose(multi_out);
			}

			return SUCCESS;
		}

/********************************************************************************
 * This function writes out a ".tal" file.
 * Its synopsis is:
 *
 * success = Btk_output_tal_file(file_name, path,
 *       consensus_name, consensus_seq, called_bases, num_bases,
 *	 Match, MisMatch, Insertion, Deletion, RepeatFraction, verbose);
 *
 * where
 *	file_name	is the name of the sample file
 *	path		is the path name of the directory in which to write
 *			the .qual file, if any (NULL means current dir)
 *	consensus_name	is the name of consensus file
 *	consensus_seq	is an array of consensus sequence
 *	called_bases	is an array of base calls
 *	num_called_bases  is the number of elements in the called_bases
 *	Match		the Match premium used for the alignment
 *	MisMatch	the MisMatch penalty used for the alignment
 *	Insertion	the Insertion penalty used for the alignment
 *	Deletion	the Deletion penalty used for the alignment
 *	RepeatFraction	the RepeatFraction parameter used for the alignment
 *	verbose		is whether to write status messages to stderr, and
 *			how verbosely
 *
 *	success		is SUCCESS or ERROR
 ********************************************************************************
 */
int
Btk_output_tal_file(
    char *file_name,
    char *path,
    char *consensus_name,
    char *consensus_seq,
    char *called_bases,
    int   num_called_bases,
    int  Match,
    int  MisMatch,
    int  Insertion,
    int  Deletion,
    float  RepeatFraction,
    int  verbose)
{
    char *seq_name, tal_file_name[MAXPATHLEN];
    FILE *tal_out;
    int  alignment_size, count_del, count_ins, count_sub, ismatch, i;
    Align_params    align_pars, align_pars_IUB;
    BtkMessage message;
    Contig          consensus;
    Contig          consensusrc;
    Contig          fragment;
    Range           align_range;
    Range           clear_range;
    Align           best_alignment;
    Align           vector_start;
    Align           vector_end;
    int             num_align;
    uint8_t        *qv;

    if (consensus_seq == NULL) {
        fprintf(stderr, 
                "No consensus sequence specified; .tal file not produced\n");
        return SUCCESS;
    }

    /* Use the name of the sample file, sans path, as the sequence name */
#ifdef __WIN32
    if ((seq_name = strrchr(file_name, '\\')) != NULL) {
#else
    if ((seq_name = strrchr(file_name, '/')) != NULL) {
#endif
	seq_name++;
    }
    else {
	seq_name = file_name;
    }

    /*
     * Create the name of phd file by appending suffix ".tal" to the sample
     * file name.
     */
    if (path != NULL) {
#ifdef __WIN32
	(void)sprintf(tal_file_name, "%s\\%s.tal", path, seq_name);
#else
	(void)sprintf(tal_file_name, "%s/%s.tal", path, seq_name);
#endif
    }
    else {
	/* put it in the current dir */
	(void)sprintf(tal_file_name, "%s.tal", seq_name);
    }

    if ((tal_out = fopen(tal_file_name, "w")) == NULL) {
	error(tal_file_name, "couldn't open", errno);
	return ERROR;
    }

    if (verbose > 2) {
	(void)fprintf(stderr, "Computing alignment \n");
    }

    /* Set alignment parameters globally. */
    if (set_alignment_parameters(&align_pars, Match, MisMatch, Insertion,
            Deletion, &message ) == ERROR) {
        fprintf(stderr, message.text);
        goto error;
    }

    if (set_alignment_parameters_IUB(&align_pars_IUB, Match, MisMatch,
            Insertion, Deletion, &message ) == ERROR ) {
        fprintf(stderr, message.text);
        goto error;
    }	

    align_init(&best_alignment, &message);
    align_init(&vector_start, &message);
    align_init(&vector_end, &message);    
    qv = CALLOC(uint8_t, strlen(consensus_seq));

    if (contig_create(&consensus, consensus_seq, 
		      strlen(consensus_seq), qv, &message)
		== ERROR) {
	fprintf(stderr, "can't create contig\n");
	goto error;
    }

    /* Create reverse complement of consensus. */
    if (contig_get_reverse_comp(&consensusrc, &consensus, 
             &message)==ERROR ) {
        fprintf(stderr, message.text);
        goto error;
    }
    FREE(qv); 

    qv = CALLOC(uint8_t, num_called_bases);
    if (contig_create(&fragment, called_bases, num_called_bases, qv, &message)
		== ERROR) {
	fprintf(stderr, "can't create contig\n");
	goto error;
     }
    FREE(qv);

    if (Btk_compute_match(&align_pars, &align_pars_IUB, &consensus,
	      &consensusrc, &fragment, &num_align,
	      &align_range, RepeatFraction,
	      &best_alignment, NULL, &vector_start,
	      &vector_end, &clear_range, 0, &message) == ERROR) {
	    goto error; 
	}


    if (verbose > 2) {
	(void)fprintf(stderr, "Writing .tal output to %s\n",
			tal_file_name);
    }
    /* Create header */
    (void)fprintf(tal_out, "#CHROMAT_FILE: %s\n", seq_name);
    if (consensus_name != NULL) {
	(void)fprintf(tal_out, "#CONSENSUS_FILE: %s\n", consensus_name);
    } else {
	(void)fprintf(tal_out, "#CONSENSUS_FILE: \n");
    }
    (void)fprintf(tal_out, "#SOFTWARE_VERSION: %s\n", TT_VERSION);
    (void)fprintf(tal_out, "#\n");
    (void)fprintf(tal_out, "#Match = %d\n",Match);
    (void)fprintf(tal_out, "#MisMatch = %d\n",MisMatch);
    (void)fprintf(tal_out, "#Insertion = %d\n",Insertion);
    (void)fprintf(tal_out, "#Deletion = %d\n",Deletion);
    (void)fprintf(tal_out, "#RepeatFraction = %f\n",RepeatFraction);

    if (num_align == 0) {
        fprintf(tal_out, "#ALIGN_STATUS: No good alignments\n");
    } else if ((num_align != 1) || (best_alignment.score == -1)) {
        fprintf(tal_out, "#ALIGN_STATUS: Possible repeats\n");
    } else {
        /* Count number of mismatches */
        count_sub = 0;
        count_del = 0;
        count_ins = 0;
        alignment_size=best_alignment.trace_len;
        for (i=0; i<best_alignment.trace_len; i++) {
            if (best_alignment.trace_dir[i] == '2') {
                count_del++;
            } else if (best_alignment.trace_dir[i] == '1') {
                count_ins++;
            } else if (best_alignment.trace_dir[i] == ' ') {
                count_sub++;
            }
        }

        /* Append the results to tal file only if the fraction of
         * mismatches is < BAD_PROCESSING_MULTIPLIER
         */
        if ((double)(count_sub+count_ins+count_del) <
                (double)alignment_size * BAD_PROCESSING_MULTIPLIER) {
            fprintf(tal_out, "#ALIGN_STATUS: OK\n");
            fprintf(tal_out, "#ALIGN_LENGTH: %d\n", alignment_size);
            fprintf(tal_out, "#%s\t\t%s\t\t%s\n", "Sample", "Cons", "is_match");
	    for (i=0; i<best_alignment.trace_len; i++) {
	        fprintf(tal_out, "%d\t%c\t%d\t%c",
			best_alignment.trace_qpos[i]+1,
			best_alignment.trace_qchar[i],
			best_alignment.trace_dpos[i]+1,
			best_alignment.trace_dchar[i]);
		ismatch = (best_alignment.trace_mchar[i]=='|')?1:0;
		fprintf(tal_out, "\t%d\n", ismatch);
#if SHOW_SUBSTITUTIONS
            if (best_alignment.trace_qchar[i] != best_alignment.trace_dchar[i] &&
                best_alignment.trace_qchar[i] != '-' &&
                best_alignment.trace_dchar[i] != '-')
            {
                fprintf(stderr, "Substitution %c/%c at pos= %d\n",
                    best_alignment.trace_dchar[i],  best_alignment.trace_qchar[i],
                    best_alignment.trace_qpos[i]);
            }
#endif
	    }
	} else {
            fprintf(tal_out, "#ALIGN_STATUS: No good alignments\n");
	}
    }

    (void)fclose(tal_out);
    release(&consensus, &consensusrc, &fragment, &best_alignment,
	    &vector_start, &vector_end,
	    &align_pars, &align_pars_IUB, &message);
    return SUCCESS;

error:
    (void)fclose(tal_out);
    release(&consensus, &consensusrc, &fragment, &best_alignment,
	    &vector_start, &vector_end,
	    &align_pars, &align_pars_IUB, &message);
    return ERROR;
}

static double
average_spacing(int base_ind, int *locs)
{
    int j=0, k;

    k = base_ind;
    if (k==0)
        return 12.;

    while ((k>0) && (j<20)) {
        k--;
        if (locs[k] < locs[k+1])
            j++;
    }
    if (j==0)
        return 12.;

    return (double)(locs[base_ind] - locs[k])/(double)j;
}

int
Btk_output_hpr_file(
    char *file_name,
    char *path,
    char *called_bases,
    int *called_locs,
    uint8_t *quality_values,
    int num_bases,
    int num_datapoints,
    int verbose)
{
    int   i, left_pos = 0, right_pos = -1, first_base_pos = -1;
    char *seq_name, hpr_file_name[MAXPATHLEN];
    FILE *hpr_out;
    float hpr_length, spacing;
    char  base = ' ';
    int   num_homopolymer_bases = 0;

    /* Use the name of the sample file, sans path, as the sequence name */
#ifdef __WIN32
    if ((seq_name = strrchr(file_name, '\\')) != NULL) {
#else
    if ((seq_name = strrchr(file_name, '/')) != NULL) {
#endif
        seq_name++;
    }
    else {
        seq_name = file_name;
    }

    /*
     * Create the name of hpr file by appending suffix ".hpr" to the sample
     * file name.
     */
    if (path != NULL) {
#ifdef __WIN32
        (void)sprintf(hpr_file_name, "%s\\%s.hpr", path, seq_name);
#else
        (void)sprintf(hpr_file_name, "%s/%s.hpr", path, seq_name);
#endif
    }
    else {
        /* put it in the current dir */
        (void)sprintf(hpr_file_name, "%s.hpr", seq_name);
    }
    if ((hpr_out = fopen(hpr_file_name, "w")) == NULL) {
        error(hpr_file_name, "couldn't open", errno);
        return ERROR;
    }
    if (verbose > 2) {
        (void)fprintf(stderr, "Writing a homopolymer runs output to %s\n",
                        hpr_file_name);
    }

    /* Create header */
    (void)fprintf(hpr_out, ">%s\n\n", seq_name);

    /* Output the called bases, their quality values, and locations */
    for (i = 0; i < num_bases; i++) {

        if (base != called_bases[i])  // new homopolymer run
        {
            if (base != ' ')
            {
                right_pos = (called_locs[i] + called_locs[i-1])/2;
                spacing = average_spacing(i, called_locs);
                hpr_length = (float)(right_pos - left_pos)/spacing;
                (void)fprintf(hpr_out, "%c %d %4.2f %d\n", base,
                        num_homopolymer_bases, hpr_length, first_base_pos);
                left_pos = right_pos;
            }
            base = called_bases[i];
            first_base_pos = called_locs[i];
            num_homopolymer_bases = 1;
        }
        else
            num_homopolymer_bases++;
    }

    (void)fclose(hpr_out);
    return SUCCESS;
}

/********************************************************************************
 * This function writes out a phred-compatible file, of the '.phd.1' variety.
 * Its synopsis is:
 *
 * success = Btk_output_phd_file(file_name, path, called_bases, called_locs,
 *					quality_values, num_bases,
 *					num_datapoints, call_method, verbose)
 *
 * where
 *	file_name	is the name of the sample file
 *	path		is the path name of the directory in which to write
 *			the .qual file, if any (NULL means current dir)
 *	called_bases	is an array of base calls
 *	called_locs	is an array of positions where called_bases were called
 *	quality_values	is an array of quality values
 *	num_bases	is the number of elements in the called_bases,
 *			called_locs, and quality_values arrays
 *	num_datapoints	is the number of data points in the read (generally
 *			10-20x the number of called bases)
 *	nocall          If true, base called are from original abi basecaller
 *                      otherwise they are from ttuner
 *	chemistry	is an optional string that describes the chemistry
 *			used (primer or terminator)
 *	verbose		is whether to write status messages to stderr, and
 *			how verbosely
 *
 *	success		is SUCCESS or ERROR
 ********************************************************************************
 */
int
Btk_output_phd_file(
    char *file_name,
    char *path,
    char *called_bases, 
    int *called_locs,
    uint8_t *quality_values,
    int num_bases,
    int num_datapoints,
    int nocall,
    char *chemistry,
    int leftTrim,
    int rightTrim,
    float trim_threshold,
    int verbose)
{
    char *seq_name, phd_file_name[MAXPATHLEN];
    FILE *phd_out;
    int qv_max, i;
    time_t current;

    /* Use the name of the sample file, sans path, as the sequence name */
#ifdef __WIN32
    if ((seq_name = strrchr(file_name, '\\')) != NULL) {
#else
    if ((seq_name = strrchr(file_name, '/')) != NULL) {
#endif
	seq_name++;
    }
    else {
	seq_name = file_name;
    }

    /*
     * Create the name of phd file by appending suffix ".phd.1" to the sample
     * file name.
     */
    if (path != NULL) {
#ifdef __WIN32
	(void)sprintf(phd_file_name, "%s\\%s.phd.1", path, seq_name);
#else
	(void)sprintf(phd_file_name, "%s/%s.phd.1", path, seq_name);
#endif
    }
    else {
	/* put it in the current dir */
	(void)sprintf(phd_file_name, "%s.phd.1", seq_name);
    }

    if ((phd_out = fopen(phd_file_name, "w")) == NULL) {
	error(phd_file_name, "couldn't open", errno);
	return ERROR;
    }
    if (verbose > 2) {
	(void)fprintf(stderr, "Writing phred-compatible output to %s\n",
			phd_file_name);
    }

    /* Create header */
    (void)fprintf(phd_out, "BEGIN_SEQUENCE %s\n\n", seq_name);
    (void)fprintf(phd_out, "BEGIN_COMMENT\n\n");
    (void)fprintf(phd_out, "CHROMAT_FILE: %s\n", seq_name);
    (void)fprintf(phd_out, "ABI_THUMBPRINT: 0\n");
    (void)fprintf(phd_out, "PHRED_VERSION: %s\n", TT_VERSION);
    if(nocall) {
	(void)fprintf(phd_out, "CALL_METHOD: abi\n");
    }
    else {
	(void)fprintf(phd_out, "CALL_METHOD: ttuner\n");
    }
    qv_max = 0;
    for (i = 0; i < num_bases; i++) {
	if (qv_max < quality_values[i]) {
	    qv_max = quality_values[i];
	}
    }
    (void)fprintf(phd_out, "QUALITY_LEVELS: %d\n", qv_max);
    (void)fprintf(phd_out, "TIME: ");
    current = time(NULL);
    (void)fputs((char*)ctime(&current), phd_out);
    (void)fprintf(phd_out, "TRACE_ARRAY_MIN_INDEX: 0\n");
    (void)fprintf(phd_out, "TRACE_ARRAY_MAX_INDEX: %d\n", num_datapoints - 1);
    (void)fprintf(phd_out, "TRIM: %d %d %f\n", leftTrim, rightTrim,
                  pow(10, -0.1 * trim_threshold));
    if ((chemistry != NULL)
    && ((memcmp(chemistry, "DP", 2) == 0)
	|| (memcmp(chemistry, "DyeP", 4) == 0)))
    {
	(void)fprintf(phd_out, "CHEM: prim\n");
    }
    else if ((chemistry != NULL)
    && ((memcmp(chemistry, "DT", 2) == 0)
	|| (memcmp(chemistry, "DyeT", 4) == 0)))
    {
	(void)fprintf(phd_out, "CHEM: term\n");
    }
    else {
	(void)fprintf(phd_out, "CHEM: unknown\n");
    }
    (void)fprintf(phd_out, "DYE: big\n\n");
    (void)fprintf(phd_out, "END_COMMENT\n\n");
    (void)fprintf(phd_out, "BEGIN_DNA\n");

    /* Output the called bases, their quality values, and locations */
    for (i = 0; i < num_bases; i++) {
        (void)fprintf(phd_out, "%c %d %d\n", tolower((int)called_bases[i]),
			quality_values[i], called_locs[i]);
    }

    /* Create footer */
    (void)fprintf(phd_out, "END_DNA\n\n");
    (void)fprintf(phd_out, "END_SEQUENCE\n");

    (void)fclose(phd_out);
    return SUCCESS;
}


int
output_four_multi_fasta_files(char *multiseqsFileName,
    char *multiqualFileName, char *multilocsFileName, char *multistatFileName,
    int num_bases, char *called_bases, uint8_t *quality_values, int *called_locs, 
    double frac_QV20_with_shoulders, char *status_code, Options options)
{
    int i, j;
    char *seq_name, ttuner_name[BUFLEN];
    FILE *seqs_out = NULL, *qual_out = NULL, 
         *locs_out = NULL, *stat_out = NULL;

    if ((seqs_out = fopen(multiseqsFileName, "a")) == NULL) {
        error(multiseqsFileName, "couldn't open", errno);
        return ERROR;
    }
    
    if ((qual_out = fopen(multiqualFileName, "a")) == NULL) {
        error(multiqualFileName, "couldn't open", errno);
        return ERROR;
    }

    if ((locs_out = fopen(multilocsFileName, "a")) == NULL) {
        error(multilocsFileName, "couldn't open", errno);
        return ERROR;
    }

    if ((stat_out = fopen(multistatFileName, "a")) == NULL) {
        error(multistatFileName, "couldn't open", errno);
        return ERROR;
    }

    /* Use the name of the sample file, sans path, as the sequence name */
#ifdef __WIN32
    if ((seq_name = strrchr(options.file_name, '\\')) != NULL) {
#else
    if ((seq_name = strrchr(options.file_name, '/')) != NULL) {
#endif
        seq_name++;
    }
    else {
        seq_name = options.file_name;
    }

    fprintf(seqs_out, ">%s \n", seq_name);
    if (strcmp(status_code, "TT_SUCCESS") == 0)
    { 
        for (i = 0; i < num_bases; ) 
        {
            for (j = 0; (i < num_bases) && (j < 60); j++, i++) 
            {
                putc(called_bases[i], seqs_out);
            }
            putc('\n', seqs_out);
        }
    }
    else
        fprintf(seqs_out, "%s\n", "CACCA");

    fprintf(qual_out, ">%s \n", seq_name);
    if (strcmp(status_code, "TT_SUCCESS") == 0)
    {
        for (i = 0; i < num_bases; i++) {
            fprintf(qual_out, "%2d ", quality_values[i]);
            if ((i%20 == 19) || (i == num_bases - 1)) {
                fprintf(qual_out, "\n");               
            }
        }
    }
    else
        fprintf(qual_out, "%s\n", "00 00 00 00 00");

    fprintf(locs_out, ">%s \n", seq_name);
    if (strcmp(status_code, "TT_SUCCESS") == 0)
    {
        for (i = 0; i < num_bases; i++) {
            fprintf(locs_out, "%5d ", called_locs[i]);

            if ((i%15 == 14) || (i == num_bases - 1)) {
                fprintf(locs_out, "\n");
            }
        }
    }
    else
        fprintf(locs_out, "%s\n", "00 00 00 00 00");
    strcpy(ttuner_name, TT_VERSION + 3);
    fprintf(stat_out, ">%s ttuner%s %s %3.2f\n", seq_name, ttuner_name, 
        status_code, frac_QV20_with_shoulders);
    
    fclose(seqs_out);
    fclose(qual_out);
    fclose(locs_out);
    fclose(stat_out);
 
    return SUCCESS;
}

/********************************************************************************
 * This function writes out a FASTA file with the (potentially) recalled
 * bases.  Its synopsis is:
 *
 * success = Btk_output_fasta_file(file_name, path, called_bases, num_bases,
 *					verbose)
 *
 * where
 *	file_name	is the name of the sample file
 *	path		is the path name of the directory in which to write
 *			the .seq file, if any (NULL means current dir)
 *	called_bases	is an array of base calls
 *	num_bases	is the number of elements in called_bases
 *	verbose		is whether to write status messages to stderr, and
 *			how verbosely
 *
 *	success		is SUCCESS or ERROR
 ********************************************************************************
 */
int
Btk_output_fasta_file(
    int FastaType,
    char *file_name,
    char *path,
    char *multiseqFileName,
    char *called_bases, 
    int num_bases,
    int left_trim_point,
    int right_trim_point,
    int verbose)
{
    int i, j;
    char *seq_name, fasta_file_name[MAXPATHLEN];
    FILE *fasta_out = NULL, *dir_out = NULL, *multi_out = NULL;

    /* Use the name of the sample file, sans path, as the sequence name. */
#ifdef __WIN32
    if ((seq_name = strrchr(file_name, '\\')) != NULL) {
#else
    if ((seq_name = strrchr(file_name, '/')) != NULL) {
#endif
        seq_name++;
    }
    else {
        seq_name = file_name;
    }

    if (FastaType & NAME_DIR) {
#ifdef __WIN32
        sprintf(fasta_file_name, "%s\\%s.seq", path, seq_name);
#else
        sprintf(fasta_file_name, "%s/%s.seq", path, seq_name);
#endif
        if ((dir_out = fopen(fasta_file_name, "w")) == NULL) {
            error(fasta_file_name, "couldn't open", errno);
            return ERROR;
        }
    }

    if (FastaType & NAME_FILES) {
        sprintf(fasta_file_name, "%s.seq", seq_name);

        if ((fasta_out = fopen(fasta_file_name, "w")) == NULL) {
            error(fasta_file_name, "couldn't open", errno);
            return ERROR;
        }
    }

    if (FastaType & NAME_MULTI) {
        if ((multi_out = fopen(multiseqFileName, "a")) == NULL) {
            error(multiseqFileName, "couldn't open", errno);
            fclose(fasta_out);
            return ERROR;
        }
    }

    if (fasta_out) {
        fprintf(fasta_out, ">%s %d %d %d\n", seq_name, num_bases, left_trim_point,
                right_trim_point - left_trim_point + 1);
    }

    if (dir_out) {
        fprintf(dir_out, ">%s %d %d %d\n", seq_name, num_bases, left_trim_point,
                right_trim_point - left_trim_point + 1);
    }

    if (multi_out) {
        fprintf(multi_out, ">%s %d %d %d\n", seq_name, num_bases, left_trim_point,
                right_trim_point - left_trim_point + 1);
    }

    if (fasta_out) {
        for (i = 0; i < num_bases; ) {
            for (j = 0; (i < num_bases) && (j < BTK_FASTA_WIDTH); j++, i++) {
                putc(called_bases[i], fasta_out);
            }

            putc('\n', fasta_out);
        }
    }

    if (dir_out) {
        for (i = 0; i < num_bases; ) {
            for (j = 0; (i < num_bases) && (j < BTK_FASTA_WIDTH); j++, i++) {
                putc(called_bases[i], dir_out);
            }

            putc('\n', dir_out);
        }
    }

    if (multi_out) {
        for (i = 0; i < num_bases; ) {
            for (j = 0; (i < num_bases) && (j < BTK_FASTA_WIDTH); j++, i++) {
                putc(called_bases[i], multi_out);
            }

            putc('\n', multi_out);
        }
    }

    if (fasta_out) {
        fclose(fasta_out);
    }

    if (dir_out) {
        fclose(dir_out);
    }

    if (multi_out) {
        fclose(multi_out);
    }

    return SUCCESS;
}

/*******************************************************************************
 *
 * This function writes out a TraceTuner's Intrinsic Peaks file,
 *     of the '.tip' variety.
 * Its synopsis is:
 *
 * success = Btk_output_tip_file(data, color2base, options),
 *
 * where
 *      data            is the pointer to the data structure
 *      color2base      is a char array of length 4 which relates
 *                         a color number to base
 *      options         is a structure of type Options
 *      success         is SUCCESS or ERROR
 */

int
Btk_output_tip_file( Data *data, char *color2base, Options options)
{
    char *seq_name, tip_file_name[MAXPATHLEN];
    int i, j, halfwidth, ipos, color, ibeg, iend, *y, num_points;
    FILE *tip_out;

    /* Use the name of the sample file, sans path, as the sequence name */
#ifdef __WIN32
    if ((seq_name = strrchr(options.file_name, '\\')) != NULL) {
#else
    if ((seq_name = strrchr(options.file_name, '/')) != NULL) {
#endif
        seq_name++;
    }
    else {
        seq_name = options.file_name;
    }

    /*
     * Create the name of TIP file by appending suffix ".tip" to the sample
     * file name.
     */
#ifdef __WIN32
    if (options.tip_dir[0] != '.')
        (void)sprintf(tip_file_name, "%s\\%s.tip", options.tip_dir, 
            seq_name);
    else
        (void)sprintf(tip_file_name, "%s.tip", seq_name);           
#else
    if (options.tip_dir[0] != '.')
        (void)sprintf(tip_file_name, "%s/%s.tip", options.tip_dir, 
            seq_name);
    else
        (void)sprintf(tip_file_name, "%s.tip", seq_name);           
#endif

    /* Open the TIP file */
    if ((tip_out = fopen(tip_file_name, "w")) == NULL) {
        error(tip_file_name, "couldn't open", errno);
        return ERROR;
    }
    if (options.Verbose > 2) {
        (void)fprintf(stderr, "Writing TIP output to %s\n", tip_file_name);
    }
    (void)fprintf(tip_out, ">%s\n", seq_name);

    y = CALLOC(int, 100);
    for (i=0; i<data->peak_list_len; i++) {

        /* Determine the beginning and end of the intrinsic peak */
        ipos = data->peak_list[i]->ipos;
        halfwidth = INT_DBL(data->peak_list[i]->orig_width +
                          data->peak_list[i]->beta * 2.);
        ibeg = ipos - halfwidth;
        iend = ipos + halfwidth;

        /* Do not count intrinsic peaks of width <4 */
        if (iend - ibeg + 1 < 4) {
            continue;
        }

        /* Output peak to the TIP file if the peak position is
         * not outside of the data range
         */
        color = data->peak_list[i]->color_index;
        if ((ipos < 0) || (ipos > data->color_data[color].length-1))
            continue;

        num_points = 0;
        for (j=ipos; j<=iend; j++) {
            y[num_points] = INT_DBL(Shape(data->peak_list[i]->C0, 
                                       data->peak_list[i]->orig_width/2.,
                                       data->peak_list[i]->beta, 
                                       j-ipos, &options));
            num_points++;
            if (y[num_points-1] == 0) break;
        }

        /* Output the data for current intrinsic peak */
        (void)fprintf(tip_out, "%c %d %d %d %d %1.3f %d\n",
            color2base[color], ipos, num_points,
            data->peak_list[i]->beg, data->peak_list[i]->end,
            data->peak_list[i]->resolution, data->peak_list[i]->type);
        for (j=0; j<num_points; j++) {
            (void)fprintf(tip_out, "%d ", y[j]);
        }
        (void)fprintf(tip_out, "\n");
    }
    FREE(y);
    return SUCCESS;
}

/*******************************************************************************
 *
 * This function writes out a TraceTuner's alternative base calls file
 *     of the '.tab' variety.
 * Its synopsis is:
 *
 * success = Btk_output_tab_file(num2, altbases, options),
 *
 * where
 *      num2            is the total number of alternative base calls
 *                      (both deleted and substituted)
 *      altbases        array of alternative bases
 *      options         is a structure of type Options
 *      success         is SUCCESS or ERROR
 */

int
Btk_output_tab_file(int num2, AltBase *altbases,
    Options options)
{
    char *seq_name, tab_file_name[MAXPATHLEN];
    int i, num_sub, num_del;
    FILE *tab_out;

    /* Use the name of the sample file, sans path, as the sequence name */
#ifdef __WIN32
    if ((seq_name = strrchr(options.file_name, '\\')) != NULL) {
#else
    if ((seq_name = strrchr(options.file_name, '/')) != NULL) {
#endif
        seq_name++;
    }
    else {
        seq_name = options.file_name;
    }

    /* Removing the suffix ".Z" or ".gz" from the name of compressed sample */
    if (seq_name[strlen(seq_name) - 3] == '.' &&
        seq_name[strlen(seq_name) - 2] == 'g' &&
        seq_name[strlen(seq_name) - 1] == 'z') {
        seq_name[strlen(seq_name) - 3] = '\0';
    }
    if (seq_name[strlen(seq_name) - 2] == '.' &&
        seq_name[strlen(seq_name) - 1] == 'Z') {
        seq_name[strlen(seq_name) - 2] = '\0';
    }

    /* Create the name of ABC file by appending suffix ".tab" to the sample
     * file name.
     */
#ifdef __WIN32
    if (options.tab_dir[0] != '.')
        (void)sprintf(tab_file_name, "%s\\%s.tab", options.tab_dir,
            seq_name);
    else
        (void)sprintf(tab_file_name, "%s.tab", seq_name);
#else
    if (options.tab_dir[0] != '.')
        (void)sprintf(tab_file_name, "%s/%s.tab", options.tab_dir,
            seq_name);
    else
        (void)sprintf(tab_file_name, "%s.tab", seq_name);
#endif

    /* Open the ABC file */
    if ((tab_out = fopen(tab_file_name, "w")) == NULL) {
        error(tab_file_name, "couldn't open", errno);
        return ERROR;
    }
    if (options.Verbose > 2) {
        (void)fprintf(stderr, "Writing ABC output to %s\n", tab_file_name);
    }

    (void)fprintf(tab_out, "# CHROMAT_FILE: %s\n", seq_name);
    (void)fprintf(tab_out, "# SOFTWARE_VERSION: %s\n",TT_VERSION);

    if (num2 == 0) {
        (void)fprintf(tab_out, "#\n");
        fprintf(tab_out, "# ABC_STATUS: No alternative base calls\n");
    } else {
        /* Count substituted and deleted base calls */
        num_sub = 0;
        for (i=0; i<num2; i++) {
            if (altbases[i].base_index >= 0)
                num_sub++;
        }
        num_del = num2 - num_sub;

        (void)fprintf(tab_out, "# NUM_ABC: %d\n", num2);
        (void)fprintf(tab_out, "# NUM_SUBSTITUTIONS: %d\n", num_sub);
        (void)fprintf(tab_out, "# NUM_DELETIONS: %d\n", num_del);
        (void)fprintf(tab_out, "# base2   qv2       pos2       ind2\n");
        for (i=0; i<num2; i++) {
            (void)fprintf(tab_out, "   %c      %2d       %5d      %4d\n",
                altbases[i].base, altbases[i].qv,
                altbases[i].pos,  altbases[i].base_index);
        }
    }

    (void)fclose(tab_out);
    return SUCCESS;

}


/********************************************************************************
 * Function: delta_out_samples1
 ********************************************************************************
 */
void
delta_out_samples1(unsigned char samples[], int num_samples)
{
    int i;
    unsigned char p_delta = 0, p_sample;

    for (i = 0; i < num_samples; i++)
    {
        p_sample = samples[i];
        samples[i] = samples[i] - p_delta;
        p_delta = p_sample;
    }

    for (p_delta = 0, i = 0; i < num_samples; i++)
    {
        p_sample = samples[i];
        samples[i] = samples[i] - p_delta;
        p_delta = p_sample;
    }
}

/********************************************************************************
 * Function: delta_out_samples2
 ********************************************************************************
 */
void
delta_out_samples2(unsigned short samples[], int num_samples)
{
    int i;
    unsigned short p_delta = 0, p_sample;

    for (i = 0; i < num_samples; i++)
    {
        p_sample = samples[i];
        samples[i] = samples[i] - p_delta;
        p_delta = p_sample;
    }

    for (p_delta = 0, i = 0; i < num_samples; i++)
    {
        p_sample = samples[i];
        samples[i] = samples[i] - p_delta;
        p_delta = p_sample;
    }
}

/********************************************************************************
 * Function: write_sample1
 ********************************************************************************
 */
static int
write_sample1(FILE *fp, TT_Samples1 *sample)
{
    if (writeuchar(fp, sample->sample_A) != SUCCESS) return ERROR;
    if (writeuchar(fp, sample->sample_C) != SUCCESS) return ERROR;
    if (writeuchar(fp, sample->sample_G) != SUCCESS) return ERROR;
    if (writeuchar(fp, sample->sample_T) != SUCCESS) return ERROR;

    return SUCCESS;
}

/********************************************************************************
 * Function: write_sample2
 ********************************************************************************
 */
static int
write_sample2(FILE *fp, TT_Samples2 *sample)
{
    if (writeushort(fp, sample->sample_A) != SUCCESS) return ERROR;
    if (writeushort(fp, sample->sample_C) != SUCCESS) return ERROR;
    if (writeushort(fp, sample->sample_G) != SUCCESS) return ERROR;
    if (writeushort(fp, sample->sample_T) != SUCCESS) return ERROR;

    return SUCCESS;
}

/********************************************************************************
 * Function: write_samples31
 ********************************************************************************
 */
static int
write_samples31(FILE *fp, TT_Samples1 *samples, int num_samples)
{
    int i, base;
    unsigned char *samples_out;

    samples_out = (unsigned char *) malloc(num_samples * sizeof(unsigned char));
    if (samples_out == NULL) return ERROR;

    for (base = 0; base < 4; base++)
    {
        for (i = 0; i < num_samples; i++)
            switch (base)
            {
                case 0:
                    samples_out[i] = samples[i].sample_A;
                    break;
                case 1:
                    samples_out[i] = samples[i].sample_C;
                    break;
                case 2:
                    samples_out[i] = samples[i].sample_G;
                    break;
                case 3:
                    samples_out[i] = samples[i].sample_T;
                    break;
            }

        delta_out_samples1(samples_out, num_samples);
        if (num_samples != (int)fwrite(samples_out, 1, num_samples, fp))
        {
            FREE(samples_out);

            return ERROR;
        }
    }

    FREE(samples_out);

    return SUCCESS;
}

/********************************************************************************
 * Function: write_samples32
 ********************************************************************************
 */
static int
write_samples32(FILE *fp, TT_Samples2 *samples, int num_samples)
{
    int i, base;
    unsigned short *samples_out;

    samples_out = (unsigned short *) malloc(num_samples * sizeof(unsigned short));
    if (samples_out == NULL) return ERROR;

    for (base = 0; base < 4; base++)
    {
        for (i = 0; i < num_samples; i++)
            switch (base)
            {
                case 0:
                    samples_out[i] = (&samples[i])->sample_A;
                    break;
                case 1:
                    samples_out[i] = (&samples[i])->sample_C;
                    break;
                case 2:
                    samples_out[i] = (&samples[i])->sample_G;
                    break;
                case 3:
                    samples_out[i] = (&samples[i])->sample_T;
                    break;
            }

        delta_out_samples2(samples_out, num_samples);

        for (i = 0; i < num_samples; i++)
            if (writeushort(fp, samples_out[i]) != SUCCESS)
            {
                FREE(samples_out);

                return ERROR;
            }
    }

    FREE(samples_out);

    return SUCCESS;
}

/********************************************************************************
 * Function: write_base
 ********************************************************************************
 */
static int
write_base(FILE *fp, SCF_Bases_Rec *base)
{
    if (writeuint(fp, base->peak_index) != SUCCESS) return ERROR;
    if (writeuchar(fp, base->prob_A)    != SUCCESS) return ERROR;
    if (writeuchar(fp, base->prob_C)    != SUCCESS) return ERROR;
    if (writeuchar(fp, base->prob_G)    != SUCCESS) return ERROR;
    if (writeuchar(fp, base->prob_T)    != SUCCESS) return ERROR;
    if (writeuchar(fp, (unsigned char) base->base) != SUCCESS) return ERROR;
    if (writeuchar(fp, base->spare[0])  != SUCCESS) return ERROR;
    if (writeuchar(fp, base->spare[1])  != SUCCESS) return ERROR;
    if (writeuchar(fp, base->spare[2])  != SUCCESS) return ERROR;

    return SUCCESS;
}

/********************************************************************************
 * Function: write_bases3
 ********************************************************************************
 */
static int
write_bases3(FILE *fp, SCF_Bases_Rec *bases, int num_bases)
{
    int i;
    unsigned char *buf1;

    buf1 = (unsigned char *) malloc(8 * num_bases);
    if (buf1 == NULL) return ERROR;

    for (i = 0; i < num_bases; i++)
        if (writeuint(fp, (&bases[i])->peak_index) != SUCCESS) return ERROR;

    for (i = 0; i < num_bases; i++)
    {
        buf1[i] = (&bases[i])->prob_A;
        buf1[i + num_bases] = (&bases[i])->prob_C;
        buf1[i + 2 * num_bases] = (&bases[i])->prob_G;
        buf1[i + 3 * num_bases] = (&bases[i])->prob_T;
        buf1[i + 4 * num_bases] = (&bases[i])->base;
        buf1[i + 5 * num_bases] = (&bases[i])->spare[0];
        buf1[i + 6 * num_bases] = (&bases[i])->spare[1];
        buf1[i + 7 * num_bases] = (&bases[i])->spare[2];
    }

    if (8 * num_bases != (int)(fwrite(buf1, 1, 8 * num_bases, fp)))
    {
        FREE(buf1);

        return ERROR;
    }

    FREE(buf1);

    return SUCCESS;
}

/********************************************************************************
 * Function: write_scf_comments
 ********************************************************************************
 */
static int
write_scf_comments(FILE *fp, char *comments, size_t comments_size)
{
    if (fwrite(comments, comments_size, 1, fp) != 1) return ERROR;

    return SUCCESS;
}

/********************************************************************************
 * Function: write_scf_header
 ********************************************************************************
 */
static int
write_scf_header(FILE *fp, SCF_Header *header)
{
    int i;

    if (writeuint(fp, header->magic_number)     != SUCCESS) return ERROR;
    if (writeuint(fp, header->samples)          != SUCCESS) return ERROR;
    if (writeuint(fp, header->samples_offset)   != SUCCESS) return ERROR;
    if (writeuint(fp, header->bases)            != SUCCESS) return ERROR;
    if (writeuint(fp, header->bases_left_clip)  != SUCCESS) return ERROR;
    if (writeuint(fp, header->bases_right_clip) != SUCCESS) return ERROR;
    if (writeuint(fp, header->bases_offset)     != SUCCESS) return ERROR;
    if (writeuint(fp, header->comments_size)    != SUCCESS) return ERROR;
    if (writeuint(fp, header->comments_offset)  != SUCCESS) return ERROR;

    if (fwrite(header->version, sizeof(header->version), 1, fp) != 1) return ERROR;

    if (writeuint(fp, header->sample_size)      != SUCCESS) return ERROR;
    if (writeuint(fp, header->code_set)         != SUCCESS) return ERROR;
    if (writeuint(fp, header->private_size)     != SUCCESS) return ERROR;
    if (writeuint(fp, header->private_offset)   != SUCCESS) return ERROR;

    for (i = 0; i < 18; i++)
        if (writeuint(fp, header->spare[i]) != SUCCESS)
            return ERROR;

    return SUCCESS;
}

/********************************************************************************
 * This function writes a SCFfile
 ********************************************************************************
 */
int
output_scf_file(
    char *path,
    char *scf_dir,
    char *called_bases,
    int  *called_peak_locs,
    uint8_t  *quality_values,
    int   num_called_bases,
    int   num_datapoints,
    int  *chromatogram0,
    int  *chromatogram1,
    int  *chromatogram2,
    int  *chromatogram3,
    char *color2base,
    char *chemistry)
{
    char *seq_name, scf_file_name[MAXPATHLEN], *suffix = NULL;
    FILE *scf_out;
    char comments[2048];
    SCF_Header header;
    SCF_Bases_Rec base;
    int i, scf_version = 2;

#ifdef __WIN32
    if ((seq_name = strrchr(path, '\\')) != NULL)
#else
    if ((seq_name = strrchr(path, '/')) != NULL)
#endif
        seq_name++;
    else
        seq_name = path;

    /* Remove suffix ab1 or abi from the seq_name */
    suffix = strrchr(seq_name, '.');
    if (suffix && (!strcmp(suffix+1, "ab1") || !strcmp(suffix+1, "abi"))) {
        suffix[0] = '\0';
    }

    /* If output directory is not current, build the full path */
    if (scf_dir[0] != '\0')
    {
#ifdef __WIN32
        sprintf(scf_file_name, "%s\\%s.scf", scf_dir, seq_name);
#else
        sprintf(scf_file_name, "%s/%s.scf", scf_dir, seq_name);
#endif
    }
    else
        sprintf(scf_file_name, "%s.scf", seq_name);

    if ((scf_out = fopen(scf_file_name, "wb")) == NULL) {
        error(scf_file_name, "couldn't open", errno);
        return ERROR;
    }

    sprintf(comments, "DYEP=%s\nCONV=%s", chemistry, TT_VERSION);

    header.magic_number      = TT_SCF_MAGIC;
    header.samples           = num_datapoints;
    header.samples_offset    = (unsigned int) sizeof(SCF_Header);
    header.bases             = num_called_bases;
    header.bases_left_clip   = 0;
    header.bases_right_clip  = 0;
    header.sample_size       = (max_colordata_value < 256) ? 1 : 2;
    header.bases_offset      = (unsigned int) (header.samples_offset + header.samples
                               * ((header.sample_size == 2) ? 8 : 4));
    header.comments_size     = (unsigned int) strlen(comments) + 1;
    header.comments_offset   = (unsigned int) (header.bases_offset + header.bases
                               * sizeof(SCF_Bases_Rec));
    strncpy(header.version, (scf_version == 2) ? "2.00" : "3.00", 4);
    header.code_set          = 0;
    header.private_size      = 0;
    header.private_offset    = 0;

    for (i = 0; i < 18; i++)
        header.spare[i] = 0;

    if (write_scf_header(scf_out, &header) != SUCCESS) return ERROR;

    if (max_colordata_value < 256)
    {
        TT_Samples1 sample, *samples;

        if (scf_version == 2)
        {
            for (i = 0; i < (int)header.samples; i++)
            {
                sample.sample_A = (unsigned char) chromatogram0[i];
                sample.sample_C = (unsigned char) chromatogram1[i];
                sample.sample_G = (unsigned char) chromatogram2[i];
                sample.sample_T = (unsigned char) chromatogram3[i];

                if (write_sample1(scf_out, &sample) != SUCCESS) return ERROR;
            }
        }
        else
        {
            samples = (TT_Samples1 *) malloc(header.samples * sizeof(TT_Samples1));
            if (samples == NULL) return ERROR;

            for (i = 0; i < (int)header.samples; i++)
            {
                samples[i].sample_A = (unsigned char) chromatogram0[i];
                samples[i].sample_C = (unsigned char) chromatogram1[i];
                samples[i].sample_G = (unsigned char) chromatogram2[i];
                samples[i].sample_T = (unsigned char) chromatogram3[i];
            }

            if (write_samples31(scf_out, samples, header.samples) != SUCCESS)
                return ERROR;

            FREE(samples);
        }
    }
    else
    {
        TT_Samples2 sample, *samples;

        if (scf_version == 2)
        {
            for (i = 0; i < (int)header.samples; i++)
            {
                sample.sample_A = (unsigned short) chromatogram0[i];
                sample.sample_C = (unsigned short) chromatogram1[i];
                sample.sample_G = (unsigned short) chromatogram2[i];
                sample.sample_T = (unsigned short) chromatogram3[i];

                if (write_sample2(scf_out, &sample) != SUCCESS) return ERROR;
            }
        }
        else
        {
            samples = (TT_Samples2 *) malloc(header.samples * sizeof(TT_Samples2));
            if (samples == NULL) return ERROR;

            for (i = 0; i < (int)header.samples; i++)
            {
                samples[i].sample_A = (unsigned short) chromatogram0[i];
                samples[i].sample_C = (unsigned short) chromatogram1[i];
                samples[i].sample_G = (unsigned short) chromatogram2[i];
                samples[i].sample_T = (unsigned short) chromatogram3[i];
            }

            if (write_samples32(scf_out, samples, header.samples) != SUCCESS) return ERROR;

            FREE(samples);
        }
    }

    if (scf_version == 2)
    {
        for (i = 0; i < (int)header.bases; i++)
        {
            base.peak_index = called_peak_locs[i];
            base.base = called_bases[i];
            base.spare[0] = 0;
            base.spare[1] = 0;
            base.spare[2] = 0;
            base.prob_A = 0;
            base.prob_C = 0;
            base.prob_G = 0;
            base.prob_T = 0;

            switch(base.base)
            {
                case 'A':
                case 'a':
                    base.prob_A = quality_values[i];
                    break;
                case 'C':
                case 'c':
                    base.prob_C = quality_values[i];
                    break;
                case 'G':
                case 'g':
                    base.prob_G = quality_values[i];
                    break;
                case 'T':
                case 't':
                    base.prob_T = quality_values[i];
                    break;
            }

            if (write_base(scf_out, &base) != SUCCESS) return ERROR;
        }
    }
    else
    {
        SCF_Bases_Rec *bases;

        bases = (SCF_Bases_Rec *) malloc(header.bases * sizeof(SCF_Bases_Rec));
        if (bases == NULL) return ERROR;

        for (i = 0; i < (int)header.bases; i++)
        {
            bases[i].peak_index = called_peak_locs[i];
            bases[i].base = called_bases[i];
            bases[i].spare[0] = 0;
            bases[i].spare[1] = 0;
            bases[i].spare[2] = 0;
            bases[i].prob_A = 0;
            bases[i].prob_C = 0;
            bases[i].prob_G = 0;
            bases[i].prob_T = 0;

            switch(bases[i].base)
            {
                case 'A':
                case 'a':
                    bases[i].prob_A = quality_values[i];
                    break;
                case 'C':
                case 'c':
                    bases[i].prob_C = quality_values[i];
                    break;
                case 'G':
                case 'g':
                    bases[i].prob_G = quality_values[i];
                    break;
                case 'T':
                case 't':
                    bases[i].prob_T = quality_values[i];
                    break;
            }
        }

        if (write_bases3(scf_out, bases, header.bases) != SUCCESS) return ERROR;
    }

    if (write_scf_comments(scf_out, comments, (size_t) header.comments_size) != SUCCESS)
        return ERROR;

    fclose(scf_out);

    return SUCCESS;
}




/* *****************************************************************************
 * File: get_phd_num_bases
 * Purpose: read called bases, their locations and quality values from .phd.1
 *          file
 *******************************************************************************/

int
get_phd_num_bases(char *phd_file_name, BtkMessage *message)
{
    FILE *phd_inp;
    int   num_bases, qv, pos;
    char  c, buffer[BUFLEN];

    if ((phd_inp = fopen(phd_file_name, "r")) == NULL) {
        return ERROR;
    }

    while ( ( fgets(buffer, BUFLEN, phd_inp) != NULL)
            && (strcmp(buffer, "BEGIN_DNA\n")!= 0 ) ) {
    }

    num_bases=0;
    while ((fgets(buffer, BUFLEN, phd_inp) != NULL) && 
           (strcmp(buffer, "END_DNA\n")!= 0 ) && 
           (num_bases < MAX_BASES_LEN) ) {
        sscanf( buffer, "%c %d %d\n", &c, &qv, &pos );
        if(islower((int)c)) c=toupper((int)c); /* always use upper */
        num_bases++;
    }
    fclose(phd_inp);
    return num_bases;
}

/* *****************************************************************************
 * File: Btk_read_phd_file
 * Purpose: read called bases, their locations and quality values from .phd.1
 *          file
 *******************************************************************************/

int
Btk_read_phd_file(char *phd_file_name, char *called_bases, uint8_t *qualities,
    int *called_locs, int *num_bases, BtkMessage *message)
{
    FILE *phd_inp;
    int   i;
    char  buffer[BUFLEN];
    char  c;
    int   qv, pos;

    if ((phd_inp = fopen(phd_file_name, "r")) == NULL) {
        fprintf(stderr, "Can not open input phd file %s\n",
        phd_file_name);
        return ERROR;
    }

    while ( ( fgets(buffer, BUFLEN, phd_inp) != NULL)
            && (strcmp(buffer, "BEGIN_DNA\n")!= 0 ) ) {
    }

    i=0;
    while ( ( fgets(buffer, BUFLEN, phd_inp) != NULL)
            && (strcmp(buffer, "END_DNA\n")!= 0 ) && (i < *num_bases) ) {
        sscanf( buffer, "%c %d %d\n", &c, &qv, &pos );
        if(islower((int)c)) c=toupper((int)c); /* always use upper */
        called_bases[i]=c;
        qualities[i]=qv;
        called_locs[i]=pos;
        i++;
    }
    *num_bases=i;

    (void)fclose(phd_inp);
    return SUCCESS;
}

/* *****************************************************************************
 * File: Btk_read_tab_file
 * Purpose: read alternative bases, their quality values, locations and base
 *          indexes from TAB file
 *******************************************************************************/

int
Btk_read_tab_file(char *tab_file_name, char *bases2, uint8_t *qvs2,
    int *called_locs2, int *base_ind2, int *num_bases, BtkMessage *message)
{
    FILE *tab_inp;
    int   i;
    char  buffer[BUFLEN];
    char  c;
    int   pos, bind;
    int   qv;

    if ((tab_inp = fopen(tab_file_name, "r")) == NULL) {
        fprintf(stderr, "Can not open input tab file %s\n", tab_file_name);
        return ERROR;
    }

    i=0;
    while ( ( fgets(buffer, BUFLEN, tab_inp) != NULL) && (i < *num_bases) ) {
        if (buffer[0] == '#')
            continue;
        sscanf( buffer, "   %c      %2d       %5d      %4d\n", &c, &qv, &pos, &bind);
        if(islower((int)c)) c=toupper((int)c); /* always use upper */
        bases2[i]=c;
        qvs2[i]=(uint8_t)qv;
        called_locs2[i]=pos;
        base_ind2[i]=bind; 
        i++;
    }
    *num_bases=i;

    (void)fclose(tab_inp);
    return SUCCESS;
}

/* *****************************************************************************
 * File:    Btk_output_poly_file
 *******************************************************************************/

int
Btk_output_poly_file(Data *data, Options *options, BtkMessage *message)
{
    int    i, j, jc, peak_index, num_peaks; 
    double max_height, ave_area = 0.;      
    double ampA, ampC, ampG, ampT;
    char  *seq_name, poly_name[MAX_FILE_NAME_LENGTH];
    Peak   peak;
    FILE  *poly_out;

    /* Get sequence name by ignoring the rest of the path name */
#ifdef __WIN32
    if ((seq_name = strrchr(options->file_name, '\\')) != NULL)
        seq_name++;
#else
    if ((seq_name = strrchr(options->file_name, '/')) != NULL)
        seq_name++;
#endif
    else
        seq_name = options->file_name;

    /* Form .poly file's full name by prefixing the seq_name by
     * poly directory name, if necessary
     */
    if (options->poly_dir[0] == '\0') {
        sprintf(poly_name, "%s.poly", seq_name);
    }
    else
#ifdef __WIN32
        sprintf(poly_name, "%s\\%s.poly",
            options->poly_dir, seq_name);
#else
        sprintf(poly_name, "%s/%s.poly",
            options->poly_dir, seq_name);
#endif

    if ((poly_out=fopen(poly_name,"w"))== NULL) {
        sprintf(message->text, "Unable to open POLY file '%s'\n",
        poly_name);
        return ERROR;
    }
    fprintf(poly_out, "%s 1.0 1.0 1.0 1.0 1.0\n", seq_name);

    /* Precompute average area of 10 firts called peaks */
    num_peaks = QVMIN(data->bases.length, 10);
    if (num_peaks <= 0)
       return SUCCESS;
 
    for (i=0; i< num_peaks; i++)
    {
        ave_area += data->bases.called_peak_list[i]->area;
    }
    ave_area /= (double)num_peaks;

    /* Process all bases in the order of their positions */
    for (i=0; i<data->bases.length; i++) {

        double area, area2=0., rel_area, rel_area2=0.;
        int    position = data->bases.called_peak_list[i]->ipos, position2 = -1;
        int    height   = data->bases.called_peak_list[i]->height;
        char   base = data->bases.called_peak_list[i]->base, base2='\0';

        if (i > 10) {
            ave_area -= data->bases.called_peak_list[i-11]->area/10.;
            ave_area += data->bases.called_peak_list[i-1 ]->area/10.;
        }
        area      = data->bases.called_peak_list[i]->area;
        rel_area  = area/ave_area;
        ampA = (double)data->color_data[0].data[position];
        ampC = (double)data->color_data[1].data[position];
        ampG = (double)data->color_data[2].data[position];
        ampT = (double)data->color_data[3].data[position]; 
        
        /* Determine uncalled base and uncalled peak */
        jc = data->bases.called_peak_list[i]->color_index;
        max_height = -1.; 
        for (j=0; j<NUM_COLORS; j++) {
            if (j == jc)
                continue;
            
            if (colordata_find_peak_index_by_location(&(data->color_data[j]),
                position, &peak, &peak_index, message) != SUCCESS) {
                return ERROR;
            }

            if (peak_index < 0)
                continue;
   
            if (is_dye_blob(peak.ipos, &peak, data, 0))
                continue;
 
            if ((i>0) && 
                (peak.pos < (position + data->bases.called_peak_list[i-1]->pos)/2))
                continue;

            if ((i<data->bases.length-1) &&
                (peak.pos >= (position +
                    data->bases.called_peak_list[i+1]->pos)/2))
                continue;
#if 0
            fprintf(stderr, "min_ratio=%f\n", options->min_ratio);
#endif
            if (peak.height < (options->min_ratio*(double)height))
                continue;

            if (peak.height < max_height)
                continue;

            max_height = peak.height;

            area2     = peak.area;
            rel_area2 = area2/ave_area;
            base2     = data->color_data[j].peak_list[peak_index].base;
            position2 = data->color_data[j].peak_list[peak_index].pos;
        }
 
        /* Output */
        if (base2 != '\0') {
            fprintf(poly_out,
                "%c  %d  %.6f  %.6f  %c  %d  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f\n",
                base,  position,  area,  rel_area,
                base2, position2, area2, rel_area2,
                ampA, ampC, ampG, ampT);
        }
        else {
            fprintf(poly_out,
                "%c  %d  %.6f  %.6f  %c  %d  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f\n",
                base, position, area, rel_area, 'N', -1, -1., -1.,
                ampA, ampC, ampG, ampT);
        }
        position2 = -1;
        base2='\0';

    }

    fclose(poly_out);
    return SUCCESS;
}

/* *****************************************************************************
 * File:    Btk_output_poly_file
 *******************************************************************************/

int
Btk_output_poly_file_old(Data *data, int *data_peak_ind1, int *data_peak_ind2,
    int *qv, Options *options, BtkMessage *message)
{
    int    i, ind1, ind2, ipos;
    char  *seq_name, poly_name[MAX_FILE_NAME_LENGTH];
    FILE  *poly_out; 
    Peak **dpl = data->peak_list; 

    /* Get sequence name by ignoring the rest of the path name */
#ifdef __WIN32
    if ((seq_name = strrchr(options->file_name, '\\')) != NULL) 
        seq_name++;
#else
    if ((seq_name = strrchr(options->file_name, '/')) != NULL) 
        seq_name++;
#endif
    else
        seq_name = options->file_name;

    /* Form .poly file's full name by prefixing the seq_name by
     * poly directory name, if necessary 
     */
    if (options->poly_dir[0] == '\0') {
        sprintf(poly_name, "%s.poly", seq_name);
    }
    else
#ifdef __WIN32
        sprintf(poly_name, "%s\\%s.poly", 
            options->poly_dir, seq_name);
#else
        sprintf(poly_name, "%s/%s.poly", 
            options->poly_dir, seq_name);
#endif
 
    if ((poly_out=fopen(poly_name,"w"))== NULL) {
        sprintf(message->text, "Unable to open POLY file '%s'\n", 
        poly_name);
        return ERROR;
    }
    fprintf(poly_out, "%s 1.0 1.0 1.0 1.0 1.0\n", seq_name);

    for (i=0; i<data->bases.length; i++) {

        double area1, area2=0., rel_area1, rel_area2=0.;

        ind1 = data_peak_ind1[i];
        ind2 = data_peak_ind2[i]; 
        ipos = dpl[ind1]->ipos;
        area1     = dpl[ind1]->area;
        rel_area1 = dpl[ind1]->relative_area;
  
        /* Increase, if needed, the area of the 2nd peak 
         * if quality of mixed base is high enough 
         */
        if (ind2 >= 0) {
             area2     = dpl[ind2]->area;
             rel_area2 = dpl[ind2]->relative_area; 
             if (ADJUST_SECONDARY_PEAK && 
                 (qv[i] >= 25) && (area2 < area1 * 0.5)) {
                 area2     = area1     * 0.5;
                 rel_area2 = rel_area1 * 0.5;
             }
             if (ADJUST_SECONDARY_PEAK &&
                 (qv[i] >= 30) && (area2 < area1 * 0.75)) {
                 area2     = area1     * 0.75;
                 rel_area2 = rel_area1 * 0.75;
             }
        }

        /* Output */
        if (ind2 >= 0) { 
            fprintf(poly_out, 
                "%c %d %.6f %.6f %c %d %.6f %.6f %.6f %.6f %.6f %.6f\n",
                dpl[ind1]->base, (int)dpl[ind1]->ipos, area1, rel_area1,
                dpl[ind2]->base, (int)dpl[ind2]->ipos, area2, rel_area2,
                (double)data->color_data[0].data[ipos],
                (double)data->color_data[1].data[ipos], 
                (double)data->color_data[2].data[ipos], 
                (double)data->color_data[3].data[ipos]); 
        }
        else {
            fprintf(poly_out,
                "%c %d %.6f %.6f %c %d %.6f %.6f %.6f %.6f %.6f %.6f\n",
                dpl[ind1]->base, (int)dpl[ind1]->ipos, dpl[ind1]->area,
                dpl[ind1]->relative_area,
                'N', -1, -1., -1.,
                (double)data->color_data[0].data[ipos],
                (double)data->color_data[1].data[ipos],
                (double)data->color_data[2].data[ipos],
                (double)data->color_data[3].data[ipos]);
        }
    }
    
    fclose(poly_out);
    return SUCCESS;
}
