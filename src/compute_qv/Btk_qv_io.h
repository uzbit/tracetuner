
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
 * Copyright (c) 1999-2003 Paracel Inc.  All rights reserved.
 */

/*
 *   Btk_qv_io.h    $Id: Btk_qv_io.h,v 1.10 2009/01/01 14:15:33 gdenisov Exp $    
 */

#define NAME_NONE 0
#define NAME_FILES 1         /* input will come from files named on the command line */
#define NAME_DIR 2           /* input will come from all files in a directory */
#define NAME_FILEOFFILES 4   /* input will come from a file with one filename per line */
#define NAME_MULTI 8

extern int 
read_consensus_from_sample_file(char **, int);

extern int
read_sequence_from_fasta(char *, char **, BtkMessage *);

extern int
read_fasta(char *, Contig *, BtkMessage *);

extern int
get_file_type(char *file_name, BtkMessage *message);

extern int
Btk_read_sample_file(
    char   *file_name,
    int    *num_bases,
    char  **called_bases,
    int     use_edited,
    int   **called_locs,
    uint8_t   **quality_values, 
    int    *num_values,
    int   **avals,
    int   **cvals,
    int   **gvals,
    int   **tvals,
    char  **call_method,
    char  **chemistry,
    char   *status_code,
    int    *file_type,
    Options options,
    BtkMessage *message);


extern void
Btk_release_file_data(
    char  *called_bases,
    int   *called_locs,
    uint8_t   *quality_values, 
    int  **chromatogram,
    char **call_method,
    char **chemistry);

extern int
find_trim_points(int , uint8_t *, int win, float thr, int *left, int *right);

extern int
Btk_output_quality_values(
    int QualType,
    char *file_name,
    char *path,
    char *multiqualFileName,
    uint8_t *quality_values,
    int num_values,
    int left_trim_point,
    int right_trim_point,
    int verbose);

extern int
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
    int  verbose);

extern int
Btk_output_hpr_file(
    char *file_name,
    char *path,
    char *called_bases,
    int *called_locs,
    uint8_t *quality_values,
    int num_bases,
    int num_datapoints,
    int verbose);

extern int
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
    int verbose);

extern int
Btk_output_fastq_file(
	int FastqType,
	char *file_name,
	char *path,
	char *multiseqFileName,
	char *called_bases,
	uint8_t *quality_values,
	int num_bases,
	int left_trim_point,
	int right_trim_point,
    int verbose);

extern int
Btk_output_fasta_file(
    int FastaType,
    char *file_name,
    char *path,
    char *multiseqFileName,
    char *called_bases, 
    int num_bases,
    int left_trim_point,
    int right_trim_point,
    int verbose);

extern int
Btk_output_tip_file(
    Data *data,
    char *color2base,
    Options options);

extern int
Btk_output_tab_file(
    int num2,
    AltBase *altbases,
    Options options);


extern int
output_scf_file(
    char *path,
    char *scf_dir,
    char *called_bases,
    int *called_locs,
    uint8_t *quality_values,
    int num_bases,
    int num_datapoints,
    int *chromatogram0,
    int *chromatogram1,
    int *chromatogram2,
    int *chromatogram3,
    char *color2base,
    char *chemistry);

extern int
get_phd_num_bases(
    char *file_name, 
    BtkMessage *message);

extern int
Btk_read_phd_file(char *full_name, 
    char *called_bases, 
    uint8_t  *quality_values,
    int  *called_locs, 
    int  *num_bases, 
    BtkMessage *message);

extern int
Btk_output_poly_file(
    Data *, 
    Options *, 
    BtkMessage *);

extern int
output_four_multi_fasta_files(char *, char *, char *, char *,
    int , char *, uint8_t *, int *, double, char *, Options );

extern int
Btk_read_tab_file(char *, char *, uint8_t *, int *, int *, int *, 
    BtkMessage *);
