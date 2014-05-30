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
 * $Id: train_data.c,v 1.15 2009/01/12 22:20:15 gdenisov Exp $                  
 *
 * This file contains functions that are shared by train, trainphred 
 * and trainttphred applications
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
#include <stdint.h>

#include "ABI_Toolkit.h"
#include "FileHandler.h"
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

/*********************************************************************************
 * Function: release1                 
 *********************************************************************************
 */
void
release1_sanger(char *bases, int *peak_locs, int **chromatogram,
        uint8_t *quality_values, int NUM_PARAMS, double **params, double *iheight)
{
    int i;

    FREE(bases);
    FREE(peak_locs);

    for (i=0; i<NUM_COLORS; i++) {
        FREE(chromatogram[i]);
    }
    if (params != NULL) {
        for (i=0; i<NUM_PARAMS; i++) {
                FREE(params[i]);
        }
    }
    FREE(quality_values);
    FREE(iheight);
}

void
release1_454(sffRead *read, int NUM_PARAMS, double **params)
{
    int i;

    FREE(read->bases);
    FREE(read->flowgram_values);
    FREE(read->quality_values);
    if (params != NULL) {
        for (i=0; i<NUM_PARAMS; i++) {
                FREE(params[i]);
        }
    }
}

/*********************************************************************************
 * Function: release2
 * Purpose:  release the resources held by the specified objects
 *********************************************************************************
 */
int   
release2(Contig *fragment, Align *best_alignment, Align *vector_start_align,
    Align *vector_end_align, BtkMessage *message) 
{

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
  
/*********************************************************************************
 * Function: is_base                  
 *********************************************************************************
 */
static int
is_base(char b)
{
    if ((b=='A') || (b=='C') || (b=='G') || (b=='T') || (b=='N'))
        return 1;
    return 0;
}

static void
get_num_errors(Align* align, int *ins, int *del, int *sub, int *tot)
{
    int j;

   *ins=0, *del=0, *sub=0, *tot=0;

    for (j=0; j< align->trace_len; j++) {
        if ( is_base(align->trace_qchar[j]) &&
             is_base(align->trace_dchar[j]) &&
            (align->trace_qchar[j] != align->trace_dchar[j]))
           (*sub)++;

        if (!is_base(align->trace_qchar[j]) &&
             is_base(align->trace_dchar[j])) 
           (*del)++;

        if ( is_base(align->trace_qchar[j]) &&
            !is_base(align->trace_dchar[j])) 
           (*ins)++;
    }
   (*tot) = (*sub) + (*ins) + (*del);
   
    return;
}

/* This function appends one Align data structure to the training output
 * file.  Its synopsis is:
 *
 * result = append_align_to_train_data(out, params, align, message)
 *
 * where
 *      out             is the (already open) output stream
 *      params          is the matrix of trace parameters
 *      align           is the address of the Align
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 */
int
append_align_to_train_data(FILE* fout, int NUM_PARAMS, double** params, 
    uint8_t *orig_qvs, double *iheight, double *iheight2, double *ave_iheight, 
    Align* align, int min_dpos, int max_dpos, Options *op, BtkMessage* message)
{
    int     i, j, ismatch;

    for (j=0; j< align->trace_len; j++) {
        if ((align->trace_dpos[j] < min_dpos) ||
            (align->trace_dpos[j] > max_dpos))
            continue;

        ismatch = 0;
        if (align->trace_dir[j] != '2') {
            ismatch = (align->trace_mchar[j]=='|')?1:0;
        }

        fprintf(fout, "%d\t%c\t%d\t%d\t%c",
                align->trace_dpos[j]+1,
                align->trace_dchar[j],
                ismatch,
                align->trace_qpos[j]+1,
                align->trace_qchar[j]);

        if (align->trace_dir[j] == '2') {
            fprintf(fout, "\n");
            continue;
        }

        i = align->trace_qpos[j];

        if (NUM_PARAMS == 4)
        {
            fprintf(fout, "\t%.6f %.6f %.6f %.6f ",
                    params[0][i], params[1][i],
                    params[2][i], params[3][i]);
        }
        else /* NUM_PARAMS == 6 */
        {
            fprintf(fout, "\t%.6f %.6f %.6f %.6f %.6f %.6f ",
                    params[0][i], params[1][i], params[2][i], 
                    params[3][i], params[4][i], params[5][i]);
        }

        if (op->dev)
        {
            fprintf(fout, " %d   ", (int)orig_qvs[i]);
            if (iheight != NULL && iheight2 != NULL && ave_iheight != NULL)
                fprintf(fout, "%.2f\t%.2f\t%.2f\n",
                    iheight[i], iheight2[i], ave_iheight[i] );
            else if (iheight != NULL)
                fprintf(fout, "%.6f\n",
                    iheight[i]);
            else
                fprintf(fout, "\n");
        }
        else 
            fprintf(fout, "\n");
    }

    return SUCCESS;
}

/*******************************************************************************
 * Function: append_to_train_data
 *******************************************************************************
 */
int
append_to_train_data(FILE* fout, int num_bases, int *peak_locs, int NUM_PARAMS,
    double *params[NUM_PARAMS], uint8_t *orig_qvs, double *iheight, 
    double *iheight2, double *ave_iheight, Range clear_range, Range align_range,
    Align* vec_start, Align* best_alignment, Align* vec_end,
    int min_dpos, int max_dpos, int Verbose, Options *op, BtkMessage* message)
{
    int align_range_len;
    int r;
    int ins, del, sub, tot, alignment_size=best_alignment->trace_len;

    if (Verbose > 0)
        fprintf(stderr, "Clear range start :%4d end:%4d\n",
            clear_range.begin+1, clear_range.end+1);
        fprintf(fout, "# Clear range start :%4d end:%4d\n",
            clear_range.begin+1, clear_range.end+1);

    align_range_len = align_range.end - align_range.begin + 1;
    if (Verbose > 0)
        fprintf(stderr, "Align range start :%4d end:%4d\n",
            align_range.begin+1, align_range.end+1);
    if (Verbose > 0)
        fprintf(stderr, "Align range length:%4d\n", align_range_len);
    fprintf(fout, "# Align range start :%4d end:%4d\n",
            align_range.begin+1, align_range.end+1);
    fprintf(fout, "# Align range length:%4d\n", align_range_len);

    if (NUM_PARAMS == 4)
    {
        if (op->dev)
            fprintf(fout,
            "# Cons./Ref.\tMatch\tSample\t\tPredictors\t\t\t    Or.QV Iheights/Flows\n");
        else
            fprintf(fout,"# Cons./Ref.\tMatch\tSample\t\tPredictors\n");
    }
    else 
    {
        if (op->dev)
            fprintf(fout,
            "# Cons./Ref.\tMatch\tSample\t\tPredictors\t\t\t\t\t      Or.QV Iheights/Flows\n");
        else
            fprintf(fout,"# Cons./Ref.\tMatch\tSample\t\tPredictors\n");
    }

    if (OUTPUT_VECTOR_ALIGNMENT && (vec_start->score > 0) ) {
        fprintf(fout, "# Vector\n");
        append_align_to_train_data(fout, NUM_PARAMS, params, orig_qvs, iheight, 
            iheight2, ave_iheight, vec_start, vec_start->trace_dpos[0], 
            vec_start->trace_dpos[vec_start->trace_len-1], op, message);

        if (Verbose > 0) {
            r = align_fprint(stderr, vec_start, 70, message);
            if ( r == ERROR ) { return ERROR; }
        }
    }

    append_align_to_train_data(fout, NUM_PARAMS, params, orig_qvs, iheight, 
        iheight2, ave_iheight, best_alignment, min_dpos, max_dpos, op, message);

    if (Verbose > 0) {
        r = align_fprint(stderr, best_alignment, 70, message);
        if ( r == ERROR ) { return ERROR; }
    }

    if (OUTPUT_VECTOR_ALIGNMENT && (vec_end->score > 0) ) {
        fprintf(fout, "# Vector\n");
        append_align_to_train_data(fout, NUM_PARAMS, params, orig_qvs,
            iheight, iheight2, ave_iheight, vec_end, vec_end->trace_dpos[0],
            vec_end->trace_dpos[vec_end->trace_len-1], op, message);

        r = align_fprint(stderr, vec_end, 70, message);
        if ( r == ERROR ) { return ERROR; }
    }

    (void)get_num_errors(best_alignment, &ins, &del, &sub, &tot);
    fprintf(fout, "\n# Num errors: Ins=%d Del=%d Sub=%d Tot=%d\n\n\n", 
        ins, del, sub, tot); 
    fprintf(stderr, "\nNum errors: Ins=%d Del=%d Sub=%d Tot=%d\n",
        ins, del, sub, tot);
    fprintf(stderr, "Fraction of errors: sub=%f, ins=%f, del=%f tot=%f\n\n",
                (double)sub/(double)alignment_size,
                (double)ins/(double)alignment_size,
                (double)del/(double)alignment_size,
                (double)(sub+ins+del)/(double)alignment_size);

    return SUCCESS;
}

/********************************************************************************
 * This function appends one Align data structure to the training output
 * file.  Its synopsis is:
 *
 * result = append_align_to_trainphred_data(out, params, align, quality_values, 
 *          message)
 *
 * where
 *      out             is the (already open) output stream
 *      params          is the matrix of trace parameters
 *      align           is the address of the Align
 *      message         is the address of a BtkMessage where information about
 *                      an error will be put, if any
 *
 *      result          is 0 on success, !0 if an error occurs
 *
 ********************************************************************************/
int
append_align_to_trainphred_data(FILE* fout, Align* align, int min_dpos,
    int max_dpos, uint8_t* qualities, int num_alt_bases, char *bases2, 
    uint8_t *qvs2, int *inds2, BtkMessage* message)
{
    int     i, j, k, ismatch;

    for (j=0; j< align->trace_len; j++) {

        if ((align->trace_dpos[j] < min_dpos) ||
            (align->trace_dpos[j] > max_dpos))
 
            continue;

        ismatch = (align->trace_mchar[j]=='|')?1:0;

        fprintf(fout, "%d\t%c\t%d\t%d\t%c",
            align->trace_dpos[j]+1,
            align->trace_dchar[j],
            ismatch,
            align->trace_qpos[j]+1,
            align->trace_qchar[j]);

        if (align->trace_dir[j] == '2')
        {
            fprintf(fout, "\n");
            continue;
        }

        i = align->trace_qpos[j];
#if 0
        if ((i < clear_range.begin) || (i > clear_range.end))
        {
            fprintf(fout, "\n");
            continue;
        }
#endif
        fprintf(fout, "\t%d", qualities[i]);
        if (num_alt_bases > 0)
        {
            for (k=0; k<num_alt_bases; k++)
            {
                if (inds2[k] == align->trace_qpos[j])
                {
                    fprintf(fout, "\t%c\t%d", bases2[k], qvs2[k]);
                } 
            }
        } 
        fprintf(fout, "\n");   
    }

    return SUCCESS;
}

/*********************************************************************************
 * File: append_to_trainphred_data
 ********************************************************************************/
int
append_to_trainphred_data(FILE* fout, int num_bases, int *called_peak_locs,
    Range clear_range, Range align_range, Align* vec_start, Align* best_alignment,
    Align* vec_end, int min_dpos, int max_dpos, uint8_t *qualities, 
    int num_alt_bases, char *bases2, uint8_t *qvs2, int *inds2,
    BtkMessage* message )
{
    int align_range_len;
    int r;
    int ins, del, sub, tot, alignment_size=best_alignment->trace_len; 

    fprintf(stderr, "Clear range start :%4d end:%4d\n",
                    clear_range.begin+1, clear_range.end+1);
    fprintf(fout, "# Clear range start :%4d end:%4d\n",
                    clear_range.begin+1, clear_range.end+1);

    align_range_len = align_range.end - align_range.begin + 1;
    fprintf(stderr, "Align range start :%4d end:%4d\n",
                    align_range.begin+1, align_range.end+1);
    fprintf(stderr, "Align range length:%4d\n", align_range_len);
    fprintf(fout, "# Align range start :%4d end:%4d\n",
                    align_range.begin+1, align_range.end+1);
    fprintf(fout, "# Align range length:%4d\n", align_range_len);

    if (num_alt_bases == 0)
        fprintf(fout, "#\n# Ref. seq\tmatch\tRead_seq\tQV\n");
    else 
        fprintf(fout, "#\n# Ref. seq\tmatch\tRead_seq\tQV\tAlt. bases and QVs\n");

    if (OUTPUT_VECTOR_ALIGNMENT && (vec_start->score > 0) ) {
        fprintf(fout, "# Vector\n");
        append_align_to_trainphred_data(fout, vec_start, 
            vec_end->trace_dpos[0], vec_end->trace_dpos[vec_end->trace_len-1], 
            qualities, num_alt_bases, bases2, qvs2, inds2, message);

        r = align_fprint(stderr, vec_start, 70, message);
        if ( r == ERROR ) { return ERROR; }

    }

    append_align_to_trainphred_data(fout, best_alignment, 
        min_dpos, max_dpos, qualities, num_alt_bases, bases2, qvs2, inds2, 
        message);

        r = align_fprint(stderr, best_alignment, 70, message);
        if ( r == ERROR ) { return ERROR; }


    if (OUTPUT_VECTOR_ALIGNMENT && (vec_end->score > 0) ) {
        fprintf(fout, "# Vector\n");
        append_align_to_trainphred_data(fout, vec_end, 
            vec_end->trace_dpos[0], vec_end->trace_dpos[vec_end->trace_len-1], 
            qualities, num_alt_bases, bases2, qvs2, inds2, message);

        r = align_fprint(stderr, vec_end, 70, message);
        if ( r == ERROR ) { return ERROR; }

    }

    (void)get_num_errors(best_alignment, &ins, &del, &sub, &tot);
    fprintf(fout, "# Num errors: Ins=%d Del=%d Sub=%d Tot=%d\n\n",
        ins, del, sub, tot);
    fprintf(stderr, "\nNum errors: Ins=%d Del=%d Sub=%d Tot=%d\n",
        ins, del, sub, tot);
    fprintf(stderr, "Fraction of errors: sub=%f, ins=%f, del=%f tot=%f\n\n",
                (double)sub/(double)alignment_size,
                (double)ins/(double)alignment_size,
                (double)del/(double)alignment_size,
                (double)(sub+ins+del)/(double)alignment_size);

    return SUCCESS;
}


