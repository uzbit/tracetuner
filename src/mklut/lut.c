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

/**     Copyright (c) 1999 Paracel Inc.  All rights reserved.
 **/

/***************************************************************************
 * lut.c
 *
 *  $Id: lut.c,v 1.14 2009/01/16 15:19:51 gdenisov Exp $ 
 *
 * Purpose: 
 * - create lookup table for quality values 
 * - count the number of bases having each QV 0 to 100 and display percentages
 *   over QV 20, 30, 40, etc.
 *
 * Calls: 
 * - get_bases
 * - get_thresholds
 * - count_number_of_correct_bases_in_each_bin
 * - create_qv_table_via_dynamic_programming
 *
 **************************************************************************
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>

#include "select.h"
#include "params.h"
#include "lut.h"
#include "util.h"
#include "Btk_atod.h"
#include "get_thresholds.h"
#include "Btk_qv.h"
#include "Btk_lookup_table.h"
#include "check_data.h"
#include <unistd.h>

int          Verbose;          /* How much status info to print, if any */
int          Compress;         /* Whether to compress thresholds */

#define BUFLEN 1000
#define MAX_NUM_THRESHOLDS 100
#define MIN_DIST 0.0000001
#define DISPLAY_BASES 0
#define DISPLAY_THRESHOLDS 0
#define MIN_INCORRECT_COUNT 3

static char OutputName[BUFLEN];    /* Name of the Output lookup table file. */
static int  OutputSpecified;       /* Whether the user has specified a name. */
static int  PARAMETER_COUNT = 0;

static void 
show_usage(int argc, char *argv[])
{
    fprintf(stderr, 
    "\nVersion: %s\n"
    "usage: %s\n"
    "     [ -Q ] [ -V ]\n"
    "     [ -o <output_file> ]\n"
    "     <num_thresholds_1 num_thresholds_2 ... num_thresholds_n "
    "(n= 4 or 6)>  <  <alignment_file>\n"
      , TT_VERSION, argv[0]);
}

static void
show_usage_dev(int argc, char *argv[])
{
    fprintf(stderr,
    "\nVersion: %s\n"
    "usage: %s\n"
    "     [ -Q ] [ -V ]\n"
    "     [ -c ] [ -C ] \n"
    "     [ -b <initialbaseroom>]\n"
    "     [ -f <fileoffiles>]\n"
    "     [ -o <lookup_table_file> ]\n"
    "     <num_thresholds_1 num_thresholds_2 ... num_thresholds_n "
    "(n= 4 or 6)>  <  <alignment_file>\n"
      , TT_VERSION, argv[0]);
}


/***************************************************************************
 * get_threshold_index
 *
 * purpose:  for a given parameter, determine which of the <threshold_count>
 * <threshold>s have a value just greater than or equal to a certain <value>.
 *
 * called by: count_number_of_correct_bases_in_each_bin
 * calls: none
 *
 ***************************************************************************/
int
get_threshold_index(double value, double threshold[], int threshold_count)
{
    int i = threshold_count - 1; /* indexed 0 to threshold_count - 1 */

    if (value > threshold[i])
    {
        fprintf(stderr, "problem: a value=%f is greater than the max threshold[%d]=%f\n",
                value, i, threshold[i]);
        exit(-1);
    }

    while (i >= 0 && value <= threshold[i])
        i--;

    return i + 1;
}

/***************************************************************************
 * get_bin
 * <index> is an array of PARAMETER_COUNT indices
 *
 * purpose: get the address of a bin.
 *
 * called by:
 * count_number_of_correct_bases_in_each_bin
 * display_number_of_correct_bases_in_each_bin
 * number_in_cut
 *
 * calls:
 *
 ***************************************************************************/
BIN *
get_bin(BIN *bin, PARAMETER *parameter, int *index)
{
    unsigned long bin_number = 0;
    int i;

    for (i = 0; i < PARAMETER_COUNT; i++)
        bin_number += index[i] * parameter[i].dimension;

    return &(bin[bin_number]);
}

/***************************************************************************
 * get_cut
 *
 * Purpose:
 * get the address of a cut.
 *
 * called by:
 * number_in_cut
 * create_table_via_dynamic_programming
 * display_cubes
 *
 * calls:
 *
 * note that <l> is just used to check if we have a boundary condition.
 * otherwise <i,j,k> define the cube in the 4D space that we are
 * currently working on.  to save memory, we don't work with the whole
 * hypercubic (4d) space, just 2 cubic (3d) portions (a CURRENT and
 * a PREVIOUS) of it at a time.
 **************************************************************************
 */
CUT *
get_cut_4(int i, int j, int k, int l, TIME time, INFO *info)
{
    unsigned long cut_number;

    if (i == -1 || j == -1 || k == -1 || l == -1)
        return &(info->boundary_cut);

    cut_number = i * info->parameter[0].dimension
               + j * info->parameter[1].dimension
               + k * info->parameter[2].dimension;

    if (time == CURRENT)
        return &(info->current_cube[cut_number]);
    else if (time == PREVIOUS)
        return &(info->previous_cube[cut_number]);
    else
    {
        fprintf(stderr, "invalid time in function get_cut\n");
        exit(-1);
    }
}

CUT *
get_cut_6(int i, int j, int k, int l, int m, int n, TIME time, INFO *info)
{
    unsigned long cut_number;

    if (i == -1 || j == -1 || k == -1 || l == -1 || m == -1 || n == -1)
        return &(info->boundary_cut);

    cut_number = i * info->parameter[0].dimension
               + j * info->parameter[1].dimension
               + k * info->parameter[2].dimension
               + l * info->parameter[3].dimension
               + m * info->parameter[4].dimension;

    if (time == CURRENT)
        return &(info->current_cube[cut_number]);
    else if (time == PREVIOUS)
        return &(info->previous_cube[cut_number]);
    else
    {
        fprintf(stderr, "invalid time in function get_cut\n");
        exit(-1);
    }
}

/***************************************************************************
 *r initialize_highest_qv_cut
 *
 * purpose:  initialize the <initialize_highest_qv_cut>.
 *
 * called by: create_qv_table_via_dynamic_programming
 * calls: none
 *
 ***************************************************************************
 */
void
initialize_highest_qv_cut(HIGHEST_QV_CUT *highest_qv_cut)
{
    int m;

    highest_qv_cut->sum_of_indices = 0;

    for (m = 0; m < PARAMETER_COUNT; m++)
        highest_qv_cut->index[m] = -1;     /* the minimum real index is 0 */

    highest_qv_cut->correct_base_call_count = 0;
    highest_qv_cut->incorrect_base_call_count = 0;
    highest_qv_cut->total_base_call_count = 0;
    highest_qv_cut->error_rate = 0;
    highest_qv_cut->quality_value = 0;
}

/***************************************************************************
 * number_in_cut
 *
 * purpose: return, via <CUT>,
 * the number of correct and incorrect bases calls in
 * a cut defined by indices (i, j, k, l).
 * we need the array of <previous_highest_cut_parameter_index>s because
 * we use this function to also eliminate bases in bins that fell beneath
 * the cut from the highest qv cut from the last pass.
 *
 * called by: create_qv_table_via_dynamic_programming
 * calls: get_cut, get_bin
 *
 ***************************************************************************/
CUT *
number_in_cut_4(int i, int j, int k, int l, INFO *info,
              int *previous_highest_cut_parameter_index)
{
    int index[PARAMETER_COUNT];
    BIN *b_1111;

    CUT *c_0000;
    CUT *c_0001;
    CUT *c_0010;
    CUT *c_0011;
    CUT *c_0100;
    CUT *c_0101;
    CUT *c_0110;
    CUT *c_0111;
    CUT *c_1000;
    CUT *c_1001;
    CUT *c_1010;
    CUT *c_1011;
    CUT *c_1100;
    CUT *c_1101;
    CUT *c_1110;
    CUT *c_1111;

    index[0] = i;
    index[1] = j;
    index[2] = k;
    index[3] = l;

  /* if statements check if a cut represents a boundary condition
     in which the number of correct and incorrect calls are 0. */

  /*
   * Assume we're in the interior (which we are, most of the time), then
   * check and override for boundaries.
   */


  /*
   * Use a Gray-like coding here so that only one dimension changes from
   * line to line, so we can use pointer addition instead of array subscript
   * multiplication.
   *
  c_1111 = get_cut(i  ,j  ,k  ,l  ,CURRENT ,info);
  c_1101 = get_cut(i  ,j  ,k-1,l  ,CURRENT ,info);
  c_1001 = get_cut(i  ,j-1,k-1,l  ,CURRENT ,info);
  c_1011 = get_cut(i  ,j-1,k  ,l  ,CURRENT ,info);
  c_0011 = get_cut(i-1,j-1,k  ,l  ,CURRENT ,info);
  c_0111 = get_cut(i-1,j  ,k  ,l  ,CURRENT ,info);
  c_0101 = get_cut(i-1,j  ,k-1,l  ,CURRENT ,info);
  c_0001 = get_cut(i-1,j-1,k-1,l  ,CURRENT ,info);
  c_0000 = get_cut(i-1,j-1,k-1,l-1,PREVIOUS,info);
  c_0010 = get_cut(i-1,j-1,k  ,l-1,PREVIOUS,info);
  c_0110 = get_cut(i-1,j  ,k  ,l-1,PREVIOUS,info);
  c_0100 = get_cut(i-1,j  ,k-1,l-1,PREVIOUS,info);
  c_1100 = get_cut(i  ,j  ,k-1,l-1,PREVIOUS,info);
  c_1110 = get_cut(i  ,j  ,k  ,l-1,PREVIOUS,info);
  c_1010 = get_cut(i  ,j-1,k  ,l-1,PREVIOUS,info);
  c_1000 = get_cut(i  ,j-1,k-1,l-1,PREVIOUS,info);
   */

/* c_1111 always gives us a non-boundary cube, as all subscripts are >= 0 */

    c_1111 = get_cut_4(i, j, k, l, CURRENT, info);
    c_1101 = c_1111 - info->dimension[2];
    c_1001 = c_1101 - info->dimension[1];
    c_1011 = c_1001 + info->dimension[2];
    c_0011 = c_1011 - info->dimension[0];
    c_0111 = c_0011 + info->dimension[1];
    c_0101 = c_0111 - info->dimension[2];
    c_0001 = c_0101 - info->dimension[1];

/*
 * Use the offset of c_0001 from CURRENT as the offset of c_0000 from
 * PREVIOUS.  We can't easily start over from a get_cut() call, since we
 * can't be sure one of the subscripts isn't 0 at this point, making one
 * of the -1 values < 0 (boundary condition).
 */

    c_0000 = info->previous_cube + (c_0001 - info->current_cube);
    c_0010 = c_0000 + info->dimension[2];
    c_0110 = c_0010 + info->dimension[1];
    c_0100 = c_0110 - info->dimension[2];
    c_1100 = c_0100 + info->dimension[0];
    c_1110 = c_1100 + info->dimension[2];
    c_1010 = c_1110 - info->dimension[1];
    c_1000 = c_1010 - info->dimension[2];

    if (i == 0)
    {
        c_0000 = &info->boundary_cut;
        c_0001 = &info->boundary_cut;
        c_0010 = &info->boundary_cut;
        c_0011 = &info->boundary_cut;
        c_0100 = &info->boundary_cut;
        c_0101 = &info->boundary_cut;
        c_0110 = &info->boundary_cut;
        c_0111 = &info->boundary_cut;
    }

    if (j == 0)
    {
        c_0000 = &info->boundary_cut;
        c_0001 = &info->boundary_cut;
        c_0010 = &info->boundary_cut;
        c_0011 = &info->boundary_cut;
        c_1000 = &info->boundary_cut;
        c_1001 = &info->boundary_cut;
        c_1010 = &info->boundary_cut;
        c_1011 = &info->boundary_cut;
    }

    if (k == 0)
    {
        c_0000 = &info->boundary_cut;
        c_0001 = &info->boundary_cut;
        c_0100 = &info->boundary_cut;
        c_0101 = &info->boundary_cut;
        c_1000 = &info->boundary_cut;
        c_1001 = &info->boundary_cut;
        c_1100 = &info->boundary_cut;
        c_1101 = &info->boundary_cut;
    }

    if (l == 0)
    {
        c_0000 = &info->boundary_cut;
        c_0010 = &info->boundary_cut;
        c_0100 = &info->boundary_cut;
        c_0110 = &info->boundary_cut;
        c_1000 = &info->boundary_cut;
        c_1010 = &info->boundary_cut;
        c_1100 = &info->boundary_cut;
        c_1110 = &info->boundary_cut;
    }

    b_1111 = get_bin(info->bin, info->parameter, index);

    if (i <= previous_highest_cut_parameter_index[0] &&
        j <= previous_highest_cut_parameter_index[1] &&
        k <= previous_highest_cut_parameter_index[2] &&
        l <= previous_highest_cut_parameter_index[3])
    {
        b_1111->correct = b_1111->incorrect = 0;
    }

    /* the layout in space of the cuts are as follows--
       each row of comments corresponds to a row of equations */
    /* same cube (3d-space), same plane (2d-space):
       same bin, cut to left, cut to right, cut diagonal */
    /* same cube (3d-space), previous plane (2d-space):
       cut below, bottom left, bottom right, bottom diagonal */
    c_1111->correct =
     b_1111->correct+c_0111->correct+c_1011->correct-c_0011->correct
    +c_1101->correct-c_0101->correct-c_1001->correct+c_0001->correct
    +c_1110->correct-c_0110->correct-c_1010->correct+c_0010->correct
    -c_1100->correct+c_0100->correct+c_1000->correct-c_0000->correct;

    /* the layout in space of the cuts are as follows--
       each row of comments corresponds to a row of equations */
    /* same cube (3D-space), same plane (2D-space):
       same bin, cut to left, cut to right, cut diagonal */
    /* same cube (3D-space), previous plane (2D-space):
       cut below, bottom left, bottom right, bottom diagonal.
       the pattern (extensible to any number of dimensions) is:
       self (1111): +
       anything with 1 zero (1110, 1101, 1011, 0111): +;
       anything with 2 zeroes (1100, 1010, 1001, 0110, 0101, 0011): -;
       anything with 3 zeroes (1000, 0100, 0010, 0001): +;
       anything with 4 zeroes (0000): -;
       the idea is to add the perpendiculars, then subtract the diagonals. */
    c_1111->incorrect =
     b_1111->incorrect+c_0111->incorrect+c_1011->incorrect-c_0011->incorrect
    +c_1101->incorrect-c_0101->incorrect-c_1001->incorrect+c_0001->incorrect
    +c_1110->incorrect-c_0110->incorrect-c_1010->incorrect+c_0010->incorrect
    -c_1100->incorrect+c_0100->incorrect+c_1000->incorrect-c_0000->incorrect;

#if 0
      fprintf(stderr, "Cut(%d,%d,%d,%d) corr= %lu incorr= %lu\n",
          i,j,k,l, c_1111->correct, c_1111->incorrect);
#endif
    return c_1111;
}


CUT *
number_in_cut_6(int i, int j, int k, int l, int m, int n, INFO *info,
              int *previous_highest_cut_parameter_index)
{
    int index[PARAMETER_COUNT];
    BIN *b_111111;

    CUT *c_111111 = get_cut_6(i, j, k, l, m, n, CURRENT, info);
    CUT *c_111101 = c_111111 - info->dimension[4];
    CUT *c_111011 = c_111111 - info->dimension[3];
    CUT *c_110111 = c_111111 - info->dimension[2];
    CUT *c_101111 = c_111111 - info->dimension[1];
    CUT *c_011111 = c_111111 - info->dimension[0];
    CUT *c_111001 = c_111101 - info->dimension[3];
    CUT *c_110101 = c_111101 - info->dimension[2];
    CUT *c_101101 = c_111101 - info->dimension[1];
    CUT *c_011101 = c_111101 - info->dimension[0];
    CUT *c_110011 = c_111011 - info->dimension[2];
    CUT *c_101011 = c_111011 - info->dimension[1];
    CUT *c_011011 = c_111011 - info->dimension[0];
    CUT *c_100111 = c_110111 - info->dimension[1];
    CUT *c_010111 = c_110111 - info->dimension[0];
    CUT *c_001111 = c_101111 - info->dimension[0];
    CUT *c_110001 = c_111001 - info->dimension[2];
    CUT *c_101001 = c_111001 - info->dimension[1];
    CUT *c_100101 = c_101101 - info->dimension[2];
    CUT *c_100011 = c_110011 - info->dimension[1];
    CUT *c_011001 = c_111001 - info->dimension[0];
    CUT *c_010101 = c_110101 - info->dimension[0];
    CUT *c_010011 = c_110011 - info->dimension[0];
    CUT *c_001101 = c_001111 - info->dimension[4];
    CUT *c_001011 = c_001111 - info->dimension[3];
    CUT *c_000111 = c_001111 - info->dimension[2];
    CUT *c_100001 = c_110001 - info->dimension[1];
    CUT *c_010001 = c_110001 - info->dimension[0];
    CUT *c_001001 = c_011001 - info->dimension[1];
    CUT *c_000101 = c_000111 - info->dimension[4];
    CUT *c_000011 = c_000111 - info->dimension[3];
    CUT *c_000001 = c_000011 - info->dimension[4];

    CUT *c_000000 = info->previous_cube + (c_000001 - info->current_cube);
    CUT *c_000010 = c_000000 + info->dimension[4];
    CUT *c_000100 = c_000000 + info->dimension[3];
    CUT *c_000110 = c_000100 + info->dimension[4];
    CUT *c_001000 = c_000000 + info->dimension[2];
    CUT *c_001010 = c_001000 + info->dimension[4];
    CUT *c_001100 = c_001000 + info->dimension[3];
    CUT *c_001110 = c_001100 + info->dimension[4];
    CUT *c_010000 = c_000000 + info->dimension[1];
    CUT *c_010010 = c_010000 + info->dimension[4];
    CUT *c_010100 = c_010000 + info->dimension[3];
    CUT *c_010110 = c_010100 + info->dimension[4];
    CUT *c_011000 = c_001000 + info->dimension[1];
    CUT *c_011010 = c_011000 + info->dimension[4];
    CUT *c_011100 = c_011000 + info->dimension[3];
    CUT *c_011110 = c_011100 + info->dimension[4];
    CUT *c_100000 = c_000000 + info->dimension[0];
    CUT *c_100010 = c_100000 + info->dimension[4];
    CUT *c_100100 = c_100000 + info->dimension[3];
    CUT *c_100110 = c_100100 + info->dimension[4];
    CUT *c_101000 = c_100000 + info->dimension[2];
    CUT *c_101010 = c_101000 + info->dimension[4];
    CUT *c_101100 = c_101000 + info->dimension[3];
    CUT *c_101110 = c_101100 + info->dimension[4];
    CUT *c_110000 = c_100000 + info->dimension[1];
    CUT *c_110010 = c_110000 + info->dimension[4];
    CUT *c_110100 = c_110000 + info->dimension[3];
    CUT *c_110110 = c_110100 + info->dimension[4];
    CUT *c_111000 = c_110000 + info->dimension[2];
    CUT *c_111010 = c_111000 + info->dimension[4];
    CUT *c_111100 = c_111000 + info->dimension[3];
    CUT *c_111110 = c_111100 + info->dimension[4];

    index[0] = i;
    index[1] = j;
    index[2] = k;
    index[3] = l;
    index[4] = m;
    index[5] = n;

    if (i == 0)
    {
        c_000000 = &info->boundary_cut;
        c_000001 = &info->boundary_cut;
        c_000010 = &info->boundary_cut;
        c_000100 = &info->boundary_cut;
        c_001000 = &info->boundary_cut;
        c_010000 = &info->boundary_cut;
        c_011000 = &info->boundary_cut;
        c_010100 = &info->boundary_cut;
        c_010010 = &info->boundary_cut;
        c_010001 = &info->boundary_cut;
        c_001100 = &info->boundary_cut;
        c_001010 = &info->boundary_cut;
        c_001001 = &info->boundary_cut;
        c_000110 = &info->boundary_cut;
        c_000101 = &info->boundary_cut;
        c_000011 = &info->boundary_cut;
        c_011100 = &info->boundary_cut;
        c_011010 = &info->boundary_cut;
        c_010110 = &info->boundary_cut;
        c_001110 = &info->boundary_cut;
        c_011001 = &info->boundary_cut;
        c_010101 = &info->boundary_cut;
        c_001101 = &info->boundary_cut;
        c_010011 = &info->boundary_cut;
        c_001011 = &info->boundary_cut;
        c_000111 = &info->boundary_cut;
        c_001111 = &info->boundary_cut;
        c_010111 = &info->boundary_cut;
        c_011011 = &info->boundary_cut;
        c_011101 = &info->boundary_cut;
        c_011110 = &info->boundary_cut;
        c_011111 = &info->boundary_cut;
    }

    if (j == 0)
    {
        c_000000 = &info->boundary_cut;
        c_100000 = &info->boundary_cut;
        c_000001 = &info->boundary_cut;
        c_100001 = &info->boundary_cut;
        c_000010 = &info->boundary_cut;
        c_100010 = &info->boundary_cut;
        c_000100 = &info->boundary_cut;
        c_100100 = &info->boundary_cut;
        c_001000 = &info->boundary_cut;
        c_101000 = &info->boundary_cut;
        c_001100 = &info->boundary_cut;
        c_101100 = &info->boundary_cut;
        c_001010 = &info->boundary_cut;
        c_101010 = &info->boundary_cut;
        c_001001 = &info->boundary_cut;
        c_101001 = &info->boundary_cut;
        c_000110 = &info->boundary_cut;
        c_100110 = &info->boundary_cut;
        c_000101 = &info->boundary_cut;
        c_100101 = &info->boundary_cut;
        c_000011 = &info->boundary_cut;
        c_100011 = &info->boundary_cut;
        c_001110 = &info->boundary_cut;
        c_101110 = &info->boundary_cut;
        c_001101 = &info->boundary_cut;
        c_101101 = &info->boundary_cut;
        c_001011 = &info->boundary_cut;
        c_101011 = &info->boundary_cut;
        c_000111 = &info->boundary_cut;
        c_100111 = &info->boundary_cut;
        c_001111 = &info->boundary_cut;
        c_101111 = &info->boundary_cut;
    }

    if (k == 0)
    {
        c_000000 = &info->boundary_cut;
        c_100000 = &info->boundary_cut;
        c_000001 = &info->boundary_cut;
        c_100001 = &info->boundary_cut;
        c_000010 = &info->boundary_cut;
        c_100010 = &info->boundary_cut;
        c_000100 = &info->boundary_cut;
        c_100100 = &info->boundary_cut;
        c_010000 = &info->boundary_cut;
        c_110000 = &info->boundary_cut;
        c_010100 = &info->boundary_cut;
        c_110100 = &info->boundary_cut;
        c_010010 = &info->boundary_cut;
        c_110010 = &info->boundary_cut;
        c_010001 = &info->boundary_cut;
        c_110001 = &info->boundary_cut;
        c_000110 = &info->boundary_cut;
        c_100110 = &info->boundary_cut;
        c_000101 = &info->boundary_cut;
        c_100101 = &info->boundary_cut;
        c_000011 = &info->boundary_cut;
        c_100011 = &info->boundary_cut;
        c_010110 = &info->boundary_cut;
        c_110110 = &info->boundary_cut;
        c_010101 = &info->boundary_cut;
        c_110101 = &info->boundary_cut;
        c_010011 = &info->boundary_cut;
        c_110011 = &info->boundary_cut;
        c_000111 = &info->boundary_cut;
        c_100111 = &info->boundary_cut;
        c_010111 = &info->boundary_cut;
        c_110111 = &info->boundary_cut;
    }

    if (l == 0)
    {
        c_000000 = &info->boundary_cut;
        c_100000 = &info->boundary_cut;
        c_000001 = &info->boundary_cut;
        c_100001 = &info->boundary_cut;
        c_000010 = &info->boundary_cut;
        c_100010 = &info->boundary_cut;
        c_001000 = &info->boundary_cut;
        c_101000 = &info->boundary_cut;
        c_010000 = &info->boundary_cut;
        c_110000 = &info->boundary_cut;
        c_011000 = &info->boundary_cut;
        c_111000 = &info->boundary_cut;
        c_010010 = &info->boundary_cut;
        c_110010 = &info->boundary_cut;
        c_010001 = &info->boundary_cut;
        c_110001 = &info->boundary_cut;
        c_001010 = &info->boundary_cut;
        c_101010 = &info->boundary_cut;
        c_001001 = &info->boundary_cut;
        c_101001 = &info->boundary_cut;
        c_000011 = &info->boundary_cut;
        c_100011 = &info->boundary_cut;
        c_011010 = &info->boundary_cut;
        c_111010 = &info->boundary_cut;
        c_011001 = &info->boundary_cut;
        c_111001 = &info->boundary_cut;
        c_010011 = &info->boundary_cut;
        c_110011 = &info->boundary_cut;
        c_001011 = &info->boundary_cut;
        c_101011 = &info->boundary_cut;
        c_011011 = &info->boundary_cut;
        c_111011 = &info->boundary_cut;
    }

    if (m == 0)
    {
        c_000000 = &info->boundary_cut;
        c_100000 = &info->boundary_cut;
        c_000001 = &info->boundary_cut;
        c_100001 = &info->boundary_cut;
        c_000100 = &info->boundary_cut;
        c_100100 = &info->boundary_cut;
        c_001000 = &info->boundary_cut;
        c_101000 = &info->boundary_cut;
        c_010000 = &info->boundary_cut;
        c_110000 = &info->boundary_cut;
        c_011000 = &info->boundary_cut;
        c_111000 = &info->boundary_cut;
        c_010100 = &info->boundary_cut;
        c_110100 = &info->boundary_cut;
        c_010001 = &info->boundary_cut;
        c_110001 = &info->boundary_cut;
        c_001100 = &info->boundary_cut;
        c_101100 = &info->boundary_cut;
        c_001001 = &info->boundary_cut;
        c_101001 = &info->boundary_cut;
        c_000101 = &info->boundary_cut;
        c_100101 = &info->boundary_cut;
        c_011100 = &info->boundary_cut;
        c_111100 = &info->boundary_cut;
        c_011001 = &info->boundary_cut;
        c_111001 = &info->boundary_cut;
        c_010101 = &info->boundary_cut;
        c_110101 = &info->boundary_cut;
        c_001101 = &info->boundary_cut;
        c_101101 = &info->boundary_cut;
        c_011101 = &info->boundary_cut;
        c_111101 = &info->boundary_cut;
    }

    if (n == 0)
    {
        c_000000 = &info->boundary_cut;
        c_100000 = &info->boundary_cut;
        c_000010 = &info->boundary_cut;
        c_100010 = &info->boundary_cut;
        c_000100 = &info->boundary_cut;
        c_100100 = &info->boundary_cut;
        c_001000 = &info->boundary_cut;
        c_101000 = &info->boundary_cut;
        c_010000 = &info->boundary_cut;
        c_110000 = &info->boundary_cut;
        c_011000 = &info->boundary_cut;
        c_111000 = &info->boundary_cut;
        c_010100 = &info->boundary_cut;
        c_110100 = &info->boundary_cut;
        c_010010 = &info->boundary_cut;
        c_110010 = &info->boundary_cut;
        c_001100 = &info->boundary_cut;
        c_101100 = &info->boundary_cut;
        c_001010 = &info->boundary_cut;
        c_101010 = &info->boundary_cut;
        c_000110 = &info->boundary_cut;
        c_100110 = &info->boundary_cut;
        c_011100 = &info->boundary_cut;
        c_111100 = &info->boundary_cut;
        c_011010 = &info->boundary_cut;
        c_111010 = &info->boundary_cut;
        c_010110 = &info->boundary_cut;
        c_110110 = &info->boundary_cut;
        c_001110 = &info->boundary_cut;
        c_101110 = &info->boundary_cut;
        c_011110 = &info->boundary_cut;
        c_111110 = &info->boundary_cut;
    }

    b_111111 = get_bin(info->bin, info->parameter, index);

    if (i <= previous_highest_cut_parameter_index[0] &&
        j <= previous_highest_cut_parameter_index[1] &&
        k <= previous_highest_cut_parameter_index[2] &&
        l <= previous_highest_cut_parameter_index[3] &&
        m <= previous_highest_cut_parameter_index[4] &&
        n <= previous_highest_cut_parameter_index[5])
    {
        b_111111->correct = b_111111->incorrect = 0;
    }

    /* the layout in space of the cuts are as follows--
       each row of comments corresponds to a row of equations */
    /* same cube (3d-space), same plane (2d-space):
       same bin, cut to left, cut to right, cut diagonal */
    /* same cube (3d-space), previous plane (2d-space):
       cut below, bottom left, bottom right, bottom diagonal */

    c_111111->correct = 
    
        b_111111->correct

      + c_011111->correct
      + c_101111->correct
      + c_110111->correct
      + c_111011->correct
      + c_111101->correct
      + c_111110->correct

      - c_001111->correct
      - c_010111->correct
      - c_011011->correct
      - c_011101->correct
      - c_011110->correct
      - c_100111->correct
      - c_101011->correct
      - c_101101->correct
      - c_101110->correct
      - c_110011->correct
      - c_110101->correct
      - c_110110->correct
      - c_111001->correct
      - c_111010->correct
      - c_111100->correct

      + c_000111->correct
      + c_001011->correct
      + c_001101->correct
      + c_001110->correct
      + c_010011->correct
      + c_010101->correct
      + c_010110->correct
      + c_011001->correct
      + c_011010->correct
      + c_011100->correct
      + c_100011->correct
      + c_100101->correct
      + c_100110->correct
      + c_101001->correct
      + c_101010->correct
      + c_101100->correct
      + c_110001->correct
      + c_110010->correct
      + c_110100->correct
      + c_111000->correct

      - c_110000->correct
      - c_101000->correct
      - c_100100->correct
      - c_100010->correct
      - c_100001->correct
      - c_011000->correct
      - c_010100->correct
      - c_010010->correct
      - c_010001->correct
      - c_001100->correct
      - c_001010->correct
      - c_001001->correct
      - c_000110->correct
      - c_000101->correct
      - c_000011->correct

      + c_100000->correct
      + c_010000->correct
      + c_001000->correct
      + c_000100->correct
      + c_000010->correct
      + c_000001->correct

      - c_000000->correct;

    /* the layout in space of the cuts are as follows--
       each row of comments corresponds to a row of equations */
    /* same cube (3D-space), same plane (2D-space):
       same bin, cut to left, cut to right, cut diagonal */
    /* same cube (3D-space), previous plane (2D-space):
       cut below, bottom left, bottom right, bottom diagonal.
       the pattern (extensible to any number of dimensions) is:
       self (1111): +
       anything with 1 zero (1110, 1101, 1011, 0111): +;
       anything with 2 zeroes (1100, 1010, 1001, 0110, 0101, 0011): -;
       anything with 3 zeroes (1000, 0100, 0010, 0001): +;
       anything with 4 zeroes (0000): -;
       the idea is to add the perpendiculars, then subtract the diagonals. */
    c_111111->incorrect = 

        b_111111->incorrect

      + c_011111->incorrect
      + c_101111->incorrect
      + c_110111->incorrect
      + c_111011->incorrect
      + c_111101->incorrect
      + c_111110->incorrect

      - c_001111->incorrect
      - c_010111->incorrect
      - c_011011->incorrect
      - c_011101->incorrect
      - c_011110->incorrect
      - c_100111->incorrect
      - c_101011->incorrect
      - c_101101->incorrect
      - c_101110->incorrect
      - c_110011->incorrect
      - c_110101->incorrect
      - c_110110->incorrect
      - c_111001->incorrect
      - c_111010->incorrect
      - c_111100->incorrect

      + c_000111->incorrect
      + c_001011->incorrect
      + c_001101->incorrect
      + c_001110->incorrect
      + c_010011->incorrect
      + c_010101->incorrect
      + c_010110->incorrect
      + c_011001->incorrect
      + c_011010->incorrect
      + c_011100->incorrect
      + c_100011->incorrect
      + c_100101->incorrect
      + c_100110->incorrect
      + c_101001->incorrect
      + c_101010->incorrect
      + c_101100->incorrect
      + c_110001->incorrect
      + c_110010->incorrect
      + c_110100->incorrect
      + c_111000->incorrect

      - c_110000->incorrect
      - c_101000->incorrect
      - c_100100->incorrect
      - c_100010->incorrect
      - c_100001->incorrect
      - c_011000->incorrect
      - c_010100->incorrect
      - c_010010->incorrect
      - c_010001->incorrect
      - c_001100->incorrect
      - c_001010->incorrect
      - c_001001->incorrect
      - c_000110->incorrect
      - c_000101->incorrect
      - c_000011->incorrect

      + c_100000->incorrect
      + c_010000->incorrect
      + c_001000->incorrect
      + c_000100->incorrect
      + c_000010->incorrect
      + c_000001->incorrect

      - c_000000->incorrect;
#if 0
      fprintf(stderr, "Cut(%d,%d,%d,%d,%d,%d) corr= %lu incorr= %lu\n", 
          i,j,k,l,m,n, c_111111->correct, c_111111->incorrect);
#endif
    return c_111111;
}

/**************************************************************************
 * Function: read_bases_and_populate_bins
 **************************************************************************
 */
int 
read_bases_and_populate_bins(char *InputName, 
    unsigned long *base_count, PARAMETER parameter[], BIN *bin)
{
    int           i, read_count=0;
#if DISPLAY_THRESHOLDS
    int           j;
#endif
    int           linenum = 0, index[PARAMETER_COUNT];
    char          train_name[BUFLEN], *s, *current;
    int           is_match;
    int           spos, cpos;   /* Sample and consensus positions. */
    char          schar, cchar; /* Sample and consensus characters. */
    char          buffer[BUFLEN];
    FILE         *fileoffiles, *trainfile;
    BASE          base;
    BIN          *current_bin;

    if ((fileoffiles=fopen(InputName,"r"))== NULL) {
        fprintf(stderr, "Unable to open file_of_files '%s'\n",
            InputName);
        return ERROR;
    }

   *base_count = 0;
    current = buffer;

#if DISPLAY_THRESHOLDS
    for (i = 0; i < PARAMETER_COUNT; i++) {
        fprintf(stderr, "\nParameter=%d, thresholds: \n", i);
        for (j=0; j<parameter[i].threshold_count; j++) {
            fprintf(stderr, "%f\n", parameter[i].threshold[j]);
        }
    }
#endif
    while (fgets(current, BUFLEN, fileoffiles) != NULL)
    {
        /* Ignore all white space lines and comments */
        if (strspn(current, " \t\r\n") == strlen(current)) {
            continue;
        }
        if ((current[0] == '#')
            || ((current[0] == '/') && (current[1] == '*'))
            || (current[0] == ';'))
        {
            continue;
        }

        /* Open an individual train file */
        sscanf(current, "%s", train_name);
        if ((trainfile=fopen(train_name,"r"))== NULL) {
            fprintf(stderr, "Unable to open train file '%s'\n",
                InputName);
            return ERROR;
        }

        linenum = 0;
        while (fgets(train_name, sizeof(train_name), trainfile) != NULL)
        {
            linenum++;

/* Ignore all white space lines and comments */
            if (strspn(train_name, " \t\r\n") == strlen(train_name))
                continue;

            if (train_name[0] == '#' || train_name[0] == ';' ||
               (train_name[0] == '/' && train_name[1] == '*'))
                continue;

/* Get rid of 1) sample position, 2) sample base,
 * 3) consensus position and 4) consensus base. */

            s = strtok(train_name, " \t\n");
            if (s == NULL) continue;
            cpos = atoi(s);

            s = strtok(NULL, " \t\n");
            if (s == NULL) continue;
            cchar = s[0];

            s = strtok(NULL, " \t\n");
            if (s == NULL) continue;
            is_match = atoi(s);

            s = strtok(NULL, " \t\n");
            if (s == NULL) continue;
            spos = atoi(s);

            s = strtok(NULL, " \t\n");
            if (s == NULL) continue;
            schar = s[0];

            s += strlen(s) + 1;

/* If sample == '-', there are no training parameters. */

            if (schar == '-') continue;

            if (Btk_atod(&s, &base.parameter[0]) != 1
            ||  Btk_atod(&s, &base.parameter[1]) != 1
            ||  Btk_atod(&s, &base.parameter[2]) != 1
            ||  Btk_atod(&s, &base.parameter[3]) != 1)
            {
                fprintf(stderr, 
                "train file %s, line %d:\n%s\nmissing/garbled base; skipping\n",
                        train_name, linenum, current);
                continue;
            }

           (*base_count)++;
            read_count++;
            if ((read_count % 100000) == 0) {
                fprintf(stderr, "\r   %d bases populated", read_count);
            }

            /* Populate bins with current base */
            for (i = 0; i < PARAMETER_COUNT; i++) {
                index[i] = get_threshold_index(base.parameter[i],
                    parameter[i].threshold, parameter[i].threshold_count);
            }
            current_bin = get_bin(bin, parameter, index);
            if (is_match)
                current_bin->correct++;
            else
                current_bin->incorrect++;  

        }   /* loop in lines of a particular trainfile */
        fclose(trainfile);
    }       /* loop in lines of fileoffiles */
    fprintf(stderr, "\n   %d bases populated\n", read_count);
    fclose(fileoffiles);

    return SUCCESS;
}

/**************************************************************************
 * count_number_of_correct_bases_in_each_bin
 * 
 * purpose: record in each bin, the number of correct and incorrect bases.
 *                                                                          
 * called by: main
 * calls:
 * get_threshold_index
 * get_bin
 *
 **************************************************************************
 */
void
count_number_of_correct_bases_in_each_bin(BASE *base_array, 
    unsigned long base_count, PARAMETER *parameter, BIN *bin, FILE *fout)
{
    unsigned long n, bin_count;
    int *index;
    BASE *base;
    int i;
    BIN *current_bin;

/*
 *  bin_count = 1;
 *  for (i=0; i< PARAMETER_COUNT; i++) {
 *      bin_count *= parameter[i].threshold_count;
 *  }
 */

    bin_count = parameter[PARAMETER_COUNT - 1].threshold_count
              * parameter[PARAMETER_COUNT - 1].dimension;

    index = (int *) malloc(sizeof(int) * PARAMETER_COUNT);
    if (index == NULL)
    {
        fprintf(fout, "couldn't malloc index in count_number_of_correct_bases_in_bin\n");
        exit(-1);
    }

    for (n = 0, current_bin = bin; n < bin_count; n++, current_bin++) {
        current_bin->correct = current_bin->incorrect = 0;
    }

    for (n = 0, base = base_array; n < base_count; n++, base++)
    {
        for (i = 0; i < PARAMETER_COUNT; i++)
            index[i] = get_threshold_index(base->parameter[i], 
                                           parameter[i].threshold, 
                                           parameter[i].threshold_count);

        current_bin = get_bin(bin, parameter, index);

        if (base->is_match)
            current_bin->correct++;
        else
            current_bin->incorrect++;
    }

    FREE(index);
}

/***************************************************************************
 * update_highest_qv_cut
 *
 * purpose: check if the <correct_base_call_count> and
 * <incorrect_base_call_count> defined by the latest cut imply
 * that this latest cut is better than the <highest_qv_cut> so far.
 * if so, update the <highest_qv_cut>.
 * <i>, <j>, <k>, <l> are indices the the threshold values (e.g., 0,..,49)
 * that define the <parameter> cuts.
 *
 * called by: create_qv_table_via_dynamic_programming
 * calls: none
 *
 ***************************************************************************/
void
update_highest_qv_cut_4(
    HIGHEST_QV_CUT *highest_qv_cut, 
    PARAMETER *parameter,
    unsigned long correct_base_call_count,
    unsigned long incorrect_base_call_count, 
    int i, 
    int j,
    int k, 
    int l, 
    int *previous_highest_cut_parameter_index)
{
    unsigned long total_base_call_count;
    double error_rate;
    int quality_value, sum_of_indices;

    total_base_call_count = correct_base_call_count + incorrect_base_call_count;

/* the error rate includes a penalty for small sample size
   by adding 1 in numerator and denominator.
   hence, 1 correct and 1 incorrect base is assigned an
   error rate of 2/3, which is equivalent to the
   error rate for 66 incorrect bases and 33 correct bases. */

    if (incorrect_base_call_count == 0)
        error_rate = ((double) (1 + incorrect_base_call_count)) /
                     ((double) (1 + total_base_call_count));
    else
        error_rate = ((double) incorrect_base_call_count) /
                     ((double) total_base_call_count);

    quality_value = (int) rint(-10 * log10(error_rate));
    sum_of_indices = i + j + k + l;

    if (((incorrect_base_call_count >= MIN_INCORRECT_COUNT) &&
         (quality_value > highest_qv_cut->quality_value))
         ||
        ((incorrect_base_call_count >= MIN_INCORRECT_COUNT) &&
         (quality_value == highest_qv_cut->quality_value)
                                                   &&
         (total_base_call_count > highest_qv_cut->total_base_call_count))
         ||
        ((incorrect_base_call_count >= MIN_INCORRECT_COUNT) &&
         (quality_value == highest_qv_cut->quality_value)
                                                   &&
         (total_base_call_count == highest_qv_cut->total_base_call_count)
                                                   &&
         (sum_of_indices > highest_qv_cut->sum_of_indices)))
    {
        highest_qv_cut->sum_of_indices = sum_of_indices;
        highest_qv_cut->index[0] = i;
        highest_qv_cut->index[1] = j;
        highest_qv_cut->index[2] = k;
        highest_qv_cut->index[3] = l;
        highest_qv_cut->correct_base_call_count = correct_base_call_count;
        highest_qv_cut->incorrect_base_call_count = incorrect_base_call_count;
        highest_qv_cut->total_base_call_count = total_base_call_count;
        highest_qv_cut->error_rate = error_rate;
        highest_qv_cut->quality_value = quality_value;
        highest_qv_cut->parameter[0] = parameter[0].threshold[i];
        highest_qv_cut->parameter[1] = parameter[1].threshold[j];
        highest_qv_cut->parameter[2] = parameter[2].threshold[k];
        highest_qv_cut->parameter[3] = parameter[3].threshold[l];
    }
}

void
update_highest_qv_cut_6(
    HIGHEST_QV_CUT *highest_qv_cut, 
    PARAMETER *parameter,
    unsigned long correct_base_call_count,
    unsigned long incorrect_base_call_count,
    int i,
    int j,
    int k,
    int l,
    int m,
    int n,
    int *previous_highest_cut_parameter_index)
{
    unsigned long total_base_call_count;
    double error_rate;
    int quality_value, sum_of_indices;

    total_base_call_count = correct_base_call_count + incorrect_base_call_count;

/* the error rate includes a penalty for small sample size
   by adding 1 in numerator and denominator.
   hence, 1 correct and 1 incorrect base is assigned an
   error rate of 2/3, which is equivalent to the
   error rate for 66 incorrect bases and 33 correct bases. */

    if (incorrect_base_call_count == 0)
        error_rate = ((double) (1 + incorrect_base_call_count)) /
                     ((double) (1 + total_base_call_count));
    else
        error_rate = ((double) incorrect_base_call_count) /
                     ((double) total_base_call_count);

    quality_value = (int) rint(-10 * log10(error_rate));
    sum_of_indices = i + j + k + l + m + n;


    if (((incorrect_base_call_count >= MIN_INCORRECT_COUNT) &&
         (quality_value > highest_qv_cut->quality_value))
         ||
        ((incorrect_base_call_count >= MIN_INCORRECT_COUNT) &&
         (quality_value == highest_qv_cut->quality_value)
                                                   &&
         (total_base_call_count > highest_qv_cut->total_base_call_count))
         ||
        ((incorrect_base_call_count >= MIN_INCORRECT_COUNT) &&
         (quality_value == highest_qv_cut->quality_value)
                                                   &&
         (total_base_call_count == highest_qv_cut->total_base_call_count)
                                                   &&
         (sum_of_indices > highest_qv_cut->sum_of_indices)))
    {
        highest_qv_cut->sum_of_indices = sum_of_indices;
        highest_qv_cut->index[0] = i;
        highest_qv_cut->index[1] = j;
        highest_qv_cut->index[2] = k;
        highest_qv_cut->index[3] = l;
        highest_qv_cut->index[4] = m;
        highest_qv_cut->index[5] = n;
        highest_qv_cut->correct_base_call_count = correct_base_call_count;
        highest_qv_cut->incorrect_base_call_count = incorrect_base_call_count;
        highest_qv_cut->total_base_call_count = total_base_call_count;
        highest_qv_cut->error_rate = error_rate;
        highest_qv_cut->quality_value = quality_value;
        highest_qv_cut->parameter[0] = parameter[0].threshold[i];
        highest_qv_cut->parameter[1] = parameter[1].threshold[j];
        highest_qv_cut->parameter[2] = parameter[2].threshold[k];
        highest_qv_cut->parameter[3] = parameter[3].threshold[l];
        highest_qv_cut->parameter[4] = parameter[4].threshold[m];
        highest_qv_cut->parameter[5] = parameter[5].threshold[n];
    }
}

/***************************************************************************
 * write_to_qv_table
 *
 * purpose: write the latest <highest_qv_cut> in a "C" programming language
 * style if statement that checks if each of the PARAMETER_COUNT parameters
 * falls below the values specified in the if statement.
 * this table can be used in a "C" function that assigns quality values.
 *
 * called by: create_qv_table_via_dynamic_programming
 * calls: none
 *
 ***************************************************************************/

void
write_to_qv_table(HIGHEST_QV_CUT *highest_qv_cut, FILE *fout, 
    int parameter_count)
{
    static int table_entry_number = 0;

    table_entry_number++;

    if (parameter_count == 4)
        fprintf(fout, "   %d   %d  %d  %d  %d \n",
           highest_qv_cut->quality_value,
           highest_qv_cut->index[0],
           highest_qv_cut->index[1],
           highest_qv_cut->index[2],
           highest_qv_cut->index[3]);
    else /* parameter_count == 6 */
        fprintf(fout, "   %d   %d  %d  %d  %d %d  %d \n",
           highest_qv_cut->quality_value,
           highest_qv_cut->index[0],
           highest_qv_cut->index[1],
           highest_qv_cut->index[2],
           highest_qv_cut->index[3],
           highest_qv_cut->index[4],
           highest_qv_cut->index[5]);
}

/***************************************************************************
 * create_qv_table_via_dynamic_programming
 *
 * called by: main
 * calls: 
 * initialize_highest_qv_cut
 * get_cut
 * number_in_cut
 * write_to_qv_table
 * 
 ***************************************************************************/
/* a cut is a set of thresholds */
int 
create_qv_table_via_dynamic_programming(BIN *bin, PARAMETER parameter[],
    int parameter_count, unsigned long base_count, FILE *fout, 
    BtkMessage *message)
{
    unsigned long countdown = base_count;
    HIGHEST_QV_CUT highest_qv_cut;
    int i, j, k, l, m, n, num_entries = 0;
    INFO info;
    CUT *cut, *temp_cube;
    int *qv_counter, *qv_decade_counter;

    int *previous_highest_cut_parameter_index;
/* keep track of threshold indices from previous highest_qv_cut added to
   lookup table for the purpose of zeroing out the number of correct and
   incorrect calls in bins less than or equal to these indices */

    qv_counter        = CALLOC(int, MAX_QV);
    qv_decade_counter = CALLOC(int, MAX_QV / 10);

    fprintf(fout, "\n#  Quality value and parameter threshold indexes:\n");

/* for each cut (set of threshold values), count the number of correct and
   incorrect base calls in the data set to determine the quality value
   associated with that cut; decide which cut has the highest quality value 
   and enter that cut into a table; repeat the procedure with the remaining
   data (the data above the last set of thresholds). */

    info.parameter = parameter;
    info.boundary_cut.correct = info.boundary_cut.incorrect = 0;
    info.bin = bin;
    if (parameter_count == 4)
        info.previous_cube = CALLOC(CUT, parameter[0].threshold_count
                                       * parameter[1].threshold_count
                                       * parameter[2].threshold_count); 
    else /* parameter_count == 6 */
        info.previous_cube = CALLOC(CUT, parameter[0].threshold_count
                                       * parameter[1].threshold_count
                                       * parameter[2].threshold_count
                                       * parameter[3].threshold_count
                                       * parameter[4].threshold_count);
    MEM_ERROR(info.previous_cube);

    if (parameter_count == 4)
        info.current_cube  = CALLOC(CUT, parameter[0].threshold_count
                                       * parameter[1].threshold_count
                                       * parameter[2].threshold_count); 

    else /* parameter_count == 6 */
        info.current_cube =  CALLOC(CUT, parameter[0].threshold_count
                                       * parameter[1].threshold_count
                                       * parameter[2].threshold_count
                                       * parameter[3].threshold_count
                                       * parameter[4].threshold_count);
    MEM_ERROR(info.current_cube);

    for (i=0; i<parameter_count; i++)
    {
        info.dimension[i] = parameter[i].dimension;
    }

    /* Initialize base counts of all cuts to zero */
    if (parameter_count == 4)
    {
        for (k = 0; k < parameter[2].threshold_count; k++)
            for (j = 0; j < parameter[1].threshold_count; j++)
                for (i = 0; i < parameter[0].threshold_count; i++)
                {
                    cut = get_cut_4(i, j, k, 0, PREVIOUS, &info);
                    cut->correct = cut->incorrect = 0;
                    cut = get_cut_4(i, j, k, 0, CURRENT, &info);
                    cut->correct = cut->incorrect = 0;
                }
    }
    else /* parameter_count == 6 */
    {  
        for (m=0; m< parameter[4].threshold_count; m++)
            for (l=0; l< parameter[3].threshold_count; l++)
                for (k = 0; k < parameter[2].threshold_count; k++)
                    for (j = 0; j < parameter[1].threshold_count; j++)
                        for (i = 0; i < parameter[0].threshold_count; i++)
                        {
                            cut = get_cut_6(i, j, k, l, m, 0, PREVIOUS, &info);
                            cut->correct = cut->incorrect = 0;
                            cut = get_cut_6(i, j, k, l, m, 0, CURRENT, &info);
                            cut->correct = cut->incorrect = 0;
                        }
    }
    previous_highest_cut_parameter_index = CALLOC(int, PARAMETER_COUNT);
    MEM_ERROR(previous_highest_cut_parameter_index);

    /* Initialize the real minimum value 0 */
    for (m = 0; m < PARAMETER_COUNT; m++)
        previous_highest_cut_parameter_index[m] = -1; 

    do  /* Dynamic prorgamming */
    {
        initialize_highest_qv_cut(&highest_qv_cut);

        if (parameter_count == 4)
        {
            for (l = 0; l < parameter[3].threshold_count; l++)
            {
              temp_cube = info.previous_cube;      
              info.previous_cube = info.current_cube;
              info.current_cube = temp_cube;

              for (k = 0; k < parameter[2].threshold_count; k++)
                for (j = 0; j < parameter[1].threshold_count; j++)
                  for (i = 0; i < parameter[0].threshold_count; i++)
                  {
                      /* The key step 
                       ***************************************
                       */
                      cut = number_in_cut_4(i, j, k, l, &info,
                          previous_highest_cut_parameter_index);
                      update_highest_qv_cut_4(&highest_qv_cut, 
                          parameter, cut->correct,
                          cut->incorrect, i, j, k, l,
                          previous_highest_cut_parameter_index);
                  }
            }
        }
        else /* parameter_count == 6 */
        {
            for (n = 0; n < parameter[5].threshold_count; n++)
            {
              temp_cube = info.previous_cube;
              info.previous_cube = info.current_cube;
              info.current_cube = temp_cube;
              for (m = 0; m < parameter[4].threshold_count; m++)
                  for (l = 0; l < parameter[3].threshold_count; l++)
                      for (k = 0; k < parameter[2].threshold_count; k++)
                          for (j = 0; j < parameter[1].threshold_count; j++)
                              for (i = 0; i < parameter[0].threshold_count; i++)
                              {
                                  /* The key step
                                   ***************************************
                                   */
                                  cut = number_in_cut_6(i, j, k, l, m, n, &info,
                                      previous_highest_cut_parameter_index);
                                  update_highest_qv_cut_6(&highest_qv_cut,
                                      parameter, cut->correct,
                                      cut->incorrect, i, j, k, l, m, n, 
                                      previous_highest_cut_parameter_index);
                              }
            }
        }

        if (highest_qv_cut.total_base_call_count != 0)
        {
            countdown -= highest_qv_cut.total_base_call_count;
            write_to_qv_table(&highest_qv_cut, fout, parameter_count);

            for (m = 0; m < PARAMETER_COUNT; m++)
                previous_highest_cut_parameter_index[m] = highest_qv_cut.index[m];

            num_entries++;

            if (Verbose)
                fprintf(stderr, "\r%lu bases to go (%d entries so far)...        ",
                        countdown, num_entries);
        }

        /* Keep track of how many bases have certain qv's */
        qv_counter[highest_qv_cut.quality_value] += 
                                    highest_qv_cut.total_base_call_count;

        if (highest_qv_cut.quality_value >= 10)
            qv_decade_counter[1] += highest_qv_cut.total_base_call_count;
        if (highest_qv_cut.quality_value >= 20)
            qv_decade_counter[2] += highest_qv_cut.total_base_call_count;
        if (highest_qv_cut.quality_value >= 30)
            qv_decade_counter[3] += highest_qv_cut.total_base_call_count;
        if (highest_qv_cut.quality_value >= 40)
            qv_decade_counter[4] += highest_qv_cut.total_base_call_count;
        if (highest_qv_cut.quality_value >= 50)
            qv_decade_counter[5] += highest_qv_cut.total_base_call_count;
        if (highest_qv_cut.quality_value >= 60)
            qv_decade_counter[6] += highest_qv_cut.total_base_call_count;
        if (highest_qv_cut.quality_value >= 70)
            qv_decade_counter[7] += highest_qv_cut.total_base_call_count;
        if (highest_qv_cut.quality_value >= 80)
            qv_decade_counter[8] += highest_qv_cut.total_base_call_count;
        if (highest_qv_cut.quality_value >= 90)
            qv_decade_counter[9] += highest_qv_cut.total_base_call_count;

    } while (highest_qv_cut.total_base_call_count != 0);

    if (Verbose)
        fprintf(stderr, "\rdone.                                      \n");

    if (countdown != 0)
        fprintf(fout, "/* warning: %lu base calls were unaccounted for */\n",
               countdown);


    fprintf(fout, "\n#  Number and percentage of bases having each quality value\n");

    for (i = 0; i < MAX_QV; i++)
        if (qv_counter[i] != 0)
            fprintf(fout, "#   %2d \t %9d \t %5.2f%%\n", i, qv_counter[i], 
                   100 * (double) qv_counter[i] / (double) base_count);

    fprintf(fout, "\n#  Number and percentage of bases having >= quality value\n");

    for (i = 1; i < (MAX_QV / 10); i++)
        if (qv_decade_counter[i] != 0)
            fprintf(fout, "#   %2d \t %9d \t %5.2f%%\n", 10 * i, qv_decade_counter[i], 
                   100 * (double) qv_decade_counter[i] / (double) base_count);

        fprintf(fout, "\n");

    FREE(info.previous_cube);
    FREE(info.current_cube);
    FREE(previous_highest_cut_parameter_index);
    FREE(qv_counter);
    FREE(qv_decade_counter);
    return SUCCESS;
error:
    FREE(info.previous_cube);
    FREE(info.current_cube);
    FREE(previous_highest_cut_parameter_index);
    FREE(qv_counter);
    FREE(qv_decade_counter);
    return ERROR;
}

/***************************************************************************
 * display_cubes
 *
 * Purpose: display the (n-1) dimensional cube from the n-dimensional space.
 * 
 ****************************************************************************
 */
void
display_cubes(int l, INFO *info)
{
    int i, j, k;
    CUT *cut;

    for (k = 0; k < info->parameter[2].threshold_count; k++)
    {
        for (j = 0; j < info->parameter[3].threshold_count; j++)
            for (i = 0; i < info->parameter[3].threshold_count; i++)
            {
                cut = get_cut_4(i, j, k, l - 1, PREVIOUS, info);
                printf("prev(%d,%d,%d,%d)=%lu,%lu ", i, j, k, l - l,
                       cut->correct, cut->incorrect);
            }

        printf("\n");
    }

    for (k = 0; k < info->parameter[2].threshold_count; k++)
    {
        for (j = 0; j < info->parameter[3].threshold_count; j++)
            for (i = 0; i < info->parameter[3].threshold_count; i++)
            {
                cut = get_cut_4(i, j, k, l, CURRENT, info);
                printf("curr(%d,%d,%d,%d)=%lu,%lu ", i, j, k, l,
                       cut->correct, cut->incorrect);
            }

        printf("\n");
    }

    printf("\n");
}

#if DISPLAY_BASES
/***************************************************************************
 * display_number_of_correct_bases_in_each_bin
 *
 * called by: main
 * calls: get_bin
 *
 * ASSUMES: PARAMETER_COUNT = 4
 * 
 ***************************************************************************/
void
display_number_of_correct_bases_in_each_bin(PARAMETER *parameter, BIN *bin)
{
    int i, j, k, l, m, n;
    int *index;
    BIN *current_bin;

    index = (int *) malloc(sizeof(int) * PARAMETER_COUNT);
    if (index == NULL)
    {
        fprintf(stderr, "couldn't malloc index in display_number_of_correct_bases_in_each_bin\n");
        exit(-1);
    }

    if (PARAMETER_COUNT == 4)
    for (i = 0; i < parameter[0].threshold_count; i++)
        for (j = 0; j < parameter[1].threshold_count; j++)
            for (k = 0; k < parameter[2].threshold_count; k++)
                for (l = 0; l < parameter[3].threshold_count; l++)
                {
                    index[0] = i;
                    index[1] = j;
                    index[2] = k;
                    index[3] = l;

                    current_bin = get_bin(bin, parameter, index);

                    printf("bin(%d,%d,%d,%d)=(%f,%f,%f,%f) correct=%lu incorrect=%lu\n",
                            i, j, k, l,
                            parameter[0].threshold[i],
                            parameter[1].threshold[j],
                            parameter[2].threshold[k],
                            parameter[3].threshold[l],
                            current_bin->correct, current_bin->incorrect);
                }
    else /* PARAMETER_COUNT == 6 */
    for (i = 0; i < parameter[0].threshold_count; i++)
        for (j = 0; j < parameter[1].threshold_count; j++)
            for (k = 0; k < parameter[2].threshold_count; k++)
                for (l = 0; l < parameter[3].threshold_count; l++)
                    for (m = 0; m < parameter[4].threshold_count; m++)
                        for (n = 0; n < parameter[5].threshold_count; n++)
                        {
                            index[0] = i;
                            index[1] = j;
                            index[2] = k;
                            index[3] = l;
                            index[4] = m;
                            index[5] = n;

                            current_bin = get_bin(bin, parameter, index);

                            printf("bin(%d,%d,%d,%d,%d,%d)=(%f,%f,%f,%f,%f,%f) correct=%lu incorrect=%lu\n",
                                    i, j, k, l, m, n,
                                    parameter[0].threshold[i],
                                    parameter[1].threshold[j],
                                    parameter[2].threshold[k],
                                    parameter[3].threshold[l],
                                    parameter[4].threshold[m],
                                    parameter[5].threshold[n],
                                    current_bin->correct, current_bin->incorrect);
                        }

    FREE(index);
}
#endif

/**************************************************************************
 * Function: compress_params_array          
 * Purpose:  compress all duplicated parameter values into unique values 
 *           and set its weight to the original number of the 
 *           duplicates
 **************************************************************************
 */
int
compress_params_array(double **A, int **w, int *len)
{
    int     i, j;
    int     max_weight=0, *wn;
    double  A_min=(*A)[0];

#if 0
    for (i=0; i<*len; i++) {
        if (max_weight<(*w)[i]) {
            max_weight = (*w)[i];
        }
    }
    fprintf(stderr, "Begin compress_params_array: max_weight=%d\n",max_weight);
#endif
    wn = (int*)calloc(*len, sizeof(int));

    j=0;
    max_weight=0;
    for (i=0; i<*len; i++) {
        if ((*w)[i] <=0) {
            fprintf(stderr, "ERROR: weight[%d]=%d<0!!!\n", i, (*w)[i]);
        }
        wn[j] += (*w)[i];
#if 0  
        if (max_weight < wn[j]) {
            max_weight = wn[j];
        }
#endif
        /* Compress into one all the values that are distanced 
         * one from anotherby no more than MIN_DIST
         */
        if (i < *len-1 && (*A)[i] < (*A)[i+1]-MIN_DIST) {
            (*A)[j]  = (A_min + (*A)[i])/2.;
            j++;
            A_min = (*A)[i+1];
        }
        else if (i == *len-1) {
            (*A)[j]  = (A_min + (*A)[i])/2.;
            j++;
        }
    }        
#if 0
    fprintf(stderr, "In compress_params_array: max_weight=%d\n", 
        max_weight);
#endif
    /* Give back unused memory */                     
   *len = j;         

    max_weight=0;
    for (i=0; i<*len; i++) {
        (*w)[i]=wn[i];
#if 0
        if (max_weight<(*w)[i]) {
            max_weight = (*w)[i];
        }
#endif
    }
#if 0
    fprintf(stderr, "End of compress_params_array: max_weight=%d\n",
        max_weight);
#endif
    FREE(wn);
    wn=NULL;

    return SUCCESS;
}


/***************************************************************************
 * Function: read_and_sort_trace_parameters
 * Purpose: read 4 arrays of trace parameters from stdin and sort each the 
 *          array. If identical parameter values are found after the 
 *          sorting, compress them them into one unique value and increment
 *          the corresponding "weight" number for this parameter value 
 * Returns: 4 parameter arrays, each having (generally) different length,
 *          4 arrays of weights for each the parameter and the integer array
 *          of length 4 which indicates the length of each the parameter
 *          array
 *
 * called by: main
 * calls: sort_parameters
 *
 ***************************************************************************/
int   
read_and_sort_trace_parameters(char *InputName,
    unsigned long *base_count, PARAMETER parameter[])
{
    /* ind[j] is the current value for the index of j-th parameter
     * array in which the value from the train file will be read in
     */
    int      i;
    int      len[6]={0,0,0,0,0,0}, j;    
    int      linenum = 0;
    char     train_name[BUFLEN], *s, *current;
    double   is_match;
    int      spos, cpos;   /* Sample and consensus positions. */
    char     schar, cchar; /* Sample and consensus characters. */
    char     buffer[BUFLEN];
    FILE    *fileoffiles=NULL, *trainfile=NULL;

    if ((fileoffiles=fopen(InputName,"r"))== NULL) {
        fprintf(stderr, "Unable to open file_of_files '%s'\n",
            InputName);
        return ERROR;
    }

    for (j=0; j<PARAMETER_COUNT; j++) {
        parameter[j].num_val = BASE_COUNT_SCALE;
        parameter[j].value = (double*)malloc(BASE_COUNT_SCALE *
            sizeof(double));
        parameter[j].weight = (int*)calloc(BASE_COUNT_SCALE,
            sizeof(int));
        for (i=0; i<BASE_COUNT_SCALE; i++) {
            parameter[j].weight[i]=1;
        }
    }

   *base_count = 0;
    current = buffer;

    /* Scan the list of all train files */
    while (fgets(current, BUFLEN, fileoffiles) != NULL)
    {
        /* Ignore all white space lines and comments */
        if (strspn(current, " \t\r\n") == strlen(current)) {
            continue;
        }
        if ((current[0] == '#')
            || ((current[0] == '/') && (current[1] == '*'))
            || (current[0] == ';'))
        {
            continue;
        }

        /* Ignore white space in the train file name */
        sscanf(current, "%s", train_name);

        /* Open an individual train file */
        if ((trainfile=fopen(train_name,"r"))== NULL) {
            fprintf(stderr, "Unable to open train file '%s'\n",
                InputName);
            return ERROR;
        }

        /* Read bases from the individual train file */
        linenum = 0;
        while (fgets(train_name, sizeof(train_name), trainfile) != NULL)
        {
            linenum++;

/* Ignore all white space lines and comments */
            if (strspn(train_name, " \t\r\n") == strlen(train_name))
                continue;

            if (train_name[0] == '#' || train_name[0] == ';' ||
               (train_name[0] == '/' && train_name[1] == '*'))
                continue;

/* Get rid of 1) sample position, 2) sample base,
 * 3) consensus position and 4) consensus base. */

            s = strtok(train_name, " \t\n");
            if (s == NULL) continue;
            spos = atoi(s);

            s = strtok(NULL, " \t\n");
            if (s == NULL) continue;
            schar = s[0];

            s = strtok(NULL, " \t\n");
            if (s == NULL) continue;
            cpos = atoi(s);

            s = strtok(NULL, " \t\n");
            if (s == NULL) continue;
            cchar = s[0];

            s += strlen(s) + 1;
    
/* If sample == '-', there are no training parameters. */

            if (schar == '-') continue;

            if (PARAMETER_COUNT == 4)
            {
                if (Btk_atod(&s, &is_match) != SUCCESS
                ||  Btk_atod(&s, &parameter[0].value[len[0]]) != SUCCESS
                ||  Btk_atod(&s, &parameter[1].value[len[1]]) != SUCCESS
                ||  Btk_atod(&s, &parameter[2].value[len[2]]) != SUCCESS
                ||  Btk_atod(&s, &parameter[3].value[len[3]]) != SUCCESS)
                {
                    fprintf(stderr, "line %d:\n%s\nmissing/garbled parameter; skipping\n",
                    linenum, current);
                    continue;
                }
           }
           else if (PARAMETER_COUNT == 6)
           {
               if (Btk_atod(&s, &is_match) != SUCCESS
                ||  Btk_atod(&s, &parameter[0].value[len[0]]) != SUCCESS
                ||  Btk_atod(&s, &parameter[1].value[len[1]]) != SUCCESS
                ||  Btk_atod(&s, &parameter[2].value[len[2]]) != SUCCESS
                ||  Btk_atod(&s, &parameter[3].value[len[3]]) != SUCCESS
                ||  Btk_atod(&s, &parameter[4].value[len[4]]) != SUCCESS
                ||  Btk_atod(&s, &parameter[5].value[len[5]]) != SUCCESS)
                {
                    fprintf(stderr, "line %d:\n%s\nmissing/garbled parameter; skipping\n",
                    linenum, current);
                    continue;
                }
           } 
           (*base_count)++; 
            if ((*base_count > 1 && *base_count % BASE_COUNT_SCALE) == 0) {
                fprintf(stderr, "\r   %lu parameters have been read       ", 
               *base_count);
            }

            for (j=0; j<PARAMETER_COUNT; j++) {
                parameter[j].weight[len[j]] = 1;
    
                /* Currently allocated room for j-th parameter exhausted;
                 * allocate more room
                 */
                if (len[j] == parameter[j].num_val-1) {
                    parameter[j].value = (double*)realloc(parameter[j].value, 
                        (parameter[j].num_val+BASE_COUNT_SCALE)*sizeof(double));
                    parameter[j].weight = (int*)realloc(parameter[j].weight, 
                        (parameter[j].num_val+BASE_COUNT_SCALE)*sizeof(int));
                    parameter[j].num_val += BASE_COUNT_SCALE;
                }
            }

            /* data portion data of size BASE_COUNT_SCALE has been read;
             * Re-sort all the parameters that have been already read;
             * then compress the parameter array
             */
            if (*base_count > 1 &&
                *base_count % BASE_COUNT_SCALE == 0)
            {
//              fprintf(stderr, "\nBase count=%lu\n", *base_count);
                for (j=0; j<PARAMETER_COUNT; j++) {
#if 0
                    fprintf(stderr, "\nParameter[%d]:\n", j);
                    max_weight=0;
                    min_weight=0;
                    for (i=0; i<len[j]+1; i++) {
                        if (max_weight<parameter[j].weight[i])
                            max_weight = parameter[j].weight[i];
                        if (min_weight>parameter[j].weight[i])
                            min_weight = parameter[j].weight[i];
                    }
                    fprintf(stderr, "    before quicksort2 num_val=%d, max_weight=%d, min_weight=%d\n", 
                        len[j], max_weight, min_weight);
#endif
                    
                    quicksort2(parameter[j].value, parameter[j].weight, 0, len[j]-1); 
#if 0
                    max_weight=0;
                    min_weight=0;
                    for (i=0; i<len[j]+1; i++) {
                        if (max_weight<parameter[j].weight[i])
                            max_weight = parameter[j].weight[i];
                        if (min_weight>parameter[j].weight[i])
                            min_weight = parameter[j].weight[i];
                    }
                    fprintf(stderr, "    after  quicksort2 num_val=%d, max_weight=%d, min_weight=%d\n", 
                        len[j], max_weight, min_weight);
#endif
                    compress_params_array(&parameter[j].value,
                                          &parameter[j].weight,
                                          &len[j]);    
#if 0
                    max_weight=0;
                    min_weight=0;
                    for (i=0; i<len[j]+1; i++) {
                        if (max_weight<parameter[j].weight[i])
                            max_weight = parameter[j].weight[i];
                        if (min_weight>parameter[j].weight[i])
                            min_weight = parameter[j].weight[i];
                    }          
                    fprintf(stderr, "    after  compress   num_val=%d, max_weight=%d, min_weight=%d\n", 
                        len[j], max_weight, min_weight);
#endif
                    parameter[j].num_val = len[j]+1 + BASE_COUNT_SCALE;
                    parameter[j].value = (double*)realloc(parameter[j].value,
                        sizeof(double)*(len[j] +1+ BASE_COUNT_SCALE));
                    parameter[j].weight = (int*)realloc(parameter[j].weight,
                        sizeof(int)*(len[j] +1+ BASE_COUNT_SCALE));
                }
            } 
            for (j=0; j<PARAMETER_COUNT; j++) {
                len[j]++;
            }
        }     /* loop in lines of a particular trainfile */
        fclose(trainfile);
    }         /* loop in lines of fileoffiles */
    fclose(fileoffiles);
    fprintf(stderr, "\n");

/* Give back any unused space */
#if 0
    fprintf(stderr, "Base count=%lu\n", *base_count);
#endif
    for (j=0; j<PARAMETER_COUNT; j++) {
#if 0
        fprintf(stderr, "\nParameter[%d]:\n", j);
#endif
#if 0
        max_weight=0;
        min_weight=0;
        for (i=0; i<len[j]; i++) {
            if (max_weight<parameter[j].weight[i])
                max_weight = parameter[j].weight[i];
            if (min_weight>parameter[j].weight[i])
                min_weight = parameter[j].weight[i];
        }
        fprintf(stderr, "    before last quicksort2 num_val=%d, max_weight=%d, min_weight=%d\n",
            len[j], max_weight, min_weight);
#endif
#if 0
        fprintf(stderr, "\nFirst 20 values and weights:\n");
        for (i=0; i<20; i++) {
            fprintf(stderr, "%6d  %3.15f %d\n", i, parameter[j].value[i],
            parameter[j].weight[i]);
        }
#endif

        parameter[j].num_val = len[j];
        quicksort2(parameter[j].value, parameter[j].weight, 0, len[j]-1);
#if 0
        fprintf(stderr, "\nFirst 50 values and weights:\n");
        for (i=0; i<50; i++) {
            fprintf(stderr, "%6d  %3.15f %d\n", i, parameter[j].value[i],
            parameter[j].weight[i]);
        }
#endif
#if 0
        max_weight=0;
        min_weight=0;
        for (i=0; i<len[j]; i++) {
            if (max_weight<parameter[j].weight[i])
                max_weight = parameter[j].weight[i];
            if (min_weight>parameter[j].weight[i])
                min_weight = parameter[j].weight[i];
        }
        fprintf(stderr, "    after  last quicksort2 num_val=%d, max_weight=%d, min_weight=%d\n",
            len[j], max_weight, min_weight);
#endif
        compress_params_array(&parameter[j].value, &parameter[j].weight,
                              &len[j]);
#if 0
        fprintf(stderr, "\nFirst 50 values and weights:\n");
        for (i=0; i<50; i++) {
            fprintf(stderr, "%6d  %3.15f %d\n", i, parameter[j].value[i],
            parameter[j].weight[i]);
        }
#endif
#if 0
        max_weight=0;
        min_weight=0;
        for (i=0; i<len[j]; i++) {
            if (max_weight<parameter[j].weight[i])
                max_weight = parameter[j].weight[i];
            if (min_weight>parameter[j].weight[i])
                min_weight = parameter[j].weight[i];
        }
        fprintf(stderr, "    after  last compress   num_val=%d, max_weight=%d, min_weight=%d\n",
            len[j], max_weight, min_weight);
#endif
        parameter[j].num_val = len[j];
#if 0 
        /* Do not realloc here!!! For some reason, it does not work !!! */
        parameter[j].value = realloc(parameter[j].value, len[j]);
        parameter[j].weight = realloc(parameter[j].weight, len[j]);
#endif  
#if 0
        fprintf(stderr, "\nFirst 50 values and weights:\n");
        for (i=0; i<50; i++) {
            fprintf(stderr, "%6d  %3.15f %d\n", i, parameter[j].value[i],
            parameter[j].weight[i]);
        }
#endif
    }

    return SUCCESS;
}

int
main(int argc, char *argv[])
{
    int           i, num_bins = 1;
    BASE         *base=NULL;
    BIN          *bin;
    PARAMETER    *parameter;
    FILE         *fout;
    unsigned int  threshold_count[6];  /* Number of thresholds used */
    unsigned long initial_base_room, base_count;
    time_t        t1, t2;
    int           development = 0;    
    BtkMessage    msg;
    BtkMessage   *message = &msg;

    /* FileOfFiles = 0 if input from stdin; =1 if from inp_file_name */
    int           FileOfFiles = 0;
    char         *InputName = "\0";
    char          buffer[BUFLEN];

    t1 = time((time_t) NULL);
    InputName     = buffer;
    InputName[0]  = '\0';
    OutputName[0] = '\0';
    OutputSpecified     = 0;

/* Set defaults */

    Verbose = 1;
    Compress = 1;

    initial_base_room = BASE_COUNT_SCALE;

    opterr = 0;
    while ((i = getopt(argc, argv, "b:f:o:cCdQV")) != EOF)
        switch (i)
        {
            case 'b':
                if (sscanf(optarg, "%lu", &initial_base_room) != 1
                ||  initial_base_room < 1)
                {
                    show_usage(argc, argv);
                    exit(2);
                }
                break;
            case 'c':
                Compress++;
                break;
            case 'C':
                Compress = 0;
                break;
            case 'd':
                development++;
                break;
            case 'f':
                FileOfFiles++;
                InputName = (char*)malloc(BUFLEN*sizeof(char));
                if (sscanf(optarg, "%s", InputName) != 1) {
                    fprintf(stderr, "Case f\n");
                    FREE(InputName);
                    show_usage(argc, argv);
                    exit(2);
                }
                break;
            case 'o':
                OutputSpecified++;
                (void)strncpy(OutputName, optarg, sizeof(OutputName));
                break;
            case 'Q':
                Verbose = 0;
                break;
            case 'V':
                Verbose++;
                break;
            case '?':
            default:
                show_usage(argc, argv);
                exit(2);
        }

//  fprintf(stderr, "optind= %d argc= %d\n", optind, argc);
    if ((argc - optind != 4 &&  argc - optind != 6))
    {
        show_usage_dev(argc, argv);
        fprintf(stderr,
            "The number of specified threshold counts should be 4 or 6\n");
        exit(2);
    }    

    while (optind + PARAMETER_COUNT < argc)
    {
        if (sscanf(argv[optind+PARAMETER_COUNT], "%u", 
            &threshold_count[PARAMETER_COUNT]) != 1)
        {
            fprintf(stderr, 
                "Error reading threshold count for %d-th parameter\n",  
                PARAMETER_COUNT+1);
            exit(-2);
        }
        num_bins *= threshold_count[PARAMETER_COUNT];
//      fprintf(stderr, "parameter= %d thresholds= %d\n", 
//          PARAMETER_COUNT, threshold_count[PARAMETER_COUNT]);
        PARAMETER_COUNT++;
    }

    fprintf(stderr, "%d bins total\n", num_bins);

    parameter = CALLOC(PARAMETER, PARAMETER_COUNT);

    if (OutputSpecified) {
        /* Set up output file. */
        if ((fout=fopen(OutputName, "w"))== NULL) {
            fprintf(stderr, "Cannot open output file '%s'\n", OutputName);
            exit(ERROR);
        }
        fprintf(stderr, "Software Version: %s\n", TT_VERSION);
    } else {
        if (Verbose > 1)
            fprintf(stderr, "No output file specified.  Using stdout.\n");
        fout=stdout;
    }


    setbuf(stderr, NULL);    /* so that any status comes out w/o delay */

    if (Verbose > 1)
        fprintf(stderr, "starting with room for %lu bases\n", initial_base_room);

    /* 1. Read bases from train file
     ****************************
     */
    if (FileOfFiles && Verbose) {
        fprintf(stderr, "\nReading trace parameters ... \n");
        read_and_sort_trace_parameters(InputName, &base_count, parameter);
    }
    else {
        if (Verbose)
            fprintf(stderr, "Reading bases ... \n");
        base = get_bases(initial_base_room, &base_count, PARAMETER_COUNT, 0);
    }
    t2 = time((time_t) NULL);

    if (FileOfFiles && Verbose)
        fprintf(stderr, "%lu parameters read in %f sec\n", base_count, 
            (double)(t2 - t1));
    else 
        if (Verbose)
            fprintf(stderr, "\n%lu bases read in %f sec\n", base_count, 
            (double)(t2 - t1));

    if (base_count == 0)
    {
        fputs("no valid bases\n", stderr);
        exit(1);
    }

    t1 = t2;

    for (i = 0; i < PARAMETER_COUNT; i++)
    {
        parameter[i].threshold_count = threshold_count[i];
        parameter[i].threshold = CALLOC(double, parameter[i].threshold_count);
        MEM_ERROR(parameter[i].threshold);
    }

#if 0
    j=2;
    fprintf(stderr, "PARAMETER%d %d values and weights:\n",
        j, parameter[j].num_val);
    for (i=0; i<parameter[j].num_val; i++) {
        fprintf(stderr, "%6d  %3.15f %d\n", i, parameter[j].value[i],
        parameter[j].weight[i]);
    }
#endif
    /* 2. Compute the trace parameter thresholds
     ****************************************
     */
    fprintf(stderr, "\nComputing trace parameter thresholds ... \n");
    if (FileOfFiles) {
        get_thresholds2(parameter, PARAMETER_COUNT);
    }
    else {
        get_thresholds(base, base_count, parameter, PARAMETER_COUNT);
    }

    for (parameter[0].dimension = i = 1; i < PARAMETER_COUNT; i++)
         parameter[i].dimension = parameter[i - 1].dimension *
                                  parameter[i - 1].threshold_count;

    t2 = time((time_t) NULL);

    if (Verbose)
        fprintf(stderr, "thresholds determined in %f sec\n", (double)(t2 - t1));

    t1 = t2;

    if (PARAMETER_COUNT == 4)
        bin = (BIN *) calloc(parameter[0].threshold_count
                           * parameter[1].threshold_count
                           * parameter[2].threshold_count
                           * parameter[3].threshold_count
                           , sizeof(BIN));
    else /* PARAMETER_COUNT == 6 */
        bin = (BIN *) calloc(parameter[0].threshold_count
                           * parameter[1].threshold_count
                           * parameter[2].threshold_count
                           * parameter[3].threshold_count
                           * parameter[4].threshold_count
                           * parameter[5].threshold_count
                           , sizeof(BIN));


    MEM_ERROR(bin);

    /* 3. Populate bins
     *******************
     */
    fprintf(stderr, "\nPopulating bins ... \n");
    if (FileOfFiles) {
        read_bases_and_populate_bins(InputName, &base_count,
            parameter, bin);
    }
    else {
        count_number_of_correct_bases_in_each_bin(base, base_count,
            parameter, bin, fout);
    }
    t2 = time((time_t) NULL);

    if (Verbose)
        fprintf(stderr, "bins populated in %f sec\n", (double)(t2 - t1));

#if DISPLAY_BASES
    display_number_of_correct_bases_in_each_bin(parameter, bin);
#endif

    /* 4. Generate a lookup table 
     ****************************
     */
    t1 = t2;
    fprintf(fout, "#\n# Version %s\n#\n", TT_VERSION);

    fprintf(fout, "\n#  Numbers of parameter thresholds specified: ");
    for (i=0; i<PARAMETER_COUNT; i++)
        fprintf(fout, " %d", threshold_count[i]);
    fprintf(fout, " \n");
    fprintf(fout, 
        "\n#  Numbers of parameter thresholds actually established: ");
    for (i=0; i<PARAMETER_COUNT; i++) 
        fprintf(fout, " %d", parameter[i].threshold_count);
    fprintf(fout, " \n");

    fprintf(fout, "\n#  Parameter thresholds:\n");

    for (i=0; i<MAX_NUM_THRESHOLDS; i++) {
        int j;
        if (PARAMETER_COUNT == 4 &&
            (i >= parameter[0].threshold_count) && 
            (i >= parameter[1].threshold_count) &&
            (i >= parameter[2].threshold_count) &&
            (i >= parameter[3].threshold_count))   
            break;   

        else if (PARAMETER_COUNT == 6 &&
            (i >= parameter[0].threshold_count) &&
            (i >= parameter[1].threshold_count) &&
            (i >= parameter[2].threshold_count) &&
            (i >= parameter[3].threshold_count) &&
            (i >= parameter[4].threshold_count) &&
            (i >= parameter[5].threshold_count))
            break;

        for (j=0; j<PARAMETER_COUNT; j++) {
            if (i < parameter[j].threshold_count)
                fprintf(fout, " %13.6f", parameter[j].threshold[i]);
            else {
                fprintf(fout, " %13.6f", 
                    parameter[j].threshold[parameter[j].threshold_count-1]);
            }
        }
        fprintf(fout, "\n");
    }    
    fprintf(stderr, "\nGenerating a lookup table ... \n");
    if (create_qv_table_via_dynamic_programming(bin, parameter, PARAMETER_COUNT,
        base_count, fout, &msg) != SUCCESS)
        goto error;
    t2 = time((time_t) NULL);

    if (Verbose)
        fprintf(stderr, "lookup table generated in %f sec\n", (double)(t2 - t1));

    FREE(bin);
    FREE(base);
    for (i = 0; i < PARAMETER_COUNT; i++) {
        FREE(parameter[i].threshold);
    }
    FREE(parameter);
    if (FileOfFiles) {
        FREE(InputName);
    }
    return SUCCESS;
error:
    FREE(bin);
    FREE(base);
    for (i = 0; i < PARAMETER_COUNT; i++) {
        FREE(parameter[i].threshold);
    }
    FREE(parameter);
    if (FileOfFiles) {
        FREE(InputName);
    }
    return ERROR;
}

