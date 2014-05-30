/**************************************************************************
 * This file is part of TraceTuner, the DNA sequencing quality value,
 * base calling and trace processing software.
 *
 * Copyright (c) 2003-2008 Gennady Denisov.  All rights reserved.
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

/***************************************************************************
 * checklutparam.c       1.12
 *
 * Purpose: 
 * - create lookup table for quality values using 4 parameters.
 * - count the number of bases having each QV 0 to 100 and display percentages
 *   over QV 20, 30, 40, etc.
 *
 * Calls: 
 * - get_bases
 * - count_number_of_correct_bases_in_each_bin
 * - create_qv_table_via_dynamic_programming
 *
 **************************************************************************
 */
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

#include "func_name.h"
#include "select.h"
#include "params.h"
#include "lut.h"
#include "Btk_atod.h"
#include "get_thresholds.h"
#include "Btk_lookup_table.h"
#include "check_data.h"
#include "Btk_qv.h"

#define BASE_COUNT_SCALE	(100000)
#define MAXLINE			(1000)


int Verbose ;		/* How much status info to print, if any */
int Compress;	/* Whether to compress thresholds */


static void 
show_usage(int argc, char *argv[])
{
    (void)fprintf(stderr, "\nVersion: %s\n", TT_VERSION);
    (void)fprintf(stderr, 
		"usage: %s [-cCqv] [-b <initialbaseroom>] <num_thresholds>\n",
		argv[0]);
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
  int i = threshold_count-1;
  /* thresholds are indexed 0 to threshold_count-1 */

  if (value > threshold[i]) {
    printf("problem: a value is greater than the max threshold\n");
    exit(-1);
  }
  while ((i >= 0) && (value <= threshold[i])) {
    i--;
  }
  return (i+1);
} /* get_threshold_index */

/***************************************************************************
 * get_bin
 * <index> is an array of <parameter_count> indices
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

  for (i = 0; i < PARAMETER_COUNT; i++) {
    bin_number = index[i] * parameter[i].dimension + bin_number;
  }
#if 0
  printf("bin_number(%d,%d,%d,%d)=%lu\n",index[0],index[1],index[2],index[3],
         bin_number);
  printf("&(bin[bin_number])=%lu\n", (unsigned long)&(bin[bin_number]));
#endif
  return &(bin[bin_number]);
} /* get_bin */

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
 * ASSUMES: 4 parameters
 *
 * note that <l> is just used to check if we have a boundary condition.
 * otherwise <i,j,k> define the cube in the 4D space that we are
 * currently working on.  to save memory, we don't work with the whole
 * hypercubic (4d) space, just 2 cubic (3d) portions (a CURRENT and
 * a PREVIOUS) of it at a time.
 **************************************************************************
 */
CUT *
get_cut(int i, int j, int k, int l, TIME time, INFO *info)
{
  unsigned long cut_number;

  if (i == -1 || j == -1 || k == -1 || l == -1)
    return &(info->boundary_cut);

  cut_number  = i * info->parameter[0].dimension
              + j * info->parameter[1].dimension
              + k * info->parameter[2].dimension;

  if (time == CURRENT) {
    return &(info->current_cube[cut_number]);
  }
  else if (time == PREVIOUS) {
    return &(info->previous_cube[cut_number]);
  }
  else {
    printf("invalid time in function get_cut\n");
    exit(-1);
  }
} /* get_cut */

/***************************************************************************
 * initialize_highest_qv_cut
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
    highest_qv_cut->index[m] = -1; /* the minimum real index is 0 */
  highest_qv_cut->correct_base_call_count = 0;
  highest_qv_cut->incorrect_base_call_count = 0;
  highest_qv_cut->total_base_call_count = 0;
  highest_qv_cut->error_rate = 0;
  highest_qv_cut->quality_value = 0;
} /* initialize_highest_qv_cut */

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
 * ASSUMES: 4 parameters
 *
 ***************************************************************************/
CUT *
number_in_cut(int i, int j, int k, int l, INFO *info,
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

  /*
   * c_1111 always gives us a non-boundary cube, as all of the subscripts
   * are >= 0.
   */
  c_1111 = get_cut(i  ,j  ,k  ,l  ,CURRENT ,info);
  c_1101 = c_1111 - info->dimension2;
  c_1001 = c_1101 - info->dimension1;
  c_1011 = c_1001 + info->dimension2;
  c_0011 = c_1011 - info->dimension0;
  c_0111 = c_0011 + info->dimension1;
  c_0101 = c_0111 - info->dimension2;
  c_0001 = c_0101 - info->dimension1;
  /*
   * Use the offset of c_0001 from CURRENT as the offset of c_0000 from
   * PREVIOUS.  We can't easily start over from a get_cut() call, since we
   * can't be sure one of the subscripts isn't 0 at this point, making one
   * of the -1 values < 0 (boundary condition).
   */
  c_0000 = info->previous_cube + (c_0001 - info->current_cube);
  c_0010 = c_0000 + info->dimension2;
  c_0110 = c_0010 + info->dimension1;
  c_0100 = c_0110 - info->dimension2;
  c_1100 = c_0100 + info->dimension0;
  c_1110 = c_1100 + info->dimension2;
  c_1010 = c_1110 - info->dimension1;
  c_1000 = c_1010 - info->dimension2;

  if (i == 0) {
    c_0000 = &info->boundary_cut;
    c_0001 = &info->boundary_cut;
    c_0010 = &info->boundary_cut;
    c_0011 = &info->boundary_cut;
    c_0100 = &info->boundary_cut;
    c_0101 = &info->boundary_cut;
    c_0110 = &info->boundary_cut;
    c_0111 = &info->boundary_cut;
  }

  if (j == 0) {
    c_0000 = &info->boundary_cut;
    c_0001 = &info->boundary_cut;
    c_0010 = &info->boundary_cut;
    c_0011 = &info->boundary_cut;
    c_1000 = &info->boundary_cut;
    c_1001 = &info->boundary_cut;
    c_1010 = &info->boundary_cut;
    c_1011 = &info->boundary_cut;
  }

  if (k == 0) {
    c_0000 = &info->boundary_cut;
    c_0001 = &info->boundary_cut;
    c_0100 = &info->boundary_cut;
    c_0101 = &info->boundary_cut;
    c_1000 = &info->boundary_cut;
    c_1001 = &info->boundary_cut;
    c_1100 = &info->boundary_cut;
    c_1101 = &info->boundary_cut;
  }

  if (l == 0) {
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
      l <= previous_highest_cut_parameter_index[3]) {
    b_1111->correct = 0;
    b_1111->incorrect = 0;
  }

  c_1111->correct =
    /* the layout in space of the cuts are as follows--
       each row of comments corresponds to a row of equations */
    /* same cube (3d-space), same plane (2d-space):
       same bin, cut to left, cut to right, cut diagonal */
    /* same cube (3d-space), previous plane (2d-space):
       cut below, bottom left, bottom right, bottom diagonal */
    +b_1111->correct+c_0111->correct+c_1011->correct-c_0011->correct
    +c_1101->correct-c_0101->correct-c_1001->correct+c_0001->correct
    +c_1110->correct-c_0110->correct-c_1010->correct+c_0010->correct
    -c_1100->correct+c_0100->correct+c_1000->correct-c_0000->correct;

  c_1111->incorrect =
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
    +b_1111->incorrect+c_0111->incorrect+c_1011->incorrect-c_0011->incorrect
    +c_1101->incorrect-c_0101->incorrect-c_1001->incorrect+c_0001->incorrect
    +c_1110->incorrect-c_0110->incorrect-c_1010->incorrect+c_0010->incorrect
    -c_1100->incorrect+c_0100->incorrect+c_1000->incorrect-c_0000->incorrect;

  return (c_1111);

} /* number_in_cut */

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
					   unsigned long base_count,
					   PARAMETER *parameter,
                                           BIN *bin)
{
  unsigned long n;
  int *index;
  unsigned long bin_count;
  BASE *base;
  int i;
  BIN *current_bin;

  bin_count = parameter[PARAMETER_COUNT-1].threshold_count
	      * parameter[PARAMETER_COUNT-1].dimension;

  index = (int *) malloc(sizeof(int)*PARAMETER_COUNT);
  if (index == NULL) {
    printf("couln't malloc index in count_number_of_correct_bases_in_bin\n");
    exit(-1);
  }
  for (n = 0, current_bin = bin; n < bin_count; n++, current_bin++) {
    current_bin->correct = 0;
    current_bin->incorrect = 0;
  }

  for (n = 0, base = base_array; n < base_count; n++, base++) {
    for (i = 0; i < PARAMETER_COUNT; i++) {
      index[i] = get_threshold_index(base->parameter[i], 
				     parameter[i].threshold, 
				     parameter[i].threshold_count);
    }
    current_bin = get_bin(bin, parameter, index);
    if (base->is_match) {
      current_bin->correct++;
    }
    else {
      current_bin->incorrect++;
    }
  }

  free(index);
} /* count_number_of_correct_bases_in_each_bin */
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
 * ASSUMPTION: 4 parameters
 *
 * called by: create_qv_table_via_dynamic_programming
 * calls: none
 *
 ***************************************************************************/
void
update_highest_qv_cut(HIGHEST_QV_CUT *highest_qv_cut, PARAMETER *parameter,
                       unsigned long correct_base_call_count,
                       unsigned long incorrect_base_call_count, int i, int j,
                       int k, int l, int *previous_highest_cut_parameter_index)
{
  unsigned long total_base_call_count;
  double error_rate;
  int quality_value;
  int sum_of_indices;

  total_base_call_count
    = correct_base_call_count + incorrect_base_call_count;
  /* the error rate includes a penalty for small sample size
     by adding 1 in numerator and denominator.
     hence, 1 correct and 1 incorrect base is assigned an
     error rate of 2/3, which is equivalent to the
     error rate for 66 incorrect bases and 33 correct bases. */
  if (incorrect_base_call_count==0)
  {
     error_rate =
       (((double) (1 + incorrect_base_call_count)) /
        ((double) (1 + total_base_call_count)));
  }
  else
  {
     error_rate =
       (((double) (incorrect_base_call_count)) /
        ((double) (total_base_call_count)));
  }
  quality_value = (int) rint(-10*log10(error_rate));
  sum_of_indices = i+j+k+l;
  if ((quality_value > highest_qv_cut->quality_value)

      || ((quality_value == highest_qv_cut->quality_value)
          && (total_base_call_count >
              highest_qv_cut->total_base_call_count))

      || ((quality_value == highest_qv_cut->quality_value)
          && (total_base_call_count ==
              highest_qv_cut->total_base_call_count)
          && (sum_of_indices > highest_qv_cut->sum_of_indices)))
  {
    highest_qv_cut->sum_of_indices = sum_of_indices;
    highest_qv_cut->index[0] = i;
    highest_qv_cut->index[1] = j;
    highest_qv_cut->index[2] = k;
    highest_qv_cut->index[3] = l;
    highest_qv_cut->correct_base_call_count =
      correct_base_call_count;
    highest_qv_cut->incorrect_base_call_count =

      incorrect_base_call_count;
    highest_qv_cut->total_base_call_count =
      total_base_call_count;
    highest_qv_cut->error_rate = error_rate;
    highest_qv_cut->quality_value = quality_value;
    highest_qv_cut->parameter[0] = parameter[0].threshold[i];
    highest_qv_cut->parameter[1] = parameter[1].threshold[j];
    highest_qv_cut->parameter[2] = parameter[2].threshold[k];
    highest_qv_cut->parameter[3] = parameter[3].threshold[l];
  } /* if (quality value > highest_qv_cut->quality_value) */
} /* update_highest_qv_cut */

/***************************************************************************
 * write_to_qv_table
 *
 * purpose: write the latest <highest_qv_cut> in a "C" programming language
 * style if statement that checks if each of the <parameter_count> parameters
 * falls below the values specified in the if statement.
 * this table can be used in a "C" function that assigns quality values.
 *
 * ASSUMES: 4 parameters
 *
 * called by: create_qv_table_via_dynamic_programming
 * calls: none
 *
 ***************************************************************************/
void
write_to_qv_table(HIGHEST_QV_CUT *highest_qv_cut)
{
  int m;
  static int table_entry_number= 0;

  table_entry_number++;
  printf("\n/* ENTRY %d: INDEX=%d(%d,%d,%d,%d) CORR= %ld INC= %ld TOT= %ld E=%f QV= %d */ \n",
         table_entry_number,
         highest_qv_cut->sum_of_indices,
         highest_qv_cut->index[0],
         highest_qv_cut->index[1],
         highest_qv_cut->index[2],
         highest_qv_cut->index[3],
         highest_qv_cut->correct_base_call_count,
         highest_qv_cut->incorrect_base_call_count,
         highest_qv_cut->total_base_call_count,
         highest_qv_cut->error_rate,
         highest_qv_cut->quality_value);
  printf("if (");
// -------------------------------------------------------------
    for (m = 0; m < (PARAMETER_COUNT-1); m++) {
        if (m != 3) {
            printf("(base->parameter[%d] <= %19.15f) && \n",
                m, highest_qv_cut->parameter[m]);
        }
            else {
                printf("(base->parameter[%d] <= %19.15f) && \n",
                    m, highest_qv_cut->parameter[m]);
        }
    }
    printf("(base->parameter[%d] <= %19.15f))\n",
        (PARAMETER_COUNT-1),
        highest_qv_cut->parameter[PARAMETER_COUNT-1]);
  printf("return (quality_value = %d);\n", highest_qv_cut->quality_value);
// -------------------------------------------------------------
  printf("\n");
} /* write_to_qv_table */

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
 * ASSUMES: 4 parameters
 *
 ***************************************************************************/
/* a cut is a set of thresholds */
void
create_qv_table_via_dynamic_programming(BIN *bin, PARAMETER parameter[],
					 unsigned long base_count)
{
  unsigned long countdown = base_count;
  HIGHEST_QV_CUT highest_qv_cut;
  int i,j,k,l, num_entries;
  INFO info;
  CUT *cut;
  CUT *temp_cube; /* temp pointer for use when switching current_cube and
		     previous_cube */
  int m; /* parameter index counter */
  int *qv_counter;
  int *qv_decade_counter;
  int *previous_highest_cut_parameter_index;
			/* keep track of threshold indices from
			   previous highest_qv_cut added to lookup table
			   for the purpose of zeroing out the number
			   of correct and incorrect calls in bins less than
			   or equal to these indices */


  qv_counter = (int *) calloc(MAX_QV,sizeof(int));
  qv_decade_counter = (int *) calloc((MAX_QV/10),sizeof(int));
  num_entries = 0;
  printf("\n/* QUALITY VALUE LOOKUP TABLE */ \n");
  printf("/* table was calibrated using %lu base calls */\n", 
	 base_count);

  /* for each cut (set of threshold values),
     count the number of correct and incorrect base calls
     in the data set to determine the quality value
     associated with that cut. 
     decide which cut has the highest quality value 
     and enter that cut into a table.
     repeat the procedure with the remaining data
     (the data above the last set of thresholds). */
  info.parameter = parameter;
  info.boundary_cut.correct = 0;
  info.boundary_cut.incorrect = 0;
  info.bin = bin;
  info.previous_cube = 
    (CUT *) malloc(sizeof(CUT)
		   *parameter[0].threshold_count
		   *parameter[1].threshold_count
		   *parameter[2].threshold_count); 
  /* cube size is not dependent on parameter[2].threshold_count */
  if (info.previous_cube == NULL) {
    printf("couldn't malloc info.previous_cube\n");
    exit(-1);
  }
  info.current_cube = 
    (CUT *) malloc(sizeof(CUT)
		   *parameter[0].threshold_count
		   *parameter[1].threshold_count
		   *parameter[2].threshold_count); 
  /* cube size is not dependent on parameter[2].threshold_count */
  if (info.current_cube == NULL) {
    printf("couldn't malloc info.current_cube\n");
    exit(-1);
  }

  info.dimension0 = parameter[0].dimension;
  info.dimension1 = parameter[1].dimension;
  info.dimension2 = parameter[2].dimension;
  info.dimension3 = parameter[3].dimension;

  for (k = 0; k < parameter[2].threshold_count; k++) {
    for (j = 0; j < parameter[1].threshold_count; j++) {
      for (i = 0; i < parameter[0].threshold_count; i++) {
	cut = get_cut(i,j,k,0,PREVIOUS,&info);
	cut->correct = 0;
	cut->incorrect = 0;
	cut = get_cut(i,j,k,0,CURRENT,&info);
	cut->correct = 0;
	cut->incorrect = 0;
      }
    }
  }
  previous_highest_cut_parameter_index = 
    (int *) malloc(sizeof(int)*PARAMETER_COUNT);
  if (previous_highest_cut_parameter_index == NULL) {
    puts("couldn't malloc previous_highest_cut_parameter_index");
    exit(-1);
  }

  for (m = 0; m < PARAMETER_COUNT; m++)
    previous_highest_cut_parameter_index[m] = -1; 
  /* initialize below the real minimum value 0 */

  do {
    initialize_highest_qv_cut(&highest_qv_cut);
    for (l = 0; l < parameter[3].threshold_count; l++) {
      temp_cube = info.previous_cube;      
      info.previous_cube = info.current_cube;
      info.current_cube = temp_cube;

      for (k = 0; k < parameter[2].threshold_count; k++) {
         for (j = 0; j < parameter[1].threshold_count; j++) {
            for (i = 0; i < parameter[0].threshold_count; i++) {

	    cut = number_in_cut(i,j,k,l, &info,
				previous_highest_cut_parameter_index);
	    update_highest_qv_cut(&highest_qv_cut, parameter, cut->correct,
				  cut->incorrect, i, j, k, l,
                                  previous_highest_cut_parameter_index);
	  } /* for parameter 3 */
	} /* for parameter 2 */
      } /* for parameter 1 */
    } /* for parameter 0 */
    if (highest_qv_cut.total_base_call_count != 0) {
      countdown -= highest_qv_cut.total_base_call_count;
      write_to_qv_table(&highest_qv_cut);
      for (m = 0; m < PARAMETER_COUNT; m++)
	previous_highest_cut_parameter_index[m] = highest_qv_cut.index[m];
      num_entries++;
      if (Verbose) {
	(void)fprintf(stderr,
			"\r%lu bases to go (%d entries so far)...        ",
			countdown, num_entries);
      }
    }
    /* keep track of how many bases have certain qv's */
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
  if (Verbose) {
    (void)fprintf(stderr, "\rdone.                         \n");
  }
  if (countdown != 0) {
    printf("/* warning: %lu base calls were unaccounted for */\n", countdown);
  }
  printf("\n/* NUMBER AND PERCENTAGE OF BASES HAVING EACH QUALITY VALUE\n");
  for (i = 0; i < MAX_QV; i++) {
    if (qv_counter[i] != 0) {
      printf("%2d \t %9d \t %5.2f%%\n", i, qv_counter[i], 
	     100*((double) qv_counter[i])/((double) base_count));
    }
  }
  printf("*/\n");
  printf("\n/* NUMBER AND PERCENTAGE OF BASES HAVING >= QUALITY VALUE\n");
  for (i = 1; i < (MAX_QV/10); i++) {
    if (qv_decade_counter[i] != 0) {
      printf("%2d \t %9d \t %5.2f%%\n", 10*i, qv_decade_counter[i], 
	     100*((double)qv_decade_counter[i])/((double)base_count));
    }
  }
  printf("*/\n");
  free(info.previous_cube);
  free(info.current_cube);
  free(previous_highest_cut_parameter_index);
  free(qv_counter);
  free(qv_decade_counter);
} /* create_qv_table_via_dynamic_programming */

/***************************************************************************
 * display_cubes
 *
 * Purpose: display the (n-1) dimensional cube from the n-dimensional space.
 * 
 * ASSUMES: 4 parameters.
 * 
 ****************************************************************************
 */
void
display_cubes(int l, INFO *info)
{
  int i, j, k;
  CUT *cut;

  for (k = 0; k < info->parameter[2].threshold_count; k++) {
    for (j = 0; j < info->parameter[3].threshold_count; j++) {
      for (i = 0; i < info->parameter[3].threshold_count; i++) {
	cut = get_cut(i,j,k,l-1,PREVIOUS,info);
	printf("prev(%d,%d,%d,%d)=%lu,%lu ", i,j,k,l-l,cut->correct, cut->incorrect);
      } /* for i */
    } /* for j */
    printf("\n");
  } /* for k */


  for (k = 0; k < info->parameter[2].threshold_count; k++) {
    for (j = 0; j < info->parameter[3].threshold_count; j++) {
      for (i = 0; i < info->parameter[3].threshold_count; i++) {
	cut = get_cut(i,j,k,l,CURRENT,info);
	printf("curr(%d,%d,%d,%d)=%lu,%lu ", i,j,k,l,cut->correct, cut->incorrect);
      } /* for i */
    } /* for j */
    printf("\n");
  } /* for k */

  printf("\n");
} /* display_cubes */

#ifdef DISPLAY_BASES
/***************************************************************************
 * display_number_of_correct_bases_in_each_bin
 *
 * called by: main
 * calls: get_bin
 *
 * ASSUMES: parameter_count = 4
 * 
 ***************************************************************************/
void
display_number_of_correct_bases_in_each_bin(PARAMETER *parameter, BIN *bin)
{
  int i, j, k, l;
  int *index;
  BIN *current_bin;

  index = (int *) malloc(sizeof(int)*PARAMETER_COUNT);
  if (index == NULL) {
    printf("couln't malloc index in display_number_of_correct_bases_in_each_bin\n");
    exit(-1);
  }
  for (i = 0; i < parameter[0].threshold_count; i++) {
    for (j = 0; j < parameter[1].threshold_count; j++) {
      for (k = 0; k < parameter[2].threshold_count; k++) {
	for (l = 0; l < parameter[3].threshold_count; l++) {
	  index[0] = i;
	  index[1] = j;
	  index[2] = k;
	  index[3] = l;
	  current_bin = get_bin(bin, parameter, index);
	  printf("bin(%d,%d,%d,%d)=(%f,%f,%f,%f) correct=%lu incorrect=%lu\n",
		 i,j,k,l,
		 parameter[0].threshold[i],
		 parameter[1].threshold[j],
		 parameter[2].threshold[k],
		 parameter[3].threshold[l],
		 current_bin->correct,
		 current_bin->incorrect);
	} /* for l */
      } /* for k */
    } /* for j */
  } /* for i */
  free(index);
} /* display_number_of_correct_bases_in_each_bin */
#endif



int
main(int argc, char *argv[])
{
    BASE         *base;
    BIN          *bin;
    PARAMETER     parameter[PARAMETER_COUNT];
    int           i, ij;
    unsigned long initial_base_room, base_count;
    unsigned int  threshold_count;
    time_t        t1, t2;
    int Num_params[4][50], Num_zeros[4], Num_psrbin[10];
    double pre, phr3, phr7, psr, psrbin[10], iheight;
    int spos, cpos, match;
    char schar, cchar;
    int corr[4][100], incorr[4][100], j,qv[4][100];
 

    /*    int max_threshold_count;*/

    t1 = time((time_t)NULL);

    /* Set defaults */
    Verbose = 4;
    Compress = 1;
    initial_base_room = BASE_COUNT_SCALE;

    opterr = 0;
    while ((i = getopt(argc, argv, "b:cCqv")) != EOF) {
	switch(i) {
	case 'b':
	    if ((sscanf(optarg, "%lu", &initial_base_room) != 1)
	    || (initial_base_room < 1))
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
	case 'q':
	    Verbose = 0;
	    break;
	case 'v':
	    Verbose++;
	    break;
	case '?':
	default:
	    show_usage(argc, argv);
	    exit(2);
	}
    }
    if ((optind + 1 != argc)
    || (sscanf(argv[optind], "%u", &threshold_count) != 1)
    || (threshold_count < 2))
    {
        show_usage(argc, argv);
        exit(2);
    }
    setbuf(stderr, NULL);	/* so that any status comes out w/o delay */

    if (Verbose > 1) {
	(void)fprintf(stderr, "starting with room for %lu bases\n",
			initial_base_room);
    }
    base = get_bases(initial_base_room, &base_count, 0);
    t2 = time((time_t)NULL);
    if (Verbose) {
	(void)fprintf(stderr, "%lu bases read in %ld sec\n", base_count,
			(long)(t2 - t1));
    }
    if (base_count == 0) {
        (void)fputs("no valid bases\n", stderr);
        exit(1);
    }

    t1 = t2;
    for (i = 0; i < PARAMETER_COUNT; i++) {
	parameter[i].threshold_count = threshold_count;
	parameter[i].threshold = (double *)malloc(
					    threshold_count * sizeof(double));
	if (parameter[i].threshold == NULL) {
	    (void)fputs("couldn't malloc thresholds\n", stderr);
	    exit(1);
	}
    }
    get_thresholds(base, base_count, parameter);
    parameter[0].dimension = 1;
    for (i = 1; i < PARAMETER_COUNT; i++) {
	parameter[i].dimension =
		    parameter[i-1].dimension * parameter[i-1].threshold_count;

    }

    /*Calculate and print into file quality values of each bin for each parameters*/

    for (i=0; i<4; i++){
      for (j=0; j<parameter[i].threshold_count; j++){
	corr[i][j]=0;
	incorr[i][j]=0;
	qv[i][j]=0;
      }}

	rewind(stdin);
	while (getbase(&spos, &schar, &cpos, &cchar, &phr3,
		 &phr7, &psr, &pre, &match, &iheight)){

	      if (phr3<=parameter[0].threshold[0]){
		if (match) corr[0][0]++; 
		else incorr[0][0]++;
	      } else {
		for (j=0; j<parameter[0].threshold_count-1; j++){
		  if (phr3<=parameter[0].threshold[j+1] && phr3>parameter[0].threshold[j]){
		    if (match) corr[0][j+1]++; 
		    else incorr[0][j+1]++;
		  }
		}
	      }
	      if (phr7<=parameter[1].threshold[0]){
		if (match) corr[1][0]++; 
		else incorr[1][0]++;
	      } else {
		for (j=0; j<parameter[1].threshold_count-1; j++){
		  if (phr7<=parameter[1].threshold[j+1] && phr7>parameter[1].threshold[j]){
		    if (match) corr[1][j+1]++; 
		    else incorr[1][j+1]++;
		  }
		}
	      }

	      if (psr<=parameter[2].threshold[0]){
		if (match) corr[2][0]++; 
		else incorr[2][0]++;
	      } else {
		for (j=0; j<parameter[2].threshold_count-1; j++){
		  if (psr<=parameter[2].threshold[j+1] && 
		      psr>parameter[2].threshold[j]){
		    if (match) corr[2][j+1]++; 
		    else incorr[2][j+1]++;
		  }
		}
	      }

	      if (pre<=parameter[3].threshold[0]){
		if (match) corr[3][0]++; 
		else incorr[3][0]++;
	      } else {
		for (j=0; j<parameter[3].threshold_count-1; j++){
		  if (pre<=parameter[3].threshold[j+1] && 
		      pre>parameter[3].threshold[j]){
		    if (match) corr[3][j+1]++; 
		    else incorr[3][j+1]++;
		  }
		}
	      }
	}



	for (i=0; i< 4; i++){
	(void)printf("Parameter # %d\n",i);
	(void)printf("%10s  %10s  %10s  %10s   %10s\n", "Bin #","Threshold", 
		     "# incorrect", "# total", "QV");
	  for (j=0; j<parameter[i].threshold_count; j++){
	    qv[i][j]=(int)(-(floor)(10.*log10((double)incorr[i][j]/((double)(corr[i][j]+incorr[i][j])))+0.5));
	    (void)printf("%10d  %10.4f  %10d  %10d   %10d\n", j, 
			 parameter[i].threshold[j], incorr[i][j], corr[i][j]+incorr[i][j], qv[i][j]);
	    
	  }
	}
	      

		  


  for (ij=0; ij<4; ij++){
    Num_zeros[ij]=0;
      for (i=0; i<parameter[0].threshold_count; i++){

	Num_params[ij][i]=0;
      }}

  for (ij=0; ij<5; ij++){
    Num_psrbin[ij]=0;
  }
  psrbin[0]=1.444444444444444;
  psrbin[1]=1.454545454545454;
  psrbin[2]=1.500000000000000;
  psrbin[3]=1.555555555555555;
  psrbin[4]=1.666666666666666;

	rewind(stdin);
	while (getbase(&spos, &schar, &cpos, &cchar, &phr3,
		 &phr7, &psr, &pre, &match, &iheight)){
	  if (phr3<=parameter[0].threshold[0]){
	    Num_params[0][0]++;
	    if (phr3==0.0) Num_zeros[0]++;
	  }
	  else {
	    for (i=1; i<parameter[0].threshold_count; i++){
	      if (phr3<=parameter[0].threshold[i] && phr3>parameter[0].threshold[i-1])
		Num_params[0][i]++;
	    }}

	  if (phr7<=parameter[1].threshold[0]){
	    Num_params[1][0]++;
	    if (phr7==0.0) Num_zeros[1]++;
	  }
	  else {
	    for (i=1; i<parameter[1].threshold_count; i++){
	      if (phr7<=parameter[1].threshold[i] && phr7>parameter[1].threshold[i-1])
		Num_params[1][i]++;
	    }}

	  if (pre<=parameter[3].threshold[0]){
	    Num_params[2][0]++;
	    if (pre==0.0) Num_zeros[2]++;
	  }
	  else {
	    for (i=1; i<parameter[3].threshold_count; i++){
	      if (pre<=parameter[3].threshold[i] && pre>parameter[3].threshold[i-1])
		Num_params[2][i]++;
	    }}

	  if (psr<=parameter[2].threshold[0]){
	    Num_params[3][0]++;
	    if (psr==0.0) Num_zeros[3]++;
	      if (pre>10.0)
	      	(void)fprintf(stderr, "\n%d %d %f %f %f %f\n", spos, cpos, phr3, phr7, psr, pre);
	  }
	  else {
	    for (i=1; i<parameter[2].threshold_count; i++){
	      if (psr<=parameter[2].threshold[i] && psr>parameter[2].threshold[i-1])
		Num_params[3][i]++;

	    }}
	  	      for (ij=0; ij<5; ij++)
		if (psr<=psrbin[ij]+1.0e-14 && psr>=psrbin[ij]-1.0e-14) Num_psrbin[ij]++;

	}
	/*
    (void)printf("%10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s\n","phr3","Num_phr3","phr7","Num_phr7","pre","Num_psr","psr","Num_pre");
(void)printf("%10.4f  %10d  %10.4f  %10d  %10.4f  %10d  %10.4f  %10d\n",
		     0.0000, Num_zeros[0], 0.0000, Num_zeros[1],
    		0.0000, Num_zeros[2], 0.0000, Num_zeros[3]);


 max_threshold_count=0;
 for (i=0; i<4; i++){
   max_threshold_count = (max_threshold_count > parameter[i].threshold_count) ? max_threshold_count:parameter[i].threshold_count;
 }
    for (i = 0; i < max_threshold_count; i++) {
	(void)printf("%10.4f  %10d  %10.4f  %10d  %10.4f  %10d  %10.4f  %10d\n",
		     parameter[0].threshold[i], Num_params[0][i], parameter[1].threshold[i], Num_params[1][i],
    		parameter[2].threshold[i], Num_params[2][i], parameter[3].threshold[i], Num_params[3][i]);
    }
    */
    for (ij=0; ij<5; ij++){
      (void)printf("%10.4f  %10d  \n", psrbin[ij], Num_psrbin[ij]);
    }

    t2 = time((time_t)NULL);
    if (Verbose) {
	(void)fprintf(stderr, "thresholds determined in %d sec\n", 
        (int)(t2 - t1));
    }

    t1 = t2;
    bin = (BIN *)malloc(parameter[0].threshold_count
			* parameter[1].threshold_count
			* parameter[2].threshold_count
			* parameter[3].threshold_count
			* sizeof(BIN));
    if (bin == NULL) {
        fputs("couldn't malloc bin\n", stderr);
        exit(1);
    }

    count_number_of_correct_bases_in_each_bin(base, base_count, parameter, bin);
    t2 = time((time_t)NULL);
    free(base);
    if (Verbose) {
	(void)fprintf(stderr, "bins populated in %d sec\n", (int)(t2 - t1));
    }
	(void)fprintf(stderr, "PARAMETER_COUNT=%d\n",PARAMETER_COUNT );

    free(bin);
    for (i = 0; i < PARAMETER_COUNT; i++) {
        free(parameter[i].threshold);
    }
    return 0;

}


