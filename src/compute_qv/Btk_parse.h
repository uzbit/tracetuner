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
 * 1.7 2003/11/06 18:18:39
 */


/*
 *  Btk_parse.h  $Revision: 1.4 $
 */

int Btk_parse_sample_file(
    char *,	/* the head name for the files; will turn into head.bases,
		   head.pos, head.a, head.c, head.g, head.t */
    int *,	/* pointer to an int that will be set to the # of bases read */
    char **,	/* pointer to an char pointer that will be set to an array
		   allocated via malloc() */
    int **, 	/* pointer to an int pointer that will be set to an array
		   allocated via malloc() */
    int *,	/* pointer to an int that will be set to the number of values
		   for each color */
    int **,	/* pointer to an int pointer that will be set to an array
		   allocated via malloc() for A values */
    int **,	/* pointer to an int pointer that will be set to an array
		   allocated via malloc() for C values */
    int **,	/* pointer to an int pointer that will be set to an array
		   allocated via malloc() for G values */
    int **	/* pointer to an int pointer that will be set to an array
		   allocated via malloc() for T values */
);
