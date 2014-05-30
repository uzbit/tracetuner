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
 */

/*
 *  $Id: Btk_atod.c,v 1.5 2008/11/27 18:48:29 gdenisov Exp $ 
 */

#include <stdio.h>
#include <ctype.h>


/*
 * This function converts an ASCII string, from the standard input, of the
 * form 123.4559... into a double.  It can handle leading white space and an
 * optional leading sign, but does NOT handle scientific notation.  Its
 * synopsis is:
 *
 * success = Btk_atod(buf, dp)
 *
 * where
 *	buf	is a pointer to the string to parse, passed by REFERENCE
 *	dp	is a pointer in which to put the result
 *
 *	success	is 1 if a number was read
 *		   0 if the end of buf was reached without reading a number
 *		  -1 if a non-numeric character was encountered
 */
int
Btk_atod(char **buf, double *dp)
{
    int negative;
    double frac;


    /* Skip leading white space, if any */
    while (isspace((int)**buf)) {
	(*buf)++;
    }
    if (**buf == '\0') {
	return 0;
    }

    /* Check for leading sign */
    if (**buf == '-') {
	negative = 1;
	(*buf)++;
	if (**buf == '\0') {
	    return 0;
	}
    }
    else {
	negative = 0;
	if (**buf == '+') {
	    (*buf)++;
	    if (**buf == '\0') {
		return 0;
	    }
	}
    }
    if (!isdigit((int)**buf)) {
	return -1;
    }

    /* Get the integral part */
    *dp = 0.0;
    while (isdigit((int)**buf)) {
	*dp *= 10.0;
	if (**buf != '0') {
	    *dp += (**buf - '0');
	}
	(*buf)++;
	if (**buf == '\0') {
	    goto done;
	}
    }

    /* Get the fraction part, if any */
    if (**buf == '.') {
	(*buf)++;
	if (**buf == '\0') {
	    goto done;
	}
	frac = 0.1;
	while (isdigit((int)**buf)) {
	    if (**buf != '0') {
		*dp += (**buf - '0') * frac;
	    }
	    (*buf)++;
	    if (**buf == '\0') {
		goto done;
	    }
	    frac /= 10.0;
	}
    }

done:
    if (negative) {
	*dp = -(*dp);
    }
    return 1;

}
