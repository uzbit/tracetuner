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

/**  $Revision: 1.6 $
 **  example.c - Sample application program to output .qual files from
 **              ABI sample file inputs using the Tracetuner library.
 **/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <errno.h>
#include <stdint.h>

#include "util.h"
#include "Btk_qv.h"
#include "Btk_qv_data.h"
#include "Btk_match_data.h"
#include "Btk_qv_io.h"

#define SUCCESS 0
#define ERROR -1

static char status_code[BUFLEN];

int
main(int argc, char *argv[])
{
    uint8_t *quality_values;
    int      edited_bases = 0, file_type = -1;
    char     file_name[256];
    char    *called_bases=NULL, *call_method=NULL, *chemistry=NULL;
    int      num_called_bases=0, *called_peak_locs, num_datapoints;
    int     *chromatogram[NUM_COLORS];
    BtkMessage *message = NULL;
    Options  options;
    FILE    *fp;

    if(argc < 2) {
        fprintf(stderr, "usage: %s <sample_file>\n", argv[0]);
        exit(0);
    }

    if ((fp = fopen(argv[1], "r")) == NULL)
    {
        sprintf(message->text,
                "could not open file for reading, errno = %d", errno);
        return ERROR;
    }

    sprintf(file_name, "%s", argv[1]);
    fclose(fp);

    if (file_name[0] == '\0') {
        return SUCCESS;
    }
    options.inp_phd = 0;
    fprintf(stderr, "File name   = %s\n", file_name);
  
    if (Btk_read_sample_file(file_name, &num_called_bases, &called_bases,
        edited_bases, &called_peak_locs, &quality_values,
        &num_datapoints, &chromatogram[0], &chromatogram[1],
        &chromatogram[2], &chromatogram[3], &call_method, &chemistry,
        status_code, &file_type, options, message) != SUCCESS)
    {
        goto error;
    }
    fprintf(stderr, "Chemistry   = %s\n", 
        strlen(chemistry)>0 ? chemistry : "unknown");
    fprintf(stderr, "Call method = %s\n", 
        strlen(call_method)>0 ? call_method : "unknown");
    fprintf(stderr, "Num called bases = %d\n", num_called_bases);
    fprintf(stderr, "Num datapoints   = %d\n", num_datapoints);
    return(0);
error:
    return(1);
}
