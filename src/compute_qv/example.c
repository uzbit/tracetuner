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
 */

/**  $Id: example.c,v 1.7 2009/01/01 14:15:33 gdenisov Exp $       
 **  example.c - Sample application program to output .qual files from
 **              ABI sample file inputs using the Tracetuner library.
 **/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdint.h>

#include "Btk_qv.h"
#include "Btk_qv_data.h"
#include "Btk_lookup_table.h"
#include "Btk_default_table.h"
#include "util.h"
#include "train.h"
#include "Btk_compute_qv.h"
#include "Btk_match_data.h"
#include "Btk_qv_io.h"
#include "ABI_Toolkit.h"

int
main(int argc, char *argv[])
{
    int i, n;
    char *lookup_table, *smp, *smptail;
    BtkLookupTable *table;
    int nbases, nvals, filetype = -1;
    char *bases, *chemistry = "";
    char status_code[100];
    int *locations, *vals[4]; 
    uint8_t *qv = NULL;

    BtkMessage  msg;
    Options     options;
    Results     results;

    if(argc < 2) {
        fprintf(stderr, "usage: %s <samplefiles...>\n", argv[0]);
        exit(0);
    }

    /*
     *  If a non-standard lookup table is specified as an environment
     *  variable use it, otherwise NULL gets the default.
     */
    lookup_table = getenv("LOOKUP_TABLE");

    if (lookup_table != NULL) {
        if (strcmp(lookup_table, "3730pop7") == 0) {
            table = Btk_get_3730pop7_table();
        }
        else if (strcmp(lookup_table, "3700pop5") == 0) {
            table = Btk_get_3700pop5_table();
        }
        else if (strcmp(lookup_table, "3700pop6") == 0) {
            table = Btk_get_3700pop6_table();
        }
        else {
            if ((table = Btk_read_lookup_table(lookup_table)) == NULL) {
                fprintf(stderr, "Couldn't read lookup table '%s'.\n", 
                        lookup_table);
                exit(1);
            }
        }
    }
    else {
        table = NULL;
    }

    for(n=1; n < argc; n++) {
        smp = argv[n];
        if ((smptail = strrchr(smp, '/')) != NULL)
            smptail++;
        else
            smptail = smp;

        /* Setting default options 
         */
        options.nocall         = 0;
        options.recalln        = 0;
        options.edited_bases   = 0;
        options.gauss          = 1;
        options.tip_dir[0]     = '\0';
        options.tal_dir[0]     = '\0';
        options.tab_dir[0]     = '\0';
        options.het            = 0;
        options.time           = 0;
        options.min_ratio      = (float)0.15;
        options.Verbose        = 3;
        options.inp_phd        = 0;
        options.inp_phd_dir[0] = '\0';
        options.file_name[0]   = '\0';
        options.lut_type       = ABI3730pop7;
       
      
        if (Btk_read_sample_file(smp, &nbases, &bases, 0, &locations, 
            &qv, &nvals, &vals[0], &vals[1], &vals[2], &vals[3],
            NULL, &chemistry, status_code, &filetype, options, &msg) != 0) 
        {
            fprintf(stderr, "%s: couldn't read sample file\n", smptail);
            continue;
        }

        strcpy(options.file_name, smptail);

        if (table == NULL) {
            if (strstr(chemistry, "3730")) {
                table = Btk_get_3730pop7_table();
            }
            else if (strstr(chemistry, "POP5")) {
                table = Btk_get_3700pop5_table();
            }
            else if (strstr(chemistry, "POP6") && 
                     strstr(chemistry, "3700")) {
                table = Btk_get_3700pop6_table();
            }
            else if (strstr(chemistry, "3100")) {
                table = Btk_get_3100pop6_table();
            }
            else {
            fprintf(stderr,
                    "Can't select the lookup table automatically. \n"
                    "Using built-in ABI 3730 Pop-7 table.\n");
            table = Btk_get_3730pop7_table();
            }
        }

        if ( Btk_compute_qv(&nbases, &bases, &locations, &nvals, vals, 
                     "ACGT", table, &qv, options, &msg, &results) != 0) {
            fprintf(stderr, "%s: %s\n", smptail, msg.text);
            goto cleanup_a_file;
        }

        fprintf(stderr, "%s: %d bases. ", smp, nbases);
        fprintf(stderr, "QVs are output to .qual file\n");

        Btk_output_quality_values(NAME_FILES, smp, NULL, "", qv, nbases, 0, 
                                  nbases - 1, 0);

    cleanup_a_file:
        if (qv != NULL) {
            free(qv);
            qv = NULL;
        }
        if (bases != NULL) {
            free(bases);
            bases = NULL;
        }
        if (locations != NULL) {
            free(locations);
            locations = NULL;
        }
        for (i=0; i<4; i++) {
            if (vals[i] != NULL) {
                free(vals[i]);
                vals[i] = NULL;
            }
        }

    }

    Btk_destroy_lookup_table(table);

    return(0);
}
