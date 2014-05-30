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
 * 1.20 2003/11/06 18:18:38
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Btk_qv.h"
#include "Btk_lookup_table.h"
#include "Btk_default_table.h"
#include "Btk_atod.h"

#define MAXLINE	(1000)
#define CHUNK	(1000)
#define MAX_NUM_TPAR_THRESHOLDS	(100)

// --------------------------------------------------------------------
/*
 * This function reads in and parses a lookup table.  Its synopsis is:
 *
 * table = Btk_read_lookup_table(path)
 *
 * where
 *	path	is the name of the file containing the lookup table
 *
 *	table	is a BtkLookupTable allocated and populated by this function
 *		on success, NULL otherwise.
 *
 * If this function returns a table, then the caller should free up the
 * resources (when through with it) by calling Btk_destroy_lookup_table().
 */
BtkLookupTable *
Btk_read_lookup_table(char *path)
{
    FILE *fp;
    BtkLookupTable *table = NULL;
    int size;
    char linebuf[MAXLINE], *s, *s1, *s2, *s3, *s4;
    int qv, p1, p2, p3, p4;
    int version_correct = 0;
    char table_version[20] = "1.0";  

    if ((fp = fopen(path, "r")) == NULL) {
	return(table);
    }

    table = CALLOC(BtkLookupTable, 1);
    size = CHUNK;
    table->tpar = (TraceParamEntry *) malloc(sizeof(TraceParamEntry)
			       	* MAX_NUM_TPAR_THRESHOLDS);
    table->num_tpar_entries = 0;
    table->entries = (BtkLookupEntry *)malloc(sizeof(BtkLookupEntry) * size);
    table->num_lut_entries = 0;
    while (fgets(linebuf, MAXLINE, fp) != NULL) {
	/* Ignore all white space lines and comments */
	if (strspn(linebuf, " \t\r\n") == strlen(linebuf)) {
	    continue;
	}
	if (linebuf[0] == '#') {
            if  (sscanf(linebuf, "# Version %s", table_version) != 1)
            {
                continue;
            }
            else {
                if (strncmp(table_version, TT_VERSION, 6) == 0) {
                    version_correct = 1;
                    continue;
                }
            }
        }

        if (((linebuf[0] == '/') && (linebuf[1] == '*'))
	  || (linebuf[0] == ';')
          || (linebuf[0] == '#'))
	{
	    continue;
	}

	s = linebuf;
	if (!version_correct) {
            fprintf(stderr, 
            "\nLookup table version %s does not match the ttuner version %s\n",
            table_version, TT_VERSION);
            goto error_return;
        }

	// The trace parameter threshold entries have 4 tokens; whereas lut
	// entires have 6 tokens.
	s1 = strtok(s, " \t\n");
	if (s1 == NULL)
	    goto error_return;

	s2 = strtok(NULL, " \t\n");
	if (s2 == NULL)
	    goto error_return;

	s3 = strtok(NULL, " \t\n");
	if (s3 == NULL)
	    goto error_return;

	s4 = strtok(NULL, " \t\n");
	if (s4 == NULL)
	    goto error_return;

	s = strtok(NULL, " \t\n");
	if (s == NULL) {
	    // only 4 tokens: this is a tpar entry
	    // read the tpar thresholds
            if ((Btk_atod(&s1, &table->tpar[table->num_tpar_entries].phr3t) != 1)
    		    || (Btk_atod(&s2, &table->tpar[table->num_tpar_entries].phr7t) != 1)
    	    	    || (Btk_atod(&s3, &table->tpar[table->num_tpar_entries].psr7t) != 1)
    	    	    || (Btk_atod(&s4, &table->tpar[table->num_tpar_entries].prest) != 1))
	    {
	        goto error_return;
	    }
	    table->num_tpar_entries++;
	} else {
	    // has 5 tokens: lut entries
	    qv = atoi(s1);
	    p1 = atoi(s2);
	    p2 = atoi(s3);
	    p3 = atoi(s4);
	    p4 = atoi(s);

            s = strtok(NULL, " \t\n");

	    if (qv < 0 || p1 < 0 || p1 >= table->num_tpar_entries
		    || p2 < 0 || p2 >= table->num_tpar_entries
		    || p3 < 0 || p3 >= table->num_tpar_entries
		    || p4 < 0 || p4 >= table->num_tpar_entries)
	    {
	    	goto error_return;
	    }

	    table->entries[table->num_lut_entries].qval = qv;
	    table->entries[table->num_lut_entries].phr3i = p1;
	    table->entries[table->num_lut_entries].phr7i = p2;
	    table->entries[table->num_lut_entries].psr7i = p3;
	    table->entries[table->num_lut_entries].presi = p4;
	    table->num_lut_entries++;
	    if (table->num_lut_entries == size) {
	    	size += CHUNK;
	    	table->entries = REALLOC(table->entries, BtkLookupEntry, size);
	    }
	}
    }
    if (ferror(fp)) {
	goto error_return;
    }

    /* Close file and give back unused space */
    (void)fclose(fp);
    table->tpar = REALLOC(table->tpar,  TraceParamEntry,
			  table->num_tpar_entries);
    table->entries = REALLOC(table->entries,  BtkLookupEntry,
			     table->num_lut_entries);

    return(table);

error_return:
    (void)fclose(fp);
	FREE(table->tpar);
	FREE(table->entries);
	FREE(table);
    return(NULL);
}


/*
 * This function reclaims the storage taken up by a lookup table previously
 * returned by Btk_read_lookup_table().  Its synopsis is:
 *
 * Btk_destroy_lookup_table(table)
 *
 * where
 *	table	is a BtkLookupTable previously returned by
 *		Btk_read_lookup_table().
 */
void
Btk_destroy_lookup_table(BtkLookupTable *table)
{
    if (table == NULL)                     {  return;  }
    if (table == Btk_get_3730pop7_table()) {  return;  }
    if (table == Btk_get_3700pop5_table()) {  return;  }
    if (table == Btk_get_3700pop6_table()) {  return;  }
    if (table == Btk_get_3100pop6_table()) {  return;  }
    if (table == Btk_get_mbace_table())    {  return;  }

    FREE(table->entries);
    FREE(table->tpar);
    FREE(table);
}
