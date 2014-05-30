/**************************************************************************
 * This file is part of TraceTuner, the DNA sequencing quality value,
 * base calling and trace processing software.
 *
 * Copyright (c) 1999-2003 Paracel, Inc.  All rights reserved.
 *
 * 2.6 2003/11/06 18:18:46
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
 *  ABI_Toolkit.c $Revision: 1.7 $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ABI_Toolkit.h"

static char   *gFile = NULL;
unsigned long  dirloc;
unsigned long  tag_count;

unsigned long get_offset(unsigned char *cptr)
{
    int byte;
    unsigned long total = 0;

    byte = *cptr;
    total = byte * 256 * 256 * 256;
    byte = *(cptr + 1);
    total += byte * 256 * 256;
    byte = *(cptr + 2);
    total += byte * 256;
    byte = *(cptr + 3);
    total += byte;

    return total;
}

ABIError ABI_Open(void *file, size_t size)
{
    ABIError error = kNoError;

    if (gFile != NULL)
        error = kFileAlreadyOpen;
    else
    {
        gFile = (char *) file;

        dirloc = get_offset((unsigned char *)((unsigned char *)file + 26));

        if (dirloc > size || dirloc < 128)
            error = kBadCatalogLocation;

        tag_count = get_offset((unsigned char *)((unsigned char *)file + 18));
    }
    return error;
}

ABIError ABI_Close(void *file)
{
    ABIError error = kNoError;

    if (gFile != file)
        error = kFileNotOpen;
    else
        gFile = NULL;

    return error;
}

char *find_dir_entry(char *tag, int id)
{
    int i;
    char curtag[4];
    int curid;
    char *tagptr = gFile + dirloc;

    for (i = 0; i < (int)tag_count; i++)
    {
        curtag[0] = *tagptr;
        curtag[1] = *(tagptr + 1);
        curtag[2] = *(tagptr + 2);
        curtag[3] = *(tagptr + 3);

        if (!strncmp(tag, curtag, 4))
        {
            curid = *((unsigned char *) tagptr + 7);
            if (curid == id)
                return tagptr;
        }
        tagptr += 28;
    }
    return NULL;
}

ABIError ABI_NumCalledBases(long *num_bases)
{
    ABIError error = kDataNotFound;
    char *location;

   *num_bases = 0;
    location = find_dir_entry("PBAS", 2);
    if (location)
    {
       *num_bases = get_offset((unsigned char *)(location + 12));
        error = kNoError;
    }
    return error;
}

ABIError ABI_NumEditedBases(long *num_bases)
{
    ABIError error = kDataNotFound;
    char *location;

   *num_bases = 0;
    location = find_dir_entry("PBAS", 1);
    if (location)
    {
       *num_bases = get_offset((unsigned char *)(location + 12));
        error = kNoError;
    }
    return error;
}

ABIError ABI_NumConsensusBases(long *num_bases)
{
    ABIError error = kDataNotFound;
    char *location;

   *num_bases = 0;
    location = find_dir_entry("aSEQ", 1);
    if (location)
    {
       *num_bases = get_offset((unsigned char *)(location + 12));
        error = kNoError;
    }
    return error;
}

ABIError ABI_NumCalledPeakLocations(long *num_locs)
{
    ABIError error = kDataNotFound;
    char *location;

   *num_locs = 0;
    location = find_dir_entry("PLOC", 2);
    if (location)
    {
       *num_locs = get_offset((unsigned char *)(location + 12));
        error = kNoError;
    }
    return error;
}

ABIError ABI_NumEditedPeakLocations(long *num_locs)
{
    ABIError error = kDataNotFound;
    char *location;

   *num_locs = 0;
    location = find_dir_entry("PLOC", 1);
    if (location)
    {
       *num_locs = get_offset((unsigned char *)(location + 12));
        error = kNoError;
    }
    return error;
}

ABIError ABI_MobilityFile(long buffer_size, char *file_name, long *data_size)
{
    ABIError error = kDataNotFound;
    char *location;

    location = find_dir_entry("PDMF", 2);
    if (location)
    {
       *data_size = get_offset((unsigned char *)(location + 12));

        if (*data_size > buffer_size)
        {
            fprintf(stderr, 
                "Warning: Mobility file name will be truncated ");
            fprintf(stderr, "- need a larger buffer.\n");
            memcpy(file_name, gFile + 
                get_offset((unsigned char *)(location + 20)) + 1, 
                buffer_size - 1);
            file_name[buffer_size] = '\0';
        }
        else
        {
            memcpy(file_name, gFile + 
                get_offset((unsigned char *)(location + 20)) + 1, *data_size-1);
	    file_name[*data_size - 1] = '\0';
	}
        error = kNoError;
    }
    return error;
}

ABIError ABI_AnalysisVersion(long buffer_size, char *software, long *data_size)
{
    ABIError error = kDataNotFound;
    char *location;

    location = find_dir_entry("SVER", 2);
    if (location)
    {
       *data_size = get_offset((unsigned char *)(location + 12));

        if (*data_size > buffer_size)
        {
            fprintf(stderr, "Warning: Analysis version will be truncated - ");
            fprintf(stderr, "need a larger buffer.\n");
            memcpy(software, gFile + 
                get_offset((unsigned char *)(location+20)) + 1, buffer_size-1);
            software[buffer_size] = '\0';
        }
        else
        {
            memcpy(software, gFile + 
                get_offset((unsigned char *)(location+20)) + 1, *data_size-1);
            software[*data_size - 1] = '\0';
        }
        error = kNoError;
    }
    return error;
}

ABIError ABI_CalledBases(char *called_bases)
{
    ABIError error = kDataNotFound;
    char *location;
    unsigned long base_count;
    unsigned long base_start;

    location = find_dir_entry("PBAS", 2);
    if (location)
    {
        base_count = get_offset((unsigned char *)(location + 12));
        base_start = get_offset((unsigned char *)(location + 20));

        memcpy(called_bases, gFile + base_start, base_count);
        /* called_bases[base_count] = '\0';*/

        error = kNoError;
    }
    return error;
}

ABIError ABI_EditedBases(char *edited_bases)
{
    ABIError error = kDataNotFound;
    char *location;
    unsigned long base_count;
    unsigned long base_start;

    location = find_dir_entry("PBAS", 1);
    if (location)
    {
        base_count = get_offset((unsigned char *)(location + 12));
        base_start = get_offset((unsigned char *)(location + 20));

        memcpy(edited_bases, gFile + base_start, base_count);
        /* edited_bases[base_count] = '\0';*/

	error = kNoError;
    }
    return error;
}

ABIError ABI_ConsensusBases(char *cons_bases)
{
    ABIError error = kDataNotFound;
    char *location;
    unsigned long base_count;
    unsigned long base_start;

    location = find_dir_entry("aSEQ", 1);
    if (location)
    {
        base_count = get_offset((unsigned char *)(location + 12));
        base_start = get_offset((unsigned char *)(location + 20));

        memcpy(cons_bases, gFile + base_start, base_count);
        /* cons_bases[base_count] = '\0';*/

        error = kNoError;
    }
    return error;
}

ABIError ABI_NumQualityValues(short *num_qvs)
{
    ABIError error = kDataNotFound;
    char *location;

    location = find_dir_entry("PCON", 1);
    if (location)
    {
       *num_qvs = get_offset((unsigned char *)(location + 12));
        error = kNoError;
    }
    return error;
}

ABIError ABI_BasecallerQualityValues(char *qv)
{
    ABIError error = kDataNotFound;
    char          *location;
    unsigned long  qv_count;
    unsigned long  qv_start;

    location = find_dir_entry("PCON", 1);
    if (location)
    {
        qv_count = get_offset((unsigned char *)(location + 12));
        qv_start = get_offset((unsigned char *)(location + 20));

        memcpy(qv, gFile + qv_start, qv_count);

        error = kNoError;
    }
    return error;
}

ABIError ABI_CalledPeakLocations(short *called_locs)
{
    ABIError error = kDataNotFound;
    char *location;
    unsigned long num_peaks;
    unsigned long peak_start;
    int i;

    location = find_dir_entry("PLOC", 2);
    if (location)
    {
        num_peaks = get_offset((unsigned char *)(location + 12));
        peak_start = get_offset((unsigned char *)(location + 20));

        for (i = 0; i < (int)num_peaks; i++)
        {
            called_locs[i] = 
               *((unsigned char *) (gFile + peak_start + (i * 2))) * 256;
            called_locs[i] += 
               *((unsigned char *) (gFile + peak_start + (i * 2) + 1));
        }

        error = kNoError;
    }
    return error;
}

ABIError ABI_EditedPeakLocations(short *edited_locs)
{
    ABIError error = kDataNotFound;
    char *location;
    unsigned long num_peaks;
    unsigned long peak_start;
    int i;

    location = find_dir_entry("PLOC", 1);
    if (location)
    {
        num_peaks = get_offset((unsigned char *)(location + 12));
        peak_start = get_offset((unsigned char *)(location + 20));

        for (i = 0; i < (int)num_peaks; i++)
        {
            edited_locs[i] = 
               *((unsigned char *) (gFile + peak_start + (i * 2))) * 256;
            edited_locs[i] += 
               *((unsigned char *) (gFile + peak_start + (i * 2) + 1));
        }
        error = kNoError;
    }
    return error;
}

ABIError ABI_DyeIndexToBase(short index, char *c)
{
    ABIError error = kDataNotFound;
    char *location;

   *c = '\0';

    location = find_dir_entry("FWO_", 1);
    if (location)
    {
        if (index >= 1 && index <= 4)
	{
           *c = *(location + 19 + index);
            error = kNoError;
        }
    }
    return error;
}

ABIError ABI_NumRawData(short lane, short dye, long *num_data_points)
{
    ABIError error = kDataNotFound;
    char *location;

   *num_data_points = 0;
    lane =0;

    location = find_dir_entry("DATA", dye <= 4 ? dye : 100 + dye);
    if (location)
    {
       *num_data_points = get_offset((unsigned char *)(location + 12));
        error = kNoError;
    }
    return error;
}

ABIError ABI_NumAnalyzedData(short lane, short dye, long *num_data_points)
{
    ABIError error = kDataNotFound;
    char *location;

   *num_data_points = 0;
    lane =0;

    location = find_dir_entry("DATA", dye <= 4 ? dye + 8 : 200 + dye);
    if (location)
    {
       *num_data_points = get_offset((unsigned char *)(location + 12));
        error = kNoError;
    }
    return error;
}

ABIError ABI_RawData(short lane, short dye, int *raw_array)
{
    ABIError error = kDataNotFound;
    char *location;
    unsigned long num_data;
    unsigned long data_start;
    int i;

    lane =0;

    location = find_dir_entry("DATA", dye <= 4 ? dye : 100 + dye);
    if (location)
    {
        num_data = get_offset((unsigned char *)(location + 12));
        data_start = get_offset((unsigned char *)(location + 20));

        for (i = 0; i < (int)num_data; i++)
        {
            raw_array[i] =
               *((unsigned char *) (gFile + data_start + (i * 2))) * 256;
            raw_array[i] +=
               *((unsigned char *) (gFile + data_start + (i * 2) + 1));
        }
        error = kNoError;
    }
    return error;
}

ABIError ABI_AnalyzedData(short lane, short dye, int *analyzed_array)
{
    ABIError error = kDataNotFound;
    char *location;
    unsigned long num_data;
    unsigned long data_start;
    int i;

    lane =0;

    location = find_dir_entry("DATA", dye <= 4 ? dye + 8 : 200 + dye);
    if (location)
    {
        num_data = get_offset((unsigned char *)(location + 12));
        data_start = get_offset((unsigned char *)(location + 20));

        for (i = 0; i < (int)num_data; i++)
        {
	    analyzed_array[i] = 
               *((unsigned char *) (gFile + data_start + (i * 2))) * 256;
	    analyzed_array[i] += 
               *((unsigned char *) (gFile + data_start + (i * 2) + 1));
	}
	error = kNoError;
    }
    return error;
}

/*
 *  Return a user readable string associated with each ABI error number.
 *  K can be between 0 and -8 in the current code.
 */
char *ABI_ErrorString(ABIError k)
{
    static char line[80];

    switch(k) {
        case kNoError:
            return "Success";
        case kCantOpenFile:
            return "Unable to open file";
        case kFileError:
            return "Error seeking or reading in file";
        case kMemoryFull:
            return "Unable to open sample file due to insufficient memory";
        case kDataNotFound:
            return "Sample file did not contain expected data elements";
        case kFileNotOpen:
            return "File is not open";
        case kFileAlreadyOpen:
            return "File has already been openned";
        case kWrongFileType:
            return "Not an ABI, SCF or SFF file";
        case kBadCatalogLocation:
            return "Sample file corrupt - bad catalog location";
        default: 
            sprintf(line,"Unknown error code %d",k);

        return line;
    }
}
