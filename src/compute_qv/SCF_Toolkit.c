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

/**    Copyright (c) 1999 Paracel Inc.  All rights reserved.
 **/

/*
 *  SCF_Toolkit.c 2.2
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ABI_Toolkit.h"    /* So SCF_Toolkit will know what ABIError is. */
#include "SCF_Toolkit.h"
#include "util.h"

extern unsigned long get_offset(unsigned char *);

static char *gFile = NULL;

void SCF_NumBases(long *num_bases)
{
     *num_bases = get_offset((unsigned char *) (gFile + 12));
}

void SCF_Bases(char *edited_bases)
{
     long num_bases, i;
     unsigned long bases_offset;
     char scf_version_string[5];
     double scf_version_number;

     SCF_SCFVersion(scf_version_string);
     scf_version_number = atof(scf_version_string);

     SCF_NumBases(&num_bases);
     bases_offset = get_offset((unsigned char *) (gFile + 24));

     if (scf_version_number < 2.9)
	  for (i = 0; i < num_bases; i++)
	       edited_bases[i] = *(gFile + bases_offset + (i * 12) + 8);
     else
     {
	  bases_offset += (num_bases * 8);

	  for (i = 0; i < num_bases; i++)
	       edited_bases[i] = *(gFile + bases_offset + i);
     }
}

void SCF_PeakLocations(short *edited_locs)
{
     long num_bases, i;
     unsigned long bases_offset;
     char scf_version_string[5];
     double scf_version_number;

     SCF_SCFVersion(scf_version_string);
     scf_version_number = atof(scf_version_string);

     SCF_NumBases(&num_bases);
     bases_offset = get_offset((unsigned char *) (gFile + 24));

     if (scf_version_number < 2.9)
	  for (i = 0; i < num_bases; i++)
	       edited_locs[i] = (short) get_offset((unsigned char *) 
					(gFile + bases_offset + (i * 12)));
     else
	  for (i = 0; i < num_bases; i++)
	       edited_locs[i] = (short) get_offset((unsigned char *) 
					(gFile + bases_offset + (i * 4)));
}

void SCF_NumAnalyzedData(long *num_data_points)
{
     *num_data_points = get_offset((unsigned char *) (gFile + 4));
}

void SCF_AnalyzedData(short dye, int *analyzed_array)
{
     long num_data_points, i;
     unsigned long samples_offset;
     unsigned long sample_size;
     char scf_version_string[5];
     double scf_version_number;
     unsigned char *buf1;
     unsigned short *buf2;

     SCF_SCFVersion(scf_version_string);
     scf_version_number = atof(scf_version_string);

     SCF_NumAnalyzedData(&num_data_points);

     samples_offset = get_offset((unsigned char *) (gFile + 8));
     sample_size = get_offset((unsigned char *) (gFile + 40));

     if (scf_version_number < 2.9) {
	  for (i = 0; i < num_data_points; i++)
	       if (sample_size == 1)
		    analyzed_array[i] = *((unsigned char *) gFile + 
					  samples_offset + (i * 4) + dye);
	       else 
	       {
		    analyzed_array[i] = *((unsigned char *) gFile + 
					  samples_offset +
					  (i * 8) + (dye * 2)) * 256;
		    analyzed_array[i] += *((unsigned char *) gFile + 
					   samples_offset +
					   (i * 8) + (dye * 2) + 1);
	       }
     } else
     {
	  if (sample_size == 1)
	  {
	       buf1 = (unsigned char *) malloc(num_data_points * 
					       sizeof(unsigned char));

	       memcpy(buf1, ((unsigned char *) gFile + 
			     samples_offset + (dye * num_data_points)),
		      num_data_points);

	       delta_samples1(buf1, num_data_points);

	       for (i = 0; i < num_data_points; i++)
		    analyzed_array[i] = buf1[i];

	       free(buf1);
	  }
	  else
	  {
	       buf2 = (unsigned short *) malloc(num_data_points * 
						sizeof(unsigned short));

	       for (i = 0; i < num_data_points; i++)
	       {
		    buf2[i] = *((unsigned char *) gFile + 
				samples_offset + (dye * num_data_points * 2)
				+ (i * 2)) * 256;
		    buf2[i] += *((unsigned char *) gFile + 
				 samples_offset + (dye * num_data_points * 2)
				 + (i * 2) + 1);
	       }

	       delta_samples2(buf2, num_data_points);

	       for (i = 0; i < num_data_points; i++)
		    analyzed_array[i] = buf2[i];

	       free(buf2);
	  }
     }
}

void delta_samples1(unsigned char samples[], long num_samples)
{
     long i;
     unsigned char p_sample = 0;

     for (i = 0; i < num_samples; i++)
     {
	  samples[i] = samples[i] + p_sample;
	  p_sample = samples[i];
     }

     p_sample = 0;

     for (i = 0; i < num_samples; i++)
     {
	  samples[i] = samples[i] + p_sample;
	  p_sample = samples[i];
     }
}

void delta_samples2(unsigned short samples[], long num_samples)
{
     long i;
     unsigned short p_sample = 0;

     for (i = 0; i < num_samples; i++)
     {
	  samples[i] = samples[i] + p_sample;
	  p_sample = samples[i];
     }

     p_sample = 0;

     for (i = 0; i < num_samples; i++)
     {
	  samples[i] = samples[i] + p_sample;
	  p_sample = samples[i];
     }
}

void SCF_SCFVersion(char *scf_version_string)
{
     int i;

     for (i = 0; i < 4; i++)
	  scf_version_string[i] = *((char *) gFile + 36 + i);

     scf_version_string[4] = '\0';
}

ABIError SCF_Open(void *file, size_t size)
{
     ABIError error = kNoError;

     if (gFile != NULL)
	  error = kFileAlreadyOpen;
     else
	  gFile = (char *) file;

     return error;
}

ABIError SCF_Close(void *file)
{
     ABIError error = kNoError;

     if (gFile != file)
	  error = kFileNotOpen;
     else
	  gFile = NULL;

     return error;
}

unsigned int swapuint(unsigned int i)
{
    unsigned int n;
    char *p;

    p = (char *) &n;
    p[0] = (i >> 24) & 0xff;
    p[1] = (i >> 16) & 0xff;
    p[2] = (i >> 8) & 0xff;
    p[3] = i & 0xff;

    /* hexadecimal number 0xff s represented in a binary form as: 
     * 00000000000000000000000011111111 
     * Thus, operation i & 0xff sets to zero first 24 bits and 
     * retains the last 8 bits
     */
    return n;
}

unsigned short swapushort(unsigned short i)
{
    unsigned short n;
    char *p;

    p = (char *) &n;
    p[0] = (i >> 8) & 0xff;
    p[1] = i & 0xff;

    return n;
}

int writeuint(FILE *fp, unsigned int i)
{
    unsigned int n;

    n = swapuint(i);

    if (fwrite((void *) &n, sizeof(unsigned int), 1, fp) != 1)
        return ERROR;

    return SUCCESS;
}

int writeushort(FILE *fp, unsigned short i)
{
    unsigned short n;

    n = swapushort(i);

    if (fwrite((char *) &n, sizeof(unsigned short), 1, fp) != 1)
        return ERROR;

    return SUCCESS;
}

int writeuchar(FILE *fp, unsigned char i)
{
    if (fwrite((char *) &i, sizeof(unsigned char), 1, fp) != 1)
        return ERROR;

    return SUCCESS;
}
