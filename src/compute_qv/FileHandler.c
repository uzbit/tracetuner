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
 * 1.8 2003/11/06 18:18:44
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "ABI_Toolkit.h"
#include "SCF_Toolkit.h"
#include "SFF_Toolkit.h"
#include "FileHandler.h"
#include "Btk_qv.h"

ABIError F_Open(char *file_name, void **ptr, long *file_size, int *file_type)
{
    ABIError error = kNoError;
    FILE          *stream;
    size_t         size;

   *ptr = NULL;
   *file_size = size = 0;

    stream = fopen(file_name, "rb");
    if (stream == NULL)
        error = kCantOpenFile;

    if (error == kNoError)
        error = fseek(stream, 0, SEEK_END) ? kFileError : kNoError;

    if (error == kNoError)
        size = ftell(stream);

    if (error == kNoError)
        error = fseek(stream, 0, SEEK_SET) ? kFileError : kNoError;

    if (error == kNoError)
    {
       *file_size = size;
       *ptr = malloc(size);
        if (*ptr == NULL)
            error = kMemoryFull;
    }

    if (error == kNoError)
        error = size != fread(*ptr, 1, size, stream) ? kFileError : kNoError;

    if (error == kNoError)
    {
        if (strncmp((char *) *ptr, "ABIF", 4) == 0)
            *file_type = ABI;
        else if (strncmp((char *) *ptr, ".scf", 4) == 0)
            *file_type = SCF;
        else if (strncmp((char *) *ptr + 1, "ZTR", 3) == 0)
        {
            *file_type = ZTR;
        }
        else // check if SFF
        {
            error = kWrongFileType;
        }
    }

    if (error == kNoError)
    {
        if (*file_type == ABI)
            error = ABI_Open(*ptr, size);
        else if (*file_type == SCF)
            error = SCF_Open(*ptr, size);
        else if (*file_type == ZTR)
            error = kNoError;              
        else if (*file_type == SFF)
            error = kNoError;   
        else 
        {
            fprintf(stderr, "Unknown file type\n");
            error = kWrongFileType;
            return error;
        }
    }

    if (stream != NULL)
        fclose(stream);

     return error;
}

ABIError F_Close(void *ptr, int file_type)
{
     ABIError error;

     if (file_type == ABI)
	  error = ABI_Close(ptr);
     else
	  error = SCF_Close(ptr);

     free(ptr);

     return error;
}
