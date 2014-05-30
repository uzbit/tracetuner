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
 * 2.2 2003/11/06 18:18:43
 */

#define TT_SCF_MAGIC 779314022        

typedef struct
{
    unsigned int magic_number;
    unsigned int samples;
    unsigned int samples_offset;
    unsigned int bases;
    unsigned int bases_left_clip;
    unsigned int bases_right_clip;
    unsigned int bases_offset;
    unsigned int comments_size;
    unsigned int comments_offset;
    char version[4];
    unsigned int sample_size;
    unsigned int code_set;
    unsigned int private_size;
    unsigned int private_offset;
    unsigned int spare[18];
} SCF_Header;

typedef struct
{
    unsigned int peak_index;
    unsigned char prob_A;
    unsigned char prob_C;
    unsigned char prob_G;
    unsigned char prob_T;
    char base;
    unsigned char spare[3];
} SCF_Bases_Rec;

typedef struct
{
    unsigned char sample_A;
    unsigned char sample_C;
    unsigned char sample_G;
    unsigned char sample_T;
} TT_Samples1;

typedef struct
{
    unsigned short sample_A;
    unsigned short sample_C;
    unsigned short sample_G;
    unsigned short sample_T;
} TT_Samples2;
