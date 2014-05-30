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
 * 2.3 2003/11/06 18:18:45
 */

ABIError SCF_Open(void *, size_t);
ABIError SCF_Close(void *);


void SCF_NumAnalyzedData(long *);
void SCF_NumBases(long *);
void SCF_Bases(char *);
void SCF_PeakLocations(short *);
void SCF_AnalyzedData(short, int *);
void SCF_SCFVersion(char *);

void delta_samples1(unsigned char *, long);
void delta_samples2(unsigned short *, long);

unsigned int swapuint(unsigned int);
int writeuint(FILE *, unsigned int);

unsigned short swapushort(unsigned short);
int writeushort(FILE *, unsigned short);

int writeuchar(FILE *, unsigned char);
