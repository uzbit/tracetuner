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
 *  ABI_Toolkit.h $Revision: 1.7 $
 */

typedef short ABIError;

#define kNoError             0
#define kCantOpenFile       -1
#define kFileError          -2
#define kMemoryFull         -3
#define kDataNotFound       -5
#define kFileNotOpen        -6
#define kFileAlreadyOpen    -7
#define kWrongFileType      -8
#define kBadCatalogLocation -9

char *ABI_ErrorString(ABIError);

ABIError ABI_Open(void *, size_t);
ABIError ABI_Close(void *);

ABIError ABI_AnalyzedData(short, short, int *);
ABIError ABI_AnalysisVersion(long, char *, long *);
ABIError ABI_CalledBases(char *);
ABIError ABI_CalledPeakLocations(short *);
ABIError ABI_DyeIndexToBase(short, char *);
ABIError ABI_EditedBases(char *);
ABIError ABI_ConsensusBases(char *);
ABIError ABI_EditedPeakLocations(short *);
ABIError ABI_MobilityFile(long, char *, long *);
ABIError ABI_NumAnalyzedData(short, short, long *);
ABIError ABI_NumCalledBases(long *);
ABIError ABI_NumCalledPeakLocations(long *);
ABIError ABI_NumEditedBases(long *);
ABIError ABI_NumConsensusBases(long *);   
ABIError ABI_NumEditedPeakLocations(long *);
ABIError ABI_NumRawData(short, short, long *);
ABIError ABI_RawData(short, short, int *);
ABIError ABI_NumQualityValues(short *num_qvs);
ABIError ABI_BasecallerQualityValues(char *qv);
