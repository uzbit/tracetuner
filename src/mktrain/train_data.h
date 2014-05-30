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
 * $Id: train_data.h,v 1.8 2009/01/12 20:20:59 gdenisov Exp $                  
 */

#ifndef TRAINDATA_H
#define TRAINDATA_H 

int append_align_to_train_data(FILE *, int, double **, uint8_t *, double *, 
    double *, double *, Align *, int, int, Options *, BtkMessage *);
int append_to_train_data(FILE *, int, int *, int, double **, uint8_t *, 
    double *, double *, double *, Range, Range, Align *, Align *,  Align *, 
    int, int, int, Options *, BtkMessage *);
int append_align_to_trainphred_data(FILE *, Align *, int, int, uint8_t *,
    int , char *, uint8_t *, int *, BtkMessage *);
int append_to_trainphred_data(FILE *, int, int *, Range, Range, Align *,
    Align *, Align *, int, int, uint8_t *,  int, char *, uint8_t *, int *, 
    BtkMessage *);
void release1_sanger(char *, int *, int **, uint8_t *, int, double **, double *);
void release1_454(sffRead *read, int, double **params);
#endif
