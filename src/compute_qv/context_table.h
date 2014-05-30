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
 * 2.3 2003/11/06 18:18:46
 */

#ifndef CONTEXT_TABLE_H__
#define CONTEXT_TABLE_H__

typedef struct _context_table {
        int     dimension;
        double *weights;
}  ContextTable;  

extern double         weight_from_reverse_context( const char *, ContextTable * );
extern double         weight_from_context( const char *, ContextTable * );
extern ContextTable* read_context_table( char * );
extern void           destroy_context_table( ContextTable * );

#endif          /* CONTEXT_TABLE_H__ */
