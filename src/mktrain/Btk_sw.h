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
 * $Id: Btk_sw.h,v 1.6 2009/01/13 20:46:47 gdenisov Exp $                 
 */

#ifndef BTK_SW_H
#define BTK_SW_H

extern int Btk_sw_alignment(Align_params *ap, Contig *query, Contig *library,
               Contig *rev_comp_library, Align *align, Range align_range,
               int read_direction, int phd, BtkMessage* message);

#endif
