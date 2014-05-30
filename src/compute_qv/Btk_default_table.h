
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
 * 2.9 2003/11/06 18:18:37
 */


/*
 *  Btk_default_table.h  $Revision: 1.7 $
 */


extern BtkLookupTable *
Btk_get_3700pop5_table(void);

extern BtkLookupTable *
Btk_get_3700pop6_table(void);

extern BtkLookupTable *
Btk_get_3100pop6_table(void);

extern BtkLookupTable *
Btk_get_3730pop7_table(void);

extern BtkLookupTable *
Btk_get_ABI_3xx_recalln_table(void);

extern BtkLookupTable *
Btk_get_ABI_KB_recalln_table(void);

extern BtkLookupTable *
Btk_get_mbace_table(void);
