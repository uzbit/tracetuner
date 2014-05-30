# *************************************************************************
#  This file is part of TraceTuner, the DNA sequencing quality value,
#  base calling and trace processing software.
# 
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received (LICENSE.txt) a copy of the GNU General Public
#  License along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
# ************************************************************************/

#
# Copyright (c) 1999-2003 Paracel, Inc.  All rights reserved.
#
# $Id: Makefile,v 1.9 2008/12/19 15:16:35 gdenisov Exp $
#

include ../include.mk

.PHONY: all lut 
 
all:  lut 

lut:		$(RELDIR)/lut

INCTOOLSDIR = ../../..
INCDIR      = ../compute_qv
INCTRAINDIR = ../mktrain
CURDIR      = .
LIBS        = -lm -lz
INCCHKDIR   = ../mkchk

LUTOBJS		= $(OBJDIR)/lut.o $(OBJDIR)/select.o $(OBJDIR)/func_name.o \
                  $(OBJDIR)/get_thresholds.o $(OBJDIR)/check_data.o
TTLIB		= $(LIBDIR)/libtt.a
CFLAGS     += -I$(INCDIR)
CFLAGS     += -I$(INCTRAINDIR)
CFLAGS     += -I$(CURDIR)
CFLAGS     += -I$(INCCHKDIR)

$(OBJDIR)/%.o: $(INCCHKDIR)/%.c
	mkdir -p $(OBJDIR)
	$(COMPILE.c) $< -o $@

$(RELDIR)/lut: $(LUTOBJS) $(TTLIB) 
	@mkdir -p $(RELDIR)
	$(LINK.c) $(LUTOBJS) -o $@ $(LIBS) $(TTLIB) 

$(OBJDIR)/lut.o:	lut.c lut.h get_thresholds.h select.h func_name.h params.h \
			$(INCDIR)/Btk_atod.h $(INCDIR)/Btk_qv.h $(INCTRAINDIR)/train.h

$(DIRS):
	mkdir -p $@

clean:
	@/bin/rm -f $(LUTOBJS)
	@/bin/rm -f $(RELDIR)/lut
