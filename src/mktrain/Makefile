#**************************************************************************
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
#*************************************************************************
#
# Copyright (c) 1999-2003 Paracel, Inc.  All rights reserved.
#
# $Id: Makefile,v 1.9 2008/12/16 12:13:51 gdenisov Exp $                  
#

include ../include.mk

.PHONY: all pure train trainphd train.pure trainphd.pure

all:  train trainphd 

train:                $(RELDIR)/train
trainphd:             $(RELDIR)/trainphd
train.pure:           $(RELDIR)/train.pure
trainphd.pure:        $(RELDIR)/trainphd.pure

INCDIR      = ../compute_qv
CURDIR      = .
TTLIB       =  $(LIBDIR)/libtt.a
LIBS        = -lm
TRAINOBJS   =  $(OBJDIR)/train.o $(OBJDIR)/train_data.o\
               $(OBJDIR)/Btk_compute_match.o \
	       $(OBJDIR)/Btk_match_data.o $(OBJDIR)/Btk_sw.o
TRAINPHREDOBJS=$(OBJDIR)/trainphd.o $(OBJDIR)/train_data.o\
               $(OBJDIR)/Btk_compute_match.o \
	       $(OBJDIR)/Btk_match_data.o $(OBJDIR)/Btk_sw.o
CFLAGS     += -I$(INCDIR)
CFLAGS     += -I$(CURDIR)

$(RELDIR)/train: $(TRAINOBJS) $(TTLIB)
	@mkdir -p $(RELDIR)
	$(LINK.c) $(TRAINOBJS) $(TTLIB) $(LIBS) -o $@

$(RELDIR)/trainphd: $(TRAINPHREDOBJS) $(TTLIB)
	@mkdir -p $(RELDIR)
	$(LINK.c) $(TRAINPHREDOBJS) $(TTLIB) $(LIBS) -o $@

$(RELDIR)/train.pure: $(TRAINOBJS) $(TTLIB)
	@mkdir -p $(RELDIR)
	$(PURIFY) $(LINK.c) $(TRAINOBJS) $(TTLIB) $(LIBS) -o $@

$(RELDIR)/trainphd.pure: $(TRAINPHREDOBJS) $(TTLIB)
	@mkdir -p $(RELDIR)
	$(PURIFY) $(LINK.c) $(TRAINPHREDOBJS) $(TTLIB) $(LIBS) -o $@

$(OBJDIR)/Btk_match_data.o: Btk_match_data.c \
	Btk_match_data.h $(INCDIR)/Btk_qv.h $(INCDIR)/util.h 

$(OBJDIR)/Btk_sw.o: Btk_sw.c \
	Btk_sw.h $(INCDIR)/Btk_qv.h $(INCDIR)/util.h \
	$(INCDIR)/Btk_qv_data.h Btk_match_data.h \
	Btk_sw.h 

$(OBJDIR)/Btk_compute_match.o: Btk_compute_match.c \
	Btk_compute_match.h $(INCDIR)/Btk_qv.h $(INCDIR)/util.h \
	$(INCDIR)/Btk_qv_data.h Btk_match_data.h \
	Btk_sw.h

$(OBJDIR)/train.o:  train.c train_data.c \
        $(INCDIR)/ABI_Toolkit.h $(INCDIR)/FileHandler.h \
	$(INCDIR)/Btk_qv.h $(INCDIR)/Btk_qv_data.h $(INCDIR)/Btk_lookup_table.h\
	$(INCDIR)/Btk_compute_tpars.h $(INCDIR)/util.h \
	Btk_compute_match.h Btk_match_data.h $(INCDIR)/Btk_get_mixed_bases.h \
	Btk_sw.h 

$(OBJDIR)/trainphd.o:  trainphd.c train_data.c \
        $(INCDIR)/ABI_Toolkit.h $(INCDIR)/FileHandler.h \
	$(INCDIR)/Btk_qv.h $(INCDIR)/Btk_qv_data.h $(INCDIR)/Btk_lookup_table.h\
	$(INCDIR)/Btk_compute_tpars.h $(INCDIR)/util.h $(INCDIR)/Btk_qv_io.h\
	Btk_compute_match.h Btk_match_data.h \
	Btk_sw.h 

$(DIRS):
	mkdir -p $@

clean:
	@/bin/rm -f $(TRAINOBJS)
	@/bin/rm -f $(TRAINPHREDOBJS)
	@/bin/rm -f $(RELDIR)/train.pure
	@/bin/rm -f $(RELDIR)/trainphd.pure
	@/bin/rm -f $(RELDIR)/train
	@/bin/rm -f $(RELDIR)/trainphd
