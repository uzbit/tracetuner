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
# $Id: Makefile,v 1.11 2009/01/01 16:21:50 gdenisov Exp $
#
#  Makefile for main TraceTuner library and executable
#

include ../include.mk

OSNAME:="$(shell uname -s)"


.PHONY: all pure ttuner example ttuner.pure showfile


all:  ttuner 
pure: ttuner.pure

ttuner:          $(RELDIR)/ttuner
ttuner.pure:     $(RELDIR)/ttuner.pure
example:         $(RELDIR)/example
showfile:        $(RELDIR)/showfile

INCDIR      = ../mktrain
CURDIR      = .
QVLIB       = $(LIBDIR)/libtt.a
LIBS        = -lm
QVOBJS      = $(OBJDIR)/main.o
QVLIBSRCS   = $(OBJDIR)/Btk_match_data.c $(OBJDIR)/Btk_compute_match.c \
	      $(OBJDIR)/Btk_sw.c $(OBJDIR)/Btk_process_indels.c        \
	      $(OBJDIR)/Btk_atod.c $(OBJDIR)/Btk_compute_qv.c          \
              $(OBJDIR)/Btk_call_bases.c                               \
              $(OBJDIR)/Btk_process_raw_data.c                          \
	      $(OBJDIR)/util.c  $(OBJDIR)/nr.c			       \
              $(OBJDIR)/Btk_process_peaks.c $(OBJDIR)/Btk_qv_io.c      \
              $(OBJDIR)/Btk_compute_tp.c $(OBJDIR)/Btk_compute_tpars.c \
              $(OBJDIR)/Btk_get_mixed_bases.c                           \
              $(OBJDIR)/Btk_lookup_table.c $(OBJDIR)/Btk_qv_funs.c     \
              $(OBJDIR)/ABI_Toolkit.c $(OBJDIR)/SFF_Toolkit.c          \
              $(OBJDIR)/Btk_default_table.c                            \
              $(OBJDIR)/FileHandler.c $(OBJDIR)/SCF_Toolkit.c          \
              $(OBJDIR)/context_table.c                                \
              $(OBJDIR)/tracepoly.c 				

QVLIBOBJS  = $(patsubst %.c,%.o,$(QVLIBSRCS))
EXAMPLEOBJS = $(OBJDIR)/example.o
SHOWFILEOBJS = $(OBJDIR)/showfile.o
CFLAGS	   += -I$(INCDIR) -DOS_NAME='$(OSNAME)'
CFLAGS	   += -I$(CURDIR)

$(RELDIR)/ttuner: $(QVOBJS) $(QVLIB)
	@mkdir -p $(RELDIR)
	$(LINK.c) $(QVOBJS) $(QVLIB) $(LIBS) -o $@

$(RELDIR)/example: $(EXAMPLEOBJS) $(QVLIB)
	@mkdir -p $(RELDIR)
	$(LINK.c) $(EXAMPLEOBJS) $(QVLIB) $(LIBS) -o $@

$(RELDIR)/ttuner.pure: $(QVOBJS) $(QVLIB)
	@mkdir -p $(RELDIR)
	$(PURIFY) $(LINK.c) $(QVOBJS) $(QVLIB) $(LIBS) -o $@

$(RELDIR)/showfile: $(SHOWFILEOBJS) $(QVLIB)
	@mkdir -p $(RELDIR)
	$(LINK.c) $(SHOWFILEOBJS) $(QVLIB) $(LIBS) -o $@

$(QVLIB): $(QVLIBOBJS)
	@mkdir -p $(LIBDIR)
	$(AR) rc $@ $(QVLIBOBJS)

$(DIRS):
	mkdir -p $@

clean:
	@/bin/rm -f $(EXAMPLEOBJS) $(QVOBJS) $(CHECKORDEROBJS)
	@/bin/rm -f $(QVLIBOBJS) $(QVLIB)
	@/bin/rm -f $(RELDIR)/qvdata.pure 
	@/bin/rm -f $(RELDIR)/ttuner.pure
	@/bin/rm -f $(RELDIR)/qvdata  $(RELDIR)/ttuner 
	@/bin/rm -f $(RELDIR)/get_default_lut.perl
	@/bin/rm -f $(RELDIR)/get_context.perl
	@/bin/rm -f $(RELDIR)/lut_to_default_lut.perl         
        

#depend:
#	$(CC) $(INCLUDE) -M $(SOURCES) > depend.mk

# The rule to make Btk_match_data.o, Btk_sw.o, Btk_compute_match.o
# It is not the same as the rule for the other object files, because
# the source code of these 3 object files are located in a different
# directory.
$(OBJDIR)/%.o: $(INCDIR)/%.c 
	mkdir -p $(OBJDIR)
	$(COMPILE.c) $< -o $@
$(OBJDIR)/Btk_match_data.o: Btk_qv.h util.h $(INCDIR)/Btk_match_data.h
$(OBJDIR)/Btk_sw.o: Btk_qv.h util.h Btk_qv_data.h
$(OBJDIR)/Btk_sw.o: $(INCDIR)/Btk_sw.h $(INCDIR)/Btk_match_data.h
$(OBJDIR)/Btk_compute_match.o: Btk_qv.h util.h Btk_qv_data.h
$(OBJDIR)/Btk_compute_match.o: $(INCDIR)/Btk_sw.h $(INCDIR)/Btk_match_data.h
$(OBJDIR)/Btk_compute_match.o: $(INCDIR)/Btk_compute_match.h

# The rule to make the other object files whose source code are located
# in the current directory.
$(OBJDIR)/%.o: %.c
	mkdir -p $(OBJDIR)
	$(COMPILE.c) $< -o $@
$(OBJDIR)/example.o: example.c Btk_lookup_table.h Btk_qv_io.h
$(OBJDIR)/showfile.o: showfile.c  Btk_qv_io.h
$(OBJDIR)/ABI_Toolkit.o: ABI_Toolkit.h
$(OBJDIR)/SCF_Toolkit.o: ABI_Toolkit.h SCF_Toolkit.h 
$(OBJDIR)/Btk_call_bases.o: Btk_qv.h util.h Btk_qv_data.h
$(OBJDIR)/Btk_call_bases.o: Btk_qv_funs.h Btk_process_peaks.h Btk_call_bases.h
$(OBJDIR)/Btk_call_bases.o: context_table.h Btk_lookup_table.h
$(OBJDIR)/Btk_call_bases.o: tracepoly.h
$(OBJDIR)/Btk_call_bases.o: Btk_qv_data.h
$(OBJDIR)/Btk_process_raw_data.o: Btk_process_raw_data.h Btk_qv.h
$(OBJDIR)/Btk_compute_qv.o: Btk_qv.h util.h Btk_lookup_table.h
$(OBJDIR)/Btk_compute_qv.o: Btk_qv_data.h Btk_compute_qv.h Btk_compute_tpars.h
$(OBJDIR)/Btk_compute_qv.o: Btk_atod.h Btk_default_table.h
$(OBJDIR)/Btk_compute_tp.o: Btk_qv.h util.h Btk_qv_data.h
$(OBJDIR)/Btk_compute_tp.o: Btk_lookup_table.h Btk_compute_qv.h 
$(OBJDIR)/Btk_compute_tp.o: Btk_compute_tp.h
$(OBJDIR)/Btk_compute_tp.o: Btk_qv_funs.h Btk_process_peaks.h Btk_call_bases.h
$(OBJDIR)/Btk_compute_tpars.o: Btk_qv.h util.h Btk_lookup_table.h
$(OBJDIR)/Btk_compute_tpars.o: Btk_qv_data.h Btk_process_peaks.h 
$(OBJDIR)/Btk_compute_tpars.o: Btk_call_bases.h
$(OBJDIR)/Btk_compute_tpars.o: Btk_compute_tp.h Btk_atod.h Btk_compute_tpars.h
$(OBJDIR)/Btk_compute_tpars.o: Btk_qv_data.h
$(OBJDIR)/Btk_default_table.o: Btk_lookup_table.h
$(OBJDIR)/Btk_lookup_table.o: Btk_qv.h Btk_lookup_table.h
$(OBJDIR)/Btk_lookup_table.o: Btk_atod.h
$(OBJDIR)/context_table.o: context_table.h
$(OBJDIR)/tracepoly.o: tracepoly.h
$(OBJDIR)/Btk_process_peaks.o: Btk_qv_funs.h Btk_process_peaks.h
$(OBJDIR)/Btk_process_peaks.o: Btk_qv_data.h
$(OBJDIR)/Btk_qv_io.o: FileHandler.h Btk_qv.h util.h Btk_qv_io.h 
$(OBJDIR)/Btk_qv_io.o: ABI_Toolkit.h SCF_Toolkit.h
$(OBJDIR)/Btk_qv_io.o: $(INCDIR)/Btk_match_data.h
$(OBJDIR)/Btk_qv_io.o: $(INCDIR)/Btk_compute_match.h
$(OBJDIR)/Btk_qv_io.o: Btk_qv_data.h
$(OBJDIR)/FileHandler.o: ABI_Toolkit.h SCF_Toolkit.h
$(OBJDIR)/FileHandler.o: FileHandler.h
$(OBJDIR)/Btk_qv_funs.o: Btk_qv_funs.h 
$(OBJDIR)/Btk_qv_funs.o: Btk_qv_data.h 
$(OBJDIR)/SFF_Toolkit.o: SFF_Toolkit.h
$(OBJDIR)/main.o: ABI_Toolkit.h FileHandler.h Btk_qv.h util.h Btk_qv_data.h
$(OBJDIR)/main.o: Btk_lookup_table.h Btk_compute_qv.h Btk_qv_io.h
