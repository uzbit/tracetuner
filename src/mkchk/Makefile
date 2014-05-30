#
# Copyright (c) 1999-2003 Paracel, Inc.  All rights reserved.
#
# $Id: Makefile,v 1.8 2009/01/05 22:45:12 gdenisov Exp $                   
#


include ../include.mk

.PHONY: all checkqv checkbc checkbcphd checkqvphd
 
all:  checkbc checkqv checkbcphd checkqvphd 

checkqv:	$(RELDIR)/checkqv
checkbc:	$(RELDIR)/checkbc
checkbcphd:	$(RELDIR)/checkbcphd
checkqvphd:	$(RELDIR)/checkqvphd

INCTOOLSDIR = ../../..
INCJTCTTDIR = $(INCTOOLSDIR)/jtc_tt
INCDIR      = ../compute_qv
INCTRAINDIR = ../mktrain
INCLUTDIR   = ../mklut
CURDIR      = .
LIBS        = -lm -lz

CHECKQVOBJS	= $(OBJDIR)/checkqv.o $(OBJDIR)/check_data.o
CHECKBCOBJS	= $(OBJDIR)/checkbc.o $(OBJDIR)/check_data.o
CHECKBCPHDOBJS	= $(OBJDIR)/checkbcphd.o $(OBJDIR)/check_data.o
CHECKQVPHDOBJS		= $(OBJDIR)/checkqvphd.o $(OBJDIR)/check_data.o
TTLIB		= $(LIBDIR)/libtt.a
CFLAGS     += -I$(INCDIR)
CFLAGS     += -I$(INCTRAINDIR)
CFLAGS     += -I$(INCLUTDIR)
CFLAGS     += -I$(CURDIR)
CFLAGS     += -I$(INCIOLIBDIR) -I$(INCJTCTTDIR) -I$(CURDIR)
CFLAGS     += -I$(INCTOOLSDIR)/io_lib
CFLAGS     += -I$(INCTOOLSDIR)/io_lib/utils

LINK.c = g++

$(RELDIR)/checkqv: $(CHECKQVOBJS) $(TTLIB)
	@mkdir -p $(RELDIR)
	$(LINK.c) $(CHECKQVOBJS) $(TTLIB) $(LIBS) -o $@ $(LIBS) $(TTLIB) $(IOLIB)

$(RELDIR)/checkbc: $(CHECKBCOBJS) $(TTLIB)
	@mkdir -p $(RELDIR)
	$(LINK.c) $(CHECKBCOBJS) $(TTLIB) $(LIBS) -o $@ $(LIBS) $(TTLIB) $(IOLIB)

$(RELDIR)/checkbcphd: $(CHECKBCPHDOBJS) $(TTLIB)
	@mkdir -p $(RELDIR)
	$(LINK.c) $(CHECKBCPHDOBJS) $(TTLIB) $(LIBS) -o $@ $(LIBS) $(TTLIB) $(IOLIB)

$(RELDIR)/checkqvphd: $(CHECKQVPHDOBJS) $(TTLIB)
	@mkdir -p $(RELDIR)
	$(LINK.c) $(CHECKQVPHDOBJS) $(TTLIB) $(LIBS) -o $@ $(LIBS) $(TTLIB) $(IOLIB)

$(OBJDIR)/check_data.o: check_data.c $(INCLUTDIR)/lut.h\
                        $(INCDIR)/Btk_atod.h $(INCDIR)/Btk_qv.h \
                        $(INCDIR)/Btk_lookup_table.h \
                        check_data.h

$(OBJDIR)/checkqv.o:	checkqv.c $(INCLUTDIR)/lut.h\
			$(INCDIR)/Btk_atod.h $(INCDIR)/Btk_qv.h\
			$(INCDIR)/Btk_lookup_table.h $(INCTRAINDIR)/train.h \
                        check_data.h

$(OBJDIR)/checkbc.o:	checkbc.c $(INCLUTDIR)/lut.h\
			$(INCDIR)/Btk_atod.h $(INCDIR)/Btk_qv.h \
			$(INCDIR)/Btk_lookup_table.h check_data.h

$(OBJDIR)/checkqvphd.o:	checkqvphd.c $(INCLUTDIR)/lut.h\
			$(INCDIR)/Btk_atod.h $(INCDIR)/Btk_qv.h \
			$(INCDIR)/Btk_lookup_table.h

$(OBJDIR)/checkbcphd.o: checkbcphd.c \
			$(INCDIR)/Btk_atod.h $(INCDIR)/Btk_qv.h \
			$(INCDIR)/Btk_lookup_table.h  check_data.h

$(DIRS):
	mkdir -p $@

clean:
	@/bin/rm -f $(CHECKQVOBJS) $(CHECKBCOBJS) $(CHECKBCPHDOBJS) $(CHECKQVPHDOBJS)
	@/bin/rm -f $(RELDIR)/checkbc $(RELDIR)/checkqv 
	@/bin/rm -f $(RELDIR)/checkbcphd $(RELDIR)/checkqvphd 
