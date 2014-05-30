#
# Copyright (c) 1999-2003 Paracel Inc.  All rights reserved.
#
# example.mk $Revision: 1.2 $

CC          = gcc
QVLIB       = libtt.a
EXAMPLEOBJS = example.o
LIBS        = -lm

example: $(EXAMPLEOBJS) $(QVLIB)
	$(LINK.c) $(EXAMPLEOBJS) $(QVLIB) $(LIBS) -o $@

example.o: example.c Btk_qv.h Btk_compute_qv.h Btk_lookup_table.h Btk_qv_io.h

clean:
	@/bin/rm -f example example.o
