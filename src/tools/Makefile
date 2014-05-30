##	Copyright (c) 1999 Paracel Inc.  All rights reserved.
##
##
include ../include.mk

.PHONY: all pure visualize visualize.pure

all:  visualize
pure: visualize.pure

visualize:            $(RELDIR)/visualize
visualize.pure:       $(RELDIR)/visualize.pure


INCDIR      = ../compute_qv
JNIDIR      = $(JAVADIR)/include
CFLAGS     += -I$(INCDIR)
CFLAGS1     += -I$(JNIDIR)
CFLAGS2     += -I$(JNIDIR)/solaris
QVLIB       = $(LIBDIR)/libtt.a
VISLIB	=$(LIBDIR)/libvis.so
LIBS        = -lm
VISOBJS      = $(OBJDIR)/vis_main.o $(OBJDIR)/vis_process.o
VISLIBOBJS	=$(OBJDIR)/Vis_process.o 

visualize.lib:	$(VISLIB)

QVLIBOBJS   = $(OBJDIR)/Btk_atod.o $(OBJDIR)/Btk_compute_qv.o \
		$(OBJDIR)/Btk_qv_io.o $(OBJDIR)/Btk_compute_tpars.o \
		$(OBJDIR)/Btk_lookup_table.o \
		$(OBJDIR)/ABI_Toolkit.o $(OBJDIR)/ABI_Utils.o \
		$(OBJDIR)/FileHandler.o $(OBJDIR)/MemoryHandler.o

$(RELDIR)/visualize: $(VISOBJS) $(QVLIB)
	@mkdir -p $(RELDIR)
	$(LINK.c) $(VISOBJS) $(QVLIB) $(LIBS) -o $@
$(RELDIR)/visualize.pure: $(VISOBJS) $(QVLIB)
	@mkdir -p $(RELDIR)
	$(PURIFY) $(LINK.c) $(VISOBJS) $(QVLIB) $(LIBS) -o $@
      
$(VISLIB): $(VISLIBOBJS) $(QVLIBOBJS)
	@mkdir -p $(LIBDIR)
	$(CC)  -G -o $@ $(VISLIBOBJS) $(QVLIB)

$(OBJDIR)/vis_main.o: vis_main.c $(INCDIR)/ABI_Toolkit.h $(INCDIR)/util.h
	$(CC) $(CFLAGS) -c  vis_main.c -o $@

$(OBJDIR)/vis_process.o: vis_process.c $(INCDIR)/ABI_Toolkit.h $(INCDIR)/util.h
	$(CC)  $(CFLAGS)   -c  vis_process.c -o $@


$(OBJDIR)/Vis_process.o: Vis_process.c $(INCDIR)/ABI_Toolkit.h $(INCDIR)/util.h $(JNIDIR)/jni.h $(JNIDIR)/solaris/jni_md.h
	$(CC)  $(CFLAGS) $(CFLAGS1) $(CFLAGS2) -c  Vis_process.c -o $@


$(DIRS):
	mkdir -p $@

clean:
	@/bin/rm -f $(VISOBJS)
	@/bin/rm -f visualize.pure
	@/bin/rm -f visualize
	@/bin/rm -f $(VISLIB)




