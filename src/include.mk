#   Copyright (c) 1998, Paracel Inc. All rights reserved.
#
#   include.mk $Revision: 1.5 $

PFLAGS      = -messages=first -leaks-at-exit=yes -thread-stack-change=2000000
PURIFY      = $(PURE_PROG) $(PFLAGS)

OPTIMIZE = -g
ifdef RELEASE
OPTIMIZE = -O3
endif

ARCH :=$(shell uname -s)-$(shell uname -m)
BASEMACHINE = Linux

ifeq ($(ARCH),IRIX-IP22)
  BASEMACHINE = IRIX
endif
ifeq ($(ARCH),IRIX-IP32)
BASEMACHINE     = IRIX
endif
ifeq ($(ARCH),SunOS-sun4c)
BASEMACHINE     = Sun-Solaris
endif
ifeq ($(ARCH),SunOS-sun4m)
BASEMACHINE     = Sun-Solaris
endif
ifeq ($(ARCH),SunOS-sun4u)
BASEMACHINE     = Sun-Solaris
endif
ifeq ($(ARCH),SunOS-sun4c)
BASEMACHINE     = Sun-Solaris
endif
ifeq ($(ARCH),HP-UX-9000/715)
BASEMACHINE     = HP-UX
endif
ifeq ($(ARCH),SunOS-i86pc)
BASEMACHINE	= x86-Solaris
endif
ifeq ($(ARCH),Linux-i686)
BASEMACHINE     = Linux
endif
ifeq ($(ARCH),Linux-x86_64)
BASEMACHINE     = Linux_64
endif
ifeq ($(ARCH),OSF1-alpha)
BASEMACHINE	= OSF-Alpha
endif

ifeq ($(BASEMACHINE),IRIX)
CC              = /usr/local/gnu/bin/gcc
CFLAGS          = $(OPTIMIZE) -Wall -Wstrict-prototypes
AR              = /bin/ar
LD              = /bin/ld
DEPENDFLAGS     = -M
LDFLAGS         =
CP              = /sbin/cp
PURE_PROG       = $(PURE_HOME)/purify-4.1-irix6/purify
LOADLIBES       = -lm
CLIENTRANLIB    =
CP              = /sbin/cp
DYNAMIC         =
PURE_PROG       = /home/blackstone13/marc/junk/purify/purify-4.1-irix6/purify
endif

ifeq ($(BASEMACHINE),Sun-Solaris)
ifeq ($(PURE), 1)
	CC              = /usr/local/gnu/bin/gcc-2.7.2.2
	CFLAGS          = $(OPTIMIZE)  -Wall  -Wstrict-prototypes 
else
	CC              = /usr/local/gnu/bin/gcc
	CFLAGS          = $(OPTIMIZE) -Wall  -Wstrict-prototypes 
endif
AR              = /usr/ccs/bin/ar
LD              = /usr/ccs/bin/ld
CFLAGS_DYNAMIC  = $(OPTIMIZE) -Wall
DEFINES         = -DSUNOS_5_7 -D__EXTENSIONS__ -D_REENTRANT -D_POSIX_PTHREAD_SEMANTICS
LDFLAGS         =
LOADLIBES       = -lnsl -lm
CLIENTRANLIB    =
CP		= /bin/cp
DYNAMIC		= dynamic
BTK_MPI         = 1
LD_FLAGS_SHARED =  -dy -G -z text 
PURE_PROG        = /home/blackstone1/purify/pure/purify-4.2-solaris2/purify
endif

ifeq ($(BASEMACHINE),x86-Solaris)
ifeq ($(PURE), 1)
    CC              = /usr/local/gnu/bin/gcc-2.7.2.2
    CFLAGS          = $(OPTIMIZE)  -Wall  -Wstrict-prototypes
else
    CC              = /usr/local/gnu/bin/gcc
    CFLAGS          = $(OPTIMIZE) -Wall  -Wstrict-prototypes
endif
AR              = /usr/ccs/bin/ar
LD              = /usr/ccs/bin/ld
CFLAGS_DYNAMIC  = $(OPTIMIZE) -Wall
DEFINES         = -DSUNOS_5_6 -D__EXTENSIONS__ -D_REENTRANT -D_POSIX_PTHREAD_SEMANTICS
LDFLAGS         =
LOADLIBES       = -lnsl -lm
CLIENTRANLIB    =
CP      = /bin/cp
DYNAMIC     = dynamic
BTK_MPI         = 1
LD_FLAGS_SHARED =  -dy -G -z text
PURE_PROG        = /home/blackstone1/purify/pure/purify-4.2-solaris2/purify
endif

ifeq ($(BASEMACHINE),Linux)
CC              = /usr/bin/gcc
AR              = /usr/bin/ar
LD              = /usr/bin/ld
CFLAGS          = $(OPTIMIZE) -Wall
CFLAGS_DYNAMIC  = $(OPTIMIZE) -Wall
DEFINES         = -D_REENTRANT -D_POSIX_PTHREAD_SEMANTICS
LDFLAGS         =
LOADLIBES       =  -lm
CLIENTRANLIB    =
CP		= /bin/cp
DYNAMIC		= dynamic
LD_FLAGS_SHARED =  -dy -G -z text
endif

ifeq ($(BASEMACHINE),Linux_64)
CC              = /usr/bin/gcc
AR              = /usr/bin/ar
LD              = /usr/bin/ld
CFLAGS          = $(OPTIMIZE) -Wall
CFLAGS_DYNAMIC  = $(OPTIMIZE) -Wall
DEFINES         = -D_REENTRANT -D_POSIX_PTHREAD_SEMANTICS
LDFLAGS         =
LOADLIBES       =  -lm
CLIENTRANLIB    =
CP              = /bin/cp
DYNAMIC         = dynamic
LD_FLAGS_SHARED =  -dy -G -z text
endif

ifeq ($(BASEMACHINE),HP-UX)
PURE_PROG        = /home/blackstone13/marc/junk/purify/purify-4.2-hpux/purify
CC              = /usr/local/gnu/bin/gcc
AR              = ar
LD              =
CFLAGS          = $(OPTIMIZE) -Wall -Wstrict-prototypes
DEFINES         = -DHPPA_1_1_HPUX -DFDF_SIMAIO
LDFLAGS         =
LOADLIBES       = -lm
CLIENTRANLIB    =
CP		= /bin/cp
DYNAMIC		=
LD_FLAGS_SHARED =  -b
endif

ifeq ($(BASEMACHINE),OSF-Alpha)
CC              = /usr/local/bin/gcc
#CC              = cc -I/usr/include
MAKE            = /usr/local/bin/gmake
AR              = ar
LD              =
CFLAGS          = $(OPTIMIZE) -Wall -Wstrict-prototypes -taso
#DEFINES         = -pthread -D_REENTRANT -D_POSIX_PTHREAD_SEMANTICS 
DEFINES         = -D_XOPEN_SOURCE_EXTENDED -OSF1
LDFLAGS         =
LOADLIBES       = -lm
CLIENTRANLIB    =
CP		= /bin/cp
DYNAMIC		=
LD_FLAGS_SHARED =  -S -shared #Maybe -S isn't necessary
endif

LIBDIR := ../../lib/$(BASEMACHINE)
OBJDIR := ../../obj/$(BASEMACHINE)
RELDIR := ../../rel/$(BASEMACHINE)
ABIDIR := ../abitoolkit/$(BASEMACHINE)
DIRS   := $(LIBDIR) $(OBJDIR) $(RELDIR)

$(OBJDIR)/%.o: %.c
	mkdir -p $(OBJDIR)
	$(COMPILE.c) $< -o $@
