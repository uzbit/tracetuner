#Top level Makefile for the tracetuner package com.paracel.tt
# Copyright (c) 2000 Paracel Inc.

BASE_DIR        = .
include $(BASE_DIR)/include.mk

SUBDIRS = io util chrom event view run
        
all :
	for dir in $(SUBDIRS); \
                do (cd $$dir; $(MAKE) $@); \
        done

	if [ -f ttuner_tools.jar ]; then \
		$(JAR) umf Manifest.txt ttuner_tools.jar com/ images/; \
	else \
		$(JAR) cmf Manifest.txt ttuner_tools.jar com/ images/; \
	fi

clean :
	rm -rf com
	rm -f ttuner_tools.jar
