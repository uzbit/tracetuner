# Top level Makefile for TraceTuner source code
# Copyright (c) 2008 J. Craig Venter Institute

SUBDIRS = compute_qv mktrain mklut mkchk 

all:
	for dir in $(SUBDIRS); \
            do (cd $$dir; $(MAKE) $@); \
        done
	
clean:
	cd compute_qv; $(MAKE) clean
	cd mktrain; $(MAKE) clean
	cd mklut; $(MAKE) clean
	cd mkchk; $(MAKE) clean
