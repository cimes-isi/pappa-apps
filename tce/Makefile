# Copyright (C) 2013 ET International, Inc.

ifeq (x$(MKLROOT),x)
BLASCFLAGS=-DUSE_OPENBLAS
BLASLIBS=-lblas
#BLASCFLAGS=-DUSE_ATLAS_BLAS
#BLASLDFLAGS=-L/usr/lib/atlas-base -Wl,-rpath,/usr/lib/atlas-base
#BLASLIBS=-lcblas
else
BLASCFLAGS=-DUSE_MKL_BLAS -I$(MKLROOT)/include
BLASFFLAGS=-I$(MKLROOT)/include
BLASLDFLAGS=-L$(MKLROOT)/lib/intel64 -Wl,-rpath,$(MKLROOT)/lib/intel64
BLASLIBS=-lmkl_intel
endif

CC=gcc
CFLAGS=-g -O3 -Wall -march=native -mtune=native $(BLASCFLAGS) $(USERCFLAGS)
LDFLAGS=$(BLASLDFLAGS) $(USERLDFLAGS)
LIBS=$(BLASLIBS) $(USERLIBS)

CC2_O3_URL=https://www.extreme-scale.com/downloads/tce/test_cc2_o3.tar.lzma
CCSD_O3_URL=https://www.extreme-scale.com/downloads/tce/test_ccsd_o3.tar.lzma
CCSD_C6H6_2EORB_URL=https://www.extreme-scale.com/downloads/tce/test_ccsd_c6h6_2eorb.tar.lzma
CC2_O3_DATAPATH=data/test_cc2_o3
CCSD_O3_DATAPATH=data/test_ccsd_o3
CCSD_C6H6_2EORB_DATAPATH=data/test_ccsd_c6h6_2eorb

argslist = $(1)/params.i1.dump $(1)/f1.i1.before.dump $(1)/v2.i1.before.dump $(1)/t1.i1.before.dump $(1)/t2.i1.before.dump $(1)/r1.i1.after.dump $(1)/r2.i1.after.dump

all: cc2_t1_t2_standalone ccsd_t1_t2_standalone $(CC2_O3_DATAPATH)/README $(CCSD_O3_DATAPATH)/README $(CCSD_C6H6_2EORB_DATAPATH)/README


run: run_cc2_O3 run_ccsd_O3 run_ccsd_C6H6_2eorb

run_cc2_O3: cc2_t1_t2_standalone $(CC2_O3_DATAPATH)/README
	time ./$< $(call argslist,$(CC2_O3_DATAPATH))

run_ccsd_O3: ccsd_t1_t2_standalone $(CCSD_O3_DATAPATH)/README
	time ./$< $(call argslist,$(CCSD_O3_DATAPATH))

run_ccsd_C6H6_2eorb: ccsd_t1_t2_standalone $(CCSD_C6H6_2EORB_DATAPATH)/README
	time ./$< $(call argslist,$(CCSD_C6H6_2EORB_DATAPATH))

cc2_t1_t2_standalone: cc2_t1.o cc2_t2.o cc2_t1_t2_standalone.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

ccsd_t1_t2_standalone: ccsd_t1.o ccsd_t2.o ccsd_t1_t2_standalone.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

cc2_t1.tex cc2_t1.F cc2_t1.f90 cc2_t1.c: cc2_t1.tt driver.py tce.py
	python driver.py cc2_t1 cc2_t1.tt >/dev/null

cc2_t2.tex cc2_t2.F cc2_t2.f90 cc2_t2.c: cc2_t2.tt driver.py tce.py
	python driver.py cc2_t2 cc2_t2.tt >/dev/null

ccsd_t1.tex ccsd_t1.F ccsd_t1.f90 ccsd_t1.c: ccsd_t1.tt driver.py tce.py
	python driver.py ccsd_t1 ccsd_t1.tt >/dev/null

ccsd_t2.tex ccsd_t2.F ccsd_t2.f90 ccsd_t2.c: ccsd_t2.tt driver.py tce.py
	python driver.py ccsd_t2 ccsd_t2.tt >/dev/null

%.pdf: %.tex Makefile
	(echo '\\documentclass{article}'; echo '\\begin{document}'; cat $<; echo '\\end{document}') > $<.tex
	pdflatex $<.tex
	mv $<.pdf $@

%.o: %.c *.h Makefile
	$(CC) $(CFLAGS) -c -o $@ $<

cc2_t1_t2_standalone.o: generic_t1_t2_standalone.c Makefile
	$(CC) $(CFLAGS) -DMODEL=cc2 -c -o $@ $<

ccsd_t1_t2_standalone.o: generic_t1_t2_standalone.c Makefile
	$(CC) $(CFLAGS) -DMODEL=ccsd -c -o $@ $<

$(CC2_O3_DATAPATH)/README: data/test_cc2_o3.tar.lzma
	cd data && tar --lzma -xvf test_cc2_o3.tar.lzma
	touch $@

$(CCSD_O3_DATAPATH)/README: data/test_ccsd_o3.tar.lzma
	cd data && tar --lzma -xvf test_ccsd_o3.tar.lzma
	touch $@

$(CCSD_C6H6_2EORB_DATAPATH)/README: data/test_ccsd_c6h6_2eorb.tar.lzma
	cd data && tar --lzma -xvf test_ccsd_c6h6_2eorb.tar.lzma
	touch $@

data/test_cc2_o3.tar.lzma:
	mkdir -p data && cd data && wget $(CC2_O3_URL) -O tmp && mv tmp test_cc2_o3.tar.lzma

data/test_ccsd_o3.tar.lzma:
	mkdir -p data && cd data && wget $(CCSD_O3_URL) -O tmp && mv tmp test_ccsd_o3.tar.lzma

data/test_ccsd_c6h6_2eorb.tar.lzma:
	mkdir -p data && cd data && wget $(CCSD_C6H6_2EORB_URL) -O tmp && mv tmp test_ccsd_c6h6_2eorb.tar.lzma

clean:
	for EXT in c F f90 tex tex.tex tex.aux tex.log pdf py; do rm -f cc2_t[012].$$EXT ccsd_t[012].$$EXT; done
	rm -f *.py~ *.pyc *.o tmp.* cc2_t1_t2_standalone ccsd_t1_t2_standalone
