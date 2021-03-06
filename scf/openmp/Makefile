## Makefile for SCF-C


CC = gcc
CFLAGS += -g -O3 -mtune=native -march=native -Wall -fopenmp
LDFLAGS += -lm 

OBJFILES = input.o integ.o output.o timer.o scf.o diagonalize.o twoel.o

run: scf.x
	time ./scf.x

vg: scf.x
	valgrind ./scf.x

MTARGET = scf.x

ifdef BLAS
CFLAGS += -DBLAS
LDFLAGS += -lblas
endif

ifeq (x$(SWARM_PATH),x)
else
CFLAGS += -I$(SWARM_PATH)/include
LDFLAGS += $(SWARM_PATH)/lib/libswarm-tracing.a 
endif

ifdef NO_PRECALC
CFLAGS += -DNO_PRECALC
endif

ifdef USE_CACHE
CFLAGS += -DUSE_CACHE
endif

ifndef NO_BOOKKEEPING
CFLAGS += -DBOOKKEEPING
endif

ifndef NO_SYMMETRY
CFLAGS += -DSYMMETRY
endif

ifdef BAD_CACHE
CFLAGS += -DBAD_CACHE_UTILIZATION
endif

ifdef LAPACKE
CFLAGS += -DUSE_LAPACK
LDFLAGS += -llapacke -llapack -lblas
endif

ifeq (x$(MKLROOT),x)
else
CFLAGS += -I$(MKLROOT)/include -DUSE_MKL -DUSE_LAPACK -DBLAS
LDFLAGS += -L$(MKLROOT)/lib/intel64 -L$(MKLROOT)/../compiler/lib/intel64 -Wl,--start-group -lmkl_intel_lp64  -lmkl_intel_thread -lmkl_core -Wl,--end-group -liomp5
endif

ifndef OMP_VERSION
CFLAGS += -DOMP_VERSION=3
endif

all: $(MTARGET)

scf.x:  $(OBJFILES)
	$(CC) $(OBJFILES) $(CFLAGS) -o $(MTARGET) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *~ *.o $(MTARGET)
