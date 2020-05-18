## Makefile for SCF-C

CC  = clang++
CXX = clang++
HALIDE_INCDIR=../../halide/build/include
HALIDE_LIBDIR=../../halide/build/lib
CFLAGS = -O3 -mtune=native -march=native -g -Wall -Wno-deprecated -I$(HALIDE_INCDIR)
CXXFLAGS = $(CFLAGS)
CLINK = -L$(HALIDE_LIBDIR) -ldl -lpthread
CXXLINK = -L$(HALIDE_LIBDIR) -lHalide

TARGET=host-x86-64
#TARGET=$(TARGET)-profile
MACHINE_PARAMS=32,16777216,40

# to use manually defined scheduler in twoel_gen.cpp, leave AUTOSCHEDULER commented.
#AUTOSCHEDULER=Mullapudi2016
#AUTOSCHEDULER=Li2018
#AUTOSCHEDULER=Adams2019

ifneq ($(AUTOSCHEDULER),)
SCHEDLIB=../../halide/build/lib/libauto_schedule.so
ifeq ($(AUTOSCHEDULER),Li2018)
SCHEDLIB=../../halide/build/lib/libgradient_autoscheduler.so
endif
endif

OBJFILES = input.o integ.o output.o timer.o scf.o diagonalize.o twoel.a

MTARGET = scf.x

all: $(MTARGET)

scf.x:  $(OBJFILES)
	$(CC) $(OBJFILES) $(CFLAGS) $(CLINK) -o $(MTARGET) -lm

run: scf.x
	time ./scf.x

trace scf.trace: scf.x be.inpt
	rm -f scf.trace
	HL_TRACE_FILE=scf.trace ./scf.x

pictures: scf.trace
	../../halide/build/bin/HalideTraceDump -i scf.trace -t png

VIDEO_WIDTH=600
VIDEO_HEIGHT=480
video: scf.trace
	rm -f scf.mp4
	cat scf.trace | ../../halide/build/bin/HalideTraceViz \
	  --size $(VIDEO_WIDTH) $(VIDEO_HEIGHT) --zoom 4 --timestep 32 --hold 100 --decay 5 5 --gray \
	  --move  32  48 --func g_fock_out               --move  32  40 --label g_fock_out               "g_fock out"     1 \
	  --move 232  48 --func g_fock_no_symmetry       --move 232  40 --label g_fock_no_symmetry       "fock nosymm"    1 \
	  --move 432  48 --func g_fock_pairwise_symmetry --move 432  40 --label g_fock_pairwise_symmetry "fock pairwise"  1 \
	  --move  32 148 --func g_fock_single_symmetry   --move  32 140 --label g_fock_single_symmetry   "fock single"    1 \
	  --move 232 148 --func g_fock_double_symmetry   --move 232 140 --label g_fock_double_symmetry   "fock double"    1 \
	  --move 432 148 --func g_fock_triple_symmetry   --move 432 140 --label g_fock_triple_symmetry   "fock triple"    1 \
	  --move  32 248 --func g_dens_in                --move  32 240 --label g_dens_in                g_dens           1 \
	  --move 232 248 --func 'r2$$0'                  --move 232 240 --label 'r2$$0'                  r2               1 \
	  --move 432 248 --func fac2                     --move 432 240 --label fac2                     fac2             1 \
	  --move  32 348 --func 'x2$$0'                  --move  32 340 --label 'x2$$0'                  x2               1 \
	  --move 232 348 --func 'y2$$0'                  --move 232 340 --label 'y2$$0'                  y2               1 \
	  --move 432 348 --func 'z2$$0'                  --move 432 340 --label 'z2$$0'                  z2               1 \
	  --move 432 450 --func rnorm_in                 --move 432 442 --label rnorm_in                 rnorm            1 \
	  --zoom 1.6 \
	  --move  32 450 --func fm_in                    --move  32 442 --label fm_in                    fm               1 \
	| ffmpeg -loglevel warning -f rawvideo -pix_fmt bgr32 -s $(VIDEO_WIDTH)x$(VIDEO_HEIGHT) -i /dev/stdin -c:v h264 scf.mp4

vg: scf.x
	valgrind ./scf.x

cg: scf.x
	valgrind --tool=callgrind ./scf.x

twoel.a: tools/loops.py
	time python3 tools/loops.py

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

scf.o: twoel.a

clean:
	rm -f *~ *.o $(MTARGET) *.a *_gen *.stmt *.halide_generated.cpp twoel.h twoel.schedule.h twoel.registration.cpp *.trace *.png *.mp4
