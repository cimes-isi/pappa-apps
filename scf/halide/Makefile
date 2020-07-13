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

VIDEO_WIDTH=1920
VIDEO_HEIGHT=1080
video: scf.trace
	rm -f scf.mp4
	cat scf.trace | HalideTraceViz \
	  --size $(VIDEO_WIDTH) $(VIDEO_HEIGHT) --zoom 4 --timestep 32 --hold 100 --decay 5 5 --gray --auto_layout \
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
