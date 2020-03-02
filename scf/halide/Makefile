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
AUTOSCHEDULER=Mullapudi2016
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

VIDEO_WIDTH=500
VIDEO_HEIGHT=400
video: scf.trace
	rm -f scf.mp4
	cat scf.trace | ../../halide/build/bin/HalideTraceViz \
	  --size $(VIDEO_WIDTH) $(VIDEO_HEIGHT) --zoom 4 --timestep 8 --hold 100 --decay 5 5 --gray \
	  --move  32  48 --func g_fock_out --move  32  40 --label g_fock_out "g_fock out" 1 --move 240  48 --func g_dens     --move 240  40 --label g_dens      g_dens      1 \
	  --move  32 248 --func x          --move  32 240 --label x           x           1 --move 240 248 --func expnt      --move 240 240 --label expnt       expnt       1 \
	  --move  32 292 --func y          --move  32 284 --label y           y           1 --move 240 292 --func rnorm      --move 240 284 --label rnorm       rnorm       1 \
	  --move  32 340 --func z          --move  32 332 --label z           z           1 \
      --zoom 2 \
      --move  32 200 --func fm         --move  32 192 --label fm          fm          1 \
	| ffmpeg -f rawvideo -pix_fmt bgr32 -s $(VIDEO_WIDTH)x$(VIDEO_HEIGHT) -i /dev/stdin -c:v h264 scf.mp4

vg: scf.x
	valgrind ./scf.x

cg: scf.x
	valgrind --tool=callgrind ./scf.x

twoel_gen: twoel_gen.cpp
	$(CXX) $(CXXFLAGS) -o $@ ../../halide/tools/GenGen.cpp $^ $(CXXLINK)

twoel.a: twoel_gen Makefile
ifeq ($(AUTOSCHEDULER),)
	LD_LIBRARY_PATH=$(HALIDE_LIBDIR) ./twoel_gen -g twoel -o . -e static_library,h,schedule,cpp,stmt,registration target=$(TARGET) machine_params=$(MACHINE_PARAMS)
else
	LD_LIBRARY_PATH=$(HALIDE_LIBDIR) ./twoel_gen -g twoel -o . -e static_library,h,schedule,cpp,stmt,registration -p $(SCHEDLIB) -s $(AUTOSCHEDULER) target=$(TARGET) auto_schedule=true machine_params=$(MACHINE_PARAMS)
endif

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

scf.o: twoel.a

clean:
	rm -f *~ *.o $(MTARGET) *.a *_gen *.stmt *.halide_generated.cpp twoel.h twoel.schedule.h twoel.registration.cpp *.trace *.png *.mp4
