## Makefile for SCF-C

N ?= 4

run:
	python3 faketwoel.py $(N)

twoel.trace: faketwoel.py ../tools/loops.py
	rm -f $@
	HL_TRACE_FILE=$@ python3 faketwoel.py

VIDEO_WIDTH=1920
VIDEO_HEIGHT=1080
video: twoel.trace
	rm -f twoel.mp4
	time cat twoel.trace | HalideTraceViz \
	  --size $(VIDEO_WIDTH) $(VIDEO_HEIGHT) --zoom 4 --timestep 2 --hold 100 --decay 5 5 --gray --auto_layout \
	| ffmpeg -loglevel warning -f rawvideo -pix_fmt bgr32 -s $(VIDEO_WIDTH)x$(VIDEO_HEIGHT) -i /dev/stdin -c:v h264 twoel.mp4

clean:
	rm -f *.trace
