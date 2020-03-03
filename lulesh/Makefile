run: original/build/lulesh2.0
	cd original/build && OMP_NUM_THREADS=32 ./lulesh2.0 -s 40

original/build/lulesh2.0: original/build/Makefile
	make -C original/build

original/build/Makefile:
	mkdir original/build
	cd original/build && cmake ..

clean:
	rm -rf original/build
