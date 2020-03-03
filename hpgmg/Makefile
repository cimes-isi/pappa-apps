run: original/build/bin/hpgmg
	OMP_NUM_THREADS=32 time original/build/bin/hpgmg-fv 8 9

original/build/Makefile:
	cd original && ./configure --CFLAGS=-fopenmp

original/build/bin/hpgmg: original/build/Makefile
	make -C original/build

clean:
	rm -rf original/build
