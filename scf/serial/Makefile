## Makefile for SCF-C

CC  = clang++
#CFLAGS = -O3 -mtune=native -march=native -Wall -Wno-deprecated -DNO_PRECALC
#CFLAGS = -O3 -mtune=native -march=native -Wall -Wno-deprecated -DNO_PRECALC -DNEGLIGIBILITY_CHECKS
#CFLAGS = -O3 -mtune=native -march=native -Wall -Wno-deprecated
#CFLAGS = -O3 -mtune=native -march=native -Wall -Wno-deprecated -DNO_PRECALC -DSYMMETRY
CFLAGS = -O3 -mtune=native -march=native -Wall -Wno-deprecated -DNEGLIGIBILITY_CHECKS -DSYMMETRY
LDFLAGS =

OBJFILES = input.o integ.o output.o timer.o scf.o diagonalize.o twoel.o

MTARGET = scf.x

all: $(MTARGET)

scf.x:  $(OBJFILES)
	$(CC) $(OBJFILES) $(CFLAGS) $(LDFLAGS) -o $(MTARGET) -lm

run: scf.x
	./scf.x

vg: scf.x
	valgrind ./scf.x

cg: scf.x
	valgrind --tool=callgrind ./scf.x

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *~ *.o $(MTARGET)
