# Intro

SCF

* Self-Consistent Field
* A.k.a. The Hartree-Fock method
* Determines the wave function and energy of a quantum system
* Part of NWChem, a computational quantum chemistry app from PNNL

# Implementations

There are several implementations here; each one lives in its own folder.

## original
This is a simple serial C implentation of the SCF module, for use as a starting
point.

## serial
This is hand-optimized serial C code.  Individual optimizations may be enabled
by modifying CFLAGS in the Makefile.

## openmp
This is based on the above hand-optimized C code, with openmp parallelism.

## halide
This uses a Halide-generated **twoel** function.  Looks for a halide build at
`../../halide/build/`.
The autoscheduler can be chosen with a Makefile variable.  If the Li2018
scheduler is selected, the
corresponding library must be built (manually).


# Basis for comparison

To keep comparisons clean, the problem size, compiler version and CFLAGS should
be kept consistent for all builds.

## Compiler & CFLAGS

The code is built with clang++ version 9, and g++ version 8.3, with parameters
`-O3 -march=native -mtune=native`.

## Baseline

This implementation of SCF is a simplified form of the SCF feature in
[nwchem](https://github.com/nwchemgit/nwchem/).  This provides our basis for
comparison, which is threaded, but no symmetry or precomputation.

## Input data

Each implementation has a data file included with it.  These datasets are all
the same, and consist of 30 Be (Beryllium) atoms, in a row, 4 units apart from
each other.

# Running
You can go into the appropriate implementation folder and run:

    make run

to build the program (if necessary) and run it.

SCF takes an input file called `be.inpt`.  This input file contains the list of
atoms to simulate; the number of atoms determines the problem size.  A different
problem size can be selected quickly by running:

    `./be.sh <num>`

and specifying the number of atoms to generate.

When it runs, it will print some initial parameters, then execute a convergence
algorithm.  This algorithm usually takes 8 to 15 iterations to complete.  Once
it completes, it will print the final energy value, the list of eigenvalues,
some execution timings and stats.

Here is an example of what the execution timings look like:

          init       onel      twoel       diag       dens       print       ncpu
          ----       ----      -----       ----       ----       ----        ----
          0.00       0.09      73.68       0.75       0.02       0.00
     elapsed time in seconds      74.54

In this example, the complete SCF module took 74.54 seconds to run, and 73.68 of
those were spent in the `twoel()` function.

# Reference stuff

* [Hartree-Fock algorithm @ wikipedia](https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method#Hartree%E2%80%93Fock_algorithm)
* [nwchem docs for SCF input directive](https://github.com/nwchemgit/nwchem/wiki/Hartree-Fock-Theory-for-Molecules)
