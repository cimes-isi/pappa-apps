# Intro

TCE

* Tensor Contraction Engine
* Implements multiple iterative methods to solve the Schrödinger equation
* Part of NWChem, a computational quantum chemistry app from PNNL

# Implementations

## original
This is a simple serial C implentation of the TCE module, for use as a starting
point.


# Basis for comparison

To keep comparisons clean, the problem size, compiler version and CFLAGS should
be kept consistent for all builds.

## Compiler & CFLAGS

The code is built with clang++ version 9, and g++ version 8.3, with parameters
`-O3 -march=native -mtune=native`.

## Baseline

This implementation of TCE is a simplified form of the TCE feature in
[nwchem](https://github.com/nwchemgit/nwchem/).  This implementation uses
thread-level parallelism within the BLAS library matrix multiply calls.

## Input data

Our baseline performance metric will be a simulation of a Methane molecule (`C6H6`).


# Reference stuff

* [TCE homepage](https://www.csc.lsu.edu/~gb/TCE/)
* [TCE code within nwchem](https://github.com/nwchemgit/nwchem/tree/master/src/tce)
* [nwchem docs for TCE feature](http://www.nwchem-sw.org/index.php/Release66:TCE)
