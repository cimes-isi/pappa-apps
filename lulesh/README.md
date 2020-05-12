# Intro

LULESH

Livermore Unstructured Legrangian Explicit Shock Hydrodynamics

# Implementations

## original
This is the original code base.


# Basis for comparison

To keep comparisons clean, the problem size, compiler version and CFLAGS should
be kept consistent for all builds.

## Compiler & CFLAGS

The code is built with clang++ version 9, and g++ version 8.3, with parameters
`-O3 -march=native -mtune=native`.

## Baseline

Our baseline performance metric will be simulation of a 40x40x40 cubical space.


# Reference stuff

* [hpgmg website](https://hpgmg.org/)
* [hpgmg code](https://bitbucket.org/hpgmg/hpgmg)
