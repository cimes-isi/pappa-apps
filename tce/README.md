# Intro

TCE

* Tensor Contraction Engine
* Implements multiple iterative methods to solve the Schrödinger equation
* Part of NWChem, a computational quantum chemistry app from PNNL

# What is this

This is a modified version of TCE which generates C code, and a
standalone test harness to run that C code and verify its output.
It is intended for use as a research proxy application, so that
language and runtime researchers can work with and improve the TCE
code generation, without having to dive too deeply into the guts of
nwchem.

Hopefully, improvements found here can be applied to nwchem as well.


# Details

[Author's note: I am not a quantum chemist or physicist, so I don't
have a good idea what's happening on that side of things.  Apologies
if my description below takes all the fun out of it.]

The files `cc2_t1.tt` and `cc2_t2.tt` define a set of tensor operations
which embody a correlation model (CC2), which provides a solution to a
quantum chemistry problem.  The TCE python library (tce.py) reads one of
these files and does some symmetry analysis, reducing the amount of
duplicate computation and deciding what order to do the resulting
operations in, to minimize some computational and memory cost functions.
It then generates Fortran (`cc2_t1.F`, `cc2_t2.F`) or C (`cc2_t1.c`,
`cc2_t2.c`) source files to implement this set of operations efficiently.
A second correlation model, CCSD, is also provided in this example, with
inputs and source files following the same naming scheme.

Many other correlation models exist and are supported by TCE; see
[here](http://www.nwchem-sw.org/index.php/Release66:TCE#Overview) for a list.

The Fortran output can be built into nwchem unchanged; the C code expects
to run in an environment similar to nwchem.  In both cases, they expect the
input data to be in the form that nwchem uses, which is a block-sparse tensor
format with associated details about index spaces, block sizes, symmetry
domains, and restrictions on the kinds of interactions, most of which is
defined by the nwchem input file.

One such input file is data/test_cc2_o3/test_cc2.nw, which defines a simple
problem (running CC2 on an ozone molecule) for testing purposes.  Nwchem
solves this problem by invoking cc2_t1 and cc2_t2 iteratively, until the
output converges.  For this input file, convergence happens in 21
iterations.  Various other pieces of input are also provided to these
functions; I have included the other data in other files in the
data/test_cc2_o3/ directory.

A more complex problem can be found in data/test_ccsd_c6h6_2eorb.  This
data runs CCSD on a benzene molecule, and enables the "2eorb" feature,
which calculates some of the larger inputs on the fly, reducing the
memory footprint.  See
[here](http://www.nwchem-sw.org/index.php/Release66:TCE#2EORB_--_alternative_storage_of_two-electron_integrals)
for a description of that.

In practical terms, cc2_t1 produces an output matrix, r1, which contains
the residuals from an input matrix, t1.  Likewise, cc2_t2 produces an
output tensor, r2, containing the residuals from another input tensor,
t2.  At the end of each iteration, r1 and r2 are added (or whatever)
to t1 and t2.  Nwchem detects convergence when r1 and r2 get
sufficiently small.

For the C version, I have included a standalone test program,
`generic_t1_t2_standalone.c`, which reads all of the input files, sets
up enough of an nwchem-like execution environment for the functions to
work, then calls them with the input and compares the output to known
good output.  If the output is off by more than a very small amount
(epsilon), it will complain and return non-zero.  Note that nwchem
iterates over these functions multiple times until convergence; this
test program only runs them once.

The Makefile is set up to run TCE to generate C sources, compile and
link the test program, and run it with the right data files.  All you
should have to do is type "make run" to verify that everything runs and
produces the correct result.


# Some Performance Notes

The generated code follows a common pattern: loop over blocks of input,
then for each block, copy the block into local memory, then call a BLAS
matrix multiply function, then add the result into an output buffer.  The
input/output "copy" operations often permute indexes, rearranging the data
as it is copied.  The indices are carefully ordered so that matrix multiply
does the right thing for higher-dimensional tensor blocks... for instance,
one "vector" of the "matrix" may be a 2-dimensional chunk of a 4d input
tensor.

Much of the heavy lifting is done with the BLAS `dgemm()` function, so
your choice of BLAS library will affect the performance.  If `$MKLROOT`
is set, MKL is used automatically.  Otherwise it attempts to use openblas.
This can be changed by modifying `Makefile` and `tce.h`.  MKL, ATLAS BLAS
and OpenBLAS are known to work.

The generated code calls some functions to perform operations like
fetching, sorting, and generating blocks of input.  There are C
implementations of these functions in generic_t1_t2_standalone.c, and most
of them are direct translations of the corresponding Fortran code.  None
of these functions are multithreaded; some of them may be good candidates
for improvement.

The copy/permute functions are called `tce_sort_2` and `tce_sort_4`.  In
particular, the 4-dimensional sort function is much simpler than its NWChem
equivalent, and can be optimized further.

For larger datasets, the `v2` input tensor becomes too large to fit in
memory, and so blocks of `v2` input are calculated on the fly.  This is
done by the `tce_get_block_ind_i` function.  This adds overhead, and can
probably be optimized further.


# Other Stuff

The data files are hosted externally, as they are inconveniently
large.  They will be downloaded and unpacked automatically, the first
time you run "make" in a new git checkout.

In tex/, you can find some tex output which I fixed up to be slightly
more readable, and pdf files built from that tex.  This provides a
mathematical view of the structure of the cc2 and ccsd operations.  In
the PDF, each line corresponds to a function in the generated source
code.


# License

This software is derived from NWChem, a computational chemistry package
maintained by PNNL.  As it is a derivative work, it is subject to the
terms of NWChem's license agreement.  See LICENSE.nwchem for those terms.

Our modifications and additions are Copyright (C) 2013-2014 ET
International, Inc., and are made available under the same license terms
as NWChem.

There is more information on TCE [here](http://www.csc.lsu.edu/~gb/TCE/).

There is more information on NWChem [here](http://www.nwchem-sw.org/).
