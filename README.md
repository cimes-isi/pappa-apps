# pappa-apps

This repo contains copies of proxy apps to consider, and build scripts
to set them up and run them.

Each app has a separate folder, with a Makefile.  You should be able to
just run "make" to run the app.  The Makefile knows how to build the
app, set up datasets if needed, and run it.  For python apps (sar,
mot), the Makefile will also set up a virtualenv folder with the
necessary third party libraries.


# HPGMG

High-performance Geometric Multigrid

https://hpgmg.org/
https://bitbucket.org/hpgmg/hpgmg

# LULESH

Laurence Unstructured Lagrangian-Eulerian Shock Hydrodynamics

https://computing.llnl.gov/projects/co-design/lulesh
https://github.com/LLNL/LULESH

# MOT

Multiple-Object Tracking is a computer vision problem where multiple
objects (such as pedestrians) are detected and their motion tracked
over time.

https://motchallenge.net/
https://github.com/conan7882/GoogLeNet-Inception
https://github.com/nwojke/deep_sort

# SAR

Synthetic Aperture RADAR is a mapping algorithm which stitches together
large amounts of raw RADAR data to assemble a high quality map of the
ground surface.

https://github.com/dm6718/RITSAR

# SCF

Self-Consistent Field, a.k.a. Hartree Fock.  This is a quantum chemistry
algorithm, from NWChem.

# TCE

Tensor Contraction Engine, generating Coupled Cluster method code.
This is a quantum chemistry algorithm, from NWChem.
