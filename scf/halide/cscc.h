/*
* Copyright (c) 2013 Battelle Memorial Institute.
* All rights reserved.
*
* The software in this package is published under the terms of the BSD
* style license a copy of which has been included with this distribution in
* the LICENSE.txt file.
*
* Created on 2013 by John Feo
*/

#include <HalideBuffer.h>

#define mxiter 30
#define maxatom 384
#define maxnbfn 15 * maxatom
#define pi 3.141592653589793
#define tol 0.006
#define tol2e 0.000001
#define EPSILON 0.0000000001

#define MAX(a,b) (((a) >= (b)) ? (a) : (b))
#define MIN(a,b) (((a) <= (b)) ? (a) : (b))

extern double enrep, q[maxatom], ax[maxatom], ay[maxatom], az[maxatom];
extern Halide::Runtime::Buffer<double, 1> x_buf, y_buf, z_buf, expnt_buf, rnorm_buf, g_schwarz_max_j_buf;
extern double *x, *y, *z, *expnt, *rnorm, *g_schwarz_max_j;
extern Halide::Runtime::Buffer<double, 2> g_dens_buf, g_fock_double_buf[2], *g_fock_buf, *g_fock_out_buf, g_schwarz_buf, fm_buf;
extern double ** g_dens, ** g_fock, ** g_fock_out, ** g_schwarz, ** g_tfock, ** g_work, ** g_ident, ** g_orbs, * eigv;

extern int iky[maxnbfn], nocc, nbfn;
extern long long int icut1,icut2,icut3,icut4;
extern int natom;

