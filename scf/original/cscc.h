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

#define mxiter 30
#define maxatom 384
#define maxnbfn 15 * maxatom
#define pi 3.141592653589793
#define tol 0.006
#define tol2e 0.000001
#define EPSILON 0.0000000001

#define MAX(a,b) (((a) >= (b)) ? (a) : (b))
#define MIN(a,b) (((a) <= (b)) ? (a) : (b))

double enrep, q[maxatom], ax[maxatom], ay[maxatom], az[maxatom];
double x[maxnbfn], y[maxnbfn], z[maxnbfn], expnt[maxnbfn], rnorm[maxnbfn];
double ** g_dens, ** g_fock, ** g_tfock, ** g_schwarz, ** g_work, ** g_ident, ** g_orbs, * eigv;

int iky[maxnbfn], nocc, nbfn;
long long int icut1,icut2,icut3,icut4; 
int natom;

