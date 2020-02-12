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

#include <math.h>
#include "cscc.h"

double fm[2001][5];
double rdelta, delta, delo2;

double exprjh(double x) {

  if (x < -37.0) return 0.0;
  else           return exp(x);
}

void setfm(void)
{
  int i, ii;
  double t[2001];
  double et[2001], rr, tt;

  delta  = 0.014;
  delo2  = delta * 0.5;
  rdelta = 1.0 / delta;

  for (i = 0; i < 2001; i++) {
    tt       = delta * (double)i;
    et[i]    = exprjh(-tt);
    t[i]     = tt * 2.0;
    fm[i][4] = 0.0;
  }

  for (i = 199; i > 3; i--) {
    rr = 1.0 / (double)(2 * i + 1);
    for (ii = 0; ii < 2001; ii++) fm[ii][4] = (et[ii] + t[ii] * fm[ii][4]) * rr;
  }

  for (i = 3; i >= 0; i--) {
    rr = 1.0 / (double) (2 * i + 1);
    for (ii = 0; ii < 2001; ii++) fm[ii][i] = (et[ii] + t[ii] * fm[ii][i+1]) * rr;
  }

  return;
}


// computes f0 to a relative accuracy of better than 4.e-13 for all t. Uses 4th order taylor
// expansion on grid out to t = 28.0 asymptotic expansion accurate for t greater than 28
void  f0(double *f0val, double t) {
  if (t >= 28.0) { *f0val = 0.88622692545276 / sqrt(t); return;}

  int    n    = (int) ((t + delo2) * rdelta);
  double x    = delta * (double) n - t;
  double *fmn = fm[n];

  *f0val   = fmn[0] + x * (fmn[1] + x * 0.5  * (fmn[2] + x * (1./3.) * (fmn[3] + x * 0.25 * fmn[4])));
  return;
}


void Dgemm(char a, char b, int M, int N, int K, double **A, double **B, double **C) {
  int i, j, l;

  if ((a == 'n') && (b == 'n')) {

     for (i = 0; i < M; i++) {
     for (j = 0; j < N; j++) {
       double sum = 0.0;
       for (l = 0; l < K; l++) sum += A[i][l] * B[l][j];
       C[i][j] = sum;
     } }

  } else if ((a == 'n') && (b == 't')) {

     for (i = 0; i < M; i++) {
     for (j = 0; j < N; j++) {
       double sum = 0.0;
       for (l = 0; l < K; l++) sum += A[i][l] * B[j][l];
       C[i][j] = sum;
     } }

  } else if ((a == 't') && (b == 'n')) {

     for (i = 0; i < M; i++) {
     for (j = 0; j < N; j++) {
       double sum = 0.0;
       for (l = 0; l < K; l++) sum += A[l][i] * B[l][j];
       C[i][j] = sum;
     } }

  } else if ((a == 't') && (b == 't')) {

     for (i = 0; i < M; i++) {
     for (j = 0; j < N; j++) {
       double sum = 0.0;
       for (l = 0; l < K; l++) sum += A[l][i] * B[j][l];
       C[i][j] = sum;
     } }

} }
