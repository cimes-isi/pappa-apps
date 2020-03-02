/*
 * twoel.c
 *
 *  Created on: Jan 16, 2013
 *      Author: adam
 */

#include "cscc.h"
#include "g.h"
#include <omp.h>
#include <string.h>
double twoel_reducer(void);
double contract_matrices(double *, double *);


#define OFF(a, b)             ((a ##_off) + (b))
#ifdef SYMMETRY
#define DENS(a, b)            (g_dens_ ## a ## _ ## b)
#else
#define DENS(a, b)            (g_dens[OFF(a,b)])
#endif
#define DENS_DECL(a, b)       double g_dens_ ## a ## _ ## b = g_dens[OFF(a,b)]

#define BAD_CACHE_UTILIZATION
#ifdef BAD_CACHE_UTILIZATION
#define DENS_DECL_SYMM(a, b)  double g_dens_ ## a ## _ ## b = g_dens[OFF(a,b)]; \
                              double g_dens_ ## b ## _ ## a = g_dens[OFF(b,a)]
#else
#define DENS_DECL_SYMM(a, b)  double g_dens_ ## a ## _ ## b = g_dens[OFF(a,b)]; \
                              double g_dens_ ## b ## _ ## a = g_dens_ ## a ## _ ## b
#endif

#define TOL_DECL(a, b)        double tol2e_over_g_schwarz_ ## a ## _ ## b = tol2e / g_schwarz[OFF(a,b)]



#define UPDATE(a, b, c, d)    fock[OFF(a,b)] += (gg * DENS(c,d)); \
				   fock[OFF(a,c)] -= (0.50 * gg * DENS(b,d));


#define UPDATE1(a, b, c, d)    fock[OFF(a,b)] += (gg * DENS(c,d));
#define UPDATE2(a, b, c, d)    fock[OFF(a,c)] -= (0.50 * gg * DENS(b,d));

#define UPDATE_COULOMB1(a, b, c, d)  fock[OFF(a,b)] += (gg * DENS(c,d))
#define UPDATE_COULOMB2(a, b, c, d)  fock[OFF(a,b)] += (gg * DENS(c,d)) + (gg * DENS(c,d))
#define UPDATE_EXCHANGE(a, b, c, d)  fock[OFF(a,c)] -= (0.50 * gg * DENS(b,d));

#ifdef BOOKKEEPING
extern long long int icut1,icut2,icut3,icut4;
#define CUT_DECLS       long long int temp_icut1 = 0; \
                        long long int temp_icut2 = 0; \
                        long long int temp_icut3 = 0; \
                        long long int temp_icut4 = 0
#define CUT_WRITEBACK   icut1 += temp_icut1; \
                        icut2 += temp_icut2; \
                        icut3 += temp_icut3; \
                        icut4 += temp_icut4

#define CUT1(N)   temp_icut1 += N
#define CUT2(N)   temp_icut2 += N
#define CUT3(N)   temp_icut3 += N
#define CUT4(N)   temp_icut4 += N

#define ICUTVARS temp_icut1,temp_icut2,temp_icut3,temp_icut4

#define REDUCE_CLAUSE reduction(+:ICUTVARS)
#else
#define CUT1(N)
#define CUT2(N)
#define CUT3(N)
#define CUT4(N)
#define CUT_DECLS
#define CUT_WRITEBACK
#define ICUTVARS
#define REDUCE_CLAUSE
#endif

#ifndef NO_PRECALC
#define COMPUTE_G(a, b, c, d) g_fast(OFF(a, b), OFF(c, d))
#define CACHED_G(a, b, c, d) g_cached_fast(OFF(a, b), OFF(c, d))
#else
#define COMPUTE_G(a, b, c, d) g(a, b, c, d)
#define CACHED_G(a, b, c, d) g_cached(a, b, c, d)
#endif

#ifdef USE_CACHE
#define G CACHED_G
#else
#define G COMPUTE_G
#endif

#ifndef OMP_SCHEDULE
#define OMP_SCHEDULE schedule(runtime)
#endif





void twoel_i_j_k_l_all_different(double tol2e_over_schwmax) {
  CUT_DECLS;
  int j, k, l;
  int i_off, j_off, k_off, l_off;


#pragma omp parallel shared(g_fock, g_thread_fock, fm, g_dens,ICUTVARS) private( j, k, l, i_off, j_off, k_off, l_off)
{


#if OMP_VERSION > 2

double * fock = &g_thread_fock[omp_get_thread_num()*nbfn*nbfn];

#else

	#ifdef BAD_CACHE_UTILIZATION
	  double* fock = g_fock;
	#else
	  double* fock = g_partial_fock;
	#endif


#endif


  // 8-way symmetry
  i_off = 0;
  #pragma omp for OMP_SCHEDULE REDUCE_CLAUSE
  for (int i = 0; i < nbfn; i++) {
    i_off = i * nbfn;

#ifdef SYMMETRY
    j_off = i_off + nbfn;
    for (j = i + 1; j < nbfn; j++, j_off += nbfn) {
      DENS_DECL_SYMM(j,i);
#else
    j_off = 0;
    for (j = 0; j < nbfn; j++, j_off += nbfn) {
#endif /* SYMMETRY */

      TOL_DECL(i,j);

      if (g_schwarz[OFF(i,j)] < tol2e_over_schwmax) {
        CUT1(8 * (-2*j + i*(3-2*nbfn) + i*i + nbfn*(nbfn - 1))/2);
        continue;
      }

#ifdef SYMMETRY
      k_off = i_off;
      for (k = i; k < nbfn; k++, k_off += nbfn) {
        DENS_DECL_SYMM(k,i);
        DENS_DECL_SYMM(k,j);
#else
      k_off = 0;
      for (k = 0; k < nbfn; k++, k_off += nbfn) {
#endif /* SYMMETRY */

        if (g_schwarz_max_j[k] < tol2e_over_g_schwarz_i_j) {
          CUT4(8 * ((k==i) ? nbfn - j - 1: nbfn - k - 1));
          continue;
        }

#ifdef SYMMETRY
        l_off = ((k == i) ? j_off : k_off);
        for (l = 1 + ((k == i) ? j : k); l < nbfn; l++) {
          l_off += nbfn;
#else
        l_off = 0;
        for (l = 0; l < nbfn; l++, l_off += nbfn) {
#endif /* SYMMETRY */

          if (g_schwarz[OFF(k,l)] < tol2e_over_g_schwarz_i_j) {
            CUT2(8);
            continue;
          }

#ifdef SYMMETRY
          DENS_DECL_SYMM(l,i);
          DENS_DECL_SYMM(l,j);
          DENS_DECL_SYMM(l,k);
#endif

          CUT3(8);
          double gg = G(i, j, k, l);


#ifdef SYMMETRY
#if OMP_VERSION==1

	#ifdef BAD_CACHE_UTILIZATION
	#pragma omp atomic
	         UPDATE1(i,j,k,l);
	#pragma omp atomic
	         UPDATE1(i,j,l,k);
	#pragma omp atomic
	         UPDATE1(j,i,k,l);
	#pragma omp atomic
	         UPDATE1(j,i,l,k);
	#pragma omp atomic
	         UPDATE1(k,l,i,j);
	#pragma omp atomic
	         UPDATE1(k,l,j,i);
	#pragma omp atomic
	         UPDATE1(l,k,i,j);
	#pragma omp atomic
	         UPDATE1(l,k,j,i);

	#pragma omp atomic
	         UPDATE2(i,j,k,l);
	#pragma omp atomic
	         UPDATE2(i,j,l,k);
	#pragma omp atomic
	         UPDATE2(j,i,k,l);
	#pragma omp atomic
	         UPDATE2(j,i,l,k);
	#pragma omp atomic
	         UPDATE2(k,l,i,j);
	#pragma omp atomic
	         UPDATE2(k,l,j,i);
	#pragma omp atomic
	         UPDATE2(l,k,i,j);
	#pragma omp atomic
	         UPDATE2(l,k,j,i);


	#else
            /*
             * Exploiting yet another symmetry to make the faster and harder to read:
             * - g_dens is a symmetric matrix
             * - g_fock is a symmetric matrix
             *
             * The normal update sequence pulls identical values from both triangles
             * of g_dens and writes identical values to both triangles of g_fock.
             *
             * We can improve cache utilization by limiting or avoiding memory
             * access to one of the triangles.
             *
             * The following code separates the Coulomb and Exchange force update
             * so that we can use the symmetry of each to reduce memory operations
             * and limit access to the lower triangle. There are still some lower
             * triangle accesses, but they are significantly reduced.
             *
             * Once computed, the partial fock matrix will need some more processing
             * to make it symmetric again, but this can be done in an N^2 loop. To
             * get the final partial matrix, values from the upper and lower triangles
             * must be added together and written back to both triangles.
             *
             */

	#pragma omp atomic
	         UPDATE_COULOMB2(i,j,k,l);
	#pragma omp atomic
	         UPDATE_COULOMB2(k,l,i,j);

	#pragma omp atomic
	         UPDATE_EXCHANGE(i,j,k,l);
	#pragma omp atomic
	        UPDATE_EXCHANGE(i,j,l,k);
	#pragma omp atomic
	        UPDATE_EXCHANGE(j,i,k,l);
	#pragma omp atomic
	        UPDATE_EXCHANGE(j,i,l,k);

	#endif

#else

	#if OMP_VERSION==2
	#pragma omp critical
	{
	#endif

	#ifdef BAD_CACHE_UTILIZATION
           UPDATE(i,j,k,l);
           UPDATE(i,j,l,k);
           UPDATE(j,i,k,l);
           UPDATE(j,i,l,k);
           UPDATE(k,l,i,j);
           UPDATE(k,l,j,i);
           UPDATE(l,k,i,j);
           UPDATE(l,k,j,i);

	#else
            /*
             * Exploiting yet another symmetry to make the faster and harder to read:
             * - g_dens is a symmetric matrix
             * - g_fock is a symmetric matrix
             *
             * The normal update sequence pulls identical values from both triangles
             * of g_dens and writes identical values to both triangles of g_fock.
             *
             * We can improve cache utilization by limiting or avoiding memory
             * access to one of the triangles.
             *
             * The following code separates the Coulomb and Exchange force update
             * so that we can use the symmetry of each to reduce memory operations
             * and limit access to the lower triangle. There are still some lower
             * triangle accesses, but they are significantly reduced.
             *
             * Once computed, the partial fock matrix will need some more processing
             * to make it symmetric again, but this can be done in an N^2 loop. To
             * get the final partial matrix, values from the upper and lower triangles
             * must be added together and written back to both triangles.
             *
             */

            UPDATE_COULOMB2(i,j,k,l);
            UPDATE_COULOMB2(k,l,i,j);

            UPDATE_EXCHANGE(i,j,k,l);
           UPDATE_EXCHANGE(i,j,l,k);
            UPDATE_EXCHANGE(j,i,k,l);
           UPDATE_EXCHANGE(j,i,l,k);
	#endif


	#if OMP_VERSION==2
	}
	#endif
#endif//the not version 1 case
#else
	         UPDATE1(i,j,k,l);
	         UPDATE2(i,j,k,l);
#endif

        } // l
      } // k
    } // j
  } // i


#if OMP_VERSION>2
{

	//reassign fock pointer for reduction
	#ifdef BAD_CACHE_UTILIZATION
		double *  fock = g_fock;
	#else
		double *  fock = g_partial_fock;
	#endif

        long long int nbfn_squared=nbfn*nbfn;
        int numthreads=omp_get_max_threads();

	//reduce all the threaded focks
	#pragma omp for private(j) schedule(static)
	for(int i=0;i<nbfn_squared;i++)
	{
		double *tmp = &g_thread_fock[i];
  		double sum = 0.0;
		for(j=0;j<numthreads;j++)
		{
			sum+=tmp[0];
			tmp[0]=0;
			tmp=&tmp[nbfn_squared];
		}

		fock[i]+=sum;
	}


}

#endif//Version 3 reduction and cleanup finished






// reconstruct the full g_fock using the partial and some symmetry
#ifndef BAD_CACHE_UTILIZATION
  #pragma omp for OMP_SCHEDULE private(j,i_off,j_off)
  for (int i = 0; i < nbfn; i++) {

    i_off=i*nbfn;
    // Treat the diagonal special, we need to double it's value but if we
    // included it in the j loop, we would end up quadrupling the value.
    double fock_val = g_partial_fock[i_off + i];

    g_fock[i_off + i] += fock_val + fock_val;
    g_partial_fock[i_off + i] = 0.0;

    j_off = i_off + nbfn;
    for (j = i+1; j < nbfn; j++) {
      fock_val = g_partial_fock[i_off + j] + g_partial_fock[j_off + i];
      g_fock[i_off + j] += fock_val;
      g_fock[j_off + i] += fock_val;
      g_partial_fock[i_off + j] = 0.0;
      g_partial_fock[j_off + i] = 0.0;


      j_off += nbfn;
    }

  }
#endif
}//end of omp parallel
  CUT_WRITEBACK;
}

#ifdef SYMMETRY
void twoel_i_eq_j(double tol2e_over_schwmax) {
  CUT_DECLS;
  int i, k, l;
  int i_off, k_off, l_off;

#pragma omp parallel shared(g_fock, g_thread_fock, fm, g_dens,ICUTVARS) private( k, l, i_off, k_off, l_off)
{


#if OMP_VERSION > 2

double * fock = &g_thread_fock[omp_get_thread_num() *nbfn*nbfn];

#else

	#ifdef BAD_CACHE_UTILIZATION
	  double* fock = g_fock;
	#else
	  double* fock = g_partial_fock;
	#endif

#endif

  // 4-way symmetry
  #pragma omp for OMP_SCHEDULE REDUCE_CLAUSE nowait
  for (i = 0; i < nbfn; i++) {
    i_off = i * nbfn;

    DENS_DECL(i,i);
    TOL_DECL(i,i);

    if (g_schwarz[OFF(i,i)] < tol2e_over_schwmax) {
      CUT1(4 * (1 + i - nbfn) * (i - nbfn) / 2);
      continue;
    }
    k_off = i_off - nbfn;

    for (k = i; k < nbfn; k++) {
      k_off += nbfn;
      DENS_DECL_SYMM(k,i);
      if (g_schwarz_max_j[k] < tol2e_over_g_schwarz_i_i) {
        CUT4(4 * (nbfn - k - 1));
        continue;
      }

      l_off = k_off;
      for (l = 1 + k; l < nbfn; l++) {
        l_off += nbfn;

        if (g_schwarz[OFF(k,l)] < tol2e_over_g_schwarz_i_i) {
          CUT2(4);
          continue;
        }

        DENS_DECL_SYMM(l,i);
        DENS_DECL_SYMM(l,k);

        CUT3(4);
        double gg = G(i, i, k, l);

#if OMP_VERSION==1

	#ifdef BAD_CACHE_UTILIZATION
	#pragma omp atomic
	         UPDATE1(i,i,k,l);
	#pragma omp atomic
	         UPDATE1(i,i,l,k);
	#pragma omp atomic
	         UPDATE1(k,l,i,i);
	#pragma omp atomic
	         UPDATE1(l,k,i,i);

	#pragma omp atomic
	         UPDATE2(i,i,k,l);
	#pragma omp atomic
	         UPDATE2(i,i,l,k);
	#pragma omp atomic
	         UPDATE2(k,l,i,i);
	#pragma omp atomic
	         UPDATE2(l,k,i,i);

	#else
	#pragma omp atomic
        UPDATE_COULOMB1(i,i,k,l);
	#pragma omp atomic
        UPDATE_COULOMB1(k,l,i,i);

	#pragma omp atomic
        UPDATE_EXCHANGE(i,i,k,l);
	#pragma omp atomic
        UPDATE_EXCHANGE(i,i,l,k);
	#endif

#else

	#if OMP_VERSION==2
	#pragma omp critical
	{
	#endif

	#ifdef BAD_CACHE_UTILIZATION
        UPDATE(i,i,k,l);
        UPDATE(i,i,l,k);
        UPDATE(k,l,i,i);
        UPDATE(l,k,i,i);
	#else
        UPDATE_COULOMB1(i,i,k,l);
        UPDATE_COULOMB1(k,l,i,i);

        UPDATE_EXCHANGE(i,i,k,l);
        UPDATE_EXCHANGE(i,i,l,k);
	#endif

	#if OMP_VERSION==2
	}
	#endif

#endif//end of version 1 ifcase
      } // l
    } // k
  } // i


}//end of parallel section (remember for these functions it is all nowait)



  CUT_WRITEBACK;
}

void twoel_k_eq_l(double tol2e_over_schwmax) {
  CUT_DECLS;
  int i, j, k;
  int i_off, j_off, k_off;



#pragma omp parallel shared(g_fock, g_thread_fock, fm, g_dens,ICUTVARS) private( j, k, i_off, j_off, k_off)
{


#if OMP_VERSION > 2

double * fock = &g_thread_fock[omp_get_thread_num()*nbfn*nbfn];

#else


	#ifdef BAD_CACHE_UTILIZATION
	  double* fock = g_fock;
	#else
	  double* fock = g_partial_fock;
	#endif
#endif
  // 4-way symmetry
  #pragma omp for OMP_SCHEDULE REDUCE_CLAUSE nowait
  for (i = 0; i < nbfn; i++) {
    i_off = i * nbfn;
    j_off = i_off;

    for (j = i + 1; j < nbfn; j++) {
      j_off += nbfn;
      DENS_DECL_SYMM(j,i);
      TOL_DECL(i,j);

      if (g_schwarz[OFF(i,j)] < tol2e_over_schwmax) {
        CUT1(4 * (nbfn - i - 1));
        continue;
      }
      k_off = i_off;

      for (k = i + 1; k < nbfn; k++) {
        k_off += nbfn;

        DENS_DECL_SYMM(k,i);
        DENS_DECL_SYMM(k,j);
        DENS_DECL(k,k);

        if (g_schwarz[OFF(k,k)] < tol2e_over_g_schwarz_i_j) {
          CUT2(4);
          continue;
        }


        CUT3(4);
        double gg = G(i, j, k, k);

#if OMP_VERSION==1

	#ifdef BAD_CACHE_UTILIZATION
	#pragma omp atomic
	UPDATE1(i,j,k,k);
	#pragma omp atomic
	UPDATE1(j,i,k,k);
	#pragma omp atomic
        UPDATE1(k,k,i,j);
	#pragma omp atomic
        UPDATE1(k,k,j,i);

	#pragma omp atomic
	UPDATE2(i,j,k,k);
	#pragma omp atomic
        UPDATE2(j,i,k,k);
	#pragma omp atomic
        UPDATE2(k,k,i,j);
	#pragma omp atomic
        UPDATE2(k,k,j,i);


	#else

	#pragma omp atomic
         UPDATE_COULOMB1(i,j,k,k);
	#pragma omp atomic
         UPDATE_COULOMB1(k,k,i,j);

	#pragma omp atomic
         UPDATE_EXCHANGE(i,j,k,k);
	#pragma omp atomic
         UPDATE_EXCHANGE(j,i,k,k);

	#endif
#else



	#if OMP_VERSION==2
	#pragma omp critical
	{
	#endif

	#ifdef BAD_CACHE_UTILIZATION
        UPDATE(i,j,k,k);
        UPDATE(j,i,k,k);
        UPDATE(k,k,i,j);
        UPDATE(k,k,j,i);
	#else
        UPDATE_COULOMB1(i,j,k,k);
        UPDATE_COULOMB1(k,k,i,j);

        UPDATE_EXCHANGE(i,j,k,k);
        UPDATE_EXCHANGE(j,i,k,k);
	#endif

	#if OMP_VERSION==2
	}
	#endif

#endif//end of version1 update cases

      } // k
    } // j
  } // i



}//omp parallel


  CUT_WRITEBACK;
}

void twoel_ij_eq_kl(double tol2e_over_schwmax) {
  CUT_DECLS;
  int i, j;
  int i_off, j_off;
  double* fock = g_fock;

  // 4-way symmetry
  i_off = -nbfn;
  for (i = 0; i < nbfn; i++) {
    i_off += nbfn;
    j_off = i_off;

    DENS_DECL(i,i);

    for (j = i + 1; j < nbfn; j++) {
      j_off += nbfn;
      DENS_DECL_SYMM(j,i);
      DENS_DECL(j,j);
      TOL_DECL(i,j);

      if (g_schwarz[OFF(i,j)] < tol2e_over_g_schwarz_i_j) {
        CUT1(4);
        continue;
      }

      CUT3(4);
      double gg = G(i, j, i, j);

      UPDATE(i,j,i,j);
      UPDATE(i,j,j,i);
      UPDATE(j,i,i,j);
      UPDATE(j,i,j,i);
    } // j
  } // i

  CUT_WRITEBACK;
}

void twoel_i_eq_j_and_k_eq_l(double tol2e_over_schwmax) {
  CUT_DECLS;
  int i, k;
  int i_off, k_off;
  double* fock = g_fock;

  // 2-way symmetry
  i_off = -nbfn;
  for (i = 0; i < nbfn; i++) {
    i_off += nbfn;

    DENS_DECL(i,i);
    TOL_DECL(i,i);

    if (g_schwarz[OFF(i,i)] < tol2e_over_schwmax) {
      CUT1(2 * (nbfn - i - 1));
      continue;
    }
    k_off = i_off;

    for (k = i + 1; k < nbfn; k++) {
      k_off += nbfn;
      DENS_DECL_SYMM(k,i);
      DENS_DECL(k,k);

      if (g_schwarz[OFF(k,k)] < tol2e_over_g_schwarz_i_i) {
        CUT2(2);
        continue;
      }

      CUT3(2);
      double gg = G(i, i, k, k);

      UPDATE(i,i,k,k);
      UPDATE(k,k,i,i);
    } // k
  } // i

  CUT_WRITEBACK;
}


void twoel_i_eq_j_eq_k_eq_l(double tol2e_over_schwmax) {
  CUT_DECLS;
  int i;
  int i_off;
  double* fock = g_fock;

  // 1-way symmetry
  i_off = -nbfn;
  for (i = 0; i < nbfn; i++) {
    i_off += nbfn;

    DENS_DECL(i,i);
    TOL_DECL(i,i);

    if (g_schwarz[OFF(i,i)] < tol2e_over_g_schwarz_i_i) {
      CUT2(1);
      continue;
    }

    CUT3(1);
    double gg = G(i, i, i, i);

    UPDATE(i,i,i,i);
  } // i

  CUT_WRITEBACK;
}
#endif /* SYMMETRY */

double twoel_fast(double schwmax) {
  double tol2e_over_schwmax = tol2e / schwmax;

#ifdef SYMMETRY
#pragma omp single nowait
{
  // do the N^2 and smaller loops which directly accumulate to g_fock
  //this is done single threaded so we can save memory allocations
  twoel_ij_eq_kl(tol2e_over_schwmax);
  twoel_i_eq_j_and_k_eq_l(tol2e_over_schwmax);
  twoel_i_eq_j_eq_k_eq_l(tol2e_over_schwmax);
}
  // do the N^3 loops first since they write to the g_partial_fock matrix

  twoel_i_eq_j(tol2e_over_schwmax);
  twoel_k_eq_l(tol2e_over_schwmax);
#endif /* SYMMETRY */

  // do the N^4 loop which will write to the g_partial_fock matrix, then
  // correctly accumulate the values back to g_fock at the end

  twoel_i_j_k_l_all_different(tol2e_over_schwmax);

#ifdef USE_CACHE
  reset_cache_index();
#endif

  return (0.50 * contract_matrices(g_fock, g_dens));
}
