/*
 * Copyright (c) 2008-2012
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rgemm_dd.cpp,v 1.4 2010/08/07 05:50:09 nakatamaho Exp $
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */

#include <stdio.h>
#include <string.h>
#include <dlfcn.h>
#include <mblas.h>
#include <mlapack.h>
#include <mpack_benchmark.h>

#define TOTALSTEPS 1000

double flops_getrf(mpackint m_i, mpackint n_i)
{
  double adds, muls, flops;
  double m, n;
  m = (double)m_i;   n = (double)n_i; 
  muls = 0.5 * m * n * n - (1./6.) * n * n * n + 0.5 * m * n - 0.5 * n * n + (2./3.) * n;
  adds = 0.5 * m * n * n - (1./6.) * n * n * n - 0.5 * m * n + (1./6.) * n;
  flops = muls + adds;
  return flops;
}

#define TOTALSTEPS 1000
int main (int argc, char *argv[])
{
  REAL alpha, beta, mtemp, dummy;
  REAL *dummywork;
  double elapsedtime, t1, t2;
  char uplo, normtype;
  mpackint N0, M0, STEPN, STEPM, info;
  mpackint lda;
  int i, j, m, n, k, ka, kb, p, q;
  int check_flag = 1;

  ___MPACK_INITIALIZE___

  const char mlapack_sym[] = SYMBOL_GCC_RGETRF;
  const char raxpy_sym[] = SYMBOL_GCC_RAXPY;
  void *handle;
  void (*mlapack_ref) (mpackint, mpackint, REAL *, mpackint, mpackint *, mpackint *);
  void (*raxpy_ref) (mpackint, REAL, REAL *, mpackint, REAL *, mpackint);
  char *error;
  REAL diff;
  double diffr;

  // initialization
  N0 = 1;
  M0 = 1;
  STEPN = 1;
  STEPM = 1;
  if (argc != 1) {
    for (i = 1; i < argc; i++) {
      if (strcmp ("-STEPN", argv[i]) == 0) {
	STEPN = atoi (argv[++i]);
      } 
      else if (strcmp ("-STEPM", argv[i]) == 0) {
	STEPM = atoi (argv[++i]);
      }
      else if (strcmp ("-N0", argv[i]) == 0) {
	N0 = atoi (argv[++i]);
      }
      else if (strcmp ("-M0", argv[i]) == 0) {
	M0 = atoi (argv[++i]);
      }
      else if (strcmp ("-NOCHECK", argv[i]) == 0) {
	check_flag = 0;
      }
    }
  }

  if (check_flag) {
    handle = dlopen (MLAPACK_REF_LIB DYLIB_SUFFIX, RTLD_LAZY);
    if (!handle) {
      printf ("dlopen: %s\n", dlerror ());
      return 1;
    }
    mlapack_ref = (void (*)(mpackint, mpackint, REAL *, mpackint, mpackint *, mpackint *)) dlsym (handle, mlapack_sym);
    if ((error = dlerror ()) != NULL) {
      fprintf (stderr, "%s\n", error);
      return 1;
    }

    handle = dlopen (MBLAS_REF_LIB DYLIB_SUFFIX, RTLD_LAZY);
    if (!handle) {
      printf ("dlopen: %s\n", dlerror ());
      return 1;
    }
    raxpy_ref = (void (*)(mpackint, REAL, REAL *, mpackint, REAL *, mpackint)) dlsym (handle, raxpy_sym);
    if ((error = dlerror ()) != NULL) {
      fprintf (stderr, "%s\n", error);
      return 1;
    }
  }

  n = N0; 
  m = M0; 
  for (p = 0; p < TOTALSTEPS; p++) {
    lda = m;
    REAL *A = new REAL[lda * n];
    REAL *Ad = new REAL[lda * n];
    mpackint *ipiv = new mpackint[min(m,n)];
    mpackint *ipivd = new mpackint[min(m,n)];
    REAL mOne = -1;
    for (i = 0; i < lda * n; i++) {
      A[i] = Ad[i] = randomnumber (dummy);
    }

    if (check_flag) {
      t1 = gettime ();
      Rgetrf (m, n, A, lda, ipiv, &info);
      t2 = gettime ();
      elapsedtime = (t2 - t1);
      (*mlapack_ref) (m, n, Ad, lda, ipivd, &info);
      (*raxpy_ref) ((mpackint) (lda * n), mOne, A, (mpackint) 1, Ad, (mpackint) 1);
      diff = Rlange (&normtype, (mpackint) lda, (mpackint) n, Ad, lda, dummywork);
      diffr = cast2double (diff);
      printf ("    n     MFLOPS     error   uplo\n");
      printf ("%5d %10.3f %5.2e      %c\n", (int) n, flops_getrf(m,n) / elapsedtime * MFLOPS, diffr, uplo);
    }
    else {
      t1 = gettime ();
      Rgetrf (m, n, A, lda, ipiv, &info);
      t2 = gettime ();
      elapsedtime = (t2 - t1);
      printf ("    n     MFLOPS   uplo\n");
      printf ("%5d %10.3f      %c\n", (int) n, flops_getrf(m,n) / elapsedtime * MFLOPS, uplo);
    }
    delete[]ipivd;
    delete[]ipiv;
    delete[]Ad;
    delete[]A;
    n = n + STEPN;
    m = m + STEPM;
  }
  if (check_flag)
    dlclose (handle);
}
