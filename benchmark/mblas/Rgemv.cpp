/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rgemv_qd.cpp,v 1.3 2010/08/07 05:50:09 nakatamaho Exp $
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

int
main (int argc, char *argv[])
{
  mpackint k, l, m, n;
  mpackint incx = 1, incy = 1, STEPN, STEPM, N0, M0;
  REAL alpha, beta, dummy, *dummywork;
  REAL mOne = -1;
  double elapsedtime, t1, t2;
  int i, p;
  int check_flag = 1;
  char trans, normtype;

  ___MPACK_INITIALIZE___

  const char mblas_sym[] = SYMBOL_GCC_RGEMV;
  const char raxpy_sym[] = SYMBOL_GCC_RAXPY;
  void *handle;
  void (*mblas_ref) (const char *, mpackint, mpackint, REAL, REAL *, mpackint, REAL *, mpackint, REAL, REAL *, mpackint);
  void (*raxpy_ref) (mpackint, REAL, REAL *, mpackint, REAL *, mpackint);
  char *error;
  REAL diff;
  double diffr;

  // initialization
  N0 = M0 = 1;
  STEPN = STEPM = 1;
  trans = 'n';
  normtype = 'm';
  if (argc != 1) {
    for (i = 1; i < argc; i++) {
      if (strcmp ("-N", argv[i]) == 0) {
	N0 = atoi (argv[++i]);
      }
      else if (strcmp ("-M", argv[i]) == 0) {
	M0 = atoi (argv[++i]);
      }
      else if (strcmp ("-STEPN", argv[i]) == 0) {
	STEPN = atoi (argv[++i]);
      }
      else if (strcmp ("-STEPM", argv[i]) == 0) {
	STEPM = atoi (argv[++i]);
      }
      else if (strcmp ("-T", argv[i]) == 0) {
	trans = 't';
      }
      else if (strcmp ("-NOCHECK", argv[i]) == 0) {
	check_flag = 0;
      }
    }
  }
  if (check_flag) {
    handle = dlopen (MBLAS_REF_LIB DYLIB_SUFFIX, RTLD_LAZY);
    if (!handle) {
      printf ("dlopen: %s\n", dlerror ());
      return 1;
    }
    mblas_ref = (void (*)(const char *, mpackint, mpackint, REAL, REAL *, mpackint, REAL *, mpackint, REAL, REAL *, mpackint)) dlsym (handle, mblas_sym);
    if ((error = dlerror ()) != NULL) {
      fprintf (stderr, "%s\n", error);
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
    if (Mlsame (&trans, "n")) {
      k = n;
      l = m;
    }
    else {
      k = m;
      l = n;
    }
    REAL *x = new REAL[k];
    REAL *y = new REAL[l];
    REAL *yd = new REAL[l];
    REAL *A = new REAL[n * m];
    if (check_flag) {
      for (i = 0; i < k; i++) {
	x[i] = randomnumber (dummy);
      }
      for (i = 0; i < l; i++) {
	y[i] = yd[i] = randomnumber (dummy);
      }
      for (i = 0; i < k * l; i++) {
	A[i] = randomnumber (dummy);
      }
      alpha = randomnumber (dummy);
      beta = randomnumber (dummy);
      t1 = gettime ();
      Rgemv (&trans, m, n, alpha, A, m, x, (mpackint) 1, beta, y, (mpackint) 1);
      t2 = gettime ();
      elapsedtime = (t2 - t1);
      (*mblas_ref) (&trans, m, n, alpha, A, m, x, (mpackint) 1, beta, yd, (mpackint) 1);
      (*raxpy_ref) (l, mOne, y, (mpackint) 1, yd, (mpackint) 1);
      diff = Rlange (&normtype, (mpackint) l, (mpackint) 1, yd, 1, dummywork);
      diffr = cast2double (diff);
      printf ("     m       n      MFLOPS      error    trans\n");
      printf ("%6d  %6d  %10.3f   %5.2e        %c\n", (int) n, (int) m, (2.0 * (double) n * (double) m) / elapsedtime * MFLOPS, diffr, trans);
    }
    else {
      for (i = 0; i < k; i++) {
	x[i] = randomnumber (dummy);
      }
      for (i = 0; i < l; i++) {
	y[i] = yd[i] = randomnumber (dummy);
      }
      for (i = 0; i < k * l; i++) {
	A[i] = randomnumber (dummy);
      }
      alpha = randomnumber (dummy);
      beta = randomnumber (dummy);
      t1 = gettime ();
      Rgemv (&trans, m, n, alpha, A, m, x, (mpackint) 1, beta, y, (mpackint) 1);
      t2 = gettime ();
      elapsedtime = (t2 - t1);
      printf ("     m       n      MFLOPS  trans\n");
      printf ("%6d  %6d  %10.3f      %c\n", (int) n, (int) m, (2.0 * (double) n * (double) m) / elapsedtime * MFLOPS, trans);

    }
    delete[]yd;
    delete[]y;
    delete[]x;
    delete[]A;
    n = n + STEPN;
    m = m + STEPM;
  }
  if (check_flag)
    dlclose (handle);
}
