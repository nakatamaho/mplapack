/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rdot_dd.cpp,v 1.3 2010/08/07 05:50:09 nakatamaho Exp $
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

int main (int argc, char *argv[])
{
  mpackint n;
  mpackint incx = 1, incy = 1, STEP, N0;
  REAL dummy, ans, ans_ref;
  REAL mOne = -1;
  double elapsedtime, t1, t2;
  int i, p; 
  int check_flag = 1;
  char normtype;
  ___MPACK_INITIALIZE___

  const char mblas_sym[] = SYMBOL_GCC_RDOT;
  void *handle;
  REAL (*mblas_ref) (mpackint, REAL *, mpackint, REAL *, mpackint);
  char *error;
  REAL diff;
  double diffr;

  // initialization
  N0 = 1;
  STEP = 1;
  normtype = 'm';
  if (argc != 1) {
    for (i = 1; i < argc; i++) {
      if (strcmp ("-N", argv[i]) == 0) {
	N0 = atoi (argv[++i]);
      }
      else if (strcmp ("-STEP", argv[i]) == 0) {
	STEP = atoi (argv[++i]);
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
    mblas_ref = (REAL (*)(mpackint, REAL *, mpackint, REAL *, mpackint)) dlsym (handle, mblas_sym);
    if ((error = dlerror ()) != NULL) {
      fprintf (stderr, "%s\n", error);
      return 1;
    }
  }

  n = N0;
  for (p = 0; p < TOTALSTEPS; p++) {
    REAL *x = new REAL[n];
    REAL *y = new REAL[n];
    if (check_flag) {
      for (i = 0; i < n; i++) {
	x[i] = randomnumber (dummy);
	y[i] = randomnumber (dummy);
      }
      t1 = gettime ();
      ans = Rdot (n, x, incx, y, incy);
      t2 = gettime ();
      elapsedtime = (t2 - t1);
      ans_ref = (*mblas_ref) (n, x, incx, y, incy);
      diff = ans - ans_ref;
      diffr = cast2double (diff);
      printf ("         n       MFLOPS      error\n");
      printf ("%10d   %10.3f   %5.2e\n", (int) n, (2.0 * (double) n ) / elapsedtime * MFLOPS, diffr);
    }
    else {
      for (i = 0; i < n; i++) {
        x[i] = randomnumber (dummy);
        y[i] = randomnumber (dummy);
      }
      t1 = gettime ();
      ans = Rdot (n, x, incx, y, incy);
      t2 = gettime ();
      elapsedtime = (t2 - t1);
      printf ("         n       MFLOPS\n");
      printf ("%10d   %10.3f\n", (int) n, (2.0 * (double) n ) / elapsedtime * MFLOPS);
    }
    delete[]y;
    delete[]x;
    n = n + STEP;
  }
  if (check_flag)
    dlclose (handle);
}
