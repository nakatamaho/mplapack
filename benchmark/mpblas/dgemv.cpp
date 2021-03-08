/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: dgemv.cpp,v 1.4 2010/08/07 05:50:09 nakatamaho Exp $
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

#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define F77_FUNC(name,NAME) name ## _
#include <blas.h>
#define ___DOUBLE_BENCH___
#include <mplapack_benchmark.h>

#define TOTALSTEPS 1000

int
main (int argc, char *argv[])
{
  int k, l, m, n;
  int incx = 1, incy = 1, STEPN, STEPM, N0, M0;
  double alpha, beta, dummy, *dummywork;
  double mOne = -1;
  double elapsedtime, t1, t2;
  int i, p;
  char trans, normtype;

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
    }
  }
  n = N0;
  m = M0;
  for (p = 0; p < TOTALSTEPS; p++) {
    if (lsame_f77 (&trans, "n")) {
      k = n;
      l = m;
    }
    else {
      k = m;
      l = n;
    }
    double *x = new double[k];
    double *y = new double[l];
    double *A = new double[n * m];
    for (i = 0; i < k; i++) {
      x[i] = randomnumber (dummy);
    }
    for (i = 0; i < l; i++) {
      y[i] = randomnumber (dummy);
    }
    for (i = 0; i < k * l; i++) {
      A[i] = randomnumber (dummy);
    }
    alpha = randomnumber (dummy);
    beta = randomnumber (dummy);
    t1 = gettime ();
    dgemv_f77 (&trans, &m, &n, &alpha, A, &m, x, &incx, &beta, y, &incy);
    t2 = gettime ();
    elapsedtime = (t2 - t1);
    printf ("     m       n      MFLOPS  trans\n");
    printf ("%6d  %6d  %10.3f      %c\n", (int) n, (int) m, (2.0 * (double) n * (double) m) / elapsedtime * MFLOPS, trans);
    delete[]y;
    delete[]x;
    delete[]A;
    n = n + STEPN;
    m = m + STEPM;
  }
}
