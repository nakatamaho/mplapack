/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: ddot.cpp,v 1.4 2010/08/07 05:50:09 nakatamaho Exp $
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
#include <mpack_benchmark.h>

#define TOTALSTEPS 1000

int
main (int argc, char *argv[])
{
  int n;
  int incx = 1, incy = 1, STEP, N0;
  double dummy, ans, ans_ref;
  double mOne = -1;
  double elapsedtime, t1, t2;
  int i, p;
  int check_flag = 1;

  // initialization
  N0 = 1;
  STEP = 1;
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

  n = N0;
  for (p = 0; p < TOTALSTEPS; p++) {
    double *x = new double[n];
    double *y = new double[n];
    for (i = 0; i < n; i++) {
      x[i] = randomnumber (dummy);
      y[i] = randomnumber (dummy);
    }
    t1 = gettime ();
    ans = ddot_f77 (&n, x, &incx, y, &incy);
    t2 = gettime ();
    elapsedtime = (t2 - t1);
    printf ("         n       MFLOPS\n");
    printf ("%10d   %10.3f\n", (int) n, (2.0 * (double) n) / elapsedtime * MFLOPS);
    delete[]y;
    delete[]x;
    n = n + STEP;
  }
}
