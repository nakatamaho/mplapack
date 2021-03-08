/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlacn2.cpp,v 1.11 2010/08/13 00:14:22 nakatamaho Exp $ 
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
/*
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

$COPYRIGHT$

Additional copyrights may follow

$HEADER$

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer. 
  
- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.
  
- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.
  
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
*/

#include <mpblas.h>
#include <mplapack.h>
#include <stdio.h>

#define ITMAX 5
void Rlancn2_finalization(INTEGER * kase, INTEGER * isave, REAL * x, INTEGER n)
{
    INTEGER i;
    REAL altsgn;
    REAL One = 1.0;
    double tmp;
    altsgn = One;
    for (i = 0; i < n; i++) {
	tmp = (double) i / (double) (n - 1);
	x[i] = altsgn * (One + tmp);
	altsgn = -altsgn;
    }
    *kase = 1;
    isave[0] = 5;
    return;
}

void Rlacn2(INTEGER n, REAL * v, REAL * x, INTEGER * isgn, REAL * est, INTEGER * kase, INTEGER * isave)
{
    INTEGER i, j;
    INTEGER jlast;
    double dtmp;
    REAL temp;
    REAL estold;
    REAL Zero = 0.0, One = 1.0, Two = 2.0;

//initialization
    if (*kase == 0) {
	for (i = 0; i < n; i++) {
	    x[i] = One / n;
	}
	*kase = 1;
	isave[0] = 1;
	return;
    }
    switch (isave[0]) {
    case 1:
//Entry isave(0)=1. 
//First iteration. X has been overwritten by A*X.
	if (n == 1) {
	    v[0] = x[0];
	    *est = abs(v[0]);
	    kase = 0;
	    return;
	}
	*est = Rasum(n, x, 1);
	for (i = 0; i < n; i++) {
	    x[i] = sign(One, x[i]);
	    isgn[i] = (INTEGER) nint(x[i]);
	}
	*kase = 2;
	isave[0] = 2;
	return;
    case 2:
//Entry isave(0)=2. 
//First iteration. X has been overwritten by A^t*X.
	isave[1] = iRamax(n, &x[0], 1);
	isave[2] = 2;
	for (i = 0; i < n; i++) {
	    x[i] = Zero;
	}
	x[isave[1] - 1] = One;
	*kase = 1;
	isave[0] = 3;
	return;
    case 3:
//Entry isave(0)=3. 
//First iteration. X has been overwritten by A*X.
	Rcopy(n, &x[0], 1, &v[0], 1);
	estold = *est;
	*est = Rasum(n, &v[0], 1);
	for (i = 0; i < n; i++) {
	    if ((INTEGER) nint(sign(One, x[i])) != isgn[i]) {
//Test for cycling.
		if (*est <= estold) {
		    Rlancn2_finalization(kase, isave, x, n);
		    return;
		} else {
		    for (j = 0; j < n; j++) {
			x[j] = sign(One, x[j]);
			isgn[j] = (INTEGER) nint(x[j]);
		    }
		    *kase = 2;
		    isave[0] = 4;
		    return;
		}
	    }
	}
//Repeated sign vector detected, hence algorithm has converged.
	Rlancn2_finalization(kase, isave, x, n);
	return;

    case 4:
//Entry isave(0)=4. 
//First iteration. X has been overwritten by A^t*X.
	jlast = isave[1];
	isave[1] = iRamax(n, &x[0], 1);
	if (x[jlast - 1] != abs(x[isave[1] - 1]) && isave[2] < ITMAX) {
	    isave[2] = isave[2] + 1;
	    for (i = 0; i < n; i++) {
		x[i] = Zero;
	    }
	    x[isave[1] - 1] = One;
	    *kase = 1;
	    isave[0] = 3;
	    return;
	}
	Rlancn2_finalization(kase, isave, x, n);
	return;

    case 5:
//Entry isave(0)=5. 
//X has been overwritten by A*X.
	//QD/DD only dirty
	dtmp = n * 3;
	temp = Two * (Rasum(n, &x[0], 1) / (REAL) dtmp);
	if (temp > *est) {
	    Rcopy(n, &x[0], 1, &v[0], 1);
	    *est = temp;
	}
	*kase = 0;
	return;
    }
    return;
}
