/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clacn2.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Clacn2(INTEGER n, COMPLEX * v, COMPLEX * x, REAL * est, INTEGER * kase, INTEGER * isave)
{
    INTEGER i;
    REAL temp, absxi;
    INTEGER jlast;
    REAL safmin, altsgn, estold;
    REAL Zero = 0.0, One = 1.0, Two = 2.0;

    safmin = Rlamch("Safe minimum");
    if (*kase == 0) {
	for (i = 0; i < n; i++) {
	    x[i] = One / (REAL) double (n);
	}
	*kase = 1;
	isave[1] = 1;
	return;
    }

    switch (isave[1]) {
    case 1:
	goto L20;
    case 2:
	goto L40;
    case 3:
	goto L70;
    case 4:
	goto L90;
    case 5:
	goto L120;
    }
//................ ENTRY   (ISAVE( 1 ) = 1)
//FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
  L20:
    if (n == 1) {
	v[1] = x[0];
	*est = abs(v[1]);
//... QUIT
	goto L130;
    }
    *est = RCsum1(n, &x[0], 1);

    for (i = 0; i < n; i++) {
	absxi = abs(x[i]);
	if (absxi > safmin) {
	    x[i] = Real2Complex(x[i].real() / absxi, x[i].imag() / absxi);
	} else {
	    x[i] = One;
	}

    }
    *kase = 2;
    isave[1] = 2;
    return;
//................ ENTRY   (ISAVE( 1 ) = 2)
//FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
  L40:
    isave[2] = iCmax1(n, &x[0], 1);
    isave[3] = 2;
//MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
  L50:
    for (i = 0; i < n; i++) {
	x[i] = Zero;
    }
    x[isave[2]] = One;
    *kase = 1;
    isave[1] = 3;
    return;
//................ ENTRY   (ISAVE( 1 ) = 3)
//X HAS BEEN OVERWRITTEN BY A*X.
  L70:
    Ccopy(n, &x[0], 1, &v[1], 1);
    estold = *est;
    *est = RCsum1(n, &v[1], 1);
//TEST FOR CYCLING.
    if (*est <= estold) {
	goto L100;
    }
    for (i = 0; i < n; i++) {
	absxi = abs(x[i]);
	if (absxi > safmin) {
	    x[i] = Real2Complex(x[i].real() / absxi, x[i].imag() / absxi);
	} else {
	    x[i] = One;
	}
    }
    *kase = 2;
    isave[1] = 4;
    return;
//................ ENTRY   (ISAVE( 1 ) = 4)
//X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
  L90:
    jlast = isave[2];
    isave[2] = iCmax1(n, &x[0], 1);
    if (abs(x[jlast]) != abs(x[isave[2]]) && isave[3] < 5) {
	isave[3]++;
	goto L50;
    }
//ITERATION COMPLETE.  FINAL STAGE.
  L100:
    altsgn = One;
    for (i = 0; i < n; i++) {
	x[i] = (altsgn * (i - 1) / (n - 1) + One);
	altsgn = -altsgn;
    }
    *kase = 1;
    isave[1] = 5;
    return;
//................ ENTRY   (ISAVE( 1 ) = 5)
//X HAS BEEN OVERWRITTEN BY A*X.
  L120:
    temp = RCsum1(n, &x[0], 1) / (n * 3) * Two;
    if (temp > *est) {
	Ccopy(n, &x[0], 1, &v[1], 1);
	*est = temp;
    }
  L130:
    kase = 0;
    return;
}
