/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rsptrd.cpp,v 1.8 2010/08/07 04:48:33 nakatamaho Exp $ 
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

#include <mblas.h>
#include <mlapack.h>

void Rsptrd(const char *uplo, INTEGER n, REAL * ap, REAL * d, REAL * e, REAL * tau, INTEGER * info)
{
    INTEGER i, i1, ii, i1i1;
    INTEGER upper;
    REAL taui;
    REAL alpha;
    REAL Half = 0.5, Zero = 0.0, One = 1.0;

    *info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	Mxerbla("Rsptrd", -(*info));
	return;
    }
//Quick return if possible
    if (n <= 0)
	return;

    if (upper) {

//Reduce the upper triangle of A.
//I1 is the index in AP of A(1,I+1).
	i1 = n * (n - 1) / 2 + 1;
	for (i = n - 1; i >= 1; i--) {

//Generate elementary reflector H(i) = I - tau * v * v'
//to annihilate A(1:i-1,i+1)

	    Rlarfg(i, &ap[i1 + i - 1], &ap[i1], 1, &taui);
	    e[i] = ap[i1 + i - 1];

	    if (taui != Zero) {
//Apply H(i) from both sides to A(1:i,1:i)
		ap[i1 + i - 1] = One;

//Compute  y := tau * A * v  storing y in TAU(1:i)
		Rspmv(uplo, i, taui, &ap[0], &ap[i1], 1, Zero, &tau[0], 1);

//Compute  w := y - 1/2 * tau * (y'*v) * v
		alpha = taui * -Half * Rdot(i, &tau[0], 1, &ap[i1], 1);
		Raxpy(i, alpha, &ap[i1], 1, &tau[0], 1);

//Apply the transformation as a rank-2 update:
//A := A - v * w' - w * v'
		Rspr2(uplo, i, -One, &ap[i1], 1, &tau[0], 1, &ap[0]);
		ap[i1 + i - 1] = e[i];
	    }
	    d[i + 1] = ap[i1 + i];
	    tau[i] = taui;
	    i1 -= i;
	}
	d[0] = ap[0];
    } else {
//Reduce the lower triangle of A. II is the index in AP of
//A(i,i) and I1I1 is the index of A(i+1,i+1).
	ii = 0;
	i1 = n - 1;
	for (i = 0; i < i1; i++) {
	    i1i1 = ii + n - i + 1;

//Generate elementary reflector H(i) = I - tau * v * v'
//to annihilate A(i+2:n,i)
	    Rlarfg(n - i, &ap[ii + 1], &ap[ii + 2], 1, &taui);
	    e[i] = ap[ii + 1];

	    if (taui != Zero) {
//Apply H(i) from both sides to A(i+1:n,i+1:n)
		ap[ii + 1] = One;

//Compute  y := tau * A * v  storing y in TAU(i:n-1)
		Rspmv(uplo, n - i, taui, &ap[i1i1], &ap[ii + 1], 1, Zero, &tau[i], 1);

//Compute  w := y - 1/2 * tau * (y'*v) * v
		alpha = taui * -Half * Rdot(n - i, &tau[i], 1, &ap[ii + 1], 1);
		Raxpy(n - i, alpha, &ap[ii + 1], 1, &tau[i], 1);

//Apply the transformation as a rank-2 update:
//A := A - v * w' - w * v'
		Rspr2(uplo, n - i, -One, &ap[ii + 1], 1, &tau[i], 1, &ap[i1i1]);

		ap[ii + 1] = e[i];
	    }
	    d[i] = ap[ii];
	    tau[i] = taui;
	    ii = i1i1;

	}
	d[n] = ap[ii];
    }
    return;
}
