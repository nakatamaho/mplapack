/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgeequ.cpp,v 1.9 2010/08/11 08:48:00 nakatamaho Exp $ 
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

void Rgeequ(INTEGER m, INTEGER n, REAL * A, INTEGER lda, REAL * r, REAL * c, REAL * rowcnd, REAL * colcnd, REAL * amax, INTEGER * info)
{
    INTEGER i, j;
    REAL rcmin, rcmax;
    REAL bignum, smlnum;
    REAL One = 1.0, Zero = 0.0;
    REAL mtemp1, mtemp2;

//Test the input parameters.
    *info = 0;
    if (m < 0) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -4;
    }
    if (*info != 0) {
	Mxerbla("Rgeequ", -(*info));
	return;
    }
//Quick return if possible
    if (m == 0 || n == 0) {
	*rowcnd = One;
	*colcnd = One;
	*amax = Zero;
	return;
    }
//Get machine constants.
    smlnum = Rlamch("S");
    bignum = One / smlnum;

//Compute row scale factors.
    for (i = 0; i < m; i++) {
	r[i] = Zero;
    }
//Find the maximum element in each row.
    for (j = 0; j < n; j++) {
	for (i = 0; i < m; i++) {
	    mtemp1 = r[i], mtemp2 = abs(A[i + j * lda]);
	    r[i] = max(mtemp1, mtemp2);
	}
    }
//Find the maximum and minimum scale factors.
    rcmin = bignum;
    rcmax = Zero;
    for (i = 0; i < m; i++) {
	mtemp1 = rcmax, mtemp2 = r[i];
	rcmax = max(mtemp1, mtemp2);
	mtemp1 = rcmin, mtemp2 = r[i];
	rcmin = min(mtemp1, mtemp2);
    }
    *amax = rcmax;

    if (rcmin == Zero) {
//Find the first zero scale factor and return an error code.
	for (i = 0; i < m; i++) {
	    if (r[i] == Zero) {
		*info = i;
		return;
	    }
	}
    } else {
//Invert the scale factors.
	for (i = 0; i < m; i++) {
	    mtemp2 = r[i];
	    mtemp1 = max(mtemp2, smlnum);
	    r[i] = One / min(mtemp1, bignum);
	}
//Compute ROWCND = min(R(I)) / max(R(I))
	*rowcnd = max(rcmin, smlnum) / min(rcmax, bignum);
    }
//Compute column scale factors
    for (j = 0; j < n; j++) {
	c[j] = Zero;
    }
//Find the maximum element in each column,
//assuming the row scaling computed above.
    for (j = 0; j < n; j++) {
	for (i = 0; i < m; i++) {
	    mtemp1 = c[j], mtemp2 = abs(A[i + j * lda]) * r[i];
	    c[j] = max(mtemp1, mtemp2);
	}
    }
//Find the maximum and minimum scale factors.
    rcmin = bignum;
    rcmax = Zero;
    for (j = 0; j < n; j++) {
	mtemp1 = rcmin, mtemp2 = c[j];
	rcmin = min(mtemp1, mtemp2);
	mtemp1 = rcmax, mtemp2 = c[j];
	rcmax = max(mtemp1, mtemp2);
    }
    if (rcmin == Zero) {
//Find the first zero scale factor and return an error code.
	for (j = 0; j < n; j++) {
	    if (c[j] == Zero) {
		*info = m + j;
		return;
	    }
	}
    } else {
//Invert the scale factors.
	for (j = 0; j < n; j++) {
	    mtemp1 = max(c[j], smlnum);
	    c[j] = One / min(mtemp1, bignum);
	}
//Compute COLCND = min(C(J)) / max(C(J))
	*colcnd = max(rcmin, smlnum) / min(rcmax, bignum);
    }
    return;
}
