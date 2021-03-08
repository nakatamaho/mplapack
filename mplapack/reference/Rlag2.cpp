/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlag2.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rlag2(REAL * A, INTEGER lda, REAL * B, INTEGER ldb, REAL safmin, REAL * scale1, REAL * scale2, REAL * wr1, REAL * wr2, REAL * wi)
{
    REAL r, c1, c2, c3, c4, c5, s1, s2, a11, a12, a21, a22,
	b11, b12, b22, pp, qq, ss, as11, as12, as22, sum, abi22, diff,
	bmin, wbig, wabs, wdet, binv11, binv22, discr, anorm, bnorm, bsize, shift, rtmin, rtmax, wsize, ascale, bscale, wscale, safmax, wsmall;
    REAL mtemp0, mtemp1, mtemp2, mtemp3, mtemp4;
    REAL Zero = 0.0, Half = 0.5, One = 1.0, Fuzzy1 = 1.0 + 0.00001;

    rtmin = sqrt(safmin);
    rtmax = One / rtmin;
    safmax = One / safmin;
//Scale A
    mtemp1 = abs(A[lda + 1]) + abs(A[lda + 2]);
    mtemp2 = abs(A[(lda << 1) + 1]) + abs(A[(lda << 1) + 2]);
    mtemp3 = max(mtemp1, mtemp2);
    anorm = max(mtemp3, safmin);
    ascale = One / anorm;
    a11 = ascale * A[lda + 1];
    a21 = ascale * A[lda + 2];
    a12 = ascale * A[(lda << 1) + 1];
    a22 = ascale * A[(lda << 1) + 2];
//Perturb B if necessary to insure non-singularity
    b11 = B[ldb + 1];
    b12 = B[(ldb << 1) + 1];
    b22 = B[(ldb << 1) + 2];

    mtemp1 = abs(b11), mtemp2 = abs(b12);
    mtemp3 = max(mtemp1, mtemp2);
    mtemp1 = abs(b22);
    mtemp2 = max(mtemp1, mtemp3);
    bmin = rtmin * max(mtemp2, rtmin);
    if (abs(b11) < bmin) {
	b11 = sign(bmin, b11);
    }
    if (abs(b22) < bmin) {
	b22 = sign(bmin, b22);
    }
//Scale B
    mtemp1 = abs(b11), mtemp2 = abs(b12) + abs(b22);
    mtemp3 = max(mtemp1, mtemp2);
    bnorm = max(mtemp3, safmin);
    mtemp1 = abs(b11), mtemp2 = abs(b22);
    bsize = max(mtemp1, mtemp2);
    bscale = One / bsize;
    b11 *= bscale;
    b12 *= bscale;
    b22 *= bscale;

//Compute larger eigenvalue by method described by C. van Loan
//( AS is A shifted by -SHIFT*B )
    binv11 = One / b11;
    binv22 = One / b22;
    s1 = a11 * binv11;
    s2 = a22 * binv22;
    if (abs(s1) <= abs(s2)) {
	as12 = a12 - s1 * b12;
	as22 = a22 - s1 * b22;
	ss = a21 * (binv11 * binv22);
	abi22 = as22 * binv22 - ss * b12;
	pp = abi22 * Half;
	shift = s1;
    } else {
	as12 = a12 - s2 * b12;
	as11 = a11 - s2 * b11;
	ss = a21 * (binv11 * binv22);
	abi22 = -ss * b12;
	pp = (as11 * binv11 + abi22) * Half;
	shift = s2;
    }
    qq = ss * as12;
    if (abs(pp * rtmin) >= One) {
	discr = (rtmin * pp) * (rtmin * pp) + qq * safmin;
	r = sqrt((abs(discr))) * rtmax;
    } else {
	if (pp * pp + abs(qq) <= safmin) {
	    discr = (rtmax * pp) * (rtmax * pp) + qq * safmax;
	    r = sqrt((abs(discr))) * rtmin;
	} else {
	    discr = pp * pp + qq;
	    r = sqrt((abs(discr)));
	}
    }
//Note: the test of R in the following IF is to cover the case when
//      DISCR is small and negative and is flushed to zero during
//      the calculation of R.  On machines which have a consistent
//      flush-to-zero threshhold and handle numbers above that
//      threshhold correctly, it would not be necessary.
    if (discr >= Zero || r == Zero) {
	sum = pp + sign(r, pp);
	diff = pp - sign(r, pp);
	wbig = shift + sum;
//Compute smaller eigenvalue
	wsmall = shift + diff;
	mtemp2 = abs(wsmall);
	mtemp1 = max(mtemp2, safmin);
	if (abs(wbig) * Half > mtemp1) {
	    wdet = (a11 * a22 - a12 * a21) * (binv11 * binv22);
	    wsmall = wdet / wbig;
	}
//Choose (real) eigenvalue closest to 2,2 element of A*B**(-1)
//for WR1
	if (pp > abi22) {
	    *wr1 = min(wbig, wsmall);
	    *wr2 = max(wbig, wsmall);
	} else {
	    *wr1 = max(wbig, wsmall);
	    *wr2 = min(wbig, wsmall);
	}
	*wi = Zero;
    } else {
//Complex eigenvalues
	*wr1 = shift + pp;
	*wr2 = *wr1;
	*wi = r;
    }

//Further scaling to avoid underflow and overflow in computing */
//SCALE1 and overflow in computing w*B.
//This scale factor (WSCALE) is bounded from above using C1 and C2,
//and from below using C3 and CFour
//   C1 implements the condition  s A  must never overflow.
//   C2 implements the condition  w B  must never overflow.
//   C3, with C2,
//      implement the condition that s A - w B must never overflow.
//   C4 implements the condition  s    should not underflow.
//   C5 implements the condition  max(s,|w|) should be at least Two

    c1 = bsize * (safmin * max(One, ascale));
    c2 = safmin * max(One, bnorm);
    c3 = bsize * safmin;
    if (ascale <= One && bsize <= One) {
	mtemp1 = One, mtemp2 = ascale / safmin * bsize;
	c4 = min(mtemp1, mtemp2);
    } else {
	c4 = One;
    }
    if (ascale <= One || bsize <= One) {
	mtemp1 = One, mtemp2 = ascale * bsize;
	c5 = min(mtemp1, mtemp2);
    } else {
	c5 = One;
    }
//Scale first eigenvalue
    wabs = abs(*wr1) + abs(*wi);
    mtemp1 = c4, mtemp2 = max(wabs, c5) * Half;
    mtemp3 = min(mtemp1, mtemp2);
    mtemp1 = max(safmin, c1), mtemp2 = (wabs * c2 + c3) * Fuzzy1;
    mtemp4 = max(mtemp1, mtemp2);
    wsize = max(mtemp3, mtemp4);
    if (wsize != One) {
	wscale = One / wsize;
	if (wsize > One) {
	    *scale1 = max(ascale, bsize) * wscale * min(ascale, bsize);
	} else {
	    *scale1 = min(ascale, bsize) * wscale * max(ascale, bsize);
	}
	*wr1 *= wscale;
	if (*wi != Zero) {
	    *wi *= wscale;
	    *wr2 = *wr1;
	    *scale2 = *scale1;
	}
    } else {
	*scale1 = ascale * bsize;
	*scale2 = *scale1;
    }
//Scale second eigenvalue (if real)
    if (*wi == Zero) {
	mtemp0 = abs(*wr2);
	mtemp1 = c4, mtemp2 = max(mtemp0, c5) * Half;
	mtemp3 = min(mtemp1, mtemp2);

	mtemp1 = max(safmin, c1), mtemp2 = (abs(*wr2) * c2 + c3) * Fuzzy1;
	mtemp4 = max(mtemp1, mtemp2);

	wsize = max(mtemp3, mtemp4);
	if (wsize != One) {
	    wscale = One / wsize;
	    if (wsize > One) {
		*scale2 = max(ascale, bsize) * wscale * min(ascale, bsize);
	    } else {
		*scale2 = min(ascale, bsize) * wscale * max(ascale, bsize);
	    }
	    *wr2 *= wscale;
	} else {
	    *scale2 = ascale * bsize;
	}
    }
    return;
}
