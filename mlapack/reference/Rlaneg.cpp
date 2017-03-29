/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaneg.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

INTEGER Rlaneg(INTEGER n, REAL * d, REAL * lld, REAL sigma, REAL pivmin, INTEGER r)
{
    INTEGER j;
    REAL p, t;
    INTEGER bj;
    REAL tmp;
    INTEGER neg1, neg2;
    REAL bsav, gamma, dplus;
    INTEGER negcnt;
    INTEGER sawnan;
    REAL dminus;
    REAL Zero = 0.0, One = 1.0;

    negcnt = 0;
//I) upper part: L D L^T - SIGMA I = L+ D+ L+^T
    t = -(sigma);
    for (bj = 0; bj <= r - 1; bj += 128) {
	neg1 = 0;
	bsav = t;
	for (j = bj; j <= min(bj + 127, r - 1); j++) {
	    dplus = d[j] + t;
	    if (dplus < Zero) {
		++neg1;
	    }
	    tmp = t / dplus;
	    t = tmp * lld[j] - sigma;
	}
	sawnan = Risnan(t);
//Run a slower version of the above loop if a NaN is detected.
//A NaN should occur only with a zero pivot after an infinite
//pivot.  In that case, substituting 1 for T/DPLUS is the
//correct limit.
	if (sawnan) {
	    neg1 = 0;
	    t = bsav;
	    for (j = bj; j <= min(bj + 127, r - 1); j++) {
		dplus = d[j] + t;
		if (dplus < Zero) {
		    ++neg1;
		}
		tmp = t / dplus;
		if (Risnan(tmp)) {
		    tmp = One;
		}
		t = tmp * lld[j] - sigma;
	    }
	}
	negcnt = negcnt + neg1;
    }
//II) lower part: L D L^T - SIGMA I = U- D- U-^T
    p = d[n] - sigma;
    for (bj = n - 1; bj >= r; bj += -128) {
	neg2 = 0;
	bsav = p;
	for (j = bj; j >= max(bj - 127, r); j--) {
	    dminus = lld[j] + p;
	    if (dminus < Zero) {
		++neg2;
	    }
	    tmp = p / dminus;
	    p = tmp * d[j] - sigma;
	}
	sawnan = Risnan(p);
//As above, run a slower version that substitutes 1 for Inf/Inf.
	if (sawnan) {
	    neg2 = 0;
	    p = bsav;
	    for (j = bj; j >= max(bj - 127, r); j--) {
		dminus = lld[j] + p;
		if (dminus < Zero) {
		    ++neg2;
		}
		tmp = p / dminus;
		if (Risnan(tmp)) {
		    tmp = One;
		}
		p = tmp * d[j] - sigma;
	    }
	}
	negcnt = negcnt + neg2;
    }
//III) Twist index
//  T was shifted by SIGMA initially.
    gamma = t + sigma + p;
    if (gamma < Zero) {
	++negcnt;
    }
    return negcnt;
}
