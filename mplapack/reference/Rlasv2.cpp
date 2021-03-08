/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasv2.cpp,v 1.6 2010/08/07 04:48:33 nakatamaho Exp $ 
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

#define MTRUE 1
#define MFALSE 0

void Rlasv2(REAL f, REAL g, REAL h, REAL * ssmin, REAL * ssmax, REAL * snr, REAL * csr, REAL * snl, REAL * csl)
{
    REAL a, d, l, m, r, s, t, fa, ga, ha, ft, gt, ht, mm, tt, clt, crt, slt, srt;
    INTEGER pmax;
    REAL temp;
    INTEGER swap;
    REAL tsign = 0.0;
    INTEGER gasmal;
    REAL One = 1.0, Two = 2.0, Zero = 0.0, Half = 0.5, Four = 4.0;
    REAL mtemp1;

    ft = f;
    fa = abs(ft);
    ht = h;
    ha = abs(h);

//PMAX points to the maximum absolute element of matrix
//PMAX = 1 if F largest in absolute values
//PMAX = 2 if G largest in absolute values
//PMAX = 3 if H largest in absolute values

    pmax = 1;
    swap = ha > fa;
    if (swap) {
	pmax = 3;
	temp = ft;
	ft = ht;
	ht = temp;
	temp = fa;
	fa = ha;
	ha = temp;
//Now FA .ge. HA
    }
    gt = g;
    ga = abs(gt);
    if (ga == Zero) {

//Diagonal matrix 
	*ssmin = ha;
	*ssmax = fa;
	clt = One;
	crt = One;
	slt = Zero;
	srt = Zero;
    } else {
	gasmal = MTRUE;
	if (ga > fa) {
	    pmax = 2;
	    if (fa / ga < Rlamch("EPS")) {
//Case of very large GA
		gasmal = MFALSE;
		*ssmax = ga;
		if (ha > One) {
		    *ssmin = fa / (ga / ha);
		} else {
		    *ssmin = fa / ga * ha;
		}
		clt = One;
		slt = ht / gt;
		srt = One;
		crt = ft / gt;
	    }
	}
	if (gasmal) {
//Normal case
	    d = fa - ha;
	    if (d == fa) {
//Copes with infinite F or H
		l = One;
	    } else {
		l = d / fa;
	    }
//Note that 0 .le. L .le. 1
	    m = gt / ft;
//Note that abs(M) .le. 1/macheps
	    t = Two - l;
//Note that T .ge. 1
	    mm = m * m;
	    tt = t * t;
	    s = sqrt(tt + mm);
//Note that 1 .le. S .le. 1 + 1/macheps
	    if (l == Zero) {
		r = abs(m);
	    } else {
		r = sqrt(l * l + mm);
	    }
//Note that 0 .le. R .le. 1 + 1/macheps
	    a = (s + r) * Half;
//Note that 1 .le. A .le. 1 + abs(M)
	    *ssmin = ha / a;
	    *ssmax = fa * a;
	    if (mm == Zero) {
//Note that M is very tiny
		if (l == Zero) {
		    t = sign(Two, ft) * sign(One, gt);
		} else {
		    t = gt / sign(d, ft) + m / t;
		}
	    } else {
		t = (m / (s + t) + m / (r + l)) * (a + One);
	    }
	    l = sqrt(t * t + Four);
	    crt = Two / l;
	    srt = t / l;
	    clt = (crt + srt * m) / a;
	    slt = ht / ft * srt / a;
	}
    }
    if (swap) {
	*csl = srt;
	*snl = crt;
	*csr = slt;
	*snr = clt;
    } else {
	*csl = clt;
	*snl = slt;
	*csr = crt;
	*snr = srt;
    }
//Correct signs of SSMAX and SSMIN
    if (pmax == 1) {
	tsign = sign(One, *csr) * sign(One, *csl) * sign(One, f);
    }
    if (pmax == 2) {
	tsign = sign(One, *snr) * sign(One, *csl) * sign(One, g);
    }
    if (pmax == 3) {
	tsign = sign(One, *snr) * sign(One, *snl) * sign(One, h);
    }
    *ssmax = sign(*ssmax, tsign);
    mtemp1 = tsign * sign(One, f) * sign(One, h);
    *ssmin = sign(*ssmin, mtemp1);
    return;
}
