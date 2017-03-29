/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clartg.cpp,v 1.7 2010/08/07 04:48:32 nakatamaho Exp $ 
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
#include <iostream>

REAL abs1(COMPLEX ff)
{
    return max(abs(ff.real()), abs(ff.imag()));
}

REAL abssq(COMPLEX ff)
{
    REAL temp;

    temp = (ff.real() * ff.real()) + (ff.imag() * ff.imag());
    return temp;
}

void Clartg(COMPLEX f, COMPLEX g, REAL * cs, COMPLEX * sn, COMPLEX * r)
{
    REAL d;
    INTEGER i;
    REAL f2, g2;
    COMPLEX ff;
    REAL di, dr;
    COMPLEX fs, gs;
    REAL f2s, g2s, eps, scale;
    INTEGER count;
    REAL safmn2;
    REAL safmx2;
    REAL safmin;
    REAL Zero = 0.0, One = 1.0;

    safmin = Rlamch("S");
    eps = Rlamch("E");
    safmn2 = sqrt(safmin / eps);
    safmx2 = One / safmn2;
    scale = max(abs1(f), abs1(g));

    fs = f;
    gs = g;
    count = 0;
    if (scale >= safmx2) {
      L10:
	count++;
	fs = fs * safmn2;
	gs = gs * safmn2;
	scale = scale * safmn2;
	if (scale >= safmx2) {
	    goto L10;
	}
    } else if (scale <= safmn2) {
	if (g == Zero) {
	    *cs = One;
	    *sn = Zero;
	    *r = f;
	    return;
	}
      L20:
	--count;
	fs = fs * safmx2;
	gs = gs * safmx2;
	scale = scale * safmx2;
	if (scale <= safmn2) {
	    goto L20;
	}
    }
    f2 = abssq(fs);
    g2 = abssq(gs);
    if (f2 <= max(g2, One) * safmin) {
//This is a rare case: F is very small.
	if (f == Zero) {
	    *cs = Zero;
	    *r = Rlapy2(g.real(), g.imag());
//Do complex/real division explicitly with two real divisions
	    d = Rlapy2(gs.real(), gs.imag());
	    (*sn) = Real2Complex(gs.real() / d, -gs.imag() / d);
	    return;
	}
	std::cout << "# XXX Clartg not very well tested 1\n";
	f2s = Rlapy2(fs.real(), fs.imag());
//G2 and G2S are accurate
//G2 is at least SAFMIN, and G2S is at least SAFMN2
	g2s = sqrt(g2);
//Error in CS from underflow in F2S is at most
//UNFL / SAFMN2 .lt. sqrt(UNFL*EPS) .lt. EPS 
//If MAX(G2,ONE)=G2, then F2 .lt. G2*SAFMIN,
//and so CS .lt. sqrt(SAFMIN)
//If MAX(G2,ONE)=ONE, then F2 .lt. SAFMIN
//and so CS .lt. sqrt(SAFMIN)/SAFMN2 = sqrt(EPS)
//Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S
	*cs = f2s / g2s;
//Make sure abs(FF) = 1
//Do complex/real division explicitly with 2 real divisions
	if (abs1(f) > One) {
	    d = Rlapy2(f.real(), f.imag());
	    ff = Real2Complex(f.real() / d, f.imag() / d);
	} else {
	    dr = safmx2 * f.real();
	    di = safmx2 * f.imag();
	    d = Rlapy2(dr, di);
	    ff = Real2Complex(dr / d, di / d);
	}
	(*sn) = ff * Real2Complex(gs.real() / g2s, -gs.imag() / g2s);
	*r = *cs * f + *sn * g;
    } else {
//This is the most common case.
//Neither F2 nor F2/G2 are less than SAFMIN
//F2S cannot overflow, and it is accurate
	f2s = sqrt(g2 / f2 + One);
//Do the F2S(real)fS(complex) multiply with two real multiplies
	*r = Real2Complex(f2s * fs.real(), f2s * fs.imag());
	*cs = One / f2s;
	d = f2 + g2;
//Do complex/real division explicitly with two real divisions
	*sn = Real2Complex(r->real() / d, r->imag() / d);
	*sn = *sn * conj(gs);
	if (count != 0) {
	    if (count > 0) {
		for (i = 0; i < count; i++) {
		    *r = *r * safmx2;
		}
	    } else {
		std::cout << "# XXX Clartg not very well tested 2\n";
		for (i = 0; i < -count; i++) {
		    *r = *r * safmn2;
		}
	    }
	}
    }
    return;
}
