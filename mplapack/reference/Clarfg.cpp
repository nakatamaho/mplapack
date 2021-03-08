/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clarfg.cpp,v 1.8 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Clarfg(INTEGER n, COMPLEX * alpha, COMPLEX * x, INTEGER incx, COMPLEX * tau)
{
    INTEGER j, knt;
    REAL beta, alphi, alphr;
    REAL xnorm;
    REAL safmin;
    REAL rsafmn;
    REAL One = 1.0, Zero = 0.0;

    if (n <= 0) {
	*tau = Zero;
	return;
    }
    xnorm = RCnrm2(n - 1, &x[0], incx);
    alphr = alpha->real();
    alphi = alpha->imag();
    if (xnorm == Zero && alphi == Zero) {
//H  =  I
	*tau = Zero;
    } else {
//general case
	beta = -sign(Rlapy3(alphr, alphi, xnorm), alphr);
	safmin = Rlamch("S") / Rlamch("E");
	rsafmn = One / safmin;
	if (abs(beta) < safmin) {
//XNORM, BETA may be inaccurate; scale X and recompute them
	    knt = 0;
	    while (1) {
		knt++;
		CRscal(n - 1, rsafmn, &x[0], incx);
		beta = beta * rsafmn;
		alphi = alphi * rsafmn;
		alphr = alphr * rsafmn;
		if (abs(beta) < safmin)
		    continue;
		break;
	    }
//New BETA is at most 1, at least SAFMIN
	    xnorm = RCnrm2(n - 1, &x[0], incx);
	    *alpha = alphr;
	    beta = -sign(Rlapy3(alphr, alphi, xnorm), alphr);
	    *tau = COMPLEX((beta - alphr) / beta, -alphi / beta);
	    *alpha = Cladiv(One, *alpha - beta);
	    Cscal(n - 1, *alpha, &x[0], incx);
//If ALPHA is subnormal, it may lose relative accuracy
	    *alpha = beta;
	    for (j = 0; j < knt; j++) {
		*alpha = *alpha * safmin;
	    }
	} else {
	    *tau = COMPLEX ((beta - alphr) / beta, -alphi / beta);
	    *alpha = Cladiv((COMPLEX) One, *alpha - beta);
	    Cscal(n - 1, *alpha, &x[0], incx);
	    *alpha = beta;
	}
    }
    return;
}
