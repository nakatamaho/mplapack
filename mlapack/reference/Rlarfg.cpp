/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlarfg.cpp,v 1.7 2010/08/07 04:48:32 nakatamaho Exp $ 
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
#include <stdio.h>		//for debugging
void Rlarfg(INTEGER N, REAL * alpha, REAL * x, INTEGER incx, REAL * tau)
{
    REAL xnorm;
    REAL One = 1.0;
    REAL beta;
    REAL safmin;
    REAL rsafmn;
    INTEGER knt;

    if (N <= 1) {
	*tau = 0.0;
	return;
    }
    xnorm = Rnrm2(N - 1, x, incx);
//H  =  I
    if (xnorm == 0.0) {
	*tau = 0.0;
    } else {
	beta = -1.0 * sign(Rlapy2(*alpha, xnorm), *alpha);
	safmin = Rlamch("S") / Rlamch("E");

//XNORM, BETA may be inaccurate; scale X and recompute them
	if (abs(beta) < safmin) {
	    fprintf(stderr, "# Rlarfg: 1: XXX not very well tested\n");
	    rsafmn = One / safmin;
	    knt = 0;
	    while (abs(beta) < safmin) {
		knt++;
		Rscal(N - 1, rsafmn, x, incx);
		beta = beta * rsafmn;
		*alpha = *alpha * rsafmn;
	    }

//New BETA is at most 1, at least SAFMIN
	    xnorm = Rnrm2(N - 1, x, incx);
	    beta = -1.0 * sign(Rlapy2(*alpha, xnorm), *alpha);
	    *tau = (beta - *alpha) / beta;
	    Rscal(N - 1, One / (*alpha - beta), x, incx);

//If ALPHA is subnormal, it may lose relative accuracy
	    *alpha = beta;
	    for (INTEGER j = 0; j < knt; j++) {
		*alpha = *alpha * safmin;
	    }
	} else {
	    *tau = (beta - *alpha) / beta;
	    Rscal(N - 1, One / (*alpha - beta), x, incx);
	    *alpha = beta;
	}
    }
}
