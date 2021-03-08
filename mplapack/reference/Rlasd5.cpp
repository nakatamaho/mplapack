/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasd5.cpp,v 1.4 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rlasd5(INTEGER i, REAL * d, REAL * z, REAL * delta, REAL rho, REAL * dsigma, REAL * work)
{
    REAL b, c, w, del, tau, delsq;
    REAL Zero = 0.0, One = 1.0, Two = 2.0, Three = 3.0, Four = 4.0;

    del = d[2] - d[1];
    delsq = del * (d[2] + d[1]);
    if (i == 1) {
	w = rho * Four * (z[2] * z[2] / (d[1] + d[2] * Three) - z[1] * z[1] / (d[1] * Three + d[2])) / del + One;
	if (w > Zero) {
	    b = delsq + rho * (z[1] * z[1] + z[2] * z[2]);
	    c = rho * z[1] * z[1] * delsq;
//B > ZERO, always
//The following TAU is DSIGMA * DSIGMA - D( 1 ) * D( 1 )
	    tau = c * Two / (b + sqrt(abs(b * b - c * Four)));

//The following TAU is DSIGMA - D( 1 )
	    tau /= d[1] + sqrt(d[1] * d[1] + tau);
	    *dsigma = d[1] + tau;
	    delta[1] = -tau;
	    delta[2] = del - tau;
	    work[1] = d[1] * Two + tau;
	    work[2] = d[1] + tau + d[2];
//DELTA( 1 ) = -Z( 1 ) / TAU
//DELTA( 2 ) = Z( 2 ) / ( DEL-TAU )
	} else {
	    b = -delsq + rho * (z[1] * z[1] + z[2] * z[2]);
	    c = rho * z[2] * z[2] * delsq;
//The following TAU is DSIGMA * DSIGMA - D( 2 ) * D( 2 )
	    if (b > Zero) {
		tau = c * -Two / (b + sqrt(b * b + c * Four));
	    } else {
		tau = (b - sqrt(b * b + c * Four)) / Two;
	    }

//The following TAU is DSIGMA - D( 2 ) */
	    tau /= d[2] + sqrt(abs(d[2] * d[2] + tau));
	    *dsigma = d[2] + tau;
	    delta[1] = -(del + tau);
	    delta[2] = -tau;
	    work[1] = d[1] + tau + d[2];
	    work[2] = d[2] * Two + tau;
//DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
//DELTA( 2 ) = -Z( 2 ) / TAU
	}
//TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
//DELTA( 1 ) = DELTA( 1 ) / TEMP
//DELTA( 2 ) = DELTA( 2 ) / TEMP
    } else {

//Now I=2
	b = -delsq + rho * (z[1] * z[1] + z[2] * z[2]);
	c = rho * z[2] * z[2] * delsq;
//The following TAU is DSIGMA * DSIGMA - D( 2 ) * D( 2 )
	if (b > Zero) {
	    tau = (b + sqrt(b * b + c * Four)) / Two;
	} else {
	    tau = c * Two / (-b + sqrt(b * b + c * Four));
	}

//The following TAU is DSIGMA - D( 2 )
	tau /= d[2] + sqrt(d[2] * d[2] + tau);
	*dsigma = d[2] + tau;
	delta[1] = -(del + tau);
	delta[2] = -tau;
	work[1] = d[1] + tau + d[2];
	work[2] = d[2] * Two + tau;
//DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
//DELTA( 2 ) = -Z( 2 ) / TAU
//TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
//DELTA( 1 ) = DELTA( 1 ) / TEMP
//DELTA( 2 ) = DELTA( 2 ) / TEMP
    }
    return;
}
