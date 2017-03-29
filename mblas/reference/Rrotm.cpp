/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rrotm.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $
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
 *
 * $Id: Rrotm.cpp,v 1.5 2010/08/07 05:50:10 nakatamaho Exp $

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

/*
Based on http://www.netlib.org/blas/drotm.f
Apply the modified givens transformation.
*/

#include <mblas.h>

void Rrotm(INTEGER n, REAL * dx, INTEGER incx, REAL * dy, INTEGER incy, REAL * dparam)
{
    INTEGER i, kx, ky;
    REAL Zero = 0.0, Two = 2.0;
    REAL w, z;
    REAL mflag = dparam[0];
    REAL dh11, dh12, dh21, dh22;

    if (n <= 0 || mflag == -Two)
	return;

    kx = 0;
    ky = 0;
    if (incx < 0)
	kx = (1 - n) * incx;
    if (incy < 0)
	ky = (1 - n) * incy;

    if (mflag == Zero) {
	dh12 = dparam[3];
	dh21 = dparam[2];
	for (i = 0; i < n; i++) {
	    w = dx[kx];
	    z = dy[ky];
	    dx[kx] = w + z * dh12;
	    dy[ky] = w * dh21 + z;
	    kx = kx + incx;
	    ky = ky + incy;
	}
    } else if (mflag > Zero) {
	dh11 = dparam[1];
	dh22 = dparam[4];
	for (i = 0; i < n; i++) {
	    w = dx[kx];
	    z = dy[ky];
	    dx[kx] = w * dh11 + z;
	    dy[ky] = -w + dh22 * z;
	    kx = kx + incx;
	    ky = ky + incy;
	}
    } else if (mflag < Zero) {
	dh11 = dparam[1];
	dh12 = dparam[3];
	dh21 = dparam[2];
	dh22 = dparam[4];
	for (i = 0; i < n; i++) {
	    w = dx[kx];
	    z = dy[ky];
	    dx[kx] = w * dh11 + z * dh12;
	    dy[ky] = w * dh21 + z * dh22;
	    kx = kx + incx;
	    ky = ky + incy;
	}
    }
}
