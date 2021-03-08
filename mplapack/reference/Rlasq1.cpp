/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasq1.cpp,v 1.7 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rlasq1(INTEGER n, REAL * d, REAL * e, REAL * work, INTEGER * info)
{
    INTEGER i;
    REAL eps;
    REAL scale;
    INTEGER iinfo;
    REAL sigmn;
    REAL sigmx;
    REAL safmin;
    REAL Zero = 0.0;
    REAL mtemp1, mtemp2;

    *info = 0;
    if (n < 0) {
	*info = -2;
	Mxerbla("Rlasq1", -(*info));
	return;
    } else if (n == 0) {
	return;
    } else if (n == 1) {
	d[1] = abs(d[1]);
	return;
    } else if (n == 2) {
	Rlas2(d[0], e[0], d[2], &sigmn, &sigmx);
	d[1] = sigmx;
	d[2] = sigmn;
	return;
    }
//Estimate the largest singular value.
    sigmx = Zero;
    for (i = 0; i < n - 1; i++) {
	d[i] = abs(d[i]);
	mtemp1 = sigmx;
	mtemp2 = abs(e[i]);
	sigmx = max(mtemp1, mtemp2);

    }
    d[n] = abs(d[n]);
//Early return if SIGMX is zero (matrix is already diagonal).
    if (sigmx == Zero) {
	Rlasrt("D", n, &d[0], &iinfo);
	return;
    }
    for (i = 0; i < n; i++) {
	mtemp1 = sigmx;
	mtemp2 = d[i];
	sigmx = max(mtemp1, mtemp2);
    }
//Copy D and E into WORK (in the Z format) and scale (squaring the
//input data makes scaling by a power of the radix pointless).
    eps = Rlamch("P");
    safmin = Rlamch("S");
    scale = sqrt(eps / safmin);
    Rcopy(n, &d[0], 1, &work[0], 2);
    Rcopy(n - 1, &e[0], 1, &work[2], 2);
    Rlascl("G", 0, 0, sigmx, scale, n * 2 - 1, 1, &work[0], n * 2 - 1, &iinfo);

//Compute the q's and e's.
    for (i = 0; i < n * 2 - 1; i++) {
	work[i] = work[i] * work[i];
    }
    work[n * 2] = Zero;
    Rlasq2(n, &work[0], info);
    if (*info == 0) {
	for (i = 0; i < n; i++) {
	    d[i] = sqrt(work[i]);

	}
	Rlascl("G", 0, 0, scale, sigmx, n, 1, &d[0], n, &iinfo);
    }
    return;
}
