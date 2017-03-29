/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgttrf.cpp,v 1.7 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rgttrf(INTEGER n, REAL * dl, REAL * d, REAL * du, REAL * du2, INTEGER * ipiv, INTEGER * info)
{
    INTEGER i;
    REAL fact, temp;
    REAL Zero = 0.0;

    *info = 0;
    if (n < 0) {
	*info = -1;
	Mxerbla("Rgttrf", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0)
	return;
//Initialize IPIV(i) = i and DU2(I) = 0
    for (i = 0; i < n; i++) {
	ipiv[i] = i;

    }
    for (i = 0; i < n - 2; i++) {
	du2[i] = Zero;
    }
    for (i = 0; i < n - 2; i++) {
	if (abs(d[i]) >= abs(dl[i])) {
//No row interchange required, eliminate DL(I)
	    if (d[i] != Zero) {
		fact = dl[i] / d[i];
		dl[i] = fact;
		d[i + 1] -= fact * du[i];
	    }
	} else {
//Interchange rows I and I+1, eliminate DL(I)
	    fact = d[i] / dl[i];
	    d[i] = dl[i];
	    dl[i] = fact;
	    temp = du[i];
	    du[i] = d[i + 1];
	    d[i + 1] = temp - fact * d[i + 1];
	    du2[i] = du[i + 1];
	    du[i + 1] = -fact * du[i + 1];
	    ipiv[i] = i + 1;
	}

    }
    if (n > 1) {
	i = n - 1;
	if (abs(d[i]) >= abs(dl[i])) {
	    if (d[i] != Zero) {
		fact = dl[i] / d[i];
		dl[i] = fact;
		d[i + 1] -= fact * du[i];
	    }
	} else {
	    fact = d[i] / dl[i];
	    d[i] = dl[i];
	    dl[i] = fact;
	    temp = du[i];
	    du[i] = d[i + 1];
	    d[i + 1] = temp - fact * d[i + 1];
	    ipiv[i] = i + 1;
	}
    }
//Check for a zero on the diagonal of U.
    for (i = 0; i < n; i++) {
	if (d[i] == Zero) {
	    *info = i;
	    return;
	}
    }
    return;
}
