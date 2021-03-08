/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Ctgexc.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Ctgexc(LOGICAL wantq, LOGICAL wantz, INTEGER n, COMPLEX * A, INTEGER lda, COMPLEX * B, INTEGER ldb, COMPLEX * q, INTEGER ldq,
	    COMPLEX * z, INTEGER ldz, INTEGER ifst, INTEGER * ilst, INTEGER * info)
{
    INTEGER here;
//Decode and test input arguments.
    *info = 0;
    if (n < 0) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -5;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -7;
    } else if (ldq < 1 || (wantq && ldq < max((INTEGER) 1, n))) {
	*info = -9;
    } else if (ldz < 1 || (wantz && ldz < max((INTEGER) 1, n))) {
	*info = -11;
    } else if (ifst < 1 || ifst > n) {
	*info = -12;
    } else if (*ilst < 1 || *ilst > n) {
	*info = -13;
    }
    if (*info != 0) {
	Mxerbla("Ctgexc", -(*info));
	return;
    }
//Quick return if possible
    if (n <= 1) {
	return;
    }
    if (ifst == *ilst) {
	return;
    }
    if (ifst < *ilst) {
	here = ifst;
      L10:
//Swap with next one below
	Ctgex2(wantq, wantz, n, &A[0], lda, &B[0], ldb, &q[0], ldq, &z[0], ldz, here, info);
	if (*info != 0) {
	    *ilst = here;
	    return;
	}
	++here;
	if (here < *ilst) {
	    goto L10;
	}
	--here;
    } else {
	here = ifst - 1;
      L20:
//Swap with next one above
	Ctgex2(wantq, wantz, n, &A[0], lda, &B[0], ldb, &q[0], ldq, &z[0], ldz, here, info);
	if (*info != 0) {
	    *ilst = here;
	    return;
	}
	--here;
	if (here >= *ilst) {
	    goto L20;
	}
	++here;
    }
    *ilst = here;
    return;
}
