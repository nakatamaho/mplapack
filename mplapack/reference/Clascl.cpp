/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clascl.cpp,v 1.10 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Clascl(const char *type, INTEGER kl, INTEGER ku, REAL cfrom, REAL cto, INTEGER m, INTEGER n, COMPLEX * A, INTEGER lda, INTEGER * info)
{
    INTEGER i, j, k1, k2, k3, k4;
    INTEGER itype;
    REAL bignum, cfrom1, cfromc, cto1, ctoc, mul, smlnum;
    REAL Zero = 0.0, One = 1.0;
    INTEGER done = MFALSE;

    *info = 0;
    if (Mlsame(type, "G")) {
	itype = 0;
    } else if (Mlsame(type, "L")) {
	itype = 1;
    } else if (Mlsame(type, "U")) {
	itype = 2;
    } else if (Mlsame(type, "H")) {
	itype = 3;
    } else if (Mlsame(type, "B")) {
	itype = 4;
    } else if (Mlsame(type, "Q")) {
	itype = 5;
    } else if (Mlsame(type, "Z")) {
	itype = 6;
    } else {
	itype = -1;
    }
    if (itype == -1) {
	*info = -1;
    } else if (cfrom == Zero) {
	*info = -4;
    } else if (m < 0) {
	*info = -6;
    } else if (n < 0 || (itype == 4 && n != m) || (itype == 5 && n != m)) {
	*info = -7;
    } else if (itype <= 3 && lda < max((INTEGER) 1, m)) {
	*info = -9;
    } else if (itype >= 4) {
	if (kl < 0 || kl > max(m - 1, (INTEGER) 0)) {
	    *info = -2;
	} else {
	    if (ku < 0 || ku > max(n - 1, (INTEGER) 0) || ((itype == 4 || itype == 5) && kl != ku)) {
		*info = -3;
	    } else if ((itype == 4 && lda < kl + 1) || (itype == 5 && lda < ku + 1)
		       || (itype == 6 && lda < (kl * 2) + ku + 1)) {
		*info = -9;
	    }
	}
    }
    if (*info != 0) {
	Mxerbla("Clascl", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0 || m == 0) {
	return;
    }
//Get machine parameters
    smlnum = Rlamch("S");
    bignum = One / smlnum;

    cfromc = cfrom;
    ctoc = cto;

    while (done == MFALSE) {
	cfrom1 = cfromc * smlnum;
	cto1 = ctoc / bignum;
	if (abs(cfrom1) > abs(ctoc) && ctoc != Zero) {
	    mul = smlnum;
	    done = MFALSE;
	    cfromc = cfrom1;
	} else if (abs(cto1) > abs(cfromc)) {
	    mul = bignum;
	    done = MFALSE;
	    ctoc = cto1;
	} else {
	    mul = ctoc / cfromc;
	    done = MTRUE;
	}
	if (itype == 0) {
//Full matrix
	    for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
		    A[i + j * lda] = A[i + j * lda] * mul;
		}
	    }
	} else if (itype == 1) {
//Lower triangular matrix
	    for (j = 0; j < n; j++) {
		for (i = j; i < m; i++) {
		    A[i + j * lda] = A[i + j * lda] * mul;
		}
	    }
	} else if (itype == 2) {
//Upper triangular matrix
	    for (j = 0; j < n; j++) {
		for (i = 0; i <= min(j, m - 1); i++) {
		    A[i + j * lda] = A[i + j * lda] * mul;
		}
	    }
	} else if (itype == 3) {
//Upper Hessenberg matrix
	    for (j = 0; j < n; j++) {
		for (i = 0; i <= min(j + 1, m - 1); i++) {
		    A[i + j * lda] = A[i + j * lda] * mul;
		}
	    }
	} else if (itype == 4) {
//Lower half of a symmetric band matrix
	    k3 = kl + 1;
	    k4 = n + 1;
	    for (j = 0; j < n; j++) {
		for (i = 0; i < min(k3, k4 - j - 1); i++) {
		    A[i + j * lda] = A[i + j * lda] * mul;
		}
	    }
	} else if (itype == 5) {
//Upper half of a symmetric band matrix
	    k1 = ku + 2;
	    k3 = ku + 1;
	    for (j = 0; j < n; j++) {
		for (i = max(k1 - j - 1, (INTEGER) 1) - 1; i < k3; i++) {
		    A[i + j * lda] = A[i + j * lda] * mul;
		}
	    }
	} else if (itype == 6) {
//Band matrix
	    k1 = kl + ku + 2;
	    k2 = kl + 1;
	    k3 = kl * 2 + ku + 1;
	    k4 = kl + ku + 1 + m;
	    for (j = 0; j < n; j++) {
		for (i = max(k1 - j - 1, k2) - 1; i < min(k3, k4 - j - 1); i++) {
		    A[i + j * lda] = A[i + j * lda] * mul;
		}
	    }
	}
    }
    return;
}
