/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaln2.cpp,v 1.7 2010/08/07 04:48:32 nakatamaho Exp $ 
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

#define MFALSE 0
#define MTRUE 1

void
Rlaln2(INTEGER ltrans, INTEGER na, INTEGER nw, REAL smin, REAL ca, REAL * A,
       INTEGER lda, REAL d1, REAL d2, REAL * B, INTEGER ldb, REAL wr, REAL wi, REAL * x, INTEGER ldx, REAL * scale, REAL * xnorm, INTEGER * info)
{
    INTEGER zswap[4] = { MFALSE, MFALSE, MTRUE, MTRUE };
    INTEGER rswap[4] = { MFALSE, MTRUE, MFALSE, MTRUE };
    INTEGER ipivot[16] = { 1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, 3, 2, 1 };
    INTEGER j;
    REAL bi1, bi2, br1, br2, xi1, xi2, xr1, xr2, ci21, ci22, cr21, cr22, li21, csi, ui11, lr21, ui12, ui22;
    REAL csr, ur11, ur12, ur22;
    REAL bbnd, cmax, ui11r, ui12s, temp, ur11r, ur12s, u22abs;
    INTEGER icmax;
    REAL bnorm, cnorm, smini;
    REAL bignum, smlnum;
    REAL equiv_0[4], equiv_1[4];
    REAL mtemp1, mtemp2;
    REAL Zero = 0.0, One = 1.0, Two = 2.0;

#define ci (equiv_0)
#define cr (equiv_1)
#define civ (equiv_0)
#define crv (equiv_1)

//Compute BIGNUM
    smlnum = Two * Rlamch("Safe minimum");
    bignum = One / smlnum;
    smini = max(smin, smlnum);
//Don't check for input errors
    *info = 0;
//Standard Initializations
    *scale = One;
    if (na == 1) {
//1 x 1  (i.e., scalar) system   C X = B
	if (nw == 1) {
//Real 1x1 system.
//C = ca A - w D
	    csr = ca * A[lda + 1] - wr * d1;
	    cnorm = abs(csr);
//If | C | < SMINI, use C = SMINI
	    if (cnorm < smini) {
		csr = smini;
		cnorm = smini;
		*info = 1;
	    }
//Check scaling for  X = B / C
	    bnorm = abs(B[ldb + 1]);
	    if (cnorm < One && bnorm > One) {
		if (bnorm > bignum * cnorm) {
		    *scale = One / bnorm;
		}
	    }
//Compute X
	    x[ldx + 1] = B[ldb + 1] * *scale / csr;
	    *xnorm = abs(x[ldx + 1]);
	} else {
//Complex 1x1 system (w is complex)
//C = ca A - w D
	    csr = ca * A[lda + 1] - wr * d1;
	    csi = -(wi) * d1;
	    cnorm = abs(csr) + abs(csi);
//If | C | < SMINI, use C = SMINI
	    if (cnorm < smini) {
		csr = smini;
		csi = Zero;
		cnorm = smini;
		*info = 1;
	    }
//Check scaling for  X = B / C
	    bnorm = abs(B[ldb + 1]) + abs(B[(ldb * 2) + 1]);
	    if (cnorm < One && bnorm > One) {
		if (bnorm > bignum * cnorm) {
		    *scale = One / bnorm;
		}
	    }
//Compute X
	    mtemp1 = *scale * B[ldb + 1];
	    mtemp2 = *scale * B[(ldb * 2) + 1];
	    Rladiv(mtemp1, mtemp2, csr, csi, &x[ldx + 1], &x[(ldx * 2) + 1]);
	    *xnorm = abs(x[ldx + 1]) + abs(x[(ldx * 2) + 1]);
	}

    } else {
//2x2 System
//Compute the real part of  C = ca A - w D  (or  ca A' - w D )
	cr[0] = ca * A[lda + 1] - wr * d1;
	cr[3] = ca * A[(lda * 2) + 2] - wr * d2;
	if (ltrans) {
	    cr[2] = ca * A[lda + 2];
	    cr[1] = ca * A[(lda * 2) + 1];
	} else {
	    cr[1] = ca * A[lda + 2];
	    cr[2] = ca * A[(lda * 2) + 1];
	}
	if (nw == 1) {
//Real 2x2 system  (w is real)
//Find the largest element in C
	    cmax = Zero;
	    icmax = 0;
	    for (j = 0; j < 4; j++) {
		if (abs(crv[j - 1]) > cmax) {
		    cmax = abs(crv[j - 1]);
		    icmax = j;
		}
	    }
//If norm(C) < SMINI, use SMINI*identity.
	    if (cmax < smini) {
		mtemp1 = abs(B[ldb + 1]), mtemp2 = abs(B[ldb + 2]);
		bnorm = max(mtemp1, mtemp2);
		if (smini < One && bnorm > One) {
		    if (bnorm > bignum * smini) {
			*scale = One / bnorm;
		    }
		}
		temp = *scale / smini;
		x[ldx + 1] = temp * B[ldb + 1];
		x[ldx + 2] = temp * B[ldb + 2];
		*xnorm = temp * bnorm;
		*info = 1;
		return;
	    }
//Gaussian elimination with complete pivoting.
	    ur11 = crv[icmax - 1];
	    cr21 = crv[ipivot[(icmax * 4) - 3] - 1];
	    ur12 = crv[ipivot[(icmax * 4) - 2] - 1];
	    cr22 = crv[ipivot[(icmax * 4) - 1] - 1];
	    ur11r = One / ur11;
	    lr21 = ur11r * cr21;
	    ur22 = cr22 - ur12 * lr21;
//If smaller pivot < SMINI, use SMINI
	    if (abs(ur22) < smini) {
		ur22 = smini;
		*info = 1;
	    }
	    if (rswap[icmax - 1]) {
		br1 = B[ldb + 2];
		br2 = B[ldb + 1];
	    } else {
		br1 = B[ldb + 1];
		br2 = B[ldb + 2];
	    }
	    br2 = br2 - lr21 * br1;
	    mtemp1 = abs(br1 * (ur22 * ur11r)), mtemp2 = abs(br2);
	    bbnd = max(mtemp1, mtemp2);
	    if (bbnd > One && abs(ur22) < One) {
		if (bbnd >= bignum * abs(ur22)) {
		    *scale = One / bbnd;
		}
	    }
	    xr2 = br2 * *scale / ur22;
	    xr1 = *scale * br1 * ur11r - xr2 * (ur11r * ur12);
	    if (zswap[icmax - 1]) {
		x[ldx + 1] = xr2;
		x[ldx + 2] = xr1;
	    } else {
		x[ldx + 1] = xr1;
		x[ldx + 2] = xr2;
	    }
	    mtemp1 = abs(xr1), mtemp2 = abs(xr2);
	    *xnorm = max(mtemp1, mtemp2);
//Further scaling if  norm(A) norm(X) > overflow
	    if (*xnorm > One && cmax > One) {
		if (*xnorm > bignum / cmax) {
		    temp = cmax / bignum;
		    x[ldx + 1] = temp * x[ldx + 1];
		    x[ldx + 2] = temp * x[ldx + 2];
		    *xnorm = temp * *xnorm;
		    *scale = temp * *scale;
		}
	    }
	} else {
//Complex 2x2 system  (w is complex)
//Find the largest element in C
	    ci[0] = -(wi) * d1;
	    ci[1] = Zero;
	    ci[2] = Zero;
	    ci[3] = -(wi) * d2;
	    cmax = Zero;
	    icmax = 0;
	    for (j = 0; j < 4; j++) {
		if (abs(crv[j - 1]) + abs(civ[j - 1]) > cmax) {
		    cmax = abs(crv[j - 1]) + abs(civ[j - 1]);
		    icmax = j;
		}
	    }
//If norm(C) < SMINI, use SMINI*identity.
	    if (cmax < smini) {
		mtemp1 = abs(B[ldb + 1]) + abs(B[(ldb * 2) + 1]);
		mtemp2 = abs(B[ldb + 2]) + abs(B[(ldb * 2) + 2]);
		bnorm = max(mtemp1, mtemp2);
		if (smini < One && bnorm > One) {
		    if (bnorm > bignum * smini) {
			*scale = One / bnorm;
		    }
		}
		temp = *scale / smini;
		x[ldx + 1] = temp * B[ldb + 1];
		x[ldx + 2] = temp * B[ldb + 2];
		x[(ldx * 2) + 1] = temp * B[(ldb * 2) + 1];
		x[(ldx * 2) + 2] = temp * B[(ldb * 2) + 2];
		*xnorm = temp * bnorm;
		*info = 1;
		return;
	    }
//Gaussian elimination with complete pivoting.
	    ur11 = crv[icmax - 1];
	    ui11 = civ[icmax - 1];
	    cr21 = crv[ipivot[(icmax * 4) - 3] - 1];
	    ci21 = civ[ipivot[(icmax * 4) - 3] - 1];
	    ur12 = crv[ipivot[(icmax * 4) - 2] - 1];
	    ui12 = civ[ipivot[(icmax * 4) - 2] - 1];
	    cr22 = crv[ipivot[(icmax * 4) - 1] - 1];
	    ci22 = civ[ipivot[(icmax * 4) - 1] - 1];
	    if (icmax == 1 || icmax == 4) {
//Code when off-diagonals of pivoted C are real
		if (abs(ur11) > abs(ui11)) {
		    temp = ui11 / ur11;
		    ur11r = One / (ur11 * (temp * temp + One));
		    ui11r = -temp * ur11r;
		} else {
		    temp = ur11 / ui11;
		    ui11r = -One / (ui11 * (temp * temp + One));
		    ur11r = -temp * ui11r;
		}
		lr21 = cr21 * ur11r;
		li21 = cr21 * ui11r;
		ur12s = ur12 * ur11r;
		ui12s = ur12 * ui11r;
		ur22 = cr22 - ur12 * lr21;
		ui22 = ci22 - ur12 * li21;
	    } else {
//Code when diagonals of pivoted C are real
		ur11r = One / ur11;
		ui11r = Zero;
		lr21 = cr21 * ur11r;
		li21 = ci21 * ur11r;
		ur12s = ur12 * ur11r;
		ui12s = ui12 * ur11r;
		ur22 = cr22 - ur12 * lr21 + ui12 * li21;
		ui22 = -ur12 * li21 - ui12 * lr21;
	    }
	    u22abs = abs(ur22) + abs(ui22);
//If smaller pivot < SMINI, use SMINI
	    if (u22abs < smini) {
		ur22 = smini;
		ui22 = Zero;
		*info = 1;
	    }
	    if (rswap[icmax - 1]) {
		br2 = B[ldb + 1];
		br1 = B[ldb + 2];
		bi2 = B[(ldb * 2) + 1];
		bi1 = B[(ldb * 2) + 2];
	    } else {
		br1 = B[ldb + 1];
		br2 = B[ldb + 2];
		bi1 = B[(ldb * 2) + 1];
		bi2 = B[(ldb * 2) + 2];
	    }
	    br2 = br2 - lr21 * br1 + li21 * bi1;
	    bi2 = bi2 - li21 * br1 - lr21 * bi1;
	    mtemp1 = (abs(br1) + abs(bi1)) * (u22abs * (abs(ur11r) + abs(ui11r)));
	    mtemp2 = abs(br2) + abs(bi2);
	    bbnd = max(mtemp1, mtemp2);
	    if (bbnd > One && u22abs < One) {
		if (bbnd >= bignum * u22abs) {
		    *scale = One / bbnd;
		    br1 = *scale * br1;
		    bi1 = *scale * bi1;
		    br2 = *scale * br2;
		    bi2 = *scale * bi2;
		}
	    }

	    Rladiv(br2, bi2, ur22, ui22, &xr2, &xi2);
	    xr1 = ur11r * br1 - ui11r * bi1 - ur12s * xr2 + ui12s * xi2;
	    xi1 = ui11r * br1 + ur11r * bi1 - ui12s * xr2 - ur12s * xi2;
	    if (zswap[icmax - 1]) {
		x[ldx + 1] = xr2;
		x[ldx + 2] = xr1;
		x[(ldx * 2) + 1] = xi2;
		x[(ldx * 2) + 2] = xi1;
	    } else {
		x[ldx + 1] = xr1;
		x[ldx + 2] = xr2;
		x[(ldx * 2) + 1] = xi1;
		x[(ldx * 2) + 2] = xi2;
	    }
	    mtemp1 = abs(xr1) + abs(xi1), mtemp2 = abs(xr2) + abs(xi2);
	    *xnorm = max(mtemp1, mtemp2);
//Further scaling if  norm(A) norm(X) > overflow
	    if (*xnorm > One && cmax > One) {
		if (*xnorm > bignum / cmax) {
		    temp = cmax / bignum;
		    x[ldx + 1] = temp * x[ldx + 1];
		    x[ldx + 2] = temp * x[ldx + 2];
		    x[(ldx * 2) + 1] = temp * x[(ldx * 2) + 1];
		    x[(ldx * 2) + 2] = temp * x[(ldx * 2) + 2];
		    *xnorm = temp * *xnorm;
		    *scale = temp * *scale;
		}
	    }
	}
    }
    return;
}

#undef crv
#undef civ
#undef cr
#undef ci
