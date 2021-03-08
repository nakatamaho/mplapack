/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgecon.cpp,v 1.10 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cgecon(const char *norm, INTEGER n, COMPLEX * A, INTEGER lda, REAL anorm, REAL * rcond, COMPLEX * work, REAL * rwork, INTEGER * info)
{
    REAL sl;
    INTEGER ix;
    REAL su;
    INTEGER kase, kase1;
    REAL scale;
    INTEGER isave[3];
    REAL ainvnm;
    INTEGER onenrm;
    char normin;
    REAL smlnum;
    REAL One = 1.0, Zero = 0.0;

    *info = 0;
    onenrm = Mlsame(norm, "1") || Mlsame(norm, "O");
    if (!onenrm && !Mlsame(norm, "I")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -4;
    } else if (anorm < Zero) {
	*info = -5;
    }
    if (*info != 0) {
	Mxerbla("Cgecon", -(*info));
	return;
    }
/*     Quick return if possible */
    *rcond = Zero;
    if (n == 0) {
	*rcond = One;
	return;
    } else if (anorm == Zero) {
	return;
    }
    smlnum = Rlamch("Safe minimum");
//Estimate the norm of inv(A).
    ainvnm = Zero;
    normin = 'N';
    if (onenrm) {
	kase1 = 1;
    } else {
	kase1 = 2;
    }
    kase = 0;
  L10:
    Clacn2(n, &work[n + 1], &work[0], &ainvnm, &kase, isave);
    if (kase != 0) {
	if (kase == kase1) {
//Multiply by inv(L).
	    Clatrs("Lower", "No transpose", "Unit", &normin, n, &A[0], lda, &work[0], &sl, &rwork[1], info);
//Multiply by inv(U).
	    Clatrs("Upper", "No transpose", "Non-unit", &normin, n, &A[0], lda, &work[0], &su, &rwork[n + 1], info);
	} else {
//Multiply by inv(U').
	    Clatrs("Upper", "Conjugate transpose", "Non-unit", &normin, n, &A[0], lda, &work[0], &su, &rwork[n + 1], info);
//Multiply by inv(L').
	    Clatrs("Lower", "Conjugate transpose", "Unit", &normin, n, &A[0], lda, &work[0], &sl, &rwork[1], info);
	}
//Divide X by 1/(SL*SU) if doing so will not cause overflow.
	scale = sl * su;
	normin = 'Y';
	if (scale != One) {
	    ix = iCamax(n, &work[0], 1);
	    if (scale < Cabs1(work[ix]) * smlnum || scale == Zero) {
		goto L20;
	    }
	    CRrscl(n, scale, &work[0], 1);
	}
	goto L10;
    }
//Compute the estimate of the reciprocal condition number.
    if (ainvnm != Zero) {
	*rcond = One / ainvnm / anorm;
    }
  L20:
    return;
}
