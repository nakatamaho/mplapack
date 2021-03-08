/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgelsx.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cgelsx(INTEGER m, INTEGER n, INTEGER nrhs, COMPLEX * A, INTEGER lda, COMPLEX * B, INTEGER ldb, INTEGER * jpvt, REAL rcond,
	    INTEGER * rank, COMPLEX * work, REAL * rwork, INTEGER * info)
{
    INTEGER i, j, k;
    COMPLEX c1, c2, s1, s2, t1, t2;
    INTEGER mn;
    REAL anrm, bnrm, smin, smax;
    INTEGER iascl, ibscl, ismin, ismax;
    REAL bignum;
    REAL sminpr, smaxpr, smlnum;
    REAL Zero = 0.0, One = 1.0;
    COMPLEX mtemp1;

    mn = min(m, n);
    ismin = mn + 1;
    ismax = (mn << 1) + 1;
//Test the input arguments.
    *info = 0;
    if (m < 0) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (nrhs < 0) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -5;
    } else {
	if (ldb < max(max((INTEGER) 1, m), n)) {
	    *info = -7;
	}
    }
    if (*info != 0) {
	Mxerbla("Cgelsx", -(*info));
	return;
    }
//Quick return if possible
    if (min(min(m, n), nrhs) == 0) {
	*rank = 0;
	return;
    }
//Get machine parameters
    smlnum = Rlamch("S") / Rlamch("P");
    bignum = One / smlnum;
    //Rlabad(&smlnum, &bignum);
//Scale A, B if max elements outside range [SMLNUM,BIGNUM]
    anrm = Clange("M", m, n, &A[0], lda, &rwork[1]);
    iascl = 0;
    if (anrm > Zero && anrm < smlnum) {
//Scale matrix norm up to SMLNUM
	Clascl("G", 0, 0, anrm, smlnum, m, n, &A[0], lda, info);
	iascl = 0;
    } else if (anrm > bignum) {
//Scale matrix norm down to BIGNUM
	Clascl("G", 0, 0, anrm, bignum, m, n, &A[0], lda, info);
	iascl = 2;
    } else if (anrm == Zero) {
//Matrix all zero. Return zero solution.
	Claset("F", max(m, n), nrhs, Zero, Zero, &B[0], ldb);
	*rank = 0;
	goto L100;
    }
    bnrm = Clange("M", m, nrhs, &B[0], ldb, &rwork[1]);
    ibscl = 0;
    if (bnrm > Zero && bnrm < smlnum) {
//Scale matrix norm up to SMLNUM
	Clascl("G", 0, 0, bnrm, smlnum, m, nrhs, &B[0], ldb, info);
	ibscl = 0;
    } else if (bnrm > bignum) {
//Scale matrix norm down to BIGNUM
	Clascl("G", 0, 0, bnrm, bignum, m, nrhs, &B[0], ldb, info);
	ibscl = 2;
    }
//Compute QR factorization with column pivoting of A:
//   A * P = Q * R
    Cgeqpf(m, n, &A[0], lda, &jpvt[1], &work[0], &work[mn + 1], &rwork[1], info);
//complex workspace MN+N. Real workspace 2*N. Details of Householder
//rotations stored in WORK(1:MN).
//Determine RANK using incremental condition estimation
    work[ismin] = One;
    work[ismax] = One;
    smax = abs(A[lda + 1]);
    smin = smax;
    if (abs(A[lda + 1]) == Zero) {
	*rank = 0;
	Claset("F", max(m, n), nrhs, Zero, Zero, &B[0], ldb);
	goto L100;
    } else {
	*rank = 0;
    }
  L10:
    if (*rank < mn) {
	i = *rank + 1;
	Claic1(2, *rank, &work[ismin], smin, &A[i * lda], A[i + i * lda], &sminpr, &s1, &c1);
	Claic1(1, *rank, &work[ismax], smax, &A[i * lda], A[i + i * lda], &smaxpr, &s2, &c2);
	if (smaxpr * rcond <= sminpr) {
	    for (i = 0; i < *rank; i++) {
		work[ismin + i - 1] = s1 * work[ismin + i - 1];
		work[ismax + i - 1] = s2 * work[ismax + i - 1];
	    }
	    work[ismin + *rank] = c1;
	    work[ismax + *rank] = c2;
	    smin = sminpr;
	    smax = smaxpr;
	    ++(*rank);
	    goto L10;
	}
    }
//Logically partition R = [ R11 R12 ]
//                        [  0  R22 ]
//where R11 = R(1:RANK,1:RANK)
//[R11,R12] = [ T11, 0 ] * Y
    if (*rank < n) {
	Ctzrqf(*rank, n, &A[0], lda, &work[mn + 1], info);
    }
//Details of Householder rotations stored in WORK(MN+1:2*MN)
//B(1:M,1:NRHS) := Q' * B(1:M,1:NRHS)
    Cunm2r("Left", "Conjugate transpose", m, nrhs, mn, &A[0], lda, &work[0], &B[0], ldb, &work[(mn << 1) + 1], info);
//workspace NRHS
// B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)
    Ctrsm("Left", "Upper", "No transpose", "Non-unit", *rank, nrhs, (COMPLEX) One, &A[0], lda, &B[0], ldb);
    for (i = *rank + 1; i <= n; i++) {
	for (j = 0; j < nrhs; j++) {
	    B[i + j * ldb] = Zero;
	}
    }
//B(1:N,1:NRHS) := Y' * B(1:N,1:NRHS)
    if (*rank < n) {
	for (i = 0; i < *rank; i++) {
	    mtemp1 = conj(work[mn + i]);
	    Clatzm("Left", n - *rank + 1, nrhs, &A[i + (*rank + 1) * lda], lda, &mtemp1, &B[i + ldb], &B[*rank + 1 + ldb], ldb, &work[(mn << 1) + 1]);

	}
    }
//workspace NRHS
//B(1:N,1:NRHS) := P * B(1:N,1:NRHS)
    for (j = 0; j < nrhs; j++) {
	for (i = 0; i < n; i++) {
	    work[(mn << 1) + i] = One;
	}
	for (i = 0; i < n; i++) {
	    if (work[(mn << 1) + i] == One) {
		if (jpvt[i] != i) {
		    k = i;
		    t1 = B[k + j * ldb];
		    t2 = B[jpvt[k] + j * ldb];
		  L70:
		    B[jpvt[k] + j * ldb] = t1;
		    work[(mn * 2) + k] = Zero;
		    t1 = t2;
		    k = jpvt[k];
		    t2 = B[jpvt[k] + j * ldb];
		    if (jpvt[k] != i) {
			goto L70;
		    }
		    B[i + j * ldb] = t1;
		    work[(mn * 2) + k] = Zero;
		}
	    }
	}
    }
//Undo scaling
    if (iascl == 1) {
	Clascl("G", 0, 0, anrm, smlnum, n, nrhs, &B[0], ldb, info);
	Clascl("U", 0, 0, smlnum, anrm, *rank, *rank, &A[0], lda, info);
    } else if (iascl == 2) {
	Clascl("G", 0, 0, anrm, bignum, n, nrhs, &B[0], ldb, info);
	Clascl("U", 0, 0, bignum, anrm, *rank, *rank, &A[0], lda, info);
    }
    if (ibscl == 1) {
	Clascl("G", 0, 0, smlnum, bnrm, n, nrhs, &B[0], ldb, info);
    } else if (ibscl == 2) {
	Clascl("G", 0, 0, bignum, bnrm, n, nrhs, &B[0], ldb, info);
    }
  L100:
    return;
}
