/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgelsy.cpp,v 1.3 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rgelsy(INTEGER m, INTEGER n, INTEGER nrhs, REAL * A, INTEGER lda, REAL * B, INTEGER ldb, INTEGER * jpvt, REAL rcond,
	    INTEGER * rank, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i, j;
    REAL c1, c2, s1, s2;
    INTEGER nb, mn, nb1, nb2, nb3, nb4;
    REAL anrm, bnrm, smin, smax;
    INTEGER iascl, ibscl;
    INTEGER ismin, ismax;
    REAL wsize;
    REAL bignum;
    INTEGER lwkmin;
    REAL sminpr, smaxpr, smlnum;
    INTEGER lwkopt;
    INTEGER lquery;
    REAL One = 1.0, Zero = 0.0;
    REAL mtemp1, mtemp2;

    mn = min(m, n);
    ismin = mn + 1;
    ismax = (mn * 2) + 1;
//Test the input arguments.
    *info = 0;
    lquery = lwork == -1;
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
//Figure out optimal block size
    if (*info == 0) {
	if (mn == 0 || nrhs == 0) {
	    lwkmin = 1;
	    lwkopt = 1;
	} else {
	    nb1 = iMlaenv(1, "Rgeqrf", " ", m, n, -1, -1);
	    nb2 = iMlaenv(1, "Rgerqf", " ", m, n, -1, -1);
	    nb3 = iMlaenv(1, "Rormqr", " ", m, n, nrhs, -1);
	    nb4 = iMlaenv(1, "Rormrq", " ", m, n, nrhs, -1);
	    nb = max(max(max(nb1, nb2), nb3), nb4);
	    lwkmin = mn + max(max(mn * 2, n + 1), mn + nrhs);
	    lwkopt = max(max(lwkmin, mn + (n * 2) + nb * (n + 1)), (mn << 1) + nb * nrhs);
	}
	work[1] = lwkopt;
	if (lwork < lwkmin && !lquery) {
	    *info = -12;
	}
    }
    if (*info != 0) {
	Mxerbla("Rgelsy", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (mn == 0 || nrhs == 0) {
	*rank = 0;
	return;
    }
//Get machine parameters
    smlnum = Rlamch("S") / Rlamch("P");
    bignum = One / smlnum;
    //Rlabad(&smlnum, &bignum);
//Scale A, B if max entries outside range [SMLNUM,BIGNUM]
    anrm = Rlange("M", m, n, &A[0], lda, &work[0]);
    iascl = 0;
    if (anrm > Zero && anrm < smlnum) {
//Scale matrix norm up to SMLNUM
	Rlascl("G", 0, 0, anrm, smlnum, m, n, &A[0], lda, info);
	iascl = 0;
    } else if (anrm > bignum) {
//Scale matrix norm down to BIGNUM
	Rlascl("G", 0, 0, anrm, bignum, m, n, &A[0], lda, info);
	iascl = 2;
    } else if (anrm == Zero) {
//Matrix all zero. Return zero solution.
	Rlaset("F", max(m, n), nrhs, Zero, Zero, &B[0], ldb);
	*rank = 0;
	goto L70;
    }
    bnrm = Rlange("M", m, nrhs, &B[0], ldb, &work[0]);
    ibscl = 0;
    if (bnrm > Zero && bnrm < smlnum) {
//Scale matrix norm up to SMLNUM
	Rlascl("G", 0, 0, bnrm, smlnum, m, nrhs, &B[0], ldb, info);
	ibscl = 0;
    } else if (bnrm > bignum) {
//Scale matrix norm down to BIGNUM
	Rlascl("G", 0, 0, bnrm, bignum, m, nrhs, &B[0], ldb, info);
	ibscl = 2;
    }
//Compute QR factorization with column pivoting of A:
//   A * P = Q * R
    Rgeqp3(m, n, &A[0], lda, &jpvt[1], &work[0], &work[mn + 1], lwork - mn, info);
    wsize = mn + work[mn + 1];
//workspace: MN+2*N+NB*(N+1).
//Details of Householder rotations stored in WORK(1:MN).
//Determine RANK using incremental condition estimation
    work[ismin] = One;
    work[ismax] = One;
    smax = abs(A[lda + 1]);
    smin = smax;
    if (abs(A[lda + 1]) == Zero) {
	*rank = 0;
	Rlaset("F", max(m, n), nrhs, Zero, Zero, &B[0], ldb);
	goto L70;
    } else {
	*rank = 0;
    }
  L10:
    if (*rank < mn) {
	i = *rank + 1;
	Rlaic1(2, *rank, &work[ismin], smin, &A[i * lda], A[i + i * lda], &sminpr, &s1, &c1);
	Rlaic1(1, *rank, &work[ismax], smax, &A[i * lda], A[i + i * lda], &smaxpr, &s2, &c2);
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
//workspace: 3*MN.
//Logically partition R = [ R11 R12 ]
//                        [  0  R22 ]
//where R11 = R(1:RANK,1:RANK)
//[R11,R12] = [ T11, 0 ] * Y
    if (*rank < n) {
	Rtzrzf(*rank, n, &A[0], lda, &work[mn + 1], &work[(mn * 2) + 1], lwork - (mn * 2), info);
    }
//workspace: 2*MN.
//Details of Householder rotations stored in WORK(MN+1:2*MN)
//B(1:M,1:NRHS) := Q' * B(1:M,1:NRHS)
    Rormqr("Left", "Transpose", m, nrhs, mn, &A[0], lda, &work[0], &B[0], ldb, &work[(mn * 2) + 1], lwork - (mn * 2), info);
    mtemp1 = wsize, mtemp2 = (mn * 2) + work[(mn * 2) + 1];
    wsize = max(mtemp1, mtemp2);
//workspace: 2*MN+NB*NRHS.
//B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS)
    Rtrsm("Left", "Upper", "No transpose", "Non-unit", *rank, nrhs, One, &A[0], lda, &B[0], ldb);
    for (j = 0; j < nrhs; j++) {
	for (i = *rank + 1; i <= n; i++) {
	    B[i + j * ldb] = Zero;
	}
    }
//B(1:N,1:NRHS) := Y' * B(1:N,1:NRHS)
    if (*rank < n) {
	Rormrz("Left", "Transpose", n, nrhs, *rank, n - *rank, &A[0], lda, &work[mn + 1], &B[0], ldb, &work[(mn << 1) + 1], lwork - (mn * 2), info);
    }
//workspace: 2*MN+NRHS.
//B(1:N,1:NRHS) := P * B(1:N,1:NRHS)
    for (j = 0; j < nrhs; j++) {
	for (i = 0; i < n; i++) {
	    work[jpvt[i]] = B[i + j * ldb];
	}
	Rcopy(n, &work[0], 1, &B[j * ldb + 1], 1);
    }
//workspace: N.
//Undo scaling
    if (iascl == 1) {
	Rlascl("G", 0, 0, anrm, smlnum, n, nrhs, &B[0], ldb, info);
	Rlascl("U", 0, 0, smlnum, anrm, *rank, *rank, &A[0], lda, info);
    } else if (iascl == 2) {
	Rlascl("G", 0, 0, anrm, bignum, n, nrhs, &B[0], ldb, info);
	Rlascl("U", 0, 0, bignum, anrm, *rank, *rank, &A[0], lda, info);
    }
    if (ibscl == 1) {
	Rlascl("G", 0, 0, smlnum, bnrm, n, nrhs, &B[0], ldb, info);
    } else if (ibscl == 2) {
	Rlascl("G", 0, 0, bignum, bnrm, n, nrhs, &B[0], ldb, info);
    }
  L70:
    work[1] = lwkopt;
    return;
}
