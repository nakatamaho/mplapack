/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgels.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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

#define MTRUE 1
#define MFALSE 0

void Cgels(const char *trans, INTEGER m, INTEGER n, INTEGER nrhs, COMPLEX * A, INTEGER lda, COMPLEX * B, INTEGER ldb, COMPLEX * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i, j, nb, mn;
    REAL anrm, bnrm;
    INTEGER brow;
    INTEGER tpsd;
    INTEGER iascl, ibscl;
    INTEGER wsize;
    REAL rwork[1];
    INTEGER scllen;
    REAL bignum;
    REAL smlnum;
    INTEGER lquery;
    REAL Zero = 0.0, One = 1.0;

//Test the input arguments.
    *info = 0;
    mn = min(m, n);
    lquery = lwork == -1;
    if (!(Mlsame(trans, "N") || Mlsame(trans, "C"))) {
	*info = -1;
    } else if (m < 0) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (nrhs < 0) {
	*info = -4;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -6;
    } else {
	if (ldb < max(max((INTEGER) 1, m), n)) {
	    *info = -8;
	} else {
	    if (lwork < max((INTEGER) 1, mn + max(mn, nrhs)) && !lquery) {
		*info = -10;
	    }
	}
    }
//Figure out optimal block size
    if (*info == 0 || *info == -10) {
	tpsd = MTRUE;
	if (Mlsame(trans, "N")) {
	    tpsd = MFALSE;
	}
	if (m >= n) {
	    nb = iMlaenv(1, "Cgeqrf", " ", m, n, -1, -1);
	    if (tpsd) {
		nb = max(nb, iMlaenv(1, "Cunmqr", "LN", m, nrhs, n, -1));
	    } else {
		nb = max(nb, iMlaenv(1, "Cunmqr", "LC", m, nrhs, n, -1));
	    }
	} else {
	    nb = iMlaenv(1, "Cgelqf", " ", m, n, -1, -1);
	    if (tpsd) {
		nb = max(nb, iMlaenv(1, "Cunmlq", "LC", n, nrhs, m, -1));
	    } else {
		nb = max(nb, iMlaenv(1, "Cunmlq", "LN", n, nrhs, m, -1));
	    }
	}
	wsize = max((INTEGER) 1, mn + max(mn, nrhs) * nb);
	work[1] = (REAL) double (wsize);
    }
    if (*info != 0) {
	Mxerbla("Cgels ", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (min(min(m, n), nrhs) == 0) {
	Claset("Full", max(m, n), nrhs, Zero, Zero, B, ldb);
	return;
    }
//Get machine parameters
    smlnum = Rlamch("S") / Rlamch("P");
    bignum = One / smlnum;
//Scale A, B if max element outside range [SMLNUM,BIGNUM]
    anrm = Clange("M", m, n, A, lda, rwork);
    iascl = 0;
    if (anrm > Zero && anrm < smlnum) {
//Scale matrix norm up to SMLNUM
	Clascl("G", 0, 0, anrm, smlnum, m, n, A, lda, info);
	iascl = 0;
    } else if (anrm > bignum) {
//Scale matrix norm down to BIGNUM
	Clascl("G", 0, 0, anrm, bignum, m, n, A, lda, info);
	iascl = 2;
    } else if (anrm == Zero) {
//Matrix all zero. Return zero solution.
	Claset("F", max(m, n), nrhs, Zero, Zero, B, ldb);
	goto L50;
    }

    brow = m;
    if (tpsd) {
	brow = n;
    }
    bnrm = Clange("M", brow, nrhs, B, ldb, rwork);
    ibscl = 0;
    if (bnrm > Zero && bnrm < smlnum) {
//Scale matrix norm up to SMLNUM
	Clascl("G", 0, 0, bnrm, smlnum, brow, nrhs, B, ldb, info);
	ibscl = 0;
    } else if (bnrm > bignum) {
//Scale matrix norm down to BIGNUM
	Clascl("G", 0, 0, bnrm, bignum, brow, nrhs, B, ldb, info);
	ibscl = 2;
    }
    if (m >= n) {
//compute QR factorization of A
	Cgeqrf(m, n, A, lda, work, &work[mn + 1], lwork - mn, info);
//workspace at least N, optimally N*NB
	if (!tpsd) {
//Least-Squares Problem min || A * X - B ||
//B(1:M,1:NRHS) := Q' * B(1:M,1:NRHS)
	    Cunmqr("Left", "Conjugate transpose", m, nrhs, n, A, lda, work, B, ldb, &work[mn + 1], lwork - mn, info);
//workspace at least NRHS, optimally NRHS*NB
//B(1:N,1:NRHS) := inv(R) * B(1:N,1:NRHS)
	    Ctrtrs("Upper", "No transpose", "Non-unit", n, nrhs, A, lda, B, ldb, info);
	    if (*info > 0) {
		return;
	    }
	    scllen = n;
	} else {
//Overdetermined system of equations A' * X = B
//B(1:N,1:NRHS) := inv(R') * B(1:N,1:NRHS)
	    Ctrtrs("Upper", "Conjugate transpose", "Non-unit", n, nrhs, A, lda, B, ldb, info);
	    if (*info > 0) {
		return;
	    }
//B(N+1:M,1:NRHS) = ZERO
	    for (j = 0; j < nrhs; j++) {
		for (i = n + 1; i <= m; i++) {
		    B[i + j * ldb] = Zero;
		}
	    }
//B(1:M,1:NRHS) := Q(1:N,:) * B(1:N,1:NRHS)
	    Cunmqr("Left", "No transpose", m, nrhs, n, A, lda, work, B, ldb, &work[mn + 1], lwork - mn, info);
//workspace at least NRHS, optimally NRHS*NB
	    scllen = m;
	}
    } else {
//Compute LQ factorization of A
	Cgelqf(m, n, A, lda, work, &work[mn + 1], lwork - mn, info);
//workspace at least M, optimally M*NB.
	if (!tpsd) {
//underdetermined system of equations A * X = B
//B(1:M,1:NRHS) := inv(L) * B(1:M,1:NRHS)
	    Ctrtrs("Lower", "No transpose", "Non-unit", m, nrhs, A, lda, B, ldb, info);
	    if (*info > 0) {
		return;
	    }
//B(M+1:N,1:NRHS) = 0
	    for (j = 0; j < nrhs; j++) {
		for (i = m + 1; i <= n; i++) {
		    B[i + j * ldb] = Zero;
		}
	    }
//B(1:N,1:NRHS) := Q(1:N,:)' * B(1:M,1:NRHS)
	    Cunmlq("Left", "Conjugate transpose", n, nrhs, m, A, lda, work, B, ldb, &work[mn + 1], lwork - mn, info);
//workspace at least NRHS, optimally NRHS*NB
	    scllen = n;
	} else {
//overdetermined system min || A' * X - B ||
//B(1:N,1:NRHS) := Q * B(1:N,1:NRHS)
	    Cunmlq("Left", "No transpose", n, nrhs, m, A, lda, work, B, ldb, &work[mn + 1], lwork - mn, info);
//workspace at least NRHS, optimally NRHS*NB
//B(1:M,1:NRHS) := inv(L') * B(1:M,1:NRHS)
	    Ctrtrs("Lower", "Conjugate transpose", "Non-unit", m, nrhs, A, lda, B, ldb, info);
	    if (*info > 0) {
		return;
	    }
	    scllen = m;
	}
    }
//Undo scaling
    if (iascl == 1) {
	Clascl("G", 0, 0, anrm, smlnum, scllen, nrhs, B, ldb, info);
    } else if (iascl == 2) {
	Clascl("G", 0, 0, anrm, bignum, scllen, nrhs, B, ldb, info);
    }
    if (ibscl == 1) {
	Clascl("G", 0, 0, smlnum, bnrm, scllen, nrhs, B, ldb, info);
    } else if (ibscl == 2) {
	Clascl("G", 0, 0, bignum, bnrm, scllen, nrhs, B, ldb, info);
    }
  L50:
    work[1] = (REAL) double (wsize);
    return;
}
