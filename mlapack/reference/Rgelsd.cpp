/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgelsd.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rgelsd(INTEGER m, INTEGER n, INTEGER nrhs,
	    REAL * A, INTEGER lda, REAL * B, INTEGER ldb, REAL * s, REAL * rcond, INTEGER * rank, REAL * work, INTEGER lwork, INTEGER * iwork, INTEGER * info)
{
    INTEGER ie, il, mm;
    REAL eps, anrm, bnrm;
    INTEGER itau, nlvl, iascl, ibscl;
    REAL sfmin;
    INTEGER minmn, maxmn, itaup, itauq, mnthr, nwork;
    REAL bignum;
    INTEGER wlalsd;
    INTEGER ldwork;
    INTEGER minwrk, maxwrk;
    REAL smlnum;
    INTEGER lquery;
    INTEGER smlsiz;
    REAL Zero = 0.0, One = 1.0;
//Test the input arguments.
    *info = 0;
    minmn = min(m, n);
    maxmn = max(m, n);
    mnthr = iMlaenv(6, "Rgelsd", " ", m, n, nrhs, -1);
    lquery = lwork == -1;
    if (m < 0) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (nrhs < 0) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -5;
    } else if (ldb < max((INTEGER) 1, maxmn)) {
	*info = -7;
    }
    smlsiz = iMlaenv(9, "Rgelsd", " ", 0, 0, 0, 0);
//Compute workspace.
//(Note: Comments in the code beginning "Workspace:" describe the
//minimal amount of workspace needed at that point in the code,
//as well as the preferred amount for good performance.
//NB refers to the optimal block size for the immediately
//following subroutine, as returned by ILAENV.)
    minwrk = 0;
    minmn = max((INTEGER) 1, minmn);
    nlvl = max((INTEGER) cast2double((log(double (minmn) / double (smlsiz + 1)) / log(2.)) + 1.), (INTEGER) 0);
    if (*info == 0) {
	maxwrk = 0;
	mm = m;
	if (m >= n && m >= mnthr) {
//Path 1a - overdetermined, with many more rows than columns.
	    mm = n;
	    maxwrk = max(maxwrk, n + n * iMlaenv(1, "DGEQRF", " ", m, n, -1, -1));
	    maxwrk = max(maxwrk, n + nrhs * iMlaenv(1, "DORMQR", "LT", m, nrhs, n, -1));
	}
	if (m >= n) {
//Path 1 - overdetermined or exactly determined.
	    maxwrk = max(maxwrk, n * 3 + (mm + n) * iMlaenv(1, "DGEBRD", " ", mm, n, -1, -1));
	    maxwrk = max(maxwrk, n * 3 + nrhs * iMlaenv(1, "DORMBR", "QLT", mm, nrhs, n, -1));
	    maxwrk = max(maxwrk, n * 3 + (n - 1) * iMlaenv(1, "DORMBR", "PLN", n, nrhs, n, -1));
	    wlalsd = n * 9 + (n * 2) * smlsiz + (n * 8) * nlvl + n * nrhs + (smlsiz + 1) * (smlsiz + 1);
	    maxwrk = max(maxwrk, n * 3 + wlalsd);
	    minwrk = max(max(n * 3 + mm, n * 3 + nrhs), n * 3 + wlalsd);
	}
	if (n > m) {
	    wlalsd = m * 9 + (m * 2) * smlsiz + (m * 8) * nlvl + m * nrhs + (smlsiz + 1) * (smlsiz + 1);
	    if (n >= mnthr) {
//Path 2a - underdetermined, with many more columns
//than rows.
		maxwrk = m + m * iMlaenv(1, "DGELQF", " ", m, n, -1, -1);
		maxwrk = max(maxwrk, m * m + (m << 2) + (m * 2) * iMlaenv(1, "DGEBRD", " ", m, m, -1, -1));
		maxwrk = max(maxwrk, m * m + (m << 2) + nrhs * iMlaenv(1, "DORMBR", "QLT", m, nrhs, m, -1));
		maxwrk = max(maxwrk, m * m + (m << 2) + (m - 1) * iMlaenv(1, "DORMBR", "PLN", m, nrhs, m, -1));
		if (nrhs > 1) {
		    maxwrk = max(maxwrk, m * m + m + m * nrhs);
		} else {
		    maxwrk = max(maxwrk, m * m + (m * 2));
		}
		maxwrk = max(maxwrk, m + nrhs * iMlaenv(1, "DORMLQ", "LT", n, nrhs, m, -1));
		maxwrk = max(maxwrk, m * m + (m << 2) + wlalsd);
	    } else {
//Path 2 - remaining underdetermined cases.
		maxwrk = m * 3 + (n + m) * iMlaenv(1, "Rgebrd", " ", m, n, -1, -1);
		maxwrk = max(maxwrk, m * 3 + nrhs * iMlaenv(1, "Rormbr", "qlt", m, nrhs, n, -1));
		maxwrk = max(maxwrk, m * 3 + m * iMlaenv(1, "Rormbr", "pln", n, nrhs, m, -1));
		maxwrk = max(maxwrk, m * 3 + wlalsd);
	    }
	    minwrk = max(max(m * 3 + nrhs, m * 3 + m), m * 3 + wlalsd);
	}
	minwrk = min(minwrk, maxwrk);
	work[1] = maxwrk;
	if (lwork < minwrk && !lquery) {
	    *info = -12;
	}
    }
    if (*info != 0) {
	Mxerbla("Rgelsd", -(*info));
	return;
    } else if (lquery) {
	goto L10;
    }
//Quick return if possible.
    if (m == 0 || n == 0) {
	*rank = 0;
	return;
    }
//Get machine parameters.
    eps = Rlamch("P");
    sfmin = Rlamch("S");
    smlnum = sfmin / eps;
    bignum = One / smlnum;
    //    Rlabad(&smlnum, &bignum);
//Scale A if max entry outside range [SMLNUM,BIGNUM].
    anrm = Rlange("M", m, n, &A[0], lda, &work[0]);
    iascl = 0;
    if (anrm > Zero && anrm < smlnum) {
//Scale matrix norm up to SMLNUM.
	Rlascl("G", 0, 0, anrm, smlnum, m, n, &A[0], lda, info);
	iascl = 0;
    } else if (anrm > bignum) {
//Scale matrix norm down to BIGNUM.
	Rlascl("G", 0, 0, anrm, bignum, m, n, &A[0], lda, info);
	iascl = 2;
    } else if (anrm == Zero) {
//Matrix all zero. Return zero solution.
	Rlaset("F", max(m, n), nrhs, Zero, Zero, &B[0], ldb);
	Rlaset("F", minmn, 1, Zero, Zero, &s[1], 1);
	*rank = 0;
	goto L10;
    }
//Scale B if max entry outside range [SMLNUM,BIGNUM].
    bnrm = Rlange("M", m, nrhs, &B[0], ldb, &work[0]);
    ibscl = 0;
    if (bnrm > Zero && bnrm < smlnum) {
//Scale matrix norm up to SMLNUM.
	Rlascl("G", 0, 0, bnrm, smlnum, m, nrhs, &B[0], ldb, info);
	ibscl = 0;
    } else if (bnrm > bignum) {
//Scale matrix norm down to BIGNUM.
	Rlascl("G", 0, 0, bnrm, bignum, m, nrhs, &B[0], ldb, info);
	ibscl = 2;
    }
//If M < N make sure certain entries of B are zero.
    if (m < n) {
	Rlaset("F", n - m, nrhs, Zero, Zero, &B[m + 1 + ldb], ldb);
    }
//Overdetermined case.
    if (m >= n) {
//Path 1 - overdetermined or exactly determined.
	mm = m;
	if (m >= mnthr) {
//Path 1a - overdetermined, with many more rows than columns.
	    mm = n;
	    itau = 1;
	    nwork = itau + n;
//Compute A=Q*R.
//(Workspace: need 2*N, prefer N+N*NB)
	    Rgeqrf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, info);
//Multiply B by transpose(Q).
//(Workspace: need N+NRHS, prefer N+NRHS*NB)
	    Rormqr("L", "T", m, nrhs, n, &A[0], lda, &work[itau], &B[0], ldb, &work[nwork], lwork - nwork + 1, info);
//Zero out below R.
	    if (n > 1) {
		Rlaset("L", n - 1, n - 1, Zero, Zero, &A[lda + 2], lda);
	    }
	}
	ie = 1;
	itauq = ie + n;
	itaup = itauq + n;
	nwork = itaup + n;
//Bidiagonalize R in A.
//(Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB)
	Rgebrd(mm, n, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, info);
//Multiply B by transpose of left bidiagonalizing vectors of R.
//(Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB)
	Rormbr("Q", "L", "T", mm, nrhs, n, &A[0], lda, &work[itauq], &B[0], ldb, &work[nwork], lwork - nwork + 1, info);
//Solve the bidiagonal least squares problem.
	Rlalsd("U", smlsiz, n, nrhs, &s[1], &work[ie], &B[0], ldb, *rcond, rank, &work[nwork], &iwork[1], info);
	if (*info != 0) {
	    goto L10;
	}
//Multiply B by right bidiagonalizing vectors of R.
	Rormbr("P", "L", "N", n, nrhs, n, &A[0], lda, &work[itaup], &B[0], ldb, &work[nwork], lwork - nwork + 1, info);
    } else {
	if (n >= mnthr && lwork >= (m << 2) + m * m + max(max(max(max(m, (m * 2) - 4), nrhs), n - m * 3), wlalsd)) {
//Path 2a - underdetermined, with many more columns than rows
//and sufficient workspace for an efficient algorithm.
	    ldwork = m;
	    if (lwork >= max(max((m << 2) + m * lda + max(max(max(m, (m * 2) - 4), nrhs), n - m * 3), m * lda + m + m * nrhs), (m << 2) + m * lda + wlalsd)) {
		ldwork = lda;
	    }
	    itau = 1;
	    nwork = m + 1;
//Compute A=L*Q.
//(Workspace: need 2*M, prefer M+M*NB)
	    Rgelqf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, info);
	    il = nwork;
//Copy L to WORK(IL), zeroing out above its diagonal.
	    Rlacpy("L", m, m, &A[0], lda, &work[il], ldwork);
	    Rlaset("U", m - 1, m - 1, Zero, Zero, &work[il + ldwork], ldwork);
	    ie = il + ldwork * m;
	    itauq = ie + m;
	    itaup = itauq + m;
	    nwork = itaup + m;
//Bidiagonalize L in WORK(IL).
//(Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB)
	    Rgebrd(m, m, &work[il], ldwork, &s[1], &work[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, info);
//Multiply B by transpose of left bidiagonalizing vectors of L.
//(Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)
	    Rormbr("Q", "L", "T", m, nrhs, m, &work[il], ldwork, &work[itauq], &B[0], ldb, &work[nwork], lwork - nwork + 1, info);
//Solve the bidiagonal least squares problem.
	    Rlalsd("U", smlsiz, m, nrhs, &s[1], &work[ie], &B[0], ldb, *rcond, rank, &work[nwork], &iwork[1], info);
	    if (*info != 0) {
		goto L10;
	    }
//Multiply B by right bidiagonalizing vectors of L.
	    Rormbr("P", "L", "N", m, nrhs, m, &work[il], ldwork, &work[itaup], &B[0], ldb, &work[nwork], lwork - nwork + 1, info);
//Zero out below first M rows of B.
	    Rlaset("F", n - m, nrhs, Zero, Zero, &B[m + 1 + ldb], ldb);
	    nwork = itau + m;
//Multiply transpose(Q) by B.
//(Workspace: need M+NRHS, prefer M+NRHS*NB)
	    Rormlq("L", "T", n, nrhs, m, &A[0], lda, &work[itau], &B[0], ldb, &work[nwork], lwork - nwork + 1, info);
	} else {
//Path 2 - remaining underdetermined cases.
	    ie = 1;
	    itauq = ie + m;
	    itaup = itauq + m;
	    nwork = itaup + m;
//Bidiagonalize A.
//(Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
	    Rgebrd(m, n, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, info);
//Multiply B by transpose of left bidiagonalizing vectors.
//(Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB)
	    Rormbr("Q", "L", "T", m, nrhs, n, &A[0], lda, &work[itauq]
		   , &B[0], ldb, &work[nwork], lwork - nwork + 1, info);
//Solve the bidiagonal least squares problem.
	    Rlalsd("L", smlsiz, m, nrhs, &s[1], &work[ie], &B[0], ldb, *rcond, rank, &work[nwork], &iwork[1], info);
	    if (*info != 0) {
		goto L10;
	    }
//Multiply B by right bidiagonalizing vectors of A.
	    Rormbr("P", "L", "N", n, nrhs, m, &A[0], lda, &work[itaup]
		   , &B[0], ldb, &work[nwork], lwork - nwork + 1, info);
	}
    }
//Undo scaling.
    if (iascl == 1) {
	Rlascl("G", 0, 0, anrm, smlnum, n, nrhs, &B[0], ldb, info);
	Rlascl("G", 0, 0, smlnum, anrm, minmn, 1, &s[1], minmn, info);
    } else if (iascl == 2) {
	Rlascl("G", 0, 0, anrm, bignum, n, nrhs, &B[0], ldb, info);
	Rlascl("G", 0, 0, bignum, anrm, minmn, 1, &s[1], minmn, info);
    }
    if (ibscl == 1) {
	Rlascl("G", 0, 0, smlnum, bnrm, n, nrhs, &B[0], ldb, info);
    } else if (ibscl == 2) {
	Rlascl("G", 0, 0, bignum, bnrm, n, nrhs, &B[0], ldb, info);
    }
  L10:
    work[1] = maxwrk;
    return;
}
