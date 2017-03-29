/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgelsd.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cgelsd(INTEGER m, INTEGER n, INTEGER nrhs,
	    COMPLEX * A, INTEGER lda, COMPLEX * B, INTEGER ldb, REAL * s, REAL rcond, INTEGER * rank, COMPLEX * work, INTEGER lwork, REAL * rwork, INTEGER * iwork, INTEGER * info)
{
    INTEGER ie, il, mm;
    REAL eps, anrm, bnrm;
    INTEGER itau, nlvl, iascl, ibscl;
    REAL sfmin;
    INTEGER minmn, maxmn, itaup, itauq, mnthr = 0, nwork;
    REAL bignum;
    INTEGER ldwork;
    INTEGER liwork, minwrk, maxwrk;
    REAL smlnum;
    INTEGER lrwork;
    LOGICAL lquery;
    INTEGER nrwork, smlsiz = 0;
    REAL Zero = 0.0, One = 1.0, Two = 2.0;

//Test the input arguments.
    *info = 0;
    minmn = min(m, n);
    maxmn = max(m, n);
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
//Compute workspace.
//(Note: Comments in the code beginning "Workspace:" describe the
//minimal amount of workspace needed at that point in the code,
//as well as the preferred amount for good performance.
//NB refers to the optimal block size for the immediately
//following subroutine, as returned by ILAENV.)
    if (*info == 0) {
	minwrk = 0;
	maxwrk = 0;
	liwork = 0;
	lrwork = 0;
	if (minmn > 0) {
	    smlsiz = iMlaenv(9, "Cgelsd", " ", 0, 0, 0, 0);
	    mnthr = iMlaenv(6, "Cgelsd", " ", m, n, nrhs, -1);
	    nlvl = max((INTEGER) cast2double((log((double) minmn / (double) (smlsiz + 1)) / log(Two)) + 1), (INTEGER) 0);
	    liwork = minmn * 3 * nlvl + minmn * 11;
	    mm = m;
	    if (m >= n && m >= mnthr) {
//Path 1a - overdetermined, with many more rows than
//          columns.
		mm = n;
		maxwrk = max(maxwrk, n * iMlaenv(1, "Cgeqrf", " ", m, n, -1, -1));
		maxwrk = max(maxwrk, nrhs * iMlaenv(1, "Cunmqr", "LC", m, nrhs, n, -1));
	    }
	    if (m >= n) {
//Path 1 - overdetermined or exactly determined.
		lrwork = n * 10 + (n << 1) * smlsiz + (n << 3) * nlvl + smlsiz * 3 * nrhs + (smlsiz + 1) * (smlsiz + 1);
		maxwrk = max(maxwrk, (n << 1) + (mm + n) * iMlaenv(1, "Cgebrd", " ", mm, n, -1, -1));
		maxwrk = max(maxwrk, (n << 1) + nrhs * iMlaenv(1, "Cunmbr", "qlc", mm, nrhs, n, -1));
		maxwrk = max(maxwrk, (n << 1) + (n - 1) * iMlaenv(1, "Cunmbr", "pln", n, nrhs, n, -1));
		maxwrk = max(maxwrk, (n << 1) + n * nrhs);
		minwrk = max((n << 1) + mm, (n << 1) + n * nrhs);
	    }
	    if (n > m) {
		lrwork = m * 10 + (m << 1) * smlsiz + (m << 3) * nlvl + smlsiz * 3 * nrhs + (smlsiz + 1) * (smlsiz + 1);
		if (n >= mnthr) {
//Path 2a - underdetermined, with many more columns
//          than rows.
		    maxwrk = m + m * iMlaenv(1, "Cgelqf", " ", m, n, -1, -1);
		    maxwrk = max(maxwrk, m * m + (m << 2) + (m << 1) * iMlaenv(1, "Cgebrd", " ", m, m, -1, -1));
		    maxwrk = max(maxwrk, m * m + (m << 2) + nrhs * iMlaenv(1, "Cunmbr", "qlc", m, nrhs, m, -1));
		    maxwrk = max(maxwrk, m * m + (m << 2) + (m - 1) * iMlaenv(1, "Cunmlq", "lc", n, nrhs, m, -1));
		    if (nrhs > 1) {
			maxwrk = max(maxwrk, m * m + m + m * nrhs);
		    } else {
			maxwrk = max(maxwrk, m * m + (m << 1));
		    }
		    maxwrk = max(maxwrk, m * m + (m << 2) + m * nrhs);
		} else {
//Path 2 - underdetermined.
		    maxwrk = (m << 1) + (n + m) * iMlaenv(1, "Cgebrd", " ", m, n, -1, -1);
		    maxwrk = max(maxwrk, (m << 1) + nrhs * iMlaenv(1, "Cunmbr", "qlc", m, nrhs, m, -1));
		    maxwrk = max(maxwrk, (m << 1) + m * iMlaenv(1, "Cunmbr", "pln", n, nrhs, m, -1));
		    maxwrk = max(maxwrk, (m << 1) + m * nrhs);
		}
		minwrk = max((m << 1) + n, (m << 1) + m * nrhs);
	    }
	}
	minwrk = min(minwrk, maxwrk);
	work[1] = maxwrk;
	iwork[1] = liwork;
	rwork[1] = lrwork;
	if (lwork < minwrk && !lquery) {
	    *info = -12;
	}
    }
    if (*info != 0) {
	Mxerbla("Cgelsd", -(*info));
	return;
    } else if (lquery) {
	return;
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
    //Rlabad(&smlnum, &bignum);
//Scale A if max entry outside range [SMLNUM,BIGNUM].
    anrm = Clange("M", m, n, &A[0], lda, &rwork[1]);
    iascl = 0;
    if (anrm > Zero && anrm < smlnum) {
//Scale matrix norm up to SMLNUM
	Clascl("G", 0, 0, anrm, smlnum, m, n, &A[0], lda, info);
	iascl = 0;
    } else if (anrm > bignum) {
//Scale matrix norm down to BIGNUM.
	Clascl("G", 0, 0, anrm, bignum, m, n, &A[0], lda, info);
	iascl = 2;
    } else if (anrm == Zero) {
//Matrix all zero. Return zero solution.
	Claset("F", max(m, n), nrhs, Zero, Zero, &B[0], ldb);
	Rlaset("F", minmn, 1, Zero, Zero, &s[1], 1);
	*rank = 0;
	goto L10;
    }
//Scale B if max entry outside range [SMLNUM,BIGNUM].
    bnrm = Clange("M", m, nrhs, &B[0], ldb, &rwork[1]);
    ibscl = 0;
    if (bnrm > Zero && bnrm < smlnum) {
//Scale matrix norm up to SMLNUM.
	Clascl("G", 0, 0, bnrm, smlnum, m, nrhs, &B[0], ldb, info);
	ibscl = 0;
    } else if (bnrm > bignum) {
//Scale matrix norm down to BIGNUM.
	Clascl("G", 0, 0, bnrm, bignum, m, nrhs, &B[0], ldb, info);
	ibscl = 2;
    }
//If M < N make sure B(M+1:N,:) = 0
    if (m < n) {
	Claset("F", n - m, nrhs, Zero, Zero, &B[m + 1 + ldb], ldb);
    }
//Overdetermined case.
    if (m >= n) {
//Path 1 - overdetermined or exactly determined.
	mm = m;
	if (m >= mnthr) {
//Path 1a - overdetermined, with many more rows than columns
	    mm = n;
	    itau = 1;
	    nwork = itau + n;
//Compute A=Q*R.
//(RWorkspace: need N)
//(CWorkspace: need N, prefer N*NB)
	    Cgeqrf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, info);
//Multiply B by transpose(Q).
//(RWorkspace: need N)
//(CWorkspace: need NRHS, prefer NRHS*NB)
	    Cunmqr("L", "C", m, nrhs, n, &A[0], lda, &work[itau], &B[0], ldb, &work[nwork], lwork - nwork + 1, info);
//Zero out below R.
	    if (n > 1) {
		Claset("L", n - 1, n - 1, Zero, Zero, &A[lda + 2], lda);
	    }
	}
	itauq = 1;
	itaup = itauq + n;
	nwork = itaup + n;
	ie = 1;
	nrwork = ie + n;
//Bidiagonalize R in A.
//(RWorkspace: need N)
//(CWorkspace: need 2*N+MM, prefer 2*N+(MM+N)*NB)
	Cgebrd(mm, n, &A[0], lda, &s[1], &rwork[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, info);
//Multiply B by transpose of left bidiagonalizing vectors of R.
//(CWorkspace: need 2*N+NRHS, prefer 2*N+NRHS*NB)
	Cunmbr("Q", "L", "C", mm, nrhs, n, &A[0], lda, &work[itauq], &B[0], ldb, &work[nwork], lwork - nwork + 1, info);
//Solve the bidiagonal least squares problem.
	Clalsd("U", smlsiz, n, nrhs, &s[1], &rwork[ie], &B[0], ldb, rcond, rank, &work[nwork], &rwork[nrwork], &iwork[1], info);
	if (*info != 0) {
	    goto L10;
	}
//Multiply B by right bidiagonalizing vectors of R.
	Cunmbr("P", "L", "N", n, nrhs, n, &A[0], lda, &work[itaup], &B[0], ldb, &work[nwork], lwork - nwork + 1, info);
    } else {
	if (n >= mnthr && lwork >= (m << 2) + m * m + max(max(max(m, m * 2 - 4), nrhs), n - m * 3)) {
//Path 2a - underdetermined, with many more columns than rows
//and sufficient workspace for an efficient algorithm.
	    ldwork = m;
	    if (lwork >= max((m << 2) + m * lda + max(max(max(m, (m << 1) - 4), nrhs), n - m * 3), m * lda + m + m * nrhs)) {
		ldwork = lda;
	    }
	    itau = 1;
	    nwork = m + 1;
//Compute A=L*Q.
//(CWorkspace: need 2*M, prefer M+M*NB)
	    Cgelqf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, info);
	    il = nwork;
//Copy L to WORK(IL), zeroing out above its diagonal.
	    Clacpy("L", m, m, &A[0], lda, &work[il], ldwork);
	    Claset("U", m - 1, m - 1, Zero, Zero, &work[il + ldwork], ldwork);
	    itauq = il + ldwork * m;
	    itaup = itauq + m;
	    nwork = itaup + m;
	    ie = 1;
	    nrwork = ie + m;
//Bidiagonalize L in WORK(IL).
//(RWorkspace: need M)
//(CWorkspace: need M*M+4*M, prefer M*M+4*M+2*M*NB)
	    Cgebrd(m, m, &work[il], ldwork, &s[1], &rwork[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, info);
//Multiply B by transpose of left bidiagonalizing vectors of L.
//(CWorkspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)
	    Cunmbr("Q", "L", "C", m, nrhs, m, &work[il], ldwork, &work[itauq], &B[0], ldb, &work[nwork], lwork - nwork + 1, info);
//Solve the bidiagonal least squares problem.
	    Clalsd("U", smlsiz, m, nrhs, &s[1], &rwork[ie], &B[0], ldb, rcond, rank, &work[nwork], &rwork[nrwork], &iwork[1], info);
	    if (*info != 0) {
		goto L10;
	    }
//Multiply B by right bidiagonalizing vectors of L.
	    Cunmbr("P", "L", "N", m, nrhs, m, &work[il], ldwork, &work[itaup], &B[0], ldb, &work[nwork], lwork - nwork + 1, info);
//Zero out below first M rows of B.
	    Claset("F", n - m, nrhs, Zero, Zero, &B[m + 1 + ldb], ldb);
	    nwork = itau + m;
//Multiply transpose(Q) by B.
//(CWorkspace: need NRHS, prefer NRHS*NB)
	    Cunmlq("L", "C", n, nrhs, m, &A[0], lda, &work[itau], &B[0], ldb, &work[nwork], lwork - nwork + 1, info);
	} else {
//Path 2 - remaining underdetermined cases.
	    itauq = 1;
	    itaup = itauq + m;
	    nwork = itaup + m;
	    ie = 1;
	    nrwork = ie + m;
//Bidiagonalize A.
//(RWorkspace: need M)
//(CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)
	    Cgebrd(m, n, &A[0], lda, &s[1], &rwork[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, info);
//Multiply B by transpose of left bidiagonalizing vectors.
//(CWorkspace: need 2*M+NRHS, prefer 2*M+NRHS*NB)
	    Cunmbr("Q", "L", "C", m, nrhs, n, &A[0], lda, &work[itauq], &B[0], ldb, &work[nwork], lwork - nwork + 1, info);
//Solve the bidiagonal least squares problem.
	    Clalsd("L", smlsiz, m, nrhs, &s[1], &rwork[ie], &B[0], ldb, rcond, rank, &work[nwork], &rwork[nrwork], &iwork[1], info);
	    if (*info != 0) {
		goto L10;
	    }
//Multiply B by right bidiagonalizing vectors of A.
	    Cunmbr("P", "L", "N", n, nrhs, m, &A[0], lda, &work[itaup]
		   , &B[0], ldb, &work[nwork], lwork - nwork + 1, info);
	}
    }
//Undo scaling.
    if (iascl == 1) {
	Clascl("G", 0, 0, anrm, smlnum, n, nrhs, &B[0], ldb, info);
	Rlascl("G", 0, 0, smlnum, anrm, minmn, 1, &s[1], minmn, info);
    } else if (iascl == 2) {
	Clascl("G", 0, 0, anrm, bignum, n, nrhs, &B[0], ldb, info);
	Rlascl("G", 0, 0, bignum, anrm, minmn, 1, &s[1], minmn, info);
    }
    if (ibscl == 1) {
	Clascl("G", 0, 0, smlnum, bnrm, n, nrhs, &B[0], ldb, info);
    } else if (ibscl == 2) {
	Clascl("G", 0, 0, bignum, bnrm, n, nrhs, &B[0], ldb, info);
    }

  L10:
    work[1] = maxwrk;
    iwork[1] = liwork;
    rwork[1] = lrwork;
    return;
}
