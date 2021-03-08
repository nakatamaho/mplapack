/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgelss.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cgelss(INTEGER m, INTEGER n, INTEGER nrhs,
	    COMPLEX * A, INTEGER lda, COMPLEX * B, INTEGER ldb, REAL * s, REAL rcond, INTEGER * rank, COMPLEX * work, INTEGER lwork, REAL * rwork, INTEGER * info)
{
    INTEGER i, bl, ie, il, mm;
    REAL eps, thr, anrm, bnrm;
    INTEGER itau;
    COMPLEX vdum;
    INTEGER iascl, ibscl, chunk;
    REAL sfmin;
    INTEGER minmn;
    INTEGER maxmn, itaup, itauq, mnthr = 0;
    INTEGER iwork;
    REAL bignum;
    INTEGER ldwork;
    INTEGER minwrk, maxwrk;
    REAL smlnum;
    INTEGER irwork;
    LOGICAL lquery;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1;

//Test the input arguments
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
//Compute workspace
// (Note: Comments in the code beginning "Workspace:" describe the
//  minimal amount of workspace needed at that point in the code,
//  as well as the preferred amount for good performance.
//  CWorkspace refers to complex workspace, and RWorkspace refers
//  to real workspace. NB refers to the optimal block size for the
//  immediately following subroutine, as returned by ILAENV.)
    if (*info == 0) {
	minwrk = 0;
	maxwrk = 0;
	if (minmn > 0) {
	    mm = m;
	    mnthr = iMlaenv(6, "Cgelss", " ", m, n, nrhs, -1);
	    if (m >= n && m >= mnthr) {
//Path 1a - overdetermined, with many more rows than
//          columns
		mm = n;
		maxwrk = max(maxwrk, n + n * iMlaenv(1, "Cgeqrf", " ", m, n, -1, -1));
		maxwrk = max(maxwrk, n + nrhs * iMlaenv(1, "Cunmqr", "lc", m, nrhs, n, -1));
	    }
	    if (m >= n) {
//Path 1 - overdetermined or exactly determined
		maxwrk = max(maxwrk, (n << 1) + (mm + n) * iMlaenv(1, "Cgebrd", " ", mm, n, -1, -1));
		maxwrk = max(maxwrk, (n << 1) + nrhs * iMlaenv(1, "Cunmbr", "qlc", mm, nrhs, n, -1));
		maxwrk = max(maxwrk, (n << 1) + (n - 1) * iMlaenv(1, "Cungbr", "p", n, n, n, -1));
		maxwrk = max(maxwrk, n * nrhs);
		minwrk = (n << 1) + max(nrhs, m);
	    }
	    if (n > m) {
		minwrk = (m << 1) + max(nrhs, n);
		if (n >= mnthr) {
//Path 2a - underdetermined, with many more columns
//than rows
		    maxwrk = m + m * iMlaenv(1, "Cgelqf", " ", m, n, -1, -1);
		    maxwrk = max(maxwrk, m * 3 + m * m + (m << 1) * iMlaenv(1, "Cgebrd", " ", m, m, -1, -1));
		    maxwrk = max(maxwrk, m * 3 + m * m + nrhs * iMlaenv(1, "Cunmbr", "qlc", m, nrhs, m, -1));
		    maxwrk = max(maxwrk, m * 3 + m * m + (m - 1) * iMlaenv(1, "Cungbr", "p", m, m, m, -1));
		    if (nrhs > 1) {
			maxwrk = max(maxwrk, m * m + m + m * nrhs);
		    } else {
			maxwrk = max(maxwrk, m * m + (m << 1));
		    }
		    maxwrk = max(maxwrk, m + nrhs * iMlaenv(1, "Cunmlq", "lc", n, nrhs, m, -1));
		} else {
//Path 2 - underdetermined
		    maxwrk = (m << 1) + (n + m) * iMlaenv(1, "ZGEBRD", " ", m, n, -1, -1);
		    maxwrk = max(maxwrk, (m << 1) + nrhs * iMlaenv(1, "ZUNMBR", "QLC", m, nrhs, m, -1));
		    maxwrk = max(maxwrk, (m << 1) + m * iMlaenv(1, "ZUNGBR", "P", m, n, m, -1));
		    maxwrk = max(maxwrk, n * nrhs);
		}
	    }
	    maxwrk = max(minwrk, maxwrk);
	}
	work[1] = maxwrk;
	if (lwork < minwrk && !lquery) {
	    *info = -12;
	}
    }
    if (*info != 0) {
	Mxerbla("Cgelss", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (m == 0 || n == 0) {
	*rank = 0;
	return;
    }
//Get machine parameters
    eps = Rlamch("P");
    sfmin = Rlamch("S");
    smlnum = sfmin / eps;
    bignum = One / smlnum;
    //Rlabad(&smlnum, &bignum);
//Scale A if max element outside range [SMLNUM,BIGNUM]
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
	Rlaset("F", minmn, 1, Zero, Zero, &s[1], minmn);
	*rank = 0;
	goto L70;
    }
//Scale B if max element outside range [SMLNUM,BIGNUM]
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
//Overdetermined case
    if (m >= n) {
//Path 1 - overdetermined or exactly determined
	mm = m;
	if (m >= mnthr) {
//Path 1a - overdetermined, with many more rows than columns
	    mm = n;
	    itau = 1;
	    iwork = itau + n;
//Compute A=Q*R
//(CWorkspace: need 2*N, prefer N+N*NB)
//(RWorkspace: none)
	    Cgeqrf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, info);
//Multiply B by transpose(Q)
//(CWorkspace: need N+NRHS, prefer N+NRHS*NB)
//(RWorkspace: none)
	    Cunmqr("L", "C", m, nrhs, n, &A[0], lda, &work[itau], &B[0], ldb, &work[iwork], lwork - iwork + 1, info);
//Zero out below R
	    if (n > 1) {
		Claset("L", n - 1, n - 1, Zero, Zero, &A[lda + 2], lda);
	    }
	}
	ie = 1;
	itauq = 1;
	itaup = itauq + n;
	iwork = itaup + n;
//Bidiagonalize R in A
//(CWorkspace: need 2*N+MM, prefer 2*N+(MM+N)*NB)
//(RWorkspace: need N)
	Cgebrd(mm, n, &A[0], lda, &s[1], &rwork[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, info);
//Multiply B by transpose of left bidiagonalizing vectors of R
//(CWorkspace: need 2*N+NRHS, prefer 2*N+NRHS*NB)
//(RWorkspace: none)
	Cunmbr("Q", "L", "C", mm, nrhs, n, &A[0], lda, &work[itauq], &B[0], ldb, &work[iwork], lwork - iwork + 1, info);
//Generate right bidiagonalizing vectors of R in A
//(CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
//(RWorkspace: none)
	Cungbr("P", n, n, n, &A[0], lda, &work[itaup], &work[iwork], lwork - iwork + 1, info);
	irwork = ie + n;
//Perform bidiagonal QR iteration
//  multiply B by transpose of left singular vectors
//  compute right singular vectors in A
//(CWorkspace: none)
//(RWorkspace: need BDSPAC)
	Cbdsqr("U", n, n, 0, nrhs, &s[1], &rwork[ie], &A[0], lda, &vdum, 1, &B[0], ldb, &rwork[irwork], info);
	if (*info != 0) {
	    goto L70;
	}
//Multiply B by reciprocals of singular values
	mtemp1 = rcond * s[1];
	thr = max(mtemp1, sfmin);
	if (rcond < Zero) {
	    mtemp1 = rcond * s[1];
	    thr = max(mtemp1, sfmin);
	}
	*rank = 0;
	for (i = 0; i < n; i++) {
	    if (s[i] > thr) {
		CRrscl(nrhs, s[i], &B[i + ldb], ldb);
		++(*rank);
	    } else {
		Claset("F", 1, nrhs, Zero, Zero, &B[i + ldb], ldb);
	    }
	}
//Multiply B by right singular vectors
//(CWorkspace: need N, prefer N*NRHS)
//(RWorkspace: none)
	if (lwork >= ldb * nrhs && nrhs > 1) {
	    Cgemm("C", "N", n, nrhs, n, One, &A[0], lda, &B[0], ldb, Zero, &work[0], ldb);
	    Clacpy("G", n, nrhs, &work[0], ldb, &B[0], ldb);
	} else if (nrhs > 1) {
	    chunk = lwork / n;
	    for (i = 1; i <= nrhs; i = i + chunk) {
		bl = min(nrhs - i + 1, chunk);
		Cgemm("C", "N", n, bl, n, One, &A[0], lda, &B[i * ldb + 1], ldb, Zero, &work[0], n);
		Clacpy("G", n, bl, &work[0], n, &B[i * ldb + 1], ldb);
	    }
	} else {
	    Cgemv("C", n, n, One, &A[0], lda, &B[0], 1, Zero, &work[0], 1);
	    Ccopy(n, &work[0], 1, &B[0], 1);
	}
    } else {
	if (n >= mnthr && lwork >= m * 3 + m * m + max(max(m, nrhs), n - (m << 1))) {
//Underdetermined case, M much less than N
//Path 2a - underdetermined, with many more columns than rows
//and sufficient workspace for an efficient algorithm
	    ldwork = m;
	    if (lwork >= m * 3 + m * lda + max(max(m, nrhs), n - (m << 1))) {
		ldwork = lda;
	    }
	    itau = 1;
	    iwork = m + 1;
//Compute A=L*Q
//(CWorkspace: need 2*M, prefer M+M*NB)
//(RWorkspace: none)
	    Cgelqf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, info);
	    il = iwork;
//Copy L to WORK(IL), zeroing out above it
	    Clacpy("L", m, m, &A[0], lda, &work[il], ldwork);
	    Claset("U", m - 1, m - 1, Zero, Zero, &work[il + ldwork], ldwork);
	    ie = 1;
	    itauq = il + ldwork * m;
	    itaup = itauq + m;
	    iwork = itaup + m;
//Bidiagonalize L in WORK(IL)
//(CWorkspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
//(RWorkspace: need M)
	    Cgebrd(m, m, &work[il], ldwork, &s[1], &rwork[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, info);
//Multiply B by transpose of left bidiagonalizing vectors of L
//(CWorkspace: need M*M+3*M+NRHS, prefer M*M+3*M+NRHS*NB)
//(RWorkspace: none)
	    Cunmbr("Q", "L", "C", m, nrhs, m, &work[il], ldwork, &work[itauq], &B[0], ldb, &work[iwork], lwork - iwork + 1, info);
//Generate right bidiagonalizing vectors of R in WORK(IL)
//(CWorkspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB)
//(RWorkspace: none)
	    Cungbr("P", m, m, m, &work[il], ldwork, &work[itaup], &work[iwork], lwork - iwork + 1, info);
	    irwork = ie + m;
//Perform bidiagonal QR iteration, computing right singular
//vectors of L in WORK(IL) and multiplying B by transpose of
//left singular vectors
//(CWorkspace: need M*M)
//(RWorkspace: need BDSPAC)
	    Cbdsqr("U", m, m, 0, nrhs, &s[1], &rwork[ie], &work[il], ldwork, &A[0], lda, &B[0], ldb, &rwork[irwork], info);
	    if (*info != 0) {
		goto L70;
	    }
//Multiply B by reciprocals of singular values
	    mtemp1 = rcond * s[1];
	    thr = max(mtemp1, sfmin);
	    if (rcond < Zero) {
		mtemp1 = rcond * s[1];
		thr = max(mtemp1, sfmin);
	    }
	    *rank = 0;
	    for (i = 0; i < m; i++) {
		if (s[i] > thr) {
		    CRrscl(nrhs, s[i], &B[i + ldb], ldb);
		    ++(*rank);
		} else {
		    Claset("F", 1, nrhs, Zero, Zero, &B[i + ldb], ldb);
		}
	    }
	    iwork = il + m * ldwork;
//Multiply B by right singular vectors of L in WORK(IL)
//(CWorkspace: need M*M+2*M, prefer M*M+M+M*NRHS)
//(RWorkspace: none)
	    if (lwork >= ldb * nrhs + iwork - 1 && nrhs > 1) {
		Cgemm("C", "N", m, nrhs, m, One, &work[il], ldwork, &B[0], ldb, Zero, &work[iwork], ldb);
		Clacpy("G", m, nrhs, &work[iwork], ldb, &B[0], ldb);
	    } else if (nrhs > 1) {
		chunk = (lwork - iwork + 1) / m;
		for (i = 1; i <= nrhs; i = i + chunk) {
		    bl = min(nrhs - i + 1, chunk);
		    Cgemm("C", "N", m, bl, m, One, &work[il], ldwork, &B[i * ldb + 1], ldb, Zero, &work[iwork], m);
		    Clacpy("G", m, bl, &work[iwork], m, &B[i * ldb + 1], ldb);
		}
	    } else {
		Cgemv("C", m, m, One, &work[il], ldwork, &B[ldb + 1], 1, Zero, &work[iwork], 1);
		Ccopy(m, &work[iwork], 1, &B[ldb + 1], 1);
	    }
//Zero out below first M rows of B
	    Claset("F", n - m, nrhs, Zero, Zero, &B[m + 1 + ldb], ldb);
	    iwork = itau + m;
//Multiply transpose(Q) by B
//(CWorkspace: need M+NRHS, prefer M+NHRS*NB)
//(RWorkspace: none)
	    Cunmlq("L", "C", n, nrhs, m, &A[0], lda, &work[itau], &B[0], ldb, &work[iwork], lwork - iwork + 1, info);
	} else {
//Path 2 - remaining underdetermined cases
	    ie = 1;
	    itauq = 1;
	    itaup = itauq + m;
	    iwork = itaup + m;
//Bidiagonalize A
//(CWorkspace: need 3*M, prefer 2*M+(M+N)*NB)
//(RWorkspace: need N)
	    Cgebrd(m, n, &A[0], lda, &s[1], &rwork[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, info);
//Multiply B by transpose of left bidiagonalizing vectors
//(CWorkspace: need 2*M+NRHS, prefer 2*M+NRHS*NB)
//(RWorkspace: none)
	    Cunmbr("Q", "L", "C", m, nrhs, n, &A[0], lda, &work[itauq], &B[0], ldb, &work[iwork], lwork - iwork + 1, info);
//Generate right bidiagonalizing vectors in A
//(CWorkspace: need 3*M, prefer 2*M+M*NB)
//(RWorkspace: none)
	    Cungbr("P", m, n, m, &A[0], lda, &work[itaup], &work[iwork], lwork - iwork + 1, info);
	    irwork = ie + m;
//Perform bidiagonal QR iteration,
//   computing right singular vectors of A in A and
//   multiplying B by transpose of left singular vectors
//(CWorkspace: none)
//(RWorkspace: need BDSPAC)
	    Cbdsqr("L", m, n, 0, nrhs, &s[1], &rwork[ie], &A[0], lda, &vdum, 1, &B[0], ldb, &rwork[irwork], info);
	    if (*info != 0) {
		goto L70;
	    }
//Multiply B by reciprocals of singular values
	    mtemp1 = rcond * s[1];
	    thr = max(mtemp1, sfmin);
	    if (rcond < Zero) {
		mtemp1 = eps * s[1];
		thr = max(mtemp1, sfmin);
	    }
	    *rank = 0;
	    for (i = 0; i < m; i++) {
		if (s[i] > thr) {
		    CRrscl(nrhs, s[i], &B[i + ldb], ldb);
		    ++(*rank);
		} else {
		    Claset("F", 1, nrhs, Zero, Zero, &B[i + ldb], ldb);
		}
	    }
//Multiply B by right singular vectors of A
//(CWorkspace: need N, prefer N*NRHS)
//(RWorkspace: none)
	    if (lwork >= ldb * nrhs && nrhs > 1) {
		Cgemm("C", "N", n, nrhs, m, One, &A[0], lda, &B[0], ldb, Zero, &work[0], ldb);
		Clacpy("G", n, nrhs, &work[0], ldb, &B[0], ldb);
	    } else if (nrhs > 1) {
		chunk = lwork / n;
		for (i = 1; i <= nrhs; i = i + chunk) {
		    bl = min(nrhs - i + 1, chunk);
		    Cgemm("C", "N", n, bl, m, One, &A[0], lda, &B[i * ldb + 1], ldb, Zero, &work[0], n);
		    Clacpy("F", n, bl, &work[0], n, &B[i * ldb + 1], ldb);

		}
	    } else {
		Cgemv("C", m, n, One, &A[0], lda, &B[0], 1, Zero, &work[0], 1);
		Ccopy(n, &work[0], 1, &B[0], 1);
	    }
	}
    }
//Undo scaling
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
  L70:
    work[1] = maxwrk;
    return;
}
