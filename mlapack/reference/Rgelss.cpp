/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgelss.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rgelss(INTEGER m, INTEGER n, INTEGER nrhs, REAL * A, INTEGER lda, REAL * B, INTEGER ldb, REAL * s, REAL rcond, INTEGER * rank, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i, bl, ie, il, mm;
    REAL eps, thr, anrm, bnrm;
    INTEGER itau;
    REAL vdum;
    INTEGER iascl, ibscl;
    INTEGER chunk;
    REAL sfmin;
    INTEGER minmn;
    INTEGER maxmn, itaup, itauq, mnthr = 0, iwork;
    INTEGER bdspac;
    REAL bignum;
    INTEGER ldwork;
    INTEGER minwrk, maxwrk;
    REAL smlnum;
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
//  NB refers to the optimal block size for the immediately
//  following subroutine, as returned by ILAENV.)
    if (*info == 0) {
	minwrk = 0;
	maxwrk = 0;
	if (minmn > 0) {
	    mm = m;
	    mnthr = iMlaenv(6, "Rgelss", " ", m, n, nrhs, -1);
	    if (m >= n && m >= mnthr) {
//Path 1a - overdetermined, with many more rows than
//          columns
		mm = n;
		maxwrk = max(maxwrk, n + n * iMlaenv(1, "Rgeqrf", " ", m, n, -1, -1));
		maxwrk = max(maxwrk, n + nrhs * iMlaenv(1, "Rormqr", "lt", m, nrhs, n, -1));
	    }
	    if (m >= n) {
//Path 1 - overdetermined or exactly determined
//Compute workspace needed for DBDSQR
		bdspac = max((INTEGER) 1, n * 5);
		maxwrk = max(maxwrk, n * 3 + (mm + n) * iMlaenv(1, "Rgebrd", " ", mm, n, -1, -1));
		maxwrk = max(maxwrk, n * 3 + nrhs * iMlaenv(1, "Rormbr", "qlt", mm, nrhs, n, -1));
		maxwrk = max(maxwrk, n * 3 + (n - 1) * iMlaenv(1, "Rormgbr", "p", n, n, n, -1));
		maxwrk = max(maxwrk, bdspac);
		maxwrk = max(maxwrk, n * nrhs);
		minwrk = max(max(n * 3 + mm, n * 3 + nrhs), bdspac);
		maxwrk = max(minwrk, maxwrk);
	    }
	    if (n > m) {
//Compute workspace needed for DBDSQR
		bdspac = max((INTEGER) 1, m * 5);
		minwrk = max(max(m * 3 + nrhs, m * 3 + 1), bdspac);
		if (n >= mnthr) {
//Path 2a - underdetermined, with many more columns
//than rows
		    maxwrk = m + m * iMlaenv(1, "Rgelqf", " ", m, n, -1, -1);
		    maxwrk = max(maxwrk, m * m + (m << 2) + (m << 1) * iMlaenv(1, "Rgebrd", " ", m, m, -1, -1));
		    maxwrk = max(maxwrk, m * m + (m << 2) + nrhs * iMlaenv(1, "Rormbr", "qlt", m, nrhs, m, -1));
		    maxwrk = max(maxwrk, m * m + (m << 2) + (m - 1) * iMlaenv(1, "Rorgbr", "p", m, m, m, -1));
		    maxwrk = max(maxwrk, m * m + m + bdspac);
		    if (nrhs > 1) {
			maxwrk = max(maxwrk, m * m + m + m * nrhs);
		    } else {
			maxwrk = max(maxwrk, m * m + (m << 1));
		    }
		    maxwrk = max(maxwrk, m + nrhs * iMlaenv(1, "Rormlq", "lt", n, nrhs, m, -1));
		} else {
//Path 2 - underdetermined
		    maxwrk = m * 3 + (n + m) * iMlaenv(1, "Rgebrd", " ", m, n, -1, -1);
		    maxwrk = max(maxwrk, m * 3 + nrhs * iMlaenv(1, "Rormbr", "qlt", m, nrhs, m, -1));
		    maxwrk = max(maxwrk, m * 3 + m * iMlaenv(1, "Rorg" "br", "p", m, n, m, -1));
		    maxwrk = max(maxwrk, bdspac);
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
	Mxerbla("Rgelss", -(*info));
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
	Rlaset("F", minmn, 1, Zero, Zero, &s[1], 1);
	*rank = 0;
	goto L70;
    }
//Scale B if max element outside range [SMLNUM,BIGNUM]
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
//(Workspace: need 2*N, prefer N+N*NB)
	    Rgeqrf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, info);
//Multiply B by transpose(Q)
//(Workspace: need N+NRHS, prefer N+NRHS*NB)
	    Rormqr("L", "T", m, nrhs, n, &A[0], lda, &work[itau], &B[0], ldb, &work[iwork], lwork - iwork + 1, info);
//Zero out below R
	    if (n > 1) {
		Rlaset("L", n - 1, n - 1, Zero, Zero, &A[lda + 2], lda);
	    }
	}
	ie = 1;
	itauq = ie + n;
	itaup = itauq + n;
	iwork = itaup + n;
//Bidiagonalize R in A
//(Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB)
	Rgebrd(mm, n, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, info);
//Multiply B by transpose of left bidiagonalizing vectors of R
//(Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB)
	Rormbr("Q", "L", "T", mm, nrhs, n, &A[0], lda, &work[itauq], &B[0], ldb, &work[iwork], lwork - iwork + 1, info);
//Generate right bidiagonalizing vectors of R in A
//(Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
	Rorgbr("P", n, n, n, &A[0], lda, &work[itaup], &work[iwork], lwork - iwork + 1, info);
	iwork = ie + n;
//Perform bidiagonal QR iteration
//  multiply B by transpose of left singular vectors
//  compute right singular vectors in A
//(Workspace: need BDSPAC)
	Rbdsqr("U", n, n, 0, nrhs, &s[1], &work[ie], &A[0], lda, &vdum, 1, &B[0], ldb, &work[iwork], info);
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
	for (i = 0; i < n; i++) {
	    if (s[i] > thr) {
		Rrscl(nrhs, s[i], &B[i + ldb], ldb);
		++(*rank);
	    } else {
		Rlaset("F", 1, nrhs, Zero, Zero, &B[i + ldb], ldb);
	    }
	}
//Multiply B by right singular vectors
//(Workspace: need N, prefer N*NRHS)
	if (lwork >= ldb * nrhs && nrhs > 1) {
	    Rgemm("T", "N", n, nrhs, n, One, &A[0], lda, &B[0], ldb, Zero, &work[0], ldb);
	    Rlacpy("G", n, nrhs, &work[0], ldb, &B[0], ldb);
	} else if (nrhs > 1) {
	    chunk = lwork / n;
	    for (i = 1; i <= nrhs; i += chunk) {
		bl = min(nrhs - i + 1, chunk);
		Rgemm("T", "N", n, bl, n, One, &A[0], lda, &B[i * ldb + 1], ldb, Zero, &work[0], n);
		Rlacpy("G", n, bl, &work[0], n, &B[i * ldb + 1], ldb);
	    }
	} else {
	    Rgemv("T", n, n, One, &A[0], lda, &B[0], 1, Zero, &work[0], 1);
	    Rcopy(n, &work[0], 1, &B[0], 1);
	}
    } else {
	if (n >= mnthr && lwork >= (m << 2) + m * m + max(max(max(m, (m << 1) - 4), nrhs), n - m * 3)) {
//Path 2a - underdetermined, with many more columns than rows
//and sufficient workspace for an efficient algorithm
	    ldwork = m;
	    if (lwork >= max((m << 2) + m * lda + max(max(max(m, (m << 1) - 4), nrhs), n - m * 3), m * lda + m + m * nrhs)) {
		ldwork = lda;
	    }
	    itau = 1;
	    iwork = m + 1;
//Compute A=L*Q
//(Workspace: need 2*M, prefer M+M*NB)
	    Rgelqf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, info);
	    il = iwork;
//Copy L to WORK(IL), zeroing out above it
	    Rlacpy("L", m, m, &A[0], lda, &work[il], ldwork);
	    Rlaset("U", m - 1, m - 1, Zero, Zero, &work[il + ldwork], ldwork);
	    ie = il + ldwork * m;
	    itauq = ie + m;
	    itaup = itauq + m;
	    iwork = itaup + m;
//Bidiagonalize L in WORK(IL)
//(Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB)
	    Rgebrd(m, m, &work[il], ldwork, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, info);
//Multiply B by transpose of left bidiagonalizing vectors of L
//(Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)
	    Rormbr("Q", "L", "T", m, nrhs, m, &work[il], ldwork, &work[itauq], &B[0], ldb, &work[iwork], lwork - iwork + 1, info);
//Generate right bidiagonalizing vectors of R in WORK(IL)
//(Workspace: need M*M+5*M-1, prefer M*M+4*M+(M-1)*NB)
	    Rorgbr("P", m, m, m, &work[il], ldwork, &work[itaup], &work[iwork], lwork - iwork + 1, info);
	    iwork = ie + m;
//Perform bidiagonal QR iteration,
//   computing right singular vectors of L in WORK(IL) and
//   multiplying B by transpose of left singular vectors
//(Workspace: need M*M+M+BDSPAC)
	    Rbdsqr("U", m, m, 0, nrhs, &s[1], &work[ie], &work[il], ldwork, &A[0], lda, &B[0], ldb, &work[iwork], info);
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
		    Rrscl(nrhs, s[i], &B[i + ldb], ldb);
		    ++(*rank);
		} else {
		    Rlaset("F", 1, nrhs, Zero, Zero, &B[i + ldb], ldb);
		}
	    }
	    iwork = ie;
//Multiply B by right singular vectors of L in WORK(IL)
//(Workspace: need M*M+2*M, prefer M*M+M+M*NRHS)
	    if (lwork >= ldb * nrhs + iwork - 1 && nrhs > 1) {
		Rgemm("T", "N", m, nrhs, m, One, &work[il], ldwork, &B[0], ldb, Zero, &work[iwork], ldb);
		Rlacpy("G", m, nrhs, &work[iwork], ldb, &B[0], ldb);
	    } else if (nrhs > 1) {
		chunk = (lwork - iwork + 1) / m;
		for (i = 1; i <= nrhs; i = i + chunk) {
		    bl = min(nrhs - i + 1, chunk);
		    Rgemm("T", "N", m, bl, m, One, &work[il], ldwork, &B[i * ldb + 1], ldb, Zero, &work[iwork], m);
		    Rlacpy("G", m, bl, &work[iwork], m, &B[i * ldb + 1], ldb);
		}
	    } else {
		Rgemv("T", m, m, One, &work[il], ldwork, &B[ldb + 1], 1, Zero, &work[iwork], 1);
		Rcopy(m, &work[iwork], 1, &B[ldb + 1], 1);
	    }
//Zero out below first M rows of B
	    Rlaset("F", n - m, nrhs, Zero, Zero, &B[m + 1 + ldb], ldb);
	    iwork = itau + m;
//Multiply transpose(Q) by B
//(Workspace: need M+NRHS, prefer M+NRHS*NB)
	    Rormlq("L", "T", n, nrhs, m, &A[0], lda, &work[itau], &B[0], ldb, &work[iwork], lwork - iwork + 1, info);
	} else {
//Path 2 - remaining underdetermined cases
	    ie = 1;
	    itauq = ie + m;
	    itaup = itauq + m;
	    iwork = itaup + m;
//Bidiagonalize A
//(Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
	    Rgebrd(m, n, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, info);
//Multiply B by transpose of left bidiagonalizing vectors
//(Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB)
	    Rormbr("Q", "L", "T", m, nrhs, n, &A[0], lda, &work[itauq]
		   , &B[0], ldb, &work[iwork], lwork - iwork + 1, info);
//Generate right bidiagonalizing vectors in A
//(Workspace: need 4*M, prefer 3*M+M*NB)
	    Rorgbr("P", m, n, m, &A[0], lda, &work[itaup], &work[iwork], lwork - iwork + 1, info);
	    iwork = ie + m;
//Perform bidiagonal QR iteration,
//   computing right singular vectors of A in A and
//   multiplying B by transpose of left singular vectors
//(Workspace: need BDSPAC)
	    Rbdsqr("L", m, n, 0, nrhs, &s[1], &work[ie], &A[0], lda, &vdum, 1, &B[0], ldb, &work[iwork], info);
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
		    Rrscl(nrhs, s[i], &B[i + ldb], ldb);
		    ++(*rank);
		} else {
		    Rlaset("F", 1, nrhs, Zero, Zero, &B[i + ldb], ldb);
		}
	    }
//Multiply B by right singular vectors of A
//(Workspace: need N, prefer N*NRHS)
	    if (lwork >= ldb * nrhs && nrhs > 1) {
		Rgemm("T", "N", n, nrhs, m, One, &A[0], lda, &B[0], ldb, Zero, &work[0], ldb);
		Rlacpy("F", n, nrhs, &work[0], ldb, &B[0], ldb);
	    } else if (nrhs > 1) {
		chunk = lwork / n;
		for (i = 1; i <= nrhs; i += chunk) {
		    bl = min(nrhs - i + 1, chunk);
		    Rgemm("T", "N", n, bl, m, One, &A[0], lda, &B[i * ldb + 1], ldb, Zero, &work[0], n);
		    Rlacpy("F", n, bl, &work[0], n, &B[i * ldb + 1], ldb);
		}
	    } else {
		Rgemv("T", m, n, One, &A[0], lda, &B[0], 1, Zero, &work[0], 1);
		Rcopy(n, &work[0], 1, &B[0], 1);
	    }
	}
    }
//Undo scaling
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
  L70:
    work[1] = maxwrk;
    return;
}
