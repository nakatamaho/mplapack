/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgesdd.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rgesdd(const char *jobz, INTEGER m, INTEGER n, REAL *
	    A, INTEGER lda, REAL * s, REAL * u, INTEGER ldu, REAL * vt, INTEGER ldvt, REAL * work, INTEGER lwork, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, ie, il = 0, ir = 0, iu, blk;
    REAL dum, eps;
    INTEGER ivt, iscl;
    REAL anrm;
    INTEGER idum, ierr, itau;
    INTEGER chunk = 0, minmn, wrkbl, itaup, itauq, mnthr = 0;
    LOGICAL wntqa;
    INTEGER nwork;
    LOGICAL wntqn, wntqo, wntqs;
    INTEGER bdspac = 0, minwrk, maxwrk, ldwrkr = 0, ldwrkl, ldwrku, ldwkvt;
    REAL bignum;
    REAL smlnum;
    LOGICAL wntqas, lquery;
    REAL Zero = 0.0, One = 1.0;

//Test the input arguments
    *info = 0;
    minmn = min(m, n);
    wntqa = Mlsame(jobz, "A");
    wntqs = Mlsame(jobz, "S");
    wntqas = wntqa || wntqs;
    wntqo = Mlsame(jobz, "O");
    wntqn = Mlsame(jobz, "N");
    lquery = lwork == -1;
    if (!(wntqa || wntqs || wntqo || wntqn)) {
	*info = -1;
    } else if (m < 0) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -5;
    } else if (ldu < 1 || (wntqas && ldu < m) || (wntqo && m < n && ldu < m)) {
	*info = -8;
    } else if (ldvt < 1 || (wntqa && ldvt < n) || (wntqs && ldvt < minmn) || (wntqo && m >= n && ldvt < n)) {
	*info = -10;
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
	if (m >= n && minmn > 0) {
//Compute space needed for DBDSDC
	    mnthr = (minmn * 11.0 / 6.);
	    if (wntqn) {
		bdspac = n * 7;
	    } else {
		bdspac = n * 3 * n + (n << 2);
	    }
	    if (m >= mnthr) {
		if (wntqn) {
//Path 1 (M much larger than N, JOBZ='N')
		    wrkbl = n + n * iMlaenv(1, "RGEQRF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, n * 3 + (n << 1) * iMlaenv(1, "RGEBRD", " ", n, n, -1, -1));
		    maxwrk = max(wrkbl, bdspac + n);
		    minwrk = bdspac + n;
		} else if (wntqo) {
//Path 2 (M much larger than N, JOBZ='O')
		    wrkbl = n + n * iMlaenv(1, "RGEQRF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, n + n * iMlaenv(1, "RORGQR", " ", m, n, n, -1));
		    wrkbl = max(wrkbl, n * 3 + (n << 1) * iMlaenv(1, "RGEBRD", " ", n, n, -1, -1));
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORMBR", "QLN", n, n, n, -1));
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORMBR", "PRT", n, n, n, -1));
		    wrkbl = max(wrkbl, bdspac + n * 3);
		    maxwrk = wrkbl + (n << 1) * n;
		    minwrk = bdspac + (n << 1) * n + n * 3;
		} else if (wntqs) {
//Path 3 (M much larger than N, JOBZ='S')
		    wrkbl = n + n * iMlaenv(1, "RGEQRF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, n + n * iMlaenv(1, "RORGQR", " ", m, n, n, -1));
		    wrkbl = max(wrkbl, n * 3 + (n << 1) * iMlaenv(1, "RGEBRD", " ", n, n, -1, -1));
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORMBR", "QLN", n, n, n, -1));
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORMBR", "PRT", n, n, n, -1));
		    wrkbl = max(wrkbl, bdspac + n * 3);
		    maxwrk = wrkbl + n * n;
		    minwrk = bdspac + n * n + n * 3;
		} else if (wntqa) {
//Path 4 (M much larger than N, JOBZ='A')
		    wrkbl = n + n * iMlaenv(1, "RGEQRF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, n + m * iMlaenv(1, "RORGQR", " ", m, m, n, -1));
		    wrkbl = max(wrkbl, n * 3 + (n << 1) * iMlaenv(1, "RGEBRD", " ", n, n, -1, -1));
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORMBR", "QLN", n, n, n, -1));
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORMBR", "PRT", n, n, n, -1));
		    wrkbl = max(wrkbl, bdspac + n * 3);
		    maxwrk = wrkbl + n * n;
		    minwrk = bdspac + n * n + n * 3;
		}
	    } else {
//Path 5 (M at least N, but not much larger)
		wrkbl = n * 3 + (m + n) * iMlaenv(1, "RGEBRD", " ", m, n, -1, -1);
		if (wntqn) {
		    maxwrk = max(wrkbl, bdspac + n * 3);
		    minwrk = n * 3 + max(m, bdspac);
		} else if (wntqo) {
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORMBR", "QLN", m, n, n, -1));
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORMBR", "PRT", n, n, n, -1));
		    wrkbl = max(wrkbl, bdspac + n * 3);
		    maxwrk = wrkbl + m * n;
		    minwrk = n * 3 + max(m, n * n + bdspac);
		} else if (wntqs) {
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORMBR", "QLN", m, n, n, -1));
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORMBR", "PRT", n, n, n, -1));
		    maxwrk = max(wrkbl, bdspac + n * 3);
		    minwrk = n * 3 + max(m, bdspac);
		} else if (wntqa) {
		    wrkbl = max(wrkbl, n * 3 + m * iMlaenv(1, "RORMBR", "QLN", m, m, n, -1));
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORMBR", "PRT", n, n, n, -1));
		    maxwrk = max(maxwrk, bdspac + n * 3);
		    minwrk = n * 3 + max(m, bdspac);
		}
	    }
	} else if (minmn > 0) {
//Compute space needed for DBDSDC
	    mnthr = (minmn * 11.0 / 6.0);
	    if (wntqn) {
		bdspac = m * 7;
	    } else {
		bdspac = m * 3 * m + (m << 2);
	    }
	    if (n >= mnthr) {
		if (wntqn) {
//Path 1t (N much larger than M, JOBZ='N')
		    wrkbl = m + m * iMlaenv(1, "RGELQF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, m * 3 + (m << 1) * iMlaenv(1, "RGEBRD", " ", m, m, -1, -1));
		    maxwrk = max(wrkbl, bdspac + m);
		    minwrk = bdspac + m;
		} else if (wntqo) {
//Path 2t (N much larger than M, JOBZ='O')
		    wrkbl = m + m * iMlaenv(1, "RGELQF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, m + m * iMlaenv(1, "RORGLQ", " ", m, n, m, -1));
		    wrkbl = max(wrkbl, m * 3 + (m << 1) * iMlaenv(1, "RGEBRD", " ", m, m, -1, -1));
		    wrkbl = max(wrkbl, m * 3 + m * iMlaenv(1, "RORMBR", "QLN", m, m, m, -1));
		    wrkbl = max(wrkbl, m * 3 + m * iMlaenv(1, "RORMBR", "PRT", m, m, m, -1));
		    wrkbl = max(wrkbl, bdspac + m * 3);
		    maxwrk = wrkbl + (m << 1) * m;
		    minwrk = bdspac + (m << 1) * m + m * 3;
		} else if (wntqs) {
//Path 3t (N much larger than M, JOBZ='S')
		    wrkbl = m + m * iMlaenv(1, "RGELQF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, m + m * iMlaenv(1, "RORGLQ", " ", m, n, m, -1));
		    wrkbl = max(wrkbl, m * 3 + (m << 1) * iMlaenv(1, "RGEBRD", " ", m, m, -1, -1));
		    wrkbl = max(wrkbl, m * 3 + m * iMlaenv(1, "RORMBR", "QLN", m, m, m, -1));
		    wrkbl = max(wrkbl, m * 3 + m * iMlaenv(1, "RORMBR", "PRT", m, m, m, -1));
		    wrkbl = max(wrkbl, bdspac + m * 3);
		    maxwrk = wrkbl + m * m;
		    minwrk = bdspac + m * m + m * 3;
		} else if (wntqa) {
//Path 4t (N much larger than M, JOBZ='A')
		    wrkbl = m + m * iMlaenv(1, "RGELQF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, m + n * iMlaenv(1, "RORGLQ", " ", n, n, m, -1));
		    wrkbl = max(wrkbl, m * 3 + (m << 1) * iMlaenv(1, "RGEBRD", " ", m, m, -1, -1));
		    wrkbl = max(wrkbl, m * 3 + m * iMlaenv(1, "RORMBR", "QLN", m, m, m, -1));
		    wrkbl = max(wrkbl, m * 3 + m * iMlaenv(1, "RORMBR", "PRT", m, m, m, -1));
		    wrkbl = max(wrkbl, bdspac + m * 3);
		    maxwrk = wrkbl + m * m;
		    minwrk = bdspac + m * m + m * 3;
		}
	    } else {
//Path 5t (N greater than M, but not much larger)
		wrkbl = m * 3 + (m + n) * iMlaenv(1, "RGEBRD", " ", m, n, -1, -1);
		if (wntqn) {
		    maxwrk = max(wrkbl, bdspac + m * 3);
		    minwrk = m * 3 + max(n, bdspac);
		} else if (wntqo) {
		    wrkbl = max(wrkbl, m * 3 + m * iMlaenv(1, "RORMBR", "QLN", m, m, n, -1));
		    wrkbl = max(wrkbl, m * 3 + m * iMlaenv(1, "RORMBR", "PRT", m, n, m, -1));
		    wrkbl = max(wrkbl, bdspac + m * 3);
		    maxwrk = wrkbl + m * n;
		    minwrk = m * 3 + max(n, m * m + bdspac);
		} else if (wntqs) {
		    wrkbl = max(wrkbl, m * 3 + m * iMlaenv(1, "RORMBR", "QLN", m, m, n, -1));
		    wrkbl = max(wrkbl, m * 3 + m * iMlaenv(1, "RORMBR", "PRT", m, n, m, -1));
		    maxwrk = max(wrkbl, bdspac + m * 3);
		    minwrk = m * 3 + max(n, bdspac);
		} else if (wntqa) {
		    wrkbl = max(wrkbl, m * 3 + m * iMlaenv(1, "RORMBR", "QLN", m, m, n, -1));
		    wrkbl = max(wrkbl, m * 3 + m * iMlaenv(1, "RORMBR", "PRT", n, n, m, -1));
		    maxwrk = max(wrkbl, bdspac + m * 3);
		    minwrk = m * 3 + max(n, bdspac);
		}
	    }
	}
	maxwrk = max(maxwrk, minwrk);
	work[1] = maxwrk;
	if (lwork < minwrk && !lquery) {
	    *info = -12;
	}
    }
    if (*info != 0) {
	Mxerbla("Rgesdd", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (m == 0 || n == 0) {
	return;
    }
//Get machine constants
    eps = Rlamch("P");
    smlnum = sqrt(Rlamch("S")) / eps;
    bignum = One / smlnum;
//Scale A if max element outside range [SMLNUM,BIGNUM]
    anrm = Rlange("M", m, n, &A[0], lda, &dum);
    iscl = 0;
    if (anrm > Zero && anrm < smlnum) {
	iscl = 0;
	Rlascl("G", 0, 0, anrm, smlnum, m, n, &A[0], lda, &ierr);
    } else if (anrm > bignum) {
	iscl = 0;
	Rlascl("G", 0, 0, anrm, bignum, m, n, &A[0], lda, &ierr);
    }
    if (m >= n) {
//A has at least as many rows as columns. If A has sufficiently
//more rows than columns, first reduce using the QR
//decomposition (if sufficient workspace available)
	if (m >= mnthr) {
	    if (wntqn) {
//Path 1 (M much larger than N, JOBZ='N')
//No singular vectors to be computed
		itau = 1;
		nwork = itau + n;
//Compute A=Q*R
//(Workspace: need 2*N, prefer N+N*NB)
		Rgeqrf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
//Zero out below R
		Rlaset("L", n - 1, n - 1, Zero, Zero, &A[lda + 2], lda);
		ie = 1;
		itauq = ie + n;
		itaup = itauq + n;
		nwork = itaup + n;
//Bidiagonalize R in A
//(Workspace: need 4*N, prefer 3*N+2*N*NB)
		Rgebrd(n, n, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
		nwork = ie + n;
//Perform bidiagonal SVD, computing singular values only
//(Workspace: need N+BDSPAC)
		Rbdsdc("U", "N", n, &s[1], &work[ie], &dum, 1, &dum, 1, &dum, &idum, &work[nwork], &iwork[1], info);
	    } else if (wntqo) {
//Path 2 (M much larger than N, JOBZ = 'O')
//N left singular vectors to be overwritten on A and
//N right singular vectors to be computed in VT
		ir = 1;
//WORK(IR) is LDWRKR by N
		if (lwork >= lda * n + n * n + n * 3 + bdspac) {
		    ldwrkr = lda;
		} else {
		    ldwrkr = (lwork - n * n - n * 3 - bdspac) / n;
		}
		itau = ir + ldwrkr * n;
		nwork = itau + n;
//Compute A=Q*R
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
		Rgeqrf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
//Copy R to WORK(IR), zeroing out below it
		Rlacpy("U", n, n, &A[0], lda, &work[ir], ldwrkr);
		Rlaset("L", n - 1, n - 1, Zero, Zero, &work[ir + 1], ldwrkr);
//Generate Q in A
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
		Rorgqr(m, n, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
		ie = itau;
		itauq = ie + n;
		itaup = itauq + n;
		nwork = itaup + n;
//Bidiagonalize R in VT, copying result to WORK(IR)
//(Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
		Rgebrd(n, n, &work[ir], ldwrkr, &s[1], &work[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
//WORK(IU) is N by N
		iu = nwork;
		nwork = iu + n * n;
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in WORK(IU) and computing right
//singular vectors of bidiagonal matrix in VT
//(Workspace: need N+N*N+BDSPAC)
		Rbdsdc("U", "I", n, &s[1], &work[ie], &work[iu], n, &vt[0], ldvt, &dum, &idum, &work[nwork], &iwork[1], info);
//Overwrite WORK(IU) by left singular vectors of R
//and VT by right singular vectors of R
//(Workspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB)
		Rormbr("Q", "L", "N", n, n, n, &work[ir], ldwrkr, &work[itauq], &work[iu], n, &work[nwork], lwork - nwork + 1, &ierr);
		Rormbr("P", "R", "T", n, n, n, &work[ir], ldwrkr, &work[itaup], &vt[0], ldvt, &work[nwork], lwork - nwork + 1, &ierr);
//Multiply Q in A by left singular vectors of R in
//WORK(IU), storing result in WORK(IR) and copying to A
//(Workspace: need 2*N*N, prefer N*N+M*N)
		for (i = 1; i <= m; i = i + ldwrkr) {
		    chunk = min(m - i + 1, ldwrkr);
		    Rgemm("N", "N", chunk, n, n, One, &A[i + lda], lda, &work[iu], n, Zero, &work[ir], ldwrkr);
		    Rlacpy("F", chunk, n, &work[ir], ldwrkr, &A[i + lda], lda);
		}
	    } else if (wntqs) {
//Path 3 (M much larger than N, JOBZ='S')
//N left singular vectors to be computed in U and
//N right singular vectors to be computed in VT
		ir = 1;
//WORK(IR) is N by N
		ldwrkr = n;
		itau = ir + ldwrkr * n;
		nwork = itau + n;
//Compute A=Q*R
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
		Rgeqrf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
//Copy R to WORK(IR), zeroing out below it
		Rlacpy("U", n, n, &A[0], lda, &work[ir], ldwrkr);
		Rlaset("L", n - 1, n - 1, Zero, Zero, &work[ir + 1], ldwrkr);
//Generate Q in A
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
		Rorgqr(m, n, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
		ie = itau;
		itauq = ie + n;
		itaup = itauq + n;
		nwork = itaup + n;
//Bidiagonalize R in WORK(IR)
//(Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
		Rgebrd(n, n, &work[ir], ldwrkr, &s[1], &work[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagoal matrix in U and computing right singular
//vectors of bidiagonal matrix in VT
//(Workspace: need N+BDSPAC)
		Rbdsdc("U", "I", n, &s[1], &work[ie], &u[0], ldu, &vt[0], ldvt, &dum, &idum, &work[nwork], &iwork[1], info);
//Overwrite U by left singular vectors of R and VT
//by right singular vectors of R
//(Workspace: need N*N+3*N, prefer N*N+2*N+N*NB)
		Rormbr("Q", "L", "N", n, n, n, &work[ir], ldwrkr, &work[itauq], &u[0], ldu, &work[nwork], lwork - nwork + 1, &ierr);
		Rormbr("P", "R", "T", n, n, n, &work[ir], ldwrkr, &work[itaup], &vt[0], ldvt, &work[nwork], lwork - nwork + 1, &ierr);
//Multiply Q in A by left singular vectors of R in
//WORK(IR), storing result in U
//(Workspace: need N*N)
		Rlacpy("F", n, n, &u[0], ldu, &work[ir], ldwrkr);
		Rgemm("N", "N", m, n, n, One, &A[0], lda, &work[ir], ldwrkr, Zero, &u[0], ldu);
	    } else if (wntqa) {
//Path 4 (M much larger than N, JOBZ='A')
//M left singular vectors to be computed in U and
//N right singular vectors to be computed in VT
		iu = 1;
//WORK(IU) is N by N
		ldwrku = n;
		itau = iu + ldwrku * n;
		nwork = itau + n;
//Compute A=Q*R, copying result to U
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
		Rgeqrf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
		Rlacpy("L", m, n, &A[0], lda, &u[0], ldu);
//Generate Q in U
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
		Rorgqr(m, m, n, &u[0], ldu, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
//Produce R in A, zeroing out other entries
		Rlaset("L", n - 1, n - 1, Zero, Zero, &A[lda + 2], lda);
		ie = itau;
		itauq = ie + n;
		itaup = itauq + n;
		nwork = itaup + n;
//Bidiagonalize R in A
//(Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
		Rgebrd(n, n, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in WORK(IU) and computing right
//singular vectors of bidiagonal matrix in VT
//(Workspace: need N+N*N+BDSPAC)
		Rbdsdc("U", "I", n, &s[1], &work[ie], &work[iu], n, &vt[0], ldvt, &dum, &idum, &work[nwork], &iwork[1], info);
//Overwrite WORK(IU) by left singular vectors of R and VT
//by right singular vectors of R
//(Workspace: need N*N+3*N, prefer N*N+2*N+N*NB)
		Rormbr("Q", "L", "N", n, n, n, &A[0], lda, &work[itauq], &work[iu], ldwrku, &work[nwork], lwork - nwork + 1, &ierr);
		Rormbr("P", "R", "T", n, n, n, &A[0], lda, &work[itaup], &vt[0], ldvt, &work[nwork], lwork - nwork + 1, &ierr);
//Multiply Q in U by left singular vectors of R in
//WORK(IU), storing result in A
//(Workspace: need N*N)
		Rgemm("N", "N", m, n, n, One, &u[0], ldu, &work[iu], ldwrku, Zero, &A[0], lda);
//Copy left singular vectors of A from A to U
		Rlacpy("F", m, n, &A[0], lda, &u[0], ldu);
	    }
	} else {
//M .LT. MNTHR
//Path 5 (M at least N, but not much larger)
//Reduce to bidiagonal form without QR decomposition
	    ie = 1;
	    itauq = ie + n;
	    itaup = itauq + n;
	    nwork = itaup + n;
//Bidiagonalize A
//(Workspace: need 3*N+M, prefer 3*N+(M+N)*NB)
	    Rgebrd(m, n, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
	    if (wntqn) {
//Perform bidiagonal SVD, only computing singular values
//(Workspace: need N+BDSPAC)
		Rbdsdc("U", "N", n, &s[1], &work[ie], &dum, 1, &dum, 1, &dum, &idum, &work[nwork], &iwork[1], info);
	    } else if (wntqo) {
		iu = nwork;
		if (lwork >= m * n + n * 3 + bdspac) {
//WORK( IU ) is M by N
		    ldwrku = m;
		    nwork = iu + ldwrku * n;
		    Rlaset("F", m, n, Zero, Zero, &work[iu], ldwrku);
		} else {
//WORK( IU ) is N by N
		    ldwrku = n;
		    nwork = iu + ldwrku * n;
//WORK(IR) is LDWRKR by N
		    ir = nwork;
		    ldwrkr = (lwork - n * n - n * 3) / n;
		}
		nwork = iu + ldwrku * n;
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in WORK(IU) and computing right
//singular vectors of bidiagonal matrix in VT
//(Workspace: need N+N*N+BDSPAC)
		Rbdsdc("U", "I", n, &s[1], &work[ie], &work[iu], ldwrku, &vt[0], ldvt, &dum, &idum, &work[nwork], &iwork[1], info);
//Overwrite VT by right singular vectors of A
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
		Rormbr("P", "R", "T", n, n, n, &A[0], lda, &work[itaup], &vt[0], ldvt, &work[nwork], lwork - nwork + 1, &ierr);
		if (lwork >= m * n + n * 3 + bdspac) {
//Overwrite WORK(IU) by left singular vectors of A
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
		    Rormbr("Q", "L", "N", m, n, n, &A[0], lda, &work[itauq], &work[iu], ldwrku, &work[nwork], lwork - nwork + 1, &ierr);
//Copy left singular vectors of A from WORK(IU) to A
		    Rlacpy("F", m, n, &work[iu], ldwrku, &A[0], lda);
		} else {
//Generate Q in A
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
		    Rorgbr("Q", m, n, n, &A[0], lda, &work[itauq], &work[nwork], lwork - nwork + 1, &ierr);
//Multiply Q in A by left singular vectors of
//bidiagonal matrix in WORK(IU), storing result in
//WORK(IR) and copying to A
//(Workspace: need 2*N*N, prefer N*N+M*N)
		    for (i = 1; i <= m; i += ldwrkr) {
			chunk = min(m - i + 1, ldwrkr);
			Rgemm("N", "N", chunk, n, n, One, &A[i + lda], lda, &work[iu], ldwrku, Zero, &work[ir], ldwrkr);
			Rlacpy("F", chunk, n, &work[ir], ldwrkr, &A[i + lda], lda);
		    }
		}
	    } else if (wntqs) {
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in U and computing right singular
//vectors of bidiagonal matrix in VT
//(Workspace: need N+BDSPAC)
		Rlaset("F", m, n, Zero, Zero, &u[0], ldu);
		Rbdsdc("U", "I", n, &s[1], &work[ie], &u[0], ldu, &vt[0], ldvt, &dum, &idum, &work[nwork], &iwork[1], info);
//Overwrite U by left singular vectors of A and VT
//by right singular vectors of A
//(Workspace: need 3*N, prefer 2*N+N*NB)

		Rormbr("Q", "L", "N", m, n, n, &A[0], lda, &work[itauq], &u[0], ldu, &work[nwork], lwork - nwork + 1, &ierr);
		Rormbr("P", "R", "T", n, n, n, &A[0], lda, &work[itaup], &vt[0], ldvt, &work[nwork], lwork - nwork + 1, &ierr);
	    } else if (wntqa) {
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in U and computing right singular
//vectors of bidiagonal matrix in VT
//(Workspace: need N+BDSPAC)
		Rlaset("F", m, m, Zero, Zero, &u[0], ldu);
		Rbdsdc("U", "I", n, &s[1], &work[ie], &u[0], ldu, &vt[0], ldvt, &dum, &idum, &work[nwork], &iwork[1], info);
//Set the right corner of U to identity matrix
		if (m > n) {
		    Rlaset("F", m - n, m - n, Zero, One, &u[n + 1 + (n + 1) * ldu], ldu);
		}
//Overwrite U by left singular vectors of A and VT
//by right singular vectors of A
//(Workspace: need N*N+2*N+M, prefer N*N+2*N+M*NB)
		Rormbr("Q", "L", "N", m, m, n, &A[0], lda, &work[itauq], &u[0], ldu, &work[nwork], lwork - nwork + 1, &ierr);
		Rormbr("P", "R", "T", n, n, m, &A[0], lda, &work[itaup], &vt[0], ldvt, &work[nwork], lwork - nwork + 1, &ierr);
	    }
	}
    } else {
//A has more columns than rows. If A has sufficiently more
//columns than rows, first reduce using the LQ decomposition (if
//sufficient workspace available)
	if (n >= mnthr) {
	    if (wntqn) {
//Path 1t (N much larger than M, JOBZ='N')
//No singular vectors to be computed
		itau = 1;
		nwork = itau + m;
//Compute A=L*Q
//(Workspace: need 2*M, prefer M+M*NB)
		Rgelqf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
//Zero out above L
		Rlaset("U", m - 1, m - 1, Zero, Zero, &A[(lda << 1) + 1], lda);
		ie = 1;
		itauq = ie + m;
		itaup = itauq + m;
		nwork = itaup + m;
//Bidiagonalize L in A
//(Workspace: need 4*M, prefer 3*M+2*M*NB)
		Rgebrd(m, m, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
		nwork = ie + m;
//Perform bidiagonal SVD, computing singular values only
//(Workspace: need M+BDSPAC)
		Rbdsdc("U", "N", m, &s[1], &work[ie], &dum, 1, &dum, 1, &dum, &idum, &work[nwork], &iwork[1], info);
	    } else if (wntqo) {
//Path 2t (N much larger than M, JOBZ='O')
//M right singular vectors to be overwritten on A and
//M left singular vectors to be computed in U
		ivt = 1;
//IVT is M by M
		il = ivt + m * m;
		if (lwork >= m * n + m * m + m * 3 + bdspac) {
//WORK(IL) is M by N
		    ldwrkl = m;
		    chunk = n;
		} else {
		    ldwrkl = m;
		    chunk = (lwork - m * m) / m;
		}
		itau = il + ldwrkl * m;
		nwork = itau + m;
//Compute A=L*Q
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
		Rgelqf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
//Copy L to WORK(IL), zeroing about above it
		Rlacpy("L", m, m, &A[0], lda, &work[il], ldwrkl);
		Rlaset("U", m - 1, m - 1, Zero, Zero, &work[il + ldwrkl], ldwrkl);
//Generate Q in A
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
		Rorglq(m, n, m, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
		ie = itau;
		itauq = ie + m;
		itaup = itauq + m;
		nwork = itaup + m;
//Bidiagonalize L in WORK(IL)
//(Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
		Rgebrd(m, m, &work[il], ldwrkl, &s[1], &work[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in U, and computing right singular
//vectors of bidiagonal matrix in WORK(IVT)
//(Workspace: need M+M*M+BDSPAC)
		Rbdsdc("U", "I", m, &s[1], &work[ie], &u[0], ldu, &work[ivt], m, &dum, &idum, &work[nwork], &iwork[1], info);
//Overwrite U by left singular vectors of L and WORK(IVT)
//by right singular vectors of L
//(Workspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB)
		Rormbr("Q", "L", "N", m, m, m, &work[il], ldwrkl, &work[itauq], &u[0], ldu, &work[nwork], lwork - nwork + 1, &ierr);
		Rormbr("P", "R", "T", m, m, m, &work[il], ldwrkl, &work[itaup], &work[ivt], m, &work[nwork], lwork - nwork + 1, &ierr);
//Multiply right singular vectors of L in WORK(IVT) by Q
//in A, storing result in WORK(IL) and copying to A
//(Workspace: need 2*M*M, prefer M*M+M*N)
		for (i = 1; i <= n; i = i + chunk) {
		    blk = min(n - i + 1, chunk);
		    Rgemm("N", "N", m, blk, m, One, &work[ivt], m, &A[i * lda + 1], lda, Zero, &work[il], ldwrkl);
		    Rlacpy("F", m, blk, &work[il], ldwrkl, &A[i * lda + 1], lda);
		}
	    } else if (wntqs) {
//Path 3t (N much larger than M, JOBZ='S')
//M right singular vectors to be computed in VT and
//M left singular vectors to be computed in U
		il = 0;
//WORK(IL) is M by M
		ldwrkl = m;
		itau = il + ldwrkl * m;
		nwork = itau + m;
//Compute A=L*Q
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
		Rgelqf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
//Copy L to WORK(IL), zeroing out above it
		Rlacpy("L", m, m, &A[0], lda, &work[il], ldwrkl);
		Rlaset("U", m - 1, m - 1, Zero, Zero, &work[il + ldwrkl], ldwrkl);
//Generate Q in A
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
		Rorglq(m, n, m, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
		ie = itau;
		itauq = ie + m;
		itaup = itauq + m;
		nwork = itaup + m;
//Bidiagonalize L in WORK(IU), copying result to U
//(Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
		Rgebrd(m, m, &work[il], ldwrkl, &s[1], &work[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in U and computing right singular
//vectors of bidiagonal matrix in VT
//(Workspace: need M+BDSPAC)
		Rbdsdc("U", "I", m, &s[1], &work[ie], &u[0], ldu, &vt[0], ldvt, &dum, &idum, &work[nwork], &iwork[1], info);
//Overwrite U by left singular vectors of L and VT
//by right singular vectors of L
//(Workspace: need M*M+3*M, prefer M*M+2*M+M*NB)
		Rormbr("Q", "L", "N", m, m, m, &work[il], ldwrkl, &work[itauq], &u[0], ldu, &work[nwork], lwork - nwork + 1, &ierr);
		Rormbr("P", "R", "T", m, m, m, &work[il], ldwrkl, &work[itaup], &vt[0], ldvt, &work[nwork], lwork - nwork + 1, &ierr);
//Multiply right singular vectors of L in WORK(IL) by
//Q in A, storing result in VT
//(Workspace: need M*M)
		Rlacpy("F", m, m, &vt[0], ldvt, &work[il], ldwrkl);
		Rgemm("N", "N", m, n, m, One, &work[il], ldwrkl, &A[0], lda, Zero, &vt[0], ldvt);
	    } else if (wntqa) {
//Path 4t (N much larger than M, JOBZ='A')
//N right singular vectors to be computed in VT and
//M left singular vectors to be computed in U
		ivt = 1;
//WORK(IVT) is M by M
		ldwkvt = m;
		itau = ivt + ldwkvt * m;
		nwork = itau + m;
//Compute A=L*Q, copying result to VT
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
		Rgelqf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
		Rlacpy("U", m, n, &A[0], lda, &vt[0], ldvt);
//Generate Q in VT
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
		Rorglq(n, n, m, &vt[0], ldvt, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
//Produce L in A, zeroing out other entries
		Rlaset("U", m - 1, m - 1, Zero, Zero, &A[(lda << 1) + 1], lda);
		ie = itau;
		itauq = ie + m;
		itaup = itauq + m;
		nwork = itaup + m;
//Bidiagonalize L in A
//(Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
		Rgebrd(m, m, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in U and computing right singular
//vectors of bidiagonal matrix in WORK(IVT)
//(Workspace: need M+M*M+BDSPAC)
		Rbdsdc("U", "I", m, &s[1], &work[ie], &u[0], ldu, &work[ivt], ldwkvt, &dum, &idum, &work[nwork], &iwork[1]
		       , info);
//Overwrite U by left singular vectors of L and WORK(IVT)
//by right singular vectors of L
//(Workspace: need M*M+3*M, prefer M*M+2*M+M*NB)
		Rormbr("Q", "L", "N", m, m, m, &A[0], lda, &work[itauq], &u[0], ldu, &work[nwork], lwork - nwork + 1, &ierr);
		Rormbr("P", "R", "T", m, m, m, &A[0], lda, &work[itaup], &work[ivt], ldwkvt, &work[nwork], lwork - nwork + 1, &ierr);
//Multiply right singular vectors of L in WORK(IVT) by
//Q in VT, storing result in A
//(Workspace: need M*M)
		Rgemm("N", "N", m, n, m, One, &work[ivt], ldwkvt, &vt[0], ldvt, Zero, &A[0], lda);
//Copy right singular vectors of A from A to VT
		Rlacpy("F", m, n, &A[0], lda, &vt[0], ldvt);
	    }
	} else {
//N .LT. MNTHR
//Path 5t (N greater than M, but not much larger)
//Reduce to bidiagonal form without LQ decomposition
	    ie = 1;
	    itauq = ie + m;
	    itaup = itauq + m;
	    nwork = itaup + m;
//Bidiagonalize A
//(Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
	    Rgebrd(m, n, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
	    if (wntqn) {
//Perform bidiagonal SVD, only computing singular values
//(Workspace: need M+BDSPAC)
		Rbdsdc("L", "N", m, &s[1], &work[ie], &dum, 1, &dum, 1, &dum, &idum, &work[nwork], &iwork[1], info);
	    } else if (wntqo) {
		ldwkvt = m;
		ivt = nwork;
		if (lwork >= m * n + m * 3 + bdspac) {
//WORK( IVT ) is M by N
		    Rlaset("F", m, n, Zero, Zero, &work[ivt], ldwkvt);
		    nwork = ivt + ldwkvt * n;
		} else {
//WORK( IVT ) is M by M
		    nwork = ivt + ldwkvt * m;
		    il = nwork;
//WORK(IL) is M by CHUNK
		    chunk = (lwork - m * m - m * 3) / m;
		}
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in U and computing right singular
//vectors of bidiagonal matrix in WORK(IVT)
//(Workspace: need M*M+BDSPAC)
		Rbdsdc("L", "I", m, &s[1], &work[ie], &u[0], ldu, &work[ivt], ldwkvt, &dum, &idum, &work[nwork], &iwork[1]
		       , info);
//Overwrite U by left singular vectors of A
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
		Rormbr("Q", "L", "N", m, m, n, &A[0], lda, &work[itauq], &u[0], ldu, &work[nwork], lwork - nwork + 1, &ierr);
		if (lwork >= m * n + m * 3 + bdspac) {
//Overwrite WORK(IVT) by left singular vectors of A
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
		    Rormbr("P", "R", "T", m, n, m, &A[0], lda, &work[itaup], &work[ivt], ldwkvt, &work[nwork], lwork - nwork + 1, &ierr);
//Copy right singular vectors of A from WORK(IVT) to A
		    Rlacpy("F", m, n, &work[ivt], ldwkvt, &A[0], lda);
		} else {
//Generate P**T in A
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
		    Rorgbr("P", m, n, m, &A[0], lda, &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
//Multiply Q in A by right singular vectors of
//bidiagonal matrix in WORK(IVT), storing result in
//WORK(IL) and copying to A
//(Workspace: need 2*M*M, prefer M*M+M*N)
		    for (i = 1; i <= n; i += chunk) {
			blk = min(n - i + 1, chunk);
			Rgemm("N", "N", m, blk, m, One, &work[ivt], ldwkvt, &A[i * lda + 1], lda, Zero, &work[il], m);
			Rlacpy("F", m, blk, &work[il], m, &A[i * lda + 1], lda);
		    }
		}
	    } else if (wntqs) {
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in U and computing right singular
//vectors of bidiagonal matrix in VT
//(Workspace: need M+BDSPAC)
		Rlaset("F", m, n, Zero, Zero, &vt[0], ldvt);
		Rbdsdc("L", "I", m, &s[1], &work[ie], &u[0], ldu, &vt[0], ldvt, &dum, &idum, &work[nwork], &iwork[1], info);
//Overwrite U by left singular vectors of A and VT
//by right singular vectors of A
//(Workspace: need 3*M, prefer 2*M+M*NB)
		Rormbr("Q", "L", "N", m, m, n, &A[0], lda, &work[itauq], &u[0], ldu, &work[nwork], lwork - nwork + 1, &ierr);
		Rormbr("P", "R", "T", m, n, m, &A[0], lda, &work[itaup], &vt[0], ldvt, &work[nwork], lwork - nwork + 1, &ierr);
	    } else if (wntqa) {
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in U and computing right singular
//vectors of bidiagonal matrix in VT
//(Workspace: need M+BDSPAC)
		Rlaset("F", n, n, Zero, Zero, &vt[0], ldvt);
		Rbdsdc("L", "I", m, &s[1], &work[ie], &u[0], ldu, &vt[0], ldvt, &dum, &idum, &work[nwork], &iwork[1], info);
//Set the right corner of VT to identity matrix
		if (n > m) {
		    Rlaset("F", n - m, n - m, Zero, One, &vt[m + 1 + (m + 1) * ldvt], ldvt);
		}
//Overwrite U by left singular vectors of A and VT
//by right singular vectors of A
//(Workspace: need 2*M+N, prefer 2*M+N*NB)
		Rormbr("Q", "L", "N", m, m, n, &A[0], lda, &work[itauq], &u[0], ldu, &work[nwork], lwork - nwork + 1, &ierr);
		Rormbr("P", "R", "T", n, n, m, &A[0], lda, &work[itaup], &vt[0], ldvt, &work[nwork], lwork - nwork + 1, &ierr);
	    }
	}
    }
//Undo scaling if necessary
    if (iscl == 1) {
	if (anrm > bignum) {
	    Rlascl("G", 0, 0, bignum, anrm, minmn, 1, &s[1], minmn, &ierr);
	}
	if (anrm < smlnum) {
	    Rlascl("G", 0, 0, smlnum, anrm, minmn, 1, &s[1], minmn, &ierr);
	}
    }
//Return optimal workspace in WORK(1)
    work[1] = maxwrk;
    return;
}
