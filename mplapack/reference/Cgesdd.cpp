/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgesdd.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cgesdd(const char *jobz, INTEGER m, INTEGER n,
	    COMPLEX * A, INTEGER lda, REAL * s, COMPLEX * u, INTEGER ldu, COMPLEX * vt, INTEGER ldvt, COMPLEX * work, INTEGER lwork, REAL * rwork, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, ie, il, ir, iu, blk;
    REAL dum, eps;
    INTEGER iru, ivt, iscl;
    REAL anrm;
    INTEGER idum, ierr, itau, irvt;
    INTEGER chunk = 0, minmn;
    INTEGER wrkbl, itaup, itauq;
    LOGICAL wntqa;
    INTEGER nwork;
    LOGICAL wntqn, wntqo, wntqs;
    INTEGER mnthr1, mnthr2;
    REAL bignum;
    INTEGER ldwrkl;
    INTEGER ldwrkr, minwrk, ldwrku, maxwrk;
    INTEGER ldwkvt;
    REAL smlnum;
    LOGICAL wntqas;
    INTEGER nrwork;
    REAL Zero = 0.0, One = 1.0;

    *info = 0;
    minmn = min(m, n);
    mnthr1 = (INTEGER) (minmn * 17. / 9.);
    mnthr2 = (INTEGER) (minmn * 5. / 3.);
    wntqa = Mlsame(jobz, "A");
    wntqs = Mlsame(jobz, "S");
    wntqas = wntqa || wntqs;
    wntqo = Mlsame(jobz, "O");
    wntqn = Mlsame(jobz, "N");
    minwrk = 0;
    maxwrk = 0;
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
//  CWorkspace refers to complex workspace, and RWorkspace to
//  real workspace. NB refers to the optimal block size for the
//  immediately following subroutine, as returned by ILAENV.)
    if (*info == 0 && m > 0 && n > 0) {
	if (m >= n) {
//There is no complex work space needed for bidiagonal SVD
//The real work space needed for bidiagonal SVD is BDSPAC
//for computing singular values and singular vectors; BDSPAN
//for computing singular values only.
//BDSPAC = 5*N*N + 7*N
//BDSPAN = MAX(7*N+4, 3*N+2+SMLSIZ*(SMLSIZ+8))
	    if (m >= mnthr1) {
		if (wntqn) {
//Path 1 (M much larger than N, JOBZ='N')
		    maxwrk = n + n * iMlaenv(1, "Cgeqrf", " ", m, n, -1, -1);
		    maxwrk = max(maxwrk, (n << 1) + (n << 1) * iMlaenv(1, "Cgebrd", " ", n, n, -1, -1));
		    minwrk = n * 3;
		} else if (wntqo) {
//Path 2 (M much larger than N, JOBZ='O')
		    wrkbl = n + n * iMlaenv(1, "CGEQRF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, n + n * iMlaenv(1, "CUNGQR", " ", m, n, n, -1));
		    wrkbl = max(wrkbl, (n << 1) + (n << 1) * iMlaenv(1, "CGEBRD", " ", n, n, -1, -1));
		    wrkbl = max(wrkbl, (n << 1) + n * iMlaenv(1, "CUNMBR", "QLN", n, n, n, -1));
		    wrkbl = max(wrkbl, (n << 1) + n * iMlaenv(1, "CUNMBR", "PRC", n, n, n, -1));
		    maxwrk = m * n + n * n + wrkbl;
		    minwrk = (n << 1) * n + n * 3;
		} else if (wntqs) {
//Path 3 (M much larger than N, JOBZ='S')
		    wrkbl = n + n * iMlaenv(1, "CGEQRF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, n + n * iMlaenv(1, "CUNGQR", " ", m, n, n, -1));
		    wrkbl = max(wrkbl, (n << 1) + (n << 1) * iMlaenv(1, "CGEBRD", " ", n, n, -1, -1));
		    wrkbl = max(wrkbl, (n << 1) + n * iMlaenv(1, "CUNMBR", "QLN", n, n, n, -1));
		    wrkbl = max(wrkbl, (n << 1) + n * iMlaenv(1, "CUNMBR", "PRC", n, n, n, -1));
		    maxwrk = n * n + wrkbl;
		    minwrk = n * n + n * 3;
		} else if (wntqa) {
//Path 4 (M much larger than N, JOBZ='A')
		    wrkbl = n + n * iMlaenv(1, "CGEQRF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, n + m * iMlaenv(1, "CUNGQR", " ", m, m, n, -1));
		    wrkbl = max(wrkbl, (n << 1) + (n << 1) * iMlaenv(1, "CGEBRD", " ", n, n, -1, -1));
		    wrkbl = max(wrkbl, (n << 1) + n * iMlaenv(1, "CUNMBR", "QLN", n, n, n, -1));
		    wrkbl = max(wrkbl, (n << 1) + n * iMlaenv(1, "CUNMBR", "PRC", n, n, n, -1));
		    maxwrk = n * n + wrkbl;
		    minwrk = n * n + (n << 1) + m;
		}
	    } else if (m >= mnthr2) {
//Path 5 (M much larger than N, but not as much as MNTHR1)
		maxwrk = (n << 1) + (m + n) * iMlaenv(1, "CGEBRD", " ", m, n, -1, -1);
		minwrk = (n << 1) + m;
		if (wntqo) {
		    maxwrk = max(maxwrk, (n << 1) + n * iMlaenv(1, "CUNGBR", "P", n, n, n, -1));
		    maxwrk = max(maxwrk, (n << 1) + n * iMlaenv(1, "CUNGBR", "Q", m, n, n, -1));
		    maxwrk = maxwrk + m * n;
		    minwrk = minwrk + n * n;
		} else if (wntqs) {
		    maxwrk = max(maxwrk, (n << 1) + n * iMlaenv(1, "CUNGBR", "P", n, n, n, -1));
		    maxwrk = max(maxwrk, (n << 1) + n * iMlaenv(1, "CUNGBR", "Q", m, n, n, -1));
		} else if (wntqa) {
		    maxwrk = max(maxwrk, (n << 1) + n * iMlaenv(1, "CUNGBR", "P", n, n, n, -1));
		    maxwrk = max(maxwrk, (n << 1) + m * iMlaenv(1, "CUNGBR", "Q", m, m, n, -1));
		}
	    } else {
//Path 6 (M at least N, but not much larger)
		maxwrk = (n << 1) + (m + n) * iMlaenv(1, "CGEBRD", " ", m, n, -1, -1);
		minwrk = (n << 1) + m;
		if (wntqo) {
		    maxwrk = max(maxwrk, (n << 1) + n * iMlaenv(1, "CUNMBR", "PRC", n, n, n, -1));
		    maxwrk = max(maxwrk, (n << 1) + n * iMlaenv(1, "CUNMBR", "QLN", m, n, n, -1));
		    maxwrk = maxwrk + m * n;
		    minwrk = minwrk + n * n;
		} else if (wntqs) {
		    maxwrk = max(maxwrk, (n << 1) + n * iMlaenv(1, "CUNMBR", "PRC", n, n, n, -1));
		    maxwrk = max(maxwrk, (n << 1) + n * iMlaenv(1, "CUNMBR", "QLN", m, n, n, -1));
		} else if (wntqa) {
		    maxwrk = max(maxwrk, (n << 1) + n * iMlaenv(1, "CUNGBR", "PRC", n, n, n, -1));
		    maxwrk = max(maxwrk, (n << 1) + m * iMlaenv(1, "CUNGBR", "QLN", m, m, n, -1));
		}
	    }
	} else {
//There is no complex work space needed for bidiagonal SVD
//The real work space needed for bidiagonal SVD is BDSPAC
//for computing singular values and singular vectors; BDSPAN
//for computing singular values only.
//BDSPAC = 5*M*M + 7*M
//BDSPAN = MAX(7*M+4, 3*M+2+SMLSIZ*(SMLSIZ+8))
	    if (n >= mnthr1) {
		if (wntqn) {
//Path 1t (N much larger than M, JOBZ='N')
		    maxwrk = m + m * iMlaenv(1, "CGELQF", " ", m, n, -1, -1);
		    maxwrk = max(maxwrk, (m << 1) + (m << 1) * iMlaenv(1, "CGEBRD", " ", m, m, -1, -1));
		    minwrk = m * 3;
		} else if (wntqo) {
//Path 2t (N much larger than M, JOBZ='O')
		    wrkbl = m + m * iMlaenv(1, "CGELQF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, m + m * iMlaenv(1, "CUNGLQ", " ", m, n, m, -1));
		    wrkbl = max(wrkbl, (m << 1) + (m << 1) * iMlaenv(1, "CGEBRD", " ", m, m, -1, -1));
		    wrkbl = max(wrkbl, (m << 1) + m * iMlaenv(1, "CUNMBR", "PRC", m, m, m, -1));
		    wrkbl = max(wrkbl, (m << 1) + m * iMlaenv(1, "CUNMBR", "QLN", m, m, m, -1));
		    maxwrk = m * n + m * m + wrkbl;
		    minwrk = (m << 1) * m + m * 3;
		} else if (wntqs) {
//Path 3t (N much larger than M, JOBZ='S')
		    wrkbl = m + m * iMlaenv(1, "CGELQF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, m + m * iMlaenv(1, "CUNGLQ", " ", m, n, m, -1));
		    wrkbl = max(wrkbl, (m << 1) + (m << 1) * iMlaenv(1, "CGEBRD", " ", m, m, -1, -1));
		    wrkbl = max(wrkbl, (m << 1) + m * iMlaenv(1, "CUNMBR", "PRC", m, m, m, -1));
		    wrkbl = max(wrkbl, (m << 1) + m * iMlaenv(1, "CUNMBR", "QLN", m, m, m, -1));
		    maxwrk = m * m + wrkbl;
		    minwrk = m * m + m * 3;
		} else if (wntqa) {
//Path 4t (N much larger than M, JOBZ='A')
		    wrkbl = m + m * iMlaenv(1, "CGELQF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, m + n * iMlaenv(1, "CUNGLQ", " ", n, n, m, -1));
		    wrkbl = max(wrkbl, (m << 1) + (m << 1) * iMlaenv(1, "CGEBRD", " ", m, m, -1, -1));
		    wrkbl = max(wrkbl, (m << 1) + m * iMlaenv(1, "CUNMBR", "PRC", m, m, m, -1));
		    wrkbl = max(wrkbl, (m << 1) + m * iMlaenv(1, "CUNMBR", "QLN", m, m, m, -1));
		    maxwrk = m * m + wrkbl;
		    minwrk = m * m + (m << 1) + n;
		}
	    } else if (n >= mnthr2) {
//Path 5t (N much larger than M, but not as much as MNTHR1)
		maxwrk = (m << 1) + (m + n) * iMlaenv(1, "CGEBRD", " ", m, n, -1, -1);
		minwrk = (m << 1) + n;
		if (wntqo) {
		    maxwrk = max(maxwrk, (m << 1) + m * iMlaenv(1, "CUNGBR", "P", m, n, m, -1));
		    maxwrk = max(maxwrk, (m << 1) + m * iMlaenv(1, "CUNGBR", "Q", m, m, n, -1));
		    maxwrk += m * n;
		    minwrk += m * m;
		} else if (wntqs) {
		    maxwrk = max(maxwrk, (m << 1) + m * iMlaenv(1, "CUNGBR", "P", m, n, m, -1));
		    maxwrk = max(maxwrk, (m << 1) + m * iMlaenv(1, "CUNGBR", "Q", m, m, n, -1));
		} else if (wntqa) {
		    maxwrk = max(maxwrk, (m << 1) + n * iMlaenv(1, "CUNGBR", "P", n, n, m, -1));
		    maxwrk = max(maxwrk, (m << 1) + m * iMlaenv(1, "CUNGBR", "Q", m, m, n, -1));
		}
	    } else {
//Path 6t (N greater than M, but not much larger)
		maxwrk = (m << 1) + (m + n) * iMlaenv(1, "CGEBRD", " ", m, n, -1, -1);
		minwrk = (m << 1) + n;
		if (wntqo) {
		    maxwrk = max(maxwrk, (m << 1) + m * iMlaenv(1, "CUNMBR", "PRC", m, n, m, -1));
		    maxwrk = max(maxwrk, (m << 1) + m * iMlaenv(1, "CUNMBR", "QLN", m, m, n, -1));
		    maxwrk += m * n;
		    minwrk += m * m;
		} else if (wntqs) {
		    maxwrk = max(maxwrk, (m << 1) + m * iMlaenv(1, "CUNGBR", "PRC", m, n, m, -1));
		    maxwrk = max(maxwrk, (m << 1) + m * iMlaenv(1, "CUNGBR", "QLN", m, m, n, -1));
		} else if (wntqa) {
		    maxwrk = max(maxwrk, (m << 1) + n * iMlaenv(1, "CUNGBR", "PRC", n, n, m, -1));
		    maxwrk = max(maxwrk, (m << 1) + m * iMlaenv(1, "CUNGBR", "QLN", m, m, n, -1));
		}
	    }
	}
	maxwrk = max(maxwrk, minwrk);
    }
    if (*info == 0) {
	work[1] = maxwrk;
	if (lwork < minwrk && lwork != -1) {
	    *info = -13;
	}
    }
//Quick returns
    if (*info != 0) {
	Mxerbla("Cgesdd", -(*info));
	return;
    }
    if (lwork == -1) {
	return;
    }
    if (m == 0 || n == 0) {
	return;
    }
//Get machine constants
    eps = Rlamch("P");
    smlnum = sqrt(Rlamch("S")) / eps;
    bignum = One / smlnum;
//Scale A if max element outside range [SMLNUM,BIGNUM]
    anrm = Clange("M", m, n, &A[0], lda, &dum);
    iscl = 0;
    if (anrm > Zero && anrm < smlnum) {
	iscl = 0;
	Clascl("G", 0, 0, anrm, smlnum, m, n, &A[0], lda, &ierr);
    } else if (anrm > bignum) {
	iscl = 0;
	Clascl("G", 0, 0, anrm, bignum, m, n, &A[0], lda, &ierr);
    }
    if (m >= n) {
//A has at least as many rows as columns. If A has sufficiently
//more rows than columns, first reduce using the QR
//decomposition (if sufficient workspace available)
	if (m >= mnthr1) {
	    if (wntqn) {
//Path 1 (M much larger than N, JOBZ='N')
//No singular vectors to be computed
		itau = 1;
		nwork = itau + n;
//Compute A=Q*R
//(CWorkspace: need 2*N, prefer N+N*NB)
//(RWorkspace: need 0)
		Cgeqrf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
//Zero out below R
		Claset("L", n - 1, n - 1, Zero, Zero, &A[lda + 2], lda);
		ie = 1;
		itauq = 1;
		itaup = itauq + n;
		nwork = itaup + n;
//Bidiagonalize R in A
//(CWorkspace: need 3*N, prefer 2*N+2*N*NB)
//(RWorkspace: need N)
		Cgebrd(n, n, &A[0], lda, &s[1], &rwork[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
		nrwork = ie + n;
//Perform bidiagonal SVD, compute singular values only
//(CWorkspace: 0)
//(RWorkspace: need BDSPAN)
		Rbdsdc("U", "N", n, &s[1], &rwork[ie], &dum, 1, &dum, 1, &dum, &idum, &rwork[nrwork], &iwork[1], info);
	    } else if (wntqo) {
//Path 2 (M much larger than N, JOBZ='O')
//N left singular vectors to be overwritten on A and
//N right singular vectors to be computed in VT
		iu = 1;
//WORK(IU) is N by N
		ldwrku = n;
		ir = iu + ldwrku * n;
		if (lwork >= m * n + n * n + n * 3) {
//WORK(IR) is M by N
		    ldwrkr = m;
		} else {
		    ldwrkr = (lwork - n * n - n * 3) / n;
		}
		itau = ir + ldwrkr * n;
		nwork = itau + n;
//Compute A=Q*R
//(CWorkspace: need N*N+2*N, prefer M*N+N+N*NB)
//(RWorkspace: 0)
		Cgeqrf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
//Copy R to WORK( IR ), zeroing out below it
		Clacpy("U", n, n, &A[0], lda, &work[ir], ldwrkr);
		Claset("L", n - 1, n - 1, Zero, Zero, &work[ir + 1], ldwrkr);
//Generate Q in A
//(CWorkspace: need 2*N, prefer N+N*NB)
//(RWorkspace: 0)
		Cungqr(m, n, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
		ie = 1;
		itauq = itau;
		itaup = itauq + n;
		nwork = itaup + n;
//Bidiagonalize R in WORK(IR)
//(CWorkspace: need N*N+3*N, prefer M*N+2*N+2*N*NB)
//(RWorkspace: need N)
		Cgebrd(n, n, &work[ir], ldwrkr, &s[1], &rwork[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
//Perform bidiagonal SVD, computing left singular vectors
//of R in WORK(IRU) and computing right singular vectors
//of R in WORK(IRVT)
//(CWorkspace: need 0)
//(RWorkspace: need BDSPAC)
		iru = ie + n;
		irvt = iru + n * n;
		nrwork = irvt + n * n;
		Rbdsdc("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &rwork[irvt], n, &dum, &idum, &rwork[nrwork], &iwork[1], info);
//Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
//Overwrite WORK(IU) by the left singular vectors of R
//(CWorkspace: need 2*N*N+3*N, prefer M*N+N*N+2*N+N*NB)
//(RWorkspace: 0)
		Clacp2("F", n, n, &rwork[iru], n, &work[iu], ldwrku);
		Cunmbr("Q", "L", "N", n, n, n, &work[ir], ldwrkr, &work[itauq], &work[iu], ldwrku, &work[nwork], lwork - nwork + 1, &ierr);
//Copy real matrix RWORK(IRVT) to complex matrix VT
//Overwrite VT by the right singular vectors of R
//(CWorkspace: need N*N+3*N, prefer M*N+2*N+N*NB)
//(RWorkspace: 0)
		Clacp2("F", n, n, &rwork[irvt], n, &vt[0], ldvt);
		Cunmbr("P", "R", "C", n, n, n, &work[ir], ldwrkr, &work[itaup], &vt[0], ldvt, &work[nwork], lwork - nwork + 1, &ierr);
//Multiply Q in A by left singular vectors of R in
//WORK(IU), storing result in WORK(IR) and copying to A
//(CWorkspace: need 2*N*N, prefer N*N+M*N)
//k(RWorkspace: 0)
		for (i = 1; i <= m; i = i + ldwrkr) {
		    chunk = min(m - i + 1, ldwrkr);
		    Cgemm("N", "N", chunk, n, n, One, &A[i + lda], lda, &work[iu], ldwrku, Zero, &work[ir], ldwrkr);
		    Clacpy("F", chunk, n, &work[ir], ldwrkr, &A[i + lda], lda);
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
//(CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
//(RWorkspace: 0)
		Cgeqrf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
//Copy R to WORK(IR), zeroing out below it
		Clacpy("U", n, n, &A[0], lda, &work[ir], ldwrkr);
		Claset("L", n - 1, n - 1, Zero, Zero, &work[ir + 1], ldwrkr);
//Generate Q in A
//(CWorkspace: need 2*N, prefer N+N*NB)
//(RWorkspace: 0)
		Cungqr(m, n, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
		ie = 1;
		itauq = itau;
		itaup = itauq + n;
		nwork = itaup + n;
//Bidiagonalize R in WORK(IR)
//(CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
//(RWorkspace: need N)
		Cgebrd(n, n, &work[ir], ldwrkr, &s[1], &rwork[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in RWORK(IRU) and computing right
//singular vectors of bidiagonal matrix in RWORK(IRVT)
//(CWorkspace: need 0)
//(RWorkspace: need BDSPAC)
		iru = ie + n;
		irvt = iru + n * n;
		nrwork = irvt + n * n;
		Rbdsdc("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &rwork[irvt], n, &dum, &idum, &rwork[nrwork], &iwork[1], info);
//Copy real matrix RWORK(IRU) to complex matrix U
//Overwrite U by left singular vectors of R
//(CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
//(RWorkspace: 0)
		Clacp2("F", n, n, &rwork[iru], n, &u[0], ldu);
		Cunmbr("Q", "L", "N", n, n, n, &work[ir], ldwrkr, &work[itauq], &u[0], ldu, &work[nwork], lwork - nwork + 1, &ierr);
//Copy real matrix RWORK(IRVT) to complex matrix VT
//Overwrite VT by right singular vectors of R
//(CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
//(RWorkspace: 0)
		Clacp2("F", n, n, &rwork[irvt], n, &vt[0], ldvt);
		Cunmbr("P", "R", "C", n, n, n, &work[ir], ldwrkr, &work[itaup], &vt[0], ldvt, &work[nwork], lwork - nwork + 1, &ierr);
//Multiply Q in A by left singular vectors of R in
//WORK(IR), storing result in U
//(CWorkspace: need N*N)
//(RWorkspace: 0)
		Clacpy("F", n, n, &u[0], ldu, &work[ir], ldwrkr);
		Cgemm("N", "N", m, n, n, One, &A[0], lda, &work[ir], ldwrkr, Zero, &u[0], ldu);
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
//(CWorkspace: need 2*N, prefer N+N*NB)
//(RWorkspace: 0)
		Cgeqrf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
		Clacpy("L", m, n, &A[0], lda, &u[0], ldu);
//Generate Q in U
//(CWorkspace: need N+M, prefer N+M*NB)
//(RWorkspace: 0)
		Cungqr(m, m, n, &u[0], ldu, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
//Produce R in A, zeroing out below it
		Claset("L", n - 1, n - 1, Zero, Zero, &A[lda + 2], lda);
		ie = 1;
		itauq = itau;
		itaup = itauq + n;
		nwork = itaup + n;
//Bidiagonalize R in A
//(CWorkspace: need 3*N, prefer 2*N+2*N*NB)
//(RWorkspace: need N)
		Cgebrd(n, n, &A[0], lda, &s[1], &rwork[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
		iru = ie + n;
		irvt = iru + n * n;
		nrwork = irvt + n * n;
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in RWORK(IRU) and computing right
//singular vectors of bidiagonal matrix in RWORK(IRVT)
//(CWorkspace: need 0)
//(RWorkspace: need BDSPAC)
		Rbdsdc("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &rwork[irvt], n, &dum, &idum, &rwork[nrwork], &iwork[1], info);
//Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
//Overwrite WORK(IU) by left singular vectors of R
//(CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
//(RWorkspace: 0)
		Clacp2("F", n, n, &rwork[iru], n, &work[iu], ldwrku);
		Cunmbr("Q", "L", "N", n, n, n, &A[0], lda, &work[itauq], &work[iu], ldwrku, &work[nwork], lwork - nwork + 1, &ierr);
//Copy real matrix RWORK(IRVT) to complex matrix VT
//Overwrite VT by right singular vectors of R
//(CWorkspace: need 3*N, prefer 2*N+N*NB)
//(RWorkspace: 0)
		Clacp2("F", n, n, &rwork[irvt], n, &vt[0], ldvt);
		Cunmbr("P", "R", "C", n, n, n, &A[0], lda, &work[itaup], &vt[0], ldvt, &work[nwork], lwork - nwork + 1, &ierr);
//Multiply Q in U by left singular vectors of R in
//WORK(IU), storing result in A
//(CWorkspace: need N*N)
//(RWorkspace: 0)
		Cgemm("N", "N", m, n, n, One, &u[0], ldu, &work[iu], ldwrku, Zero, &A[0], lda);
//Copy left singular vectors of A from A to U
		Clacpy("F", m, n, &A[0], lda, &u[0], ldu);
	    }
	} else if (m >= mnthr2) {
//MNTHR2 <= M < MNTHR1
//Path 5 (M much larger than N, but not as much as MNTHR1)
//Reduce to bidiagonal form without QR decomposition, use
//ZUNGBR and matrix multiplication to compute singular vectors
	    ie = 1;
	    nrwork = ie + n;
	    itauq = 1;
	    itaup = itauq + n;
	    nwork = itaup + n;
//Bidiagonalize A
//(CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB)
//(RWorkspace: need N)
	    Cgebrd(m, n, &A[0], lda, &s[1], &rwork[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
	    if (wntqn) {
//Compute singular values only
//(Cworkspace: 0)
//(Rworkspace: need BDSPAN)
		Rbdsdc("U", "N", n, &s[1], &rwork[ie], &dum, 1, &dum, 1, &dum, &idum, &rwork[nrwork], &iwork[1], info);
	    } else if (wntqo) {
		iu = nwork;
		iru = nrwork;
		irvt = iru + n * n;
		nrwork = irvt + n * n;
//Copy A to VT, generate P**H
//(Cworkspace: need 2*N, prefer N+N*NB)
//(Rworkspace: 0)
		Clacpy("U", n, n, &A[0], lda, &vt[0], ldvt);
		Cungbr("P", n, n, n, &vt[0], ldvt, &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
//Generate Q in A
//(CWorkspace: need 2*N, prefer N+N*NB)
//(RWorkspace: 0)
		Cungbr("Q", m, n, n, &A[0], lda, &work[itauq], &work[nwork], lwork - nwork + 1, &ierr);
		if (lwork >= m * n + n * 3) {
//WORK( IU ) is M by N
		    ldwrku = m;
		} else {
//WORK(IU) is LDWRKU by N
		    ldwrku = (lwork - n * 3) / n;
		}
		nwork = iu + ldwrku * n;
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in RWORK(IRU) and computing right
//singular vectors of bidiagonal matrix in RWORK(IRVT)
//(CWorkspace: need 0)
//(RWorkspace: need BDSPAC)
		Rbdsdc("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &rwork[irvt], n, &dum, &idum, &rwork[nrwork], &iwork[1], info);
//Multiply real matrix RWORK(IRVT) by P**H in VT,
//storing the result in WORK(IU), copying to VT
//(Cworkspace: need 0)
//(Rworkspace: need 3*N*N)
		Clarcm(n, n, &rwork[irvt], n, &vt[0], ldvt, &work[iu]
		       , ldwrku, &rwork[nrwork]);
		Clacpy("F", n, n, &work[iu], ldwrku, &vt[0], ldvt);
//Multiply Q in A by real matrix RWORK(IRU), storing the
//result in WORK(IU), copying to A
//(CWorkspace: need N*N, prefer M*N)
//(Rworkspace: need 3*N*N, prefer N*N+2*M*N)
		nrwork = irvt;
		for (i = 1; i <= m; i += ldwrku) {
		    chunk = min(m - i + 1, ldwrku);
		    Clacrm(chunk, n, &A[i + lda], lda, &rwork[iru], n, &work[iu], ldwrku, &rwork[nrwork]);
		    Clacpy("F", chunk, n, &work[iu], ldwrku, &A[i + lda], lda);
		}
	    } else if (wntqs) {
//Copy A to VT, generate P**H
//(Cworkspace: need 2*N, prefer N+N*NB)
//(Rworkspace: 0)
		Clacpy("U", n, n, &A[0], lda, &vt[0], ldvt);
		Cungbr("P", n, n, n, &vt[0], ldvt, &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
//Copy A to U, generate Q
//(Cworkspace: need 2*N, prefer N+N*NB)
//(Rworkspace: 0)
		Clacpy("L", m, n, &A[0], lda, &u[0], ldu);
		Cungbr("Q", m, n, n, &u[0], ldu, &work[itauq], &work[nwork], lwork - nwork + 1, &ierr);
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in RWORK(IRU) and computing right
//singular vectors of bidiagonal matrix in RWORK(IRVT)
//(CWorkspace: need 0)
//(RWorkspace: need BDSPAC)
		iru = nrwork;
		irvt = iru + n * n;
		nrwork = irvt + n * n;
		Rbdsdc("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &rwork[irvt], n, &dum, &idum, &rwork[nrwork], &iwork[1], info);
//Multiply real matrix RWORK(IRVT) by P**H in VT,
//storing the result in A, copying to VT
//(Cworkspace: need 0)
//(Rworkspace: need 3*N*N)
		Clarcm(n, n, &rwork[irvt], n, &vt[0], ldvt, &A[0], lda, &rwork[nrwork]);
		Clacpy("F", n, n, &A[0], lda, &vt[0], ldvt);
//Multiply Q in U by real matrix RWORK(IRU), storing the
//result in A, copying to U
//(CWorkspace: need 0)
//(Rworkspace: need N*N+2*M*N)
		nrwork = irvt;
		Clacrm(m, n, &u[0], ldu, &rwork[iru], n, &A[0], lda, &rwork[nrwork]);
		Clacpy("F", m, n, &A[0], lda, &u[0], ldu);
	    } else {
//Copy A to VT, generate P**H
//(Cworkspace: need 2*N, prefer N+N*NB)
//(Rworkspace: 0)
		Clacpy("U", n, n, &A[0], lda, &vt[0], ldvt);
		Cungbr("P", n, n, n, &vt[0], ldvt, &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
//Copy A to U, generate Q
//(Cworkspace: need 2*N, prefer N+N*NB)
//(Rworkspace: 0)
		Clacpy("L", m, n, &A[0], lda, &u[0], ldu);
		Cungbr("Q", m, m, n, &u[0], ldu, &work[itauq], &work[nwork], lwork - nwork + 1, &ierr);
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in RWORK(IRU) and computing right
//singular vectors of bidiagonal matrix in RWORK(IRVT)
//(CWorkspace: need 0)
//(RWorkspace: need BDSPAC)
		iru = nrwork;
		irvt = iru + n * n;
		nrwork = irvt + n * n;
		Rbdsdc("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &rwork[irvt], n, &dum, &idum, &rwork[nrwork], &iwork[1], info);
//Multiply real matrix RWORK(IRVT) by P**H in VT,
//storing the result in A, copying to VT
//(Cworkspace: need 0)
//(Rworkspace: need 3*N*N)
		Clarcm(n, n, &rwork[irvt], n, &vt[0], ldvt, &A[0], lda, &rwork[nrwork]);
		Clacpy("F", n, n, &A[0], lda, &vt[0], ldvt);
//Multiply Q in U by real matrix RWORK(IRU), storing the
//result in A, copying to U
//(CWorkspace: 0)
//(Rworkspace: need 3*N*N)
		nrwork = irvt;
		Clacrm(m, n, &u[0], ldu, &rwork[iru], n, &A[0], lda, &rwork[nrwork]);
		Clacpy("F", m, n, &A[0], lda, &u[0], ldu);
	    }
	} else {
//M .LT. MNTHR2
//Path 6 (M at least N, but not much larger)
//Reduce to bidiagonal form without QR decomposition
//Use ZUNMBR to compute singular vectors
	    ie = 1;
	    nrwork = ie + n;
	    itauq = 1;
	    itaup = itauq + n;
	    nwork = itaup + n;
//Bidiagonalize A
//(CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB)
//(RWorkspace: need N)
	    Cgebrd(m, n, &A[0], lda, &s[1], &rwork[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
	    if (wntqn) {
//Compute singular values only
//(Cworkspace: 0)
//(Rworkspace: need BDSPAN)
		Rbdsdc("U", "N", n, &s[1], &rwork[ie], &dum, 1, &dum, 1, &dum, &idum, &rwork[nrwork], &iwork[1], info);
	    } else if (wntqo) {
		iu = nwork;
		iru = nrwork;
		irvt = iru + n * n;
		nrwork = irvt + n * n;
		if (lwork >= m * n + n * 3) {
//WORK( IU ) is M by N
		    ldwrku = m;
		} else {
//WORK( IU ) is LDWRKU by N
		    ldwrku = (lwork - n * 3) / n;
		}
		nwork = iu + ldwrku * n;
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in RWORK(IRU) and computing right
//singular vectors of bidiagonal matrix in RWORK(IRVT)
//(CWorkspace: need 0)
//(RWorkspace: need BDSPAC)
		Rbdsdc("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &rwork[irvt], n, &dum, &idum, &rwork[nrwork], &iwork[1], info);
//Copy real matrix RWORK(IRVT) to complex matrix VT
//Overwrite VT by right singular vectors of A
//(Cworkspace: need 2*N, prefer N+N*NB)
//(Rworkspace: need 0)
		Clacp2("F", n, n, &rwork[irvt], n, &vt[0], ldvt);
		Cunmbr("P", "R", "C", n, n, n, &A[0], lda, &work[itaup], &vt[0], ldvt, &work[nwork], lwork - nwork + 1, &ierr);
		if (lwork >= m * n + n * 3) {
//Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
//Overwrite WORK(IU) by left singular vectors of A, copying
//to A
//(Cworkspace: need M*N+2*N, prefer M*N+N+N*NB)
//(Rworkspace: need 0)
		    Claset("F", m, n, Zero, Zero, &work[iu], ldwrku);
		    Clacp2("F", n, n, &rwork[iru], n, &work[iu], ldwrku);
		    Cunmbr("Q", "L", "N", m, n, n, &A[0], lda, &work[itauq], &work[iu], ldwrku, &work[nwork], lwork - nwork + 1, &ierr);
		    Clacpy("F", m, n, &work[iu], ldwrku, &A[0], lda);
		} else {
//Generate Q in A
//(Cworkspace: need 2*N, prefer N+N*NB)
//(Rworkspace: need 0)
		    Cungbr("Q", m, n, n, &A[0], lda, &work[itauq], &work[nwork], lwork - nwork + 1, &ierr);
//Multiply Q in A by real matrix RWORK(IRU), storing the
//result in WORK(IU), copying to A
//(CWorkspace: need N*N, prefer M*N)
//(Rworkspace: need 3*N*N, prefer N*N+2*M*N)
		    nrwork = irvt;
		    for (i = 1; i <= m; i = i + ldwrku) {
			chunk = min(m - i + 1, ldwrku);
			Clacrm(chunk, n, &A[i + lda], lda, &rwork[iru], n, &work[iu], ldwrku, &rwork[nrwork]);
			Clacpy("F", chunk, n, &work[iu], ldwrku, &A[i + lda], lda);

		    }
		}
	    } else if (wntqs) {
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in RWORK(IRU) and computing right
//singular vectors of bidiagonal matrix in RWORK(IRVT)
//(CWorkspace: need 0)
//(RWorkspace: need BDSPAC)
		iru = nrwork;
		irvt = iru + n * n;
		nrwork = irvt + n * n;
		Rbdsdc("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &rwork[irvt], n, &dum, &idum, &rwork[nrwork], &iwork[1], info);
//Copy real matrix RWORK(IRU) to complex matrix U
//Overwrite U by left singular vectors of A
//(CWorkspace: need 3*N, prefer 2*N+N*NB)
//(RWorkspace: 0)
		Claset("F", m, n, Zero, Zero, &u[0], ldu);
		Clacp2("F", n, n, &rwork[iru], n, &u[0], ldu);
		Cunmbr("Q", "L", "N", m, n, n, &A[0], lda, &work[itauq], &u[0], ldu, &work[nwork], lwork - nwork + 1, &ierr);
//Copy real matrix RWORK(IRVT) to complex matrix VT
//Overwrite VT by right singular vectors of A
//(CWorkspace: need 3*N, prefer 2*N+N*NB)
//(RWorkspace: 0)
		Clacp2("F", n, n, &rwork[irvt], n, &vt[0], ldvt);
		Cunmbr("P", "R", "C", n, n, n, &A[0], lda, &work[itaup], &vt[0], ldvt, &work[nwork], lwork - nwork + 1, &ierr);
	    } else {
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in RWORK(IRU) and computing right
//singular vectors of bidiagonal matrix in RWORK(IRVT)
//(CWorkspace: need 0)
//(RWorkspace: need BDSPAC)
		iru = nrwork;
		irvt = iru + n * n;
		nrwork = irvt + n * n;
		Rbdsdc("U", "I", n, &s[1], &rwork[ie], &rwork[iru], n, &rwork[irvt], n, &dum, &idum, &rwork[nrwork], &iwork[1], info);
//Set the right corner of U to identity matrix
		Claset("F", m, m, Zero, Zero, &u[0], ldu);
		if (m > n) {
		    Claset("F", m - n, m - n, Zero, One, &u[n + 1 + (n + 1) * ldu], ldu);
		}
//Copy real matrix RWORK(IRU) to complex matrix U
//Overwrite U by left singular vectors of A
//(CWorkspace: need 2*N+M, prefer 2*N+M*NB)
//(RWorkspace: 0)
		Clacp2("F", n, n, &rwork[iru], n, &u[0], ldu);
		Cunmbr("Q", "L", "N", m, m, n, &A[0], lda, &work[itauq], &u[0], ldu, &work[nwork], lwork - nwork + 1, &ierr);
//Copy real matrix RWORK(IRVT) to complex matrix VT
//Overwrite VT by right singular vectors of A
//(CWorkspace: need 3*N, prefer 2*N+N*NB)
//(RWorkspace: 0)
		Clacp2("F", n, n, &rwork[irvt], n, &vt[0], ldvt);
		Cunmbr("P", "R", "C", n, n, n, &A[0], lda, &work[itaup], &vt[0], ldvt, &work[nwork], lwork - nwork + 1, &ierr);
	    }
	}
    } else {
//A has more columns than rows. If A has sufficiently more
//columns than rows, first reduce using the LQ decomposition (if
//sufficient workspace available)
	if (n >= mnthr1) {
	    if (wntqn) {
//Path 1t (N much larger than M, JOBZ='N')
//No singular vectors to be computed
		itau = 1;
		nwork = itau + m;
//Compute A=L*Q
//(CWorkspace: need 2*M, prefer M+M*NB)
//(RWorkspace: 0)
		Cgelqf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
//Zero out above L
		Claset("U", m - 1, m - 1, Zero, Zero, &A[(lda << 1) + 1], lda);
		ie = 1;
		itauq = 1;
		itaup = itauq + m;
		nwork = itaup + m;
//Bidiagonalize L in A
//(CWorkspace: need 3*M, prefer 2*M+2*M*NB)
//(RWorkspace: need M)
		Cgebrd(m, m, &A[0], lda, &s[1], &rwork[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
		nrwork = ie + m;
//Perform bidiagonal SVD, compute singular values only
//(CWorkspace: 0)
//(RWorkspace: need BDSPAN)
		Rbdsdc("U", "N", m, &s[1], &rwork[ie], &dum, 1, &dum, 1, &dum, &idum, &rwork[nrwork], &iwork[1], info);
	    } else if (wntqo) {
//Path 2t (N much larger than M, JOBZ='O')
//M right singular vectors to be overwritten on A and
//M left singular vectors to be computed in U
		ivt = 1;
		ldwkvt = m;
//WORK(IVT) is M by M
		il = ivt + ldwkvt * m;
		if (lwork >= m * n + m * m + m * 3) {
//WORK(IL) M by N
		    ldwrkl = m;
		    chunk = n;
		} else {
//WORK(IL) is M by CHUNK
		    ldwrkl = m;
		    chunk = (lwork - m * m - m * 3) / m;
		}
		itau = il + ldwrkl * chunk;
		nwork = itau + m;
//Compute A=L*Q
//(CWorkspace: need 2*M, prefer M+M*NB)
//(RWorkspace: 0)
		Cgelqf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
//Copy L to WORK(IL), zeroing about above it
		Clacpy("L", m, m, &A[0], lda, &work[il], ldwrkl);
		Claset("U", m - 1, m - 1, Zero, Zero, &work[il + ldwrkl], ldwrkl);
//Generate Q in A
//(CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
//(RWorkspace: 0)
		Cunglq(m, n, m, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
		ie = 1;
		itauq = itau;
		itaup = itauq + m;
		nwork = itaup + m;
//Bidiagonalize L in WORK(IL)
//(CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
//(RWorkspace: need M)
		Cgebrd(m, m, &work[il], ldwrkl, &s[1], &rwork[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in RWORK(IRU) and computing right
//singular vectors of bidiagonal matrix in RWORK(IRVT)
//(CWorkspace: need 0)
//(RWorkspace: need BDSPAC)
		iru = ie + m;
		irvt = iru + m * m;
		nrwork = irvt + m * m;
		Rbdsdc("U", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &rwork[irvt], m, &dum, &idum, &rwork[nrwork], &iwork[1], info);
//Copy real matrix RWORK(IRU) to complex matrix WORK(IU)
//Overwrite WORK(IU) by the left singular vectors of L
//(CWorkspace: need N*N+3*N, prefer M*N+2*N+N*NB)
//(RWorkspace: 0)
		Clacp2("F", m, m, &rwork[iru], m, &u[0], ldu);
		Cunmbr("Q", "L", "N", m, m, m, &work[il], ldwrkl, &work[itauq], &u[0], ldu, &work[nwork], lwork - nwork + 1, &ierr);
//Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT)
//Overwrite WORK(IVT) by the right singular vectors of L
//(CWorkspace: need N*N+3*N, prefer M*N+2*N+N*NB)
//(RWorkspace: 0)
		Clacp2("F", m, m, &rwork[irvt], m, &work[ivt], ldwkvt);
		Cunmbr("P", "R", "C", m, m, m, &work[il], ldwrkl, &work[itaup], &work[ivt], ldwkvt, &work[nwork], lwork - nwork + 1, &ierr);
//Multiply right singular vectors of L in WORK(IL) by Q
//in A, storing result in WORK(IL) and copying to A
//(CWorkspace: need 2*M*M, prefer M*M+M*N))
//(RWorkspace: 0)
		for (i = 1; i <= n; i += chunk) {
		    blk = min(n - i + 1, chunk);
		    Cgemm("N", "N", m, blk, m, One, &work[ivt], m, &A[i * lda + 1], lda, Zero, &work[il], ldwrkl);
		    Clacpy("F", m, blk, &work[il], ldwrkl, &A[i * lda + 1], lda);
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
//(CWorkspace: need 2*M, prefer M+M*NB)
//(RWorkspace: 0)
		Cgelqf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
//Copy L to WORK(IL), zeroing out above it
		Clacpy("L", m, m, &A[0], lda, &work[il], ldwrkl);
		Claset("U", m - 1, m - 1, Zero, Zero, &work[il + ldwrkl], ldwrkl);
//Generate Q in A
//(CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
//(RWorkspace: 0)
		Cunglq(m, n, m, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
		ie = 1;
		itauq = itau;
		itaup = itauq + m;
		nwork = itaup + m;
//Bidiagonalize L in WORK(IL)
//(CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
//(RWorkspace: need M)
		Cgebrd(m, m, &work[il], ldwrkl, &s[1], &rwork[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in RWORK(IRU) and computing right
//singular vectors of bidiagonal matrix in RWORK(IRVT)
//(CWorkspace: need 0)
//(RWorkspace: need BDSPAC)
		iru = ie + m;
		irvt = iru + m * m;
		nrwork = irvt + m * m;
		Rbdsdc("U", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &rwork[irvt], m, &dum, &idum, &rwork[nrwork], &iwork[1], info);
//Copy real matrix RWORK(IRU) to complex matrix U
//Overwrite U by left singular vectors of L
//(CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB)
//(RWorkspace: 0)
		Clacp2("F", m, m, &rwork[iru], m, &u[0], ldu);
		Cunmbr("Q", "L", "N", m, m, m, &work[il], ldwrkl, &work[itauq], &u[0], ldu, &work[nwork], lwork - nwork + 1, &ierr);
//Copy real matrix RWORK(IRVT) to complex matrix VT
//Overwrite VT by left singular vectors of L
//(CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB)
//(RWorkspace: 0)
		Clacp2("F", m, m, &rwork[irvt], m, &vt[0], ldvt);
		Cunmbr("P", "R", "C", m, m, m, &work[il], ldwrkl, &work[itaup], &vt[0], ldvt, &work[nwork], lwork - nwork + 1, &ierr);
//Copy VT to WORK(IL), multiply right singular vectors of L
//in WORK(IL) by Q in A, storing result in VT
//(CWorkspace: need M*M)
//(RWorkspace: 0)
		Clacpy("F", m, m, &vt[0], ldvt, &work[il], ldwrkl);
		Cgemm("N", "N", m, n, m, One, &work[il], ldwrkl, &A[0], lda, Zero, &vt[0], ldvt);
	    } else if (wntqa) {
//Path 9t (N much larger than M, JOBZ='A')
//N right singular vectors to be computed in VT and
//M left singular vectors to be computed in U
		ivt = 1;
//WORK(IVT) is M by M
		ldwkvt = m;
		itau = ivt + ldwkvt * m;
		nwork = itau + m;
//Compute A=L*Q, copying result to VT
//(CWorkspace: need 2*M, prefer M+M*NB)
//(RWorkspace: 0)
		Cgelqf(m, n, &A[0], lda, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
		Clacpy("U", m, n, &A[0], lda, &vt[0], ldvt);
//Generate Q in VT
//(CWorkspace: need M+N, prefer M+N*NB)
//(RWorkspace: 0)
		Cunglq(n, n, m, &vt[0], ldvt, &work[itau], &work[nwork], lwork - nwork + 1, &ierr);
//Produce L in A, zeroing out above it
		Claset("U", m - 1, m - 1, Zero, Zero, &A[(lda << 1) + 1], lda);
		ie = 1;
		itauq = itau;
		itaup = itauq + m;
		nwork = itaup + m;
//Bidiagonalize L in A
//(CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
//(RWorkspace: need M)
		Cgebrd(m, m, &A[0], lda, &s[1], &rwork[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in RWORK(IRU) and computing right
//singular vectors of bidiagonal matrix in RWORK(IRVT)
//(CWorkspace: need 0)
//(RWorkspace: need BDSPAC)
		iru = ie + m;
		irvt = iru + m * m;
		nrwork = irvt + m * m;
		Rbdsdc("U", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &rwork[irvt], m, &dum, &idum, &rwork[nrwork], &iwork[1], info);
//Copy real matrix RWORK(IRU) to complex matrix U
//Overwrite U by left singular vectors of L
//(CWorkspace: need 3*M, prefer 2*M+M*NB)
//(RWorkspace: 0)
		Clacp2("F", m, m, &rwork[iru], m, &u[0], ldu);
		Cunmbr("Q", "L", "N", m, m, m, &A[0], lda, &work[itauq], &u[0], ldu, &work[nwork], lwork - nwork + 1, &ierr);
//Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT)
//Overwrite WORK(IVT) by right singular vectors of L
//(CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB)
//(RWorkspace: 0)
		Clacp2("F", m, m, &rwork[irvt], m, &work[ivt], ldwkvt);
		Cunmbr("P", "R", "C", m, m, m, &A[0], lda, &work[itaup], &work[ivt], ldwkvt, &work[nwork], lwork - nwork + 1, &ierr);
//Multiply right singular vectors of L in WORK(IVT) by
//Q in VT, storing result in A
//(CWorkspace: need M*M)
//(RWorkspace: 0)
		Cgemm("N", "N", m, n, m, One, &work[ivt], ldwkvt, &vt[0], ldvt, Zero, &A[0], lda);
//Copy right singular vectors of A from A to VT
		Clacpy("F", m, n, &A[0], lda, &vt[0], ldvt);
	    }
	} else if (n >= mnthr2) {
//MNTHR2 <= N < MNTHR1
//Path 5t (N much larger than M, but not as much as MNTHR1)
//Reduce to bidiagonal form without QR decomposition, use
//ZUNGBR and matrix multiplication to compute singular vectors
	    ie = 1;
	    nrwork = ie + m;
	    itauq = 1;
	    itaup = itauq + m;
	    nwork = itaup + m;
//Bidiagonalize A
//(CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)
//(RWorkspace: M)
	    Cgebrd(m, n, &A[0], lda, &s[1], &rwork[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
	    if (wntqn) {
//Compute singular values only
//(Cworkspace: 0)
//(Rworkspace: need BDSPAN)
		Rbdsdc("L", "N", m, &s[1], &rwork[ie], &dum, 1, &dum, 1, &dum, &idum, &rwork[nrwork], &iwork[1], info);
	    } else if (wntqo) {
		irvt = nrwork;
		iru = irvt + m * m;
		nrwork = iru + m * m;
		ivt = nwork;
//Copy A to U, generate Q
//(Cworkspace: need 2*M, prefer M+M*NB)
//(Rworkspace: 0)
		Clacpy("L", m, m, &A[0], lda, &u[0], ldu);
		Cungbr("Q", m, m, n, &u[0], ldu, &work[itauq], &work[nwork], lwork - nwork + 1, &ierr);
//Generate P**H in A
//(Cworkspace: need 2*M, prefer M+M*NB)
//(Rworkspace: 0)
		Cungbr("P", m, n, m, &A[0], lda, &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
		ldwkvt = m;
		if (lwork >= m * n + m * 3) {
//WORK( IVT ) is M by N
		    nwork = ivt + ldwkvt * n;
		    chunk = n;
		} else {
//WORK( IVT ) is M by CHUNK
		    chunk = (lwork - m * 3) / m;
		    nwork = ivt + ldwkvt * chunk;
		}
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in RWORK(IRU) and computing right
//singular vectors of bidiagonal matrix in RWORK(IRVT)
//(CWorkspace: need 0)
//(RWorkspace: need BDSPAC)
		Rbdsdc("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &rwork[irvt], m, &dum, &idum, &rwork[nrwork], &iwork[1], info);
//Multiply Q in U by real matrix RWORK(IRVT)
//storing the result in WORK(IVT), copying to U
//(Cworkspace: need 0)
//(Rworkspace: need 2*M*M)
		Clacrm(m, m, &u[0], ldu, &rwork[iru], m, &work[ivt], ldwkvt, &rwork[nrwork]);
		Clacpy("F", m, m, &work[ivt], ldwkvt, &u[0], ldu);
//Multiply RWORK(IRVT) by P**H in A, storing the
//result in WORK(IVT), copying to A
//(CWorkspace: need M*M, prefer M*N)
//(Rworkspace: need 2*M*M, prefer 2*M*N)
		nrwork = iru;
		for (i = 1; i <= n; i += chunk) {
		    blk = min(n - i + 1, chunk);
		    Clarcm(m, blk, &rwork[irvt], m, &A[i * lda + 1], lda, &work[ivt], ldwkvt, &rwork[nrwork]);
		    Clacpy("F", m, blk, &work[ivt], ldwkvt, &A[i * lda + 1], lda);
		}
	    } else if (wntqs) {
//Copy A to U, generate Q
//(Cworkspace: need 2*M, prefer M+M*NB)
//(Rworkspace: 0)
		Clacpy("L", m, m, &A[0], lda, &u[0], ldu);
		Cungbr("Q", m, m, n, &u[0], ldu, &work[itauq], &work[nwork], lwork - nwork + 1, &ierr);
//Copy A to VT, generate P**H
//(Cworkspace: need 2*M, prefer M+M*NB)
//(Rworkspace: 0)
		Clacpy("U", m, n, &A[0], lda, &vt[0], ldvt);
		Cungbr("P", m, n, m, &vt[0], ldvt, &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in RWORK(IRU) and computing right
//singular vectors of bidiagonal matrix in RWORK(IRVT)
//(CWorkspace: need 0)
//(RWorkspace: need BDSPAC)
		irvt = nrwork;
		iru = irvt + m * m;
		nrwork = iru + m * m;
		Rbdsdc("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &rwork[irvt], m, &dum, &idum, &rwork[nrwork], &iwork[1], info);
//Multiply Q in U by real matrix RWORK(IRU), storing the
//result in A, copying to U
//(CWorkspace: need 0)
//(Rworkspace: need 3*M*M)
		Clacrm(m, m, &u[0], ldu, &rwork[iru], m, &A[0], lda, &rwork[nrwork]);
		Clacpy("F", m, m, &A[0], lda, &u[0], ldu);
//Multiply real matrix RWORK(IRVT) by P**H in VT,
//storing the result in A, copying to VT
//(Cworkspace: need 0)
//(Rworkspace: need M*M+2*M*N)
		nrwork = iru;
		Clarcm(m, n, &rwork[irvt], m, &vt[0], ldvt, &A[0], lda, &rwork[nrwork]);
		Clacpy("F", m, n, &A[0], lda, &vt[0], ldvt);
	    } else {
//Copy A to U, generate Q
//(Cworkspace: need 2*M, prefer M+M*NB)
//(Rworkspace: 0)
		Clacpy("L", m, m, &A[0], lda, &u[0], ldu);
		Cungbr("Q", m, m, n, &u[0], ldu, &work[itauq], &work[nwork], lwork - nwork + 1, &ierr);
//Copy A to VT, generate P**H
//(Cworkspace: need 2*M, prefer M+M*NB)
//(Rworkspace: 0)
		Clacpy("U", m, n, &A[0], lda, &vt[0], ldvt);
		Cungbr("P", n, n, m, &vt[0], ldvt, &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in RWORK(IRU) and computing right
//singular vectors of bidiagonal matrix in RWORK(IRVT)
//(CWorkspace: need 0)
//(RWorkspace: need BDSPAC)
		irvt = nrwork;
		iru = irvt + m * m;
		nrwork = iru + m * m;
		Rbdsdc("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &rwork[irvt], m, &dum, &idum, &rwork[nrwork], &iwork[1], info);
//Multiply Q in U by real matrix RWORK(IRU), storing the
//result in A, copying to U
//(CWorkspace: need 0)
//(Rworkspace: need 3*M*M)
		Clacrm(m, m, &u[0], ldu, &rwork[iru], m, &A[0], lda, &rwork[nrwork]);
		Clacpy("F", m, m, &A[0], lda, &u[0], ldu);
//Multiply real matrix RWORK(IRVT) by P**H in VT,
//storing the result in A, copying to VT
//(Cworkspace: need 0)
//(Rworkspace: need M*M+2*M*N)
		Clarcm(m, n, &rwork[irvt], m, &vt[0], ldvt, &A[0], lda, &rwork[nrwork]);
		Clacpy("F", m, n, &A[0], lda, &vt[0], ldvt);
	    }
	} else {
//N .LT. MNTHR2
//Path 6t (N greater than M, but not much larger)
//Reduce to bidiagonal form without LQ decomposition
//Use ZUNMBR to compute singular vectors
	    ie = 1;
	    nrwork = ie + m;
	    itauq = 1;
	    itaup = itauq + m;
	    nwork = itaup + m;
//Bidiagonalize A
//(CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)
//(RWorkspace: M)
	    Cgebrd(m, n, &A[0], lda, &s[1], &rwork[ie], &work[itauq], &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
	    if (wntqn) {
//Compute singular values only
//(Cworkspace: 0)
//(Rworkspace: need BDSPAN)
		Rbdsdc("L", "N", m, &s[1], &rwork[ie], &dum, 1, &dum, 1, &dum, &idum, &rwork[nrwork], &iwork[1], info);
	    } else if (wntqo) {
		ldwkvt = m;
		ivt = nwork;
		if (lwork >= m * n + m * 3) {
//WORK( IVT ) is M by N
		    Claset("F", m, n, Zero, Zero, &work[ivt], ldwkvt);
		    nwork = ivt + ldwkvt * n;
		} else {
//WORK( IVT ) is M by CHUNK
		    chunk = (lwork - m * 3) / m;
		    nwork = ivt + ldwkvt * chunk;
		}
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in RWORK(IRU) and computing right
//singular vectors of bidiagonal matrix in RWORK(IRVT)
//(CWorkspace: need 0)
//(RWorkspace: need BDSPAC)
		irvt = nrwork;
		iru = irvt + m * m;
		nrwork = iru + m * m;
		Rbdsdc("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &rwork[irvt], m, &dum, &idum, &rwork[nrwork], &iwork[1], info);
//Copy real matrix RWORK(IRU) to complex matrix U
//Overwrite U by left singular vectors of A
//(Cworkspace: need 2*M, prefer M+M*NB)
//(Rworkspace: need 0)
		Clacp2("F", m, m, &rwork[iru], m, &u[0], ldu);
		Cunmbr("Q", "L", "N", m, m, n, &A[0], lda, &work[itauq], &u[0], ldu, &work[nwork], lwork - nwork + 1, &ierr);
		if (lwork >= m * n + m * 3) {
//Copy real matrix RWORK(IRVT) to complex matrix WORK(IVT)
//Overwrite WORK(IVT) by right singular vectors of A,
//copying to A
//(Cworkspace: need M*N+2*M, prefer M*N+M+M*NB)
//(Rworkspace: need 0)
		    Clacp2("F", m, m, &rwork[irvt], m, &work[ivt], ldwkvt);
		    Cunmbr("P", "R", "C", m, n, m, &A[0], lda, &work[itaup], &work[ivt], ldwkvt, &work[nwork], lwork - nwork + 1, &ierr);
		    Clacpy("F", m, n, &work[ivt], ldwkvt, &A[0], lda);
		} else {
//Generate P**H in A
//(Cworkspace: need 2*M, prefer M+M*NB)
//(Rworkspace: need 0)
		    Cungbr("P", m, n, m, &A[0], lda, &work[itaup], &work[nwork], lwork - nwork + 1, &ierr);
//Multiply Q in A by real matrix RWORK(IRU), storing the
//result in WORK(IU), copying to A
//(CWorkspace: need M*M, prefer M*N)
//(Rworkspace: need 3*M*M, prefer M*M+2*M*N)
		    nrwork = iru;
		    for (i = 1; i <= n; i += chunk) {
			blk = min(n - i + 1, chunk);
			Clarcm(m, blk, &rwork[irvt], m, &A[i * lda + 1], lda, &work[ivt], ldwkvt, &rwork[nrwork]);
			Clacpy("F", m, blk, &work[ivt], ldwkvt, &A[i * lda + 1], lda);
		    }
		}
	    } else if (wntqs) {
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in RWORK(IRU) and computing right
//singular vectors of bidiagonal matrix in RWORK(IRVT)
//(CWorkspace: need 0)
//(RWorkspace: need BDSPAC)
		irvt = nrwork;
		iru = irvt + m * m;
		nrwork = iru + m * m;
		Rbdsdc("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &rwork[irvt], m, &dum, &idum, &rwork[nrwork], &iwork[1], info);
//Copy real matrix RWORK(IRU) to complex matrix U
//Overwrite U by left singular vectors of A
//(CWorkspace: need 3*M, prefer 2*M+M*NB)
//(RWorkspace: M*M)
		Clacp2("F", m, m, &rwork[iru], m, &u[0], ldu);
		Cunmbr("Q", "L", "N", m, m, n, &A[0], lda, &work[itauq], &u[0], ldu, &work[nwork], lwork - nwork + 1, &ierr);
//Copy real matrix RWORK(IRVT) to complex matrix VT
//Overwrite VT by right singular vectors of A
//(CWorkspace: need 3*M, prefer 2*M+M*NB)
//(RWorkspace: M*M)
		Claset("F", m, n, Zero, Zero, &vt[0], ldvt);
		Clacp2("F", m, m, &rwork[irvt], m, &vt[0], ldvt);
		Cunmbr("P", "R", "C", m, n, m, &A[0], lda, &work[itaup], &vt[0], ldvt, &work[nwork], lwork - nwork + 1, &ierr);
	    } else {
//Perform bidiagonal SVD, computing left singular vectors
//of bidiagonal matrix in RWORK(IRU) and computing right
//singular vectors of bidiagonal matrix in RWORK(IRVT)
//(CWorkspace: need 0)
//(RWorkspace: need BDSPAC)
		irvt = nrwork;
		iru = irvt + m * m;
		nrwork = iru + m * m;
		Rbdsdc("L", "I", m, &s[1], &rwork[ie], &rwork[iru], m, &rwork[irvt], m, &dum, &idum, &rwork[nrwork], &iwork[1], info);
//Copy real matrix RWORK(IRU) to complex matrix U
//Overwrite U by left singular vectors of A
//(CWorkspace: need 3*M, prefer 2*M+M*NB)
//(RWorkspace: M*M)
		Clacp2("F", m, m, &rwork[iru], m, &u[0], ldu);
		Cunmbr("Q", "L", "N", m, m, n, &A[0], lda, &work[itauq], &u[0], ldu, &work[nwork], lwork - nwork + 1, &ierr);
//Set all of VT to identity matrix
		Claset("F", n, n, Zero, One, &vt[0], ldvt);
//Copy real matrix RWORK(IRVT) to complex matrix VT
//Overwrite VT by right singular vectors of A
//(CWorkspace: need 2*M+N, prefer 2*M+N*NB)
//(RWorkspace: M*M)
		Clacp2("F", m, m, &rwork[irvt], m, &vt[0], ldvt);
		Cunmbr("P", "R", "C", n, n, m, &A[0], lda, &work[itaup], &vt[0], ldvt, &work[nwork], lwork - nwork + 1, &ierr);
	    }
	}
    }
//Undo scaling if necessary
    if (iscl == 1) {
	if (anrm > bignum) {
	    Rlascl("G", 0, 0, bignum, anrm, minmn, 1, &s[1], minmn, &ierr);
	}
	if (*info != 0 && anrm > bignum) {
	    Rlascl("G", 0, 0, bignum, anrm, minmn - 1, 1, &rwork[ie], minmn, &ierr);
	}
	if (anrm < smlnum) {
	    Rlascl("G", 0, 0, smlnum, anrm, minmn, 1, &s[1], minmn, &ierr);
	}
	if (*info != 0 && anrm < smlnum) {
	    Rlascl("G", 0, 0, smlnum, anrm, minmn - 1, 1, &rwork[ie], minmn, &ierr);
	}
    }
//Return optimal workspace in WORK(1)
    work[1] = maxwrk;
    return;
}
