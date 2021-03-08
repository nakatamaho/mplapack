/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cggev.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

#define MTRUE 1
#define MFALSE 0

void Cggev(const char *jobvl, const char *jobvr, INTEGER n,
	   COMPLEX * A, INTEGER lda, COMPLEX * B, INTEGER ldb,
	   COMPLEX * alpha, COMPLEX * beta, COMPLEX * vl, INTEGER ldvl, COMPLEX * vr, INTEGER ldvr, COMPLEX * work, INTEGER lwork, REAL * rwork, INTEGER * info)
{
    INTEGER jc, in, jr, ihi, ilo;
    REAL eps;
    LOGICAL ilv;
    REAL anrm, bnrm;
    INTEGER ierr, itau;
    REAL temp;
    LOGICAL ilvl, ilvr;
    INTEGER iwrk;
    INTEGER ileft, icols, irwrk, irows;
    LOGICAL ilascl, ilbscl;
    LOGICAL ldummy;
    char chtemp;
    REAL bignum;
    INTEGER ijobvl, iright;
    INTEGER ijobvr;
    REAL anrmto = 0.0;
    INTEGER lwkmin;
    REAL bnrmto = 0.0;
    REAL smlnum;
    INTEGER lwkopt;
    LOGICAL lquery;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1, mtemp2;
    if (Mlsame(jobvl, "N")) {
	ijobvl = 0;
	ilvl = MFALSE;
    } else if (Mlsame(jobvl, "V")) {
	ijobvl = 2;
	ilvl = MTRUE;
    } else {
	ijobvl = -1;
	ilvl = MFALSE;
    }
    if (Mlsame(jobvr, "N")) {
	ijobvr = 1;
	ilvr = MFALSE;
    } else if (Mlsame(jobvr, "V")) {
	ijobvr = 2;
	ilvr = MTRUE;
    } else {
	ijobvr = -1;
	ilvr = MFALSE;
    }
    ilv = ilvl || ilvr;
//Test the input arguments
    *info = 0;
    lquery = lwork == -1;
    if (ijobvl <= 0) {
	*info = -1;
    } else if (ijobvr <= 0) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -5;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -7;
    } else if (ldvl < 1 || (ilvl && ldvl < n)) {
	*info = -11;
    } else if (ldvr < 1 || (ilvr && ldvr < n)) {
	*info = -13;
    }
//Compute workspace
// (Note: Comments in the code beginning "Workspace:" describe the
//  minimal amount of workspace needed at that point in the code,
//  as well as the preferred amount for good performance.
//  NB refers to the optimal block size for the immediately
//  following subroutine, as returned by ILAENV. The workspace is
//  computed assuming ILO = 1 and IHI = N, the worst case.)
    if (*info == 0) {
	lwkmin = max((INTEGER) 1, n << 1);
	lwkopt = max((INTEGER) 1, n + n * iMlaenv(1, "Cgeqrf", " ", n, 1, n, 0));
	lwkopt = max(lwkopt, n + n * iMlaenv(1, "Cunmqr", " ", n, 1, n, 0));
	if (ilvl) {
	    lwkopt = max(lwkopt, n + n * iMlaenv(1, "CZungqr", " ", n, 1, n, -1));
	}
	work[1] = lwkopt;
	if (lwork < lwkmin && !lquery) {
	    *info = -15;
	}
    }
    if (*info != 0) {
	Mxerbla("Cggev ", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//Get machine constants
    eps = Rlamch("E") * Rlamch("B");
    smlnum = Rlamch("S");
    bignum = One / smlnum;
    //Rlabad(&smlnum, &bignum);
    smlnum = sqrt(smlnum) / eps;
    bignum = One / smlnum;
//Scale A if max element outside range [SMLNUM,BIGNUM]
    anrm = Clange("M", n, n, &A[0], lda, &rwork[1]);
    ilascl = MFALSE;
    if (anrm > Zero && anrm < smlnum) {
	anrmto = smlnum;
	ilascl = MTRUE;
    } else if (anrm > bignum) {
	anrmto = bignum;
	ilascl = MTRUE;
    }
    if (ilascl) {
	Clascl("G", 0, 0, anrm, anrmto, n, n, &A[0], lda, &ierr);
    }
//Scale B if max element outside range [SMLNUM,BIGNUM]
    bnrm = Clange("M", n, n, &B[0], ldb, &rwork[1]);
    ilbscl = MFALSE;
    if (bnrm > Zero && bnrm < smlnum) {
	bnrmto = smlnum;
	ilbscl = MTRUE;
    } else if (bnrm > bignum) {
	bnrmto = bignum;
	ilbscl = MTRUE;
    }
    if (ilbscl) {
	Clascl("G", 0, 0, bnrm, bnrmto, n, n, &B[0], ldb, &ierr);
    }
//Permute the matrices A, B to isolate eigenvalues if possible
//(Real Workspace: need 6*N)
    ileft = 1;
    iright = n + 1;
    irwrk = iright + n;
    Cggbal("P", n, &A[0], lda, &B[0], ldb, &ilo, &ihi, &rwork[ileft], &rwork[iright], &rwork[irwrk], &ierr);
//Reduce B to triangular form (QR decomposition of B)
//(Complex Workspace: need N, prefer N*NB)
    irows = ihi + 1 - ilo;
    if (ilv) {
	icols = n + 1 - ilo;
    } else {
	icols = irows;
    }
    itau = 1;
    iwrk = itau + irows;
    Cgeqrf(irows, icols, &B[ilo + ilo * ldb], ldb, &work[itau], &work[iwrk], lwork + 1 - iwrk, &ierr);
//Apply the orthogonal transformation to matrix A
//(Complex Workspace: need N, prefer N*NB)
    Cunmqr("L", "C", irows, icols, irows, &B[ilo + ilo * ldb], ldb, &work[itau], &A[ilo + ilo * lda], lda, &work[iwrk], lwork + 1 - iwrk, &ierr);
//Initialize VL
//(Complex Workspace: need N, prefer N*NB)
    if (ilvl) {
	Claset("Full", n, n, Zero, One, &vl[0], ldvl);
	if (irows > 1) {
	    Clacpy("L", irows - 1, irows - 1, &B[ilo + 1 + ilo * ldb], ldb, &vl[ilo + 1 + ilo * ldvl], ldvl);
	}
	Cungqr(irows, irows, irows, &vl[ilo + ilo * ldvl], ldvl, &work[itau], &work[iwrk], lwork + 1 - iwrk, &ierr);
    }
//Initialize VR
    if (ilvr) {
	Claset("Full", n, n, Zero, One, &vr[0], ldvr);
    }
//Reduce to generalized Hessenberg form
    if (ilv) {
//Eigenvectors requested -- work on whole matrix.
	Cgghrd(jobvl, jobvr, n, ilo, ihi, &A[0], lda, &B[0], ldb, &vl[0], ldvl, &vr[0], ldvr, &ierr);
    } else {
	Cgghrd("N", "N", irows, 1, irows, &A[ilo + ilo * lda], lda, &B[ilo + ilo * ldb], ldb, &vl[0], ldvl, &vr[0], ldvr, &ierr);
    }
//Perform QZ algorithm (Compute eigenvalues, and optionally, the
//Schur form and Schur vectors)
//(Complex Workspace: need N)
//(Real Workspace: need N)
    iwrk = itau;
    if (ilv) {
	chtemp = 'S';
    } else {
	chtemp = 'E';
    }
    Chgeqz((const char *) chtemp, jobvl, jobvr, n, ilo, ihi, &A[0], lda, &B[0], ldb, &alpha[1], &beta[1], &vl[0], ldvl, &vr[0],
	   ldvr, &work[iwrk], lwork + 1 - iwrk, &rwork[irwrk], &ierr);
    if (ierr != 0) {
	if (ierr > 0 && ierr <= n) {
	    *info = ierr;
	} else if (ierr > n && ierr <= n << 1) {
	    *info = ierr - n;
	} else {
	    *info = n + 1;
	}
	goto L70;
    }
//Compute Eigenvectors
//(Real Workspace: need 2*N)
//(Complex Workspace: need 2*N)
    if (ilv) {
	if (ilvl) {
	    if (ilvr) {
		chtemp = 'B';
	    } else {
		chtemp = 'L';
	    }
	} else {
	    chtemp = 'R';
	}
	Ctgevc((const char *) chtemp, "B", &ldummy, n, &A[0], lda, &B[0], ldb, &vl[0], ldvl, &vr[0], ldvr, n, &in, &work[iwrk], &rwork[irwrk], &ierr);
	if (ierr != 0) {
	    *info = n + 2;
	    goto L70;
	}
//Undo balancing on VL and VR and normalization
//(Workspace: none needed)
	if (ilvl) {
	    Cggbak("P", "L", n, ilo, ihi, &rwork[ileft], &rwork[iright], n, &vl[0], ldvl, &ierr);
	    for (jc = 1; jc <= n; jc++) {
		temp = Zero;
		for (jr = 1; jr <= n; jr++) {
		    mtemp1 = temp, mtemp2 = Cabs1(vl[jr + jc * ldvl]);
		    temp = max(mtemp1, mtemp2);
		}
		if (temp < smlnum) {
		    goto L30;
		}
		temp = One / temp;
		for (jr = 1; jr <= n; jr++) {
		    vl[jr + jc * ldvl] = vl[jr + jc * ldvl] * temp;
		}
	      L30:
		;
	    }
	}
	if (ilvr) {
	    Cggbak("P", "R", n, ilo, ihi, &rwork[ileft], &rwork[iright], n, &vr[0], ldvr, &ierr);
	    for (jc = 1; jc <= n; jc++) {
		temp = Zero;
		for (jr = 1; jr <= n; jr++) {
		    mtemp1 = temp, mtemp2 = Cabs1(vr[jr + jc * ldvr]);
		    temp = max(mtemp1, mtemp2);
		}
		if (temp < smlnum) {
		    goto L60;
		}
		temp = One / temp;
		for (jr = 1; jr <= n; jr++) {
		    vr[jr + jc * ldvr] = vr[jr + jc * ldvr] * temp;
		}
	      L60:
		;
	    }
	}
    }
//Undo scaling if necessary
    if (ilascl) {
	Clascl("G", 0, 0, anrmto, anrm, n, 1, &alpha[1], n, &ierr);
    }
    if (ilbscl) {
	Clascl("G", 0, 0, bnrmto, bnrm, n, 1, &beta[1], n, &ierr);
    }
  L70:
    work[1] = lwkopt;
    return;
}
