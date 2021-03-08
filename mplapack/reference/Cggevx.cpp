/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cggevx.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cggevx(const char *balanc, const char *jobvl, const char *jobvr, const char *sense, INTEGER n, COMPLEX * A, INTEGER lda,
	    COMPLEX * B, INTEGER ldb, COMPLEX * alpha, COMPLEX * beta, COMPLEX * vl, INTEGER ldvl, COMPLEX * vr, INTEGER ldvr,
	    INTEGER * ilo, INTEGER * ihi, REAL * lscale, REAL * rscale, REAL * abnrm, REAL * bbnrm, REAL * rconde, REAL * rcondv,
	    COMPLEX * work, INTEGER lwork, REAL * rwork, INTEGER * iwork, LOGICAL * bwork, INTEGER * info)
{
    INTEGER i, j, m, jc, in, jr;
    REAL eps;
    LOGICAL ilv;
    REAL anrm, bnrm;
    INTEGER ierr, itau;
    REAL temp;
    LOGICAL ilvl, ilvr;
    INTEGER iwrk, iwrk1;
    INTEGER icols;
    LOGICAL noscl;
    INTEGER irows;
    LOGICAL ilascl, ilbscl;
    LOGICAL ldummy;
    char chtemp;
    REAL bignum;
    INTEGER ijobvl;
    INTEGER ijobvr;
    LOGICAL wantsb;
    REAL anrmto = 0.0;
    LOGICAL wantse;
    REAL bnrmto = 0.0;
    INTEGER minwrk;
    INTEGER maxwrk;
    LOGICAL wantsn;
    REAL smlnum;
    LOGICAL lquery, wantsv;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1, mtemp2;

//Decode the input arguments
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
    noscl = Mlsame(balanc, "N") || Mlsame(balanc, "P");
    wantsn = Mlsame(sense, "N");
    wantse = Mlsame(sense, "E");
    wantsv = Mlsame(sense, "V");
    wantsb = Mlsame(sense, "B");
//Test the input arguments
    *info = 0;
    lquery = lwork == -1;
    if (!(noscl || Mlsame(balanc, "S") || Mlsame(balanc, "B"))) {
	*info = -1;
    } else if (ijobvl <= 0) {
	*info = -2;
    } else if (ijobvr <= 0) {
	*info = -3;
    } else if (!(wantsn || wantse || wantsb || wantsv)) {
	*info = -4;
    } else if (n < 0) {
	*info = -5;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -7;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -9;
    } else if (ldvl < 1 || (ilvl && ldvl < n)) {
	*info = -13;
    } else if (ldvr < 1 || (ilvr && ldvr < n)) {
	*info = -15;
    }
//Compute workspace
// (Note: Comments in the code beginning "Workspace:" describe the
//  minimal amount of workspace needed at that point in the code,
//  as well as the preferred amount for good performance.
//  NB refers to the optimal block size for the immediately
//  following subroutine, as returned by ILAENV. The workspace is
//  computed assuming ILO = 1 and IHI = N, the worst case.)
    if (*info == 0) {
	if (n == 0) {
	    minwrk = 0;
	    maxwrk = 0;
	} else {
	    minwrk = n << 1;
	    if (wantse) {
		minwrk = n << 2;
	    } else if (wantsv || wantsb) {
		minwrk = (n << 1) * (n + 1);
	    }
	    maxwrk = minwrk;
	    maxwrk = max(maxwrk, n + n * iMlaenv(1, "Cgeqrf", " ", n, 1, n, 0));
	    maxwrk = max(maxwrk, n + n * iMlaenv(1, "Cunmqr", " ", n, 1, n, 0));
	    if (ilvl) {
		maxwrk = max(maxwrk, n + n * iMlaenv(1, "Cungqr", " ", n, 1, n, 0));
	    }
	}
	work[1] = maxwrk;
	if (lwork < minwrk && !lquery) {
	    *info = -25;
	}
    }
    if (*info != 0) {
	Mxerbla("Cggevx", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//Get machine constants
    eps = Rlamch("P");
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
//Permute and/or balance the matrix pair (A,B)
//(Real Workspace: need 6*N if BALANC = 'S' or 'B', 1 otherwise)
    Cggbal(balanc, n, &A[0], lda, &B[0], ldb, ilo, ihi, &lscale[1], &rscale[1], &rwork[1], &ierr);
//Compute ABNRM and BBNRM
    *abnrm = Clange("1", n, n, &A[0], lda, &rwork[1]);
    if (ilascl) {
	rwork[0] = *abnrm;
	Rlascl("G", 0, 0, anrmto, anrm, 1, 1, &rwork[1], 1, &ierr);
	*abnrm = rwork[0];
    }
    *bbnrm = Clange("1", n, n, &B[0], ldb, &rwork[1]);
    if (ilbscl) {
	rwork[0] = *bbnrm;
	Rlascl("G", 0, 0, bnrmto, bnrm, 1, 1, &rwork[1], 1, &ierr);
	*bbnrm = rwork[0];
    }
//Reduce B to triangular form (QR decomposition of B)
//(Complex Workspace: need N, prefer N*NB )
    irows = *ihi + 1 - *ilo;
    if (ilv || !wantsn) {
	icols = n + 1 - *ilo;
    } else {
	icols = irows;
    }
    itau = 1;
    iwrk = itau + irows;
    Cgeqrf(irows, icols, &B[*ilo + *ilo * ldb], ldb, &work[itau], &work[iwrk], lwork + 1 - iwrk, &ierr);
//Apply the unitary transformation to A
//(Complex Workspace: need N, prefer N*NB)
    Cunmqr("L", "C", irows, icols, irows, &B[*ilo + *ilo * ldb], ldb, &work[itau], &A[*ilo + *ilo * lda], lda, &work[iwrk], lwork + 1 - iwrk, &ierr);
//Initialize VL and/or VR
//(Workspace: need N, prefer N*NB)
    if (ilvl) {
	Claset("Full", n, n, Zero, One, &vl[0], ldvl);
	if (irows > 1) {
	    Clacpy("L", irows - 1, irows - 1, &B[*ilo + 1 + *ilo * ldb], ldb, &vl[*ilo + 1 + *ilo * ldvl], ldvl);
	}
	Cungqr(irows, irows, irows, &vl[*ilo + *ilo * ldvl], ldvl, &work[itau], &work[iwrk], lwork + 1 - iwrk, &ierr);
    }
    if (ilvr) {
	Claset("Full", n, n, Zero, One, &vr[0], ldvr);
    }
//Reduce to generalized Hessenberg form
//(Workspace: none needed)
    if (ilv || !wantsn) {
//Eigenvectors requested -- work on whole matrix.
	Cgghrd(jobvl, jobvr, n, *ilo, *ihi, &A[0], lda, &B[0], ldb, &vl[0], ldvl, &vr[0], ldvr, &ierr);
    } else {
	Cgghrd("N", "N", irows, 1, irows, &A[*ilo + *ilo * lda], lda, &B[*ilo + *ilo * ldb], ldb, &vl[0], ldvl, &vr[0], ldvr, &ierr);
    }
//Perform QZ algorithm (Compute eigenvalues, and optionally, the
//Schur forms and Schur vectors)
//(Complex Workspace: need N)
//(Real Workspace: need N)
    iwrk = itau;
    if (ilv || !wantsn) {
	chtemp = 'S';
    } else {
	chtemp = 'E';
    }
    Chgeqz((const char *) chtemp, jobvl, jobvr, n, *ilo, *ihi, &A[0], lda, &B[0]
	   , ldb, &alpha[1], &beta[1], &vl[0], ldvl, &vr[0], ldvr, &work[iwrk], lwork + 1 - iwrk, &rwork[1], &ierr);
    if (ierr != 0) {
	if (ierr > 0 && ierr <= n) {
	    *info = ierr;
	} else if (ierr > n && ierr <= n << 1) {
	    *info = ierr - n;
	} else {
	    *info = n + 1;
	}
	goto L90;
    }
//Compute Eigenvectors and estimate condition numbers if desired
//ZTGEVC: (Complex Workspace: need 2*N )
//        (Real Workspace:    need 2*N )
//ZTGSNA: (Complex Workspace: need 2*N*N if SENSE='V' or 'B')
//        (Integer Workspace: need N+2 )
    if (ilv || !wantsn) {
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
	    Ctgevc((const char *) chtemp, "B", &ldummy, n, &A[0], lda, &B[0], ldb, &vl[0], ldvl, &vr[0], ldvr, n, &in, &work[iwrk], &rwork[1], &ierr);
	    if (ierr != 0) {
		*info = n + 2;
		goto L90;
	    }
	}
	if (!wantsn) {
//compute eigenvectors (DTGEVC) and estimate condition
//numbers (DTGSNA). Note that the definition of the condition
//number is not invariant under transformation (u,v) to
//(Q*u, Z*v), where (u,v) are eigenvectors of the generalized
//Schur form (S,T), Q and Z are orthogonal matrices. In order
//to avoid using extra 2*N*N workspace, we have to
//re-calculate eigenvectors and estimate the condition numbers
//one at a time.
	    for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
		    bwork[j] = MFALSE;
		}
		bwork[i] = MTRUE;
		iwrk = n + 1;
		iwrk1 = iwrk + n;
		if (wantse || wantsb) {
		    Ctgevc("B", "S", &bwork[1], n, &A[0], lda, &B[0], ldb, &work[0], n, &work[iwrk], n, 1, &m, &work[iwrk1], &rwork[1], &ierr);
		    if (ierr != 0) {
			*info = n + 2;
			goto L90;
		    }
		}
		Ctgsna(sense, "S", &bwork[1], n, &A[0], lda, &B[0], ldb, &work[0], n, &work[iwrk], n, &rconde[i], &rcondv[i], 1, &m,
		       &work[iwrk1], lwork - iwrk1 + 1, &iwork[1], &ierr);
	    }
	}
    }
//Undo balancing on VL and VR and normalization
//(Workspace: none needed)
    if (ilvl) {
	Cggbak(balanc, "L", n, *ilo, *ihi, &lscale[1], &rscale[1], n, &vl[0], ldvl, &ierr);
	for (jc = 1; jc <= n; jc++) {
	    temp = Zero;
	    for (jr = 1; jr <= n; jr++) {
		mtemp1 = temp, mtemp2 = Cabs1(vl[jr + jc * ldvl]);
		temp = max(mtemp1, mtemp2);
	    }
	    if (temp < smlnum) {
		goto L50;
	    }
	    temp = One / temp;
	    for (jr = 1; jr <= n; jr++) {
		vl[jr + jc * ldvl] = vl[jr + jc * ldvl] * temp;
	    }
	  L50:
	    ;
	}
    }
    if (ilvr) {
	Cggbak(balanc, "R", n, *ilo, *ihi, &lscale[1], &rscale[1], n, &vr[0], ldvr, &ierr);
	for (jc = 1; jc <= n; jc++) {
	    temp = Zero;
	    for (jr = 1; jr <= n; jr++) {
		mtemp1 = temp, mtemp2 = Cabs1(vr[jr + jc * ldvr]);
		temp = max(mtemp1, mtemp2);
	    }
	    if (temp < smlnum) {
		goto L80;
	    }
	    temp = One / temp;
	    for (jr = 1; jr <= n; jr++) {
		vr[jr + jc * ldvr] = vr[jr + jc * ldvr] * temp;
	    }
	  L80:
	    ;
	}
    }
//Undo scaling if necessary
    if (ilascl) {
	Clascl("G", 0, 0, anrmto, anrm, n, 1, &alpha[1], n, &ierr);
    }
    if (ilbscl) {
	Clascl("G", 0, 0, bnrmto, bnrm, n, 1, &beta[1], n, &ierr);
    }
  L90:
    work[1] = maxwrk;
    return;
}
