/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rggevx.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void
Rggevx(const char *balanc, const char *jobvl, const char *jobvr, const char *sense, INTEGER n, REAL * A,
       INTEGER lda, REAL * B, INTEGER ldb, REAL * alphar, REAL * alphai, REAL * beta,
       REAL * vl, INTEGER ldvl, REAL * vr, INTEGER ldvr,
       INTEGER * ilo, INTEGER * ihi, REAL * lscale, REAL * rscale,
       REAL * abnrm, REAL * bbnrm, REAL * rconde, REAL * rcondv, REAL * work, INTEGER lwork, INTEGER * iwork, LOGICAL * bwork, INTEGER * info)
{
    INTEGER i, j, m, jc, in, mm, jr;
    REAL eps;
    LOGICAL ilv, pair;
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
    INTEGER minwrk, maxwrk;
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
    if (!(Mlsame(balanc, "N") || Mlsame(balanc, "S") || Mlsame(balanc, "P")
	  || Mlsame(balanc, "B"))) {
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
	*info = -14;
    } else if (ldvr < 1 || (ilvr && ldvr < n)) {
	*info = -16;
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
	    if (noscl && !ilv) {
		minwrk = n << 1;
	    } else {
		minwrk = n * 6;
	    }
	    if (wantse || wantsb) {
		minwrk = n * 10;
	    }
	    if (wantsv || wantsb) {
		minwrk = max(minwrk, (n << 1) * (n + 4) + 16);
	    }
	    maxwrk = minwrk;
	    maxwrk = max(maxwrk, n + n * iMlaenv(1, "Rgeqrf", " ", n, 1, n, 0));
	    maxwrk = max(maxwrk, n + n * iMlaenv(1, "Rormqr", " ", n, 1, n, 0));
	    if (ilvl) {
		maxwrk = max(maxwrk, n + n * iMlaenv(1, "Rorgqr", " ", n, 1, n, 0));
	    }
	}
	work[1] = maxwrk;
	if (lwork < minwrk && !lquery) {
	    *info = -26;
	}
    }
    if (*info != 0) {
	Mxerbla("Rggevx", -(*info));
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
    anrm = Rlange("M", n, n, &A[0], lda, &work[0]);
    ilascl = MFALSE;
    if (anrm > Zero && anrm < smlnum) {
	anrmto = smlnum;
	ilascl = MTRUE;
    } else if (anrm > bignum) {
	anrmto = bignum;
	ilascl = MTRUE;
    }
    if (ilascl) {
	Rlascl("G", 0, 0, anrm, anrmto, n, n, &A[0], lda, &ierr);
    }
//Scale B if max element outside range [SMLNUM,BIGNUM]
    bnrm = Rlange("M", n, n, &B[0], ldb, &work[0]);
    ilbscl = MFALSE;
    if (bnrm > Zero && bnrm < smlnum) {
	bnrmto = smlnum;
	ilbscl = MTRUE;
    } else if (bnrm > bignum) {
	bnrmto = bignum;
	ilbscl = MTRUE;
    }
    if (ilbscl) {
	Rlascl("G", 0, 0, bnrm, bnrmto, n, n, &B[0], ldb, &ierr);
    }
//Permute and/or balance the matrix pair (A,B)
//(Workspace: need 6*N if BALANC = 'S' or 'B', 1 otherwise)
    Rggbal(balanc, n, &A[0], lda, &B[0], ldb, ilo, ihi, &lscale[1], &rscale[1], &work[0], &ierr);
//Compute ABNRM and BBNRM
    *abnrm = Rlange("1", n, n, &A[0], lda, &work[0]);
    if (ilascl) {
	work[1] = *abnrm;
	Rlascl("G", 0, 0, anrmto, anrm, 1, 1, &work[0], 1, &ierr);
	*abnrm = work[1];
    }
    *bbnrm = Rlange("1", n, n, &B[0], ldb, &work[0]);
    if (ilbscl) {
	work[1] = *bbnrm;
	Rlascl("G", 0, 0, bnrmto, bnrm, 1, 1, &work[0], 1, &ierr);
	*bbnrm = work[1];
    }
//Reduce B to triangular form (QR decomposition of B)
//(Workspace: need N, prefer N*NB )
    irows = *ihi + 1 - *ilo;
    if (ilv || !wantsn) {
	icols = n + 1 - *ilo;
    } else {
	icols = irows;
    }
    itau = 1;
    iwrk = itau + irows;
    Rgeqrf(irows, icols, &B[*ilo + *ilo * ldb], ldb, &work[itau], &work[iwrk], lwork + 1 - iwrk, &ierr);
//Apply the orthogonal transformation to A
//(Workspace: need N, prefer N*NB)
    Rormqr("L", "T", irows, icols, irows, &B[*ilo + *ilo * ldb], ldb, &work[itau], &A[*ilo + *ilo * lda], lda, &work[iwrk], lwork + 1 - iwrk, &ierr);
//Initialize VL and/or VR
//(Workspace: need N, prefer N*NB)
    if (ilvl) {
	Rlaset("Full", n, n, Zero, One, &vl[0], ldvl);
	if (irows > 1) {
	    Rlacpy("L", irows - 1, irows - 1, &B[*ilo + 1 + *ilo * ldb], ldb, &vl[*ilo + 1 + *ilo * ldvl], ldvl);
	}
	Rorgqr(irows, irows, irows, &vl[*ilo + *ilo * ldvl], ldvl, &work[itau], &work[iwrk], lwork + 1 - iwrk, &ierr);
    }
    if (ilvr) {
	Rlaset("Full", n, n, Zero, One, &vr[0], ldvr);
    }
//Reduce to generalized Hessenberg form
//(Workspace: none needed)
    if (ilv || !wantsn) {
//Eigenvectors requested -- work on whole matrix.
	Rgghrd(jobvl, jobvr, n, *ilo, *ihi, &A[0], lda, &B[0], ldb, &vl[0], ldvl, &vr[0], ldvr, &ierr);
    } else {
	Rgghrd("N", "N", irows, 1, irows, &A[*ilo + *ilo * lda], lda, &B[*ilo + *ilo * ldb], ldb, &vl[0], ldvl, &vr[0], ldvr, &ierr);
    }
//Perform QZ algorithm (Compute eigenvalues, and optionally, the
//Schur forms and Schur vectors)
//(Workspace: need N)
    if (ilv || !wantsn) {
	chtemp = 'S';
    } else {
	chtemp = 'E';
    }
    Rhgeqz((const char *) chtemp, jobvl, jobvr, n, *ilo, *ihi, &A[0], lda, &B[0]
	   , ldb, &alphar[1], &alphai[1], &beta[1], &vl[0], ldvl, &vr[0], ldvr, &work[0], lwork, &ierr);
    if (ierr != 0) {
	if (ierr > 0 && ierr <= n) {
	    *info = ierr;
	} else if (ierr > n && ierr <= n << 1) {
	    *info = ierr - n;
	} else {
	    *info = n + 1;
	}
	goto L130;
    }
//Compute Eigenvectors and estimate condition numbers if desired
//(Workspace: DTGEVC: need 6*N
//            DTGSNA: need 2*N*(N+2)+16 if SENSE = 'V' or 'B',
//                    need N otherwise )
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
	    Rtgevc((const char *) chtemp, "B", &ldummy, n, &A[0], lda, &B[0], ldb, &vl[0], ldvl, &vr[0], ldvr, n, &in, &work[0], &ierr);
	    if (ierr != 0) {
		*info = n + 2;
		goto L130;
	    }
	}
	if (!wantsn) {
//compute eigenvectors (DTGEVC) and estimate condition
//numbers (DTGSNA). Note that the definition of the condition
//number is not invariant under transformation (u,v) to
//(Q*u, Z*v), where (u,v) are eigenvectors of the generalized
//Schur form (S,T), Q and Z are orthogonal matrices. In order
//to avoid using extra 2*N*N workspace, we have to recalculate
//eigenvectors and estimate one condition numbers at a time.
	    pair = MFALSE;
	    for (i = 0; i < n; i++) {
		if (pair) {
		    pair = MFALSE;
		    goto L20;
		}
		mm = 1;
		if (i < n) {
		    if (A[i + 1 + i * lda] != Zero) {
			pair = MTRUE;
			mm = 2;
		    }
		}
		for (j = 0; j < n; j++) {
		    bwork[j] = MFALSE;
		}
		if (mm == 1) {
		    bwork[i] = MTRUE;
		} else if (mm == 2) {
		    bwork[i] = MTRUE;
		    bwork[i + 1] = MTRUE;
		}
		iwrk = mm * n + 1;
		iwrk1 = iwrk + mm * n;
//Compute a pair of left and right eigenvectors.
//(compute workspace: need up to 4*N + 6*N)
		if (wantse || wantsb) {
		    Rtgevc("B", "S", &bwork[1], n, &A[0], lda, &B[0], ldb, &work[0], n, &work[iwrk], n, mm, &m, &work[iwrk1], &ierr);
		    if (ierr != 0) {
			*info = n + 2;
			goto L130;
		    }
		}
		Rtgsna(sense, "S", &bwork[1], n, &A[0], lda, &B[0], ldb, &work[0], n, &work[iwrk], n, &rconde[i], &rcondv[i], mm,
		       &m, &work[iwrk1], lwork - iwrk1 + 1, &iwork[1], &ierr);
	      L20:
		;
	    }
	}
    }
//Undo balancing on VL and VR and normalization
//(Workspace: none needed)
    if (ilvl) {
	Rggbak(balanc, "L", n, *ilo, *ihi, &lscale[1], &rscale[1], n, &vl[0], ldvl, &ierr);
	for (jc = 1; jc <= n; jc++) {
	    if (alphai[jc] < Zero) {
		goto L70;
	    }
	    temp = Zero;
	    if (alphai[jc] == Zero) {
		for (jr = 1; jr <= n; jr++) {
		    mtemp1 = temp, mtemp2 = abs(vl[jr + jc * ldvl]);
		    temp = max(mtemp1, mtemp2);
		}
	    } else {
		for (jr = 1; jr <= n; jr++) {
		    mtemp1 = temp, mtemp2 = abs(vl[jr + jc * ldvl]) + abs(vl[jr + (jc + 1) * ldvl]);
		    temp = max(mtemp1, mtemp2);
		}
	    }
	    if (temp < smlnum) {
		goto L70;
	    }
	    temp = One / temp;
	    if (alphai[jc] == Zero) {
		for (jr = 1; jr <= n; jr++) {
		    vl[jr + jc * ldvl] = vl[jr + jc * ldvl] * temp;
		}
	    } else {
		for (jr = 1; jr <= n; jr++) {
		    vl[jr + jc * ldvl] = vl[jr + jc * ldvl] * temp;
		    vl[jr + (jc + 1) * ldvl] = vl[jr + (jc + 1) * ldvl] * temp;
		}
	    }
	  L70:
	    ;
	}
    }
    if (ilvr) {
	Rggbak(balanc, "R", n, *ilo, *ihi, &lscale[1], &rscale[1], n, &vr[0], ldvr, &ierr);
	for (jc = 1; jc <= n; jc++) {
	    if (alphai[jc] < Zero) {
		goto L120;
	    }
	    temp = Zero;
	    if (alphai[jc] == Zero) {
		for (jr = 1; jr <= n; jr++) {
		    mtemp1 = temp, mtemp2 = abs(vr[jr + jc * ldvr]);
		    temp = max(mtemp1, mtemp2);
		}
	    } else {
		for (jr = 1; jr <= n; jr++) {
		    mtemp1 = temp, mtemp2 = abs(vr[jr + jc * ldvr]) + abs(vr[jr + (jc + 1) * ldvr]);
		    temp = max(mtemp1, mtemp2);
		}
	    }
	    if (temp < smlnum) {
		goto L120;
	    }
	    temp = One / temp;
	    if (alphai[jc] == Zero) {
		for (jr = 1; jr <= n; jr++) {
		    vr[jr + jc * ldvr] = vr[jr + jc * ldvr] * temp;
		}
	    } else {
		for (jr = 1; jr <= n; jr++) {
		    vr[jr + jc * ldvr] = vr[jr + jc * ldvr] * temp;
		    vr[jr + (jc + 1) * ldvr] = vr[jr + (jc + 1) * ldvr] * temp;
		}
	    }
	  L120:
	    ;
	}
    }
//Undo scaling if necessary
    if (ilascl) {
	Rlascl("G", 0, 0, anrmto, anrm, n, 1, &alphar[1], n, &ierr);
	Rlascl("G", 0, 0, anrmto, anrm, n, 1, &alphai[1], n, &ierr);
    }
    if (ilbscl) {
	Rlascl("G", 0, 0, bnrmto, bnrm, n, 1, &beta[1], n, &ierr);
    }
  L130:
    work[1] = maxwrk;
    return;
}
