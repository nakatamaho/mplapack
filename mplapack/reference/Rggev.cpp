/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rggev.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rggev(const char *jobvl, const char *jobvr, INTEGER n, REAL *
	   A, INTEGER lda, REAL * B, INTEGER ldb, REAL * alphar,
	   REAL * alphai, REAL * beta, REAL * vl, INTEGER ldvl, REAL * vr, INTEGER ldvr, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER jc, in, jr, ihi, ilo;
    REAL eps;
    LOGICAL ilv;
    REAL anrm, bnrm;
    INTEGER ierr, itau;
    REAL temp;
    LOGICAL ilvl, ilvr;
    INTEGER iwrk;
    INTEGER ileft, icols, irows;
    LOGICAL ilascl, ilbscl;
    LOGICAL ldummy;
    char chtemp;
    REAL bignum;
    INTEGER ijobvl, iright, ijobvr;
    REAL anrmto = 0.0, bnrmto = 0.0;
    INTEGER minwrk, maxwrk;
    REAL smlnum;
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
	*info = -12;
    } else if (ldvr < 1 || (ilvr && ldvr < n)) {
	*info = -14;
    }
//Compute workspace
// (Note: Comments in the code beginning "Workspace:" describe the
//  minimal amount of workspace needed at that point in the code,
//  as well as the preferred amount for good performance.
//  NB refers to the optimal block size for the immediately
//  following subroutine, as returned by ILAENV. The workspace is
//  computed assuming ILO = 1 and IHI = N, the worst case.)
    if (*info == 0) {
	minwrk = max((INTEGER) 1, n * 8);
	maxwrk = max((INTEGER) 1, n * (iMlaenv(1, "Rgeqrf", " ", n, 1, n, 0) + 7));
	maxwrk = max(maxwrk, n * (iMlaenv(1, "Rormqr", " ", n, 1, n, 0) + 7));
	if (ilvl) {
	    maxwrk = max(maxwrk, n * (iMlaenv(1, "Rorgqr", " ", n, 1, n, -1) + 7));
	}
	work[0] = maxwrk;
	if (lwork < minwrk && !lquery) {
	    *info = -16;
	}
    }
    if (*info != 0) {
	Mxerbla("Rggev ", -(*info));
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
    //    Rlabad(&smlnum, &bignum);
    smlnum = sqrt(smlnum) / eps;
    bignum = One / smlnum;
//Scale A if max element outside range [SMLNUM,BIGNUM]
    anrm = Rlange("M", n, n, A, lda, work);
    ilascl = MFALSE;
    if (anrm > Zero && anrm < smlnum) {
	anrmto = smlnum;
	ilascl = MTRUE;
    } else if (anrm > bignum) {
	anrmto = bignum;
	ilascl = MTRUE;
    }
    if (ilascl) {
	Rlascl("G", 0, 0, anrm, anrmto, n, n, A, lda, &ierr);
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
	Rlascl("G", 0, 0, bnrm, bnrmto, n, n, B, ldb, &ierr);
    }
//Permute the matrices A, B to isolate eigenvalues if possible
//(Workspace: need 6*N)
    ileft = 1;
    iright = n + 1;
    iwrk = iright + n;
    Rggbal("P", n, A, lda, B, ldb, &ilo, &ihi, &work[ileft], &work[iright], &work[iwrk], &ierr);
//Reduce B to triangular form (QR decomposition of B)
//(Workspace: need N, prefer N*NB)
    irows = ihi + 1 - ilo;
    if (ilv) {
	icols = n + 1 - ilo;
    } else {
	icols = irows;
    }
    itau = iwrk;
    iwrk = itau + irows;
    Rgeqrf(irows, icols, &B[ilo + ilo * ldb], ldb, &work[itau], &work[iwrk], lwork + 1 - iwrk, &ierr);
//Apply the orthogonal transformation to matrix A
//(Workspace: need N, prefer N*NB)
    Rormqr("L", "T", irows, icols, irows, &B[ilo + ilo * ldb], ldb, &work[itau], &A[ilo + ilo * lda], lda, &work[iwrk], lwork + 1 - iwrk, &ierr);
//Initialize VL
//(Workspace: need N, prefer N*NB)
    if (ilvl) {
	Rlaset("Full", n, n, Zero, One, vl, ldvl);
	if (irows > 1) {
	    Rlacpy("L", irows - 1, irows - 1, &B[ilo + 1 + ilo * ldb], ldb, &vl[ilo + 1 + ilo * ldvl], ldvl);
	}
	Rorgqr(irows, irows, irows, &vl[ilo + ilo * ldvl], ldvl, &work[itau], &work[iwrk], lwork + 1 - iwrk, &ierr);
    }
//Initialize VR
    if (ilvr) {
	Rlaset("Full", n, n, Zero, One, vr, ldvr);
    }
//Reduce to generalized Hessenberg form
//(Workspace: none needed)
    if (ilv) {
//Eigenvectors requested -- work on whole matrix.
	Rgghrd(jobvl, jobvr, n, ilo, ihi, A, lda, B, ldb, vl, ldvl, vr, ldvr, &ierr);
    } else {
	Rgghrd("N", "N", irows, (INTEGER) 1, irows, &A[ilo + ilo * lda], lda, &B[ilo + ilo * ldb], ldb, vl, ldvl, vr, ldvr, &ierr);
    }
//Perform QZ algorithm (Compute eigenvalues, and optionally, the
//Schur forms and Schur vectors)
//(Workspace: need N)
    iwrk = itau;
    if (ilv) {
	chtemp = 'S';
    } else {
	chtemp = 'E';
    }
    Rhgeqz(&chtemp, jobvl, jobvr, n, ilo, ihi, A, lda, B, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, &work[iwrk], lwork + 1 - iwrk, &ierr);
    if (ierr != 0) {
	if (ierr > 0 && ierr <= n) {
	    *info = ierr;
	} else if (ierr > n && ierr <= n << 1) {
	    *info = ierr - n;
	} else {
	    *info = n + 1;
	}
	goto L110;
    }
//Compute Eigenvectors
//(Workspace: need 6*N)
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
	Rtgevc(&chtemp, "B", &ldummy, n, A, lda, B, ldb, &vl[0], ldvl, &vr[0], ldvr, n, &in, &work[iwrk], &ierr);
	if (ierr != 0) {
	    *info = n + 2;
	    goto L110;
	}
//Undo balancing on VL and VR and normalization
//(Workspace: none needed)
	if (ilvl) {
	    Rggbak("P", "L", n, ilo, ihi, &work[ileft], &work[iright], n, &vl[0], ldvl, &ierr);
	    for (jc = 1; jc <= n; jc++) {
		if (alphai[jc] < Zero) {
		    goto L50;
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
		    goto L50;
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
	      L50:
		;
	    }
	}
	if (ilvr) {
	    Rggbak("P", "R", n, ilo, ihi, &work[ileft], &work[iright], n, vr, ldvr, &ierr);
	    for (jc = 1; jc <= n; jc++) {
		if (alphai[jc] < Zero) {
		    goto L100;
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
		    goto L100;
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
	      L100:
		;
	    }
	}
//End of eigenvector calculation
    }
//Undo scaling if necessary
    if (ilascl) {
	Rlascl("G", 0, 0, anrmto, anrm, n, 1, alphar, n, &ierr);
	Rlascl("G", 0, 0, anrmto, anrm, n, 1, alphai, n, &ierr);
    }
    if (ilbscl) {
	Rlascl("G", 0, 0, bnrmto, bnrm, n, 1, beta, n, &ierr);
    }
  L110:
    work[0] = maxwrk;
    return;
}
