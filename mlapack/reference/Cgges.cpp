/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgges.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

#define MTRUE 1
#define MFALSE 0

void Cgges(const char *jobvsl, const char *jobvsr, const char *sort, LFP selctg, INTEGER n, COMPLEX * A, INTEGER lda, COMPLEX * B,
	   INTEGER ldb, INTEGER * sdim, COMPLEX * alpha, COMPLEX * beta, COMPLEX * vsl, INTEGER ldvsl, COMPLEX * vsr, INTEGER ldvsr,
	   COMPLEX * work, INTEGER lwork, REAL * rwork, LOGICAL * bwork, INTEGER * info)
{
    INTEGER i;
    REAL dif[2];
    INTEGER ihi, ilo;
    REAL eps, anrm, bnrm;
    INTEGER idum[1], ierr, itau, iwrk;
    REAL pvsl, pvsr;
    INTEGER ileft, icols;
    LOGICAL cursl, ilvsl, ilvsr;
    INTEGER irwrk, irows;
    LOGICAL ilascl, ilbscl;
    REAL bignum;
    INTEGER ijobvl, iright;
    INTEGER ijobvr;
    REAL anrmto = 0.0;
    INTEGER lwkmin;
    LOGICAL lastsl;
    REAL bnrmto = 0.0;
    REAL smlnum;
    LOGICAL wantst, lquery;
    INTEGER lwkopt;
    REAL Zero = 0.0, One = 1.0;

//Decode the input arguments
    if (Mlsame(jobvsl, "N")) {
	ijobvl = 0;
	ilvsl = MFALSE;
    } else if (Mlsame(jobvsl, "V")) {
	ijobvl = 2;
	ilvsl = MTRUE;
    } else {
	ijobvl = -1;
	ilvsl = MFALSE;
    }
    if (Mlsame(jobvsr, "N")) {
	ijobvr = 1;
	ilvsr = MFALSE;
    } else if (Mlsame(jobvsr, "V")) {
	ijobvr = 2;
	ilvsr = MTRUE;
    } else {
	ijobvr = -1;
	ilvsr = MFALSE;
    }
    wantst = Mlsame(sort, "S");
//Test the input arguments
    *info = 0;
    lquery = lwork == -1;
    if (ijobvl <= 0) {
	*info = -1;
    } else if (ijobvr <= 0) {
	*info = -2;
    } else if (!wantst && !Mlsame(sort, "N")) {
	*info = -3;
    } else if (n < 0) {
	*info = -5;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -7;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -9;
    } else if (ldvsl < 1 || (ilvsl && ldvsl < n)) {
	*info = -14;
    } else if (ldvsr < 1 || (ilvsr && ldvsr < n)) {
	*info = -16;
    }
//Compute workspace
// (Note: Comments in the code beginning "Workspace:" describe the
//  minimal amount of workspace needed at that point in the code,
//  as well as the preferred amount for good performance.
//  NB refers to the optimal block size for the immediately
//  following subroutine, as returned by ILAENV.)
    if (*info == 0) {
	lwkmin = max((INTEGER) 1, n << 1);
	lwkopt = max((INTEGER) 1, n + n * iMlaenv(1, "ZGEQRF", " ", n, 1, n, 0));
	lwkopt = max(lwkopt, n + n * iMlaenv(1, "ZUNMQR", " ", n, 1, n, -1));
	if (ilvsl) {
	    lwkopt = max(lwkopt, n + n * iMlaenv(1, "ZUNGQR", " ", n, 1, n, -1));
	}
	work[1] = lwkopt;
	if (lwork < lwkmin && !lquery) {
	    *info = -18;
	}
    }
    if (*info != 0) {
	Mxerbla("Cgges ", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	*sdim = 0;
	return;
    }
//Get machine constants
    eps = Rlamch("P");
    smlnum = Rlamch("S");
    bignum = One / smlnum;
    //  Rlabad(&smlnum, &bignum);
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
//Permute the matrix to make it more nearly triangular
//(Real Workspace: need 6*N)
    ileft = 1;
    iright = n + 1;
    irwrk = iright + n;
    Cggbal("P", n, &A[0], lda, &B[0], ldb, &ilo, &ihi, &rwork[ileft], &rwork[iright], &rwork[irwrk], &ierr);
//Reduce B to triangular form (QR decomposition of B)
//(Complex Workspace: need N, prefer N*NB)
    irows = ihi + 1 - ilo;
    icols = n + 1 - ilo;
    itau = 1;
    iwrk = itau + irows;
    Cgeqrf(irows, icols, &B[ilo + ilo * ldb], ldb, &work[itau], &work[iwrk], lwork + 1 - iwrk, &ierr);
//Apply the orthogonal transformation to matrix A
//(Complex Workspace: need N, prefer N*NB)
    Cunmqr("L", "C", irows, icols, irows, &B[ilo + ilo * ldb], ldb, &work[itau], &A[ilo + ilo * lda], lda, &work[iwrk], lwork + 1 - iwrk, &ierr);
//Initialize VSL
//(Complex Workspace: need N, prefer N*NB)
    if (ilvsl) {
	Claset("Full", n, n, Zero, One, &vsl[0], ldvsl);
	if (irows > 1) {
	    Clacpy("L", irows - 1, irows - 1, &B[ilo + 1 + ilo * ldb], ldb, &vsl[ilo + 1 + ilo * ldvsl], ldvsl);
	}
	Cungqr(irows, irows, irows, &vsl[ilo + ilo * ldvsl], ldvsl, &work[itau], &work[iwrk], lwork + 1 - iwrk, &ierr);
    }
//Initialize VSR
    if (ilvsr) {
	Claset("Full", n, n, Zero, One, &vsr[0], ldvsr);
    }
//Reduce to generalized Hessenberg form
//(Workspace: none needed)
    Cgghrd(jobvsl, jobvsr, n, ilo, ihi, &A[0], lda, &B[0], ldb, &vsl[0], ldvsl, &vsr[0], ldvsr, &ierr);
    *sdim = 0;
//Perform QZ algorithm, computing Schur vectors if desired
//(Complex Workspace: need N)
//(Real Workspace: need N)
    iwrk = itau;
    Chgeqz("S", jobvsl, jobvsr, n, ilo, ihi, &A[0], lda, &B[0], ldb, &alpha[1], &beta[1], &vsl[0], ldvsl, &vsr[0], ldvsr, &work[iwrk], lwork + 1 - iwrk, &rwork[irwrk], &ierr);
    if (ierr != 0) {
	if (ierr > 0 && ierr <= n) {
	    *info = ierr;
	} else if (ierr > n && ierr <= n << 1) {
	    *info = ierr - n;
	} else {
	    *info = n + 1;
	}
	goto L30;
    }
//Sort eigenvalues ALPHA/BETA if desired
//(Workspace: none needed)
    if (wantst) {
//Undo scaling on eigenvalues before selecting
	if (ilascl) {
	    Clascl("G", 0, 0, anrm, anrmto, n, 1, &alpha[1], n, &ierr);
	}
	if (ilbscl) {
	    Clascl("G", 0, 0, bnrm, bnrmto, n, 1, &beta[1], n, &ierr);
	}
//Select eigenvalues
	for (i = 0; i < n; i++) {
	    bwork[i] = (*selctg) (&alpha[i], &beta[i]);
	}
	Ctgsen(0, ilvsl, ilvsr, &bwork[1], n, &A[0], lda, &B[0], ldb, &alpha[1], &beta[1], &vsl[0], ldvsl,
	       &vsr[0], ldvsr, sdim, &pvsl, &pvsr, dif, &work[iwrk], lwork - iwrk + 1, idum, 1, &ierr);
	if (ierr == 1) {
	    *info = n + 3;
	}
    }
//Apply back-permutation to VSL and VSR
//(Workspace: none needed)
    if (ilvsl) {
	Cggbak("P", "L", n, ilo, ihi, &rwork[ileft], &rwork[iright], n, &vsl[0], ldvsl, &ierr);
    }
    if (ilvsr) {
	Cggbak("P", "R", n, ilo, ihi, &rwork[ileft], &rwork[iright], n, &vsr[0], ldvsr, &ierr);
    }
//Undo scaling
    if (ilascl) {
	Clascl("U", 0, 0, anrmto, anrm, n, n, &A[0], lda, &ierr);
	Clascl("G", 0, 0, anrmto, anrm, n, 1, &alpha[1], n, &ierr);
    }
    if (ilbscl) {
	Clascl("U", 0, 0, bnrmto, bnrm, n, n, &B[0], ldb, &ierr);
	Clascl("G", 0, 0, bnrmto, bnrm, n, 1, &beta[1], n, &ierr);
    }
    if (wantst) {
//Check if reordering is correct
	lastsl = MTRUE;
	*sdim = 0;
	for (i = 0; i < n; i++) {
	    cursl = (*selctg) (&alpha[i], &beta[i]);
	    if (cursl) {
		++(*sdim);
	    }
	    if (cursl && !lastsl) {
		*info = n + 2;
	    }
	    lastsl = cursl;
	}
    }
  L30:
    work[1] = lwkopt;
    return;
}
