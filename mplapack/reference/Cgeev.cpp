/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgeev.cpp,v 1.7 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cgeev(const char *jobvl, const char *jobvr, INTEGER n,
	   COMPLEX * A, INTEGER lda, COMPLEX * w, COMPLEX * vl, INTEGER ldvl, COMPLEX * vr, INTEGER ldvr, COMPLEX * work, INTEGER lwork, REAL * rwork, INTEGER * info)
{
    INTEGER i, k, ihi;
    REAL scl;
    INTEGER ilo;
    REAL dum[1], eps;
    COMPLEX tmp;
    INTEGER ibal;
    char side = 'A';
    REAL anrm;
    INTEGER ierr, itau, iwrk, nout;
    LOGICAL scalea;
    REAL cscale = 0.0;
    LOGICAL select[1];
    REAL bignum;
    INTEGER minwrk, maxwrk;
    LOGICAL wantvl;
    REAL smlnum;
    INTEGER hswork, irwork;
    LOGICAL lquery, wantvr;
    REAL Zero = 0.0, One = 1.0;

//Test the input arguments
    *info = 0;
    lquery = lwork == -1;
    wantvl = Mlsame(jobvl, "V");
    wantvr = Mlsame(jobvr, "V");
    if (!wantvl && !Mlsame(jobvl, "N")) {
	*info = -1;
    } else if (!wantvr && !Mlsame(jobvr, "N")) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -5;
    } else if (ldvl < 1 || (wantvl && ldvl < n)) {
	*info = -8;
    } else if (ldvr < 1 || (wantvr && ldvr < n)) {
	*info = -10;
    }
//Compute workspace
// (Note: Comments in the code beginning "Workspace:" describe the
//  minimal amount of workspace needed at that point in the code,
//  as well as the preferred amount for good performance.
//  CWorkspace refers to complex workspace, and RWorkspace to real
//  workspace. NB refers to the optimal block size for the
//  immediately following subroutine, as returned by ILAENV.
//  HSWORK refers to the workspace preferred by ZHSEQR, as
//  calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
//  the worst case.)
    if (*info == 0) {
	if (n == 0) {
	    minwrk = 0;
	    maxwrk = 0;
	} else {
	    maxwrk = n + n * iMlaenv(1, "Cgehrd", " ", n, 1, n, 0);
	    minwrk = n << 1;
	    if (wantvl) {
		maxwrk = max(maxwrk, n + (n - 1) * iMlaenv(1, "Cunghr", " ", n, 1, n, -1));
		Chseqr("S", "V", n, 1, n, &A[0], lda, &w[1], &vl[0], ldvl, &work[0], -1, info);
	    } else if (wantvr) {
		maxwrk = max(maxwrk, n + (n - 1) * iMlaenv(1, "ZUNGHR", " ", n, 1, n, -1));
		Chseqr("S", "V", n, 1, n, &A[0], lda, &w[1], &vr[0], ldvr, &work[0], -1, info);
	    } else {
		Chseqr("E", "N", n, 1, n, &A[0], lda, &w[1], &vr[0], ldvr, &work[0], -1, info);
	    }
	    hswork = (INTEGER) cast2double(work[1].real());
	    maxwrk = max(max(maxwrk, hswork), minwrk);
	}
	work[1] = maxwrk;
	if (lwork < minwrk && !lquery) {
	    *info = -12;
	}
    }
    if (*info != 0) {
	Mxerbla("Cgeev ", -(*info));
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
    anrm = Clange("M", n, n, &A[0], lda, dum);
    scalea = MFALSE;
    if (anrm > Zero && anrm < smlnum) {
	scalea = MTRUE;
	cscale = smlnum;
    } else if (anrm > bignum) {
	scalea = MTRUE;
	cscale = bignum;
    }
    if (scalea) {
	Clascl("G", 0, 0, anrm, cscale, n, n, &A[0], lda, &ierr);
    }
//Balance the matrix
//(CWorkspace: none)
//(RWorkspace: need N)
    ibal = 0;
    Cgebal("B", n, &A[0], lda, &ilo, &ihi, &rwork[ibal], &ierr);
//Reduce to upper Hessenberg form
//(CWorkspace: need 2*N, prefer N+N*NB)
//(RWorkspace: none)
    itau = 1;
    iwrk = itau + n;
    Cgehrd(n, ilo, ihi, &A[0], lda, &work[itau], &work[iwrk], lwork - iwrk + 1, &ierr);
    if (wantvl) {
//Want left eigenvectors
//Copy Householder vectors to VL
	side = 'L';
	Clacpy("L", n, n, &A[0], lda, &vl[0], ldvl);
//Generate unitary matrix in VL
//(CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
//(RWorkspace: none)
	Cunghr(n, ilo, ihi, &vl[0], ldvl, &work[itau], &work[iwrk], lwork - iwrk + 1, &ierr);
//Perform QR iteration, accumulating Schur vectors in VL
//(CWorkspace: need 1, prefer HSWORK (see comments) )
//(RWorkspace: none)
	iwrk = itau;
	Chseqr("S", "V", n, ilo, ihi, &A[0], lda, &w[1], &vl[0], ldvl, &work[iwrk], lwork - iwrk + 1, info);
	if (wantvr) {
//Want left and right eigenvectors
//Copy Schur vectors to VR
	    side = 'B';
	    Clacpy("F", n, n, &vl[0], ldvl, &vr[0], ldvr);
	}
    } else if (wantvr) {
//Want right eigenvectors
//Copy Householder vectors to VR
	side = 'R';
	Clacpy("L", n, n, &A[0], lda, &vr[0], ldvr);
//Generate unitary matrix in VR
//(CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
//(RWorkspace: none)
	Cunghr(n, ilo, ihi, &vr[0], ldvr, &work[itau], &work[iwrk], lwork - iwrk + 1, &ierr);
//Perform QR iteration, accumulating Schur vectors in VR
//(CWorkspace: need 1, prefer HSWORK (see comments) )
//(RWorkspace: none)
	iwrk = itau;
	Chseqr("S", "V", n, ilo, ihi, &A[0], lda, &w[1], &vr[0], ldvr, &work[iwrk], lwork - iwrk + 1, info);
    } else {
//Compute eigenvalues only
//(CWorkspace: need 1, prefer HSWORK (see comments) )
//(RWorkspace: none)
	iwrk = itau;
	Chseqr("E", "N", n, ilo, ihi, &A[0], lda, &w[1], &vr[0], ldvr, &work[iwrk], lwork - iwrk + 1, info);
    }
// If INFO > 0 from ZHSEQR, then quit
    if (*info > 0) {
	goto L50;
    }
    if (wantvl || wantvr) {
//Compute left and/or right eigenvectors
//(CWorkspace: need 2*N)
//(RWorkspace: need 2*N)
	irwork = ibal + n;
	Ctrevc((const char *) side, "B", select, n, &A[0], lda, &vl[0], ldvl, &vr[0], ldvr, n, &nout, &work[iwrk], &rwork[irwork], &ierr);
    }
    if (wantvl) {
//Undo balancing of left eigenvectors
//(CWorkspace: none)
//(RWorkspace: need N)
	Cgebak("B", "L", n, ilo, ihi, &rwork[ibal], n, &vl[0], ldvl, &ierr);
//Normalize left eigenvectors and make largest component real
	for (i = 0; i < n; i++) {
	    scl = Zero / RCnrm2(n, &vl[i * ldvl + 1], 1);
	    CRscal(n, scl, &vl[i * ldvl + 1], 1);
	    for (k = 0; k < n; k++) {
		rwork[irwork + k - 1] = vl[k + i * ldvl].real() * vl[k + i * ldvl].real() + vl[k + i * ldvl].imag() * vl[k + i * ldvl].imag();
	    }
	    k = iRamax(n, &rwork[irwork], 1);
	    //tmp = conj(vl[k + i * ldvl]) / (COMPLEX)sqrt(rwork[irwork + k - 1]);
	    Cscal(n, tmp, &vl[i * ldvl + 1], 1);
	    vl[k + i * ldvl] = vl[k + i * ldvl].real();
	}
    }
    if (wantvr) {
//Undo balancing of right eigenvectors
//(CWorkspace: none)
//(RWorkspace: need N)
	Cgebak("B", "R", n, ilo, ihi, &rwork[ibal], n, &vr[0], ldvr, &ierr);
//Normalize right eigenvectors and make largest component real
	for (i = 0; i < n; i++) {
	    scl = Zero / RCnrm2(n, &vr[i * ldvr + 1], 1);
	    CRscal(n, scl, &vr[i * ldvr + 1], 1);
	    for (k = 0; k < n; k++) {
		rwork[irwork + k - 1] = vr[k + i * ldvr].real() * vr[k + i * ldvr].real() + vr[k + i * ldvr].imag() * vr[k + i * ldvr].imag();
	    }
	    k = iRamax(n, &rwork[irwork], 1);
        //	    tmp = conj(vr[k + i * ldvr]) / sqrt(rwork[irwork + k - 1]);
	    Cscal(n, tmp, &vr[i * ldvr + 1], 1);
	    vr[k + i * ldvr] = vr[k + i * ldvr].real();
	}
    }
//Undo scaling if necessary
  L50:
    if (scalea) {
	Clascl("G", 0, 0, cscale, anrm, n - *info, 1, &w[*info + 1], max(n - *info, (INTEGER) 1), &ierr);
	if (*info > 0) {
	    Clascl("G", 0, 0, cscale, anrm, ilo - 1, 1, &w[1], n, &ierr);
	}
    }
    work[1] = maxwrk;
    return;
}
