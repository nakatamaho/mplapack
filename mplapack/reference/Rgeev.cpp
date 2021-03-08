/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgeev.cpp,v 1.11 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Rgeev(const char *jobvl, const char *jobvr, INTEGER n, REAL *
      A, INTEGER lda, REAL * wr, REAL * wi, REAL * vl, INTEGER ldvl, REAL * vr, INTEGER ldvr, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i, k;
    REAL r, cs, sn;
    INTEGER ihi;
    REAL scl;
    INTEGER ilo;
    REAL dum[1], eps;
    INTEGER ibal;
    char side[1];
    REAL anrm;
    INTEGER ierr, itau;
    INTEGER iwrk, nout;
    INTEGER scalea;
    REAL cscale = 0.0;
    REAL bignum;
    INTEGER select[1];
    INTEGER minwrk, maxwrk;
    INTEGER wantvl;
    REAL smlnum;
    INTEGER hswork;
    INTEGER lquery, wantvr;

    REAL mtemp1, mtemp2;
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
	*info = -9;
    } else if (ldvr < 1 || (wantvr && ldvr < n)) {
	*info = -11;
    }
//Compute workspace
// (Note: Comments in the code beginning "Workspace:" describe the
//  minimal amount of workspace needed at that point in the code,
//  as well as the preferred amount for good performance.
//  NB refers to the optimal block size for the immediately
//  following subroutine, as returned by ILAENV.
//  HSWORK refers to the workspace preferred by DHSEQR, as
//  calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
//  the worst case.)
    if (*info == 0) {
	if (n == 0) {
	    minwrk = 0;
	    maxwrk = 0;
	} else {
	    maxwrk = (n << 1) + n * iMlaenv(1, "Rgehrd", " ", n, 1, n, 0);
	    if (wantvl) {
		minwrk = n * 2;
		maxwrk = max(maxwrk, (n * 1) + (n - 1) * iMlaenv(1, "Rorghr", " ", n, 1, n, -1));
		Rhseqr("S", "V", n, 1, n, A, lda, wr, wi, vl, ldvl, work, -1, info);
		hswork = (INTEGER) cast2double(work[1]);
		maxwrk = max(max(maxwrk, n + 1), n + hswork);
		maxwrk = max(maxwrk, n * 4);
	    } else if (wantvr) {
		minwrk = n * 4;
		maxwrk = max(maxwrk, (n * 2) + (n - 1) * iMlaenv(1, "Rorghr", " ", n, 1, n, -1));
		Rhseqr("S", "V", n, 1, n, A, lda, wr, wi, vr, ldvr, work, -1, info);
		hswork = (INTEGER) cast2double(work[1]);
		maxwrk = max(max(maxwrk, n + 1), n + hswork);
		maxwrk = max(maxwrk, n * 4);
	    } else {
		minwrk = n * 3;
		Rhseqr("E", "N", n, 1, n, A, lda, wr, wi, vr, ldvr, work, -1, info);
		hswork = (INTEGER) cast2double(work[1]);
		maxwrk = max(max(maxwrk, n + 1), n + hswork);
	    }
	    maxwrk = max(maxwrk, minwrk);
	}
	work[1] = maxwrk;
	if (lwork < minwrk && !lquery) {
	    *info = -13;
	}
    }
    if (*info != 0) {
	Mxerbla("Rgeev ", -(*info));
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
    smlnum = sqrt(smlnum) / eps;
    bignum = One / smlnum;
//Scale A if max element outside range [SMLNUM,BIGNUM]
    anrm = Rlange("M", n, n, A, lda, dum);
    scalea = MFALSE;
    if (anrm > Zero && anrm < smlnum) {
	scalea = MTRUE;
	cscale = smlnum;
    } else if (anrm > bignum) {
	scalea = MTRUE;
	cscale = bignum;
    }
    if (scalea) {
	Rlascl("G", 0, 0, anrm, cscale, n, n, A, lda, &ierr);
    }
//Balance the matrix
//(Workspace: need N)
    ibal = 0;
    Rgebal("B", n, A, lda, &ilo, &ihi, &work[ibal], &ierr);
//Reduce to upper Hessenberg form
//(Workspace: need 3*N, prefer 2*N+N*NB)
    itau = ibal + n;
    iwrk = itau + n;
    Rgehrd(n, ilo, ihi, A, lda, &work[itau], &work[iwrk], lwork - iwrk + 1, &ierr);
    if (wantvl) {
//Want left eigenvectors
//Copy Householder vectors to VL
	side[0] = 'L';
	Rlacpy("L", n, n, A, lda, vl, ldvl);
//Generate orthogonal matrix in VL
//(Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
	Rorghr(n, ilo, ihi, vl, ldvl, &work[itau], &work[iwrk], lwork - iwrk + 1, &ierr);
//Perform QR iteration, accumulating Schur vectors in VL
//(Workspace: need N+1, prefer N+HSWORK (see comments) )
	iwrk = itau;
	Rhseqr("S", "V", n, ilo, ihi, A, lda, wr, wi, vl, ldvl, &work[iwrk], lwork - iwrk + 1, info);
	if (wantvr) {
//Want left and right eigenvectors
//Copy Schur vectors to VR
	    side[0] = 'B';
	    Rlacpy("F", n, n, vl, ldvl, vr, ldvr);
	}
    } else if (wantvr) {
//Want right eigenvectors
//Copy Householder vectors to VR
	side[0] = 'R';
	Rlacpy("L", n, n, A, lda, vr, ldvr);
//Generate orthogonal matrix in VR
//(Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
	Rorghr(n, ilo, ihi, vr, ldvr, &work[itau], &work[iwrk], lwork - iwrk + 1, &ierr);
//Perform QR iteration, accumulating Schur vectors in VR
//(Workspace: need N+1, prefer N+HSWORK (see comments) )
	iwrk = itau;
	Rhseqr("S", "V", n, ilo, ihi, A, lda, wr, wi, vr, ldvr, &work[iwrk], lwork - iwrk + 1, info);
    } else {
//Compute eigenvalues only
//(Workspace: need N+1, prefer N+HSWORK (see comments) )
	iwrk = itau;
	Rhseqr("E", "N", n, ilo, ihi, A, lda, wr, wi, vr, ldvr, &work[iwrk], lwork - iwrk + 1, info);
    }
//If INFO > 0 from DHSEQR, then quit
    if (*info > 0) {
	goto L50;
    }
    if (wantvl || wantvr) {
//Compute left and/or right eigenvectors
//(Workspace: need 4*N)
	Rtrevc(side, "B", select, n, A, lda, vl, ldvl, vr, ldvr, n, &nout, &work[iwrk], &ierr);
    }
    if (wantvl) {
//Undo balancing of left eigenvectors
//(Workspace: need N)
	Rgebak("B", "L", n, ilo, ihi, &work[ibal], n, vl, ldvl, &ierr);
//Normalize left eigenvectors and make largest component real
	for (i = 0; i < n; i++) {
	    if (wi[i] == Zero) {
		scl = Zero / Rnrm2(n, &vl[i * ldvl + 1], 1);
		Rscal(n, scl, &vl[i * ldvl + 1], 1);
	    } else if (wi[i] > Zero) {
		mtemp1 = Rnrm2(n, &vl[i * ldvl + 1], 1);
		mtemp2 = Rnrm2(n, &vl[(i + 1) * ldvl + 1], 1);
		scl = Zero / Rlapy2(mtemp1, mtemp2);
		Rscal(n, scl, &vl[i * ldvl + 1], 1);
		Rscal(n, scl, &vl[(i + 1) * ldvl + 1], 1);
		for (k = 0; k < n; k++) {
		    work[iwrk + k - 1] = vl[k + i * ldvl] * vl[k + i * ldvl] + vl[k + (i + 1) * ldvl] * vl[k + (i + 1) * ldvl];
		}
		k = iRamax(n, &work[iwrk], 1);
		Rlartg(vl[k + i * ldvl], vl[k + (i + 1) * ldvl], &cs, &sn, &r);
		Rrot(n, &vl[i * ldvl + 1], 1, &vl[(i + 1) * ldvl + 1], 1, cs, sn);
		vl[k + (i + 1) * ldvl] = Zero;
	    }
	}
    }
    if (wantvr) {
//Undo balancing of right eigenvectors
//(Workspace: need N)
	Rgebak("B", "R", n, ilo, ihi, &work[ibal], n, vr, ldvr, &ierr);
//Normalize right eigenvectors and make largest component real
	for (i = 0; i < n; i++) {
	    if (wi[i] == Zero) {
		scl = Zero / Rnrm2(n, &vr[i * ldvr + 1], 1);
		Rscal(n, scl, &vr[i * ldvr + 1], 1);
	    } else if (wi[i] > Zero) {
		mtemp1 = Rnrm2(n, &vr[i * ldvr + 1], 1);
		mtemp2 = Rnrm2(n, &vr[(i + 1) * ldvr + 1], 1);
		scl = Zero / Rlapy2(mtemp1, mtemp2);
		Rscal(n, scl, &vr[i * ldvr + 1], 1);
		Rscal(n, scl, &vr[(i + 1) * ldvr + 1], 1);
		for (k = 0; k < n; k++) {
		    work[iwrk + k - 1] = vr[k + i * ldvr] * vr[k + i * ldvr] + vr[k + (i + 1) * ldvr] * vr[k + (i + 1) * ldvr];
		}
		k = iRamax(n, &work[iwrk], 1);
		Rlartg(vr[k + i * ldvr], vr[k + (i + 1) * ldvr], &cs, &sn, &r);
		Rrot(n, &vr[i * ldvr + 1], 1, &vr[(i + 1) * ldvr + 1], 1, cs, sn);
		vr[k + (i + 1) * ldvr] = Zero;
	    }
	}
    }
//Undo scaling if necessary
  L50:
    if (scalea) {
	Rlascl("G", 0, 0, cscale, anrm, n - *info, 1, &wr[*info + 1], max(n - *info, (INTEGER) 1), &ierr);
	Rlascl("G", 0, 0, cscale, anrm, n - *info, 1, &wi[*info + 1], max(n - *info, (INTEGER) 1), &ierr);
	if (*info > 0) {
	    Rlascl("G", 0, 0, cscale, anrm, ilo - 1, 1, wr, n, &ierr);
	    Rlascl("G", 0, 0, cscale, anrm, ilo - 1, 1, wi, n, &ierr);
	}
    }
    work[1] = (REAL) double (maxwrk);
    return;
}
