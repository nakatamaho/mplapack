/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgeevx.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rgeevx(const char *balanc, const char *jobvl, const char *jobvr, const char *sense, INTEGER n, REAL * A, INTEGER lda,
	    REAL * wr, REAL * wi, REAL * vl, INTEGER ldvl, REAL * vr, INTEGER ldvr, INTEGER * ilo, INTEGER * ihi, REAL * scale,
	    REAL * abnrm, REAL * rconde, REAL * rcondv, REAL * work, INTEGER lwork, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, k;
    REAL r, cs, sn;
    char job;
    REAL scl, dum, eps;
    char side = 'A';
    REAL anrm;
    INTEGER ierr, itau;
    INTEGER iwrk, nout;
    INTEGER icond;
    LOGICAL scalea;
    REAL cscale = 0.0;
    LOGICAL select;
    REAL bignum;
    INTEGER minwrk, maxwrk;
    LOGICAL wantvl, wntsnb;
    INTEGER hswork;
    LOGICAL wntsne;
    REAL smlnum;
    LOGICAL lquery, wantvr, wntsnn, wntsnv;
    REAL Zero = 0.0, One = 1.0;

//Test the input arguments
    *info = 0;
    lquery = lwork == -1;
    wantvl = Mlsame(jobvl, "V");
    wantvr = Mlsame(jobvr, "V");
    wntsnn = Mlsame(sense, "N");
    wntsne = Mlsame(sense, "E");
    wntsnv = Mlsame(sense, "V");
    wntsnb = Mlsame(sense, "B");
    if (!(Mlsame(balanc, "N") || Mlsame(balanc, "S") || Mlsame(balanc, "P")
	  || Mlsame(balanc, "B"))) {
	*info = -1;
    } else if (!wantvl && !Mlsame(jobvl, "N")) {
	*info = -2;
    } else if (!wantvr && !Mlsame(jobvr, "N")) {
	*info = -3;
    } else if (!(wntsnn || wntsne || wntsnb || wntsnv) || ((wntsne || wntsnb) && !(wantvl && wantvr))) {
	*info = -4;
    } else if (n < 0) {
	*info = -5;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -7;
    } else if (ldvl < 1 || (wantvl && ldvl < n)) {
	*info = -11;
    } else if (ldvr < 1 || (wantvr && ldvr < n)) {
	*info = -13;
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
	    maxwrk = n + n * iMlaenv(1, "Rgehrd", " ", n, 1, n, 0);
	    if (wantvl) {
		Rhseqr("S", "V", n, 1, n, &A[0], lda, &wr[1], &wi[1], &vl[0], ldvl, &work[0], -1, info);
	    } else if (wantvr) {
		Rhseqr("S", "V", n, 1, n, &A[0], lda, &wr[1], &wi[1], &vr[0], ldvr, &work[0], -1, info);
	    } else {
		if (wntsnn) {
		    Rhseqr("E", "N", n, 1, n, &A[0], lda, &wr[1], &wi[1], &vr[0], ldvr, &work[0], -1, info);
		} else {
		    Rhseqr("S", "N", n, 1, n, &A[0], lda, &wr[1], &wi[1], &vr[0], ldvr, &work[0], -1, info);
		}
	    }
	    hswork = (INTEGER) cast2double(work[1]);
	    if (!wantvl && !wantvr) {
		minwrk = n << 1;
		if (!wntsnn) {
		    minwrk = max(minwrk, n * n + n * 6);
		}
		maxwrk = max(maxwrk, hswork);
		if (!wntsnn) {
		    maxwrk = max(maxwrk, n * n + n * 6);
		}
	    } else {
		minwrk = n * 3;
		if (!wntsnn && !wntsne) {
		    minwrk = max(minwrk, n * n + n * 6);
		}
		maxwrk = max(maxwrk, hswork);
		maxwrk = max(maxwrk, n + (n - 1) * iMlaenv(1, "Rorghr", " ", n, 1, n, -1));
		if (!wntsnn && !wntsne) {
		    maxwrk = max(maxwrk, n * n + n * 6);
		}
		maxwrk = max(maxwrk, n * 3);
	    }
	    maxwrk = max(maxwrk, minwrk);
	}
	work[1] = maxwrk;
	if (lwork < minwrk && !lquery) {
	    *info = -21;
	}
    }
    if (*info != 0) {
	Mxerbla("Rgeevx", -(*info));
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
    icond = 0;
    anrm = Rlange("M", n, n, &A[0], lda, &dum);
    scalea = MFALSE;
    if (anrm > Zero && anrm < smlnum) {
	scalea = MTRUE;
	cscale = smlnum;
    } else if (anrm > bignum) {
	scalea = MTRUE;
	cscale = bignum;
    }
    if (scalea) {
	Rlascl("G", 0, 0, anrm, cscale, n, n, &A[0], lda, &ierr);
    }
//Balance the matrix and compute ABNRM
    Rgebal(balanc, n, &A[0], lda, ilo, ihi, &scale[1], &ierr);
    *abnrm = Rlange("1", n, n, &A[0], lda, &dum);
    if (scalea) {
	dum = *abnrm;
	Rlascl("G", 0, 0, cscale, anrm, 1, 1, &dum, 1, &ierr);
	*abnrm = dum;
    }
//Reduce to upper Hessenberg form
//(Workspace: need 2*N, prefer N+N*NB)
    itau = 1;
    iwrk = itau + n;
    Rgehrd(n, *ilo, *ihi, &A[0], lda, &work[itau], &work[iwrk], lwork - iwrk + 1, &ierr);
    if (wantvl) {
//Want left eigenvectors
//Copy Householder vectors to VL
	side = 'L';
	Rlacpy("L", n, n, &A[0], lda, &vl[0], ldvl);
//Generate orthogonal matrix in VL
//(Workspace: need 2*N-1, prefer N+(N-1)*NB)
	Rorghr(n, *ilo, *ihi, &vl[0], ldvl, &work[itau], &work[iwrk], lwork - iwrk + 1, &ierr);
//Perform QR iteration, accumulating Schur vectors in VL
//(Workspace: need 1, prefer HSWORK (see comments) )
	iwrk = itau;
	Rhseqr("S", "V", n, *ilo, *ihi, &A[0], lda, &wr[1], &wi[1], &vl[0], ldvl, &work[iwrk], lwork - iwrk + 1, info);
	if (wantvr) {
//Want left and right eigenvectors
//Copy Schur vectors to VR
	    side = 'B';
	    Rlacpy("F", n, n, &vl[0], ldvl, &vr[0], ldvr);
	}
    } else if (wantvr) {
//Want right eigenvectors
//Copy Householder vectors to VR
	side = 'R';
	Rlacpy("L", n, n, &A[0], lda, &vr[0], ldvr);
//Generate orthogonal matrix in VR
//(Workspace: need 2*N-1, prefer N+(N-1)*NB)
	Rorghr(n, *ilo, *ihi, &vr[0], ldvr, &work[itau], &work[iwrk], lwork - iwrk + 1, &ierr);
//Perform QR iteration, accumulating Schur vectors in VR
//(Workspace: need 1, prefer HSWORK (see comments) )
	iwrk = itau;
	Rhseqr("S", "V", n, *ilo, *ihi, &A[0], lda, &wr[1], &wi[1], &vr[0], ldvr, &work[iwrk], lwork - iwrk + 1, info);
    } else {
//Compute eigenvalues only
//If condition numbers desired, compute Schur form
	if (wntsnn) {
	    job = 'E';
	} else {
	    job = 'S';
	}
//(Workspace: need 1, prefer HSWORK (see comments) )
	iwrk = itau;
	Rhseqr((const char *) job, "N", n, *ilo, *ihi, &A[0], lda, &wr[1], &wi[1], &vr[0], ldvr, &work[iwrk], lwork - iwrk + 1, info);
    }
//If INFO > 0 from DHSEQR, then quit
    if (*info > 0) {
	goto L50;
    }
    if (wantvl || wantvr) {
//Compute left and/or right eigenvectors
//(Workspace: need 3*N)
	Rtrevc((const char *) side, "B", &select, n, &A[0], lda, &vl[0], ldvl, &vr[0], ldvr, n, &nout, &work[iwrk], &ierr);
    }
//Compute condition numbers if desired
//(Workspace: need N*N+6*N unless SENSE = 'E')
    if (!wntsnn) {
	Rtrsna(sense, "A", &select, n, &A[0], lda, &vl[0], ldvl, &vr[0], ldvr, &rconde[1], &rcondv[1], n, &nout, &work[iwrk], n, &iwork[1], &icond);
    }
    if (wantvl) {
//Undo balancing of left eigenvectors
	Rgebak(balanc, "L", n, *ilo, *ihi, &scale[1], n, &vl[0], ldvl, &ierr);
//Normalize left eigenvectors and make largest component real
	for (i = 0; i < n; i++) {
	    if (wi[i] == Zero) {
		scl = Zero / Rnrm2(n, &vl[i * ldvl + 1], 1);
		Rscal(n, scl, &vl[i * ldvl + 1], 1);
	    } else if (wi[i] > Zero) {
		scl = Zero / Rlapy2(Rnrm2(n, &vl[i * ldvl + 1], 1), Rnrm2(n, &vl[(i + 1) * ldvl + 1], 1));
		Rscal(n, scl, &vl[i * ldvl + 1], 1);
		Rscal(n, scl, &vl[(i + 1) * ldvl + 1], 1);
		for (k = 0; k < n; k++) {
		    work[k] = vl[k + i * ldvl] * vl[k + i * ldvl] + vl[k + (i + 1) * ldvl] * vl[k + (i + 1) * ldvl];
		}
		k = iRamax(n, &work[0], 1);
		Rlartg(vl[k + i * ldvl], vl[k + (i + 1) * ldvl], &cs, &sn, &r);
		Rrot(n, &vl[i * ldvl + 1], 1, &vl[(i + 1) * ldvl + 1], 1, cs, sn);
		vl[k + (i + 1) * ldvl] = Zero;
	    }
	}
    }
    if (wantvr) {
//Undo balancing of right eigenvectors
	Rgebak(balanc, "R", n, *ilo, *ihi, &scale[1], n, &vr[0], ldvr, &ierr);
//Normalize right eigenvectors and make largest component real
	for (i = 0; i < n; i++) {
	    if (wi[i] == Zero) {
		scl = Zero / Rnrm2(n, &vr[i * ldvr + 1], 1);
		Rscal(n, scl, &vr[i * ldvr + 1], 1);
	    } else if (wi[i] > Zero) {
		scl = Zero / Rlapy2(Rnrm2(n, &vr[i * ldvr + 1], 1), Rnrm2(n, &vr[(i + 1) * ldvr + 1], 1));
		Rscal(n, scl, &vr[i * ldvr + 1], 1);
		Rscal(n, scl, &vr[(i + 1) * ldvr + 1], 1);
		for (k = 0; k < n; k++) {
		    work[k] = vr[k + i * ldvr] * vr[k + i * ldvr] + vr[k + (i + 1) * ldvr] * vr[k + (i + 1) * ldvr];
		}
		k = iRamax(n, &work[0], 1);
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
	if (*info == 0) {
	    if ((wntsnv || wntsnb) && icond == 0) {
		Rlascl("G", 0, 0, cscale, anrm, n, 1, &rcondv[1], n, &ierr);
	    }
	} else {
	    Rlascl("G", 0, 0, cscale, anrm, *ilo - 1, 1, &wr[1], n, &ierr);
	    Rlascl("G", 0, 0, cscale, anrm, *ilo - 1, 1, &wi[1], n, &ierr);
	}
    }
    work[1] = maxwrk;
    return;
}
