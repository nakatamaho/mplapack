/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgeesx.cpp,v 1.7 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cgeesx(const char *jobvs, const char *sort, LFP select, const char *sense, INTEGER n, COMPLEX * A, INTEGER lda, INTEGER * sdim,
	    COMPLEX * w, COMPLEX * vs, INTEGER ldvs, REAL * rconde, REAL * rcondv, COMPLEX * work, INTEGER lwork, REAL * rwork, LOGICAL * bwork, INTEGER * info)
{
    INTEGER i, ihi, ilo;
    REAL dum[1], eps;
    INTEGER ibal;
    REAL anrm;
    INTEGER ierr, itau, iwrk, lwrk, icond, ieval;
    INTEGER scalea;
    REAL cscale = 0.0;
    REAL bignum;
    INTEGER wantsb, wantse;
    INTEGER minwrk, maxwrk = 0;
    INTEGER wantsn;
    REAL smlnum;
    INTEGER hswork;
    INTEGER wantst, wantsv, wantvs;
    REAL Zero = 0.0, One = 1.0;

//Test the input arguments
    *info = 0;
    wantvs = Mlsame(jobvs, "V");
    wantst = Mlsame(sort, "S");
    wantsn = Mlsame(sense, "N");
    wantse = Mlsame(sense, "E");
    wantsv = Mlsame(sense, "V");
    wantsb = Mlsame(sense, "B");
    if (!wantvs && !Mlsame(jobvs, "N")) {
	*info = -1;
    } else if (!wantst && !Mlsame(sort, "N")) {
	*info = -2;
    } else if (!(wantsn || wantse || wantsv || wantsb) || (!wantst && !wantsn)) {
	*info = -4;
    } else if (n < 0) {
	*info = -5;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -7;
    } else if (ldvs < 1 || (wantvs && ldvs < n)) {
	*info = -11;
    }
//Compute workspace
// (Note: Comments in the code beginning "Workspace:" describe the
//  minimal amount of real workspace needed at that point in the
//  code, as well as the preferred amount for good performance.
//  CWorkspace refers to complex workspace, and RWorkspace to real
//  workspace. NB refers to the optimal block size for the
//  immediately following subroutine, as returned by ILAENV.
//  HSWORK refers to the workspace preferred by ZHSEQR, as
//  calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
//  the worst case.
//  If SENSE = 'E', 'V' or 'B', then the amount of workspace needed
//  depends on SDIM, which is computed by the routine ZTRSEN later
//  in the code.)
    if (*info == 0) {
	if (n == 0) {
	    minwrk = 0;
	    lwrk = 0;
	} else {
	    maxwrk = n + n * iMlaenv(1, "Cgehrd", " ", n, 1, n, 0);
	    minwrk = n << 1;
	    Chseqr("S", jobvs, n, 1, n, &A[0], lda, &w[1], &vs[0], ldvs, &work[0], -1, &ieval);
	    hswork = (INTEGER) cast2double(work[1].real());
	    if (!wantvs) {
		maxwrk = max(maxwrk, hswork);
	    } else {
		maxwrk = max(maxwrk, n + (n - 1) * iMlaenv(1, "Cunghr", " ", n, 1, n, -1));
		maxwrk = max(maxwrk, hswork);
	    }
	    lwrk = maxwrk;
	    if (!wantsn) {
		lwrk = max(lwrk, n * n / 2);
	    }
	}
	work[1] = lwrk;
	if (lwork < minwrk) {
	    *info = -15;
	}
    }
    if (*info != 0) {
	Mxerbla("Cgeesx", -(*info));
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
//Permute the matrix to make it more nearly triangular
//(CWorkspace: none)
//(RWorkspace: need N)
    ibal = 0;
    Cgebal("P", n, &A[0], lda, &ilo, &ihi, &rwork[ibal], &ierr);
//Reduce to upper Hessenberg form
//(CWorkspace: need 2*N, prefer N+N*NB)
//(RWorkspace: none)
    itau = 1;
    iwrk = n + itau;
    Cgehrd(n, ilo, ihi, &A[0], lda, &work[itau], &work[iwrk], lwork - iwrk + 1, &ierr);
    if (wantvs) {
//Copy Householder vectors to VS
	Clacpy("L", n, n, &A[0], lda, &vs[0], ldvs);
//Generate unitary matrix in VS
//(CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
//(RWorkspace: none)
	Cunghr(n, ilo, ihi, &vs[0], ldvs, &work[itau], &work[iwrk], lwork - iwrk + 1, &ierr);
    }
    *sdim = 0;
//Perform QR iteration, accumulating Schur vectors in VS if desired
//(CWorkspace: need 1, prefer HSWORK (see comments) )
//(RWorkspace: none) */
    iwrk = itau;
    Chseqr("S", jobvs, n, ilo, ihi, &A[0], lda, &w[1], &vs[0], ldvs, &work[iwrk], lwork - iwrk + 1, &ieval);
    if (ieval > 0) {
	*info = ieval;
    }
//Sort eigenvalues if desired
    if (wantst && *info == 0) {
	if (scalea) {
	    Clascl("G", 0, 0, cscale, anrm, n, 1, &w[1], n, &ierr);
	}
	for (i = 0; i < n; i++) {
	    bwork[i] = (*select) (&w[i]);
	}
//Reorder eigenvalues, transform Schur vectors, and compute
//reciprocal condition numbers
//(CWorkspace: if SENSE is not 'N', need 2*SDIM*(N-SDIM)
//             otherwise, need none )
//(RWorkspace: none)
	Ctrsen(sense, jobvs, &bwork[1], n, &A[0], lda, &vs[0], ldvs, &w[1], *sdim, rconde, rcondv, &work[iwrk], lwork - iwrk + 1, &icond);
	if (!wantsn) {
	    maxwrk = max(maxwrk, (*sdim << 1) * (n - *sdim));
	}
	if (icond == -14) {
//Not enough complex workspace
	    *info = -15;
	}
    }
    if (wantvs) {
//Undo balancing
//(CWorkspace: none)
//(RWorkspace: need N)
	Cgebak("P", "R", n, ilo, ihi, &rwork[ibal], n, &vs[0], ldvs, &ierr);
    }
    if (scalea) {
//Undo scaling for the Schur form of A
	Clascl("U", 0, 0, cscale, anrm, n, n, &A[0], lda, &ierr);
	Ccopy(n, &A[0], lda + 1, &w[1], 1);
	if ((wantsv || wantsb) && *info == 0) {
	    dum[0] = *rcondv;
	    Rlascl("G", 0, 0, cscale, anrm, 1, 1, dum, 1, &ierr);
	    *rcondv = dum[0];
	}
    }
    work[1] = maxwrk;
    return;
}
