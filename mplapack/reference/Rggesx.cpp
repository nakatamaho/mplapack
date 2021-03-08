/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rggesx.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rggesx(const char *jobvsl, const char *jobvsr, const char *sort, LFP selctg,
	    const char *sense, INTEGER n, REAL * A, INTEGER lda,
	    REAL * B, INTEGER ldb, INTEGER * sdim, REAL * alphar,
	    REAL * alphai, REAL * beta, REAL * vsl, INTEGER ldvsl,
	    REAL * vsr, INTEGER ldvsr, REAL * rconde, REAL * rcondv, REAL * work, INTEGER lwork, INTEGER * iwork, INTEGER liwork, LOGICAL * bwork, INTEGER * info)
{
    INTEGER i, ip;
    REAL pl, pr, dif[2];
    INTEGER ihi, ilo;
    REAL eps;
    INTEGER ijob = 0;
    REAL anrm, bnrm;
    INTEGER ierr, itau, iwrk, lwrk;
    INTEGER ileft, icols;
    LOGICAL cursl, ilvsl, ilvsr;
    INTEGER irows;
    LOGICAL lst2sl;
    LOGICAL ilascl, ilbscl;
    REAL safmin;
    REAL safmax;
    REAL bignum;
    INTEGER ijobvl, iright;
    INTEGER ijobvr;
    LOGICAL wantsb;
    INTEGER liwmin;
    LOGICAL wantse, lastsl;
    REAL anrmto = 0.0, bnrmto = 0.0;
    INTEGER minwrk, maxwrk;
    LOGICAL wantsn;
    REAL smlnum;
    LOGICAL wantst, lquery, wantsv;
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
    wantsn = Mlsame(sense, "N");
    wantse = Mlsame(sense, "E");
    wantsv = Mlsame(sense, "V");
    wantsb = Mlsame(sense, "B");
    lquery = lwork == -1 || liwork == -1;
    if (wantsn) {
	ijob = 0;
    } else if (wantse) {
	ijob = 1;
    } else if (wantsv) {
	ijob = 2;
    } else if (wantsb) {
	ijob = 4;
    }
//Test the input arguments
    *info = 0;
    if (ijobvl <= 0) {
	*info = -1;
    } else if (ijobvr <= 0) {
	*info = -2;
    } else if (!wantst && !Mlsame(sort, "N")) {
	*info = -3;
    } else if (!(wantsn || wantse || wantsv || wantsb) || (!wantst && !wantsn)) {
	*info = -5;
    } else if (n < 0) {
	*info = -6;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -8;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -10;
    } else if (ldvsl < 1 || (ilvsl && ldvsl < n)) {
	*info = -16;
    } else if (ldvsr < 1 || (ilvsr && ldvsr < n)) {
	*info = -18;
    }
//Compute workspace
// (Note: Comments in the code beginning "Workspace:" describe the
//  minimal amount of workspace needed at that point in the code,
//  as well as the preferred amount for good performance.
//  NB refers to the optimal block size for the immediately
//  following subroutine, as returned by ILAENV.)
    if (*info == 0) {
	if (n > 0) {
	    minwrk = max(n << 3, n * 6 + 16);
	    maxwrk = minwrk - n + n * iMlaenv(1, "DGEQRF", " ", n, 1, n, 0);
	    maxwrk = max(maxwrk, minwrk - n + n * iMlaenv(1, "DORMQR", " ", n, 1, n, -1));
	    if (ilvsl) {
		maxwrk = max(maxwrk, minwrk - n + n * iMlaenv(1, "DOR" "GQR", " ", n, 1, n, -1));
	    }
	    lwrk = maxwrk;
	    if (ijob >= 1) {
		lwrk = max(lwrk, n * n / 2);
	    }
	} else {
	    minwrk = 0;
	    maxwrk = 0;
	    lwrk = 0;
	}
	work[1] = lwrk;
	if (wantsn || n == 0) {
	    liwmin = 1;
	} else {
	    liwmin = n + 6;
	}
	iwork[1] = liwmin;

	if (lwork < minwrk && !lquery) {
	    *info = -22;
	} else if (liwork < liwmin && !lquery) {
	    *info = -24;
	}
    }
    if (*info != 0) {
	Mxerbla("Rggesx", -(*info));
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
    safmin = Rlamch("S");
    safmax = One / safmin;
    //Rlabad(&safmin, &safmax);
    smlnum = sqrt(safmin) / eps;
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
//Permute the matrix to make it more nearly triangular
//(Workspace: need 6*N + 2*N for permutation parameters)
    ileft = 1;
    iright = n + 1;
    iwrk = iright + n;
    Rggbal("P", n, &A[0], lda, &B[0], ldb, &ilo, &ihi, &work[ileft], &work[iright], &work[iwrk], &ierr);
//Reduce B to triangular form (QR decomposition of B)
//(Workspace: need N, prefer N*NB)
    irows = ihi + 1 - ilo;
    icols = n + 1 - ilo;
    itau = iwrk;
    iwrk = itau + irows;
    Rgeqrf(irows, icols, &B[ilo + ilo * ldb], ldb, &work[itau], &work[iwrk], lwork + 1 - iwrk, &ierr);
//Apply the orthogonal transformation to matrix A
//(Workspace: need N, prefer N*NB)
    Rormqr("L", "T", irows, icols, irows, &B[ilo + ilo * ldb], ldb, &work[itau], &A[ilo + ilo * lda], lda, &work[iwrk], lwork + 1 - iwrk, &ierr);
//Initialize VSL
//(Workspace: need N, prefer N*NB)
    if (ilvsl) {
	Rlaset("Full", n, n, Zero, One, &vsl[0], ldvsl);
	if (irows > 1) {
	    Rlacpy("L", irows - 1, irows - 1, &B[ilo + 1 + ilo * ldb], ldb, &vsl[ilo + 1 + ilo * ldvsl], ldvsl);
	}
	Rorgqr(irows, irows, irows, &vsl[ilo + ilo * ldvsl], ldvsl, &work[itau], &work[iwrk], lwork + 1 - iwrk, &ierr);
    }
//Initialize VSR
    if (ilvsr) {
	Rlaset("Full", n, n, Zero, One, &vsr[0], ldvsr);
    }
//Reduce to generalized Hessenberg form
//(Workspace: none needed)
    Rgghrd(jobvsl, jobvsr, n, ilo, ihi, &A[0], lda, &B[0], ldb, &vsl[0], ldvsl, &vsr[0], ldvsr, &ierr);
    *sdim = 0;
//Perform QZ algorithm, computing Schur vectors if desired
//(Workspace: need N)
    iwrk = itau;
    Rhgeqz("S", jobvsl, jobvsr, n, ilo, ihi, &A[0], lda, &B[0], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[0]
	   , ldvsl, &vsr[0], ldvsr, &work[iwrk], lwork + 1 - iwrk, &ierr);
    if (ierr != 0) {
	if (ierr > 0 && ierr <= n) {
	    *info = ierr;
	} else if (ierr > n && ierr <= n << 1) {
	    *info = ierr - n;
	} else {
	    *info = n + 1;
	}
	goto L60;
    }
//Sort eigenvalues ALPHA/BETA and compute the reciprocal of
//condition number(s)
//(Workspace: If IJOB >= 1, need MAX( 8*(N+1), 2*SDIM*(N-SDIM) )
//            otherwise, need 8*(N+1) )
    if (wantst) {
//Undo scaling on eigenvalues before SELCTGing
	if (ilascl) {
	    Rlascl("G", 0, 0, anrmto, anrm, n, 1, &alphar[1], n, &ierr);
	    Rlascl("G", 0, 0, anrmto, anrm, n, 1, &alphai[1], n, &ierr);
	}
	if (ilbscl) {
	    Rlascl("G", 0, 0, bnrmto, bnrm, n, 1, &beta[1], n, &ierr);
	}
//Select eigenvalues
	for (i = 0; i < n; i++) {
	    bwork[i] = (*selctg) (&alphar[i], &alphai[i], &beta[i]);
	}
//Reorder eigenvalues, transform Generalized Schur vectors, and
//compute reciprocal condition numbers
	Rtgsen(ijob, ilvsl, ilvsr, &bwork[1], n, &A[0], lda, &B[0], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[0], ldvsl, &vsr[0],
	       ldvsr, sdim, &pl, &pr, dif, &work[iwrk], lwork - iwrk + 1, &iwork[1], liwork, &ierr);
	if (ijob >= 1) {
	    maxwrk = max(maxwrk, (*sdim << 1) * (n - *sdim));
	}
	if (ierr == -22) {
//not enough real workspace
	    *info = -22;
	} else {
	    if (ijob == 1 || ijob == 4) {
		rconde[1] = pl;
		rconde[2] = pr;
	    }
	    if (ijob == 2 || ijob == 4) {
		rcondv[1] = dif[0];
		rcondv[2] = dif[1];
	    }
	    if (ierr == 1) {
		*info = n + 3;
	    }
	}
    }
//Apply permutation to VSL and VSR
//(Workspace: none needed)
    if (ilvsl) {
	Rggbak("P", "L", n, ilo, ihi, &work[ileft], &work[iright], n, &vsl[0], ldvsl, &ierr);
    }
    if (ilvsr) {
	Rggbak("P", "R", n, ilo, ihi, &work[ileft], &work[iright], n, &vsr[0], ldvsr, &ierr);
    }
//Check if unscaling would cause over/underflow, if so, rescale
//(ALPHAR(I),ALPHAI(I),BETA(I)) so BETA(I) is on the order of
//B(I,I) and ALPHAR(I) and ALPHAI(I) are on the order of A(I,I)
    if (ilascl) {
	for (i = 0; i < n; i++) {
	    if (alphai[i] != Zero) {
		if (alphar[i] / safmax > anrmto / anrm || safmin / alphar[i] > anrm / anrmto) {
		    work[1] = abs(A[i + i * lda] / alphar[i]);
		    beta[i] = beta[i] * work[1];
		    alphar[i] = alphar[i] * work[1];
		    alphai[i] = alphai[i] * work[1];
		} else if (alphai[i] / safmax > anrmto / anrm || safmin / alphai[i] > anrm / anrmto) {
		    work[1] = abs(A[i + (i + 1) * lda] / alphai[i]);
		    beta[i] = beta[i] * work[1];
		    alphar[i] = alphar[i] * work[1];
		    alphai[i] = alphai[i] * work[1];
		}
	    }
	}
    }
    if (ilbscl) {
	for (i = 0; i < n; i++) {
	    if (alphai[i] != Zero) {
		if (beta[i] / safmax > bnrmto / bnrm || safmin / beta[i] > bnrm / bnrmto) {
		    work[1] = abs(B[i + i * ldb] / beta[i]);
		    beta[i] = beta[i] * work[1];
		    alphar[i] = alphar[i] * work[1];
		    alphai[i] = alphai[i] * work[1];
		}
	    }
	}
    }
//Undo scaling
    if (ilascl) {
	Rlascl("H", 0, 0, anrmto, anrm, n, n, &A[0], lda, &ierr);
	Rlascl("G", 0, 0, anrmto, anrm, n, 1, &alphar[1], n, &ierr);
	Rlascl("G", 0, 0, anrmto, anrm, n, 1, &alphai[1], n, &ierr);
    }
    if (ilbscl) {
	Rlascl("U", 0, 0, bnrmto, bnrm, n, n, &B[0], ldb, &ierr);
	Rlascl("G", 0, 0, bnrmto, bnrm, n, 1, &beta[1], n, &ierr);
    }
    if (wantst) {
//Check if reordering is correct
	lastsl = MTRUE;
	lst2sl = MTRUE;
	*sdim = 0;
	ip = 0;
	for (i = 0; i < n; i++) {
	    cursl = (*selctg) (&alphar[i], &alphai[i], &beta[i]);
	    if (alphai[i] == Zero) {
		if (cursl) {
		    ++(*sdim);
		}
		ip = 0;
		if (cursl && !lastsl) {
		    *info = n + 2;
		}
	    } else {
		if (ip == 1) {
//Last eigenvalue of conjugate pair
		    cursl = cursl || lastsl;
		    lastsl = cursl;
		    if (cursl) {
			*sdim += 2;
		    }
		    ip = -1;
		    if (cursl && !lst2sl) {
			*info = n + 2;
		    }
		} else {
//First eigenvalue of conjugate pair
		    ip = 1;
		}
	    }
	    lst2sl = lastsl;
	    lastsl = cursl;
	}
    }
  L60:
    work[1] = maxwrk;
    iwork[1] = liwmin;
    return;
}
