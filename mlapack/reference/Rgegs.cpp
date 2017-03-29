/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgegs.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rgegs(const char *jobvsl, const char *jobvsr, INTEGER n,
	   REAL * A, INTEGER lda, REAL * B, INTEGER ldb, REAL *
	   alphar, REAL * alphai, REAL * beta, REAL * vsl, INTEGER ldvsl, REAL * vsr, INTEGER ldvsr, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER nb1, nb2, nb3, ihi, ilo;
    REAL eps, anrm, bnrm;
    INTEGER itau, lopt;
    INTEGER ileft, iright, iinfo, icols;
    INTEGER ilvsl;
    INTEGER iwork;
    INTEGER ilvsr;
    INTEGER irows;
    INTEGER ilascl, ilbscl;
    REAL safmin;
    REAL bignum;
    REAL anrmto = 0.0;
    INTEGER lwkmin;
    REAL bnrmto = 0.0;
    REAL smlnum;
    INTEGER lwkopt, ijobvr, ijobvl;
    INTEGER lquery;
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
//Test the input arguments
    lwkmin = max(n << 2, (INTEGER) 1);
    lwkopt = lwkmin;
    work[1] = lwkopt;
    lquery = lwork == -1;
    *info = 0;
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
    } else if (ldvsl < 1 || (ilvsl && ldvsl < n)) {
	*info = -12;
    } else if (ldvsr < 1 || (ilvsr && ldvsr < n)) {
	*info = -14;
    } else if (lwork < lwkmin && !lquery) {
	*info = -16;
    }
    if (*info == 0) {
	nb1 = iMlaenv(1, "Rgeqrf", " ", n, n, -1, -1);
	nb2 = iMlaenv(1, "Rormqr", " ", n, n, n, -1);
	nb3 = iMlaenv(1, "Rorgqr", " ", n, n, n, -1);
	lopt = (n << 1) + n * (max(max(nb1, nb2), nb3) + 1);
	work[1] = lopt;
    }
    if (*info != 0) {
	Mxerbla("Rgegs ", -(*info));
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
    safmin = Rlamch("S");
    smlnum = n * safmin / eps;
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
	Rlascl("G", -1, -1, anrm, anrmto, n, n, &A[0], lda, &iinfo);
	if (iinfo != 0) {
	    *info = n + 9;
	    return;
	}
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
	Rlascl("G", -1, -1, bnrm, bnrmto, n, n, &B[0], ldb, &iinfo);
	if (iinfo != 0) {
	    *info = n + 9;
	    return;
	}
    }
//Permute the matrix to make it more nearly triangular
//Workspace layout:  (2*N words -- "work..." not actually used)
//   left_permutation, right_permutation, work...
    ileft = 1;
    iright = n + 1;
    iwork = iright + n;
    Rggbal("P", n, &A[0], lda, &B[0], ldb, &ilo, &ihi, &work[ileft], &work[iright], &work[iwork], &iinfo);
    if (iinfo != 0) {
	*info = n + 1;
	goto L10;
    }
//Reduce B to triangular form, and initialize VSL and/or VSR
//Workspace layout:  ("work..." must have at least N words)
//   left_permutation, right_permutation, tau, work...
    irows = ihi + 1 - ilo;
    icols = n + 1 - ilo;
    itau = iwork;
    iwork = itau + irows;
    Rgeqrf(irows, icols, &B[ilo + ilo * ldb], ldb, &work[itau], &work[iwork], lwork + 1 - iwork, &iinfo);
    if (iinfo >= 0) {
	lwkopt = max(lwkopt, (INTEGER) cast2double(work[iwork]) + iwork - 1);
    }
    if (iinfo != 0) {
	*info = n + 2;
	goto L10;
    }
    Rormqr("L", "T", irows, icols, irows, &B[ilo + ilo * ldb], ldb, &work[itau], &A[ilo + ilo * lda], lda, &work[iwork], lwork + 1 - iwork, &iinfo);
    if (iinfo >= 0) {
	lwkopt = max(lwkopt, (INTEGER) cast2double(work[iwork]) + iwork - 1);
    }
    if (iinfo != 0) {
	*info = n + 3;
	goto L10;
    }
    if (ilvsl) {
	Rlaset("Full", n, n, Zero, One, &vsl[0], ldvsl);
	Rlacpy("L", irows - 1, irows - 1, &B[ilo + 1 + ilo * ldb], ldb, &vsl[ilo + 1 + ilo * ldvsl], ldvsl);
	Rorgqr(irows, irows, irows, &vsl[ilo + ilo * ldvsl], ldvsl, &work[itau], &work[iwork], lwork + 1 - iwork, &iinfo);
	if (iinfo >= 0) {
	    lwkopt = max(lwkopt, (INTEGER) cast2double(work[iwork]) + iwork - 1);
	}
	if (iinfo != 0) {
	    *info = n + 4;
	    goto L10;
	}
    }
    if (ilvsr) {
	Rlaset("Full", n, n, Zero, One, &vsr[0], ldvsr);
    }
//Reduce to generalized Hessenberg form
    Rgghrd(jobvsl, jobvsr, n, ilo, ihi, &A[0], lda, &B[0], ldb, &vsl[0], ldvsl, &vsr[0], ldvsr, &iinfo);
    if (iinfo != 0) {
	*info = n + 5;
	goto L10;
    }
//Perform QZ algorithm, computing Schur vectors if desired
//Workspace layout:  ("work..." must have at least 1 word)
//  left_permutation, right_permutation, work...
    iwork = itau;
    Rhgeqz("S", jobvsl, jobvsr, n, ilo, ihi, &A[0], lda, &B[0], ldb, &alphar[1], &alphai[1], &beta[1], &vsl[0]
	   , ldvsl, &vsr[0], ldvsr, &work[iwork], lwork + 1 - iwork, &iinfo);
    if (iinfo >= 0) {
	lwkopt = max(lwkopt, (INTEGER) cast2double(work[iwork] + iwork - 1));
    }
    if (iinfo != 0) {
	if (iinfo > 0 && iinfo <= n) {
	    *info = iinfo;
	} else if (iinfo > n && iinfo <= n << 1) {
	    *info = iinfo - n;
	} else {
	    *info = n + 6;
	}
	goto L10;
    }
//Apply permutation to VSL and VSR
    if (ilvsl) {
	Rggbak("P", "L", n, ilo, ihi, &work[ileft], &work[iright], n, &vsl[0], ldvsl, &iinfo);
	if (iinfo != 0) {
	    *info = n + 7;
	    goto L10;
	}
    }
    if (ilvsr) {
	Rggbak("P", "R", n, ilo, ihi, &work[ileft], &work[iright], n, &vsr[0], ldvsr, &iinfo);
	if (iinfo != 0) {
	    *info = n + 8;
	    goto L10;
	}
    }
//Undo scaling
    if (ilascl) {
	Rlascl("H", -1, -1, anrmto, anrm, n, n, &A[0], lda, &iinfo);
	if (iinfo != 0) {
	    *info = n + 9;
	    return;
	}
	Rlascl("G", -1, -1, anrmto, anrm, n, 1, &alphar[1], n, &iinfo);
	if (iinfo != 0) {
	    *info = n + 9;
	    return;
	}
	Rlascl("G", -1, -1, anrmto, anrm, n, 1, &alphai[1], n, &iinfo);
	if (iinfo != 0) {
	    *info = n + 9;
	    return;
	}
    }
    if (ilbscl) {
	Rlascl("U", -1, -1, bnrmto, bnrm, n, n, &B[0], ldb, &iinfo);
	if (iinfo != 0) {
	    *info = n + 9;
	    return;
	}
	Rlascl("G", -1, -1, bnrmto, bnrm, n, 1, &beta[1], n, &iinfo);
	if (iinfo != 0) {
	    *info = n + 9;
	    return;
	}
    }
  L10:
    work[1] = lwkopt;
    return;
}
