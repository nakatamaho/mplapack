/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgegs.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cgegs(const char *jobvsl, const char *jobvsr, INTEGER n,
	   COMPLEX * A, INTEGER lda, COMPLEX * B, INTEGER ldb,
	   COMPLEX * alpha, COMPLEX * beta, COMPLEX * vsl, INTEGER ldvsl, COMPLEX * vsr, INTEGER ldvsr, COMPLEX * work, INTEGER lwork, REAL * rwork, INTEGER * info)
{
    INTEGER nb, nb1, nb2, nb3, ihi, ilo;
    REAL eps, anrm, bnrm;
    INTEGER itau, lopt;
    INTEGER ileft, iinfo, icols;
    LOGICAL ilvsl;
    INTEGER iwork;
    LOGICAL ilvsr;
    INTEGER irows;
    LOGICAL ilascl, ilbscl;
    REAL safmin;
    REAL bignum;
    INTEGER ijobvl, iright;
    INTEGER ijobvr;
    REAL anrmto = 0.0;
    INTEGER lwkmin;
    REAL bnrmto = 0.0;
    REAL smlnum;
    INTEGER irwork, lwkopt;
    LOGICAL lquery;
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
    lwkmin = max(n << 1, (INTEGER) 1);
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
	*info = -11;
    } else if (ldvsr < 1 || (ilvsr && ldvsr < n)) {
	*info = -13;
    } else if (lwork < lwkmin && !lquery) {
	*info = -15;
    }
    if (*info == 0) {
	nb1 = iMlaenv(1, "Cgeqrf", " ", n, n, -1, -1);
	nb2 = iMlaenv(1, "Cunmqr", " ", n, n, n, -1);
	nb3 = iMlaenv(1, "Cungq", " ", n, n, n, -1);
	nb = max(max(nb1, nb2), nb3);
	lopt = n * (nb + 1);
	work[1] = lopt;
    }
    if (*info != 0) {
	Mxerbla("Cgegs ", -(*info));
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
	Clascl("G", -1, -1, anrm, anrmto, n, n, &A[0], lda, &iinfo);
	if (iinfo != 0) {
	    *info = n + 9;
	    return;
	}
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
	Clascl("G", -1, -1, bnrm, bnrmto, n, n, &B[0], ldb, &iinfo);
	if (iinfo != 0) {
	    *info = n + 9;
	    return;
	}
    }
//Permute the matrix to make it more nearly triangular
    ileft = 1;
    iright = n + 1;
    irwork = iright + n;
    iwork = 0;
    Cggbal("P", n, &A[0], lda, &B[0], ldb, &ilo, &ihi, &rwork[ileft], &rwork[iright], &rwork[irwork], &iinfo);
    if (iinfo != 0) {
	*info = n + 1;
	goto L10;
    }
//Reduce B to triangular form, and initialize VSL and/or VSR
    irows = ihi + 1 - ilo;
    icols = n + 1 - ilo;
    itau = iwork;
    iwork = itau + irows;
    Cgeqrf(irows, icols, &B[ilo + ilo * ldb], ldb, &work[itau], &work[iwork], lwork + 1 - iwork, &iinfo);
    if (iinfo >= 0) {
	lwkopt = max(lwkopt, (INTEGER) cast2double(work[iwork].real()) + iwork - 1);
    }
    if (iinfo != 0) {
	*info = n + 2;
	goto L10;
    }
    Cunmqr("L", "C", irows, icols, irows, &B[ilo + ilo * ldb], ldb, &work[itau], &A[ilo + ilo * lda], lda, &work[iwork], lwork + 1 - iwork, &iinfo);
    if (iinfo >= 0) {
	lwkopt = max(lwkopt, (INTEGER) cast2double(work[iwork].real()) + iwork - 1);
    }
    if (iinfo != 0) {
	*info = n + 3;
	goto L10;
    }
    if (ilvsl) {
	Claset("Full", n, n, Zero, One, &vsl[0], ldvsl);
	Clacpy("L", irows - 1, irows - 1, &B[ilo + 1 + ilo * ldb], ldb, &vsl[ilo + 1 + ilo * ldvsl], ldvsl);
	Cungqr(irows, irows, irows, &vsl[ilo + ilo * ldvsl], ldvsl, &work[itau], &work[iwork], lwork + 1 - iwork, &iinfo);
	if (iinfo >= 0) {
	    lwkopt = max(lwkopt, (INTEGER) cast2double(work[iwork].real()) + iwork - 1);
	}
	if (iinfo != 0) {
	    *info = n + 4;
	    goto L10;
	}
    }
    if (ilvsr) {
	Claset("Full", n, n, Zero, One, &vsr[0], ldvsr);
    }
//Reduce to generalized Hessenberg form
    Cgghrd(jobvsl, jobvsr, n, ilo, ihi, &A[0], lda, &B[0], ldb, &vsl[0], ldvsl, &vsr[0], ldvsr, &iinfo);
    if (iinfo != 0) {
	*info = n + 5;
	goto L10;
    }
//Perform QZ algorithm, computing Schur vectors if desired
    iwork = itau;
    Chgeqz("S", jobvsl, jobvsr, n, ilo, ihi, &A[0], lda, &B[0], ldb, &alpha[1], &beta[1], &vsl[0], ldvsl, &vsr[0], ldvsr, &work[iwork], lwork + 1 - iwork, &rwork[irwork], &iinfo);
    if (iinfo >= 0) {
	lwkopt = max(lwkopt, (INTEGER) cast2double(work[iwork].real()) + iwork - 1);
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
	Cggbak("P", "L", n, ilo, ihi, &rwork[ileft], &rwork[iright], n, &vsl[0], ldvsl, &iinfo);
	if (iinfo != 0) {
	    *info = n + 7;
	    goto L10;
	}
    }
    if (ilvsr) {
	Cggbak("P", "R", n, ilo, ihi, &rwork[ileft], &rwork[iright], n, &vsr[0], ldvsr, &iinfo);
	if (iinfo != 0) {
	    *info = n + 8;
	    goto L10;
	}
    }
//Undo scaling
    if (ilascl) {
	Clascl("U", -1, -1, anrmto, anrm, n, n, &A[0], lda, &iinfo);
	if (iinfo != 0) {
	    *info = n + 9;
	    return;
	}
	Clascl("G", -1, -1, anrmto, anrm, n, 1, &alpha[1], n, &iinfo);
	if (iinfo != 0) {
	    *info = n + 9;
	    return;
	}
    }
    if (ilbscl) {
	Clascl("U", -1, -1, bnrmto, bnrm, n, n, &B[0], ldb, &iinfo);
	if (iinfo != 0) {
	    *info = n + 9;
	    return;
	}
	Clascl("G", -1, -1, bnrmto, bnrm, n, 1, &beta[1], n, &iinfo);
	if (iinfo != 0) {
	    *info = n + 9;
	    return;
	}
    }
  L10:
    work[1] = lwkopt;
    return;
}
