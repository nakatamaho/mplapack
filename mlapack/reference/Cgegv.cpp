/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgegv.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cgegv(const char *jobvl, const char *jobvr, INTEGER n,
	   COMPLEX * A, INTEGER lda, COMPLEX * B, INTEGER ldb,
	   COMPLEX * alpha, COMPLEX * beta, COMPLEX * vl, INTEGER ldvl, COMPLEX * vr, INTEGER ldvr, COMPLEX * work, INTEGER lwork, REAL * rwork, INTEGER * info)
{
    INTEGER jc, nb, in, jr, nb1, nb2, nb3, ihi, ilo;
    REAL eps;
    LOGICAL ilv;
    REAL absb, anrm, bnrm;
    INTEGER itau;
    REAL temp;
    LOGICAL ilvl, ilvr;
    INTEGER lopt;
    REAL anrm1, anrm2, bnrm1, bnrm2, absai, scale, absar, sbeta;
    INTEGER ileft, iinfo, icols, iwork, irows;
    REAL salfai;
    REAL salfar, safmin;
    REAL safmax;
    char chtemp;
    LOGICAL ldumma;
    INTEGER ijobvl, iright;
    LOGICAL ilimit;
    INTEGER ijobvr;
    INTEGER lwkmin;
    INTEGER irwork, lwkopt;
    LOGICAL lquery;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1, mtemp2;
    REAL mtemp3, mtemp4;

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
    } else if (ldvl < 1 || (ilvl && ldvl < n)) {
	*info = -11;
    } else if (ldvr < 1 || (ilvr && ldvr < n)) {
	*info = -13;
    } else if (lwork < lwkmin && !lquery) {
	*info = -15;
    }
    if (*info == 0) {
	nb1 = iMlaenv(1, "Cgeqrf", " ", n, n, -1, -1);
	nb2 = iMlaenv(1, "Cunmqr", " ", n, n, n, -1);
	nb3 = iMlaenv(1, "Cungqr", " ", n, n, n, -1);
	nb = max(max(nb1, nb2), nb3);
	lopt = max(n << 1, n * (nb + 1));
	work[1] = lopt;
    }
    if (*info != 0) {
	Mxerbla("Cgegv ", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
// Get machine constants
    eps = Rlamch("E") * Rlamch("B");
    safmin = Rlamch("S");
    safmin = safmin + safmin;
    safmax = One / safmin;
//Scale A
    anrm = Clange("M", n, n, &A[0], lda, &rwork[1]);
    anrm1 = anrm;
    anrm2 = One;
    if (anrm < One) {
	if (safmax * anrm < One) {
	    anrm1 = safmin;
	    anrm2 = safmax * anrm;
	}
    }
    if (anrm > Zero) {
	Clascl("G", -1, -1, anrm, One, n, n, &A[0], lda, &iinfo);
	if (iinfo != 0) {
	    *info = n + 10;
	    return;
	}
    }
//Scale B
    bnrm = Clange("M", n, n, &B[0], ldb, &rwork[1]);
    bnrm1 = bnrm;
    bnrm2 = One;
    if (bnrm < One) {
	if (safmax * bnrm < One) {
	    bnrm1 = safmin;
	    bnrm2 = safmax * bnrm;
	}
    }
    if (bnrm > Zero) {
	Clascl("G", -1, -1, bnrm, One, n, n, &B[0], ldb, &iinfo);
	if (iinfo != 0) {
	    *info = n + 10;
	    return;
	}
    }
//Permute the matrix to make it more nearly triangular
//Also "balance" the matrix.
    ileft = 1;
    iright = n + 1;
    irwork = iright + n;
    Cggbal("P", n, &A[0], lda, &B[0], ldb, &ilo, &ihi, &rwork[ileft], &rwork[iright], &rwork[irwork], &iinfo);
    if (iinfo != 0) {
	*info = n + 1;
	goto L80;
    }
//Reduce B to triangular form, and initialize VL and/or VR
    irows = ihi + 1 - ilo;
    if (ilv) {
	icols = n + 1 - ilo;
    } else {
	icols = irows;
    }
    itau = 1;
    iwork = itau + irows;
    Cgeqrf(irows, icols, &B[ilo + ilo * ldb], ldb, &work[itau], &work[iwork], lwork + 1 - iwork, &iinfo);
    if (iinfo >= 0) {
	lwkopt = max(lwkopt, (INTEGER) cast2double(work[iwork].real()) + iwork - 1);
    }
    if (iinfo != 0) {
	*info = n + 2;
	goto L80;
    }
    Cunmqr("L", "C", irows, icols, irows, &B[ilo + ilo * ldb], ldb, &work[itau], &A[ilo + ilo * lda], lda, &work[iwork], lwork + 1 - iwork, &iinfo);
    if (iinfo >= 0) {
	lwkopt = max(lwkopt, (INTEGER) cast2double(work[iwork].real()) + iwork - 1);
    }
    if (iinfo != 0) {
	*info = n + 3;
	goto L80;
    }
    if (ilvl) {
	Claset("Full", n, n, Zero, One, &vl[0], ldvl);
	Clacpy("L", irows - 1, irows - 1, &B[ilo + 1 + ilo * ldb], ldb, &vl[ilo + 1 + ilo * ldvl], ldvl);
	Cungqr(irows, irows, irows, &vl[ilo + ilo * ldvl], ldvl, &work[itau], &work[iwork], lwork + 1 - iwork, &iinfo);
	if (iinfo >= 0) {
	    lwkopt = max(lwkopt, (INTEGER) cast2double(work[iwork].real()) + iwork - 1);
	}
	if (iinfo != 0) {
	    *info = n + 4;
	    goto L80;
	}
    }
    if (ilvr) {
	Claset("Full", n, n, Zero, One, &vr[0], ldvr);
    }
//Reduce to generalized Hessenberg form
    if (ilv) {
//Eigenvectors requested -- work on whole matrix.
	Cgghrd(jobvl, jobvr, n, ilo, ihi, &A[0], lda, &B[0], ldb, &vl[0], ldvl, &vr[0], ldvr, &iinfo);
    } else {
	Cgghrd("N", "N", irows, 1, irows, &A[ilo + ilo * lda], lda, &B[ilo + ilo * ldb], ldb, &vl[0], ldvl, &vr[0], ldvr, &iinfo);
    }
    if (iinfo != 0) {
	*info = n + 5;
	goto L80;
    }
//Perform QZ algorithm
    iwork = itau;
    if (ilv) {
	chtemp = 'S';
    } else {
	chtemp = 'E';
    }
    Chgeqz((const char *) chtemp, jobvl, jobvr, n, ilo, ihi, &A[0], lda, &B[0], ldb, &alpha[1], &beta[1], &vl[0], ldvl, &vr[0],
	   ldvr, &work[iwork], lwork + 1 - iwork, &rwork[irwork], &iinfo);
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
	goto L80;
    }
    if (ilv) {
//Compute Eigenvectors
	if (ilvl) {
	    if (ilvr) {
		chtemp = 'B';
	    } else {
		chtemp = 'L';
	    }
	} else {
	    chtemp = 'R';
	}
	Ctgevc((const char *) chtemp, "B", &ldumma, n, &A[0], lda, &B[0], ldb, &vl[0], ldvl, &vr[0], ldvr, n, &in, &work[iwork], &rwork[irwork], &iinfo);
	if (iinfo != 0) {
	    *info = n + 7;
	    goto L80;
	}
//Undo balancing on VL and VR, rescale
	if (ilvl) {
	    Cggbak("P", "L", n, ilo, ihi, &rwork[ileft], &rwork[iright], n, &vl[0], ldvl, &iinfo);
	    if (iinfo != 0) {
		*info = n + 8;
		goto L80;
	    }
	    for (jc = 1; jc <= n; jc++) {
		temp = Zero;
		for (jr = 1; jr <= n; jr++) {
		    mtemp1 = temp, mtemp2 = Cabs1(vl[jr + jc * ldvl]);
		    temp = max(mtemp1, mtemp2);
		}
		if (temp < safmin) {
		    goto L30;
		}
		temp = One / temp;
		for (jr = 1; jr <= n; jr++) {
		    vl[jr + jc * ldvl] = vl[jr + jc * ldvl] * temp;
		}
	      L30:
		;
	    }
	}
	if (ilvr) {
	    Cggbak("P", "R", n, ilo, ihi, &rwork[ileft], &rwork[iright], n, &vr[0], ldvr, &iinfo);
	    if (iinfo != 0) {
		*info = n + 9;
		goto L80;
	    }
	    for (jc = 1; jc <= n; jc++) {
		temp = Zero;
		for (jr = 1; jr <= n; jr++) {
		    mtemp1 = temp, mtemp2 = Cabs1(vr[jr + jc * ldvr]);
		    temp = max(mtemp1, mtemp2);
		}
		if (temp < safmin) {
		    goto L60;
		}
		temp = One / temp;
		for (jr = 1; jr <= n; jr++) {
		    vr[jr + jc * ldvr] = vr[jr + jc * ldvr] * temp;
		}
	      L60:
		;
	    }
	}
//End of eigenvector calculation
    }
//Undo scaling in alpha, beta
//Note: this does not give the alpha and beta for the unscaled
//problem.
//Un-scaling is limited to avoid underflow in alpha and beta
//if they are significant.
    for (jc = 1; jc <= n; jc++) {
	absar = abs(alpha[jc].real());
	absai = abs(alpha[jc].imag());
	absb = abs(beta[jc].real());
	salfar = anrm * alpha[jc].real();
	salfai = anrm * alpha[jc].imag();
	sbeta = bnrm * beta[jc].real();
	ilimit = MFALSE;
	scale = One;
//Check for significant underflow in imaginary part of ALPHA
	mtemp1 = safmin, mtemp2 = eps * absar;
	mtemp3 = max(mtemp1, mtemp2), mtemp4 = eps * absb;
	if (abs(salfai) < safmin && absai >= max(mtemp3, mtemp4)) {
	    ilimit = MTRUE;
	    mtemp1 = safmin, mtemp2 = anrm2 * absai;
	    scale = safmin / anrm1 / max(mtemp1, mtemp2);
	}
//Check for significant underflow in real part of ALPHA
	mtemp1 = safmin, mtemp2 = eps * absai;
	mtemp3 = max(mtemp1, mtemp2), mtemp4 = eps * absb;
	if (abs(salfar) < safmin && absar >= max(mtemp3, mtemp4)) {
	    ilimit = MTRUE;
	    mtemp1 = safmin, mtemp2 = anrm2 * absar;
	    mtemp3 = scale, mtemp4 = safmin / anrm1 / max(mtemp1, mtemp2);
	    scale = max(mtemp3, mtemp4);
	}
//Check for significant underflow in BETA
	mtemp1 = safmin, mtemp2 = eps * absar, mtemp3 = max(mtemp1, mtemp2), mtemp4 = eps * absai;
	if (abs(sbeta) < safmin && absb >= max(mtemp3, mtemp4)) {
	    ilimit = MTRUE;
	    mtemp1 = safmin, mtemp2 = bnrm2 * absb;
	    mtemp3 = scale, mtemp4 = safmin / bnrm1 / max(mtemp1, mtemp2);
	    scale = max(mtemp3, mtemp4);
	}
//Check for possible overflow when limiting scaling
	if (ilimit) {
	    mtemp1 = abs(salfar), mtemp2 = abs(salfai);
	    mtemp3 = max(mtemp1, mtemp2), mtemp4 = abs(sbeta);
	    temp = scale * safmin * max(mtemp3, mtemp4);
	    if (temp > One) {
		scale /= temp;
	    }
	    if (scale < One) {
		ilimit = MFALSE;
	    }
	}
//Recompute un-scaled ALPHA, BETA if necessary.
	if (ilimit) {
	    salfar = scale * alpha[jc].real() * anrm;
	    salfai = scale * alpha[jc].imag() * anrm;
	    sbeta = scale * (beta[jc] * bnrm).real();
	}
	alpha[jc] = salfar;
	beta[jc] = sbeta;
    }
  L80:
    work[1] = lwkopt;
    return;
}
