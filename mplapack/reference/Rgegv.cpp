/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgegv.cpp,v 1.3 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rgegv(const char *jobvl, const char *jobvr, INTEGER n, REAL *
	   A, INTEGER lda, REAL * B, INTEGER ldb, REAL * alphar,
	   REAL * alphai, REAL * beta, REAL * vl, INTEGER ldvl, REAL * vr, INTEGER ldvr, REAL * work, INTEGER lwork, INTEGER * info)
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
    REAL salfar;
    REAL safmin;
    REAL safmax;
    char chtemp;
    LOGICAL ldumma;
    INTEGER ijobvl, iright;
    LOGICAL ilimit;
    INTEGER ijobvr;
    REAL onepls;
    INTEGER lwkmin;
    INTEGER lwkopt;
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
    lwkmin = max(n << 3, (INTEGER) 1);
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
	*info = -12;
    } else if (ldvr < 1 || (ilvr && ldvr < n)) {
	*info = -14;
    } else if (lwork < lwkmin && !lquery) {
	*info = -16;
    }
    if (*info == 0) {
	nb1 = iMlaenv(1, "Rgeqrf", " ", n, n, -1, -1);
	nb2 = iMlaenv(1, "Rormqr", " ", n, n, n, -1);
	nb3 = iMlaenv(1, "Rorgqr", " ", n, n, n, -1);
	nb = max(max(nb1, nb2), nb3);
	lopt = (n << 1) + max(n * 6, n * (nb + 1));
	work[1] = lopt;
    }
    if (*info != 0) {
	Mxerbla("Rgegv ", -(*info));
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
    safmin = safmin + safmin;
    safmax = One / safmin;
    onepls = eps * 4 + One;
//Scale A
    anrm = Rlange("M", n, n, &A[0], lda, &work[0]);
    anrm1 = anrm;
    anrm2 = One;
    if (anrm < One) {
	if (safmax * anrm < One) {
	    anrm1 = safmin;
	    anrm2 = safmax * anrm;
	}
    }
    if (anrm > Zero) {
	Rlascl("G", -1, -1, anrm, One, n, n, &A[0], lda, &iinfo);
	if (iinfo != 0) {
	    *info = n + 10;
	    return;
	}
    }
//Scale B
    bnrm = Rlange("M", n, n, &B[0], ldb, &work[0]);
    bnrm1 = bnrm;
    bnrm2 = One;
    if (bnrm < One) {
	if (safmax * bnrm < One) {
	    bnrm1 = safmin;
	    bnrm2 = safmax * bnrm;
	}
    }
    if (bnrm > Zero) {
	Rlascl("G", -1, -1, bnrm, One, n, n, &B[0], ldb, &iinfo);
	if (iinfo != 0) {
	    *info = n + 10;
	    return;
	}
    }
//Permute the matrix to make it more nearly triangular
//Workspace layout:  (8*N words -- "work" requires 6*N words)
//   left_permutation, right_permutation, work...
    ileft = 1;
    iright = n + 1;
    iwork = iright + n;
    Rggbal("P", n, &A[0], lda, &B[0], ldb, &ilo, &ihi, &work[ileft], &work[iright], &work[iwork], &iinfo);
    if (iinfo != 0) {
	*info = n + 1;
	goto L120;
    }
//Reduce B to triangular form, and initialize VL and/or VR
//Workspace layout:  ("work..." must have at least N words)
//   left_permutation, right_permutation, tau, work...
    irows = ihi + 1 - ilo;
    if (ilv) {
	icols = n + 1 - ilo;
    } else {
	icols = irows;
    }
    itau = iwork;
    iwork = itau + irows;
    Rgeqrf(irows, icols, &B[ilo + ilo * ldb], ldb, &work[itau], &work[iwork], lwork + 1 - iwork, &iinfo);
    if (iinfo >= 0) {
	lwkopt = max(lwkopt, (INTEGER) cast2double(work[iwork]) + iwork - 1);
    }
    if (iinfo != 0) {
	*info = n + 2;
	goto L120;
    }
    Rormqr("L", "T", irows, icols, irows, &B[ilo + ilo * ldb], ldb, &work[itau], &A[ilo + ilo * lda], lda, &work[iwork], lwork + 1 - iwork, &iinfo);
    if (iinfo >= 0) {
	lwkopt = max(lwkopt, (INTEGER) cast2double(work[iwork]) + iwork - 1);
    }
    if (iinfo != 0) {
	*info = n + 3;
	goto L120;
    }
    if (ilvl) {
	Rlaset("Full", n, n, Zero, One, &vl[0], ldvl);
	Rlacpy("L", irows - 1, irows - 1, &B[ilo + 1 + ilo * ldb], ldb, &vl[ilo + 1 + ilo * ldvl], ldvl);
	Rorgqr(irows, irows, irows, &vl[ilo + ilo * ldvl], ldvl, &work[itau], &work[iwork], lwork + 1 - iwork, &iinfo);
	if (iinfo >= 0) {
	    lwkopt = max(lwkopt, (INTEGER) cast2double(work[iwork]) + iwork - 1);
	}
	if (iinfo != 0) {
	    *info = n + 4;
	    goto L120;
	}
    }
    if (ilvr) {
	Rlaset("Full", n, n, Zero, One, &vr[0], ldvr);
    }
//Reduce to generalized Hessenberg form
    if (ilv) {
//Eigenvectors requested -- work on whole matrix.
	Rgghrd(jobvl, jobvr, n, ilo, ihi, &A[0], lda, &B[0], ldb, &vl[0], ldvl, &vr[0], ldvr, &iinfo);
    } else {
	Rgghrd("N", "N", irows, 1, irows, &A[ilo + ilo * lda], lda, &B[ilo + ilo * ldb], ldb, &vl[0], ldvl, &vr[0], ldvr, &iinfo);
    }
    if (iinfo != 0) {
	*info = n + 5;
	goto L120;
    }
//Perform QZ algorithm
//Workspace layout:  ("work..." must have at least 1 word)
//   left_permutation, right_permutation, work...
    iwork = itau;
    if (ilv) {
	chtemp = 'S';
    } else {
	chtemp = 'E';
    }
    Rhgeqz((const char *) chtemp, jobvl, jobvr, n, ilo, ihi, &A[0], lda, &B[0], ldb, &alphar[1], &alphai[1], &beta[1], &vl[0], ldvl,
	   &vr[0], ldvr, &work[iwork], lwork + 1 - iwork, &iinfo);
    if (iinfo >= 0) {
	lwkopt = max(lwkopt, (INTEGER) cast2double(work[iwork]) + iwork - 1);
    }
    if (iinfo != 0) {
	if (iinfo > 0 && iinfo <= n) {
	    *info = iinfo;
	} else if (iinfo > n && iinfo <= n << 1) {
	    *info = iinfo - n;
	} else {
	    *info = n + 6;
	}
	goto L120;
    }
    if (ilv) {
//Compute Eigenvectors  (DTGEVC requires 6*N words of workspace)
	if (ilvl) {
	    if (ilvr) {
		chtemp = 'B';
	    } else {
		chtemp = 'L';
	    }
	} else {
	    chtemp = 'R';
	}
	Rtgevc((const char *) chtemp, "B", &ldumma, n, &A[0], lda, &B[0], ldb, &vl[0], ldvl, &vr[0], ldvr, n, &in, &work[iwork], &iinfo);
	if (iinfo != 0) {
	    *info = n + 7;
	    goto L120;
	}
//Undo balancing on VL and VR, rescale
	if (ilvl) {
	    Rggbak("P", "L", n, ilo, ihi, &work[ileft], &work[iright], n, &vl[0], ldvl, &iinfo);
	    if (iinfo != 0) {
		*info = n + 8;
		goto L120;
	    }
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
		if (temp < safmin) {
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
	    Rggbak("P", "R", n, ilo, ihi, &work[ileft], &work[iright], n, &vr[0], ldvr, &iinfo);
	    if (iinfo != 0) {
		*info = n + 9;
		goto L120;
	    }
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
		if (temp < safmin) {
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
/*     Undo scaling in alpha, beta */
/*     Note: this does not give the alpha and beta for the unscaled */
/*     problem. */
/*     Un-scaling is limited to avoid underflow in alpha and beta */
/*     if they are significant. */
    for (jc = 1; jc <= n; jc++) {
	absar = abs(alphar[jc]);
	absai = abs(alphai[jc]);
	absb = abs(beta[jc]);
	salfar = anrm * alphar[jc];
	salfai = anrm * alphai[jc];
	sbeta = bnrm * beta[jc];
	ilimit = MFALSE;
	scale = One;
//Check for significant underflow in ALPHAI
	mtemp1 = safmin, mtemp2 = eps * absar;
	mtemp3 = max(mtemp1, mtemp2), mtemp4 = eps * absb;
	if (abs(salfai) < safmin && absai >= max(mtemp3, mtemp4)) {
	    ilimit = MTRUE;
	    mtemp1 = onepls * safmin, mtemp2 = anrm2 * absai;
	    scale = onepls * safmin / anrm1 / max(mtemp1, mtemp2);
	} else if (salfai == Zero) {
//If insignificant underflow in ALPHAI, then make the
//conjugate eigenvalue real.
	    if (alphai[jc] < Zero && jc > 1) {
		alphai[jc - 1] = Zero;
	    } else if (alphai[jc] > Zero && jc < n) {
		alphai[jc + 1] = Zero;
	    }
	}
//Check for significant underflow in ALPHAR
	mtemp1 = safmin, mtemp2 = eps * absai;
	mtemp3 = max(mtemp1, mtemp2), mtemp4 = eps * absb;
	if (abs(salfar) < safmin && absar >= max(mtemp3, mtemp4)) {
	    ilimit = MTRUE;
	    mtemp1 = onepls * safmin, mtemp2 = anrm2 * absar;
	    mtemp3 = scale, mtemp4 = onepls * safmin / anrm1 / max(mtemp3, mtemp4);
	    scale = max(mtemp3, mtemp4);
	}
	mtemp1 = safmin, mtemp2 = eps * absar;
	mtemp3 = max(mtemp1, mtemp2), mtemp4 = eps * absai;
	if (abs(sbeta) < safmin && absb >= max(mtemp3, mtemp4)) {
	    ilimit = MTRUE;
	    mtemp1 = onepls * safmin, mtemp2 = bnrm2 * absb;
	    mtemp3 = scale, mtemp4 = onepls * safmin / bnrm1 / max(mtemp1, mtemp2);
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
//Recompute un-scaled ALPHAR, ALPHAI, BETA if necessary.
	if (ilimit) {
	    salfar = scale * alphar[jc] * anrm;
	    salfai = scale * alphai[jc] * anrm;
	    sbeta = scale * beta[jc] * bnrm;
	}
	alphar[jc] = salfar;
	alphai[jc] = salfai;
	beta[jc] = sbeta;
    }
  L120:
    work[1] = lwkopt;
    return;
}
