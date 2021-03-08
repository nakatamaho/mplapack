/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cggbal.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cggbal(const char *job, INTEGER n, COMPLEX * A, INTEGER lda, COMPLEX * B, INTEGER ldb, INTEGER * ilo, INTEGER * ihi, REAL * lscale, REAL * rscale, REAL * work, INTEGER * info)
{
    INTEGER i, j, k, l, m;
    REAL t;
    INTEGER jc;
    REAL ta, tb, tc;
    INTEGER ir;
    REAL ew;
    INTEGER it, nr, ip1, jp1, lm1;
    REAL cab, rab, ewc, cor, sum;
    INTEGER nrp2, icab, lcab;
    REAL beta, coef;
    INTEGER irab, lrab;
    REAL basl, cmax;
    REAL coef2, coef5, gamma, alpha;
    REAL sfmin, sfmax;
    INTEGER iflow;
    INTEGER kount;
    REAL pgamma = 0.0;
    INTEGER lsfmin;
    INTEGER lsfmax;
    REAL Zero = 0.0, Half = .5, One = 1.0, Ten = 10.0;
    REAL mtemp1, mtemp2;

//Test the input parameters
    *info = 0;
    if (!Mlsame(job, "N") && !Mlsame(job, "P") && !Mlsame(job, "S")
	&& !Mlsame(job, "B")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -4;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -6;
    }
    if (*info != 0) {
	Mxerbla("Cggbal", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	*ilo = 1;
	*ihi = n;
	return;
    }
    if (n == 1) {
	*ilo = 1;
	*ihi = n;
	lscale[1] = One;
	rscale[1] = One;
	return;
    }
    if (Mlsame(job, "N")) {
	*ilo = 1;
	*ihi = n;
	for (i = 0; i < n; i++) {
	    lscale[i] = One;
	    rscale[i] = One;
	}
	return;
    }
    k = 0;
    l = n;
    if (Mlsame(job, "S")) {
	goto L190;
    }
    goto L30;
//Permute the matrices A and B to isolate the eigenvalues.
//Find row with one nonzero in columns 1 through L
  L20:
    l = lm1;
    if (l != 1) {
	goto L30;
    }
    rscale[1] = One;
    lscale[1] = One;
    goto L190;

  L30:
    lm1 = l - 1;
    for (i = l; i >= 1; i--) {
	for (j = 0; j < lm1; j++) {
	    jp1 = j + 1;
	    if (A[i + j * lda] != Zero || B[i + j * ldb] != Zero) {
		goto L50;
	    }
	}
	j = l;
	goto L70;
      L50:
	for (j = jp1; j <= l; j++) {
	    if (A[i + j * lda] != Zero || B[i + j * ldb] != Zero) {
		goto L80;
	    }
	}
	j = jp1 - 1;
      L70:
	m = l;
	iflow = 1;
	goto L160;
      L80:
	;
    }
    goto L100;
//Find column with one nonzero in rows K through N
  L90:
    k++;
  L100:
    for (j = k; j <= l; j++) {
	for (i = k; i <= lm1; i++) {
	    ip1 = i + 1;
	    if (A[i + j * lda] != Zero || B[i + j * ldb] != Zero) {
		goto L120;
	    }
	}
	i = l;
	goto L140;
      L120:
	for (i = ip1; i <= l; i++) {
	    if (A[i + j * lda] != Zero || B[i + j * ldb] != Zero) {
		goto L150;
	    }
	}
	i = ip1 - 1;
      L140:
	m = k;
	iflow = 2;
	goto L160;
      L150:
	;
    }
    goto L190;
//Permute rows M and I
  L160:
    lscale[m] = i;
    if (i == m) {
	goto L170;
    }
    Cswap(n - k + 1, &A[i + k * lda], lda, &A[m + k * lda], lda);
    Cswap(n - k + 1, &B[i + k * ldb], ldb, &B[m + k * ldb], ldb);
//Permute columns M and J
  L170:
    rscale[m] = j;
    if (j == m) {
	goto L180;
    }
    Cswap(l, &A[j * lda], 1, &A[m * lda], 1);
    Cswap(l, &B[j * ldb + 1], 1, &B[m * ldb + 1], 1);
  L180:
    switch (iflow) {
    case 1:
	goto L20;
    case 2:
	goto L90;
    }

  L190:
    *ilo = k;
    *ihi = l;
    if (Mlsame(job, "P")) {
	for (i = *ilo; i <= *ihi; i++) {
	    lscale[i] = One;
	    rscale[i] = One;
	}
	return;
    }
    if (*ilo == *ihi) {
	return;
    }
//Balance the submatrix in rows ILO to IHI.
    nr = *ihi - *ilo + 1;
    for (i = *ilo; i <= *ihi; i++) {
	rscale[i] = Zero;
	lscale[i] = Zero;
	work[i] = Zero;
	work[i + n] = Zero;
	work[i + n * 2] = Zero;
	work[i + n * 3] = Zero;
	work[i + n * 4] = Zero;
	work[i + n * 5] = Zero;
    }
//Compute right side vector in resulting linear equations
    basl = log10(Ten);
    for (i = *ilo; i <= *ihi; i++) {
	for (j = *ilo; j <= *ihi; j++) {
	    if (A[i + j * lda] == Zero) {
		ta = Zero;
		goto L210;
	    }
	    mtemp1 = Cabs1(A[i + j * lda]);
	    ta = log10(mtemp1) / basl;

	  L210:
	    if (B[i + j * ldb] == Zero) {
		tb = Zero;
		goto L220;
	    }
	    tb = log10(Cabs1(B[i + j * ldb])) / basl;
	  L220:
	    work[i + (n * 4)] = work[i + (n * 4)] - ta - tb;
	    work[j + n * 5] = work[j + n * 5] - ta - tb;
	}
    }
    coef = One / (double) (nr * 2);
    coef2 = coef * coef;
    coef5 = coef2 * Half;
    nrp2 = nr + 2;
    beta = Zero;
    it = 1;
//Start generalized conjugate gradient iteration
  L250:
    gamma = Rdot(nr, &work[*ilo + (n * 4)], 1, &work[*ilo + (n * 4)], 1) + Rdot(nr, &work[*ilo + n * 5], 1, &work[*ilo + n * 5], 1);
    ew = Zero;
    ewc = Zero;
    for (i = *ilo; i <= *ihi; i++) {
	ew = ew + work[i + (n * 4)];
	ewc = ewc + work[i + n * 5];
    }
    gamma = coef * gamma - coef2 * (ew * ew + ewc * ewc) - coef5 * ((ew - ewc) * (ew - ewc));
    if (gamma == Zero) {
	goto L350;
    }
    if (it != 1) {
	beta = gamma / pgamma;
    }
    t = coef5 * (ewc - ew * 3.);
    tc = coef5 * (ew - ewc * 3.);
    Rscal(nr, beta, &work[*ilo], 1);
    Rscal(nr, beta, &work[*ilo + n], 1);
    Raxpy(nr, coef, &work[*ilo + (n * 4)], 1, &work[*ilo + n], 1);
    Raxpy(nr, coef, &work[*ilo + n * 5], 1, &work[*ilo], 1);
    for (i = *ilo; i <= *ihi; i++) {
	work[i] = work[i] + tc;
	work[i + n] = work[i + n] + t;
    }
//Apply matrix to vector
    for (i = *ilo; i <= *ihi; i++) {
	kount = 0;
	sum = Zero;
	for (j = *ilo; j <= *ihi; j++) {
	    if (A[i + j * lda] == Zero) {
		goto L280;
	    }
	    kount++;
	    sum = sum + work[j];
	  L280:
	    if (B[i + j * ldb] == Zero) {
		goto L290;
	    }
	    kount++;
	    sum = sum + work[j];
	  L290:;
	}
	work[i + (n * 2)] = kount * work[i + n] + sum;
    }
    for (j = *ilo; j <= *ihi; j++) {
	kount = 0;
	sum = Zero;
	for (i = *ilo; i <= *ihi; i++) {
	    if (A[i + j * lda] == Zero) {
		goto L310;
	    }
	    kount++;
	    sum = sum + work[i + n];
	  L310:
	    if (B[i + j * ldb] == Zero) {
		goto L320;
	    }
	    kount++;
	    sum = sum + work[i + n];
	  L320:
	    ;
	}
	work[j + n * 3] = kount * work[j] + sum;
    }
    sum = Rdot(nr, &work[*ilo + n], 1, &work[*ilo + (n * 2)], 1)
	+ Rdot(nr, &work[*ilo], 1, &work[*ilo + n * 3], 1);
    alpha = gamma / sum;
//Determine correction to current iteration
    cmax = Zero;
    for (i = *ilo; i <= *ihi; i++) {
	cor = alpha * work[i + n];
	if (abs(cor) > cmax) {
	    cmax = abs(cor);
	}
	lscale[i] = lscale[i] + cor;
	cor = alpha * work[i];
	if (abs(cor) > cmax) {
	    cmax = abs(cor);
	}
	rscale[i] = rscale[i] + cor;
    }
    if (cmax < Half) {
	goto L350;
    }
    Raxpy(nr, -alpha, &work[*ilo + (n * 2)], 1, &work[*ilo + (n * 4)], 1);
    Raxpy(nr, -alpha, &work[*ilo + n * 3], 1, &work[*ilo + n * 5], 1);
    pgamma = gamma;
    it++;
    if (it <= nrp2) {
	goto L250;
    }
//End generalized conjugate gradient iteration
  L350:
    sfmin = Rlamch("S");
    sfmax = One / sfmin;
    lsfmin = (INTEGER) cast2double(log10(sfmin) / basl + One);
    lsfmax = (INTEGER) cast2double(log10(sfmax) / basl);
    for (i = *ilo; i <= *ihi; i++) {
	irab = iCamax(n - *ilo + 1, &A[i + *ilo * lda], lda);
	rab = abs(A[i + (irab + *ilo - 1) * lda]);
	irab = iCamax(n - *ilo + 1, &B[i + *ilo * ldb], ldb);
	rab = max(rab, abs(B[i + (irab + *ilo - 1) * ldb]));
	lrab = (INTEGER) cast2double(log10(rab + sfmin) / basl + One);
	ir = (INTEGER) cast2double(lscale[i] + sign(Half, lscale[i]));
	ir = min(min(max(ir, lsfmin), lsfmax), lsfmax - lrab);
	lscale[i] = 10 ^ ir;
	icab = iCamax(*ihi, &A[i * lda], 1);
	cab = abs(A[icab + i * lda]);
	icab = iCamax(*ihi, &B[i * ldb + 1], 1);
	mtemp1 = cab, mtemp2 = abs(B[icab + i * ldb]);
	cab = max(mtemp1, mtemp2);
	lcab = (INTEGER) cast2double(log10(cab + sfmin) / basl + One);
	jc = (INTEGER) cast2double(rscale[i] + sign(Half, rscale[i]));
	jc = min(min(max(jc, lsfmin), lsfmax), lsfmax - lcab);
	rscale[i] = 10 ^ jc;
    }
//Row scaling of matrices A and B
    for (i = *ilo; i <= *ihi; i++) {
	CRscal(n - *ilo + 1, lscale[i], &A[i + *ilo * lda], lda);
	CRscal(n - *ilo + 1, lscale[i], &B[i + *ilo * ldb], ldb);
    }
//Column scaling of matrices A and B
    for (j = *ilo; j <= *ihi; j++) {
	CRscal(*ihi, rscale[j], &A[j * lda], 1);
	CRscal(*ihi, rscale[j], &B[j * ldb + 1], 1);
    }
    return;
}
