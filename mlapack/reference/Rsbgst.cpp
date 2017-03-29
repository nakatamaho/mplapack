/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rsbgst.cpp,v 1.12 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void
Rsbgst(const char *vect, const char *uplo, INTEGER n, INTEGER ka, INTEGER kb, REAL * AB, INTEGER ldab, REAL * bb, INTEGER ldbb, REAL * x, INTEGER ldx, REAL * work, INTEGER * info)
{
    INTEGER i, j, k, l, m;
    REAL t;
    INTEGER i0 = 0, i1 = 0, i2 = 0, j1, j2;
    REAL ra;
    INTEGER nr, nx, ka1, kb1;
    REAL ra1 = 0;
    INTEGER j1t, j2t;
    REAL bii;
    INTEGER kbt = 0, nrt, inca;
    INTEGER upper, wantx;
    INTEGER update;
    REAL Zero = 0.0, One = 1.0;

//Test the input parameters
    wantx = Mlsame(vect, "V");
    upper = Mlsame(uplo, "U");
    ka1 = ka + 1;
    kb1 = kb + 1;
    *info = 0;
    if (!wantx && !Mlsame(vect, "N")) {
	*info = -1;
    } else if (!upper && !Mlsame(uplo, "L")) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (ka < 0) {
	*info = -4;
    } else if (kb < 0 || kb > ka) {
	*info = -5;
    } else if (ldab < ka + 1) {
	*info = -7;
    } else if (ldbb < kb + 1) {
	*info = -9;
    } else if (ldx < 1 || (wantx && ldx < max((INTEGER) 1, n))) {
	*info = -11;
    }
    if (*info != 0) {
	Mxerbla("Rsbgst", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0)
	return;
    inca = ldab * ka1;
//Initialize X to the unit matrix, if needed
    if (wantx) {
	Rlaset("Full", n, n, Zero, One, &x[0], ldx);
    }
//Set M to the splitting point m. It must be the same value as is
//used in DPBSTF. The chosen value allows the arrays WORK and RWORK
//to be of dimension (N).

    m = (n + kb) / 2;

/*     The routine works in two phases, corresponding to the two halves */
/*     of the split Cholesky factorization of B as S**T*S where */
/*     S = ( U    ) */
/*         ( M  L ) */
/*     with U upper triangular of order m, and L lower triangular of */
/*     order n-m. S has the same bandwidth as B. */
/*     S is treated as a product of elementary matrices: */
/*     S = S(m)*S(m-1)*...*S(2)*S(1)*S(m+1)*S(m+2)*...*S(n-1)*S(n) */
/*     where S(i) is determined by the i-th row of S. */
/*     In phase 1, the index i takes the values n, n-1, ... , m+1; */
/*     in phase 2, it takes the values 1, 2, ... , m. */
/*     For each value of i, the current matrix A is updated by forming */
/*     inv(S(i))**T*A*inv(S(i)). This creates a triangular bulge outside */
/*     the band of A. The bulge is then pushed down toward the bottom of */
/*     A in phase 1, and up toward the top of A in phase 2, by applying */
/*     plane rotations. */
/*     There are kb*(kb+1)/2 elements in the bulge, but at most 2kb-1 */
/*     of them are linearly independent, so annihilating a bulge requires */
/*     only 2kb-1 plane rotations. The rotations are divided into a 1st */
/*     set of kb-1 rotations, and a 2nd set of kb rotations. */
/*     Wherever possible, rotations are generated and applied in vector */
/*     operations of length NR between the indices J1 and J2 (sometimes */
/*     replaced by modified values NRT, J1T or J2T). */
/*     The cosines and sines of the rotations are stored in the array */
/*     WORK. The cosines of the 1st set of rotations are stored in */
/*     elements n+2:n+m-kb-1 and the sines of the 1st set in elements */
/*     2:m-kb-1; the cosines of the 2nd set are stored in elements */
/*     n+m-kb+1:2n and the sines of the second set in elements m-kb+1:n. */
/*     The bulges are not formed explicitly; nonzero elements outside the */
/*     band are created only when they are required for generating new */
/*     rotations; they are stored in the array WORK, in positions where */
/*     they are later overwritten by the sines of the rotations which */
/*     annihilate them. */

/*     **************************** Phase 1 ***************************** */

/*     The logical structure of this phase is: */

/*     UPDATE = .TRUE. */
/*     DO I = N, M + 1, -1 */
/*        use S(i) to update A and create a new bulge */
/*        apply rotations to push all bulges KA positions downward */
/*     END DO */
/*     UPDATE = .FALSE. */
/*     DO I = M + KA + 1, N - 1 */
/*        apply rotations to push all bulges KA positions downward */
/*     END DO */

//To avoid duplicating code, the two loops are merged.

    update = MTRUE;
    i = n + 1;
  L10:
    if (update) {
	i--;
	kbt = min(kb, i - 1);
	i0 = i - 1;
	i1 = min(n, i + ka);
	i2 = i - kbt + ka1;
	if (i < m + 1) {
	    update = MFALSE;
	    i++;
	    i0 = m;
	    if (ka == 0) {
		goto L480;
	    }
	    goto L10;
	}
    } else {
	i += ka;
	if (i > n - 1) {
	    goto L480;
	}
    }

    if (upper) {
//Transform A, working with the upper triangle
	if (update) {
//Form  inv(S(i))**T * A * inv(S(i))

	    bii = bb[kb1 + i * ldbb];
	    for (j = i; j <= i1; j++) {
		AB[i - j + ka1 + j * ldab] /= bii;
	    }
	    for (j = max((INTEGER) 1, i - ka); j <= i; j++) {
		AB[j - i + ka1 + i * ldab] /= bii;
	    }
	    for (k = i - kbt; k <= i - 1; k++) {
		for (j = i - kbt; j <= k; j++) {
		    AB[j - k + ka1 + k * ldab] = AB[j - k + ka1 + k *
						    ldab] - bb[j - i + kb1 + i * ldbb] * AB[k - i + ka1 +
											    i * ldab] - bb[k - i + kb1 +
													   i * ldbb] * AB[j - i +
															  ka1 +
															  i *
															  ldab] +
			AB[ka1 + i * ldab] * bb[j - i + kb1 + i * ldbb] * bb[k - i + kb1 + i * ldbb];
		}
		for (j = max((INTEGER) 1, i - ka); j <= i - kbt - 1; j++) {
		    AB[j - k + ka1 + k * ldab] -= bb[k - i + kb1 + i * ldbb] * AB[j - i + ka1 + i * ldab];
		}
	    }
	    for (j = i; j <= i1; j++) {
		for (k = max(j - ka, i - kbt); k <= i - 1; k++) {
		    AB[k - j + ka1 + j * ldab] -= bb[k - i + kb1 + i * ldbb] * AB[i - j + ka1 + j * ldab];
		}
	    }
	    if (wantx) {
//post-multiply X by inv(S(i))
		Rscal(n - m, One / bii, &x[m + 1 + i * ldx], 1);
		if (kbt > 0) {
		    Rger(n - m, kbt, -One, &x[m + 1 + i * ldx], 1, &bb[kb1 - kbt + i * ldbb], 1, &x[m + 1 + (i - kbt) * ldx], ldx);
		}
	    }
//store a(i,i1) in RA1 for use in next loop over K
	    ra1 = AB[i - i1 + ka1 + i1 * ldab];
	}
//Generate and apply vectors of rotations to chase all the
//existing bulges KA positions down toward the bottom of the
//band
	for (k = 0; k < kb - 1; k++) {
	    if (update) {
//Determine the rotations which would annihilate the bulge
//which has in theory just been created
		if (i - k + ka < n && i - k > 1) {
//generate rotation to annihilate a(i,i-k+ka+1)
		    Rlartg(AB[k + 1 + (i - k + ka) * ldab], ra1, &work[n + i - k + ka - m], &work[i - k + ka - m], &ra);
//create nonzero element a(i-k,i-k+ka+1) outside the
//band and store it in WORK(i-k)
		    t = -bb[kb1 - k + i * ldbb] * ra1;
		    work[i - k] = work[n + i - k + ka - m] * t - work[i - k + ka - m] * AB[(i - k + ka) * ldab + 1];
		    AB[(i - k + ka) * ldab + 1] = work[i - k + ka - m] * t + work[n + i - k + ka - m] * AB[(i - k + ka) * ldab + 1];
		    ra1 = ra;
		}
	    }
	    j2 = i - k - 1 + max((INTEGER) 1, k - i0 + 2) * ka1;
	    nr = (n - j2 + ka) / ka1;
	    j1 = j2 + (nr - 1) * ka1;
	    if (update) {
		j2t = max(j2, i + (ka << 1) - k + 1);
	    } else {
		j2t = j2;
	    }
	    nrt = (n - j2t + ka) / ka1;
	    for (j = j2t; j < j1; j += ka1) {
//create nonzero element a(j-ka,j+1) outside the band
//and store it in WORK(j-m)
		work[j - m] *= AB[(j + 1) * ldab + 1];
		AB[(j + 1) * ldab + 1] = work[n + j - m] * AB[(j + 1) * ldab + 1];
	    }
//generate rotations in 1st set to annihilate elements which
//have been created outside the band
	    if (nrt > 0) {
		Rlargv(nrt, &AB[j2t * ldab + 1], inca, &work[j2t - m], ka1, &work[n + j2t - m], ka1);
	    }
	    if (nr > 0) {
//apply rotations in 1st set from the right
		for (l = 0; l < ka - 1; l++) {
		    Rlartv(nr, &AB[ka1 - l + j2 * ldab], inca, &AB[ka - l + (j2 + 1) * ldab], inca, &work[n + j2 - m], &work[j2 - m], ka1);
		}
//apply rotations in 1st set from both sides to diagonal
//blocks
		Rlar2v(nr, &AB[ka1 + j2 * ldab], &AB[ka1 + (j2 + 1) * ldab], &AB[ka + (j2 + 1) * ldab], inca, &work[n + j2 - m], &work[j2 - m], ka1);
	    }
//start applying rotations in 1st set from the left
	    for (l = ka - 1; l >= kb - k + 1; l--) {
		nrt = (n - j2 + l) / ka1;
		if (nrt > 0) {
		    Rlartv(nrt, &AB[l + (j2 + ka1 - l) * ldab], inca, &AB[l + 1 + (j2 + ka1 - l) * ldab], inca, &work[n + j2 - m], &work[j2 - m], ka1);
		}
	    }
	    if (wantx) {
//post-multiply X by product of rotations in 1st set
		for (j = j2; j < j1; j += ka1) {
		    Rrot(n - m, &x[m + 1 + j * ldx], 1, &x[m + 1 + (j + 1) * ldx], 1, work[n + j - m], work[j - m]);
		}
	    }
	}
	if (update) {
	    if (i2 <= n && kbt > 0) {
//create nonzero element a(i-kbt,i-kbt+ka+1) outside the
//band and store it in WORK(i-kbt)
		work[i - kbt] = -bb[kb1 - kbt + i * ldbb] * ra1;
	    }
	}
	for (k = kb; k >= 1; k--) {
	    if (update) {
		j2 = i - k - 1 + max((INTEGER) 2, k - i0 + 1) * ka1;
	    } else {
		j2 = i - k - 1 + max((INTEGER) 1, k - i0 + 1) * ka1;
	    }
//finish applying rotations in 2nd set from the left
	    for (l = kb - k; l >= 1; l--) {
		nrt = (n - j2 + ka + l) / ka1;
		if (nrt > 0) {
		    Rlartv(nrt, &AB[l + (j2 - l + 1) * ldab], inca, &AB[l + 1 + (j2 - l + 1) * ldab], inca, &work[n + j2 - ka], &work[j2 - ka], ka1);
		}
	    }
	    nr = (n - j2 + ka) / ka1;
	    j1 = j2 + (nr - 1) * ka1;
	    for (j = j1; j < j2; j = j - ka1) {
		work[j] = work[j - ka];
		work[n + j] = work[n + j - ka];

	    }
	    for (j = j2; j < j1; j += ka1) {
//create nonzero element a(j-ka,j+1) outside the band
//and store it in WORK(j)
		work[j] *= AB[(j + 1) * ldab + 1];
		AB[(j + 1) * ldab + 1] = work[n + j] * AB[(j + 1) * ldab + 1];
	    }
	    if (update) {
		if (i - k < n - ka && k <= kbt) {
		    work[i - k + ka] = work[i - k];
		}
	    }
	}
	for (k = kb; k >= 1; k--) {
	    j2 = i - k - 1 + max((INTEGER) 1, k - i0 + 1) * ka1;
	    nr = (n - j2 + ka) / ka1;
	    j1 = j2 + (nr - 1) * ka1;
	    if (nr > 0) {
//generate rotations in 2nd set to annihilate elements
//which have been created outside the band
		Rlargv(nr, &AB[j2 * ldab + 1], inca, &work[j2], ka1, &work[n + j2], ka1);
//apply rotations in 2nd set from the right
		for (l = 0; l < ka - 1; l++) {
		    Rlartv(nr, &AB[ka1 - l + j2 * ldab], inca, &AB[ka - l + (j2 + 1) * ldab], inca, &work[n + j2], &work[j2], ka1);
		}
//apply rotations in 2nd set from both sides to diagonal
//blocks
		Rlar2v(nr, &AB[ka1 + j2 * ldab], &AB[ka1 + (j2 + 1) * ldab], &AB[ka + (j2 + 1) * ldab], inca, &work[n + j2], &work[j2], ka1);
	    }
//start applying rotations in 2nd set from the left
	    for (l = ka - 1; l >= kb - k + 1; l--) {
		nrt = (n - j2 + l) / ka1;
		if (nrt > 0) {
		    Rlartv(nrt, &AB[l + (j2 + ka1 - l) * ldab], inca, &AB[l + 1 + (j2 + ka1 - l) * ldab], inca, &work[n + j2], &work[j2], ka1);
		}
	    }
	    if (wantx) {
//post-multiply X by product of rotations in 2nd set
		for (j = j2; j < j1; j += ka1) {
		    Rrot(n - m, &x[m + 1 + j * ldx], 1, &x[m + 1 + (j + 1) * ldx], 1, work[n + j], work[j]);
		}
	    }
	}
	for (k = 0; k < kb - 1; k++) {
	    j2 = i - k - 1 + max((INTEGER) 1, k - i0 + 2) * ka1;
//finish applying rotations in 1st set from the left
	    for (l = kb - k; l >= 1; l--) {
		nrt = (n - j2 + l) / ka1;
		if (nrt > 0) {
		    Rlartv(nrt, &AB[l + (j2 + ka1 - l) * ldab], inca, &AB[l + 1 + (j2 + ka1 - l) * ldab], inca, &work[n + j2 - m], &work[j2 - m], ka1);
		}
	    }
	}
	if (kb > 1) {
	    for (j = n - 1; j >= i - kb + (ka << 1) + 1; j--) {
		work[n + j - m] = work[n + j - ka - m];
		work[j - m] = work[j - ka - m];
	    }
	}
    } else {
//Transform A, working with the lower triangle
	if (update) {
//Form  inv(S(i))**T * A * inv(S(i))
	    bii = bb[i * ldbb + 1];
	    for (j = i; j <= i1; j++) {
		AB[j - i + 1 + i * ldab] /= bii;
	    }
	    for (j = max((INTEGER) 1, i - ka); j <= i; j++) {
		AB[i - j + 1 + j * ldab] /= bii;

	    }
	    for (k = i - kbt; k <= i - 1; k++) {
		for (j = i - kbt; j <= k; j++) {
		    AB[k - j + 1 + j * ldab] = AB[k - j + 1 + j * ldab]
			- bb[i - j + 1 + j * ldbb] * AB[i - k + 1
							+ k * ldab] - bb[i - k + 1 + k * ldbb] * AB[i - j + 1 +
												    j * ldab] + AB[i * ldab +
														   1] * bb[i - j + 1 + j * ldbb] * bb[i - k + 1 + k * ldbb];

		}
		for (j = max((INTEGER) 1, i - ka); j <= i - kbt - 1; j++) {
		    AB[k - j + 1 + j * ldab] -= bb[i - k + 1 + k * ldbb] * AB[i - j + 1 + j * ldab];
		}
	    }
	    for (j = i; j <= i1; j++) {
		for (k = max(j - ka, i - kbt); k <= i - 1; k++) {
		    AB[j - k + 1 + k * ldab] -= bb[i - k + 1 + k * ldbb] * AB[j - i + 1 + i * ldab];
		}
	    }
	    if (wantx) {
//post-multiply X by inv(S(i))
		Rscal(n - m, One / bii, &x[m + 1 + i * ldx], 1);
		if (kbt > 0) {
		    Rger(n - m, kbt, -One, &x[m + 1 + i * ldx], 1, &bb[kbt + 1 + (i - kbt) * ldbb], ldbb - 1, &x[m + 1 + (i - kbt) * ldx], ldx);
		}
	    }
//store a(i1,i) in RA1 for use in next loop over K
	    ra1 = AB[i1 - i + 1 + i * ldab];
	}
//Generate and apply vectors of rotations to chase all the
//existing bulges KA positions down toward the bottom of the
//band
	for (k = 0; k < kb - 1; k++) {
	    if (update) {
//Determine the rotations which would annihilate the bulge
//which has in theory just been created
		if (i - k + ka < n && i - k > 1) {
//generate rotation to annihilate a(i-k+ka+1,i)
		    Rlartg(AB[ka1 - k + i * ldab], ra1, &work[n + i - k + ka - m], &work[i - k + ka - m], &ra);

//create nonzero element a(i-k+ka+1,i-k) outside the
//band and store it in WORK(i-k)
		    t = -bb[k + 1 + (i - k) * ldbb] * ra1;
		    work[i - k] = work[n + i - k + ka - m] * t - work[i - k + ka - m] * AB[ka1 + (i - k) * ldab];
		    AB[ka1 + (i - k) * ldab] = work[i - k + ka - m] * t + work[n + i - k + ka - m] * AB[ka1 + (i - k) * ldab];
		    ra1 = ra;
		}
	    }
	    j2 = i - k - 1 + max((INTEGER) 1, k - i0 + 2) * ka1;
	    nr = (n - j2 + ka) / ka1;
	    j1 = j2 + (nr - 1) * ka1;
	    if (update) {
		j2t = max(j2, i + (ka << 1) - k + 1);
	    } else {
		j2t = j2;
	    }
	    nrt = (n - j2t + ka) / ka1;
	    for (j = j2t; j < j1; j += ka1) {
//create nonzero element a(j+1,j-ka) outside the band
//and store it in WORK(j-m)
		work[j - m] *= AB[ka1 + (j - ka + 1) * ldab];
		AB[ka1 + (j - ka + 1) * ldab] = work[n + j - m] * AB[ka1 + (j - ka + 1) * ldab];
	    }
//generate rotations in 1st set to annihilate elements which
//have been created outside the band
	    if (nrt > 0) {
		Rlargv(nrt, &AB[ka1 + (j2t - ka) * ldab], inca, &work[j2t - m], ka1, &work[n + j2t - m], ka1);
	    }
	    if (nr > 0) {
//apply rotations in 1st set from the left
		for (l = 0; l < ka - 1; l++) {
		    Rlartv(nr, &AB[l + 1 + (j2 - l) * ldab], inca, &AB[l + 2 + (j2 - l) * ldab], inca, &work[n + j2 - m], &work[j2 - m], ka1);
		}
//apply rotations in 1st set from both sides to diagonal
//blocks
		Rlar2v(nr, &AB[j2 * ldab + 1], &AB[(j2 + 1) * ldab + 1], &AB[j2 * ldab + 2], inca, &work[n + j2 - m], &work[j2 - m], ka1);
	    }
//start applying rotations in 1st set from the right
	    for (l = ka - 1; l >= kb - k + 1; l--) {
		nrt = (n - j2 + l) / ka1;
		if (nrt > 0) {
		    Rlartv(nrt, &AB[ka1 - l + 1 + j2 * ldab], inca, &AB[ka1 - l + (j2 + 1) * ldab], inca, &work[n + j2 - m], &work[j2 - m], ka1);
		}
	    }
	    if (wantx) {
//post-multiply X by product of rotations in 1st set
		for (j = j2; j <= j1; j += ka1) {
		    Rrot(n - m, &x[m + 1 + j * ldx], 1, &x[m + 1 + (j + 1) * ldx], 1, work[n + j - m], work[j - m]);
		}
	    }
	}
	if (update) {
	    if (i2 <= n && kbt > 0) {
//create nonzero element a(i-kbt+ka+1,i-kbt) outside the
//band and store it in WORK(i-kbt)
		work[i - kbt] = -bb[kbt + 1 + (i - kbt) * ldbb] * ra1;
	    }
	}
	for (k = kb; k >= 1; k--) {
	    if (update) {
		j2 = i - k - 1 + max((INTEGER) 2, k - i0 + 1) * ka1;
	    } else {
		j2 = i - k - 1 + max((INTEGER) 1, k - i0 + 1) * ka1;
	    }
//finish applying rotations in 2nd set from the right
	    for (l = kb - k; l >= 1; l--) {
		nrt = (n - j2 + ka + l) / ka1;
		if (nrt > 0) {
		    Rlartv(nrt, &AB[ka1 - l + 1 + (j2 - ka) * ldab], inca, &AB[ka1 - l + (j2 - ka + 1) * ldab], inca, &work[n + j2 - ka], &work[j2 - ka], ka1);
		}
	    }
	    nr = (n - j2 + ka) / ka1;
	    j1 = j2 + (nr - 1) * ka1;
	    for (j = j1; j <= j2; j = j - ka1) {
		work[j] = work[j - ka];
		work[n + j] = work[n + j - ka];
	    }
	    for (j = j2; j < j1; j += ka1) {
//create nonzero element a(j+1,j-ka) outside the band
//and store it in WORK(j)
		work[j] *= AB[ka1 + (j - ka + 1) * ldab];
		AB[ka1 + (j - ka + 1) * ldab] = work[n + j] * AB[ka1 + (j - ka + 1) * ldab];
	    }
	    if (update) {
		if (i - k < n - ka && k <= kbt) {
		    work[i - k + ka] = work[i - k];
		}
	    }
	}
	for (k = kb; k >= 1; k--) {
	    j2 = i - k - 1 + max((INTEGER) 1, k - i0 + 1) * ka1;
	    nr = (n - j2 + ka) / ka1;
	    j1 = j2 + (nr - 1) * ka1;
	    if (nr > 0) {
//generate rotations in 2nd set to annihilate elements
//which have been created outside the band
		Rlargv(nr, &AB[ka1 + (j2 - ka) * ldab], inca, &work[j2], ka1, &work[n + j2], ka1);
//apply rotations in 2nd set from the left
		for (l = 0; l < ka - 1; l++) {
		    Rlartv(nr, &AB[l + 1 + (j2 - l) * ldab], inca, &AB[l + 2 + (j2 - l) * ldab], inca, &work[n + j2], &work[j2], ka1);
		}
//apply rotations in 2nd set from both sides to diagonal
//blocks
		Rlar2v(nr, &AB[j2 * ldab + 1], &AB[(j2 + 1) * ldab + 1], &AB[j2 * ldab + 2], inca, &work[n + j2], &work[j2], ka1);

	    }
//start applying rotations in 2nd set from the right
	    for (l = ka - 1; l >= kb - k + 1; l--) {
		nrt = (n - j2 + l) / ka1;
		if (nrt > 0) {
		    Rlartv(nrt, &AB[ka1 - l + 1 + j2 * ldab], inca, &AB[ka1 - l + (j2 + 1) * ldab], inca, &work[n + j2], &work[j2], ka1);
		}
	    }
	    if (wantx) {
//post-multiply X by product of rotations in 2nd set
		for (j = j2; j <= j1; j += ka1) {
		    Rrot(n - m, &x[m + 1 + j * ldx], 1, &x[m + 1 + (j + 1) * ldx], 1, work[n + j], work[j]);
		}
	    }
	}
	for (k = 0; k < kb - 1; k++) {
	    j2 = i - k - 1 + max((INTEGER) 1, k - i0 + 2) * ka1;
//finish applying rotations in 1st set from the right
	    for (l = kb - k; l >= 1; l--) {
		nrt = (n - j2 + l) / ka1;
		if (nrt > 0) {
		    Rlartv(nrt, &AB[ka1 - l + 1 + j2 * ldab], inca, &AB[ka1 - l + (j2 + 1) * ldab], inca, &work[n + j2 - m], &work[j2 - m], ka1);
		}
	    }
	}
	if (kb > 1) {
	    for (j = n - 1; j >= i - kb + (ka << 1) + 1; j--) {
		work[n + j - m] = work[n + j - ka - m];
		work[j - m] = work[j - ka - m];
	    }
	}
    }
    goto L10;

  L480:
/*     **************************** Phase 2 ***************************** */
/*     The logical structure of this phase is: */
/*     UPDATE = .TRUE. */
/*     DO I = 1, M */
/*        use S(i) to update A and create a new bulge */
/*        apply rotations to push all bulges KA positions upward */
/*     END DO */
/*     UPDATE = .FALSE. */
/*     DO I = M - KA - 1, 2, -1 */
/*        apply rotations to push all bulges KA positions upward */
/*     END DO */
/*     To avoid duplicating code, the two loops are merged. */

    update = MTRUE;
    i = 0;
  L490:
    if (update) {
	i++;
	kbt = min(kb, m - i);
	i0 = i + 1;
	i1 = max((INTEGER) 1, i - ka);
	i2 = i + kbt - ka1;
	if (i > m) {
	    update = MFALSE;
	    i--;
	    i0 = m + 1;
	    if (ka == 0) {
		return;
	    }
	    goto L490;
	}
    } else {
	i -= ka;
	if (i < 2) {
	    return;
	}
    }

    if (i < m - kbt) {
	nx = m;
    } else {
	nx = n;
    }

    if (upper) {
//Transform A, working with the upper triangle
	if (update) {
//Form  inv(S(i))**T * A * inv(S(i))
	    bii = bb[kb1 + i * ldbb];
	    for (j = i1; j <= i; j++) {
		AB[j - i + ka1 + i * ldab] /= bii;
	    }
	    for (j = i; j <= min(n, i + ka1); j++) {
		AB[i - j + ka1 + j * ldab] /= bii;

	    }
	    for (k = i + 1; k <= i + kbt; k++) {
		for (j = k; j <= i + kbt; j++) {
		    AB[k - j + ka1 + j * ldab] = AB[k - j + ka1 + j *
						    ldab] - bb[i - j + kb1 + j * ldbb] * AB[i - k + ka1 +
											    k * ldab] - bb[i - k + kb1 +
													   k * ldbb] * AB[i - j +
															  ka1 +
															  j *
															  ldab] +
			AB[ka1 + i * ldab] * bb[i - j + kb1 + j * ldbb] * bb[i - k + kb1 + k * ldbb];

		}
		for (j = i + kbt + 1; j <= min(n, i + ka); j++) {
		    AB[k - j + ka1 + j * ldab] -= bb[i - k + kb1 + k * ldbb] * AB[i - j + ka1 + j * ldab];
		}
	    }
	    for (j = i1; j <= i; j++) {
		for (k = i + 1; k <= min(j + ka, i + kbt); k++) {
		    AB[j - k + ka1 + k * ldab] -= bb[i - k + kb1 + k * ldbb] * AB[j - i + ka1 + i * ldab];

		}
	    }
	    if (wantx) {
//post-multiply X by inv(S(i))
		Rscal(nx, One / bii, &x[i * ldx + 1], 1);
		if (kbt > 0) {
		    Rger(nx, kbt, -One, &x[i * ldx + 1], 1, &bb[kb + (i + 1) * ldbb], ldbb - 1, &x[(i + 1) * ldx + 1], ldx);
		}
	    }
//store a(i1,i) in RA1 for use in next loop over K
	    ra1 = AB[i1 - i + ka1 + i * ldab];
	}
//Generate and apply vectors of rotations to chase all the
//existing bulges KA positions up toward the top of the band
	for (k = 0; k < kb - 1; k++) {
	    if (update) {
//Dtermine the rotations which would annihilate the bulge
//which has in theory just been created
		if (i + k - ka1 > 0 && i + k < m) {
//generate rotation to annihilate a(i+k-ka-1,i)
		    Rlartg(AB[k + 1 + i * ldab], ra1, &work[n + i + k - ka], &work[i + k - ka], &ra);
//create nonzero element a(i+k-ka-1,i+k) outside the
//band and store it in WORK(m-kb+i+k)
		    t = -bb[kb1 - k + (i + k) * ldbb] * ra1;
		    work[m - kb + i + k] = work[n + i + k - ka] * t - work[i + k - ka] * AB[(i + k) * ldab + 1];
		    AB[(i + k) * ldab + 1] = work[i + k - ka] * t + work[n + i + k - ka] * AB[(i + k) * ldab + 1];
		    ra1 = ra;
		}
	    }
	    j2 = i + k + 1 - max((INTEGER) 1, k + i0 - m + 1) * ka1;
	    nr = (j2 + ka - 1) / ka1;
	    j1 = j2 - (nr - 1) * ka1;
	    if (update) {
		j2t = min(j2, i - (ka << 1) + k - 1);
	    } else {
		j2t = j2;
	    }
	    nrt = (j2t + ka - 1) / ka1;
	    for (j = j1; j <= j2t; j += ka1) {
//create nonzero element a(j-1,j+ka) outside the band
//and store it in WORK(j)
		work[j] *= AB[(j + ka - 1) * ldab + 1];
		AB[(j + ka - 1) * ldab + 1] = work[n + j] * AB[(j + ka - 1) * ldab + 1];
	    }
//generate rotations in 1st set to annihilate elements which
//have been created outside the band

	    if (nrt > 0) {
		Rlargv(nrt, &AB[(j1 + ka) * ldab + 1], inca, &work[j1], ka1, &work[n + j1], ka1);
	    }
	    if (nr > 0) {
//apply rotations in 1st set from the left
		for (l = 0; l < ka - 1; l++) {
		    Rlartv(nr, &AB[ka1 - l + (j1 + l) * ldab], inca, &AB[ka - l + (j1 + l) * ldab], inca, &work[n + j1], &work[j1], ka1);

		}
//apply rotations in 1st set from both sides to diagonal
//blocks
		Rlar2v(nr, &AB[ka1 + j1 * ldab], &AB[ka1 + (j1 - 1) * ldab], &AB[ka + j1 * ldab], inca, &work[n + j1], &work[j1], ka1);
	    }
//start applying rotations in 1st set from the right
	    for (l = ka - 1; l >= kb - k + 1; l--) {
		nrt = (j2 + l - 1) / ka1;
		j1t = j2 - (nrt - 1) * ka1;
		if (nrt > 0) {
		    Rlartv(nrt, &AB[l + j1t * ldab], inca, &AB[l + 1 + (j1t - 1) * ldab], inca, &work[n + j1t], &work[j1t], ka1);
		}
	    }
	    if (wantx) {
//post-multiply X by product of rotations in 1st set
		for (j = j1; j <= j2; j += ka1) {
		    Rrot(nx, &x[j * ldx + 1], 1, &x[(j - 1) * ldx + 1], 1, work[n + j], work[j]);
		}
	    }
	}
	if (update) {
	    if (i2 > 0 && kbt > 0) {
//create nonzero element a(i+kbt-ka-1,i+kbt) outside the
//band and store it in WORK(m-kb+i+kbt)
		work[m - kb + i + kbt] = -bb[kb1 - kbt + (i + kbt) * ldbb] * ra1;
	    }
	}
	for (k = kb; k >= 1; k--) {
	    if (update) {
		j2 = i + k + 1 - max((INTEGER) 2, k + i0 - m) * ka1;
	    } else {
		j2 = i + k + 1 - max((INTEGER) 1, k + i0 - m) * ka1;
	    }
//finish applying rotations in 2nd set from the right
	    for (l = kb - k; l >= 1; l--) {
		nrt = (j2 + ka + l - 1) / ka1;
		j1t = j2 - (nrt - 1) * ka1;
		if (nrt > 0) {
		    Rlartv(nrt, &AB[l + (j1t + ka) * ldab], inca, &AB[l + 1 + (j1t + ka - 1) * ldab], inca, &work[n + m - kb + j1t + ka], &work[m - kb + j1t + ka], ka1);
		}

	    }
	    nr = (j2 + ka - 1) / ka1;
	    j1 = j2 - (nr - 1) * ka1;
	    for (j = j1; j <= j2; j += ka1) {
		work[m - kb + j] = work[m - kb + j + ka];
		work[n + m - kb + j] = work[n + m - kb + j + ka];

	    }
	    for (j = j1; j <= j2; j += ka1) {
//create nonzero element a(j-1,j+ka) outside the band
//and store it in WORK(m-kb+j)
		work[m - kb + j] *= AB[(j + ka - 1) * ldab + 1];
		AB[(j + ka - 1) * ldab + 1] = work[n + m - kb + j] * AB[(j + ka - 1) * ldab + 1];
	    }
	    if (update) {
		if (i + k > ka1 && k <= kbt) {
		    work[m - kb + i + k - ka] = work[m - kb + i + k];
		}
	    }
	}
	for (k = kb; k >= 1; k--) {
	    j2 = i + k + 1 - max((INTEGER) 1, k + i0 - m) * ka1;
	    nr = (j2 + ka - 1) / ka1;
	    j1 = j2 - (nr - 1) * ka1;
	    if (nr > 0) {

//generate rotations in 2nd set to annihilate elements
//which have been created outside the band
		Rlargv(nr, &AB[(j1 + ka) * ldab + 1], inca, &work[m - kb + j1], ka1, &work[n + m - kb + j1], ka1);
//apply rotations in 2nd set from the left
		for (l = 0; l < ka - 1; l++) {
		    Rlartv(nr, &AB[ka1 - l + (j1 + l) * ldab], inca, &AB[ka - l + (j1 + l) * ldab], inca, &work[n + m - kb + j1], &work[m - kb + j1], ka1);

		}
//apply rotations in 2nd set from both sides to diagonal
//blocks
		Rlar2v(nr, &AB[ka1 + j1 * ldab], &AB[ka1 + (j1 - 1) * ldab], &AB[ka + j1 * ldab], inca, &work[n + m - kb + j1], &work[m - kb + j1], ka1);

	    }
//start applying rotations in 2nd set from the right
	    for (l = ka - 1; l >= kb - k + 1; l--) {
		nrt = (j2 + l - 1) / ka1;
		j1t = j2 - (nrt - 1) * ka1;
		if (nrt > 0) {
		    Rlartv(nrt, &AB[l + j1t * ldab], inca, &AB[l + 1 + (j1t - 1) * ldab], inca, &work[n + m - kb + j1t], &work[m - kb + j1t], ka1);
		}
	    }
	    if (wantx) {
//post-multiply X by product of rotations in 2nd set
		for (j = j1; j <= j2; j += ka1) {
		    Rrot(nx, &x[j * ldx + 1], 1, &x[(j - 1) * ldx + 1], 1, work[n + m - kb + j], work[m - kb + j]);
		}
	    }
	}
	for (k = 0; k < kb - 1; k++) {
	    j2 = i + k + 1 - max((INTEGER) 1, k + i0 - m + 1) * ka1;
//finish applying rotations in 1st set from the right
	    for (l = kb - k; l >= 1; l--) {
		nrt = (j2 + l - 1) / ka1;
		j1t = j2 - (nrt - 1) * ka1;
		if (nrt > 0) {
		    Rlartv(nrt, &AB[l + j1t * ldab], inca, &AB[l + 1 + (j1t - 1) * ldab], inca, &work[n + j1t], &work[j1t], ka1);
		}
	    }
	}

	if (kb > 1) {
	    for (j = 2; j <= min(i + kb, m) - (ka << 1) - 1; j++) {
		work[n + j] = work[n + j + ka];
		work[j] = work[j + ka];
	    }
	}
    } else {
//Transform A, working with the lower triangle
	if (update) {
//Form  inv(S(i))**T * A * inv(S(i))
	    bii = bb[i * ldbb + 1];
	    for (j = i1; j <= i; j++) {
		AB[i - j + 1 + j * ldab] /= bii;

	    }
	    for (j = i; j <= min(n, i + ka); j++) {
		AB[j - i + 1 + i * ldab] /= bii;

	    }
	    for (k = i + 1; k <= i + kbt; k++) {
		for (j = k; j <= i + kbt; j++) {
		    AB[j - k + 1 + k * ldab] = AB[j - k + 1 + k * ldab]
			- bb[j - i + 1 + i * ldbb] * AB[k - i + 1 + i * ldab] - bb[k - i + 1 + i * ldbb] * AB[j - i + 1 + i * ldab] + AB[i * ldab + 1] * bb[j - i + 1 + i * ldbb]
			* bb[k - i + 1 + i * ldbb];

		}
		for (j = i + kbt + 1; j <= min(n, i + ka); j++) {
		    AB[j - k + 1 + k * ldab] -= bb[k - i + 1 + i * ldbb] * AB[j - i + 1 + i * ldab];

		}
	    }
	    for (j = i1; j <= i; j++) {
		for (k = i + 1; k <= min(j + ka, i + kbt); k++) {
		    AB[k - j + 1 + j * ldab] -= bb[k - i + 1 + i * ldbb] * AB[i - j + 1 + j * ldab];
		}
	    }
	    if (wantx) {
//post-multiply X by inv(S(i))
		Rscal(nx, One / bii, &x[i * ldx + 1], 1);
		if (kbt > 0) {
		    Rger(nx, kbt, -One, &x[i * ldx + 1], 1, &bb[i * ldbb + 2], 1, &x[(i + 1) * ldx + 1], ldx);
		}
	    }
//store a(i,i1) in RA1 for use in next loop over K
	    ra1 = AB[i - i1 + 1 + i1 * ldab];
	}
//Generate and apply vectors of rotations to chase all the
//existing bulges KA positions up toward the top of the band
	for (k = 0; k < kb - 1; k++) {
	    if (update) {
//Determine the rotations which would annihilate the bulge
//which has in theory just been created
		if (i + k - ka1 > 0 && i + k < m) {
//generate rotation to annihilate a(i,i+k-ka-1)
		    Rlartg(AB[ka1 - k + (i + k - ka) * ldab], ra1, &work[n + i + k - ka], &work[i + k - ka], &ra);

//create nonzero element a(i+k,i+k-ka-1) outside the
//band and store it in WORK(m-kb+i+k)
		    t = -bb[k + 1 + i * ldbb] * ra1;
		    work[m - kb + i + k] = work[n + i + k - ka] * t - work[i + k - ka] * AB[ka1 + (i + k - ka) * ldab];
		    AB[ka1 + (i + k - ka) * ldab] = work[i + k - ka]
			* t + work[n + i + k - ka] * AB[ka1 + (i + k - ka) * ldab];
		    ra1 = ra;
		}
	    }
	    j2 = i + k + 1 - max((INTEGER) 1, k + i0 - m + 1) * ka1;
	    nr = (j2 + ka - 1) / ka1;
	    j1 = j2 - (nr - 1) * ka1;
	    if (update) {
		j2t = min(j2, i - (ka << 1) + k - 1);
	    } else {
		j2t = j2;
	    }
	    nrt = (j2t + ka - 1) / ka1;
	    for (j = j1; j <= j2; j += ka1) {
//create nonzero element a(j+ka,j-1) outside the band
//and store it in WORK(j)
		work[j] *= AB[ka1 + (j - 1) * ldab];
		AB[ka1 + (j - 1) * ldab] = work[n + j] * AB[ka1 + (j - 1) * ldab];
	    }
//generate rotations in 1st set to annihilate elements which
//have been created outside the band
	    if (nrt > 0) {
		Rlargv(nrt, &AB[ka1 + j1 * ldab], inca, &work[j1], ka1, &work[n + j1], ka1);
	    }
	    if (nr > 0) {
//apply rotations in 1st set from the right
		for (l = 0; l < ka - 1; l++) {
		    Rlartv(nr, &AB[l + 1 + j1 * ldab], inca, &AB[l + 2 + (j1 - 1) * ldab], inca, &work[n + j1], &work[j1], ka1);
		}
//apply rotations in 1st set from both sides to diagonal
//blocks
		Rlar2v(nr, &AB[j1 * ldab + 1], &AB[(j1 - 1) * ldab + 1], &AB[(j1 - 1) * ldab + 2], inca, &work[n + j1]
		       , &work[j1], ka1);
	    }
//start applying rotations in 1st set from the left
	    for (l = ka - 1; l >= kb - k + 1; l--) {
		nrt = (j2 + l - 1) / ka1;
		j1t = j2 - (nrt - 1) * ka1;
		if (nrt > 0) {
		    Rlartv(nrt, &AB[ka1 - l + 1 + (j1t - ka1 + l) * ldab]
			   , inca, &AB[ka1 - l + (j1t - ka1 + l) * ldab], inca, &work[n + j1t], &work[j1t], ka1);
		}
	    }
	    if (wantx) {
//post-multiply X by product of rotations in 1st set
		for (j = j1; j <= j2; j += ka1) {
		    Rrot(nx, &x[j * ldx + 1], 1, &x[(j - 1) * ldx + 1], 1, work[n + j], work[j]);
		}
	    }
	}
	if (update) {
	    if (i2 > 0 && kbt > 0) {
//create nonzero element a(i+kbt,i+kbt-ka-1) outside the
//band and store it in WORK(m-kb+i+kbt)
		work[m - kb + i + kbt] = -bb[kbt + 1 + i * ldbb] * ra1;
	    }
	}
	for (k = kb; k >= 1; k--) {
	    if (update) {
		j2 = i + k + 1 - max((INTEGER) 2, k + i0 - m) * ka1;
	    } else {
		j2 = i + k + 1 - max((INTEGER) 1, k + i0 - m) * ka1;
	    }
//finish applying rotations in 2nd set from the left
	    for (l = kb - k; l >= 1; l--) {
		nrt = (j2 + ka + l - 1) / ka1;
		j1t = j2 - (nrt - 1) * ka1;
		if (nrt > 0) {
		    Rlartv(nrt, &AB[ka1 - l + 1 + (j1t + l - 1) * ldab],
			   inca, &AB[ka1 - l + (j1t + l - 1) * ldab], inca, &work[n + m - kb + j1t + ka], &work[m - kb + j1t + ka], ka1);
		}

	    }
	    nr = (j2 + ka - 1) / ka1;
	    j1 = j2 - (nr - 1) * ka1;
	    for (j = j1; j <= j2; j += ka1) {
		work[m - kb + j] = work[m - kb + j + ka];
		work[n + m - kb + j] = work[n + m - kb + j + ka];

	    }
	    for (j = j1; j < j2; j += ka1) {
//create nonzero element a(j+ka,j-1) outside the band
//and store it in WORK(m-kb+j)
		work[m - kb + j] *= AB[ka1 + (j - 1) * ldab];
		AB[ka1 + (j - 1) * ldab] = work[n + m - kb + j] * AB[ka1 + (j - 1) * ldab];
	    }
	    if (update) {
		if (i + k > ka1 && k <= kbt) {
		    work[m - kb + i + k - ka] = work[m - kb + i + k];
		}
	    }
	}
	for (k = kb; k >= 1; k--) {
	    j2 = i + k + 1 - max((INTEGER) 1, k + i0 - m) * ka1;
	    nr = (j2 + ka - 1) / ka1;
	    j1 = j2 - (nr - 1) * ka1;
	    if (nr > 0) {
//generate rotations in 2nd set to annihilate elements
//which have been created outside the band
		Rlargv(nr, &AB[ka1 + j1 * ldab], inca, &work[m - kb + j1], ka1, &work[n + m - kb + j1], ka1);
//apply rotations in 2nd set from the right
		for (l = 0; l < ka - 1; l++) {
		    Rlartv(nr, &AB[l + 1 + j1 * ldab], inca, &AB[l + 2 + (j1 - 1) * ldab], inca, &work[n + m - kb + j1], &work[m - kb + j1], ka1);
		}

//apply rotations in 2nd set from both sides to diagonal
//blocks
		Rlar2v(nr, &AB[j1 * ldab + 1], &AB[(j1 - 1) * ldab + 1], &AB[(j1 - 1) * ldab + 2], inca, &work[n + m - kb + j1], &work[m - kb + j1], ka1);
	    }
//start applying rotations in 2nd set from the left
	    for (l = ka - 1; l >= kb - k + 1; l--) {
		nrt = (j2 + l - 1) / ka1;
		j1t = j2 - (nrt - 1) * ka1;
		if (nrt > 0) {
		    Rlartv(nrt, &AB[ka1 - l + 1 + (j1t - ka1 + l) * ldab]
			   , inca, &AB[ka1 - l + (j1t - ka1 + l) * ldab], inca, &work[n + m - kb + j1t], &work[m - kb + j1t], ka1);
		}
	    }
	    if (wantx) {
//post-multiply X by product of rotations in 2nd set
		for (j = j1; j < j2; j += ka1) {
		    Rrot(nx, &x[j * ldx + 1], 1, &x[(j - 1) * ldx + 1], 1, work[n + m - kb + j], work[m - kb + j]);
		}
	    }
	}
	for (k = 0; k < kb - 1; k++) {
	    j2 = i + k + 1 - max((INTEGER) 1, k + i0 - m + 1) * ka1;
//finish applying rotations in 1st set from the left
	    for (l = kb - k; l >= 1; l--) {
		nrt = (j2 + l - 1) / ka1;
		j1t = j2 - (nrt - 1) * ka1;
		if (nrt > 0) {
		    Rlartv(nrt, &AB[ka1 - l + 1 + (j1t - ka1 + l) * ldab]
			   , inca, &AB[ka1 - l + (j1t - ka1 + l) * ldab], inca, &work[n + j1t], &work[j1t], ka1);
		}
	    }
	}
	if (kb > 1) {
	    for (j = 2; j <= min(i + kb, m) - (ka << 1) - 1; j++) {
		work[n + j] = work[n + j + ka];
		work[j] = work[j + ka];
	    }
	}
    }
    goto L490;
}
