/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rsbtrd.cpp,v 1.9 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rsbtrd(const char *vect, const char *uplo, INTEGER n, INTEGER kd, REAL * AB, INTEGER ldab, REAL * d, REAL * e, REAL * q, INTEGER ldq, REAL * work, INTEGER * info)
{
    INTEGER i, j, k, l, i2 = 0, j1, j2, nq, nr, kd1, ibl, iqb, kdn, jin, nrt, kdm1, inca, jend, lend, jinc, incx, last;
    REAL temp;
    INTEGER j1end, j1inc, iqend;
    INTEGER initq, wantq, upper;
    INTEGER iqaend;
    REAL Zero = 0.0, One = 1.0;

//Test the input parameters
    initq = Mlsame(vect, "V");
    wantq = initq || Mlsame(vect, "U");
    upper = Mlsame(uplo, "U");
    kd1 = kd + 1;
    kdm1 = kd - 1;
    incx = ldab - 1;
    iqend = 1;
    *info = 0;
    if (!wantq && !Mlsame(vect, "N")) {
	*info = -1;
    } else if (!upper && !Mlsame(uplo, "L")) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (kd < 0) {
	*info = -4;
    } else if (ldab < kd1) {
	*info = -6;
    } else if (ldq < max((INTEGER) 1, n) && wantq) {
	*info = -10;
    }
    if (*info != 0) {
	Mxerbla("Rsbtrd", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0)
	return;

//Initialize Q to the unit matrix, if needed
    if (initq) {
	Rlaset("Full", n, n, Zero, One, &q[0], ldq);
    }
//Wherever possible, plane rotations are generated and applied in
//vector operations of length NR over the index set J1:J2:KD1
//The cosines and sines of the plane rotations are stored in the
// arrays D and WORK.

    inca = kd1 * ldab;
    kdn = min(n - 1, kd);
    if (upper) {
	if (kd > 1) {
//Reduce to tridiagonal form, working with upper triangle
	    nr = 0;
	    j1 = kdn + 2;
	    j2 = 1;
	    for (i = 0; i < n - 2; i++) {
//Reduce i-th row of matrix to tridiagonal form
		for (k = kdn + 1; k >= 2; k--) {
		    j1 += kdn;
		    j2 += kdn;
		    if (nr > 0) {
//generate plane rotations to annihilate nonzero
//elements which have been created outside the band
			Rlargv(nr, &AB[(j1 - 1) * ldab + 1], inca, &work[j1], kd1, &d[j1], kd1);
//apply rotations from the right
//Dependent on the the number of diagonals either
//DLARTV or DROT is used */

			if (nr >= (kd << 1) - 1) {
			    for (l = 0; l < kd - 1; l++) {
				Rlartv(nr, &AB[l + 1 + (j1 - 1) * ldab], inca, &AB[l + j1 * ldab], inca, &d[j1], &work[j1], kd1);
			    }
			} else {
			    jend = j1 + (nr - 1) * kd1;
			    for (jinc = j1; jinc < jend; jinc = jinc + kd1) {
				Rrot(kdm1, &AB[(jinc - 1) * ldab + 2], 1, &AB[jinc * ldab + 1], 1, d[jinc], work[jinc]);
			    }
			}
		    }
		    if (k > 2) {
			if (k <= n - i + 1) {
//generate plane rotation to annihilate a(i,i+k-1)
//within the band

			    Rlartg(AB[kd - k + 3 + (i + k - 2) * ldab], AB[kd - k + 2 + (i + k - 1) * ldab], &d[i + k - 1], &work[i + k - 1], &temp);
			    AB[kd - k + 3 + (i + k - 2) * ldab] = temp;

//apply rotation from the right
			    Rrot(k - 3, &AB[kd - k + 4 + (i + k - 2) * ldab], 1, &AB[kd - k + 3 + (i + k - 1) * ldab], 1, d[i + k - 1], work[i + k - 1]);
			}
			++nr;
			j1 = j1 - kdn - 1;
		    }
//apply plane rotations from both sides to diagonal
//blocks
		    if (nr > 0) {
			Rlar2v(nr, &AB[kd1 + (j1 - 1) * ldab], &AB[kd1 + j1 * ldab], &AB[kd + j1 * ldab], inca, &d[j1], &work[j1], kd1);
		    }
//apply plane rotations from the left
		    if (nr > 0) {
			if ((kd << 1) - 1 < nr) {
//Dependent on the the number of diagonals either
//DLARTV or DROT is used
			    for (l = 0; l < kd - 1; l++) {
				if (j2 + l > n) {
				    nrt = nr - 1;
				} else {
				    nrt = nr;
				}
				if (nrt > 0) {
				    Rlartv(nrt, &AB[kd - l + (j1 + l) * ldab], inca, &AB[kd - l + 1 + (j1 + l) * ldab], inca, &d[j1], &work[j1], kd1);
				}
			    }
			} else {
			    j1end = j1 + kd1 * (nr - 2);
			    if (j1end >= j1) {
				for (jin = j1; jin < j1end; jin = jin + kd1) {
				    Rrot(kd - 1, &AB[kd - 1 + (jin + 1) * ldab], incx, &AB[kd + (jin + 1) * ldab], incx, d[jin], work[jin]);
				}
			    }
			    lend = min(kdm1, n - j2);
			    last = j1end + kd1;
			    if (lend > 0) {
				Rrot(lend, &AB[kd - 1 + (last + 1) * ldab], incx, &AB[kd + (last + 1)
										      * ldab], incx, d[last], work[last]);
			    }
			}
		    }

		    if (wantq) {
//accumulate product of plane rotations in Q
			if (initq) {
//take advantage of the fact that Q was
//initially the Identity matrix
			    iqend = max(iqend, j2);
			    i2 = max((INTEGER) 0, k - 3);
			    iqaend = i * kd + 1;
			    if (k == 2) {
				iqaend += kd;
			    }
			    iqaend = min(iqaend, iqend);
			    for (j = j1; j <= j2; j += kd1) {
				ibl = i - i2 / kdm1;
				i2++;
				iqb = max((INTEGER) 1, j - ibl);
				nq = iqaend + 1 - iqb;
				Rrot(iqaend = min(iqaend + kd, iqend), &q[iqb + (j - 1) * ldq], 1, &q[iqb + j * ldq], 1, d[j], work[j]);
			    }
			} else {
			    for (j = j1; j <= j2; j += kd1) {
				Rrot(n, &q[(j - 1) * ldq + 1], 1, &q[j * ldq + 1], 1, d[j], work[j]);
			    }
			}
		    }

		    if (j2 + kdn > n) {
//adjust J2 to keep within the bounds of the matrix
			--nr;
			j2 = j2 - kdn - 1;
		    }
		    for (j = j1; j <= j2; j += kd1) {
//create nonzero element a(j-1,j+kd) outside the band
//and store it in WORK
			work[j + kd] = work[j] * AB[(j + kd) * ldab + 1];
			AB[(j + kd) * ldab + 1] = d[j] * AB[(j + kd) * ldab + 1];

		    }

		}

	    }
	}

	if (kd > 0) {
//copy off-diagonal elements to E
	    for (i = 0; i < n - 1; i++) {
		e[i] = AB[kd + (i + 1) * ldab];
	    }
	} else {
//set E to zero if original matrix was diagonal
	    for (i = 0; i < n - 1; i++) {
		e[i] = Zero;

	    }
	}
//copy diagonal elements to D
	for (i = 0; i < n; i++) {
	    d[i] = AB[kd1 + i * ldab];
	}
    } else {
	if (kd > 1) {
//Reduce to tridiagonal form, working with lower triangle
	    nr = 0;
	    j1 = kdn + 2;
	    j2 = 1;

	    for (i = 0; i < n - 2; i++) {
//Reduce i-th column of matrix to tridiagonal form
		for (k = kdn + 1; k >= 2; k--) {
		    j1 += kdn;
		    j2 += kdn;

		    if (nr > 0) {
//generate plane rotations to annihilate nonzero
//elements which have been created outside the band
			Rlargv(nr, &AB[kd1 + (j1 - kd1) * ldab], inca, &work[j1], kd1, &d[j1], kd1);

//apply plane rotations from one side
//Dependent on the the number of diagonals either
//DLARTV or DROT is used
			if (nr > (kd << 1) - 1) {
			    for (l = 0; l < kd - 1; l++) {
				Rlartv(nr, &AB[kd1 - l + (j1 - kd1 + l) * ldab], inca, &AB[kd1 - l + 1 + (j1 - kd1 + l) * ldab], inca, &d[j1], &work[j1], kd1);

			    }
			} else {
			    jend = j1 + kd1 * (nr - 1);
			    for (jinc = j1; jinc <= jend; jinc += kd1) {
				Rrot(kdm1, &AB[kd + (jinc - kd) * ldab]
				     , incx, &AB[kd1 + (jinc - kd) * ldab], incx, d[jinc], work[jinc]);

			    }
			}
		    }
		    if (k > 2) {
			if (k <= n - i + 1) {
//generate plane rotation to annihilate a(i+k-1,i)
//within the band
			    Rlartg(AB[k - 1 + i * ldab], AB[k + i * ldab], &d[i + k - 1], &work[i + k - 1], &temp);
			    AB[k - 1 + i * ldab] = temp;
			    Rrot(k - 3, &AB[k - 2 + (i + 1) * ldab], ldab - 1, &AB[k - 1 + (i + 1) * ldab], ldab - 1, d[i + k - 1], work[i + k - 1]);
			}
			++nr;
			j1 = j1 - kdn - 1;
		    }
//apply plane rotations from both sides to diagonal
//blocks
		    if (nr > 0) {
			Rlar2v(nr, &AB[(j1 - 1) * ldab + 1], &AB[j1 * ldab + 1], &AB[(j1 - 1) * ldab + 2], inca, &d[j1], &work[j1], kd1);
		    }
//apply plane rotations from the right
//Dependent on the the number of diagonals either
//DLARTV or DROT is used
		    if (nr > 0) {
			if (nr > (kd << 1) - 1) {
			    for (l = 0; l < kd - 1; l++) {
				if (j2 + l > n) {
				    nrt = nr - 1;
				} else {
				    nrt = nr;
				}
				if (nrt > 0) {
				    Rlartv(nrt, &AB[l + 2 + (j1 - 1) * ldab], inca, &AB[l + 1 + j1 * ldab], inca, &d[j1], &work[j1], kd1);
				}

			    }
			} else {
			    j1end = j1 + kd1 * (nr - 2);
			    if (j1end >= j1) {
				for (j1inc = j1; j1inc <= j1end; j1inc += kd1) {
				    Rrot(kdm1, &AB[(j1inc - 1) * ldab + 3], 1, &AB[j1inc * ldab + 2], 1, d[j1inc], work[j1inc]);
				}
			    }
			    lend = min(kdm1, n - j2);
			    last = j1end + kd1;
			    if (lend > 0) {
				Rrot(lend, &AB[(last - 1) * ldab + 3], 1, &AB[last * ldab + 2], 1, d[last], work[last]);
			    }
			}
		    }
		    if (wantq) {
//accumulate product of plane rotations in Q
			if (initq) {
//take advantage of the fact that Q was
//initially the Identity matrix
			    iqend = max(iqend, j2);
			    i2 = max((INTEGER) 0, k - 32);
			    iqaend = i * kd + 1;
			    if (k == 2) {
				iqaend += kd;
			    }
			    iqaend = min(iqaend, iqend);
			    for (j = j1; j < j2; j += kd1) {
				ibl = i - i2 / kdm1;
				i2++;
				iqb = max((INTEGER) 1, j - ibl);
				nq = iqaend + 1 - iqb;
				iqaend = min(iqaend + kd, iqend);
				Rrot(nq, &q[iqb + (j - 1) * ldq], 1, &q[iqb + j * ldq], 1, d[j], work[j]);
			    }
			} else {
			    for (j = j1; j < i2; j = j + kd1) {
				Rrot(n, &q[(j - 1) * ldq + 1], 1, &q[j * ldq + 1], 1, d[j], work[j]);
			    }
			}
		    }

		    if (j2 + kdn > n) {
//adjust J2 to keep within the bounds of the matrix
			--nr;
			j2 = j2 - kdn - 1;
		    }
		    for (j = j1; j < j2; j += kd1) {
//create nonzero element a(j+kd,j-1) outside the
//band and store it in WORK
			work[j + kd] = work[j] * AB[kd1 + j * ldab];
			AB[kd1 + j * ldab] = d[j] * AB[kd1 + j * ldab];
		    }
		}
	    }
	}

	if (kd > 0) {
//copy off-diagonal elements to E
	    for (i = 0; i < n - 1; i++) {
		e[i] = AB[i * ldab + 2];
	    }
	} else {
//set E to zero if original matrix was diagonal
	    for (i = 0; i < n - 1; i++) {
		e[i] = Zero;
	    }
	}

//copy diagonal elements to D
	for (i = 0; i < n; i++) {
	    d[i] = AB[i * ldab + 1];

	}
    }
    return;
}
