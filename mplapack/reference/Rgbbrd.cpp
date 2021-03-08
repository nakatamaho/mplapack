/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgbbrd.cpp,v 1.3 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rgbbrd(const char *vect, INTEGER m, INTEGER n, INTEGER ncc,
	    INTEGER kl, INTEGER ku, REAL * AB, INTEGER ldab, REAL * d, REAL * e, REAL * q, INTEGER ldq, REAL * pt, INTEGER ldpt, REAL * c, INTEGER ldc, REAL * work, INTEGER * info)
{
    INTEGER i, j, l, j1, j2, kb;
    REAL ra, rb, rc;
    INTEGER kk, ml, mn, nr, mu;
    REAL rs;
    INTEGER kb1, ml0, mu0, klm, kun, nrt, klu1, inca;
    INTEGER wantb, wantc;
    INTEGER minmn;
    INTEGER wantq;
    INTEGER wantpt;
    REAL One = 1.0, Zero = 0.0;

//Test the input parameters
    wantb = Mlsame(vect, "B");
    wantq = Mlsame(vect, "Q") || wantb;
    wantpt = Mlsame(vect, "P") || wantb;
    wantc = ncc > 0;
    klu1 = kl + ku + 1;
    *info = 0;
    if (!wantq && !wantpt && !Mlsame(vect, "N")) {
	*info = -1;
    } else if (m < 0) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (ncc < 0) {
	*info = -4;
    } else if (kl < 0) {
	*info = -5;
    } else if (ku < 0) {
	*info = -6;
    } else if (ldab < klu1) {
	*info = -8;
    } else if (ldq < 1 || (wantq && ldq < max((INTEGER) 1, m))) {
	*info = -12;
    } else if (ldpt < 1 || (wantpt && ldpt < max((INTEGER) 1, n))) {
	*info = -14;
    } else if (ldc < 1 || (wantc && ldc < max((INTEGER) 1, m))) {
	*info = -16;
    }
    if (*info != 0) {
	Mxerbla("Rgbbrd", -(*info));
	return;
    }
//Initialize Q and P' to the unit matrix, if needed
    if (wantq) {
	Rlaset("Full", m, m, Zero, One, &q[0], ldq);
    }
    if (wantpt) {
	Rlaset("Full", n, n, Zero, One, &pt[0], ldpt);
    }
//Quick return if possible.
    if (m == 0 || n == 0) {
	return;
    }
    minmn = min(m, n);
    if (kl + ku > 1) {
//Reduce to upper bidiagonal form if KU > 0; if KU = 0, reduce
//first to lower bidiagonal form and then transform to upper
//bidiagonal
	if (ku > 0) {
	    ml0 = 1;
	    mu0 = 2;
	} else {
	    ml0 = 2;
	    mu0 = 1;
	}
//Wherever possible, plane rotations are generated and applied in
//vector operations of length NR over the index set J1:J2:KLU1
//The sines of the plane rotations are stored in WORK(1:max(m,n))
//and the cosines in WORK(max(m,n)+1:2max(m,n)).
	mn = max(m, n);
	klm = min(m - 1, kl);
	kun = min(n - 1, ku);
	kb = klm + kun;
	kb1 = kb + 1;
	inca = kb1 * ldab;
	nr = 0;
	j1 = klm + 2;
	j2 = 1 - kun;

	for (i = 0; i < minmn; i++) {
//Reduce i-th column and i-th row of matrix to bidiagonal form
	    ml = klm + 1;
	    mu = kun + 1;
	    for (kk = 0; kk <= kb; kk++) {
		j1 = j1 + kb;
		j2 = j2 + kb;
//generate plane rotations to annihilate nonzero elements
//which have been created below the band
		if (nr > 0) {
		    Rlargv(nr, &AB[klu1 + (j1 - klm - 1) * ldab], inca, &work[j1], kb1, &work[mn + j1], kb1);
		}
//apply plane rotations from the left
		for (l = 0; l < kb; l++) {
		    if (j2 - klm + l - 1 > n) {
			nrt = nr - 1;
		    } else {
			nrt = nr;
		    }
		    if (nrt > 0) {
			Rlartv(nrt, &AB[klu1 - l + (j1 - klm + l - 1) * ldab], inca, &AB[klu1 - l + 1 + (j1 - klm + l - 1) * ldab], inca, &work[mn + j1], &work[j1], kb1);
		    }
		}
		if (ml > ml0) {
		    if (ml <= m - i + 1) {
//generate plane rotation to annihilate a(i+ml-1,i)
//within the band, and apply rotation from the left
			Rlartg(AB[ku + ml - 1 + i * ldab], AB[ku + ml + i * ldab], &work[mn + i + ml - 1], &work[i + ml - 1], &ra);
			AB[ku + ml - 1 + i * ldab] = ra;
			if (i < n) {
			    Rrot(min(ku + ml - 2, n - i), &AB[ku + ml - 2 + (i + 1) * ldab], ldab - 1,
				 &AB[ku + ml - 1 + (i + 1) * ldab], ldab - 1, work[mn + i + ml - 1], work[i + ml - 1]);
			}
		    }
		    ++nr;
		    j1 = j1 - kb1;
		}
		if (wantq) {
//accumulate product of plane rotations in Q
		    for (j = j1; j <= j2; j = j + kb1) {
			Rrot(m, &q[(j - 1) * ldq + 1], 1, &q[j * ldq + 1], 1, work[mn + j], work[j]);
		    }
		}
		if (wantc) {
//apply plane rotations to C
		    for (j = j1; j <= j2; j = j + kb1) {
			Rrot(ncc, &c[j - 1 + ldc], ldc, &c[j + ldc], ldc, work[mn + j], work[j]);
		    }
		}
		if (j2 + kun > n) {
//adjust J2 to keep within the bounds of the matrix
		    --nr;
		    j2 = j2 - kb1;
		}
		for (j = j1; j <= j2; j += kb1) {
//create nonzero element a(j-1,j+ku) above the band
//and store it in WORK(n+1:2n)
		    work[j + kun] = work[j] * AB[(j + kun) * ldab + 1];
		    AB[(j + kun) * ldab + 1] = work[mn + j] * AB[(j + kun)
								 * ldab + 1];
		}
//generate plane rotations to annihilate nonzero elements
//which have been generated above the band
		if (nr > 0) {
		    Rlargv(nr, &AB[(j1 + kun - 1) * ldab + 1], inca, &work[j1 + kun], kb1, &work[mn + j1 + kun], kb1);
		}
//apply plane rotations from the right
		for (l = 0; l < kb; l++) {
		    if (j2 + l - 1 > m) {
			nrt = nr - 1;
		    } else {
			nrt = nr;
		    }
		    if (nrt > 0) {
			Rlartv(nrt, &AB[l + 1 + (j1 + kun - 1) * ldab], inca, &AB[l + (j1 + kun) * ldab], inca, &work[mn + j1 + kun], &work[j1 + kun], kb1);
		    }
		}
		if (ml == ml0 && mu > mu0) {
		    if (mu <= n - i + 1) {
//generate plane rotation to annihilate a(i,i+mu-1)
//within the band, and apply rotation from the right
			Rlartg(AB[ku - mu + 3 + (i + mu - 2) * ldab], AB[ku - mu + 2 + (i + mu - 1) * ldab], &work[mn + i + mu - 1], &work[i + mu - 1], &ra);
			AB[ku - mu + 3 + (i + mu - 2) * ldab] = ra;
			Rrot(min(kl + mu - 2, m - i), &AB[ku - mu + 4 + (i + mu - 2) * ldab], 1,
			     &AB[ku - mu + 3 + (i + mu - 1) * ldab], 1, work[mn + i + mu - 1], work[i + mu - 1]);
		    }
		    ++nr;
		    j1 = j1 - kb1;
		}
		if (wantpt) {
//accumulate product of plane rotations in P'
		    for (j = j1; j <= j2; j += kb1) {
			Rrot(n, &pt[j + kun - 1 + ldpt], ldpt, &pt[j + kun + ldpt], ldpt, work[mn + j + kun], work[j + kun]);
		    }
		}
		if (j2 + kb > m) {
//adjust J2 to keep within the bounds of the matrix
		    --nr;
		    j2 = j2 - kb1;
		}
		for (j = j1; j <= j2; j = j + kb1) {
//create nonzero element a(j+kl+ku,j+ku-1) below the
//band and store it in WORK(1:n)
		    work[j + kb] = work[j + kun] * AB[klu1 + (j + kun) * ldab];
		    AB[klu1 + (j + kun) * ldab] = work[mn + j + kun] * AB[klu1 + (j + kun) * ldab];
		}
		if (ml > ml0) {
		    ml--;
		} else {
		    mu--;
		}
	    }
	}
    }
    if (ku == 0 && kl > 0) {
//A has been reduced to lower bidiagonal form
//Transform lower bidiagonal form to upper bidiagonal by applying
//plane rotations from the left, storing diagonal elements in D
//and off-diagonal elements in E
	for (i = 0; i < min(m - 1, n); i++) {
	    Rlartg(AB[i * ldab + 1], AB[i * ldab + 2], &rc, &rs, &ra);
	    d[i] = ra;
	    if (i < n) {
		e[i] = rs * AB[(i + 1) * ldab + 1];
		AB[(i + 1) * ldab + 1] = rc * AB[(i + 1) * ldab + 1];
	    }
	    if (wantq) {
		Rrot(m, &q[i * ldq + 1], 1, &q[(i + 1) * ldq + 1], 1, rc, rs);
	    }
	    if (wantc) {
		Rrot(ncc, &c[i + ldc], ldc, &c[i + 1 + ldc], ldc, rc, rs);
	    }
	}
	if (m <= n) {
	    d[m] = AB[m * ldab + 1];
	}
    } else if (ku > 0) {
//A has been reduced to upper bidiagonal form
	if (m < n) {
//Annihilate a(m,m+1) by applying plane rotations from the
//right, storing diagonal elements in D and off-diagonal
//elements in E
	    rb = AB[ku + (m + 1) * ldab];
	    for (i = m; i >= 1; i--) {
		Rlartg(AB[ku + 1 + i * ldab], rb, &rc, &rs, &ra);
		d[i] = ra;
		if (i > 1) {
		    rb = -rs * AB[ku + i * ldab];
		    e[i - 1] = rc * AB[ku + i * ldab];
		}
		if (wantpt) {
		    Rrot(n, &pt[i + ldpt], ldpt, &pt[m + 1 + ldpt], ldpt, rc, rs);
		}
	    }
	} else {
//Copy off-diagonal elements to E and diagonal elements to D
	    for (i = 0; i < minmn - 1; i++) {
		e[i] = AB[ku + (i + 1) * ldab];
	    }
	    for (i = 0; i < minmn; i++) {
		d[i] = AB[ku + 1 + i * ldab];
	    }
	}
    } else {
//A is diagonal. Set elements of E to zero and copy diagonal
//elements to D.
	for (i = 0; i < minmn - 1; i++) {
	    e[i] = Zero;
	}
	for (i = 0; i < minmn; i++) {
	    d[i] = AB[i * ldab + 1];
	}
    }
    return;
}
