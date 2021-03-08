/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rtgsyl.cpp,v 1.5 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rtgsyl(const char *trans, INTEGER ijob, INTEGER m, INTEGER n, REAL * A, INTEGER lda, REAL * B, INTEGER ldb,
	    REAL * c, INTEGER ldc, REAL * d, INTEGER ldd, REAL * e, INTEGER lde, REAL * f, INTEGER ldf, REAL * scale, REAL * dif,
	    REAL * work, INTEGER lwork, INTEGER * iwork, INTEGER * info)
{
    INTEGER i, j, k, p, q, ie, je, mb, nb, is, js, pq;
    REAL dsum;
    INTEGER ppqq;
    INTEGER ifunc, linfo, lwmin;
    REAL scale2 = 0;
    REAL dscale, scaloc;
    INTEGER iround;
    LOGICAL notran;
    INTEGER isolve;
    LOGICAL lquery;
    REAL Zero = 0.0, One = 1.0;

//Decode and test input parameters
    *info = 0;
    notran = Mlsame(trans, "N");
    lquery = lwork == -1;

    if (!notran && !Mlsame(trans, "T")) {
	*info = -1;
    } else if (notran) {
	if (ijob < 0 || ijob > 4) {
	    *info = -2;
	}
    }
    if (*info == 0) {
	if (m <= 0) {
	    *info = -3;
	} else if (n <= 0) {
	    *info = -4;
	} else if (lda < max((INTEGER) 1, m)) {
	    *info = -6;
	} else if (ldb < max((INTEGER) 1, n)) {
	    *info = -8;
	} else if (ldc < max((INTEGER) 1, m)) {
	    *info = -10;
	} else if (ldd < max((INTEGER) 1, m)) {
	    *info = -12;
	} else if (lde < max((INTEGER) 1, n)) {
	    *info = -14;
	} else if (ldf < max((INTEGER) 1, m)) {
	    *info = -16;
	}
    }
    if (*info == 0) {
	if (notran) {
	    if (ijob == 1 || ijob == 2) {
		lwmin = max((INTEGER) 1, (m << 1) * n);
	    } else {
		lwmin = 1;
	    }
	} else {
	    lwmin = 1;
	}
	work[1] = lwmin;
	if (lwork < lwmin && !lquery) {
	    *info = -20;
	}
    }
    if (*info != 0) {
	Mxerbla("Rtgsyl", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (m == 0 || n == 0) {
	*scale = One;
	if (notran) {
	    if (ijob != 0) {
		*dif = Zero;
	    }
	}
	return;
    }
//Determine optimal block sizes MB and NB
    mb = iMlaenv(2, "Rtgsyl", trans, m, n, -1, -1);
    nb = iMlaenv(5, "Rtgsyl", trans, m, n, -1, -1);
    isolve = 1;
    ifunc = 0;
    if (notran) {
	if (ijob >= 3) {
	    ifunc = ijob - 2;
	    Rlaset("F", m, n, Zero, Zero, &c[0], ldc);
	    Rlaset("F", m, n, Zero, Zero, &f[0], ldf);
	} else if (ijob >= 1) {
	    isolve = 2;
	}
    }
    if ((mb <= 1 && nb <= 1) || (mb >= m && nb >= n)) {
	for (iround = 1; iround <= isolve; iround++) {
//Use unblocked Level 2 solver
	    dscale = Zero;
	    dsum = One;
	    pq = 0;
	    Rtgsy2(trans, ifunc, m, n, &A[0], lda, &B[0], ldb, &c[0], ldc, &d[0], ldd, &e[0], lde, &f[0], ldf, scale, &dsum, &dscale, &iwork[1], &pq, info);
	    if (dscale != Zero) {
		if (ijob == 1 || ijob == 3) {
//		    *dif = sqrt((REAL) ((m << 1) * n)) / (dscale * sqrt(dsum));
		} else {
//		    *dif = sqrt((REAL) pq) / (dscale * sqrt(dsum));
		}
	    }
	    if (isolve == 2 && iround == 1) {
		if (notran) {
		    ifunc = ijob;
		}
		scale2 = *scale;
		Rlacpy("F", m, n, &c[0], ldc, &work[0], m);
		Rlacpy("F", m, n, &f[0], ldf, &work[m * n + 1], m);
		Rlaset("F", m, n, Zero, Zero, &c[0], ldc);
		Rlaset("F", m, n, Zero, Zero, &f[0], ldf);
	    } else if (isolve == 2 && iround == 2) {
		Rlacpy("F", m, n, &work[0], m, &c[0], ldc);
		Rlacpy("F", m, n, &work[m * n + 1], m, &f[0], ldf);
		*scale = scale2;
	    }
	}
	return;
    }
//Determine block structure of A
    p = 0;
    i = 1;
  L40:
    if (i > m) {
	goto L50;
    }
    ++p;
    iwork[p] = i;
    i = i + mb;
    if (i >= m) {
	goto L50;
    }
    if (A[i + (i - 1) * lda] != Zero) {
	i++;
    }
    goto L40;
  L50:
    iwork[p + 1] = m + 1;
    if (iwork[p] == iwork[p + 1]) {
	--p;
    }
//Determine block structure of B
    q = p + 1;
    j = 0;
  L60:
    if (j > n) {
	goto L70;
    }
    ++q;
    iwork[q] = j;
    j = j + nb;
    if (j >= n) {
	goto L70;
    }
    if (B[j + (j - 1) * ldb] != Zero) {
	j++;
    }
    goto L60;
  L70:
    iwork[q + 1] = n + 1;
    if (iwork[q] == iwork[q + 1]) {
	--q;
    }
    if (notran) {
	for (iround = 1; iround <= isolve; iround++) {
//Solve (I, J)-subsystem
//    A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
//    D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
//for I = P, P - 1,..., 1; J = 1, 2,..., Q
	    dscale = Zero;
	    dsum = One;
	    pq = 0;
	    *scale = One;
	    for (j = p + 2; j <= q; j++) {
		js = iwork[j];
		je = iwork[j + 1] - 1;
		nb = je - js + 1;
		for (i = p; i >= 1; i--) {
		    is = iwork[i];
		    ie = iwork[i + 1] - 1;
		    mb = ie - is + 1;
		    ppqq = 0;
		    Rtgsy2(trans, ifunc, mb, nb, &A[is + is * lda],
			   lda, &B[js + js * ldb], ldb, &c[is + js * ldc], ldc, &d[is + is * ldd], ldd, &e[js + js * lde], lde,
			   &f[is + js * ldf], ldf, &scaloc, &dsum, &dscale, &iwork[q + 2], &ppqq, &linfo);
		    if (linfo > 0) {
			*info = linfo;
		    }
		    pq = pq + ppqq;
		    if (scaloc != One) {
			for (k = 0; k < js - 1; k++) {
			    Rscal(m, scaloc, &c[k * ldc + 1], 1);
			    Rscal(m, scaloc, &f[k * ldf + 1], 1);
			}
			for (k = js; k <= je; k++) {
			    Rscal(is - 1, scaloc, &c[k * ldc + 1], 1);
			    Rscal(is - 1, scaloc, &f[k * ldf + 1], 1);
			}
			for (k = js; k <= je; k++) {
			    Rscal(m - ie, scaloc, &c[ie + 1 + k * ldc], 1);
			    Rscal(m - ie, scaloc, &f[ie + 1 + k * ldf], 1);
			}
			for (k = je + 1; k <= n; k++) {
			    Rscal(m, scaloc, &c[k * ldc + 1], 1);
			    Rscal(m, scaloc, &f[k * ldf + 1], 1);
			}
			*scale = *scale * scaloc;
		    }
//Substitute R(I, J) and L(I, J) into remaining
//equation.
		    if (i > 1) {
			Rgemm("N", "N", is - 1, nb, mb, -One, &A[is * lda + 1], lda, &c[is + js * ldc], ldc, One, &c[js * ldc + 1], ldc);
			Rgemm("N", "N", is - 1, nb, mb, -One, &d[is * ldd + 1], ldd, &c[is + js * ldc], ldc, One, &f[js * ldf + 1], ldf);
		    }
		    if (j < q) {
			Rgemm("N", "N", mb, n - je, nb, One, &f[is + js * ldf], ldf, &B[js + (je + 1) * ldb], ldb, One, &c[is + (je + 1) * ldc], ldc);
			Rgemm("N", "N", mb, n - je, nb, One, &f[is + js * ldf], ldf, &e[js + (je + 1) * lde], lde, One, &f[is + (je + 1) * ldf], ldf);
		    }
		}
	    }
	    if (dscale != Zero) {
		if (ijob == 1 || ijob == 3) {
//		    *dif = (REAL) (sqrt(((m << 1) * n)) / ((double)dscale * sqrt(double(dsum))));
		} else {
//		    *dif = (REAL) (sqrt((double)pq) / ((double)dscale * sqrt((double)dsum)));
		}
	    }
	    if (isolve == 2 && iround == 1) {
		if (notran) {
		    ifunc = ijob;
		}
		scale2 = *scale;
		Rlacpy("F", m, n, &c[0], ldc, &work[0], m);
		Rlacpy("F", m, n, &f[0], ldf, &work[m * n + 1], m);
		Rlaset("F", m, n, Zero, Zero, &c[0], ldc);
		Rlaset("F", m, n, Zero, Zero, &f[0], ldf);
	    } else if (isolve == 2 && iround == 2) {
		Rlacpy("F", m, n, &work[0], m, &c[0], ldc);
		Rlacpy("F", m, n, &work[m * n + 1], m, &f[0], ldf);
		*scale = scale2;
	    }
	}
    } else {
//Solve transposed (I, J)-subsystem
//     A(I, I)' * R(I, J)  + D(I, I)' * L(I, J)  =  C(I, J)
//     R(I, J)  * B(J, J)' + L(I, J)  * E(J, J)' = -F(I, J)
//for I = 1,2,..., P; J = Q, Q-1,..., 1
	*scale = One;
	for (i = 0; i < p; i++) {
	    is = iwork[i];
	    ie = iwork[i + 1] - 1;
	    mb = ie - is + 1;
	    for (j = q; j >= p + 2; j--) {
		js = iwork[j];
		je = iwork[j + 1] - 1;
		nb = je - js + 1;
		Rtgsy2(trans, ifunc, mb, nb, &A[is + is * lda], lda, &B[js + js * ldb], ldb, &c[is + js * ldc], ldc,
		       &d[is + is * ldd], ldd, &e[js + js * lde], lde, &f[is + js * ldf], ldf, &scaloc, &dsum, &dscale, &iwork[q + 2], &ppqq, &linfo);
		if (linfo > 0) {
		    *info = linfo;
		}
		if (scaloc != One) {
		    for (k = 0; k < js - 1; k++) {
			Rscal(m, scaloc, &c[k * ldc + 1], 1);
			Rscal(m, scaloc, &f[k * ldf + 1], 1);
		    }
		    for (k = js; k <= je; k++) {
			Rscal(is - 1, scaloc, &c[k * ldc + 1], 1);
			Rscal(is - 1, scaloc, &f[k * ldf + 1], 1);
		    }
		    for (k = js; k <= je; k++) {
			Rscal(m - ie, scaloc, &c[ie + 1 + k * ldc], 1);
			Rscal(m - ie, scaloc, &f[ie + 1 + k * ldf], 1);
		    }
		    for (k = je + 1; k <= n; k++) {
			Rscal(m, scaloc, &c[k * ldc + 1], 1);
			Rscal(m, scaloc, &f[k * ldf + 1], 1);
		    }
		    *scale *= scaloc;
		}
//Substitute R(I, J) and L(I, J) into remaining equation.
		if (j > p + 2) {
		    Rgemm("N", "T", mb, js - 1, nb, One, &c[is + js * ldc], ldc, &B[js * ldb + 1], ldb, One, &f[is + ldf], ldf);
		    Rgemm("N", "T", mb, js - 1, nb, One, &f[is + js * ldf], ldf, &e[js * lde + 1], lde, One, &f[is + ldf], ldf);
		}
		if (i < p) {
		    Rgemm("T", "N", m - ie, nb, mb, -One, &A[is + (ie + 1) * lda], lda, &c[is + js * ldc], ldc, One, &c[ie + 1 + js * ldc], ldc);
		    Rgemm("T", "N", m - ie, nb, mb, -One, &d[is + (ie + 1) * ldd], ldd, &f[is + js * ldf], ldf, One, &c[ie + 1 + js * ldc], ldc);
		}
	    }
	}
    }
    work[1] = lwmin;
    return;
}
