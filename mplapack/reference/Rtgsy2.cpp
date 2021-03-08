/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rtgsy2.cpp,v 1.8 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void
Rtgsy2(const char *trans, INTEGER ijob, INTEGER m, INTEGER n, REAL * A, INTEGER lda,
       REAL * B, INTEGER ldb, REAL * c, INTEGER ldc, REAL * d, INTEGER ldd,
       REAL * e, INTEGER lde, REAL * f, INTEGER ldf, REAL * scale, REAL * rdsum, REAL * rdscal, INTEGER * iwork, INTEGER * pq, INTEGER * info)
{
    INTEGER i, j, k, p, q;
    REAL z[64] /* was [8][8] */ ;
    INTEGER ie, je, mb, nb, ii, jj, is, js;
    REAL rhs[8];
    INTEGER isp1, jsp1;
    INTEGER ierr, zdim, ipiv[8], jpiv[8];
    REAL alpha;
    REAL scaloc;
    INTEGER notran;
    REAL Zero = 0.0, One = 1.0;

    *info = 0;
    ierr = 0;
    notran = Mlsame(trans, "N");
    if (!notran && !Mlsame(trans, "T")) {
	*info = -1;
    } else if (notran) {
	if (ijob < 0 || ijob > 2) {
	    *info = -2;
	}
    }
    if (*info == 0) {
	if (m <= 0) {
	    *info = -3;
	} else if (n <= 0) {
	    *info = -4;
	} else if (lda < max((INTEGER) 1, m)) {
	    *info = -5;
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
    if (*info != 0) {
	Mxerbla("Rtgsy2", -(*info));
	return;
    }
//Determine block structure of A
    *pq = 0;
    p = 0;
    i = 1;
  L10:
    if (i > m) {
	goto L20;
    }
    ++p;
    iwork[p] = i;
    if (i == m) {
	goto L20;
    }
    if (A[i + 1 + i * lda] != Zero) {
	i += 2;
    } else {
	i++;
    }
    goto L10;
  L20:
    iwork[p + 1] = m + 1;
//Determine block structure of B
    q = p + 1;
    j = 0;
  L30:
    if (j > n) {
	goto L40;
    }
    ++q;
    iwork[q] = j;
    if (j == n) {
	goto L40;
    }
    if (B[j + 1 + j * ldb] != Zero) {
	j += 2;
    } else {
	j++;
    }
    goto L30;
  L40:
    iwork[q + 1] = n + 1;
    *pq = p * (q - p - 1);

    if (notran) {
//Solve (I, J) - subsystem
//   A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
//   D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
//for I = P, P - 1, ..., 1; J = 1, 2, ..., Q
	*scale = One;
	scaloc = One;
	for (j = p + 2; j <= q; j++) {
	    js = iwork[j];
	    jsp1 = js + 1;
	    je = iwork[j + 1] - 1;
	    nb = je - js + 1;
	    for (i = p; i >= 1; i--) {
		is = iwork[i];
		isp1 = is + 1;
		ie = iwork[i + 1] - 1;
		mb = ie - is + 1;
		zdim = mb * nb << 1;
		if (mb == 1 && nb == 1) {
//Build a 2-by-2 system Z * x = RHS
		    z[0] = A[is + is * lda];
		    z[1] = d[is + is * ldd];
		    z[8] = -B[js + js * ldb];
		    z[9] = -e[js + js * lde];
//Set up right hand side(s)
		    rhs[0] = c[is + js * ldc];
		    rhs[1] = f[is + js * ldf];
//Solve Z * x = RHS
		    Rgetc2(zdim, z, 8, ipiv, jpiv, &ierr);
		    if (ierr > 0) {
			*info = ierr;
		    }
		    if (ijob == 0) {
			Rgesc2(zdim, z, 8, rhs, ipiv, jpiv, &scaloc);
			if (scaloc != One) {
			    for (k = 0; k < n; k++) {
				Rscal(m, scaloc, &c[k * ldc + 1], 1);
				Rscal(m, scaloc, &f[k * ldf + 1], 1);
			    }
			    *scale = *scale * scaloc;
			}
		    } else {
			Rlatdf(ijob, zdim, z, 8, rhs, rdsum, rdscal, ipiv, jpiv);
		    }
//Unpack solution vector(s)
		    c[is + js * ldc] = rhs[0];
		    f[is + js * ldf] = rhs[1];
//Substitute R(I, J) and L(I, J) into remaining
//equation.
		    if (i > 1) {
			alpha = -rhs[0];
			Raxpy(is - 1, alpha, &A[is * lda], 1, &c[js * ldc + 1], 1);
			Raxpy(is - 1, alpha, &d[is * ldd + 1], 1, &f[js * ldf + 1], 1);
		    }
		    if (j < q) {
			Raxpy(n - je, rhs[1], &B[js + (je + 1) * ldb], ldb, &c[is + (je + 1) * ldc], ldc);
			Raxpy(n - je, rhs[1], &e[js + (je + 1) * lde], lde, &f[is + (je + 1) * ldf], ldf);
		    }
		} else if (mb == 1 && nb == 2) {
//Build a 4-by-4 system Z * x = RHS
		    z[0] = A[is + is * lda];
		    z[1] = Zero;
		    z[2] = d[is + is * ldd];
		    z[3] = Zero;
		    z[8] = Zero;
		    z[9] = A[is + is * lda];
		    z[10] = Zero;
		    z[11] = d[is + is * ldd];
		    z[16] = -B[js + js * ldb];
		    z[17] = -B[js + jsp1 * ldb];
		    z[18] = -e[js + js * lde];
		    z[19] = -e[js + jsp1 * lde];
		    z[24] = -B[jsp1 + js * ldb];
		    z[25] = -B[jsp1 + jsp1 * ldb];
		    z[26] = Zero;
		    z[27] = -e[jsp1 + jsp1 * lde];
//Set up right hand side(s)
		    rhs[0] = c[is + js * ldc];
		    rhs[1] = c[is + jsp1 * ldc];
		    rhs[2] = f[is + js * ldf];
		    rhs[3] = f[is + jsp1 * ldf];
//Solve Z * x = RHS
		    Rgetc2(zdim, z, 8, ipiv, jpiv, &ierr);
		    if (ierr > 0) {
			*info = ierr;
		    }
		    if (ijob == 0) {
			Rgesc2(zdim, z, 8, rhs, ipiv, jpiv, &scaloc);
			if (scaloc != One) {
			    for (k = 0; k < n; k++) {
				Rscal(m, scaloc, &c[k * ldc + 1], 1);
				Rscal(m, scaloc, &f[k * ldf + 1], 1);
			    }
			    *scale = *scale * scaloc;
			}
		    } else {
			Rlatdf(ijob, zdim, z, 8, rhs, rdsum, rdscal, ipiv, jpiv);
		    }
//Unpack solution vector(s)
		    c[is + js * ldc] = rhs[0];
		    c[is + jsp1 * ldc] = rhs[1];
		    f[is + js * ldf] = rhs[2];
		    f[is + jsp1 * ldf] = rhs[3];
//Substitute R(I, J) and L(I, J) into remaining
//equation.
		    if (i > 1) {
			Rger(is - 1, nb, -One, &A[is * lda], 1, rhs, 1, &c[js * ldc + 1], ldc);
			Rger(is - 1, nb, -One, &d[is * ldd + 1], 1, rhs, 1, &f[js * ldf + 1], ldf);
		    }
		    if (j < q) {
			Raxpy(n - je, rhs[2], &B[js + (je + 1) * ldb], ldb, &c[is + (je + 1) * ldc], ldc);
			Raxpy(n - je, rhs[2], &e[js + (je + 1) * lde], lde, &f[is + (je + 1) * ldf], ldf);
			Raxpy(n - je, rhs[3], &B[jsp1 + (je + 1) * ldb], ldb, &c[is + (je + 1) * ldc], ldc);
			Raxpy(n - je, rhs[3], &e[jsp1 + (je + 1) * lde], lde, &f[is + (je + 1) * ldf], ldf);
		    }
		} else if (mb == 2 && nb == 1) {
//Build a 4-by-4 system Z * x = RHS
		    z[0] = A[is + is * lda];
		    z[1] = A[isp1 + is * lda];
		    z[2] = d[is + is * ldd];
		    z[3] = Zero;

		    z[8] = A[is + isp1 * lda];
		    z[9] = A[isp1 + isp1 * lda];
		    z[10] = d[is + isp1 * ldd];
		    z[11] = d[isp1 + isp1 * ldd];

		    z[16] = -B[js + js * ldb];
		    z[17] = Zero;
		    z[18] = -e[js + js * lde];
		    z[19] = Zero;

		    z[24] = Zero;
		    z[25] = -B[js + js * ldb];
		    z[26] = Zero;
		    z[27] = -e[js + js * lde];
//Set up right hand side(s)
		    rhs[0] = c[is + js * ldc];
		    rhs[1] = c[isp1 + js * ldc];
		    rhs[2] = f[is + js * ldf];
		    rhs[3] = f[isp1 + js * ldf];
//Solve Z * x = RHS
		    Rgetc2(zdim, z, 8, ipiv, jpiv, &ierr);
		    if (ierr > 0) {
			*info = ierr;
		    }
		    if (ijob == 0) {
			Rgesc2(zdim, z, 8, rhs, ipiv, jpiv, &scaloc);
			if (scaloc != One) {
			    for (k = 0; k < n; k++) {
				Rscal(m, scaloc, &c[k * ldc + 1], 1);
				Rscal(m, scaloc, &f[k * ldf + 1], 1);
			    }
			    *scale = *scale * scaloc;
			}
		    } else {
			Rlatdf(ijob, zdim, z, 8, rhs, rdsum, rdscal, ipiv, jpiv);
		    }
//Unpack solution vector(s)
		    c[is + js * ldc] = rhs[0];
		    c[isp1 + js * ldc] = rhs[1];
		    f[is + js * ldf] = rhs[2];
		    f[isp1 + js * ldf] = rhs[3];
//Substitute R(I, J) and L(I, J) into remaining
//equation.
		    if (i > 1) {
			Rgemv("N", is - 1, mb, -One, &A[is * lda], lda, rhs, 1, One, &c[js * ldc + 1], 1);
			Rgemv("N", is - 1, mb, -One, &d[is * ldd + 1], ldd, rhs, 1, One, &f[js * ldf + 1], 1);
		    }
		    if (j < q) {
			Rger(mb, n - je, One, &rhs[2], 1, &B[js + (je + 1) * ldb], ldb, &c[is + (je + 1) * ldc], ldc);
			Rger(mb, n - je, One, &rhs[2], 1, &e[js + (je + 1) * lde], lde, &f[is + (je + 1) * ldf], ldf);
		    }
		} else if (mb == 2 && nb == 2) {
//Build an 8-by-8 system Z * x = RHS
		    Rlaset("F", 8, 8, Zero, Zero, z, 8);
		    z[0] = A[is + is * lda];
		    z[1] = A[isp1 + is * lda];
		    z[4] = d[is + is * ldd];
		    z[8] = A[is + isp1 * lda];
		    z[9] = A[isp1 + isp1 * lda];
		    z[12] = d[is + isp1 * ldd];
		    z[13] = d[isp1 + isp1 * ldd];
		    z[18] = A[is + is * lda];
		    z[19] = A[isp1 + is * lda];
		    z[22] = d[is + is * ldd];
		    z[26] = A[is + isp1 * lda];
		    z[27] = A[isp1 + isp1 * lda];
		    z[30] = d[is + isp1 * ldd];
		    z[31] = d[isp1 + isp1 * ldd];
		    z[32] = -B[js + js * ldb];
		    z[34] = -B[js + jsp1 * ldb];
		    z[36] = -e[js + js * lde];
		    z[38] = -e[js + jsp1 * lde];
		    z[41] = -B[js + js * ldb];
		    z[43] = -B[js + jsp1 * ldb];
		    z[45] = -e[js + js * lde];
		    z[47] = -e[js + jsp1 * lde];
		    z[48] = -B[jsp1 + js * ldb];
		    z[50] = -B[jsp1 + jsp1 * ldb];
		    z[54] = -e[jsp1 + jsp1 * lde];
		    z[57] = -B[jsp1 + js * ldb];
		    z[59] = -B[jsp1 + jsp1 * ldb];
		    z[63] = -e[jsp1 + jsp1 * lde];
//Set up right hand side(s)
		    k = 0;
		    ii = mb * nb + 1;
		    for (jj = 0; jj <= nb - 1; jj++) {
			Rcopy(mb, &c[is + (js + jj) * ldc], 1, &rhs[k - 1], 1);
			Rcopy(mb, &f[is + (js + jj) * ldf], 1, &rhs[ii - 1], 1);
			k = k + mb;
			ii = ii + mb;
		    }
//Solve Z * x = RHS
		    Rgetc2(zdim, z, 8, ipiv, jpiv, &ierr);
		    if (ierr > 0) {
			*info = ierr;
		    }
		    if (ijob == 0) {
			Rgesc2(zdim, z, 8, rhs, ipiv, jpiv, &scaloc);
			if (scaloc != One) {
			    for (k = 0; k < n; k++) {
				Rscal(m, scaloc, &c[k * ldc + 1], 1);
				Rscal(m, scaloc, &f[k * ldf + 1], 1);
			    }
			    *scale *= scaloc;
			}
		    } else {
			Rlatdf(ijob, zdim, z, 8, rhs, rdsum, rdscal, ipiv, jpiv);
		    }
//Unpack solution vector(s)
		    k = 0;
		    ii = mb * nb + 1;
		    for (jj = 0; jj <= nb - 1; jj++) {
			Rcopy(mb, &rhs[k - 1], 1, &c[is + (js + jj) * ldc], 1);
			Rcopy(mb, &rhs[ii - 1], 1, &f[is + (js + jj) * ldf], 1);
			k = k + mb;
			ii = ii + mb;
		    }
//Substitute R(I, J) and L(I, J) into remaining
//equation.
		    if (i > 1) {
			Rgemm("N", "N", is - 1, nb, mb, -One, &A[is * lda], lda, rhs, mb, One, &c[js * ldc + 1], ldc);
			Rgemm("N", "N", is - 1, nb, mb, -One, &d[is * ldd + 1], ldd, rhs, mb, One, &f[js * ldf + 1], ldf);
		    }
		    if (j < q) {
			k = mb * nb + 1;
			Rgemm("N", "N", mb, n - je, nb, One, &rhs[k - 1], mb, &B[js + (je + 1) * ldb], ldb, One, &c[is + (je + 1) * ldc], ldc);
			Rgemm("N", "N", mb, n - je, nb, One, &rhs[k - 1], mb, &e[js + (je + 1) * lde], lde, One, &f[is + (je + 1) * ldf], ldf);
		    }
		}
	    }
	}
    } else {
//Solve (I, J) - subsystem 
//     A(I, I)' * R(I, J) + D(I, I)' * L(J, J)  =  C(I, J)
//     R(I, I)  * B(J, J) + L(I, J)  * E(J, J)  = -F(I, J)
//for I = 1, 2, ..., P, J = Q, Q - 1, ..., 1
	*scale = One;
	scaloc = One;
	for (i = 0; i < p; i++) {
	    is = iwork[i];
	    isp1 = is + 1;
	    ie = i;
	    mb = ie - is + 1;
	    for (j = q; j >= p + 2; j--) {
		js = iwork[j];
		jsp1 = js + 1;
		je = iwork[j + 1] - 1;
		nb = je - js + 1;
		zdim = mb * nb << 1;
		if (mb == 1 && nb == 1) {
//Build a 2-by-2 system Z' * x = RHS
		    z[0] = A[is + is * lda];
		    z[1] = -B[js + js * ldb];
		    z[8] = d[is + is * ldd];
		    z[9] = -e[js + js * lde];
//Set up right hand side(s)
		    rhs[0] = c[is + js * ldc];
		    rhs[1] = f[is + js * ldf];
//Solve Z' * x = RHS
		    Rgetc2(zdim, z, 8, ipiv, jpiv, &ierr);
		    if (ierr > 0) {
			*info = ierr;
		    }
		    Rgesc2(zdim, z, 8, rhs, ipiv, jpiv, &scaloc);
		    if (scaloc != One) {
			for (k = 0; k < n; k++) {
			    Rscal(m, scaloc, &c[k * ldc + 1], 1);
			    Rscal(m, scaloc, &f[k * ldf + 1], 1);
			}
			*scale *= scaloc;
		    }
//Unpack solution vector(s)
		    c[is + js * ldc] = rhs[0];
		    f[is + js * ldf] = rhs[1];
//Substitute R(I, J) and L(I, J) into remaining
//equation.
		    if (j > p + 2) {
			alpha = rhs[0];
			Raxpy(js - 1, alpha, &B[js * ldb + 1], 1, &f[is + ldf], ldf);
			alpha = rhs[1];
			Raxpy(js - 1, alpha, &e[js * lde + 1], 1, &f[is + ldf], ldf);
		    }
		    if (i < p) {
			alpha = -rhs[0];
			Raxpy(m - ie, alpha, &A[is + (ie + 1) * lda], lda, &c[ie + 1 + js * ldc], 1);
			alpha = -rhs[1];
			Raxpy(m - ie, alpha, &d[is + (ie + 1) * ldd], ldd, &c[ie + 1 + js * ldc], 1);
		    }

		} else if (mb == 1 && nb == 2) {
//Build a 4-by-4 system Z' * x = RHS
		    z[0] = A[is + is * lda];
		    z[1] = Zero;
		    z[2] = -B[js + js * ldb];
		    z[3] = -B[jsp1 + js * ldb];

		    z[8] = Zero;
		    z[9] = A[is + is * lda];
		    z[10] = -B[js + jsp1 * ldb];
		    z[11] = -B[jsp1 + jsp1 * ldb];

		    z[16] = d[is + is * ldd];
		    z[17] = Zero;
		    z[18] = -e[js + js * lde];
		    z[19] = Zero;

		    z[24] = Zero;
		    z[25] = d[is + is * ldd];
		    z[26] = -e[js + jsp1 * lde];
		    z[27] = -e[jsp1 + jsp1 * lde];
//Set up right hand side(s)
		    rhs[0] = c[is + js * ldc];
		    rhs[1] = c[is + jsp1 * ldc];
		    rhs[2] = f[is + js * ldf];
		    rhs[3] = f[is + jsp1 * ldf];
//Solve Z' * x = RHS
		    Rgetc2(zdim, z, 8, ipiv, jpiv, &ierr);
		    if (ierr > 0) {
			*info = ierr;
		    }
		    Rgesc2(zdim, z, 8, rhs, ipiv, jpiv, &scaloc);
		    if (scaloc != One) {
			for (k = 0; k < n; k++) {
			    Rscal(m, scaloc, &c[k * ldc + 1], 1);
			    Rscal(m, scaloc, &f[k * ldf + 1], 1);

			}
			*scale *= scaloc;
		    }
//Unpack solution vector(s)
		    c[is + js * ldc] = rhs[0];
		    c[is + jsp1 * ldc] = rhs[1];
		    f[is + js * ldf] = rhs[2];
		    f[is + jsp1 * ldf] = rhs[3];
//Substitute R(I, J) and L(I, J) into remaining
//equation.
		    if (j > p + 2) {
			Raxpy(js - 1, rhs[0], &B[js * ldb + 1], 1, &f[is + ldf], ldf);
			Raxpy(js - 1, rhs[1], &B[jsp1 * ldb + 1], 1, &f[is + ldf], ldf);
			Raxpy(js - 1, rhs[2], &e[js * lde + 1], 1, &f[is + ldf], ldf);
			Raxpy(js - 1, rhs[3], &e[jsp1 * lde + 1], 1, &f[is + ldf], ldf);
		    }
		    if (i < p) {
			Rger(m - ie, nb, -One, &A[is + (ie + 1) * lda], lda, rhs, 1, &c[ie + 1 + js * ldc], ldc);
			Rger(m - ie, nb, -One, &d[is + (ie + 1) * ldd], ldd, &rhs[2], 1, &c[ie + 1 + js * ldc], ldc);
		    }
		} else if (mb == 2 && nb == 1) {
//Build a 4-by-4 system Z' * x = RHS
		    z[0] = A[is + is * lda];
		    z[1] = A[is + isp1 * lda];
		    z[2] = -B[js + js * ldb];
		    z[3] = Zero;
		    z[8] = A[isp1 + is * lda];
		    z[9] = A[isp1 + isp1 * lda];
		    z[10] = Zero;
		    z[11] = -B[js + js * ldb];
		    z[16] = d[is + is * ldd];
		    z[17] = d[is + isp1 * ldd];
		    z[18] = -e[js + js * lde];
		    z[19] = Zero;
		    z[24] = Zero;
		    z[25] = d[isp1 + isp1 * ldd];
		    z[26] = Zero;
		    z[27] = -e[js + js * lde];
//Set up right hand side(s)
		    rhs[0] = c[is + js * ldc];
		    rhs[1] = c[isp1 + js * ldc];
		    rhs[2] = f[is + js * ldf];
		    rhs[3] = f[isp1 + js * ldf];
//Solve Z' * x = RHS
		    Rgetc2(zdim, z, 8, ipiv, jpiv, &ierr);
		    if (ierr > 0) {
			*info = ierr;
		    }
		    Rgesc2(zdim, z, 8, rhs, ipiv, jpiv, &scaloc);
		    if (scaloc != One) {
			for (k = 0; k < n; k++) {
			    Rscal(m, scaloc, &c[k * ldc + 1], 1);
			    Rscal(m, scaloc, &f[k * ldf + 1], 1);
			}
			*scale = *scale * scaloc;
		    }
//Unpack solution vector(s)
		    c[is + js * ldc] = rhs[0];
		    c[isp1 + js * ldc] = rhs[1];
		    f[is + js * ldf] = rhs[2];
		    f[isp1 + js * ldf] = rhs[3];
//Substitute R(I, J) and L(I, J) into remaining
//equation.
		    if (j > p + 2) {
			Rger(mb, js - 1, One, rhs, 1, &B[js * ldb + 1], 1, &f[is + ldf], ldf);
			Rger(mb, js - 1, One, &rhs[2], 1, &e[js * lde + 1], 1, &f[is + ldf], ldf);
		    }
		    if (i < p) {
			Rgemv("T", mb, m - ie, -One, &A[is + (ie + 1) * lda], lda, rhs, 1, One, &c[ie + 1 + js * ldc], 1);
			Rgemv("T", mb, m - ie, -One, &d[is + (ie + 1) * ldd], ldd, &rhs[2], 1, One, &c[ie + 1 + js * ldc], 1);
		    }
		} else if (mb == 2 && nb == 2) {
//Build an 8-by-8 system Z' * x = RHS
		    Rlaset("F", 8, 8, Zero, Zero, z, 8);
		    z[0] = A[is + is * lda];
		    z[1] = A[is + isp1 * lda];
		    z[4] = -B[js + js * ldb];
		    z[6] = -B[jsp1 + js * ldb];

		    z[8] = A[isp1 + is * lda];
		    z[9] = A[isp1 + isp1 * lda];
		    z[13] = -B[js + js * ldb];
		    z[15] = -B[jsp1 + js * ldb];

		    z[18] = A[is + is * lda];
		    z[19] = A[is + isp1 * lda];
		    z[20] = -B[js + jsp1 * ldb];
		    z[22] = -B[jsp1 + jsp1 * ldb];

		    z[26] = A[isp1 + is * lda];
		    z[27] = A[isp1 + isp1 * lda];
		    z[29] = -B[js + jsp1 * ldb];
		    z[31] = -B[jsp1 + jsp1 * ldb];

		    z[32] = d[is + is * ldd];
		    z[33] = d[is + isp1 * ldd];
		    z[36] = -e[js + js * lde];

		    z[41] = d[isp1 + isp1 * ldd];
		    z[45] = -e[js + js * lde];

		    z[50] = d[is + is * ldd];
		    z[51] = d[is + isp1 * ldd];
		    z[52] = -e[js + jsp1 * lde];
		    z[54] = -e[jsp1 + jsp1 * lde];

		    z[59] = d[isp1 + isp1 * ldd];
		    z[61] = -e[js + jsp1 * lde];
		    z[63] = -e[jsp1 + jsp1 * lde];
//Set up right hand side(s)
		    k = 0;
		    ii = mb * nb + 1;
		    for (jj = 0; jj <= nb - 1; jj++) {
			Rcopy(mb, &c[is + (js + jj) * ldc], 1, &rhs[k - 1], 1);
			Rcopy(mb, &f[is + (js + jj) * ldf], 1, &rhs[ii - 1], 1);
			k = k + mb;
			ii = ii + mb;
		    }
//Solve Z' * x = RHS
		    Rgetc2(zdim, z, 8, ipiv, jpiv, &ierr);
		    if (ierr > 0) {
			*info = ierr;
		    }
		    Rgesc2(zdim, z, 8, rhs, ipiv, jpiv, &scaloc);
		    if (scaloc != One) {
			for (k = 0; k < n; k++) {
			    Rscal(m, scaloc, &c[k * ldc + 1], 1);
			    Rscal(m, scaloc, &f[k * ldf + 1], 1);
			}
			*scale = *scale * scaloc;
		    }
//Unpack solution vector(s)
		    k = 0;
		    ii = mb * nb + 1;
		    for (jj = 0; jj <= nb - 1; jj++) {
			Rcopy(mb, &rhs[k - 1], 1, &c[is + (js + jj) * ldc], 1);
			Rcopy(mb, &rhs[ii - 1], 1, &f[is + (js + jj) * ldf], 1);
			k = k + mb;
			ii = ii + mb;
		    }
//Substitute R(I, J) and L(I, J) into remaining
//equation.
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
    }
    return;
}
