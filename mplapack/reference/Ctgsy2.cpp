/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Ctgsy2.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Ctgsy2(const char *trans, INTEGER ijob, INTEGER m, INTEGER n, COMPLEX * A, INTEGER lda,
       COMPLEX * B, INTEGER ldb, COMPLEX * C, INTEGER ldc, COMPLEX * d, INTEGER ldd,
       COMPLEX * e, INTEGER lde, COMPLEX * f, INTEGER ldf, REAL * scale, REAL * rdsum, REAL * rdscal, INTEGER * info)
{
    INTEGER i, j, k;
    COMPLEX z[4], rhs[2];
    INTEGER ierr, ipiv[2], jpiv[2];
    COMPLEX alpha;
    REAL scaloc;
    INTEGER notran;
    REAL One = 1.0;

    *info = 0;
    ierr = 0;
    notran = Mlsame(trans, "N");
    if (!notran && !Mlsame(trans, "C")) {
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
	Mxerbla("Ctgsy2", -(*info));
	return;
    }
    if (notran) {
//Solve (I, J) - system
//   A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
//   D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
//for I = M, M - 1, ..., 1; J = 1, 2, ..., N
	*scale = One;
	scaloc = One;
	for (j = 0; j < n; j++) {
	    for (i = m; i >= 1; i--) {
//Build 2 by 2 system
		z[0] = A[i + i * lda];
		z[1] = d[i + i * ldd];
		z[2] = -B[j + j * ldb];
		z[3] = -e[j + j * lde];
//Set up right hand side(s)
		rhs[0] = C[i + j * ldc];
		rhs[1] = f[i + j * ldf];
//Solve Z * x = RHS
		Cgetc2(2, z, 2, ipiv, jpiv, &ierr);
		if (ierr > 0) {
		    *info = ierr;
		}
		if (ijob == 0) {
		    Cgesc2(2, z, 2, rhs, ipiv, jpiv, &scaloc);
		    if (scaloc != One) {
			for (k = 0; k < n; k++) {
			    Cscal(m, (COMPLEX) scaloc, &C[k * ldc + 1], 1);
			    Cscal(m, (COMPLEX) scaloc, &f[k * ldf + 1], 1);
			}
			*scale = *scale * scaloc;
		    }
		} else {
		    Clatdf(ijob, 2, z, 2, rhs, rdsum, rdscal, ipiv, jpiv);
		}
//Unpack solution vector(s)
		C[i + j * ldc] = rhs[0];
		f[i + j * ldf] = rhs[1];
//Substitute R(I, J) and L(I, J) into remaining equation.
		if (i > 1) {
		    alpha = -rhs[0];
		    Caxpy(i - 1, alpha, &A[i * lda], 1, &C[j * ldc + 1], 1);
		    Caxpy(i - 1, alpha, &d[i * ldd + 1], 1, &f[j * ldf + 1], 1);
		}
		if (j < n) {
		    Caxpy(n - j, rhs[1], &B[j + (j + 1) * ldb], ldb, &C[i + (j + 1) * ldc], ldc);
		    Caxpy(n - j, rhs[1], &e[j + (j + 1) * lde], lde, &f[i + (j + 1) * ldf], ldf);
		}
	    }
	}
    } else {
//Solve transposed (I, J) - system:
//   A(I, I)' * R(I, J) + D(I, I)' * L(J, J) = C(I, J)
//   R(I, I) * B(J, J) + L(I, J) * E(J, J)   = -F(I, J)
//for I = 1, 2, ..., M, J = N, N - 1, ..., 1
	*scale = One;
	scaloc = One;
	for (i = 0; i < m; i++) {
	    for (j = n; j >= 1; j--) {
//Build 2 by 2 system Z'
		z[0] = conj(A[i + i * lda]);
		z[1] = conj(B[j + j * ldb]);
		z[2] = conj(d[i + i * ldd]);
		z[3] = conj(e[j + j * lde]);
//Set up right hand side(s)
		rhs[0] = C[i + j * ldc];
		rhs[1] = f[i + j * ldf];
//Solve Z' * x = RHS
		Cgetc2(2, z, 2, ipiv, jpiv, &ierr);
		if (ierr > 0) {
		    *info = ierr;
		}
		Cgesc2(2, z, 2, rhs, ipiv, jpiv, &scaloc);
		if (scaloc != One) {
		    for (k = 0; k < n; k++) {
			Cscal(m, (COMPLEX) scaloc, &C[k * ldc + 1], 1);
			Cscal(m, (COMPLEX) scaloc, &f[k * ldf + 1], 1);
		    }
		    *scale = *scale * scaloc;
		}
//Unpack solution vector(s)
		C[i + j * ldc] = rhs[0];
		f[i + j * ldf] = rhs[1];
//Substitute R(I, J) and L(I, J) into remaining equation.
		for (k = 0; k < j - 1; k++) {
		    f[i + k * ldf] = f[i + k * ldf] + rhs[0] * conj(B[k + j * ldb]) + rhs[1] * conj(e[k + j * lde]);
		}
		for (k = i + 1; k <= m; k++) {
		    C[k + j * ldc] = C[k + j * ldc] - conj(A[i + k * lda]) * rhs[0] - conj(d[i + k * ldd]) * rhs[1];
		}
	    }
	}
    }
    return;
}
