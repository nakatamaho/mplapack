/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Ctrsyl.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Ctrsyl(const char *trana, const char *tranb, INTEGER isgn, INTEGER
	    m, INTEGER n, COMPLEX * A, INTEGER lda, COMPLEX * B, INTEGER ldb, COMPLEX * c, INTEGER ldc, REAL * scale, INTEGER * info)
{
    INTEGER j, k, l;
    COMPLEX a11;
    REAL db;
    COMPLEX x11;
    REAL da11;
    COMPLEX vec;
    REAL dum[1], eps, sgn, smin;
    COMPLEX suml, sumr;
    REAL scaloc;
    REAL bignum;
    INTEGER notrna, notrnb;
    REAL smlnum;
    REAL One = 1.0;
    REAL mtemp1, mtemp2, mtemp3;

//Decode and Test input parameters
    notrna = Mlsame(trana, "N");
    notrnb = Mlsame(tranb, "N");
    *info = 0;
    if (!notrna && !Mlsame(trana, "C")) {
	*info = -1;
    } else if (!notrnb && !Mlsame(tranb, "C")) {
	*info = -2;
    } else if (isgn != 1 && isgn != -1) {
	*info = -3;
    } else if (m < 0) {
	*info = -4;
    } else if (n < 0) {
	*info = -5;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -7;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -9;
    } else if (ldc < max((INTEGER) 1, m)) {
	*info = -11;
    }
    if (*info != 0) {
	Mxerbla("Ctrsyl", -(*info));
	return;
    }
//Quick return if possible
    if (m == 0 || n == 0) {
	return;
    }
//Set constants to control overflow
    eps = Rlamch("P");
    smlnum = Rlamch("S");
    bignum = One / smlnum;
    smlnum = smlnum * (m * n) / eps;
    bignum = One / smlnum;
    mtemp1 = smlnum, mtemp2 = eps * Clange("M", m, m, &A[0], lda, dum), mtemp3 = max(mtemp1, mtemp2), mtemp1 = eps * Clange("M", n, n, &B[0], ldb, dum);
    smin = max(mtemp3, mtemp1);
    *scale = One;
    sgn = (isgn);
    if (notrna && notrnb) {
//Solve    A*X + ISGN*X*B = scale*C.
//The (K,L)th block of X is determined starting from
//bottom-left corner column by column by
//    A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
//Where
//            M                        L-1
//  R(K,L) = SUM [A(K,I)*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)].
//          I=K+1                      J=1
	for (l = 0; l < n; l++) {
	    for (k = m; k >= 1; k--) {
		suml = Cdotu(m - k, &A[k + min(k + 1, m) * lda], lda, &c[min(k + 1, m) + l * ldc], 1);
		sumr = Cdotu(l - 1, &c[k + ldc], ldc, &B[l * ldb + 1], 1);
		vec = c[k + l * ldc] - (suml + sgn * sumr);
		scaloc = One;
		a11 = A[k + k * lda] + sgn * B[l + l * ldb];
		da11 = abs(a11.real()) + abs(a11.imag());
		if (da11 <= smin) {
		    a11 = smin;
		    da11 = smin;
		    *info = 1;
		}
		db = abs(vec.real()) + abs(vec.imag());
		if (da11 < One && db > One) {
		    if (db > bignum * da11) {
			scaloc = One / db;
		    }
		}
		x11 = Cladiv(vec * scaloc, a11);
		if (scaloc != One) {
		    for (j = 0; j < n; j++) {
			CRscal(m, scaloc, &c[j * ldc + 1], 1);
		    }
		    (*scale) = (*scale) * scaloc;
		}
		c[k + l * ldc] = x11;
	    }
	}
    } else if (!notrna && notrnb) {
//Solve    A' *X + ISGN*X*B = scale*C.
//The (K,L)th block of X is determined starting from
//upper-left corner column by column by
//    A'(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
//Where
//           K-1                         L-1
//  R(K,L) = SUM [A'(I,K)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)]
//           I=1                         J=1
	for (l = 0; l < n; l++) {
	    for (k = 0; k < m; k++) {
		suml = Cdotc(k - 1, &A[k * lda], 1, &c[l * ldc + 1], 1);
		sumr = Cdotu(l - 1, &c[k + ldc], ldc, &B[l * ldb + 1], 1);
		vec = c[k + l * ldc] - (suml + sgn * sumr);
		scaloc = One;
		a11 = conj(A[k + k * lda]) + sgn * B[l + l * ldb];
		da11 = abs(a11.real()) + abs(a11.imag());
		if (da11 <= smin) {
		    a11 = smin;
		    da11 = smin;
		    *info = 1;
		}
		db = abs(vec.real()) + abs(vec.imag());
		if (da11 < One && db > One) {
		    if (db > bignum * da11) {
			scaloc = One / db;
		    }
		}
		x11 = Cladiv(vec * scaloc, a11);
		if (scaloc != One) {
		    for (j = 0; j < n; j++) {
			CRscal(m, scaloc, &c[j * ldc + 1], 1);
		    }
		    (*scale) = (*scale) * scaloc;
		}
		c[k + l * ldc] = x11;
	    }
	}
    } else if (!notrna && !notrnb) {
//Solve    A'*X + ISGN*X*B' = C.
//The (K,L)th block of X is determined starting from
//upper-right corner column by column by
//    A'(K,K)*X(K,L) + ISGN*X(K,L)*B'(L,L) = C(K,L) - R(K,L)
//Where
//            K-1
//   R(K,L) = SUM [A'(I,K)*X(I,L)] +
//            I=1
//                   N
//             ISGN*SUM [X(K,J)*B'(L,J)].
//                  J=L+1
	for (l = n; l >= 1; l--) {
	    for (k = 0; k < m; k++) {
		suml = Cdotc(k - 1, &A[k * lda], 1, &c[l * ldc + 1], 1);
		sumr = Cdotc(n - l, &c[k + min(l + 1, n) * ldc], ldc, &B[l + min(l + 1, n) * ldb], ldb);
		vec = c[k + l * ldc] - (suml + sgn * conj(sumr));
		scaloc = One;
		a11 = conj(A[k + k * lda] + sgn * B[l + l * ldb]);
		da11 = abs(a11.real()) + abs(a11.imag());
		if (da11 <= smin) {
		    a11 = smin;
		    da11 = smin;
		    *info = 1;
		}
		db = abs(vec.real()) + abs(vec.imag());
		if (da11 < One && db > One) {
		    if (db > bignum * da11) {
			scaloc = One / db;
		    }
		}
		x11 = Cladiv(vec * scaloc, a11);
		if (scaloc != One) {
		    for (j = 0; j < n; j++) {
			CRscal(m, scaloc, &c[j * ldc + 1], 1);
		    }
		    (*scale) = (*scale) * scaloc;
		}
		c[k + l * ldc] = x11;
	    }
	}
    } else if (notrna && !notrnb) {
//Solve    A*X + ISGN*X*B' = C.
//The (K,L)th block of X is determined starting from
//bottom-left corner column by column by
//   A(K,K)*X(K,L) + ISGN*X(K,L)*B'(L,L) = C(K,L) - R(K,L)
//Where
//            M                          N
//  R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B'(L,J)]
//          I=K+1                      J=L+1
	for (l = n; l >= 1; l--) {
	    for (k = m; k >= 1; k--) {
		suml = Cdotu(m - k, &A[k + min(k + 1, m) * lda], lda, &c[min(k + 1, m) + l * ldc], 1);
		sumr = Cdotc(n - l, &c[k + min(l + 1, n) * ldc], ldc, &B[l + min(l + 1, n) * ldb], ldb);
		vec = c[k + l * ldc] - (suml + sgn * conj(sumr));
		scaloc = One;
		a11 = A[k + k * lda] + sgn * conj(B[l + l * ldb]);
		da11 = abs(a11.real()) + abs(a11.imag());
		if (da11 <= smin) {
		    a11 = smin;
		    da11 = smin;
		    *info = 1;
		}
		db = abs(vec.real()) + abs(vec.imag());
		if (da11 < One && db > One) {
		    if (db > bignum * da11) {
			scaloc = One / db;
		    }
		}
		x11 = Cladiv(vec * scaloc, a11);
		if (scaloc != One) {
		    for (j = 0; j < n; j++) {
			CRscal(m, scaloc, &c[j * ldc + 1], 1);
		    }
		    (*scale) = (*scale) * scaloc;
		}
		c[k + l * ldc] = x11;
	    }
	}
    }
    return;
}
