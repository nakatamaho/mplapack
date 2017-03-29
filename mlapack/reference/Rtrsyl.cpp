/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rtrsyl.cpp,v 1.9 2010/08/07 04:48:33 nakatamaho Exp $ 
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
Rtrsyl(const char *trana, const char *tranb, INTEGER isgn, INTEGER m, INTEGER n, REAL * A, INTEGER lda, REAL * B, INTEGER ldb, REAL * c, INTEGER ldc, REAL * scale, INTEGER * info)
{
    INTEGER j, k, l;
    REAL x[4];
    INTEGER k1, k2, l1, l2;
    REAL a11, db, da11, vec[4], dum[1], eps, sgn;
    INTEGER ierr;
    REAL smin, suml, sumr;
    INTEGER knext, lnext;
    REAL xnorm;
    REAL scaloc;
    REAL bignum;
    INTEGER notrna, notrnb;
    REAL smlnum;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1, mtemp2;
    REAL mtemp3, mtemp4;

//Decode and Test input parameters
    notrna = Mlsame(trana, "N");
    notrnb = Mlsame(tranb, "N");
    *info = 0;
    if (!notrna && !Mlsame(trana, "T") && !Mlsame(trana, "C")) {
	*info = -1;
    } else if (!notrnb && !Mlsame(tranb, "T") && !Mlsame(tranb, "C")) {
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
	Mxerbla("Rtrsyl", -(*info));
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
    smlnum = smlnum * (REAL) double (m * n) / eps;
    bignum = One / smlnum;
    mtemp1 = smlnum, mtemp2 = eps * Rlange("M", m, m, A, lda, dum);
    mtemp3 = max(mtemp1, mtemp2), mtemp4 = eps * Rlange("M", n, n, B, ldb, dum);
    smin = max(mtemp3, mtemp4);
    *scale = One;
    sgn = (isgn);
    if (notrna && notrnb) {
//Solve    A*X + ISGN*X*B = scale*C.
//The (K,L)th block of X is determined starting from
//bottom-left corner column by column by
// A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
//Where
//          M                         L-1
//R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)].
//        I=K+1                       J=1
//Start column loop (index = L)
//L1 (L2) : column index of the first (first) row of X(K,L).
	lnext = 1;
	for (l = 0; l < n; l++) {
	    if (l < lnext) {
		goto L60;
	    }
	    if (l == n) {
		l1 = l;
		l2 = l;
	    } else {
		if (B[l + 1 + l * ldb] != Zero) {
		    l1 = l;
		    l2 = l + 1;
		    lnext = l + 2;
		} else {
		    l1 = l;
		    l2 = l;
		    lnext = l + 1;
		}
	    }
//Start row loop (index = K)
//K1 (K2): row index of the first (last) row of X(K,L).
	    knext = m;
	    for (k = m; k >= 1; k--) {
		if (k > knext) {
		    goto L50;
		}
		if (k == 1) {
		    k1 = k;
		    k2 = k;
		} else {
		    if (A[k + (k - 1) * lda] != Zero) {
			k1 = k - 1;
			k2 = k;
			knext = k - 2;
		    } else {
			k1 = k;
			k2 = k;
			knext = k - 1;
		    }
		}
		if (l1 == l2 && k1 == k2) {
		    suml = Rdot(m - k1, &A[k1 + min(k1 + 1, m) * lda], lda, &c[min(k1 + 1, m) + l1 * ldc], 1);
		    sumr = Rdot(l1 - 1, &c[k1 + ldc], ldc, &B[l1 * ldb + 1], 1);
		    vec[0] = c[k1 + l1 * ldc] - (suml + sgn * sumr);
		    scaloc = One;
		    a11 = A[k1 + k1 * lda] + sgn * B[l1 + l1 * ldb];
		    da11 = abs(a11);
		    if (da11 <= smin) {
			a11 = smin;
			da11 = smin;
			*info = 1;
		    }
		    db = abs(vec[0]);
		    if (da11 < One && db > One) {
			if (db > bignum * da11) {
			    scaloc = One / db;
			}
		    }
		    x[0] = vec[0] * scaloc / a11;
		    if (scaloc != One) {
			for (j = 0; j < n; j++) {
			    Rscal(m, scaloc, &c[j * ldc + 1], 1);
			}
			*scale *= scaloc;
		    }
		    c[k1 + l1 * ldc] = x[0];

		} else if (l1 == l2 && k1 != k2) {
		    suml = Rdot(m - k2, &A[k1 + min(k2 + 1, m) * lda], lda, &c[min(k2 + 1, m) + l1 * ldc], 1);
		    sumr = Rdot(l1 - 1, &c[k1 + ldc], ldc, &B[l1 * ldb + 1], 1);
		    vec[0] = c[k1 + l1 * ldc] - (suml + sgn * sumr);
		    suml = Rdot(m - k2, &A[k2 + min(k2 + 1, m) * lda], lda, &c[min(k2 + 1, m) + l1 * ldc], 1);
		    sumr = Rdot(l1 - 1, &c[k2 + ldc], ldc, &B[l1 * ldb + 1], 1);
		    vec[1] = c[k2 + l1 * ldc] - (suml + sgn * sumr);
		    Rlaln2(MFALSE, 2, 1, smin, One, &A[k1 + k1 * lda], lda, One, One, vec, 2, -sgn * B[l1 + l1 * ldb], Zero, x, 2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }
		    if (scaloc != One) {
			for (j = 0; j < n; j++) {
			    Rscal(m, scaloc, &c[j * ldc + 1], 1);
			}
			*scale *= scaloc;
		    }
		    c[k1 + l1 * ldc] = x[0];
		    c[k2 + l1 * ldc] = x[0];

		} else if (l1 != l2 && k1 == k2) {
		    suml = Rdot(m - k1, &A[k1 + min(k1 + 1, m) * lda], lda, &c[min(k1 + 1, m) + l1 * ldc], 1);
		    sumr = Rdot(l1 - 1, &c[k1 + ldc], ldc, &B[l1 * ldb + 1], 1);
		    vec[0] = sgn * (c[k1 + l1 * ldc] - (suml + sgn * sumr));
		    suml = Rdot(m - k1, &A[k1 + min(k1 + 1, m) * lda], lda, &c[min(k1 + 1, m) + l2 * ldc], 1);
		    sumr = Rdot(l1 - 1, &c[k1 + ldc], ldc, &B[l2 * ldb + 1], 1);
		    vec[1] = sgn * (c[k1 + l2 * ldc] - (suml + sgn * sumr));
		    Rlaln2(MTRUE, 2, 1, smin, One, &B[l1 + l1 * ldb], ldb, One, One, vec, 2, -sgn * A[k1 + k1 * lda], Zero, x, 2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }
		    if (scaloc != One) {
			for (j = 0; j < n; j++) {
			    Rscal(m, scaloc, &c[j * ldc + 1], 1);
			}
			*scale = *scale * scaloc;
		    }
		    c[k1 + l1 * ldc] = x[0];
		    c[k1 + l2 * ldc] = x[0];
		} else if (l1 != l2 && k1 != k2) {
		    suml = Rdot(m - k2, &A[k1 + min(k2 + 1, m) * lda], lda, &c[min(k2 + 1, m) + l1 * ldc], 1);
		    sumr = Rdot(l1 - 1, &c[k1 + ldc], ldc, &B[l1 * ldb + 1], 1);
		    vec[0] = c[k1 + l1 * ldc] - (suml + sgn * sumr);
		    suml = Rdot(m - k2, &A[k1 + min(k2 + 1, m) * lda], lda, &c[min(k2 + 1, m) + l2 * ldc], 1);
		    sumr = Rdot(l1 - 1, &c[k1 + ldc], ldc, &B[l2 * ldb + 1], 1);
		    vec[2] = c[k1 + l2 * ldc] - (suml + sgn * sumr);
		    suml = Rdot(m - k2, &A[k2 + min(k2 + 1, m) * lda], lda, &c[min(k2 + 1, m) + l1 * ldc], 1);
		    sumr = Rdot(l1 - 1, &c[k2 + ldc], ldc, &B[l1 * ldb + 1], 1);
		    vec[1] = c[k2 + l1 * ldc] - (suml + sgn * sumr);
		    suml = Rdot(m - k2, &A[k2 + min(k2 + 1, m) * lda], lda, &c[min(k2 + 1, m) + l2 * ldc], 1);
		    sumr = Rdot(l1 - 1, &c[k2 + ldc], ldc, &B[l2 * ldb + 1], 1);
		    vec[3] = c[k2 + l2 * ldc] - (suml + sgn * sumr);
		    Rlasy2(MFALSE, MFALSE, isgn, 2, 2, &A[k1 + k1 * lda], lda, &B[l1 + l1 * ldb], ldb, vec, 2, &scaloc, x, 2, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }
		    if (scaloc != One) {
			for (j = 0; j < n; j++) {
			    Rscal(m, scaloc, &c[j * ldc + 1], 1);
			}
			*scale *= scaloc;
		    }
		    c[k1 + l1 * ldc] = x[0];
		    c[k1 + l2 * ldc] = x[2];
		    c[k2 + l1 * ldc] = x[0];
		    c[k2 + l2 * ldc] = x[3];
		}
	      L50:
		;
	    }
	  L60:
	    ;
	}
    } else if (!notrna && notrnb) {
//Solve    A' *X + ISGN*X*B = scale*C.
//The (K,L)th block of X is determined starting from
//upper-left corner column by column by
//  A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
//Where
//           K-1                        L-1
//  R(K,L) = SUM [A(I,K)'*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]
//           I=1                        J=1
//Start column loop (index = L)
//L1 (L2): column index of the first (last) row of X(K,L)
	lnext = 1;
	for (l = 0; l < n; l++) {
	    if (l < lnext) {
		goto L120;
	    }
	    if (l == n) {
		l1 = l;
		l2 = l;
	    } else {
		if (B[l + 1 + l * ldb] != Zero) {
		    l1 = l;
		    l2 = l + 1;
		    lnext = l + 2;
		} else {
		    l1 = l;
		    l2 = l;
		    lnext = l + 1;
		}
	    }
//Start row loop (index = K)
//K1 (K2): row index of the first (last) row of X(K,L)
	    knext = 1;
	    for (k = 0; k < m; k++) {
		if (k < knext) {
		    goto L110;
		}
		if (k == m) {
		    k1 = k;
		    k2 = k;
		} else {
		    if (A[k + 1 + k * lda] != Zero) {
			k1 = k;
			k2 = k + 1;
			knext = k + 2;
		    } else {
			k1 = k;
			k2 = k;
			knext = k + 1;
		    }
		}
		if (l1 == l2 && k1 == k2) {
		    suml = Rdot(k1 - 1, &A[k1 * lda], 1, &c[l1 * ldc + 1], 1);
		    sumr = Rdot(l1 - 1, &c[k1 + ldc], ldc, &B[l1 * ldb + 1], 1);
		    vec[0] = c[k1 + l1 * ldc] - (suml + sgn * sumr);
		    scaloc = One;
		    a11 = A[k1 + k1 * lda] + sgn * B[l1 + l1 * ldb];
		    da11 = abs(a11);
		    if (da11 <= smin) {
			a11 = smin;
			da11 = smin;
			*info = 1;
		    }
		    db = abs(vec[0]);
		    if (da11 < One && db > One) {
			if (db > bignum * da11) {
			    scaloc = One / db;
			}
		    }
		    x[0] = vec[0] * scaloc / a11;
		    if (scaloc != One) {
			for (j = 0; j < n; j++) {
			    Rscal(m, scaloc, &c[j * ldc + 1], 1);
			}
			*scale *= scaloc;
		    }
		    c[k1 + l1 * ldc] = x[0];
		} else if (l1 == l2 && k1 != k2) {
		    suml = Rdot(k1 - 1, &A[k1 * lda], 1, &c[l1 * ldc + 1], 1);
		    sumr = Rdot(l1 - 1, &c[k1 + ldc], ldc, &B[l1 * ldb + 1], 1);
		    vec[0] = c[k1 + l1 * ldc] - (suml + sgn * sumr);
		    suml = Rdot(k1 - 1, &A[k2 * lda], 1, &c[l1 * ldc + 1], 1);
		    sumr = Rdot(l1 - 1, &c[k2 + ldc], ldc, &B[l1 * ldb + 1], 1);
		    vec[1] = c[k2 + l1 * ldc] - (suml + sgn * sumr);
		    Rlaln2(MTRUE, 2, 1, smin, One, &A[k1 + k1 * lda], lda, One, One, vec, 2, -sgn * B[l1 + l1 * ldb], Zero, x, 2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }
		    if (scaloc != One) {
			for (j = 0; j < n; j++) {
			    Rscal(m, scaloc, &c[j * ldc + 1], 1);
			}
			*scale = *scale * scaloc;
		    }
		    c[k1 + l1 * ldc] = x[0];
		    c[k2 + l1 * ldc] = x[0];
		} else if (l1 != l2 && k1 == k2) {
		    suml = Rdot(k1 - 1, &A[k1 * lda], 1, &c[l1 * ldc + 1], 1);
		    sumr = Rdot(l1 - 1, &c[k1 + ldc], ldc, &B[l1 * ldb + 1], 1);
		    vec[0] = sgn * (c[k1 + l1 * ldc] - (suml + sgn * sumr));
		    suml = Rdot(k1 - 1, &A[k1 * lda], 1, &c[l2 * ldc + 1], 1);
		    sumr = Rdot(l1 - 1, &c[k1 + ldc], ldc, &B[l2 * ldb + 1], 1);
		    vec[1] = sgn * (c[k1 + l2 * ldc] - (suml + sgn * sumr));
		    Rlaln2(MTRUE, 2, 1, smin, One, &B[l1 + l1 * ldb], ldb, One, One, vec, 2, -sgn * A[k1 + k1 * lda], Zero, x, 2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }
		    if (scaloc != One) {
			for (j = 0; j < n; j++) {
			    Rscal(m, scaloc, &c[j * ldc + 1], 1);
			}
			*scale = *scale * scaloc;
		    }
		    c[k1 + l1 * ldc] = x[0];
		    c[k1 + l2 * ldc] = x[0];

		} else if (l1 != l2 && k1 != k2) {
		    suml = Rdot(k1 - 1, &A[k1 * lda], 1, &c[l1 * ldc + 1], 1);
		    sumr = Rdot(l1 - 1, &c[k1 + ldc], ldc, &B[l1 * ldb + 1], 1);
		    vec[0] = c[k1 + l1 * ldc] - (suml + sgn * sumr);
		    suml = Rdot(k1 - 1, &A[k1 * lda], 1, &c[l2 * ldc + 1], 1);
		    sumr = Rdot(l1 - 1, &c[k1 + ldc], ldc, &B[l2 * ldb + 1], 1);
		    vec[2] = c[k1 + l2 * ldc] - (suml + sgn * sumr);
		    suml = Rdot(k1 - 1, &A[k2 * lda], 1, &c[l1 * ldc + 1], 1);
		    sumr = Rdot(l1 - 1, &c[k2 + ldc], ldc, &B[l1 * ldb + 1], 1);
		    vec[1] = c[k2 + l1 * ldc] - (suml + sgn * sumr);
		    suml = Rdot(k1 - 1, &A[k2 * lda], 1, &c[l2 * ldc + 1], 1);
		    sumr = Rdot(l1 - 1, &c[k2 + ldc], ldc, &B[l2 * ldb + 1], 1);
		    vec[3] = c[k2 + l2 * ldc] - (suml + sgn * sumr);
		    Rlasy2(MTRUE, MFALSE, isgn, 2, 2, &A[k1 + k1 * lda], lda, &B[l1 + l1 * ldb], ldb, vec, 2, &scaloc, x, 2, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }
		    if (scaloc != One) {
			for (j = 0; j < n; j++) {
			    Rscal(m, scaloc, &c[j * ldc + 1], 1);
			}
			*scale = *scale * scaloc;
		    }
		    c[k1 + l1 * ldc] = x[0];
		    c[k1 + l2 * ldc] = x[2];
		    c[k2 + l1 * ldc] = x[0];
		    c[k2 + l2 * ldc] = x[3];
		}
	      L110:
		;
	    }
	  L120:
	    ;
	}
    } else if (!notrna && !notrnb) {
//Solve    A'*X + ISGN*X*B' = scale*C.
//The (K,L)th block of X is determined starting from
//top-right corner column by column by
//   A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L)' = C(K,L) - R(K,L)
//Where
//             K-1                          N
//    R(K,L) = SUM [A(I,K)'*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)'].
//             I=1                        J=L+1
//Start column loop (index = L)
//L1 (L2): column index of the first (last) row of X(K,L)
	lnext = n;
	for (l = n; l >= 1; l--) {
	    if (l > lnext) {
		goto L180;
	    }
	    if (l == 1) {
		l1 = l;
		l2 = l;
	    } else {
		if (B[l + (l - 1) * ldb] != Zero) {
		    l1 = l - 1;
		    l2 = l;
		    lnext = l - 2;
		} else {
		    l1 = l;
		    l2 = l;
		    lnext = l - 1;
		}
	    }
//Start row loop (index = K)
//K1 (K2): row index of the first (last) row of X(K,L)
	    knext = 1;
	    for (k = 0; k < m; k++) {
		if (k < knext) {
		    goto L170;
		}
		if (k == m) {
		    k1 = k;
		    k2 = k;
		} else {
		    if (A[k + 1 + k * lda] != Zero) {
			k1 = k;
			k2 = k + 1;
			knext = k + 2;
		    } else {
			k1 = k;
			k2 = k;
			knext = k + 1;
		    }
		}
		if (l1 == l2 && k1 == k2) {
		    suml = Rdot(k1 - 1, &A[k1 * lda], 1, &c[l1 * ldc + 1], 1);
		    sumr = Rdot(n - l1, &c[k1 + min(l1 + 1, n) * ldc], ldc, &B[l1 + min(l1 + 1, n) * ldb], ldb);
		    vec[0] = c[k1 + l1 * ldc] - (suml + sgn * sumr);
		    scaloc = One;
		    a11 = A[k1 + k1 * lda] + sgn * B[l1 + l1 * ldb];
		    da11 = abs(a11);
		    if (da11 <= smin) {
			a11 = smin;
			da11 = smin;
			*info = 1;
		    }
		    db = abs(vec[0]);
		    if (da11 < One && db > One) {
			if (db > bignum * da11) {
			    scaloc = One / db;
			}
		    }
		    x[0] = vec[0] * scaloc / a11;
		    if (scaloc != One) {
			for (j = 0; j < n; j++) {
			    Rscal(m, scaloc, &c[j * ldc + 1], 1);
			}
			*scale = *scale * scaloc;
		    }
		    c[k1 + l1 * ldc] = x[0];
		} else if (l1 == l2 && k1 != k2) {
		    suml = Rdot(k1 - 1, &A[k1 * lda], 1, &c[l1 * ldc + 1], 1);
		    sumr = Rdot(n - l2, &c[k1 + min(l2 + 1, n) * ldc], ldc, &B[l1 + min(l2 + 1, n) * ldb], ldb);
		    vec[0] = c[k1 + l1 * ldc] - (suml + sgn * sumr);
		    suml = Rdot(k1 - 1, &A[k2 * lda], 1, &c[l1 * ldc + 1], 1);
		    sumr = Rdot(n - l2, &c[k2 + min(l2 + 1, n) * ldc], ldc, &B[l1 + min(l2 + 1, n) * ldb], ldb);
		    vec[1] = c[k2 + l1 * ldc] - (suml + sgn * sumr);
		    Rlaln2(MTRUE, 2, 1, smin, One, &A[k1 + k1 * lda], lda, One, One, vec, 2, -sgn * B[l1 + l1 * ldb], Zero, x, 2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }
		    if (scaloc != One) {
			for (j = 0; j < n; j++) {
			    Rscal(m, scaloc, &c[j * ldc + 1], 1);
			}
			*scale = *scale * scaloc;
		    }
		    c[k1 + l1 * ldc] = x[0];
		    c[k2 + l1 * ldc] = x[0];
		} else if (l1 != l2 && k1 == k2) {
		    suml = Rdot(k1 - 1, &A[k1 * lda], 1, &c[l1 * ldc + 1], 1);
		    sumr = Rdot(n - l2, &c[k1 + min(l2 + 1, n) * ldc], ldc, &B[l1 + min(l2 + 1, n) * ldb], ldb);
		    vec[0] = sgn * (c[k1 + l1 * ldc] - (suml + sgn * sumr));
		    suml = Rdot(k1 - 1, &A[k1 * lda], 1, &c[l2 * ldc + 1], 1);
		    sumr = Rdot(n - l2, &c[k1 + min(l2 + 1, n) * ldc], ldc, &B[l2 + min(l2 + 1, n) * ldb], ldb);
		    vec[1] = sgn * (c[k1 + l2 * ldc] - (suml + sgn * sumr));
		    Rlaln2(MFALSE, 2, 1, smin, One, &B[l1 + l1 * ldb], ldb, One, One, vec, 2, -sgn * A[k1 + k1 * lda], Zero, x, 2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }
		    if (scaloc != One) {
			for (j = 0; j < n; j++) {
			    Rscal(m, scaloc, &c[j * ldc + 1], 1);
			}
			*scale = *scale * scaloc;
		    }
		    c[k1 + l1 * ldc] = x[0];
		    c[k1 + l2 * ldc] = x[0];
		} else if (l1 != l2 && k1 != k2) {
		    suml = Rdot(k1 - 1, &A[k1 * lda], 1, &c[l1 * ldc + 1], 1);
		    sumr = Rdot(n - l2, &c[k1 + min(l2 + 1, n) * ldc], ldc, &B[l1 + min(l2 + 1, n) * ldb], ldb);
		    vec[0] = c[k1 + l1 * ldc] - (suml + sgn * sumr);
		    suml = Rdot(k1 - 1, &A[k1 * lda], 1, &c[l2 * ldc + 1], 1);
		    sumr = Rdot(n - l2, &c[k1 + min(l2 + 1, n) * ldc], ldc, &B[l2 + min(l2 + 1, n) * ldb], ldb);
		    vec[2] = c[k1 + l2 * ldc] - (suml + sgn * sumr);
		    suml = Rdot(k1 - 1, &A[k2 * lda], 1, &c[l1 * ldc + 1], 1);
		    sumr = Rdot(n - l2, &c[k2 + min(l2 + 1, n) * ldc], ldc, &B[l1 + min(l2 + 1, n) * ldb], ldb);
		    vec[1] = c[k2 + l1 * ldc] - (suml + sgn * sumr);
		    suml = Rdot(k1 - 1, &A[k2 * lda], 1, &c[l2 * ldc + 1], 1);
		    sumr = Rdot(n - l2, &c[k2 + min(l2 + 1, n) * ldc], ldc, &B[l2 + min(l2 + 1, n) * ldb], ldb);
		    vec[3] = c[k2 + l2 * ldc] - (suml + sgn * sumr);
		    Rlasy2(MTRUE, MTRUE, isgn, 2, 2, &A[k1 + k1 * lda], lda, &B[l1 + l1 * ldb], ldb, vec, 2, &scaloc, x, 2, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }
		    if (scaloc != One) {
			for (j = 0; j < n; j++) {
			    Rscal(m, scaloc, &c[j * ldc + 1], 1);
			}
			*scale = *scale * scaloc;
		    }
		    c[k1 + l1 * ldc] = x[0];
		    c[k1 + l2 * ldc] = x[2];
		    c[k2 + l1 * ldc] = x[0];
		    c[k2 + l2 * ldc] = x[3];
		}
	      L170:
		;
	    }
	  L180:
	    ;
	}
    } else if (notrna && !notrnb) {
//Solve    A*X + ISGN*X*B' = scale*C.
//The (K,L)th block of X is determined starting from
//bottom-right corner column by column by

//    A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)' = C(K,L) - R(K,L)

//Where
//              M                          N
//    R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)'].
//            I=K+1                      J=L+1

//Start column loop (index = L)
//L1 (L2): column index of the first (last) row of X(K,L)
	lnext = n;
	for (l = n; l >= 1; l--) {
	    if (l > lnext) {
		goto L240;
	    }
	    if (l == 1) {
		l1 = l;
		l2 = l;
	    } else {
		if (B[l + (l - 1) * ldb] != Zero) {
		    l1 = l - 1;
		    l2 = l;
		    lnext = l - 2;
		} else {
		    l1 = l;
		    l2 = l;
		    lnext = l - 1;
		}
	    }
//Start row loop (index = K)
//K1 (K2): row index of the first (last) row of X(K,L)
	    knext = m;
	    for (k = m; k >= 1; k--) {
		if (k > knext) {
		    goto L230;
		}
		if (k == 1) {
		    k1 = k;
		    k2 = k;
		} else {
		    if (A[k + (k - 1) * lda] != Zero) {
			k1 = k - 1;
			k2 = k;
			knext = k - 2;
		    } else {
			k1 = k;
			k2 = k;
			knext = k - 1;
		    }
		}

		if (l1 == l2 && k1 == k2) {
		    suml = Rdot(m - k1, &A[k1 + min(k1 + 1, m) * lda], lda, &c[min(k1 + 1, m) + l1 * ldc], 1);
		    sumr = Rdot(n - l1, &c[k1 + min(l1 + 1, n) * ldc], ldc, &B[l1 + min(l1 + 1, n) * ldb], ldb);
		    vec[0] = c[k1 + l1 * ldc] - (suml + sgn * sumr);
		    scaloc = One;
		    a11 = A[k1 + k1 * lda] + sgn * B[l1 + l1 * ldb];
		    da11 = abs(a11);
		    if (da11 <= smin) {
			a11 = smin;
			da11 = smin;
			*info = 1;
		    }
		    db = abs(vec[0]);
		    if (da11 < One && db > One) {
			if (db > bignum * da11) {
			    scaloc = One / db;
			}
		    }
		    x[0] = vec[0] * scaloc / a11;
		    if (scaloc != One) {
			for (j = 0; j < n; j++) {
			    Rscal(m, scaloc, &c[j * ldc + 1], 1);
			}
			*scale = *scale * scaloc;
		    }
		    c[k1 + l1 * ldc] = x[0];
		} else if (l1 == l2 && k1 != k2) {
		    suml = Rdot(m - k2, &A[k1 + min(k2 + 1, m) * lda], lda, &c[min(k2 + 1, m) + l1 * ldc], 1);
		    sumr = Rdot(n - l2, &c[k1 + min(l2 + 1, n) * ldc], ldc, &B[l1 + min(l2 + 1, n) * ldb], ldb);
		    vec[0] = c[k1 + l1 * ldc] - (suml + sgn * sumr);
		    suml = Rdot(m - k2, &A[k2 + min(k2 + 1, m) * lda], lda, &c[min(k2 + 1, m) + l1 * ldc], 1);
		    sumr = Rdot(n - l2, &c[k2 + min(l2 + 1, n) * ldc], ldc, &B[l1 + min(l2 + 1, n) * ldb], ldb);
		    vec[1] = c[k2 + l1 * ldc] - (suml + sgn * sumr);
		    Rlaln2(MFALSE, 2, 1, smin, One, &A[k1 + k1 * lda], lda, One, One, vec, 2, -sgn * B[l1 + l1 * ldb], Zero, x, 2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }
		    if (scaloc != One) {
			for (j = 0; j < n; j++) {
			    Rscal(m, scaloc, &c[j * ldc + 1], 1);
			}
			*scale *= scaloc;
		    }
		    c[k1 + l1 * ldc] = x[0];
		    c[k2 + l1 * ldc] = x[0];
		} else if (l1 != l2 && k1 == k2) {
		    suml = Rdot(m - k1, &A[k1 + min(k1 + 1, m) * lda], lda, &c[min(k1 + 1, m) + l1 * ldc], 1);
		    sumr = Rdot(n - l2, &c[k1 + min(l2 + 1, n) * ldc], ldc, &B[l1 + min(l2 + 1, n) * ldb], ldb);
		    vec[0] = sgn * (c[k1 + l1 * ldc] - (suml + sgn * sumr));
		    suml = Rdot(m - k1, &A[k1 + min(k1 + 1, m) * lda], lda, &c[min(k1 + 1, m) + l2 * ldc], 1);
		    sumr = Rdot(n - l2, &c[k1 + min(l2 + 1, n) * ldc], ldc, &B[l2 + min(l2 + 1, n) * ldb], ldb);
		    vec[1] = sgn * (c[k1 + l2 * ldc] - (suml + sgn * sumr));
		    Rlaln2(MFALSE, 2, 1, smin, One, &B[l1 + l1 * ldb], ldb, One, One, vec, 2, -sgn * A[k1 + k1 * lda], Zero, x, 2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }
		    if (scaloc != One) {
			for (j = 0; j < n; j++) {
			    Rscal(m, scaloc, &c[j * ldc + 1], 1);
			}
			*scale = *scale * scaloc;
		    }
		    c[k1 + l1 * ldc] = x[0];
		    c[k1 + l2 * ldc] = x[0];
		} else if (l1 != l2 && k1 != k2) {
		    suml = Rdot(m - k2, &A[k1 + min(k2 + 1, m) * lda], lda, &c[min(k2 + 1, m) + l1 * ldc], 1);
		    sumr = Rdot(n - l2, &c[k1 + min(l2 + 1, n) * ldc], ldc, &B[l1 + min(l2 + 1, n) * ldb], ldb);
		    vec[0] = c[k1 + l1 * ldc] - (suml + sgn * sumr);
		    suml = Rdot(m - k2, &A[k1 + min(k2 + 1, m) * lda], lda, &c[min(k2 + 1, m) + l2 * ldc], 1);
		    sumr = Rdot(n - l2, &c[k1 + min(l2 + 1, n) * ldc], ldc, &B[l2 + min(l2 + 1, n) * ldb], ldb);
		    vec[2] = c[k1 + l2 * ldc] - (suml + sgn * sumr);
		    suml = Rdot(m - k2, &A[k2 + min(k2 + 1, m) * lda], lda, &c[min(k2 + 1, m) + l1 * ldc], 1);
		    sumr = Rdot(n - l2, &c[k2 + min(l2 + 1, n) * ldc], ldc, &B[l1 + min(l2 + 1, n) * ldb], ldb);
		    vec[1] = c[k2 + l1 * ldc] - (suml + sgn * sumr);
		    suml = Rdot(m - k2, &A[k2 + min(k2 + 1, m) * lda], lda, &c[min(k2 + 1, m) + l2 * ldc], 1);
		    sumr = Rdot(n - l2, &c[k2 + min(l2 + 1, n) * ldc], ldc, &B[l2 + min(l2 + 1, n) * ldb], ldb);
		    vec[3] = c[k2 + l2 * ldc] - (suml + sgn * sumr);
		    Rlasy2(MFALSE, MTRUE, isgn, 2, 2, &A[k1 + k1 * lda], lda, &B[l1 + l1 * ldb], ldb, vec, 2, &scaloc, x, 2, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 1;
		    }
		    if (scaloc != One) {
			for (j = 0; j < n; j++) {
			    Rscal(m, scaloc, &c[j * ldc + 1], 1);
			}
			*scale = *scale * scaloc;
		    }
		    c[k1 + l1 * ldc] = x[0];
		    c[k1 + l2 * ldc] = x[2];
		    c[k2 + l1 * ldc] = x[0];
		    c[k2 + l2 * ldc] = x[3];
		}
	      L230:
		;
	    }
	  L240:
	    ;
	}
    }
    return;
}
