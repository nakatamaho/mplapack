/*
 * Copyright (c) 2008-2021
 *      Nakata, Maho
 *      All rights reserved.
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

#include <mpblas.h>
#include <mplapack.h>

void Ctgsja(const char *jobu, const char *jobv, const char *jobq, INTEGER const m, INTEGER const p, INTEGER const n, INTEGER const k, INTEGER const l, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, REAL const tola, REAL const tolb, REAL *alpha, REAL *beta, COMPLEX *u, INTEGER const ldu, COMPLEX *v, INTEGER const ldv, COMPLEX *q, INTEGER const ldq, COMPLEX *work, INTEGER &ncycle, INTEGER &info) {
    bool initu = false;
    bool wantu = false;
    bool initv = false;
    bool wantv = false;
    bool initq = false;
    bool wantq = false;
    const COMPLEX czero = (0.0, 0.0);
    const COMPLEX cone = (1.0, 0.0);
    bool upper = false;
    INTEGER kcycle = 0;
    const INTEGER maxit = 40;
    INTEGER i = 0;
    INTEGER j = 0;
    const REAL zero = 0.0;
    REAL a1 = 0.0;
    COMPLEX a2 = 0.0;
    REAL a3 = 0.0;
    REAL b1 = 0.0;
    REAL b3 = 0.0;
    COMPLEX b2 = 0.0;
    REAL csu = 0.0;
    COMPLEX snu = 0.0;
    REAL csv = 0.0;
    COMPLEX snv = 0.0;
    REAL csq = 0.0;
    COMPLEX snq = 0.0;
    REAL error = 0.0;
    REAL ssmin = 0.0;
    const REAL one = 1.0;
    REAL gamma = 0.0;
    REAL rwk = 0.0;
    //
    //  -- LAPACK computational routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Decode and test the input parameters
    //
    initu = Mlsame(jobu, "I");
    wantu = initu || Mlsame(jobu, "U");
    //
    initv = Mlsame(jobv, "I");
    wantv = initv || Mlsame(jobv, "V");
    //
    initq = Mlsame(jobq, "I");
    wantq = initq || Mlsame(jobq, "Q");
    //
    info = 0;
    if (!(initu || wantu || Mlsame(jobu, "N"))) {
        info = -1;
    } else if (!(initv || wantv || Mlsame(jobv, "N"))) {
        info = -2;
    } else if (!(initq || wantq || Mlsame(jobq, "N"))) {
        info = -3;
    } else if (m < 0) {
        info = -4;
    } else if (p < 0) {
        info = -5;
    } else if (n < 0) {
        info = -6;
    } else if (lda < max((INTEGER)1, m)) {
        info = -10;
    } else if (ldb < max((INTEGER)1, p)) {
        info = -12;
    } else if (ldu < 1 || (wantu && ldu < m)) {
        info = -18;
    } else if (ldv < 1 || (wantv && ldv < p)) {
        info = -20;
    } else if (ldq < 1 || (wantq && ldq < n)) {
        info = -22;
    }
    if (info != 0) {
        Mxerbla("Ctgsja", -info);
        return;
    }
    //
    //     Initialize U, V and Q, if necessary
    //
    if (initu) {
        Claset("Full", m, m, czero, cone, u, ldu);
    }
    if (initv) {
        Claset("Full", p, p, czero, cone, v, ldv);
    }
    if (initq) {
        Claset("Full", n, n, czero, cone, q, ldq);
    }
    //
    //     Loop until convergence
    //
    upper = false;
    for (kcycle = 1; kcycle <= maxit; kcycle = kcycle + 1) {
        //
        upper = !upper;
        //
        for (i = 1; i <= l - 1; i = i + 1) {
            for (j = i + 1; j <= l; j = j + 1) {
                //
                a1 = zero;
                a2 = czero;
                a3 = zero;
                if (k + i <= m) {
                    a1 = a[((k + i) - 1) + ((n - l + i) - 1) * lda].real();
                }
                if (k + j <= m) {
                    a3 = a[((k + j) - 1) + ((n - l + j) - 1) * lda].real();
                }
                //
                b1 = b[(i - 1) + ((n - l + i) - 1) * ldb].real();
                b3 = b[(j - 1) + ((n - l + j) - 1) * ldb].real();
                //
                if (upper) {
                    if (k + i <= m) {
                        a2 = a[((k + i) - 1) + ((n - l + j) - 1) * lda];
                    }
                    b2 = b[(i - 1) + ((n - l + j) - 1) * ldb];
                } else {
                    if (k + j <= m) {
                        a2 = a[((k + j) - 1) + ((n - l + i) - 1) * lda];
                    }
                    b2 = b[(j - 1) + ((n - l + i) - 1) * ldb];
                }
                //
                Clags2(upper, a1, a2, a3, b1, b2, b3, csu, snu, csv, snv, csq, snq);
                //
                //              Update (K+I)-th and (K+J)-th rows of matrix A: U**H *A
                //
                if (k + j <= m) {
                    Crot(l, &a[((k + j) - 1) + ((n - l + 1) - 1) * lda], lda, &a[((k + i) - 1) + ((n - l + 1) - 1) * lda], lda, csu, conj(snu));
                }
                //
                //              Update I-th and J-th rows of matrix B: V**H *B
                //
                Crot(l, &b[(j - 1) + ((n - l + 1) - 1) * ldb], ldb, &b[(i - 1) + ((n - l + 1) - 1) * ldb], ldb, csv, conj(snv));
                //
                //              Update (N-L+I)-th and (N-L+J)-th columns of matrices
                //              A and B: A*Q and B*Q
                //
                Crot(min(k + l, m), &a[((n - l + j) - 1) * lda], 1, &a[((n - l + i) - 1) * lda], 1, csq, snq);
                //
                Crot(l, &b[((n - l + j) - 1) * ldb], 1, &b[((n - l + i) - 1) * ldb], 1, csq, snq);
                //
                if (upper) {
                    if (k + i <= m) {
                        a[((k + i) - 1) + ((n - l + j) - 1) * lda] = czero;
                    }
                    b[(i - 1) + ((n - l + j) - 1) * ldb] = czero;
                } else {
                    if (k + j <= m) {
                        a[((k + j) - 1) + ((n - l + i) - 1) * lda] = czero;
                    }
                    b[(j - 1) + ((n - l + i) - 1) * ldb] = czero;
                }
                //
                //              Ensure that the diagonal elements of A and B are real.
                //
                if (k + i <= m) {
                    a[((k + i) - 1) + ((n - l + i) - 1) * lda] = a[((k + i) - 1) + ((n - l + i) - 1) * lda].real();
                }
                if (k + j <= m) {
                    a[((k + j) - 1) + ((n - l + j) - 1) * lda] = a[((k + j) - 1) + ((n - l + j) - 1) * lda].real();
                }
                b[(i - 1) + ((n - l + i) - 1) * ldb] = b[(i - 1) + ((n - l + i) - 1) * ldb].real();
                b[(j - 1) + ((n - l + j) - 1) * ldb] = b[(j - 1) + ((n - l + j) - 1) * ldb].real();
                //
                //              Update unitary matrices U, V, Q, if desired.
                //
                if (wantu && k + j <= m) {
                    Crot(m, &u[((k + j) - 1) * ldu], 1, &u[((k + i) - 1) * ldu], 1, csu, snu);
                }
                //
                if (wantv) {
                    Crot(p, &v[(j - 1) * ldv], 1, &v[(i - 1) * ldv], 1, csv, snv);
                }
                //
                if (wantq) {
                    Crot(n, &q[((n - l + j) - 1) * ldq], 1, &q[((n - l + i) - 1) * ldq], 1, csq, snq);
                }
                //
            }
        }
        //
        if (!upper) {
            //
            //           The matrices A13 and B13 were lower triangular at the start
            //           of the cycle, and are now upper triangular.
            //
            //           Convergence test: test the parallelism of the corresponding
            //           rows of A and B.
            //
            error = zero;
            for (i = 1; i <= min(l, m - k); i = i + 1) {
                Ccopy(l - i + 1, &a[((k + i) - 1) + ((n - l + i) - 1) * lda], lda, work, 1);
                Ccopy(l - i + 1, &b[(i - 1) + ((n - l + i) - 1) * ldb], ldb, &work[(l + 1) - 1], 1);
                Clapll(l - i + 1, work, 1, &work[(l + 1) - 1], 1, ssmin);
                error = max(error, ssmin);
            }
            //
            if (abs(error) <= min(tola, tolb)) {
                goto statement_50;
            }
        }
        //
        //        End of cycle loop
        //
    }
    //
    //     The algorithm has not converged after MAXIT cycles.
    //
    info = 1;
    goto statement_100;
//
statement_50:
    //
    //     If ERROR <= MIN(TOLA,TOLB), then the algorithm has converged.
    //     Compute the generalized singular value pairs (ALPHA, BETA), and
    //     set the triangular matrix R to array A.
    //
    for (i = 1; i <= k; i = i + 1) {
        alpha[i - 1] = one;
        beta[i - 1] = zero;
    }
    //
    for (i = 1; i <= min(l, m - k); i = i + 1) {
        //
        a1 = a[((k + i) - 1) + ((n - l + i) - 1) * lda].real();
        b1 = b[(i - 1) + ((n - l + i) - 1) * ldb].real();
        //
        if (a1 != zero) {
            gamma = b1 / a1;
            //
            if (gamma < zero) {
                CRscal(l - i + 1, -one, &b[(i - 1) + ((n - l + i) - 1) * ldb], ldb);
                if (wantv) {
                    CRscal(p, -one, &v[(i - 1) * ldv], 1);
                }
            }
            //
            Rlartg(abs(gamma), one, beta[(k + i) - 1], alpha[(k + i) - 1], rwk);
            //
            if (alpha[(k + i) - 1] >= beta[(k + i) - 1]) {
                CRscal(l - i + 1, one / alpha[(k + i) - 1], &a[((k + i) - 1) + ((n - l + i) - 1) * lda], lda);
            } else {
                CRscal(l - i + 1, one / beta[(k + i) - 1], &b[(i - 1) + ((n - l + i) - 1) * ldb], ldb);
                Ccopy(l - i + 1, &b[(i - 1) + ((n - l + i) - 1) * ldb], ldb, &a[((k + i) - 1) + ((n - l + i) - 1) * lda], lda);
            }
            //
        } else {
            //
            alpha[(k + i) - 1] = zero;
            beta[(k + i) - 1] = one;
            Ccopy(l - i + 1, &b[(i - 1) + ((n - l + i) - 1) * ldb], ldb, &a[((k + i) - 1) + ((n - l + i) - 1) * lda], lda);
        }
    }
    //
    //     Post-assignment
    //
    for (i = m + 1; i <= k + l; i = i + 1) {
        alpha[i - 1] = zero;
        beta[i - 1] = one;
    }
    //
    if (k + l < n) {
        for (i = k + l + 1; i <= n; i = i + 1) {
            alpha[i - 1] = zero;
            beta[i - 1] = zero;
        }
    }
//
statement_100:
    ncycle = kcycle;
    //
    //     End of Ctgsja
    //
}
