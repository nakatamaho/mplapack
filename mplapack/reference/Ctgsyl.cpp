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

void Ctgsyl(const char *trans, INTEGER const ijob, INTEGER const m, INTEGER const n, COMPLEX *a, INTEGER const lda, COMPLEX *b, INTEGER const ldb, COMPLEX *c, INTEGER const ldc, COMPLEX *d, INTEGER const ldd, COMPLEX *e, INTEGER const lde, COMPLEX *f, INTEGER const ldf, REAL &scale, REAL &dif, COMPLEX *work, INTEGER const lwork, INTEGER *iwork, INTEGER &info) {
    bool notran = false;
    bool lquery = false;
    INTEGER lwmin = 0;
    INTEGER mb = 0;
    INTEGER nb = 0;
    INTEGER isolve = 0;
    INTEGER ifunc = 0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    INTEGER iround = 0;
    const REAL one = 1.0;
    const REAL zero = 0.0;
    REAL Rscale = 0.0;
    REAL dsum = 0.0;
    INTEGER pq = 0;
    REAL scale2 = 0.0;
    INTEGER p = 0;
    INTEGER i = 0;
    INTEGER q = 0;
    INTEGER j = 0;
    INTEGER js = 0;
    INTEGER je = 0;
    INTEGER is = 0;
    INTEGER ie = 0;
    REAL scaloc = 0.0;
    INTEGER linfo = 0;
    INTEGER k = 0;
    //
    //     Decode and test input parameters
    //
    info = 0;
    notran = Mlsame(trans, "N");
    lquery = (lwork == -1);
    //
    if (!notran && !Mlsame(trans, "C")) {
        info = -1;
    } else if (notran) {
        if ((ijob < 0) || (ijob > 4)) {
            info = -2;
        }
    }
    if (info == 0) {
        if (m <= 0) {
            info = -3;
        } else if (n <= 0) {
            info = -4;
        } else if (lda < max((INTEGER)1, m)) {
            info = -6;
        } else if (ldb < max((INTEGER)1, n)) {
            info = -8;
        } else if (ldc < max((INTEGER)1, m)) {
            info = -10;
        } else if (ldd < max((INTEGER)1, m)) {
            info = -12;
        } else if (lde < max((INTEGER)1, n)) {
            info = -14;
        } else if (ldf < max((INTEGER)1, m)) {
            info = -16;
        }
    }
    //
    if (info == 0) {
        if (notran) {
            if (ijob == 1 || ijob == 2) {
                lwmin = max((INTEGER)1, 2 * m * n);
            } else {
                lwmin = 1;
            }
        } else {
            lwmin = 1;
        }
        work[1 - 1] = lwmin;
        //
        if (lwork < lwmin && !lquery) {
            info = -20;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Ctgsyl", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0) {
        scale = 1;
        if (notran) {
            if (ijob != 0) {
                dif = 0.0;
            }
        }
        return;
    }
    //
    //     Determine  optimal block sizes MB and NB
    //
    mb = iMlaenv(2, "Ctgsyl", trans, m, n, -1, -1);
    nb = iMlaenv(5, "Ctgsyl", trans, m, n, -1, -1);
    //
    isolve = 1;
    ifunc = 0;
    if (notran) {
        if (ijob >= 3) {
            ifunc = ijob - 2;
            Claset("F", m, n, czero, czero, c, ldc);
            Claset("F", m, n, czero, czero, f, ldf);
        } else if (ijob >= 1 && notran) {
            isolve = 2;
        }
    }
    //
    if ((mb <= 1 && nb <= 1) || (mb >= m && nb >= n)) {
        //
        //        Use unblocked Level 2 solver
        //
        for (iround = 1; iround <= isolve; iround = iround + 1) {
            //
            scale = one;
            Rscale = zero;
            dsum = one;
            pq = m * n;
            Ctgsy2(trans, ifunc, m, n, a, lda, b, ldb, c, ldc, d, ldd, e, lde, f, ldf, scale, dsum, Rscale, info);
            if (Rscale != zero) {
                if (ijob == 1 || ijob == 3) {
                    dif = sqrt(castREAL(2 * m * n)) / (Rscale * sqrt(dsum));
                } else {
                    dif = sqrt(castREAL(pq)) / (Rscale * sqrt(dsum));
                }
            }
            if (isolve == 2 && iround == 1) {
                if (notran) {
                    ifunc = ijob;
                }
                scale2 = scale;
                Clacpy("F", m, n, c, ldc, work, m);
                Clacpy("F", m, n, f, ldf, &work[(m * n + 1) - 1], m);
                Claset("F", m, n, czero, czero, c, ldc);
                Claset("F", m, n, czero, czero, f, ldf);
            } else if (isolve == 2 && iround == 2) {
                Clacpy("F", m, n, work, m, c, ldc);
                Clacpy("F", m, n, &work[(m * n + 1) - 1], m, f, ldf);
                scale = scale2;
            }
        }
        //
        return;
        //
    }
    //
    //     Determine block structure of A
    //
    p = 0;
    i = 1;
statement_40:
    if (i > m) {
        goto statement_50;
    }
    p++;
    iwork[p - 1] = i;
    i += mb;
    if (i >= m) {
        goto statement_50;
    }
    goto statement_40;
statement_50:
    iwork[(p + 1) - 1] = m + 1;
    if (iwork[p - 1] == iwork[(p + 1) - 1]) {
        p = p - 1;
    }
    //
    //     Determine block structure of B
    //
    q = p + 1;
    j = 1;
statement_60:
    if (j > n) {
        goto statement_70;
    }
    //
    q++;
    iwork[q - 1] = j;
    j += nb;
    if (j >= n) {
        goto statement_70;
    }
    goto statement_60;
//
statement_70:
    iwork[(q + 1) - 1] = n + 1;
    if (iwork[q - 1] == iwork[(q + 1) - 1]) {
        q = q - 1;
    }
    //
    if (notran) {
        for (iround = 1; iround <= isolve; iround = iround + 1) {
            //
            //           Solve (I, J) - subsystem
            //               A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
            //               D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
            //           for I = P, P - 1, ..., 1; J = 1, 2, ..., Q
            //
            pq = 0;
            scale = one;
            Rscale = zero;
            dsum = one;
            for (j = p + 2; j <= q; j = j + 1) {
                js = iwork[j - 1];
                je = iwork[(j + 1) - 1] - 1;
                nb = je - js + 1;
                for (i = p; i >= 1; i = i - 1) {
                    is = iwork[i - 1];
                    ie = iwork[(i + 1) - 1] - 1;
                    mb = ie - is + 1;
                    Ctgsy2(trans, ifunc, mb, nb, &a[(is - 1) + (is - 1) * lda], lda, &b[(js - 1) + (js - 1) * ldb], ldb, &c[(is - 1) + (js - 1) * ldc], ldc, &d[(is - 1) + (is - 1) * ldd], ldd, &e[(js - 1) + (js - 1) * lde], lde, &f[(is - 1) + (js - 1) * ldf], ldf, scaloc, dsum, Rscale, linfo);
                    if (linfo > 0) {
                        info = linfo;
                    }
                    pq += mb * nb;
                    if (scaloc != one) {
                        for (k = 1; k <= js - 1; k = k + 1) {
                            Cscal(m, COMPLEX(scaloc, zero), &c[(k - 1) * ldc], 1);
                            Cscal(m, COMPLEX(scaloc, zero), &f[(k - 1) * ldf], 1);
                        }
                        for (k = js; k <= je; k = k + 1) {
                            Cscal(is - 1, COMPLEX(scaloc, zero), &c[(k - 1) * ldc], 1);
                            Cscal(is - 1, COMPLEX(scaloc, zero), &f[(k - 1) * ldf], 1);
                        }
                        for (k = js; k <= je; k = k + 1) {
                            Cscal(m - ie, COMPLEX(scaloc, zero), &c[((ie + 1) - 1) + (k - 1) * ldc], 1);
                            Cscal(m - ie, COMPLEX(scaloc, zero), &f[((ie + 1) - 1) + (k - 1) * ldf], 1);
                        }
                        for (k = je + 1; k <= n; k = k + 1) {
                            Cscal(m, COMPLEX(scaloc, zero), &c[(k - 1) * ldc], 1);
                            Cscal(m, COMPLEX(scaloc, zero), &f[(k - 1) * ldf], 1);
                        }
                        scale = scale * scaloc;
                    }
                    //
                    //                 Substitute R(I,J) and L(I,J) into remaining equation.
                    //
                    if (i > 1) {
                        Cgemm("N", "N", is - 1, nb, mb, COMPLEX(-one, zero), &a[(is - 1) * lda], lda, &c[(is - 1) + (js - 1) * ldc], ldc, COMPLEX(one, zero), &c[(js - 1) * ldc], ldc);
                        Cgemm("N", "N", is - 1, nb, mb, COMPLEX(-one, zero), &d[(is - 1) * ldd], ldd, &c[(is - 1) + (js - 1) * ldc], ldc, COMPLEX(one, zero), &f[(js - 1) * ldf], ldf);
                    }
                    if (j < q) {
                        Cgemm("N", "N", mb, n - je, nb, COMPLEX(one, zero), &f[(is - 1) + (js - 1) * ldf], ldf, &b[(js - 1) + ((je + 1) - 1) * ldb], ldb, COMPLEX(one, zero), &c[(is - 1) + ((je + 1) - 1) * ldc], ldc);
                        Cgemm("N", "N", mb, n - je, nb, COMPLEX(one, zero), &f[(is - 1) + (js - 1) * ldf], ldf, &e[(js - 1) + ((je + 1) - 1) * lde], lde, COMPLEX(one, zero), &f[(is - 1) + ((je + 1) - 1) * ldf], ldf);
                    }
                }
            }
            if (Rscale != zero) {
                if (ijob == 1 || ijob == 3) {
                    dif = sqrt(castREAL(2 * m * n)) / (Rscale * sqrt(dsum));
                } else {
                    dif = sqrt(castREAL(pq)) / (Rscale * sqrt(dsum));
                }
            }
            if (isolve == 2 && iround == 1) {
                if (notran) {
                    ifunc = ijob;
                }
                scale2 = scale;
                Clacpy("F", m, n, c, ldc, work, m);
                Clacpy("F", m, n, f, ldf, &work[(m * n + 1) - 1], m);
                Claset("F", m, n, czero, czero, c, ldc);
                Claset("F", m, n, czero, czero, f, ldf);
            } else if (isolve == 2 && iround == 2) {
                Clacpy("F", m, n, work, m, c, ldc);
                Clacpy("F", m, n, &work[(m * n + 1) - 1], m, f, ldf);
                scale = scale2;
            }
        }
    } else {
        //
        //        Solve transposed (I, J)-subsystem
        //            A(I, I)**H * R(I, J) + D(I, I)**H * L(I, J) = C(I, J)
        //            R(I, J) * B(J, J)  + L(I, J) * E(J, J) = -F(I, J)
        //        for I = 1,2,..., P; J = Q, Q-1,..., 1
        //
        scale = one;
        for (i = 1; i <= p; i = i + 1) {
            is = iwork[i - 1];
            ie = iwork[(i + 1) - 1] - 1;
            mb = ie - is + 1;
            for (j = q; j >= p + 2; j = j - 1) {
                js = iwork[j - 1];
                je = iwork[(j + 1) - 1] - 1;
                nb = je - js + 1;
                Ctgsy2(trans, ifunc, mb, nb, &a[(is - 1) + (is - 1) * lda], lda, &b[(js - 1) + (js - 1) * ldb], ldb, &c[(is - 1) + (js - 1) * ldc], ldc, &d[(is - 1) + (is - 1) * ldd], ldd, &e[(js - 1) + (js - 1) * lde], lde, &f[(is - 1) + (js - 1) * ldf], ldf, scaloc, dsum, Rscale, linfo);
                if (linfo > 0) {
                    info = linfo;
                }
                if (scaloc != one) {
                    for (k = 1; k <= js - 1; k = k + 1) {
                        Cscal(m, COMPLEX(scaloc, zero), &c[(k - 1) * ldc], 1);
                        Cscal(m, COMPLEX(scaloc, zero), &f[(k - 1) * ldf], 1);
                    }
                    for (k = js; k <= je; k = k + 1) {
                        Cscal(is - 1, COMPLEX(scaloc, zero), &c[(k - 1) * ldc], 1);
                        Cscal(is - 1, COMPLEX(scaloc, zero), &f[(k - 1) * ldf], 1);
                    }
                    for (k = js; k <= je; k = k + 1) {
                        Cscal(m - ie, COMPLEX(scaloc, zero), &c[((ie + 1) - 1) + (k - 1) * ldc], 1);
                        Cscal(m - ie, COMPLEX(scaloc, zero), &f[((ie + 1) - 1) + (k - 1) * ldf], 1);
                    }
                    for (k = je + 1; k <= n; k = k + 1) {
                        Cscal(m, COMPLEX(scaloc, zero), &c[(k - 1) * ldc], 1);
                        Cscal(m, COMPLEX(scaloc, zero), &f[(k - 1) * ldf], 1);
                    }
                    scale = scale * scaloc;
                }
                //
                //              Substitute R(I,J) and L(I,J) into remaining equation.
                //
                if (j > p + 2) {
                    Cgemm("N", "C", mb, js - 1, nb, COMPLEX(one, zero), &c[(is - 1) + (js - 1) * ldc], ldc, &b[(js - 1) * ldb], ldb, COMPLEX(one, zero), &f[(is - 1)], ldf);
                    Cgemm("N", "C", mb, js - 1, nb, COMPLEX(one, zero), &f[(is - 1) + (js - 1) * ldf], ldf, &e[(js - 1) * lde], lde, COMPLEX(one, zero), &f[(is - 1)], ldf);
                }
                if (i < p) {
                    Cgemm("C", "N", m - ie, nb, mb, COMPLEX(-one, zero), &a[(is - 1) + ((ie + 1) - 1) * lda], lda, &c[(is - 1) + (js - 1) * ldc], ldc, COMPLEX(one, zero), &c[((ie + 1) - 1) + (js - 1) * ldc], ldc);
                    Cgemm("C", "N", m - ie, nb, mb, COMPLEX(-one, zero), &d[(is - 1) + ((ie + 1) - 1) * ldd], ldd, &f[(is - 1) + (js - 1) * ldf], ldf, COMPLEX(one, zero), &c[((ie + 1) - 1) + (js - 1) * ldc], ldc);
                }
            }
        }
    }
    //
    work[1 - 1] = lwmin;
    //
    //     End of Ctgsyl
    //
}
