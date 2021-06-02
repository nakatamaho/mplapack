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

void Ctfsm(const char *transr, const char *side, const char *uplo, const char *trans, const char *diag, INTEGER const m, INTEGER const n, COMPLEX const alpha, COMPLEX *a, COMPLEX *b, INTEGER const ldb) {
    //
    //     Test the input parameters.
    //
    INTEGER info = 0;
    bool normaltransr = Mlsame(transr, "N");
    bool lside = Mlsame(side, "L");
    bool lower = Mlsame(uplo, "L");
    bool notrans = Mlsame(trans, "N");
    if (!normaltransr && !Mlsame(transr, "C")) {
        info = -1;
    } else if (!lside && !Mlsame(side, "R")) {
        info = -2;
    } else if (!lower && !Mlsame(uplo, "U")) {
        info = -3;
    } else if (!notrans && !Mlsame(trans, "C")) {
        info = -4;
    } else if (!Mlsame(diag, "N") && !Mlsame(diag, "U")) {
        info = -5;
    } else if (m < 0) {
        info = -6;
    } else if (n < 0) {
        info = -7;
    } else if (ldb < max((INTEGER)1, m)) {
        info = -11;
    }
    if (info != 0) {
        Mxerbla("Ctfsm", -info);
        return;
    }
    //
    //     Quick return when ( (N.EQ.0).OR.(M.EQ.0) )
    //
    if ((m == 0) || (n == 0)) {
        return;
    }
    //
    //     Quick return when ALPHA.EQ.(0D+0,0D+0)
    //
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    INTEGER j = 0;
    INTEGER i = 0;
    if (alpha == czero) {
        for (j = 0; j <= n - 1; j = j + 1) {
            for (i = 0; i <= m - 1; i = i + 1) {
                b[(i - 1) + (j - 1) * ldb] = czero;
            }
        }
        return;
    }
    //
    bool misodd = false;
    INTEGER k = 0;
    INTEGER m2 = 0;
    INTEGER m1 = 0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    bool nisodd = false;
    INTEGER n2 = 0;
    INTEGER n1 = 0;
    if (lside) {
        //
        //        SIDE = 'L'
        //
        //        A is M-by-M.
        //        If M is odd, set NISODD = .TRUE., and M1 and M2.
        //        If M is even, NISODD = .FALSE., and M.
        //
        if (mod(m, 2) == 0) {
            misodd = false;
            k = m / 2;
        } else {
            misodd = true;
            if (lower) {
                m2 = m / 2;
                m1 = m - m2;
            } else {
                m1 = m / 2;
                m2 = m - m1;
            }
        }
        //
        if (misodd) {
            //
            //           SIDE = 'L' and N is odd
            //
            if (normaltransr) {
                //
                //              SIDE = 'L', N is odd, and TRANSR = 'N'
                //
                if (lower) {
                    //
                    //                 SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'L'
                    //
                    if (notrans) {
                        //
                        //                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and
                        //                    TRANS = 'N'
                        //
                        if (m == 1) {
                            Ctrsm("L", "L", "N", diag, m1, n, alpha, a, m, b, ldb);
                        } else {
                            Ctrsm("L", "L", "N", diag, m1, n, alpha, &a[0 - 1], m, b, ldb);
                            Cgemm("N", "N", m2, n, m1, -cone, &a[m1 - 1], m, b, ldb, alpha, &b[(m1 - 1) + (0 - 1) * ldb], ldb);
                            Ctrsm("L", "U", "C", diag, m2, n, cone, &a[m - 1], m, &b[(m1 - 1) + (0 - 1) * ldb], ldb);
                        }
                        //
                    } else {
                        //
                        //                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and
                        //                    TRANS = 'C'
                        //
                        if (m == 1) {
                            Ctrsm("L", "L", "C", diag, m1, n, alpha, &a[0 - 1], m, b, ldb);
                        } else {
                            Ctrsm("L", "U", "N", diag, m2, n, alpha, &a[m - 1], m, &b[(m1 - 1) + (0 - 1) * ldb], ldb);
                            Cgemm("C", "N", m1, n, m2, -cone, &a[m1 - 1], m, &b[(m1 - 1) + (0 - 1) * ldb], ldb, alpha, b, ldb);
                            Ctrsm("L", "L", "C", diag, m1, n, cone, &a[0 - 1], m, b, ldb);
                        }
                        //
                    }
                    //
                } else {
                    //
                    //                 SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'U'
                    //
                    if (!notrans) {
                        //
                        //                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and
                        //                    TRANS = 'N'
                        //
                        Ctrsm("L", "L", "N", diag, m1, n, alpha, &a[m2 - 1], m, b, ldb);
                        Cgemm("C", "N", m2, n, m1, -cone, &a[0 - 1], m, b, ldb, alpha, &b[(m1 - 1) + (0 - 1) * ldb], ldb);
                        Ctrsm("L", "U", "C", diag, m2, n, cone, &a[m1 - 1], m, &b[(m1 - 1) + (0 - 1) * ldb], ldb);
                        //
                    } else {
                        //
                        //                    SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and
                        //                    TRANS = 'C'
                        //
                        Ctrsm("L", "U", "N", diag, m2, n, alpha, &a[m1 - 1], m, &b[(m1 - 1) + (0 - 1) * ldb], ldb);
                        Cgemm("N", "N", m1, n, m2, -cone, &a[0 - 1], m, &b[(m1 - 1) + (0 - 1) * ldb], ldb, alpha, b, ldb);
                        Ctrsm("L", "L", "C", diag, m1, n, cone, &a[m2 - 1], m, b, ldb);
                        //
                    }
                    //
                }
                //
            } else {
                //
                //              SIDE = 'L', N is odd, and TRANSR = 'C'
                //
                if (lower) {
                    //
                    //                 SIDE  ='L', N is odd, TRANSR = 'C', and UPLO = 'L'
                    //
                    if (notrans) {
                        //
                        //                    SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'L', and
                        //                    TRANS = 'N'
                        //
                        if (m == 1) {
                            Ctrsm("L", "U", "C", diag, m1, n, alpha, &a[0 - 1], m1, b, ldb);
                        } else {
                            Ctrsm("L", "U", "C", diag, m1, n, alpha, &a[0 - 1], m1, b, ldb);
                            Cgemm("C", "N", m2, n, m1, -cone, &a[(m1 * m1) - 1], m1, b, ldb, alpha, &b[(m1 - 1) + (0 - 1) * ldb], ldb);
                            Ctrsm("L", "L", "N", diag, m2, n, cone, &a[1 - 1], m1, &b[(m1 - 1) + (0 - 1) * ldb], ldb);
                        }
                        //
                    } else {
                        //
                        //                    SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'L', and
                        //                    TRANS = 'C'
                        //
                        if (m == 1) {
                            Ctrsm("L", "U", "N", diag, m1, n, alpha, &a[0 - 1], m1, b, ldb);
                        } else {
                            Ctrsm("L", "L", "C", diag, m2, n, alpha, &a[1 - 1], m1, &b[(m1 - 1) + (0 - 1) * ldb], ldb);
                            Cgemm("N", "N", m1, n, m2, -cone, &a[(m1 * m1) - 1], m1, &b[(m1 - 1) + (0 - 1) * ldb], ldb, alpha, b, ldb);
                            Ctrsm("L", "U", "N", diag, m1, n, cone, &a[0 - 1], m1, b, ldb);
                        }
                        //
                    }
                    //
                } else {
                    //
                    //                 SIDE  ='L', N is odd, TRANSR = 'C', and UPLO = 'U'
                    //
                    if (!notrans) {
                        //
                        //                    SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'U', and
                        //                    TRANS = 'N'
                        //
                        Ctrsm("L", "U", "C", diag, m1, n, alpha, &a[(m2 * m2) - 1], m2, b, ldb);
                        Cgemm("N", "N", m2, n, m1, -cone, &a[0 - 1], m2, b, ldb, alpha, &b[(m1 - 1) + (0 - 1) * ldb], ldb);
                        Ctrsm("L", "L", "N", diag, m2, n, cone, &a[(m1 * m2) - 1], m2, &b[(m1 - 1) + (0 - 1) * ldb], ldb);
                        //
                    } else {
                        //
                        //                    SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'U', and
                        //                    TRANS = 'C'
                        //
                        Ctrsm("L", "L", "C", diag, m2, n, alpha, &a[(m1 * m2) - 1], m2, &b[(m1 - 1) + (0 - 1) * ldb], ldb);
                        Cgemm("C", "N", m1, n, m2, -cone, &a[0 - 1], m2, &b[(m1 - 1) + (0 - 1) * ldb], ldb, alpha, b, ldb);
                        Ctrsm("L", "U", "N", diag, m1, n, cone, &a[(m2 * m2) - 1], m2, b, ldb);
                        //
                    }
                    //
                }
                //
            }
            //
        } else {
            //
            //           SIDE = 'L' and N is even
            //
            if (normaltransr) {
                //
                //              SIDE = 'L', N is even, and TRANSR = 'N'
                //
                if (lower) {
                    //
                    //                 SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'L'
                    //
                    if (notrans) {
                        //
                        //                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L',
                        //                    and TRANS = 'N'
                        //
                        Ctrsm("L", "L", "N", diag, k, n, alpha, &a[1 - 1], m + 1, b, ldb);
                        Cgemm("N", "N", k, n, k, -cone, &a[(k + 1) - 1], m + 1, b, ldb, alpha, &b[(k - 1) + (0 - 1) * ldb], ldb);
                        Ctrsm("L", "U", "C", diag, k, n, cone, &a[0 - 1], m + 1, &b[(k - 1) + (0 - 1) * ldb], ldb);
                        //
                    } else {
                        //
                        //                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L',
                        //                    and TRANS = 'C'
                        //
                        Ctrsm("L", "U", "N", diag, k, n, alpha, &a[0 - 1], m + 1, &b[(k - 1) + (0 - 1) * ldb], ldb);
                        Cgemm("C", "N", k, n, k, -cone, &a[(k + 1) - 1], m + 1, &b[(k - 1) + (0 - 1) * ldb], ldb, alpha, b, ldb);
                        Ctrsm("L", "L", "C", diag, k, n, cone, &a[1 - 1], m + 1, b, ldb);
                        //
                    }
                    //
                } else {
                    //
                    //                 SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'U'
                    //
                    if (!notrans) {
                        //
                        //                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U',
                        //                    and TRANS = 'N'
                        //
                        Ctrsm("L", "L", "N", diag, k, n, alpha, &a[(k + 1) - 1], m + 1, b, ldb);
                        Cgemm("C", "N", k, n, k, -cone, &a[0 - 1], m + 1, b, ldb, alpha, &b[(k - 1) + (0 - 1) * ldb], ldb);
                        Ctrsm("L", "U", "C", diag, k, n, cone, &a[k - 1], m + 1, &b[(k - 1) + (0 - 1) * ldb], ldb);
                        //
                    } else {
                        //
                        //                    SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U',
                        //                    and TRANS = 'C'
                        Ctrsm("L", "U", "N", diag, k, n, alpha, &a[k - 1], m + 1, &b[(k - 1) + (0 - 1) * ldb], ldb);
                        Cgemm("N", "N", k, n, k, -cone, &a[0 - 1], m + 1, &b[(k - 1) + (0 - 1) * ldb], ldb, alpha, b, ldb);
                        Ctrsm("L", "L", "C", diag, k, n, cone, &a[(k + 1) - 1], m + 1, b, ldb);
                        //
                    }
                    //
                }
                //
            } else {
                //
                //              SIDE = 'L', N is even, and TRANSR = 'C'
                //
                if (lower) {
                    //
                    //                 SIDE  ='L', N is even, TRANSR = 'C', and UPLO = 'L'
                    //
                    if (notrans) {
                        //
                        //                    SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'L',
                        //                    and TRANS = 'N'
                        //
                        Ctrsm("L", "U", "C", diag, k, n, alpha, &a[k - 1], k, b, ldb);
                        Cgemm("C", "N", k, n, k, -cone, &a[(k * (k + 1)) - 1], k, b, ldb, alpha, &b[(k - 1) + (0 - 1) * ldb], ldb);
                        Ctrsm("L", "L", "N", diag, k, n, cone, &a[0 - 1], k, &b[(k - 1) + (0 - 1) * ldb], ldb);
                        //
                    } else {
                        //
                        //                    SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'L',
                        //                    and TRANS = 'C'
                        //
                        Ctrsm("L", "L", "C", diag, k, n, alpha, &a[0 - 1], k, &b[(k - 1) + (0 - 1) * ldb], ldb);
                        Cgemm("N", "N", k, n, k, -cone, &a[(k * (k + 1)) - 1], k, &b[(k - 1) + (0 - 1) * ldb], ldb, alpha, b, ldb);
                        Ctrsm("L", "U", "N", diag, k, n, cone, &a[k - 1], k, b, ldb);
                        //
                    }
                    //
                } else {
                    //
                    //                 SIDE  ='L', N is even, TRANSR = 'C', and UPLO = 'U'
                    //
                    if (!notrans) {
                        //
                        //                    SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'U',
                        //                    and TRANS = 'N'
                        //
                        Ctrsm("L", "U", "C", diag, k, n, alpha, &a[(k * (k + 1)) - 1], k, b, ldb);
                        Cgemm("N", "N", k, n, k, -cone, &a[0 - 1], k, b, ldb, alpha, &b[(k - 1) + (0 - 1) * ldb], ldb);
                        Ctrsm("L", "L", "N", diag, k, n, cone, &a[(k * k) - 1], k, &b[(k - 1) + (0 - 1) * ldb], ldb);
                        //
                    } else {
                        //
                        //                    SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'U',
                        //                    and TRANS = 'C'
                        //
                        Ctrsm("L", "L", "C", diag, k, n, alpha, &a[(k * k) - 1], k, &b[(k - 1) + (0 - 1) * ldb], ldb);
                        Cgemm("C", "N", k, n, k, -cone, &a[0 - 1], k, &b[(k - 1) + (0 - 1) * ldb], ldb, alpha, b, ldb);
                        Ctrsm("L", "U", "N", diag, k, n, cone, &a[(k * (k + 1)) - 1], k, b, ldb);
                        //
                    }
                    //
                }
                //
            }
            //
        }
        //
    } else {
        //
        //        SIDE = 'R'
        //
        //        A is N-by-N.
        //        If N is odd, set NISODD = .TRUE., and N1 and N2.
        //        If N is even, NISODD = .FALSE., and K.
        //
        if (mod(n, 2) == 0) {
            nisodd = false;
            k = n / 2;
        } else {
            nisodd = true;
            if (lower) {
                n2 = n / 2;
                n1 = n - n2;
            } else {
                n1 = n / 2;
                n2 = n - n1;
            }
        }
        //
        if (nisodd) {
            //
            //           SIDE = 'R' and N is odd
            //
            if (normaltransr) {
                //
                //              SIDE = 'R', N is odd, and TRANSR = 'N'
                //
                if (lower) {
                    //
                    //                 SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'L'
                    //
                    if (notrans) {
                        //
                        //                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and
                        //                    TRANS = 'N'
                        //
                        Ctrsm("R", "U", "C", diag, m, n2, alpha, &a[n - 1], n, &b[(0 - 1) + (n1 - 1) * ldb], ldb);
                        Cgemm("N", "N", m, n1, n2, -cone, &b[(0 - 1) + (n1 - 1) * ldb], ldb, &a[n1 - 1], n, alpha, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        Ctrsm("R", "L", "N", diag, m, n1, cone, &a[0 - 1], n, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        //
                    } else {
                        //
                        //                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and
                        //                    TRANS = 'C'
                        //
                        Ctrsm("R", "L", "C", diag, m, n1, alpha, &a[0 - 1], n, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        Cgemm("N", "C", m, n2, n1, -cone, &b[(0 - 1) + (0 - 1) * ldb], ldb, &a[n1 - 1], n, alpha, &b[(0 - 1) + (n1 - 1) * ldb], ldb);
                        Ctrsm("R", "U", "N", diag, m, n2, cone, &a[n - 1], n, &b[(0 - 1) + (n1 - 1) * ldb], ldb);
                        //
                    }
                    //
                } else {
                    //
                    //                 SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'U'
                    //
                    if (notrans) {
                        //
                        //                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and
                        //                    TRANS = 'N'
                        //
                        Ctrsm("R", "L", "C", diag, m, n1, alpha, &a[n2 - 1], n, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        Cgemm("N", "N", m, n2, n1, -cone, &b[(0 - 1) + (0 - 1) * ldb], ldb, &a[0 - 1], n, alpha, &b[(0 - 1) + (n1 - 1) * ldb], ldb);
                        Ctrsm("R", "U", "N", diag, m, n2, cone, &a[n1 - 1], n, &b[(0 - 1) + (n1 - 1) * ldb], ldb);
                        //
                    } else {
                        //
                        //                    SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and
                        //                    TRANS = 'C'
                        //
                        Ctrsm("R", "U", "C", diag, m, n2, alpha, &a[n1 - 1], n, &b[(0 - 1) + (n1 - 1) * ldb], ldb);
                        Cgemm("N", "C", m, n1, n2, -cone, &b[(0 - 1) + (n1 - 1) * ldb], ldb, &a[0 - 1], n, alpha, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        Ctrsm("R", "L", "N", diag, m, n1, cone, &a[n2 - 1], n, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        //
                    }
                    //
                }
                //
            } else {
                //
                //              SIDE = 'R', N is odd, and TRANSR = 'C'
                //
                if (lower) {
                    //
                    //                 SIDE  ='R', N is odd, TRANSR = 'C', and UPLO = 'L'
                    //
                    if (notrans) {
                        //
                        //                    SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'L', and
                        //                    TRANS = 'N'
                        //
                        Ctrsm("R", "L", "N", diag, m, n2, alpha, &a[1 - 1], n1, &b[(0 - 1) + (n1 - 1) * ldb], ldb);
                        Cgemm("N", "C", m, n1, n2, -cone, &b[(0 - 1) + (n1 - 1) * ldb], ldb, &a[(n1 * n1) - 1], n1, alpha, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        Ctrsm("R", "U", "C", diag, m, n1, cone, &a[0 - 1], n1, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        //
                    } else {
                        //
                        //                    SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'L', and
                        //                    TRANS = 'C'
                        //
                        Ctrsm("R", "U", "N", diag, m, n1, alpha, &a[0 - 1], n1, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        Cgemm("N", "N", m, n2, n1, -cone, &b[(0 - 1) + (0 - 1) * ldb], ldb, &a[(n1 * n1) - 1], n1, alpha, &b[(0 - 1) + (n1 - 1) * ldb], ldb);
                        Ctrsm("R", "L", "C", diag, m, n2, cone, &a[1 - 1], n1, &b[(0 - 1) + (n1 - 1) * ldb], ldb);
                        //
                    }
                    //
                } else {
                    //
                    //                 SIDE  ='R', N is odd, TRANSR = 'C', and UPLO = 'U'
                    //
                    if (notrans) {
                        //
                        //                    SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'U', and
                        //                    TRANS = 'N'
                        //
                        Ctrsm("R", "U", "N", diag, m, n1, alpha, &a[(n2 * n2) - 1], n2, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        Cgemm("N", "C", m, n2, n1, -cone, &b[(0 - 1) + (0 - 1) * ldb], ldb, &a[0 - 1], n2, alpha, &b[(0 - 1) + (n1 - 1) * ldb], ldb);
                        Ctrsm("R", "L", "C", diag, m, n2, cone, &a[(n1 * n2) - 1], n2, &b[(0 - 1) + (n1 - 1) * ldb], ldb);
                        //
                    } else {
                        //
                        //                    SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'U', and
                        //                    TRANS = 'C'
                        //
                        Ctrsm("R", "L", "N", diag, m, n2, alpha, &a[(n1 * n2) - 1], n2, &b[(0 - 1) + (n1 - 1) * ldb], ldb);
                        Cgemm("N", "N", m, n1, n2, -cone, &b[(0 - 1) + (n1 - 1) * ldb], ldb, &a[0 - 1], n2, alpha, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        Ctrsm("R", "U", "C", diag, m, n1, cone, &a[(n2 * n2) - 1], n2, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        //
                    }
                    //
                }
                //
            }
            //
        } else {
            //
            //           SIDE = 'R' and N is even
            //
            if (normaltransr) {
                //
                //              SIDE = 'R', N is even, and TRANSR = 'N'
                //
                if (lower) {
                    //
                    //                 SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'L'
                    //
                    if (notrans) {
                        //
                        //                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L',
                        //                    and TRANS = 'N'
                        //
                        Ctrsm("R", "U", "C", diag, m, k, alpha, &a[0 - 1], n + 1, &b[(0 - 1) + (k - 1) * ldb], ldb);
                        Cgemm("N", "N", m, k, k, -cone, &b[(0 - 1) + (k - 1) * ldb], ldb, &a[(k + 1) - 1], n + 1, alpha, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        Ctrsm("R", "L", "N", diag, m, k, cone, &a[1 - 1], n + 1, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        //
                    } else {
                        //
                        //                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L',
                        //                    and TRANS = 'C'
                        //
                        Ctrsm("R", "L", "C", diag, m, k, alpha, &a[1 - 1], n + 1, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        Cgemm("N", "C", m, k, k, -cone, &b[(0 - 1) + (0 - 1) * ldb], ldb, &a[(k + 1) - 1], n + 1, alpha, &b[(0 - 1) + (k - 1) * ldb], ldb);
                        Ctrsm("R", "U", "N", diag, m, k, cone, &a[0 - 1], n + 1, &b[(0 - 1) + (k - 1) * ldb], ldb);
                        //
                    }
                    //
                } else {
                    //
                    //                 SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'U'
                    //
                    if (notrans) {
                        //
                        //                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U',
                        //                    and TRANS = 'N'
                        //
                        Ctrsm("R", "L", "C", diag, m, k, alpha, &a[(k + 1) - 1], n + 1, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        Cgemm("N", "N", m, k, k, -cone, &b[(0 - 1) + (0 - 1) * ldb], ldb, &a[0 - 1], n + 1, alpha, &b[(0 - 1) + (k - 1) * ldb], ldb);
                        Ctrsm("R", "U", "N", diag, m, k, cone, &a[k - 1], n + 1, &b[(0 - 1) + (k - 1) * ldb], ldb);
                        //
                    } else {
                        //
                        //                    SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U',
                        //                    and TRANS = 'C'
                        //
                        Ctrsm("R", "U", "C", diag, m, k, alpha, &a[k - 1], n + 1, &b[(0 - 1) + (k - 1) * ldb], ldb);
                        Cgemm("N", "C", m, k, k, -cone, &b[(0 - 1) + (k - 1) * ldb], ldb, &a[0 - 1], n + 1, alpha, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        Ctrsm("R", "L", "N", diag, m, k, cone, &a[(k + 1) - 1], n + 1, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        //
                    }
                    //
                }
                //
            } else {
                //
                //              SIDE = 'R', N is even, and TRANSR = 'C'
                //
                if (lower) {
                    //
                    //                 SIDE  ='R', N is even, TRANSR = 'C', and UPLO = 'L'
                    //
                    if (notrans) {
                        //
                        //                    SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'L',
                        //                    and TRANS = 'N'
                        //
                        Ctrsm("R", "L", "N", diag, m, k, alpha, &a[0 - 1], k, &b[(0 - 1) + (k - 1) * ldb], ldb);
                        Cgemm("N", "C", m, k, k, -cone, &b[(0 - 1) + (k - 1) * ldb], ldb, &a[((k + 1) * k) - 1], k, alpha, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        Ctrsm("R", "U", "C", diag, m, k, cone, &a[k - 1], k, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        //
                    } else {
                        //
                        //                    SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'L',
                        //                    and TRANS = 'C'
                        //
                        Ctrsm("R", "U", "N", diag, m, k, alpha, &a[k - 1], k, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        Cgemm("N", "N", m, k, k, -cone, &b[(0 - 1) + (0 - 1) * ldb], ldb, &a[((k + 1) * k) - 1], k, alpha, &b[(0 - 1) + (k - 1) * ldb], ldb);
                        Ctrsm("R", "L", "C", diag, m, k, cone, &a[0 - 1], k, &b[(0 - 1) + (k - 1) * ldb], ldb);
                        //
                    }
                    //
                } else {
                    //
                    //                 SIDE  ='R', N is even, TRANSR = 'C', and UPLO = 'U'
                    //
                    if (notrans) {
                        //
                        //                    SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'U',
                        //                    and TRANS = 'N'
                        //
                        Ctrsm("R", "U", "N", diag, m, k, alpha, &a[((k + 1) * k) - 1], k, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        Cgemm("N", "C", m, k, k, -cone, &b[(0 - 1) + (0 - 1) * ldb], ldb, &a[0 - 1], k, alpha, &b[(0 - 1) + (k - 1) * ldb], ldb);
                        Ctrsm("R", "L", "C", diag, m, k, cone, &a[(k * k) - 1], k, &b[(0 - 1) + (k - 1) * ldb], ldb);
                        //
                    } else {
                        //
                        //                    SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'U',
                        //                    and TRANS = 'C'
                        //
                        Ctrsm("R", "L", "N", diag, m, k, alpha, &a[(k * k) - 1], k, &b[(0 - 1) + (k - 1) * ldb], ldb);
                        Cgemm("N", "N", m, k, k, -cone, &b[(0 - 1) + (k - 1) * ldb], ldb, &a[0 - 1], k, alpha, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        Ctrsm("R", "U", "C", diag, m, k, cone, &a[((k + 1) * k) - 1], k, &b[(0 - 1) + (0 - 1) * ldb], ldb);
                        //
                    }
                    //
                }
                //
            }
            //
        }
    }
    //
    //     End of Ctfsm
    //
}
