/*
 * Copyright (c) 2021
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

void Cgejsv(const char *joba, const char *jobu, const char *jobv, const char *jobr, const char *jobt, const char *jobp, INTEGER const &m, INTEGER const &n, COMPLEX *a, INTEGER const &lda, REAL *sva, COMPLEX *u, INTEGER const &ldu, COMPLEX *v, INTEGER const &ldv, COMPLEX *cwork, INTEGER const &lwork, REAL *rwork, INTEGER const &lrwork, arr_ref<INTEGER> iwork, INTEGER &info) {
    bool lsvec = false;
    bool jracc = false;
    bool rsvec = false;
    bool rowpiv = false;
    bool l2rank = false;
    bool l2aber = false;
    bool errest = false;
    bool l2tran = false;
    bool l2kill = false;
    bool defr = false;
    bool l2pert = false;
    bool lquery = false;
    INTEGER lwqp3 = 0;
    INTEGER lwqrf = 0;
    INTEGER lwlqf = 0;
    INTEGER lwunmlq = 0;
    INTEGER lwunmqr = 0;
    INTEGER lwunmqrm = 0;
    INTEGER lwcon = 0;
    INTEGER lwsvdj = 0;
    INTEGER lwsvdjv = 0;
    INTEGER lrwqp3 = 0;
    INTEGER lrwcon = 0;
    INTEGER lrwsvdj = 0;
    arr_1d<1, COMPLEX> cdummy(fill0);
    arr_1d<1, REAL> rdummy(fill0);
    INTEGER ierr = 0;
    INTEGER lwrk_Cgeqp3 = 0;
    INTEGER lwrk_Cgeqrf = 0;
    INTEGER lwrk_Cgelqf = 0;
    INTEGER minwrk = 0;
    INTEGER optwrk = 0;
    INTEGER miniwrk = 0;
    INTEGER lwrk_Cgesvj = 0;
    INTEGER minrwrk = 0;
    INTEGER lwrk_Cunmlq = 0;
    INTEGER lwrk_Cunmqrm = 0;
    INTEGER lwrk_Cunmqr = 0;
    INTEGER lwrk_Cgeqp3n = 0;
    INTEGER lwrk_Cgesvju = 0;
    INTEGER lwrk_Cgesvjv = 0;
    INTEGER n1 = 0;
    REAL epsln = 0.0;
    REAL sfmin = 0.0;
    REAL small = 0.0;
    REAL big = 0.0;
    const REAL one = 1.0;
    REAL scalem = 0.0;
    bool noscal = false;
    bool goscal = false;
    INTEGER p = 0;
    const REAL zero = 0.0;
    REAL aapp = 0.0;
    REAL aaqq = 0.0;
    const COMPLEX czero = (0.0, 0.0);
    const COMPLEX cone = (1.0, 0.0);
    INTEGER warning = 0;
    bool transp = false;
    REAL aatmax = 0.0;
    REAL aatmin = 0.0;
    REAL xsc = 0.0;
    REAL temp1 = 0.0;
    REAL entra = 0.0;
    REAL entrat = 0.0;
    REAL big1 = 0.0;
    INTEGER q = 0;
    COMPLEX ctemp = 0.0;
    bool kill = false;
    REAL uscal1 = 0.0;
    REAL uscal2 = 0.0;
    INTEGER iwoff = 0;
    INTEGER nr = 0;
    bool almort = false;
    REAL maxprj = 0.0;
    REAL sconda = 0.0;
    REAL condr1 = 0.0;
    REAL condr2 = 0.0;
    INTEGER numrank = 0;
    REAL cond_ok = 0.0;
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
    //  ===========================================================================
    //
    //     .. Local Parameters ..
    //     ..
    //     .. Local Scalars ..
    //
    //     ..
    //     .. Local Arrays
    //
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //
    //     ..
    //
    //     Test the input arguments
    //
    lsvec = Mlsame(jobu, "U") || Mlsame(jobu, "F");
    jracc = Mlsame(jobv, "J");
    rsvec = Mlsame(jobv, "V") || jracc;
    rowpiv = Mlsame(joba, "F") || Mlsame(joba, "G");
    l2rank = Mlsame(joba, "R");
    l2aber = Mlsame(joba, "A");
    errest = Mlsame(joba, "E") || Mlsame(joba, "G");
    l2tran = Mlsame(jobt, "T") && (m == n);
    l2kill = Mlsame(jobr, "R");
    defr = Mlsame(jobr, "N");
    l2pert = Mlsame(jobp, "P");
    //
    lquery = (lwork == -1) || (lrwork == -1);
    //
    if (!(rowpiv || l2rank || l2aber || errest || Mlsame(joba, "C"))) {
        info = -1;
    } else if (!(lsvec || Mlsame(jobu, "N") || (Mlsame(jobu, "W") && rsvec && l2tran))) {
        info = -2;
    } else if (!(rsvec || Mlsame(jobv, "N") || (Mlsame(jobv, "W") && lsvec && l2tran))) {
        info = -3;
    } else if (!(l2kill || defr)) {
        info = -4;
    } else if (!(Mlsame(jobt, "T") || Mlsame(jobt, "N"))) {
        info = -5;
    } else if (!(l2pert || Mlsame(jobp, "N"))) {
        info = -6;
    } else if (m < 0) {
        info = -7;
    } else if ((n < 0) || (n > m)) {
        info = -8;
    } else if (lda < m) {
        info = -10;
    } else if (lsvec && (ldu < m)) {
        info = -13;
    } else if (rsvec && (ldv < n)) {
        info = -15;
    } else {
        //        #:)
        info = 0;
    }
    //
    if (info == 0) {
        //         .. compute the minimal and the optimal workspace lengths
        //         [[The expressions for computing the minimal and the optimal
        //         values of LCWORK, LRWORK are written with a lot of redundancy and
        //         can be simplified. However, this verbose form is useful for
        //         maINTEGERenance and modifications of the code.]]
        //
        //        .. minimal workspace length for Cgeqp3 of an M x N matrix,
        //         Cgeqrf of an N x N matrix, Cgelqf of an N x N matrix,
        //         Cunmlq for computing N x N matrix, Cunmqr for computing N x N
        //         matrix, Cunmqr for computing M x N matrix, respectively.
        lwqp3 = n + 1;
        lwqrf = max((INTEGER)1, n);
        lwlqf = max((INTEGER)1, n);
        lwunmlq = max((INTEGER)1, n);
        lwunmqr = max((INTEGER)1, n);
        lwunmqrm = max((INTEGER)1, m);
        //        .. minimal workspace length for Cpocon of an N x N matrix
        lwcon = 2 * n;
        //        .. minimal workspace length for Cgesvj of an N x N matrix,
        //         without and with explicit accumulation of Jacobi rotations
        lwsvdj = max(2 * n, 1);
        lwsvdjv = max(2 * n, 1);
        //         .. minimal REAL workspace length for Cgeqp3, Cpocon, Cgesvj
        lrwqp3 = 2 * n;
        lrwcon = n;
        lrwsvdj = n;
        if (lquery) {
            Cgeqp3(m, n, a, lda, iwork, cdummy, cdummy, -1, rdummy, ierr);
            lwrk_Cgeqp3 = cdummy[1 - 1];
            Cgeqrf(n, n, a, lda, cdummy, cdummy, -1, ierr);
            lwrk_Cgeqrf = cdummy[1 - 1];
            Cgelqf(n, n, a, lda, cdummy, cdummy, -1, ierr);
            lwrk_Cgelqf = cdummy[1 - 1];
        }
        minwrk = 2;
        optwrk = 2;
        miniwrk = n;
        if (!(lsvec || rsvec)) {
            //             .. minimal and optimal sizes of the complex workspace if
            //             only the singular values are requested
            if (errest) {
                minwrk = max(n + lwqp3, pow2(n) + lwcon, n + lwqrf, lwsvdj);
            } else {
                minwrk = max(n + lwqp3, n + lwqrf, lwsvdj);
            }
            if (lquery) {
                Cgesvj("L", "N", "N", n, n, a, lda, sva, n, v, ldv, cdummy, -1, rdummy, -1, ierr);
                lwrk_Cgesvj = cdummy[1 - 1];
                if (errest) {
                    optwrk = max(n + lwrk_Cgeqp3, pow2(n) + lwcon, n + lwrk_Cgeqrf, lwrk_Cgesvj);
                } else {
                    optwrk = max(n + lwrk_Cgeqp3, n + lwrk_Cgeqrf, lwrk_Cgesvj);
                }
            }
            if (l2tran || rowpiv) {
                if (errest) {
                    minrwrk = max(7, 2 * m, lrwqp3, lrwcon, lrwsvdj);
                } else {
                    minrwrk = max(7, 2 * m, lrwqp3, lrwsvdj);
                }
            } else {
                if (errest) {
                    minrwrk = max(7, lrwqp3, lrwcon, lrwsvdj);
                } else {
                    minrwrk = max(7, lrwqp3, lrwsvdj);
                }
            }
            if (rowpiv || l2tran) {
                miniwrk += m;
            }
        } else if (rsvec && (!lsvec)) {
            //            .. minimal and optimal sizes of the complex workspace if the
            //            singular values and the right singular vectors are requested
            if (errest) {
                minwrk = max(n + lwqp3, lwcon, lwsvdj, n + lwlqf, 2 * n + lwqrf, n + lwsvdj, n + lwunmlq);
            } else {
                minwrk = max(n + lwqp3, lwsvdj, n + lwlqf, 2 * n + lwqrf, n + lwsvdj, n + lwunmlq);
            }
            if (lquery) {
                Cgesvj("L", "U", "N", n, n, u, ldu, sva, n, a, lda, cdummy, -1, rdummy, -1, ierr);
                lwrk_Cgesvj = cdummy[1 - 1];
                Cunmlq("L", "C", n, n, n, a, lda, cdummy, v, ldv, cdummy, -1, ierr);
                lwrk_Cunmlq = cdummy[1 - 1];
                if (errest) {
                    optwrk = max(n + lwrk_Cgeqp3, lwcon, lwrk_Cgesvj, n + lwrk_Cgelqf, 2 * n + lwrk_Cgeqrf, n + lwrk_Cgesvj, n + lwrk_Cunmlq);
                } else {
                    optwrk = max(n + lwrk_Cgeqp3, lwrk_Cgesvj, n + lwrk_Cgelqf, 2 * n + lwrk_Cgeqrf, n + lwrk_Cgesvj, n + lwrk_Cunmlq);
                }
            }
            if (l2tran || rowpiv) {
                if (errest) {
                    minrwrk = max(7, 2 * m, lrwqp3, lrwsvdj, lrwcon);
                } else {
                    minrwrk = max(7, 2 * m, lrwqp3, lrwsvdj);
                }
            } else {
                if (errest) {
                    minrwrk = max(7, lrwqp3, lrwsvdj, lrwcon);
                } else {
                    minrwrk = max(7, lrwqp3, lrwsvdj);
                }
            }
            if (rowpiv || l2tran) {
                miniwrk += m;
            }
        } else if (lsvec && (!rsvec)) {
            //            .. minimal and optimal sizes of the complex workspace if the
            //            singular values and the left singular vectors are requested
            if (errest) {
                minwrk = n + max(lwqp3, lwcon, n + lwqrf, lwsvdj, lwunmqrm);
            } else {
                minwrk = n + max(lwqp3, n + lwqrf, lwsvdj, lwunmqrm);
            }
            if (lquery) {
                Cgesvj("L", "U", "N", n, n, u, ldu, sva, n, a, lda, cdummy, -1, rdummy, -1, ierr);
                lwrk_Cgesvj = cdummy[1 - 1];
                Cunmqr("L", "N", m, n, n, a, lda, cdummy, u, ldu, cdummy, -1, ierr);
                lwrk_Cunmqrm = cdummy[1 - 1];
                if (errest) {
                    optwrk = n + max(lwrk_Cgeqp3, lwcon, n + lwrk_Cgeqrf, lwrk_Cgesvj, lwrk_Cunmqrm);
                } else {
                    optwrk = n + max(lwrk_Cgeqp3, n + lwrk_Cgeqrf, lwrk_Cgesvj, lwrk_Cunmqrm);
                }
            }
            if (l2tran || rowpiv) {
                if (errest) {
                    minrwrk = max(7, 2 * m, lrwqp3, lrwsvdj, lrwcon);
                } else {
                    minrwrk = max(7, 2 * m, lrwqp3, lrwsvdj);
                }
            } else {
                if (errest) {
                    minrwrk = max(7, lrwqp3, lrwsvdj, lrwcon);
                } else {
                    minrwrk = max(7, lrwqp3, lrwsvdj);
                }
            }
            if (rowpiv || l2tran) {
                miniwrk += m;
            }
        } else {
            //            .. minimal and optimal sizes of the complex workspace if the
            //            full SVD is requested
            if (!jracc) {
                if (errest) {
                    minwrk = max(n + lwqp3, n + lwcon, 2 * n + pow2(n) + lwcon, 2 * n + lwqrf, 2 * n + lwqp3, 2 * n + pow2(n) + n + lwlqf, 2 * n + pow2(n) + n + pow2(n) + lwcon, 2 * n + pow2(n) + n + lwsvdj, 2 * n + pow2(n) + n + lwsvdjv, 2 * n + pow2(n) + n + lwunmqr, 2 * n + pow2(n) + n + lwunmlq, n + pow2(n) + lwsvdj, n + lwunmqrm);
                } else {
                    minwrk = max(n + lwqp3, 2 * n + pow2(n) + lwcon, 2 * n + lwqrf, 2 * n + lwqp3, 2 * n + pow2(n) + n + lwlqf, 2 * n + pow2(n) + n + pow2(n) + lwcon, 2 * n + pow2(n) + n + lwsvdj, 2 * n + pow2(n) + n + lwsvdjv, 2 * n + pow2(n) + n + lwunmqr, 2 * n + pow2(n) + n + lwunmlq, n + pow2(n) + lwsvdj, n + lwunmqrm);
                }
                miniwrk += n;
                if (rowpiv || l2tran) {
                    miniwrk += m;
                }
            } else {
                if (errest) {
                    minwrk = max(n + lwqp3, n + lwcon, 2 * n + lwqrf, 2 * n + pow2(n) + lwsvdjv, 2 * n + pow2(n) + n + lwunmqr, n + lwunmqrm);
                } else {
                    minwrk = max(n + lwqp3, 2 * n + lwqrf, 2 * n + pow2(n) + lwsvdjv, 2 * n + pow2(n) + n + lwunmqr, n + lwunmqrm);
                }
                if (rowpiv || l2tran) {
                    miniwrk += m;
                }
            }
            if (lquery) {
                Cunmqr("L", "N", m, n, n, a, lda, cdummy, u, ldu, cdummy, -1, ierr);
                lwrk_Cunmqrm = cdummy[1 - 1];
                Cunmqr("L", "N", n, n, n, a, lda, cdummy, u, ldu, cdummy, -1, ierr);
                lwrk_Cunmqr = cdummy[1 - 1];
                if (!jracc) {
                    Cgeqp3(n, n, a, lda, iwork, cdummy, cdummy, -1, rdummy, ierr);
                    lwrk_Cgeqp3n = cdummy[1 - 1];
                    Cgesvj("L", "U", "N", n, n, u, ldu, sva, n, v, ldv, cdummy, -1, rdummy, -1, ierr);
                    lwrk_Cgesvj = cdummy[1 - 1];
                    Cgesvj("U", "U", "N", n, n, u, ldu, sva, n, v, ldv, cdummy, -1, rdummy, -1, ierr);
                    lwrk_Cgesvju = cdummy[1 - 1];
                    Cgesvj("L", "U", "V", n, n, u, ldu, sva, n, v, ldv, cdummy, -1, rdummy, -1, ierr);
                    lwrk_Cgesvjv = cdummy[1 - 1];
                    Cunmlq("L", "C", n, n, n, a, lda, cdummy, v, ldv, cdummy, -1, ierr);
                    lwrk_Cunmlq = cdummy[1 - 1];
                    if (errest) {
                        optwrk = max(n + lwrk_Cgeqp3, n + lwcon, 2 * n + pow2(n) + lwcon, 2 * n + lwrk_Cgeqrf, 2 * n + lwrk_Cgeqp3n, 2 * n + pow2(n) + n + lwrk_Cgelqf, 2 * n + pow2(n) + n + pow2(n) + lwcon, 2 * n + pow2(n) + n + lwrk_Cgesvj, 2 * n + pow2(n) + n + lwrk_Cgesvjv, 2 * n + pow2(n) + n + lwrk_Cunmqr, 2 * n + pow2(n) + n + lwrk_Cunmlq, n + pow2(n) + lwrk_Cgesvju, n + lwrk_Cunmqrm);
                    } else {
                        optwrk = max(n + lwrk_Cgeqp3, 2 * n + pow2(n) + lwcon, 2 * n + lwrk_Cgeqrf, 2 * n + lwrk_Cgeqp3n, 2 * n + pow2(n) + n + lwrk_Cgelqf, 2 * n + pow2(n) + n + pow2(n) + lwcon, 2 * n + pow2(n) + n + lwrk_Cgesvj, 2 * n + pow2(n) + n + lwrk_Cgesvjv, 2 * n + pow2(n) + n + lwrk_Cunmqr, 2 * n + pow2(n) + n + lwrk_Cunmlq, n + pow2(n) + lwrk_Cgesvju, n + lwrk_Cunmqrm);
                    }
                } else {
                    Cgesvj("L", "U", "V", n, n, u, ldu, sva, n, v, ldv, cdummy, -1, rdummy, -1, ierr);
                    lwrk_Cgesvjv = cdummy[1 - 1];
                    Cunmqr("L", "N", n, n, n, cdummy, n, cdummy, v, ldv, cdummy, -1, ierr);
                    lwrk_Cunmqr = cdummy[1 - 1];
                    Cunmqr("L", "N", m, n, n, a, lda, cdummy, u, ldu, cdummy, -1, ierr);
                    lwrk_Cunmqrm = cdummy[1 - 1];
                    if (errest) {
                        optwrk = max(n + lwrk_Cgeqp3, n + lwcon, 2 * n + lwrk_Cgeqrf, 2 * n + pow2(n), 2 * n + pow2(n) + lwrk_Cgesvjv, 2 * n + pow2(n) + n + lwrk_Cunmqr, n + lwrk_Cunmqrm);
                    } else {
                        optwrk = max(n + lwrk_Cgeqp3, 2 * n + lwrk_Cgeqrf, 2 * n + pow2(n), 2 * n + pow2(n) + lwrk_Cgesvjv, 2 * n + pow2(n) + n + lwrk_Cunmqr, n + lwrk_Cunmqrm);
                    }
                }
            }
            if (l2tran || rowpiv) {
                minrwrk = max(7, 2 * m, lrwqp3, lrwsvdj, lrwcon);
            } else {
                minrwrk = max(7, lrwqp3, lrwsvdj, lrwcon);
            }
        }
        minwrk = max(2, minwrk);
        optwrk = max(minwrk, optwrk);
        if (lwork < minwrk && (!lquery)) {
            info = -17;
        }
        if (lrwork < minrwrk && (!lquery)) {
            info = -19;
        }
    }
    //
    if (info != 0) {
        //       #:(
        Mxerbla("Cgejsv", -info);
        return;
    } else if (lquery) {
        cwork[1 - 1] = optwrk;
        cwork[2 - 1] = minwrk;
        rwork[1 - 1] = minrwrk;
        iwork[1 - 1] = max(4, miniwrk);
        return;
    }
    //
    //     Quick return for void matrix (Y3K safe)
    // #:)
    if ((m == 0) || (n == 0)) {
        iwork[(4 - 1) * ldiwork] = 0;
        rwork[(7 - 1) * ldrwork] = 0;
        return;
    }
    //
    //     Determine whether the matrix U should be M x N or M x M
    //
    if (lsvec) {
        n1 = n;
        if (Mlsame(jobu, "F")) {
            n1 = m;
        }
    }
    //
    //     Set numerical parameters
    //
    //!    NOTE: Make sure DLAMCH() does not fail on the target architecture.
    //
    epsln = dlamch("Epsilon");
    sfmin = dlamch("SafeMinimum");
    small = sfmin / epsln;
    big = dlamch("O");
    //     BIG   = ONE / SFMIN
    //
    //     Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N
    //
    //(!)  If necessary, scale SVA() to protect the largest norm from
    //     overflow. It is possible that this scaling pushes the smallest
    //     column norm left from the underflow threshold (extreme case).
    //
    scalem = one / sqrt(m.real() * n.real());
    noscal = true;
    goscal = true;
    for (p = 1; p <= n; p = p + 1) {
        aapp = zero;
        aaqq = one;
        Classq(m, a[(p - 1) * lda], 1, aapp, aaqq);
        if (aapp > big) {
            info = -9;
            Mxerbla("Cgejsv", -info);
            return;
        }
        aaqq = sqrt(aaqq);
        if ((aapp < (big / aaqq)) && noscal) {
            sva[p - 1] = aapp * aaqq;
        } else {
            noscal = false;
            sva[p - 1] = aapp * (aaqq * scalem);
            if (goscal) {
                goscal = false;
                Rscal(p - 1, scalem, sva, 1);
            }
        }
    }
    //
    if (noscal) {
        scalem = one;
    }
    //
    aapp = zero;
    aaqq = big;
    for (p = 1; p <= n; p = p + 1) {
        aapp = max(aapp, sva[p - 1]);
        if (sva[p - 1] != zero) {
            aaqq = min(aaqq, sva[p - 1]);
        }
    }
    //
    //     Quick return for zero M x N matrix
    // #:)
    if (aapp == zero) {
        if (lsvec) {
            Claset("G", m, n1, czero, cone, u, ldu);
        }
        if (rsvec) {
            Claset("G", n, n, czero, cone, v, ldv);
        }
        rwork[1 - 1] = one;
        rwork[2 - 1] = one;
        if (errest) {
            rwork[3 - 1] = one;
        }
        if (lsvec && rsvec) {
            rwork[4 - 1] = one;
            rwork[5 - 1] = one;
        }
        if (l2tran) {
            rwork[6 - 1] = zero;
            rwork[7 - 1] = zero;
        }
        iwork[1 - 1] = 0;
        iwork[2 - 1] = 0;
        iwork[3 - 1] = 0;
        iwork[4 - 1] = -1;
        return;
    }
    //
    //     Issue warning if denormalized column norms detected. Override the
    //     high relative accuracy request. Issue licence to kill nonzero columns
    //     (set them to zero) whose norm is less than sigma_max / BIG (roughly).
    // #:(
    warning = 0;
    if (aaqq <= sfmin) {
        l2rank = true;
        l2kill = true;
        warning = 1;
    }
    //
    //     Quick return for one-column matrix
    // #:)
    if (n == 1) {
        //
        if (lsvec) {
            Clascl("G", 0, 0, sva[1 - 1], scalem, m, 1, a[(1 - 1)], lda, ierr);
            Clacpy("A", m, 1, a, lda, u, ldu);
            //           computing all M left singular vectors of the M x 1 matrix
            if (n1 != n) {
                Cgeqrf(m, n, u, ldu, cwork, cwork[(n + 1) - 1], lwork - n, ierr);
                Cungqr(m, n1, 1, u, ldu, cwork, cwork[(n + 1) - 1], lwork - n, ierr);
                Ccopy(m, a[(1 - 1)], 1, u[(1 - 1)], 1);
            }
        }
        if (rsvec) {
            v[(1 - 1)] = cone;
        }
        if (sva[1 - 1] < (big * scalem)) {
            sva[1 - 1] = sva[1 - 1] / scalem;
            scalem = one;
        }
        rwork[1 - 1] = one / scalem;
        rwork[2 - 1] = one;
        if (sva[1 - 1] != zero) {
            iwork[1 - 1] = 1;
            if ((sva[1 - 1] / scalem) >= sfmin) {
                iwork[2 - 1] = 1;
            } else {
                iwork[2 - 1] = 0;
            }
        } else {
            iwork[1 - 1] = 0;
            iwork[2 - 1] = 0;
        }
        iwork[3 - 1] = 0;
        iwork[4 - 1] = -1;
        if (errest) {
            rwork[3 - 1] = one;
        }
        if (lsvec && rsvec) {
            rwork[4 - 1] = one;
            rwork[5 - 1] = one;
        }
        if (l2tran) {
            rwork[6 - 1] = zero;
            rwork[7 - 1] = zero;
        }
        return;
        //
    }
    //
    transp = false;
    //
    aatmax = -one;
    aatmin = big;
    if (rowpiv || l2tran) {
        //
        //     Compute the row norms, needed to determine row pivoting sequence
        //     (in the case of heavily row weighted A, row pivoting is strongly
        //     advised) and to collect information needed to compare the
        //     structures of A * A^* and A^* * A (in the case L2TRAN.EQ..TRUE.).
        //
        if (l2tran) {
            for (p = 1; p <= m; p = p + 1) {
                xsc = zero;
                temp1 = one;
                Classq(n, a[(p - 1)], lda, xsc, temp1);
                //              Classq gets both the ell_2 and the ell_infinity norm
                //              in one pass through the vector
                rwork[(m + p) - 1] = xsc * scalem;
                rwork[p - 1] = xsc * (scalem * sqrt(temp1));
                aatmax = max(aatmax, rwork[p - 1]);
                if (rwork[p - 1] != zero) {
                    aatmin = min(aatmin, rwork[p - 1]);
                }
            }
        } else {
            for (p = 1; p <= m; p = p + 1) {
                rwork[(m + p) - 1] = scalem * abs(a[(p - 1) + (iCamax[(n - 1) + (a[(p - 1)] - 1) * ldiCamax] - 1) * lda]);
                aatmax = max(aatmax, rwork[(m + p) - 1]);
                aatmin = min(aatmin, rwork[(m + p) - 1]);
            }
        }
        //
    }
    //
    //     For square matrix A try to determine whether A^*  would be better
    //     input for the preconditioned Jacobi SVD, with faster convergence.
    //     The decision is based on an O(N) function of the vector of column
    //     and row norms of A, based on the Shannon entropy. This should give
    //     the right choice in most cases when the difference actually matters.
    //     It may fail and pick the slower converging side.
    //
    entra = zero;
    entrat = zero;
    if (l2tran) {
        //
        xsc = zero;
        temp1 = one;
        Rlassq(n, sva, 1, xsc, temp1);
        temp1 = one / temp1;
        //
        entra = zero;
        for (p = 1; p <= n; p = p + 1) {
            big1 = (pow2((sva[p - 1] / xsc))) * temp1;
            if (big1 != zero) {
                entra += big1 * dlog[big1 - 1];
            }
        }
        entra = -entra / dlog[n.real() - 1];
        //
        //        Now, SVA().^2/Trace(A^* * A) is a poINTEGER in the probability simplex.
        //        It is derived from the diagonal of  A^* * A.  Do the same with the
        //        diagonal of A * A^*, compute the entropy of the corresponding
        //        probability distribution. Note that A * A^* and A^* * A have the
        //        same trace.
        //
        entrat = zero;
        for (p = 1; p <= m; p = p + 1) {
            big1 = (pow2((rwork[p - 1] / xsc))) * temp1;
            if (big1 != zero) {
                entrat += big1 * dlog[big1 - 1];
            }
        }
        entrat = -entrat / dlog[m.real() - 1];
        //
        //        Analyze the entropies and decide A or A^*. Smaller entropy
        //        usually means better input for the algorithm.
        //
        transp = (entrat < entra);
        //
        //        If A^* is better than A, take the adjoINTEGER of A. This is allowed
        //        only for square matrices, M=N.
        if (transp) {
            //           In an optimal implementation, this trivial transpose
            //           should be replaced with faster transpose.
            for (p = 1; p <= n - 1; p = p + 1) {
                a[(p - 1) + (p - 1) * lda] = conjg[a[(p - 1) + (p - 1) * lda] - 1];
                for (q = p + 1; q <= n; q = q + 1) {
                    ctemp = conjg[a[(q - 1) + (p - 1) * lda] - 1];
                    a[(q - 1) + (p - 1) * lda] = conjg[a[(p - 1) + (q - 1) * lda] - 1];
                    a[(p - 1) + (q - 1) * lda] = ctemp;
                }
            }
            a[(n - 1) + (n - 1) * lda] = conjg[a[(n - 1) + (n - 1) * lda] - 1];
            for (p = 1; p <= n; p = p + 1) {
                rwork[(m + p) - 1] = sva[p - 1];
                sva[p - 1] = rwork[p - 1];
                //              previously computed row 2-norms are now column 2-norms
                //              of the transposed matrix
            }
            temp1 = aapp;
            aapp = aatmax;
            aatmax = temp1;
            temp1 = aaqq;
            aaqq = aatmin;
            aatmin = temp1;
            kill = lsvec;
            lsvec = rsvec;
            rsvec = kill;
            if (lsvec) {
                n1 = n;
            }
            //
            rowpiv = true;
        }
        //
    }
    //     END IF L2TRAN
    //
    //     Scale the matrix so that its maximal singular value remains less
    //     than SQRT(BIG) -- the matrix is scaled so that its maximal column
    //     has Euclidean norm equal to SQRT(BIG/N). The only reason to keep
    //     SQRT(BIG) instead of BIG is the fact that Cgejsv uses LAPACK and
    //     BLAS routines that, in some implementations, are not capable of
    //     working in the full INTEGERerval [SFMIN,BIG] and that they may provoke
    //     overflows in the INTEGERermediate results. If the singular values spread
    //     from SFMIN to BIG, then Cgesvj will compute them. So, in that case,
    //     one should use Cgesvj instead of Cgejsv.
    //     >> change in the April 2016 update: allow bigger range, i.e. the
    //     largest column is allowed up to BIG/N and Cgesvj will do the rest.
    big1 = sqrt(big);
    temp1 = sqrt(big / n.real());
    //      TEMP1  = BIG/DBLE(N)
    //
    Rlascl("G", 0, 0, aapp, temp1, n, 1, sva, n, ierr);
    if (aaqq > (aapp * sfmin)) {
        aaqq = (aaqq / aapp) * temp1;
    } else {
        aaqq = (aaqq * temp1) / aapp;
    }
    temp1 = temp1 * scalem;
    Clascl("G", 0, 0, aapp, temp1, m, n, a, lda, ierr);
    //
    //     To undo scaling at the end of this procedure, multiply the
    //     computed singular values with USCAL2 / USCAL1.
    //
    uscal1 = temp1;
    uscal2 = aapp;
    //
    if (l2kill) {
        //        L2KILL enforces computation of nonzero singular values in
        //        the restricted range of condition number of the initial A,
        //        sigma_max(A) / sigma_min(A) approx. SQRT(BIG)/SQRT(SFMIN).
        xsc = sqrt(sfmin);
    } else {
        xsc = small;
        //
        //        Now, if the condition number of A is too big,
        //        sigma_max(A) / sigma_min(A) .GT. SQRT(BIG/N) * EPSLN / SFMIN,
        //        as a precaution measure, the full SVD is computed using Cgesvj
        //        with accumulated Jacobi rotations. This provides numerically
        //        more robust computation, at the cost of slightly increased run
        //        time. Depending on the concrete implementation of BLAS and LAPACK
        //        (i.e. how they behave in presence of extreme ill-conditioning) the
        //        implementor may decide to remove this switch.
        if ((aaqq < sqrt(sfmin)) && lsvec && rsvec) {
            jracc = true;
        }
        //
    }
    if (aaqq < xsc) {
        for (p = 1; p <= n; p = p + 1) {
            if (sva[p - 1] < xsc) {
                Claset("A", m, 1, czero, czero, a[(p - 1) * lda], lda);
                sva[p - 1] = zero;
            }
        }
    }
    //
    //     Preconditioning using QR factorization with pivoting
    //
    if (rowpiv) {
        //        Optional row permutation (Bjoerck row pivoting):
        //        A result by Cox and Higham shows that the Bjoerck's
        //        row pivoting combined with standard column pivoting
        //        has similar effect as Powell-Reid complete pivoting.
        //        The ell-infinity norms of A are made nonincreasing.
        if ((lsvec && rsvec) && !(jracc)) {
            iwoff = 2 * n;
        } else {
            iwoff = n;
        }
        for (p = 1; p <= m - 1; p = p + 1) {
            q = iRamax[((m - p + 1) - 1) + ((rwork[(m + p) - 1]) - 1) * ldiRamax] + p - 1;
            iwork[(iwoff + p) - 1] = q;
            if (p != q) {
                temp1 = rwork[(m + p) - 1];
                rwork[(m + p) - 1] = rwork[(m + q) - 1];
                rwork[(m + q) - 1] = temp1;
            }
        }
        Claswp(n, a, lda, 1, m - 1, iwork[(iwoff + 1) - 1], 1);
    }
    //
    //     End of the preparation phase (scaling, optional sorting and
    //     transposing, optional flushing of small columns).
    //
    //     Preconditioning
    //
    //     If the full SVD is needed, the right singular vectors are computed
    //     from a matrix equation, and for that we need theoretical analysis
    //     of the Businger-Golub pivoting. So we use Cgeqp3 as the first RR QRF.
    //     In all other cases the first RR QRF can be chosen by other criteria
    //     (eg speed by replacing global with restricted window pivoting, such
    //     as in xGEQPX from TOMS # 782). Good results will be obtained using
    //     xGEQPX with properly (!) chosen numerical parameters.
    //     Any improvement of Cgeqp3 improves overall performance of Cgejsv.
    //
    //     A * P1 = Q1 * [ R1^* 0]^*:
    for (p = 1; p <= n; p = p + 1) {
        //        .. all columns are free columns
        iwork[p - 1] = 0;
    }
    Cgeqp3(m, n, a, lda, iwork, cwork, cwork[(n + 1) - 1], lwork - n, rwork, ierr);
    //
    //     The upper triangular matrix R1 from the first QRF is inspected for
    //     rank deficiency and possibilities for deflation, or possible
    //     ill-conditioning. Depending on the user specified flag L2RANK,
    //     the procedure explores possibilities to reduce the numerical
    //     rank by inspecting the computed upper triangular factor. If
    //     L2RANK or L2ABER are up, then Cgejsv will compute the SVD of
    //     A + dA, where ||dA|| <= f(M,N)*EPSLN.
    //
    nr = 1;
    if (l2aber) {
        //        Standard absolute error bound suffices. All sigma_i with
        //        sigma_i < N*EPSLN*||A|| are flushed to zero. This is an
        //        aggressive enforcement of lower numerical rank by INTEGERroducing a
        //        backward error of the order of N*EPSLN*||A||.
        temp1 = sqrt(n.real()) * epsln;
        for (p = 2; p <= n; p = p + 1) {
            if (abs(a[(p - 1) + (p - 1) * lda]) >= (temp1 * abs(a[(1 - 1)]))) {
                nr++;
            } else {
                goto statement_3002;
            }
        }
    statement_3002:;
    } else if (l2rank) {
        //        .. similarly as above, only slightly more gentle (less aggressive).
        //        Sudden drop on the diagonal of R1 is used as the criterion for
        //        close-to-rank-deficient.
        temp1 = sqrt(sfmin);
        for (p = 2; p <= n; p = p + 1) {
            if ((abs(a[(p - 1) + (p - 1) * lda]) < (epsln * abs(a[((p - 1) - 1) + ((p - 1) - 1) * lda]))) || (abs(a[(p - 1) + (p - 1) * lda]) < small) || (l2kill && (abs(a[(p - 1) + (p - 1) * lda]) < temp1))) {
                goto statement_3402;
            }
            nr++;
        }
    statement_3402:;
        //
    } else {
        //        The goal is high relative accuracy. However, if the matrix
        //        has high scaled condition number the relative accuracy is in
        //        general not feasible. Later on, a condition number estimator
        //        will be deployed to estimate the scaled condition number.
        //        Here we just remove the underflowed part of the triangular
        //        factor. This prevents the situation in which the code is
        //        working hard to get the accuracy not warranted by the data.
        temp1 = sqrt(sfmin);
        for (p = 2; p <= n; p = p + 1) {
            if ((abs(a[(p - 1) + (p - 1) * lda]) < small) || (l2kill && (abs(a[(p - 1) + (p - 1) * lda]) < temp1))) {
                goto statement_3302;
            }
            nr++;
        }
    statement_3302:;
        //
    }
    //
    almort = false;
    if (nr == n) {
        maxprj = one;
        for (p = 2; p <= n; p = p + 1) {
            temp1 = abs(a[(p - 1) + (p - 1) * lda]) / sva[iwork[p - 1] - 1];
            maxprj = min(maxprj, temp1);
        }
        if (pow2(maxprj) >= one - n.real() * epsln) {
            almort = true;
        }
    }
    //
    sconda = -one;
    condr1 = -one;
    condr2 = -one;
    //
    if (errest) {
        if (n == nr) {
            if (rsvec) {
                //              .. V is available as workspace
                Clacpy("U", n, n, a, lda, v, ldv);
                for (p = 1; p <= n; p = p + 1) {
                    temp1 = sva[iwork[p - 1] - 1];
                    CRscal(p, one / temp1, v[(p - 1) * ldv], 1);
                }
                if (lsvec) {
                    Cpocon("U", n, v, ldv, one, temp1, cwork[(n + 1) - 1], rwork, ierr);
                } else {
                    Cpocon("U", n, v, ldv, one, temp1, cwork, rwork, ierr);
                }
                //
            } else if (lsvec) {
                //              .. U is available as workspace
                Clacpy("U", n, n, a, lda, u, ldu);
                for (p = 1; p <= n; p = p + 1) {
                    temp1 = sva[iwork[p - 1] - 1];
                    CRscal(p, one / temp1, u[(p - 1) * ldu], 1);
                }
                Cpocon("U", n, u, ldu, one, temp1, cwork[(n + 1) - 1], rwork, ierr);
            } else {
                Clacpy("U", n, n, a, lda, cwork, n);
                //[]            CALL Clacpy( 'U', N, N, A, LDA, CWORK(N+1), N )
                //              Change: here index shifted by N to the left, CWORK(1:N)
                //              not needed for SIGMA only computation
                for (p = 1; p <= n; p = p + 1) {
                    temp1 = sva[iwork[p - 1] - 1];
                    //[]               CALL CRscal( p, ONE/TEMP1, CWORK(N+(p-1)*N+1), 1 )
                    CRscal(p, one / temp1, cwork[((p - 1) * n + 1) - 1], 1);
                }
                //           .. the columns of R are scaled to have unit Euclidean lengths.
                //[]               CALL Cpocon( 'U', N, CWORK(N+1), N, ONE, TEMP1,
                //[]     $              CWORK(N+N*N+1), RWORK, IERR )
                Cpocon("U", n, cwork, n, one, temp1, cwork[(n * n + 1) - 1], rwork, ierr);
                //
            }
            if (temp1 != zero) {
                sconda = one / sqrt(temp1);
            } else {
                sconda = -one;
            }
            //           SCONDA is an estimate of SQRT(||(R^* * R)^(-1)||_1).
            //           N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
        } else {
            sconda = -one;
        }
    }
    //
    l2pert = l2pert && (abs(a[(1 - 1)] / a[(nr - 1) + (nr - 1) * lda]) > sqrt(big1));
    //     If there is no violent scaling, artificial perturbation is not needed.
    //
    //     Phase 3:
    //
    if (!(rsvec || lsvec)) {
        //
        //         Singular Values only
        //
        //         .. transpose A(1:NR,1:N)
        for (p = 1; p <= min(n - 1, nr); p = p + 1) {
            Ccopy(n - p, a[(p - 1) + ((p + 1) - 1) * lda], lda, a[((p + 1) - 1) + (p - 1) * lda], 1);
            Clacgv(n - p + 1, a[(p - 1) + (p - 1) * lda], 1);
        }
        if (nr == n) {
            a[(n - 1) + (n - 1) * lda] = conjg[a[(n - 1) + (n - 1) * lda] - 1];
        }
        //
        //        The following two DO-loops INTEGERroduce small relative perturbation
        //        INTEGERo the strict upper triangle of the lower triangular matrix.
        //        Small entries below the main diagonal are also changed.
        //        This modification is useful if the computing environment does not
        //        provide/allow FLUSH TO ZERO underflow, for it prevents many
        //        annoying denormalized numbers in case of strongly scaled matrices.
        //        The perturbation is structured so that it does not INTEGERroduce any
        //        new perturbation of the singular values, and it does not destroy
        //        the job done by the preconditioner.
        //        The licence for this perturbation is in the variable L2PERT, which
        //        should be .FALSE. if FLUSH TO ZERO underflow is active.
        //
        if (!almort) {
            //
            if (l2pert) {
                //              XSC = SQRT(SMALL)
                xsc = epsln / n.real();
                for (q = 1; q <= nr; q = q + 1) {
                    ctemp = COMPLEX(xsc * abs(a[(q - 1) + (q - 1) * lda]), zero);
                    for (p = 1; p <= n; p = p + 1) {
                        //     $                     A(p,q) = TEMP1 * ( A(p,q) / ABS(A(p,q)) )
                        if (((p > q) && (abs(a[(p - 1) + (q - 1) * lda]) <= temp1)) || (p < q)) {
                            a[(p - 1) + (q - 1) * lda] = ctemp;
                        }
                    }
                }
            } else {
                Claset("U", nr - 1, nr - 1, czero, czero, a[(2 - 1) * lda], lda);
            }
            //
            //            .. second preconditioning using the QR factorization
            //
            Cgeqrf(n, nr, a, lda, cwork, cwork[(n + 1) - 1], lwork - n, ierr);
            //
            //           .. and transpose upper to lower triangular
            for (p = 1; p <= nr - 1; p = p + 1) {
                Ccopy(nr - p, a[(p - 1) + ((p + 1) - 1) * lda], lda, a[((p + 1) - 1) + (p - 1) * lda], 1);
                Clacgv(nr - p + 1, a[(p - 1) + (p - 1) * lda], 1);
            }
            //
        }
        //
        //           Row-cyclic Jacobi SVD algorithm with column pivoting
        //
        //           .. again some perturbation (a "background noise") is added
        //           to drown denormals
        if (l2pert) {
            //              XSC = SQRT(SMALL)
            xsc = epsln / n.real();
            for (q = 1; q <= nr; q = q + 1) {
                ctemp = COMPLEX(xsc * abs(a[(q - 1) + (q - 1) * lda]), zero);
                for (p = 1; p <= nr; p = p + 1) {
                    //     $                   A(p,q) = TEMP1 * ( A(p,q) / ABS(A(p,q)) )
                    if (((p > q) && (abs(a[(p - 1) + (q - 1) * lda]) <= temp1)) || (p < q)) {
                        a[(p - 1) + (q - 1) * lda] = ctemp;
                    }
                }
            }
        } else {
            Claset("U", nr - 1, nr - 1, czero, czero, a[(2 - 1) * lda], lda);
        }
        //
        //           .. and one-sided Jacobi rotations are started on a lower
        //           triangular matrix (plus perturbation which is ignored in
        //           the part which destroys triangular form (confusing?!))
        //
        Cgesvj("L", "N", "N", nr, nr, a, lda, sva, n, v, ldv, cwork, lwork, rwork, lrwork, info);
        //
        scalem = rwork[1 - 1];
        numrank = nINTEGER[rwork[2 - 1] - 1];
        //
    } else if ((rsvec && (!lsvec) && (!jracc)) || (jracc && (!lsvec) && (nr != n))) {
        //
        //        -> Singular Values and Right Singular Vectors <-
        //
        if (almort) {
            //
            //           .. in this case NR equals N
            for (p = 1; p <= nr; p = p + 1) {
                Ccopy(n - p + 1, a[(p - 1) + (p - 1) * lda], lda, v[(p - 1) + (p - 1) * ldv], 1);
                Clacgv(n - p + 1, v[(p - 1) + (p - 1) * ldv], 1);
            }
            Claset("U", nr - 1, nr - 1, czero, czero, v[(2 - 1) * ldv], ldv);
            //
            Cgesvj("L", "U", "N", n, nr, v, ldv, sva, nr, a, lda, cwork, lwork, rwork, lrwork, info);
            scalem = rwork[1 - 1];
            numrank = nINTEGER[rwork[2 - 1] - 1];
            //
        } else {
            //
            //        .. two more QR factorizations ( one QRF is not enough, two require
            //        accumulated product of Jacobi rotations, three are perfect )
            //
            Claset("L", nr - 1, nr - 1, czero, czero, a[(2 - 1)], lda);
            Cgelqf(nr, n, a, lda, cwork, cwork[(n + 1) - 1], lwork - n, ierr);
            Clacpy("L", nr, nr, a, lda, v, ldv);
            Claset("U", nr - 1, nr - 1, czero, czero, v[(2 - 1) * ldv], ldv);
            Cgeqrf(nr, nr, v, ldv, cwork[(n + 1) - 1], cwork[(2 * n + 1) - 1], lwork - 2 * n, ierr);
            for (p = 1; p <= nr; p = p + 1) {
                Ccopy(nr - p + 1, v[(p - 1) + (p - 1) * ldv], ldv, v[(p - 1) + (p - 1) * ldv], 1);
                Clacgv(nr - p + 1, v[(p - 1) + (p - 1) * ldv], 1);
            }
            Claset("U", nr - 1, nr - 1, czero, czero, v[(2 - 1) * ldv], ldv);
            //
            Cgesvj("L", "U", "N", nr, nr, v, ldv, sva, nr, u, ldu, cwork[(n + 1) - 1], lwork - n, rwork, lrwork, info);
            scalem = rwork[1 - 1];
            numrank = nINTEGER[rwork[2 - 1] - 1];
            if (nr < n) {
                Claset("A", n - nr, nr, czero, czero, v[((nr + 1) - 1)], ldv);
                Claset("A", nr, n - nr, czero, czero, v[((nr + 1) - 1) * ldv], ldv);
                Claset("A", n - nr, n - nr, czero, cone, v[((nr + 1) - 1) + ((nr + 1) - 1) * ldv], ldv);
            }
            //
            Cunmlq("L", "C", n, n, nr, a, lda, cwork, v, ldv, cwork[(n + 1) - 1], lwork - n, ierr);
            //
        }
        //         .. permute the rows of V
        //         DO 8991 p = 1, N
        //            CALL Ccopy( N, V(p,1), LDV, A(IWORK(p),1), LDA )
        // 8991    CONTINUE
        //         CALL Clacpy( 'All', N, N, A, LDA, V, LDV )
        Clapmr(false, n, n, v, ldv, iwork);
        //
        if (transp) {
            Clacpy("A", n, n, v, ldv, u, ldu);
        }
        //
    } else if (jracc && (!lsvec) && (nr == n)) {
        //
        Claset("L", n - 1, n - 1, czero, czero, a[(2 - 1)], lda);
        //
        Cgesvj("U", "N", "V", n, n, a, lda, sva, n, v, ldv, cwork, lwork, rwork, lrwork, info);
        scalem = rwork[1 - 1];
        numrank = nINTEGER[rwork[2 - 1] - 1];
        Clapmr(false, n, n, v, ldv, iwork);
        //
    } else if (lsvec && (!rsvec)) {
        //
        //        .. Singular Values and Left Singular Vectors                 ..
        //
        //        .. second preconditioning step to avoid need to accumulate
        //        Jacobi rotations in the Jacobi iterations.
        for (p = 1; p <= nr; p = p + 1) {
            Ccopy(n - p + 1, a[(p - 1) + (p - 1) * lda], lda, u[(p - 1) + (p - 1) * ldu], 1);
            Clacgv(n - p + 1, u[(p - 1) + (p - 1) * ldu], 1);
        }
        Claset("U", nr - 1, nr - 1, czero, czero, u[(2 - 1) * ldu], ldu);
        //
        Cgeqrf(n, nr, u, ldu, cwork[(n + 1) - 1], cwork[(2 * n + 1) - 1], lwork - 2 * n, ierr);
        //
        for (p = 1; p <= nr - 1; p = p + 1) {
            Ccopy(nr - p, u[(p - 1) + ((p + 1) - 1) * ldu], ldu, u[((p + 1) - 1) + (p - 1) * ldu], 1);
            Clacgv(n - p + 1, u[(p - 1) + (p - 1) * ldu], 1);
        }
        Claset("U", nr - 1, nr - 1, czero, czero, u[(2 - 1) * ldu], ldu);
        //
        Cgesvj("L", "U", "N", nr, nr, u, ldu, sva, nr, a, lda, cwork[(n + 1) - 1], lwork - n, rwork, lrwork, info);
        scalem = rwork[1 - 1];
        numrank = nINTEGER[rwork[2 - 1] - 1];
        //
        if (nr < m) {
            Claset("A", m - nr, nr, czero, czero, u[((nr + 1) - 1)], ldu);
            if (nr < n1) {
                Claset("A", nr, n1 - nr, czero, czero, u[((nr + 1) - 1) * ldu], ldu);
                Claset("A", m - nr, n1 - nr, czero, cone, u[((nr + 1) - 1) + ((nr + 1) - 1) * ldu], ldu);
            }
        }
        //
        Cunmqr("L", "N", m, n1, n, a, lda, cwork, u, ldu, cwork[(n + 1) - 1], lwork - n, ierr);
        //
        if (rowpiv) {
            Claswp(n1, u, ldu, 1, m - 1, iwork[(iwoff + 1) - 1], -1);
        }
        //
        for (p = 1; p <= n1; p = p + 1) {
            xsc = one / RCnrm2[(m - 1) + (u[(p - 1) * ldu] - 1) * ldRCnrm2];
            CRscal(m, xsc, u[(p - 1) * ldu], 1);
        }
        //
        if (transp) {
            Clacpy("A", n, n, u, ldu, v, ldv);
        }
        //
    } else {
        //
        //        .. Full SVD ..
        //
        if (!jracc) {
            //
            if (!almort) {
                //
                //           Second Preconditioning Step (QRF [with pivoting])
                //           Note that the composition of TRANSPOSE, QRF and TRANSPOSE is
                //           equivalent to an LQF CALL. Since in many libraries the QRF
                //           seems to be better optimized than the LQF, we do explicit
                //           transpose and use the QRF. This is subject to changes in an
                //           optimized implementation of Cgejsv.
                //
                for (p = 1; p <= nr; p = p + 1) {
                    Ccopy(n - p + 1, a[(p - 1) + (p - 1) * lda], lda, v[(p - 1) + (p - 1) * ldv], 1);
                    Clacgv(n - p + 1, v[(p - 1) + (p - 1) * ldv], 1);
                }
                //
                //           .. the following two loops perturb small entries to avoid
                //           denormals in the second QR factorization, where they are
                //           as good as zeros. This is done to avoid painfully slow
                //           computation with denormals. The relative size of the perturbation
                //           is a parameter that can be changed by the implementer.
                //           This perturbation device will be obsolete on machines with
                //           properly implemented arithmetic.
                //           To switch it off, set L2PERT=.FALSE. To remove it from  the
                //           code, remove the action under L2PERT=.TRUE., leave the ELSE part.
                //           The following two loops should be blocked and fused with the
                //           transposed copy above.
                //
                if (l2pert) {
                    xsc = sqrt(small);
                    for (q = 1; q <= nr; q = q + 1) {
                        ctemp = COMPLEX(xsc * abs(v[(q - 1) + (q - 1) * ldv]), zero);
                        for (p = 1; p <= n; p = p + 1) {
                            //     $                   V(p,q) = TEMP1 * ( V(p,q) / ABS(V(p,q)) )
                            if ((p > q) && (abs(v[(p - 1) + (q - 1) * ldv]) <= temp1) || (p < q)) {
                                v[(p - 1) + (q - 1) * ldv] = ctemp;
                            }
                            if (p < q) {
                                v[(p - 1) + (q - 1) * ldv] = -v[(p - 1) + (q - 1) * ldv];
                            }
                        }
                    }
                } else {
                    Claset("U", nr - 1, nr - 1, czero, czero, v[(2 - 1) * ldv], ldv);
                }
                //
                //           Estimate the row scaled condition number of R1
                //           (If R1 is rectangular, N > NR, then the condition number
                //           of the leading NR x NR submatrix is estimated.)
                //
                Clacpy("L", nr, nr, v, ldv, cwork[(2 * n + 1) - 1], nr);
                for (p = 1; p <= nr; p = p + 1) {
                    temp1 = RCnrm2[((nr - p + 1) - 1) + ((cwork[(2 * n + (p - 1) * nr + p) - 1]) - 1) * ldRCnrm2];
                    CRscal(nr - p + 1, one / temp1, cwork[(2 * n + (p - 1) * nr + p) - 1], 1);
                }
                Cpocon("L", nr, cwork[(2 * n + 1) - 1], nr, one, temp1, cwork[(2 * n + nr * nr + 1) - 1], rwork, ierr);
                condr1 = one / sqrt(temp1);
                //           .. here need a second opinion on the condition number
                //           .. then assume worst case scenario
                //           R1 is OK for inverse <=> CONDR1 .LT. DBLE(N)
                //           more conservative    <=> CONDR1 .LT. SQRT(DBLE(N))
                //
                cond_ok = sqrt(sqrt(nr.real()));
                //[TP]       COND_OK is a tuning parameter.
                //
                if (condr1 < cond_ok) {
                    //              .. the second QRF without pivoting. Note: in an optimized
                    //              implementation, this QRF should be implemented as the QRF
                    //              of a lower triangular matrix.
                    //              R1^* = Q2 * R2
                    Cgeqrf(n, nr, v, ldv, cwork[(n + 1) - 1], cwork[(2 * n + 1) - 1], lwork - 2 * n, ierr);
                    //
                    if (l2pert) {
                        xsc = sqrt(small) / epsln;
                        for (p = 2; p <= nr; p = p + 1) {
                            for (q = 1; q <= p - 1; q = q + 1) {
                                ctemp = COMPLEX(xsc * min(abs(v[(p - 1) + (p - 1) * ldv]), abs(v[(q - 1) + (q - 1) * ldv])), zero);
                                //     $                     V(q,p) = TEMP1 * ( V(q,p) / ABS(V(q,p)) )
                                if (abs(v[(q - 1) + (p - 1) * ldv]) <= temp1) {
                                    v[(q - 1) + (p - 1) * ldv] = ctemp;
                                }
                            }
                        }
                    }
                    //
                    if (nr != n) {
                        Clacpy("A", n, nr, v, ldv, cwork[(2 * n + 1) - 1], n);
                    }
                    //              .. save ...
                    //
                    //           .. this transposed copy should be better than naive
                    for (p = 1; p <= nr - 1; p = p + 1) {
                        Ccopy(nr - p, v[(p - 1) + ((p + 1) - 1) * ldv], ldv, v[((p + 1) - 1) + (p - 1) * ldv], 1);
                        Clacgv(nr - p + 1, v[(p - 1) + (p - 1) * ldv], 1);
                    }
                    v[(nr - 1) + (nr - 1) * ldv] = conjg[v[(nr - 1) + (nr - 1) * ldv] - 1];
                    //
                    condr2 = condr1;
                    //
                } else {
                    //
                    //              .. ill-conditioned case: second QRF with pivoting
                    //              Note that windowed pivoting would be equally good
                    //              numerically, and more run-time efficient. So, in
                    //              an optimal implementation, the next call to Cgeqp3
                    //              should be replaced with eg. CALL ZGEQPX (ACM TOMS #782)
                    //              with properly (carefully) chosen parameters.
                    //
                    //              R1^* * P2 = Q2 * R2
                    for (p = 1; p <= nr; p = p + 1) {
                        iwork[(n + p) - 1] = 0;
                    }
                    Cgeqp3(n, nr, v, ldv, iwork[(n + 1) - 1], cwork[(n + 1) - 1], cwork[(2 * n + 1) - 1], lwork - 2 * n, rwork, ierr);
                    //*               CALL Cgeqrf( N, NR, V, LDV, CWORK(N+1), CWORK(2*N+1),
                    //*     $              LWORK-2*N, IERR )
                    if (l2pert) {
                        xsc = sqrt(small);
                        for (p = 2; p <= nr; p = p + 1) {
                            for (q = 1; q <= p - 1; q = q + 1) {
                                ctemp = COMPLEX(xsc * min(abs(v[(p - 1) + (p - 1) * ldv]), abs(v[(q - 1) + (q - 1) * ldv])), zero);
                                //     $                     V(q,p) = TEMP1 * ( V(q,p) / ABS(V(q,p)) )
                                if (abs(v[(q - 1) + (p - 1) * ldv]) <= temp1) {
                                    v[(q - 1) + (p - 1) * ldv] = ctemp;
                                }
                            }
                        }
                    }
                    //
                    Clacpy("A", n, nr, v, ldv, cwork[(2 * n + 1) - 1], n);
                    //
                    if (l2pert) {
                        xsc = sqrt(small);
                        for (p = 2; p <= nr; p = p + 1) {
                            for (q = 1; q <= p - 1; q = q + 1) {
                                ctemp = COMPLEX(xsc * min(abs(v[(p - 1) + (p - 1) * ldv]), abs(v[(q - 1) + (q - 1) * ldv])), zero);
                                //                        V(p,q) = - TEMP1*( V(q,p) / ABS(V(q,p)) )
                                v[(p - 1) + (q - 1) * ldv] = -ctemp;
                            }
                        }
                    } else {
                        Claset("L", nr - 1, nr - 1, czero, czero, v[(2 - 1)], ldv);
                    }
                    //              Now, compute R2 = L3 * Q3, the LQ factorization.
                    Cgelqf(nr, nr, v, ldv, cwork[(2 * n + n * nr + 1) - 1], cwork[(2 * n + n * nr + nr + 1) - 1], lwork - 2 * n - n * nr - nr, ierr);
                    //              .. and estimate the condition number
                    Clacpy("L", nr, nr, v, ldv, cwork[(2 * n + n * nr + nr + 1) - 1], nr);
                    for (p = 1; p <= nr; p = p + 1) {
                        temp1 = RCnrm2[(p - 1) + ((cwork[(2 * n + n * nr + nr + p) - 1]) - 1) * ldRCnrm2];
                        CRscal(p, one / temp1, cwork[(2 * n + n * nr + nr + p) - 1], nr);
                    }
                    Cpocon("L", nr, cwork[(2 * n + n * nr + nr + 1) - 1], nr, one, temp1, cwork[(2 * n + n * nr + nr + nr * nr + 1) - 1], rwork, ierr);
                    condr2 = one / sqrt(temp1);
                    //
                    if (condr2 >= cond_ok) {
                        //                 .. save the Householder vectors used for Q3
                        //                 (this overwrites the copy of R2, as it will not be
                        //                 needed in this branch, but it does not overwritte the
                        //                 Huseholder vectors of Q2.).
                        Clacpy("U", nr, nr, v, ldv, cwork[(2 * n + 1) - 1], n);
                        //                 .. and the rest of the information on Q3 is in
                        //                 WORK(2*N+N*NR+1:2*N+N*NR+N)
                    }
                    //
                }
                //
                if (l2pert) {
                    xsc = sqrt(small);
                    for (q = 2; q <= nr; q = q + 1) {
                        ctemp = xsc * v[(q - 1) + (q - 1) * ldv];
                        for (p = 1; p <= q - 1; p = p + 1) {
                            //                     V(p,q) = - TEMP1*( V(p,q) / ABS(V(p,q)) )
                            v[(p - 1) + (q - 1) * ldv] = -ctemp;
                        }
                    }
                } else {
                    Claset("U", nr - 1, nr - 1, czero, czero, v[(2 - 1) * ldv], ldv);
                }
                //
                //        Second preconditioning finished; continue with Jacobi SVD
                //        The input matrix is lower trinagular.
                //
                //        Recover the right singular vectors as solution of a well
                //        conditioned triangular matrix equation.
                //
                if (condr1 < cond_ok) {
                    //
                    Cgesvj("L", "U", "N", nr, nr, v, ldv, sva, nr, u, ldu, cwork[(2 * n + n * nr + nr + 1) - 1], lwork - 2 * n - n * nr - nr, rwork, lrwork, info);
                    scalem = rwork[1 - 1];
                    numrank = nINTEGER[rwork[2 - 1] - 1];
                    for (p = 1; p <= nr; p = p + 1) {
                        Ccopy(nr, v[(p - 1) * ldv], 1, u[(p - 1) * ldu], 1);
                        CRscal(nr, sva[p - 1], v[(p - 1) * ldv], 1);
                    }
                    //
                    //        .. pick the right matrix equation and solve it
                    //
                    if (nr == n) {
                        // :))             .. best case, R1 is inverted. The solution of this matrix
                        //                 equation is Q2*V2 = the product of the Jacobi rotations
                        //                 used in Cgesvj, premultiplied with the orthogonal matrix
                        //                 from the second QR factorization.
                        Ctrsm("L", "U", "N", "N", nr, nr, cone, a, lda, v, ldv);
                    } else {
                        //                 .. R1 is well conditioned, but non-square. AdjoINTEGER of R2
                        //                 is inverted to get the product of the Jacobi rotations
                        //                 used in Cgesvj. The Q-factor from the second QR
                        //                 factorization is then built in explicitly.
                        Ctrsm("L", "U", "C", "N", nr, nr, cone, cwork[(2 * n + 1) - 1], n, v, ldv);
                        if (nr < n) {
                            Claset("A", n - nr, nr, czero, czero, v[((nr + 1) - 1)], ldv);
                            Claset("A", nr, n - nr, czero, czero, v[((nr + 1) - 1) * ldv], ldv);
                            Claset("A", n - nr, n - nr, czero, cone, v[((nr + 1) - 1) + ((nr + 1) - 1) * ldv], ldv);
                        }
                        Cunmqr("L", "N", n, n, nr, cwork[(2 * n + 1) - 1], n, cwork[(n + 1) - 1], v, ldv, cwork[(2 * n + n * nr + nr + 1) - 1], lwork - 2 * n - n * nr - nr, ierr);
                    }
                    //
                } else if (condr2 < cond_ok) {
                    //
                    //              The matrix R2 is inverted. The solution of the matrix equation
                    //              is Q3^* * V3 = the product of the Jacobi rotations (appplied to
                    //              the lower triangular L3 from the LQ factorization of
                    //              R2=L3*Q3), pre-multiplied with the transposed Q3.
                    Cgesvj("L", "U", "N", nr, nr, v, ldv, sva, nr, u, ldu, cwork[(2 * n + n * nr + nr + 1) - 1], lwork - 2 * n - n * nr - nr, rwork, lrwork, info);
                    scalem = rwork[1 - 1];
                    numrank = nINTEGER[rwork[2 - 1] - 1];
                    for (p = 1; p <= nr; p = p + 1) {
                        Ccopy(nr, v[(p - 1) * ldv], 1, u[(p - 1) * ldu], 1);
                        CRscal(nr, sva[p - 1], u[(p - 1) * ldu], 1);
                    }
                    Ctrsm("L", "U", "N", "N", nr, nr, cone, cwork[(2 * n + 1) - 1], n, u, ldu);
                    //              .. apply the permutation from the second QR factorization
                    for (q = 1; q <= nr; q = q + 1) {
                        for (p = 1; p <= nr; p = p + 1) {
                            cwork[(2 * n + n * nr + nr + iwork[(n + p) - 1]) - 1] = u[(p - 1) + (q - 1) * ldu];
                        }
                        for (p = 1; p <= nr; p = p + 1) {
                            u[(p - 1) + (q - 1) * ldu] = cwork[(2 * n + n * nr + nr + p) - 1];
                        }
                    }
                    if (nr < n) {
                        Claset("A", n - nr, nr, czero, czero, v[((nr + 1) - 1)], ldv);
                        Claset("A", nr, n - nr, czero, czero, v[((nr + 1) - 1) * ldv], ldv);
                        Claset("A", n - nr, n - nr, czero, cone, v[((nr + 1) - 1) + ((nr + 1) - 1) * ldv], ldv);
                    }
                    Cunmqr("L", "N", n, n, nr, cwork[(2 * n + 1) - 1], n, cwork[(n + 1) - 1], v, ldv, cwork[(2 * n + n * nr + nr + 1) - 1], lwork - 2 * n - n * nr - nr, ierr);
                } else {
                    //              Last line of defense.
                    // #:(          This is a rather pathological case: no scaled condition
                    //              improvement after two pivoted QR factorizations. Other
                    //              possibility is that the rank revealing QR factorization
                    //              or the condition estimator has failed, or the COND_OK
                    //              is set very close to ONE (which is unnecessary). Normally,
                    //              this branch should never be executed, but in rare cases of
                    //              failure of the RRQR or condition estimator, the last line of
                    //              defense ensures that Cgejsv completes the task.
                    //              Compute the full SVD of L3 using Cgesvj with explicit
                    //              accumulation of Jacobi rotations.
                    Cgesvj("L", "U", "V", nr, nr, v, ldv, sva, nr, u, ldu, cwork[(2 * n + n * nr + nr + 1) - 1], lwork - 2 * n - n * nr - nr, rwork, lrwork, info);
                    scalem = rwork[1 - 1];
                    numrank = nINTEGER[rwork[2 - 1] - 1];
                    if (nr < n) {
                        Claset("A", n - nr, nr, czero, czero, v[((nr + 1) - 1)], ldv);
                        Claset("A", nr, n - nr, czero, czero, v[((nr + 1) - 1) * ldv], ldv);
                        Claset("A", n - nr, n - nr, czero, cone, v[((nr + 1) - 1) + ((nr + 1) - 1) * ldv], ldv);
                    }
                    Cunmqr("L", "N", n, n, nr, cwork[(2 * n + 1) - 1], n, cwork[(n + 1) - 1], v, ldv, cwork[(2 * n + n * nr + nr + 1) - 1], lwork - 2 * n - n * nr - nr, ierr);
                    //
                    Cunmlq("L", "C", nr, nr, nr, cwork[(2 * n + 1) - 1], n, cwork[(2 * n + n * nr + 1) - 1], u, ldu, cwork[(2 * n + n * nr + nr + 1) - 1], lwork - 2 * n - n * nr - nr, ierr);
                    for (q = 1; q <= nr; q = q + 1) {
                        for (p = 1; p <= nr; p = p + 1) {
                            cwork[(2 * n + n * nr + nr + iwork[(n + p) - 1]) - 1] = u[(p - 1) + (q - 1) * ldu];
                        }
                        for (p = 1; p <= nr; p = p + 1) {
                            u[(p - 1) + (q - 1) * ldu] = cwork[(2 * n + n * nr + nr + p) - 1];
                        }
                    }
                    //
                }
                //
                //           Permute the rows of V using the (column) permutation from the
                //           first QRF. Also, scale the columns to make them unit in
                //           Euclidean norm. This applies to all cases.
                //
                temp1 = sqrt(n.real()) * epsln;
                for (q = 1; q <= n; q = q + 1) {
                    for (p = 1; p <= n; p = p + 1) {
                        cwork[(2 * n + n * nr + nr + iwork[p - 1]) - 1] = v[(p - 1) + (q - 1) * ldv];
                    }
                    for (p = 1; p <= n; p = p + 1) {
                        v[(p - 1) + (q - 1) * ldv] = cwork[(2 * n + n * nr + nr + p) - 1];
                    }
                    xsc = one / RCnrm2[(n - 1) + (v[(q - 1) * ldv] - 1) * ldRCnrm2];
                    if ((xsc < (one - temp1)) || (xsc > (one + temp1))) {
                        CRscal(n, xsc, v[(q - 1) * ldv], 1);
                    }
                }
                //           At this moment, V contains the right singular vectors of A.
                //           Next, assemble the left singular vector matrix U (M x N).
                if (nr < m) {
                    Claset("A", m - nr, nr, czero, czero, u[((nr + 1) - 1)], ldu);
                    if (nr < n1) {
                        Claset("A", nr, n1 - nr, czero, czero, u[((nr + 1) - 1) * ldu], ldu);
                        Claset("A", m - nr, n1 - nr, czero, cone, u[((nr + 1) - 1) + ((nr + 1) - 1) * ldu], ldu);
                    }
                }
                //
                //           The Q matrix from the first QRF is built INTEGERo the left singular
                //           matrix U. This applies to all cases.
                //
                Cunmqr("L", "N", m, n1, n, a, lda, cwork, u, ldu, cwork[(n + 1) - 1], lwork - n, ierr);
                //
                //           The columns of U are normalized. The cost is O(M*N) flops.
                temp1 = sqrt(m.real()) * epsln;
                for (p = 1; p <= nr; p = p + 1) {
                    xsc = one / RCnrm2[(m - 1) + (u[(p - 1) * ldu] - 1) * ldRCnrm2];
                    if ((xsc < (one - temp1)) || (xsc > (one + temp1))) {
                        CRscal(m, xsc, u[(p - 1) * ldu], 1);
                    }
                }
                //
                //           If the initial QRF is computed with row pivoting, the left
                //           singular vectors must be adjusted.
                //
                if (rowpiv) {
                    Claswp(n1, u, ldu, 1, m - 1, iwork[(iwoff + 1) - 1], -1);
                }
                //
            } else {
                //
                //        .. the initial matrix A has almost orthogonal columns and
                //        the second QRF is not needed
                //
                Clacpy("U", n, n, a, lda, cwork[(n + 1) - 1], n);
                if (l2pert) {
                    xsc = sqrt(small);
                    for (p = 2; p <= n; p = p + 1) {
                        ctemp = xsc * cwork[(n + (p - 1) * n + p) - 1];
                        for (q = 1; q <= p - 1; q = q + 1) {
                            //                     CWORK(N+(q-1)*N+p)=-TEMP1 * ( CWORK(N+(p-1)*N+q) /
                            //     $                                        ABS(CWORK(N+(p-1)*N+q)) )
                            cwork[(n + (q - 1) * n + p) - 1] = -ctemp;
                        }
                    }
                } else {
                    Claset("L", n - 1, n - 1, czero, czero, cwork[(n + 2) - 1], n);
                }
                //
                Cgesvj("U", "U", "N", n, n, cwork[(n + 1) - 1], n, sva, n, u, ldu, cwork[(n + n * n + 1) - 1], lwork - n - n * n, rwork, lrwork, info);
                //
                scalem = rwork[1 - 1];
                numrank = nINTEGER[rwork[2 - 1] - 1];
                for (p = 1; p <= n; p = p + 1) {
                    Ccopy(n, cwork[(n + (p - 1) * n + 1) - 1], 1, u[(p - 1) * ldu], 1);
                    CRscal(n, sva[p - 1], cwork[(n + (p - 1) * n + 1) - 1], 1);
                }
                //
                Ctrsm("L", "U", "N", "N", n, n, cone, a, lda, cwork[(n + 1) - 1], n);
                for (p = 1; p <= n; p = p + 1) {
                    Ccopy(n, cwork[(n + p) - 1], n, v[(iwork[p - 1] - 1)], ldv);
                }
                temp1 = sqrt(n.real()) * epsln;
                for (p = 1; p <= n; p = p + 1) {
                    xsc = one / RCnrm2[(n - 1) + (v[(p - 1) * ldv] - 1) * ldRCnrm2];
                    if ((xsc < (one - temp1)) || (xsc > (one + temp1))) {
                        CRscal(n, xsc, v[(p - 1) * ldv], 1);
                    }
                }
                //
                //           Assemble the left singular vector matrix U (M x N).
                //
                if (n < m) {
                    Claset("A", m - n, n, czero, czero, u[((n + 1) - 1)], ldu);
                    if (n < n1) {
                        Claset("A", n, n1 - n, czero, czero, u[((n + 1) - 1) * ldu], ldu);
                        Claset("A", m - n, n1 - n, czero, cone, u[((n + 1) - 1) + ((n + 1) - 1) * ldu], ldu);
                    }
                }
                Cunmqr("L", "N", m, n1, n, a, lda, cwork, u, ldu, cwork[(n + 1) - 1], lwork - n, ierr);
                temp1 = sqrt(m.real()) * epsln;
                for (p = 1; p <= n1; p = p + 1) {
                    xsc = one / RCnrm2[(m - 1) + (u[(p - 1) * ldu] - 1) * ldRCnrm2];
                    if ((xsc < (one - temp1)) || (xsc > (one + temp1))) {
                        CRscal(m, xsc, u[(p - 1) * ldu], 1);
                    }
                }
                //
                if (rowpiv) {
                    Claswp(n1, u, ldu, 1, m - 1, iwork[(iwoff + 1) - 1], -1);
                }
                //
            }
            //
            //        end of the  >> almost orthogonal case <<  in the full SVD
            //
        } else {
            //
            //        This branch deploys a preconditioned Jacobi SVD with explicitly
            //        accumulated rotations. It is included as optional, mainly for
            //        experimental purposes. It does perform well, and can also be used.
            //        In this implementation, this branch will be automatically activated
            //        if the  condition number sigma_max(A) / sigma_min(A) is predicted
            //        to be greater than the overflow threshold. This is because the
            //        a posteriori computation of the singular vectors assumes robust
            //        implementation of BLAS and some LAPACK procedures, capable of working
            //        in presence of extreme values, e.g. when the singular values spread from
            //        the underflow to the overflow threshold.
            //
            for (p = 1; p <= nr; p = p + 1) {
                Ccopy(n - p + 1, a[(p - 1) + (p - 1) * lda], lda, v[(p - 1) + (p - 1) * ldv], 1);
                Clacgv(n - p + 1, v[(p - 1) + (p - 1) * ldv], 1);
            }
            //
            if (l2pert) {
                xsc = sqrt(small / epsln);
                for (q = 1; q <= nr; q = q + 1) {
                    ctemp = COMPLEX(xsc * abs(v[(q - 1) + (q - 1) * ldv]), zero);
                    for (p = 1; p <= n; p = p + 1) {
                        //     $                V(p,q) = TEMP1 * ( V(p,q) / ABS(V(p,q)) )
                        if ((p > q) && (abs(v[(p - 1) + (q - 1) * ldv]) <= temp1) || (p < q)) {
                            v[(p - 1) + (q - 1) * ldv] = ctemp;
                        }
                        if (p < q) {
                            v[(p - 1) + (q - 1) * ldv] = -v[(p - 1) + (q - 1) * ldv];
                        }
                    }
                }
            } else {
                Claset("U", nr - 1, nr - 1, czero, czero, v[(2 - 1) * ldv], ldv);
            }
            //
            Cgeqrf(n, nr, v, ldv, cwork[(n + 1) - 1], cwork[(2 * n + 1) - 1], lwork - 2 * n, ierr);
            Clacpy("L", n, nr, v, ldv, cwork[(2 * n + 1) - 1], n);
            //
            for (p = 1; p <= nr; p = p + 1) {
                Ccopy(nr - p + 1, v[(p - 1) + (p - 1) * ldv], ldv, u[(p - 1) + (p - 1) * ldu], 1);
                Clacgv(nr - p + 1, u[(p - 1) + (p - 1) * ldu], 1);
            }
            //
            if (l2pert) {
                xsc = sqrt(small / epsln);
                for (q = 2; q <= nr; q = q + 1) {
                    for (p = 1; p <= q - 1; p = p + 1) {
                        ctemp = COMPLEX(xsc * min(abs(u[(p - 1) + (p - 1) * ldu]), abs(u[(q - 1) + (q - 1) * ldu])), zero);
                        //                  U(p,q) = - TEMP1 * ( U(q,p) / ABS(U(q,p)) )
                        u[(p - 1) + (q - 1) * ldu] = -ctemp;
                    }
                }
            } else {
                Claset("U", nr - 1, nr - 1, czero, czero, u[(2 - 1) * ldu], ldu);
            }
            //
            Cgesvj("L", "U", "V", nr, nr, u, ldu, sva, n, v, ldv, cwork[(2 * n + n * nr + 1) - 1], lwork - 2 * n - n * nr, rwork, lrwork, info);
            scalem = rwork[1 - 1];
            numrank = nINTEGER[rwork[2 - 1] - 1];
            //
            if (nr < n) {
                Claset("A", n - nr, nr, czero, czero, v[((nr + 1) - 1)], ldv);
                Claset("A", nr, n - nr, czero, czero, v[((nr + 1) - 1) * ldv], ldv);
                Claset("A", n - nr, n - nr, czero, cone, v[((nr + 1) - 1) + ((nr + 1) - 1) * ldv], ldv);
            }
            //
            Cunmqr("L", "N", n, n, nr, cwork[(2 * n + 1) - 1], n, cwork[(n + 1) - 1], v, ldv, cwork[(2 * n + n * nr + nr + 1) - 1], lwork - 2 * n - n * nr - nr, ierr);
            //
            //           Permute the rows of V using the (column) permutation from the
            //           first QRF. Also, scale the columns to make them unit in
            //           Euclidean norm. This applies to all cases.
            //
            temp1 = sqrt(n.real()) * epsln;
            for (q = 1; q <= n; q = q + 1) {
                for (p = 1; p <= n; p = p + 1) {
                    cwork[(2 * n + n * nr + nr + iwork[p - 1]) - 1] = v[(p - 1) + (q - 1) * ldv];
                }
                for (p = 1; p <= n; p = p + 1) {
                    v[(p - 1) + (q - 1) * ldv] = cwork[(2 * n + n * nr + nr + p) - 1];
                }
                xsc = one / RCnrm2[(n - 1) + (v[(q - 1) * ldv] - 1) * ldRCnrm2];
                if ((xsc < (one - temp1)) || (xsc > (one + temp1))) {
                    CRscal(n, xsc, v[(q - 1) * ldv], 1);
                }
            }
            //
            //           At this moment, V contains the right singular vectors of A.
            //           Next, assemble the left singular vector matrix U (M x N).
            //
            if (nr < m) {
                Claset("A", m - nr, nr, czero, czero, u[((nr + 1) - 1)], ldu);
                if (nr < n1) {
                    Claset("A", nr, n1 - nr, czero, czero, u[((nr + 1) - 1) * ldu], ldu);
                    Claset("A", m - nr, n1 - nr, czero, cone, u[((nr + 1) - 1) + ((nr + 1) - 1) * ldu], ldu);
                }
            }
            //
            Cunmqr("L", "N", m, n1, n, a, lda, cwork, u, ldu, cwork[(n + 1) - 1], lwork - n, ierr);
            //
            if (rowpiv) {
                Claswp(n1, u, ldu, 1, m - 1, iwork[(iwoff + 1) - 1], -1);
            }
            //
        }
        if (transp) {
            //           .. swap U and V because the procedure worked on A^*
            for (p = 1; p <= n; p = p + 1) {
                Cswap(n, u[(p - 1) * ldu], 1, v[(p - 1) * ldv], 1);
            }
        }
        //
    }
    //     end of the full SVD
    //
    //     Undo scaling, if necessary (and possible)
    //
    if (uscal2 <= (big / sva[1 - 1]) * uscal1) {
        Rlascl("G", 0, 0, uscal1, uscal2, nr, 1, sva, n, ierr);
        uscal1 = one;
        uscal2 = one;
    }
    //
    if (nr < n) {
        for (p = nr + 1; p <= n; p = p + 1) {
            sva[p - 1] = zero;
        }
    }
    //
    rwork[1 - 1] = uscal2 * scalem;
    rwork[2 - 1] = uscal1;
    if (errest) {
        rwork[3 - 1] = sconda;
    }
    if (lsvec && rsvec) {
        rwork[4 - 1] = condr1;
        rwork[5 - 1] = condr2;
    }
    if (l2tran) {
        rwork[6 - 1] = entra;
        rwork[7 - 1] = entrat;
    }
    //
    iwork[1 - 1] = nr;
    iwork[2 - 1] = numrank;
    iwork[3 - 1] = warning;
    if (transp) {
        iwork[4 - 1] = 1;
    } else {
        iwork[4 - 1] = -1;
    }
    //
    //     ..
    //     .. END OF Cgejsv
    //     ..
}
