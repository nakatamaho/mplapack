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

void Chseqr(const char *job, const char *compz, INTEGER const n, INTEGER const ilo, INTEGER const ihi, COMPLEX *h, INTEGER const ldh, COMPLEX *w, COMPLEX *z, INTEGER const ldz, COMPLEX *work, INTEGER const lwork, INTEGER &info) {
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
    //
    //     ==== Matrices of order NTINY or smaller must be processed by
    //     .    Clahqr because of insufficient subdiagonal scratch space.
    //     .    (This is a hard limit.) ====
    //
    //     ==== NL allocates some local workspace to help small matrices
    //     .    through a rare Clahqr failure.  NL > NTINY = 15 is
    //     .    required and NL <= NMIN = iMlaenv(ISPEC=12,...) is recom-
    //     .    mended.  (The default value of NMIN is 75.)  Using NL = 49
    //     .    allows up to six simultaneous shifts and a 16-by-16
    //     .    deflation window.  ====
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     ==== Decode and check the input parameters. ====
    //
    bool wantt = Mlsame(job, "S");
    bool initz = Mlsame(compz, "I");
    bool wantz = initz || Mlsame(compz, "V");
    const REAL rzero = 0.0;
    work[1 - 1] = COMPLEX((max((INTEGER)1, n)).real(), rzero);
    bool lquery = lwork == -1;
    //
    info = 0;
    if (!Mlsame(job, "E") && !wantt) {
        info = -1;
    } else if (!Mlsame(compz, "N") && !wantz) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (ilo < 1 || ilo > max((INTEGER)1, n)) {
        info = -4;
    } else if (ihi < min(ilo, n) || ihi > n) {
        info = -5;
    } else if (ldh < max((INTEGER)1, n)) {
        info = -7;
    } else if (ldz < 1 || (wantz && ldz < max((INTEGER)1, n))) {
        info = -10;
    } else if (lwork < max((INTEGER)1, n) && !lquery) {
        info = -12;
    }
    //
    const COMPLEX zero = (0.0, 0.0);
    const COMPLEX one = (1.0, 0.0);
    INTEGER nmin = 0;
    const INTEGER ntiny = 15;
    INTEGER kbot = 0;
    const INTEGER nl = 49;
    arr_1d<nl, COMPLEX> workl(fill0);
    if (info != 0) {
        //
        //        ==== Quick return in case of invalid argument. ====
        //
        Mxerbla("Chseqr", -info);
        return;
        //
    } else if (n == 0) {
        //
        //        ==== Quick return in case N = 0; nothing to do. ====
        //
        return;
        //
    } else if (lquery) {
        //
        //        ==== Quick return in case of a workspace query ====
        //
        Claqr0(wantt, wantz, n, ilo, ihi, h, ldh, w, ilo, ihi, z, ldz, work, lwork, info);
        //        ==== Ensure reported workspace size is backward-compatible with
        //        .    previous LAPACK versions. ====
        work[1 - 1] = COMPLEX(max(work[1 - 1].real(), (max((INTEGER)1, n)).real()), rzero);
        return;
        //
    } else {
        //
        //        ==== copy eigenvalues isolated by Cgebal ====
        //
        if (ilo > 1) {
            Ccopy(ilo - 1, h, ldh + 1, w, 1);
        }
        if (ihi < n) {
            Ccopy(n - ihi, &h[((ihi + 1) - 1) + ((ihi + 1) - 1) * ldh], ldh + 1, &w[(ihi + 1) - 1], 1);
        }
        //
        //        ==== Initialize Z, if requested ====
        //
        if (initz) {
            Claset("A", n, n, zero, one, z, ldz);
        }
        //
        //        ==== Quick return if possible ====
        //
        if (ilo == ihi) {
            w[ilo - 1] = h[(ilo - 1) + (ilo - 1) * ldh];
            return;
        }
        //
        //        ==== Clahqr/Claqr0 crossover poINTEGER ====
        //
        nmin = iMlaenv(12, "Chseqr", job[(1 - 1)] + compz[(1 - 1)], n, ilo, ihi, lwork);
        nmin = max(ntiny, nmin);
        //
        //        ==== Claqr0 for big matrices; Clahqr for small ones ====
        //
        if (n > nmin) {
            Claqr0(wantt, wantz, n, ilo, ihi, h, ldh, w, ilo, ihi, z, ldz, work, lwork, info);
        } else {
            //
            //           ==== Small matrix ====
            //
            Clahqr(wantt, wantz, n, ilo, ihi, h, ldh, w, ilo, ihi, z, ldz, info);
            //
            if (info > 0) {
                //
                //              ==== A rare Clahqr failure!  Claqr0 sometimes succeeds
                //              .    when Clahqr fails. ====
                //
                kbot = info;
                //
                if (n >= nl) {
                    //
                    //                 ==== Larger matrices have enough subdiagonal scratch
                    //                 .    space to call Claqr0 directly. ====
                    //
                    Claqr0(wantt, wantz, n, ilo, kbot, h, ldh, w, ilo, ihi, z, ldz, work, lwork, info);
                    //
                } else {
                    //
                    //                 ==== Tiny matrices don't have enough subdiagonal
                    //                 .    scratch space to benefit from Claqr0.  Hence,
                    //                 .    tiny matrices must be copied into a larger
                    //                 .    array before calling Claqr0. ====
                    //
                    Clacpy("A", n, n, h, ldh, hl, nl);
                    hl[((n + 1) - 1) + (n - 1) * ldhl] = zero;
                    Claset("A", nl, nl - n, zero, zero, hl[((n + 1) - 1) * ldhl], nl);
                    Claqr0(wantt, wantz, nl, ilo, kbot, hl, nl, w, ilo, ihi, z, ldz, workl, nl, info);
                    if (wantt || info != 0) {
                        Clacpy("A", n, n, hl, nl, h, ldh);
                    }
                }
            }
        }
        //
        //        ==== Clear out the trash, if necessary. ====
        //
        if ((wantt || info != 0) && n > 2) {
            Claset("L", n - 2, n - 2, zero, zero, &h[(3 - 1)], ldh);
        }
        //
        //        ==== Ensure reported workspace size is backward-compatible with
        //        .    previous LAPACK versions. ====
        //
        work[1 - 1] = COMPLEX(max((max((INTEGER)1, n)).real(), &work[1 - 1].real()), rzero);
    }
    //
    //     ==== End of Chseqr ====
    //
}
