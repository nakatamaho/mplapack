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
    // .    through a rare Clahqr failure.  NL > NTINY = 15 is
    // .    required and NL <= NMIN = iMlaenv(ISPEC=12,...) is recom-
    // .    mended.  (The default value of NMIN is 75.)  Using NL = 49
    // .    allows up to six simultaneous shifts and a 16-by-16
    bool wantt = Mlsame(job, "S");
    bool initz = Mlsame(compz, "I");
    bool wantz = initz || Mlsame(compz, "V");
    const REAL rzero = 0.0;
    work[0] = COMPLEX(castREAL(max((INTEGER)1, n)), rzero);
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
    const COMPLEX zero = COMPLEX(0.0, 0.0);
    const COMPLEX one = COMPLEX(1.0, 0.0);
    INTEGER nmin = 0;
    const INTEGER ntiny = 15;
    INTEGER kbot = 0;
    const INTEGER nl = 49;
    COMPLEX hl[nl * nl];
    COMPLEX workl[nl];
    if (info != 0) {
        Mxerbla("Chseqr", -info);
        return;
    } else if (n == 0) {
        return;
    } else if (lquery) {
        Claqr0(wantt, wantz, n, ilo, ihi, h, ldh, w, ilo, ihi, z, ldz, work, lwork, info);
        work[0] = COMPLEX(max(work[0].real(), castREAL(max((INTEGER)1, n))), rzero);
        return;
    } else {
        if (ilo > 1) {
            Ccopy(ilo - 1, h, ldh + 1, w, 1);
        }
        if (ihi < n) {
            Ccopy(n - ihi, &h[((ihi + 1) - 1) + ((ihi + 1) - 1) * ldh], ldh + 1, &w[(ihi + 1) - 1], 1);
        }
        if (initz) {
            Claset("A", n, n, zero, one, z, ldz);
        }
        if (ilo == ihi) {
            w[ilo - 1] = h[(ilo - 1) + (ilo - 1) * ldh];
            return;
        }
        nmin = iMlaenv(12, "Chseqr", CHAR2(job, compz), n, ilo, ihi, lwork);
        nmin = max(ntiny, nmin);
        if (n > nmin) {
            Claqr0(wantt, wantz, n, ilo, ihi, h, ldh, w, ilo, ihi, z, ldz, work, lwork, info);
        } else {
            Clahqr(wantt, wantz, n, ilo, ihi, h, ldh, w, ilo, ihi, z, ldz, info);
            if (info > 0) {
                kbot = info;
                if (n >= nl) {
                    Claqr0(wantt, wantz, n, ilo, kbot, h, ldh, w, ilo, ihi, z, ldz, work, lwork, info);
                } else {
                    // .    scratch space to benefit from Claqr0.  Hence,
                    // .    tiny matrices must be copied into a larger
                    //
                    Clacpy("A", n, n, h, ldh, hl, nl);
                    hl[((n + 1) - 1) + (n - 1) * ldhl] = zero;
                    Claset("A", nl, nl - n, zero, zero, &hl[((n + 1) - 1) * ldhl], nl);
                    Claqr0(wantt, wantz, nl, ilo, kbot, hl, nl, w, ilo, ihi, z, ldz, workl, nl, info);
                    if (wantt || info != 0) {
                        Clacpy("A", n, n, hl, nl, h, ldh);
                    }
                }
            }
        }
        if ((wantt || info != 0) && n > 2) {
            Claset("L", n - 2, n - 2, zero, zero, &h[(3 - 1)], ldh);
        }
        work[0] = COMPLEX(max(castREAL(max((INTEGER)1, n)), work[0].real()), rzero);
    }
}
