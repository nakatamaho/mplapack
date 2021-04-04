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

void Chbgst(const char *vect, const char *uplo, INTEGER const &n, INTEGER const &ka, INTEGER const &kb, COMPLEX *ab, INTEGER const &ldab, COMPLEX *bb, INTEGER const &ldbb, COMPLEX *x, INTEGER const &ldx, COMPLEX *work, REAL *rwork, INTEGER &info) {
    bool wantx = false;
    bool upper = false;
    INTEGER ka1 = 0;
    INTEGER kb1 = 0;
    INTEGER inca = 0;
    const COMPLEX czero = (0.0, 0.0);
    const COMPLEX cone = (1.0, 0.0);
    INTEGER m = 0;
    bool update = false;
    INTEGER i = 0;
    INTEGER kbt = 0;
    INTEGER i0 = 0;
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    REAL bii = 0.0;
    INTEGER j = 0;
    INTEGER k = 0;
    const REAL one = 1.0;
    COMPLEX ra1 = 0.0;
    COMPLEX ra = 0.0;
    COMPLEX t = 0.0;
    INTEGER j2 = 0;
    INTEGER nr = 0;
    INTEGER j1 = 0;
    INTEGER j2t = 0;
    INTEGER nrt = 0;
    INTEGER l = 0;
    INTEGER nx = 0;
    INTEGER j1t = 0;
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
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input parameters
    //
    wantx = Mlsame(vect, "V");
    upper = Mlsame(uplo, "U");
    ka1 = ka + 1;
    kb1 = kb + 1;
    info = 0;
    if (!wantx && !Mlsame(vect, "N")) {
        info = -1;
    } else if (!upper && !Mlsame(uplo, "L")) {
        info = -2;
    } else if (n < 0) {
        info = -3;
    } else if (ka < 0) {
        info = -4;
    } else if (kb < 0 || kb > ka) {
        info = -5;
    } else if (ldab < ka + 1) {
        info = -7;
    } else if (ldbb < kb + 1) {
        info = -9;
    } else if (ldx < 1 || wantx && ldx < max((INTEGER)1, n)) {
        info = -11;
    }
    if (info != 0) {
        Mxerbla("Chbgst", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n == 0) {
        return;
    }
    //
    inca = ldab * ka1;
    //
    //     Initialize X to the unit matrix, if needed
    //
    if (wantx) {
        Claset("Full", n, n, czero, cone, x, ldx);
    }
    //
    //     Set M to the splitting poINTEGER m. It must be the same value as is
    //     used in Cpbstf. The chosen value allows the arrays WORK and RWORK
    //
    m = (n + kb) / 2;
    //
    //     The routine works in two phases, corresponding to the two halves
    //     of the split Cholesky factorization of B as S**H*S where
    //
    //     S = ( U    )
    //         ( M  L )
    //
    //     with U upper triangular of order m, and L lower triangular of
    //     order n-m. S has the same bandwidth as B.
    //
    //     S is treated as a product of elementary matrices:
    //
    //     S = S(m)*S(m-1)*...*S(2)*S(1)*S(m+1)*S(m+2)*...*S(n-1)*S(n)
    //
    //     where S(i) is determined by the i-th row of S.
    //
    //     In phase 1, the index i takes the values n, n-1, ... , m+1;
    //     in phase 2, it takes the values 1, 2, ... , m.
    //
    //     For each value of i, the current matrix A is updated by forming
    //     inv(S(i))**H*A*inv(S(i)). This creates a triangular bulge outside
    //     the band of A. The bulge is then pushed down toward the bottom of
    //     A in phase 1, and up toward the top of A in phase 2, by applying
    //     plane rotations.
    //
    //     There are kb*(kb+1)/2 elements in the bulge, but at most 2*kb-1
    //     of them are linearly independent, so annihilating a bulge requires
    //     only 2*kb-1 plane rotations. The rotations are divided INTEGERo a 1st
    //     set of kb-1 rotations, and a 2nd set of kb rotations.
    //
    //     Wherever possible, rotations are generated and applied in vector
    //     operations of length NR between the indices J1 and J2 (sometimes
    //     replaced by modified values NRT, J1T or J2T).
    //
    //     The real cosines and complex sines of the rotations are stored in
    //     the arrays RWORK and WORK, those of the 1st set in elements
    //     2:m-kb-1, and those of the 2nd set in elements m-kb+1:n.
    //
    //     The bulges are not formed explicitly; nonzero elements outside the
    //     band are created only when they are required for generating new
    //     rotations; they are stored in the array WORK, in positions where
    //     they are later overwritten by the sines of the rotations which
    //     annihilate them.
    //
    //     **************************** Phase 1 *****************************
    //
    //     The logical structure of this phase is:
    //
    //     UPDATE = .TRUE.
    //     DO I = N, M + 1, -1
    //        use S(i) to update A and create a new bulge
    //        apply rotations to push all bulges KA positions downward
    //     END DO
    //     UPDATE = .FALSE.
    //     DO I = M + KA + 1, N - 1
    //        apply rotations to push all bulges KA positions downward
    //     END DO
    //
    //     To avoid duplicating code, the two loops are merged.
    //
    update = true;
    i = n + 1;
statement_10:
    if (update) {
        i = i - 1;
        kbt = min(kb, i - 1);
        i0 = i - 1;
        i1 = min(n, i + ka);
        i2 = i - kbt + ka1;
        if (i < m + 1) {
            update = false;
            i++;
            i0 = m;
            if (ka == 0) {
                goto statement_480;
            }
            goto statement_10;
        }
    } else {
        i += ka;
        if (i > n - 1) {
            goto statement_480;
        }
    }
    //
    if (upper) {
        //
        //        Transform A, working with the upper triangle
        //
        if (update) {
            //
            //           Form  inv(S(i))**H * A * inv(S(i))
            //
            bii = bb[(kb1 - 1) + (i - 1) * ldbb].real();
            ab[(ka1 - 1) + (i - 1) * ldab] = (ab[(ka1 - 1) + (i - 1) * ldab].real() / bii) / bii;
            for (j = i + 1; j <= i1; j = j + 1) {
                ab[((i - j + ka1) - 1) + (j - 1) * ldab] = ab[((i - j + ka1) - 1) + (j - 1) * ldab] / bii;
            }
            for (j = max((INTEGER)1, i - ka); j <= i - 1; j = j + 1) {
                ab[((j - i + ka1) - 1) + (i - 1) * ldab] = ab[((j - i + ka1) - 1) + (i - 1) * ldab] / bii;
            }
            for (k = i - kbt; k <= i - 1; k = k + 1) {
                for (j = i - kbt; j <= k; j = j + 1) {
                    ab[((j - k + ka1) - 1) + (k - 1) * ldab] = ab[((j - k + ka1) - 1) + (k - 1) * ldab] - bb[((j - i + kb1) - 1) + (i - 1) * ldbb] * conj(ab[((k - i + ka1) - 1) + (i - 1) * ldab]) - conj(bb[((k - i + kb1) - 1) + (i - 1) * ldbb]) * ab[((j - i + ka1) - 1) + (i - 1) * ldab] + ab[(ka1 - 1) + (i - 1) * ldab].real() * bb[((j - i + kb1) - 1) + (i - 1) * ldbb] * conj(bb[((k - i + kb1) - 1) + (i - 1) * ldbb]);
                }
                for (j = max((INTEGER)1, i - ka); j <= i - kbt - 1; j = j + 1) {
                    ab[((j - k + ka1) - 1) + (k - 1) * ldab] = ab[((j - k + ka1) - 1) + (k - 1) * ldab] - conj(bb[((k - i + kb1) - 1) + (i - 1) * ldbb]) * ab[((j - i + ka1) - 1) + (i - 1) * ldab];
                }
            }
            for (j = i; j <= i1; j = j + 1) {
                for (k = max(j - ka, i - kbt); k <= i - 1; k = k + 1) {
                    ab[((k - j + ka1) - 1) + (j - 1) * ldab] = ab[((k - j + ka1) - 1) + (j - 1) * ldab] - bb[((k - i + kb1) - 1) + (i - 1) * ldbb] * ab[((i - j + ka1) - 1) + (j - 1) * ldab];
                }
            }
            //
            if (wantx) {
                //
                //              post-multiply X by inv(S(i))
                //
                CRscal(n - m, one / bii, x[((m + 1) - 1) + (i - 1) * ldx], 1);
                if (kbt > 0) {
                    Cgerc(n - m, kbt, -cone, x[((m + 1) - 1) + (i - 1) * ldx], 1, bb[((kb1 - kbt) - 1) + (i - 1) * ldbb], 1, x[((m + 1) - 1) + ((i - kbt) - 1) * ldx], ldx);
                }
            }
            //
            //           store a(i,i1) in RA1 for use in next loop over K
            //
            ra1 = ab[((i - i1 + ka1) - 1) + (i1 - 1) * ldab];
        }
        //
        //        Generate and apply vectors of rotations to chase all the
        //        existing bulges KA positions down toward the bottom of the
        //        band
        //
        for (k = 1; k <= kb - 1; k = k + 1) {
            if (update) {
                //
                //              Determine the rotations which would annihilate the bulge
                //              which has in theory just been created
                //
                if (i - k + ka < n && i - k > 1) {
                    //
                    //                 generate rotation to annihilate a(i,i-k+ka+1)
                    //
                    Clartg(ab[((k + 1) - 1) + ((i - k + ka) - 1) * ldab], ra1, rwork[(i - k + ka - m) - 1], work[(i - k + ka - m) - 1], ra);
                    //
                    //                 create nonzero element a(i-k,i-k+ka+1) outside the
                    //                 band and store it in WORK(i-k)
                    //
                    t = -bb[((kb1 - k) - 1) + (i - 1) * ldbb] * ra1;
                    work[(i - k) - 1] = rwork[(i - k + ka - m) - 1] * t - conj(work[(i - k + ka - m) - 1]) * ab[((i - k + ka) - 1) * ldab];
                    ab[((i - k + ka) - 1) * ldab] = work[(i - k + ka - m) - 1] * t + rwork[(i - k + ka - m) - 1] * ab[((i - k + ka) - 1) * ldab];
                    ra1 = ra;
                }
            }
            j2 = i - k - 1 + max((INTEGER)1, k - i0 + 2) * ka1;
            nr = (n - j2 + ka) / ka1;
            j1 = j2 + (nr - 1) * ka1;
            if (update) {
                j2t = max(j2, i + 2 * ka - k + 1);
            } else {
                j2t = j2;
            }
            nrt = (n - j2t + ka) / ka1;
            for (j = j2t; j <= j1; j = j + ka1) {
                //
                //              create nonzero element a(j-ka,j+1) outside the band
                //              and store it in WORK(j-m)
                //
                work[(j - m) - 1] = work[(j - m) - 1] * ab[((j + 1) - 1) * ldab];
                ab[((j + 1) - 1) * ldab] = rwork[(j - m) - 1] * ab[((j + 1) - 1) * ldab];
            }
            //
            //           generate rotations in 1st set to annihilate elements which
            //           have been created outside the band
            //
            if (nrt > 0) {
                Clargv(nrt, ab[(j2t - 1) * ldab], inca, work[(j2t - m) - 1], ka1, rwork[(j2t - m) - 1], ka1);
            }
            if (nr > 0) {
                //
                //              apply rotations in 1st set from the right
                //
                for (l = 1; l <= ka - 1; l = l + 1) {
                    Clartv(nr, ab[((ka1 - l) - 1) + (j2 - 1) * ldab], inca, ab[((ka - l) - 1) + ((j2 + 1) - 1) * ldab], inca, rwork[(j2 - m) - 1], work[(j2 - m) - 1], ka1);
                }
                //
                //              apply rotations in 1st set from both sides to diagonal
                //              blocks
                //
                Clar2v(nr, ab[(ka1 - 1) + (j2 - 1) * ldab], ab[(ka1 - 1) + ((j2 + 1) - 1) * ldab], ab[(ka - 1) + ((j2 + 1) - 1) * ldab], inca, rwork[(j2 - m) - 1], work[(j2 - m) - 1], ka1);
                //
                Clacgv(nr, work[(j2 - m) - 1], ka1);
            }
            //
            //           start applying rotations in 1st set from the left
            //
            for (l = ka - 1; l >= kb - k + 1; l = l - 1) {
                nrt = (n - j2 + l) / ka1;
                if (nrt > 0) {
                    Clartv(nrt, ab[(l - 1) + ((j2 + ka1 - l) - 1) * ldab], inca, ab[((l + 1) - 1) + ((j2 + ka1 - l) - 1) * ldab], inca, rwork[(j2 - m) - 1], work[(j2 - m) - 1], ka1);
                }
            }
            //
            if (wantx) {
                //
                //              post-multiply X by product of rotations in 1st set
                //
                for (j = j2; j <= j1; j = j + ka1) {
                    Crot(n - m, x[((m + 1) - 1) + (j - 1) * ldx], 1, x[((m + 1) - 1) + ((j + 1) - 1) * ldx], 1, rwork[(j - m) - 1], conj(work[(j - m) - 1]));
                }
            }
        }
        //
        if (update) {
            if (i2 <= n && kbt > 0) {
                //
                //              create nonzero element a(i-kbt,i-kbt+ka+1) outside the
                //              band and store it in WORK(i-kbt)
                //
                work[(i - kbt) - 1] = -bb[((kb1 - kbt) - 1) + (i - 1) * ldbb] * ra1;
            }
        }
        //
        for (k = kb; k >= 1; k = k - 1) {
            if (update) {
                j2 = i - k - 1 + max(2, k - i0 + 1) * ka1;
            } else {
                j2 = i - k - 1 + max((INTEGER)1, k - i0 + 1) * ka1;
            }
            //
            //           finish applying rotations in 2nd set from the left
            //
            for (l = kb - k; l >= 1; l = l - 1) {
                nrt = (n - j2 + ka + l) / ka1;
                if (nrt > 0) {
                    Clartv(nrt, ab[(l - 1) + ((j2 - l + 1) - 1) * ldab], inca, ab[((l + 1) - 1) + ((j2 - l + 1) - 1) * ldab], inca, rwork[(j2 - ka) - 1], work[(j2 - ka) - 1], ka1);
                }
            }
            nr = (n - j2 + ka) / ka1;
            j1 = j2 + (nr - 1) * ka1;
            for (j = j1; j <= j2; j = j + -ka1) {
                work[j - 1] = work[(j - ka) - 1];
                rwork[j - 1] = rwork[(j - ka) - 1];
            }
            for (j = j2; j <= j1; j = j + ka1) {
                //
                //              create nonzero element a(j-ka,j+1) outside the band
                //              and store it in WORK(j)
                //
                work[j - 1] = work[j - 1] * ab[((j + 1) - 1) * ldab];
                ab[((j + 1) - 1) * ldab] = rwork[j - 1] * ab[((j + 1) - 1) * ldab];
            }
            if (update) {
                if (i - k < n - ka && k <= kbt) {
                    work[(i - k + ka) - 1] = work[(i - k) - 1];
                }
            }
        }
        //
        for (k = kb; k >= 1; k = k - 1) {
            j2 = i - k - 1 + max((INTEGER)1, k - i0 + 1) * ka1;
            nr = (n - j2 + ka) / ka1;
            j1 = j2 + (nr - 1) * ka1;
            if (nr > 0) {
                //
                //              generate rotations in 2nd set to annihilate elements
                //              which have been created outside the band
                //
                Clargv(nr, ab[(j2 - 1) * ldab], inca, work[j2 - 1], ka1, rwork[j2 - 1], ka1);
                //
                //              apply rotations in 2nd set from the right
                //
                for (l = 1; l <= ka - 1; l = l + 1) {
                    Clartv(nr, ab[((ka1 - l) - 1) + (j2 - 1) * ldab], inca, ab[((ka - l) - 1) + ((j2 + 1) - 1) * ldab], inca, rwork[j2 - 1], work[j2 - 1], ka1);
                }
                //
                //              apply rotations in 2nd set from both sides to diagonal
                //              blocks
                //
                Clar2v(nr, ab[(ka1 - 1) + (j2 - 1) * ldab], ab[(ka1 - 1) + ((j2 + 1) - 1) * ldab], ab[(ka - 1) + ((j2 + 1) - 1) * ldab], inca, rwork[j2 - 1], work[j2 - 1], ka1);
                //
                Clacgv(nr, work[j2 - 1], ka1);
            }
            //
            //           start applying rotations in 2nd set from the left
            //
            for (l = ka - 1; l >= kb - k + 1; l = l - 1) {
                nrt = (n - j2 + l) / ka1;
                if (nrt > 0) {
                    Clartv(nrt, ab[(l - 1) + ((j2 + ka1 - l) - 1) * ldab], inca, ab[((l + 1) - 1) + ((j2 + ka1 - l) - 1) * ldab], inca, rwork[j2 - 1], work[j2 - 1], ka1);
                }
            }
            //
            if (wantx) {
                //
                //              post-multiply X by product of rotations in 2nd set
                //
                for (j = j2; j <= j1; j = j + ka1) {
                    Crot(n - m, x[((m + 1) - 1) + (j - 1) * ldx], 1, x[((m + 1) - 1) + ((j + 1) - 1) * ldx], 1, rwork[j - 1], conj(work[j - 1]));
                }
            }
        }
        //
        for (k = 1; k <= kb - 1; k = k + 1) {
            j2 = i - k - 1 + max((INTEGER)1, k - i0 + 2) * ka1;
            //
            //           finish applying rotations in 1st set from the left
            //
            for (l = kb - k; l >= 1; l = l - 1) {
                nrt = (n - j2 + l) / ka1;
                if (nrt > 0) {
                    Clartv(nrt, ab[(l - 1) + ((j2 + ka1 - l) - 1) * ldab], inca, ab[((l + 1) - 1) + ((j2 + ka1 - l) - 1) * ldab], inca, rwork[(j2 - m) - 1], work[(j2 - m) - 1], ka1);
                }
            }
        }
        //
        if (kb > 1) {
            for (j = n - 1; j >= j2 + ka; j = j - 1) {
                rwork[(j - m) - 1] = rwork[(j - ka - m) - 1];
                work[(j - m) - 1] = work[(j - ka - m) - 1];
            }
        }
        //
    } else {
        //
        //        Transform A, working with the lower triangle
        //
        if (update) {
            //
            //           Form  inv(S(i))**H * A * inv(S(i))
            //
            bii = bb[(i - 1) * ldbb].real();
            ab[(i - 1) * ldab] = (ab[(i - 1) * ldab].real() / bii) / bii;
            for (j = i + 1; j <= i1; j = j + 1) {
                ab[((j - i + 1) - 1) + (i - 1) * ldab] = ab[((j - i + 1) - 1) + (i - 1) * ldab] / bii;
            }
            for (j = max((INTEGER)1, i - ka); j <= i - 1; j = j + 1) {
                ab[((i - j + 1) - 1) + (j - 1) * ldab] = ab[((i - j + 1) - 1) + (j - 1) * ldab] / bii;
            }
            for (k = i - kbt; k <= i - 1; k = k + 1) {
                for (j = i - kbt; j <= k; j = j + 1) {
                    ab[((k - j + 1) - 1) + (j - 1) * ldab] = ab[((k - j + 1) - 1) + (j - 1) * ldab] - bb[((i - j + 1) - 1) + (j - 1) * ldbb] * conj(ab[((i - k + 1) - 1) + (k - 1) * ldab]) - conj(bb[((i - k + 1) - 1) + (k - 1) * ldbb]) * ab[((i - j + 1) - 1) + (j - 1) * ldab] + ab[(i - 1) * ldab].real() * bb[((i - j + 1) - 1) + (j - 1) * ldbb] * conj(bb[((i - k + 1) - 1) + (k - 1) * ldbb]);
                }
                for (j = max((INTEGER)1, i - ka); j <= i - kbt - 1; j = j + 1) {
                    ab[((k - j + 1) - 1) + (j - 1) * ldab] = ab[((k - j + 1) - 1) + (j - 1) * ldab] - conj(bb[((i - k + 1) - 1) + (k - 1) * ldbb]) * ab[((i - j + 1) - 1) + (j - 1) * ldab];
                }
            }
            for (j = i; j <= i1; j = j + 1) {
                for (k = max(j - ka, i - kbt); k <= i - 1; k = k + 1) {
                    ab[((j - k + 1) - 1) + (k - 1) * ldab] = ab[((j - k + 1) - 1) + (k - 1) * ldab] - bb[((i - k + 1) - 1) + (k - 1) * ldbb] * ab[((j - i + 1) - 1) + (i - 1) * ldab];
                }
            }
            //
            if (wantx) {
                //
                //              post-multiply X by inv(S(i))
                //
                CRscal(n - m, one / bii, x[((m + 1) - 1) + (i - 1) * ldx], 1);
                if (kbt > 0) {
                    Cgeru(n - m, kbt, -cone, x[((m + 1) - 1) + (i - 1) * ldx], 1, bb[((kbt + 1) - 1) + ((i - kbt) - 1) * ldbb], ldbb - 1, x[((m + 1) - 1) + ((i - kbt) - 1) * ldx], ldx);
                }
            }
            //
            //           store a(i1,i) in RA1 for use in next loop over K
            //
            ra1 = ab[((i1 - i + 1) - 1) + (i - 1) * ldab];
        }
        //
        //        Generate and apply vectors of rotations to chase all the
        //        existing bulges KA positions down toward the bottom of the
        //        band
        //
        for (k = 1; k <= kb - 1; k = k + 1) {
            if (update) {
                //
                //              Determine the rotations which would annihilate the bulge
                //              which has in theory just been created
                //
                if (i - k + ka < n && i - k > 1) {
                    //
                    //                 generate rotation to annihilate a(i-k+ka+1,i)
                    //
                    Clartg(ab[((ka1 - k) - 1) + (i - 1) * ldab], ra1, rwork[(i - k + ka - m) - 1], work[(i - k + ka - m) - 1], ra);
                    //
                    //                 create nonzero element a(i-k+ka+1,i-k) outside the
                    //                 band and store it in WORK(i-k)
                    //
                    t = -bb[((k + 1) - 1) + ((i - k) - 1) * ldbb] * ra1;
                    work[(i - k) - 1] = rwork[(i - k + ka - m) - 1] * t - conj(work[(i - k + ka - m) - 1]) * ab[(ka1 - 1) + ((i - k) - 1) * ldab];
                    ab[(ka1 - 1) + ((i - k) - 1) * ldab] = work[(i - k + ka - m) - 1] * t + rwork[(i - k + ka - m) - 1] * ab[(ka1 - 1) + ((i - k) - 1) * ldab];
                    ra1 = ra;
                }
            }
            j2 = i - k - 1 + max((INTEGER)1, k - i0 + 2) * ka1;
            nr = (n - j2 + ka) / ka1;
            j1 = j2 + (nr - 1) * ka1;
            if (update) {
                j2t = max(j2, i + 2 * ka - k + 1);
            } else {
                j2t = j2;
            }
            nrt = (n - j2t + ka) / ka1;
            for (j = j2t; j <= j1; j = j + ka1) {
                //
                //              create nonzero element a(j+1,j-ka) outside the band
                //              and store it in WORK(j-m)
                //
                work[(j - m) - 1] = work[(j - m) - 1] * ab[(ka1 - 1) + ((j - ka + 1) - 1) * ldab];
                ab[(ka1 - 1) + ((j - ka + 1) - 1) * ldab] = rwork[(j - m) - 1] * ab[(ka1 - 1) + ((j - ka + 1) - 1) * ldab];
            }
            //
            //           generate rotations in 1st set to annihilate elements which
            //           have been created outside the band
            //
            if (nrt > 0) {
                Clargv(nrt, ab[(ka1 - 1) + ((j2t - ka) - 1) * ldab], inca, work[(j2t - m) - 1], ka1, rwork[(j2t - m) - 1], ka1);
            }
            if (nr > 0) {
                //
                //              apply rotations in 1st set from the left
                //
                for (l = 1; l <= ka - 1; l = l + 1) {
                    Clartv(nr, ab[((l + 1) - 1) + ((j2 - l) - 1) * ldab], inca, ab[((l + 2) - 1) + ((j2 - l) - 1) * ldab], inca, rwork[(j2 - m) - 1], work[(j2 - m) - 1], ka1);
                }
                //
                //              apply rotations in 1st set from both sides to diagonal
                //              blocks
                //
                Clar2v(nr, ab[(j2 - 1) * ldab], ab[((j2 + 1) - 1) * ldab], ab[(2 - 1) + (j2 - 1) * ldab], inca, rwork[(j2 - m) - 1], work[(j2 - m) - 1], ka1);
                //
                Clacgv(nr, work[(j2 - m) - 1], ka1);
            }
            //
            //           start applying rotations in 1st set from the right
            //
            for (l = ka - 1; l >= kb - k + 1; l = l - 1) {
                nrt = (n - j2 + l) / ka1;
                if (nrt > 0) {
                    Clartv(nrt, ab[((ka1 - l + 1) - 1) + (j2 - 1) * ldab], inca, ab[((ka1 - l) - 1) + ((j2 + 1) - 1) * ldab], inca, rwork[(j2 - m) - 1], work[(j2 - m) - 1], ka1);
                }
            }
            //
            if (wantx) {
                //
                //              post-multiply X by product of rotations in 1st set
                //
                for (j = j2; j <= j1; j = j + ka1) {
                    Crot(n - m, x[((m + 1) - 1) + (j - 1) * ldx], 1, x[((m + 1) - 1) + ((j + 1) - 1) * ldx], 1, rwork[(j - m) - 1], work[(j - m) - 1]);
                }
            }
        }
        //
        if (update) {
            if (i2 <= n && kbt > 0) {
                //
                //              create nonzero element a(i-kbt+ka+1,i-kbt) outside the
                //              band and store it in WORK(i-kbt)
                //
                work[(i - kbt) - 1] = -bb[((kbt + 1) - 1) + ((i - kbt) - 1) * ldbb] * ra1;
            }
        }
        //
        for (k = kb; k >= 1; k = k - 1) {
            if (update) {
                j2 = i - k - 1 + max(2, k - i0 + 1) * ka1;
            } else {
                j2 = i - k - 1 + max((INTEGER)1, k - i0 + 1) * ka1;
            }
            //
            //           finish applying rotations in 2nd set from the right
            //
            for (l = kb - k; l >= 1; l = l - 1) {
                nrt = (n - j2 + ka + l) / ka1;
                if (nrt > 0) {
                    Clartv(nrt, ab[((ka1 - l + 1) - 1) + ((j2 - ka) - 1) * ldab], inca, ab[((ka1 - l) - 1) + ((j2 - ka + 1) - 1) * ldab], inca, rwork[(j2 - ka) - 1], work[(j2 - ka) - 1], ka1);
                }
            }
            nr = (n - j2 + ka) / ka1;
            j1 = j2 + (nr - 1) * ka1;
            for (j = j1; j <= j2; j = j + -ka1) {
                work[j - 1] = work[(j - ka) - 1];
                rwork[j - 1] = rwork[(j - ka) - 1];
            }
            for (j = j2; j <= j1; j = j + ka1) {
                //
                //              create nonzero element a(j+1,j-ka) outside the band
                //              and store it in WORK(j)
                //
                work[j - 1] = work[j - 1] * ab[(ka1 - 1) + ((j - ka + 1) - 1) * ldab];
                ab[(ka1 - 1) + ((j - ka + 1) - 1) * ldab] = rwork[j - 1] * ab[(ka1 - 1) + ((j - ka + 1) - 1) * ldab];
            }
            if (update) {
                if (i - k < n - ka && k <= kbt) {
                    work[(i - k + ka) - 1] = work[(i - k) - 1];
                }
            }
        }
        //
        for (k = kb; k >= 1; k = k - 1) {
            j2 = i - k - 1 + max((INTEGER)1, k - i0 + 1) * ka1;
            nr = (n - j2 + ka) / ka1;
            j1 = j2 + (nr - 1) * ka1;
            if (nr > 0) {
                //
                //              generate rotations in 2nd set to annihilate elements
                //              which have been created outside the band
                //
                Clargv(nr, ab[(ka1 - 1) + ((j2 - ka) - 1) * ldab], inca, work[j2 - 1], ka1, rwork[j2 - 1], ka1);
                //
                //              apply rotations in 2nd set from the left
                //
                for (l = 1; l <= ka - 1; l = l + 1) {
                    Clartv(nr, ab[((l + 1) - 1) + ((j2 - l) - 1) * ldab], inca, ab[((l + 2) - 1) + ((j2 - l) - 1) * ldab], inca, rwork[j2 - 1], work[j2 - 1], ka1);
                }
                //
                //              apply rotations in 2nd set from both sides to diagonal
                //              blocks
                //
                Clar2v(nr, ab[(j2 - 1) * ldab], ab[((j2 + 1) - 1) * ldab], ab[(2 - 1) + (j2 - 1) * ldab], inca, rwork[j2 - 1], work[j2 - 1], ka1);
                //
                Clacgv(nr, work[j2 - 1], ka1);
            }
            //
            //           start applying rotations in 2nd set from the right
            //
            for (l = ka - 1; l >= kb - k + 1; l = l - 1) {
                nrt = (n - j2 + l) / ka1;
                if (nrt > 0) {
                    Clartv(nrt, ab[((ka1 - l + 1) - 1) + (j2 - 1) * ldab], inca, ab[((ka1 - l) - 1) + ((j2 + 1) - 1) * ldab], inca, rwork[j2 - 1], work[j2 - 1], ka1);
                }
            }
            //
            if (wantx) {
                //
                //              post-multiply X by product of rotations in 2nd set
                //
                for (j = j2; j <= j1; j = j + ka1) {
                    Crot(n - m, x[((m + 1) - 1) + (j - 1) * ldx], 1, x[((m + 1) - 1) + ((j + 1) - 1) * ldx], 1, rwork[j - 1], work[j - 1]);
                }
            }
        }
        //
        for (k = 1; k <= kb - 1; k = k + 1) {
            j2 = i - k - 1 + max((INTEGER)1, k - i0 + 2) * ka1;
            //
            //           finish applying rotations in 1st set from the right
            //
            for (l = kb - k; l >= 1; l = l - 1) {
                nrt = (n - j2 + l) / ka1;
                if (nrt > 0) {
                    Clartv(nrt, ab[((ka1 - l + 1) - 1) + (j2 - 1) * ldab], inca, ab[((ka1 - l) - 1) + ((j2 + 1) - 1) * ldab], inca, rwork[(j2 - m) - 1], work[(j2 - m) - 1], ka1);
                }
            }
        }
        //
        if (kb > 1) {
            for (j = n - 1; j >= j2 + ka; j = j - 1) {
                rwork[(j - m) - 1] = rwork[(j - ka - m) - 1];
                work[(j - m) - 1] = work[(j - ka - m) - 1];
            }
        }
        //
    }
    //
    goto statement_10;
//
statement_480:
    //
    //     **************************** Phase 2 *****************************
    //
    //     The logical structure of this phase is:
    //
    //     UPDATE = .TRUE.
    //     DO I = 1, M
    //        use S(i) to update A and create a new bulge
    //        apply rotations to push all bulges KA positions upward
    //     END DO
    //     UPDATE = .FALSE.
    //     DO I = M - KA - 1, 2, -1
    //        apply rotations to push all bulges KA positions upward
    //     END DO
    //
    //     To avoid duplicating code, the two loops are merged.
    //
    update = true;
    i = 0;
statement_490:
    if (update) {
        i++;
        kbt = min(kb, m - i);
        i0 = i + 1;
        i1 = max((INTEGER)1, i - ka);
        i2 = i + kbt - ka1;
        if (i > m) {
            update = false;
            i = i - 1;
            i0 = m + 1;
            if (ka == 0) {
                return;
            }
            goto statement_490;
        }
    } else {
        i = i - ka;
        if (i < 2) {
            return;
        }
    }
    //
    if (i < m - kbt) {
        nx = m;
    } else {
        nx = n;
    }
    //
    if (upper) {
        //
        //        Transform A, working with the upper triangle
        //
        if (update) {
            //
            //           Form  inv(S(i))**H * A * inv(S(i))
            //
            bii = bb[(kb1 - 1) + (i - 1) * ldbb].real();
            ab[(ka1 - 1) + (i - 1) * ldab] = (ab[(ka1 - 1) + (i - 1) * ldab].real() / bii) / bii;
            for (j = i1; j <= i - 1; j = j + 1) {
                ab[((j - i + ka1) - 1) + (i - 1) * ldab] = ab[((j - i + ka1) - 1) + (i - 1) * ldab] / bii;
            }
            for (j = i + 1; j <= min(n, i + ka); j = j + 1) {
                ab[((i - j + ka1) - 1) + (j - 1) * ldab] = ab[((i - j + ka1) - 1) + (j - 1) * ldab] / bii;
            }
            for (k = i + 1; k <= i + kbt; k = k + 1) {
                for (j = k; j <= i + kbt; j = j + 1) {
                    ab[((k - j + ka1) - 1) + (j - 1) * ldab] = ab[((k - j + ka1) - 1) + (j - 1) * ldab] - bb[((i - j + kb1) - 1) + (j - 1) * ldbb] * conj(ab[((i - k + ka1) - 1) + (k - 1) * ldab]) - conj(bb[((i - k + kb1) - 1) + (k - 1) * ldbb]) * ab[((i - j + ka1) - 1) + (j - 1) * ldab] + ab[(ka1 - 1) + (i - 1) * ldab].real() * bb[((i - j + kb1) - 1) + (j - 1) * ldbb] * conj(bb[((i - k + kb1) - 1) + (k - 1) * ldbb]);
                }
                for (j = i + kbt + 1; j <= min(n, i + ka); j = j + 1) {
                    ab[((k - j + ka1) - 1) + (j - 1) * ldab] = ab[((k - j + ka1) - 1) + (j - 1) * ldab] - conj(bb[((i - k + kb1) - 1) + (k - 1) * ldbb]) * ab[((i - j + ka1) - 1) + (j - 1) * ldab];
                }
            }
            for (j = i1; j <= i; j = j + 1) {
                for (k = i + 1; k <= min(j + ka, i + kbt); k = k + 1) {
                    ab[((j - k + ka1) - 1) + (k - 1) * ldab] = ab[((j - k + ka1) - 1) + (k - 1) * ldab] - bb[((i - k + kb1) - 1) + (k - 1) * ldbb] * ab[((j - i + ka1) - 1) + (i - 1) * ldab];
                }
            }
            //
            if (wantx) {
                //
                //              post-multiply X by inv(S(i))
                //
                CRscal(nx, one / bii, x[(i - 1) * ldx], 1);
                if (kbt > 0) {
                    Cgeru(nx, kbt, -cone, x[(i - 1) * ldx], 1, bb[(kb - 1) + ((i + 1) - 1) * ldbb], ldbb - 1, x[((i + 1) - 1) * ldx], ldx);
                }
            }
            //
            //           store a(i1,i) in RA1 for use in next loop over K
            //
            ra1 = ab[((i1 - i + ka1) - 1) + (i - 1) * ldab];
        }
        //
        //        Generate and apply vectors of rotations to chase all the
        //        existing bulges KA positions up toward the top of the band
        //
        for (k = 1; k <= kb - 1; k = k + 1) {
            if (update) {
                //
                //              Determine the rotations which would annihilate the bulge
                //              which has in theory just been created
                //
                if (i + k - ka1 > 0 && i + k < m) {
                    //
                    //                 generate rotation to annihilate a(i+k-ka-1,i)
                    //
                    Clartg(ab[((k + 1) - 1) + (i - 1) * ldab], ra1, rwork[(i + k - ka) - 1], work[(i + k - ka) - 1], ra);
                    //
                    //                 create nonzero element a(i+k-ka-1,i+k) outside the
                    //                 band and store it in WORK(m-kb+i+k)
                    //
                    t = -bb[((kb1 - k) - 1) + ((i + k) - 1) * ldbb] * ra1;
                    work[(m - kb + i + k) - 1] = rwork[(i + k - ka) - 1] * t - conj(work[(i + k - ka) - 1]) * ab[((i + k) - 1) * ldab];
                    ab[((i + k) - 1) * ldab] = work[(i + k - ka) - 1] * t + rwork[(i + k - ka) - 1] * ab[((i + k) - 1) * ldab];
                    ra1 = ra;
                }
            }
            j2 = i + k + 1 - max((INTEGER)1, k + i0 - m + 1) * ka1;
            nr = (j2 + ka - 1) / ka1;
            j1 = j2 - (nr - 1) * ka1;
            if (update) {
                j2t = min(j2, i - 2 * ka + k - 1);
            } else {
                j2t = j2;
            }
            nrt = (j2t + ka - 1) / ka1;
            for (j = j1; j <= j2t; j = j + ka1) {
                //
                //              create nonzero element a(j-1,j+ka) outside the band
                //              and store it in WORK(j)
                //
                work[j - 1] = work[j - 1] * ab[((j + ka - 1) - 1) * ldab];
                ab[((j + ka - 1) - 1) * ldab] = rwork[j - 1] * ab[((j + ka - 1) - 1) * ldab];
            }
            //
            //           generate rotations in 1st set to annihilate elements which
            //           have been created outside the band
            //
            if (nrt > 0) {
                Clargv(nrt, ab[((j1 + ka) - 1) * ldab], inca, work[j1 - 1], ka1, rwork[j1 - 1], ka1);
            }
            if (nr > 0) {
                //
                //              apply rotations in 1st set from the left
                //
                for (l = 1; l <= ka - 1; l = l + 1) {
                    Clartv(nr, ab[((ka1 - l) - 1) + ((j1 + l) - 1) * ldab], inca, ab[((ka - l) - 1) + ((j1 + l) - 1) * ldab], inca, rwork[j1 - 1], work[j1 - 1], ka1);
                }
                //
                //              apply rotations in 1st set from both sides to diagonal
                //              blocks
                //
                Clar2v(nr, ab[(ka1 - 1) + (j1 - 1) * ldab], ab[(ka1 - 1) + ((j1 - 1) - 1) * ldab], ab[(ka - 1) + (j1 - 1) * ldab], inca, rwork[j1 - 1], work[j1 - 1], ka1);
                //
                Clacgv(nr, work[j1 - 1], ka1);
            }
            //
            //           start applying rotations in 1st set from the right
            //
            for (l = ka - 1; l >= kb - k + 1; l = l - 1) {
                nrt = (j2 + l - 1) / ka1;
                j1t = j2 - (nrt - 1) * ka1;
                if (nrt > 0) {
                    Clartv(nrt, ab[(l - 1) + (j1t - 1) * ldab], inca, ab[((l + 1) - 1) + ((j1t - 1) - 1) * ldab], inca, rwork[j1t - 1], work[j1t - 1], ka1);
                }
            }
            //
            if (wantx) {
                //
                //              post-multiply X by product of rotations in 1st set
                //
                for (j = j1; j <= j2; j = j + ka1) {
                    Crot(nx, x[(j - 1) * ldx], 1, x[((j - 1) - 1) * ldx], 1, rwork[j - 1], work[j - 1]);
                }
            }
        }
        //
        if (update) {
            if (i2 > 0 && kbt > 0) {
                //
                //              create nonzero element a(i+kbt-ka-1,i+kbt) outside the
                //              band and store it in WORK(m-kb+i+kbt)
                //
                work[(m - kb + i + kbt) - 1] = -bb[((kb1 - kbt) - 1) + ((i + kbt) - 1) * ldbb] * ra1;
            }
        }
        //
        for (k = kb; k >= 1; k = k - 1) {
            if (update) {
                j2 = i + k + 1 - max(2, k + i0 - m) * ka1;
            } else {
                j2 = i + k + 1 - max((INTEGER)1, k + i0 - m) * ka1;
            }
            //
            //           finish applying rotations in 2nd set from the right
            //
            for (l = kb - k; l >= 1; l = l - 1) {
                nrt = (j2 + ka + l - 1) / ka1;
                j1t = j2 - (nrt - 1) * ka1;
                if (nrt > 0) {
                    Clartv(nrt, ab[(l - 1) + ((j1t + ka) - 1) * ldab], inca, ab[((l + 1) - 1) + ((j1t + ka - 1) - 1) * ldab], inca, rwork[(m - kb + j1t + ka) - 1], work[(m - kb + j1t + ka) - 1], ka1);
                }
            }
            nr = (j2 + ka - 1) / ka1;
            j1 = j2 - (nr - 1) * ka1;
            for (j = j1; j <= j2; j = j + ka1) {
                work[(m - kb + j) - 1] = work[(m - kb + j + ka) - 1];
                rwork[(m - kb + j) - 1] = rwork[(m - kb + j + ka) - 1];
            }
            for (j = j1; j <= j2; j = j + ka1) {
                //
                //              create nonzero element a(j-1,j+ka) outside the band
                //              and store it in WORK(m-kb+j)
                //
                work[(m - kb + j) - 1] = work[(m - kb + j) - 1] * ab[((j + ka - 1) - 1) * ldab];
                ab[((j + ka - 1) - 1) * ldab] = rwork[(m - kb + j) - 1] * ab[((j + ka - 1) - 1) * ldab];
            }
            if (update) {
                if (i + k > ka1 && k <= kbt) {
                    work[(m - kb + i + k - ka) - 1] = work[(m - kb + i + k) - 1];
                }
            }
        }
        //
        for (k = kb; k >= 1; k = k - 1) {
            j2 = i + k + 1 - max((INTEGER)1, k + i0 - m) * ka1;
            nr = (j2 + ka - 1) / ka1;
            j1 = j2 - (nr - 1) * ka1;
            if (nr > 0) {
                //
                //              generate rotations in 2nd set to annihilate elements
                //              which have been created outside the band
                //
                Clargv(nr, ab[((j1 + ka) - 1) * ldab], inca, work[(m - kb + j1) - 1], ka1, rwork[(m - kb + j1) - 1], ka1);
                //
                //              apply rotations in 2nd set from the left
                //
                for (l = 1; l <= ka - 1; l = l + 1) {
                    Clartv(nr, ab[((ka1 - l) - 1) + ((j1 + l) - 1) * ldab], inca, ab[((ka - l) - 1) + ((j1 + l) - 1) * ldab], inca, rwork[(m - kb + j1) - 1], work[(m - kb + j1) - 1], ka1);
                }
                //
                //              apply rotations in 2nd set from both sides to diagonal
                //              blocks
                //
                Clar2v(nr, ab[(ka1 - 1) + (j1 - 1) * ldab], ab[(ka1 - 1) + ((j1 - 1) - 1) * ldab], ab[(ka - 1) + (j1 - 1) * ldab], inca, rwork[(m - kb + j1) - 1], work[(m - kb + j1) - 1], ka1);
                //
                Clacgv(nr, work[(m - kb + j1) - 1], ka1);
            }
            //
            //           start applying rotations in 2nd set from the right
            //
            for (l = ka - 1; l >= kb - k + 1; l = l - 1) {
                nrt = (j2 + l - 1) / ka1;
                j1t = j2 - (nrt - 1) * ka1;
                if (nrt > 0) {
                    Clartv(nrt, ab[(l - 1) + (j1t - 1) * ldab], inca, ab[((l + 1) - 1) + ((j1t - 1) - 1) * ldab], inca, rwork[(m - kb + j1t) - 1], work[(m - kb + j1t) - 1], ka1);
                }
            }
            //
            if (wantx) {
                //
                //              post-multiply X by product of rotations in 2nd set
                //
                for (j = j1; j <= j2; j = j + ka1) {
                    Crot(nx, x[(j - 1) * ldx], 1, x[((j - 1) - 1) * ldx], 1, rwork[(m - kb + j) - 1], work[(m - kb + j) - 1]);
                }
            }
        }
        //
        for (k = 1; k <= kb - 1; k = k + 1) {
            j2 = i + k + 1 - max((INTEGER)1, k + i0 - m + 1) * ka1;
            //
            //           finish applying rotations in 1st set from the right
            //
            for (l = kb - k; l >= 1; l = l - 1) {
                nrt = (j2 + l - 1) / ka1;
                j1t = j2 - (nrt - 1) * ka1;
                if (nrt > 0) {
                    Clartv(nrt, ab[(l - 1) + (j1t - 1) * ldab], inca, ab[((l + 1) - 1) + ((j1t - 1) - 1) * ldab], inca, rwork[j1t - 1], work[j1t - 1], ka1);
                }
            }
        }
        //
        if (kb > 1) {
            for (j = 2; j <= i2 - ka; j = j + 1) {
                rwork[j - 1] = rwork[(j + ka) - 1];
                work[j - 1] = work[(j + ka) - 1];
            }
        }
        //
    } else {
        //
        //        Transform A, working with the lower triangle
        //
        if (update) {
            //
            //           Form  inv(S(i))**H * A * inv(S(i))
            //
            bii = bb[(i - 1) * ldbb].real();
            ab[(i - 1) * ldab] = (ab[(i - 1) * ldab].real() / bii) / bii;
            for (j = i1; j <= i - 1; j = j + 1) {
                ab[((i - j + 1) - 1) + (j - 1) * ldab] = ab[((i - j + 1) - 1) + (j - 1) * ldab] / bii;
            }
            for (j = i + 1; j <= min(n, i + ka); j = j + 1) {
                ab[((j - i + 1) - 1) + (i - 1) * ldab] = ab[((j - i + 1) - 1) + (i - 1) * ldab] / bii;
            }
            for (k = i + 1; k <= i + kbt; k = k + 1) {
                for (j = k; j <= i + kbt; j = j + 1) {
                    ab[((j - k + 1) - 1) + (k - 1) * ldab] = ab[((j - k + 1) - 1) + (k - 1) * ldab] - bb[((j - i + 1) - 1) + (i - 1) * ldbb] * conj(ab[((k - i + 1) - 1) + (i - 1) * ldab]) - conj(bb[((k - i + 1) - 1) + (i - 1) * ldbb]) * ab[((j - i + 1) - 1) + (i - 1) * ldab] + ab[(i - 1) * ldab].real() * bb[((j - i + 1) - 1) + (i - 1) * ldbb] * conj(bb[((k - i + 1) - 1) + (i - 1) * ldbb]);
                }
                for (j = i + kbt + 1; j <= min(n, i + ka); j = j + 1) {
                    ab[((j - k + 1) - 1) + (k - 1) * ldab] = ab[((j - k + 1) - 1) + (k - 1) * ldab] - conj(bb[((k - i + 1) - 1) + (i - 1) * ldbb]) * ab[((j - i + 1) - 1) + (i - 1) * ldab];
                }
            }
            for (j = i1; j <= i; j = j + 1) {
                for (k = i + 1; k <= min(j + ka, i + kbt); k = k + 1) {
                    ab[((k - j + 1) - 1) + (j - 1) * ldab] = ab[((k - j + 1) - 1) + (j - 1) * ldab] - bb[((k - i + 1) - 1) + (i - 1) * ldbb] * ab[((i - j + 1) - 1) + (j - 1) * ldab];
                }
            }
            //
            if (wantx) {
                //
                //              post-multiply X by inv(S(i))
                //
                CRscal(nx, one / bii, x[(i - 1) * ldx], 1);
                if (kbt > 0) {
                    Cgerc(nx, kbt, -cone, x[(i - 1) * ldx], 1, bb[(2 - 1) + (i - 1) * ldbb], 1, x[((i + 1) - 1) * ldx], ldx);
                }
            }
            //
            //           store a(i,i1) in RA1 for use in next loop over K
            //
            ra1 = ab[((i - i1 + 1) - 1) + (i1 - 1) * ldab];
        }
        //
        //        Generate and apply vectors of rotations to chase all the
        //        existing bulges KA positions up toward the top of the band
        //
        for (k = 1; k <= kb - 1; k = k + 1) {
            if (update) {
                //
                //              Determine the rotations which would annihilate the bulge
                //              which has in theory just been created
                //
                if (i + k - ka1 > 0 && i + k < m) {
                    //
                    //                 generate rotation to annihilate a(i,i+k-ka-1)
                    //
                    Clartg(ab[((ka1 - k) - 1) + ((i + k - ka) - 1) * ldab], ra1, rwork[(i + k - ka) - 1], work[(i + k - ka) - 1], ra);
                    //
                    //                 create nonzero element a(i+k,i+k-ka-1) outside the
                    //                 band and store it in WORK(m-kb+i+k)
                    //
                    t = -bb[((k + 1) - 1) + (i - 1) * ldbb] * ra1;
                    work[(m - kb + i + k) - 1] = rwork[(i + k - ka) - 1] * t - conj(work[(i + k - ka) - 1]) * ab[(ka1 - 1) + ((i + k - ka) - 1) * ldab];
                    ab[(ka1 - 1) + ((i + k - ka) - 1) * ldab] = work[(i + k - ka) - 1] * t + rwork[(i + k - ka) - 1] * ab[(ka1 - 1) + ((i + k - ka) - 1) * ldab];
                    ra1 = ra;
                }
            }
            j2 = i + k + 1 - max((INTEGER)1, k + i0 - m + 1) * ka1;
            nr = (j2 + ka - 1) / ka1;
            j1 = j2 - (nr - 1) * ka1;
            if (update) {
                j2t = min(j2, i - 2 * ka + k - 1);
            } else {
                j2t = j2;
            }
            nrt = (j2t + ka - 1) / ka1;
            for (j = j1; j <= j2t; j = j + ka1) {
                //
                //              create nonzero element a(j+ka,j-1) outside the band
                //              and store it in WORK(j)
                //
                work[j - 1] = work[j - 1] * ab[(ka1 - 1) + ((j - 1) - 1) * ldab];
                ab[(ka1 - 1) + ((j - 1) - 1) * ldab] = rwork[j - 1] * ab[(ka1 - 1) + ((j - 1) - 1) * ldab];
            }
            //
            //           generate rotations in 1st set to annihilate elements which
            //           have been created outside the band
            //
            if (nrt > 0) {
                Clargv(nrt, ab[(ka1 - 1) + (j1 - 1) * ldab], inca, work[j1 - 1], ka1, rwork[j1 - 1], ka1);
            }
            if (nr > 0) {
                //
                //              apply rotations in 1st set from the right
                //
                for (l = 1; l <= ka - 1; l = l + 1) {
                    Clartv(nr, ab[((l + 1) - 1) + (j1 - 1) * ldab], inca, ab[((l + 2) - 1) + ((j1 - 1) - 1) * ldab], inca, rwork[j1 - 1], work[j1 - 1], ka1);
                }
                //
                //              apply rotations in 1st set from both sides to diagonal
                //              blocks
                //
                Clar2v(nr, ab[(j1 - 1) * ldab], ab[((j1 - 1) - 1) * ldab], ab[(2 - 1) + ((j1 - 1) - 1) * ldab], inca, rwork[j1 - 1], work[j1 - 1], ka1);
                //
                Clacgv(nr, work[j1 - 1], ka1);
            }
            //
            //           start applying rotations in 1st set from the left
            //
            for (l = ka - 1; l >= kb - k + 1; l = l - 1) {
                nrt = (j2 + l - 1) / ka1;
                j1t = j2 - (nrt - 1) * ka1;
                if (nrt > 0) {
                    Clartv(nrt, ab[((ka1 - l + 1) - 1) + ((j1t - ka1 + l) - 1) * ldab], inca, ab[((ka1 - l) - 1) + ((j1t - ka1 + l) - 1) * ldab], inca, rwork[j1t - 1], work[j1t - 1], ka1);
                }
            }
            //
            if (wantx) {
                //
                //              post-multiply X by product of rotations in 1st set
                //
                for (j = j1; j <= j2; j = j + ka1) {
                    Crot(nx, x[(j - 1) * ldx], 1, x[((j - 1) - 1) * ldx], 1, rwork[j - 1], conj(work[j - 1]));
                }
            }
        }
        //
        if (update) {
            if (i2 > 0 && kbt > 0) {
                //
                //              create nonzero element a(i+kbt,i+kbt-ka-1) outside the
                //              band and store it in WORK(m-kb+i+kbt)
                //
                work[(m - kb + i + kbt) - 1] = -bb[((kbt + 1) - 1) + (i - 1) * ldbb] * ra1;
            }
        }
        //
        for (k = kb; k >= 1; k = k - 1) {
            if (update) {
                j2 = i + k + 1 - max(2, k + i0 - m) * ka1;
            } else {
                j2 = i + k + 1 - max((INTEGER)1, k + i0 - m) * ka1;
            }
            //
            //           finish applying rotations in 2nd set from the left
            //
            for (l = kb - k; l >= 1; l = l - 1) {
                nrt = (j2 + ka + l - 1) / ka1;
                j1t = j2 - (nrt - 1) * ka1;
                if (nrt > 0) {
                    Clartv(nrt, ab[((ka1 - l + 1) - 1) + ((j1t + l - 1) - 1) * ldab], inca, ab[((ka1 - l) - 1) + ((j1t + l - 1) - 1) * ldab], inca, rwork[(m - kb + j1t + ka) - 1], work[(m - kb + j1t + ka) - 1], ka1);
                }
            }
            nr = (j2 + ka - 1) / ka1;
            j1 = j2 - (nr - 1) * ka1;
            for (j = j1; j <= j2; j = j + ka1) {
                work[(m - kb + j) - 1] = work[(m - kb + j + ka) - 1];
                rwork[(m - kb + j) - 1] = rwork[(m - kb + j + ka) - 1];
            }
            for (j = j1; j <= j2; j = j + ka1) {
                //
                //              create nonzero element a(j+ka,j-1) outside the band
                //              and store it in WORK(m-kb+j)
                //
                work[(m - kb + j) - 1] = work[(m - kb + j) - 1] * ab[(ka1 - 1) + ((j - 1) - 1) * ldab];
                ab[(ka1 - 1) + ((j - 1) - 1) * ldab] = rwork[(m - kb + j) - 1] * ab[(ka1 - 1) + ((j - 1) - 1) * ldab];
            }
            if (update) {
                if (i + k > ka1 && k <= kbt) {
                    work[(m - kb + i + k - ka) - 1] = work[(m - kb + i + k) - 1];
                }
            }
        }
        //
        for (k = kb; k >= 1; k = k - 1) {
            j2 = i + k + 1 - max((INTEGER)1, k + i0 - m) * ka1;
            nr = (j2 + ka - 1) / ka1;
            j1 = j2 - (nr - 1) * ka1;
            if (nr > 0) {
                //
                //              generate rotations in 2nd set to annihilate elements
                //              which have been created outside the band
                //
                Clargv(nr, ab[(ka1 - 1) + (j1 - 1) * ldab], inca, work[(m - kb + j1) - 1], ka1, rwork[(m - kb + j1) - 1], ka1);
                //
                //              apply rotations in 2nd set from the right
                //
                for (l = 1; l <= ka - 1; l = l + 1) {
                    Clartv(nr, ab[((l + 1) - 1) + (j1 - 1) * ldab], inca, ab[((l + 2) - 1) + ((j1 - 1) - 1) * ldab], inca, rwork[(m - kb + j1) - 1], work[(m - kb + j1) - 1], ka1);
                }
                //
                //              apply rotations in 2nd set from both sides to diagonal
                //              blocks
                //
                Clar2v(nr, ab[(j1 - 1) * ldab], ab[((j1 - 1) - 1) * ldab], ab[(2 - 1) + ((j1 - 1) - 1) * ldab], inca, rwork[(m - kb + j1) - 1], work[(m - kb + j1) - 1], ka1);
                //
                Clacgv(nr, work[(m - kb + j1) - 1], ka1);
            }
            //
            //           start applying rotations in 2nd set from the left
            //
            for (l = ka - 1; l >= kb - k + 1; l = l - 1) {
                nrt = (j2 + l - 1) / ka1;
                j1t = j2 - (nrt - 1) * ka1;
                if (nrt > 0) {
                    Clartv(nrt, ab[((ka1 - l + 1) - 1) + ((j1t - ka1 + l) - 1) * ldab], inca, ab[((ka1 - l) - 1) + ((j1t - ka1 + l) - 1) * ldab], inca, rwork[(m - kb + j1t) - 1], work[(m - kb + j1t) - 1], ka1);
                }
            }
            //
            if (wantx) {
                //
                //              post-multiply X by product of rotations in 2nd set
                //
                for (j = j1; j <= j2; j = j + ka1) {
                    Crot(nx, x[(j - 1) * ldx], 1, x[((j - 1) - 1) * ldx], 1, rwork[(m - kb + j) - 1], conj(work[(m - kb + j) - 1]));
                }
            }
        }
        //
        for (k = 1; k <= kb - 1; k = k + 1) {
            j2 = i + k + 1 - max((INTEGER)1, k + i0 - m + 1) * ka1;
            //
            //           finish applying rotations in 1st set from the left
            //
            for (l = kb - k; l >= 1; l = l - 1) {
                nrt = (j2 + l - 1) / ka1;
                j1t = j2 - (nrt - 1) * ka1;
                if (nrt > 0) {
                    Clartv(nrt, ab[((ka1 - l + 1) - 1) + ((j1t - ka1 + l) - 1) * ldab], inca, ab[((ka1 - l) - 1) + ((j1t - ka1 + l) - 1) * ldab], inca, rwork[j1t - 1], work[j1t - 1], ka1);
                }
            }
        }
        //
        if (kb > 1) {
            for (j = 2; j <= i2 - ka; j = j + 1) {
                rwork[j - 1] = rwork[(j + ka) - 1];
                work[j - 1] = work[(j + ka) - 1];
            }
        }
        //
    }
    //
    goto statement_490;
    //
    //     End of Chbgst
    //
}
