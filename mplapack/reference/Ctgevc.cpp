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

void Ctgevc(const char *side, const char *howmny, arr_cref<bool> select, INTEGER const &n, COMPLEX *s, INTEGER const &lds, COMPLEX *p, INTEGER const &ldp, COMPLEX *vl, INTEGER const &ldvl, COMPLEX *vr, INTEGER const &ldvr, INTEGER const &mm, INTEGER &m, COMPLEX *work, REAL *rwork, INTEGER &info) {
    COMPLEX x = 0.0;
    INTEGER ihwmny = 0;
    bool ilall = false;
    bool ilback = false;
    INTEGER iside = 0;
    bool identifier_compl = false;
    bool compr = false;
    INTEGER im = 0;
    INTEGER j = 0;
    bool ilbbad = false;
    const REAL zero = 0.0;
    REAL safmin = 0.0;
    const REAL one = 1.0;
    REAL big = 0.0;
    REAL ulp = 0.0;
    REAL small = 0.0;
    REAL bignum = 0.0;
    REAL anorm = 0.0;
    REAL bnorm = 0.0;
    INTEGER i = 0;
    REAL ascale = 0.0;
    REAL bscale = 0.0;
    INTEGER ieig = 0;
    INTEGER je = 0;
    bool ilcomp = false;
    INTEGER jr = 0;
    const COMPLEX czero = (0.0, 0.0);
    const COMPLEX cone = (1.0, 0.0);
    REAL temp = 0.0;
    COMPLEX salpha = 0.0;
    REAL sbeta = 0.0;
    REAL acoeff = 0.0;
    COMPLEX bcoeff = 0.0;
    bool lsa = false;
    bool lsb = false;
    REAL scale = 0.0;
    REAL acoefa = 0.0;
    REAL bcoefa = 0.0;
    REAL xmax = 0.0;
    REAL dmin = 0.0;
    COMPLEX suma = 0.0;
    COMPLEX sumb = 0.0;
    COMPLEX sum = 0.0;
    COMPLEX d = 0.0;
    INTEGER isrc = 0;
    INTEGER ibeg = 0;
    COMPLEX ca = 0.0;
    COMPLEX cb = 0.0;
    INTEGER iend = 0;
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
    //     .. Statement Functions ..
    //     ..
    //     .. Statement Function definitions ..
    abs1[x - 1] = abs(x.real()) + abs(x.imag());
    //     ..
    //     .. Executable Statements ..
    //
    //     Decode and Test the input parameters
    //
    if (Mlsame(howmny, "A")) {
        ihwmny = 1;
        ilall = true;
        ilback = false;
    } else if (Mlsame(howmny, "S")) {
        ihwmny = 2;
        ilall = false;
        ilback = false;
    } else if (Mlsame(howmny, "B")) {
        ihwmny = 3;
        ilall = true;
        ilback = true;
    } else {
        ihwmny = -1;
    }
    //
    if (Mlsame(side, "R")) {
        iside = 1;
        identifier_compl = false;
        compr = true;
    } else if (Mlsame(side, "L")) {
        iside = 2;
        identifier_compl = true;
        compr = false;
    } else if (Mlsame(side, "B")) {
        iside = 3;
        identifier_compl = true;
        compr = true;
    } else {
        iside = -1;
    }
    //
    info = 0;
    if (iside < 0) {
        info = -1;
    } else if (ihwmny < 0) {
        info = -2;
    } else if (n < 0) {
        info = -4;
    } else if (lds < max((INTEGER)1, n)) {
        info = -6;
    } else if (ldp < max((INTEGER)1, n)) {
        info = -8;
    }
    if (info != 0) {
        Mxerbla("Ctgevc", -info);
        return;
    }
    //
    //     Count the number of eigenvectors
    //
    if (!ilall) {
        im = 0;
        for (j = 1; j <= n; j = j + 1) {
            if (select[j - 1]) {
                im++;
            }
        }
    } else {
        im = n;
    }
    //
    //     Check diagonal of B
    //
    ilbbad = false;
    for (j = 1; j <= n; j = j + 1) {
        if (p[(j - 1) + (j - 1) * ldp].imag() != zero) {
            ilbbad = true;
        }
    }
    //
    if (ilbbad) {
        info = -7;
    } else if (identifier_compl && ldvl < n || ldvl < 1) {
        info = -10;
    } else if (compr && ldvr < n || ldvr < 1) {
        info = -12;
    } else if (mm < im) {
        info = -13;
    }
    if (info != 0) {
        Mxerbla("Ctgevc", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    m = im;
    if (n == 0) {
        return;
    }
    //
    //     Machine Constants
    //
    safmin = dlamch("Safe minimum");
    big = one / safmin;
    Rlabad(safmin, big);
    ulp = dlamch("Epsilon") * dlamch("Base");
    small = safmin * n / ulp;
    big = one / small;
    bignum = one / (safmin * n);
    //
    //     Compute the 1-norm of each column of the strictly upper triangular
    //     part of A and B to check for possible overflow in the triangular
    //     solver.
    //
    anorm = abs1[s[(1 - 1)] - 1];
    bnorm = abs1[p[(1 - 1)] - 1];
    rwork[1 - 1] = zero;
    rwork[(n + 1) - 1] = zero;
    for (j = 2; j <= n; j = j + 1) {
        rwork[j - 1] = zero;
        rwork[(n + j) - 1] = zero;
        for (i = 1; i <= j - 1; i = i + 1) {
            rwork[j - 1] += abs1[s[(i - 1) + (j - 1) * lds] - 1];
            rwork[(n + j) - 1] += abs1[p[(i - 1) + (j - 1) * ldp] - 1];
        }
        anorm = max(anorm, rwork[j - 1] + abs1[s[(j - 1) + (j - 1) * lds] - 1]);
        bnorm = max(bnorm, rwork[(n + j) - 1] + abs1[p[(j - 1) + (j - 1) * ldp] - 1]);
    }
    //
    ascale = one / max(anorm, safmin);
    bscale = one / max(bnorm, safmin);
    //
    //     Left eigenvectors
    //
    if (identifier_compl) {
        ieig = 0;
        //
        //        Main loop over eigenvalues
        //
        for (je = 1; je <= n; je = je + 1) {
            if (ilall) {
                ilcomp = true;
            } else {
                ilcomp = select[je - 1];
            }
            if (ilcomp) {
                ieig++;
                //
                if (abs1[s[(je - 1) + (je - 1) * lds] - 1] <= safmin && abs(p[(je - 1) + (je - 1) * ldp].real()) <= safmin) {
                    //
                    //                 Singular matrix pencil -- return unit eigenvector
                    //
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        vl[(jr - 1) + (ieig - 1) * ldvl] = czero;
                    }
                    vl[(ieig - 1) + (ieig - 1) * ldvl] = cone;
                    goto statement_140;
                }
                //
                //              Non-singular eigenvalue:
                //              Compute coefficients  a  and  b  in
                //                   H
                //                 y  ( a A - b B ) = 0
                //
                temp = one / max(abs1[s[(je - 1) + (je - 1) * lds] - 1] * ascale, abs(p[(je - 1) + (je - 1) * ldp].real()) * bscale, safmin);
                salpha = (temp * s[(je - 1) + (je - 1) * lds]) * ascale;
                sbeta = (temp * p[(je - 1) + (je - 1) * ldp].real()) * bscale;
                acoeff = sbeta * ascale;
                bcoeff = salpha * bscale;
                //
                //              Scale to avoid underflow
                //
                lsa = abs(sbeta) >= safmin && abs(acoeff) < small;
                lsb = abs1[salpha - 1] >= safmin && abs1[bcoeff - 1] < small;
                //
                scale = one;
                if (lsa) {
                    scale = (small / abs(sbeta)) * min(anorm, big);
                }
                if (lsb) {
                    scale = max(scale, (small / abs1[salpha - 1]) * min(bnorm, big));
                }
                if (lsa || lsb) {
                    scale = min(scale, one / (safmin * max(one, abs(acoeff), abs1[bcoeff - 1])));
                    if (lsa) {
                        acoeff = ascale * (scale * sbeta);
                    } else {
                        acoeff = scale * acoeff;
                    }
                    if (lsb) {
                        bcoeff = bscale * (scale * salpha);
                    } else {
                        bcoeff = scale * bcoeff;
                    }
                }
                //
                acoefa = abs(acoeff);
                bcoefa = abs1[bcoeff - 1];
                xmax = one;
                for (jr = 1; jr <= n; jr = jr + 1) {
                    work[jr - 1] = czero;
                }
                work[je - 1] = cone;
                dmin = max(ulp * acoefa * anorm, ulp * bcoefa * bnorm, safmin);
                //
                //                                              H
                //              Triangular solve of  (a A - b B)  y = 0
                //
                //                                      H
                //              (rowwise in  (a A - b B) , or columnwise in a A - b B)
                //
                for (j = je + 1; j <= n; j = j + 1) {
                    //
                    //                 Compute
                    //                       j-1
                    //                 SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k)
                    //                       k=je
                    //                 (Scale if necessary)
                    //
                    temp = one / xmax;
                    if (acoefa * rwork[j - 1] + bcoefa * rwork[(n + j) - 1] > bignum * temp) {
                        for (jr = je; jr <= j - 1; jr = jr + 1) {
                            work[jr - 1] = temp * work[jr - 1];
                        }
                        xmax = one;
                    }
                    suma = czero;
                    sumb = czero;
                    //
                    for (jr = je; jr <= j - 1; jr = jr + 1) {
                        suma += conj(s[(jr - 1) + (j - 1) * lds]) * work[jr - 1];
                        sumb += conj(p[(jr - 1) + (j - 1) * ldp]) * work[jr - 1];
                    }
                    sum = acoeff * suma - conj(bcoeff) * sumb;
                    //
                    //                 Form x(j) = - SUM / conjg( a*S(j,j) - b*P(j,j) )
                    //
                    //                 with scaling and perturbation of the denominator
                    //
                    d = conj(acoeff * s[(j - 1) + (j - 1) * lds] - bcoeff * p[(j - 1) + (j - 1) * ldp]);
                    if (abs1[d - 1] <= dmin) {
                        d = COMPLEX(dmin);
                    }
                    //
                    if (abs1[d - 1] < one) {
                        if (abs1[sum - 1] >= bignum * abs1[d - 1]) {
                            temp = one / abs1[sum - 1];
                            for (jr = je; jr <= j - 1; jr = jr + 1) {
                                work[jr - 1] = temp * work[jr - 1];
                            }
                            xmax = temp * xmax;
                            sum = temp * sum;
                        }
                    }
                    work[j - 1] = Cladiv[(-sum - 1) + (d - 1) * ldCladiv];
                    xmax = max(xmax, abs1[work[j - 1] - 1]);
                }
                //
                //              Back transform eigenvector if HOWMNY='B'.
                //
                if (ilback) {
                    Cgemv("N", n, n + 1 - je, cone, vl[(je - 1) * ldvl], ldvl, work[je - 1], 1, czero, work[(n + 1) - 1], 1);
                    isrc = 2;
                    ibeg = 1;
                } else {
                    isrc = 1;
                    ibeg = je;
                }
                //
                //              Copy and scale eigenvector INTEGERo column of VL
                //
                xmax = zero;
                for (jr = ibeg; jr <= n; jr = jr + 1) {
                    xmax = max(xmax, abs1[(work[((isrc - 1) * n + jr) - 1]) - 1]);
                }
                //
                if (xmax > safmin) {
                    temp = one / xmax;
                    for (jr = ibeg; jr <= n; jr = jr + 1) {
                        vl[(jr - 1) + (ieig - 1) * ldvl] = temp * work[((isrc - 1) * n + jr) - 1];
                    }
                } else {
                    ibeg = n + 1;
                }
                //
                for (jr = 1; jr <= ibeg - 1; jr = jr + 1) {
                    vl[(jr - 1) + (ieig - 1) * ldvl] = czero;
                }
                //
            }
        statement_140:;
        }
    }
    //
    //     Right eigenvectors
    //
    if (compr) {
        ieig = im + 1;
        //
        //        Main loop over eigenvalues
        //
        for (je = n; je >= 1; je = je - 1) {
            if (ilall) {
                ilcomp = true;
            } else {
                ilcomp = select[je - 1];
            }
            if (ilcomp) {
                ieig = ieig - 1;
                //
                if (abs1[s[(je - 1) + (je - 1) * lds] - 1] <= safmin && abs(p[(je - 1) + (je - 1) * ldp].real()) <= safmin) {
                    //
                    //                 Singular matrix pencil -- return unit eigenvector
                    //
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        vr[(jr - 1) + (ieig - 1) * ldvr] = czero;
                    }
                    vr[(ieig - 1) + (ieig - 1) * ldvr] = cone;
                    goto statement_250;
                }
                //
                //              Non-singular eigenvalue:
                //              Compute coefficients  a  and  b  in
                //
                //              ( a A - b B ) x  = 0
                //
                temp = one / max(abs1[s[(je - 1) + (je - 1) * lds] - 1] * ascale, abs(p[(je - 1) + (je - 1) * ldp].real()) * bscale, safmin);
                salpha = (temp * s[(je - 1) + (je - 1) * lds]) * ascale;
                sbeta = (temp * p[(je - 1) + (je - 1) * ldp].real()) * bscale;
                acoeff = sbeta * ascale;
                bcoeff = salpha * bscale;
                //
                //              Scale to avoid underflow
                //
                lsa = abs(sbeta) >= safmin && abs(acoeff) < small;
                lsb = abs1[salpha - 1] >= safmin && abs1[bcoeff - 1] < small;
                //
                scale = one;
                if (lsa) {
                    scale = (small / abs(sbeta)) * min(anorm, big);
                }
                if (lsb) {
                    scale = max(scale, (small / abs1[salpha - 1]) * min(bnorm, big));
                }
                if (lsa || lsb) {
                    scale = min(scale, one / (safmin * max(one, abs(acoeff), abs1[bcoeff - 1])));
                    if (lsa) {
                        acoeff = ascale * (scale * sbeta);
                    } else {
                        acoeff = scale * acoeff;
                    }
                    if (lsb) {
                        bcoeff = bscale * (scale * salpha);
                    } else {
                        bcoeff = scale * bcoeff;
                    }
                }
                //
                acoefa = abs(acoeff);
                bcoefa = abs1[bcoeff - 1];
                xmax = one;
                for (jr = 1; jr <= n; jr = jr + 1) {
                    work[jr - 1] = czero;
                }
                work[je - 1] = cone;
                dmin = max(ulp * acoefa * anorm, ulp * bcoefa * bnorm, safmin);
                //
                //              Triangular solve of  (a A - b B) x = 0  (columnwise)
                //
                //              WORK(1:j-1) contains sums w,
                //              WORK(j+1:JE) contains x
                //
                for (jr = 1; jr <= je - 1; jr = jr + 1) {
                    work[jr - 1] = acoeff * s[(jr - 1) + (je - 1) * lds] - bcoeff * p[(jr - 1) + (je - 1) * ldp];
                }
                work[je - 1] = cone;
                //
                for (j = je - 1; j >= 1; j = j - 1) {
                    //
                    //                 Form x(j) := - w(j) / d
                    //                 with scaling and perturbation of the denominator
                    //
                    d = acoeff * s[(j - 1) + (j - 1) * lds] - bcoeff * p[(j - 1) + (j - 1) * ldp];
                    if (abs1[d - 1] <= dmin) {
                        d = COMPLEX(dmin);
                    }
                    //
                    if (abs1[d - 1] < one) {
                        if (abs1[work[j - 1] - 1] >= bignum * abs1[d - 1]) {
                            temp = one / abs1[work[j - 1] - 1];
                            for (jr = 1; jr <= je; jr = jr + 1) {
                                work[jr - 1] = temp * work[jr - 1];
                            }
                        }
                    }
                    //
                    work[j - 1] = Cladiv[(-work[j - 1] - 1) + (d - 1) * ldCladiv];
                    //
                    if (j > 1) {
                        //
                        //                    w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling
                        //
                        if (abs1[work[j - 1] - 1] > one) {
                            temp = one / abs1[work[j - 1] - 1];
                            if (acoefa * rwork[j - 1] + bcoefa * rwork[(n + j) - 1] >= bignum * temp) {
                                for (jr = 1; jr <= je; jr = jr + 1) {
                                    work[jr - 1] = temp * work[jr - 1];
                                }
                            }
                        }
                        //
                        ca = acoeff * work[j - 1];
                        cb = bcoeff * work[j - 1];
                        for (jr = 1; jr <= j - 1; jr = jr + 1) {
                            work[jr - 1] += ca * s[(jr - 1) + (j - 1) * lds] - cb * p[(jr - 1) + (j - 1) * ldp];
                        }
                    }
                }
                //
                //              Back transform eigenvector if HOWMNY='B'.
                //
                if (ilback) {
                    Cgemv("N", n, je, cone, vr, ldvr, work, 1, czero, work[(n + 1) - 1], 1);
                    isrc = 2;
                    iend = n;
                } else {
                    isrc = 1;
                    iend = je;
                }
                //
                //              Copy and scale eigenvector INTEGERo column of VR
                //
                xmax = zero;
                for (jr = 1; jr <= iend; jr = jr + 1) {
                    xmax = max(xmax, abs1[(work[((isrc - 1) * n + jr) - 1]) - 1]);
                }
                //
                if (xmax > safmin) {
                    temp = one / xmax;
                    for (jr = 1; jr <= iend; jr = jr + 1) {
                        vr[(jr - 1) + (ieig - 1) * ldvr] = temp * work[((isrc - 1) * n + jr) - 1];
                    }
                } else {
                    iend = 0;
                }
                //
                for (jr = iend + 1; jr <= n; jr = jr + 1) {
                    vr[(jr - 1) + (ieig - 1) * ldvr] = czero;
                }
                //
            }
        statement_250:;
        }
    }
    //
    //     End of Ctgevc
    //
}
