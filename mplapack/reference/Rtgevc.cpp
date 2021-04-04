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

void Rtgevc(const char *side, const char *howmny, arr_cref<bool> select, INTEGER const &n, REAL *s, INTEGER const &lds, REAL *p, INTEGER const &ldp, REAL *vl, INTEGER const &ldvl, REAL *vr, INTEGER const &ldvr, INTEGER const &mm, INTEGER &m, REAL *work, INTEGER &info) {
    INTEGER ihwmny = 0;
    bool ilall = false;
    bool ilback = false;
    INTEGER iside = 0;
    bool identifier_compl = false;
    bool compr = false;
    INTEGER im = 0;
    bool ilcplx = false;
    INTEGER j = 0;
    const REAL zero = 0.0;
    bool ilabad = false;
    bool ilbbad = false;
    REAL safmin = 0.0;
    const REAL one = 1.0;
    REAL big = 0.0;
    REAL ulp = 0.0;
    REAL small = 0.0;
    REAL bignum = 0.0;
    REAL anorm = 0.0;
    REAL bnorm = 0.0;
    REAL temp = 0.0;
    REAL temp2 = 0.0;
    INTEGER iend = 0;
    INTEGER i = 0;
    REAL ascale = 0.0;
    REAL bscale = 0.0;
    INTEGER ieig = 0;
    INTEGER je = 0;
    INTEGER nw = 0;
    bool ilcomp = false;
    INTEGER jr = 0;
    REAL salfar = 0.0;
    REAL sbeta = 0.0;
    REAL acoef = 0.0;
    REAL bcoefr = 0.0;
    REAL bcoefi = 0.0;
    REAL scale = 0.0;
    bool lsa = false;
    bool lsb = false;
    REAL acoefa = 0.0;
    REAL bcoefa = 0.0;
    REAL xmax = 0.0;
    const REAL safety = 1.0e+2;
    REAL temp2r = 0.0;
    REAL temp2i = 0.0;
    REAL dmin = 0.0;
    bool il2by2 = false;
    INTEGER na = 0;
    arr_1d<2, REAL> bdiag(fill0);
    REAL xscale = 0.0;
    INTEGER jw = 0;
    INTEGER ja = 0;
    arr_2d<2, 2, REAL> sums(fill0);
    arr_2d<2, 2, REAL> sump(fill0);
    arr_2d<2, 2, REAL> sum(fill0);
    INTEGER iinfo = 0;
    INTEGER ibeg = 0;
    REAL creala = 0.0;
    REAL cimaga = 0.0;
    REAL crealb = 0.0;
    REAL cimagb = 0.0;
    REAL cre2a = 0.0;
    REAL cim2a = 0.0;
    REAL cre2b = 0.0;
    REAL cim2b = 0.0;
    INTEGER jc = 0;
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
    //     .. Local Arrays ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
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
        ilall = true;
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
        Mxerbla("Rtgevc", -info);
        return;
    }
    //
    //     Count the number of eigenvectors to be computed
    //
    if (!ilall) {
        im = 0;
        ilcplx = false;
        for (j = 1; j <= n; j = j + 1) {
            if (ilcplx) {
                ilcplx = false;
                goto statement_10;
            }
            if (j < n) {
                if (s[((j + 1) - 1) + (j - 1) * lds] != zero) {
                    ilcplx = true;
                }
            }
            if (ilcplx) {
                if (select[j - 1] || select[(j + 1) - 1]) {
                    im += 2;
                }
            } else {
                if (select[j - 1]) {
                    im++;
                }
            }
        statement_10:;
        }
    } else {
        im = n;
    }
    //
    //     Check 2-by-2 diagonal blocks of A, B
    //
    ilabad = false;
    ilbbad = false;
    for (j = 1; j <= n - 1; j = j + 1) {
        if (s[((j + 1) - 1) + (j - 1) * lds] != zero) {
            if (p[(j - 1) + (j - 1) * ldp] == zero || p[((j + 1) - 1) + ((j + 1) - 1) * ldp] == zero || p[(j - 1) + ((j + 1) - 1) * ldp] != zero) {
                ilbbad = true;
            }
            if (j < n - 1) {
                if (s[((j + 2) - 1) + ((j + 1) - 1) * lds] != zero) {
                    ilabad = true;
                }
            }
        }
    }
    //
    if (ilabad) {
        info = -5;
    } else if (ilbbad) {
        info = -7;
    } else if (identifier_compl && ldvl < n || ldvl < 1) {
        info = -10;
    } else if (compr && ldvr < n || ldvr < 1) {
        info = -12;
    } else if (mm < im) {
        info = -13;
    }
    if (info != 0) {
        Mxerbla("Rtgevc", -info);
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
    //     part (i.e., excluding all elements belonging to the diagonal
    //     blocks) of A and B to check for possible overflow in the
    //     triangular solver.
    //
    anorm = abs(s[(1 - 1)]);
    if (n > 1) {
        anorm += abs(s[(2 - 1)]);
    }
    bnorm = abs(p[(1 - 1)]);
    work[1 - 1] = zero;
    work[(n + 1) - 1] = zero;
    //
    for (j = 2; j <= n; j = j + 1) {
        temp = zero;
        temp2 = zero;
        if (s[(j - 1) + ((j - 1) - 1) * lds] == zero) {
            iend = j - 1;
        } else {
            iend = j - 2;
        }
        for (i = 1; i <= iend; i = i + 1) {
            temp += abs(s[(i - 1) + (j - 1) * lds]);
            temp2 += abs(p[(i - 1) + (j - 1) * ldp]);
        }
        work[j - 1] = temp;
        work[(n + j) - 1] = temp2;
        for (i = iend + 1; i <= min(j + 1, n); i = i + 1) {
            temp += abs(s[(i - 1) + (j - 1) * lds]);
            temp2 += abs(p[(i - 1) + (j - 1) * ldp]);
        }
        anorm = max(anorm, temp);
        bnorm = max(bnorm, temp2);
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
        ilcplx = false;
        for (je = 1; je <= n; je = je + 1) {
            //
            //           Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or
            //           (b) this would be the second of a complex pair.
            //           Check for complex eigenvalue, so as to be sure of which
            //           entry(-ies) of SELECT to look at.
            //
            if (ilcplx) {
                ilcplx = false;
                goto statement_220;
            }
            nw = 1;
            if (je < n) {
                if (s[((je + 1) - 1) + (je - 1) * lds] != zero) {
                    ilcplx = true;
                    nw = 2;
                }
            }
            if (ilall) {
                ilcomp = true;
            } else if (ilcplx) {
                ilcomp = select[je - 1] || select[(je + 1) - 1];
            } else {
                ilcomp = select[je - 1];
            }
            if (!ilcomp) {
                goto statement_220;
            }
            //
            //           Decide if (a) singular pencil, (b) real eigenvalue, or
            //           (c) complex eigenvalue.
            //
            if (!ilcplx) {
                if (abs(s[(je - 1) + (je - 1) * lds]) <= safmin && abs(p[(je - 1) + (je - 1) * ldp]) <= safmin) {
                    //
                    //                 Singular matrix pencil -- return unit eigenvector
                    //
                    ieig++;
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        vl[(jr - 1) + (ieig - 1) * ldvl] = zero;
                    }
                    vl[(ieig - 1) + (ieig - 1) * ldvl] = one;
                    goto statement_220;
                }
            }
            //
            //           Clear vector
            //
            for (jr = 1; jr <= nw * n; jr = jr + 1) {
                work[(2 * n + jr) - 1] = zero;
            }
            //                                                 T
            //           Compute coefficients in  ( a A - b B )  y = 0
            //              a  is  ACOEF
            //              b  is  BCOEFR + i*BCOEFI
            //
            if (!ilcplx) {
                //
                //              Real eigenvalue
                //
                temp = one / max(abs(s[(je - 1) + (je - 1) * lds]) * ascale, abs(p[(je - 1) + (je - 1) * ldp]) * bscale, safmin);
                salfar = (temp * s[(je - 1) + (je - 1) * lds]) * ascale;
                sbeta = (temp * p[(je - 1) + (je - 1) * ldp]) * bscale;
                acoef = sbeta * ascale;
                bcoefr = salfar * bscale;
                bcoefi = zero;
                //
                //              Scale to avoid underflow
                //
                scale = one;
                lsa = abs(sbeta) >= safmin && abs(acoef) < small;
                lsb = abs(salfar) >= safmin && abs(bcoefr) < small;
                if (lsa) {
                    scale = (small / abs(sbeta)) * min(anorm, big);
                }
                if (lsb) {
                    scale = max(scale, (small / abs(salfar)) * min(bnorm, big));
                }
                if (lsa || lsb) {
                    scale = min(scale, one / (safmin * max(one, abs(acoef), abs(bcoefr))));
                    if (lsa) {
                        acoef = ascale * (scale * sbeta);
                    } else {
                        acoef = scale * acoef;
                    }
                    if (lsb) {
                        bcoefr = bscale * (scale * salfar);
                    } else {
                        bcoefr = scale * bcoefr;
                    }
                }
                acoefa = abs(acoef);
                bcoefa = abs(bcoefr);
                //
                //              First component is 1
                //
                work[(2 * n + je) - 1] = one;
                xmax = one;
            } else {
                //
                //              Complex eigenvalue
                //
                Rlag2(s[(je - 1) + (je - 1) * lds], lds, p[(je - 1) + (je - 1) * ldp], ldp, safmin * safety, acoef, temp, bcoefr, temp2, bcoefi);
                bcoefi = -bcoefi;
                if (bcoefi == zero) {
                    info = je;
                    return;
                }
                //
                //              Scale to avoid over/underflow
                //
                acoefa = abs(acoef);
                bcoefa = abs(bcoefr) + abs(bcoefi);
                scale = one;
                if (acoefa * ulp < safmin && acoefa >= safmin) {
                    scale = (safmin / ulp) / acoefa;
                }
                if (bcoefa * ulp < safmin && bcoefa >= safmin) {
                    scale = max(scale, (safmin / ulp) / bcoefa);
                }
                if (safmin * acoefa > ascale) {
                    scale = ascale / (safmin * acoefa);
                }
                if (safmin * bcoefa > bscale) {
                    scale = min(scale, bscale / (safmin * bcoefa));
                }
                if (scale != one) {
                    acoef = scale * acoef;
                    acoefa = abs(acoef);
                    bcoefr = scale * bcoefr;
                    bcoefi = scale * bcoefi;
                    bcoefa = abs(bcoefr) + abs(bcoefi);
                }
                //
                //              Compute first two components of eigenvector
                //
                temp = acoef * s[((je + 1) - 1) + (je - 1) * lds];
                temp2r = acoef * s[(je - 1) + (je - 1) * lds] - bcoefr * p[(je - 1) + (je - 1) * ldp];
                temp2i = -bcoefi * p[(je - 1) + (je - 1) * ldp];
                if (abs(temp) > abs(temp2r) + abs(temp2i)) {
                    work[(2 * n + je) - 1] = one;
                    work[(3 * n + je) - 1] = zero;
                    work[(2 * n + je + 1) - 1] = -temp2r / temp;
                    work[(3 * n + je + 1) - 1] = -temp2i / temp;
                } else {
                    work[(2 * n + je + 1) - 1] = one;
                    work[(3 * n + je + 1) - 1] = zero;
                    temp = acoef * s[(je - 1) + ((je + 1) - 1) * lds];
                    work[(2 * n + je) - 1] = (bcoefr * p[((je + 1) - 1) + ((je + 1) - 1) * ldp] - acoef * s[((je + 1) - 1) + ((je + 1) - 1) * lds]) / temp;
                    work[(3 * n + je) - 1] = bcoefi * p[((je + 1) - 1) + ((je + 1) - 1) * ldp] / temp;
                }
                xmax = max(abs(work[(2 * n + je) - 1]) + abs(work[(3 * n + je) - 1]), abs(work[(2 * n + je + 1) - 1]) + abs(work[(3 * n + je + 1) - 1]));
            }
            //
            dmin = max(ulp * acoefa * anorm, ulp * bcoefa * bnorm, safmin);
            //
            //                                           T
            //           Triangular solve of  (a A - b B)  y = 0
            //
            //                                   T
            //           (rowwise in  (a A - b B) , or columnwise in (a A - b B) )
            //
            il2by2 = false;
            //
            for (j = je + nw; j <= n; j = j + 1) {
                if (il2by2) {
                    il2by2 = false;
                    goto statement_160;
                }
                //
                na = 1;
                bdiag[1 - 1] = p[(j - 1) + (j - 1) * ldp];
                if (j < n) {
                    if (s[((j + 1) - 1) + (j - 1) * lds] != zero) {
                        il2by2 = true;
                        bdiag[2 - 1] = p[((j + 1) - 1) + ((j + 1) - 1) * ldp];
                        na = 2;
                    }
                }
                //
                //              Check whether scaling is necessary for dot products
                //
                xscale = one / max(one, xmax);
                temp = max(work[j - 1], work[(n + j) - 1], acoefa * work[j - 1] + bcoefa * work[(n + j) - 1]);
                if (il2by2) {
                    temp = max(temp, work[(j + 1) - 1], work[(n + j + 1) - 1], acoefa * work[(j + 1) - 1] + bcoefa * work[(n + j + 1) - 1]);
                }
                if (temp > bignum * xscale) {
                    for (jw = 0; jw <= nw - 1; jw = jw + 1) {
                        for (jr = je; jr <= j - 1; jr = jr + 1) {
                            work[((jw + 2) * n + jr) - 1] = xscale * work[((jw + 2) * n + jr) - 1];
                        }
                    }
                    xmax = xmax * xscale;
                }
                //
                //              Compute dot products
                //
                //                    j-1
                //              SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k)
                //                    k=je
                //
                //              To reduce the op count, this is done as
                //
                //              _        j-1                  _        j-1
                //              a*conjg( sum  S(k,j)*x(k) ) - b*conjg( sum  P(k,j)*x(k) )
                //                       k=je                          k=je
                //
                //              which may cause underflow problems if A or B are close
                //              to underflow.  (E.g., less than SMALL.)
                //
                for (jw = 1; jw <= nw; jw = jw + 1) {
                    for (ja = 1; ja <= na; ja = ja + 1) {
                        sums[(ja - 1) + (jw - 1) * ldsums] = zero;
                        sump[(ja - 1) + (jw - 1) * ldsump] = zero;
                        //
                        for (jr = je; jr <= j - 1; jr = jr + 1) {
                            sums[(ja - 1) + (jw - 1) * ldsums] += s[(jr - 1) + ((j + ja - 1) - 1) * lds] * work[((jw + 1) * n + jr) - 1];
                            sump[(ja - 1) + (jw - 1) * ldsump] += p[(jr - 1) + ((j + ja - 1) - 1) * ldp] * work[((jw + 1) * n + jr) - 1];
                        }
                    }
                }
                //
                for (ja = 1; ja <= na; ja = ja + 1) {
                    if (ilcplx) {
                        sum[(ja - 1)] = -acoef * sums[(ja - 1)] + bcoefr * sump[(ja - 1)] - bcoefi * sump[(ja - 1) + (2 - 1) * ldsump];
                        sum[(ja - 1) + (2 - 1) * ldsum] = -acoef * sums[(ja - 1) + (2 - 1) * ldsums] + bcoefr * sump[(ja - 1) + (2 - 1) * ldsump] + bcoefi * sump[(ja - 1)];
                    } else {
                        sum[(ja - 1)] = -acoef * sums[(ja - 1)] + bcoefr * sump[(ja - 1)];
                    }
                }
                //
                //                                  T
                //              Solve  ( a A - b B )  y = SUM(,)
                //              with scaling and perturbation of the denominator
                //
                Rlaln2(true, na, nw, dmin, acoef, s[(j - 1) + (j - 1) * lds], lds, bdiag[1 - 1], bdiag[2 - 1], sum, 2, bcoefr, bcoefi, work[(2 * n + j) - 1], n, scale, temp, iinfo);
                if (scale < one) {
                    for (jw = 0; jw <= nw - 1; jw = jw + 1) {
                        for (jr = je; jr <= j - 1; jr = jr + 1) {
                            work[((jw + 2) * n + jr) - 1] = scale * work[((jw + 2) * n + jr) - 1];
                        }
                    }
                    xmax = scale * xmax;
                }
                xmax = max(xmax, temp);
            statement_160:;
            }
            //
            //           Copy eigenvector to VL, back transforming if
            //           HOWMNY='B'.
            //
            ieig++;
            if (ilback) {
                for (jw = 0; jw <= nw - 1; jw = jw + 1) {
                    Rgemv("N", n, n + 1 - je, one, vl[(je - 1) * ldvl], ldvl, work[((jw + 2) * n + je) - 1], 1, zero, work[((jw + 4) * n + 1) - 1], 1);
                }
                Rlacpy(" ", n, nw, work[(4 * n + 1) - 1], n, vl[(je - 1) * ldvl], ldvl);
                ibeg = 1;
            } else {
                Rlacpy(" ", n, nw, work[(2 * n + 1) - 1], n, vl[(ieig - 1) * ldvl], ldvl);
                ibeg = je;
            }
            //
            //           Scale eigenvector
            //
            xmax = zero;
            if (ilcplx) {
                for (j = ibeg; j <= n; j = j + 1) {
                    xmax = max(xmax, abs(vl[(j - 1) + (ieig - 1) * ldvl]) + abs(vl[(j - 1) + ((ieig + 1) - 1) * ldvl]));
                }
            } else {
                for (j = ibeg; j <= n; j = j + 1) {
                    xmax = max(xmax, abs(vl[(j - 1) + (ieig - 1) * ldvl]));
                }
            }
            //
            if (xmax > safmin) {
                xscale = one / xmax;
                //
                for (jw = 0; jw <= nw - 1; jw = jw + 1) {
                    for (jr = ibeg; jr <= n; jr = jr + 1) {
                        vl[(jr - 1) + ((ieig + jw) - 1) * ldvl] = xscale * vl[(jr - 1) + ((ieig + jw) - 1) * ldvl];
                    }
                }
            }
            ieig += nw - 1;
        //
        statement_220:;
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
        ilcplx = false;
        for (je = n; je >= 1; je = je - 1) {
            //
            //           Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or
            //           (b) this would be the second of a complex pair.
            //           Check for complex eigenvalue, so as to be sure of which
            //           entry(-ies) of SELECT to look at -- if complex, SELECT(JE)
            //           or SELECT(JE-1).
            //           If this is a complex pair, the 2-by-2 diagonal block
            //           corresponding to the eigenvalue is in rows/columns JE-1:JE
            //
            if (ilcplx) {
                ilcplx = false;
                goto statement_500;
            }
            nw = 1;
            if (je > 1) {
                if (s[(je - 1) + ((je - 1) - 1) * lds] != zero) {
                    ilcplx = true;
                    nw = 2;
                }
            }
            if (ilall) {
                ilcomp = true;
            } else if (ilcplx) {
                ilcomp = select[je - 1] || select[(je - 1) - 1];
            } else {
                ilcomp = select[je - 1];
            }
            if (!ilcomp) {
                goto statement_500;
            }
            //
            //           Decide if (a) singular pencil, (b) real eigenvalue, or
            //           (c) complex eigenvalue.
            //
            if (!ilcplx) {
                if (abs(s[(je - 1) + (je - 1) * lds]) <= safmin && abs(p[(je - 1) + (je - 1) * ldp]) <= safmin) {
                    //
                    //                 Singular matrix pencil -- unit eigenvector
                    //
                    ieig = ieig - 1;
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        vr[(jr - 1) + (ieig - 1) * ldvr] = zero;
                    }
                    vr[(ieig - 1) + (ieig - 1) * ldvr] = one;
                    goto statement_500;
                }
            }
            //
            //           Clear vector
            //
            for (jw = 0; jw <= nw - 1; jw = jw + 1) {
                for (jr = 1; jr <= n; jr = jr + 1) {
                    work[((jw + 2) * n + jr) - 1] = zero;
                }
            }
            //
            //           Compute coefficients in  ( a A - b B ) x = 0
            //              a  is  ACOEF
            //              b  is  BCOEFR + i*BCOEFI
            //
            if (!ilcplx) {
                //
                //              Real eigenvalue
                //
                temp = one / max(abs(s[(je - 1) + (je - 1) * lds]) * ascale, abs(p[(je - 1) + (je - 1) * ldp]) * bscale, safmin);
                salfar = (temp * s[(je - 1) + (je - 1) * lds]) * ascale;
                sbeta = (temp * p[(je - 1) + (je - 1) * ldp]) * bscale;
                acoef = sbeta * ascale;
                bcoefr = salfar * bscale;
                bcoefi = zero;
                //
                //              Scale to avoid underflow
                //
                scale = one;
                lsa = abs(sbeta) >= safmin && abs(acoef) < small;
                lsb = abs(salfar) >= safmin && abs(bcoefr) < small;
                if (lsa) {
                    scale = (small / abs(sbeta)) * min(anorm, big);
                }
                if (lsb) {
                    scale = max(scale, (small / abs(salfar)) * min(bnorm, big));
                }
                if (lsa || lsb) {
                    scale = min(scale, one / (safmin * max(one, abs(acoef), abs(bcoefr))));
                    if (lsa) {
                        acoef = ascale * (scale * sbeta);
                    } else {
                        acoef = scale * acoef;
                    }
                    if (lsb) {
                        bcoefr = bscale * (scale * salfar);
                    } else {
                        bcoefr = scale * bcoefr;
                    }
                }
                acoefa = abs(acoef);
                bcoefa = abs(bcoefr);
                //
                //              First component is 1
                //
                work[(2 * n + je) - 1] = one;
                xmax = one;
                //
                //              Compute contribution from column JE of A and B to sum
                //              (See "Further Details", above.)
                //
                for (jr = 1; jr <= je - 1; jr = jr + 1) {
                    work[(2 * n + jr) - 1] = bcoefr * p[(jr - 1) + (je - 1) * ldp] - acoef * s[(jr - 1) + (je - 1) * lds];
                }
            } else {
                //
                //              Complex eigenvalue
                //
                Rlag2(s[((je - 1) - 1) + ((je - 1) - 1) * lds], lds, p[((je - 1) - 1) + ((je - 1) - 1) * ldp], ldp, safmin * safety, acoef, temp, bcoefr, temp2, bcoefi);
                if (bcoefi == zero) {
                    info = je - 1;
                    return;
                }
                //
                //              Scale to avoid over/underflow
                //
                acoefa = abs(acoef);
                bcoefa = abs(bcoefr) + abs(bcoefi);
                scale = one;
                if (acoefa * ulp < safmin && acoefa >= safmin) {
                    scale = (safmin / ulp) / acoefa;
                }
                if (bcoefa * ulp < safmin && bcoefa >= safmin) {
                    scale = max(scale, (safmin / ulp) / bcoefa);
                }
                if (safmin * acoefa > ascale) {
                    scale = ascale / (safmin * acoefa);
                }
                if (safmin * bcoefa > bscale) {
                    scale = min(scale, bscale / (safmin * bcoefa));
                }
                if (scale != one) {
                    acoef = scale * acoef;
                    acoefa = abs(acoef);
                    bcoefr = scale * bcoefr;
                    bcoefi = scale * bcoefi;
                    bcoefa = abs(bcoefr) + abs(bcoefi);
                }
                //
                //              Compute first two components of eigenvector
                //              and contribution to sums
                //
                temp = acoef * s[(je - 1) + ((je - 1) - 1) * lds];
                temp2r = acoef * s[(je - 1) + (je - 1) * lds] - bcoefr * p[(je - 1) + (je - 1) * ldp];
                temp2i = -bcoefi * p[(je - 1) + (je - 1) * ldp];
                if (abs(temp) >= abs(temp2r) + abs(temp2i)) {
                    work[(2 * n + je) - 1] = one;
                    work[(3 * n + je) - 1] = zero;
                    work[(2 * n + je - 1) - 1] = -temp2r / temp;
                    work[(3 * n + je - 1) - 1] = -temp2i / temp;
                } else {
                    work[(2 * n + je - 1) - 1] = one;
                    work[(3 * n + je - 1) - 1] = zero;
                    temp = acoef * s[((je - 1) - 1) + (je - 1) * lds];
                    work[(2 * n + je) - 1] = (bcoefr * p[((je - 1) - 1) + ((je - 1) - 1) * ldp] - acoef * s[((je - 1) - 1) + ((je - 1) - 1) * lds]) / temp;
                    work[(3 * n + je) - 1] = bcoefi * p[((je - 1) - 1) + ((je - 1) - 1) * ldp] / temp;
                }
                //
                xmax = max(abs(work[(2 * n + je) - 1]) + abs(work[(3 * n + je) - 1]), abs(work[(2 * n + je - 1) - 1]) + abs(work[(3 * n + je - 1) - 1]));
                //
                //              Compute contribution from columns JE and JE-1
                //              of A and B to the sums.
                //
                creala = acoef * work[(2 * n + je - 1) - 1];
                cimaga = acoef * work[(3 * n + je - 1) - 1];
                crealb = bcoefr * work[(2 * n + je - 1) - 1] - bcoefi * work[(3 * n + je - 1) - 1];
                cimagb = bcoefi * work[(2 * n + je - 1) - 1] + bcoefr * work[(3 * n + je - 1) - 1];
                cre2a = acoef * work[(2 * n + je) - 1];
                cim2a = acoef * work[(3 * n + je) - 1];
                cre2b = bcoefr * work[(2 * n + je) - 1] - bcoefi * work[(3 * n + je) - 1];
                cim2b = bcoefi * work[(2 * n + je) - 1] + bcoefr * work[(3 * n + je) - 1];
                for (jr = 1; jr <= je - 2; jr = jr + 1) {
                    work[(2 * n + jr) - 1] = -creala * s[(jr - 1) + ((je - 1) - 1) * lds] + crealb * p[(jr - 1) + ((je - 1) - 1) * ldp] - cre2a * s[(jr - 1) + (je - 1) * lds] + cre2b * p[(jr - 1) + (je - 1) * ldp];
                    work[(3 * n + jr) - 1] = -cimaga * s[(jr - 1) + ((je - 1) - 1) * lds] + cimagb * p[(jr - 1) + ((je - 1) - 1) * ldp] - cim2a * s[(jr - 1) + (je - 1) * lds] + cim2b * p[(jr - 1) + (je - 1) * ldp];
                }
            }
            //
            dmin = max(ulp * acoefa * anorm, ulp * bcoefa * bnorm, safmin);
            //
            //           Columnwise triangular solve of  (a A - b B)  x = 0
            //
            il2by2 = false;
            for (j = je - nw; j >= 1; j = j - 1) {
                //
                //              If a 2-by-2 block, is in position j-1:j, wait until
                //              next iteration to process it (when it will be j:j+1)
                //
                if (!il2by2 && j > 1) {
                    if (s[(j - 1) + ((j - 1) - 1) * lds] != zero) {
                        il2by2 = true;
                        goto statement_370;
                    }
                }
                bdiag[1 - 1] = p[(j - 1) + (j - 1) * ldp];
                if (il2by2) {
                    na = 2;
                    bdiag[2 - 1] = p[((j + 1) - 1) + ((j + 1) - 1) * ldp];
                } else {
                    na = 1;
                }
                //
                //              Compute x(j) (and x(j+1), if 2-by-2 block)
                //
                Rlaln2(false, na, nw, dmin, acoef, s[(j - 1) + (j - 1) * lds], lds, bdiag[1 - 1], bdiag[2 - 1], work[(2 * n + j) - 1], n, bcoefr, bcoefi, sum, 2, scale, temp, iinfo);
                if (scale < one) {
                    //
                    for (jw = 0; jw <= nw - 1; jw = jw + 1) {
                        for (jr = 1; jr <= je; jr = jr + 1) {
                            work[((jw + 2) * n + jr) - 1] = scale * work[((jw + 2) * n + jr) - 1];
                        }
                    }
                }
                xmax = max(scale * xmax, temp);
                //
                for (jw = 1; jw <= nw; jw = jw + 1) {
                    for (ja = 1; ja <= na; ja = ja + 1) {
                        work[((jw + 1) * n + j + ja - 1) - 1] = sum[(ja - 1) + (jw - 1) * ldsum];
                    }
                }
                //
                //              w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling
                //
                if (j > 1) {
                    //
                    //                 Check whether scaling is necessary for sum.
                    //
                    xscale = one / max(one, xmax);
                    temp = acoefa * work[j - 1] + bcoefa * work[(n + j) - 1];
                    if (il2by2) {
                        temp = max(temp, acoefa * work[(j + 1) - 1] + bcoefa * work[(n + j + 1) - 1]);
                    }
                    temp = max(temp, acoefa, bcoefa);
                    if (temp > bignum * xscale) {
                        //
                        for (jw = 0; jw <= nw - 1; jw = jw + 1) {
                            for (jr = 1; jr <= je; jr = jr + 1) {
                                work[((jw + 2) * n + jr) - 1] = xscale * work[((jw + 2) * n + jr) - 1];
                            }
                        }
                        xmax = xmax * xscale;
                    }
                    //
                    //                 Compute the contributions of the off-diagonals of
                    //                 column j (and j+1, if 2-by-2 block) of A and B to the
                    //                 sums.
                    //
                    for (ja = 1; ja <= na; ja = ja + 1) {
                        if (ilcplx) {
                            creala = acoef * work[(2 * n + j + ja - 1) - 1];
                            cimaga = acoef * work[(3 * n + j + ja - 1) - 1];
                            crealb = bcoefr * work[(2 * n + j + ja - 1) - 1] - bcoefi * work[(3 * n + j + ja - 1) - 1];
                            cimagb = bcoefi * work[(2 * n + j + ja - 1) - 1] + bcoefr * work[(3 * n + j + ja - 1) - 1];
                            for (jr = 1; jr <= j - 1; jr = jr + 1) {
                                work[(2 * n + jr) - 1] = work[(2 * n + jr) - 1] - creala * s[(jr - 1) + ((j + ja - 1) - 1) * lds] + crealb * p[(jr - 1) + ((j + ja - 1) - 1) * ldp];
                                work[(3 * n + jr) - 1] = work[(3 * n + jr) - 1] - cimaga * s[(jr - 1) + ((j + ja - 1) - 1) * lds] + cimagb * p[(jr - 1) + ((j + ja - 1) - 1) * ldp];
                            }
                        } else {
                            creala = acoef * work[(2 * n + j + ja - 1) - 1];
                            crealb = bcoefr * work[(2 * n + j + ja - 1) - 1];
                            for (jr = 1; jr <= j - 1; jr = jr + 1) {
                                work[(2 * n + jr) - 1] = work[(2 * n + jr) - 1] - creala * s[(jr - 1) + ((j + ja - 1) - 1) * lds] + crealb * p[(jr - 1) + ((j + ja - 1) - 1) * ldp];
                            }
                        }
                    }
                }
                //
                il2by2 = false;
            statement_370:;
            }
            //
            //           Copy eigenvector to VR, back transforming if
            //           HOWMNY='B'.
            //
            ieig = ieig - nw;
            if (ilback) {
                //
                for (jw = 0; jw <= nw - 1; jw = jw + 1) {
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        work[((jw + 4) * n + jr) - 1] = work[((jw + 2) * n + 1) - 1] * vr[(jr - 1)];
                    }
                    //
                    //                 A series of compiler directives to defeat
                    //                 vectorization for the next loop
                    //
                    for (jc = 2; jc <= je; jc = jc + 1) {
                        for (jr = 1; jr <= n; jr = jr + 1) {
                            work[((jw + 4) * n + jr) - 1] += work[((jw + 2) * n + jc) - 1] * vr[(jr - 1) + (jc - 1) * ldvr];
                        }
                    }
                }
                //
                for (jw = 0; jw <= nw - 1; jw = jw + 1) {
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        vr[(jr - 1) + ((ieig + jw) - 1) * ldvr] = work[((jw + 4) * n + jr) - 1];
                    }
                }
                //
                iend = n;
            } else {
                for (jw = 0; jw <= nw - 1; jw = jw + 1) {
                    for (jr = 1; jr <= n; jr = jr + 1) {
                        vr[(jr - 1) + ((ieig + jw) - 1) * ldvr] = work[((jw + 2) * n + jr) - 1];
                    }
                }
                //
                iend = je;
            }
            //
            //           Scale eigenvector
            //
            xmax = zero;
            if (ilcplx) {
                for (j = 1; j <= iend; j = j + 1) {
                    xmax = max(xmax, abs(vr[(j - 1) + (ieig - 1) * ldvr]) + abs(vr[(j - 1) + ((ieig + 1) - 1) * ldvr]));
                }
            } else {
                for (j = 1; j <= iend; j = j + 1) {
                    xmax = max(xmax, abs(vr[(j - 1) + (ieig - 1) * ldvr]));
                }
            }
            //
            if (xmax > safmin) {
                xscale = one / xmax;
                for (jw = 0; jw <= nw - 1; jw = jw + 1) {
                    for (jr = 1; jr <= iend; jr = jr + 1) {
                        vr[(jr - 1) + ((ieig + jw) - 1) * ldvr] = xscale * vr[(jr - 1) + ((ieig + jw) - 1) * ldvr];
                    }
                }
            }
        statement_500:;
        }
    }
    //
    //     End of Rtgevc
    //
}
