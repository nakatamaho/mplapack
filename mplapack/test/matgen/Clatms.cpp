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

#include <mplapack_matgen.h>

void Clatms(INTEGER const m, INTEGER const n, const char *dist, INTEGER *iseed, const char *sym, REAL *d, INTEGER const mode, REAL const cond, REAL const dmax, INTEGER const kl, INTEGER const ku, const char *pack, COMPLEX *a, INTEGER const lda, COMPLEX *work, INTEGER &info) {
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
    //     1)      Decode and Test the input parameters.
    //             Initialize flags & seed.
    //
    info = 0;
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0) {
        return;
    }
    //
    //     Decode DIST
    //
    INTEGER idist = 0;
    if (Mlsame(dist, "U")) {
        idist = 1;
    } else if (Mlsame(dist, "S")) {
        idist = 2;
    } else if (Mlsame(dist, "N")) {
        idist = 3;
    } else {
        idist = -1;
    }
    //
    //     Decode SYM
    //
    INTEGER isym = 0;
    INTEGER irsign = 0;
    bool zsym = false;
    if (Mlsame(sym, "N")) {
        isym = 1;
        irsign = 0;
        zsym = false;
    } else if (Mlsame(sym, "P")) {
        isym = 2;
        irsign = 0;
        zsym = false;
    } else if (Mlsame(sym, "S")) {
        isym = 2;
        irsign = 0;
        zsym = true;
    } else if (Mlsame(sym, "H")) {
        isym = 2;
        irsign = 1;
        zsym = false;
    } else {
        isym = -1;
    }
    //
    //     Decode PACK
    //
    INTEGER isympk = 0;
    INTEGER ipack = 0;
    if (Mlsame(pack, "N")) {
        ipack = 0;
    } else if (Mlsame(pack, "U")) {
        ipack = 1;
        isympk = 1;
    } else if (Mlsame(pack, "L")) {
        ipack = 2;
        isympk = 1;
    } else if (Mlsame(pack, "C")) {
        ipack = 3;
        isympk = 2;
    } else if (Mlsame(pack, "R")) {
        ipack = 4;
        isympk = 3;
    } else if (Mlsame(pack, "B")) {
        ipack = 5;
        isympk = 3;
    } else if (Mlsame(pack, "Q")) {
        ipack = 6;
        isympk = 2;
    } else if (Mlsame(pack, "Z")) {
        ipack = 7;
    } else {
        ipack = -1;
    }
    //
    //     Set certain internal parameters
    //
    INTEGER mnmin = min(m, n);
    INTEGER llb = min(kl, m - 1);
    INTEGER uub = min(ku, n - 1);
    INTEGER mr = min(m, n + llb);
    INTEGER nc = min(n, m + uub);
    //
    INTEGER minlda = 0;
    if (ipack == 5 || ipack == 6) {
        minlda = uub + 1;
    } else if (ipack == 7) {
        minlda = llb + uub + 1;
    } else {
        minlda = m;
    }
    //
    //     Use Givens rotation method if bandwidth small enough,
    //     or if LDA is too small to store the matrix unpacked.
    //
    bool givens = false;
    if (isym == 1) {
        if (castREAL(llb + uub) < 0.3 * castREAL(max((INTEGER)1, mr + nc))) {
            givens = true;
        }
    } else {
        if ((INTEGER)2 * llb < m) {
            givens = true;
        }
    }
    if (lda < m && lda >= minlda) {
        givens = true;
    }
    //
    //     Set INFO if an error
    //
    const REAL one = 1.0;
    if (m < 0) {
        info = -1;
    } else if (m != n && isym != 1) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (idist == -1) {
        info = -3;
    } else if (isym == -1) {
        info = -5;
    } else if (abs(mode) > 6) {
        info = -7;
    } else if ((mode != 0 && abs(mode) != 6) && cond < one) {
        info = -8;
    } else if (kl < 0) {
        info = -10;
    } else if (ku < 0 || (isym != 1 && kl != ku)) {
        info = -11;
    } else if (ipack == -1 || (isympk == 1 && isym == 1) || (isympk == 2 && isym == 1 && kl > 0) || (isympk == 3 && isym == 1 && ku > 0) || (isympk != 0 && m != n)) {
        info = -12;
    } else if (lda < max((INTEGER)1, minlda)) {
        info = -14;
    }
    //
    if (info != 0) {
        Mxerbla("Clatms", -info);
        return;
    }
    //
    //     Initialize random number generator
    //
    INTEGER i = 0;
    for (i = 1; i <= 4; i = i + 1) {
        iseed[i - 1] = mod(abs(iseed[i - 1]), 4096);
    }
    //
    //     2)      Set up D  if indicated.
    //
    //             Compute D according to COND and MODE
    //
    INTEGER iinfo = 0;
    Rlatm1(mode, cond, irsign, idist, iseed, d, mnmin, iinfo);
    if (iinfo != 0) {
        info = 1;
        return;
    }
    //
    //     Choose Top-Down if D is (apparently) increasing,
    //     Bottom-Up if D is (apparently) decreasing.
    //
    bool topdwn = false;
    if (abs(d[1 - 1]) <= abs(d[mnmin - 1])) {
        topdwn = true;
    } else {
        topdwn = false;
    }
    //
    REAL temp = 0.0;
    const REAL zero = 0.0;
    REAL alpha = 0.0;
    if (mode != 0 && abs(mode) != 6) {
        //
        //        Scale by DMAX
        //
        temp = abs(d[1 - 1]);
        for (i = 2; i <= mnmin; i = i + 1) {
            temp = max(temp, abs(d[i - 1]));
        }
        //
        if (temp > zero) {
            alpha = dmax / temp;
        } else {
            info = 2;
            return;
        }
        //
        Rscal(mnmin, alpha, d, 1);
        //
    }
    //
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    Claset("Full", lda, n, czero, czero, a, lda);
    //
    //     3)      Generate Banded Matrix using Givens rotations.
    //             Also the special case of UUB=LLB=0
    //
    //               Compute Addressing constants to cover all
    //               storage formats.  Whether GE, HE, SY, GB, HB, or SB,
    //               upper or lower triangle or both,
    //               the (i,j)-th element is in
    //               A( i - ISKEW*j + IOFFST, j )
    //
    INTEGER ilda = 0;
    INTEGER iskew = 0;
    INTEGER ioffst = 0;
    if (ipack > 4) {
        ilda = lda - 1;
        iskew = 1;
        if (ipack > 5) {
            ioffst = uub + 1;
        } else {
            ioffst = 1;
        }
    } else {
        ilda = lda;
        iskew = 0;
        ioffst = 0;
    }
    //
    //     IPACKG is the format that the matrix is generated in. If this is
    //     different from IPACK, then the matrix must be repacked at the
    //     end.  It also signals how to compute the norm, for scaling.
    //
    INTEGER ipackg = 0;
    //
    //     Diagonal Matrix -- We are done, unless it
    //     is to be stored HP/SP/PP/TP (PACK='R' or 'C')
    //
    INTEGER j = 0;
    INTEGER jkl = 0;
    INTEGER jku = 0;
    INTEGER jr = 0;
    COMPLEX extra = 0.0;
    const REAL twopi = 6.28318530717958647692528676655900576839e+0;
    REAL angle = 0.0;
    COMPLEX c = 0.0;
    COMPLEX s = 0.0;
    INTEGER icol = 0;
    INTEGER il = 0;
    COMPLEX dummy = 0.0;
    INTEGER ir = 0;
    INTEGER ic = 0;
    INTEGER jch = 0;
    REAL realc = 0.0;
    INTEGER irow = 0;
    COMPLEX ctemp = 0.0;
    bool iltemp = false;
    INTEGER jc = 0;
    INTEGER iendch = 0;
    bool ilextr = false;
    INTEGER ioffg = 0;
    INTEGER k = 0;
    COMPLEX ct = 0.0;
    COMPLEX st = 0.0;
    if (llb == 0 && uub == 0) {
        for (j = 1; j <= mnmin; j = j + 1) {
            a[(((1 - iskew) * j + ioffst) - 1) + (j - 1) * lda] = COMPLEX(d[j - 1]);
        }
        //
        if (ipack <= 2 || ipack >= 5) {
            ipackg = ipack;
        }
        //
    } else if (givens) {
        //
        //        Check whether to use Givens rotations,
        //        Householder transformations, or nothing.
        //
        if (isym == 1) {
            //
            //           Non-symmetric -- A = U D V
            //
            if (ipack > 4) {
                ipackg = ipack;
            } else {
                ipackg = 0;
            }
            //
            for (j = 1; j <= mnmin; j = j + 1) {
                a[(((1 - iskew) * j + ioffst) - 1) + (j - 1) * lda] = COMPLEX(d[j - 1]);
            }
            //
            if (topdwn) {
                jkl = 0;
                for (jku = 1; jku <= uub; jku = jku + 1) {
                    //
                    //                 Transform from bandwidth JKL, JKU-1 to JKL, JKU
                    //
                    //                 Last row actually rotated is M
                    //                 Last column actually rotated is MIN( M+JKU, N )
                    //
                    for (jr = 1; jr <= min(m + jku, n) + jkl - 1; jr = jr + 1) {
                        extra = czero;
                        angle = twopi * Rlarnd(1, iseed);
                        c = cos(angle) * Clarnd(5, iseed);
                        s = sin(angle) * Clarnd(5, iseed);
                        icol = max((INTEGER)1, jr - jkl);
                        if (jr < m) {
                            il = min(n, jr + jku) + 1 - icol;
                            Clarot(true, jr > jkl, false, il, c, s, &a[((jr - iskew * icol + ioffst) - 1) + (icol - 1) * lda], ilda, extra, dummy);
                        }
                        //
                        //                    Chase "EXTRA" back up
                        //
                        ir = jr;
                        ic = icol;
                        for (jch = jr - jkl; jch >= 1; jch = jch - jkl - jku) {
                            if (ir < m) {
                                Clartg(a[((ir + 1 - iskew * (ic + 1) + ioffst) - 1) + ((ic + 1) - 1) * lda], extra, realc, s, dummy);
                                dummy = Clarnd(5, iseed);
                                c = conj(realc * dummy);
                                s = conj(-s * dummy);
                            }
                            irow = max((INTEGER)1, jch - jku);
                            il = ir + 2 - irow;
                            ctemp = czero;
                            iltemp = jch > jku;
                            Clarot(false, iltemp, true, il, c, s, &a[((irow - iskew * ic + ioffst) - 1) + (ic - 1) * lda], ilda, ctemp, extra);
                            if (iltemp) {
                                Clartg(a[((irow + 1 - iskew * (ic + 1) + ioffst) - 1) + ((ic + 1) - 1) * lda], ctemp, realc, s, dummy);
                                dummy = Clarnd(5, iseed);
                                c = conj(realc * dummy);
                                s = conj(-s * dummy);
                                //
                                icol = max((INTEGER)1, jch - jku - jkl);
                                il = ic + 2 - icol;
                                extra = czero;
                                Clarot(true, jch > jku + jkl, true, il, c, s, &a[((irow - iskew * icol + ioffst) - 1) + (icol - 1) * lda], ilda, extra, ctemp);
                                ic = icol;
                                ir = irow;
                            }
                        }
                    }
                }
                //
                jku = uub;
                for (jkl = 1; jkl <= llb; jkl = jkl + 1) {
                    //
                    //                 Transform from bandwidth JKL-1, JKU to JKL, JKU
                    //
                    for (jc = 1; jc <= min(n + jkl, m) + jku - 1; jc = jc + 1) {
                        extra = czero;
                        angle = twopi * Rlarnd(1, iseed);
                        c = cos(angle) * Clarnd(5, iseed);
                        s = sin(angle) * Clarnd(5, iseed);
                        irow = max((INTEGER)1, jc - jku);
                        if (jc < n) {
                            il = min(m, jc + jkl) + 1 - irow;
                            Clarot(false, jc > jku, false, il, c, s, &a[((irow - iskew * jc + ioffst) - 1) + (jc - 1) * lda], ilda, extra, dummy);
                        }
                        //
                        //                    Chase "EXTRA" back up
                        //
                        ic = jc;
                        ir = irow;
                        for (jch = jc - jku; jch >= 1; jch = jch - jkl - jku) {
                            if (ic < n) {
                                Clartg(a[((ir + 1 - iskew * (ic + 1) + ioffst) - 1) + ((ic + 1) - 1) * lda], extra, realc, s, dummy);
                                dummy = Clarnd(5, iseed);
                                c = conj(realc * dummy);
                                s = conj(-s * dummy);
                            }
                            icol = max((INTEGER)1, jch - jkl);
                            il = ic + 2 - icol;
                            ctemp = czero;
                            iltemp = jch > jkl;
                            Clarot(true, iltemp, true, il, c, s, &a[((ir - iskew * icol + ioffst) - 1) + (icol - 1) * lda], ilda, ctemp, extra);
                            if (iltemp) {
                                Clartg(a[((ir + 1 - iskew * (icol + 1) + ioffst) - 1) + ((icol + 1) - 1) * lda], ctemp, realc, s, dummy);
                                dummy = Clarnd(5, iseed);
                                c = conj(realc * dummy);
                                s = conj(-s * dummy);
                                irow = max((INTEGER)1, jch - jkl - jku);
                                il = ir + 2 - irow;
                                extra = czero;
                                Clarot(false, jch > jkl + jku, true, il, c, s, &a[((irow - iskew * icol + ioffst) - 1) + (icol - 1) * lda], ilda, extra, ctemp);
                                ic = icol;
                                ir = irow;
                            }
                        }
                    }
                }
                //
            } else {
                //
                //              Bottom-Up -- Start at the bottom right.
                //
                jkl = 0;
                for (jku = 1; jku <= uub; jku = jku + 1) {
                    //
                    //                 Transform from bandwidth JKL, JKU-1 to JKL, JKU
                    //
                    //                 First row actually rotated is M
                    //                 First column actually rotated is MIN( M+JKU, N )
                    //
                    iendch = min(m, n + jkl) - 1;
                    for (jc = min(m + jku, n) - 1; jc >= 1 - jkl; jc = jc - 1) {
                        extra = czero;
                        angle = twopi * Rlarnd(1, iseed);
                        c = cos(angle) * Clarnd(5, iseed);
                        s = sin(angle) * Clarnd(5, iseed);
                        irow = max((INTEGER)1, jc - jku + 1);
                        if (jc > 0) {
                            il = min(m, jc + jkl + 1) + 1 - irow;
                            Clarot(false, false, jc + jkl < m, il, c, s, &a[((irow - iskew * jc + ioffst) - 1) + (jc - 1) * lda], ilda, dummy, extra);
                        }
                        //
                        //                    Chase "EXTRA" back down
                        //
                        ic = jc;
                        for (jch = jc + jkl; jch <= iendch; jch = jch + jkl + jku) {
                            ilextr = ic > 0;
                            if (ilextr) {
                                Clartg(a[((jch - iskew * ic + ioffst) - 1) + (ic - 1) * lda], extra, realc, s, dummy);
                                dummy = Clarnd(5, iseed);
                                c = realc * dummy;
                                s = s * dummy;
                            }
                            ic = max((INTEGER)1, ic);
                            icol = min(n - 1, jch + jku);
                            iltemp = jch + jku < n;
                            ctemp = czero;
                            Clarot(true, ilextr, iltemp, icol + 2 - ic, c, s, &a[((jch - iskew * ic + ioffst) - 1) + (ic - 1) * lda], ilda, extra, ctemp);
                            if (iltemp) {
                                Clartg(a[((jch - iskew * icol + ioffst) - 1) + (icol - 1) * lda], ctemp, realc, s, dummy);
                                dummy = Clarnd(5, iseed);
                                c = realc * dummy;
                                s = s * dummy;
                                il = min(iendch, jch + jkl + jku) + 2 - jch;
                                extra = czero;
                                Clarot(false, true, jch + jkl + jku <= iendch, il, c, s, &a[((jch - iskew * icol + ioffst) - 1) + (icol - 1) * lda], ilda, ctemp, extra);
                                ic = icol;
                            }
                        }
                    }
                }
                //
                jku = uub;
                for (jkl = 1; jkl <= llb; jkl = jkl + 1) {
                    //
                    //                 Transform from bandwidth JKL-1, JKU to JKL, JKU
                    //
                    //                 First row actually rotated is MIN( N+JKL, M )
                    //                 First column actually rotated is N
                    //
                    iendch = min(n, m + jku) - 1;
                    for (jr = min(n + jkl, m) - 1; jr >= 1 - jku; jr = jr - 1) {
                        extra = czero;
                        angle = twopi * Rlarnd(1, iseed);
                        c = cos(angle) * Clarnd(5, iseed);
                        s = sin(angle) * Clarnd(5, iseed);
                        icol = max((INTEGER)1, jr - jkl + 1);
                        if (jr > 0) {
                            il = min(n, jr + jku + 1) + 1 - icol;
                            Clarot(true, false, jr + jku < n, il, c, s, &a[((jr - iskew * icol + ioffst) - 1) + (icol - 1) * lda], ilda, dummy, extra);
                        }
                        //
                        //                    Chase "EXTRA" back down
                        //
                        ir = jr;
                        for (jch = jr + jku; jch <= iendch; jch = jch + jkl + jku) {
                            ilextr = ir > 0;
                            if (ilextr) {
                                Clartg(a[((ir - iskew * jch + ioffst) - 1) + (jch - 1) * lda], extra, realc, s, dummy);
                                dummy = Clarnd(5, iseed);
                                c = realc * dummy;
                                s = s * dummy;
                            }
                            ir = max((INTEGER)1, ir);
                            irow = min(m - 1, jch + jkl);
                            iltemp = jch + jkl < m;
                            ctemp = czero;
                            Clarot(false, ilextr, iltemp, irow + 2 - ir, c, s, &a[((ir - iskew * jch + ioffst) - 1) + (jch - 1) * lda], ilda, extra, ctemp);
                            if (iltemp) {
                                Clartg(a[((irow - iskew * jch + ioffst) - 1) + (jch - 1) * lda], ctemp, realc, s, dummy);
                                dummy = Clarnd(5, iseed);
                                c = realc * dummy;
                                s = s * dummy;
                                il = min(iendch, jch + jkl + jku) + 2 - jch;
                                extra = czero;
                                Clarot(true, true, jch + jkl + jku <= iendch, il, c, s, &a[((irow - iskew * jch + ioffst) - 1) + (jch - 1) * lda], ilda, ctemp, extra);
                                ir = irow;
                            }
                        }
                    }
                }
                //
            }
            //
        } else {
            //
            //           Symmetric -- A = U D U'
            //           Hermitian -- A = U D U*
            //
            ipackg = ipack;
            ioffg = ioffst;
            //
            if (topdwn) {
                //
                //              Top-Down -- Generate Upper triangle only
                //
                if (ipack >= 5) {
                    ipackg = 6;
                    ioffg = uub + 1;
                } else {
                    ipackg = 1;
                }
                //
                for (j = 1; j <= mnmin; j = j + 1) {
                    a[(((1 - iskew) * j + ioffg) - 1) + (j - 1) * lda] = COMPLEX(d[j - 1]);
                }
                //
                for (k = 1; k <= uub; k = k + 1) {
                    for (jc = 1; jc <= n - 1; jc = jc + 1) {
                        irow = max((INTEGER)1, jc - k);
                        il = min(jc + 1, k + 2);
                        extra = czero;
                        ctemp = a[((jc - iskew * (jc + 1) + ioffg) - 1) + ((jc + 1) - 1) * lda];
                        angle = twopi * Rlarnd(1, iseed);
                        c = cos(angle) * Clarnd(5, iseed);
                        s = sin(angle) * Clarnd(5, iseed);
                        if (zsym) {
                            ct = c;
                            st = s;
                        } else {
                            ctemp = conj(ctemp);
                            ct = conj(c);
                            st = conj(s);
                        }
                        Clarot(false, jc > k, true, il, c, s, &a[((irow - iskew * jc + ioffg) - 1) + (jc - 1) * lda], ilda, extra, ctemp);
                        Clarot(true, true, false, min(k, n - jc) + 1, ct, st, &a[(((1 - iskew) * jc + ioffg) - 1) + (jc - 1) * lda], ilda, ctemp, dummy);
                        //
                        //                    Chase EXTRA back up the matrix
                        //
                        icol = jc;
                        for (jch = jc - k; jch >= 1; jch = jch - k) {
                            Clartg(a[((jch + 1 - iskew * (icol + 1) + ioffg) - 1) + ((icol + 1) - 1) * lda], extra, realc, s, dummy);
                            dummy = Clarnd(5, iseed);
                            c = conj(realc * dummy);
                            s = conj(-s * dummy);
                            ctemp = a[((jch - iskew * (jch + 1) + ioffg) - 1) + ((jch + 1) - 1) * lda];
                            if (zsym) {
                                ct = c;
                                st = s;
                            } else {
                                ctemp = conj(ctemp);
                                ct = conj(c);
                                st = conj(s);
                            }
                            Clarot(true, true, true, k + 2, c, s, &a[(((1 - iskew) * jch + ioffg) - 1) + (jch - 1) * lda], ilda, ctemp, extra);
                            irow = max((INTEGER)1, jch - k);
                            il = min(jch + 1, k + 2);
                            extra = czero;
                            Clarot(false, jch > k, true, il, ct, st, &a[((irow - iskew * jch + ioffg) - 1) + (jch - 1) * lda], ilda, extra, ctemp);
                            icol = jch;
                        }
                    }
                }
                //
                //              If we need lower triangle, copy from upper. Note that
                //              the order of copying is chosen to work for 'q' -> 'b'
                //
                if (ipack != ipackg && ipack != 3) {
                    for (jc = 1; jc <= n; jc = jc + 1) {
                        irow = ioffst - iskew * jc;
                        if (zsym) {
                            for (jr = jc; jr <= min(n, jc + uub); jr = jr + 1) {
                                a[((jr + irow) - 1) + (jc - 1) * lda] = a[((jc - iskew * jr + ioffg) - 1) + (jr - 1) * lda];
                            }
                        } else {
                            for (jr = jc; jr <= min(n, jc + uub); jr = jr + 1) {
                                a[((jr + irow) - 1) + (jc - 1) * lda] = conj(a[((jc - iskew * jr + ioffg) - 1) + (jr - 1) * lda]);
                            }
                        }
                    }
                    if (ipack == 5) {
                        for (jc = n - uub + 1; jc <= n; jc = jc + 1) {
                            for (jr = n + 2 - jc; jr <= uub + 1; jr = jr + 1) {
                                a[(jr - 1) + (jc - 1) * lda] = czero;
                            }
                        }
                    }
                    if (ipackg == 6) {
                        ipackg = ipack;
                    } else {
                        ipackg = 0;
                    }
                }
            } else {
                //
                //              Bottom-Up -- Generate Lower triangle only
                //
                if (ipack >= 5) {
                    ipackg = 5;
                    if (ipack == 6) {
                        ioffg = 1;
                    }
                } else {
                    ipackg = 2;
                }
                //
                for (j = 1; j <= mnmin; j = j + 1) {
                    a[(((1 - iskew) * j + ioffg) - 1) + (j - 1) * lda] = COMPLEX(d[j - 1]);
                }
                //
                for (k = 1; k <= uub; k = k + 1) {
                    for (jc = n - 1; jc >= 1; jc = jc - 1) {
                        il = min(n + 1 - jc, k + 2);
                        extra = czero;
                        ctemp = a[((1 + (1 - iskew) * jc + ioffg) - 1) + (jc - 1) * lda];
                        angle = twopi * Rlarnd(1, iseed);
                        c = cos(angle) * Clarnd(5, iseed);
                        s = sin(angle) * Clarnd(5, iseed);
                        if (zsym) {
                            ct = c;
                            st = s;
                        } else {
                            ctemp = conj(ctemp);
                            ct = conj(c);
                            st = conj(s);
                        }
                        Clarot(false, true, n - jc > k, il, c, s, &a[(((1 - iskew) * jc + ioffg) - 1) + (jc - 1) * lda], ilda, ctemp, extra);
                        icol = max((INTEGER)1, jc - k + 1);
                        Clarot(true, false, true, jc + 2 - icol, ct, st, &a[((jc - iskew * icol + ioffg) - 1) + (icol - 1) * lda], ilda, dummy, ctemp);
                        //
                        //                    Chase EXTRA back down the matrix
                        //
                        icol = jc;
                        for (jch = jc + k; jch <= n - 1; jch = jch + k) {
                            Clartg(a[((jch - iskew * icol + ioffg) - 1) + (icol - 1) * lda], extra, realc, s, dummy);
                            dummy = Clarnd(5, iseed);
                            c = realc * dummy;
                            s = s * dummy;
                            ctemp = a[((1 + (1 - iskew) * jch + ioffg) - 1) + (jch - 1) * lda];
                            if (zsym) {
                                ct = c;
                                st = s;
                            } else {
                                ctemp = conj(ctemp);
                                ct = conj(c);
                                st = conj(s);
                            }
                            Clarot(true, true, true, k + 2, c, s, &a[((jch - iskew * icol + ioffg) - 1) + (icol - 1) * lda], ilda, extra, ctemp);
                            il = min(n + 1 - jch, k + 2);
                            extra = czero;
                            Clarot(false, true, n - jch > k, il, ct, st, &a[(((1 - iskew) * jch + ioffg) - 1) + (jch - 1) * lda], ilda, ctemp, extra);
                            icol = jch;
                        }
                    }
                }
                //
                //              If we need upper triangle, copy from lower. Note that
                //              the order of copying is chosen to work for 'b' -> 'q'
                //
                if (ipack != ipackg && ipack != 4) {
                    for (jc = n; jc >= 1; jc = jc - 1) {
                        irow = ioffst - iskew * jc;
                        if (zsym) {
                            for (jr = jc; jr >= max((INTEGER)1, jc - uub); jr = jr - 1) {
                                a[((jr + irow) - 1) + (jc - 1) * lda] = a[((jc - iskew * jr + ioffg) - 1) + (jr - 1) * lda];
                            }
                        } else {
                            for (jr = jc; jr >= max((INTEGER)1, jc - uub); jr = jr - 1) {
                                a[((jr + irow) - 1) + (jc - 1) * lda] = conj(a[((jc - iskew * jr + ioffg) - 1) + (jr - 1) * lda]);
                            }
                        }
                    }
                    if (ipack == 6) {
                        for (jc = 1; jc <= uub; jc = jc + 1) {
                            for (jr = 1; jr <= uub + 1 - jc; jr = jr + 1) {
                                a[(jr - 1) + (jc - 1) * lda] = czero;
                            }
                        }
                    }
                    if (ipackg == 5) {
                        ipackg = ipack;
                    } else {
                        ipackg = 0;
                    }
                }
            }
            //
            //           Ensure that the diagonal is real if Hermitian
            //
            if (!zsym) {
                for (jc = 1; jc <= n; jc = jc + 1) {
                    irow = ioffst + (1 - iskew) * jc;
                    a[(irow - 1) + (jc - 1) * lda] = COMPLEX(a[(irow - 1) + (jc - 1) * lda].real());
                }
            }
            //
        }
        //
    } else {
        //
        //        4)      Generate Banded Matrix by first
        //                Rotating by random Unitary matrices,
        //                then reducing the bandwidth using Householder
        //                transformations.
        //
        //                Note: we should get here only if LDA .ge. N
        //
        if (isym == 1) {
            //
            //           Non-symmetric -- A = U D V
            //
            Clagge(mr, nc, llb, uub, d, a, lda, iseed, work, iinfo);
        } else {
            //
            //           Symmetric -- A = U D U' or
            //           Hermitian -- A = U D U*
            //
            if (zsym) {
                Clagsy(m, llb, d, a, lda, iseed, work, iinfo);
            } else {
                Claghe(m, llb, d, a, lda, iseed, work, iinfo);
            }
        }
        //
        if (iinfo != 0) {
            info = 3;
            return;
        }
    }
    //
    //     5)      Pack the matrix
    //
    INTEGER ir1 = 0;
    INTEGER ir2 = 0;
    if (ipack != ipackg) {
        if (ipack == 1) {
            //
            //           'U' -- Upper triangular, not packed
            //
            for (j = 1; j <= m; j = j + 1) {
                for (i = j + 1; i <= m; i = i + 1) {
                    a[(i - 1) + (j - 1) * lda] = czero;
                }
            }
            //
        } else if (ipack == 2) {
            //
            //           'L' -- Lower triangular, not packed
            //
            for (j = 2; j <= m; j = j + 1) {
                for (i = 1; i <= j - 1; i = i + 1) {
                    a[(i - 1) + (j - 1) * lda] = czero;
                }
            }
            //
        } else if (ipack == 3) {
            //
            //           'C' -- Upper triangle packed Columnwise.
            //
            icol = 1;
            irow = 0;
            for (j = 1; j <= m; j = j + 1) {
                for (i = 1; i <= j; i = i + 1) {
                    irow++;
                    if (irow > lda) {
                        irow = 1;
                        icol++;
                    }
                    a[(irow - 1) + (icol - 1) * lda] = a[(i - 1) + (j - 1) * lda];
                }
            }
            //
        } else if (ipack == 4) {
            //
            //           'R' -- Lower triangle packed Columnwise.
            //
            icol = 1;
            irow = 0;
            for (j = 1; j <= m; j = j + 1) {
                for (i = j; i <= m; i = i + 1) {
                    irow++;
                    if (irow > lda) {
                        irow = 1;
                        icol++;
                    }
                    a[(irow - 1) + (icol - 1) * lda] = a[(i - 1) + (j - 1) * lda];
                }
            }
            //
        } else if (ipack >= 5) {
            //
            //           'B' -- The lower triangle is packed as a band matrix.
            //           'Q' -- The upper triangle is packed as a band matrix.
            //           'Z' -- The whole matrix is packed as a band matrix.
            //
            if (ipack == 5) {
                uub = 0;
            }
            if (ipack == 6) {
                llb = 0;
            }
            //
            for (j = 1; j <= uub; j = j + 1) {
                for (i = min(j + llb, m); i >= 1; i = i - 1) {
                    a[((i - j + uub + 1) - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda];
                }
            }
            //
            for (j = uub + 2; j <= n; j = j + 1) {
                for (i = j - uub; i <= min(j + llb, m); i = i + 1) {
                    a[((i - j + uub + 1) - 1) + (j - 1) * lda] = a[(i - 1) + (j - 1) * lda];
                }
            }
        }
        //
        //        If packed, zero out extraneous elements.
        //
        //        Symmetric/Triangular Packed --
        //        zero out everything after A(IROW,ICOL)
        //
        if (ipack == 3 || ipack == 4) {
            for (jc = icol; jc <= m; jc = jc + 1) {
                for (jr = irow + 1; jr <= lda; jr = jr + 1) {
                    a[(jr - 1) + (jc - 1) * lda] = czero;
                }
                irow = 0;
            }
            //
        } else if (ipack >= 5) {
            //
            //           Packed Band --
            //              1st row is now in A( UUB+2-j, j), zero above it
            //              m-th row is now in A( M+UUB-j,j), zero below it
            //              last non-zero diagonal is now in A( UUB+LLB+1,j ),
            //                 zero below it, too.
            //
            ir1 = uub + llb + 2;
            ir2 = uub + m + 2;
            for (jc = 1; jc <= n; jc = jc + 1) {
                for (jr = 1; jr <= uub + 1 - jc; jr = jr + 1) {
                    a[(jr - 1) + (jc - 1) * lda] = czero;
                }
                for (jr = max((INTEGER)1, min(ir1, ir2 - jc)); jr <= lda; jr = jr + 1) {
                    a[(jr - 1) + (jc - 1) * lda] = czero;
                }
            }
        }
    }
    //
    //     End of Clatms
    //
}
