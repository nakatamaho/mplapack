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

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_lin.h>

void Rlattb(INTEGER const imat, const char *uplo, const char *trans, char *diag, INTEGER *iseed, INTEGER const n, INTEGER const kd, REAL *ab, INTEGER const ldab, REAL *b, REAL *work, INTEGER &info) {
    //
    //  -- LAPACK test routine --
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
    char path[3];
    path[0] = 'R';
    path[1] = 'T';
    path[2] = 'B';
    REAL unfl = Rlamch("Safe minimum");
    REAL ulp = Rlamch("Epsilon") * Rlamch("Base");
    REAL smlnum = unfl;
    const REAL one = 1.0;
    REAL bignum = (one - ulp) / smlnum;
    Rlabad(smlnum, bignum);
    if ((imat >= 6 && imat <= 9) || imat == 17) {
        *diag = 'U';
    } else {
        *diag = 'N';
    }
    info = 0;
    //
    //     Quick return if N.LE.0.
    //
    if (n <= 0) {
        return;
    }
    //
    //     Call Rlatb4 to set parameters for SLATMS.
    //
    bool upper = Mlsame(uplo, "U");
    char type;
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cndnum = 0.0;
    char dist;
    INTEGER ioff = 0;
    char packit;
    if (upper) {
        Rlatb4(path, imat, n, n, &type, kl, ku, anorm, mode, cndnum, &dist);
        ku = kd;
        ioff = 1 + max((INTEGER)0, kd - n + 1);
        kl = 0;
        packit = 'Q';
    } else {
        Rlatb4(path, -imat, n, n, &type, kl, ku, anorm, mode, cndnum, &dist);
        kl = kd;
        ioff = 1;
        ku = 0;
        packit = 'B';
    }
    //
    //     IMAT <= 5:  Non-unit triangular matrix
    //
    INTEGER j = 0;
    INTEGER i = 0;
    const REAL zero = 0.0;
    REAL tnorm = 0.0;
    INTEGER lenj = 0;
    REAL star1 = 0.0;
    REAL sfac = 0.0;
    REAL plus1 = 0.0;
    REAL plus2 = 0.0;
    REAL rexp = 0.0;
    const REAL two = 2.0e+0;
    INTEGER iy = 0;
    REAL bnorm = 0.0;
    REAL bscal = 0.0;
    REAL tscal = 0.0;
    INTEGER jcount = 0;
    REAL texp = 0.0;
    REAL tleft = 0.0;
    if (imat <= 5) {
        Rlatms(n, n, &dist, iseed, &type, b, mode, cndnum, anorm, kl, ku, &packit, &ab[(ioff - 1)], ldab, work, info);
        //
        //     IMAT > 5:  Unit triangular matrix
        //     The diagonal is deliberately set to something other than 1.
        //
        //     IMAT = 6:  Matrix is the identity
        //
    } else if (imat == 6) {
        if (upper) {
            for (j = 1; j <= n; j = j + 1) {
                for (i = max((INTEGER)1, kd + 2 - j); i <= kd; i = i + 1) {
                    ab[(i - 1) + (j - 1) * ldab] = zero;
                }
                ab[((kd + 1) - 1) + (j - 1) * ldab] = j;
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                ab[(j - 1) * ldab] = j;
                for (i = 2; i <= min(kd + 1, n - j + 1); i = i + 1) {
                    ab[(i - 1) + (j - 1) * ldab] = zero;
                }
            }
        }
        //
        //     IMAT > 6:  Non-trivial unit triangular matrix
        //
        //     A unit triangular matrix T with condition CNDNUM is formed.
        //     In this version, T only has bandwidth 2, the rest of it is zero.
        //
    } else if (imat <= 9) {
        tnorm = sqrt(cndnum);
        //
        //        Initialize AB to zero.
        //
        if (upper) {
            for (j = 1; j <= n; j = j + 1) {
                for (i = max((INTEGER)1, kd + 2 - j); i <= kd; i = i + 1) {
                    ab[(i - 1) + (j - 1) * ldab] = zero;
                }
                ab[((kd + 1) - 1) + (j - 1) * ldab] = castREAL(j);
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                for (i = 2; i <= min(kd + 1, n - j + 1); i = i + 1) {
                    ab[(i - 1) + (j - 1) * ldab] = zero;
                }
                ab[(j - 1) * ldab] = castREAL(j);
            }
        }
        //
        //        Special case:  T is tridiagonal.  Set every other offdiagonal
        //        so that the matrix has norm TNORM+1.
        //
        if (kd == 1) {
            if (upper) {
                ab[(1 - 1) + (2 - 1) * ldab] = sign(tnorm, Rlarnd(2, iseed));
                lenj = (n - 3) / 2;
                Rlarnv(2, iseed, lenj, work);
                for (j = 1; j <= lenj; j = j + 1) {
                    ab[((2 * (j + 1)) - 1) * ldab] = tnorm * work[j - 1];
                }
            } else {
                ab[(2 - 1)] = sign(tnorm, Rlarnd(2, iseed));
                lenj = (n - 3) / 2;
                Rlarnv(2, iseed, lenj, work);
                for (j = 1; j <= lenj; j = j + 1) {
                    ab[(2 - 1) + ((2 * j + 1) - 1) * ldab] = tnorm * work[j - 1];
                }
            }
        } else if (kd > 1) {
            //
            //           Form a unit triangular matrix T with condition CNDNUM.  T is
            //           given by
            //                   | 1   +   *                      |
            //                   |     1   +                      |
            //               T = |         1   +   *              |
            //                   |             1   +              |
            //                   |                 1   +   *      |
            //                   |                     1   +      |
            //                   |                          . . . |
            //        Each element marked with a '*' is formed by taking the product
            //        of the adjacent elements marked with '+'.  The '*'s can be
            //        chosen freely, and the '+'s are chosen so that the inverse of
            //        T will have elements of the same magnitude as T.
            //
            //        The two offdiagonals of T are stored in WORK.
            //
            star1 = sign(tnorm, Rlarnd(2, iseed));
            sfac = sqrt(tnorm);
            plus1 = sign(sfac, Rlarnd(2, iseed));
            for (j = 1; j <= n; j = j + 2) {
                plus2 = star1 / plus1;
                work[j - 1] = plus1;
                work[(n + j) - 1] = star1;
                if (j + 1 <= n) {
                    work[(j + 1) - 1] = plus2;
                    work[(n + j + 1) - 1] = zero;
                    plus1 = star1 / plus2;
                    //
                    //                 Generate a new *-value with norm between sqrt(TNORM)
                    //                 and TNORM.
                    //
                    rexp = Rlarnd(2, iseed);
                    if (rexp < zero) {
                        star1 = -pow(sfac, (one - rexp));
                    } else {
                        star1 = pow(sfac, (one + rexp));
                    }
                }
            }
            //
            //           Copy the tridiagonal T to AB.
            //
            if (upper) {
                Rcopy(n - 1, work, 1, &ab[(kd - 1) + (2 - 1) * ldab], ldab);
                Rcopy(n - 2, &work[(n + 1) - 1], 1, &ab[((kd - 1) - 1) + (3 - 1) * ldab], ldab);
            } else {
                Rcopy(n - 1, work, 1, &ab[(2 - 1)], ldab);
                Rcopy(n - 2, &work[(n + 1) - 1], 1, &ab[(3 - 1)], ldab);
            }
        }
        //
        //     IMAT > 9:  Pathological test cases.  These triangular matrices
        //     are badly scaled or badly conditioned, so when used in solving a
        //     triangular system they may cause overflow in the solution vector.
        //
    } else if (imat == 10) {
        //
        //        Type 10:  Generate a triangular matrix with elements between
        //        -1 and 1. Give the diagonal norm 2 to make it well-conditioned.
        //        Make the right hand side large so that it requires scaling.
        //
        if (upper) {
            for (j = 1; j <= n; j = j + 1) {
                lenj = min(j, kd + 1);
                Rlarnv(2, iseed, lenj, &ab[((kd + 2 - lenj) - 1) + (j - 1) * ldab]);
                ab[((kd + 1) - 1) + (j - 1) * ldab] = sign(two, ab[((kd + 1) - 1) + (j - 1) * ldab]);
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                lenj = min(n - j + 1, kd + 1);
                if (lenj > 0) {
                    Rlarnv(2, iseed, lenj, &ab[(j - 1) * ldab]);
                }
                ab[(j - 1) * ldab] = sign(two, ab[(j - 1) * ldab]);
            }
        }
        //
        //        Set the right hand side so that the largest value is BIGNUM.
        //
        Rlarnv(2, iseed, n, b);
        iy = iRamax(n, b, 1);
        bnorm = abs(b[iy - 1]);
        bscal = bignum / max(one, bnorm);
        Rscal(n, bscal, b, 1);
        //
    } else if (imat == 11) {
        //
        //        Type 11:  Make the first diagonal element in the solve small to
        //        cause immediate overflow when dividing by T(j,j).
        //        In type 11, the offdiagonal elements are small (CNORM(j) < 1).
        //
        Rlarnv(2, iseed, n, b);
        tscal = one / castREAL(kd + 1);
        if (upper) {
            for (j = 1; j <= n; j = j + 1) {
                lenj = min(j, kd + 1);
                Rlarnv(2, iseed, lenj, &ab[((kd + 2 - lenj) - 1) + (j - 1) * ldab]);
                Rscal(lenj - 1, tscal, &ab[((kd + 2 - lenj) - 1) + (j - 1) * ldab], 1);
                ab[((kd + 1) - 1) + (j - 1) * ldab] = sign(one, ab[((kd + 1) - 1) + (j - 1) * ldab]);
            }
            ab[((kd + 1) - 1) + (n - 1) * ldab] = smlnum * ab[((kd + 1) - 1) + (n - 1) * ldab];
        } else {
            for (j = 1; j <= n; j = j + 1) {
                lenj = min(n - j + 1, kd + 1);
                Rlarnv(2, iseed, lenj, &ab[(j - 1) * ldab]);
                if (lenj > 1) {
                    Rscal(lenj - 1, tscal, &ab[(2 - 1) + (j - 1) * ldab], 1);
                }
                ab[(j - 1) * ldab] = sign(one, ab[(j - 1) * ldab]);
            }
            ab[(1 - 1) + (1 - 1) * ldab] = smlnum * ab[(1 - 1) + (1 - 1) * ldab];
        }
        //
    } else if (imat == 12) {
        //
        //        Type 12:  Make the first diagonal element in the solve small to
        //        cause immediate overflow when dividing by T(j,j).
        //        In type 12, the offdiagonal elements are O(1) (CNORM(j) > 1).
        //
        Rlarnv(2, iseed, n, b);
        if (upper) {
            for (j = 1; j <= n; j = j + 1) {
                lenj = min(j, kd + 1);
                Rlarnv(2, iseed, lenj, &ab[((kd + 2 - lenj) - 1) + (j - 1) * ldab]);
                ab[((kd + 1) - 1) + (j - 1) * ldab] = sign(one, ab[((kd + 1) - 1) + (j - 1) * ldab]);
            }
            ab[((kd + 1) - 1) + (n - 1) * ldab] = smlnum * ab[((kd + 1) - 1) + (n - 1) * ldab];
        } else {
            for (j = 1; j <= n; j = j + 1) {
                lenj = min(n - j + 1, kd + 1);
                Rlarnv(2, iseed, lenj, &ab[(j - 1) * ldab]);
                ab[(j - 1) * ldab] = sign(one, ab[(j - 1) * ldab]);
            }
            ab[(1 - 1) + (1 - 1) * ldab] = smlnum * ab[(1 - 1) + (1 - 1) * ldab];
        }
        //
    } else if (imat == 13) {
        //
        //        Type 13:  T is diagonal with small numbers on the diagonal to
        //        make the growth factor underflow, but a small right hand side
        //        chosen so that the solution does not overflow.
        //
        if (upper) {
            jcount = 1;
            for (j = n; j >= 1; j = j - 1) {
                for (i = max((INTEGER)1, kd + 1 - (j - 1)); i <= kd; i = i + 1) {
                    ab[(i - 1) + (j - 1) * ldab] = zero;
                }
                if (jcount <= 2) {
                    ab[((kd + 1) - 1) + (j - 1) * ldab] = smlnum;
                } else {
                    ab[((kd + 1) - 1) + (j - 1) * ldab] = one;
                }
                jcount++;
                if (jcount > 4) {
                    jcount = 1;
                }
            }
        } else {
            jcount = 1;
            for (j = 1; j <= n; j = j + 1) {
                for (i = 2; i <= min(n - j + 1, kd + 1); i = i + 1) {
                    ab[(i - 1) + (j - 1) * ldab] = zero;
                }
                if (jcount <= 2) {
                    ab[(j - 1) * ldab] = smlnum;
                } else {
                    ab[(j - 1) * ldab] = one;
                }
                jcount++;
                if (jcount > 4) {
                    jcount = 1;
                }
            }
        }
        //
        //        Set the right hand side alternately zero and small.
        //
        if (upper) {
            b[1 - 1] = zero;
            for (i = n; i >= 2; i = i - 2) {
                b[i - 1] = zero;
                b[(i - 1) - 1] = smlnum;
            }
        } else {
            b[n - 1] = zero;
            for (i = 1; i <= n - 1; i = i + 2) {
                b[i - 1] = zero;
                b[(i + 1) - 1] = smlnum;
            }
        }
        //
    } else if (imat == 14) {
        //
        //        Type 14:  Make the diagonal elements small to cause gradual
        //        overflow when dividing by T(j,j).  To control the amount of
        //        scaling needed, the matrix is bidiagonal.
        //
        texp = one / castREAL(kd + 1);
        tscal = pow(smlnum, texp);
        Rlarnv(2, iseed, n, b);
        if (upper) {
            for (j = 1; j <= n; j = j + 1) {
                for (i = max((INTEGER)1, kd + 2 - j); i <= kd; i = i + 1) {
                    ab[(i - 1) + (j - 1) * ldab] = zero;
                }
                if (j > 1 && kd > 0) {
                    ab[(kd - 1) + (j - 1) * ldab] = -one;
                }
                ab[((kd + 1) - 1) + (j - 1) * ldab] = tscal;
            }
            b[n - 1] = one;
        } else {
            for (j = 1; j <= n; j = j + 1) {
                for (i = 3; i <= min(n - j + 1, kd + 1); i = i + 1) {
                    ab[(i - 1) + (j - 1) * ldab] = zero;
                }
                if (j < n && kd > 0) {
                    ab[(2 - 1) + (j - 1) * ldab] = -one;
                }
                ab[(j - 1) * ldab] = tscal;
            }
            b[1 - 1] = one;
        }
        //
    } else if (imat == 15) {
        //
        //        Type 15:  One zero diagonal element.
        //
        iy = n / 2 + 1;
        if (upper) {
            for (j = 1; j <= n; j = j + 1) {
                lenj = min(j, kd + 1);
                Rlarnv(2, iseed, lenj, &ab[((kd + 2 - lenj) - 1) + (j - 1) * ldab]);
                if (j != iy) {
                    ab[((kd + 1) - 1) + (j - 1) * ldab] = sign(two, ab[((kd + 1) - 1) + (j - 1) * ldab]);
                } else {
                    ab[((kd + 1) - 1) + (j - 1) * ldab] = zero;
                }
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                lenj = min(n - j + 1, kd + 1);
                Rlarnv(2, iseed, lenj, &ab[(j - 1) * ldab]);
                if (j != iy) {
                    ab[(j - 1) * ldab] = sign(two, ab[(j - 1) * ldab]);
                } else {
                    ab[(j - 1) * ldab] = zero;
                }
            }
        }
        Rlarnv(2, iseed, n, b);
        Rscal(n, two, b, 1);
        //
    } else if (imat == 16) {
        //
        //        Type 16:  Make the offdiagonal elements large to cause overflow
        //        when adding a column of T.  In the non-transposed case, the
        //        matrix is constructed to cause overflow when adding a column in
        //        every other step.
        //
        tscal = unfl / ulp;
        tscal = (one - ulp) / tscal;
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= kd + 1; i = i + 1) {
                ab[(i - 1) + (j - 1) * ldab] = zero;
            }
        }
        texp = one;
        if (kd > 0) {
            if (upper) {
                for (j = n; j >= 1; j = j - kd) {
                    for (i = j; i >= max((INTEGER)1, j - kd + 1); i = i - 2) {
                        ab[((1 + (j - i)) - 1) + (i - 1) * ldab] = -tscal / castREAL(kd + 2);
                        ab[((kd + 1) - 1) + (i - 1) * ldab] = one;
                        b[i - 1] = texp * (one - ulp);
                        if (i > max((INTEGER)1, j - kd + 1)) {
                            ab[((2 + (j - i)) - 1) + ((i - 1) - 1) * ldab] = -(tscal / castREAL(kd + 2)) / castREAL(kd + 3);
                            ab[((kd + 1) - 1) + ((i - 1) - 1) * ldab] = one;
                            b[(i - 1) - 1] = texp * castREAL((kd + 1) * (kd + 1) + kd);
                        }
                        texp = texp * two;
                    }
                    b[max(n, j - kd + 1) - 1] = (castREAL(kd + 2) / castREAL(kd + 3)) * tscal;
                }
            } else {
                for (j = 1; j <= n; j = j + kd) {
                    texp = one;
                    lenj = min(kd + 1, n - j + 1);
                    for (i = j; i <= min(n, j + kd - 1); i = i + 2) {
                        ab[((lenj - (i - j)) - 1) + (j - 1) * ldab] = -tscal / castREAL(kd + 2);
                        ab[(j - 1) * ldab] = one;
                        b[j - 1] = texp * (one - ulp);
                        if (i < min(n, j + kd - 1)) {
                            ab[((lenj - (i - j + 1)) - 1) + ((i + 1) - 1) * ldab] = -(tscal / castREAL(kd + 2)) / castREAL(kd + 3);
                            ab[((i + 1) - 1) * ldab] = one;
                            b[(i + 1) - 1] = texp * castREAL((kd + 1) * (kd + 1) + kd);
                        }
                        texp = texp * two;
                    }
                    b[min(n, j + kd - 1) - 1] = (castREAL(kd + 2) / castREAL(kd + 3)) * tscal;
                }
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                ab[(j - 1) * ldab] = one;
                b[j - 1] = castREAL(j);
            }
        }
        //
    } else if (imat == 17) {
        //
        //        Type 17:  Generate a unit triangular matrix with elements
        //        between -1 and 1, and make the right hand side large so that it
        //        requires scaling.
        //
        if (upper) {
            for (j = 1; j <= n; j = j + 1) {
                lenj = min(j - 1, kd);
                Rlarnv(2, iseed, lenj, &ab[((kd + 1 - lenj) - 1) + (j - 1) * ldab]);
                ab[((kd + 1) - 1) + (j - 1) * ldab] = castREAL(j);
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                lenj = min(n - j, kd);
                if (lenj > 0) {
                    Rlarnv(2, iseed, lenj, &ab[(2 - 1) + (j - 1) * ldab]);
                }
                ab[(j - 1) * ldab] = castREAL(j);
            }
        }
        //
        //        Set the right hand side so that the largest value is BIGNUM.
        //
        Rlarnv(2, iseed, n, b);
        iy = iRamax(n, b, 1);
        bnorm = abs(b[iy - 1]);
        bscal = bignum / max(one, bnorm);
        Rscal(n, bscal, b, 1);
        //
    } else if (imat == 18) {
        //
        //        Type 18:  Generate a triangular matrix with elements between
        //        BIGNUM/KD and BIGNUM so that at least one of the column
        //        norms will exceed BIGNUM.
        //
        tleft = bignum / max(one, castREAL(kd));
        tscal = bignum * (castREAL(kd) / castREAL(kd + 1));
        if (upper) {
            for (j = 1; j <= n; j = j + 1) {
                lenj = min(j, kd + 1);
                Rlarnv(2, iseed, lenj, &ab[((kd + 2 - lenj) - 1) + (j - 1) * ldab]);
                for (i = kd + 2 - lenj; i <= kd + 1; i = i + 1) {
                    ab[(i - 1) + (j - 1) * ldab] = sign(tleft, ab[(i - 1) + (j - 1) * ldab]) + tscal * ab[(i - 1) + (j - 1) * ldab];
                }
            }
        } else {
            for (j = 1; j <= n; j = j + 1) {
                lenj = min(n - j + 1, kd + 1);
                Rlarnv(2, iseed, lenj, &ab[(j - 1) * ldab]);
                for (i = 1; i <= lenj; i = i + 1) {
                    ab[(i - 1) + (j - 1) * ldab] = sign(tleft, ab[(i - 1) + (j - 1) * ldab]) + tscal * ab[(i - 1) + (j - 1) * ldab];
                }
            }
        }
        Rlarnv(2, iseed, n, b);
        Rscal(n, two, b, 1);
    }
    //
    //     Flip the matrix if the transpose will be used.
    //
    if (!Mlsame(trans, "N")) {
        if (upper) {
            for (j = 1; j <= n / 2; j = j + 1) {
                lenj = min(n - 2 * j + 1, kd + 1);
                Rswap(lenj, &ab[((kd + 1) - 1) + (j - 1) * ldab], ldab - 1, &ab[((kd + 2 - lenj) - 1) + ((n - j + 1) - 1) * ldab], -1);
            }
        } else {
            for (j = 1; j <= n / 2; j = j + 1) {
                lenj = min(n - 2 * j + 1, kd + 1);
                Rswap(lenj, &ab[(j - 1) * ldab], 1, &ab[(lenj - 1) + ((n - j + 2 - lenj) - 1) * ldab], -ldab + 1);
            }
        }
    }
    //
    //     End of Rlattb
    //
}
