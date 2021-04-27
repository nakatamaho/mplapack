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

void Clattp(INTEGER const imat, const char *uplo, const char *trans, char *diag, INTEGER *iseed, INTEGER const n, COMPLEX *ap, COMPLEX *b, COMPLEX *work, REAL *rwork, INTEGER &info) {
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
    str<3> path = "Zomplex precision";
    path[(2 - 1) + (3 - 1) * ldpath] = "TP";
    REAL unfl = Rlamch("Safe minimum");
    REAL ulp = Rlamch("Epsilon") * Rlamch("Base");
    REAL smlnum = unfl;
    const REAL one = 1.0;
    REAL bignum = (one - ulp) / smlnum;
    Rlabad(smlnum, bignum);
    if ((imat >= 7 && imat <= 10) || imat == 18) {
        diag = "U";
    } else {
        diag = "N";
    }
    info = 0;
    //
    //     Quick return if N.LE.0.
    //
    if (n <= 0) {
        return;
    }
    //
    //     Call Clatb4 to set parameters for CLATMS.
    //
    bool upper = Mlsame(uplo, "U");
    char type = char0;
    INTEGER kl = 0;
    INTEGER ku = 0;
    REAL anorm = 0.0;
    INTEGER mode = 0;
    REAL cndnum = 0.0;
    char dist = char0;
    char packit = char0;
    if (upper) {
        Clatb4(path, imat, n, n, type, kl, ku, anorm, mode, cndnum, dist);
        packit = "C";
    } else {
        Clatb4(path, -imat, n, n, type, kl, ku, anorm, mode, cndnum, dist);
        packit = "R";
    }
    //
    //     IMAT <= 6:  Non-unit triangular matrix
    //
    INTEGER jc = 0;
    INTEGER j = 0;
    INTEGER i = 0;
    const REAL zero = 0.0;
    COMPLEX star1 = 0.0;
    REAL sfac = 0.0;
    COMPLEX plus1 = 0.0;
    COMPLEX plus2 = 0.0;
    REAL rexp = 0.0;
    REAL x = 0.0;
    const REAL two = 2.0e+0;
    REAL y = 0.0;
    REAL z = 0.0;
    INTEGER jcnext = 0;
    COMPLEX ra = 0.0;
    COMPLEX rb = 0.0;
    REAL c = 0.0;
    COMPLEX s = 0.0;
    INTEGER jx = 0;
    COMPLEX ctemp = 0.0;
    INTEGER iy = 0;
    REAL bnorm = 0.0;
    REAL bscal = 0.0;
    REAL tscal = 0.0;
    INTEGER jcount = 0;
    REAL texp = 0.0;
    REAL tleft = 0.0;
    if (imat <= 6) {
        zlatms(n, n, dist, iseed, type, rwork, mode, cndnum, anorm, kl, ku, packit, ap, n, work, info);
        //
        //     IMAT > 6:  Unit triangular matrix
        //     The diagonal is deliberately set to something other than 1.
        //
        //     IMAT = 7:  Matrix is the identity
        //
    } else if (imat == 7) {
        if (upper) {
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= j - 1; i = i + 1) {
                    ap[(jc + i - 1) - 1] = zero;
                }
                ap[(jc + j - 1) - 1] = j;
                jc += j;
            }
        } else {
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                ap[jc - 1] = j;
                for (i = j + 1; i <= n; i = i + 1) {
                    ap[(jc + i - j) - 1] = zero;
                }
                jc += n - j + 1;
            }
        }
        //
        //     IMAT > 7:  Non-trivial unit triangular matrix
        //
        //     Generate a unit triangular matrix T with condition CNDNUM by
        //     forming a triangular matrix with known singular values and
        //     filling in the zero entries with Givens rotations.
        //
    } else if (imat <= 10) {
        if (upper) {
            jc = 0;
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= j - 1; i = i + 1) {
                    ap[(jc + i) - 1] = zero;
                }
                ap[(jc + j) - 1] = j;
                jc += j;
            }
        } else {
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                ap[jc - 1] = j;
                for (i = j + 1; i <= n; i = i + 1) {
                    ap[(jc + i - j) - 1] = zero;
                }
                jc += n - j + 1;
            }
        }
        //
        //        Since the trace of a unit triangular matrix is 1, the product
        //        of its singular values must be 1.  Let s = sqrt(CNDNUM),
        //        x = sqrt(s) - 1/sqrt(s), y = sqrt(2/(n-2))*x, and z = x**2.
        //        The following triangular matrix has singular values s, 1, 1,
        //        ..., 1, 1/s:
        //
        //        1  y  y  y  ...  y  y  z
        //           1  0  0  ...  0  0  y
        //              1  0  ...  0  0  y
        //                 .  ...  .  .  .
        //                     .   .  .  .
        //                         1  0  y
        //                            1  y
        //                               1
        //
        //        To fill in the zeros, we first multiply by a matrix with small
        //        condition number of the form
        //
        //        1  0  0  0  0  ...
        //           1  +  *  0  0  ...
        //              1  +  0  0  0
        //                 1  +  *  0  0
        //                    1  +  0  0
        //                       ...
        //                          1  +  0
        //                             1  0
        //                                1
        //
        //        Each element marked with a '*' is formed by taking the product
        //        of the adjacent elements marked with '+'.  The '*'s can be
        //        chosen freely, and the '+'s are chosen so that the inverse of
        //        T will have elements of the same magnitude as T.  If the *'s in
        //        both T and inv(T) have small magnitude, T is well conditioned.
        //        The two offdiagonals of T are stored in WORK.
        //
        //        The product of these two matrices has the form
        //
        //        1  y  y  y  y  y  .  y  y  z
        //           1  +  *  0  0  .  0  0  y
        //              1  +  0  0  .  0  0  y
        //                 1  +  *  .  .  .  .
        //                    1  +  .  .  .  .
        //                       .  .  .  .  .
        //                          .  .  .  .
        //                             1  +  y
        //                                1  y
        //                                   1
        //
        //        Now we multiply by Givens rotations, using the fact that
        //
        //              [  c   s ] [  1   w ] [ -c  -s ] =  [  1  -w ]
        //              [ -s   c ] [  0   1 ] [  s  -c ]    [  0   1 ]
        //        and
        //              [ -c  -s ] [  1   0 ] [  c   s ] =  [  1   0 ]
        //              [  s  -c ] [  w   1 ] [ -s   c ]    [ -w   1 ]
        //
        //        where c = w / sqrt(w**2+4) and s = 2 / sqrt(w**2+4).
        //
        star1 = 0.25e0 * zlarnd(5, iseed);
        sfac = 0.5e0;
        plus1 = sfac * zlarnd(5, iseed);
        for (j = 1; j <= n; j = j + 2) {
            plus2 = star1 / plus1;
            work[j - 1] = plus1;
            work[(n + j) - 1] = star1;
            if (j + 1 <= n) {
                work[(j + 1) - 1] = plus2;
                work[(n + j + 1) - 1] = zero;
                plus1 = star1 / plus2;
                rexp = zlarnd(2, iseed);
                if (rexp < zero) {
                    star1 = -pow(sfac, [(one - rexp) - 1]) * zlarnd(5, iseed);
                } else {
                    star1 = pow(sfac, [(one + rexp) - 1]) * zlarnd(5, iseed);
                }
            }
        }
        //
        x = sqrt(cndnum) - one / sqrt(cndnum);
        if (n > 2) {
            y = sqrt(two / (n - 2).real()) * x;
        } else {
            y = zero;
        }
        z = x * x;
        //
        if (upper) {
            //
            //           Set the upper triangle of A with a unit triangular matrix
            //           of known condition number.
            //
            jc = 1;
            for (j = 2; j <= n; j = j + 1) {
                ap[(jc + 1) - 1] = y;
                if (j > 2) {
                    ap[(jc + j - 1) - 1] = work[(j - 2) - 1];
                }
                if (j > 3) {
                    ap[(jc + j - 2) - 1] = work[(n + j - 3) - 1];
                }
                jc += j;
            }
            jc = jc - n;
            ap[(jc + 1) - 1] = z;
            for (j = 2; j <= n - 1; j = j + 1) {
                ap[(jc + j) - 1] = y;
            }
        } else {
            //
            //           Set the lower triangle of A with a unit triangular matrix
            //           of known condition number.
            //
            for (i = 2; i <= n - 1; i = i + 1) {
                ap[i - 1] = y;
            }
            ap[n - 1] = z;
            jc = n + 1;
            for (j = 2; j <= n - 1; j = j + 1) {
                ap[(jc + 1) - 1] = work[(j - 1) - 1];
                if (j < n - 1) {
                    ap[(jc + 2) - 1] = work[(n + j - 1) - 1];
                }
                ap[(jc + n - j) - 1] = y;
                jc += n - j + 1;
            }
        }
        //
        //        Fill in the zeros using Givens rotations
        //
        if (upper) {
            jc = 1;
            for (j = 1; j <= n - 1; j = j + 1) {
                jcnext = jc + j;
                ra = ap[(jcnext + j - 1) - 1];
                rb = two;
                Crotg(ra, rb, c, s);
                //
                //              Multiply by [ c  s; -conj(s)  c] on the left.
                //
                if (n > j + 1) {
                    jx = jcnext + j;
                    for (i = j + 2; i <= n; i = i + 1) {
                        ctemp = c * ap[(jx + j) - 1] + s * ap[(jx + j + 1) - 1];
                        ap[(jx + j + 1) - 1] = -conj(s) * ap[(jx + j) - 1] + c * ap[(jx + j + 1) - 1];
                        ap[(jx + j) - 1] = ctemp;
                        jx += i;
                    }
                }
                //
                //              Multiply by [-c -s;  conj(s) -c] on the right.
                //
                if (j > 1) {
                    Crot(j - 1, &ap[jcnext - 1], 1, &ap[jc - 1], 1, -c, -s);
                }
                //
                //              Negate A(J,J+1).
                //
                ap[(jcnext + j - 1) - 1] = -ap[(jcnext + j - 1) - 1];
                jc = jcnext;
            }
        } else {
            jc = 1;
            for (j = 1; j <= n - 1; j = j + 1) {
                jcnext = jc + n - j + 1;
                ra = ap[(jc + 1) - 1];
                rb = two;
                Crotg(ra, rb, c, s);
                s = conj(s);
                //
                //              Multiply by [ c -s;  conj(s) c] on the right.
                //
                if (n > j + 1) {
                    Crot(n - j - 1, &ap[(jcnext + 1) - 1], 1, &ap[(jc + 2) - 1], 1, c, -s);
                }
                //
                //              Multiply by [-c  s; -conj(s) -c] on the left.
                //
                if (j > 1) {
                    jx = 1;
                    for (i = 1; i <= j - 1; i = i + 1) {
                        ctemp = -c * ap[(jx + j - i) - 1] + s * ap[(jx + j - i + 1) - 1];
                        ap[(jx + j - i + 1) - 1] = -conj(s) * ap[(jx + j - i) - 1] - c * ap[(jx + j - i + 1) - 1];
                        ap[(jx + j - i) - 1] = ctemp;
                        jx += n - i + 1;
                    }
                }
                //
                //              Negate A(J+1,J).
                //
                ap[(jc + 1) - 1] = -ap[(jc + 1) - 1];
                jc = jcnext;
            }
        }
        //
        //     IMAT > 10:  Pathological test cases.  These triangular matrices
        //     are badly scaled or badly conditioned, so when used in solving a
        //     triangular system they may cause overflow in the solution vector.
        //
    } else if (imat == 11) {
        //
        //        Type 11:  Generate a triangular matrix with elements between
        //        -1 and 1. Give the diagonal norm 2 to make it well-conditioned.
        //        Make the right hand side large so that it requires scaling.
        //
        if (upper) {
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                Clarnv(4, iseed, j - 1, &ap[jc - 1]);
                ap[(jc + j - 1) - 1] = zlarnd(5, iseed) * two;
                jc += j;
            }
        } else {
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                if (j < n) {
                    Clarnv(4, iseed, n - j, &ap[(jc + 1) - 1]);
                }
                ap[jc - 1] = zlarnd(5, iseed) * two;
                jc += n - j + 1;
            }
        }
        //
        //        Set the right hand side so that the largest value is BIGNUM.
        //
        Clarnv(2, iseed, n, b);
        iy = iCamax(n, b, 1);
        bnorm = abs(b[iy - 1]);
        bscal = bignum / max(one, bnorm);
        CRscal(n, bscal, b, 1);
        //
    } else if (imat == 12) {
        //
        //        Type 12:  Make the first diagonal element in the solve small to
        //        cause immediate overflow when dividing by T(j,j).
        //        In type 12, the offdiagonal elements are small (CNORM(j) < 1).
        //
        Clarnv(2, iseed, n, b);
        tscal = one / max(one, (n - 1).real());
        if (upper) {
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                Clarnv(4, iseed, j - 1, &ap[jc - 1]);
                CRscal(j - 1, tscal, &ap[jc - 1], 1);
                ap[(jc + j - 1) - 1] = zlarnd(5, iseed);
                jc += j;
            }
            ap[(n * (n + 1) / 2) - 1] = smlnum * ap[(n * (n + 1) / 2) - 1];
        } else {
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                Clarnv(2, iseed, n - j, &ap[(jc + 1) - 1]);
                CRscal(n - j, tscal, &ap[(jc + 1) - 1], 1);
                ap[jc - 1] = zlarnd(5, iseed);
                jc += n - j + 1;
            }
            ap[1 - 1] = smlnum * ap[1 - 1];
        }
        //
    } else if (imat == 13) {
        //
        //        Type 13:  Make the first diagonal element in the solve small to
        //        cause immediate overflow when dividing by T(j,j).
        //        In type 13, the offdiagonal elements are O(1) (CNORM(j) > 1).
        //
        Clarnv(2, iseed, n, b);
        if (upper) {
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                Clarnv(4, iseed, j - 1, &ap[jc - 1]);
                ap[(jc + j - 1) - 1] = zlarnd(5, iseed);
                jc += j;
            }
            ap[(n * (n + 1) / 2) - 1] = smlnum * ap[(n * (n + 1) / 2) - 1];
        } else {
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                Clarnv(4, iseed, n - j, &ap[(jc + 1) - 1]);
                ap[jc - 1] = zlarnd(5, iseed);
                jc += n - j + 1;
            }
            ap[1 - 1] = smlnum * ap[1 - 1];
        }
        //
    } else if (imat == 14) {
        //
        //        Type 14:  T is diagonal with small numbers on the diagonal to
        //        make the growth factor underflow, but a small right hand side
        //        chosen so that the solution does not overflow.
        //
        if (upper) {
            jcount = 1;
            jc = (n - 1) * n / 2 + 1;
            for (j = n; j >= 1; j = j - 1) {
                for (i = 1; i <= j - 1; i = i + 1) {
                    ap[(jc + i - 1) - 1] = zero;
                }
                if (jcount <= 2) {
                    ap[(jc + j - 1) - 1] = smlnum * zlarnd(5, iseed);
                } else {
                    ap[(jc + j - 1) - 1] = zlarnd(5, iseed);
                }
                jcount++;
                if (jcount > 4) {
                    jcount = 1;
                }
                jc = jc - j + 1;
            }
        } else {
            jcount = 1;
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                for (i = j + 1; i <= n; i = i + 1) {
                    ap[(jc + i - j) - 1] = zero;
                }
                if (jcount <= 2) {
                    ap[jc - 1] = smlnum * zlarnd(5, iseed);
                } else {
                    ap[jc - 1] = zlarnd(5, iseed);
                }
                jcount++;
                if (jcount > 4) {
                    jcount = 1;
                }
                jc += n - j + 1;
            }
        }
        //
        //        Set the right hand side alternately zero and small.
        //
        if (upper) {
            b[1 - 1] = zero;
            for (i = n; i >= 2; i = i - 2) {
                b[i - 1] = zero;
                b[(i - 1) - 1] = smlnum * zlarnd(5, iseed);
            }
        } else {
            b[n - 1] = zero;
            for (i = 1; i <= n - 1; i = i + 2) {
                b[i - 1] = zero;
                b[(i + 1) - 1] = smlnum * zlarnd(5, iseed);
            }
        }
        //
    } else if (imat == 15) {
        //
        //        Type 15:  Make the diagonal elements small to cause gradual
        //        overflow when dividing by T(j,j).  To control the amount of
        //        scaling needed, the matrix is bidiagonal.
        //
        texp = one / max(one, (n - 1).real());
        tscal = pow(smlnum, texp);
        Clarnv(4, iseed, n, b);
        if (upper) {
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                for (i = 1; i <= j - 2; i = i + 1) {
                    ap[(jc + i - 1) - 1] = zero;
                }
                if (j > 1) {
                    ap[(jc + j - 2) - 1] = COMPLEX(-one, -one);
                }
                ap[(jc + j - 1) - 1] = tscal * zlarnd(5, iseed);
                jc += j;
            }
            b[n - 1] = COMPLEX(one, one);
        } else {
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                for (i = j + 2; i <= n; i = i + 1) {
                    ap[(jc + i - j) - 1] = zero;
                }
                if (j < n) {
                    ap[(jc + 1) - 1] = COMPLEX(-one, -one);
                }
                ap[jc - 1] = tscal * zlarnd(5, iseed);
                jc += n - j + 1;
            }
            b[1 - 1] = COMPLEX(one, one);
        }
        //
    } else if (imat == 16) {
        //
        //        Type 16:  One zero diagonal element.
        //
        iy = n / 2 + 1;
        if (upper) {
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                Clarnv(4, iseed, j, &ap[jc - 1]);
                if (j != iy) {
                    ap[(jc + j - 1) - 1] = zlarnd(5, iseed) * two;
                } else {
                    ap[(jc + j - 1) - 1] = zero;
                }
                jc += j;
            }
        } else {
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                Clarnv(4, iseed, n - j + 1, &ap[jc - 1]);
                if (j != iy) {
                    ap[jc - 1] = zlarnd(5, iseed) * two;
                } else {
                    ap[jc - 1] = zero;
                }
                jc += n - j + 1;
            }
        }
        Clarnv(2, iseed, n, b);
        CRscal(n, two, b, 1);
        //
    } else if (imat == 17) {
        //
        //        Type 17:  Make the offdiagonal elements large to cause overflow
        //        when adding a column of T.  In the non-transposed case, the
        //        matrix is constructed to cause overflow when adding a column in
        //        every other step.
        //
        tscal = unfl / ulp;
        tscal = (one - ulp) / tscal;
        for (j = 1; j <= n * (n + 1) / 2; j = j + 1) {
            ap[j - 1] = zero;
        }
        texp = one;
        if (upper) {
            jc = (n - 1) * n / 2 + 1;
            for (j = n; j >= 2; j = j - 2) {
                ap[jc - 1] = -tscal / (n + 1).real();
                ap[(jc + j - 1) - 1] = one;
                b[j - 1] = texp * (one - ulp);
                jc = jc - j + 1;
                ap[jc - 1] = -(tscal / (n + 1).real()) / (n + 2).real();
                ap[(jc + j - 2) - 1] = one;
                b[(j - 1) - 1] = texp * (n * n + n - 1).real();
                texp = texp * two;
                jc = jc - j + 2;
            }
            b[1 - 1] = ((n + 1).real() / (n + 2).real()) * tscal;
        } else {
            jc = 1;
            for (j = 1; j <= n - 1; j = j + 2) {
                ap[(jc + n - j) - 1] = -tscal / (n + 1).real();
                ap[jc - 1] = one;
                b[j - 1] = texp * (one - ulp);
                jc += n - j + 1;
                ap[(jc + n - j - 1) - 1] = -(tscal / (n + 1).real()) / (n + 2).real();
                ap[jc - 1] = one;
                b[(j + 1) - 1] = texp * (n * n + n - 1).real();
                texp = texp * two;
                jc += n - j;
            }
            b[n - 1] = ((n + 1).real() / (n + 2).real()) * tscal;
        }
        //
    } else if (imat == 18) {
        //
        //        Type 18:  Generate a unit triangular matrix with elements
        //        between -1 and 1, and make the right hand side large so that it
        //        requires scaling.
        //
        if (upper) {
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                Clarnv(4, iseed, j - 1, &ap[jc - 1]);
                ap[(jc + j - 1) - 1] = zero;
                jc += j;
            }
        } else {
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                if (j < n) {
                    Clarnv(4, iseed, n - j, &ap[(jc + 1) - 1]);
                }
                ap[jc - 1] = zero;
                jc += n - j + 1;
            }
        }
        //
        //        Set the right hand side so that the largest value is BIGNUM.
        //
        Clarnv(2, iseed, n, b);
        iy = iCamax(n, b, 1);
        bnorm = abs(b[iy - 1]);
        bscal = bignum / max(one, bnorm);
        CRscal(n, bscal, b, 1);
        //
    } else if (imat == 19) {
        //
        //        Type 19:  Generate a triangular matrix with elements between
        //        BIGNUM/(n-1) and BIGNUM so that at least one of the column
        //        norms will exceed BIGNUM.
        //        1/3/91:  Clatps no longer can handle this case
        //
        tleft = bignum / max(one, (n - 1).real());
        tscal = bignum * ((n - 1).real() / max(one, n.real()));
        if (upper) {
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                Clarnv(5, iseed, j, &ap[jc - 1]);
                Rlarnv(1, iseed, j, rwork);
                for (i = 1; i <= j; i = i + 1) {
                    ap[(jc + i - 1) - 1] = ap[(jc + i - 1) - 1] * (tleft + rwork[i - 1] * tscal);
                }
                jc += j;
            }
        } else {
            jc = 1;
            for (j = 1; j <= n; j = j + 1) {
                Clarnv(5, iseed, n - j + 1, &ap[jc - 1]);
                Rlarnv(1, iseed, n - j + 1, rwork);
                for (i = j; i <= n; i = i + 1) {
                    ap[(jc + i - j) - 1] = ap[(jc + i - j) - 1] * (tleft + rwork[(i - j + 1) - 1] * tscal);
                }
                jc += n - j + 1;
            }
        }
        Clarnv(2, iseed, n, b);
        CRscal(n, two, b, 1);
    }
    //
    //     Flip the matrix across its counter-diagonal if the transpose will
    //     be used.
    //
    INTEGER jj = 0;
    INTEGER jr = 0;
    INTEGER jl = 0;
    REAL t = 0.0;
    if (!Mlsame(trans, "N")) {
        if (upper) {
            jj = 1;
            jr = n * (n + 1) / 2;
            for (j = 1; j <= n / 2; j = j + 1) {
                jl = jj;
                for (i = j; i <= n - j; i = i + 1) {
                    t = ap[(jr - i + j) - 1];
                    ap[(jr - i + j) - 1] = ap[jl - 1];
                    ap[jl - 1] = t;
                    jl += i;
                }
                jj += j + 1;
                jr = jr - (n - j + 1);
            }
        } else {
            jl = 1;
            jj = n * (n + 1) / 2;
            for (j = 1; j <= n / 2; j = j + 1) {
                jr = jj;
                for (i = j; i <= n - j; i = i + 1) {
                    t = ap[(jl + i - j) - 1];
                    ap[(jl + i - j) - 1] = ap[jr - 1];
                    ap[jr - 1] = t;
                    jr = jr - i;
                }
                jl += n - j + 1;
                jj = jj - j - 1;
            }
        }
    }
    //
    //     End of Clattp
    //
}
