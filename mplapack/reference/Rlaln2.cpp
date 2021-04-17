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

void Rlaln2(bool const ltrans, INTEGER const na, INTEGER const nw, REAL const smin, REAL const ca, REAL *a, INTEGER const lda, REAL const d1, REAL const d2, REAL *b, INTEGER const ldb, REAL const wr, REAL const wi, REAL *x, INTEGER const ldx, REAL &scale, REAL &xnorm, INTEGER &info) {
    static bool zswap[] = {false, false, true, true};
    static bool rswap[] = {false, true, false, true};
    static INTEGER ipivot[] = {1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, 3, 2, 1};
    const REAL two = 2.0;
    REAL smlnum = two * Rlamch("Safe minimum");
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    REAL smini = max(smin, smlnum);
    //
    //     Don't check for input errors
    //
    info = 0;
    //
    //     Standard Initializations
    //
    scale = one;
    //
    REAL csr = 0.0;
    REAL cnorm = 0.0;
    REAL bnorm = 0.0;
    REAL csi = 0.0;
    const REAL zero = 0.0;
    REAL cmax = 0.0;
    INTEGER icmax = 0;
    INTEGER j = 0;
    REAL temp = 0.0;
    REAL ur11 = 0.0;
    REAL cr21 = 0.0;
    REAL ur12 = 0.0;
    REAL cr22 = 0.0;
    REAL ur11r = 0.0;
    REAL lr21 = 0.0;
    REAL ur22 = 0.0;
    REAL br1 = 0.0;
    REAL br2 = 0.0;
    REAL bbnd = 0.0;
    REAL xr2 = 0.0;
    REAL xr1 = 0.0;
    REAL ui11 = 0.0;
    REAL ci21 = 0.0;
    REAL ui12 = 0.0;
    REAL ci22 = 0.0;
    REAL ui11r = 0.0;
    REAL li21 = 0.0;
    REAL ur12s = 0.0;
    REAL ui12s = 0.0;
    REAL ui22 = 0.0;
    REAL u22abs = 0.0;
    REAL bi2 = 0.0;
    REAL bi1 = 0.0;
    REAL xi2 = 0.0;
    REAL xi1 = 0.0;
    REAL ci[4];
    REAL civ[4];
    REAL cr[4];
    REAL crv[4];
    INTEGER ldci = 2;
    INTEGER ldcr = 2;
    INTEGER ldipivot = 4;
    if (na == 1) {
        //
        //        1 x 1  (i.e., scalar) system   C X = B
        //
        if (nw == 1) {
            //
            //           Real 1x1 system.
            //
            //           C = ca A - w D
            //
            csr = ca * a[(1 - 1)] - wr * d1;
            cnorm = abs(csr);
            //
            //           If | C | < SMINI, use C = SMINI
            //
            if (cnorm < smini) {
                csr = smini;
                cnorm = smini;
                info = 1;
            }
            //
            //           Check scaling for  X = B / C
            //
            bnorm = abs(b[(1 - 1)]);
            if (cnorm < one && bnorm > one) {
                if (bnorm > bignum * cnorm) {
                    scale = one / bnorm;
                }
            }
            //
            //           Compute X
            //
            x[(1 - 1)] = (b[(1 - 1)] * scale) / csr;
            xnorm = abs(x[(1 - 1)]);
        } else {
            //
            //           Complex 1x1 system (w is complex)
            //
            //           C = ca A - w D
            //
            csr = ca * a[(1 - 1)] - wr * d1;
            csi = -wi * d1;
            cnorm = abs(csr) + abs(csi);
            //
            //           If | C | < SMINI, use C = SMINI
            //
            if (cnorm < smini) {
                csr = smini;
                csi = zero;
                cnorm = smini;
                info = 1;
            }
            //
            //           Check scaling for  X = B / C
            //
            bnorm = abs(b[(1 - 1)]) + abs(b[(2 - 1) * ldb]);
            if (cnorm < one && bnorm > one) {
                if (bnorm > bignum * cnorm) {
                    scale = one / bnorm;
                }
            }
            //
            //           Compute X
            //
            Rladiv(scale * b[(1 - 1)], scale * b[(2 - 1) * ldb], csr, csi, x[(1 - 1)], x[(2 - 1) * ldx]);
            xnorm = abs(x[(1 - 1)]) + abs(x[(2 - 1) * ldx]);
        }
        //
    } else {
        //
        //        2x2 System
        //
        //        Compute the real part of  C = ca A - w D  (or  ca A**T - w D )
        //
        cr[(1 - 1)] = ca * a[(1 - 1)] - wr * d1;
        cr[(2 - 1) + (2 - 1) * ldcr] = ca * a[(2 - 1) + (2 - 1) * lda] - wr * d2;
        if (ltrans) {
            cr[(2 - 1) * ldcr] = ca * a[(2 - 1)];
            cr[(2 - 1)] = ca * a[(2 - 1) * lda];
        } else {
            cr[(2 - 1)] = ca * a[(2 - 1)];
            cr[(2 - 1) * ldcr] = ca * a[(2 - 1) * lda];
        }
        //
        if (nw == 1) {
            //
            //           Real 2x2 system  (w is real)
            //
            //           Find the largest element in C
            //
            cmax = zero;
            icmax = 0;
            //
            for (j = 1; j <= 4; j = j + 1) {
                if (abs(crv[j - 1]) > cmax) {
                    cmax = abs(crv[j - 1]);
                    icmax = j;
                }
            }
            //
            //           If norm(C) < SMINI, use SMINI*identity.
            //
            if (cmax < smini) {
                bnorm = max(abs(b[(1 - 1)]), abs(b[(2 - 1)]));
                if (smini < one && bnorm > one) {
                    if (bnorm > bignum * smini) {
                        scale = one / bnorm;
                    }
                }
                temp = scale / smini;
                x[(1 - 1)] = temp * b[(1 - 1)];
                x[(2 - 1)] = temp * b[(2 - 1)];
                xnorm = temp * bnorm;
                info = 1;
                return;
            }
            //
            //           Gaussian elimination with complete pivoting.
            //
            ur11 = crv[icmax - 1];
            cr21 = crv[ipivot[(2 - 1) + (icmax - 1) * ldipivot] - 1];
            ur12 = crv[ipivot[(3 - 1) + (icmax - 1) * ldipivot] - 1];
            cr22 = crv[ipivot[(4 - 1) + (icmax - 1) * ldipivot] - 1];
            ur11r = one / ur11;
            lr21 = ur11r * cr21;
            ur22 = cr22 - ur12 * lr21;
            //
            //           If smaller pivot < SMINI, use SMINI
            //
            if (abs(ur22) < smini) {
                ur22 = smini;
                info = 1;
            }
            if (rswap[icmax - 1]) {
                br1 = b[(2 - 1)];
                br2 = b[(1 - 1)];
            } else {
                br1 = b[(1 - 1)];
                br2 = b[(2 - 1)];
            }
            br2 = br2 - lr21 * br1;
            bbnd = max(REAL(abs(br1 * (ur22 * ur11r))), REAL(abs(br2)));
            if (bbnd > one && abs(ur22) < one) {
                if (bbnd >= bignum * abs(ur22)) {
                    scale = one / bbnd;
                }
            }
            //
            xr2 = (br2 * scale) / ur22;
            xr1 = (scale * br1) * ur11r - xr2 * (ur11r * ur12);
            if (zswap[icmax - 1]) {
                x[(1 - 1)] = xr2;
                x[(2 - 1)] = xr1;
            } else {
                x[(1 - 1)] = xr1;
                x[(2 - 1)] = xr2;
            }
            xnorm = max(abs(xr1), abs(xr2));
            //
            //           Further scaling if  norm(A) norm(X) > overflow
            //
            if (xnorm > one && cmax > one) {
                if (xnorm > bignum / cmax) {
                    temp = cmax / bignum;
                    x[(1 - 1)] = temp * x[(1 - 1)];
                    x[(2 - 1)] = temp * x[(2 - 1)];
                    xnorm = temp * xnorm;
                    scale = temp * scale;
                }
            }
        } else {
            //
            //           Complex 2x2 system  (w is complex)
            //
            //           Find the largest element in C
            //
            ci[(1 - 1)] = -wi * d1;
            ci[(2 - 1)] = zero;
            ci[(2 - 1) * ldci] = zero;
            ci[(2 - 1) + (2 - 1) * ldci] = -wi * d2;
            cmax = zero;
            icmax = 0;
            //
            for (j = 1; j <= 4; j = j + 1) {
                if (abs(crv[j - 1]) + abs(civ[j - 1]) > cmax) {
                    cmax = abs(crv[j - 1]) + abs(civ[j - 1]);
                    icmax = j;
                }
            }
            //
            //           If norm(C) < SMINI, use SMINI*identity.
            //
            if (cmax < smini) {
                bnorm = max(abs(b[(1 - 1)]) + abs(b[(2 - 1) * ldb]), abs(b[(2 - 1)]) + abs(b[(2 - 1) + (2 - 1) * ldb]));
                if (smini < one && bnorm > one) {
                    if (bnorm > bignum * smini) {
                        scale = one / bnorm;
                    }
                }
                temp = scale / smini;
                x[(1 - 1)] = temp * b[(1 - 1)];
                x[(2 - 1)] = temp * b[(2 - 1)];
                x[(2 - 1) * ldx] = temp * b[(2 - 1) * ldb];
                x[(2 - 1) + (2 - 1) * ldx] = temp * b[(2 - 1) + (2 - 1) * ldb];
                xnorm = temp * bnorm;
                info = 1;
                return;
            }
            //
            //           Gaussian elimination with complete pivoting.
            //
            ur11 = crv[icmax - 1];
            ui11 = civ[icmax - 1];
            cr21 = crv[ipivot[(2 - 1) + (icmax - 1) * ldipivot] - 1];
            ci21 = civ[ipivot[(2 - 1) + (icmax - 1) * ldipivot] - 1];
            ur12 = crv[ipivot[(3 - 1) + (icmax - 1) * ldipivot] - 1];
            ui12 = civ[ipivot[(3 - 1) + (icmax - 1) * ldipivot] - 1];
            cr22 = crv[ipivot[(4 - 1) + (icmax - 1) * ldipivot] - 1];
            ci22 = civ[ipivot[(4 - 1) + (icmax - 1) * ldipivot] - 1];
            if (icmax == 1 || icmax == 4) {
                //
                //              Code when off-diagonals of pivoted C are real
                //
                if (abs(ur11) > abs(ui11)) {
                    temp = ui11 / ur11;
                    ur11r = one / (ur11 * (one + pow2(temp)));
                    ui11r = -temp * ur11r;
                } else {
                    temp = ur11 / ui11;
                    ui11r = -one / (ui11 * (one + pow2(temp)));
                    ur11r = -temp * ui11r;
                }
                lr21 = cr21 * ur11r;
                li21 = cr21 * ui11r;
                ur12s = ur12 * ur11r;
                ui12s = ur12 * ui11r;
                ur22 = cr22 - ur12 * lr21;
                ui22 = ci22 - ur12 * li21;
            } else {
                //
                //              Code when diagonals of pivoted C are real
                //
                ur11r = one / ur11;
                ui11r = zero;
                lr21 = cr21 * ur11r;
                li21 = ci21 * ur11r;
                ur12s = ur12 * ur11r;
                ui12s = ui12 * ur11r;
                ur22 = cr22 - ur12 * lr21 + ui12 * li21;
                ui22 = -ur12 * li21 - ui12 * lr21;
            }
            u22abs = abs(ur22) + abs(ui22);
            //
            //           If smaller pivot < SMINI, use SMINI
            //
            if (u22abs < smini) {
                ur22 = smini;
                ui22 = zero;
                info = 1;
            }
            if (rswap[icmax - 1]) {
                br2 = b[(1 - 1)];
                br1 = b[(2 - 1)];
                bi2 = b[(2 - 1) * ldb];
                bi1 = b[(2 - 1) + (2 - 1) * ldb];
            } else {
                br1 = b[(1 - 1)];
                br2 = b[(2 - 1)];
                bi1 = b[(2 - 1) * ldb];
                bi2 = b[(2 - 1) + (2 - 1) * ldb];
            }
            br2 = br2 - lr21 * br1 + li21 * bi1;
            bi2 = bi2 - li21 * br1 - lr21 * bi1;
            bbnd = max(REAL((abs(br1) + abs(bi1)) * (u22abs * (abs(ur11r) + abs(ui11r)))), REAL(abs(br2) + abs(bi2)));
            if (bbnd > one && u22abs < one) {
                if (bbnd >= bignum * u22abs) {
                    scale = one / bbnd;
                    br1 = scale * br1;
                    bi1 = scale * bi1;
                    br2 = scale * br2;
                    bi2 = scale * bi2;
                }
            }
            //
            Rladiv(br2, bi2, ur22, ui22, xr2, xi2);
            xr1 = ur11r * br1 - ui11r * bi1 - ur12s * xr2 + ui12s * xi2;
            xi1 = ui11r * br1 + ur11r * bi1 - ui12s * xr2 - ur12s * xi2;
            if (zswap[icmax - 1]) {
                x[(1 - 1)] = xr2;
                x[(2 - 1)] = xr1;
                x[(2 - 1) * ldx] = xi2;
                x[(2 - 1) + (2 - 1) * ldx] = xi1;
            } else {
                x[(1 - 1)] = xr1;
                x[(2 - 1)] = xr2;
                x[(2 - 1) * ldx] = xi1;
                x[(2 - 1) + (2 - 1) * ldx] = xi2;
            }
            xnorm = max(abs(xr1) + abs(xi1), abs(xr2) + abs(xi2));
            //
            //           Further scaling if  norm(A) norm(X) > overflow
            //
            if (xnorm > one && cmax > one) {
                if (xnorm > bignum / cmax) {
                    temp = cmax / bignum;
                    x[(1 - 1)] = temp * x[(1 - 1)];
                    x[(2 - 1)] = temp * x[(2 - 1)];
                    x[(2 - 1) * ldx] = temp * x[(2 - 1) * ldx];
                    x[(2 - 1) + (2 - 1) * ldx] = temp * x[(2 - 1) + (2 - 1) * ldx];
                    xnorm = temp * xnorm;
                    scale = temp * scale;
                }
            }
        }
    }
    //
    //     End of Rlaln2
    //
}
