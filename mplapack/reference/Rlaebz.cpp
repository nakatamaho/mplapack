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

void Rlaebz(INTEGER const ijob, INTEGER const nitmax, INTEGER const n, INTEGER const mmax, INTEGER const minp, INTEGER const nbmin, REAL const abstol, REAL const reltol, REAL const pivmin, REAL *d, REAL *e, REAL *e2, INTEGER *nval, REAL *ab, REAL *c, INTEGER &mout, INTEGER *nab, REAL *work, INTEGER *iwork, INTEGER &info) {
    INTEGER ji = 0;
    INTEGER jp = 0;
    REAL tmp1 = 0.0;
    const REAL zero = 0.0;
    INTEGER j = 0;
    INTEGER kf = 0;
    INTEGER kl = 0;
    const REAL two = 2.0;
    const REAL half = 1.0 / two;
    INTEGER jit = 0;
    INTEGER klnew = 0;
    REAL tmp2 = 0.0;
    INTEGER itmp1 = 0;
    INTEGER kfnew = 0;
    INTEGER itmp2 = 0;
    INTEGER ldab = mmax;
    INTEGER ldnab = mmax;
    //
    //     Check for Errors
    //
    info = 0;
    if (ijob < 1 || ijob > 3) {
        info = -1;
        return;
    }
    //
    //     Initialize NAB
    //
    if (ijob == 1) {
        //
        //        Compute the number of eigenvalues in the initial intervals.
        //
        mout = 0;
        for (ji = 1; ji <= minp; ji = ji + 1) {
            for (jp = 1; jp <= 2; jp = jp + 1) {
                tmp1 = d[1 - 1] - ab[(ji - 1) + (jp - 1) * ldab];
                if (abs(tmp1) < pivmin) {
                    tmp1 = -pivmin;
                }
                nab[(ji - 1) + (jp - 1) * ldnab] = 0;
                if (tmp1 <= zero) {
                    nab[(ji - 1) + (jp - 1) * ldnab] = 1;
                }
                //
                for (j = 2; j <= n; j = j + 1) {
                    tmp1 = d[j - 1] - e2[(j - 1) - 1] / tmp1 - ab[(ji - 1) + (jp - 1) * ldab];
                    if (abs(tmp1) < pivmin) {
                        tmp1 = -pivmin;
                    }
                    if (tmp1 <= zero) {
                        nab[(ji - 1) + (jp - 1) * ldnab]++;
                    }
                }
            }
            mout += nab[(ji - 1) + (2 - 1) * ldnab] - nab[(ji - 1)];
        }
        return;
    }
    //
    //     Initialize for loop
    //
    //     KF and KL have the following meaning:
    //        Intervals 1,...,KF-1 have converged.
    //        Intervals KF,...,KL  still need to be refined.
    //
    kf = 1;
    kl = minp;
    //
    //     If IJOB=2, initialize C.
    //     If IJOB=3, use the user-supplied starting point.
    //
    if (ijob == 2) {
        for (ji = 1; ji <= minp; ji = ji + 1) {
            c[ji - 1] = half * (ab[(ji - 1)] + ab[(ji - 1) + (2 - 1) * ldab]);
        }
    }
    //
    //     Iteration loop
    //
    for (jit = 1; jit <= nitmax; jit = jit + 1) {
        //
        //        Loop over intervals
        //
        if (kl - kf + 1 >= nbmin && nbmin > 0) {
            //
            //           Begin of Parallel Version of the loop
            //
            for (ji = kf; ji <= kl; ji = ji + 1) {
                //
                //              Compute N(c), the number of eigenvalues less than c
                //
                work[ji - 1] = d[1 - 1] - c[ji - 1];
                iwork[ji - 1] = 0;
                if (work[ji - 1] <= pivmin) {
                    iwork[ji - 1] = 1;
                    work[ji - 1] = min(work[ji - 1], -pivmin);
                }
                //
                for (j = 2; j <= n; j = j + 1) {
                    work[ji - 1] = d[j - 1] - e2[(j - 1) - 1] / work[ji - 1] - c[ji - 1];
                    if (work[ji - 1] <= pivmin) {
                        iwork[ji - 1]++;
                        work[ji - 1] = min(work[ji - 1], -pivmin);
                    }
                }
            }
            //
            if (ijob <= 2) {
                //
                //              IJOB=2: Choose all intervals containing eigenvalues.
                //
                klnew = kl;
                for (ji = kf; ji <= kl; ji = ji + 1) {
                    //
                    //                 Insure that N(w) is monotone
                    //
                    iwork[ji - 1] = min({nab[(ji - 1) + (2 - 1) * ldnab], max(nab[(ji - 1)], iwork[ji - 1])});
                    //
                    //                 Update the Queue -- add intervals if both halves
                    //                 contain eigenvalues.
                    //
                    if (iwork[ji - 1] == nab[(ji - 1) + (2 - 1) * ldnab]) {
                        //
                        //                    No eigenvalue in the upper interval:
                        //                    just use the lower interval.
                        //
                        ab[(ji - 1) + (2 - 1) * ldab] = c[ji - 1];
                        //
                    } else if (iwork[ji - 1] == nab[(ji - 1)]) {
                        //
                        //                    No eigenvalue in the lower interval:
                        //                    just use the upper interval.
                        //
                        ab[(ji - 1)] = c[ji - 1];
                    } else {
                        klnew++;
                        if (klnew <= mmax) {
                            //
                            //                       Eigenvalue in both intervals -- add upper to
                            //                       queue.
                            //
                            ab[(klnew - 1) + (2 - 1) * ldab] = ab[(ji - 1) + (2 - 1) * ldab];
                            nab[(klnew - 1) + (2 - 1) * ldnab] = nab[(ji - 1) + (2 - 1) * ldnab];
                            ab[(klnew - 1)] = c[ji - 1];
                            nab[(klnew - 1)] = iwork[ji - 1];
                            ab[(ji - 1) + (2 - 1) * ldab] = c[ji - 1];
                            nab[(ji - 1) + (2 - 1) * ldnab] = iwork[ji - 1];
                        } else {
                            info = mmax + 1;
                        }
                    }
                }
                if (info != 0) {
                    return;
                }
                kl = klnew;
            } else {
                //
                //              IJOB=3: Binary search.  Keep only the interval containing
                //                      w   s.t. N(w) = NVAL
                //
                for (ji = kf; ji <= kl; ji = ji + 1) {
                    if (iwork[ji - 1] <= nval[ji - 1]) {
                        ab[(ji - 1)] = c[ji - 1];
                        nab[(ji - 1)] = iwork[ji - 1];
                    }
                    if (iwork[ji - 1] >= nval[ji - 1]) {
                        ab[(ji - 1) + (2 - 1) * ldab] = c[ji - 1];
                        nab[(ji - 1) + (2 - 1) * ldnab] = iwork[ji - 1];
                    }
                }
            }
            //
        } else {
            //
            //           End of Parallel Version of the loop
            //
            //           Begin of Serial Version of the loop
            //
            klnew = kl;
            for (ji = kf; ji <= kl; ji = ji + 1) {
                //
                //              Compute N(w), the number of eigenvalues less than w
                //
                tmp1 = c[ji - 1];
                tmp2 = d[1 - 1] - tmp1;
                itmp1 = 0;
                if (tmp2 <= pivmin) {
                    itmp1 = 1;
                    tmp2 = min(tmp2, -pivmin);
                }
                //
                for (j = 2; j <= n; j = j + 1) {
                    tmp2 = d[j - 1] - e2[(j - 1) - 1] / tmp2 - tmp1;
                    if (tmp2 <= pivmin) {
                        itmp1++;
                        tmp2 = min(tmp2, -pivmin);
                    }
                }
                //
                if (ijob <= 2) {
                    //
                    //                 IJOB=2: Choose all intervals containing eigenvalues.
                    //
                    //                 Insure that N(w) is monotone
                    //
                    itmp1 = min({nab[(ji - 1) + (2 - 1) * ldnab], max(nab[(ji - 1)], itmp1)});
                    //
                    //                 Update the Queue -- add intervals if both halves
                    //                 contain eigenvalues.
                    //
                    if (itmp1 == nab[(ji - 1) + (2 - 1) * ldnab]) {
                        //
                        //                    No eigenvalue in the upper interval:
                        //                    just use the lower interval.
                        //
                        ab[(ji - 1) + (2 - 1) * ldab] = tmp1;
                        //
                    } else if (itmp1 == nab[(ji - 1)]) {
                        //
                        //                    No eigenvalue in the lower interval:
                        //                    just use the upper interval.
                        //
                        ab[(ji - 1)] = tmp1;
                    } else if (klnew < mmax) {
                        //
                        //                    Eigenvalue in both intervals -- add upper to queue.
                        //
                        klnew++;
                        ab[(klnew - 1) + (2 - 1) * ldab] = ab[(ji - 1) + (2 - 1) * ldab];
                        nab[(klnew - 1) + (2 - 1) * ldnab] = nab[(ji - 1) + (2 - 1) * ldnab];
                        ab[(klnew - 1)] = tmp1;
                        nab[(klnew - 1)] = itmp1;
                        ab[(ji - 1) + (2 - 1) * ldab] = tmp1;
                        nab[(ji - 1) + (2 - 1) * ldnab] = itmp1;
                    } else {
                        info = mmax + 1;
                        return;
                    }
                } else {
                    //
                    //                 IJOB=3: Binary search.  Keep only the interval
                    //                         containing  w  s.t. N(w) = NVAL
                    //
                    if (itmp1 <= nval[ji - 1]) {
                        ab[(ji - 1)] = tmp1;
                        nab[(ji - 1)] = itmp1;
                    }
                    if (itmp1 >= nval[ji - 1]) {
                        ab[(ji - 1) + (2 - 1) * ldab] = tmp1;
                        nab[(ji - 1) + (2 - 1) * ldnab] = itmp1;
                    }
                }
            }
            kl = klnew;
            //
        }
        //
        //        Check for convergence
        //
        kfnew = kf;
        for (ji = kf; ji <= kl; ji = ji + 1) {
            tmp1 = abs(ab[(ji - 1) + (2 - 1) * ldab] - ab[(ji - 1)]);
            tmp2 = max(abs(ab[(ji - 1) + (2 - 1) * ldab]), abs(ab[(ji - 1)]));
            if (tmp1 < max({abstol, pivmin, reltol * tmp2}) || nab[(ji - 1)] >= nab[(ji - 1) + (2 - 1) * ldnab]) {
                //
                //              Converged -- Swap with position KFNEW,
                //                           then increment KFNEW
                //
                if (ji > kfnew) {
                    tmp1 = ab[(ji - 1)];
                    tmp2 = ab[(ji - 1) + (2 - 1) * ldab];
                    itmp1 = nab[(ji - 1)];
                    itmp2 = nab[(ji - 1) + (2 - 1) * ldnab];
                    ab[(ji - 1)] = ab[(kfnew - 1)];
                    ab[(ji - 1) + (2 - 1) * ldab] = ab[(kfnew - 1) + (2 - 1) * ldab];
                    nab[(ji - 1)] = nab[(kfnew - 1)];
                    nab[(ji - 1) + (2 - 1) * ldnab] = nab[(kfnew - 1) + (2 - 1) * ldnab];
                    ab[(kfnew - 1)] = tmp1;
                    ab[(kfnew - 1) + (2 - 1) * ldab] = tmp2;
                    nab[(kfnew - 1)] = itmp1;
                    nab[(kfnew - 1) + (2 - 1) * ldnab] = itmp2;
                    if (ijob == 3) {
                        itmp1 = nval[ji - 1];
                        nval[ji - 1] = nval[kfnew - 1];
                        nval[kfnew - 1] = itmp1;
                    }
                }
                kfnew++;
            }
        }
        kf = kfnew;
        //
        //        Choose Midpoints
        //
        for (ji = kf; ji <= kl; ji = ji + 1) {
            c[ji - 1] = half * (ab[(ji - 1)] + ab[(ji - 1) + (2 - 1) * ldab]);
        }
        //
        //        If no more intervals to refine, quit.
        //
        if (kf > kl) {
            goto statement_140;
        }
    }
//
//     Converged
//
statement_140:
    info = max(kl + 1 - kf, (INTEGER)0);
    mout = kl;
    //
    //     End of Rlaebz
    //
}
