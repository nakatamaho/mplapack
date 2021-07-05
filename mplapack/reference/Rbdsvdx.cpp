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

#include <mplapack_debug.h>

void Rbdsvdx(const char *uplo, const char *jobz, const char *range, INTEGER const n, REAL *d, REAL *e, REAL const vl, REAL const vu, INTEGER const il, INTEGER const iu, INTEGER &ns, REAL *s, REAL *z, INTEGER const ldz, REAL *work, INTEGER *iwork, INTEGER &info) {
    //
    //     Test the input parameters.
    //
    bool allsv = Mlsame(range, "A");
    bool valsv = Mlsame(range, "V");
    bool indsv = Mlsame(range, "I");
    bool wantz = Mlsame(jobz, "V");
    bool lower = Mlsame(uplo, "L");
    //
    info = 0;
    const REAL zero = 0.0;
    if (!Mlsame(uplo, "U") && !lower) {
        info = -1;
    } else if (!(wantz || Mlsame(jobz, "N"))) {
        info = -2;
    } else if (!(allsv || valsv || indsv)) {
        info = -3;
    } else if (n < 0) {
        info = -4;
    } else if (n > 0) {
        if (valsv) {
            if (vl < zero) {
                info = -7;
            } else if (vu <= vl) {
                info = -8;
            }
        } else if (indsv) {
            if (il < 1 || il > max((INTEGER)1, n)) {
                info = -9;
            } else if (iu < min(n, il) || iu > n) {
                info = -10;
            }
        }
    }
    if (info == 0) {
        if (ldz < 1 || (wantz && ldz < n * 2)) {
            info = -14;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rbdsvdx", -info);
        return;
    }
    //
    //     Quick return if possible (N.LE.1)
    //
    ns = 0;
    if (n == 0) {
        return;
    }
    //
    const REAL one = 1.0;
    if (n == 1) {
        if (allsv || indsv) {
            ns = 1;
            s[1 - 1] = abs(d[1 - 1]);
        } else {
            if (vl < abs(d[1 - 1]) && vu >= abs(d[1 - 1])) {
                ns = 1;
                s[1 - 1] = abs(d[1 - 1]);
            }
        }
        if (wantz) {
            z[(1 - 1)] = sign(one, d[1 - 1]);
            z[(2 - 1)] = one;
        }
        return;
    }
    //
    REAL two = 2.0;
    REAL abstol = two * Rlamch("Safe Minimum");
    REAL ulp = Rlamch("Precision");
    REAL eps = Rlamch("Epsilon");
    REAL sqrt2 = sqrt(two);
    REAL ortol = sqrt(ulp);
    //
    //     Criterion for splitting is taken from Rbdsqr when singular
    //     values are computed to relative accuracy TOL. (See J. Demmel and
    //     W. Kahan, Accurate singular values of bidiagonal matrices, SIAM
    //     J. Sci. and Stat. Comput., 11:873912, 1990.)
    //
    const REAL ten = 10.0;
    const REAL hndrd = 100.0;
    const REAL meigth = -0.1250e0;
    REAL tol = max(ten, min(hndrd, pow(eps, meigth))) * eps;
    //
    //     Compute approximate maximum, minimum singular values.
    //
    INTEGER i = iRamax(n, d, 1);
    REAL smax = abs(d[i - 1]);
    i = iRamax(n - 1, e, 1);
    smax = max(smax, abs(e[i - 1]));
    //
    //     Compute threshold for neglecting D's and E's.
    //
    REAL smin = abs(d[1 - 1]);
    REAL mu = 0.0;
    if (smin != zero) {
        mu = smin;
        for (i = 2; i <= n; i = i + 1) {
            mu = abs(d[i - 1]) * (mu / (mu + abs(e[(i - 1) - 1])));
            smin = min(smin, mu);
            if (smin == zero) {
                break;
            }
        }
    }
    smin = smin / sqrt(castREAL(n));
    REAL thresh = tol * smin;
    //
    //     Check for zeros in D and E (splits), i.e. submatrices.
    //
    for (i = 1; i <= n - 1; i = i + 1) {
        if (abs(d[i - 1]) <= thresh) {
            d[i - 1] = zero;
        }
        if (abs(e[i - 1]) <= thresh) {
            e[i - 1] = zero;
        }
    }
    if (abs(d[n - 1]) <= thresh) {
        d[n - 1] = zero;
    }
    //
    //     Pointers for arrays used by Rstevx.
    //
    INTEGER idtgk = 1;
    INTEGER ietgk = idtgk + n * 2;
    INTEGER itemp = ietgk + n * 2;
    INTEGER iifail = 1;
    INTEGER iiwork = iifail + n * 2;
    //
    //     Set RNGVX, which corresponds to RANGE for Rstevx in TGK mode.
    //     VL,VU or IL,IU are redefined to conform to implementation a)
    //     described in the leading comments.
    //
    INTEGER iltgk = 0;
    INTEGER iutgk = 0;
    REAL vltgk = zero;
    REAL vutgk = zero;
    //
    char rngvx;
    const REAL fudge = 2.0;
    if (allsv) {
        //
        //        All singular values will be found. We aim at -s (see
        //        leading comments) with RNGVX = 'I'. IL and IU are set
        //        of the active submatrix.
        //
        rngvx = 'I';
        if (wantz) {
            Rlaset("F", n * 2, n + 1, zero, zero, z, ldz);
        }
    } else if (valsv) {
        //
        //        Find singular values in a half-open interval. We aim
        //        at -s (see leading comments) and we swap VL and VU
        //        (as VUTGK and VLTGK), changing their signs.
        //
        rngvx = 'V';
        vltgk = -vu;
        vutgk = -vl;
        for (int l = idtgk; l <= idtgk + 2 * n - 1; l++)
            work[l - 1] = zero;
        Rcopy(n, d, 1, &work[ietgk - 1], 2);
        Rcopy(n - 1, e, 1, &work[(ietgk + 1) - 1], 2);
        Rstevx("N", "V", n * 2, &work[idtgk - 1], &work[ietgk - 1], vltgk, vutgk, iltgk, iltgk, abstol, ns, s, z, ldz, &work[itemp - 1], &iwork[iiwork - 1], &iwork[iifail - 1], info);
        if (ns == 0) {
            return;
        } else {
            if (wantz) {
                Rlaset("F", n * 2, ns, zero, zero, z, ldz);
            }
        }
    } else if (indsv) {
        //
        //        Find the IL-th through the IU-th singular values. We aim
        //        at -s (see leading comments) and indices are mapped into
        //        values, therefore mimicking Rstebz, where
        //
        //        GL = GL - FUDGE*TNORM*ULP*N - FUDGE*TWO*PIVMIN
        //        GU = GU + FUDGE*TNORM*ULP*N + FUDGE*PIVMIN
        //
        iltgk = il;
        iutgk = iu;
        rngvx = 'V';
        for (int l = idtgk; l <= idtgk + 2 * n - 1; l++)
            work[l - 1] = zero;
        Rcopy(n, d, 1, &work[ietgk - 1], 2);
        Rcopy(n - 1, e, 1, &work[(ietgk + 1) - 1], 2);
        Rstevx("N", "I", n * 2, &work[idtgk - 1], &work[ietgk - 1], vltgk, vltgk, iltgk, iltgk, abstol, ns, s, z, ldz, &work[itemp - 1], &iwork[iiwork - 1], &iwork[iifail - 1], info);
        vltgk = s[1 - 1] - fudge * smax * ulp * n;
        for (int l = idtgk; l <= idtgk + 2 * n - 1; l++)
            work[l - 1] = zero;
        Rcopy(n, d, 1, &work[ietgk - 1], 2);
        Rcopy(n - 1, e, 1, &work[(ietgk + 1) - 1], 2);
        Rstevx("N", "I", n * 2, &work[idtgk - 1], &work[ietgk - 1], vutgk, vutgk, iutgk, iutgk, abstol, ns, s, z, ldz, &work[itemp - 1], &iwork[iiwork - 1], &iwork[iifail - 1], info);
        vutgk = s[1 - 1] + fudge * smax * ulp * n;
        vutgk = min(vutgk, zero);
        //
        //        If VLTGK=VUTGK, Rstevx returns an error message,
        //        so if needed we change VUTGK slightly.
        //
        if (vltgk == vutgk) {
            vltgk = vltgk - tol;
        }
        //
        if (wantz) {
            Rlaset("F", n * 2, iu - il + 1, zero, zero, z, ldz);
        }
    }
    //
    //     Initialize variables and pointers for S, Z, and WORK.
    //
    //     NRU, NRV: number of rows in U and V for the active submatrix
    //     IDBEG, ISBEG: offsets for the entries of D and S
    //     IROWZ, ICOLZ: offsets for the rows and columns of Z
    //     IROWU, IROWV: offsets for the rows of U and V
    //
    ns = 0;
    INTEGER nru = 0;
    INTEGER nrv = 0;
    INTEGER idbeg = 1;
    INTEGER isbeg = 1;
    INTEGER irowz = 1;
    INTEGER icolz = 1;
    INTEGER irowu = 2;
    INTEGER irowv = 1;
    bool split = false;
    bool sveq0 = false;
    //
    //     Form the tridiagonal TGK matrix.
    //
    for (int i = 1; i <= n; i++)
        s[i - 1] = zero;
    work[(ietgk + 2 * n - 1) - 1] = zero;
    for (int l = idtgk; l <= idtgk + 2 * n - 1; l++)
        work[l - 1] = zero;
    Rcopy(n, d, 1, &work[ietgk - 1], 2);
    Rcopy(n - 1, e, 1, &work[(ietgk + 1) - 1], 2);
    //
    //     Check for splits in two levels, outer level
    //     in E and inner level in D.
    //
    INTEGER ieptr = 0;
    INTEGER isplt = 0;
    INTEGER idend = 0;
    INTEGER idptr = 0;
    INTEGER ntgk = 0;
    INTEGER nsl = 0;
    REAL emin = 0.0;
    REAL nrmu = 0.0;
    INTEGER j = 0;
    REAL zjtji = 0.0;
    REAL nrmv = 0.0;
    for (ieptr = 2; ieptr <= n * 2; ieptr = ieptr + 2) {
        if (work[(ietgk + ieptr - 1) - 1] == zero) {
            //
            //           Split in E (this piece of B is square) or bottom
            //           of the (input bidiagonal) matrix.
            //
            isplt = idbeg;
            idend = ieptr - 1;
            for (idptr = idbeg; idptr <= idend; idptr = idptr + 2) {
                if (work[(ietgk + idptr - 1) - 1] == zero) {
                    //
                    //                 Split in D (rectangular submatrix). Set the number
                    //                 of rows in U and V (NRU and NRV) accordingly.
                    //
                    if (idptr == idbeg) {
                        //
                        //                    D=0 at the top.
                        //
                        sveq0 = true;
                        if (idbeg == idend) {
                            nru = 1;
                            nrv = 1;
                        }
                    } else if (idptr == idend) {
                        //
                        //                    D=0 at the bottom.
                        //
                        sveq0 = true;
                        nru = (idend - isplt) / 2 + 1;
                        nrv = nru;
                        if (isplt != idbeg) {
                            nru++;
                        }
                    } else {
                        if (isplt == idbeg) {
                            //
                            //                       Split: top rectangular submatrix.
                            //
                            nru = (idptr - idbeg) / 2;
                            nrv = nru + 1;
                        } else {
                            //
                            //                       Split: middle square submatrix.
                            //
                            nru = (idptr - isplt) / 2 + 1;
                            nrv = nru;
                        }
                    }
                } else if (idptr == idend) {
                    //
                    //                 Last entry of D in the active submatrix.
                    //
                    if (isplt == idbeg) {
                        //
                        //                    No split (trivial case).
                        //
                        nru = (idend - idbeg) / 2 + 1;
                        nrv = nru;
                    } else {
                        //
                        //                    Split: bottom rectangular submatrix.
                        //
                        nrv = (idend - isplt) / 2 + 1;
                        nru = nrv + 1;
                    }
                }
                //
                ntgk = nru + nrv;
                //
                if (ntgk > 0) {
                    //
                    //                 Compute eigenvalues/vectors of the active
                    //                 submatrix according to RANGE:
                    //                 if RANGE='A' (ALLSV) then RNGVX = 'I'
                    //                 if RANGE='V' (VALSV) then RNGVX = 'V'
                    //                 if RANGE='I' (INDSV) then RNGVX = 'V'
                    //
                    iltgk = 1;
                    iutgk = ntgk / 2;
                    if (allsv || vutgk == zero) {
                        if (sveq0 || smin < eps || mod(ntgk, 2) > 0) {
                            //                        Special case: eigenvalue equal to zero or very
                            //                        small, additional eigenvector is needed.
                            iutgk++;
                        }
                    }
                    //
                    //                 Workspace needed by Rstevx:
                    //                 WORK( ITEMP: ): 2*5*NTGK
                    //                 IWORK( 1: ): 2*6*NTGK
                    //
                    Rstevx(jobz, &rngvx, ntgk, &work[(idtgk + isplt - 1) - 1], &work[(ietgk + isplt - 1) - 1], vltgk, vutgk, iltgk, iutgk, abstol, nsl, &s[isbeg - 1], &z[(irowz - 1) + (icolz - 1) * ldz], ldz, &work[itemp - 1], &iwork[iiwork - 1], &iwork[iifail - 1], info);
                    if (info != 0) {
                        //                    Exit with the error code from Rstevx.
                        return;
                    }
                    emin = abs(Mmaxval(s, isbeg, isbeg + nsl - 1, 1));
                    //
                    if (nsl > 0 && wantz) {
                        //
                        //                    Normalize u=Z([2,4,...],:) and v=Z([1,3,...],:),
                        //                    changing the sign of v as discussed in the leading
                        //                    comments. The norms of u and v may be (slightly)
                        //                    different from 1/sqrt(2) if the corresponding
                        //                    eigenvalues are very small or too close. We check
                        //                    those norms and, if needed, reorthogonalize the
                        //                    vectors.
                        //
                        if (nsl > 1 && vutgk == zero && mod(ntgk, 2) == 0 && emin == 0 && !split) {
                            //
                            //                       D=0 at the top or bottom of the active submatrix:
                            //                       one eigenvalue is equal to zero; concatenate the
                            //                       eigenvectors corresponding to the two smallest
                            //                       eigenvalues.
                            //
                            for (int l = irowz; l <= irowz + ntgk - 1; l++)
                                z[(l - 1) + ((icolz + nsl - 2) - 1) * ldz] = z[(l - 1) + ((icolz + nsl - 2) - 1) * ldz] + z[(l - 1) + ((icolz + nsl - 1) - 1) * ldz];
                            for (int l = irowz; l <= irowz + ntgk - 1; l++)
                                z[(l - 1) + ((icolz + nsl - 1) - 1) * ldz] = 0.0;
                            //                       IF( IUTGK*2.GT.NTGK ) THEN
                            //                          Eigenvalue equal to zero or very small.
                            //                          NSL = NSL - 1
                            //                       END IF
                        }
                        //
                        for (i = 0; i <= min(nsl - 1, nru - 1); i = i + 1) {
                            nrmu = Rnrm2(nru, &z[(irowu - 1) + ((icolz + i) - 1) * ldz], 2);
                            if (nrmu == zero) {
                                info = n * 2 + 1;
                                return;
                            }
                            Rscal(nru, one / nrmu, &z[(irowu - 1) + ((icolz + i) - 1) * ldz], 2);
                            if (nrmu != one && abs(nrmu - ortol) * sqrt2 > one) {
                                for (j = 0; j <= i - 1; j = j + 1) {
                                    zjtji = -Rdot(nru, &z[(irowu - 1) + ((icolz + j) - 1) * ldz], 2, &z[(irowu - 1) + ((icolz + i) - 1) * ldz], 2);
                                    Raxpy(nru, zjtji, &z[(irowu - 1) + ((icolz + j) - 1) * ldz], 2, &z[(irowu - 1) + ((icolz + i) - 1) * ldz], 2);
                                }
                                nrmu = Rnrm2(nru, &z[(irowu - 1) + ((icolz + i) - 1) * ldz], 2);
                                Rscal(nru, one / nrmu, &z[(irowu - 1) + ((icolz + i) - 1) * ldz], 2);
                            }
                        }
                        for (i = 0; i <= min(nsl - 1, nrv - 1); i = i + 1) {
                            nrmv = Rnrm2(nrv, &z[(irowv - 1) + ((icolz + i) - 1) * ldz], 2);
                            if (nrmv == zero) {
                                info = n * 2 + 1;
                                return;
                            }
                            Rscal(nrv, -one / nrmv, &z[(irowv - 1) + ((icolz + i) - 1) * ldz], 2);
                            if (nrmv != one && abs(nrmv - ortol) * sqrt2 > one) {
                                for (j = 0; j <= i - 1; j = j + 1) {
                                    zjtji = -Rdot(nrv, &z[(irowv - 1) + ((icolz + j) - 1) * ldz], 2, &z[(irowv - 1) + ((icolz + i) - 1) * ldz], 2);
                                    Raxpy(nru, zjtji, &z[(irowv - 1) + ((icolz + j) - 1) * ldz], 2, &z[(irowv - 1) + ((icolz + i) - 1) * ldz], 2);
                                }
                                nrmv = Rnrm2(nrv, &z[(irowv - 1) + ((icolz + i) - 1) * ldz], 2);
                                Rscal(nrv, one / nrmv, &z[(irowv - 1) + ((icolz + i) - 1) * ldz], 2);
                            }
                        }
                        if (vutgk == zero && idptr < idend && mod(ntgk, 2) > 0) {
                            //
                            //                       D=0 in the middle of the active submatrix (one
                            //                       eigenvalue is equal to zero): save the corresponding
                            //                       eigenvector for later use (when bottom of the
                            //                       active submatrix is reached).
                            //
                            split = true;

                            for (int l = irowz; l <= irowz + ntgk - 1; l++)
                                z[l + ((n + 1) - 1) * ldz] = z[l + ((ns + nsl) - 1) * ldz];
                            for (int l = irowz; l <= irowz + ntgk - 1; l++)
                                z[l + ((ns + nsl) - 1) * ldz] = 0.0;
                        }
                        //** WANTZ **!
                    }
                    //
                    nsl = min(nsl, nru);
                    sveq0 = false;
                    //
                    //                 Absolute values of the eigenvalues of TGK.
                    //
                    for (i = 0; i <= nsl - 1; i = i + 1) {
                        s[(isbeg + i) - 1] = abs(s[(isbeg + i) - 1]);
                    }
                    //
                    //                 Update pointers for TGK, S and Z.
                    //
                    isbeg += nsl;
                    irowz += ntgk;
                    icolz += nsl;
                    irowu = irowz;
                    irowv = irowz + 1;
                    isplt = idptr + 1;
                    ns += nsl;
                    nru = 0;
                    nrv = 0;
                    //** NTGK.GT.0 **!
                }
                if (irowz < n * 2 && wantz) {
                    for (int l = 1; l <= irowz - 1; l++)
                        z[(l - 1) + (icolz - 1) * ldz] = zero;
                }
                //** IDPTR loop **!
            }
            if (split && wantz) {
                //
                //              Bring back eigenvector corresponding
                //              to eigenvalue equal to zero.
                //
                for (int l = idbeg; l <= idend - ntgk + 1; l++)
                    z[(l - 1) + (isbeg - 1) * ldz] = z[(l - 1) + (isbeg - 1) * ldz] + z[(l - 1) + ((n + 1) - 1) * ldz];
                for (int l = idbeg; l <= idend - ntgk + 1; l++)
                    z[(l - 1) + ((n + 1) - 1) * ldz] = 0.0;
            }
            irowv = irowv - 1;
            irowu++;
            idbeg = ieptr + 1;
            sveq0 = false;
            split = false;
            //** Check for split in E **!
        }
        //** IEPTR loop **!
    }
    //
    //     Sort the singular values into decreasing order (insertion sort on
    //     singular values, but only one transposition per singular vector)
    //
    INTEGER k = 0;
    for (i = 1; i <= ns - 1; i = i + 1) {
        k = 1;
        smin = s[1 - 1];
        for (j = 2; j <= ns + 1 - i; j = j + 1) {
            if (s[j - 1] <= smin) {
                k = j;
                smin = s[j - 1];
            }
        }
        if (k != ns + 1 - i) {
            s[k - 1] = s[(ns + 1 - i) - 1];
            s[(ns + 1 - i) - 1] = smin;
            if (wantz) {
                Rswap(n * 2, &z[(k - 1) * ldz], 1, &z[((ns + 1 - i) - 1) * ldz], 1);
            }
        }
    }
    //
    //     If RANGE=I, check for singular values/vectors to be discarded.
    //
    if (indsv) {
        k = iu - il + 1;
        if (k < ns) {
            for (int l = k + 1; l <= ns; l++)
                s[l - 1] = zero;
            if (wantz) {
                for (int l = 1; l <= n * 2; l++)
                    for (int m = k + 1; m <= ns; m++)
                        z[(l - 1) + (m - 1) * ldz] = zero;
            }
            ns = k;
        }
    }
    //
    //     Reorder Z: U = Z( 1:N,1:NS ), V = Z( N+1:N*2,1:NS ).
    //     If B is a lower diagonal, swap U and V.
    //
    if (wantz) {
        for (i = 1; i <= ns; i = i + 1) {
            Rcopy(n * 2, &z[(i - 1) * ldz], 1, work, 1);
            if (lower) {
                Rcopy(n, &work[2 - 1], 2, &z[((n + 1) - 1) + (i - 1) * ldz], 1);
                Rcopy(n, &work[1 - 1], 2, &z[(i - 1) * ldz], 1);
            } else {
                Rcopy(n, &work[2 - 1], 2, &z[(i - 1) * ldz], 1);
                Rcopy(n, &work[1 - 1], 2, &z[((n + 1) - 1) + (i - 1) * ldz], 1);
            }
        }
    }
    //
    //     End of Rbdsvdx
    //
}
