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

void Cgesvdx(const char *jobu, const char *jobvt, const char *range, INTEGER const &m, INTEGER const &n, COMPLEX *a, INTEGER const &lda, REAL const &vl, REAL const &vu, INTEGER const &il, INTEGER const &iu, INTEGER &ns, REAL *s, COMPLEX *u, INTEGER const &ldu, COMPLEX *vt, INTEGER const &ldvt, COMPLEX *work, INTEGER const &lwork, REAL *rwork, INTEGER *iwork, INTEGER &info) {
    //
    //  -- LAPACK driver routine --
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
    //     .. External Subroutines ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Test the input arguments.
    //
    ns = 0;
    info = 0;
    REAL abstol = 2 * dlamch("S");
    bool lquery = (lwork == -1);
    INTEGER minmn = min(m, n);
    //
    bool wantu = Mlsame(jobu, "V");
    bool wantvt = Mlsame(jobvt, "V");
    str<1> jobz = char0;
    if (wantu || wantvt) {
        jobz = "V";
    } else {
        jobz = "N";
    }
    bool alls = Mlsame(range, "A");
    bool vals = Mlsame(range, "V");
    bool inds = Mlsame(range, "I");
    //
    info = 0;
    const REAL zero = 0.0;
    if (!Mlsame(jobu, "V") && !Mlsame(jobu, "N")) {
        info = -1;
    } else if (!Mlsame(jobvt, "V") && !Mlsame(jobvt, "N")) {
        info = -2;
    } else if (!(alls || vals || inds)) {
        info = -3;
    } else if (m < 0) {
        info = -4;
    } else if (n < 0) {
        info = -5;
    } else if (m > lda) {
        info = -7;
    } else if (minmn > 0) {
        if (vals) {
            if (vl < zero) {
                info = -8;
            } else if (vu <= vl) {
                info = -9;
            }
        } else if (inds) {
            if (il < 1 || il > max((INTEGER)1, minmn)) {
                info = -10;
            } else if (iu < min(minmn, il) || iu > minmn) {
                info = -11;
            }
        }
        if (info == 0) {
            if (wantu && ldu < m) {
                info = -15;
            } else if (wantvt) {
                if (inds) {
                    if (ldvt < iu - il + 1) {
                        info = -17;
                    }
                } else if (ldvt < minmn) {
                    info = -17;
                }
            }
        }
    }
    //
    //     Compute workspace
    //     (Note: Comments in the code beginning "Workspace:" describe the
    //     minimal amount of workspace needed at that poINTEGER in the code,
    //     as well as the preferred amount for good performance.
    //     NB refers to the optimal block size for the immediately
    //     following subroutine, as returned by iMlaenv.)
    //
    INTEGER minwrk = 0;
    INTEGER maxwrk = 0;
    INTEGER mnthr = 0;
    if (info == 0) {
        minwrk = 1;
        maxwrk = 1;
        if (minmn > 0) {
            if (m >= n) {
                mnthr = iMlaenv[(6 - 1) + ("Cgesvd" - 1) * ldiMlaenv];
                if (m >= mnthr) {
                    //
                    //                 Path 1 (M much larger than N)
                    //
                    minwrk = n * (n + 5);
                    maxwrk = n + n * iMlaenv[("Cgeqrf" - 1) * ldiMlaenv];
                    maxwrk = max(maxwrk, n * n + 2 * n + 2 * n * iMlaenv[("Cgebrd" - 1) * ldiMlaenv]);
                    if (wantu || wantvt) {
                        maxwrk = max(maxwrk, n * n + 2 * n + n * iMlaenv[("Cunmqr" - 1) * ldiMlaenv]);
                    }
                } else {
                    //
                    //                 Path 2 (M at least N, but not much larger)
                    //
                    minwrk = 3 * n + m;
                    maxwrk = 2 * n + (m + n) * iMlaenv[("Cgebrd" - 1) * ldiMlaenv];
                    if (wantu || wantvt) {
                        maxwrk = max(maxwrk, 2 * n + n * iMlaenv[("Cunmqr" - 1) * ldiMlaenv]);
                    }
                }
            } else {
                mnthr = iMlaenv[(6 - 1) + ("Cgesvd" - 1) * ldiMlaenv];
                if (n >= mnthr) {
                    //
                    //                 Path 1t (N much larger than M)
                    //
                    minwrk = m * (m + 5);
                    maxwrk = m + m * iMlaenv[("Cgelqf" - 1) * ldiMlaenv];
                    maxwrk = max(maxwrk, m * m + 2 * m + 2 * m * iMlaenv[("Cgebrd" - 1) * ldiMlaenv]);
                    if (wantu || wantvt) {
                        maxwrk = max(maxwrk, m * m + 2 * m + m * iMlaenv[("Cunmqr" - 1) * ldiMlaenv]);
                    }
                } else {
                    //
                    //                 Path 2t (N greater than M, but not much larger)
                    //
                    minwrk = 3 * m + n;
                    maxwrk = 2 * m + (m + n) * iMlaenv[("Cgebrd" - 1) * ldiMlaenv];
                    if (wantu || wantvt) {
                        maxwrk = max(maxwrk, 2 * m + m * iMlaenv[("Cunmqr" - 1) * ldiMlaenv]);
                    }
                }
            }
        }
        maxwrk = max(maxwrk, minwrk);
        work[1 - 1] = COMPLEX(maxwrk.real(), zero);
        //
        if (lwork < minwrk && !lquery) {
            info = -19;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Cgesvdx", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (m == 0 || n == 0) {
        return;
    }
    //
    //     Set singular values indices accord to RANGE='A'.
    //
    str<1> rngtgk = char0;
    INTEGER iltgk = 0;
    INTEGER iutgk = 0;
    if (alls) {
        rngtgk = "I";
        iltgk = 1;
        iutgk = min(m, n);
    } else if (inds) {
        rngtgk = "I";
        iltgk = il;
        iutgk = iu;
    } else {
        rngtgk = "V";
        iltgk = 0;
        iutgk = 0;
    }
    //
    //     Get machine constants
    //
    REAL eps = dlamch("P");
    REAL smlnum = sqrt(dlamch("S")) / eps;
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    //
    //     Scale A if max element outside range [SMLNUM,BIGNUM]
    //
    arr_1d<1, REAL> dum(fill0);
    REAL anrm = Clange[("M" - 1) + (m - 1) * ldClange];
    INTEGER iscl = 0;
    if (anrm > zero && anrm < smlnum) {
        iscl = 1;
        Clascl("G", 0, 0, anrm, smlnum, m, n, a, lda, info);
    } else if (anrm > bignum) {
        iscl = 1;
        Clascl("G", 0, 0, anrm, bignum, m, n, a, lda, info);
    }
    //
    INTEGER itau = 0;
    INTEGER itemp = 0;
    INTEGER iqrf = 0;
    INTEGER itauq = 0;
    INTEGER itaup = 0;
    INTEGER id = 0;
    INTEGER ie = 0;
    INTEGER itgkz = 0;
    const COMPLEX czero = (0.0, 0.0);
    INTEGER itempr = 0;
    INTEGER k = 0;
    INTEGER i = 0;
    INTEGER j = 0;
    INTEGER ierr = 0;
    INTEGER ilqf = 0;
    if (m >= n) {
        //
        //        A has at least as many rows as columns. If A has sufficiently
        //        more rows than columns, first reduce A using the QR
        //        decomposition.
        //
        if (m >= mnthr) {
            //
            //           Path 1 (M much larger than N):
            //           A = Q * R = Q * ( QB * B * PB**T )
            //                     = Q * ( QB * ( UB * S * VB**T ) * PB**T )
            //           U = Q * QB * UB; V**T = VB**T * PB**T
            //
            //           Compute A=Q*R
            //           (Workspace: need 2*N, prefer N+N*NB)
            //
            itau = 1;
            itemp = itau + n;
            Cgeqrf(m, n, a, lda, work[itau - 1], work[itemp - 1], lwork - itemp + 1, info);
            //
            //           Copy R INTEGERo WORK and bidiagonalize it:
            //           (Workspace: need N*N+3*N, prefer N*N+N+2*N*NB)
            //
            iqrf = itemp;
            itauq = itemp + n * n;
            itaup = itauq + n;
            itemp = itaup + n;
            id = 1;
            ie = id + n;
            itgkz = ie + n;
            Clacpy("U", n, n, a, lda, work[iqrf - 1], n);
            Claset("L", n - 1, n - 1, czero, czero, work[(iqrf + 1) - 1], n);
            Cgebrd(n, n, work[iqrf - 1], n, rwork[id - 1], rwork[ie - 1], work[itauq - 1], work[itaup - 1], work[itemp - 1], lwork - itemp + 1, info);
            itempr = itgkz + n * (n * 2 + 1);
            //
            //           Solve eigenvalue problem TGK*Z=Z*S.
            //           (Workspace: need 2*N*N+14*N)
            //
            Rbdsvdx("U", jobz, rngtgk, n, rwork[id - 1], rwork[ie - 1], vl, vu, iltgk, iutgk, ns, s, rwork[itgkz - 1], n * 2, rwork[itempr - 1], iwork, info);
            //
            //           If needed, compute left singular vectors.
            //
            if (wantu) {
                k = itgkz;
                for (i = 1; i <= ns; i = i + 1) {
                    for (j = 1; j <= n; j = j + 1) {
                        u[(j - 1) + (i - 1) * ldu] = COMPLEX(rwork[k - 1], zero);
                        k++;
                    }
                    k += n;
                }
                Claset("A", m - n, ns, czero, czero, u[((n + 1) - 1)], ldu);
                //
                //              Call Cunmbr to compute QB*UB.
                //              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
                //
                Cunmbr("Q", "L", "N", n, ns, n, work[iqrf - 1], n, work[itauq - 1], u, ldu, work[itemp - 1], lwork - itemp + 1, info);
                //
                //              Call Cunmqr to compute Q*(QB*UB).
                //              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
                //
                Cunmqr("L", "N", m, ns, n, a, lda, work[itau - 1], u, ldu, work[itemp - 1], lwork - itemp + 1, info);
            }
            //
            //           If needed, compute right singular vectors.
            //
            if (wantvt) {
                k = itgkz + n;
                for (i = 1; i <= ns; i = i + 1) {
                    for (j = 1; j <= n; j = j + 1) {
                        vt[(i - 1) + (j - 1) * ldvt] = COMPLEX(rwork[k - 1], zero);
                        k++;
                    }
                    k += n;
                }
                //
                //              Call Cunmbr to compute VB**T * PB**T
                //              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
                //
                Cunmbr("P", "R", "C", ns, n, n, work[iqrf - 1], n, work[itaup - 1], vt, ldvt, work[itemp - 1], lwork - itemp + 1, info);
            }
        } else {
            //
            //           Path 2 (M at least N, but not much larger)
            //           Reduce A to bidiagonal form without QR decomposition
            //           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
            //           U = QB * UB; V**T = VB**T * PB**T
            //
            //           Bidiagonalize A
            //           (Workspace: need 2*N+M, prefer 2*N+(M+N)*NB)
            //
            itauq = 1;
            itaup = itauq + n;
            itemp = itaup + n;
            id = 1;
            ie = id + n;
            itgkz = ie + n;
            Cgebrd(m, n, a, lda, rwork[id - 1], rwork[ie - 1], work[itauq - 1], work[itaup - 1], work[itemp - 1], lwork - itemp + 1, info);
            itempr = itgkz + n * (n * 2 + 1);
            //
            //           Solve eigenvalue problem TGK*Z=Z*S.
            //           (Workspace: need 2*N*N+14*N)
            //
            Rbdsvdx("U", jobz, rngtgk, n, rwork[id - 1], rwork[ie - 1], vl, vu, iltgk, iutgk, ns, s, rwork[itgkz - 1], n * 2, rwork[itempr - 1], iwork, info);
            //
            //           If needed, compute left singular vectors.
            //
            if (wantu) {
                k = itgkz;
                for (i = 1; i <= ns; i = i + 1) {
                    for (j = 1; j <= n; j = j + 1) {
                        u[(j - 1) + (i - 1) * ldu] = COMPLEX(rwork[k - 1], zero);
                        k++;
                    }
                    k += n;
                }
                Claset("A", m - n, ns, czero, czero, u[((n + 1) - 1)], ldu);
                //
                //              Call Cunmbr to compute QB*UB.
                //              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
                //
                Cunmbr("Q", "L", "N", m, ns, n, a, lda, work[itauq - 1], u, ldu, work[itemp - 1], lwork - itemp + 1, ierr);
            }
            //
            //           If needed, compute right singular vectors.
            //
            if (wantvt) {
                k = itgkz + n;
                for (i = 1; i <= ns; i = i + 1) {
                    for (j = 1; j <= n; j = j + 1) {
                        vt[(i - 1) + (j - 1) * ldvt] = COMPLEX(rwork[k - 1], zero);
                        k++;
                    }
                    k += n;
                }
                //
                //              Call Cunmbr to compute VB**T * PB**T
                //              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
                //
                Cunmbr("P", "R", "C", ns, n, n, a, lda, work[itaup - 1], vt, ldvt, work[itemp - 1], lwork - itemp + 1, ierr);
            }
        }
    } else {
        //
        //        A has more columns than rows. If A has sufficiently more
        //        columns than rows, first reduce A using the LQ decomposition.
        //
        if (n >= mnthr) {
            //
            //           Path 1t (N much larger than M):
            //           A = L * Q = ( QB * B * PB**T ) * Q
            //                     = ( QB * ( UB * S * VB**T ) * PB**T ) * Q
            //           U = QB * UB ; V**T = VB**T * PB**T * Q
            //
            //           Compute A=L*Q
            //           (Workspace: need 2*M, prefer M+M*NB)
            //
            itau = 1;
            itemp = itau + m;
            Cgelqf(m, n, a, lda, work[itau - 1], work[itemp - 1], lwork - itemp + 1, info);
            //
            //           Copy L INTEGERo WORK and bidiagonalize it:
            //           (Workspace in WORK( ITEMP ): need M*M+3*M, prefer M*M+M+2*M*NB)
            //
            ilqf = itemp;
            itauq = ilqf + m * m;
            itaup = itauq + m;
            itemp = itaup + m;
            id = 1;
            ie = id + m;
            itgkz = ie + m;
            Clacpy("L", m, m, a, lda, work[ilqf - 1], m);
            Claset("U", m - 1, m - 1, czero, czero, work[(ilqf + m) - 1], m);
            Cgebrd(m, m, work[ilqf - 1], m, rwork[id - 1], rwork[ie - 1], work[itauq - 1], work[itaup - 1], work[itemp - 1], lwork - itemp + 1, info);
            itempr = itgkz + m * (m * 2 + 1);
            //
            //           Solve eigenvalue problem TGK*Z=Z*S.
            //           (Workspace: need 2*M*M+14*M)
            //
            Rbdsvdx("U", jobz, rngtgk, m, rwork[id - 1], rwork[ie - 1], vl, vu, iltgk, iutgk, ns, s, rwork[itgkz - 1], m * 2, rwork[itempr - 1], iwork, info);
            //
            //           If needed, compute left singular vectors.
            //
            if (wantu) {
                k = itgkz;
                for (i = 1; i <= ns; i = i + 1) {
                    for (j = 1; j <= m; j = j + 1) {
                        u[(j - 1) + (i - 1) * ldu] = COMPLEX(rwork[k - 1], zero);
                        k++;
                    }
                    k += m;
                }
                //
                //              Call Cunmbr to compute QB*UB.
                //              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
                //
                Cunmbr("Q", "L", "N", m, ns, m, work[ilqf - 1], m, work[itauq - 1], u, ldu, work[itemp - 1], lwork - itemp + 1, info);
            }
            //
            //           If needed, compute right singular vectors.
            //
            if (wantvt) {
                k = itgkz + m;
                for (i = 1; i <= ns; i = i + 1) {
                    for (j = 1; j <= m; j = j + 1) {
                        vt[(i - 1) + (j - 1) * ldvt] = COMPLEX(rwork[k - 1], zero);
                        k++;
                    }
                    k += m;
                }
                Claset("A", ns, n - m, czero, czero, vt[((m + 1) - 1) * ldvt], ldvt);
                //
                //              Call Cunmbr to compute (VB**T)*(PB**T)
                //              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
                //
                Cunmbr("P", "R", "C", ns, m, m, work[ilqf - 1], m, work[itaup - 1], vt, ldvt, work[itemp - 1], lwork - itemp + 1, info);
                //
                //              Call Cunmlq to compute ((VB**T)*(PB**T))*Q.
                //              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
                //
                Cunmlq("R", "N", ns, n, m, a, lda, work[itau - 1], vt, ldvt, work[itemp - 1], lwork - itemp + 1, info);
            }
        } else {
            //
            //           Path 2t (N greater than M, but not much larger)
            //           Reduce to bidiagonal form without LQ decomposition
            //           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
            //           U = QB * UB; V**T = VB**T * PB**T
            //
            //           Bidiagonalize A
            //           (Workspace: need 2*M+N, prefer 2*M+(M+N)*NB)
            //
            itauq = 1;
            itaup = itauq + m;
            itemp = itaup + m;
            id = 1;
            ie = id + m;
            itgkz = ie + m;
            Cgebrd(m, n, a, lda, rwork[id - 1], rwork[ie - 1], work[itauq - 1], work[itaup - 1], work[itemp - 1], lwork - itemp + 1, info);
            itempr = itgkz + m * (m * 2 + 1);
            //
            //           Solve eigenvalue problem TGK*Z=Z*S.
            //           (Workspace: need 2*M*M+14*M)
            //
            Rbdsvdx("L", jobz, rngtgk, m, rwork[id - 1], rwork[ie - 1], vl, vu, iltgk, iutgk, ns, s, rwork[itgkz - 1], m * 2, rwork[itempr - 1], iwork, info);
            //
            //           If needed, compute left singular vectors.
            //
            if (wantu) {
                k = itgkz;
                for (i = 1; i <= ns; i = i + 1) {
                    for (j = 1; j <= m; j = j + 1) {
                        u[(j - 1) + (i - 1) * ldu] = COMPLEX(rwork[k - 1], zero);
                        k++;
                    }
                    k += m;
                }
                //
                //              Call Cunmbr to compute QB*UB.
                //              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
                //
                Cunmbr("Q", "L", "N", m, ns, n, a, lda, work[itauq - 1], u, ldu, work[itemp - 1], lwork - itemp + 1, info);
            }
            //
            //           If needed, compute right singular vectors.
            //
            if (wantvt) {
                k = itgkz + m;
                for (i = 1; i <= ns; i = i + 1) {
                    for (j = 1; j <= m; j = j + 1) {
                        vt[(i - 1) + (j - 1) * ldvt] = COMPLEX(rwork[k - 1], zero);
                        k++;
                    }
                    k += m;
                }
                Claset("A", ns, n - m, czero, czero, vt[((m + 1) - 1) * ldvt], ldvt);
                //
                //              Call Cunmbr to compute VB**T * PB**T
                //              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
                //
                Cunmbr("P", "R", "C", ns, n, m, a, lda, work[itaup - 1], vt, ldvt, work[itemp - 1], lwork - itemp + 1, info);
            }
        }
    }
    //
    //     Undo scaling if necessary
    //
    if (iscl == 1) {
        if (anrm > bignum) {
            Rlascl("G", 0, 0, bignum, anrm, minmn, 1, s, minmn, info);
        }
        if (anrm < smlnum) {
            Rlascl("G", 0, 0, smlnum, anrm, minmn, 1, s, minmn, info);
        }
    }
    //
    //     Return optimal workspace in WORK(1)
    //
    work[1 - 1] = COMPLEX(maxwrk.real(), zero);
    //
    //     End of Cgesvdx
    //
}
