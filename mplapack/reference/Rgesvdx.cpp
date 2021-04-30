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

void Rgesvdx(const char *jobu, const char *jobvt, const char *range, INTEGER const m, INTEGER const n, REAL *a, INTEGER const lda, REAL const vl, REAL const vu, INTEGER const il, INTEGER const iu, INTEGER &ns, REAL *s, REAL *u, INTEGER const ldu, REAL *vt, INTEGER const ldvt, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER &info) {
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
    REAL abstol = 2 * Rlamch("S");
    bool lquery = (lwork == -1);
    INTEGER minmn = min(m, n);
    //
    bool wantu = Mlsame(jobu, "V");
    bool wantvt = Mlsame(jobvt, "V");
    char jobz;
    if (wantu || wantvt) {
        jobz = 'V';
    } else {
        jobz = 'N';
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
    //     minimal amount of workspace needed at that point in the code,
    //     as well as the preferred amount for good performance.
    //     NB refers to the optimal block size for the immediately
    //     following subroutine, as returned by iMlaenv.)
    //
    INTEGER minwrk = 0;
    INTEGER maxwrk = 0;
    INTEGER mnthr = 0;
    char jobu_jobvt[3];
    jobu_jobvt[0] = jobu[0];
    jobu_jobvt[1] = jobvt[0];
    jobu_jobvt[2] = '\0';
    if (info == 0) {
        minwrk = 1;
        maxwrk = 1;
        if (minmn > 0) {
            if (m >= n) {
                mnthr = iMlaenv(6, "Rgesvd", jobu_jobvt, m, n, 0, 0);
                if (m >= mnthr) {
                    //
                    //                 Path 1 (M much larger than N)
                    //
                    maxwrk = n + n * iMlaenv(1, "Rgeqrf", " ", m, n, -1, -1);
                    maxwrk = max({maxwrk, n * (n + 5) + 2 * n * iMlaenv(1, "Rgebrd", " ", n, n, -1, -1)});
                    if (wantu) {
                        maxwrk = max({maxwrk, n * (n * 3 + 6) + n * iMlaenv(1, "Rormqr", " ", n, n, -1, -1)});
                    }
                    if (wantvt) {
                        maxwrk = max({maxwrk, n * (n * 3 + 6) + n * iMlaenv(1, "Rormlq", " ", n, n, -1, -1)});
                    }
                    minwrk = n * (n * 3 + 20);
                } else {
                    //
                    //                 Path 2 (M at least N, but not much larger)
                    //
                    maxwrk = 4 * n + (m + n) * iMlaenv(1, "Rgebrd", " ", m, n, -1, -1);
                    if (wantu) {
                        maxwrk = max({maxwrk, n * (n * 2 + 5) + n * iMlaenv(1, "Rormqr", " ", n, n, -1, -1)});
                    }
                    if (wantvt) {
                        maxwrk = max({maxwrk, n * (n * 2 + 5) + n * iMlaenv(1, "Rormlq", " ", n, n, -1, -1)});
                    }
                    minwrk = max(n * (n * 2 + 19), 4 * n + m);
                }
            } else {
                mnthr = iMlaenv(6, "Rgesvd", jobu_jobvt, m, n, 0, 0);
                if (n >= mnthr) {
                    //
                    //                 Path 1t (N much larger than M)
                    //
                    maxwrk = m + m * iMlaenv(1, "Rgelqf", " ", m, n, -1, -1);
                    maxwrk = max({maxwrk, m * (m + 5) + 2 * m * iMlaenv(1, "Rgebrd", " ", m, m, -1, -1)});
                    if (wantu) {
                        maxwrk = max({maxwrk, m * (m * 3 + 6) + m * iMlaenv(1, "Rormqr", " ", m, m, -1, -1)});
                    }
                    if (wantvt) {
                        maxwrk = max({maxwrk, m * (m * 3 + 6) + m * iMlaenv(1, "Rormlq", " ", m, m, -1, -1)});
                    }
                    minwrk = m * (m * 3 + 20);
                } else {
                    //
                    //                 Path 2t (N at least M, but not much larger)
                    //
                    maxwrk = 4 * m + (m + n) * iMlaenv(1, "Rgebrd", " ", m, n, -1, -1);
                    if (wantu) {
                        maxwrk = max({maxwrk, m * (m * 2 + 5) + m * iMlaenv(1, "Rormqr", " ", m, m, -1, -1)});
                    }
                    if (wantvt) {
                        maxwrk = max({maxwrk, m * (m * 2 + 5) + m * iMlaenv(1, "Rormlq", " ", m, m, -1, -1)});
                    }
                    minwrk = max(m * (m * 2 + 19), 4 * m + n);
                }
            }
        }
        maxwrk = max(maxwrk, minwrk);
        work[1 - 1] = castREAL(maxwrk);
        //
        if (lwork < minwrk && !lquery) {
            info = -19;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rgesvdx", -info);
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
    //     Set singular values indices accord to RANGE.
    //
    char rngtgk;
    INTEGER iltgk = 0;
    INTEGER iutgk = 0;
    if (alls) {
        rngtgk = 'I';
        iltgk = 1;
        iutgk = min(m, n);
    } else if (inds) {
        rngtgk = 'I';
        iltgk = il;
        iutgk = iu;
    } else {
        rngtgk = 'V';
        iltgk = 0;
        iutgk = 0;
    }
    //
    //     Get machine constants
    //
    REAL eps = Rlamch("P");
    REAL smlnum = sqrt(Rlamch("S")) / eps;
    const REAL one = 1.0;
    REAL bignum = one / smlnum;
    //
    //     Scale A if max element outside range [SMLNUM,BIGNUM]
    //
    REAL dum[1];
    REAL anrm = Rlange("M", m, n, a, lda, dum);
    INTEGER iscl = 0;
    if (anrm > zero && anrm < smlnum) {
        iscl = 1;
        Rlascl("G", 0, 0, anrm, smlnum, m, n, a, lda, info);
    } else if (anrm > bignum) {
        iscl = 1;
        Rlascl("G", 0, 0, anrm, bignum, m, n, a, lda, info);
    }
    //
    INTEGER itau = 0;
    INTEGER itemp = 0;
    INTEGER iqrf = 0;
    INTEGER id = 0;
    INTEGER ie = 0;
    INTEGER itauq = 0;
    INTEGER itaup = 0;
    INTEGER itgkz = 0;
    INTEGER j = 0;
    INTEGER i = 0;
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
            Rgeqrf(m, n, a, lda, &work[itau - 1], &work[itemp - 1], lwork - itemp + 1, info);
            //
            //           Copy R into WORK and bidiagonalize it:
            //           (Workspace: need N*N+5*N, prefer N*N+4*N+2*N*NB)
            //
            iqrf = itemp;
            id = iqrf + n * n;
            ie = id + n;
            itauq = ie + n;
            itaup = itauq + n;
            itemp = itaup + n;
            Rlacpy("U", n, n, a, lda, &work[iqrf - 1], n);
            Rlaset("L", n - 1, n - 1, zero, zero, &work[(iqrf + 1) - 1], n);
            Rgebrd(n, n, &work[iqrf - 1], n, &work[id - 1], &work[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[itemp - 1], lwork - itemp + 1, info);
            //
            //           Solve eigenvalue problem TGK*Z=Z*S.
            //           (Workspace: need 14*N + 2*N*(N+1))
            //
            itgkz = itemp;
            itemp = itgkz + n * (n * 2 + 1);
            Rbdsvdx("U", &jobz, &rngtgk, n, &work[id - 1], &work[ie - 1], vl, vu, iltgk, iutgk, ns, s, &work[itgkz - 1], n * 2, &work[itemp - 1], iwork, info);
            //
            //           If needed, compute left singular vectors.
            //
            if (wantu) {
                j = itgkz;
                for (i = 1; i <= ns; i = i + 1) {
                    Rcopy(n, &work[j - 1], 1, &u[(i - 1) * ldu], 1);
                    j += n * 2;
                }
                Rlaset("A", m - n, ns, zero, zero, &u[((n + 1) - 1)], ldu);
                //
                //              Call Rormbr to compute QB*UB.
                //              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
                //
                Rormbr("Q", "L", "N", n, ns, n, &work[iqrf - 1], n, &work[itauq - 1], u, ldu, &work[itemp - 1], lwork - itemp + 1, info);
                //
                //              Call Rormqr to compute Q*(QB*UB).
                //              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
                //
                Rormqr("L", "N", m, ns, n, a, lda, &work[itau - 1], u, ldu, &work[itemp - 1], lwork - itemp + 1, info);
            }
            //
            //           If needed, compute right singular vectors.
            //
            if (wantvt) {
                j = itgkz + n;
                for (i = 1; i <= ns; i = i + 1) {
                    Rcopy(n, &work[j - 1], 1, &vt[(i - 1)], ldvt);
                    j += n * 2;
                }
                //
                //              Call Rormbr to compute VB**T * PB**T
                //              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
                //
                Rormbr("P", "R", "T", ns, n, n, &work[iqrf - 1], n, &work[itaup - 1], vt, ldvt, &work[itemp - 1], lwork - itemp + 1, info);
            }
        } else {
            //
            //           Path 2 (M at least N, but not much larger)
            //           Reduce A to bidiagonal form without QR decomposition
            //           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
            //           U = QB * UB; V**T = VB**T * PB**T
            //
            //           Bidiagonalize A
            //           (Workspace: need 4*N+M, prefer 4*N+(M+N)*NB)
            //
            id = 1;
            ie = id + n;
            itauq = ie + n;
            itaup = itauq + n;
            itemp = itaup + n;
            Rgebrd(m, n, a, lda, &work[id - 1], &work[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[itemp - 1], lwork - itemp + 1, info);
            //
            //           Solve eigenvalue problem TGK*Z=Z*S.
            //           (Workspace: need 14*N + 2*N*(N+1))
            //
            itgkz = itemp;
            itemp = itgkz + n * (n * 2 + 1);
            Rbdsvdx("U", &jobz, &rngtgk, n, &work[id - 1], &work[ie - 1], vl, vu, iltgk, iutgk, ns, s, &work[itgkz - 1], n * 2, &work[itemp - 1], iwork, info);
            //
            //           If needed, compute left singular vectors.
            //
            if (wantu) {
                j = itgkz;
                for (i = 1; i <= ns; i = i + 1) {
                    Rcopy(n, &work[j - 1], 1, &u[(i - 1) * ldu], 1);
                    j += n * 2;
                }
                Rlaset("A", m - n, ns, zero, zero, &u[((n + 1) - 1)], ldu);
                //
                //              Call Rormbr to compute QB*UB.
                //              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
                //
                Rormbr("Q", "L", "N", m, ns, n, a, lda, &work[itauq - 1], u, ldu, &work[itemp - 1], lwork - itemp + 1, ierr);
            }
            //
            //           If needed, compute right singular vectors.
            //
            if (wantvt) {
                j = itgkz + n;
                for (i = 1; i <= ns; i = i + 1) {
                    Rcopy(n, &work[j - 1], 1, &vt[(i - 1)], ldvt);
                    j += n * 2;
                }
                //
                //              Call Rormbr to compute VB**T * PB**T
                //              (Workspace in WORK( ITEMP ): need N, prefer N*NB)
                //
                Rormbr("P", "R", "T", ns, n, n, a, lda, &work[itaup - 1], vt, ldvt, &work[itemp - 1], lwork - itemp + 1, ierr);
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
            Rgelqf(m, n, a, lda, &work[itau - 1], &work[itemp - 1], lwork - itemp + 1, info);
            //
            //           Copy L into WORK and bidiagonalize it:
            //           (Workspace in WORK( ITEMP ): need M*M+5*N, prefer M*M+4*M+2*M*NB)
            //
            ilqf = itemp;
            id = ilqf + m * m;
            ie = id + m;
            itauq = ie + m;
            itaup = itauq + m;
            itemp = itaup + m;
            Rlacpy("L", m, m, a, lda, &work[ilqf - 1], m);
            Rlaset("U", m - 1, m - 1, zero, zero, &work[(ilqf + m) - 1], m);
            Rgebrd(m, m, &work[ilqf - 1], m, &work[id - 1], &work[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[itemp - 1], lwork - itemp + 1, info);
            //
            //           Solve eigenvalue problem TGK*Z=Z*S.
            //           (Workspace: need 2*M*M+14*M)
            //
            itgkz = itemp;
            itemp = itgkz + m * (m * 2 + 1);
            Rbdsvdx("U", &jobz, &rngtgk, m, &work[id - 1], &work[ie - 1], vl, vu, iltgk, iutgk, ns, s, &work[itgkz - 1], m * 2, &work[itemp - 1], iwork, info);
            //
            //           If needed, compute left singular vectors.
            //
            if (wantu) {
                j = itgkz;
                for (i = 1; i <= ns; i = i + 1) {
                    Rcopy(m, &work[j - 1], 1, &u[(i - 1) * ldu], 1);
                    j += m * 2;
                }
                //
                //              Call Rormbr to compute QB*UB.
                //              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
                //
                Rormbr("Q", "L", "N", m, ns, m, &work[ilqf - 1], m, &work[itauq - 1], u, ldu, &work[itemp - 1], lwork - itemp + 1, info);
            }
            //
            //           If needed, compute right singular vectors.
            //
            if (wantvt) {
                j = itgkz + m;
                for (i = 1; i <= ns; i = i + 1) {
                    Rcopy(m, &work[j - 1], 1, &vt[(i - 1)], ldvt);
                    j += m * 2;
                }
                Rlaset("A", ns, n - m, zero, zero, &vt[((m + 1) - 1) * ldvt], ldvt);
                //
                //              Call Rormbr to compute (VB**T)*(PB**T)
                //              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
                //
                Rormbr("P", "R", "T", ns, m, m, &work[ilqf - 1], m, &work[itaup - 1], vt, ldvt, &work[itemp - 1], lwork - itemp + 1, info);
                //
                //              Call Rormlq to compute ((VB**T)*(PB**T))*Q.
                //              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
                //
                Rormlq("R", "N", ns, n, m, a, lda, &work[itau - 1], vt, ldvt, &work[itemp - 1], lwork - itemp + 1, info);
            }
        } else {
            //
            //           Path 2t (N greater than M, but not much larger)
            //           Reduce to bidiagonal form without LQ decomposition
            //           A = QB * B * PB**T = QB * ( UB * S * VB**T ) * PB**T
            //           U = QB * UB; V**T = VB**T * PB**T
            //
            //           Bidiagonalize A
            //           (Workspace: need 4*M+N, prefer 4*M+(M+N)*NB)
            //
            id = 1;
            ie = id + m;
            itauq = ie + m;
            itaup = itauq + m;
            itemp = itaup + m;
            Rgebrd(m, n, a, lda, &work[id - 1], &work[ie - 1], &work[itauq - 1], &work[itaup - 1], &work[itemp - 1], lwork - itemp + 1, info);
            //
            //           Solve eigenvalue problem TGK*Z=Z*S.
            //           (Workspace: need 2*M*M+14*M)
            //
            itgkz = itemp;
            itemp = itgkz + m * (m * 2 + 1);
            Rbdsvdx("L", &jobz, &rngtgk, m, &work[id - 1], &work[ie - 1], vl, vu, iltgk, iutgk, ns, s, &work[itgkz - 1], m * 2, &work[itemp - 1], iwork, info);
            //
            //           If needed, compute left singular vectors.
            //
            if (wantu) {
                j = itgkz;
                for (i = 1; i <= ns; i = i + 1) {
                    Rcopy(m, &work[j - 1], 1, &u[(i - 1) * ldu], 1);
                    j += m * 2;
                }
                //
                //              Call Rormbr to compute QB*UB.
                //              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
                //
                Rormbr("Q", "L", "N", m, ns, n, a, lda, &work[itauq - 1], u, ldu, &work[itemp - 1], lwork - itemp + 1, info);
            }
            //
            //           If needed, compute right singular vectors.
            //
            if (wantvt) {
                j = itgkz + m;
                for (i = 1; i <= ns; i = i + 1) {
                    Rcopy(m, &work[j - 1], 1, &vt[(i - 1)], ldvt);
                    j += m * 2;
                }
                Rlaset("A", ns, n - m, zero, zero, &vt[((m + 1) - 1) * ldvt], ldvt);
                //
                //              Call Rormbr to compute VB**T * PB**T
                //              (Workspace in WORK( ITEMP ): need M, prefer M*NB)
                //
                Rormbr("P", "R", "T", ns, n, m, a, lda, &work[itaup - 1], vt, ldvt, &work[itemp - 1], lwork - itemp + 1, info);
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
    work[1 - 1] = castREAL(maxwrk);
    //
    //     End of Rgesvdx
    //
}
