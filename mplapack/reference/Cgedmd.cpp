/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine ZGEDMD.
// Original LAPACK authors:
//   Univ. of Tennessee,
//   Univ. of California Berkeley,
//   Univ. of Colorado Denver,
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

void Cgedmd(const char *jobs, const char *jobz, const char *jobr, const char *jobf, INTEGER const whtsvd, INTEGER const m, INTEGER const n, COMPLEX *x, INTEGER const ldx, COMPLEX *y, INTEGER const ldy, INTEGER const nrnk, REAL const tol, INTEGER &k, COMPLEX *eigs, COMPLEX *z, INTEGER const ldz, REAL *res, COMPLEX *b, INTEGER const ldb, COMPLEX *w, INTEGER const ldw, COMPLEX *s, INTEGER const lds, COMPLEX *zwork, INTEGER const lzwork, REAL *rwork, INTEGER const lrwork, INTEGER *iwork, INTEGER const liwork, INTEGER &info) {
    //
    //
    // Test the input arguments
    //
    bool wntres = Mlsame(jobr, "R");
    bool sccolx = Mlsame(jobs, "S") || Mlsame(jobs, "C");
    bool sccoly = Mlsame(jobs, "Y");
    bool wntvec = Mlsame(jobz, "V");
    bool wntref = Mlsame(jobf, "R");
    bool wntex = Mlsame(jobf, "E");
    info = 0;
    bool lquery = ((lzwork == -1) || (liwork == -1) || (lrwork == -1));
    //
    const REAL zero = 0.0;
    const REAL one = 1.0;
    if (!(sccolx || sccoly || Mlsame(jobs, "N"))) {
        info = -1;
    } else if (!(wntvec || Mlsame(jobz, "N") || Mlsame(jobz, "F"))) {
        info = -2;
    } else if (!(wntres || Mlsame(jobr, "N")) || (wntres && (!wntvec))) {
        info = -3;
    } else if (!(wntref || wntex || Mlsame(jobf, "N"))) {
        info = -4;
    } else if (!((whtsvd == 1) || (whtsvd == 2) || (whtsvd == 3) || (whtsvd == 4))) {
        info = -5;
    } else if (m < 0) {
        info = -6;
    } else if ((n < 0) || (n > m)) {
        info = -7;
    } else if (ldx < m) {
        info = -9;
    } else if (ldy < m) {
        info = -11;
    } else if (!((nrnk == -2) || (nrnk == -1) || ((nrnk >= 1) && (nrnk <= n)))) {
        info = -12;
    } else if ((tol < zero) || (tol >= one)) {
        info = -13;
    } else if (ldz < m) {
        info = -17;
    } else if ((wntref || wntex) && (ldb < m)) {
        info = -20;
    } else if (ldw < n) {
        info = -22;
    } else if (lds < n) {
        info = -24;
    }
    //
    INTEGER iminwr = 0;
    INTEGER mlrwrk = 0;
    INTEGER mlwork = 0;
    INTEGER olwork = 0;
    INTEGER mwrsvd = 0;
    REAL rdummy[2];
    INTEGER info1 = 0;
    INTEGER lwrsvd = 0;
    INTEGER mwrsdd = 0;
    INTEGER lwrsdd = 0;
    INTEGER numrnk = 0;
    INTEGER mwrsvq = 0;
    INTEGER lwrsvq = 0;
    char jsvopt;
    INTEGER mwrsvj = 0;
    INTEGER lwrsvj = 0;
    char jobzl;
    INTEGER mwrkev = 0;
    INTEGER lwrkev = 0;
    if (info == 0) {
        // Compute the minimal and the optimal workspace
        // requirements. Simulate running the code and
        // determine minimal and optimal sizes of the
        // workspace at any moment of the run.
        if (n == 0) {
            // Quick return. All output except K is void.
            // INFO=1 signals the void input.
            // In case of a workspace query, the default
            // minimal workspace lengths are returned.
            if (lquery) {
                iwork[1 - 1] = 1;
                rwork[1 - 1] = 1.0;
                zwork[1 - 1] = 2.0;
                zwork[2 - 1] = 2.0;
            } else {
                k = 0;
            }
            info = 1;
            return;
        }
        //
        iminwr = 1;
        mlrwrk = max((INTEGER)1, n);
        mlwork = 2;
        olwork = 2;
        if ((whtsvd == 1)) {
            // The following is specified as the minimal
            // length of WORK in the definition of Cgesvd:
            // MWRSVD = MAX(1,2*MIN(M,N)+MAX(M,N))
            mwrsvd = max((INTEGER)1, 2 * min(m, n) + max(m, n));
            mlwork = max(mlwork, mwrsvd);
            mlrwrk = max(mlrwrk, n + 5 * min(m, n));
            if (lquery) {
                Cgesvd("O", "S", m, n, x, ldx, rwork, b, ldb, w, ldw, zwork, -1, rdummy, info1);
                lwrsvd = castINTEGER(zwork[1 - 1].real());
                olwork = max(olwork, lwrsvd);
            }
        } else if ((whtsvd == 2)) {
            // The following is specified as the minimal
            // length of WORK in the definition of Cgesdd:
            // MWRSDD = 2*min(M,N)*min(M,N)+2*min(M,N)+max(M,N).
            // RWORK length: 5*MIN(M,N)*MIN(M,N)+7*MIN(M,N)
            // In LAPACK 3.10.1 RWORK is defined differently.
            // Below we take max over the two versions.
            // IMINWR = 8*MIN(M,N)
            mwrsdd = 2 * min(m, n) * min(m, n) + 2 * min(m, n) + max(m, n);
            mlwork = max(mlwork, mwrsdd);
            iminwr = 8 * min(m, n);
            mlrwrk = max(mlrwrk, n + max(5 * min(m, n) * min(m, n) + 7 * min(m, n), 5 * min(m, n) * min(m, n) + 5 * min(m, n), 2 * max(m, n) * min(m, n) + 2 * min(m, n) * min(m, n) + min(m, n)));
            if (lquery) {
                Cgesdd("O", m, n, x, ldx, rwork, b, ldb, w, ldw, zwork, -1, rdummy, iwork, info1);
                lwrsdd = max(mwrsdd, castINTEGER(zwork[1 - 1].real()));
                // Possible bug in Cgesdd optimal workspace size.
                olwork = max(olwork, lwrsdd);
            }
        } else if ((whtsvd == 3)) {
            Cgesvdq("H", "P", "N", "R", "R", m, n, x, ldx, rwork, z, ldz, w, ldw, numrnk, iwork, -1, zwork, -1, rdummy, -1, info1);
            iminwr = iwork[1 - 1];
            mwrsvq = castINTEGER(zwork[2 - 1].real());
            mlwork = max(mlwork, mwrsvq);
            mlrwrk = max(mlrwrk, n + castINTEGER(rdummy[1 - 1]));
            if (lquery) {
                lwrsvq = castINTEGER(zwork[1 - 1].real());
                olwork = max(olwork, lwrsvq);
            }
        } else if ((whtsvd == 4)) {
            jsvopt = 'J';
            Cgejsv("F", "U", &jsvopt, "R", "N", "P", m, n, x, ldx, rwork, z, ldz, w, ldw, zwork, -1, rdummy, -1, iwork, info1);
            iminwr = iwork[1 - 1];
            mwrsvj = castINTEGER(zwork[2 - 1].real());
            mlwork = max(mlwork, mwrsvj);
            mlrwrk = max(mlrwrk, n + max((INTEGER)7, castINTEGER(rdummy[1 - 1])));
            if (lquery) {
                lwrsvj = castINTEGER(zwork[1 - 1].real());
                olwork = max(olwork, lwrsvj);
            }
        }
        if (wntvec || wntex || Mlsame(jobz, "F")) {
            jobzl = 'V';
        } else {
            jobzl = 'N';
        }
        // Workspace calculation to the Cgeev call
        mwrkev = max((INTEGER)1, 2 * n);
        mlwork = max(mlwork, mwrkev);
        mlrwrk = max(mlrwrk, n + 2 * n);
        if (lquery) {
            Cgeev("N", &jobzl, n, s, lds, eigs, w, ldw, w, ldw, zwork, -1, rwork, info1);
            lwrkev = castINTEGER(zwork[1 - 1].real());
            olwork = max(olwork, lwrkev);
        }
        //
        if (liwork < iminwr && (!lquery)) {
            info = -30;
        }
        if (lrwork < mlrwrk && (!lquery)) {
            info = -28;
        }
        if (lzwork < mlwork && (!lquery)) {
            info = -26;
        }
        //
    }
    //
    if (info != 0) {
        Mxerbla("Cgedmd", -info);
        return;
    } else if (lquery) {
        // Return minimal and optimal workspace sizes
        iwork[1 - 1] = iminwr;
        rwork[1 - 1] = mlrwrk;
        zwork[1 - 1] = mlwork;
        zwork[2 - 1] = olwork;
        return;
    }
    // ............................................................
    //
    REAL ofl = Rlamch("O");
    REAL small = Rlamch("S");
    bool badxy = false;
    //
    // <1> Optional scaling of the snapshots (columns of X, Y)
    // ==========================================================
    INTEGER i = 0;
    REAL ssum = 0.0;
    REAL scale = 0.0;
    REAL rootsc = 0.0;
    INTEGER info2 = 0;
    if (sccolx) {
        // The columns of X will be normalized.
        // To prevent overflows, the column norms of X are
        // carefully computed using Classq.
        k = 0;
        for (i = 1; i <= n; i = i + 1) {
            // WORK(i) = RCnrm2( M, X(1,i), 1 )
            ssum = one;
            scale = zero;
            Classq(m, &x[(i - 1) * ldx], 1, scale, ssum);
            if (Risnan(scale) || Risnan(ssum)) {
                k = 0;
                info = -8;
                Mxerbla("Cgedmd", -info);
            }
            if ((scale != zero) && (ssum != zero)) {
                rootsc = sqrt(ssum);
                if (scale >= (ofl / rootsc)) {
                    // Norm of X(:,i) overflows. First, X(:,i)
                    // is scaled by
                    // ( ONE / ROOTSC ) / SCALE = 1/||X(:,i)||_2.
                    // Next, the norm of X(:,i) is stored without
                    // overflow as RWORK(i) = - SCALE * (ROOTSC/M),
                    // the minus sign indicating the 1/M factor.
                    // Scaling is performed without overflow, and
                    // underflow may occur in the smallest entries
                    // of X(:,i). The relative backward and forward
                    // errors are small in the ell_2 norm.
                    Clascl("G", 0, 0, scale, one / rootsc, m, 1, &x[(i - 1) * ldx], ldx, info2);
                    rwork[i - 1] = -scale * (rootsc / castREAL(m));
                } else {
                    // X(:,i) will be scaled to unit 2-norm
                    rwork[i - 1] = scale * rootsc;
                    // LAPACK CALL
                    Clascl("G", 0, 0, rwork[i - 1], one, m, 1, &x[(i - 1) * ldx], ldx, info2);
                    // X(1:M,i) = (ONE/RWORK(i)) * X(1:M,i)   ! INTRINSIC
                }
            } else {
                rwork[i - 1] = zero;
                k++;
            }
        }
        if (k == n) {
            // All columns of X are zero. Return error code -8.
            // (the 8th input variable had an illegal value)
            k = 0;
            info = -8;
            Mxerbla("Cgedmd", -info);
            return;
        }
        for (i = 1; i <= n; i = i + 1) {
            // Now, apply the same scaling to the columns of Y.
            if (rwork[i - 1] > zero) {
                // BLAS CALL
                CRscal(m, one / rwork[i - 1], &y[(i - 1) * ldy], 1);
                // Y(1:M,i) = (ONE/RWORK(i)) * Y(1:M,i)       ! INTRINSIC
            } else if (rwork[i - 1] < zero) {
                // LAPACK CALL
                Clascl("G", 0, 0, -rwork[i - 1], one / castREAL(m), m, 1, &y[(i - 1) * ldy], ldy, info2);
            } else if (abs(y[((iCamax(m, &y[(i - 1) * ldy], 1)) - 1) + (i - 1) * ldy]) != zero) {
                // X(:,i) is zero vector. For consistency,
                // Y(:,i) should also be zero. If Y(:,i) is not
                // zero, then the data might be inconsistent or
                // corrupted. If JOBS == 'C', Y(:,i) is set to
                // zero and a warning flag is raised.
                // The computation continues but the
                // situation will be reported in the output.
                badxy = true;
                // BLAS CALL
                if (Mlsame(jobs, "C")) {
                    CRscal(m, zero, &y[(i - 1) * ldy], 1);
                }
            }
        }
    }
    //
    if (sccoly) {
        // The columns of Y will be normalized.
        // To prevent overflows, the column norms of Y are
        // carefully computed using Classq.
        for (i = 1; i <= n; i = i + 1) {
            // RWORK(i) = RCnrm2( M, Y(1,i), 1 )
            ssum = one;
            scale = zero;
            Classq(m, &y[(i - 1) * ldy], 1, scale, ssum);
            if (Risnan(scale) || Risnan(ssum)) {
                k = 0;
                info = -10;
                Mxerbla("Cgedmd", -info);
            }
            if (scale != zero && (ssum != zero)) {
                rootsc = sqrt(ssum);
                if (scale >= (ofl / rootsc)) {
                    // Norm of Y(:,i) overflows. First, Y(:,i)
                    // is scaled by
                    // ( ONE / ROOTSC ) / SCALE = 1/||Y(:,i)||_2.
                    // Next, the norm of Y(:,i) is stored without
                    // overflow as RWORK(i) = - SCALE * (ROOTSC/M),
                    // the minus sign indicating the 1/M factor.
                    // Scaling is performed without overflow, and
                    // underflow may occur in the smallest entries
                    // of Y(:,i). The relative backward and forward
                    // errors are small in the ell_2 norm.
                    Clascl("G", 0, 0, scale, one / rootsc, m, 1, &y[(i - 1) * ldy], ldy, info2);
                    rwork[i - 1] = -scale * (rootsc / castREAL(m));
                } else {
                    // Y(:,i) will be scaled to unit 2-norm
                    rwork[i - 1] = scale * rootsc;
                    // LAPACK CALL
                    Clascl("G", 0, 0, rwork[i - 1], one, m, 1, &y[(i - 1) * ldy], ldy, info2);
                    // Y(1:M,i) = (ONE/RWORK(i)) * Y(1:M,i)          ! INTRINSIC
                }
            } else {
                rwork[i - 1] = zero;
            }
        }
        for (i = 1; i <= n; i = i + 1) {
            // Now, apply the same scaling to the columns of X.
            if (rwork[i - 1] > zero) {
                // BLAS CALL
                CRscal(m, one / rwork[i - 1], &x[(i - 1) * ldx], 1);
                // X(1:M,i) = (ONE/RWORK(i)) * X(1:M,i)      ! INTRINSIC
            } else if (rwork[i - 1] < zero) {
                // LAPACK CALL
                Clascl("G", 0, 0, -rwork[i - 1], one / castREAL(m), m, 1, &x[(i - 1) * ldx], ldx, info2);
            } else if (abs(x[((iCamax(m, &x[(i - 1) * ldx], 1)) - 1) + (i - 1) * ldx]) != zero) {
                // Y(:,i) is zero vector.  If X(:,i) is not
                // zero, then a warning flag is raised.
                // The computation continues but the
                // situation will be reported in the output.
                badxy = true;
            }
        }
    }
    //
    // <2> SVD of the data snapshot matrix X.
    // =====================================
    // The left singular vectors are stored in the array X.
    // The right singular vectors are in the array W.
    // The array W will later on contain the eigenvectors
    // of a Rayleigh quotient.
    numrnk = n;
    char t_or_n;
    REAL xscl1 = 0.0;
    REAL xscl2 = 0.0;
    if ((whtsvd == 1)) {
        // LAPACK CALL
        Cgesvd("O", "S", m, n, x, ldx, rwork, b, ldb, w, ldw, zwork, lzwork, &rwork[(n + 1) - 1], info1);
        t_or_n = 'C';
    } else if ((whtsvd == 2)) {
        // LAPACK CALL
        Cgesdd("O", m, n, x, ldx, rwork, b, ldb, w, ldw, zwork, lzwork, &rwork[(n + 1) - 1], iwork, info1);
        t_or_n = 'C';
    } else if ((whtsvd == 3)) {
        // LAPACK CALL
        Cgesvdq("H", "P", "N", "R", "R", m, n, x, ldx, rwork, z, ldz, w, ldw, numrnk, iwork, liwork, zwork, lzwork, &rwork[(n + 1) - 1], lrwork - n, info1);
        // LAPACK CALL
        Clacpy("A", m, numrnk, z, ldz, x, ldx);
        t_or_n = 'C';
    } else if ((whtsvd == 4)) {
        // LAPACK CALL
        Cgejsv("F", "U", &jsvopt, "R", "N", "P", m, n, x, ldx, rwork, z, ldz, w, ldw, zwork, lzwork, &rwork[(n + 1) - 1], lrwork - n, iwork, info1);
        // LAPACK CALL
        Clacpy("A", m, n, z, ldz, x, ldx);
        t_or_n = 'N';
        xscl1 = rwork[(n + 1) - 1];
        xscl2 = rwork[(n + 2) - 1];
        if (xscl1 != xscl2) {
            // This is an exceptional situation. If the
            // data matrices are not scaled and the
            // largest singular value of X overflows.
            // In that case Cgejsv can return the SVD
            // in scaled form. The scaling factor can be used
            // to rescale the data (X and Y).
            Clascl("G", 0, 0, xscl1, xscl2, m, n, y, ldy, info2);
        }
    }
    //
    if (info1 > 0) {
        // The SVD selected subroutine did not converge.
        // Return with an error code.
        info = 2;
        return;
    }
    //
    if (rwork[1 - 1] == zero) {
        // The largest computed singular value of (scaled)
        // X is zero. Return error code -8
        // (the 8th input variable had an illegal value).
        k = 0;
        info = -8;
        Mxerbla("Cgedmd", -info);
        return;
    }
    //
    // <3> Determine the numerical rank of the data
    // snapshots matrix X. This depends on the
    // parameters NRNK and TOL.
    //
    if ((nrnk == -1)) {
        k = 1;
        for (i = 2; i <= numrnk; i = i + 1) {
            if ((rwork[i - 1] <= rwork[1 - 1] * tol) || (rwork[i - 1] <= small)) {
                break;
            }
            k++;
        }
    } else if ((nrnk == -2)) {
        k = 1;
        for (i = 1; i <= numrnk - 1; i = i + 1) {
            if ((rwork[(i + 1) - 1] <= rwork[i - 1] * tol) || (rwork[i - 1] <= small)) {
                break;
            }
            k++;
        }
    } else {
        k = 1;
        for (i = 2; i <= nrnk; i = i + 1) {
            if (rwork[i - 1] <= small) {
                break;
            }
            k++;
        }
    }
    // Now, U = X(1:M,1:K) is the SVD/POD basis for the
    // snapshot data in the input matrix X.
    //
    // <4> Compute the Rayleigh quotient S = U^H * A * U.
    // Depending on the requested outputs, the computation
    // is organized to compute additional auxiliary
    // matrices (for the residuals and refinements).
    //
    // In all formulas below, we need V_k*Sigma_k^(-1)
    // where either V_k is in W(1:N,1:K), or V_k^H is in
    // W(1:K,1:N). Here Sigma_k=diag(WORK(1:K)).
    INTEGER j = 0;
    if (Mlsame(&t_or_n, "N")) {
        for (i = 1; i <= k; i = i + 1) {
            // BLAS CALL
            CRscal(n, one / rwork[i - 1], &w[(i - 1) * ldw], 1);
            // W(1:N,i) = (ONE/RWORK(i)) * W(1:N,i)      ! INTRINSIC
        }
    } else {
        // This non-unit stride access is due to the fact
        // that Cgesvd, Cgesvdq and Cgesdd return the
        // adjoint matrix of the right singular vectors.
        // DO i = 1, K
        // CALL CRscal( N, ONE/RWORK(i), W(i,1), LDW )    ! BLAS CALL
        // ! W(i,1:N) = (ONE/RWORK(i)) * W(i,1:N)      ! INTRINSIC
        // END DO
        for (i = 1; i <= k; i = i + 1) {
            rwork[(n + i) - 1] = one / rwork[i - 1];
        }
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                w[(i - 1) + (j - 1) * ldw] = COMPLEX(rwork[(n + i) - 1], zero) * w[(i - 1) + (j - 1) * ldw];
            }
        }
    }
    //
    const COMPLEX zone = COMPLEX(1.0, 0.0);
    const COMPLEX zzero = COMPLEX(0.0, 0.0);
    if (wntref) {
        //
        // Need A*U(:,1:K)=Y*V_k*inv(diag(WORK(1:K)))
        // for computing the refined Ritz vectors
        // (optionally, outside Cgedmd).
        // BLAS CALL
        Cgemm("N", &t_or_n, m, k, n, zone, y, ldy, w, ldw, zzero, z, ldz);
        // Z(1:M,1:K)=MATMUL(Y(1:M,1:N),TRANSPOSE(CONJG(W(1:K,1:N)))) ! INTRINSIC, for T_OR_N=='C'
        // Z(1:M,1:K)=MATMUL(Y(1:M,1:N),W(1:N,1:K))                   ! INTRINSIC, for T_OR_N=='N'
        //
        // At this point Z contains
        // A * U(:,1:K) = Y * V_k * Sigma_k^(-1), and
        // this is needed for computing the residuals.
        // This matrix is  returned in the array B and
        // it can be used to compute refined Ritz vectors.
        // BLAS CALL
        Clacpy("A", m, k, z, ldz, b, ldb);
        // B(1:M,1:K) = Z(1:M,1:K)                  ! INTRINSIC
        //
        // BLAS CALL
        Cgemm("C", "N", k, k, m, zone, x, ldx, z, ldz, zzero, s, lds);
        // S(1:K,1:K) = MATMUL(TRANSPOSE(CONJG(X(1:M,1:K))),Z(1:M,1:K)) ! INTRINSIC
        // At this point S = U^H * A * U is the Rayleigh quotient.
    } else {
        // A * U(:,1:K) is not explicitly needed and the
        // computation is organized differently. The Rayleigh
        // quotient is computed more efficiently.
        // BLAS CALL
        Cgemm("C", "N", k, n, m, zone, x, ldx, y, ldy, zzero, z, ldz);
        // Z(1:K,1:N) = MATMUL( TRANSPOSE(CONJG(X(1:M,1:K))), Y(1:M,1:N) )  ! INTRINSIC
        //
        // BLAS CALL
        Cgemm("N", &t_or_n, k, k, n, zone, z, ldz, w, ldw, zzero, s, lds);
        // S(1:K,1:K) = MATMUL(Z(1:K,1:N),TRANSPOSE(CONJG(W(1:K,1:N)))) ! INTRINSIC, for T_OR_N=='T'
        // S(1:K,1:K) = MATMUL(Z(1:K,1:N),(W(1:N,1:K)))                 ! INTRINSIC, for T_OR_N=='N'
        // At this point S = U^H * A * U is the Rayleigh quotient.
        // If the residuals are requested, save scaled V_k into Z.
        // Recall that V_k or V_k^H is stored in W.
        if (wntres || wntex) {
            if (Mlsame(&t_or_n, "N")) {
                Clacpy("A", n, k, w, ldw, z, ldz);
            } else {
                Clacpy("A", k, n, w, ldw, z, ldz);
            }
        }
    }
    //
    // <5> Compute the Ritz values and (if requested) the
    // right eigenvectors of the Rayleigh quotient.
    //
    // LAPACK CALL
    Cgeev("N", &jobzl, k, s, lds, eigs, w, ldw, w, ldw, zwork, lzwork, &rwork[(n + 1) - 1], info1);
    //
    // W(1:K,1:K) contains the eigenvectors of the Rayleigh
    // quotient.  See the description of Z.
    // Also, see the description of Cgeev.
    if (info1 > 0) {
        // Cgeev failed to compute the eigenvalues and
        // eigenvectors of the Rayleigh quotient.
        info = 3;
        return;
    }
    //
    // <6> Compute the eigenvectors (if requested) and,
    // the residuals (if requested).
    //
    if (wntvec || wntex) {
        if (wntres) {
            if (wntref) {
                // Here, if the refinement is requested, we have
                // A*U(:,1:K) already computed and stored in Z.
                // For the residuals, need Y = A * U(:,1;K) * W.
                // BLAS CALL
                Cgemm("N", "N", m, k, k, zone, z, ldz, w, ldw, zzero, y, ldy);
                // Y(1:M,1:K) = Z(1:M,1:K) * W(1:K,1:K)        ! INTRINSIC
                // This frees Z; Y contains A * U(:,1:K) * W.
            } else {
                // Compute S = V_k * Sigma_k^(-1) * W, where
                // V_k * Sigma_k^(-1) (or its adjoint) is stored in Z
                Cgemm(&t_or_n, "N", n, k, k, zone, z, ldz, w, ldw, zzero, s, lds);
                // Then, compute Z = Y * S =
                // = Y * V_k * Sigma_k^(-1) * W(1:K,1:K) =
                // = A * U(:,1:K) * W(1:K,1:K)
                Cgemm("N", "N", m, k, n, zone, y, ldy, s, lds, zzero, z, ldz);
                // Save a copy of Z into Y and free Z for holding
                // the Ritz vectors.
                Clacpy("A", m, k, z, ldz, y, ldy);
                if (wntex) {
                    Clacpy("A", m, k, z, ldz, b, ldb);
                }
            }
        } else if (wntex) {
            // Compute S = V_k * Sigma_k^(-1) * W, where
            // V_k * Sigma_k^(-1) is stored in Z
            Cgemm(&t_or_n, "N", n, k, k, zone, z, ldz, w, ldw, zzero, s, lds);
            // Then, compute Z = Y * S =
            // = Y * V_k * Sigma_k^(-1) * W(1:K,1:K) =
            // = A * U(:,1:K) * W(1:K,1:K)
            Cgemm("N", "N", m, k, n, zone, y, ldy, s, lds, zzero, b, ldb);
            // The above call replaces the following two calls
            // that were used in the developing-testing phase.
            // CALL Cgemm( 'N', 'N', M, K, N, ZONE, Y, LDY, S, &
            // LDS, ZZERO, Z, LDZ)
            // Save a copy of Z into B and free Z for holding
            // the Ritz vectors.
            // CALL Clacpy( 'A', M, K, Z, LDZ, B, LDB )
        }
        //
        // Compute the Ritz vectors
        // BLAS CALL
        if (wntvec) {
            Cgemm("N", "N", m, k, k, zone, x, ldx, w, ldw, zzero, z, ldz);
        }
        // Z(1:M,1:K) = MATMUL(X(1:M,1:K), W(1:K,1:K))         ! INTRINSIC
        //
        if (wntres) {
            for (i = 1; i <= k; i = i + 1) {
                // BLAS CALL
                Caxpy(m, -eigs[i - 1], &z[(i - 1) * ldz], 1, &y[(i - 1) * ldy], 1);
                // Y(1:M,i) = Y(1:M,i) - EIGS(i) * Z(1:M,i)            ! INTRINSIC
                // BLAS CALL
                res[i - 1] = RCnrm2(m, &y[(i - 1) * ldy], 1);
            }
        }
    }
    //
    if (whtsvd == 4) {
        rwork[(n + 1) - 1] = xscl1;
        rwork[(n + 2) - 1] = xscl2;
    }
    //
    // Successful exit.
    if (!badxy) {
        info = 0;
    } else {
        // A warning on possible data inconsistency.
        // This should be a rare event.
        info = 4;
    }
    // ............................................................
    // ......
}
