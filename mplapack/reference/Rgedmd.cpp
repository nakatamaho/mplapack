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

// Derived from LAPACK routine DGEDMD.
// Original LAPACK authors:
//   Univ. of Tennessee,
//   Univ. of California Berkeley,
//   Univ. of Colorado Denver,
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

void Rgedmd(const char *jobs, const char *jobz, const char *jobr, const char *jobf, INTEGER const whtsvd, INTEGER const m, INTEGER const n, REAL *x, INTEGER const ldx, REAL *y, INTEGER const ldy, INTEGER const nrnk, REAL const tol, INTEGER &k, REAL *reig, REAL *imeig, REAL *z, INTEGER const ldz, REAL *res, REAL *b, INTEGER const ldb, REAL *w, INTEGER const ldw, REAL *s, INTEGER const lds, REAL *work, INTEGER const lwork, INTEGER *iwork, INTEGER const liwork, INTEGER &info) {
    INTEGER ldab = 2;
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
    bool lquery = ((lwork == -1) || (liwork == -1));
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
        info = -18;
    } else if ((wntref || wntex) && (ldb < m)) {
        info = -21;
    } else if (ldw < n) {
        info = -23;
    } else if (lds < n) {
        info = -25;
    }
    //
    INTEGER mlwork = 0;
    INTEGER olwork = 0;
    INTEGER iminwr = 0;
    INTEGER mwrsvd = 0;
    REAL rdummy[2];
    INTEGER info1 = 0;
    INTEGER lwrsvd = 0;
    INTEGER mwrsdd = 0;
    INTEGER lwrsdd = 0;
    INTEGER numrnk = 0;
    REAL rdummy2[2];
    INTEGER mwrsvq = 0;
    INTEGER lwrsvq = 0;
    char jsvopt;
    INTEGER mwrsvj = 0;
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
                work[1 - 1] = 2.0;
                work[2 - 1] = 2.0;
            } else {
                k = 0;
            }
            info = 1;
            return;
        }
        mlwork = max((INTEGER)2, n);
        olwork = max((INTEGER)2, n);
        iminwr = 1;
        if ((whtsvd == 1)) {
            // The following is specified as the minimal
            // length of WORK in the definition of Rgesvd:
            // MWRSVD = MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
            mwrsvd = max((INTEGER)1, 3 * min(m, n) + max(m, n), 5 * min(m, n));
            mlwork = max(mlwork, n + mwrsvd);
            if (lquery) {
                Rgesvd("O", "S", m, n, x, ldx, work, b, ldb, w, ldw, rdummy, -1, info1);
                lwrsvd = max(mwrsvd, castINTEGER(rdummy[1 - 1]));
                olwork = max(olwork, n + lwrsvd);
            }
        } else if ((whtsvd == 2)) {
            // The following is specified as the minimal
            // length of WORK in the definition of Rgesdd:
            // MWRSDD = 3*MIN(M,N)*MIN(M,N) +
            // MAX( MAX(M,N),5*MIN(M,N)*MIN(M,N)+4*MIN(M,N) )
            // IMINWR = 8*MIN(M,N)
            mwrsdd = 3 * min(m, n) * min(m, n) + max(max(m, n), 5 * min(m, n) * min(m, n) + 4 * min(m, n));
            mlwork = max(mlwork, n + mwrsdd);
            iminwr = 8 * min(m, n);
            if (lquery) {
                Rgesdd("O", m, n, x, ldx, work, b, ldb, w, ldw, rdummy, -1, iwork, info1);
                lwrsdd = max(mwrsdd, castINTEGER(rdummy[1 - 1]));
                olwork = max(olwork, n + lwrsdd);
            }
        } else if ((whtsvd == 3)) {
            // LWQP3 = 3*N+1
            // LWORQ = MAX(N, 1)
            // MWRSVD = MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
            // MWRSVQ = N + MAX( LWQP3, MWRSVD, LWORQ ) + MAX(M,2)
            // MLWORK = N +  MWRSVQ
            // IMINWR = M+N-1
            Rgesvdq("H", "P", "N", "R", "R", m, n, x, ldx, work, z, ldz, w, ldw, numrnk, iwork, liwork, rdummy, -1, rdummy2, -1, info1);
            iminwr = iwork[1 - 1];
            mwrsvq = castINTEGER(rdummy[2 - 1]);
            mlwork = max(mlwork, n + mwrsvq + castINTEGER(rdummy2[1 - 1]));
            if (lquery) {
                lwrsvq = max(mwrsvq, castINTEGER(rdummy[1 - 1]));
                olwork = max(olwork, n + lwrsvq + castINTEGER(rdummy2[1 - 1]));
            }
        } else if ((whtsvd == 4)) {
            jsvopt = 'J';
            // MWRSVJ = MAX( 7, 2*M+N, 6*N+2*N*N ) ! for JSVOPT='V'
            mwrsvj = max((INTEGER)7, 2 * m + n, 4 * n + n * n, 2 * n + n * n + 6);
            mlwork = max(mlwork, n + mwrsvj);
            iminwr = max((INTEGER)3, m + 3 * n);
            if (lquery) {
                olwork = max(olwork, n + mwrsvj);
            }
        }
        if (wntvec || wntex || Mlsame(jobz, "F")) {
            jobzl = 'V';
        } else {
            jobzl = 'N';
        }
        // Workspace calculation to the Rgeev call
        if (Mlsame(&jobzl, "V")) {
            mwrkev = max((INTEGER)1, 4 * n);
        } else {
            mwrkev = max((INTEGER)1, 3 * n);
        }
        mlwork = max(mlwork, n + mwrkev);
        if (lquery) {
            Rgeev("N", &jobzl, n, s, lds, reig, imeig, w, ldw, w, ldw, rdummy, -1, info1);
            lwrkev = max(mwrkev, castINTEGER(rdummy[1 - 1]));
            olwork = max(olwork, n + lwrkev);
        }
        //
        if (liwork < iminwr && (!lquery)) {
            info = -29;
        }
        if (lwork < mlwork && (!lquery)) {
            info = -27;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rgedmd", -info);
        return;
    } else if (lquery) {
        // Return minimal and optimal workspace sizes
        iwork[1 - 1] = iminwr;
        work[1 - 1] = mlwork;
        work[2 - 1] = olwork;
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
        // carefully computed using Rlassq.
        k = 0;
        for (i = 1; i <= n; i = i + 1) {
            // WORK(i) = Rnrm2( M, X(1,i), 1 )
            ssum = one;
            scale = zero;
            Rlassq(m, &x[(i - 1) * ldx], 1, scale, ssum);
            if (Risnan(scale) || Risnan(ssum)) {
                k = 0;
                info = -8;
                Mxerbla("Rgedmd", -info);
            }
            if ((scale != zero) && (ssum != zero)) {
                rootsc = sqrt(ssum);
                if (scale >= (ofl / rootsc)) {
                    // Norm of X(:,i) overflows. First, X(:,i)
                    // is scaled by
                    // ( ONE / ROOTSC ) / SCALE = 1/||X(:,i)||_2.
                    // Next, the norm of X(:,i) is stored without
                    // overflow as WORK(i) = - SCALE * (ROOTSC/M),
                    // the minus sign indicating the 1/M factor.
                    // Scaling is performed without overflow, and
                    // underflow may occur in the smallest entries
                    // of X(:,i). The relative backward and forward
                    // errors are small in the ell_2 norm.
                    Rlascl("G", 0, 0, scale, one / rootsc, m, 1, &x[(i - 1) * ldx], m, info2);
                    work[i - 1] = -scale * (rootsc / castREAL(m));
                } else {
                    // X(:,i) will be scaled to unit 2-norm
                    work[i - 1] = scale * rootsc;
                    // LAPACK CALL
                    Rlascl("G", 0, 0, work[i - 1], one, m, 1, &x[(i - 1) * ldx], m, info2);
                    // X(1:M,i) = (ONE/WORK(i)) * X(1:M,i)          ! INTRINSIC
                }
            } else {
                work[i - 1] = zero;
                k++;
            }
        }
        if (k == n) {
            // All columns of X are zero. Return error code -8.
            // (the 8th input variable had an illegal value)
            k = 0;
            info = -8;
            Mxerbla("Rgedmd", -info);
            return;
        }
        for (i = 1; i <= n; i = i + 1) {
            // Now, apply the same scaling to the columns of Y.
            if (work[i - 1] > zero) {
                // BLAS CALL
                Rscal(m, one / work[i - 1], &y[(i - 1) * ldy], 1);
                // Y(1:M,i) = (ONE/WORK(i)) * Y(1:M,i)      ! INTRINSIC
            } else if (work[i - 1] < zero) {
                // LAPACK CALL
                Rlascl("G", 0, 0, -work[i - 1], one / castREAL(m), m, 1, &y[(i - 1) * ldy], m, info2);
            } else if (y[((iRamax(m, &y[(i - 1) * ldy], 1)) - 1) + (i - 1) * ldy] != zero) {
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
                    Rscal(m, zero, &y[(i - 1) * ldy], 1);
                }
            }
        }
    }
    //
    if (sccoly) {
        // The columns of Y will be normalized.
        // To prevent overflows, the column norms of Y are
        // carefully computed using Rlassq.
        for (i = 1; i <= n; i = i + 1) {
            // WORK(i) = Rnrm2( M, Y(1,i), 1 )
            ssum = one;
            scale = zero;
            Rlassq(m, &y[(i - 1) * ldy], 1, scale, ssum);
            if (Risnan(scale) || Risnan(ssum)) {
                k = 0;
                info = -10;
                Mxerbla("Rgedmd", -info);
            }
            if (scale != zero && (ssum != zero)) {
                rootsc = sqrt(ssum);
                if (scale >= (ofl / rootsc)) {
                    // Norm of Y(:,i) overflows. First, Y(:,i)
                    // is scaled by
                    // ( ONE / ROOTSC ) / SCALE = 1/||Y(:,i)||_2.
                    // Next, the norm of Y(:,i) is stored without
                    // overflow as WORK(i) = - SCALE * (ROOTSC/M),
                    // the minus sign indicating the 1/M factor.
                    // Scaling is performed without overflow, and
                    // underflow may occur in the smallest entries
                    // of Y(:,i). The relative backward and forward
                    // errors are small in the ell_2 norm.
                    Rlascl("G", 0, 0, scale, one / rootsc, m, 1, &y[(i - 1) * ldy], m, info2);
                    work[i - 1] = -scale * (rootsc / castREAL(m));
                } else {
                    // X(:,i) will be scaled to unit 2-norm
                    work[i - 1] = scale * rootsc;
                    // LAPACK CALL
                    Rlascl("G", 0, 0, work[i - 1], one, m, 1, &y[(i - 1) * ldy], m, info2);
                    // Y(1:M,i) = (ONE/WORK(i)) * Y(1:M,i)          ! INTRINSIC
                }
            } else {
                work[i - 1] = zero;
            }
        }
        for (i = 1; i <= n; i = i + 1) {
            // Now, apply the same scaling to the columns of X.
            if (work[i - 1] > zero) {
                // BLAS CALL
                Rscal(m, one / work[i - 1], &x[(i - 1) * ldx], 1);
                // X(1:M,i) = (ONE/WORK(i)) * X(1:M,i)      ! INTRINSIC
            } else if (work[i - 1] < zero) {
                // LAPACK CALL
                Rlascl("G", 0, 0, -work[i - 1], one / castREAL(m), m, 1, &x[(i - 1) * ldx], m, info2);
            } else if (x[((iRamax(m, &x[(i - 1) * ldx], 1)) - 1) + (i - 1) * ldx] != zero) {
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
        Rgesvd("O", "S", m, n, x, ldx, work, b, ldb, w, ldw, &work[(n + 1) - 1], lwork - n, info1);
        t_or_n = 'T';
    } else if ((whtsvd == 2)) {
        // LAPACK CALL
        Rgesdd("O", m, n, x, ldx, work, b, ldb, w, ldw, &work[(n + 1) - 1], lwork - n, iwork, info1);
        t_or_n = 'T';
    } else if ((whtsvd == 3)) {
        // LAPACK CALL
        Rgesvdq("H", "P", "N", "R", "R", m, n, x, ldx, work, z, ldz, w, ldw, numrnk, iwork, liwork, &work[(n + max((INTEGER)2, m) + 1) - 1], lwork - n - max((INTEGER)2, m), &work[(n + 1) - 1], max((INTEGER)2, m), info1);
        // LAPACK CALL
        Rlacpy("A", m, numrnk, z, ldz, x, ldx);
        t_or_n = 'T';
    } else if ((whtsvd == 4)) {
        // LAPACK CALL
        Rgejsv("F", "U", &jsvopt, "N", "N", "P", m, n, x, ldx, work, z, ldz, w, ldw, &work[(n + 1) - 1], lwork - n, iwork, info1);
        // LAPACK CALL
        Rlacpy("A", m, n, z, ldz, x, ldx);
        t_or_n = 'N';
        xscl1 = work[(n + 1) - 1];
        xscl2 = work[(n + 2) - 1];
        if (xscl1 != xscl2) {
            // This is an exceptional situation. If the
            // data matrices are not scaled and the
            // largest singular value of X overflows.
            // In that case Rgejsv can return the SVD
            // in scaled form. The scaling factor can be used
            // to rescale the data (X and Y).
            Rlascl("G", 0, 0, xscl1, xscl2, m, n, y, ldy, info2);
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
    if (work[1 - 1] == zero) {
        // The largest computed singular value of (scaled)
        // X is zero. Return error code -8
        // (the 8th input variable had an illegal value).
        k = 0;
        info = -8;
        Mxerbla("Rgedmd", -info);
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
            if ((work[i - 1] <= work[1 - 1] * tol) || (work[i - 1] <= small)) {
                break;
            }
            k++;
        }
    } else if ((nrnk == -2)) {
        k = 1;
        for (i = 1; i <= numrnk - 1; i = i + 1) {
            if ((work[(i + 1) - 1] <= work[i - 1] * tol) || (work[i - 1] <= small)) {
                break;
            }
            k++;
        }
    } else {
        k = 1;
        for (i = 2; i <= nrnk; i = i + 1) {
            if (work[i - 1] <= small) {
                break;
            }
            k++;
        }
    }
    // Now, U = X(1:M,1:K) is the SVD/POD basis for the
    // snapshot data in the input matrix X.
    //
    // <4> Compute the Rayleigh quotient S = U^T * A * U.
    // Depending on the requested outputs, the computation
    // is organized to compute additional auxiliary
    // matrices (for the residuals and refinements).
    //
    // In all formulas below, we need V_k*Sigma_k^(-1)
    // where either V_k is in W(1:N,1:K), or V_k^T is in
    // W(1:K,1:N). Here Sigma_k=diag(WORK(1:K)).
    INTEGER j = 0;
    if (Mlsame(&t_or_n, "N")) {
        for (i = 1; i <= k; i = i + 1) {
            // BLAS CALL
            Rscal(n, one / work[i - 1], &w[(i - 1) * ldw], 1);
            // W(1:N,i) = (ONE/WORK(i)) * W(1:N,i)      ! INTRINSIC
        }
    } else {
        // This non-unit stride access is due to the fact
        // that Rgesvd, Rgesvdq and Rgesdd return the
        // transposed matrix of the right singular vectors.
        // DO i = 1, K
        // CALL Rscal( N, ONE/WORK(i), W(i,1), LDW )    ! BLAS CALL
        // ! W(i,1:N) = (ONE/WORK(i)) * W(i,1:N)      ! INTRINSIC
        // END DO
        for (i = 1; i <= k; i = i + 1) {
            work[(n + i) - 1] = one / work[i - 1];
        }
        for (j = 1; j <= n; j = j + 1) {
            for (i = 1; i <= k; i = i + 1) {
                w[(i - 1) + (j - 1) * ldw] = (work[(n + i) - 1]) * w[(i - 1) + (j - 1) * ldw];
            }
        }
    }
    //
    if (wntref) {
        //
        // Need A*U(:,1:K)=Y*V_k*inv(diag(WORK(1:K)))
        // for computing the refined Ritz vectors
        // (optionally, outside Rgedmd).
        // BLAS CALL
        Rgemm("N", &t_or_n, m, k, n, one, y, ldy, w, ldw, zero, z, ldz);
        // Z(1:M,1:K)=MATMUL(Y(1:M,1:N),TRANSPOSE(W(1:K,1:N)))  ! INTRINSIC, for T_OR_N=='T'
        // Z(1:M,1:K)=MATMUL(Y(1:M,1:N),W(1:N,1:K))             ! INTRINSIC, for T_OR_N=='N'
        //
        // At this point Z contains
        // A * U(:,1:K) = Y * V_k * Sigma_k^(-1), and
        // this is needed for computing the residuals.
        // This matrix is  returned in the array B and
        // it can be used to compute refined Ritz vectors.
        // BLAS CALL
        Rlacpy("A", m, k, z, ldz, b, ldb);
        // B(1:M,1:K) = Z(1:M,1:K)                  ! INTRINSIC
        //
        // BLAS CALL
        Rgemm("T", "N", k, k, m, one, x, ldx, z, ldz, zero, s, lds);
        // S(1:K,1:K) = MATMUL(TANSPOSE(X(1:M,1:K)),Z(1:M,1:K)) ! INTRINSIC
        // At this point S = U^T * A * U is the Rayleigh quotient.
    } else {
        // A * U(:,1:K) is not explicitly needed and the
        // computation is organized differently. The Rayleigh
        // quotient is computed more efficiently.
        // BLAS CALL
        Rgemm("T", "N", k, n, m, one, x, ldx, y, ldy, zero, z, ldz);
        // Z(1:K,1:N) = MATMUL( TRANSPOSE(X(1:M,1:K)), Y(1:M,1:N) )  ! INTRINSIC
        // In the two Rgemm calls here, can use K for LDZ.
        // BLAS CALL
        Rgemm("N", &t_or_n, k, k, n, one, z, ldz, w, ldw, zero, s, lds);
        // S(1:K,1:K) = MATMUL(Z(1:K,1:N),TRANSPOSE(W(1:K,1:N))) ! INTRINSIC, for T_OR_N=='T'
        // S(1:K,1:K) = MATMUL(Z(1:K,1:N),(W(1:N,1:K)))          ! INTRINSIC, for T_OR_N=='N'
        // At this point S = U^T * A * U is the Rayleigh quotient.
        // If the residuals are requested, save scaled V_k into Z.
        // Recall that V_k or V_k^T is stored in W.
        if (wntres || wntex) {
            if (Mlsame(&t_or_n, "N")) {
                Rlacpy("A", n, k, w, ldw, z, ldz);
            } else {
                Rlacpy("A", k, n, w, ldw, z, ldz);
            }
        }
    }
    //
    // <5> Compute the Ritz values and (if requested) the
    // right eigenvectors of the Rayleigh quotient.
    //
    // LAPACK CALL
    Rgeev("N", &jobzl, k, s, lds, reig, imeig, w, ldw, w, ldw, &work[(n + 1) - 1], lwork - n, info1);
    //
    // W(1:K,1:K) contains the eigenvectors of the Rayleigh
    // quotient. Even in the case of complex spectrum, all
    // computation is done in real arithmetic. REIG and
    // IMEIG are the real and the imaginary parts of the
    // eigenvalues, so that the spectrum is given as
    // REIG(:) + sqrt(-1)*IMEIG(:). Complex conjugate pairs
    // are listed at consecutive positions. For such a
    // complex conjugate pair of the eigenvalues, the
    // corresponding eigenvectors are also a complex
    // conjugate pair with the real and imaginary parts
    // stored column-wise in W at the corresponding
    // consecutive column indices. See the description of Z.
    // Also, see the description of Rgeev.
    if (info1 > 0) {
        // Rgeev failed to compute the eigenvalues and
        // eigenvectors of the Rayleigh quotient.
        info = 3;
        return;
    }
    //
    // <6> Compute the eigenvectors (if requested) and,
    // the residuals (if requested).
    //
    REAL ab[2 * 2];
    if (wntvec || wntex) {
        if (wntres) {
            if (wntref) {
                // Here, if the refinement is requested, we have
                // A*U(:,1:K) already computed and stored in Z.
                // For the residuals, need Y = A * U(:,1;K) * W.
                // BLAS CALL
                Rgemm("N", "N", m, k, k, one, z, ldz, w, ldw, zero, y, ldy);
                // Y(1:M,1:K) = Z(1:M,1:K) * W(1:K,1:K)       ! INTRINSIC
                // This frees Z; Y contains A * U(:,1:K) * W.
            } else {
                // Compute S = V_k * Sigma_k^(-1) * W, where
                // V_k * Sigma_k^(-1) is stored in Z
                Rgemm(&t_or_n, "N", n, k, k, one, z, ldz, w, ldw, zero, s, lds);
                // Then, compute Z = Y * S =
                // = Y * V_k * Sigma_k^(-1) * W(1:K,1:K) =
                // = A * U(:,1:K) * W(1:K,1:K)
                Rgemm("N", "N", m, k, n, one, y, ldy, s, lds, zero, z, ldz);
                // Save a copy of Z into Y and free Z for holding
                // the Ritz vectors.
                Rlacpy("A", m, k, z, ldz, y, ldy);
                if (wntex) {
                    Rlacpy("A", m, k, z, ldz, b, ldb);
                }
            }
        } else if (wntex) {
            // Compute S = V_k * Sigma_k^(-1) * W, where
            // V_k * Sigma_k^(-1) is stored in Z
            Rgemm(&t_or_n, "N", n, k, k, one, z, ldz, w, ldw, zero, s, lds);
            // Then, compute Z = Y * S =
            // = Y * V_k * Sigma_k^(-1) * W(1:K,1:K) =
            // = A * U(:,1:K) * W(1:K,1:K)
            Rgemm("N", "N", m, k, n, one, y, ldy, s, lds, zero, b, ldb);
            // The above call replaces the following two calls
            // that were used in the developing-testing phase.
            // CALL Rgemm( 'N', 'N', M, K, N, ONE, Y, LDY, S, &
            // LDS, ZERO, Z, LDZ)
            // Save a copy of Z into B and free Z for holding
            // the Ritz vectors.
            // CALL Rlacpy( 'A', M, K, Z, LDZ, B, LDB )
        }
        //
        // Compute the real form of the Ritz vectors
        // BLAS CALL
        if (wntvec) {
            Rgemm("N", "N", m, k, k, one, x, ldx, w, ldw, zero, z, ldz);
        }
        // Z(1:M,1:K) = MATMUL(X(1:M,1:K), W(1:K,1:K))         ! INTRINSIC
        //
        if (wntres) {
            i = 1;
            while (i <= k) {
                if (imeig[i - 1] == zero) {
                    // have a real eigenvalue with real eigenvector
                    // BLAS CALL
                    Raxpy(m, -reig[i - 1], &z[(i - 1) * ldz], 1, &y[(i - 1) * ldy], 1);
                    // Y(1:M,i) = Y(1:M,i) - REIG(i) * Z(1:M,i)            ! INTRINSIC
                    // BLAS CALL
                    res[i - 1] = Rnrm2(m, &y[(i - 1) * ldy], 1);
                    i++;
                } else {
                    // Have a complex conjugate pair
                    // REIG(i) +- sqrt(-1)*IMEIG(i).
                    // Since all computation is done in real
                    // arithmetic, the formula for the residual
                    // is recast for real representation of the
                    // complex conjugate eigenpair. See the
                    // description of RES.
                    ab[0] = reig[i - 1];
                    ab[(2 - 1)] = -imeig[i - 1];
                    ab[(2 - 1) * ldab] = imeig[i - 1];
                    ab[(2 - 1) + (2 - 1) * ldab] = reig[i - 1];
                    // BLAS CALL
                    Rgemm("N", "N", m, 2, 2, -one, &z[(i - 1) * ldz], ldz, ab, 2, one, &y[(i - 1) * ldy], ldy);
                    // Y(1:M,i:i+1) = Y(1:M,i:i+1) - Z(1:M,i:i+1) * AB   ! INTRINSIC
                    // LAPACK CALL
                    res[i - 1] = Rlange("F", m, 2, &y[(i - 1) * ldy], ldy, &work[(n + 1) - 1]);
                    res[(i + 1) - 1] = res[i - 1];
                    i += 2;
                }
            }
        }
    }
    //
    if (whtsvd == 4) {
        work[(n + 1) - 1] = xscl1;
        work[(n + 2) - 1] = xscl2;
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
