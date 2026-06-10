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

// Derived from LAPACK routine DCHKDMD.
// Original LAPACK authors:
//   Univ. of Tennessee,
//   Univ. of California Berkeley,
//   Univ. of Colorado Denver,
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_eig.h>
#include <memory>

template <class T>
static inline void copy_matrix_block(T *dst, INTEGER ldd, INTEGER dst_col, const T *src, INTEGER lds, INTEGER src_col, INTEGER rows, INTEGER cols) {
    for (INTEGER j = 0; j < cols; j = j + 1) {
        for (INTEGER i = 0; i < rows; i = i + 1) {
            dst[i + (dst_col - 1 + j) * ldd] = src[i + (src_col - 1 + j) * lds];
        }
    }
}


// This is a test program for checking the implementations of
// the implementations of the following subroutines
//
// Rgedmd  for computation of the
// Dynamic Mode Decomposition (DMD)
// Rgedmdq for computation of a
// QR factorization based compressed DMD
//
// Developed and supported by:
// ===========================
// Developed and coded by Zlatko Drmac, Faculty of Science,
// University of Zagreb;  drmac@math.hr
// In cooperation with
// AIMdyn Inc., Santa Barbara, CA.
// ========================================================
// How to run the code (compiler, link info)
// ========================================================
// Compile as FORTRAN 90 (or later) and link with BLAS and
// LAPACK libraries.
// NOTE: The code is developed and tested on top of the
// Intel MKL library (versions 2022.0.3 and 2022.2.0),
// using the Intel Fortran compiler.
//
// For developers of the C++ implementation
// ========================================================
// See the LAPACK++ and Template Numerical Toolkit (TNT)
//
// Note on a development of the GPU HP implementation
// ========================================================
// Work in progress. See CUDA, MAGMA, SLATE.
// NOTE: The four SVD subroutines used in this code are
// included as a part of R&D and for the completeness.
// This was also an opportunity to test those SVD codes.
// If the scaling option is used all four are essentially
// equally good. For implementations on HP platforms,
// one can use whichever SVD is available.
// ...  .........................................................
// NOTE:
// When using the Intel MKL 2022.0.3 the subroutine xGESVDQ
// (optionally used in xGEDMD) may cause access violation
// error for x = S, D, C, Z, but only if called with the
// work space query. (At least in our Windows 10 MSVS 2019.)
// The problem can be mitigated by downloading the source
// code of xGESVDQ from the LAPACK repository and use it
// localy instead of the one in the MKL. This seems to
// indicate that the problem is indeed in the MKL.
// This problem did not appear whith Intel MKL 2022.2.0.
//
// NOTE:
// xGESDD seems to have a problem with workspace. In some
// cases the length of the optimal workspace is returned
// smaller than the minimal workspace, as specified in the
// code. As a precaution, all optimal workspaces are
// set as MAX(minimal, optimal).
// Latest implementations of complex xGESDD have different
// length of the real worksapce. We use max value over
// two versions.
// ............................................................
// ............................................................
//
void program_dmd_test(int argc, char const *argv[]) {
    common cmn(argc, argv);
    common_read read(cmn);
    common_write write(cmn);
    INTEGER ldab = 2;
    //
    // use iso_fortran_env, only: real64
    // integer, parameter :: WP = real64
    //
    // ............................................................
    // ............................................................
    // ............................................................
    //
    // ..... external subroutines (BLAS and LAPACK)
    // .....external subroutines DMD package, part 1
    // subroutines under test
    //
    // ..... external functions (BLAS and LAPACK)
    //
    // ............................................................
    //
    // The test is always in pairs : ( Rgedmd and Rgedmdq )
    // because the test includes comparing the results (in pairs).
    // .....................................................................................
    bool test_qrdmd = true;
    // Since the QR factorizations based algorithm is designed for
    // single trajectory data, only single trajectory tests will
    // be performed with xGEDMDQ.
    fem::str<1> wantq = "Q";
    fem::str<1> wantr = "R";
    // .................................................................................
    //
    // machine precision DP
    REAL eps = Rlamch("P");
    //
    // Global counters of failures of some particular tests
    INTEGER nfail = 0;
    INTEGER nfail_rez = 0;
    INTEGER nfail_rezq = 0;
    INTEGER nfail_z_xv = 0;
    INTEGER nfail_f_qr = 0;
    INTEGER nfail_au = 0;
    INTEGER kdiff = 0;
    INTEGER nfail_svdiff = 0;
    INTEGER nfail_total = 0;
    INTEGER nfailq_total = 0;
    //
    INTEGER lloop = 0;
    INTEGER m = 0;
    INTEGER n = 0;
    INTEGER iseed[4];
    INTEGER lda = 0;
    INTEGER ldf = 0;
    INTEGER ldx = 0;
    INTEGER ldy = 0;
    INTEGER ldw = 0;
    INTEGER ldz = 0;
    INTEGER ldau = 0;
    INTEGER lds = 0;
    const REAL zero = 0.0;
    REAL tmp_zxw = 0.0;
    REAL tmp_au = 0.0;
    REAL tmp_rez = 0.0;
    REAL tmp_rezq = 0.0;
    REAL svdiff = 0.0;
    REAL tmp_ex = 0.0;
    std::unique_ptr<REAL[]> a_storage;
    REAL *a = nullptr;
    std::unique_ptr<REAL[]> ac_storage;
    REAL *ac = nullptr;
    std::unique_ptr<REAL[]> da_storage;
    REAL *da = nullptr;
    std::unique_ptr<REAL[]> dl_storage;
    REAL *dl = nullptr;
    std::unique_ptr<REAL[]> f_storage;
    REAL *f = nullptr;
    std::unique_ptr<REAL[]> f1_storage;
    REAL *f1 = nullptr;
    std::unique_ptr<REAL[]> f2_storage;
    REAL *f2 = nullptr;
    std::unique_ptr<REAL[]> x_storage;
    REAL *x = nullptr;
    std::unique_ptr<REAL[]> x0_storage;
    REAL *x0 = nullptr;
    std::unique_ptr<REAL[]> singvx_storage;
    REAL *singvx = nullptr;
    std::unique_ptr<REAL[]> singvqx_storage;
    REAL *singvqx = nullptr;
    std::unique_ptr<REAL[]> y_storage;
    REAL *y = nullptr;
    std::unique_ptr<REAL[]> y0_storage;
    REAL *y0 = nullptr;
    std::unique_ptr<REAL[]> y1_storage;
    REAL *y1 = nullptr;
    std::unique_ptr<REAL[]> z_storage;
    REAL *z = nullptr;
    std::unique_ptr<REAL[]> z1_storage;
    REAL *z1 = nullptr;
    std::unique_ptr<REAL[]> res_storage;
    REAL *res = nullptr;
    std::unique_ptr<REAL[]> res1_storage;
    REAL *res1 = nullptr;
    std::unique_ptr<REAL[]> resex_storage;
    REAL *resex = nullptr;
    std::unique_ptr<REAL[]> reig_storage;
    REAL *reig = nullptr;
    std::unique_ptr<REAL[]> ieig_storage;
    REAL *ieig = nullptr;
    std::unique_ptr<REAL[]> reigq_storage;
    REAL *reigq = nullptr;
    std::unique_ptr<REAL[]> ieigq_storage;
    REAL *ieigq = nullptr;
    std::unique_ptr<REAL[]> reiga_storage;
    REAL *reiga = nullptr;
    std::unique_ptr<REAL[]> ieiga_storage;
    REAL *ieiga = nullptr;
    std::unique_ptr<REAL[]> va_storage;
    REAL *va = nullptr;
    std::unique_ptr<REAL[]> lambda_storage;
    REAL *lambda = nullptr;
    std::unique_ptr<REAL[]> lambdaq_storage;
    REAL *lambdaq = nullptr;
    std::unique_ptr<REAL[]> eiga_storage;
    REAL *eiga = nullptr;
    std::unique_ptr<REAL[]> w_storage;
    REAL *w = nullptr;
    std::unique_ptr<REAL[]> au_storage;
    REAL *au = nullptr;
    std::unique_ptr<REAL[]> s_storage;
    REAL *s = nullptr;
    REAL tol = 0.0;
    REAL tol2 = 0.0;
    INTEGER k_traj = 0;
    REAL cond = 0.0;
    REAL dmax = 0.0;
    fem::str<1> rsign;
    fem::str<1> grade;
    INTEGER model = 0;
    REAL condl = 0.0;
    INTEGER moder = 0;
    REAL condr = 0.0;
    fem::str<1> pivtng;
    INTEGER mode = 0;
    std::unique_ptr<INTEGER[]> iwork_storage;
    INTEGER *iwork = nullptr;
    std::unique_ptr<REAL[]> dr_storage;
    REAL *dr = nullptr;
    const REAL one = 1.0;
    INTEGER info = 0;
    INTEGER lwork = 0;
    std::unique_ptr<REAL[]> work_storage;
    REAL *work = nullptr;
    REAL tmp = 0.0;
    INTEGER i = 0;
    REAL wdummy[2];
    REAL anorm = 0.0;
    REAL xnorm = 0.0;
    REAL ynorm = 0.0;
    INTEGER ijobz = 0;
    fem::str<1> jobz;
    fem::str<1> resids;
    INTEGER ijobref = 0;
    fem::str<1> jobref;
    INTEGER iscale = 0;
    fem::str<1> scale;
    INTEGER inrnk = 0;
    INTEGER nrnk = 0;
    INTEGER iwhtsvd = 0;
    INTEGER whtsvd = 0;
    INTEGER lwminopt = 0;
    INTEGER k = 0;
    INTEGER idummy[2];
    INTEGER liwork = 0;
    REAL ab[2 * 2];
    INTEGER rjobdata[8];
    INTEGER kq = 0;
    INTEGER j = 0;
    REAL tmp_fqr = 0.0;
    for (lloop = 1; lloop <= 4; lloop = lloop + 1) {
        //
        write(6, star), "L Loop Index = ", lloop;
        //
        // Set the dimensions of the problem ...
        write(6, star), "M = ";
        read(5, star), m;
        write(6, star), m;
        // ... and the number of snapshots.
        write(6, star), "N = ";
        read(5, star), n;
        write(6, star), n;
        //
        // ... Test the dimensions
        if ((min(m, n) == 0) || (m < n)) {
            write(6, star), "Bad dimensions. Required: M >= N > 0.";
            FEM_STOP(0);
        }
        // .............
        // The seed inside the LLOOP so that each pass can be reproduced easily.
        //
        iseed[1 - 1] = 4;
        iseed[2 - 1] = 3;
        iseed[3 - 1] = 2;
        iseed[4 - 1] = 1;
        //
        lda = m;
        ldf = m;
        ldx = max(m, n + 1);
        ldy = max(m, n + 1);
        ldw = n;
        ldz = m;
        ldau = max(m, n + 1);
        lds = n;
        //
        tmp_zxw = zero;
        tmp_au = zero;
        tmp_rez = zero;
        tmp_rezq = zero;
        svdiff = zero;
        tmp_ex = zero;
        //
        // Test the subroutines on real data snapshots. All
        // computation is done in real arithmetic, even when
        // Koopman eigenvalues and modes are real.
        //
        // Allocate memory space
        a_storage = std::make_unique<REAL[]>(max((INTEGER)1, lda * m));
        a = a_storage.get();
        ac_storage = std::make_unique<REAL[]>(max((INTEGER)1, lda * m));
        ac = ac_storage.get();
        da_storage = std::make_unique<REAL[]>(max((INTEGER)1, m));
        da = da_storage.get();
        dl_storage = std::make_unique<REAL[]>(max((INTEGER)1, m));
        dl = dl_storage.get();
        f_storage = std::make_unique<REAL[]>(max((INTEGER)1, ldf * (n + 1)));
        f = f_storage.get();
        f1_storage = std::make_unique<REAL[]>(max((INTEGER)1, ldf * (n + 1)));
        f1 = f1_storage.get();
        f2_storage = std::make_unique<REAL[]>(max((INTEGER)1, ldf * (n + 1)));
        f2 = f2_storage.get();
        x_storage = std::make_unique<REAL[]>(max((INTEGER)1, ldx * n));
        x = x_storage.get();
        x0_storage = std::make_unique<REAL[]>(max((INTEGER)1, ldx * n));
        x0 = x0_storage.get();
        singvx_storage = std::make_unique<REAL[]>(max((INTEGER)1, n));
        singvx = singvx_storage.get();
        singvqx_storage = std::make_unique<REAL[]>(max((INTEGER)1, n));
        singvqx = singvqx_storage.get();
        y_storage = std::make_unique<REAL[]>(max((INTEGER)1, ldy * (n + 1)));
        y = y_storage.get();
        y0_storage = std::make_unique<REAL[]>(max((INTEGER)1, ldy * (n + 1)));
        y0 = y0_storage.get();
        y1_storage = std::make_unique<REAL[]>(max((INTEGER)1, m * (n + 1)));
        y1 = y1_storage.get();
        z_storage = std::make_unique<REAL[]>(max((INTEGER)1, ldz * n));
        z = z_storage.get();
        z1_storage = std::make_unique<REAL[]>(max((INTEGER)1, ldz * n));
        z1 = z1_storage.get();
        res_storage = std::make_unique<REAL[]>(max((INTEGER)1, n));
        res = res_storage.get();
        res1_storage = std::make_unique<REAL[]>(max((INTEGER)1, n));
        res1 = res1_storage.get();
        resex_storage = std::make_unique<REAL[]>(max((INTEGER)1, n));
        resex = resex_storage.get();
        reig_storage = std::make_unique<REAL[]>(max((INTEGER)1, n));
        reig = reig_storage.get();
        ieig_storage = std::make_unique<REAL[]>(max((INTEGER)1, n));
        ieig = ieig_storage.get();
        reigq_storage = std::make_unique<REAL[]>(max((INTEGER)1, n));
        reigq = reigq_storage.get();
        ieigq_storage = std::make_unique<REAL[]>(max((INTEGER)1, n));
        ieigq = ieigq_storage.get();
        reiga_storage = std::make_unique<REAL[]>(max((INTEGER)1, m));
        reiga = reiga_storage.get();
        ieiga_storage = std::make_unique<REAL[]>(max((INTEGER)1, m));
        ieiga = ieiga_storage.get();
        va_storage = std::make_unique<REAL[]>(max((INTEGER)1, lda * m));
        va = va_storage.get();
        lambda_storage = std::make_unique<REAL[]>(max((INTEGER)1, n * 2));
        lambda = lambda_storage.get();
        lambdaq_storage = std::make_unique<REAL[]>(max((INTEGER)1, n * 2));
        lambdaq = lambdaq_storage.get();
        eiga_storage = std::make_unique<REAL[]>(max((INTEGER)1, m * 2));
        eiga = eiga_storage.get();
        w_storage = std::make_unique<REAL[]>(max((INTEGER)1, ldw * n));
        w = w_storage.get();
        au_storage = std::make_unique<REAL[]>(max((INTEGER)1, ldau * n));
        au = au_storage.get();
        s_storage = std::make_unique<REAL[]>(max((INTEGER)1, n * n));
        s = s_storage.get();
        //
        tol = m * eps;
        // This mimics O(M*N)*EPS bound for accumulated roundoff error.
        // The factor 10 is somewhat arbitrary.
        tol2 = 10 * m * n * eps;
        //
        // .............
        //
        for (k_traj = 1; k_traj <= 2; k_traj = k_traj + 1) {
            // Number of intial conditions in the simulation/trajectories (1 or 2)
            //
            cond = 100000000.0;
            dmax = 100.0;
            rsign = "F";
            grade = "N";
            model = 6;
            condl = 100.0;
            moder = 6;
            condr = 100.0;
            pivtng = "N";
            //
            // Loop over all parameter MODE values for Clatmr (+1,..,+6)
            for (mode = 1; mode <= 6; mode = mode + 1) {
                //
                iwork_storage = std::make_unique<INTEGER[]>(max((INTEGER)1, (2 * m)));
                iwork = iwork_storage.get();
                dr_storage = std::make_unique<REAL[]>(max((INTEGER)1, n));
                dr = dr_storage.get();
                Rlatmr(m, m, "S", iseed, "N", da, mode, cond, dmax, rsign, grade, dl, model, condl, dr, moder, condr, pivtng, iwork, m, m, zero, -one, "N", a, lda, &iwork[(m + 1) - 1], info);
                //
                lwork = 4 * m + 1;
                work_storage = std::make_unique<REAL[]>(max((INTEGER)1, lwork));
                work = work_storage.get();
                copy_matrix_block(ac, lda, 1, a, lda, 1, m, m);
                // LAPACK CALL
                Rgeev("N", "V", m, ac, m, reiga, ieiga, va, m, va, m, work, lwork, info);
                tmp = zero;
                for (i = 1; i <= m; i = i + 1) {
                    eiga[(i - 1)] = reiga[i - 1];
                    eiga[(i - 1) + (2 - 1) * m] = ieiga[i - 1];
                    tmp = max(tmp, sqrt(pow2(reiga[i - 1]) + pow2(ieiga[i - 1])));
                }
                //
                // Scale A to have the desirable spectral radius.
                Rlascl("G", 0, 0, tmp, one, m, m, a, m, info);
                Rlascl("G", 0, 0, tmp, one, m, 2, eiga, m, info);
                //
                // Compute the norm of A
                anorm = Rlange("F", n, n, a, m, wdummy);
                //
                if (k_traj == 2) {
                    // generate data with two inital conditions
                    Rlarnv(2, iseed, m, &f1[0]);
                    for (INTEGER i_ = 1; i_ <= m; i_++) {
                        f1[(i_ - 1)] = 0.0000000001 * f1[(i_ - 1)];
                    }
                    for (i = 1; i <= n / 2; i = i + 1) {
                        Rgemv("N", m, m, one, a, m, &f1[(i - 1) * ldf], 1, zero, &f1[((i + 1) - 1) * ldf], 1);
                    }
                    copy_matrix_block(x0, ldx, 1, f1, ldf, 1, m, n / 2);
                    copy_matrix_block(y0, ldy, 1, f1, ldf, 2, m, n / 2);
                    //
                    Rlarnv(2, iseed, m, &f1[0]);
                    for (i = 1; i <= n - n / 2; i = i + 1) {
                        Rgemv("N", m, m, one, a, m, &f1[(i - 1) * ldf], 1, zero, &f1[((i + 1) - 1) * ldf], 1);
                    }
                    copy_matrix_block(x0, ldx, n / 2 + 1, f1, ldf, 1, m, n - n / 2);
                    copy_matrix_block(y0, ldy, n / 2 + 1, f1, ldf, 2, m, n - n / 2);
                } else {
                    Rlarnv(2, iseed, m, &f[0]);
                    for (i = 1; i <= n; i = i + 1) {
                        Rgemv("N", m, m, one, a, m, &f[(i - 1) * ldf], 1, zero, &f[((i + 1) - 1) * ldf], 1);
                    }
                    copy_matrix_block(x0, ldx, 1, f, ldf, 1, m, n);
                    copy_matrix_block(y0, ldy, 1, f, ldf, 2, m, n);
                }
                //
                xnorm = Rlange("F", m, n, x0, ldx, wdummy);
                ynorm = Rlange("F", m, n, y0, ldx, wdummy);
                // ............................................................
                //
                for (ijobz = 1; ijobz <= 4; ijobz = ijobz + 1) {
                    //
                    if ((ijobz == 1)) {
                        // Ritz vectors will be computed
                        jobz = "V";
                        // Residuals will be computed
                        resids = "R";
                    } else if ((ijobz == 2)) {
                        jobz = "V";
                        resids = "N";
                    } else if ((ijobz == 3)) {
                        // Ritz vectors in factored form
                        jobz = "F";
                        resids = "N";
                    } else if ((ijobz == 4)) {
                        jobz = "N";
                        resids = "N";
                    }
                    //
                    for (ijobref = 1; ijobref <= 3; ijobref = ijobref + 1) {
                        //
                        if ((ijobref == 1)) {
                            // Data for refined Ritz vectors
                            jobref = "R";
                        } else if ((ijobref == 2)) {
                            // Exact DMD vectors
                            jobref = "E";
                        } else if ((ijobref == 3)) {
                            jobref = "N";
                        }
                        //
                        for (iscale = 1; iscale <= 4; iscale = iscale + 1) {
                            //
                            if ((iscale == 1)) {
                                // X data normalized
                                scale = "S";
                            } else if ((iscale == 2)) {
                                // X normalized, consist. check
                                scale = "C";
                            } else if ((iscale == 3)) {
                                // Y data normalized
                                scale = "Y";
                            } else if ((iscale == 4)) {
                                scale = "N";
                            }
                            //
                            for (inrnk = -1; inrnk >= -2; inrnk = inrnk - 1) {
                                // Two truncation strategies. The "-2" case for R&D
                                // purposes only - it uses possibly low accuracy small
                                // singular values, in which case the formulas used in
                                // the DMD are highly sensitive.
                                nrnk = inrnk;
                                //
                                for (iwhtsvd = 1; iwhtsvd <= 4; iwhtsvd = iwhtsvd + 1) {
                                    // Check all four options to compute the POD basis
                                    // via the SVD.
                                    whtsvd = iwhtsvd;
                                    //
                                    for (lwminopt = 1; lwminopt <= 2; lwminopt = lwminopt + 1) {
                                        // Workspace query for the minimal (1) and for the optimal
                                        // (2) workspace lengths determined by workspace query.
                                        //
                                        copy_matrix_block(x, ldx, 1, x0, ldx, 1, m, n);
                                        copy_matrix_block(y, ldy, 1, y0, ldy, 1, m, n);
                                        //
                                        // Rgedmd: Workspace query and workspace allocation
                                        Rgedmd(scale.elems, jobz.elems, resids.elems, jobref.elems, whtsvd, m, n, x, ldx, y, ldy, nrnk, tol, k, reig, ieig, z, ldz, res, au, ldau, w, ldw, s, lds, wdummy, -1, idummy, -1, info);
                                        //
                                        liwork = idummy[1 - 1];
                                        iwork_storage = std::make_unique<INTEGER[]>(max((INTEGER)1, liwork));
                                        iwork = iwork_storage.get();
                                        lwork = castINTEGER(wdummy[lwminopt - 1]);
                                        work_storage = std::make_unique<REAL[]>(max((INTEGER)1, lwork));
                                        work = work_storage.get();
                                        //
                                        // Rgedmd test: CALL Rgedmd
                                        Rgedmd(scale.elems, jobz.elems, resids.elems, jobref.elems, whtsvd, m, n, x, ldx, y, ldy, nrnk, tol, k, reig, ieig, z, ldz, res, au, ldau, w, ldw, s, lds, work, lwork, iwork, liwork, info);
                                        //
                                        for (INTEGER i_ = 1; i_ <= n; i_++) {
                                            singvx[i_ - 1] = work[i_ - 1];
                                        }
                                        //
                                        // ...... Rgedmd check point
                                        if (Mlsame(jobz.elems, "V")) {
                                            // Check that Z = X*W, on return from Rgedmd
                                            // This checks that the returned aigenvectors in Z are
                                            // the product of the SVD'POD basis returned in X
                                            // and the eigenvectors of the rayleigh quotient
                                            // returned in W
                                            Rgemm("N", "N", m, k, k, one, x, ldx, w, ldw, zero, z1, ldz);
                                            tmp = zero;
                                            for (i = 1; i <= k; i = i + 1) {
                                                Raxpy(m, -one, &z[(i - 1) * ldz], 1, &z1[(i - 1) * ldz], 1);
                                                tmp = max(tmp, Rnrm2(m, &z1[(i - 1) * ldz], 1));
                                            }
                                            tmp_zxw = max(tmp_zxw, tmp);
                                            //
                                            if (tmp_zxw > 10 * m * eps) {
                                                nfail_z_xv++;
                                                write(6, star), ":( .................Rgedmd FAILED!", "Check the code for implementation errors.";
                                                write(6, star), "The input parameters were ", scale, jobz, resids, jobref, whtsvd, m, n, ldx, ldy, nrnk, tol;
                                            }
                                            //
                                        }
                                        //
                                        // ...... Rgedmd check point
                                        if (Mlsame(jobref.elems, "R")) {
                                            // The matrix A*U is returned for computing refined Ritz vectors.
                                            // Check that A*U is computed correctly using the formula
                                            // A*U = Y * V * inv(SIGMA). This depends on the
                                            // accuracy in the computed singular values and vectors of X.
                                            // See the paper for an error analysis.
                                            // Note that the left singular vectors of the input matrix X
                                            // are returned in the array X.
                                            Rgemm("N", "N", m, k, m, one, a, lda, x, ldx, zero, z1, ldz);
                                            tmp = zero;
                                            for (i = 1; i <= k; i = i + 1) {
                                                Raxpy(m, -one, &au[(i - 1) * ldau], 1, &z1[(i - 1) * ldz], 1);
                                                tmp = max(tmp, Rnrm2(m, &z1[(i - 1) * ldz], 1) * singvx[k - 1] / (anorm * singvx[1 - 1]));
                                            }
                                            tmp_au = max(tmp_au, tmp);
                                            //
                                            if (tmp > tol2) {
                                                nfail_au++;
                                                write(6, star), ":( .................Rgedmd FAILED!", "Check the code for implementation errors.";
                                                write(6, star), "The input parameters were ", scale, jobz, resids, jobref, whtsvd, m, n, ldx, ldy, nrnk, tol;
                                            }
                                            //
                                        } else if (Mlsame(jobref.elems, "E")) {
                                            // The unscaled vectors of the Exact DMD are computed.
                                            // This option is included for the sake of completeness,
                                            // for users who prefer the Exact DMD vectors. The
                                            // returned vectors are in the real form, in the same way
                                            // as the Ritz vectors. Here we just save the vectors
                                            // and test them separately using a Matlab script.
                                            //
                                            Rgemm("N", "N", m, k, m, one, a, lda, au, ldau, zero, y1, m);
                                            i = 1;
                                            while (i <= k) {
                                                if (ieig[i - 1] == zero) {
                                                    // have a real eigenvalue with real eigenvector
                                                    Raxpy(m, -reig[i - 1], &au[(i - 1) * ldau], 1, &y1[(i - 1) * m], 1);
                                                    resex[i - 1] = Rnrm2(m, &y1[(i - 1) * m], 1) / Rnrm2(m, &au[(i - 1) * ldau], 1);
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
                                                    ab[(2 - 1)] = -ieig[i - 1];
                                                    ab[(2 - 1) * ldab] = ieig[i - 1];
                                                    ab[(2 - 1) + (2 - 1) * ldab] = reig[i - 1];
                                                    Rgemm("N", "N", m, 2, 2, -one, &au[(i - 1) * ldau], m, ab, 2, one, &y1[(i - 1) * m], m);
                                                    resex[i - 1] = Rlange("F", m, 2, &y1[(i - 1) * m], m, work) / Rlange("F", m, 2, &au[(i - 1) * ldau], m, work);
                                                    resex[(i + 1) - 1] = resex[i - 1];
                                                    i += 2;
                                                }
                                            }
                                            //
                                        }
                                        //
                                        // ...... Rgedmd check point
                                        if (Mlsame(resids.elems, "R")) {
                                            // Compare the residuals returned by Rgedmd with the
                                            // explicitly computed residuals using the matrix A.
                                            // Compute explicitly Y1 = A*Z
                                            Rgemm("N", "N", m, k, m, one, a, lda, z, ldz, zero, y1, m);
                                            // ... and then A*Z(:,i) - LAMBDA(i)*Z(:,i), using the real forms
                                            // of the invariant subspaces that correspond to complex conjugate
                                            // pairs of eigencalues. (See the description of Z in Rgedmd,)
                                            i = 1;
                                            while (i <= k) {
                                                if (ieig[i - 1] == zero) {
                                                    // have a real eigenvalue with real eigenvector
                                                    Raxpy(m, -reig[i - 1], &z[(i - 1) * ldz], 1, &y1[(i - 1) * m], 1);
                                                    res1[i - 1] = Rnrm2(m, &y1[(i - 1) * m], 1);
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
                                                    ab[(2 - 1)] = -ieig[i - 1];
                                                    ab[(2 - 1) * ldab] = ieig[i - 1];
                                                    ab[(2 - 1) + (2 - 1) * ldab] = reig[i - 1];
                                                    Rgemm("N", "N", m, 2, 2, -one, &z[(i - 1) * ldz], m, ab, 2, one, &y1[(i - 1) * m], m);
                                                    res1[i - 1] = Rlange("F", m, 2, &y1[(i - 1) * m], m, work);
                                                    res1[(i + 1) - 1] = res1[i - 1];
                                                    i += 2;
                                                }
                                            }
                                            tmp = zero;
                                            for (i = 1; i <= k; i = i + 1) {
                                                tmp = max(tmp, abs(res[i - 1] - res1[i - 1]) * singvx[k - 1] / (anorm * singvx[1 - 1]));
                                            }
                                            tmp_rez = max(tmp_rez, tmp);
                                            //
                                            if (tmp > tol2) {
                                                nfail_rez++;
                                                write(6, star), ":( ..................Rgedmd FAILED!", "Check the code for implementation errors.";
                                                write(6, star), "The input parameters were ", scale, jobz, resids, jobref, whtsvd, m, n, ldx, ldy, nrnk, tol;
                                            }
                                            //
                                            if (Mlsame(jobref.elems, "E")) {
                                                tmp = zero;
                                                for (i = 1; i <= k; i = i + 1) {
                                                    tmp = max(tmp, abs(res1[i - 1] - resex[i - 1]) / (res1[i - 1] + resex[i - 1]));
                                                }
                                                tmp_ex = max(tmp_ex, tmp);
                                            }
                                            //
                                        }
                                        //
                                        // ..... store the results for inspection
                                        for (i = 1; i <= k; i = i + 1) {
                                            lambda[(i - 1)] = reig[i - 1];
                                            lambda[(i - 1) + (2 - 1) * n] = ieig[i - 1];
                                        }
                                        //
                                        // ======================================================================
                                        // Now test the Rgedmdq
                                        // ======================================================================
                                        if (test_qrdmd && (k_traj == 1)) {
                                            rjobdata[2 - 1] = 1;
                                            copy_matrix_block(f1, ldf, 1, f, ldf, 1, m, n + 1);
                                            //
                                            // Rgedmdq test: Workspace query and workspace allocation
                                            Rgedmdq(scale.elems, jobz.elems, resids.elems, wantq.elems, wantr.elems, jobref.elems, whtsvd, m, n + 1, f1, ldf, x, ldx, y, ldy, nrnk, tol, kq, reigq, ieigq, z, ldz, res, au, ldau, w, ldw, s, lds, wdummy, -1, idummy, -1, info);
                                            liwork = idummy[1 - 1];
                                            iwork_storage = std::make_unique<INTEGER[]>(max((INTEGER)1, liwork));
                                            iwork = iwork_storage.get();
                                            lwork = castINTEGER(wdummy[lwminopt - 1]);
                                            work_storage = std::make_unique<REAL[]>(max((INTEGER)1, lwork));
                                            work = work_storage.get();
                                            // Rgedmdq test: CALL Rgedmdq
                                            Rgedmdq(scale.elems, jobz.elems, resids.elems, wantq.elems, wantr.elems, jobref.elems, whtsvd, m, n + 1, f1, ldf, x, ldx, y, ldy, nrnk, tol, kq, reigq, ieigq, z, ldz, res, au, ldau, w, ldw, s, lds, work, lwork, iwork, liwork, info);
                                            //
                                            for (INTEGER i_ = 1; i_ <= kq; i_++) {
                                                singvqx[i_ - 1] = work[min(m, n + 1) + i_ - 1];
                                            }
                                            //
                                            // ..... Rgedmdq check point
                                            if (kq != k) {
                                                kdiff++;
                                            }
                                            //
                                            tmp = zero;
                                            for (i = 1; i <= min(k, kq); i = i + 1) {
                                                tmp = max(tmp, abs(singvx[i - 1] - singvqx[i - 1]) / singvx[1 - 1]);
                                            }
                                            svdiff = max(svdiff, tmp);
                                            if (tmp > m * n * eps) {
                                                write(6, star), "FAILED! Something was wrong with the run.";
                                                nfail_svdiff++;
                                                for (j = 1; j <= 3; j = j + 1) {
                                                    write(6, star), j, singvx[j - 1], singvqx[j - 1];
                                                    read(5, star);
                                                }
                                            }
                                            //
                                            // ..... Rgedmdq check point
                                            if (Mlsame(wantq.elems, "Q") && Mlsame(wantr.elems, "R")) {
                                                // Check that the QR factors are computed and returned
                                                // as requested. The residual ||F-Q*R||_F / ||F||_F
                                                // is compared to M*N*EPS.
                                                copy_matrix_block(f2, ldf, 1, f, ldf, 1, m, n + 1);
                                                Rgemm("N", "N", m, n + 1, min(m, n + 1), -one, f1, ldf, y, ldy, one, f2, ldf);
                                                tmp_fqr = Rlange("F", m, n + 1, f2, ldf, work) / Rlange("F", m, n + 1, f, ldf, work);
                                                if (tmp_fqr > tol2) {
                                                    write(6, star), "FAILED! Something was wrong with the run.";
                                                    nfail_f_qr++;
                                                }
                                            }
                                            //
                                            // ..... Rgedmdq check point
                                            if (Mlsame(resids.elems, "R")) {
                                                // Compare the residuals returned by Rgedmdq with the
                                                // explicitly computed residuals using the matrix A.
                                                // Compute explicitly Y1 = A*Z
                                                Rgemm("N", "N", m, kq, m, one, a, m, z, m, zero, y1, m);
                                                // ... and then A*Z(:,i) - LAMBDA(i)*Z(:,i), using the real forms
                                                // of the invariant subspaces that correspond to complex conjugate
                                                // pairs of eigencalues. (See the description of Z in Rgedmdq)
                                                i = 1;
                                                while (i <= kq) {
                                                    if (ieigq[i - 1] == zero) {
                                                        // have a real eigenvalue with real eigenvector
                                                        Raxpy(m, -reigq[i - 1], &z[(i - 1) * ldz], 1, &y1[(i - 1) * m], 1);
                                                        // Y(1:M,i) = Y(1:M,i) - REIG(i)*Z(1:M,i)
                                                        res1[i - 1] = Rnrm2(m, &y1[(i - 1) * m], 1);
                                                        i++;
                                                    } else {
                                                        // Have a complex conjugate pair
                                                        // REIG(i) +- sqrt(-1)*IMEIG(i).
                                                        // Since all computation is done in real
                                                        // arithmetic, the formula for the residual
                                                        // is recast for real representation of the
                                                        // complex conjugate eigenpair. See the
                                                        // description of RES.
                                                        ab[0] = reigq[i - 1];
                                                        ab[(2 - 1)] = -ieigq[i - 1];
                                                        ab[(2 - 1) * ldab] = ieigq[i - 1];
                                                        ab[(2 - 1) + (2 - 1) * ldab] = reigq[i - 1];
                                                        // BLAS CALL
                                                        Rgemm("N", "N", m, 2, 2, -one, &z[(i - 1) * ldz], m, ab, 2, one, &y1[(i - 1) * m], m);
                                                        // Y(1:M,i:i+1) = Y(1:M,i:i+1) - Z(1:M,i:i+1) * AB   ! INTRINSIC
                                                        // LAPACK CALL
                                                        res1[i - 1] = Rlange("F", m, 2, &y1[(i - 1) * m], m, work);
                                                        res1[(i + 1) - 1] = res1[i - 1];
                                                        i += 2;
                                                    }
                                                }
                                                tmp = zero;
                                                for (i = 1; i <= kq; i = i + 1) {
                                                    tmp = max(tmp, abs(res[i - 1] - res1[i - 1]) * singvqx[k - 1] / (anorm * singvqx[1 - 1]));
                                                }
                                                tmp_rezq = max(tmp_rezq, tmp);
                                                if (tmp > tol2) {
                                                    nfail_rezq++;
                                                    write(6, star), "................ Rgedmdq FAILED!", "Check the code for implementation errors.";
                                                    FEM_STOP(0);
                                                }
                                                //
                                            }
                                            //
                                            for (i = 1; i <= kq; i = i + 1) {
                                                lambdaq[(i - 1)] = reigq[i - 1];
                                                lambdaq[(i - 1) + (2 - 1) * n] = ieigq[i - 1];
                                            }
                                            //
                                            // TEST_QRDMD
                                        }
                                        // ======================================================================
                                        //
                                        // LWMINOPT
                                    }
                                    // write(*,*) 'LWMINOPT loop completed'
                                    // WHTSVD LOOP
                                }
                                // write(*,*) 'WHTSVD loop completed'
                                // NRNK LOOP
                            }
                            // write(*,*) 'NRNK loop completed'
                            // SCALE LOOP
                        }
                        // write(*,*) 'SCALE loop completed'
                        // JOBF LOOP
                    }
                    // write(*,*) 'JOBREF loop completed'
                    // JOBZ LOOP
                }
                // write(*,*) 'JOBZ loop completed'
                //
                // MODE -6:6
            }
            // write(*,*) 'MODE loop completed'
            // 1 or 2 trajectories
        }
        // write(*,*) 'trajectories  loop completed'
        //
        // ............................................................
        // Generate random M-by-M matrix A. Use Rlatmr from
        // LLOOP
    }
    //
    write(6, star), ">>>>>>>>>>>>>>>>>>>>>>>>>>";
    write(6, star), " Test summary for Rgedmd :";
    write(6, star), ">>>>>>>>>>>>>>>>>>>>>>>>>>";
    write(6, star);
    if (nfail_z_xv == 0) {
        write(6, star), ">>>> Z - U*V test PASSED.";
    } else {
        write(6, star), "Z - U*V test FAILED ", nfail_z_xv, " time(s)";
        write(6, star), "Max error ||Z-U*V||_F was ", tmp_zxw;
        nfail_total += nfail_z_xv;
    }
    if (nfail_au == 0) {
        write(6, star), ">>>> A*U test PASSED. ";
    } else {
        write(6, star), "A*U test FAILED ", nfail_au, " time(s)";
        write(6, star), "Max A*U test adjusted error measure was ", tmp_au;
        write(6, star), "It should be up to O(M*N) times EPS, EPS = ", eps;
        nfail_total += nfail_au;
    }
    //
    if (nfail_rez == 0) {
        write(6, star), ">>>> Rezidual computation test PASSED.";
    } else {
        write(6, star), "Rezidual computation test FAILED ", nfail_rez, "time(s)";
        write(6, star), "Max residual computing test adjusted error measure was ", tmp_rez;
        write(6, star), "It should be up to O(M*N) times EPS, EPS = ", eps;
        nfail_total += nfail_rez;
    }
    //
    if (nfail_total == 0) {
        write(6, star), ">>>> Rgedmd :: ALL TESTS PASSED.";
    } else {
        write(6, star), nfail_total, "FAILURES!";
        write(6, star), ">>>>>>>>>>>>>> Rgedmd :: TESTS FAILED. CHECK THE IMPLEMENTATION.";
    }
    //
    if (test_qrdmd) {
        write(6, star);
        write(6, star), ">>>>>>>>>>>>>>>>>>>>>>>>>>";
        write(6, star), " Test summary for Rgedmdq :";
        write(6, star), ">>>>>>>>>>>>>>>>>>>>>>>>>>";
        write(6, star);
        //
        if (nfail_svdiff == 0) {
            write(6, star), ">>>> Rgedmd and Rgedmdq computed singularvalues test PASSED.";
        } else {
            write(6, star), "Rgedmd and Rgedmdq discrepancies inthe singular values unacceptable ", nfail_svdiff, " times. Test FAILED.";
            write(6, star),
                "The maximal discrepancy in the singular values (relative to the norm)"
                " was ",
                svdiff;
            write(6, star), "It should be up to O(M*N) times EPS, EPS = ", eps;
            nfailq_total += nfail_svdiff;
        }
        //
        if (nfail_f_qr == 0) {
            write(6, star), ">>>> F - Q*R test PASSED.";
        } else {
            write(6, star), "F - Q*R test FAILED ", nfail_f_qr, " time(s)";
            write(6, star), "The largest relative residual was ", tmp_fqr;
            write(6, star), "It should be up to O(M*N) times EPS, EPS = ", eps;
            nfailq_total += nfail_f_qr;
        }
        //
        if (nfail_rezq == 0) {
            write(6, star), ">>>> Rezidual computation test PASSED.";
        } else {
            write(6, star), "Rezidual computation test FAILED ", nfail_rezq, "time(s)";
            write(6, star), "Max residual computing test adjusted error measure was ", tmp_rezq;
            write(6, star), "It should be up to O(M*N) times EPS, EPS = ", eps;
            nfailq_total += nfail_rezq;
        }
        //
        if (nfailq_total == 0) {
            write(6, star), ">>>>>>> Rgedmdq :: ALL TESTS PASSED.";
        } else {
            write(6, star), nfailq_total, "FAILURES!";
            write(6, star), ">>>>>>> Rgedmdq :: TESTS FAILED. CHECK THE IMPLEMENTATION.";
        }
        //
    }
    //
    write(6, star);
    write(6, star), "Test completed.";
}

int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_dmd_test); }
