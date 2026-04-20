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

// Derived from LAPACK routine ZCHKDMD.
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
// Cgedmd,  for computation of the
// Dynamic Mode Decomposition (DMD)
// Cgedmdq, for computation of a
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
// ............................................................
//
// ............................................................
// ............................................................
//
void program_dmd_test(int argc, char const *argv[]) {
    common cmn(argc, argv);
    common_read read(cmn);
    common_write write(cmn);
    //
    // use iso_fortran_env, only: real64
    // integer, parameter :: WP = real64
    //
    // ............................................................
    //
    // ............................................................
    //
    // ............................................................
    // ............................................................
    //
    // .....external subroutines (BLAS and LAPACK)
    // .....external subroutines DMD package, part 1
    // subroutines under test
    // .....external functions (BLAS and LAPACK)
    //
    // ............................................................
    //
    // The test is always in pairs : ( Cgedmd and Cgedmdq )
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
    std::unique_ptr<COMPLEX[]> za_storage;
    COMPLEX *za = nullptr;
    std::unique_ptr<COMPLEX[]> zac_storage;
    COMPLEX *zac = nullptr;
    std::unique_ptr<COMPLEX[]> zf_storage;
    COMPLEX *zf = nullptr;
    std::unique_ptr<COMPLEX[]> zf0_storage;
    COMPLEX *zf0 = nullptr;
    std::unique_ptr<COMPLEX[]> zf1_storage;
    COMPLEX *zf1 = nullptr;
    std::unique_ptr<COMPLEX[]> zx_storage;
    COMPLEX *zx = nullptr;
    std::unique_ptr<COMPLEX[]> zx0_storage;
    COMPLEX *zx0 = nullptr;
    std::unique_ptr<COMPLEX[]> zy_storage;
    COMPLEX *zy = nullptr;
    std::unique_ptr<COMPLEX[]> zy0_storage;
    COMPLEX *zy0 = nullptr;
    std::unique_ptr<COMPLEX[]> zy1_storage;
    COMPLEX *zy1 = nullptr;
    std::unique_ptr<COMPLEX[]> zau_storage;
    COMPLEX *zau = nullptr;
    std::unique_ptr<COMPLEX[]> zw_storage;
    COMPLEX *zw = nullptr;
    std::unique_ptr<COMPLEX[]> zs_storage;
    COMPLEX *zs = nullptr;
    std::unique_ptr<COMPLEX[]> zz_storage;
    COMPLEX *zz = nullptr;
    std::unique_ptr<COMPLEX[]> zz1_storage;
    COMPLEX *zz1 = nullptr;
    std::unique_ptr<REAL[]> res_storage;
    REAL *res = nullptr;
    std::unique_ptr<REAL[]> res1_storage;
    REAL *res1 = nullptr;
    std::unique_ptr<REAL[]> resex_storage;
    REAL *resex = nullptr;
    std::unique_ptr<COMPLEX[]> zeigs_storage;
    COMPLEX *zeigs = nullptr;
    std::unique_ptr<REAL[]> singvx_storage;
    REAL *singvx = nullptr;
    std::unique_ptr<REAL[]> singvqx_storage;
    REAL *singvqx = nullptr;
    REAL tol = 0.0;
    REAL tol2 = 0.0;
    INTEGER k_traj = 0;
    REAL cond = 0.0;
    COMPLEX zmax = 0.0;
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
    std::unique_ptr<COMPLEX[]> zda_storage;
    COMPLEX *zda = nullptr;
    std::unique_ptr<COMPLEX[]> zdl_storage;
    COMPLEX *zdl = nullptr;
    std::unique_ptr<COMPLEX[]> zdr_storage;
    COMPLEX *zdr = nullptr;
    const REAL one = 1.0;
    INTEGER info = 0;
    INTEGER lzwork = 0;
    std::unique_ptr<COMPLEX[]> zeigsa_storage;
    COMPLEX *zeigsa = nullptr;
    std::unique_ptr<COMPLEX[]> zwork_storage;
    COMPLEX *zwork = nullptr;
    std::unique_ptr<REAL[]> work_storage;
    REAL *work = nullptr;
    COMPLEX zdum2x2[2 * 2];
    REAL tmp = 0.0;
    REAL wdummy[2];
    REAL anorm = 0.0;
    INTEGER i = 0;
    const COMPLEX zone = COMPLEX(1.0, 0.0);
    const COMPLEX zzero = COMPLEX(0.0, 0.0);
    INTEGER ijobz = 0;
    fem::str<1> jobz;
    fem::str<1> resids;
    INTEGER ijobref = 0;
    fem::str<1> jobref;
    INTEGER iscale = 0;
    fem::str<1> scale;
    INTEGER inrnk = 0;
    INTEGER nrnk = 0;
    INTEGER nrnksp = 0;
    INTEGER iwhtsvd = 0;
    INTEGER whtsvd = 0;
    INTEGER whtsvdsp = 0;
    INTEGER lwminopt = 0;
    INTEGER k = 0;
    COMPLEX zdummy[22];
    INTEGER idummy[4];
    INTEGER lwork = 0;
    INTEGER liwork = 0;
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
        iseed[1 - 1] = 4;
        iseed[2 - 1] = 3;
        iseed[3 - 1] = 2;
        iseed[4 - 1] = 1;
        //
        lda = m;
        ldf = m;
        ldx = m;
        ldy = m;
        ldw = n;
        ldz = m;
        ldau = m;
        lds = n;
        //
        tmp_zxw = zero;
        tmp_au = zero;
        tmp_rez = zero;
        tmp_rezq = zero;
        svdiff = zero;
        tmp_ex = zero;
        //
        za_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, lda * m));
        za = za_storage.get();
        zac_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, lda * m));
        zac = zac_storage.get();
        zf_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, ldf * (n + 1)));
        zf = zf_storage.get();
        zf0_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, ldf * (n + 1)));
        zf0 = zf0_storage.get();
        zf1_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, ldf * (n + 1)));
        zf1 = zf1_storage.get();
        zx_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, ldx * n));
        zx = zx_storage.get();
        zx0_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, ldx * n));
        zx0 = zx0_storage.get();
        zy_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, ldy * (n + 1)));
        zy = zy_storage.get();
        zy0_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, ldy * (n + 1)));
        zy0 = zy0_storage.get();
        zy1_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, ldy * (n + 1)));
        zy1 = zy1_storage.get();
        zau_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, ldau * n));
        zau = zau_storage.get();
        zw_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, ldw * n));
        zw = zw_storage.get();
        zs_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, lds * n));
        zs = zs_storage.get();
        zz_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, ldz * n));
        zz = zz_storage.get();
        zz1_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, ldz * n));
        zz1 = zz1_storage.get();
        res_storage = std::make_unique<REAL[]>(max((INTEGER)1, n));
        res = res_storage.get();
        res1_storage = std::make_unique<REAL[]>(max((INTEGER)1, n));
        res1 = res1_storage.get();
        resex_storage = std::make_unique<REAL[]>(max((INTEGER)1, n));
        resex = resex_storage.get();
        zeigs_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, n));
        zeigs = zeigs_storage.get();
        singvx_storage = std::make_unique<REAL[]>(max((INTEGER)1, n));
        singvx = singvx_storage.get();
        singvqx_storage = std::make_unique<REAL[]>(max((INTEGER)1, n));
        singvqx = singvqx_storage.get();
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
            cond = 10000.0;
            zmax = COMPLEX(10.0, 10.0);
            rsign = "F";
            grade = "N";
            model = 6;
            condl = 10.0;
            moder = 6;
            condr = 10.0;
            pivtng = "N";
            //
            // Loop over all parameter MODE values for Clatmr (+1,..,+6)
            for (mode = 1; mode <= 6; mode = mode + 1) {
                //
                iwork_storage = std::make_unique<INTEGER[]>(max((INTEGER)1, (2 * m)));
                iwork = iwork_storage.get();
                zda_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, m));
                zda = zda_storage.get();
                zdl_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, m));
                zdl = zdl_storage.get();
                zdr_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, m));
                zdr = zdr_storage.get();
                //
                Clatmr(m, m, "N", iseed, "N", zda, mode, cond, zmax, rsign, grade, zdl, model, condl, zdr, moder, condr, pivtng, iwork, m, m, zero, -one, "N", za, lda, &iwork[(m + 1) - 1], info);
                //
                lzwork = max((INTEGER)1, 2 * m);
                zeigsa_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, m));
                zeigsa = zeigsa_storage.get();
                zwork_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, lzwork));
                zwork = zwork_storage.get();
                work_storage = std::make_unique<REAL[]>(max((INTEGER)1, (2 * m)));
                work = work_storage.get();
                copy_matrix_block(zac, lda, 1, za, lda, 1, m, m);
                // LAPACK CALL
                Cgeev("N", "N", m, zac, lda, zeigsa, zdum2x2, 2, zdum2x2, 2, zwork, lzwork, work, info);
                //
                // The spectral radius of ZA
                tmp = abs(zeigsa[iCamax(m, zeigsa, 1) - 1]);
                // Scale the matrix ZA to have unit spectral radius.
                Clascl("G", 0, 0, tmp, one, m, m, za, lda, info);
                Clascl("G", 0, 0, tmp, one, m, 1, zeigsa, m, info);
                anorm = Clange("F", m, m, za, lda, wdummy);
                //
                if (k_traj == 2) {
                    // generate data as two trajectories
                    // with two inital conditions
                    Clarnv(2, iseed, m, &zf[0]);
                    for (i = 1; i <= n / 2; i = i + 1) {
                        Cgemv("N", m, m, zone, za, lda, &zf[(i - 1) * ldf], 1, zzero, &zf[((i + 1) - 1) * ldf], 1);
                    }
                    copy_matrix_block(zx0, ldx, 1, zf, ldf, 1, m, n / 2);
                    copy_matrix_block(zy0, ldy, 1, zf, ldf, 2, m, n / 2);
                    //
                    Clarnv(2, iseed, m, &zf[0]);
                    for (i = 1; i <= n - n / 2; i = i + 1) {
                        Cgemv("N", m, m, zone, za, lda, &zf[(i - 1) * ldf], 1, zzero, &zf[((i + 1) - 1) * ldf], 1);
                    }
                    copy_matrix_block(zx0, ldx, n / 2 + 1, zf, ldf, 1, m, n - n / 2);
                    copy_matrix_block(zy0, ldy, n / 2 + 1, zf, ldf, 2, m, n - n / 2);
                } else {
                    Clarnv(2, iseed, m, &zf[0]);
                    for (i = 1; i <= n; i = i + 1) {
                        Cgemv("N", m, m, zone, za, m, &zf[(i - 1) * ldf], 1, zzero, &zf[((i + 1) - 1) * ldf], 1);
                    }
                    copy_matrix_block(zf0, ldf, 1, zf, ldf, 1, m, n + 1);
                    copy_matrix_block(zx0, ldx, 1, zf0, ldf, 1, m, n);
                    copy_matrix_block(zy0, ldy, 1, zf0, ldf, 2, m, n);
                }
                //
                // ........................................................................
                //
                for (ijobz = 1; ijobz <= 4; ijobz = ijobz + 1) {
                    //
                    if ((ijobz == 1)) {
                        jobz = "V";
                        resids = "R";
                    } else if ((ijobz == 2)) {
                        jobz = "V";
                        resids = "N";
                    } else if ((ijobz == 3)) {
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
                            jobref = "R";
                        } else if ((ijobref == 2)) {
                            jobref = "E";
                        } else if ((ijobref == 3)) {
                            jobref = "N";
                        }
                        //
                        for (iscale = 1; iscale <= 4; iscale = iscale + 1) {
                            //
                            if ((iscale == 1)) {
                                scale = "S";
                            } else if ((iscale == 2)) {
                                scale = "C";
                            } else if ((iscale == 3)) {
                                scale = "Y";
                            } else if ((iscale == 4)) {
                                scale = "N";
                            }
                            //
                            for (inrnk = -1; inrnk >= -2; inrnk = inrnk - 1) {
                                nrnk = inrnk;
                                nrnksp = inrnk;
                                //
                                for (iwhtsvd = 1; iwhtsvd <= 3; iwhtsvd = iwhtsvd + 1) {
                                    // Check all four options to compute the POD basis
                                    // via the SVD.
                                    whtsvd = iwhtsvd;
                                    whtsvdsp = iwhtsvd;
                                    //
                                    for (lwminopt = 1; lwminopt <= 2; lwminopt = lwminopt + 1) {
                                        // Workspace query for the minimal (1) and for the optimal
                                        // (2) workspace lengths determined by workspace query.
                                        //
                                        // Cgedmd is always tested and its results are also used for
                                        // comparisons with Cgedmdq.
                                        //
                                        copy_matrix_block(zx, ldx, 1, zx0, ldx, 1, m, n);
                                        copy_matrix_block(zy, ldy, 1, zy0, ldy, 1, m, n);
                                        //
                                        Cgedmd(scale.elems, jobz.elems, resids.elems, jobref.elems, whtsvd, m, n, zx, ldx, zy, ldy, nrnk, tol, k, zeigs, zz, ldz, res, zau, ldau, zw, ldw, zs, lds, zdummy, -1, wdummy, -1, idummy, -1, info);
                                        if ((info == 2) || (info == 3) || (info < 0)) {
                                            write(6, star), "Call to Cgedmd workspace query failed.Check the calli"
                                                            "ng sequence and the code.";
                                            write(6, star), "The error code is ", info;
                                            write(6, star), "The input parameters were ", scale, jobz, resids, jobref, whtsvd, m, n, ldx, ldy, nrnk, tol, ldz, ldau, ldw, lds;
                                            FEM_STOP(0);
                                        }
                                        //
                                        lzwork = castINTEGER(zdummy[lwminopt - 1].real());
                                        lwork = castINTEGER(wdummy[1 - 1]);
                                        liwork = idummy[1 - 1];
                                        //
                                        zwork_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, lzwork));
                                        zwork = zwork_storage.get();
                                        work_storage = std::make_unique<REAL[]>(max((INTEGER)1, lwork));
                                        work = work_storage.get();
                                        iwork_storage = std::make_unique<INTEGER[]>(max((INTEGER)1, liwork));
                                        iwork = iwork_storage.get();
                                        //
                                        Cgedmd(scale.elems, jobz.elems, resids.elems, jobref.elems, whtsvd, m, n, zx, ldx, zy, ldy, nrnk, tol, k, zeigs, zz, ldz, res, zau, ldau, zw, ldw, zs, lds, zwork, lzwork, work, lwork, iwork, liwork, info);
                                        //
                                        if (info != 0) {
                                            write(6, star), "Call to Cgedmd failed.Check the calling sequence and "
                                                            "the code.";
                                            write(6, star), "The error code is ", info;
                                            write(6, star), "The input parameters were ", scale, jobz, resids, jobref, whtsvd, m, n, ldx, ldy, nrnk, tol;
                                            FEM_STOP(0);
                                        }
                                        //
                                        for (INTEGER i_ = 1; i_ <= n; i_++) {
                                            singvx[i_ - 1] = work[i_ - 1];
                                        }
                                        //
                                        // ...... Cgedmd check point
                                        if (Mlsame(jobz.elems, "V")) {
                                            // Check that Z = X*W, on return from Cgedmd
                                            // This checks that the returned eigenvectors in Z are
                                            // the product of the SVD'POD basis returned in X
                                            // and the eigenvectors of the rayleigh quotient
                                            // returned in W
                                            Cgemm("N", "N", m, k, k, zone, zx, ldx, zw, ldw, zzero, zz1, ldz);
                                            tmp = zero;
                                            for (i = 1; i <= k; i = i + 1) {
                                                Caxpy(m, -zone, &zz[(i - 1) * ldz], 1, &zz1[(i - 1) * ldz], 1);
                                                tmp = max(tmp, RCnrm2(m, &zz1[(i - 1) * ldz], 1));
                                            }
                                            tmp_zxw = max(tmp_zxw, tmp);
                                            if (tmp_zxw <= 10 * m * eps) {
                                                // WRITE(*,*) ' :) .... OK .........Cgedmd PASSED.'
                                            } else {
                                                nfail_z_xv++;
                                                write(6, star), ":( .................Cgedmd FAILED!", "Check the code for implementation errors.";
                                                write(6, star), "The input parameters were ", scale, jobz, resids, jobref, whtsvd, m, n, ldx, ldy, nrnk, tol;
                                            }
                                        }
                                        //
                                        // ...... Cgedmd check point
                                        if (Mlsame(jobref.elems, "R")) {
                                            // The matrix A*U is returned for computing refined Ritz vectors.
                                            // Check that A*U is computed correctly using the formula
                                            // A*U = Y * V * inv(SIGMA). This depends on the
                                            // accuracy in the computed singular values and vectors of X.
                                            // See the paper for an error analysis.
                                            // Note that the left singular vectors of the input matrix X
                                            // are returned in the array X.
                                            Cgemm("N", "N", m, k, m, zone, za, lda, zx, ldx, zzero, zz1, ldz);
                                            tmp = zero;
                                            for (i = 1; i <= k; i = i + 1) {
                                                Caxpy(m, -zone, &zau[(i - 1) * ldau], 1, &zz1[(i - 1) * ldz], 1);
                                                tmp = max(tmp, RCnrm2(m, &zz1[(i - 1) * ldz], 1) * singvx[k - 1] / (anorm * singvx[1 - 1]));
                                            }
                                            tmp_au = max(tmp_au, tmp);
                                            if (tmp <= tol2) {
                                                // WRITE(*,*) ':) .... OK .........Cgedmd PASSED.'
                                            } else {
                                                nfail_au++;
                                                write(6, star), ":( .................Cgedmd FAILED!", "Check the code for implementation errors.";
                                                write(6, star), "The input parameters were ", scale, jobz, resids, jobref, whtsvd, m, n, ldx, ldy, nrnk, tol;
                                            }
                                        } else if (Mlsame(jobref.elems, "E")) {
                                            // The unscaled vectors of the Exact DMD are computed.
                                            // This option is included for the sake of completeness,
                                            // for users who prefer the Exact DMD vectors. The
                                            // returned vectors are in the real form, in the same way
                                            // as the Ritz vectors. Here we just save the vectors
                                            // and test them separately using a Matlab script.
                                            //
                                            Cgemm("N", "N", m, k, m, zone, za, lda, zau, ldau, zzero, zy1, ldy);
                                            //
                                            for (i = 1; i <= k; i = i + 1) {
                                                // have a real eigenvalue with real eigenvector
                                                Caxpy(m, -zeigs[i - 1], &zau[(i - 1) * ldau], 1, &zy1[(i - 1) * ldy], 1);
                                                resex[i - 1] = RCnrm2(m, &zy1[(i - 1) * ldy], 1) / RCnrm2(m, &zau[(i - 1) * ldau], 1);
                                            }
                                        }
                                        // ...... Cgedmd check point
                                        //
                                        if (Mlsame(resids.elems, "R")) {
                                            // Compare the residuals returned by Cgedmd with the
                                            // explicitly computed residuals using the matrix A.
                                            // Compute explicitly Y1 = A*Z
                                            Cgemm("N", "N", m, k, m, zone, za, lda, zz, ldz, zzero, zy1, ldy);
                                            // ... and then A*Z(:,i) - LAMBDA(i)*Z(:,i), using the real forms
                                            // of the invariant subspaces that correspond to complex conjugate
                                            // pairs of eigencalues. (See the description of Z in Cgedmd,)
                                            //
                                            for (i = 1; i <= k; i = i + 1) {
                                                // have a real eigenvalue with real eigenvector
                                                Caxpy(m, -zeigs[i - 1], &zz[(i - 1) * ldz], 1, &zy1[(i - 1) * ldy], 1);
                                                res1[i - 1] = RCnrm2(m, &zy1[(i - 1) * ldy], 1);
                                            }
                                            tmp = zero;
                                            for (i = 1; i <= k; i = i + 1) {
                                                tmp = max(tmp, abs(res[i - 1] - res1[i - 1]) * singvx[k - 1] / (anorm * singvx[1 - 1]));
                                            }
                                            tmp_rez = max(tmp_rez, tmp);
                                            if (tmp <= tol2) {
                                                // WRITE(*,*) ':) .... OK ..........Cgedmd PASSED.'
                                            } else {
                                                nfail_rez++;
                                                write(6, star), ":( ..................Cgedmd FAILED!", "Check the code for implementation errors.";
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
                                        if (test_qrdmd && (k_traj == 1)) {
                                            //
                                            copy_matrix_block(zf, ldf, 1, zf0, ldf, 1, m, n + 1);
                                            //
                                            Cgedmdq(scale.elems, jobz.elems, resids.elems, wantq.elems, wantr.elems, jobref.elems, whtsvd, m, n + 1, zf, ldf, zx, ldx, zy, ldy, nrnk, tol, k, zeigs, zz, ldz, res, zau, ldau, zw, ldw, zs, lds, zdummy, -1, wdummy, -1, idummy, -1, info);
                                            //
                                            lzwork = castINTEGER(zdummy[lwminopt - 1].real());
                                            zwork_storage = std::make_unique<COMPLEX[]>(max((INTEGER)1, lzwork));
                                            zwork = zwork_storage.get();
                                            liwork = idummy[1 - 1];
                                            iwork_storage = std::make_unique<INTEGER[]>(max((INTEGER)1, liwork));
                                            iwork = iwork_storage.get();
                                            lwork = castINTEGER(wdummy[1 - 1]);
                                            work_storage = std::make_unique<REAL[]>(max((INTEGER)1, lwork));
                                            work = work_storage.get();
                                            //
                                            Cgedmdq(scale.elems, jobz.elems, resids.elems, wantq.elems, wantr.elems, jobref.elems, whtsvd, m, n + 1, zf, ldf, zx, ldx, zy, ldy, nrnk, tol, kq, zeigs, zz, ldz, res, zau, ldau, zw, ldw, zs, lds, zwork, lzwork, work, lwork, iwork, liwork, info);
                                            //
                                            if (info != 0) {
                                                write(6, star), "Call to Cgedmdq failed.Check the calling sequence a"
                                                                "nd the code.";
                                                write(6, star), "The error code is ", info;
                                                write(6, star), "The input parameters were ", scale, jobz, resids, wantq, wantr, whtsvd, m, n, ldx, ldy, nrnk, tol;
                                                FEM_STOP(0);
                                            }
                                            for (INTEGER i_ = 1; i_ <= n; i_++) {
                                                singvqx[i_ - 1] = work[i_ - 1];
                                            }
                                            //
                                            // ..... Cgedmdq check point
                                            //
                                            if (1 == 0) {
                                                // Comparison of Cgedmd and Cgedmdq singular values disabled
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
                                                    //
                                                }
                                            }
                                            //
                                            // ..... Cgedmdq check point
                                            if (Mlsame(wantq.elems, "Q") && Mlsame(wantr.elems, "R")) {
                                                // Check that the QR factors are computed and returned
                                                // as requested. The residual ||F-Q*R||_F / ||F||_F
                                                // is compared to M*N*EPS.
                                                copy_matrix_block(zf1, ldf, 1, zf0, ldf, 1, m, n + 1);
                                                Cgemm("N", "N", m, n + 1, min(m, n + 1), -zone, zf, ldf, zy, ldy, zone, zf1, ldf);
                                                tmp_fqr = Clange("F", m, n + 1, zf1, ldf, work) / Clange("F", m, n + 1, zf0, ldf, work);
                                                if (tmp_fqr > tol2) {
                                                    write(6, star), "FAILED! Something was wrong with the run.";
                                                    nfail_f_qr++;
                                                } else {
                                                    // WRITE(*,*) '........ PASSED.'
                                                }
                                            }
                                            //
                                            // ..... Cgedmdq check point
                                            if (Mlsame(resids.elems, "R")) {
                                                // Compare the residuals returned by Cgedmdq with the
                                                // explicitly computed residuals using the matrix A.
                                                // Compute explicitly Y1 = A*Z
                                                Cgemm("N", "N", m, kq, m, zone, za, lda, zz, ldz, zzero, zy1, ldy);
                                                // ... and then A*Z(:,i) - LAMBDA(i)*Z(:,i), using the real forms
                                                // of the invariant subspaces that correspond to complex conjugate
                                                // pairs of eigencalues. (See the description of Z in Cgedmdq)
                                                //
                                                for (i = 1; i <= kq; i = i + 1) {
                                                    // have a real eigenvalue with real eigenvector
                                                    Caxpy(m, -zeigs[i - 1], &zz[(i - 1) * ldz], 1, &zy1[(i - 1) * ldy], 1);
                                                    // Y(1:M,i) = Y(1:M,i) - REIG(i)*Z(1:M,i)
                                                    res1[i - 1] = RCnrm2(m, &zy1[(i - 1) * ldy], 1);
                                                }
                                                tmp = zero;
                                                for (i = 1; i <= kq; i = i + 1) {
                                                    tmp = max(tmp, abs(res[i - 1] - res1[i - 1]) * singvqx[kq - 1] / (anorm * singvqx[1 - 1]));
                                                }
                                                tmp_rezq = max(tmp_rezq, tmp);
                                                if (tmp <= tol2) {
                                                    // WRITE(*,*) '.... OK ........ Cgedmdq PASSED.'
                                                } else {
                                                    nfail_rezq++;
                                                    write(6, star), "................ Cgedmdq FAILED!", "Check the code for implementation errors.";
                                                    FEM_STOP(0);
                                                }
                                                //
                                            }
                                            //
                                            // Cgedmdq
                                        }
                                        //
                                        // .......................................................................................................
                                        //
                                        // LWMINOPT
                                    }
                                    // write(*,*) 'LWMINOPT loop completed'
                                    // iWHTSVD
                                }
                                // write(*,*) 'WHTSVD loop completed'
                                // iNRNK  -2:-1
                            }
                            // write(*,*) 'NRNK loop completed'
                            // iSCALE  1:4
                        }
                        // write(*,*) 'SCALE loop completed'
                    }
                    // write(*,*) 'JOBREF loop completed'
                    // iJOBZ
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
        // LLOOP
    }
    //
    write(6, star), ">>>>>>>>>>>>>>>>>>>>>>>>>>";
    write(6, star), " Test summary for Cgedmd :";
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
        write(6, star), ">>>> Cgedmd :: ALL TESTS PASSED.";
    } else {
        write(6, star), nfail_total, "FAILURES!";
        write(6, star), ">>>>>>>>>>>>>> Cgedmd :: TESTS FAILED. CHECK THE IMPLEMENTATION.";
    }
    //
    if (test_qrdmd) {
        write(6, star);
        write(6, star), ">>>>>>>>>>>>>>>>>>>>>>>>>>";
        write(6, star), " Test summary for Cgedmdq :";
        write(6, star), ">>>>>>>>>>>>>>>>>>>>>>>>>>";
        write(6, star);
        //
        if (nfail_svdiff == 0) {
            write(6, star), ">>>> Cgedmd and Cgedmdq computed singularvalues test PASSED.";
        } else {
            write(6, star), "Cgedmd and Cgedmdq discrepancies inthe singular values unacceptable ", nfail_svdiff, " times. Test FAILED.";
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
            write(6, star), ">>>>>>> Cgedmdq :: ALL TESTS PASSED.";
        } else {
            write(6, star), nfailq_total, "FAILURES!";
            write(6, star), ">>>>>>> Cgedmdq :: TESTS FAILED. CHECK THE IMPLEMENTATION.";
        }
        //
    }
    //
    write(6, star);
    write(6, star), "Test completed.";
}

int main(int argc, char const *argv[]) { return fem::main_with_catch(argc, argv, program_dmd_test); }
