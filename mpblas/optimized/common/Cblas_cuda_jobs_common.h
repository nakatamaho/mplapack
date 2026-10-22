/*
 * Copyright (c) 2026
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

// The jobs of Cblas_cuda_common.cuh: one per routine, with the arguments of
// the call and the device addresses of its operands (A, B, C in this order);
// job(u) runs the reference routine on part u of the output.  Shared by the
// device code and its host side (blas_cuda/Cblas_bridge_cuda.cpp).

#ifndef _MPLAPACK_CBLAS_CUDA_JOBS_COMMON_H_
#define _MPLAPACK_CBLAS_CUDA_JOBS_COMMON_H_

#include "blas_dev/blas_dev.h"
#include "Rblas_gpu_operand_common.h"

namespace mplapack_blas_dev {

template <class Ops> struct job_Cgemm {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    char ta[2], tb[2]; long m, n, k, lda, ldb, ldc; DC alpha, beta;
    E *A, *B, *C;
    MPLAPACK_DEV_HD long units() const { return m * n; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; B = d[1]; C = d[2]; }
    MPLAPACK_DEV_HD void operator()(long u) const { const long i = u % m, j = u / m; Cgemm_dev<DR, DC>(ta, tb, 1, 1, k, alpha, A + (ta[0] == 'N' || ta[0] == 'n' ? i : i * lda), lda, B + (tb[0] == 'N' || tb[0] == 'n' ? j * ldb : j), ldb, beta, C + i + j * ldc, ldc); }
};

template <class Ops> struct job_Csymm {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    char side[2], uplo[2]; long m, n, lda, ldb, ldc; DC alpha, beta;
    E *A, *B, *C;
    MPLAPACK_DEV_HD long units() const { return lsame(side, "L") ? n : m; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; B = d[1]; C = d[2]; }
    MPLAPACK_DEV_HD void operator()(long u) const { if (lsame(side, "L")) Csymm_dev<DR, DC>(side, uplo, m, 1, alpha, A, lda, B + u * ldb, ldb, beta, C + u * ldc, ldc); else Csymm_dev<DR, DC>(side, uplo, 1, n, alpha, A, lda, B + u, ldb, beta, C + u, ldc); }
};

template <class Ops> struct job_Chemm {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    char side[2], uplo[2]; long m, n, lda, ldb, ldc; DC alpha, beta;
    E *A, *B, *C;
    MPLAPACK_DEV_HD long units() const { return lsame(side, "L") ? n : m; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; B = d[1]; C = d[2]; }
    MPLAPACK_DEV_HD void operator()(long u) const { if (lsame(side, "L")) Chemm_dev<DR, DC>(side, uplo, m, 1, alpha, A, lda, B + u * ldb, ldb, beta, C + u * ldc, ldc); else Chemm_dev<DR, DC>(side, uplo, 1, n, alpha, A, lda, B + u, ldb, beta, C + u, ldc); }
};

template <class Ops> struct job_Ctrmm {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    char side[2], uplo[2], ta[2], diag[2]; long m, n, lda, ldb; DC alpha;
    E *A, *B;
    MPLAPACK_DEV_HD long units() const { return lsame(side, "L") ? n : m; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; B = d[1]; }
    MPLAPACK_DEV_HD void operator()(long u) const { if (lsame(side, "L")) Ctrmm_dev<DR, DC>(side, uplo, ta, diag, m, 1, alpha, A, lda, B + u * ldb, ldb); else Ctrmm_dev<DR, DC>(side, uplo, ta, diag, 1, n, alpha, A, lda, B + u, ldb); }
};

template <class Ops> struct job_Ctrsm {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    char side[2], uplo[2], ta[2], diag[2]; long m, n, lda, ldb; DC alpha;
    E *A, *B;
    MPLAPACK_DEV_HD long units() const { return lsame(side, "L") ? n : m; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; B = d[1]; }
    MPLAPACK_DEV_HD void operator()(long u) const { if (lsame(side, "L")) Ctrsm_dev<DR, DC>(side, uplo, ta, diag, m, 1, alpha, A, lda, B + u * ldb, ldb); else Ctrsm_dev<DR, DC>(side, uplo, ta, diag, 1, n, alpha, A, lda, B + u, ldb); }
};

template <class Ops> struct job_Cgemv {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    char tr[2]; long m, n, lda, incx, incy; DC alpha, beta;
    E *A, *B, *C;
    MPLAPACK_DEV_HD long units() const { return lsame(tr, "N") ? m : n; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; B = d[1]; C = d[2]; }
    MPLAPACK_DEV_HD void operator()(long u) const { if (lsame(tr, "N")) Cgemv_dev<DR, DC>(tr, 1, n, alpha, A + u, lda, B, incx, beta, C + u * incy, incy); else Cgemv_dev<DR, DC>(tr, m, 1, alpha, A + u * lda, lda, B, incx, beta, C + u * incy, incy); }
};

template <class Ops> struct job_Cgeru {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    long m, n, incx, incy, lda; DC alpha;
    E *A, *B, *C;
    MPLAPACK_DEV_HD long units() const { return m * n; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; B = d[1]; C = d[2]; }
    MPLAPACK_DEV_HD void operator()(long u) const { const long i = u % m, j = u / m; Cgeru_dev<DR, DC>(1, 1, alpha, A + i * incx, incx, B + j * incy, incy, C + i + j * lda, lda); }
};

template <class Ops> struct job_Cgerc {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    long m, n, incx, incy, lda; DC alpha;
    E *A, *B, *C;
    MPLAPACK_DEV_HD long units() const { return m * n; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; B = d[1]; C = d[2]; }
    MPLAPACK_DEV_HD void operator()(long u) const { const long i = u % m, j = u / m; Cgerc_dev<DR, DC>(1, 1, alpha, A + i * incx, incx, B + j * incy, incy, C + i + j * lda, lda); }
};

template <class Ops> struct job_Caxpy {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    long n, incx, incy; DC alpha;
    E *A, *C;
    MPLAPACK_DEV_HD long units() const { return n; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; C = d[1]; }
    MPLAPACK_DEV_HD void operator()(long u) const { Caxpy_dev<DR, DC>(1, alpha, A + u * incx, incx, C + u * incy, incy); }
};

template <class Ops> struct job_Ccopy {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    long n, incx, incy;
    E *A, *C;
    MPLAPACK_DEV_HD long units() const { return n; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; C = d[1]; }
    MPLAPACK_DEV_HD void operator()(long u) const { Ccopy_dev<DR, DC>(1, A + u * incx, incx, C + u * incy, incy); }
};

template <class Ops> struct job_Cswap {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    long n, incx, incy;
    E *A, *C;
    MPLAPACK_DEV_HD long units() const { return n; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; C = d[1]; }
    MPLAPACK_DEV_HD void operator()(long u) const { Cswap_dev<DR, DC>(1, A + u * incx, incx, C + u * incy, incy); }
};

template <class Ops> struct job_CRrot {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    long n, incx, incy; DR c, s;
    E *A, *C;
    MPLAPACK_DEV_HD long units() const { return n; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; C = d[1]; }
    MPLAPACK_DEV_HD void operator()(long u) const { CRrot_dev<DR, DC>(1, A + u * incx, incx, C + u * incy, incy, c, s); }
};

template <class Ops> struct job_Cscal {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    long n, incx; DC alpha;
    E *C;
    MPLAPACK_DEV_HD long units() const { return n; }
    MPLAPACK_DEV_HD void bind(E **d) { C = d[0]; }
    MPLAPACK_DEV_HD void operator()(long u) const { Cscal_dev<DR, DC>(1, alpha, C + u * incx, incx); }
};

template <class Ops> struct job_CRscal {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    long n, incx; DR da;
    E *C;
    MPLAPACK_DEV_HD long units() const { return n; }
    MPLAPACK_DEV_HD void bind(E **d) { C = d[0]; }
    MPLAPACK_DEV_HD void operator()(long u) const { CRscal_dev<DR, DC>(1, da, C + u * incx, incx); }
};

template <class Ops> struct job_Csyrk {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    char uplo[2], tr[2]; long n, k, lda, ldc; DC alpha, beta;
    E *A, *C;
    MPLAPACK_DEV_HD long units() const { return n; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; C = d[1]; }
    MPLAPACK_DEV_HD void operator()(long u) const { Csyrk_columns_dev<DR, DC>(uplo, tr, n, k, alpha, A, lda, beta, C, ldc, u, u + 1); }
};

template <class Ops> struct job_Cherk {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    char uplo[2], tr[2]; long n, k, lda, ldc; DR alpha, beta;
    E *A, *C;
    MPLAPACK_DEV_HD long units() const { return n; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; C = d[1]; }
    MPLAPACK_DEV_HD void operator()(long u) const { Cherk_columns_dev<DR, DC>(uplo, tr, n, k, alpha, A, lda, beta, C, ldc, u, u + 1); }
};

template <class Ops> struct job_Csyr2k {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    char uplo[2], tr[2]; long n, k, lda, ldb, ldc; DC alpha, beta;
    E *A, *B, *C;
    MPLAPACK_DEV_HD long units() const { return n; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; B = d[1]; C = d[2]; }
    MPLAPACK_DEV_HD void operator()(long u) const { Csyr2k_columns_dev<DR, DC>(uplo, tr, n, k, alpha, A, lda, B, ldb, beta, C, ldc, u, u + 1); }
};

template <class Ops> struct job_Cher2k {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    char uplo[2], tr[2]; long n, k, lda, ldb, ldc; DC alpha; DR beta;
    E *A, *B, *C;
    MPLAPACK_DEV_HD long units() const { return n; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; B = d[1]; C = d[2]; }
    MPLAPACK_DEV_HD void operator()(long u) const { Cher2k_columns_dev<DR, DC>(uplo, tr, n, k, alpha, A, lda, B, ldb, beta, C, ldc, u, u + 1); }
};

template <class Ops> struct job_Cgemmtr {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    char uplo[2], ta[2], tb[2]; long n, k, lda, ldb, ldc; DC alpha, beta;
    E *A, *B, *C;
    MPLAPACK_DEV_HD long units() const { return n; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; B = d[1]; C = d[2]; }
    MPLAPACK_DEV_HD void operator()(long u) const { Cgemmtr_columns_dev<DR, DC>(uplo, ta, tb, n, k, alpha, A, lda, B, ldb, beta, C, ldc, u, u + 1); }
};

template <class Ops> struct job_Rgemmtr {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DR E;
    char uplo[2], ta[2], tb[2]; long n, k, lda, ldb, ldc; DR alpha, beta;
    E *A, *B, *C;
    MPLAPACK_DEV_HD long units() const { return n; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; B = d[1]; C = d[2]; }
    MPLAPACK_DEV_HD void operator()(long u) const { Rgemmtr_columns_dev<DR, DC>(uplo, ta, tb, n, k, alpha, A, lda, B, ldb, beta, C, ldc, u, u + 1); }
};

template <class Ops> struct job_Cher {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    char uplo[2]; long n, incx, lda; DR alpha;
    E *A, *C;
    MPLAPACK_DEV_HD long units() const { return n; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; C = d[1]; }
    MPLAPACK_DEV_HD void operator()(long u) const { Cher_columns_dev<DR, DC>(uplo, n, alpha, A, incx, C, lda, u, u + 1); }
};

template <class Ops> struct job_Cher2 {
    typedef typename cgpu_types<Ops>::DR DR;
    typedef typename cgpu_types<Ops>::DC DC;
    typedef DC E;
    char uplo[2]; long n, incx, incy, lda; DC alpha;
    E *A, *B, *C;
    MPLAPACK_DEV_HD long units() const { return n; }
    MPLAPACK_DEV_HD void bind(E **d) { A = d[0]; B = d[1]; C = d[2]; }
    MPLAPACK_DEV_HD void operator()(long u) const { Cher2_columns_dev<DR, DC>(uplo, n, alpha, A, incx, B, incy, C, lda, u, u + 1); }
};

// Runs the job on the GPU: 0, or a CUDA error (the call is then left to the CPU).
template <class Job> int cgpu_run(Job job, const mplapack_blas::gpu_operand<typename Job::E> *ops, int nops);

} // namespace mplapack_blas_dev

#endif
