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

/*
 * The top level of the Winograd recursion (Rgemm_winograd_common.h,
 * ../mpfr/cuda/Rgemm_winograd_cuda.h) on several devices: its seven products
 * M1..M7 are dealt out cyclically (Mt on device (t - 1) % G), each device
 * forming the S terms from its own copy of the operands and computing its
 * products with the one-device recursion; then the products are copied to
 * device 0, which combines them and adds the peeled row, column and step of k.
 * Every element goes through the same operations in the same order as in
 * winograd() on one device, so the result is the same.
 *
 * be[d]: a backend (Rgemm_winograd_common.h interface) on device d; A[d] and
 * B[d]: op(A) (m x k, ld m) and op(B) (k x n, ld k) on device d; C (m x n,
 * ld m) on device 0.  rec(backend, m, n, k, A, lda, B, ldb, C, ldc) runs the
 * one-device recursion with the caller's cutoff.  Only for m, n, k whose
 * even parts exceed the cutoff (winograd_multi_applies).  Returns 0 or a
 * CUDA error code; the backends' failed() is checked by the caller as well.
 */

#ifndef _MPLAPACK_RGEMM_WINOGRAD_MULTI_GPU_CUH_
#define _MPLAPACK_RGEMM_WINOGRAD_MULTI_GPU_CUH_

#include <cuda_runtime.h>
#include <thread>
#include <vector>

namespace mplapack_gemm {

inline bool winograd_multi_applies(long m, long n, long k, long cutoff) {
    return cutoff > 0 && (m & ~1L) > cutoff && (n & ~1L) > cutoff && (k & ~1L) > cutoff;
}

template <class F, class Backend, class Rec> int winograd_multi(int G, Backend *const *be, long m, long n, long k, const F *const *A, const F *const *B, F *C, Rec rec) {
    const long me = m & ~1L, ne = n & ~1L, ke = k & ~1L;
    const long mh = me / 2, nh = ne / 2, kh = ke / 2, lda = m, ldb = k, ldc = m;
    Backend &b0 = *be[0];
    if (cudaSetDevice(0) != cudaSuccess)
        return 1;
    // The products cross devices: plain device memory (the backends may
    // allocate from a pool that other devices cannot reach).
    const size_t bytesM = sizeof(F) * (size_t)mh * nh;
    F *M[7] = {}, *T[2];
    int err = 0;
    for (int t = 0; t < 7 && !err; t++)
        err = (int)cudaMalloc((void **)&M[t], bytesM);
    for (int t = 0; t < 7 && err; t++) {
        cudaFree(M[t]);
        M[t] = NULL;
    }
    if (err)
        return err;
    for (int t = 0; t < 2; t++)
        T[t] = b0.alloc((size_t)mh * nh);
    std::vector<int> rc((size_t)G, 0);
    if (!b0.failed()) {
        std::vector<std::thread> workers;
        for (int d = 0; d < G; d++)
            workers.emplace_back([&, d]() {
                cudaError_t e = cudaSetDevice(d);
                if (e != cudaSuccess) {
                    rc[(size_t)d] = (int)e;
                    return;
                }
                Backend &b = *be[d];
                const F *a = A[d], *bb = B[d];
                const F *A11 = a, *A21 = a + mh, *A12 = a + kh * lda, *A22 = a + mh + kh * lda;
                const F *B11 = bb, *B21 = bb + kh, *B12 = bb + nh * ldb, *B22 = bb + kh + nh * ldb;
                F *S[8];
                for (int t = 0; t < 4; t++)
                    S[t] = b.alloc((size_t)mh * kh);
                for (int t = 4; t < 8; t++)
                    S[t] = b.alloc((size_t)kh * nh);
                if (!b.failed()) {
                    b.addsub(mh, kh, S[0], mh, A21, lda, A22, lda, 0); // S1 = A21 + A22
                    b.addsub(mh, kh, S[1], mh, S[0], mh, A11, lda, 1); // S2 = S1 - A11
                    b.addsub(mh, kh, S[2], mh, A11, lda, A21, lda, 1); // S3 = A11 - A21
                    b.addsub(mh, kh, S[3], mh, A12, lda, S[1], mh, 1); // S4 = A12 - S2
                    b.addsub(kh, nh, S[4], kh, B12, ldb, B11, ldb, 1); // S5 = B12 - B11
                    b.addsub(kh, nh, S[5], kh, B22, ldb, S[4], kh, 1); // S6 = B22 - S5
                    b.addsub(kh, nh, S[6], kh, B22, ldb, B12, ldb, 1); // S7 = B22 - B12
                    b.addsub(kh, nh, S[7], kh, S[5], kh, B21, ldb, 1); // S8 = S6 - B21
                }
                // M1 = S2 S6, M2 = A11 B11, M3 = A12 B21, M4 = S3 S7,
                // M5 = S1 S5, M6 = S4 B22, M7 = A22 S8
                const F *const opA[7] = {S[1], A11, A12, S[2], S[0], S[3], A22};
                const long ldA[7] = {mh, lda, lda, mh, mh, mh, lda};
                const F *const opB[7] = {S[5], B11, B21, S[6], S[4], B22, S[7]};
                const long ldB[7] = {kh, ldb, ldb, kh, kh, ldb, kh};
                // All of this device's products first, then the copies to
                // device 0 (a peer copy waits for the work queued on both
                // devices, so copying each product at once would stall
                // this device behind device 0).
                F *out[7] = {};
                for (int t = d; t < 7 && !b.failed(); t += G) {
                    if (d == 0)
                        out[t] = M[t];
                    else if ((e = cudaMalloc((void **)&out[t], bytesM)) != cudaSuccess)
                        break;
                    if (!b.failed())
                        rec(b, mh, nh, kh, opA[t], ldA[t], opB[t], ldB[t], out[t], mh);
                }
                for (int t = 0; t < 8; t++)
                    b.release(S[t]);
                if (d != 0)
                    for (int t = d; t < 7; t += G) {
                        if (e == cudaSuccess && !b.failed() && out[t])
                            e = cudaMemcpyPeer(M[t], 0, out[t], d, bytesM);
                        if (e == cudaSuccess)
                            e = cudaDeviceSynchronize();
                        cudaFree(out[t]);
                    }
                if (e == cudaSuccess)
                    e = cudaDeviceSynchronize();
                rc[(size_t)d] = e != cudaSuccess ? (int)e : b.failed() ? 1 : 0;
            });
        for (auto &w : workers)
            w.join();
    }
    err = cudaSetDevice(0);
    for (int r : rc)
        if (!err && r)
            err = r;
    F *C11 = C, *C21 = C + mh, *C12 = C + nh * ldc, *C22 = C + mh + nh * ldc;
    if (!err && !b0.failed()) {
        b0.addsub(mh, nh, T[0], mh, M[0], mh, M[1], mh, 0); // T1 = M1 + M2
        b0.addsub(mh, nh, T[1], mh, T[0], mh, M[3], mh, 0); // T2 = T1 + M4
        b0.addsub(mh, nh, C11, ldc, M[1], mh, M[2], mh, 0); // C11 = M2 + M3
        b0.addsub(mh, nh, C12, ldc, T[0], mh, M[4], mh, 0); // C12 = T1 + M5
        b0.addsub(mh, nh, C12, ldc, C12, ldc, M[5], mh, 0); //     + M6
        b0.addsub(mh, nh, C21, ldc, T[1], mh, M[6], mh, 1); // C21 = T2 - M7
        b0.addsub(mh, nh, C22, ldc, T[1], mh, M[4], mh, 0); // C22 = T2 + M5
    }
    if (!err)
        err = (int)cudaDeviceSynchronize();
    for (int t = 0; t < 7; t++)
        cudaFree(M[t]);
    for (int t = 0; t < 2; t++)
        b0.release(T[t]);
    // The peeled parts, as winograd() after its even part.
    const F *a = A[0], *bb = B[0];
    if (!err && !b0.failed() && (k & 1)) {
        F *Tp = b0.alloc((size_t)me * ne);
        if (!b0.failed()) {
            b0.leaf(me, ne, 1, Tp, me, a + (k - 1) * lda, lda, bb + (k - 1), ldb);
            b0.addsub(me, ne, C, ldc, C, ldc, Tp, me, 0);
        }
        b0.release(Tp);
    }
    if (!err && !b0.failed() && (n & 1))
        b0.leaf(m, 1, k, C + (n - 1) * ldc, ldc, a, lda, bb + (n - 1) * ldb, ldb);
    if (!err && !b0.failed() && (m & 1))
        b0.leaf(1, ne, k, C + (m - 1), ldc, a + (m - 1), lda, bb, ldb);
    return err ? err : b0.failed() ? 1 : 0;
}

} // namespace mplapack_gemm

#endif
