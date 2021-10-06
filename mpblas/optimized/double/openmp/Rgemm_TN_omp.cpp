/*
 * Copyright (c) 2008-2021
 *	Nakata, Maho
 * 	All rights reserved.
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

#include <mpblas_double.h>
#include <omp.h>

void Rgemm_TN_omp(mplapackint m, mplapackint n, mplapackint k, double alpha, double *A, mplapackint lda, double *B, mplapackint ldb, double beta, double *C, mplapackint ldc) {
    mplapackint i, j, l;
    double temp;

    // Form C := alpha*A'*B + beta*C.
    for (j = 0; j < n; j++) {
        if (beta == 0.0) {
            for (i = 0; i < m; i++) {
                C[i + j * ldc] = 0.0;
            }
        } else if (beta != 1.0) {
            for (i = 0; i < m; i++) {
                C[i + j * ldc] = beta * C[i + j * ldc];
            }
        }
    }
    // main loop
    mplapackint p, q, r;
    mplapackint qq, rr;
    mplapackint Bq = 16, Br = 16;

    for (qq = 0; qq < n; qq = qq + Bq) {
        for (rr = 0; rr < k; rr = rr + Br) {
            for (p = 0; p < m; p++) {
                for (q = qq; q < std::min(qq + Bq, n); q++) {
                    temp = 0.0;
                    for (r = rr; r < std::min(rr + Br, k); r++) {
                        temp = temp + A[r + p * lda] * B[r + q * ldb];
                    }
                    C[p + q * ldc] += alpha * temp;
                }
            }
        }
    }
}
