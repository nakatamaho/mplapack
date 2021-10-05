/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rgemm.cpp,v 1.1 2010/12/28 06:13:53 nakatamaho Exp $
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

void Rgemm_NN(mplapackint m, mplapackint n, mplapackint k, double alpha, double *A, mplapackint lda, double *B, mplapackint ldb, double beta, double *C, mplapackint ldc);
void Rgemm_TN(mplapackint m, mplapackint n, mplapackint k, double alpha, double *A, mplapackint lda, double *B, mplapackint ldb, double beta, double *C, mplapackint ldc);
void Rgemm_NT(mplapackint m, mplapackint n, mplapackint k, double alpha, double *A, mplapackint lda, double *B, mplapackint ldb, double beta, double *C, mplapackint ldc);
void Rgemm_TT(mplapackint m, mplapackint n, mplapackint k, double alpha, double *A, mplapackint lda, double *B, mplapackint ldb, double beta, double *C, mplapackint ldc);

void Rgemm(const char *transa, const char *transb, mplapackint m, mplapackint n, mplapackint k, double alpha, double *A, mplapackint lda, double *B, mplapackint ldb, double beta, double *C, mplapackint ldc) {
    mplapackint i, j, l, nota, notb, nrowa, ncola, nrowb, info;
    double temp;
    double Zero = 0.0, One = 1.0;

    nota = Mlsame_double(transa, "N");
    notb = Mlsame_double(transb, "N");
    if (nota) {
        nrowa = m;
        ncola = k;
    } else {
        nrowa = k;
        ncola = m;
    }
    if (notb) {
        nrowb = k;
    } else {
        nrowb = n;
    }
    // Test the input parameters.
    info = 0;
    if (!nota && (!Mlsame_double(transa, "C")) && (!Mlsame_double(transa, "T")))
        info = 1;
    else if (!notb && (!Mlsame_double(transb, "C")) && (!Mlsame_double(transb, "T")))
        info = 2;
    else if (m < 0)
        info = 3;
    else if (n < 0)
        info = 4;
    else if (k < 0)
        info = 5;
    else if (lda < std::max((mplapackint)1, nrowa))
        info = 8;
    else if (ldb < std::max((mplapackint)1, nrowb))
        info = 10;
    else if (ldc < std::max((mplapackint)1, m))
        info = 13;
    if (info != 0) {
        Mxerbla_double("Rgemm ", info);
        return;
    }
    // Quick return if possible.
    if ((m == 0) || (n == 0) || (((alpha == Zero) || (k == 0)) && (beta == One)))
        return;
    // And when alpha == 0.0
    if (alpha == Zero) {
        if (beta == Zero) {
            for (j = 0; j < n; j++) {
                for (i = 0; i < m; i++) {
                    C[i + j * ldc] = Zero;
                }
            }
        } else {
            for (j = 0; j < n; j++) {
                for (i = 0; i < m; i++) {
                    C[i + j * ldc] = beta * C[i + j * ldc];
                }
            }
        }
        return;
    }
    // Start the operations.
    if (notb) {
        if (nota) {
            // Form C := alpha*A*B + beta*C.
            Rgemm_NN(m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
        } else {
            // Form  C := alpha*A'*B + beta*C.
            Rgemm_TN(m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
        }
    } else {
        if (nota) {
            // Form  C := alpha*A*B' + beta*C.
            Rgemm_NT(m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
        } else {
            // Form  C := alpha*A'*B' + beta*C.
            Rgemm_TT(m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
        }
    }
    return;
}
