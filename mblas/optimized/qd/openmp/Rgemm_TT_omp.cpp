/*
 * Copyright (c) 2008-2012
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Rgemm_TT.cpp,v 1.1 2010/12/28 06:13:53 nakatamaho Exp $
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
#include <mblas_qd.h>

void Rgemm_TT_omp(mpackint m, mpackint n, mpackint k, qd_real alpha, qd_real *A, mpackint lda, qd_real *B, mpackint ldb, qd_real beta,
	      qd_real *C, mpackint ldc)
{
//Form  C := alpha*A'*B' + beta*C.
    mpackint i, j, l;
    qd_real temp;
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
//main loop
#ifdef _OPENMP
#pragma omp parallel for private(i, j, l, temp)
#endif
    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            temp = 0.0;
            for (l = 0; l < k; l++) {
                temp += A[l + i * lda] * B[j + l * ldb];
	    }
	    C[i + j * ldc] += alpha * temp;
	}
    }
    return;
}
