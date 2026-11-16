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
 * Argument checks of the Rgemm entry points <backend>_cuda / _opencl
 * (AGENTS.md, "Accelerator entry points"): the checks of Rgemm, without
 * Mxerbla.  gemm_entry_args() returns false for invalid arguments and for
 * the quick returns of Rgemm (m == 0, n == 0, or (alpha == 0 or k == 0) and
 * beta == 1), which the entry points leave to the caller.
 */

#ifndef _MPLAPACK_RGEMM_ENTRY_COMMON_H_
#define _MPLAPACK_RGEMM_ENTRY_COMMON_H_

namespace mplapack_gemm {

// 'N', 'T' or 'C' (either case), as Mlsame; 0 otherwise.
inline char gemm_entry_trans(const char *t) {
    if (!t)
        return 0;
    const char c = (char)(*t >= 'a' && *t <= 'z' ? *t - 'a' + 'A' : *t);
    return (c == 'N' || c == 'T' || c == 'C') ? c : 0;
}

inline bool gemm_entry_args(const char *transa, const char *transb, long m, long n, long k, long lda, long ldb, long ldc, bool alpha_zero, bool beta_one, bool &nota, bool &notb) {
    const char ta = gemm_entry_trans(transa), tb = gemm_entry_trans(transb);
    if (!ta || !tb || m < 0 || n < 0 || k < 0)
        return false;
    nota = ta == 'N';
    notb = tb == 'N';
    const long nrowa = nota ? m : k, nrowb = notb ? k : n;
    if (lda < (nrowa > 1 ? nrowa : 1) || ldb < (nrowb > 1 ? nrowb : 1) || ldc < (m > 1 ? m : 1))
        return false;
    if (m == 0 || n == 0 || ((alpha_zero || k == 0) && beta_one))
        return false;
    return true;
}

} // namespace mplapack_gemm

#endif
