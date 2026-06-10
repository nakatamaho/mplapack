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

// Derived from LAPACK routine ZLARFT.
// Original LAPACK authors:
//   Univ. of Tennessee
//   Univ. of California Berkeley
//   Univ. of Colorado Denver
//   NAG Ltd.

#include <mpblas.h>
#include <mplapack.h>

void Clarft(const char *direct, const char *storev, INTEGER const n, INTEGER const k, COMPLEX *v, INTEGER const ldv, COMPLEX *tau, COMPLEX *t, INTEGER const ldt) {
    //
    // Quick return if possible
    //
    if (n == 0 || k == 0) {
        return;
    }
    //
    // Base case
    //
    if (n == 1 || k == 1) {
        t[0] = tau[1 - 1];
        return;
    }
    //
    // Beginning of executable statements
    //
    INTEGER l = k / 2;
    //
    // Determine what kind of Q we need to compute
    // We assume that if the user doesn't provide 'F' for DIRECT,
    // then they meant to provide 'B' and if they don't provide
    // 'C' for STOREV, then they meant to provide 'R'
    //
    bool dirf = Mlsame(direct, "F");
    bool colv = Mlsame(storev, "C");
    //
    // QR happens when we have forward direction in column storage
    //
    bool qr = dirf && colv;
    //
    // LQ happens when we have forward direction in row storage
    //
    bool lq = dirf && (!colv);
    //
    // QL happens when we have backward direction in column storage
    //
    bool ql = (!dirf) && colv;
    //
    // The last case is RQ. Due to how we structured this, if the
    // above 3 are false, then RQ must be true, so we never store
    // this
    // RQ happens when we have backward direction in row storage
    // RQ = (.NOT.DIRF).AND.(.NOT.COLV)
    //
    INTEGER j = 0;
    INTEGER i = 0;
    const COMPLEX one = 1.0;
    const COMPLEX neg_one = -1.0;
    if (qr) {
        //
        // Break V apart into 6 components
        //
        // V = |---------------|
        // |V_{1,1} 0      |
        // |V_{2,1} V_{2,2}|
        // |V_{3,1} V_{3,2}|
        // |---------------|
        //
        // V_{1,1}\in\C^{l,l}      unit lower triangular
        // V_{2,1}\in\C^{k-l,l}    rectangular
        // V_{3,1}\in\C^{n-k,l}    rectangular
        //
        // V_{2,2}\in\C^{k-l,k-l}  unit lower triangular
        // V_{3,2}\in\C^{n-k,k-l}  rectangular
        //
        // We will construct the T matrix
        // T = |---------------|
        // |T_{1,1} T_{1,2}|
        // |0       T_{2,2}|
        // |---------------|
        //
        // T is the triangular factor obtained from block reflectors.
        // To motivate the structure, assume we have already computed T_{1,1}
        // and T_{2,2}. Then collect the associated reflectors in V_1 and V_2
        //
        // T_{1,1}\in\C^{l, l}     upper triangular
        // T_{2,2}\in\C^{k-l, k-l} upper triangular
        // T_{1,2}\in\C^{l, k-l}   rectangular
        //
        // Where l = floor(k/2)
        //
        // Then, consider the product:
        //
        // (I - V_1*T_{1,1}*V_1')*(I - V_2*T_{2,2}*V_2')
        // = I - V_1*T_{1,1}*V_1' - V_2*T_{2,2}*V_2' + V_1*T_{1,1}*V_1'*V_2*T_{2,2}*V_2'
        //
        // Define T_{1,2} = -T_{1,1}*V_1'*V_2*T_{2,2}
        //
        // Then, we can define the matrix V as
        // V = |-------|
        // |V_1 V_2|
        // |-------|
        //
        // So, our product is equivalent to the matrix product
        // I - V*T*V'
        // This means, we can compute T_{1,1} and T_{2,2}, then use this information
        // to compute T_{1,2}
        //
        // Compute T_{1,1} recursively
        //
        Clarft(direct, storev, n, l, v, ldv, tau, t, ldt);
        //
        // Compute T_{2,2} recursively
        //
        Clarft(direct, storev, n - l, k - l, &v[((l + 1) - 1) + ((l + 1) - 1) * ldv], ldv, &tau[(l + 1) - 1], &t[((l + 1) - 1) + ((l + 1) - 1) * ldt], ldt);
        //
        // Compute T_{1,2}
        // T_{1,2} = V_{2,1}'
        //
        for (j = 1; j <= l; j = j + 1) {
            for (i = 1; i <= k - l; i = i + 1) {
                t[(j - 1) + ((l + i) - 1) * ldt] = conj(v[((l + i) - 1) + (j - 1) * ldv]);
            }
        }
        //
        // T_{1,2} = T_{1,2}*V_{2,2}
        //
        Ctrmm("Right", "Lower", "No transpose", "Unit", l, k - l, one, &v[((l + 1) - 1) + ((l + 1) - 1) * ldv], ldv, &t[((l + 1) - 1) * ldt], ldt);
        //
        // T_{1,2} = V_{3,1}'*V_{3,2} + T_{1,2}
        // Note: We assume K <= N, and GEMM will do nothing if N=K
        //
        Cgemm("Conjugate", "No transpose", l, k - l, n - k, one, &v[((k + 1) - 1)], ldv, &v[((k + 1) - 1) + ((l + 1) - 1) * ldv], ldv, one, &t[((l + 1) - 1) * ldt], ldt);
        //
        // At this point, we have that T_{1,2} = V_1'*V_2
        // All that is left is to pre and post multiply by -T_{1,1} and T_{2,2}
        // respectively.
        //
        // T_{1,2} = -T_{1,1}*T_{1,2}
        //
        Ctrmm("Left", "Upper", "No transpose", "Non-unit", l, k - l, neg_one, t, ldt, &t[((l + 1) - 1) * ldt], ldt);
        //
        // T_{1,2} = T_{1,2}*T_{2,2}
        //
        Ctrmm("Right", "Upper", "No transpose", "Non-unit", l, k - l, one, &t[((l + 1) - 1) + ((l + 1) - 1) * ldt], ldt, &t[((l + 1) - 1) * ldt], ldt);
        //
    } else if (lq) {
        //
        // Break V apart into 6 components
        //
        // V = |----------------------|
        // |V_{1,1} V_{1,2} V{1,3}|
        // |0       V_{2,2} V{2,3}|
        // |----------------------|
        //
        // V_{1,1}\in\C^{l,l}      unit upper triangular
        // V_{1,2}\in\C^{l,k-l}    rectangular
        // V_{1,3}\in\C^{l,n-k}    rectangular
        //
        // V_{2,2}\in\C^{k-l,k-l}  unit upper triangular
        // V_{2,3}\in\C^{k-l,n-k}  rectangular
        //
        // Where l = floor(k/2)
        //
        // We will construct the T matrix
        // T = |---------------|
        // |T_{1,1} T_{1,2}|
        // |0       T_{2,2}|
        // |---------------|
        //
        // T is the triangular factor obtained from block reflectors.
        // To motivate the structure, assume we have already computed T_{1,1}
        // and T_{2,2}. Then collect the associated reflectors in V_1 and V_2
        //
        // T_{1,1}\in\C^{l, l}         upper triangular
        // T_{2,2}\in\C^{k-l, k-l}     upper triangular
        // T_{1,2}\in\C^{l, k-l}       rectangular
        //
        // Then, consider the product:
        //
        // (I - V_1'*T_{1,1}*V_1)*(I - V_2'*T_{2,2}*V_2)
        // = I - V_1'*T_{1,1}*V_1 - V_2'*T_{2,2}*V_2 + V_1'*T_{1,1}*V_1*V_2'*T_{2,2}*V_2
        //
        // Define T_{1,2} = -T_{1,1}*V_1*V_2'*T_{2,2}
        //
        // Then, we can define the matrix V as
        // V = |---|
        // |V_1|
        // |V_2|
        // |---|
        //
        // So, our product is equivalent to the matrix product
        // I - V'*T*V
        // This means, we can compute T_{1,1} and T_{2,2}, then use this information
        // to compute T_{1,2}
        //
        // Compute T_{1,1} recursively
        //
        Clarft(direct, storev, n, l, v, ldv, tau, t, ldt);
        //
        // Compute T_{2,2} recursively
        //
        Clarft(direct, storev, n - l, k - l, &v[((l + 1) - 1) + ((l + 1) - 1) * ldv], ldv, &tau[(l + 1) - 1], &t[((l + 1) - 1) + ((l + 1) - 1) * ldt], ldt);
        //
        // Compute T_{1,2}
        // T_{1,2} = V_{1,2}
        //
        Clacpy("All", l, k - l, &v[((l + 1) - 1) * ldv], ldv, &t[((l + 1) - 1) * ldt], ldt);
        //
        // T_{1,2} = T_{1,2}*V_{2,2}'
        //
        Ctrmm("Right", "Upper", "Conjugate", "Unit", l, k - l, one, &v[((l + 1) - 1) + ((l + 1) - 1) * ldv], ldv, &t[((l + 1) - 1) * ldt], ldt);
        //
        // T_{1,2} = V_{1,3}*V_{2,3}' + T_{1,2}
        // Note: We assume K <= N, and GEMM will do nothing if N=K
        //
        Cgemm("No transpose", "Conjugate", l, k - l, n - k, one, &v[((k + 1) - 1) * ldv], ldv, &v[((l + 1) - 1) + ((k + 1) - 1) * ldv], ldv, one, &t[((l + 1) - 1) * ldt], ldt);
        //
        // At this point, we have that T_{1,2} = V_1*V_2'
        // All that is left is to pre and post multiply by -T_{1,1} and T_{2,2}
        // respectively.
        //
        // T_{1,2} = -T_{1,1}*T_{1,2}
        //
        Ctrmm("Left", "Upper", "No transpose", "Non-unit", l, k - l, neg_one, t, ldt, &t[((l + 1) - 1) * ldt], ldt);
        //
        // T_{1,2} = T_{1,2}*T_{2,2}
        //
        Ctrmm("Right", "Upper", "No transpose", "Non-unit", l, k - l, one, &t[((l + 1) - 1) + ((l + 1) - 1) * ldt], ldt, &t[((l + 1) - 1) * ldt], ldt);
    } else if (ql) {
        //
        // Break V apart into 6 components
        //
        // V = |---------------|
        // |V_{1,1} V_{1,2}|
        // |V_{2,1} V_{2,2}|
        // |0       V_{3,2}|
        // |---------------|
        //
        // V_{1,1}\in\C^{n-k,k-l}  rectangular
        // V_{2,1}\in\C^{k-l,k-l}  unit upper triangular
        //
        // V_{1,2}\in\C^{n-k,l}    rectangular
        // V_{2,2}\in\C^{k-l,l}    rectangular
        // V_{3,2}\in\C^{l,l}      unit upper triangular
        //
        // We will construct the T matrix
        // T = |---------------|
        // |T_{1,1} 0      |
        // |T_{2,1} T_{2,2}|
        // |---------------|
        //
        // T is the triangular factor obtained from block reflectors.
        // To motivate the structure, assume we have already computed T_{1,1}
        // and T_{2,2}. Then collect the associated reflectors in V_1 and V_2
        //
        // T_{1,1}\in\C^{k-l, k-l} non-unit lower triangular
        // T_{2,2}\in\C^{l, l}     non-unit lower triangular
        // T_{2,1}\in\C^{k-l, l}   rectangular
        //
        // Where l = floor(k/2)
        //
        // Then, consider the product:
        //
        // (I - V_2*T_{2,2}*V_2')*(I - V_1*T_{1,1}*V_1')
        // = I - V_2*T_{2,2}*V_2' - V_1*T_{1,1}*V_1' + V_2*T_{2,2}*V_2'*V_1*T_{1,1}*V_1'
        //
        // Define T_{2,1} = -T_{2,2}*V_2'*V_1*T_{1,1}
        //
        // Then, we can define the matrix V as
        // V = |-------|
        // |V_1 V_2|
        // |-------|
        //
        // So, our product is equivalent to the matrix product
        // I - V*T*V'
        // This means, we can compute T_{1,1} and T_{2,2}, then use this information
        // to compute T_{2,1}
        //
        // Compute T_{1,1} recursively
        //
        Clarft(direct, storev, n - l, k - l, v, ldv, tau, t, ldt);
        //
        // Compute T_{2,2} recursively
        //
        Clarft(direct, storev, n, l, &v[((k - l + 1) - 1) * ldv], ldv, &tau[(k - l + 1) - 1], &t[((k - l + 1) - 1) + ((k - l + 1) - 1) * ldt], ldt);
        //
        // Compute T_{2,1}
        // T_{2,1} = V_{2,2}'
        //
        for (j = 1; j <= k - l; j = j + 1) {
            for (i = 1; i <= l; i = i + 1) {
                t[((k - l + i) - 1) + (j - 1) * ldt] = conj(v[((n - k + j) - 1) + ((k - l + i) - 1) * ldv]);
            }
        }
        //
        // T_{2,1} = T_{2,1}*V_{2,1}
        //
        Ctrmm("Right", "Upper", "No transpose", "Unit", l, k - l, one, &v[((n - k + 1) - 1)], ldv, &t[((k - l + 1) - 1)], ldt);
        //
        // T_{2,1} = V_{2,2}'*V_{2,1} + T_{2,1}
        // Note: We assume K <= N, and GEMM will do nothing if N=K
        //
        Cgemm("Conjugate", "No transpose", l, k - l, n - k, one, &v[((k - l + 1) - 1) * ldv], ldv, v, ldv, one, &t[((k - l + 1) - 1)], ldt);
        //
        // At this point, we have that T_{2,1} = V_2'*V_1
        // All that is left is to pre and post multiply by -T_{2,2} and T_{1,1}
        // respectively.
        //
        // T_{2,1} = -T_{2,2}*T_{2,1}
        //
        Ctrmm("Left", "Lower", "No transpose", "Non-unit", l, k - l, neg_one, &t[((k - l + 1) - 1) + ((k - l + 1) - 1) * ldt], ldt, &t[((k - l + 1) - 1)], ldt);
        //
        // T_{2,1} = T_{2,1}*T_{1,1}
        //
        Ctrmm("Right", "Lower", "No transpose", "Non-unit", l, k - l, one, t, ldt, &t[((k - l + 1) - 1)], ldt);
    } else {
        //
        // Else means RQ case
        //
        // Break V apart into 6 components
        //
        // V = |-----------------------|
        // |V_{1,1} V_{1,2} 0      |
        // |V_{2,1} V_{2,2} V_{2,3}|
        // |-----------------------|
        //
        // V_{1,1}\in\C^{k-l,n-k}  rectangular
        // V_{1,2}\in\C^{k-l,k-l}  unit lower triangular
        //
        // V_{2,1}\in\C^{l,n-k}    rectangular
        // V_{2,2}\in\C^{l,k-l}    rectangular
        // V_{2,3}\in\C^{l,l}      unit lower triangular
        //
        // We will construct the T matrix
        // T = |---------------|
        // |T_{1,1} 0      |
        // |T_{2,1} T_{2,2}|
        // |---------------|
        //
        // T is the triangular factor obtained from block reflectors.
        // To motivate the structure, assume we have already computed T_{1,1}
        // and T_{2,2}. Then collect the associated reflectors in V_1 and V_2
        //
        // T_{1,1}\in\C^{k-l, k-l} non-unit lower triangular
        // T_{2,2}\in\C^{l, l}     non-unit lower triangular
        // T_{2,1}\in\C^{k-l, l}   rectangular
        //
        // Where l = floor(k/2)
        //
        // Then, consider the product:
        //
        // (I - V_2'*T_{2,2}*V_2)*(I - V_1'*T_{1,1}*V_1)
        // = I - V_2'*T_{2,2}*V_2 - V_1'*T_{1,1}*V_1 + V_2'*T_{2,2}*V_2*V_1'*T_{1,1}*V_1
        //
        // Define T_{2,1} = -T_{2,2}*V_2*V_1'*T_{1,1}
        //
        // Then, we can define the matrix V as
        // V = |---|
        // |V_1|
        // |V_2|
        // |---|
        //
        // So, our product is equivalent to the matrix product
        // I - V'*T*V
        // This means, we can compute T_{1,1} and T_{2,2}, then use this information
        // to compute T_{2,1}
        //
        // Compute T_{1,1} recursively
        //
        Clarft(direct, storev, n - l, k - l, v, ldv, tau, t, ldt);
        //
        // Compute T_{2,2} recursively
        //
        Clarft(direct, storev, n, l, &v[((k - l + 1) - 1)], ldv, &tau[(k - l + 1) - 1], &t[((k - l + 1) - 1) + ((k - l + 1) - 1) * ldt], ldt);
        //
        // Compute T_{2,1}
        // T_{2,1} = V_{2,2}
        //
        Clacpy("All", l, k - l, &v[((k - l + 1) - 1) + ((n - k + 1) - 1) * ldv], ldv, &t[((k - l + 1) - 1)], ldt);
        //
        // T_{2,1} = T_{2,1}*V_{1,2}'
        //
        Ctrmm("Right", "Lower", "Conjugate", "Unit", l, k - l, one, &v[((n - k + 1) - 1) * ldv], ldv, &t[((k - l + 1) - 1)], ldt);
        //
        // T_{2,1} = V_{2,1}*V_{1,1}' + T_{2,1}
        // Note: We assume K <= N, and GEMM will do nothing if N=K
        //
        Cgemm("No transpose", "Conjugate", l, k - l, n - k, one, &v[((k - l + 1) - 1)], ldv, v, ldv, one, &t[((k - l + 1) - 1)], ldt);
        //
        // At this point, we have that T_{2,1} = V_2*V_1'
        // All that is left is to pre and post multiply by -T_{2,2} and T_{1,1}
        // respectively.
        //
        // T_{2,1} = -T_{2,2}*T_{2,1}
        //
        Ctrmm("Left", "Lower", "No tranpose", "Non-unit", l, k - l, neg_one, &t[((k - l + 1) - 1) + ((k - l + 1) - 1) * ldt], ldt, &t[((k - l + 1) - 1)], ldt);
        //
        // T_{2,1} = T_{2,1}*T_{1,1}
        //
        Ctrmm("Right", "Lower", "No tranpose", "Non-unit", l, k - l, one, t, ldt, &t[((k - l + 1) - 1)], ldt);
    }
}
