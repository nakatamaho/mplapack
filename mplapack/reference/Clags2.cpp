/*
 * Copyright (c) 2008-2021
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

#include <mpblas.h>
#include <mplapack.h>

inline REAL abs1(COMPLEX t) { return abs(t.real()) + abs(t.imag()); }

void Clags2(bool const upper, REAL const a1, COMPLEX const a2, REAL const a3, REAL const b1, COMPLEX const b2, REAL const b3, REAL &csu, COMPLEX &snu, REAL &csv, COMPLEX &snv, REAL &csq, COMPLEX &snq) {
    //
    COMPLEX t = 0.0;
    //
    REAL a = 0.0;
    REAL d = 0.0;
    COMPLEX b = 0.0;
    REAL fb = 0.0;
    const REAL one = 1.0;
    COMPLEX d1 = 0.0;
    const REAL zero = 0.0;
    REAL s1 = 0.0;
    REAL s2 = 0.0;
    REAL snr = 0.0;
    REAL csr = 0.0;
    REAL snl = 0.0;
    REAL csl = 0.0;
    REAL ua11r = 0.0;
    COMPLEX ua12 = 0.0;
    REAL vb11r = 0.0;
    COMPLEX vb12 = 0.0;
    REAL aua12 = 0.0;
    REAL avb12 = 0.0;
    COMPLEX r = 0.0;
    COMPLEX ua21 = 0.0;
    COMPLEX ua22 = 0.0;
    COMPLEX vb21 = 0.0;
    COMPLEX vb22 = 0.0;
    REAL aua22 = 0.0;
    REAL avb22 = 0.0;
    COMPLEX c = 0.0;
    REAL fc = 0.0;
    REAL ua22r = 0.0;
    REAL vb22r = 0.0;
    REAL aua21 = 0.0;
    REAL avb21 = 0.0;
    COMPLEX ua11 = 0.0;
    COMPLEX vb11 = 0.0;
    REAL aua11 = 0.0;
    REAL avb11 = 0.0;
    if (upper) {
        //
        //        Input matrices A and B are upper triangular matrices
        //
        //        Form matrix C = A*adj(B) = ( a b )
        //                                   ( 0 d )
        //
        a = a1 * b3;
        d = a3 * b1;
        b = a2 * b1 - a1 * b2;
        fb = abs(b);
        //
        //        Transform complex 2-by-2 matrix C to real matrix by unitary
        //        diagonal matrix diag(1,D1).
        //
        d1 = one;
        if (fb != zero) {
            d1 = b / fb;
        }
        //
        //        The SVD of real 2 by 2 triangular C
        //
        //         ( CSL -SNL )*( A B )*(  CSR  SNR ) = ( R 0 )
        //         ( SNL  CSL ) ( 0 D ) ( -SNR  CSR )   ( 0 T )
        //
        Rlasv2(a, fb, d, s1, s2, snr, csr, snl, csl);
        //
        if (abs(csl) >= abs(snl) || abs(csr) >= abs(snr)) {
            //
            //           Compute the (1,1) and (1,2) elements of U**H *A and V**H *B,
            //           and (1,2) element of |U|**H *|A| and |V|**H *|B|.
            //
            ua11r = csl * a1;
            ua12 = csl * a2 + d1 * snl * a3;
            //
            vb11r = csr * b1;
            vb12 = csr * b2 + d1 * snr * b3;
            //
            aua12 = abs(csl) * abs1(a2) + abs(snl) * abs(a3);
            avb12 = abs(csr) * abs1(b2) + abs(snr) * abs(b3);
            //
            //           zero (1,2) elements of U**H *A and V**H *B
            //
            if ((abs(ua11r) + abs1(ua12)) == zero) {
                Clartg(-COMPLEX(vb11r), conj(vb12), csq, snq, r);
            } else if ((abs(vb11r) + abs1(vb12)) == zero) {
                Clartg(-COMPLEX(ua11r), conj(ua12), csq, snq, r);
            } else if (aua12 / (abs(ua11r) + abs1(ua12)) <= avb12 / (abs(vb11r) + abs1(vb12))) {
                Clartg(-COMPLEX(ua11r), conj(ua12), csq, snq, r);
            } else {
                Clartg(-COMPLEX(vb11r), conj(vb12), csq, snq, r);
            }
            //
            csu = csl;
            snu = -d1 * snl;
            csv = csr;
            snv = -d1 * snr;
            //
        } else {
            //
            //           Compute the (2,1) and (2,2) elements of U**H *A and V**H *B,
            //           and (2,2) element of |U|**H *|A| and |V|**H *|B|.
            //
            ua21 = -conj(d1) * snl * a1;
            ua22 = -conj(d1) * snl * a2 + csl * a3;
            //
            vb21 = -conj(d1) * snr * b1;
            vb22 = -conj(d1) * snr * b2 + csr * b3;
            //
            aua22 = abs(snl) * abs1(a2) + abs(csl) * abs(a3);
            avb22 = abs(snr) * abs1(b2) + abs(csr) * abs(b3);
            //
            //           zero (2,2) elements of U**H *A and V**H *B, and then swap.
            //
            if ((abs1(ua21) + abs1(ua22)) == zero) {
                Clartg(-conj(vb21), conj(vb22), csq, snq, r);
            } else if ((abs1(vb21) + abs(vb22)) == zero) {
                Clartg(-conj(ua21), conj(ua22), csq, snq, r);
            } else if (aua22 / (abs1(ua21) + abs1(ua22)) <= avb22 / (abs1(vb21) + abs1(vb22))) {
                Clartg(-conj(ua21), conj(ua22), csq, snq, r);
            } else {
                Clartg(-conj(vb21), conj(vb22), csq, snq, r);
            }
            //
            csu = snl;
            snu = d1 * csl;
            csv = snr;
            snv = d1 * csr;
            //
        }
        //
    } else {
        //
        //        Input matrices A and B are lower triangular matrices
        //
        //        Form matrix C = A*adj(B) = ( a 0 )
        //                                   ( c d )
        //
        a = a1 * b3;
        d = a3 * b1;
        c = a2 * b3 - a3 * b2;
        fc = abs(c);
        //
        //        Transform complex 2-by-2 matrix C to real matrix by unitary
        //        diagonal matrix diag(d1,1).
        //
        d1 = one;
        if (fc != zero) {
            d1 = c / fc;
        }
        //
        //        The SVD of real 2 by 2 triangular C
        //
        //         ( CSL -SNL )*( A 0 )*(  CSR  SNR ) = ( R 0 )
        //         ( SNL  CSL ) ( C D ) ( -SNR  CSR )   ( 0 T )
        //
        Rlasv2(a, fc, d, s1, s2, snr, csr, snl, csl);
        //
        if (abs(csr) >= abs(snr) || abs(csl) >= abs(snl)) {
            //
            //           Compute the (2,1) and (2,2) elements of U**H *A and V**H *B,
            //           and (2,1) element of |U|**H *|A| and |V|**H *|B|.
            //
            ua21 = -d1 * snr * a1 + csr * a2;
            ua22r = csr * a3;
            //
            vb21 = -d1 * snl * b1 + csl * b2;
            vb22r = csl * b3;
            //
            aua21 = abs(snr) * abs(a1) + abs(csr) * abs1(a2);
            avb21 = abs(snl) * abs(b1) + abs(csl) * abs1(b2);
            //
            //           zero (2,1) elements of U**H *A and V**H *B.
            //
            if ((abs1(ua21) + abs(ua22r)) == zero) {
                Clartg(COMPLEX(vb22r), vb21, csq, snq, r);
            } else if ((abs1(vb21) + abs(vb22r)) == zero) {
                Clartg(COMPLEX(ua22r), ua21, csq, snq, r);
            } else if (aua21 / (abs1(ua21) + abs(ua22r)) <= avb21 / (abs1(vb21) + abs(vb22r))) {
                Clartg(COMPLEX(ua22r), ua21, csq, snq, r);
            } else {
                Clartg(COMPLEX(vb22r), vb21, csq, snq, r);
            }
            //
            csu = csr;
            snu = -conj(d1) * snr;
            csv = csl;
            snv = -conj(d1) * snl;
            //
        } else {
            //
            //           Compute the (1,1) and (1,2) elements of U**H *A and V**H *B,
            //           and (1,1) element of |U|**H *|A| and |V|**H *|B|.
            //
            ua11 = csr * a1 + conj(d1) * snr * a2;
            ua12 = conj(d1) * snr * a3;
            //
            vb11 = csl * b1 + conj(d1) * snl * b2;
            vb12 = conj(d1) * snl * b3;
            //
            aua11 = abs(csr) * abs(a1) + abs(snr) * abs1(a2);
            avb11 = abs(csl) * abs(b1) + abs(snl) * abs1(b2);
            //
            //           zero (1,1) elements of U**H *A and V**H *B, and then swap.
            //
            if ((abs1(ua11) + abs1(ua12)) == zero) {
                Clartg(vb12, vb11, csq, snq, r);
            } else if ((abs1(vb11) + abs1(vb12)) == zero) {
                Clartg(ua12, ua11, csq, snq, r);
            } else if (aua11 / (abs1(ua11) + abs1(ua12)) <= avb11 / (abs1(vb11) + abs1(vb12))) {
                Clartg(ua12, ua11, csq, snq, r);
            } else {
                Clartg(vb12, vb11, csq, snq, r);
            }
            //
            csu = snr;
            snu = conj(d1) * csr;
            csv = snl;
            snv = conj(d1) * csl;
            //
        }
        //
    }
    //
    //     End of Clags2
    //
}
