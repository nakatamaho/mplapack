/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clags2.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

$COPYRIGHT$

Additional copyrights may follow

$HEADER$

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer. 
  
- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.
  
- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.
  
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
*/

#include <mblas.h>
#include <mlapack.h>

void Clags2(LOGICAL * upper, REAL a1, COMPLEX a2, REAL a3, REAL b1, COMPLEX b2, REAL b3, REAL * csu, COMPLEX * snu, REAL * csv, COMPLEX * snv, REAL * csq, COMPLEX * snq)
{
    REAL a;
    COMPLEX b, c;
    REAL d;
    COMPLEX r, d1;
    REAL s1, s2, fb, fc;
    COMPLEX ua11, ua12, ua21, ua22, vb11, vb12, vb21, vb22;
    REAL csl, csr, snl, snr, aua11, aua12, aua21, aua22, avb12, avb11, avb21, avb22, ua11r, ua22r, vb11r, vb22r;
    REAL Zero = 0.0, One = 1.0;

    if (*upper) {
//Input matrices A and B are upper triangular matrices
//Form matrix C = A*adj(B) = ( a b )
//                           ( 0 d )
	a = a1 * b3;
	d = a3 * b1;
	b = a2 * b1 - a1 * b2;
	fb = abs(b);
//Transform complex 2-by-2 matrix C to real matrix by unitary
//diagonal matrix diag(1,D1).
	if (fb != Zero) {
	    d1 = One / fb;
	}
//The SVD of real 2 by 2 triangular C
// ( CSL -SNL )*( A B )*(  CSR  SNR ) = ( R 0 )
// ( SNL  CSL ) ( 0 D ) ( -SNR  CSR )   ( 0 T )
	Rlasv2(a, fb, d, &s1, &s2, &snr, &csr, &snl, &csl);
	if (abs(csl) >= abs(snl) || abs(csr) >= abs(snr)) {
//Compute the (1,1) and (1,2) elements of U'*A and V'*B,
//and (1,2) element of |U|'*|A| and |V|'*|B|.
	    ua11r = csl * a1;
	    ua12 = csl * a2 + d1 * snl * a3;
	    vb11r = csr * b1;
	    vb12 = csr * b2 + d1 * snr * b3;
	    aua12 = abs(csl) * Cabs1(a2) + abs(snl) * abs(a3);
	    avb12 = abs(csr) * Cabs1(b2) + abs(snr) * abs(b3);
//zero (1,2) elements of U'*A and V'*B
	    if (abs(ua11r) + Cabs1(ua12) == Zero) {
		Clartg(-vb11, conj(vb12), csq, snq, &r);
	    } else if (abs(vb11r) + Cabs1(vb12) == Zero) {
		Clartg((COMPLEX) - ua11r, conj(ua12), csq, snq, &r);
	    } else if (aua12 / (abs(ua11r) + Cabs1(ua12)) <= avb12 / (abs(vb11r) + Cabs1(vb12))) {
		Clartg((COMPLEX) - ua11r, conj(ua12), csq, snq, &r);
	    } else {
		Clartg((COMPLEX) - vb11r, conj(vb12), csq, snq, &r);
	    }
	    *csu = csl;
	    *snu = -d1 * snl;
	    *csv = csr;
	    *snv = -d1 * snr;
	} else {
//Compute the (2,1) and (2,2) elements of U'*A and V'*B,
//and (2,2) element of |U|'*|A| and |V|'*|B|.
	    ua21 = -conj(d1) * snl * a1;
	    ua22 = -conj(d1) * snl * a2 + csl * a3;
	    vb21 = -conj(d1) * snr * b1;
	    vb22 = -conj(d1) * snr * b2 + csr * b3;
	    aua22 = abs(snl) * Cabs1(a2) + abs(csl) * abs(a3);
	    avb22 = abs(snr) * Cabs1(b2) + abs(csr) * abs(b3);
//zero (2,2) elements of U'*A and V'*B, and then swap.
	    if (Cabs1(ua21) + Cabs1(ua22) == Zero) {
		Clartg(-conj(vb21), conj(vb22), csq, snq, &r);
	    } else if (Cabs1(vb21) + abs(vb22) == Zero) {
		Clartg(-conj(ua21), conj(ua22), csq, snq, &r);
	    } else if (aua22 / (Cabs1(ua21) + Cabs1(ua22))
		       <= avb22 / Cabs1(vb21) + Cabs1(vb22)) {
		Clartg(-conj(ua21), conj(ua22), csq, snq, &r);
	    } else {
		Clartg(-conj(vb21), conj(vb22), csq, snq, &r);
	    }
	    *csu = snl;
	    *snu = d1 * csl;
	    *csv = snr;
	    *snv = d1 * csr;
	}
    } else {
//Input matrices A and B are lower triangular matrices
//Form matrix C = A*adj(B) = ( a 0 )
//                           ( c d )
	a = a1 * b3;
	d = a3 * b1;
	c = a2 * b3 - a3 * b2;
	fc = abs(c);
//Transform complex 2-by-2 matrix C to real matrix by unitary
//diagonal matrix diag(d1,1).
	d1 = One;
	if (fc != Zero) {
	    d1 = c / fc;
	}
//The SVD of real 2 by 2 triangular C
// ( CSL -SNL )*( A 0 )*(  CSR  SNR ) = ( R 0 )
// ( SNL  CSL ) ( C D ) ( -SNR  CSR )   ( 0 T )
	Rlasv2(a, fc, d, &s1, &s2, &snr, &csr, &snl, &csl);
	if (abs(csr) >= abs(snr) || abs(csl) >= abs(snl)) {
//Compute the (2,1) and (2,2) elements of U'*A and V'*B,
//and (2,1) element of |U|'*|A| and |V|'*|B|.
	    ua21 = -d1 * snr * a1 + csr * a2;
	    ua22r = csr * a3;
	    vb21 = -d1 * snl * b1 + csl * b2;
	    vb22r = csl * b3;
	    aua21 = abs(snr) * abs(a1) + abs(csr) * Cabs1(a2);
	    avb21 = abs(snl) * abs(b1) + abs(csl) * Cabs1(b2);
//zero (2,1) elements of U'*A and V'*B.
	    if (Cabs1(ua21) + abs(ua22r) == Zero) {
		Clartg(vb22, vb21, csq, snq, &r);
	    } else if (Cabs1(vb21) + abs(vb22r) == Zero) {
		Clartg(ua22r, ua21, csq, snq, &r);
	    } else if (aua21 / (Cabs1(ua21) + abs(ua22r)) <= avb21 / (Cabs1(vb21) + abs(vb22r))) {
		Clartg(ua22r, ua21, csq, snq, &r);
	    } else {
		Clartg(vb22r, vb21, csq, snq, &r);
	    }
	    *csu = csr;
	    *snu = -conj(d1) * snr;
	    *csv = csl;
	    *snv = -conj(d1) * snl;
	} else {
//Compute the (1,1) and (1,2) elements of U'*A and V'*B,
//and (1,1) element of |U|'*|A| and |V|'*|B|.
	    ua11 = csr * a1 + conj(d1) * snr * a2;
	    ua12 = conj(d1) * snr * a3;
	    vb11 = csl * b1 + conj(d1) * snl * b2;
	    vb12 = conj(d1) * snl * b3;
	    aua11 = abs(csr) * abs(a1) + abs(snr) * Cabs1(a2);
	    avb11 = abs(csl) * abs(b1) + abs(snl) * Cabs1(b2);
//zero (1,1) elements of U'*A and V'*B, and then swap.
	    if ((Cabs1(ua11) + Cabs1(ua12)) == Zero) {
		Clartg(vb12, vb11, csq, snq, &r);
	    } else if ((Cabs1(vb11) + Cabs1(vb12)) == Zero) {
		Clartg(ua12, ua11, csq, snq, &r);
	    } else if (aua11 / (Cabs1(ua11) + Cabs1(ua12)) <= avb11 / (Cabs1(vb11) + Cabs1(vb12))) {
		Clartg(ua12, ua11, csq, snq, &r);
	    } else {
		Clartg(vb12, vb11, csq, snq, &r);
	    }
	    *csu = snr;
	    *snu = conj(d1) * csr;
	    *csv = snl;
	    *snv = conj(d1) * csl;
	}
    }
    return;
}
