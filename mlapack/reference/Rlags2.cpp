/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlags2.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rlags2(INTEGER upper, REAL a1, REAL a2, REAL a3, REAL b1, REAL b2, REAL b3, REAL * csu, REAL * snu, REAL * csv, REAL * snv, REAL * csq, REAL * snq)
{
    REAL a, b, c, d, r, s1, s2;
    REAL ua11, ua12, ua21, ua22, vb11, vb12, vb21, vb22;
    REAL csl, csr, snl, snr;
    REAL aua11, aua12, aua21, aua22;
    REAL avb11, avb12, avb21, avb22;
    REAL ua11r, ua22r, vb11r, vb22r;
    REAL Zero = 0.0;

    if (upper) {
//Input matrices A and B are upper triangular matrices
//Form matrix C = A*adj(B) = ( a b )
//                           ( 0 d )
	a = a1 * b3;
	d = a3 * b1;
	b = a2 * b1 - a1 * b2;
//The SVD of real 2-by-2 triangular C
// ( CSL -SNL )*( A B )*(  CSR  SNR ) = ( R 0 )
// ( SNL  CSL ) ( 0 D ) ( -SNR  CSR )   ( 0 T )

	Rlasv2(a, b, d, &s1, &s2, &snr, &csr, &snl, &csl);
	if (abs(csl) >= abs(snl) || abs(csr) >= abs(snr)) {

//Compute the (1,1) and (1,2) elements of U'*A and V'*B,
//and (1,2) element of |U|'*|A| and |V|'*|B|.
	    ua11r = csl * a1;
	    ua12 = csl * a2 + snl * a3;

	    vb11r = csr * b1;
	    vb12 = csr * b2 + snr * b3;

	    aua12 = abs(csl) * abs(a2) + abs(snl) * abs(a3);
	    avb12 = abs(csr) * abs(b2) + abs(snr) * abs(b3);
//zero (1,2) elements of U'*A and V'*B 
	    if (abs(ua11r) + abs(ua12) != Zero) {
		if (aua12 / (abs(ua11r) + abs(ua12)) <= avb12 / (abs(vb11r) + abs(vb12))) {
		    Rlartg(-ua11r, ua12, csq, snq, &r);
		} else {
		    Rlartg(-vb11r, vb12, csq, snq, &r);
		}
	    } else {
		Rlartg(-vb11r, vb12, csq, snq, &r);
	    }
	    *csu = csl;
	    *snu = -snl;
	    *csv = csr;
	    *snv = -snr;
	} else {
//Compute the (2,1) and (2,2) elements of U'*A and V'*B,
//and (2,2) element of |U|'*|A| and |V|'*|B|.
	    ua21 = -snl * a1;
	    ua22 = -snl * a2 + csl * a3;

	    vb21 = -snr * b1;
	    vb22 = -snr * b2 + csr * b3;

	    aua22 = abs(snl) * abs(a2) + abs(csl) * abs(a3);
	    avb22 = abs(snr) * abs(b2) + abs(csr) * abs(b3);

//zero (2,2) elements of U'*A and V'*B, and then swap.
	    if (abs(ua21) + abs(ua22) != Zero) {
		if (aua22 / (abs(ua21) + abs(ua22)) <= avb22 / (abs(vb21) + abs(vb22))) {
		    Rlartg(-ua21, ua22, csq, snq, &r);
		} else {
		    Rlartg(-vb21, vb22, csq, snq, &r);
		}
	    } else {
		Rlartg(-vb21, vb22, csq, snq, &r);
	    }
	    *csu = snl;
	    *snu = csl;
	    *csv = snr;
	    *snv = csr;
	}
    } else {
//Input matrices A and B are lower triangular matrices
//Form matrix C = A*adj(B) = ( a 0 )
//                           ( c d )
	a = a1 * b3;
	d = a3 * b1;
	c = a2 * b3 - a3 * b2;
//The SVD of real 2-by-2 triangular C
// ( CSL -SNL )*( A 0 )*(  CSR  SNR ) = ( R 0 )
// ( SNL  CSL ) ( C D ) ( -SNR  CSR )   ( 0 T )
	Rlasv2(a, c, d, &s1, &s2, &snr, &csr, &snl, &csl);
	if (abs(csr) >= abs(snr) || abs(csl) >= abs(snl)) {
//Compute the (2,1) and (2,2) elements of U'*A and V'*B,
//and (2,1) element of |U|'*|A| and |V|'*|B|.
	    ua21 = -snr * a1 + csr * a2;
	    ua22r = csr * a3;
	    vb21 = -snl * b1 + csl * b2;
	    vb22r = csl * b3;

	    aua21 = abs(snr) * abs(a1) + abs(csr) * abs(a2);
	    avb21 = abs(snl) * abs(b1) + abs(csl) * abs(b2);
//zero (2,1) elements of U'*A and V'*B.
	    if (abs(ua21) + abs(ua22r) != Zero) {
		if (aua21 / (abs(ua21) + abs(ua22r)) <= avb21 / (abs(vb21) + abs(vb22r))) {
		    Rlartg(ua22r, ua21, csq, snq, &r);
		} else {
		    Rlartg(vb22r, vb21, csq, snq, &r);
		}
	    } else {
		Rlartg(vb22r, vb21, csq, snq, &r);
	    }
	    *csu = csr;
	    *snu = -snr;
	    *csv = csl;
	    *snv = -snl;

	} else {

//Compute the (1,1) and (1,2) elements of U'*A and V'*B,
//and (1,1) element of |U|'*|A| and |V|'*|B|.
	    ua11 = csr * a1 + snr * a2;
	    ua12 = snr * a3;

	    vb11 = csl * b1 + snl * b2;
	    vb12 = snl * b3;

	    aua11 = abs(csr) * abs(a1) + abs(snr) * abs(a2);
	    avb11 = abs(csl) * abs(b1) + abs(snl) * abs(b2);
//zero (1,1) elements of U'*A and V'*B, and then swap.
	    if (abs(ua11) + abs(ua12) != Zero) {
		if (aua11 / (abs(ua11) + abs(ua12)) <= avb11 / (abs(vb11) + abs(vb12))) {
		    Rlartg(ua12, ua11, csq, snq, &r);
		} else {
		    Rlartg(vb12, vb11, csq, snq, &r);
		}
	    } else {
		Rlartg(vb12, vb11, csq, snq, &r);
	    }

	    *csu = snr;
	    *snu = csr;
	    *csv = snl;
	    *snv = csl;
	}
    }
    return;
}
