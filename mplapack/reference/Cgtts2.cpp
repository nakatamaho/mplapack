/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgtts2.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

#include <mpblas.h>
#include <mplapack.h>

void Cgtts2(INTEGER itrans, INTEGER n, INTEGER nrhs, COMPLEX * dl, COMPLEX * d, COMPLEX * du, COMPLEX * du2, INTEGER * ipiv, COMPLEX * B, INTEGER ldb)
{
    INTEGER i, j;
    COMPLEX temp;

    if (n == 0 || nrhs == 0) {
	return;
    }
    if (itrans == 0) {
//Solve A*X = B using the LU factorization of A,
//overwriting each right hand side vector with its solution.
	if (nrhs <= 1) {
	    j = 0;
	  L10:
//Solve L*x = b.
	    for (i = 0; i < n - 1; i++) {
		if (ipiv[i] == i) {
		    B[i + 1 + j * ldb] = B[i + 1 + j * ldb] - dl[i] * B[i + j * ldb];
		} else {
		    temp = B[i + j * ldb];
		    B[i + j * ldb] = B[i + 1 + j * ldb];
		    B[i + 1 + j * ldb] = temp - dl[i] * B[i + j * ldb];
		}
	    }
//Solve U*x = b.
	    B[n + j * ldb] = B[n + j * ldb] / d[n];
	    if (n > 1) {
		B[n - 1 + j * ldb] = (B[n - 1 + j * ldb] - du[n - 1]
				      * B[n + j * ldb]) / d[n - 1];
	    }
	    for (i = n - 2; i >= 1; i--) {
		B[i + j * ldb] = (B[i + j * ldb] - du[i] * B[i + 1 + j * ldb] - du2[i] * B[i + 2 + j * ldb]
		    ) / d[i];
	    }
	    if (j < nrhs) {
		j++;
		goto L10;
	    }
	} else {
	    for (j = 0; j < nrhs; j++) {
//Solve L*x = b.
		for (i = 0; i < n - 1; i++) {
		    if (ipiv[i] == i) {
			B[i + 1 + j * ldb] = B[i + 1 + j * ldb] - dl[i] * B[i + j * ldb];
		    } else {
			temp = B[i + j * ldb];
			B[i + j * ldb] = B[i + 1 + j * ldb];
			B[i + 1 + j * ldb] = temp - dl[i] * B[i + j * ldb];
		    }
		}
//Solve U*x = b.
		B[n + j * ldb] = B[n + j * ldb] / d[n];
		if (n > 1) {
		    B[n - 1 + j * ldb] = (B[n - 1 + j * ldb] - du[n - 1] * B[n + j * ldb]) / d[n - 1];
		}
		for (i = n - 2; i >= 1; i--) {
		    B[i + j * ldb] = (B[i + j * ldb] - du[i] * B[i + 1 + j * ldb] - du2[i] * B[i + 2 + j * ldb]) / d[i];
		}
	    }
	}
    } else if (itrans == 1) {
//Solve A**T * X = B.
	if (nrhs <= 1) {
	    j = 0;
	  L70:
//Solve U**T * x = b.
	    B[j * ldb + 1] = B[j * ldb + 1] / d[1];
	    if (n > 1) {
		B[j * ldb + 2] = (B[j * ldb + 2] - du[1] * B[j * ldb + 1]) / d[2];
	    }
	    for (i = 3; i <= n; i++) {
		B[i + j * ldb] = (B[i + j * ldb] - du[i - 1] * B[i - 1 + j * ldb] - du2[i - 2] * B[i - 2 + j * ldb]) / d[i];
	    }
//Solve L**T * x = b.
	    for (i = n - 1; i >= 1; i--) {
		if (ipiv[i] == i) {
		    B[i + j * ldb] = B[i + j * ldb] - dl[i] * B[i + 1 + j * ldb];
		} else {
		    temp = B[i + 1 + j * ldb];
		    B[i + 1 + j * ldb] = B[i + j * ldb] - dl[i] * temp;
		    B[i + j * ldb] = temp;
		}
	    }
	    if (j < nrhs) {
		j++;
		goto L70;
	    }
	} else {
	    for (j = 0; j < nrhs; j++) {
//Solve U**T * x = b.
		B[j * ldb + 1] = B[j * ldb + 1] / d[1];
		if (n > 1) {
		    B[j * ldb + 2] = (B[j * ldb + 2] - du[1] * B[j * ldb + 1]) / d[2];
		}
		for (i = 3; i <= n; i++) {
		    B[i + j * ldb] = (B[i + j * ldb] - du[i - 1] * B[i - 1 + j * ldb] - du2[i - 2] * B[i - 2 + j * ldb]) / d[i];
		}
//Solve L**T * x = b.
		for (i = n - 1; i >= 1; i--) {
		    if (ipiv[i] == i) {
			B[i + j * ldb] = B[i + j * ldb] - dl[i] * B[i + 1 + j * ldb];
		    } else {
			temp = B[i + 1 + j * ldb];
			B[i + 1 + j * ldb] = B[i + j * ldb] - dl[i] * temp;
			B[i + j * ldb] = temp;
		    }

		}
	    }
	}
    } else {
//Solve A**H * X = B.
	if (nrhs <= 1) {
	    j = 0;
	  L130:
//Solve U**H * x = b.
	    B[j * ldb + 1] = B[j * ldb + 1] / conj(d[0]);
	    if (n > 1) {
		B[j * ldb + 2] = (B[j * ldb + 2] - conj(du[1]) * B[j * ldb + 1]) / conj(d[2]);
	    }
	    for (i = 3; i <= n; i++) {
		B[i + j * ldb] = (B[i + j * ldb] - conj(du[i - 1]) * B[i - 1 + j * ldb] - conj(du2[i - 2]) * B[i - 2 + j * ldb]) / conj(d[i]);

	    }
//Solve L**H * x = b.
	    for (i = n - 1; i >= 1; i--) {
		if (ipiv[i] == i) {
		    B[i + j * ldb] = B[i + j * ldb] - conj(dl[i]) * B[i + 1 + j * ldb];
		} else {
		    temp = B[i + 1 + j * ldb];
		    B[i + 1 + j * ldb] = B[i + j * ldb] - conj(dl[i]) * temp;
		    B[i + j * ldb] = temp;
		}

	    }
	    if (j < nrhs) {
		j++;
		goto L130;
	    }
	} else {
	    for (j = 0; j < nrhs; j++) {
//Solve U**H * x = b.
		B[j * ldb + 1] = B[j * ldb + 1] / conj(d[0]);
		if (n > 1) {
		    B[j * ldb + 2] = (B[j * ldb + 2] - conj(du[1]) * B[j * ldb + 1]) / conj(d[2]);
		}
		for (i = 3; i <= n; i++) {
		    B[i + j * ldb] = (B[i + j * ldb] - conj(du[i - 1]) * B[i - 1 + j * ldb] - conj(du2[i - 2]) * B[i - 2 + j * ldb]) / conj(d[i]);
		}
//Solve L**H * x = b.
		for (i = n - 1; i >= 1; i--) {
		    if (ipiv[i] == i) {
			B[i + j * ldb] = B[i + j * ldb] - conj(dl[i]) * B[i + 1 + j * ldb];
		    } else {
			temp = B[i + 1 + j * ldb];
			B[i + 1 + j * ldb] = B[i + j * ldb] - conj(dl[i]) * temp;
			B[i + j * ldb] = temp;
		    }
		}
	    }
	}
    }
    return;
}
