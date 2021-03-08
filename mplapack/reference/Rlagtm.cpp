/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlagtm.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rlagtm(const char *trans, INTEGER n, INTEGER nrhs, REAL alpha, REAL * dl, REAL * d, REAL * du, REAL * x, INTEGER ldx, REAL * beta, REAL * B, INTEGER ldb)
{
    INTEGER i, j;
    REAL One = 1.0, Zero = 0.0;

    if (n == 0) {
	return;
    }
//Multiply B by BETA if BETA!=1

    if (*beta == Zero) {
	for (j = 0; j < nrhs; j++) {
	    for (i = 0; i < n; i++) {
		B[i + j * ldb] = Zero;
	    }
	}
    } else if (*beta == -One) {
	for (j = 0; j < nrhs; j++) {
	    for (i = 0; i < n; i++) {
		B[i + j * ldb] = -B[i + j * ldb];
	    }
	}
    }
    if (alpha == One) {
	if (Mlsame(trans, "N")) {
//Compute B := B + A*X
	    for (j = 0; j < nrhs; j++) {
		if (n == 1) {
		    B[j * ldb + 1] += d[1] * x[j * ldx + 1];
		} else {
		    B[j * ldb + 1] = B[j * ldb + 1] + d[1] * x[j * ldx + 1] + du[1] * x[j * ldx + 2];
		    B[n + j * ldb] = B[n + j * ldb] + dl[n - 1] * x[n - 1 + j * ldx] + d[n] * x[n + j * ldx];
		    for (i = 1; i < n - 1; i++) {
			B[i + j * ldb] = B[i + j * ldb] + dl[i - 1] * x[i - 1 + j * ldx] + d[i] * x[i + j * ldx] + du[i] * x[i + 1 + j * ldx];
		    }
		}
	    }
	} else {
//Compute B := B + A'*X
	    for (j = 0; j < nrhs; j++) {
		if (n == 1) {
		    B[j * ldb + 1] += d[1] * x[j * ldx + 1];
		} else {
		    B[j * ldb + 1] = B[j * ldb + 1] + d[1] * x[j * ldx + 1] + dl[1] * x[j * ldx + 2];
		    B[n + j * ldb] = B[n + j * ldb] + du[n - 1] * x[n - 1 + j * ldx] + d[n] * x[n + j * ldx];
		    for (i = 1; i < n - 1; i++) {
			B[i + j * ldb] = B[i + j * ldb] + du[i - 1] * x[i - 1 + j * ldx] + d[i] * x[i + j * ldx] + dl[i] * x[i + 1 + j * ldx];

		    }
		}
	    }
	}
    } else if (alpha == -One) {
	if (Mlsame(trans, "N")) {
//Compute B := B - A*X
	    for (j = 0; j < nrhs; j++) {
		if (n == 1) {
		    B[j * ldb + 1] -= d[1] * x[j * ldx + 1];
		} else {
		    B[j * ldb + 1] = B[j * ldb + 1] - d[1] * x[j * ldx + 1] - du[1] * x[j * ldx + 2];
		    B[n + j * ldb] = B[n + j * ldb] - dl[n - 1] * x[n - 1 + j * ldx] - d[n] * x[n + j * ldx];
		    for (i = 1; i < n - 1; i++) {
			B[i + j * ldb] = B[i + j * ldb] - dl[i - 1] * x[i - 1 + j * ldx] - d[i] * x[i + j * ldx] - du[i] * x[i + 1 + j * ldx];

		    }
		}
	    }
	} else {
//Compute B := B - A'*X
	    for (j = 0; j < nrhs; j++) {
		if (n == 1) {
		    B[j * ldb + 1] -= d[1] * x[j * ldx + 1];
		} else {
		    B[j * ldb + 1] = B[j * ldb + 1] - d[1] * x[j * ldx + 1] - dl[1] * x[j * ldx + 2];
		    B[n + j * ldb] = B[n + j * ldb] - du[n - 1] * x[n - 1 + j * ldx] - d[n] * x[n + j * ldx];
		    for (i = 1; i < n - 1; i++) {
			B[i + j * ldb] = B[i + j * ldb] - du[i - 1] * x[i - 1 + j * ldx] - d[i] * x[i + j * ldx] - dl[i] * x[i + 1 + j * ldx];

		    }
		}
	    }
	}
    }
    return;
}
