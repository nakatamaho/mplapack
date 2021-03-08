/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Ctptri.cpp,v 1.8 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Ctptri(const char *uplo, const char *diag, INTEGER n, COMPLEX * ap, INTEGER * info)
{
    INTEGER j, jc, jj;
    COMPLEX ajj;
    INTEGER upper;
    INTEGER jclast = 0;
    INTEGER nounit;
    REAL Zero = 0.0, One = 1.0;

//Test the input parameters.
    *info = 0;
    upper = Mlsame(uplo, "U");
    nounit = Mlsame(diag, "N");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (!nounit && !Mlsame(diag, "U")) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    }
    if (*info != 0) {
	Mxerbla("Ctptri", -(*info));
	return;
    }
//Check for singularity if non-unit.
    if (nounit) {
	if (upper) {
	    jj = 0;
	    for (*info = 1; *info <= n; ++(*info)) {
		jj = jj + *info;
		if (ap[jj] == Zero) {
		    return;
		}
	    }
	} else {
	    jj = 0;
	    for (*info = 1; *info <= n; ++(*info)) {
		if (ap[jj] == Zero) {
		    return;
		}
		jj = jj + n - *info + 1;
	    }
	}
	*info = 0;
    }
    if (upper) {
//Compute inverse of upper triangular matrix.
	jc = 1;
	for (j = 0; j < n; j++) {
	    if (nounit) {
		ap[jc + j - 1] = One / ap[jc + j - 1];
		ajj = -ap[jc + j - 1];
	    } else {
		ajj = -One;
	    }
//Compute elements 1:j-1 of j-th column.
	    Ctpmv("Upper", "No transpose", diag, j - 1, &ap[1], &ap[jc], 1);
	    Cscal(j - 1, ajj, &ap[jc], 1);
	    jc = jc + j;
	}
    } else {
//Compute inverse of lower triangular matrix.
	jc = n * (n + 1) / 2;
	for (j = n; j >= 1; j--) {
	    if (nounit) {
		ap[jc] = One / ap[jc];
		ajj = -ap[jc];
	    } else {
		ajj = -One;
	    }
	    if (j < n) {
//Compute elements j+1:n of j-th column.
		Ctpmv("Lower", "No transpose", diag, n - j, &ap[jclast], &ap[jc + 1], 12);
		Cscal(n - j, ajj, &ap[jc + 1], 1);
	    }
	    jclast = jc;
	    jc = jc - n + j - 2;
	}
    }
    return;
}
