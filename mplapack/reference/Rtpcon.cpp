/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rtpcon.cpp,v 1.8 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rtpcon(const char *norm, const char *uplo, const char *diag, INTEGER n, REAL * ap, REAL * rcond, REAL * work, INTEGER * iwork, INTEGER * info)
{
    INTEGER ix, kase, kase1;
    REAL scale;
    INTEGER isave[3];
    REAL anorm;
    INTEGER upper;
    REAL xnorm;
    REAL ainvnm;
    INTEGER onenrm;
    char normin;
    REAL smlnum;
    INTEGER nounit;
    REAL One = 1.0, Zero = 0.0;

//Test the input parameters.
    *info = 0;
    upper = Mlsame(uplo, "U");
    onenrm = Mlsame(norm, "1") || Mlsame(norm, "O");
    nounit = Mlsame(diag, "N");

    if (!onenrm && !Mlsame(norm, "I")) {
	*info = -1;
    } else if (!upper && !Mlsame(uplo, "L")) {
	*info = -2;
    } else if (!nounit && !Mlsame(diag, "U")) {
	*info = -3;
    } else if (n < 0) {
	*info = -4;
    }
    if (*info != 0) {
	Mxerbla("Rtpcon", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	*rcond = One;
	return;
    }

    *rcond = Zero;
    smlnum = Rlamch("Safe minimum") * max((INTEGER) 1, n);

//Compute the norm of the triangular matrix A.
    anorm = Rlantp(norm, uplo, diag, n, &ap[1], &work[0]);

//Continue only if ANORM > Zero

    if (anorm > Zero) {
//Estimate the norm of the inverse of A.
	ainvnm = Zero;
	normin = 'N';
	if (onenrm) {
	    kase1 = 1;
	} else {
	    kase1 = 2;
	}
	kase = 0;
	while (1) {
	    Rlacn2(n, &work[n + 1], &work[0], &iwork[0], &ainvnm, &kase, isave);
	    if (kase != 0) {
		if (kase == kase1) {
//Multiply by inv(A).
		    Rlatps(uplo, "No transpose", diag, &normin, n, &ap[0], &work[0], &scale, &work[(n << 1) + 1], info);
		} else {
//Multiply by inv(A').
		    Rlatps(uplo, "Transpose", diag, &normin, n, &ap[0], &work[0], &scale, &work[(n << 1) + 1], info);
		}
		normin = 'Y';

//Multiply by 1/SCALE if doing so will not cause overflow.
		if (scale != One) {
		    ix = iRamax(n, &work[0], 1);
		    xnorm = abs(work[ix]);
		    if (scale < xnorm * smlnum || scale == Zero) {
			return;
		    }
		    Rrscl(n, scale, &work[0], 1);
		}
	    }
	}
//Compute the estimate of the reciprocal condition number.
	if (ainvnm != Zero) {
	    *rcond = One / anorm / ainvnm;
	}
    }
    return;
}
