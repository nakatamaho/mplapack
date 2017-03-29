/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgbcon.cpp,v 1.7 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rgbcon(const char *norm, INTEGER n, INTEGER kl, INTEGER ku, REAL * AB, INTEGER ldab, INTEGER * ipiv, REAL anorm, REAL * rcond, REAL * work, INTEGER * iwork, INTEGER * info)
{
    INTEGER j;
    REAL t;
    INTEGER kd, lm, jp, ix, kase;
    INTEGER kase1;
    REAL scale;
    INTEGER isave[3];
    INTEGER lnoti;
    REAL ainvnm;
    INTEGER onenrm;
    char normin;
    REAL smlnum;
    REAL Zero = 0.0, One = 1.0;

//Test the input parameters.
    *info = 0;
    onenrm = Mlsame(norm, "1") || Mlsame(norm, "O");
    if (!onenrm && !Mlsame(norm, "I")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (kl < 0) {
	*info = -3;
    } else if (ku < 0) {
	*info = -4;
    } else if (ldab < (kl << 1) + ku + 1) {
	*info = -6;
    } else if (anorm < Zero) {
	*info = -8;
    }
    if (*info != 0) {
	Mxerbla("Rgbcon", -(*info));
	return;
    }
//Quick return if possible
    *rcond = Zero;
    if (n == 0) {
	*rcond = One;
	return;
    } else if (anorm == Zero) {
	return;
    }
    smlnum = Rlamch("Safe minimum");

//Estimate the norm of inv(A).
    ainvnm = Zero;
    normin = 'N';
    if (onenrm) {
	kase1 = 1;
    } else {
	kase1 = 2;
    }
    kd = kl + ku + 1;
    lnoti = kl > 0;
    kase = 0;
    while (1) {
	Rlacn2(n, &work[n + 1], &work[0], &iwork[1], &ainvnm, &kase, isave);
	if (kase != 0) {
	    if (kase == kase1) {
//Multiply by inv(L).
		if (lnoti) {
		    for (j = 0; j < n - 1; j++) {
			lm = min(kl, n - j);
			jp = ipiv[j];
			t = work[jp];
			if (jp != j) {
			    work[jp] = work[j];
			    work[j] = t;
			}
			Raxpy(lm, -t, &AB[kd + 1 + j * ldab], 1, &work[j + 1], 1);
		    }
		}
//Multiply by inv(U).
		Rlatbs("Upper", "No transpose", "Non-unit", &normin, n, kl + ku, &AB[0], ldab, &work[0], &scale, &work[(n << 1) + 1], info);
	    } else {
//Multiply by inv(U').
		Rlatbs("Upper", "Transpose", "Non-unit", &normin, n, kl + ku, &AB[0], ldab, &work[0], &scale, &work[(n << 1) + 1], info);
//Multiply by inv(L').
		if (lnoti) {
		    for (j = n - 1; j >= 1; j--) {
			lm = min(kl, n - j);
			work[j] -= Rdot(lm, &AB[kd + 1 + j * ldab], 1, &work[j + 1], 1);
			jp = ipiv[j];
			if (jp != j) {
			    t = work[jp];
			    work[jp] = work[j];
			    work[j] = t;
			}

		    }
		}
	    }

//Divide X by 1/SCALE if doing so will not cause overflow.
	    normin = 'Y';
	    if (scale != One) {
		ix = iRamax(n, &work[0], 1);
		if (scale < abs(work[ix]) * smlnum || scale == Zero) {
		    return;
		}
		Rrscl(n, scale, &work[0], 1);
	    }
	}
    }

//Compute the estimate of the reciprocal condition number.
    if (ainvnm != Zero) {
	*rcond = One / ainvnm / anorm;
    }
    return;
}
