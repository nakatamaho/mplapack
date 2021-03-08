/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rppcon.cpp,v 1.7 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rppcon(const char *uplo, INTEGER n, REAL * ap, REAL anorm, REAL * rcond, REAL * work, INTEGER * iwork, INTEGER * info)
{
    INTEGER ix, kase;
    REAL scale;
    INTEGER isave[3];
    INTEGER upper;
    REAL scalel;
    REAL scaleu;
    REAL ainvnm;
    char normin;
    REAL smlnum;
    REAL Zero = 0.0, One = 1.0;

    *info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    } else if (anorm < Zero) {
	*info = -4;
    }
    if (*info != 0) {
	Mxerbla("Rppcon", -(*info));
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

//Estimate the 1-norm of the inverse.
    kase = 0;
    normin = 'N';
    Rlacn2(n, &work[n + 1], &work[0], &iwork[1], &ainvnm, &kase, isave);
    while (1) {
	if (kase != 0) {
	    if (upper) {

//Multiply by inv(U').
		Rlatps("Upper", "Transpose", "Non-unit", &normin, n, &ap[0], &work[0], &scalel, &work[(n << 1) + 1], info);
		normin = 'Y';

//Multiply by inv(U).
		Rlatps("Upper", "No transpose", "Non-unit", &normin, n, &ap[0], &work[0], &scaleu, &work[(n << 1) + 1], info);
	    } else {

//Multiply by inv(L).
		Rlatps("Lower", "No transpose", "Non-unit", &normin, n, &ap[0], &work[0], &scalel, &work[(n << 1) + 1], info);
		normin = 'Y';

//Multiply by inv(L').
		Rlatps("Lower", "Transpose", "Non-unit", &normin, n, &ap[0], &work[0], &scaleu, &work[(n << 1) + 1], info);
	    }

//Multiply by 1/SCALE if doing so will not cause overflow.
	    scale = scalel * scaleu;
	    if (scale != One) {
		ix = iRamax(n, &work[0], 1);
		if (scale < abs(work[ix]) * smlnum || scale == Zero) {
		    return;
		}
		Rrscl(n, scale, &work[0], 1);
	    }
	}
	break;
    }
// Compute the estimate of the reciprocal condition number.
    if (ainvnm != Zero) {
	*rcond = One / ainvnm / anorm;
    }

    return;
}
