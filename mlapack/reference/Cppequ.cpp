/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cppequ.cpp,v 1.8 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Cppequ(const char *uplo, INTEGER n, COMPLEX * ap, REAL * s, REAL * scond, REAL * amax, INTEGER * info)
{
    INTEGER i, jj;
    REAL smin;
    INTEGER upper;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1, mtemp2;

//Test the input parameters.
    *info = 0;
    upper = Mlsame(uplo, "U");
    if (!upper && !Mlsame(uplo, "L")) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	Mxerbla("Cppequ", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	*scond = One;
	*amax = Zero;
	return;
    }
//Initialize SMIN and AMAX.
    s[0] = ap[1].real();
    smin = s[0];
    *amax = s[0];
    if (upper) {
//UPLO = 'U':  Upper triangle of A is stored.
//Find the minimum and maximum diagonal elements.
	jj = 0;
	for (i = 1; i < n; i++) {
	    jj = jj + i;
	    s[i] = ap[jj].real();
	    mtemp1 = smin, mtemp2 = s[i];
	    smin = min(mtemp1, mtemp2);
	    mtemp1 = *amax, mtemp2 = s[i];
	    *amax = max(mtemp1, mtemp2);
	}
    } else {
//UPLO = 'L':  Lower triangle of A is stored.
//Find the minimum and maximum diagonal elements.
	jj = 0;
	for (i = 1; i < n; i++) {
	    jj = jj + n - i + 2;
	    s[i] = ap[jj].real();
	    mtemp1 = smin, mtemp2 = s[i];
	    smin = min(mtemp1, mtemp2);
	    mtemp1 = *amax, mtemp2 = s[i];
	    *amax = max(mtemp1, mtemp2);
	}
    }
    if (smin <= Zero) {
//Find the first non-positive diagonal element and return.
	for (i = 0; i < n; i++) {
	    if (s[i] <= Zero) {
		*info = i;
		return;
	    }
	}
    } else {
//Set the scale factors to the reciprocals
//of the diagonal elements.
	for (i = 0; i < n; i++) {
	    s[i] = One / sqrt(s[i]);
	}
//Compute SCOND = min(S(I)) / max(S(I))
	*scond = sqrt(smin) / sqrt(*amax);
    }
    return;
}
