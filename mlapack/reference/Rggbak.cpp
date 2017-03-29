/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rggbak.cpp,v 1.9 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rggbak(const char *job, const char *side, INTEGER n, INTEGER ilo, INTEGER ihi, REAL * lscale, REAL * rscale, INTEGER m, REAL * v, INTEGER ldv, INTEGER * info)
{
    INTEGER i, k;
    INTEGER leftv;
    INTEGER rightv;

    rightv = Mlsame(side, "R");
    leftv = Mlsame(side, "L");

    *info = 0;
    if (!Mlsame(job, "N") && !Mlsame(job, "P") && !Mlsame(job, "S")
	&& !Mlsame(job, "B")) {
	*info = -1;
    } else if (!rightv && !leftv) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (ilo < 1) {
	*info = -4;
    } else if (n == 0 && ihi == 0 && ilo != 1) {
	*info = -4;
    } else if (n > 0 && (ihi < ilo || ihi > max((INTEGER) 1, n))) {
	*info = -5;
    } else if (n == 0 && ilo == 1 && ihi != 0) {
	*info = -5;
    } else if (m < 0) {
	*info = -8;
    } else if (ldv < max((INTEGER) 1, n)) {
	*info = -10;
    }
    if (*info != 0) {
	Mxerbla("Rggbak", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0)
	return;
    if (m == 0)
	return;
    if (Mlsame(job, "N"))
	return;
    if (ilo != ihi) {
//Backward balance
	if (Mlsame(job, "S") || Mlsame(job, "B")) {
//Backward transformation on right eigenvectors
	    if (rightv) {
		for (i = ilo; i <= ihi; i++) {
		    Rscal(m, rscale[i], &v[i + ldv], ldv);
		}
	    }
//Backward transformation on left eigenvectors
	    if (leftv) {
		for (i = ilo; i <= ihi; i++) {
		    Rscal(m, lscale[i], &v[i + ldv], ldv);
		}
	    }
	}
    }
//Backward permutation
    if (Mlsame(job, "P") || Mlsame(job, "B")) {
//Backward permutation on right eigenvectors
	if (rightv) {
	    if (ilo != 1) {
		for (i = ilo - 1; i >= 1; i--) {
		    k = (INTEGER) cast2double(rscale[i]);
		    if (k != i) {
			Rswap(m, &v[i + ldv], ldv, &v[k + ldv], ldv);
		    }
		    if (ihi == n)
			break;
		    for (i = ihi + 1; i <= n; i++) {
			k = (INTEGER) cast2double(rscale[i]);
			if (k != i) {
			    Rswap(m, &v[i + ldv], ldv, &v[k + ldv], ldv);
			}
		    }
		}
	    }
	}
//Backward permutation on left eigenvectors
	if (leftv) {
	    if (ilo != 1) {
		for (i = ilo - 1; i >= 1; i--) {
		    k = (INTEGER) cast2double(lscale[i]);
		    if (k != i) {
			Rswap(m, &v[i + ldv], ldv, &v[k + ldv], ldv);
		    }
		}
	    }
	    if (ihi == n)
		return;
	    for (i = ihi + 1; i <= n; i++) {
		k = (INTEGER) cast2double(lscale[i]);
		if (k != i) {
		    Rswap(m, &v[i + ldv], ldv, &v[k + ldv], ldv);
		}
	    }
	}
    }
    return;
}
