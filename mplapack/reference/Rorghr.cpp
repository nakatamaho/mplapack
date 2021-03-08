/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rorghr.cpp,v 1.8 2010/08/07 04:48:33 nakatamaho Exp $ 
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


void Rorghr(INTEGER n, INTEGER ilo, INTEGER ihi, REAL * A, INTEGER lda, REAL * tau, REAL * work, INTEGER lwork, INTEGER * info)
{

    INTEGER i, j, nb, nh, iinfo;
    INTEGER lwkopt;
    INTEGER lquery;
    REAL One = 1.0, Zero = 0.0;

    *info = 0;
    nh = ihi - ilo;
    lquery = lwork == -1;
    if (n < 0) {
	*info = -1;
    } else if (ilo < 1 || ilo > max((INTEGER) 1, n)) {
	*info = -2;
    } else if (ihi < min(ilo, n) || ihi > n) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -5;
    } else if (lwork < max((INTEGER) 1, nh) && !lquery) {
	*info = -8;
    }

    if (*info == 0) {
	nb = iMlaenv(1, "Rorgqr", " ", nh, nh, nh, -1);
	lwkopt = max((INTEGER) 1, nh) * nb;
	work[0] = lwkopt;
    }
    if (*info != 0) {
	Mxerbla("Rorghr", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	work[0] = One;
	return;
    }
//Shift the vectors which define the elementary reflectors one
//column to the right, and set the first ilo and the last n-ihi
//rows and columns to those of the unit matrix
    for (j = ihi; j >= ilo + 1; j--) {
	for (i = 0; i < j - 1; i++) {
	    A[i + j * lda] = Zero;

	}
	for (i = j + 1; i < ihi; i++) {
	    A[i + j * lda] = A[i + (j - 1) * lda];

	}
	for (i = ihi + 1; i < n; i++) {
	    A[i + j * lda] = Zero;
	}
    }

    for (j = 0; j < ilo; j++) {
	for (i = 0; i < n; i++) {
	    A[i + j * lda] = Zero;
	}
	A[j + j * lda] = One;

    }
    for (j = ihi + 1; j < n; j++) {
	for (i = 0; i < n; i++) {
	    A[i + j * lda] = Zero;
	}
	A[j + j * lda] = One;
    }

    if (nh > 0) {
//Generate Q(ilo+1:ihi,ilo+1:ihi)
	Rorgqr(nh, nh, nh, &A[ilo + 1 + (ilo + 1) * lda], lda, &tau[ilo], &work[0], lwork, &iinfo);
    }
    work[0] = lwkopt;
    return;
}
