/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Cgghrd.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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

#define MTRUE 1
#define MFALSE 0

void Cgghrd(const char *compq, const char *compz, INTEGER n, INTEGER ilo, INTEGER ihi, COMPLEX * A, INTEGER lda, COMPLEX * B,
	    INTEGER ldb, COMPLEX * q, INTEGER ldq, COMPLEX * z, INTEGER ldz, INTEGER * info)
{
    REAL c;
    COMPLEX s;
    INTEGER ilq = 0, ilz = 0;
    INTEGER jcol, jrow;
    COMPLEX ctemp;
    INTEGER icompq, icompz;
    REAL One = 1.0, Zero = 0.0;

//Decode COMPQ
    if (Mlsame(compq, "N")) {
	ilq = MFALSE;
	icompq = 1;
    } else if (Mlsame(compq, "V")) {
	ilq = MTRUE;
	icompq = 2;
    } else if (Mlsame(compq, "I")) {
	ilq = MTRUE;
	icompq = 3;
    } else {
	icompq = 0;
    }
//Decode COMPZ
    if (Mlsame(compz, "N")) {
	ilz = MFALSE;
	icompz = 1;
    } else if (Mlsame(compz, "V")) {
	ilz = MTRUE;
	icompz = 2;
    } else if (Mlsame(compz, "I")) {
	ilz = MTRUE;
	icompz = 3;
    } else {
	icompz = 0;
    }
//Test the input parameters.
    *info = 0;
    if (icompq <= 0) {
	*info = -1;
    } else if (icompz <= 0) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (ilo < 1) {
	*info = -4;
    } else if (ihi > n || ihi < ilo - 1) {
	*info = -5;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -7;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -9;
    } else if ((ilq && ldq < n) || ldq < 1) {
	*info = -11;
    } else if ((ilz && ldz < n) || ldz < 1) {
	*info = -13;
    }
    if (*info != 0) {
	Mxerbla("Cgghrd", -(*info));
	return;
    }
//Initialize Q and Z if desired.
    if (icompq == 3) {
	Claset("Full", n, n, Zero, One, &q[0], ldq);
    }
    if (icompz == 3) {
	Claset("Full", n, n, Zero, One, &z[0], ldz);
    }
//Quick return if possible
    if (n <= 1) {
	return;
    }
//Zero out lower triangle of B
    for (jcol = 0; jcol <= n - 1; jcol++) {
	for (jrow = jcol + 1; jrow <= n; jrow++) {
	    B[jrow + jcol * ldb] = Zero;
	}
    }
//Reduce A and B
    for (jcol = ilo; jcol <= ihi - 2; jcol++) {
	for (jrow = ihi; jrow >= jcol + 2; jrow--) {
//Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL)
	    Clartg(A[jrow - 1 + jcol * lda], A[jrow + jcol * lda], &c, &s, &A[jrow - 1 + jcol * lda]);
	    A[jrow + jcol * lda] = Zero;
	    Crot(n - jcol, &A[jrow - 1 + (jcol + 1) * lda], lda, &A[jrow + (jcol + 1) * lda], lda, c, s);
	    Crot(n + 2 - jrow, &B[jrow - 1 + (jrow - 1) * ldb], ldb, &B[jrow + (jrow - 1) * ldb], ldb, c, s);
	    if (ilq) {
		Crot(n, &q[(jrow - 1) * ldq + 1], 1, &q[jrow * ldq + 1], 1, c, conj(s));
	    }
//Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1)
	    Clartg(B[jrow + jrow * ldb], B[jrow + (jrow - 1) * ldb], &c, &s, &B[jrow + jrow * ldb]);
	    B[jrow + (jrow - 1) * ldb] = Zero;
	    Crot(ihi, &A[jrow * lda], 1, &A[(jrow - 1) * lda], 1, c, s);
	    Crot(jrow - 1, &B[jrow * ldb + 1], 1, &B[(jrow - 1) * ldb + 1], 1, c, s);
	    if (ilz) {
		Crot(n, &z[jrow * ldz + 1], 1, &z[(jrow - 1) * ldz + 1], 1, c, s);
	    }
	}
    }
    return;
}
