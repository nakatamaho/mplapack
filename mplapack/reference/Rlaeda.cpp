/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaeda.cpp,v 1.8 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void
Rlaeda(INTEGER n, INTEGER tlvls, INTEGER curlvl,
       INTEGER curpbm, INTEGER * prmptr, INTEGER * perm, INTEGER * givptr, INTEGER * givcol, REAL * givnum, REAL * q, INTEGER * qptr, REAL * z, REAL * ztemp, INTEGER * info)
{
    INTEGER i, k, mid, ptr;
    INTEGER curr, bsiz1, bsiz2, psiz1, psiz2, zptr1;
    REAL One = 1.0, Zero = 0.0;
    double Half = 0.5;

//Test the input parameters.
    *info = 0;
    if (n < 0) {
	*info = -1;
    }
    if (*info != 0) {
	Mxerbla("Rlaeda", -(*info));
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//Determine location of first number in second half.
    mid = n / 2 + 1;
//Gather last/first rows of appropriate eigenblocks into center of Z
    ptr = 1;
//Determine location of lowest level subproblem in the full storage
//scheme
    curr = ptr + curpbm * (2 ^ (curlvl)) + (2 ^ (curlvl - 1)) - 1;
//Determine size of these matrices.  We add HALF to the value of
//the SQRT in case the machine underestimates one of these square
//roots.

    bsiz1 = (INTEGER) (sqrt((qptr[curr + 1] - qptr[curr])) + Half);
    bsiz2 = (INTEGER) (sqrt((qptr[curr + 2] - qptr[curr + 1])) + Half);
    for (k = 0; k < mid - bsiz1 - 1; k++) {
	z[k] = Zero;
    }
    Rcopy(bsiz1, &q[qptr[curr] + bsiz1 - 1], bsiz1, &z[mid - bsiz1], 1);
    Rcopy(bsiz2, &q[qptr[curr + 1]], bsiz2, &z[mid], 1);
    for (k = mid + bsiz2; k <= n; k++) {
	z[k] = Zero;
    }

//Loop thru remaining levels 1 -> CURLVL applying the Givens
//rotations and permutation and then multiplying the center matrices
//against the current Z.
    ptr = (2 ^ tlvls) + 1;
    for (k = 0; k < curlvl - 1; k++) {
	curr = ptr + curpbm * (2 ^ (curlvl - k)) + (2 ^ (curlvl - k - 1)) - 1;
	psiz1 = prmptr[curr + 1] - prmptr[curr];
	psiz2 = prmptr[curr + 2] - prmptr[curr + 1];
	zptr1 = mid - psiz1;
//Apply Givens at CURR and CURR+1
	for (i = givptr[curr]; i <= givptr[curr + 1] - 1; i++) {
	    Rrot(1, &z[zptr1 + givcol[(i << 1) + 1] - 1], 1, &z[zptr1 + givcol[(i << 1) + 2] - 1], 1, givnum[(i << 1) + 1], givnum[(i << 1) + 2]);

	}
	for (i = givptr[curr + 1]; i <= givptr[curr + 2] - 1; i++) {
	    Rrot(1, &z[mid - 1 + givcol[(i << 1) + 1]], 1, &z[mid - 1 + givcol[(i << 1) + 2]], 1, givnum[(i << 1) + 1], givnum[(i << 1) + 2]);

	}
	psiz1 = prmptr[curr + 1] - prmptr[curr];
	psiz2 = prmptr[curr + 2] - prmptr[curr + 1];
	for (i = 0; i < psiz1 - 1; i++) {
	    ztemp[i + 1] = z[zptr1 + perm[prmptr[curr] + i] - 1];

	}
	for (i = 0; i < psiz2 - 1; i++) {
	    ztemp[psiz1 + i + 1] = z[mid + perm[prmptr[curr + 1] + i] - 1];
	}

//Multiply Blocks at CURR and CURR+1
//Determine size of these matrices.  We add HALF to the value of
//the SQRT in case the machine underestimates one of these
//square roots.
	bsiz1 = (INTEGER) (sqrt((qptr[curr + 1] - qptr[curr])) + Half);
	bsiz2 = (INTEGER) (sqrt((qptr[curr + 2] - qptr[curr + 1])) + Half);
	if (bsiz1 > 0) {
	    Rgemv("T", bsiz1, bsiz1, One, &q[qptr[curr]], bsiz1, &ztemp[1], 1, Zero, &z[zptr1], 1);
	}
	Rcopy(psiz1 - bsiz1, &ztemp[bsiz1 + 1], 1, &z[zptr1 + bsiz1], 1);
	if (bsiz2 > 0) {
	    Rgemv("T", bsiz2, bsiz2, One, &q[qptr[curr + 1]], bsiz2, &ztemp[psiz1 + 1], 1, Zero, &z[mid], 1);
	}
	Rcopy(psiz2 - bsiz2, &ztemp[psiz1 + bsiz2 + 1], 1, &z[mid + bsiz2], 1);
	ptr += 2 ^ (tlvls - k);
    }
    return;
}
