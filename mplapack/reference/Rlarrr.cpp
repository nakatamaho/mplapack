/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlarrr.cpp,v 1.3 2010/08/07 04:48:33 nakatamaho Exp $ 
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

#define MTRUE 1
#define MFALSE 0

void Rlarrr(INTEGER n, REAL * d, REAL * e, INTEGER * info)
{
    INTEGER i;
    REAL eps, tmp, tmp2, rmin;
    REAL offdig, safmin;
    INTEGER yesrel;
    REAL smlnum, offdig2;
    REAL Zero = 0.0;

    *info = 1;
    safmin = Rlamch("Safe minimum");
    eps = Rlamch("Precision");
    smlnum = safmin / eps;
    rmin = sqrt(smlnum);
//Tests for relative accuracy 
//Test for scaled diagonal dominance 
//Scale the diagonal entries to one and check whether the sum of the 
//off-diagonals is less than one 
//The sdd relative error bounds have a 1/(1- 2*x) factor in them, 
//x = max(OFFDIG + OFFDIG2), so when x is close to 1/2, no relative 
//accuracy is promised.  In the notation of the code fragment below, 
//1/(1 - (OFFDIG + OFFDIG2)) is the condition number. 
//We don't think it is worth going into "sdd mode" unless the relative 
//condition number is reasonable, not 1/macheps. 
//The threshold should be compatible with other thresholds used in the 
//code. We set  OFFDIG + OFFDIG2 <= .999 =: RELCOND, it corresponds 
//to losing at most 3 decimal digits: 1 / (1 - (OFFDIG + OFFDIG2)) <= 1000 
//instead of the current OFFDIG + OFFDIG2 < 1 
    yesrel = MTRUE;
    offdig = Zero;
    tmp = sqrt(abs(d[1]));
    if (tmp < rmin) {
	yesrel = MFALSE;
    }
    if (!yesrel) {
	goto L11;
    }
    for (i = 1; i < n; i++) {
	tmp2 = sqrt(abs(d[i]));
	if (tmp2 < rmin) {
	    yesrel = MFALSE;
	}
	if (!yesrel) {
	    goto L11;
	}
	offdig2 = abs(e[i - 1]) / (tmp * tmp2);
	if (offdig + offdig2 >= .999) {
	    yesrel = MFALSE;
	}
	if (!yesrel) {
	    goto L11;
	}
	tmp = tmp2;
	offdig = offdig2;

    }
  L11:
    if (yesrel) {
	*info = 0;
	return;
    } else {
    }
//*** MORE TO BE IMPLEMENTED *** 
//Test if the lower bidiagonal matrix L from T = L D L^T 
//(zero shift facto) is well conditioned 
//Test if the upper bidiagonal matrix U from T = U D U^T 
//(zero shift facto) is well conditioned. 
//In this case, the matrix needs to be flipped and, at the end 
//of the eigenvector computation, the flip needs to be applied 
//to the computed eigenvectors (and the support) 
    return;
}
