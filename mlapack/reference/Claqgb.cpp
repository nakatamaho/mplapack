/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Claqgb.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Claqgb(INTEGER m, INTEGER n, INTEGER kl, INTEGER ku, COMPLEX * AB, INTEGER ldab, REAL * r, REAL * c, REAL rowcnd, REAL colcnd, REAL amax, char *equed)
{
    INTEGER i, j;
    REAL cj, large, small;
    REAL One = 1.0;

//Quick return if possible
    if (m <= 0 || n <= 0) {
	*equed = 'N';
	return;
    }
//Initialize LARGE and SMALL.
    small = Rlamch("Safe minimum") / Rlamch("Precision");
    large = One / small;

    if (rowcnd >= .1 && amax >= small && amax <= large) {
//No row scaling
	if (colcnd >= .1) {
//No column scaling
	    *equed = 'N';
	} else {
//Column scaling
	    for (j = 0; j < n; j++) {
		cj = c[j];
		for (i = max((INTEGER) 1, j - ku); i <= min(m, j + kl); i++) {
		    AB[ku + 1 + i - j + j * ldab] = cj * AB[ku + 1 + i - j + j * ldab];
		}

	    }
	    *equed = 'C';
	}
    } else if (colcnd >= .1) {
//Row scaling, no column scaling
	for (j = 0; j < n; j++) {
	    for (i = max((INTEGER) 1, j - ku); i <= min(m, j + kl); i++) {
		AB[ku + 1 + i - j + j * ldab] = r[i] * AB[ku + 1 + i - j + j * ldab];

	    }

	}
	*equed = 'R';
    } else {
//Row and column scaling
	for (j = 0; j < n; j++) {
	    cj = c[j];
	    for (i = max((INTEGER) 1, j - ku); i <= min(m, j + kl); i++) {
		AB[ku + 1 + i - j + j * ldab] = cj * r[i] * AB[ku + 1 + i - j + j * ldab];
	    }
	}
	*equed = 'B';
    }
    return;
}
