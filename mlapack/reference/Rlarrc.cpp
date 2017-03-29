/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlarrc.cpp,v 1.7 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rlarrc(const char *jobt, INTEGER n, REAL vl, REAL vu, REAL * d, REAL * e, REAL pivmin, INTEGER * eigcnt, INTEGER * lcnt, INTEGER * rcnt, INTEGER * info)
{
    INTEGER i;
    REAL sl, su, tmp, tmp2;
    INTEGER matt;
    REAL lpivot, rpivot;
    REAL Zero = 0.0;

    *info = 0;
    lcnt = 0;
    *rcnt = 0;
    *eigcnt = 0;
    matt = Mlsame(jobt, "T");
    if (matt) {

//Sturm sequence count on T
	lpivot = d[0] - vl;
	rpivot = d[0] - vu;
	if (lpivot <= Zero) {
	    ++(lcnt);
	}
	if (rpivot <= Zero) {
	    ++(*rcnt);
	}
	for (i = 0; i < n - 1; i++) {
	    tmp = e[i] * e[i];
	    lpivot = d[i + 1] - vl - tmp / lpivot;
	    rpivot = d[i + 1] - vu - tmp / rpivot;
	    if (lpivot <= Zero) {
		++(lcnt);
	    }
	    if (rpivot <= Zero) {
		++(*rcnt);
	    }

	}
    } else {
//Sturm sequence count on L D L^T
	sl = -vl;
	su = -vu;
	for (i = 0; i < n - 1; i++) {
	    lpivot = d[i] + sl;
	    rpivot = d[i] + su;
	    if (lpivot <= Zero) {
		++(lcnt);
	    }
	    if (rpivot <= Zero) {
		++(*rcnt);
	    }
	    tmp = e[i] * d[i] * e[i];

	    tmp2 = tmp / lpivot;
	    if (tmp2 == Zero) {
		sl = tmp - vl;
	    } else {
		sl = sl * tmp2 - vl;
	    }

	    tmp2 = tmp / rpivot;
	    if (tmp2 == Zero) {
		su = tmp - vu;
	    } else {
		su = su * tmp2 - vu;
	    }

	}
	lpivot = d[n] + sl;
	rpivot = d[n] + su;
	if (lpivot <= Zero) {
	    ++(lcnt);
	}
	if (rpivot <= Zero) {
	    ++(*rcnt);
	}
    }
    *eigcnt = *rcnt - *lcnt;
    return;
}
