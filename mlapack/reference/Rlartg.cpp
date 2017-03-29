/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlartg.cpp,v 1.8 2010/08/07 04:48:33 nakatamaho Exp $ 
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
#include <stdio.h>		//for printf

void Rlartg(REAL f, REAL g, REAL * cs, REAL * sn, REAL * r)
{
    REAL f1, g1;
    INTEGER i, count;
    REAL safmin;
    REAL safmn2;
    REAL safmx2, eps, scale;
    REAL Zero = 0.0, One = 1.0;

    safmin = Rlamch("S");
    eps = Rlamch("E");
// SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) / LOG( DLAMCH( 'B' ) ) / TWO );
//        ~ 2^(ln(safmin/eps) / 2ln2 ) (dlamchB=2)  = sqrt(safmin/eps).
    safmn2 = sqrt(safmin / eps);
    safmx2 = 1.0 / safmn2;

    if (g == Zero) {
	*cs = One;
	*sn = Zero;
	*r = f;
    } else if (f == Zero) {
	*cs = Zero;
	*sn = One;
	*r = g;
    } else {
	f1 = f;
	g1 = g;
	scale = max(abs(f1), abs(g1));
	count = 0;
	if (scale >= safmx2) {
	    printf("#XXX Rlartg :1: not yet implemented.\n");
	    while (1) {
		count++;
		f1 = f1 * safmn2;
		g1 = g1 * safmn2;
		scale = max(abs(f1), abs(g1));
		if (scale >= safmx2)
		    continue;

		*r = sqrt(f1 * f1 + g1 * g1);
		*cs = f1 / (*r);
		*sn = g1 / (*r);
		for (i = 0; i < count; i++) {
		    *r = (*r) * safmx2;
		}
		break;
	    }
	} else if (scale <= safmn2) {
	    printf("#XXX Rlartg :3:very well tested. \n");
	    while (1) {
		count++;
		f1 = f1 * safmx2;
		g1 = g1 * safmn2;
		scale = max(abs(f1), abs(g1));
		if (scale >= safmx2)
		    continue;
		*r = sqrt(f1 * f1 + g1 * g1);
		*cs = f1 / (*r);
		*sn = g1 / (*r);
		for (i = 0; i < count; i++) {
		    *r = (*r) * safmx2;
		}
		break;
	    }
	} else {
	    *r = sqrt(f1 * f1 + g1 * g1);
	    *cs = f1 / (*r);
	    *sn = g1 / (*r);
	}
	if (abs(f) > abs(g) && (*cs) < Zero) {
	    *cs = -(*cs);
	    *sn = -(*sn);
	    *r = -(*r);
	}
    }
    return;
}
