/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaqr1.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rlaqr1(INTEGER n, REAL * h, INTEGER ldh, REAL sr1, REAL si1, REAL sr2, REAL si2, REAL * v)
{
    REAL s, h21s, h31s;
    REAL Zero = 0.0;

    if (n == 2) {
	s = (abs(h[0] - sr2)) + abs(si2) + (abs(h[ldh + 2]));
	if (s == Zero) {
	    v[0] = Zero;
	    v[1] = Zero;
	} else {
	    h21s = h[ldh + 2] / s;
	    v[1] = h21s * h[(ldh * 2) + 1] + (h[ldh + 1] - sr1) * ((h[ldh + 1] - sr2) / s) - si1 * (si2 / s);
	    v[2] = h21s * (h[ldh + 1] + h[(ldh * 2) + 2] - sr1 - sr2);
	}
    } else {
	s = (abs(h[ldh + 1] - sr2)) + abs(si2) + (abs(h[ldh + 2])) + (abs(h[ldh + 3]));
	if (s == Zero) {
	    v[1] = Zero;
	    v[2] = Zero;
	    v[3] = Zero;
	} else {
	    h21s = h[ldh + 2] / s;
	    h31s = h[ldh + 3] / s;
	    v[1] = (h[ldh + 1] - sr1) * ((h[ldh + 1] - sr2) / s)
		- si1 * (si2 / s) + h[(ldh * 2) + 1] * h21s + h[ldh * 3 + 1] * h31s;
	    v[2] = h21s * (h[ldh + 1] + h[(ldh * 2) + 2] - sr1 - sr2) + h[ldh * 3 + 2] * h31s;
	    v[3] = h31s * (h[ldh + 1] + h[ldh * 3 + 3] - sr1 - sr2) + h21s * h[ldh * 2 + 3];
	}
    }
    return;
}
