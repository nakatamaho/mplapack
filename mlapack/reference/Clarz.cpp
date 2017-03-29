/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Clarz.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Clarz(const char *side, INTEGER m, INTEGER n, INTEGER l, COMPLEX * v, INTEGER incv, COMPLEX tau, COMPLEX * c, INTEGER ldc, COMPLEX * work)
{
    REAL Zero = 0.0, One = 1.0;

    if (Mlsame(side, "L")) {
//Form  H * C
	if (tau != Zero) {
//w( 1:n ) = conjg( C( 1, 1:n ) )
	    Ccopy(n, &c[0], ldc, &work[0], 1);
	    Clacgv(n, &work[0], 1);
//w( 1:n ) = conjg( w( 1:n ) + C( m-l+1:m, 1:n )' * v( 1:l ) )
	    Cgemv("Conjugate transpose", l, n, One, &c[m - l + 1 + ldc], ldc, &v[1], incv, One, &work[0], 1);
	    Clacgv(n, &work[0], 1);
//C( 1, 1:n ) = C( 1, 1:n ) - tau * w( 1:n )
	    Caxpy(n, -tau, &work[0], 1, &c[0], ldc);
//C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
//                    tau * v( 1:l ) * conjg( w( 1:n )' )
	    Cgeru(l, n, -tau, &v[1], incv, &work[0], 1, &c[m - l + 1 + ldc], ldc);
	}
    } else {
//Form  C * H
	if (tau != Zero) {
//w( 1:m ) = C( 1:m, 1 )
	    Ccopy(m, &c[0], 1, &work[0], 1);
//w( 1:m ) = w( 1:m ) + C( 1:m, n-l+1:n, 1:n ) * v( 1:l )
	    Cgemv("No transpose", m, l, One, &c[(n - l + 1) * ldc + 1], ldc, &v[1], incv, One, &work[0], 1);
//C( 1:m, 1 ) = C( 1:m, 1 ) - tau * w( 1:m )
	    Caxpy(m, -(tau), &work[0], 1, &c[0], 1);
//C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
//                    tau * w( 1:m ) * v( 1:l )'
	    Cgerc(m, l, -tau, &work[0], 1, &v[1], incv, &c[(n - l + 1) * ldc + 1], ldc);
	}
    }
    return;
}
