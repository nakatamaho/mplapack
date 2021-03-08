/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlarfx.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void Rlarfx(const char *side, INTEGER m, INTEGER n, REAL * v, REAL tau, REAL * c, INTEGER ldc, REAL * work)
{
    INTEGER j;
    REAL t1, t2, t3, t4, t5, t6, t7, t8, t9, v1, v2, v3, v4, v5, v6, v7, v8, v9, t10, v10, sum;
    REAL Zero = 0.0, One = 1.0;

    if (tau == Zero) {
	return;
    }
    if (Mlsame(side, "L")) {
//Form  H * C, where H has order m.
	switch (m) {
	case 1:
	    goto L10;
	case 2:
	    goto L30;
	case 3:
	    goto L50;
	case 4:
	    goto L70;
	case 5:
	    goto L90;
	case 6:
	    goto L110;
	case 7:
	    goto L130;
	case 8:
	    goto L150;
	case 9:
	    goto L170;
	case 10:
	    goto L190;
	}
//Code for general M
//w := C'*v
	Rgemv("Transpose", m, n, One, &c[0], ldc, &v[1], 1, Zero, &work[0], 1);
//C := C - tau * v * w'
	Rger(m, n, -(tau), &v[0], 1, &work[0], 1, &c[0], ldc);
	goto L410;
      L10:
//Special code for 1 x 1 Householder

	t1 = One - tau * v[1] * v[1];
	for (j = 0; j < n; j++) {
	    c[j * ldc + 1] = t1 * c[j * ldc + 1];
	}
	goto L410;
      L30:
//Special code for 2 x 2 Householder
	v1 = v[1];
	t1 = tau * v1;
	v2 = v[2];
	t2 = tau * v2;
	for (j = 0; j < n; j++) {
	    sum = v1 * c[j * ldc + 1] + v2 * c[j * ldc + 2];
	    c[j * ldc + 1] -= sum * t1;
	    c[j * ldc + 2] -= sum * t2;
	}
	goto L410;
      L50:
//Special code for 3 x 3 Householder
	v1 = v[1];
	t1 = tau * v1;
	v2 = v[2];
	t2 = tau * v2;
	v3 = v[3];
	t3 = tau * v3;
	for (j = 0; j < n; j++) {
	    sum = v1 * c[j * ldc + 1] + v2 * c[j * ldc + 2] + v3 * c[j * ldc + 3];
	    c[j * ldc + 1] = c[j * ldc + 1] - sum * t1;
	    c[j * ldc + 2] = c[j * ldc + 2] - sum * t2;
	    c[j * ldc + 3] = c[j * ldc + 3] - sum * t3;
	}
	goto L410;
      L70:
//Special code for 4 x 4 Householder
	v1 = v[1];
	t1 = tau * v1;
	v2 = v[2];
	t2 = tau * v2;
	v3 = v[3];
	t3 = tau * v3;
	v4 = v[4];
	t4 = tau * v4;
	for (j = 0; j < n; j++) {
	    sum = v1 * c[j * ldc + 1] + v2 * c[j * ldc + 2] + v3 * c[j * ldc + 3] + v4 * c[j * ldc + 4];
	    c[j * ldc + 1] = c[j * ldc + 1] - sum * t1;
	    c[j * ldc + 2] = c[j * ldc + 2] - sum * t2;
	    c[j * ldc + 3] = c[j * ldc + 3] - sum * t3;
	    c[j * ldc + 4] = c[j * ldc + 4] - sum * t4;
	}
	goto L410;
      L90:
//Special code for 5 x 5 Householder
	v1 = v[1];
	t1 = tau * v1;
	v2 = v[2];
	t2 = tau * v2;
	v3 = v[3];
	t3 = tau * v3;
	v4 = v[4];
	t4 = tau * v4;
	v5 = v[5];
	t5 = tau * v5;
	for (j = 0; j < n; j++) {
	    sum = v1 * c[j * ldc + 1] + v2 * c[j * ldc + 2] + v3 * c[j * ldc + 3] + v4 * c[j * ldc + 4] + v5 * c[j * ldc + 5];
	    c[j * ldc + 1] = c[j * ldc + 1] - sum * t1;
	    c[j * ldc + 2] = c[j * ldc + 2] - sum * t2;
	    c[j * ldc + 3] = c[j * ldc + 3] - sum * t3;
	    c[j * ldc + 4] = c[j * ldc + 4] - sum * t4;
	    c[j * ldc + 5] = c[j * ldc + 5] - sum * t5;

	}
	goto L410;
      L110:
//Special code for 6 x 6 Householder
	v1 = v[1];
	t1 = tau * v1;
	v2 = v[2];
	t2 = tau * v2;
	v3 = v[3];
	t3 = tau * v3;
	v4 = v[4];
	t4 = tau * v4;
	v5 = v[5];
	t5 = tau * v5;
	v6 = v[6];
	t6 = tau * v6;
	for (j = 0; j < n; j++) {
	    sum = v1 * c[j * ldc + 1] + v2 * c[j * ldc + 2] + v3 * c[j * ldc + 3] + v4 * c[j * ldc + 4] + v5 * c[j * ldc + 5] + v6 * c[j * ldc + 6];
	    c[j * ldc + 1] = c[j * ldc + 1] - sum * t1;
	    c[j * ldc + 2] = c[j * ldc + 2] - sum * t2;
	    c[j * ldc + 3] = c[j * ldc + 3] - sum * t3;
	    c[j * ldc + 4] = c[j * ldc + 4] - sum * t4;
	    c[j * ldc + 5] = c[j * ldc + 5] - sum * t5;
	    c[j * ldc + 6] = c[j * ldc + 6] - sum * t6;
	}
	goto L410;
      L130:
//Special code for 7 x 7 Householder
	v1 = v[1];
	t1 = tau * v1;
	v2 = v[2];
	t2 = tau * v2;
	v3 = v[3];
	t3 = tau * v3;
	v4 = v[4];
	t4 = tau * v4;
	v5 = v[5];
	t5 = tau * v5;
	v6 = v[6];
	t6 = tau * v6;
	v7 = v[7];
	t7 = tau * v7;
	for (j = 0; j < n; j++) {
	    sum = v1 * c[j * ldc + 1] + v2 * c[j * ldc + 2] + v3 * c[j * ldc + 3] + v4 * c[j * ldc + 4] + v5 * c[j * ldc + 5] + v6 * c[j * ldc + 6] + v7 * c[j * ldc + 7];
	    c[j * ldc + 1] = c[j * ldc + 1] - sum * t1;
	    c[j * ldc + 2] = c[j * ldc + 2] - sum * t2;
	    c[j * ldc + 3] = c[j * ldc + 3] - sum * t3;
	    c[j * ldc + 4] = c[j * ldc + 4] - sum * t4;
	    c[j * ldc + 5] = c[j * ldc + 5] - sum * t5;
	    c[j * ldc + 6] = c[j * ldc + 6] - sum * t6;
	    c[j * ldc + 7] = c[j * ldc + 7] - sum * t7;
	}
	goto L410;
      L150:
//Special code for 8 x 8 Householder
	v1 = v[1];
	t1 = tau * v1;
	v2 = v[2];
	t2 = tau * v2;
	v3 = v[3];
	t3 = tau * v3;
	v4 = v[4];
	t4 = tau * v4;
	v5 = v[5];
	t5 = tau * v5;
	v6 = v[6];
	t6 = tau * v6;
	v7 = v[7];
	t7 = tau * v7;
	v8 = v[8];
	t8 = tau * v8;
	for (j = 0; j < n; j++) {
	    sum = v1 * c[j * ldc + 1] + v2 * c[j * ldc + 2] + v3 *
		c[j * ldc + 3] + v4 * c[j * ldc + 4] + v5 * c[j * ldc + 5] + v6 * c[j * ldc + 6] + v7 * c[j * ldc + 7] + v8 * c[j * ldc + 8];
	    c[j * ldc + 1] = c[j * ldc + 1] - sum * t1;
	    c[j * ldc + 2] = c[j * ldc + 2] - sum * t2;
	    c[j * ldc + 3] = c[j * ldc + 3] - sum * t3;
	    c[j * ldc + 4] = c[j * ldc + 4] - sum * t4;
	    c[j * ldc + 5] = c[j * ldc + 5] - sum * t5;
	    c[j * ldc + 6] = c[j * ldc + 6] - sum * t6;
	    c[j * ldc + 7] = c[j * ldc + 7] - sum * t7;
	    c[j * ldc + 8] = c[j * ldc + 8] - sum * t8;
	}
	goto L410;
      L170:
//Special code for 9 x 9 Householder
	v1 = v[1];
	t1 = tau * v1;
	v2 = v[2];
	t2 = tau * v2;
	v3 = v[3];
	t3 = tau * v3;
	v4 = v[4];
	t4 = tau * v4;
	v5 = v[5];
	t5 = tau * v5;
	v6 = v[6];
	t6 = tau * v6;
	v7 = v[7];
	t7 = tau * v7;
	v8 = v[8];
	t8 = tau * v8;
	v9 = v[9];
	t9 = tau * v9;
	for (j = 0; j < n; j++) {
	    sum = v1 * c[j * ldc + 1] + v2 * c[j * ldc + 2] + v3 *
		c[j * ldc + 3] + v4 * c[j * ldc + 4] + v5 * c[j * ldc + 5] + v6 * c[j * ldc + 6] + v7 * c[j * ldc + 7] + v8 * c[j * ldc + 8] + v9 * c[j * ldc + 9];
	    c[j * ldc + 1] = c[j * ldc + 1] - sum * t1;
	    c[j * ldc + 2] = c[j * ldc + 2] - sum * t2;
	    c[j * ldc + 3] = c[j * ldc + 3] - sum * t3;
	    c[j * ldc + 4] = c[j * ldc + 4] - sum * t4;
	    c[j * ldc + 5] = c[j * ldc + 5] - sum * t5;
	    c[j * ldc + 6] = c[j * ldc + 6] - sum * t6;
	    c[j * ldc + 7] = c[j * ldc + 7] - sum * t7;
	    c[j * ldc + 8] = c[j * ldc + 8] - sum * t8;
	    c[j * ldc + 9] = c[j * ldc + 9] - sum * t9;
	}
	goto L410;
      L190:
//Special code for 10 x 10 Householder
	v1 = v[1];
	t1 = tau * v1;
	v2 = v[2];
	t2 = tau * v2;
	v3 = v[3];
	t3 = tau * v3;
	v4 = v[4];
	t4 = tau * v4;
	v5 = v[5];
	t5 = tau * v5;
	v6 = v[6];
	t6 = tau * v6;
	v7 = v[7];
	t7 = tau * v7;
	v8 = v[8];
	t8 = tau * v8;
	v9 = v[9];
	t9 = tau * v9;
	v10 = v[10];
	t10 = tau * v10;
	for (j = 0; j < n; j++) {
	    sum = v1 * c[j * ldc + 1] + v2 * c[j * ldc + 2] + v3 *
		c[j * ldc + 3] + v4 * c[j * ldc + 4] + v5 * c[j * ldc + 5] +
		v6 * c[j * ldc + 6] + v7 * c[j * ldc + 7] + v8 * c[j * ldc + 8] + v9 * c[j * ldc + 9] + v10 * c[j * ldc + 10];
	    c[j * ldc + 1] = c[j * ldc + 1] - sum * t1;
	    c[j * ldc + 2] = c[j * ldc + 2] - sum * t2;
	    c[j * ldc + 3] = c[j * ldc + 3] - sum * t3;
	    c[j * ldc + 4] = c[j * ldc + 4] - sum * t4;
	    c[j * ldc + 5] = c[j * ldc + 5] - sum * t5;
	    c[j * ldc + 6] = c[j * ldc + 6] - sum * t6;
	    c[j * ldc + 7] = c[j * ldc + 7] - sum * t7;
	    c[j * ldc + 8] = c[j * ldc + 8] - sum * t8;
	    c[j * ldc + 9] = c[j * ldc + 9] - sum * t9;
	    c[j * ldc + 10] = c[j * ldc + 10] - sum * t10;
	}
	goto L410;
    } else {
//Form  C * H, where H has order n.
	switch (n) {
	case 1:
	    goto L210;
	case 2:
	    goto L230;
	case 3:
	    goto L250;
	case 4:
	    goto L270;
	case 5:
	    goto L290;
	case 6:
	    goto L310;
	case 7:
	    goto L330;
	case 8:
	    goto L350;
	case 9:
	    goto L370;
	case 10:
	    goto L390;
	}
//Code for general N
//w := C * v
	Rgemv("No transpose", m, n, One, &c[0], ldc, &v[1], 1, Zero, &work[0], 1);
//C := C - tau * w * v'
	Rger(m, n, -(tau), &work[0], 1, &v[1], 1, &c[0], ldc);
	goto L410;
      L210:
//Special code for 1 x 1 Householder
	t1 = One - tau * v[1] * v[1];
	for (j = 0; j < m; j++) {
	    c[j + ldc] = t1 * c[j + ldc];
	}
	goto L410;
      L230:
//Special code for 2 x 2 Householder
	v1 = v[1];
	t1 = tau * v1;
	v2 = v[2];
	t2 = tau * v2;
	for (j = 0; j < m; j++) {
	    sum = v1 * c[j + ldc] + v2 * c[j + ldc * 2];
	    c[j + ldc] = c[j + ldc] - sum * t1;
	    c[j + ldc * 2] = c[j + ldc * 2] - sum * t2;
	}
	goto L410;
      L250:
//Special code for 3 x 3 Householder
	v1 = v[1];
	t1 = tau * v1;
	v2 = v[2];
	t2 = tau * v2;
	v3 = v[3];
	t3 = tau * v3;
	for (j = 0; j < m; j++) {
	    sum = v1 * c[j + ldc] + v2 * c[j + ldc * 2] + v3 * c[j + ldc * 3];
	    c[j + ldc] = c[j + ldc] - sum * t1;
	    c[j + ldc * 2] = c[j + ldc * 2] - sum * t2;
	    c[j + ldc * 3] = c[j + ldc * 3] - sum * t3;
	}
	goto L410;
      L270:
//Special code for 4 x 4 Householder
	v1 = v[1];
	t1 = tau * v1;
	v2 = v[2];
	t2 = tau * v2;
	v3 = v[3];
	t3 = tau * v3;
	v4 = v[4];
	t4 = tau * v4;
	for (j = 0; j < m; j++) {
	    sum = v1 * c[j + ldc] + v2 * c[j + ldc * 2] + v3 * c[j + ldc * 3] + v4 * c[j + ldc * 4];
	    c[j + ldc] = c[j + ldc] - sum * t1;
	    c[j + ldc * 2] = c[j + ldc * 2] - sum * t2;
	    c[j + ldc * 3] = c[j + ldc * 3] - sum * t3;
	    c[j + ldc * 4] = c[j + ldc * 4] - sum * t4;

	}
	goto L410;
      L290:

/*        Special code for 5 x 5 Householder */

	v1 = v[1];
	t1 = tau * v1;
	v2 = v[2];
	t2 = tau * v2;
	v3 = v[3];
	t3 = tau * v3;
	v4 = v[4];
	t4 = tau * v4;
	v5 = v[5];
	t5 = tau * v5;
	for (j = 0; j < m; j++) {
	    sum = v1 * c[j + ldc] + v2 * c[j + ldc * 2] + v3 * c[j + ldc * 3] + v4 * c[j + ldc * 4] + v5 * c[j + ldc * 5];
	    c[j + ldc] = c[j + ldc] - sum * t1;
	    c[j + ldc * 2] = c[j + ldc * 2] - sum * t2;
	    c[j + ldc * 3] = c[j + ldc * 3] - sum * t3;
	    c[j + ldc * 4] = c[j + ldc * 4] - sum * t4;
	    c[j + ldc * 5] = c[j + ldc * 5] - sum * t5;
	}
	goto L410;
      L310:
//Special code for 6 x 6 Householder
	v1 = v[1];
	t1 = tau * v1;
	v2 = v[2];
	t2 = tau * v2;
	v3 = v[3];
	t3 = tau * v3;
	v4 = v[4];
	t4 = tau * v4;
	v5 = v[5];
	t5 = tau * v5;
	v6 = v[6];
	t6 = tau * v6;
	for (j = 0; j < m; j++) {
	    sum = v1 * c[j + ldc] + v2 * c[j + ldc * 2] + v3 * c[j + ldc * 3] + v4 * c[j + ldc * 4] + v5 * c[j + ldc * 5] + v6 * c[j + ldc * 6];
	    c[j + ldc] = c[j + ldc] - sum * t1;
	    c[j + ldc * 2] = c[j + ldc * 2] - sum * t2;
	    c[j + ldc * 3] = c[j + ldc * 3] - sum * t3;
	    c[j + ldc * 4] = c[j + ldc * 4] - sum * t4;
	    c[j + ldc * 5] = c[j + ldc * 5] - sum * t5;
	    c[j + ldc * 6] = c[j + ldc * 6] - sum * t6;
	}
	goto L410;
      L330:
//Special code for 7 x 7 Householder
	v1 = v[1];
	t1 = tau * v1;
	v2 = v[2];
	t2 = tau * v2;
	v3 = v[3];
	t3 = tau * v3;
	v4 = v[4];
	t4 = tau * v4;
	v5 = v[5];
	t5 = tau * v5;
	v6 = v[6];
	t6 = tau * v6;
	v7 = v[7];
	t7 = tau * v7;
	for (j = 0; j < m; j++) {
	    sum = v1 * c[j + ldc] + v2 * c[j + ldc * 2] + v3 * c[j + ldc * 3] + v4 * c[j + ldc * 4] + v5 * c[j + ldc * 5] + v6 * c[j + ldc * 6] + v7 * c[j + ldc * 7];
	    c[j + ldc] = c[j + ldc] - sum * t1;
	    c[j + ldc * 2] = c[j + ldc * 2] - sum * t2;
	    c[j + ldc * 3] = c[j + ldc * 3] - sum * t3;
	    c[j + ldc * 4] = c[j + ldc * 4] - sum * t4;
	    c[j + ldc * 5] = c[j + ldc * 5] - sum * t5;
	    c[j + ldc * 6] = c[j + ldc * 6] - sum * t6;
	    c[j + ldc * 7] = c[j + ldc * 7] - sum * t7;
	}
	goto L410;
      L350:
//Special code for 8 x 8 Householder
	v1 = v[1];
	t1 = tau * v1;
	v2 = v[2];
	t2 = tau * v2;
	v3 = v[3];
	t3 = tau * v3;
	v4 = v[4];
	t4 = tau * v4;
	v5 = v[5];
	t5 = tau * v5;
	v6 = v[6];
	t6 = tau * v6;
	v7 = v[7];
	t7 = tau * v7;
	v8 = v[8];
	t8 = tau * v8;
	for (j = 0; j < m; j++) {
	    sum = v1 * c[j + ldc] + v2 * c[j + ldc * 2] + v3 *
		c[j + ldc * 3] + v4 * c[j + ldc * 4] + v5 * c[j + ldc * 5] + v6 * c[j + ldc * 6] + v7 * c[j + ldc * 7] + v8 * c[j + ldc * 8];
	    c[j + ldc] = c[j + ldc] - sum * t1;
	    c[j + ldc * 2] = c[j + ldc * 2] - sum * t2;
	    c[j + ldc * 3] = c[j + ldc * 3] - sum * t3;
	    c[j + ldc * 4] = c[j + ldc * 4] - sum * t4;
	    c[j + ldc * 5] = c[j + ldc * 5] - sum * t5;
	    c[j + ldc * 6] = c[j + ldc * 6] - sum * t6;
	    c[j + ldc * 7] = c[j + ldc * 7] - sum * t7;
	    c[j + ldc * 8] = c[j + ldc * 8] - sum * t8;
	}
	goto L410;
      L370:
//Special code for 9 x 9 Householder
	v1 = v[1];
	t1 = tau * v1;
	v2 = v[2];
	t2 = tau * v2;
	v3 = v[3];
	t3 = tau * v3;
	v4 = v[4];
	t4 = tau * v4;
	v5 = v[5];
	t5 = tau * v5;
	v6 = v[6];
	t6 = tau * v6;
	v7 = v[7];
	t7 = tau * v7;
	v8 = v[8];
	t8 = tau * v8;
	v9 = v[9];
	t9 = tau * v9;
	for (j = 0; j < m; j++) {
	    sum = v1 * c[j + ldc] + v2 * c[j + ldc * 2] + v3 *
		c[j + ldc * 3] + v4 * c[j + ldc * 4] + v5 * c[j + ldc * 5] + v6 * c[j + ldc * 6] + v7 * c[j + ldc * 7] + v8 * c[j + ldc * 8] + v9 * c[j + ldc * 9];
	    c[j + ldc] = c[j + ldc] - sum * t1;
	    c[j + ldc * 2] = c[j + ldc * 2] - sum * t2;
	    c[j + ldc * 3] = c[j + ldc * 3] - sum * t3;
	    c[j + ldc * 4] = c[j + ldc * 4] - sum * t4;
	    c[j + ldc * 5] = c[j + ldc * 5] - sum * t5;
	    c[j + ldc * 6] = c[j + ldc * 6] - sum * t6;
	    c[j + ldc * 7] = c[j + ldc * 7] - sum * t7;
	    c[j + ldc * 8] = c[j + ldc * 8] - sum * t8;
	    c[j + ldc * 9] = c[j + ldc * 9] - sum * t9;

	}
	goto L410;
      L390:
//Special code for 10 x 10 Householder
	v1 = v[1];
	t1 = tau * v1;
	v2 = v[2];
	t2 = tau * v2;
	v3 = v[3];
	t3 = tau * v3;
	v4 = v[4];
	t4 = tau * v4;
	v5 = v[5];
	t5 = tau * v5;
	v6 = v[6];
	t6 = tau * v6;
	v7 = v[7];
	t7 = tau * v7;
	v8 = v[8];
	t8 = tau * v8;
	v9 = v[9];
	t9 = tau * v9;
	v10 = v[10];
	t10 = tau * v10;
	for (j = 0; j < m; j++) {
	    sum = v1 * c[j + ldc] + v2 * c[j + ldc * 2] + v3 *
		c[j + ldc * 3] + v4 * c[j + ldc * 4] + v5 *
		c[j + ldc * 5] + v6 * c[j + ldc * 6] + v7 * c[j + ldc * 7] + v8 * c[j + ldc * 8] + v9 * c[j + ldc * 9] + v10 * c[j + ldc * 10];
	    c[j + ldc] = c[j + ldc] - sum * t1;
	    c[j + ldc * 2] = c[j + ldc * 2] - sum * t2;
	    c[j + ldc * 3] = c[j + ldc * 3] - sum * t3;
	    c[j + ldc * 4] = c[j + ldc * 4] - sum * t4;
	    c[j + ldc * 5] = c[j + ldc * 5] - sum * t5;
	    c[j + ldc * 6] = c[j + ldc * 6] - sum * t6;
	    c[j + ldc * 7] = c[j + ldc * 7] - sum * t7;
	    c[j + ldc * 8] = c[j + ldc * 8] - sum * t8;
	    c[j + ldc * 9] = c[j + ldc * 9] - sum * t9;
	    c[j + ldc * 10] = c[j + ldc * 10] - sum * t10;

	}
	goto L410;
    }
  L410:
    return;
}
