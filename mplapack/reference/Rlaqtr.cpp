/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaqtr.cpp,v 1.7 2010/08/07 04:48:32 nakatamaho Exp $ 
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

#define MFALSE 0
#define MTRUE  1

void Rlaqtr(INTEGER ltran, INTEGER lreal, INTEGER n, REAL * t, INTEGER ldt, REAL * B, REAL w, REAL * scale, REAL * x, REAL * work, INTEGER * info)
{
    REAL d[4];
    INTEGER i, j, k;
    REAL v[4], z;
    INTEGER j1, j2, n1, n2;
    REAL si, xj, sr, rec, eps, tjj, tmp;
    INTEGER ierr;
    REAL smin, xmax;
    INTEGER jnext;
    REAL sminw, xnorm;
    REAL scaloc;
    REAL bignum;
    INTEGER notran;
    REAL smlnum;
    REAL One = 1.0, Zero = 0.0;
    REAL mtemp1, mtemp2, mtemp3;

    notran = !(ltran);
    *info = 0;
//Quick return if possible
    if (n == 0) {
	return;
    }
//Set constants to control overflow
    eps = Rlamch("P");
    smlnum = Rlamch("S") / eps;
    bignum = One / smlnum;
    xnorm = Rlange("M", n, n, &t[0], ldt, d);
    if (!(lreal)) {
	mtemp1 = xnorm, mtemp2 = abs(w);
	mtemp3 = max(mtemp1, mtemp2);
	mtemp2 = Rlange("M", n, 1, &B[1], n, d);
	xnorm = max(mtemp3, mtemp2);
    }
    mtemp1 = smlnum;
    mtemp2 = eps * xnorm;
    smin = max(mtemp1, mtemp2);

//Compute 1-norm of each column of strictly upper triangular
//part of T to control overflow in triangular solver.
    work[1] = Zero;
    for (j = 2; j <= n; j++) {
	work[j] = Rasum(j - 1, &t[j * ldt + 1], 1);
    }
    if (!(lreal)) {
	for (i = 1; i < n; i++) {
	    work[i] = work[i] + abs(B[i]);
	}
    }
    n2 = n << 1;
    n1 = n;
    if (!(lreal)) {
	n1 = n2;
    }
    k = iRamax(n1, &x[0], 1);
    xmax = abs(x[k]);
    *scale = One;

    if (xmax > bignum) {
	*scale = bignum / xmax;
	Rscal(n1, *scale, &x[0], 1);
	xmax = bignum;
    }

    if (lreal) {
	if (notran) {
//Solve Tp = scale*c
	    jnext = n;
	    for (j = n; j >= 1; j--) {
		if (j > jnext) {
		    goto L30;
		}
		j1 = j;
		j2 = j;
		jnext = j - 1;
		if (j > 1) {
		    if (t[j + (j - 1) * ldt] != Zero) {
			j1 = j - 1;
			jnext = j - 2;
		    }
		}
		if (j1 == j2) {
//Meet 1 by 1 diagonal block
//Scale to avoid overflow when computing
//    x(j) = b(j)/T(j,j)
		    xj = abs(x[j1]);
		    tjj = abs(t[j1 + j1 * ldt]);
		    tmp = t[j1 + j1 * ldt];
		    if (tjj < smin) {
			tmp = smin;
			tjj = smin;
			*info = 1;
		    }
		    if (xj == Zero) {
			goto L30;
		    }
		    if (tjj < One) {
			if (xj > bignum * tjj) {
			    rec = One / xj;
			    Rscal(n, rec, &x[0], 1);
			    *scale = *scale * rec;
			    xmax = xmax * rec;
			}
		    }
		    x[j1] = x[j1] / tmp;
		    xj = abs(x[j1]);
//Scale x if necessary to avoid overflow when adding a
//multiple of column j1 of T.
		    if (xj > One) {
			rec = One / xj;
			if (work[j1] > (bignum - xmax) * rec) {
			    Rscal(n, rec, &x[0], 1);
			    *scale = *scale * rec;
			}
		    }
		    if (j1 > 1) {
			Raxpy(j1 - 1, -x[j1], &t[j1 * ldt + 1], 1, &x[0], 1);
			k = iRamax(j1 - 1, &x[0], 1);
			xmax = abs(x[k]);
		    }
		} else {
//Meet 2 by 2 diagonal block
//Call 2 by 2 linear system solve, to take
//care of possible overflow by scaling factor.
		    d[0] = x[j1];
		    d[1] = x[j2];
		    Rlaln2(MFALSE, 2, 1, smin, One, &t[j1 + j1 * ldt], ldt, One, One, d, 2, One, One, v, 2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 2;
		    }
		    if (scaloc != One) {
			Rscal(n, scaloc, &x[0], 1);
			*scale *= scaloc;
		    }
		    x[j1] = v[0];
		    x[j2] = v[1];
//Scale V(1,1) (= X(J1)) and/or V(2,1) (=X(J2))
//to avoid overflow in updating right-hand side.
		    mtemp1 = abs(v[0]), mtemp2 = abs(v[1]);
		    xj = max(mtemp1, mtemp2);
		    if (xj > One) {
			rec = One / xj;
			mtemp1 = work[j1], mtemp2 = work[j2];
			if (max(mtemp1, mtemp2) > (bignum - xmax) * rec) {
			    Rscal(n, rec, &x[0], 1);
			    *scale = *scale * rec;
			}
		    }
//Update right-hand side
		    if (j1 > 1) {
			Raxpy(j1 - 1, -x[j1], &t[j1 * ldt + 1], 1, &x[0], 1);
			Raxpy(j1 - 1, -x[j2], &t[j2 * ldt + 1], 1, &x[0], 1);
			k = iRamax(j1 - 1, &x[0], 1);
			xmax = abs(x[k]);
		    }
		}
	      L30:
		;
	    }
	} else {
//Solve T'p = scale*c
	    jnext = 1;
	    for (j = 0; j < n; j++) {
		if (j < jnext) {
		    goto L40;
		}
		j1 = j;
		j2 = j;
		jnext = j + 1;
		if (j < n) {
		    if (t[j + 1 + j * ldt] != Zero) {
			j2 = j + 1;
			jnext = j + 2;
		    }
		}
		if (j1 == j2) {
//1 by 1 diagonal block
//Scale if necessary to avoid overflow in forming the
//right-hand side element by inner product.
		    xj = abs(x[j1]);
		    if (xmax > One) {
			rec = One / xmax;
			if (work[j1] > (bignum - xj) * rec) {
			    Rscal(n, rec, &x[0], 1);
			    *scale = *scale * rec;
			    xmax = xmax * rec;
			}
		    }
		    x[j1] = x[j1] - Rdot(j1 - 1, &t[j1 * ldt + 1], 1, &x[0], 1);
		    xj = abs(x[j1]);
		    tjj = abs(t[j1 + j1 * ldt]);
		    tmp = t[j1 + j1 * ldt];
		    if (tjj < smin) {
			tmp = smin;
			tjj = smin;
			*info = 1;
		    }

		    if (tjj < One) {
			if (xj > bignum * tjj) {
			    rec = One / xj;
			    Rscal(n, rec, &x[0], 1);
			    *scale = *scale * rec;
			    xmax = xmax * rec;
			}
		    }
		    x[j1] = x[j1] / tmp;
		    mtemp1 = xmax, mtemp2 = abs(x[j1]);
		    xmax = max(mtemp1, mtemp2);
		} else {
//2 by 2 diagonal block
//Scale if necessary to avoid overflow in forming the
//right-hand side elements by inner product.
/* Computing MAX */
		    mtemp1 = abs(x[j1]), mtemp2 = abs(x[j2]);
		    xj = max(mtemp1, mtemp2);
		    if (xmax > One) {
			rec = One / xmax;
			mtemp1 = work[j2], mtemp2 = work[j1];
			if (max(mtemp1, mtemp2) > (bignum - xj) * rec) {
			    Rscal(n, rec, &x[0], 1);
			    *scale = *scale * rec;
			    xmax = xmax * rec;
			}
		    }
		    d[0] = x[j1] - Rdot(j1 - 1, &t[j1 * ldt + 1], 1, &x[0], 1);
		    d[1] = x[j2] - Rdot(j1 - 1, &t[j2 * ldt + 1], 1, &x[0], 1);
		    Rlaln2(MTRUE, 2, 1, smin, Zero, &t[j1 + j1 * ldt], ldt, Zero, Zero, d, 2, One, One, v, 2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 2;
		    }
		    if (scaloc != One) {
			Rscal(n, scaloc, &x[0], 1);
			*scale = *scale * scaloc;
		    }
		    x[j1] = v[0];
		    x[j2] = v[1];
		    mtemp1 = abs(x[j1]), mtemp2 = abs(x[j2]), mtemp3 = max(mtemp1, mtemp2);
		    xmax = max(mtemp3, xmax);
		}
	      L40:
		;
	    }
	}
    } else {
	mtemp1 = eps * abs(w);
	sminw = max(mtemp1, smin);
	if (notran) {
//Solve (T + iB)*(p+iq) = c+id
	    jnext = n;
	    for (j = n; j >= 1; j--) {
		if (j > jnext) {
		    goto L70;
		}
		j1 = j;
		j2 = j;
		jnext = j - 1;
		if (j > 1) {
		    if (t[j + (j - 1) * ldt] != Zero) {
			j1 = j - 1;
			jnext = j - 2;
		    }
		}

		if (j1 == j2) {
//1 by 1 diagonal block
//Scale if necessary to avoid overflow in division
		    z = w;
		    if (j1 == 1) {
			z = B[1];
		    }
		    xj = abs(x[j1]) + abs(x[n + j1]);
		    tjj = abs(t[j1 + j1 * ldt]) + abs(z);
		    tmp = t[j1 + j1 * ldt];
		    if (tjj < sminw) {
			tmp = sminw;
			tjj = sminw;
			*info = 1;
		    }
		    if (xj == Zero) {
			goto L70;
		    }
		    if (tjj < One) {
			if (xj > bignum * tjj) {
			    rec = One / xj;
			    Rscal(n2, rec, &x[0], 1);
			    *scale = *scale * rec;
			    xmax = xmax * rec;
			}
		    }
		    Rladiv(x[j1], x[n + j1], tmp, z, &sr, &si);
		    x[j1] = sr;
		    x[n + j1] = si;
		    xj = abs(x[j1]) + abs(x[n + j1]);
//Scale x if necessary to avoid overflow when adding a
//multiple of column j1 of T.
		    if (xj > One) {
			rec = One / xj;
			if (work[j1] > (bignum - xmax) * rec) {
			    Rscal(n2, rec, &x[0], 1);
			    *scale = *scale * rec;
			}
		    }
		    if (j1 > 1) {
			Raxpy(j1 - 1, -x[j1], &t[j1 * ldt + 1], 1, &x[0], 1);
			Raxpy(j1 - 1, -x[n + j1], &t[j1 * ldt + 1], 1, &x[n + 1], 1);

			x[0] = x[0] + B[j1] * x[n + j1];
			x[n + 1] = x[n + 1] - B[j1] * x[j1];

			xmax = Zero;
			for (k = 0; k < j1 - 1; k++) {
			    mtemp1 = xmax, mtemp2 = abs(x[k]) + abs(x[k + n]);
			    xmax = max(mtemp1, mtemp2);
			}
		    }
		} else {
//Meet 2 by 2 diagonal block
		    d[0] = x[j1];
		    d[1] = x[j2];
		    d[2] = x[n + j1];
		    d[3] = x[n + j2];
		    Rlaln2(MFALSE, 2, 2, sminw, Zero, &t[j1 + j1 * ldt], ldt, Zero, Zero, d, 2, One, -w, v, 2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 2;
		    }
		    if (scaloc != One) {
			Rscal(n * 2, scaloc, &x[0], 1);
			*scale = scaloc * *scale;
		    }
		    x[j1] = v[0];
		    x[j2] = v[1];
		    x[n + j1] = v[2];
		    x[n + j2] = v[3];
//Scale X(J1), .... to avoid overflow in
//updating right hand side.
		    mtemp1 = abs(v[0]) + abs(v[2]), mtemp2 = abs(v[1]) + abs(v[3]);
		    xj = max(mtemp1, mtemp2);
		    if (xj > One) {
			rec = One / xj;
			mtemp1 = work[j1], mtemp2 = work[j2];
			if (max(mtemp1, mtemp2) > (bignum - xmax) * rec) {
			    Rscal(n2, rec, &x[0], 1);
			    *scale = *scale * rec;
			}
		    }
//Update the right-hand side.
		    if (j1 > 1) {
			Raxpy(j1 - 1, -x[j1], &t[j1 * ldt + 1], 1, &x[0], 1);
			Raxpy(j1 - 1, -x[j2], &t[j2 * ldt + 1], 1, &x[0], 1);
			Raxpy(j1 - 1, -x[n + j1], &t[j1 * ldt + 1], 1, &x[n + 1], 1);
			Raxpy(j1 - 1, -x[n + j2], &t[j2 * ldt + 1], 1, &x[n + 1], 1);
			x[0] = x[0] + B[j1] * x[n + j1] + B[j2] * x[n + j2];
			x[n + 1] = x[n + 1] - B[j1] * x[j1] - B[j2] * x[j2];
			xmax = Zero;
			for (k = 0; k < j1 - 1; k++) {
			    mtemp1 = abs(x[k]) + abs(x[k + n]);
			    xmax = max(mtemp1, xmax);
			}
		    }
		}
	      L70:
		;
	    }
	} else {
//Solve (T + iB)'*(p+iq) = c+id
	    jnext = 1;
	    for (j = 0; j < n; j++) {
		if (j < jnext) {
		    goto L80;
		}
		j1 = j;
		j2 = j;
		jnext = j + 1;
		if (j < n) {
		    if (t[j + 1 + j * ldt] != Zero) {
			j2 = j + 1;
			jnext = j + 2;
		    }
		}
		if (j1 == j2) {
//1 by 1 diagonal block
//Scale if necessary to avoid overflow in forming the
//right-hand side element by inner product.
		    xj = abs(x[j1]) + abs(x[j1 + n]);
		    if (xmax > One) {
			rec = One / xmax;
			if (work[j1] > (bignum - xj) * rec) {
			    Rscal(n2, rec, &x[0], 1);
			    *scale = *scale * rec;
			    xmax = xmax * rec;
			}
		    }
		    x[j1] = x[j1] - Rdot(j1 - 1, &t[j1 * ldt + 1], 1, &x[0], 1);
		    x[n + j1] = x[n + j1] - Rdot(j1 - 1, &t[j1 * ldt + 1], 1, &x[n + 1], 1);
		    if (j1 > 1) {
			x[j1] = x[j1] - B[j1] * x[n + 1];
			x[n + j1] = x[n + j1] + B[j1] * x[0];
		    }
		    xj = abs(x[j1]) + abs(x[j1 + n]);
		    z = w;
		    if (j1 == 1) {
			z = B[1];
		    }
//Scale if necessary to avoid overflow in
//complex division
		    tjj = abs(t[j1 + j1 * ldt]) + abs(z);
		    tmp = t[j1 + j1 * ldt];
		    if (tjj < sminw) {
			tmp = sminw;
			tjj = sminw;
			*info = 1;
		    }

		    if (tjj < One) {
			if (xj > bignum * tjj) {
			    rec = One / xj;
			    Rscal(n2, rec, &x[0], 1);
			    *scale = *scale * rec;
			    xmax = xmax * rec;
			}
		    }
		    Rladiv(x[j1], x[n + j1], tmp, -z, &sr, &si);
		    x[j1] = sr;
		    x[j1 + n] = si;
		    mtemp1 = abs(x[j1]) + abs(x[j1 + n]);
		    xmax = max(mtemp1, xmax);
		} else {
//2 by 2 diagonal block
//Scale if necessary to avoid overflow in forming the
//right-hand side element by inner product.
		    mtemp1 = abs(x[j1]) + abs(x[n + j1]), mtemp2 = abs(x[j2]) + abs(x[n + j2]);
		    xj = max(mtemp1, mtemp2);
		    if (xmax > One) {
			rec = One / xmax;
			mtemp1 = work[j1], mtemp2 = work[j2];
			if (max(mtemp1, mtemp2) > (bignum - xj) / xmax) {
			    Rscal(n2, rec, &x[0], 1);
			    *scale = *scale * rec;
			    xmax = xmax * rec;
			}
		    }
		    d[0] = x[j1] - Rdot(j1 - 1, &t[j1 * ldt + 1], 1, &x[0], 1);
		    d[1] = x[j2] - Rdot(j1 - 1, &t[j2 * ldt + 1], 1, &x[0], 1);
		    d[2] = x[n + j1] - Rdot(j1 - 1, &t[j1 * ldt + 1], 1, &x[n + 1], 1);
		    d[3] = x[n + j2] - Rdot(j1 - 1, &t[j2 * ldt + 1], 1, &x[n + 1], 1);
		    d[0] = d[0] - B[j1] * x[n + 1];
		    d[1] = d[1] - B[j2] * x[n + 1];
		    d[2] = d[2] + B[j1] * x[0];
		    d[3] = d[3] + B[j2] * x[0];
		    Rlaln2(MTRUE, 2, 2, sminw, Zero, &t[j1 + j1 * ldt], ldt, Zero, Zero, d, 2, One, w, v, 2, &scaloc, &xnorm, &ierr);
		    if (ierr != 0) {
			*info = 2;
		    }
		    if (scaloc != One) {
			Rscal(n2, scaloc, &x[0], 1);
			*scale = scaloc * *scale;
		    }
		    x[j1] = v[0];
		    x[j2] = v[1];
		    x[n + j1] = v[2];
		    x[n + j2] = v[3];
		    mtemp1 = abs(x[j1]) + abs(x[n + j1]), mtemp2 = abs(x[j2]) + abs(x[n + j2]);
		    mtemp3 = max(mtemp1, mtemp2);
		    xmax = max(mtemp3, xmax);
		}
	      L80:
		;
	    }
	}
    }
    return;
}
