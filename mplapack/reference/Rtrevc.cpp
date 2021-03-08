/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rtrevc.cpp,v 1.11 2010/08/07 04:48:33 nakatamaho Exp $ 
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
#define MTRUE 1

void
Rtrevc(const char *side, const char *howmny, LOGICAL * select,
       INTEGER n, REAL * t, INTEGER ldt, REAL * vl, INTEGER ldvl, REAL * vr, INTEGER ldvr, INTEGER mm, INTEGER * m, REAL * work, INTEGER * info)
{
    INTEGER i, j, k;
    REAL x[4];
    INTEGER j1, j2, n2, ii, ki, ip, is;
    REAL wi, wr, rec, ulp, beta, emax;
    INTEGER pair;
    INTEGER allv;
    INTEGER ierr;
    REAL unfl, ovfl, smin;
    INTEGER over;
    REAL vmax;
    INTEGER jnxt;
    REAL scale;
    REAL remax;
    INTEGER leftv, bothv;
    REAL vcrit;
    INTEGER somev;
    REAL xnorm;
    REAL bignum;
    INTEGER rightv;
    REAL smlnum;
    REAL mtemp1, mtemp2;
    REAL mtemp3, mtemp4;
    REAL Zero = 0.0, One = 1.0;

//Decode and test the input parameters
    bothv = Mlsame(side, "B");
    rightv = Mlsame(side, "R") || bothv;
    leftv = Mlsame(side, "L") || bothv;

    allv = Mlsame(howmny, "A");
    over = Mlsame(howmny, "B");
    somev = Mlsame(howmny, "S");

    *info = 0;
    if (!rightv && !leftv) {
	*info = -1;
    } else if (!allv && !over && !somev) {
	*info = -2;
    } else if (n < 0) {
	*info = -4;
    } else if (ldt < max((INTEGER) 1, n)) {
	*info = -6;
    } else if (ldvl < 1 || (leftv && ldvl < n)) {
	*info = -8;
    } else if (ldvr < 1 || (rightv && ldvr < n)) {
	*info = -10;
    } else {
//Set M to the number of columns required to store the selected
//eigenvectors, standardize the array SELECT if necessary, and
//test MM.
	if (somev) {
	    m = 0;
	    pair = MFALSE;
	    for (j = 0; j < n; j++) {
		if (pair) {
		    pair = MFALSE;
		    select[j] = MFALSE;
		} else {
		    if (j < n) {
			if (t[j + 1 + j * ldt] == Zero) {
			    if (select[j]) {
				++(m);
			    }
			} else {
			    pair = MTRUE;
			    if (select[j] || select[j + 1]) {
				select[j] = MTRUE;
				m += 2;
			    }
			}
		    } else {
			if (select[n]) {
			    ++(m);
			}
		    }
		}
	    }
	} else {
	    *m = n;
	}
	if (mm < *m) {
	    *info = -11;
	}
    }
    if (*info != 0) {
	Mxerbla("Rtrevc", -(*info));
	return;
    }
//Quick return if possible.
    if (n == 0) {
	return;
    }
//Set the constants to control overflow.
    unfl = Rlamch("Safe minimum");
    ovfl = Zero / unfl;
    ulp = Rlamch("Precision");
    smlnum = unfl * (n / ulp);
    bignum = (One - ulp) / smlnum;
//Compute 1-norm of each column of strictly upper triangular
//part of T to control overflow in triangular solver.
    work[1] = Zero;
    for (j = 2; j <= n; j++) {
	work[j] = Zero;
	for (i = 0; i < j - 1; i++) {
	    work[j] = work[j] + abs(t[i + j * ldt]);
	}
    }
//Index IP is used to specify the real or complex eigenvalue:
//  IP = 0, real eigenvalue,
//       1, first of conjugate complex pair: (wr,wi)
//      -1, second of conjugate complex pair: (wr,wi)
    n2 = n * 2;
    if (rightv) {
//Compute right eigenvectors.
	ip = 0;
	is = *m;
	for (ki = n; ki >= 1; ki--) {
	    if (ip == 1) {
		goto L130;
	    }
	    if (ki == 1) {
		goto L40;
	    }
	    if (t[ki + (ki - 1) * ldt] == Zero) {
		goto L40;
	    }
	    ip = -1;
	  L40:
	    if (somev) {
		if (ip == 0) {
		    if (!select[ki]) {
			goto L130;
		    }
		} else {
		    if (!select[ki - 1]) {
			goto L130;
		    }
		}
	    }
//Compute the KI-th eigenvalue (WR,WI).
	    wr = t[ki + ki * ldt];
	    wi = Zero;
	    if (ip != 0) {
		wi = sqrt(abs(t[ki + (ki - 1) * ldt])) * sqrt(abs(t[ki - 1 + ki * ldt]));
	    }
	    mtemp1 = ulp * (abs(wr) + abs(wi));
	    smin = max(mtemp1, smlnum);
	    if (ip == 0) {
//Real right eigenvector
		work[ki + n] = One;
//Form right-hand side
		for (k = 0; k < ki - 1; k++) {
		    work[k + n] = -t[k + ki * ldt];
		}
//Solve the upper quasi-triangular system:
//   (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK.
		jnxt = ki - 1;
		for (j = ki - 1; j >= 1; j--) {
		    if (j > jnxt) {
			goto L60;
		    }
		    j1 = j;
		    j2 = j;
		    jnxt = j - 1;
		    if (j > 1) {
			if (t[j + (j - 1) * ldt] != Zero) {
			    j1 = j - 1;
			    jnxt = j - 2;
			}
		    }
		    if (j1 == j2) {
//1-by-1 diagonal block
			Rlaln2(MFALSE, 1, 1, smin, One, &t[j + j * ldt], ldt, One, One, &work[j + n], n, wr, Zero, x, 2, &scale, &xnorm, &ierr);
//Scale X(1,1) to avoid overflow when updating
//the right-hand side.
			if (xnorm > One) {
			    if (work[j] > bignum / xnorm) {
				x[0] = x[0] / xnorm;
				scale = scale / xnorm;
			    }
			}
//Scale if necessary
			if (scale != One) {
			    Rscal(ki, scale, &work[n + 1], 1);
			}
			work[j + n] = x[0];
//Update right-hand side
			Raxpy(j - 1, -x[0], &t[j * ldt + 1], 1, &work[n + 1], 1);
		    } else {
//2-by-2 diagonal block
			Rlaln2(MFALSE, 2, 1, smin, One, &t[j - 1 + (j - 1) * ldt], ldt, One, One, &work[j - 1 + n], n, wr, Zero, x, 2, &scale, &xnorm, &ierr);
//Scale X(1,1) and X(2,1) to avoid overflow when
//updating the right-hand side.
			if (xnorm > One) {
			    beta = max(work[j - 1], work[j]);
			    if (beta > bignum / xnorm) {
				x[0] = x[0] / xnorm;
				x[0] = x[0] / xnorm;
				scale = scale / xnorm;
			    }
			}
//Scale if necessary
			if (scale != One) {
			    Rscal(ki, scale, &work[n + 1], 1);
			}
			work[j - 1 + n] = x[0];
			work[j + n] = x[0];
//Update right-hand side
			Raxpy(j - 2, -x[0], &t[(j - 1) * ldt + 1], 1, &work[n + 1], 1);
			Raxpy(j - 2, -x[0], &t[j * ldt + 1], 1, &work[n + 1], 1);
		    }
		  L60:
		    ;
		}
//Copy the vector x or Q*x to VR and normalize.
		if (!over) {
		    Rcopy(ki, &work[n + 1], 1, &vr[is * ldvr + 1], 1);
		    ii = iRamax(ki, &vr[is * ldvr + 1], 1);
		    remax = One / abs(vr[ii + is * ldvr]);
		    Rscal(ki, remax, &vr[is * ldvr + 1], 1);
		    for (k = ki + 1; k <= n; k++) {
			vr[k + is * ldvr] = Zero;
		    }
		} else {
		    if (ki > 1) {
			Rgemv("N", n, ki - 1, One, vr, ldvr, &work[n + 1], 1, work[ki + n], &vr[ki * ldvr + 1], 1);
		    }
		    ii = iRamax(n, &vr[ki * ldvr + 1], 1);
		    remax = One / abs(vr[ii + ki * ldvr]);
		    Rscal(n, remax, &vr[ki * ldvr + 1], 1);
		}
	    } else {
//Complex right eigenvector.
//Initial solve
//  [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = Zero
//  [ (T(KI,KI-1)   T(KI,KI)   )               ]
		if (abs(t[ki - 1 + ki * ldt]) >= abs(t[ki + (ki - 1) * ldt])) {
		    work[ki - 1 + n] = One;
		    work[ki + n2] = wi / t[ki - 1 + ki * ldt];
		} else {
		    work[ki - 1 + n] = -wi / t[ki + (ki - 1) * ldt];
		    work[ki + n2] = One;
		}
		work[ki + n] = Zero;
		work[ki - 1 + n2] = Zero;
//Form right-hand side
		for (k = 0; k < ki - 2; k++) {
		    work[k + n] = -work[ki - 1 + n] * t[k + (ki - 1) * ldt];
		    work[k + n2] = -work[ki + n2] * t[k + ki * ldt];
		}
//Solve upper quasi-triangular system:
//(T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2)
		jnxt = ki - 2;
		for (j = ki - 2; j >= 1; j--) {
		    if (j > jnxt) {
			goto L90;
		    }
		    j1 = j;
		    j2 = j;
		    jnxt = j - 1;
		    if (j > 1) {
			if (t[j + (j - 1) * ldt] != Zero) {
			    j1 = j - 1;
			    jnxt = j - 2;
			}
		    }
		    if (j1 == j2) {
//1-by-1 diagonal block
			Rlaln2(MFALSE, 1, 2, smin, One, &t[j + j * ldt], ldt, One, One, &work[j + n], n, wr, wi, x, 2, &scale, &xnorm, &ierr);
//Scale X(1,1) and X(1,2) to avoid overflow when
//updating the right-hand side.
			if (xnorm > One) {
			    if (work[j] > bignum / xnorm) {
				x[0] = x[0] / xnorm;
				x[2] = x[2] / xnorm;
				scale = scale / xnorm;
			    }
			}
//Scale if necessary
			if (scale != One) {
			    Rscal(ki, scale, &work[n + 1], 1);
			    Rscal(ki, scale, &work[n2 + 1], 1);
			}
			work[j + n] = x[0];
			work[j + n2] = x[2];
//Update the right-hand side
			Raxpy(j - 1, -x[0], &t[j * ldt + 1], 1, &work[n + 1], 1);
			Raxpy(j - 1, -x[2], &t[j * ldt + 1], 1, &work[n2 + 1], 1);
		    } else {
//2-by-2 diagonal block
			Rlaln2(MFALSE, 2, 2, smin, One, &t[j - 1 + (j - 1) * ldt], ldt, One, One, &work[j - 1 + n], n, wr, wi, x, 2, &scale, &xnorm, &ierr);
//Scale X to avoid overflow when updating
//the right-hand side.
			if (xnorm > One) {
			    beta = max(work[j - 1], work[j]);
			    if (beta > bignum / xnorm) {
				rec = One / xnorm;
				x[0] = x[0] * rec;
				x[1] = x[1] * rec;
				x[2] = x[2] * rec;
				x[3] = x[3] * rec;
				scale = scale * rec;
			    }
			}
//Scale if necessary
			if (scale != One) {
			    Rscal(ki, scale, &work[n + 1], 1);
			    Rscal(ki, scale, &work[n2 + 1], 1);
			}
			work[j - 1 + n] = x[0];
			work[j + n] = x[0];
			work[j - 1 + n2] = x[2];
			work[j + n2] = x[3];
//Update the right-hand side
			Raxpy(j - 2, -x[0], &t[(j - 1) * ldt + 1], 1, &work[n + 1], 1);
			Raxpy(j - 2, -x[0], &t[j * ldt + 1], 1, &work[n + 1], 1);
			Raxpy(j - 2, -x[2], &t[(j - 1) * ldt + 1], 1, &work[n2 + 1], 1);
			Raxpy(j - 2, -x[3], &t[j * ldt + 1], 1, &work[n2 + 1], 1);
		    }
		  L90:
		    ;
		}
//Copy the vector x or Q*x to VR and normalize.
		if (!over) {
		    Rcopy(ki, &work[n + 1], 1, &vr[(is - 1) * ldvr + 1], 1);
		    Rcopy(ki, &work[n2 + 1], 1, &vr[is * ldvr + 1], 1);
		    emax = Zero;
		    for (k = 0; k < ki; k++) {
			mtemp1 = emax, mtemp2 = abs(vr[k + (is - 1) * ldvr]) + abs(vr[k + is * ldvr]);
			emax = max(mtemp1, mtemp2);
		    }
		    remax = One / emax;
		    Rscal(ki, remax, &vr[(is - 1) * ldvr + 1], 1);
		    Rscal(ki, remax, &vr[is * ldvr + 1], 1);
		    for (k = ki + 1; k <= n; k++) {
			vr[k + (is - 1) * ldvr] = Zero;
			vr[k + is * ldvr] = Zero;
		    }
		} else {
		    if (ki > 2) {
			Rgemv("N", n, ki - 2, One, vr, ldvr, &work[n + 1], 1, work[ki - 1 + n], &vr[(ki - 1) * ldvr + 1], 1);
			Rgemv("N", n, ki - 2, One, vr, ldvr, &work[n2 + 1], 1, work[ki + n2], &vr[ki * ldvr + 1], 1);
		    } else {
			Rscal(n, work[ki - 1 + n], &vr[(ki - 1) * ldvr + 1], 1);
			Rscal(n, work[ki + n2], &vr[ki * ldvr + 1], 1);
		    }
		    emax = Zero;
		    for (k = 0; k < n; k++) {
			mtemp1 = emax, mtemp2 = abs(vr[k + (ki - 1) * ldvr]) + abs(vr[k + ki * ldvr]);
			emax = max(mtemp1, mtemp2);
		    }
		    remax = One / emax;
		    Rscal(n, remax, &vr[(ki - 1) * ldvr + 1], 1);
		    Rscal(n, remax, &vr[ki * ldvr + 1], 1);
		}
	    }
	    is--;
	    if (ip != 0) {
		is--;
	    }
	  L130:
	    if (ip == 1) {
		ip = 0;
	    }
	    if (ip == -1) {
		ip = 1;
	    }
	}
    }
    if (leftv) {
//Compute left eigenvectors.
	ip = 0;
	is = 1;
	for (ki = 0; ki <= n; ki++) {
	    if (ip == -1) {
		goto L250;
	    }
	    if (ki == n) {
		goto L150;
	    }
	    if (t[ki + 1 + ki * ldt] == Zero) {
		goto L150;
	    }
	    ip = 1;
	  L150:
	    if (somev) {
		if (!select[ki]) {
		    goto L250;
		}
	    }
//Compute the KI-th eigenvalue (WR,WI).
	    wr = t[ki + ki * ldt];
	    wi = Zero;
	    if (ip != 0) {
		wi = sqrt(abs(t[ki + (ki + 1) * ldt])) * sqrt(abs(t[ki + 1 + ki * ldt]));
	    }
	    mtemp1 = ulp * (abs(wr) + abs(wi));
	    smin = max(mtemp1, smlnum);
	    if (ip == 0) {
//Real left eigenvector.
		work[ki + n] = One;
//Form right-hand side
		for (k = ki + 1; k <= n; k++) {
		    work[k + n] = -t[ki + k * ldt];
		}
//Solve the quasi-triangular system:
//  (T(KI+1:N,KI+1:N) - WR)'*X = SCALE*WORK
		vmax = One;
		vcrit = bignum;
		jnxt = ki + 1;
		for (j = ki + 1; j <= n; j++) {
		    if (j < jnxt) {
			goto L170;
		    }
		    j1 = j;
		    j2 = j;
		    jnxt = j + 1;
		    if (j < n) {
			if (t[j + 1 + j * ldt] != Zero) {
			    j2 = j + 1;
			    jnxt = j + 2;
			}
		    }
		    if (j1 == j2) {
//1-by-1 diagonal block
//Scale if necessary to avoid overflow when forming
//the right-hand side.
			if (work[j] > vcrit) {
			    rec = One / vmax;
			    Rscal(n - ki + 1, rec, &work[ki + n], 1);
			    vmax = One;
			    vcrit = bignum;
			}
			work[j + n] = work[j + n] - Rdot(j - ki - 1, &t[ki + 1 + j * ldt], 1, &work[ki + 1 + n], 1);
//Solve (T(J,J)-WR)'*X = WORK
			Rlaln2(MFALSE, 1, 1, smin, One, &t[j + j * ldt], ldt, One, One, &work[j + n], n, wr, Zero, x, 2, &scale, &xnorm, &ierr);
//Scale if necessary
			if (scale != One) {
			    Rscal(n - ki + 1, scale, &work[ki + n], 1);
			}
			work[j + n] = x[0];
			mtemp1 = abs(work[j + n]);
			vmax = max(mtemp1, vmax);
			vcrit = bignum / vmax;
		    } else {
//2-by-2 diagonal block
//Scale if necessary to avoid overflow when forming
//the right-hand side.
			beta = max(work[j], work[j + 1]);
			if (beta > vcrit) {
			    rec = One / vmax;
			    Rscal(n - ki + 1, rec, &work[ki + n], 1);
			    vmax = One;
			    vcrit = bignum;
			}
			work[j + n] = work[j + n] - Rdot(j - ki - 1, &t[ki + 1 + j * ldt], 1, &work[ki + 1 + n], 1);
			work[j + 1 + n] = work[j + 1 + n] - Rdot(j - ki - 1, &t[ki + 1 + (j + 1) * ldt], 1, &work[ki + 1 + n], 1);
//Solve
//  [T(J,J)-WR   T(J,J+1)     ]'* X = SCALE*( WORK1 )
//  [T(J+1,J)    T(J+1,J+1)-WR]             ( WORK2 )
			Rlaln2(MTRUE, 2, 1, smin, One, &t[j + j * ldt], ldt, One, One, &work[j + n], n, wr, Zero, x, 2, &scale, &xnorm, &ierr);
//Scale if necessary
			if (scale != One) {
			    Rscal(n - ki + 1, scale, &work[ki + n], 1);
			}
			work[j + n] = x[0];
			work[j + 1 + n] = x[0];

			mtemp1 = abs(work[j + n]), mtemp2 = abs(work[j + 1 + n]);
			mtemp3 = max(mtemp1, mtemp2);
			vmax = max(mtemp3, vmax);
			vcrit = bignum / vmax;
		    }
		  L170:
		    ;
		}
//Copy the vector x or Q*x to VL and normalize.
		if (!over) {
		    Rcopy(n - ki + 1, &work[ki + n], 1, &vl[ki + is * ldvl], 1);
		    ii = iRamax(n - ki + 1, &vl[ki + is * ldvl], 1) + ki - 1;
		    remax = One / abs(vl[ii + is * ldvl]);
		    Rscal(n - ki + 1, remax, &vl[ki + is * ldvl], 1);
		    for (k = 0; k < ki - 1; k++) {
			vl[k + is * ldvl] = Zero;
		    }
		} else {
		    if (ki < n) {
			Rgemv("N", n, n - ki, One, &vl[(ki + 1) * ldvl + 1], ldvl, &work[ki + 1 + n], 1, work[ki + n], &vl[ki * ldvl + 1], 1);
		    }
		    ii = iRamax(n, &vl[ki * ldvl + 1], 1);
		    remax = One / abs(vl[ii + ki * ldvl]);
		    Rscal(n, remax, &vl[ki * ldvl + 1], 1);
		}
	    } else {
//Complex left eigenvector.
// Initial solve:
//   ((T(KI,KI)    T(KI,KI+1) )' - (WR - I* WI))*X = Zero
//   ((T(KI+1,KI) T(KI+1,KI+1))                )
		if (abs(t[ki + (ki + 1) * ldt]) >= abs(t[ki + 1 + ki * ldt])) {
		    work[ki + n] = wi / t[ki + (ki + 1) * ldt];
		    work[ki + 1 + n2] = One;
		} else {
		    work[ki + n] = One;
		    work[ki + 1 + n2] = -wi / t[ki + 1 + ki * ldt];
		}
		work[ki + 1 + n] = Zero;
		work[ki + n2] = Zero;
//Form right-hand side
		for (k = ki + 2; k <= n; k++) {
		    work[k + n] = -work[ki + n] * t[ki + k * ldt];
		    work[k + n2] = -work[ki + 1 + n2] * t[ki + 1 + k * ldt];
		}
//Solve complex quasi-triangular system:
//( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2
		vmax = One;
		vcrit = bignum;
		jnxt = ki + 2;
		for (j = ki + 2; j <= n; j++) {
		    if (j < jnxt) {
			goto L200;
		    }
		    j1 = j;
		    j2 = j;
		    jnxt = j + 1;
		    if (j < n) {
			if (t[j + 1 + j * ldt] != Zero) {
			    j2 = j + 1;
			    jnxt = j + 2;
			}
		    }
		    if (j1 == j2) {
//1-by-1 diagonal block
//Scale if necessary to avoid overflow when
//forming the right-hand side elements.
			if (work[j] > vcrit) {
			    rec = One / vmax;
			    Rscal(n - ki + 1, rec, &work[ki + n], 1);
			    Rscal(n - ki + 1, rec, &work[ki + n2], 1);
			    vmax = One;
			    vcrit = bignum;
			}
			work[j + n] = work[j + n] - Rdot(j - ki - 2, &t[ki + 2 + j * ldt], 1, &work[ki + 2 + n], 1);
			work[j + n2] = work[j + n2] - Rdot(j - ki - 2, &t[ki + 2 + j * ldt], 1, &work[ki + 2 + n2], 1);
//Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2
			Rlaln2(MFALSE, 1, 2, smin, One, &t[j + j * ldt], ldt, One, One, &work[j + n], n, wr, -wi, x, 2, &scale, &xnorm, &ierr);
//Scale if necessary
			if (scale != One) {
			    Rscal(n - ki + 1, scale, &work[ki + n], 1);
			    Rscal(n - ki + 1, scale, &work[ki + n2], 1);
			}
			work[j + n] = x[0];
			work[j + n2] = x[2];
			mtemp1 = abs(work[j + n]), mtemp2 = abs(work[j + n2]);
			mtemp3 = max(mtemp1, mtemp2);
			vmax = max(mtemp3, vmax);
			vcrit = bignum / vmax;
		    } else {
//2-by-2 diagonal block
//Scale if necessary to avoid overflow when forming
//the right-hand side elements.
			beta = max(work[j], work[j + 1]);
			if (beta > vcrit) {
			    rec = One / vmax;
			    Rscal(n - ki + 1, rec, &work[ki + n], 1);
			    Rscal(n - ki + 1, rec, &work[ki + n2], 1);
			    vmax = One;
			    vcrit = bignum;
			}
			work[j + n] = work[j + n] - Rdot(j - ki - 2, &t[ki + 2 + j * ldt], 1, &work[ki + 2 + n], 1);
			work[j + n2] = work[j + n2] - Rdot(j - ki - 2, &t[ki + 2 + j * ldt], 1, &work[ki + 2 + n2], 1);
			work[j + 1 + n] = work[j + 1 + n] - Rdot(j - ki - 2, &t[ki + 2 + (j + 1) * ldt], 1, &work[ki + 2 + n], 1);
			work[j + 1 + n2] = work[j + 1 + n2] - Rdot(j - ki - 2, &t[ki + 2 + (j + 1) * ldt], 1, &work[ki + 2 + n2], 1);
//Solve 2-by-2 complex linear equation
//  ([T(j,j)   T(j,j+1)  ]'-(wr-i*wi)*I)*X = SCALE*B
//  ([T(j+1,j) T(j+1,j+1)]             )
			Rlaln2(MTRUE, 2, 2, smin, One, &t[j + j * ldt], ldt, One, One, &work[j + n], n, wr, -wi, x, 2, &scale, &xnorm, &ierr);
//Scale if necessary
			if (scale != One) {
			    Rscal(n - ki + 1, scale, &work[ki + n], 1);
			    Rscal(n - ki + 1, scale, &work[ki + n2], 1);
			}
			work[j + n] = x[0];
			work[j + n2] = x[1];
			work[j + 1 + n] = x[2];
			work[j + 1 + n2] = x[3];
			mtemp3 = abs(x[0]), mtemp4 = abs(x[1]);
			mtemp1 = max(mtemp3, mtemp4);
			mtemp3 = abs(x[2]), mtemp4 = abs(x[3]);
			mtemp2 = max(mtemp3, mtemp4);
			mtemp3 = max(mtemp1, mtemp2);
			vmax = max(mtemp3, vmax);
			vcrit = bignum / vmax;
		    }
		  L200:
		    ;
		}
//Copy the vector x or Q*x to VL and normalize.
		if (!over) {
		    Rcopy(n - ki + 1, &work[ki + n], 1, &vl[ki + is * ldvl], 1);
		    Rcopy(n - ki + 1, &work[ki + n2], 1, &vl[ki + (is + 1) * ldvl], 1);
		    emax = Zero;
		    for (k = ki; k <= n; k++) {
			mtemp1 = emax, mtemp2 = abs(vl[k + is * ldvl]) + abs(vl[k + (is + 1) * ldvl]);
			emax = max(mtemp1, mtemp2);
		    }
		    remax = One / emax;
		    Rscal(n - ki + 1, remax, &vl[ki + is * ldvl], 1);
		    Rscal(n - ki + 1, remax, &vl[ki + (is + 1) * ldvl], 1);
		    for (k = 0; k < ki - 1; k++) {
			vl[k + is * ldvl] = Zero;
			vl[k + (is + 1) * ldvl] = Zero;
		    }
		} else {
		    if (ki < n - 1) {
			Rgemv("N", n, n - ki - 1, One, &vl[(ki + 2) * ldvl + 1], ldvl, &work[ki + 2 + n], 1, work[ki + n], &vl[ki * ldvl + 1], 1);
			Rgemv("N", n, n - ki - 1, One, &vl[(ki + 2) * ldvl + 1], ldvl, &work[ki + 2 + n2], 1, work[ki + 1 + n2], &vl[(ki + 1) * ldvl + 1], 1);
		    } else {
			Rscal(n, work[ki + n], &vl[ki * ldvl + 1], 1);
			Rscal(n, work[ki + 1 + n2], &vl[(ki + 1) * ldvl + 1], 1);
		    }
		    emax = Zero;
		    for (k = 0; k < n; k++) {
			mtemp1 = emax, mtemp2 = abs(vl[k + ki * ldvl]) + abs(vl[k + (ki + 1) * ldvl]);
			emax = max(mtemp1, mtemp2);
		    }
		    remax = One / emax;
		    Rscal(n, remax, &vl[ki * ldvl + 1], 1);
		    Rscal(n, remax, &vl[(ki + 1) * ldvl + 1], 1);
		}
	    }
	    is++;
	    if (ip != 0) {
		is++;
	    }
	  L250:
	    if (ip == -1) {
		ip = 0;
	    }
	    if (ip == 1) {
		ip = -1;
	    }
	}
    }
    return;
}
