/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasd4.cpp,v 1.7 2010/08/07 04:48:33 nakatamaho Exp $ 
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

#define MTRUE 1
#define MFALSE 0

void Rlasd4(INTEGER n, INTEGER i, REAL * d, REAL * z, REAL * delta, REAL rho, REAL * sigma, REAL * work, INTEGER * info)
{
    REAL a, b, c;
    INTEGER j;
    REAL w, dd[3];
    INTEGER ii;
    REAL dw, zz[3];
    INTEGER ip1;
    REAL eta, phi, eps, tau, psi;
    INTEGER iim1, iip1;
    REAL dphi, dpsi;
    INTEGER iter;
    REAL temp, prew, sg2lb, sg2ub, temp1, temp2, dtiim, delsq, dtiip;
    INTEGER niter;
    REAL dtisq;
    INTEGER swtch;
    REAL dtnsq;
    REAL delsq2, dtnsq1;
    INTEGER swtch3;
    INTEGER orgati;
    REAL erretm, dtipsq, rhoinv;
    REAL Zero = 0.0, One = 1.0, Two = 2.0, Three = 3.0, Four = 4.0, Eight = 8.0, Ten = 10.0;

//Since this routine is called in an inner loop, we do no argument
//checking.
//Quick return for N=1 and Two

    *info = 0;
    if (n == 1) {
//Presumably, I=1 upon entry
	*sigma = sqrt(d[1] * d[1] + rho * z[1] * z[1]);
	delta[1] = One;
	work[1] = One;
	return;
    }
    if (n == 2) {
	Rlasd5(i, &d[0], &z[1], &delta[1], rho, sigma, &work[0]);
	return;
    }
//Compute machine epsilon
    eps = Rlamch("Epsilon");
    rhoinv = One / rho;

//The case I = N
    if (i == n) {
//Initialize some basic variables
	ii = n - 1;
	niter = 1;
//Calculate initial guess
	temp = rho / Two;

//If ||Z||_2 is not one, then TEMP should be set to
//RHO * ||Z||_2^2 / TWO
	temp1 = temp / (d[n] + sqrt(d[n] * d[n] + temp));
	for (j = 0; j < n; j++) {
	    work[j] = d[j] + d[n] + temp1;
	    delta[j] = d[j] - d[n] - temp1;
	}
	psi = Zero;
	for (j = 0; j < n - 2; j++) {
	    psi += z[j] * z[j] / (delta[j] * work[j]);
	}
	c = rhoinv + psi;
	w = c + z[ii] * z[ii] / (delta[ii] * work[ii]) + z[n] * z[n] / (delta[n] * work[n]);
	if (w <= Zero) {
	    temp1 = sqrt(d[n] * d[n] + rho);
	    temp = z[n - 1] * z[n - 1] / ((d[n - 1] + temp1) * (d[n] - d[n - 1] + rho / (d[n] + temp1))) + z[n] * z[n] / rho;
//The following TAU is to approximate
//SIGMA_n^2 - D( N )*D( N )
	    if (c <= temp) {
		tau = rho;
	    } else {
		delsq = (d[n] - d[n - 1]) * (d[n] + d[n - 1]);
		a = -c * delsq + z[n - 1] * z[n - 1] + z[n] * z[n];
		b = z[n] * z[n] * delsq;
		if (a < Zero) {
		    tau = b * Two / (sqrt(a * a + b * Four * c) - a);
		} else {
		    tau = (a + sqrt(a * a + b * Four * c)) / (c * Two);
		}
	    }
//It can be proved that
//D(N)^2+RHO/2 <= SIGMA_n^2 < D(N)^2+TAU <= D(N)^2+RHO
	} else {
	    delsq = (d[n] - d[n - 1]) * (d[n] + d[n - 1]);
	    a = -c * delsq + z[n - 1] * z[n - 1] + z[n] * z[n];
	    b = z[n] * z[n] * delsq;
//The following TAU is to approximate
//SIGMA_n^2 - D( N )*D( N )
	    if (a < Zero) {
		tau = b * Two / (sqrt(a * a + b * Four * c) - a);
	    } else {
		tau = (a + sqrt(a * a + b * Four * c)) / (c * Two);
	    }
//It can be proved that
//D(N)^2 < D(N)^2+TAU < SIGMA(N)^2 < D(N)^2+RHO/2
	}
//The following ETA is to approximate SIGMA_n - D( N )
	eta = tau / (d[n] + sqrt(d[n] * d[n] + tau));
	*sigma = d[n] + eta;
	for (j = 0; j < n; j++) {
	    delta[j] = d[j] - d[i] - eta;
	    work[j] = d[j] + d[i] + eta;
	}
//Evaluate PSI and the derivative DPSI
	dpsi = Zero;
	psi = Zero;
	erretm = Zero;
	for (j = 0; j < ii; j++) {
	    temp = z[j] / (delta[j] * work[j]);
	    psi += z[j] * temp;
	    dpsi += temp * temp;
	    erretm += psi;
	}
	erretm = abs(erretm);
//Evaluate PHI and the derivative DPHI
	temp = z[n] / (delta[n] * work[n]);
	phi = z[n] * temp;
	dphi = temp * temp;
	erretm = (-phi - psi) * Eight + erretm - phi + rhoinv + abs(tau) * (dpsi + dphi);
	w = rhoinv + phi + psi;

//Test for convergence
	if (abs(w) <= eps * erretm) {
	    goto L240;
	}
//Calculate the new step
	++niter;
	dtnsq1 = work[n - 1] * delta[n - 1];
	dtnsq = work[n] * delta[n];
	c = w - dtnsq1 * dpsi - dtnsq * dphi;
	a = (dtnsq + dtnsq1) * w - dtnsq * dtnsq1 * (dpsi + dphi);
	b = dtnsq * dtnsq1 * w;
	if (c < Zero) {
	    c = abs(c);
	}
	if (c == Zero) {
	    eta = rho - *sigma * *sigma;
	} else if (a >= Zero) {
	    eta = (a + sqrt(abs(a * a - b * Four * c))) / (c * Two);
	} else {
	    eta = b * Two / (a - sqrt(abs(a * a - b * Four * c)));
	}

//Note, eta should be positive if w is negative, and
//eta should be negative otherwise. However,
//if for some reason caused by roundoff, eta*w > 0,
//we simply use one Newton step instead. This way
//will guarantee eta*w < Zero
	if (w * eta > Zero) {
	    eta = -w / (dpsi + dphi);
	}
	temp = eta - dtnsq;
	if (temp > rho) {
	    eta = rho + dtnsq;
	}
	tau += eta;
	eta /= *sigma + sqrt(eta + *sigma * *sigma);
	for (j = 0; j < n; j++) {
	    delta[j] -= eta;
	    work[j] += eta;
	}
	*sigma += eta;
//Evaluate PSI and the derivative DPSI
	dpsi = Zero;
	psi = Zero;
	erretm = Zero;
	for (j = 0; j < ii; j++) {
	    temp = z[j] / (work[j] * delta[j]);
	    psi += z[j] * temp;
	    dpsi += temp * temp;
	    erretm += psi;

	}
	erretm = abs(erretm);
//Evaluate PHI and the derivative DPHI
	temp = z[n] / (work[n] * delta[n]);
	phi = z[n] * temp;
	dphi = temp * temp;
	erretm = (-phi - psi) * Eight + erretm - phi + rhoinv + abs(tau) * (dpsi + dphi);
	w = rhoinv + phi + psi;

//Main loop to update the values of the array   DELTA
	iter = niter + 1;
	for (niter = iter; niter <= 20; niter++) {

//Test for convergenc
	    if (abs(w) <= eps * erretm) {
		goto L240;
	    }
//Calculate the new step
	    dtnsq1 = work[n - 1] * delta[n - 1];
	    dtnsq = work[n] * delta[n];
	    c = w - dtnsq1 * dpsi - dtnsq * dphi;
	    a = (dtnsq + dtnsq1) * w - dtnsq1 * dtnsq * (dpsi + dphi);
	    b = dtnsq1 * dtnsq * w;
	    if (a >= Zero) {
		eta = (a + sqrt(abs(a * a - b * Four * c))) / (c * Two);
	    } else {
		eta = b * Two / (a - sqrt(abs(a * a - b * Four * c)));
	    }
//Note, eta should be positive if w is negative, and
//eta should be negative otherwise. However,
//if for some reason caused by roundoff, eta*w > 0,
//we simply use one Newton step instead. This way
//will guarantee eta*w < Zero
	    if (w * eta > Zero) {
		eta = -w / (dpsi + dphi);
	    }
	    temp = eta - dtnsq;
	    if (temp <= Zero) {
		eta /= Two;
	    }
	    tau += eta;
	    eta /= *sigma + sqrt(eta + *sigma * *sigma);
	    for (j = 0; j < n; j++) {
		delta[j] -= eta;
		work[j] += eta;
	    }
	    *sigma += eta;
//Evaluate PSI and the derivative DPSI
	    dpsi = Zero;
	    psi = Zero;
	    erretm = Zero;
	    for (j = 0; j < ii; j++) {
		temp = z[j] / (work[j] * delta[j]);
		psi += z[j] * temp;
		dpsi += temp * temp;
		erretm += psi;
	    }
	    erretm = abs(erretm);
//Evaluate PHI and the derivative DPHI
	    temp = z[n] / (work[n] * delta[n]);
	    phi = z[n] * temp;
	    dphi = temp * temp;
	    erretm = (-phi - psi) * Eight + erretm - phi + rhoinv + abs(tau) * (dpsi + dphi);
	    w = rhoinv + phi + psi;
	}
//Return with INFO = 1, NITER = MAXIT and not converged
	*info = 1;
	goto L240;
//End for the case I = N
    } else {
//The case for I < N
	niter = 1;
	ip1 = i + 1;
//Calculate initial guess
	delsq = (d[ip1] - d[i]) * (d[ip1] + d[i]);
	delsq2 = delsq / Two;
	temp = delsq2 / (d[i] + sqrt(d[i] * d[i] + delsq2));
	for (j = 0; j < n; j++) {
	    work[j] = d[j] + d[i] + temp;
	    delta[j] = d[j] - d[i] - temp;
	}
	psi = Zero;
	for (j = 0; j < i - 1; j++) {
	    psi += z[j] * z[j] / (work[j] * delta[j]);
	}

	phi = Zero;
	for (j = n; j >= i + 2; j--) {
	    phi += z[j] * z[j] / (work[j] * delta[j]);

	}
	c = rhoinv + psi + phi;
	w = c + z[i] * z[i] / (work[i] * delta[i]) + z[ip1] * z[ip1] / (work[ip1] * delta[ip1]);
	if (w > Zero) {
//d(i)^2 < the ith sigma^2 < (d(i)^2+d(i+1)^2)/2
//We choose d(i) as origin.
	    orgati = MTRUE;
	    sg2lb = Zero;
	    sg2ub = delsq2;
	    a = c * delsq + z[i] * z[i] + z[ip1] * z[ip1];
	    b = z[i] * z[i] * delsq;
	    if (a > Zero) {
		tau = b * Two / (a + sqrt(abs(a * a - b * Four * c)));
	    } else {
		tau = (a - sqrt(abs(a * a - b * Four * c))) / (c * Two);
	    }
//TAU now is an estimation of SIGMA^2 - D( I )^Two The
//following, however, is the corresponding estimation of
//SIGMA - D( I ).
	    eta = tau / (d[i] + sqrt(d[i] * d[i] + tau));
	} else {
//(d(i)^2+d(i+1)^2)/2 <= the ith sigma^2 < d(i+1)^2/2
//We choose d(i+1) as origin.
	    orgati = MFALSE;
	    sg2lb = -delsq2;
	    sg2ub = Zero;
	    a = c * delsq - z[i] * z[i] - z[ip1] * z[ip1];
	    b = z[ip1] * z[ip1] * delsq;
	    if (a < Zero) {
		tau = b * Two / (a - sqrt(abs(a * a + b * Four * c)));
	    } else {
		tau = -(a + sqrt(abs(a * a + b * Four * c))) / (c * Two);
	    }
//TAU now is an estimation of SIGMA^2 - D( IP1 )^Two The
//following, however, is the corresponding estimation of
//SIGMA - D( IP1 ).
	    eta = tau / (d[ip1] + sqrt(abs(d[ip1] * d[ip1] + tau)));
	}

	if (orgati) {
	    ii = i;
	    *sigma = d[i] + eta;
	    for (j = 0; j < n; j++) {
		work[j] = d[j] + d[i] + eta;
		delta[j] = d[j] - d[i] - eta;
	    }
	} else {
	    ii = i + 1;
	    *sigma = d[ip1] + eta;
	    for (j = 0; j < n; j++) {
		work[j] = d[j] + d[ip1] + eta;
		delta[j] = d[j] - d[ip1] - eta;
	    }
	}
	iim1 = ii - 1;
	iip1 = ii + 1;
//Evaluate PSI and the derivative DPSI
	dpsi = Zero;
	psi = Zero;
	erretm = Zero;
	for (j = 0; j < iim1; j++) {
	    temp = z[j] / (work[j] * delta[j]);
	    psi += z[j] * temp;
	    dpsi += temp * temp;
	    erretm += psi;

	}
	erretm = abs(erretm);
//Evaluate PHI and the derivative DPHI
	dphi = Zero;
	phi = Zero;
	for (j = n; j >= iip1; j--) {
	    temp = z[j] / (work[j] * delta[j]);
	    phi += z[j] * temp;
	    dphi += temp * temp;
	    erretm += phi;
	}
	w = rhoinv + phi + psi;
//W is the value of the secular function with
//its ii-th element removed.
	swtch3 = MFALSE;
	if (orgati) {
	    if (w < Zero) {
		swtch3 = MTRUE;
	    }
	} else {
	    if (w > Zero) {
		swtch3 = MTRUE;
	    }
	}
	if (ii == 1 || ii == n) {
	    swtch3 = MFALSE;
	}
	temp = z[ii] / (work[ii] * delta[ii]);
	dw = dpsi + dphi + temp * temp;
	temp = z[ii] * temp;
	w += temp;
	erretm = (phi - psi) * Eight + erretm + rhoinv * Two + abs(temp) * 3. + abs(tau) * dw;
//Test for convergence
	if (abs(w) <= eps * erretm) {
	    goto L240;
	}
	if (w <= Zero) {
	    sg2lb = max(sg2lb, tau);
	} else {
	    sg2ub = min(sg2ub, tau);
	}
//Calculate the new step
	++niter;
	if (!swtch3) {
	    dtipsq = work[ip1] * delta[ip1];
	    dtisq = work[i] * delta[i];
	    if (orgati) {
		c = w - dtipsq * dw + delsq * (z[i] / dtisq * z[i] / dtisq);
	    } else {
		c = w - dtisq * dw - delsq * (z[ip1] / dtipsq * z[ip1] / dtipsq);
	    }
	    a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
	    b = dtipsq * dtisq * w;
	    if (c == Zero) {
		if (a == Zero) {
		    if (orgati) {
			a = z[i] * z[i] + dtipsq * dtipsq * (dpsi + dphi);
		    } else {
			a = z[ip1] * z[ip1] + dtisq * dtisq * (dpsi + dphi);
		    }
		}
		eta = b / a;
	    } else if (a <= Zero) {
		eta = (a - sqrt(abs(a * a - b * Four * c))) / (c * Two);
	    } else {
		eta = b * Two / (a + sqrt(abs(a * a - b * Four * c)));
	    }
	} else {
//Interpolation using THREE most relevant poles
	    dtiim = work[iim1] * delta[iim1];
	    dtiip = work[iip1] * delta[iip1];
	    temp = rhoinv + psi + phi;
	    if (orgati) {
		temp1 = z[iim1] / dtiim;
		temp1 *= temp1;
		c = temp - dtiip * (dpsi + dphi) - (d[iim1] - d[iip1]) * (d[iim1] + d[iip1]) * temp1;
		zz[0] = z[iim1] * z[iim1];
		if (dpsi < temp1) {
		    zz[2] = dtiip * dtiip * dphi;
		} else {
		    zz[2] = dtiip * dtiip * (dpsi - temp1 + dphi);
		}
	    } else {
		temp1 = z[iip1] / dtiip;
		temp1 *= temp1;
		c = temp - dtiim * (dpsi + dphi) - (d[iip1] - d[iim1]) * (d[iim1] + d[iip1]) * temp1;
		if (dphi < temp1) {
		    zz[0] = dtiim * dtiim * dpsi;
		} else {
		    zz[0] = dtiim * dtiim * (dpsi + (dphi - temp1));
		}
		zz[2] = z[iip1] * z[iip1];
	    }
	    zz[1] = z[ii] * z[ii];
	    dd[0] = dtiim;
	    dd[1] = delta[ii] * work[ii];
	    dd[2] = dtiip;
	    Rlaed6(niter, orgati, c, dd, zz, &w, &eta, info);
	    if (*info != 0) {
		goto L240;
	    }
	}
//Note, eta should be positive if w is negative, and
//eta should be negative otherwise. However,
//if for some reason caused by roundoff, eta*w > 0,
//we simply use one Newton step instead. This way
//will guarantee eta*w < Zero
	if (w * eta >= Zero) {
	    eta = -w / dw;
	}
	if (orgati) {
	    temp1 = work[i] * delta[i];
	    temp = eta - temp1;
	} else {
	    temp1 = work[ip1] * delta[ip1];
	    temp = eta - temp1;
	}
	if (temp > sg2ub || temp < sg2lb) {
	    if (w < Zero) {
		eta = (sg2ub - tau) / Two;
	    } else {
		eta = (sg2lb - tau) / Two;
	    }
	}
	tau += eta;
	eta /= *sigma + sqrt(*sigma * *sigma + eta);
	prew = w;
	*sigma += eta;
	for (j = 0; j < n; j++) {
	    work[j] += eta;
	    delta[j] -= eta;
	}
//Evaluate PSI and the derivative DPSI
	dpsi = Zero;
	psi = Zero;
	erretm = Zero;
	for (j = 0; j < iim1; j++) {
	    temp = z[j] / (work[j] * delta[j]);
	    psi += z[j] * temp;
	    dpsi += temp * temp;
	    erretm += psi;
	}
	erretm = abs(erretm);
//Evaluate PHI and the derivative DPHI
	dphi = Zero;
	phi = Zero;
	for (j = n; j >= iip1; j--) {
	    temp = z[j] / (work[j] * delta[j]);
	    phi += z[j] * temp;
	    dphi += temp * temp;
	    erretm += phi;
	}
	temp = z[ii] / (work[ii] * delta[ii]);
	dw = dpsi + dphi + temp * temp;
	temp = z[ii] * temp;
	w = rhoinv + phi + psi + temp;
	erretm = (phi - psi) * Eight + erretm + rhoinv * Two + abs(temp) * Three + abs(tau) * dw;
	if (w <= Zero) {
	    sg2lb = max(sg2lb, tau);
	} else {
	    sg2ub = min(sg2ub, tau);
	}

	swtch = MFALSE;
	if (orgati) {
	    if (-w > abs(prew) / Ten) {
		swtch = MTRUE;
	    }
	} else {
	    if (w > abs(prew) / Ten) {
		swtch = MTRUE;
	    }
	}
//Main loop to update the values of the array   DELTA and WORK
	iter = niter + 1;
	for (niter = iter; niter <= 20; ++niter) {
//Test for convergence
	    if (abs(w) <= eps * erretm) {
		goto L240;
	    }
//Calculate the new step
	    if (!swtch3) {
		dtipsq = work[ip1] * delta[ip1];
		dtisq = work[i] * delta[i];
		if (!swtch) {
		    if (orgati) {
			c = w - dtipsq * dw + delsq * (z[i] / dtisq * z[i] / dtisq);
		    } else {
			c = w - dtisq * dw - delsq * (z[ip1] / dtipsq * z[ip1] / dtipsq);
		    }
		} else {
		    temp = z[ii] / (work[ii] * delta[ii]);
		    if (orgati) {
			dpsi += temp * temp;
		    } else {
			dphi += temp * temp;
		    }
		    c = w - dtisq * dpsi - dtipsq * dphi;
		}
		a = (dtipsq + dtisq) * w - dtipsq * dtisq * dw;
		b = dtipsq * dtisq * w;
		if (c == Zero) {
		    if (a == Zero) {
			if (!swtch) {
			    if (orgati) {
				a = z[i] * z[i] + dtipsq * dtipsq * (dpsi + dphi);
			    } else {
				a = z[ip1] * z[ip1] + dtisq * dtisq * (dpsi + dphi);
			    }
			} else {
			    a = dtisq * dtisq * dpsi + dtipsq * dtipsq * dphi;
			}
		    }
		    eta = b / a;
		} else if (a <= Zero) {
		    eta = (a - sqrt(abs(a * a - b * Four * c))) / (c * Two);
		} else {
		    eta = b * Two / (a + sqrt(a * a - b * Four * c));
		}
	    } else {
//Interpolation using THREE most relevant poles
		dtiim = work[iim1] * delta[iim1];
		dtiip = work[iip1] * delta[iip1];
		temp = rhoinv + psi + phi;
		if (swtch) {
		    c = temp - dtiim * dpsi - dtiip * dphi;
		    zz[0] = dtiim * dtiim * dpsi;
		    zz[2] = dtiip * dtiip * dphi;
		} else {
		    if (orgati) {
			temp1 = z[iim1] / dtiim;
			temp1 *= temp1;
			temp2 = (d[iim1] - d[iip1]) * (d[iim1] + d[iip1]) * temp1;
			c = temp - dtiip * (dpsi + dphi) - temp2;
			zz[0] = z[iim1] * z[iim1];
			if (dpsi < temp1) {
			    zz[2] = dtiip * dtiip * dphi;
			} else {
			    zz[2] = dtiip * dtiip * (dpsi - temp1 + dphi);
			}
		    } else {
			temp1 = z[iip1] / dtiip;
			temp1 *= temp1;
			temp2 = (d[iip1] - d[iim1]) * (d[iim1] + d[iip1]) * temp1;
			c = temp - dtiim * (dpsi + dphi) - temp2;
			if (dphi < temp1) {
			    zz[0] = dtiim * dtiim * dpsi;
			} else {
			    zz[0] = dtiim * dtiim * (dpsi + (dphi - temp1));
			}
			zz[2] = z[iip1] * z[iip1];
		    }
		}
		dd[0] = dtiim;
		dd[1] = delta[ii] * work[ii];
		dd[2] = dtiip;
		Rlaed6(niter, orgati, c, dd, zz, &w, &eta, info);
		if (*info != 0) {
		    goto L240;
		}
	    }
//Note, eta should be positive if w is negative, and
//eta should be negative otherwise. However,
//if for some reason caused by roundoff, eta*w > 0,
//we simply use one Newton step instead. This way
//will guarantee eta*w < Zero
	    if (w * eta >= Zero) {
		eta = -w / dw;
	    }
	    if (orgati) {
		temp1 = work[i] * delta[i];
		temp = eta - temp1;
	    } else {
		temp1 = work[ip1] * delta[ip1];
		temp = eta - temp1;
	    }
	    if (temp > sg2ub || temp < sg2lb) {
		if (w < Zero) {
		    eta = (sg2ub - tau) / Two;
		} else {
		    eta = (sg2lb - tau) / Two;
		}
	    }
	    tau += eta;
	    eta /= *sigma + sqrt(*sigma * *sigma + eta);
	    *sigma += eta;
	    for (j = 0; j < n; j++) {
		work[j] += eta;
		delta[j] -= eta;

	    }

	    prew = w;
//Evaluate PSI and the derivative DPSI
	    dpsi = Zero;
	    psi = Zero;
	    erretm = Zero;
	    for (j = 0; j < iim1; j++) {
		temp = z[j] / (work[j] * delta[j]);
		psi += z[j] * temp;
		dpsi += temp * temp;
		erretm += psi;

	    }
	    erretm = abs(erretm);
//Evaluate PHI and the derivative DPHI
	    dphi = Zero;
	    phi = Zero;
	    for (j = n; j >= iip1; j--) {
		temp = z[j] / (work[j] * delta[j]);
		phi += z[j] * temp;
		dphi += temp * temp;
		erretm += phi;

	    }
	    temp = z[ii] / (work[ii] * delta[ii]);
	    dw = dpsi + dphi + temp * temp;
	    temp = z[ii] * temp;
	    w = rhoinv + phi + psi + temp;
	    erretm = (phi - psi) * Eight + erretm + rhoinv * Two + abs(temp) * Three + abs(tau) * dw;
	    if (w * prew > Zero && abs(w) > abs(prew) / Ten) {
		swtch = !swtch;
	    }

	    if (w <= Zero) {
		sg2lb = max(sg2lb, tau);
	    } else {
		sg2ub = min(sg2ub, tau);
	    }


	}

//Return with INFO = 1, NITER = MAXIT and not converged
	*info = 1;
    }

  L240:
    return;
}
