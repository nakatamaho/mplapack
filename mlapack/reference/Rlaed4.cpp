/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlaed4.cpp,v 1.7 2010/08/07 04:48:32 nakatamaho Exp $ 
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

#define MTRUE 0
#define MFALSE 1
void Rlaed4(INTEGER n, INTEGER i, REAL * d, REAL * z, REAL * delta, REAL rho, REAL * dlam, INTEGER * info)
{
    REAL a, b, c;
    INTEGER j;
    REAL w;
    INTEGER ii;
    REAL dw, zz[3];
    INTEGER ip1;
    REAL del, eta, phi, eps, tau, psi;
    INTEGER iim1, iip1;
    REAL dphi, dpsi;
    INTEGER iter;
    REAL temp, prew, temp1, dltlb, dltub, midpt;
    INTEGER niter;
    INTEGER swtch, swtch3;
    INTEGER orgati;
    REAL erretm, rhoinv;
    REAL Two = 2.0, Three = 3.0, Four = 4.0, One = 1.0, Zero = 0.0, Eight = 8.0, Ten = 10.0;

//Since this routine is called in an inner loop, we do no argument
//checking.
//Quick return for N=1 and 2
    *info = 0;
    if (n == 1) {

//Presumably, I=1 upon entry
	*dlam = d[1] + rho * z[1] * z[1];
	delta[1] = One;
	return;
    }
    if (n == 2) {
	Rlaed5(i, &d[0], &z[1], &delta[1], rho, dlam);
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
	midpt = rho / Two;

//If ||Z||_2 is not one, then TEMP should be set to
//RHO * ||Z||_2^2 / TWO
	for (j = 0; j < n; j++) {
	    delta[j] = d[j] - d[i] - midpt;
	}

	psi = Zero;
	for (j = 0; j < n - 2; j++) {
	    psi += z[j] * z[j] / delta[j];
	}

	c = rhoinv + psi;
	w = c + z[ii] * z[ii] / delta[ii] + z[n] * z[n] / delta[n];

	if (w <= Zero) {
	    temp = z[n - 1] * z[n - 1] / (d[n] - d[n - 1] + rho)
		+ z[n] * z[n] / rho;
	    if (c <= temp) {
		tau = rho;
	    } else {
		del = d[n] - d[n - 1];
		a = -c * del + z[n - 1] * z[n - 1] + z[n] * z[n];
		b = z[n] * z[n] * del;
		if (a < Zero) {
		    tau = b * Two / (sqrt(a * a + b * Four * c) - a);
		} else {
		    tau = (a + sqrt(a * a + b * Four * c)) / (c * Two);
		}
	    }
//It can be proved that
//D(N)+RHO/2 <= LAMBDA(N) < D(N)+TAU <= D(N)+RHO
	    dltlb = midpt;
	    dltub = rho;
	} else {
	    del = d[n] - d[n - 1];
	    a = -c * del + z[n - 1] * z[n - 1] + z[n] * z[n];
	    b = z[n] * z[n] * del;
	    if (a < Zero) {
		tau = b * Two / (sqrt(a * a + b * Four * c) - a);
	    } else {
		tau = (a + sqrt(a * a + b * Four * c)) / (c * Two);
	    }
//It can be proved that
//D(N) < D(N)+TAU < LAMBDA(N) < D(N)+RHO/2

	    dltlb = Zero;
	    dltub = midpt;
	}
	for (j = 0; j < n; j++) {
	    delta[j] = d[j] - d[i] - tau;
	}
//Evaluate PSI and the derivative DPSI
	dpsi = Zero;
	psi = Zero;
	erretm = Zero;
	for (j = 0; j < ii; j++) {
	    temp = z[j] / delta[j];
	    psi += z[j] * temp;
	    dpsi += temp * temp;
	    erretm += psi;

	}
	erretm = abs(erretm);
//Evaluate PHI and the derivative DPHI
	temp = z[n] / delta[n];
	phi = z[n] * temp;
	dphi = temp * temp;
	erretm = (-phi - psi) * Eight + erretm - phi + rhoinv + abs(tau) * (dpsi + dphi);

	w = rhoinv + phi + psi;

//Test for convergence
	if (abs(w) <= eps * erretm) {
	    *dlam = d[i] + tau;
	    goto L250;
	}
	if (w <= Zero) {
	    dltlb = max(dltlb, tau);
	} else {
	    dltub = min(dltub, tau);
	}
//Calculate the new step
	++niter;
	c = w - delta[n - 1] * dpsi - delta[n] * dphi;
	a = (delta[n - 1] + delta[n]) * w - delta[n - 1] * delta[n] * (dpsi + dphi);
	b = delta[n - 1] * delta[n] * w;
	if (c < Zero) {
	    c = abs(c);
	}
	if (c == Zero) {
//ETA = B/A
//ETA = RHO - TAU
	    eta = dltub - tau;
	} else if (a >= Zero) {
	    eta = (a + sqrt(abs(a * a - b * Four * c))) / (c * Two);
	} else {
	    eta = b * Two / (a - sqrt((abs(a * a - b * Four * c))));
	}

//Note, eta should be positive if w is negative, and
//eta should be negative otherwise. However,
//for some reason caused by roundoff, eta*w > 0,
//we simply use one Newton step instead. This way
//will guarantee eta*w < Zero
	if (w * eta > Zero) {
	    eta = -w / (dpsi + dphi);
	}
	temp = tau + eta;
	if (temp > dltub || temp < dltlb) {
	    if (w < Zero) {
		eta = (dltub - tau) / Two;
	    } else {
		eta = (dltlb - tau) / Two;
	    }
	}
	for (j = 0; j < n; j++) {
	    delta[j] -= eta;
	}
	tau += eta;
//Evaluate PSI and the derivative DPSI
	dpsi = Zero;
	psi = Zero;
	erretm = Zero;
	for (j = 0; j < ii; j++) {
	    temp = z[j] / delta[j];
	    psi += z[j] * temp;
	    dpsi += temp * temp;
	    erretm += psi;
	}
	erretm = abs(erretm);

//Evaluate PHI and the derivative DPHI
	temp = z[n] / delta[n];
	phi = z[n] * temp;
	dphi = temp * temp;
	erretm = (-phi - psi) * Eight + erretm - phi + rhoinv + abs(tau) * (dpsi + dphi);

	w = rhoinv + phi + psi;
//Main loop to update the values of the array   DELTA
	iter = niter + 1;
	for (niter = iter; niter <= 30; ++niter) {
//Test for convergence
	    if (abs(w) <= eps * erretm) {
		*dlam = d[i] + tau;
		goto L250;
	    }
	    if (w <= Zero) {
		dltlb = max(dltlb, tau);
	    } else {
		dltub = min(dltub, tau);
	    }
//Calculate the new step
	    c = w - delta[n - 1] * dpsi - delta[n] * dphi;
	    a = (delta[n - 1] + delta[n]) * w - delta[n - 1] * delta[n] * (dpsi + dphi);
	    b = delta[n - 1] * delta[n] * w;
	    if (a >= Zero) {
		eta = (a + sqrt(abs(a * a - b * Four * c))) / (c * Two);
	    } else {
		eta = b * Two / (a - sqrt(abs(a * a - b * Four * c)));
	    }
/*           Note, eta should be positive if w is negative, and */
/*           eta should be negative otherwise. However, */
/*           if for some reason caused by roundoff, eta*w > 0, */
/*           we simply use one Newton step instead. This way */
/*           will guarantee eta*w < Zero */
	    if (w * eta > Zero) {
		eta = -w / (dpsi + dphi);
	    }
	    temp = tau + eta;
	    if (temp > dltub || temp < dltlb) {
		if (w < Zero) {
		    eta = (dltub - tau) / Two;
		} else {
		    eta = (dltlb - tau) / Two;
		}
	    }
	    for (j = 0; j < n; j++) {
		delta[j] -= eta;
	    }
	    tau += eta;
//Evaluate PSI and the derivative DPSI
	    dpsi = Zero;
	    psi = Zero;
	    erretm = Zero;
	    for (j = 0; j < ii; j++) {
		temp = z[j] / delta[j];
		psi += z[j] * temp;
		dpsi += temp * temp;
		erretm += psi;
	    }
	    erretm = abs(erretm);
//Evaluate PHI and the derivative DPHI
	    temp = z[n] / delta[n];
	    phi = z[n] * temp;
	    dphi = temp * temp;
	    erretm = (-phi - psi) * 8. + erretm - phi + rhoinv + abs(tau) * (dpsi + dphi);
	    w = rhoinv + phi + psi;
	}
//Return with INFO = 1, NITER = MAXIT and not converged
	*info = 1;
	*dlam = d[i] + tau;
	goto L250;
//End for the case I = N
    } else {
//The case for I < N
	niter = 1;
	ip1 = i + 1;
//Calculate initial guess
	del = d[ip1] - d[i];
	midpt = del / Two;
	for (j = 0; j < n; j++) {
	    delta[j] = d[j] - d[i] - midpt;
	}
	psi = Zero;
	for (j = 0; j < i - 1; j++) {
	    psi += z[j] * z[j] / delta[j];
	}
	phi = Zero;
	for (j = n; j >= i + 2; j--) {
	    phi += z[j] * z[j] / delta[j];
	}
	c = rhoinv + psi + phi;
	w = c + z[i] * z[i] / delta[i] + z[ip1] * z[ip1] / delta[ip1];

	if (w > Zero) {
//d(i)< the ith eigenvalue < (d(i)+d(i+1))/2
//We choose d(i) as origin.
	    orgati = MTRUE;
	    a = c * del + z[i] * z[i] + z[ip1] * z[ip1];
	    b = z[i] * z[i] * del;
	    if (a > Zero) {
		tau = b * Two / (a + sqrt(abs(a * a - b * Four * c)));
	    } else {
		tau = (a - sqrt(abs(a * a - b * Four * c))) / (c * Two);
	    }
	    dltlb = Zero;
	    dltub = midpt;
	} else {
//(d(i)+d(i+1))/2 <= the ith eigenvalue < d(i+1)
//We choose d(i+1) as origin.
	    orgati = MFALSE;
	    a = c * del - z[i] * z[i] - z[ip1] * z[ip1];
	    b = z[ip1] * z[ip1] * del;
	    if (a < Zero) {
		tau = b * Two / (a - sqrt(abs(a * a + b * Four * c)));
	    } else {
		tau = -(a + sqrt(abs(a * a + b * Four * c))) / (c * Two);
	    }
	    dltlb = -midpt;
	    dltub = Zero;
	}
	if (orgati) {
	    for (j = 0; j < n; j++) {
		delta[j] = d[j] - d[i] - tau;

	    }
	} else {
	    for (j = 0; j < n; j++) {
		delta[j] = d[j] - d[ip1] - tau;
	    }
	}
	if (orgati) {
	    ii = i;
	} else {
	    ii = i + 1;
	}
	iim1 = ii - 1;
	iip1 = ii + 1;
//Evaluate PSI and the derivative DPSI
	dpsi = Zero;
	psi = Zero;
	erretm = Zero;
	for (j = 0; j < iim1; j++) {
	    temp = z[j] / delta[j];
	    psi += z[j] * temp;
	    dpsi += temp * temp;
	    erretm += psi;
	}
	erretm = abs(erretm);
//Evaluate PHI and the derivative DPHI
	dphi = Zero;
	phi = Zero;
	for (j = n; j >= iip1; j--) {
	    temp = z[j] / delta[j];
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
	temp = z[ii] / delta[ii];
	dw = dpsi + dphi + temp * temp;
	temp = z[ii] * temp;
	w += temp;
	erretm = (phi - psi) * Eight + erretm + rhoinv * Two + abs(temp) * Three + abs(tau) * dw;
//Test for convergence
	if (abs(w) <= eps * erretm) {
	    if (orgati) {
		*dlam = d[i] + tau;
	    } else {
		*dlam = d[ip1] + tau;
	    }
	    goto L250;
	}
	if (w <= Zero) {
	    dltlb = max(dltlb, tau);
	} else {
	    dltub = min(dltub, tau);
	}
//Calculate the new step
	++niter;
	if (!swtch3) {
	    if (orgati) {
		c = w - delta[ip1] * dw - (d[i] - d[ip1]) * (z[i] / delta[i] * z[i] / delta[i]);
	    } else {
		c = w - delta[i] * dw - (d[ip1] - d[i]) * (z[ip1] / delta[ip1] * z[ip1] / delta[ip1]);
	    }
	    a = (delta[i] + delta[ip1]) * w - delta[i] * delta[ip1] * dw;
	    b = delta[i] * delta[ip1] * w;
	    if (c == Zero) {
		if (a == Zero) {
		    if (orgati) {
			a = z[i] * z[i] + delta[ip1] * delta[ip1] * (dpsi + dphi);
		    } else {
			a = z[ip1] * z[ip1] + delta[i] * delta[i] * (dpsi + dphi);
		    }
		}
		eta = b / a;
	    } else if (a <= Zero) {
		eta = (a - sqrt(abs(a * a - b * Four * c))) / (c * Two);
	    } else {
		eta = b * Two / (a + sqrt(abs(a * a - b * Four * c)));
	    }
	} else {

/*           Interpolation using THREE most relevant poles */

	    temp = rhoinv + psi + phi;
	    if (orgati) {
		temp1 = z[iim1] / delta[iim1];
		temp1 *= temp1;
		c = temp - delta[iip1] * (dpsi + dphi) - (d[iim1] - d[iip1]) * temp1;
		zz[0] = z[iim1] * z[iim1];
		zz[2] = delta[iip1] * delta[iip1] * (dpsi - temp1 + dphi);
	    } else {
		temp1 = z[iip1] / delta[iip1];
		temp1 *= temp1;
		c = temp - delta[iim1] * (dpsi + dphi) - (d[iip1] - d[iim1]) * temp1;
		zz[0] = delta[iim1] * delta[iim1] * (dpsi + (dphi - temp1));
		zz[2] = z[iip1] * z[iip1];
	    }
	    zz[1] = z[ii] * z[ii];
	    Rlaed6(niter, orgati, c, &delta[iim1], zz, &w, &eta, info);
	    if (*info != 0) {
		goto L250;
	    }
	}

/*        Note, eta should be positive if w is negative, and */
/*        eta should be negative otherwise. However, */
/*        if for some reason caused by roundoff, eta*w > 0, */
/*        we simply use one Newton step instead. This way */
/*        will guarantee eta*w < Zero */

	if (w * eta >= Zero) {
	    eta = -w / dw;
	}
	temp = tau + eta;
	if (temp > dltub || temp < dltlb) {
	    if (w < Zero) {
		eta = (dltub - tau) / Two;
	    } else {
		eta = (dltlb - tau) / Two;
	    }
	}
	prew = w;
	for (j = 0; j < n; j++) {
	    delta[j] -= eta;
	}
//Evaluate PSI and the derivative DPSI
	dpsi = Zero;
	psi = Zero;
	erretm = Zero;
	for (j = 0; j < iim1; j++) {
	    temp = z[j] / delta[j];
	    psi += z[j] * temp;
	    dpsi += temp * temp;
	    erretm += psi;

	}
	erretm = abs(erretm);

//Evaluate PHI and the derivative DPHI
	dphi = Zero;
	phi = Zero;
	for (j = n; j >= iip1; j--) {
	    temp = z[j] / delta[j];
	    phi += z[j] * temp;
	    dphi += temp * temp;
	    erretm += phi;
	}

	temp = z[ii] / delta[ii];
	dw = dpsi + dphi + temp * temp;
	temp = z[ii] * temp;
	w = rhoinv + phi + psi + temp;
	erretm = (phi - psi) * Eight + erretm + rhoinv * Two + abs(temp) * Three + abs(tau + eta) * dw;
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
	tau += eta;
//Main loop to update the values of the array   DELTA
	iter = niter + 1;
	for (niter = iter; niter <= 30; ++niter) {
//Test for convergence
	    if (abs(w) <= eps * erretm) {
		if (orgati) {
		    *dlam = d[i] + tau;
		} else {
		    *dlam = d[ip1] + tau;
		}
		goto L250;
	    }
	    if (w <= Zero) {
		dltlb = max(dltlb, tau);
	    } else {
		dltub = min(dltub, tau);
	    }
//Calculate the new step
	    if (!swtch3) {
		if (!swtch) {
		    if (orgati) {
			c = w - delta[ip1] * dw - (d[i] - d[ip1]) * ((z[i] / delta[i]) * (z[i] / delta[i]));
		    } else {
			c = w - delta[i] * dw - (d[ip1] - d[i]) * ((z[ip1] / delta[ip1]) * (z[ip1] / delta[ip1]));
		    }
		} else {
		    temp = z[ii] / delta[ii];
		    if (orgati) {
			dpsi += temp * temp;
		    } else {
			dphi += temp * temp;
		    }
		    c = w - delta[i] * dpsi - delta[ip1] * dphi;
		}
		a = (delta[i] + delta[ip1]) * w - delta[i] * delta[ip1]
		    * dw;
		b = delta[i] * delta[ip1] * w;
		if (c == Zero) {
		    if (a == Zero) {
			if (!swtch) {
			    if (orgati) {
				a = z[i] * z[i] + delta[ip1] * delta[ip1] * (dpsi + dphi);
			    } else {
				a = z[ip1] * z[ip1] + delta[i] * delta[i] * (dpsi + dphi);
			    }
			} else {
			    a = delta[i] * delta[i] * dpsi + delta[ip1]
				* delta[ip1] * dphi;
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
		temp = rhoinv + psi + phi;
		if (swtch) {
		    c = temp - delta[iim1] * dpsi - delta[iip1] * dphi;
		    zz[0] = delta[iim1] * delta[iim1] * dpsi;
		    zz[2] = delta[iip1] * delta[iip1] * dphi;
		} else {
		    if (orgati) {
			temp1 = z[iim1] / delta[iim1];
			temp1 *= temp1;
			c = temp - delta[iip1] * (dpsi + dphi) - (d[iim1]
								  - d[iip1]) * temp1;
			zz[0] = z[iim1] * z[iim1];
			zz[2] = delta[iip1] * delta[iip1] * (dpsi - temp1 + dphi);
		    } else {
			temp1 = z[iip1] / delta[iip1];
			temp1 *= temp1;
			c = temp - delta[iim1] * (dpsi + dphi) - (d[iip1]
								  - d[iim1]) * temp1;
			zz[0] = delta[iim1] * delta[iim1] * (dpsi + (dphi - temp1));
			zz[2] = z[iip1] * z[iip1];
		    }
		}
		Rlaed6(niter, orgati, c, &delta[iim1], zz, &w, &eta, info);
		if (*info != 0) {
		    goto L250;
		}
	    }

/*           Note, eta should be positive if w is negative, and */
/*           eta should be negative otherwise. However, */
/*           if for some reason caused by roundoff, eta*w > 0, */
/*           we simply use one Newton step instead. This way */
/*           will guarantee eta*w < Zero */
	    if (w * eta >= Zero) {
		eta = -w / dw;
	    }
	    temp = tau + eta;
	    if (temp > dltub || temp < dltlb) {
		if (w < Zero) {
		    eta = (dltub - tau) / Two;
		} else {
		    eta = (dltlb - tau) / Two;
		}
	    }
	    for (j = 0; j < n; j++) {
		delta[j] -= eta;
	    }
	    tau += eta;
	    prew = w;
//Evaluate PSI and the derivative DPSI
	    dpsi = Zero;
	    psi = Zero;
	    erretm = Zero;
	    for (j = 0; j < iim1; j++) {
		temp = z[j] / delta[j];
		psi += z[j] * temp;
		dpsi += temp * temp;
		erretm += psi;

	    }
	    erretm = abs(erretm);
//Evaluate PHI and the derivative DPHI
	    dphi = Zero;
	    phi = Zero;
	    for (j = n; j >= iip1; j--) {
		temp = z[j] / delta[j];
		phi += z[j] * temp;
		dphi += temp * temp;
		erretm += phi;
	    }
	    temp = z[ii] / delta[ii];
	    dw = dpsi + dphi + temp * temp;
	    temp = z[ii] * temp;
	    w = rhoinv + phi + psi + temp;
	    erretm = (phi - psi) * Eight + erretm + rhoinv * Two + abs(temp) * Three + abs(tau) * dw;
	    if (w * prew > Zero && abs(w) > abs(prew) / Ten) {
		swtch = !swtch;
	    }
	}
//Return with INFO = 1, NITER = MAXIT and not converged
	*info = 1;
	if (orgati) {
	    *dlam = d[i] + tau;
	} else {
	    *dlam = d[ip1] + tau;
	}

    }

  L250:
    return;
}
