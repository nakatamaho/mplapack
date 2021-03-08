/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlazq4.cpp,v 1.5 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rlazq4(INTEGER i0, INTEGER n0, REAL * z, INTEGER pp, INTEGER n0in, REAL dmin, REAL dmin1, REAL dmin2, REAL dn, REAL dn1, REAL dn2, REAL * tau, INTEGER * ttype, REAL * g)
{
    REAL s = 0.0, a2, b1, b2;
    INTEGER i4, nn, np;
    REAL gam, gap1, gap2;
    REAL Cnst1 = 9 / 16, Cnst2 = 1.01, Cnst3 = 1.05;
    REAL Qurtr = 0.25, Third = 1 / 3, Half = 0.5, Zero = 0.0, One = 1.0;
    REAL Two = 2.0, Hundrd = 100.0;
    REAL mtemp1, mtemp2;

//A negative DMIN forces the shift to take that absolute value
//TTYPE records the type of shift.
    if (dmin <= Zero) {
	*tau = -(dmin);
	*ttype = -1;
	return;
    }
    nn = (n0 << 2) + pp;
    if (n0in == n0) {
//No eigenvalues deflated.
	if (dmin == dn || dmin == dn1) {
	    b1 = sqrt(z[nn - 3]) * sqrt(z[nn - 5]);
	    b2 = sqrt(z[nn - 7]) * sqrt(z[nn - 9]);
	    a2 = z[nn - 7] + z[nn - 5];
//Cases 2 and 3.
	    if (dmin == dn && dmin1 == dn1) {
		gap2 = dmin2 - a2 - dmin2 * Qurtr;
		if (gap2 > Zero && gap2 > b2) {
		    gap1 = a2 - dn - b2 / gap2 * b2;
		} else {
		    gap1 = a2 - dn - (b1 + b2);
		}
		if (gap1 > Zero && gap1 > b1) {
		    mtemp1 = dn - b1 / gap1 * b1, mtemp2 = dmin * Half;
		    s = max(mtemp1, mtemp2);
		    *ttype = -2;
		} else {
		    s = Zero;
		    if (dn > b1) {
			s = dn - b1;
		    }
		    if (a2 > b1 + b2) {
			mtemp1 = s, mtemp2 = a2 - (b1 + b2);
			s = min(mtemp1, mtemp2);
		    }
		    mtemp1 = s, mtemp2 = dmin * Third;
		    s = max(mtemp1, mtemp2);
		    *ttype = -3;
		}
	    } else {
//Case Four 
		*ttype = -4;
		s = dmin * Qurtr;
		if (dmin == dn) {
		    gam = dn;
		    a2 = Zero;
		    if (z[nn - 5] > z[nn - 7]) {
			return;
		    }
		    b2 = z[nn - 5] / z[nn - 7];
		    np = nn - 9;
		} else {
		    np = nn - (pp << 1);
		    b2 = z[np - 2];
		    gam = dn1;
		    if (z[np - 4] > z[np - 2]) {
			return;
		    }
		    a2 = z[np - 4] / z[np - 2];
		    if (z[nn - 9] > z[nn - 11]) {
			return;
		    }
		    b2 = z[nn - 9] / z[nn - 11];
		    np = nn - 13;
		}
//Approximate contribution to norm squared from I < NN-One
		a2 += b2;
		for (i4 = np; i4 >= (i0 << 2) - 1 + pp; i4 += -4) {
		    if (b2 == Zero) {
			goto L20;
		    }
		    b1 = b2;
		    if (z[i4] > z[i4 - 2]) {
			return;
		    }
		    b2 *= z[i4] / z[i4 - 2];
		    a2 += b2;
		    if (max(b2, b1) * Hundrd < a2 || Cnst1 < a2) {
			goto L20;
		    }

		}
	      L20:
		a2 *= Cnst3;
//Rayleigh quotient residual bound.
		if (a2 < Cnst1) {
		    s = gam * (One - sqrt(a2)) / (a2 + One);
		}
	    }
	} else if (dmin == dn2) {
//Case 5.
	    *ttype = -5;
	    s = dmin * Qurtr;
//Compute contribution to norm squared from I > NN-Two
	    np = nn - (pp << 1);
	    b1 = z[np - 2];
	    b2 = z[np - 6];
	    gam = dn2;
	    if (z[np - 8] > b2 || z[np - 4] > b1) {
		return;
	    }
	    a2 = z[np - 8] / b2 * (z[np - 4] / b1 + One);
//Approximate contribution to norm squared from I < NN-Two
	    if (n0 - i0 > 2) {
		b2 = z[nn - 13] / z[nn - 15];
		a2 += b2;
		for (i4 = nn - 17; i4 >= (i0 << 2) - 1 + pp; i4 += -4) {
		    if (b2 == Zero) {
			goto L40;
		    }
		    b1 = b2;
		    if (z[i4] > z[i4 - 2]) {
			return;
		    }
		    b2 *= z[i4] / z[i4 - 2];
		    a2 += b2;
		    if (max(b2, b1) * Hundrd < a2 || Cnst1 < a2) {
			goto L40;
		    }
		}
	      L40:
		a2 *= Cnst3;
	    }
	    if (a2 < Cnst1) {
		s = gam * (One - sqrt(a2)) / (a2 + One);
	    }
	} else {
// Case 6, no information to guide us.
	    if (*ttype == -6) {
		*g += (One - *g) * Third;
	    } else if (*ttype == -18) {
		*g = Qurtr * Third;
	    } else {
		*g = Qurtr;
	    }
	    s = *g * dmin;
	    *ttype = -6;
	}
    } else if (n0in == n0 + 1) {
//One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN.
	if (dmin1 == dn1 && dmin2 == dn2) {
//Cases 7 and 8.
	    *ttype = -7;
	    s = dmin1 * Third;
	    if (z[nn - 5] > z[nn - 7]) {
		return;
	    }
	    b1 = z[nn - 5] / z[nn - 7];
	    b2 = b1;
	    if (b2 == Zero) {
		goto L60;
	    }
	    for (i4 = (n0 << 2) - 9 + pp; i4 >= (i0 << 2) - 1 + pp; i4 += -4) {
		a2 = b1;
		if (z[i4] > z[i4 - 2]) {
		    return;
		}
		b1 *= z[i4] / z[i4 - 2];
		b2 += b1;
		if (max(b1, a2) * Hundrd < b2) {
		    goto L60;
		}
	    }
	  L60:
	    b2 = sqrt(b2 * Cnst3);
	    a2 = dmin1 / (b2 * b2 + One);
	    gap2 = dmin2 * Half - a2;
	    if (gap2 > Zero && gap2 > b2 * a2) {
		mtemp1 = s, mtemp2 = a2 * (One - a2 * Cnst2 * (b2 / gap2) * b2);
		s = max(mtemp1, mtemp2);
	    } else {
		mtemp1 = s, mtemp2 = a2 * (One - b2 * Cnst2);
		s = max(mtemp1, mtemp2);
		*ttype = -8;
	    }
	} else {
//Case 9.
	    s = dmin1 * Qurtr;
	    if (dmin1 == dn1) {
		s = dmin1 * Half;
	    }
	    *ttype = -9;
	}
    } else if (n0in == n0 + 2) {
//Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN.
//Cases 10 and 1O
	if (dmin2 == dn2 && z[nn - 5] * Two < z[nn - 7]) {
	    *ttype = -10;
	    s = dmin2 * Third;
	    if (z[nn - 5] > z[nn - 7]) {
		return;
	    }
	    b1 = z[nn - 5] / z[nn - 7];
	    b2 = b1;
	    if (b2 == Zero) {
		goto L80;
	    }
	    for (i4 = (n0 << 2) - 9 + pp; i4 >= (i0 << 2) - 1 + pp; i4 += -4) {
		if (z[i4] > z[i4 - 2]) {
		    return;
		}
		b1 *= z[i4] / z[i4 - 2];
		b2 += b1;
		if (b1 * Hundrd < b2) {
		    goto L80;
		}

	    }
	  L80:
	    b2 = sqrt(b2 * Cnst3);
	    a2 = dmin2 / (b2 * b2 + One);
	    gap2 = z[nn - 7] + z[nn - 9] - sqrt(z[nn - 11]) * sqrt(z[nn - 9]) - a2;
	    if (gap2 > Zero && gap2 > b2 * a2) {
		mtemp1 = s, mtemp2 = a2 * (One - a2 * Cnst2 * (b2 / gap2) * b2);
		s = max(mtemp1, mtemp2);
	    } else {
		mtemp1 = s, mtemp2 = a2 * (One - b2 * Cnst2);
		s = max(mtemp1, mtemp2);
	    }
	} else {
	    s = dmin2 * Qurtr;
	    *ttype = -11;
	}
    } else if (n0in > n0 + 2) {
//Case 12, more than two eigenvalues deflated. No information.
	s = Zero;
	*ttype = -12;
    }

    *tau = s;
    return;
}
