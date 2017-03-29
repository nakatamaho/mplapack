/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlagts.cpp,v 1.10 2010/08/07 04:48:32 nakatamaho Exp $ 
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
#include <stdlib.h>		//abs

void Rlagts(INTEGER job, INTEGER n, REAL * a, REAL * b, REAL * c, REAL * d, INTEGER * in, REAL * y, REAL * tol, INTEGER * info)
{
    INTEGER k;
    REAL ak, eps, temp, pert, absak, sfmin;
    REAL bignum;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1, mtemp2, mtemp3;

    *info = 0;

    if (abs(job) > 2 || job == 0) {
	*info = -1;
    } else if (n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	Mxerbla("Rlagts", -(*info));
	return;
    }

    if (n == 0) {
	return;
    }

    eps = Rlamch("Epsilon");
    sfmin = Rlamch("Safe minimum");
    bignum = One / sfmin;

    if (job < 0) {
	if (*tol <= Zero) {
	    *tol = abs(a[1]);
	    if (n > 1) {
		mtemp1 = *tol, mtemp2 = abs(a[2]);
		mtemp3 = max(mtemp1, mtemp2), mtemp1 = abs(b[1]);
		*tol = max(mtemp3, mtemp1);
	    }
	    for (k = 3; k <= n; k++) {
		mtemp1 = *tol, mtemp2 = abs(a[k]);
		mtemp3 = max(mtemp1, mtemp2);
		mtemp2 = abs(b[k - 1]);
		mtemp1 = max(mtemp2, mtemp3);
		mtemp2 = abs(d[k - 2]);
		*tol = max(mtemp1, mtemp2);
	    }
	    *tol = *tol * eps;
	    if (*tol == Zero) {
		*tol = eps;
	    }
	}
    }

    if (abs(job) == 1) {
	for (k = 2; k <= n; k++) {
	    if (in[k - 1] == 0) {
		y[k] = y[k] - c[k - 1] * y[k - 1];
	    } else {
		temp = y[k - 1];
		y[k - 1] = y[k];
		y[k] = temp - c[k - 1] * y[k];
	    }

	}
	if (job == 1) {
	    for (k = n; k >= 1; k--) {
		if (k <= n - 2) {
		    temp = y[k] - b[k] * y[k + 1] - d[k] * y[k + 2];
		} else if (k == n - 1) {
		    temp = y[k] - b[k] * y[k + 1];
		} else {
		    temp = y[k];
		}
		ak = a[k];
		absak = abs(ak);
		if (absak < One) {
		    if (absak < sfmin) {
			if (absak == Zero || abs(temp) * sfmin > absak) {
			    *info = k;
			    return;
			} else {
			    temp *= bignum;
			    ak *= bignum;
			}
		    } else if (abs(temp) > absak * bignum) {
			*info = k;
			return;
		    }
		}
		y[k] = temp / ak;

	    }
	} else {
	    for (k = n; k >= 1; k--) {
		if (k <= n - 2) {
		    temp = y[k] - b[k] * y[k + 1] - d[k] * y[k + 2];
		} else if (k == n - 1) {
		    temp = y[k] - b[k] * y[k + 1];
		} else {
		    temp = y[k];
		}
		ak = a[k];
		pert = sign(*tol, ak);
	      L40:
		absak = abs(ak);
		if (absak < One) {
		    if (absak < sfmin) {
			if (absak == Zero || abs(temp) * sfmin > absak) {
			    ak = ak + pert;
			    pert = pert * 2;
			    goto L40;
			} else {
			    temp = temp * bignum;
			    ak = ak * bignum;
			}
		    } else if (abs(temp) > absak * bignum) {
			ak = ak + pert;
			pert = pert * 2;
			goto L40;
		    }
		}
		y[k] = temp / ak;

	    }
	}
    } else {
//Come to here if  JOB = 2 or -2
	if (job == 2) {
	    for (k = 0; k < n; k++) {
		if (k >= 3) {
		    temp = y[k] - b[k - 1] * y[k - 1] - d[k - 2] * y[k - 2];
		} else if (k == 2) {
		    temp = y[k] - b[k - 1] * y[k - 1];
		} else {
		    temp = y[k];
		}
		ak = a[k];
		absak = abs(ak);
		if (absak < One) {
		    if (absak < sfmin) {
			if (absak == Zero || abs(temp) * sfmin > absak) {
			    *info = k;
			    return;
			} else {
			    temp = temp * bignum;
			    ak = ak * bignum;
			}
		    } else if (abs(temp) > absak * bignum) {
			*info = k;
			return;
		    }
		}
		y[k] = temp / ak;

	    }
	} else {
	    for (k = 0; k < n; k++) {
		if (k >= 3) {
		    temp = y[k] - b[k - 1] * y[k - 1] - d[k - 2] * y[k - 2];
		} else if (k == 2) {
		    temp = y[k] - b[k - 1] * y[k - 1];
		} else {
		    temp = y[k];
		}
		ak = a[k];
		pert = sign(*tol, ak);
	      L70:
		absak = abs(ak);
		if (absak < One) {
		    if (absak < sfmin) {
			if (absak == Zero || abs(temp) * sfmin > absak) {
			    ak += pert;
			    pert *= 2;
			    goto L70;
			} else {
			    temp *= bignum;
			    ak *= bignum;
			}
		    } else if (abs(temp) > absak * bignum) {
			ak += pert;
			pert *= 2;
			goto L70;
		    }
		}
		y[k] = temp / ak;

	    }
	}

	for (k = n; k >= 2; k--) {
	    if (in[k - 1] == 0) {
		y[k - 1] = y[k - 1] - c[k - 1] * y[k];
	    } else {
		temp = y[k - 1];
		y[k - 1] = y[k];
		y[k] = temp - c[k - 1] * y[k];
	    }

	}
    }
    return;
}
