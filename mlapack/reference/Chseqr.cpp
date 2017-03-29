/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Chseqr.cpp,v 1.5 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void
Chseqr(const char *job, const char *compz, INTEGER n, INTEGER ilo, INTEGER ihi, COMPLEX * h, INTEGER ldh, COMPLEX * w, COMPLEX * z,
       INTEGER ldz, COMPLEX * work, INTEGER lwork, INTEGER * info)
{
    COMPLEX hl[2401];
    INTEGER kbot, nmin;
    INTEGER initz;
    COMPLEX workl[49];
    INTEGER wantt, wantz;
    INTEGER lquery;
    REAL Zero = 0.0, One = 1.0;
    REAL mtemp1, mtemp2;
    char job_compz[3];

//==== Decode and check the input parameters. ====
    wantt = Mlsame(job, "S");
    initz = Mlsame(compz, "I");
    wantz = initz || Mlsame(compz, "V");
    work[1] = max((INTEGER) 1, n);
    lquery = lwork == -1;
    *info = 0;
    if (!Mlsame(job, "E") && !wantt) {
	*info = -1;
    } else if (!Mlsame(compz, "N") && !wantz) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (ilo < 1 || ilo > max((INTEGER) 1, n)) {
	*info = -4;
    } else if (ihi < min(ilo, n) || ihi > n) {
	*info = -5;
    } else if (ldh < max((INTEGER) 1, n)) {
	*info = -7;
    } else if (ldz < 1 || (wantz && ldz < max((INTEGER) 1, n))) {
	*info = -10;
    } else if (lwork < max((INTEGER) 1, n) && !lquery) {
	*info = -12;
    }
    if (*info != 0) {
//==== Quick return in case of invalid argument. ====
	Mxerbla("Chseqr", -(*info));
	return;
    } else if (n == 0) {
//==== Quick return in case N = 0; nothing to do. ====
	return;
    } else if (lquery) {
//==== Quick return in case of a workspace query ====
	Claqr0(wantt, wantz, n, ilo, ihi, &h[0], ldh, &w[1], &ilo, &ihi, &z[0], ldz, &work[0], lwork, info);
//==== Ensure reported workspace size is backward-compatible with
//.    previous LAPACK versions. ====
	mtemp1 = work[1].real(), mtemp2 = max((INTEGER) 1, n);
	work[1] = max(mtemp1, mtemp2);
	return;
    } else {
//==== copy eigenvalues isolated by ZGEBAL ====
	if (ilo > 1) {
	    Ccopy(ilo - 1, &h[0], ldh + 1, &w[1], 1);
	}
	if (ihi < n) {
	    Ccopy(n - ihi, &h[ihi + 1 + (ihi + 1) * ldh], ldh + 1, &w[ihi + 1], 1);
	}
//==== Initialize Z, if requested ====
	if (initz) {
	    Claset("A", n, n, Zero, One, &z[0], ldz);
	}
//==== Quick return if possible ====
	if (ilo == ihi) {
	    w[ilo] = h[ilo + ilo * ldh];
	    return;
	}
//==== ZLAHQR/ZLAQR0 crossover point ====
	job_compz[0] = job[0];
	job_compz[1] = compz[0];
	job_compz[2] = '\0';
	nmin = iMlaenv(1, "Chseqr", job_compz, n, ilo, ihi, lwork);
	nmin = max((INTEGER) 1, nmin);
//==== ZLAQR0 for big matrices; ZLAHQR for small ones ====
	if (n > nmin) {
	    Claqr0(wantt, wantz, n, ilo, ihi, &h[0], ldh, &w[1], &ilo, &ihi, &z[0], ldz, &work[0], lwork, info);
	} else {
//==== Small matrix ====
	    Clahqr(wantt, wantz, n, ilo, ihi, &h[0], ldh, &w[1], ilo, ihi, &z[0], ldz, info);
	    if (*info > 0) {
//==== A rare ZLAHQR failure!  ZLAQR0 sometimes succeeds
//.    when ZLAHQR fails. ====
		kbot = *info;
		if (n >= 49) {
//==== Larger matrices have enough subdiagonal scratch
//.    space to call ZLAQR0 directly. ====
		    Claqr0(wantt, wantz, n, ilo, kbot, &h[0], ldh, &w[1], &ilo, &ihi, &z[0], ldz, &work[0], lwork, info);
		} else {

//==== Tiny matrices don't have enough subdiagonal
//.    scratch space to benefit from ZLAQR0  Hence,
//.    tiny matrices must be copied into a larger
//.    array before calling ZLAQR0 ====
		    Clacpy("A", n, n, &h[0], ldh, hl, 49);
		    hl[n + 1 + n * 49 - 50] = Zero;
		    Claset("A", 49, 49 - n, Zero, Zero, &hl[(n + 1) * 49 - 49], 49);
		    Claqr0(wantt, wantz, 49, ilo, kbot, hl, 49, &w[1], &ilo, &ihi, &z[0], ldz, workl, 49, info);
		    if (wantt || *info != 0) {
			Clacpy("A", n, n, hl, 49, &h[0], ldh);
		    }
		}
	    }
	}
//==== Clear out the trash, if necessary. ====
	if ((wantt || *info != 0) && n > 2) {
	    Claset("L", n - 2, n - 2, Zero, Zero, &h[ldh + 3], ldh);
	}
//==== Ensure reported workspace size is backward-compatible with
//.    previous LAPACK versions. ====
	mtemp1 = max((INTEGER) 1, n), mtemp2 = work[1].real();
	work[1] = max(mtemp1, mtemp2);
    }
    return;
}
