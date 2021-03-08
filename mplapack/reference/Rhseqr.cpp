/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rhseqr.cpp,v 1.12 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void
Rhseqr(const char *job, const char *compz, INTEGER n, INTEGER ilo, INTEGER ihi,
       REAL * h, INTEGER ldh, REAL * wr, REAL * wi, REAL * z, INTEGER ldz, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i;
    REAL hl[2401];
    INTEGER kbot, nmin;
    INTEGER initz;
    REAL workl[49];
    LOGICAL wantt, wantz;
    INTEGER lquery;
    char ch[3];
    REAL Zero = 0.0, One = 1.0;

//==== Decode and check the input parameters. ====
    wantt = Mlsame(job, "S");
    initz = Mlsame(compz, "I");
    wantz = initz || Mlsame(compz, "V");
    work[1] = (REAL) double (max((INTEGER) 1, n));
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
	*info = -11;
    } else if (lwork < max((INTEGER) 1, n) && !lquery) {
	*info = -13;
    }
    if (*info != 0) {
//==== Quick return in case of invalid argument. ====
	Mxerbla("Rhseqr", -(*info));
	return;
    } else if (n == 0) {
//==== Quick return in case N = 0; nothing to do. ====
	return;
    } else if (lquery) {
//==== Quick return in case of a workspace query ====
	Rlaqr0(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, ilo, ihi, z, ldz, work, lwork, info);
//==== Ensure reported workspace size is backward-compatible with
//.    previous LAPACK versions. ====
	work[1] = max(max((INTEGER) 1, n), (INTEGER) cast2double(work[1]));
	return;
    } else {
//==== copy eigenvalues isolated by DGEBAL ====
	for (i = 0; i < ilo - 1; i++) {
	    wr[i] = h[i + i * ldh];
	    wi[i] = Zero;
	}
	for (i = ihi + 1; i <= n; i++) {
	    wr[i] = h[i + i * ldh];
	    wi[i] = Zero;
	}
//==== Initialize Z, if requested ====
	if (initz) {
	    Rlaset("A", n, n, Zero, One, &z[0], ldz);
	}
//==== Quick return if possible ====
	if (ilo == ihi) {
	    wr[ilo] = h[ilo + ilo * ldh];
	    wi[ilo] = Zero;
	    return;
	}
//==== DLAHQR/DLAQR0 crossover point ====
	ch[0] = (*job);
	ch[1] = (*compz);
	ch[2] = '\0';
	nmin = iMlaenv(12, "Rhseqr", ch, n, ilo, ihi, lwork);
	nmin = max((INTEGER) 11, nmin);
//==== DLAQR0 for big matrices; DLAHQR for small ones ====
	if (n > nmin) {
	    Rlaqr0(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, ilo, ihi, z, ldz, work, lwork, info);
	} else {
//==== Small matrix ====
	    Rlahqr(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, ilo, ihi, z, ldz, info);
	    if (*info > 0) {
//==== A rare DLAHQR failure!  DLAQR0 sometimes succeeds
//.    when DLAHQR fails. ====
		kbot = *info;
		if (n >= 49) {
//==== Larger matrices have enough subdiagonal scratch
//.    space to call DLAQR0 directly. ====
		    Rlaqr0(wantt, wantz, n, ilo, kbot, h, ldh, wr, wi, ilo, ihi, z, ldz, work, lwork, info);
		} else {
//==== Tiny matrices don't have enough subdiagonal
//.    scratch space to benefit from DLAQRZero  Hence,
//.    tiny matrices must be copied into a larger
//.    array before calling DLAQR0 ====
		    Rlacpy("A", n, n, h, ldh, hl, 49);
		    hl[n + 1 + n * 49 - 50] = Zero;
		    Rlaset("A", 49, 49 - n, Zero, Zero, &hl[(n + 1) * 49 - 49], 49);
		    Rlaqr0(wantt, wantz, 49, ilo, kbot, hl, 49, wr, wi, ilo, ihi, z, ldz, workl, 49, info);
		    if (wantt || *info != 0) {
			Rlacpy("A", n, n, hl, 49, h, ldh);
		    }
		}
	    }
	}
//==== Clear out the trash, if necessary. ====
	if ((wantt || *info != 0) && n > 2) {
	    Rlaset("L", n - 2, n - 2, Zero, Zero, &h[ldh + 3], ldh);
	}
//==== Ensure reported workspace size is backward-compatible with
//.    previous LAPACK versions. ====
	work[1] = max(max((INTEGER) 1, n), (INTEGER) cast2double(work[1]));
    }
    return;
}
