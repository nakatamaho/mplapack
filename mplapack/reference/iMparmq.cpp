/*
 * Copyright (c) 2008-2021
 *      Nakata, Maho
 *      All rights reserved.
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

#include <mpblas.h>
#include <mplapack.h>

INTEGER iMparmq(INTEGER const ispec, const char *name, const char *opts, INTEGER const n, INTEGER const ilo, INTEGER const ihi, INTEGER const lwork) {
    INTEGER return_value = 0;
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //
    //  ================================================================
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    const INTEGER ishfts = 15;
    const INTEGER inwin = 13;
    const INTEGER iacc22 = 16;
    INTEGER nh = 0;
    INTEGER ns = 0;
    const float two = 2.0f;
    if ((ispec == ishfts) || (ispec == inwin) || (ispec == iacc22)) {
        //
        //        ==== Set the number simultaneous shifts ====
        //
        nh = ihi - ilo + 1;
        ns = 2;
        if (nh >= 30) {
            ns = 4;
        }
        if (nh >= 60) {
            ns = 10;
        }
        if (nh >= 150) {
            ns = max((INTEGER)10, nh / nint[(log(castREAL(nh - 1)) / log(two)) - 1]);
        }
        if (nh >= 590) {
            ns = 64;
        }
        if (nh >= 3000) {
            ns = 128;
        }
        if (nh >= 6000) {
            ns = 256;
        }
        ns = max((INTEGER)2, ns - mod(ns, 2));
    }
    //
    const INTEGER inmin = 12;
    const INTEGER nmin = 75;
    const INTEGER inibl = 14;
    const INTEGER nibble = 14;
    const INTEGER knwswp = 500;
    str<6> subnam = char0;
    INTEGER ic = 0;
    INTEGER iz = 0;
    INTEGER i = 0;
    const INTEGER k22min = 14;
    const INTEGER kacmin = 14;
    if (ispec == inmin) {
        //
        //        ===== Matrices of order smaller than NMIN get sent
        //        .     to xLAHQR, the classic REAL shift algorithm.
        //        .     This must be at least 11. ====
        //
        return_value = nmin;
        //
    } else if (ispec == inibl) {
        //
        //        ==== INIBL: skip a multi-shift qr iteration and
        //        .    whenever aggressive early deflation finds
        //        .    at least (NIBBLE*(window size)/100) deflations. ====
        //
        return_value = nibble;
        //
    } else if (ispec == ishfts) {
        //
        //        ==== NSHFTS: The number of simultaneous shifts =====
        //
        return_value = ns;
        //
    } else if (ispec == inwin) {
        //
        //        ==== NW: deflation window size.  ====
        //
        if (nh <= knwswp) {
            return_value = ns;
        } else {
            return_value = 3 * ns / 2;
        }
        //
    } else if (ispec == iacc22) {
        //
        //        ==== IACC22: Whether to accumulate reflections
        //        .     before updating the far-from-diagonal elements
        //        .     and whether to use 2-by-2 block structure while
        //        .     doing it.  A small amount of work could be saved
        //        .     by making this choice dependent also upon the
        //        .     NH=IHI-ILO+1.
        //
        //        Convert NAME to upper case if the first character is lower case.
        //
        return_value = 0;
        subnam = name;
        ic = *subnam;
        iz = 'Z';
        if (iz == 90 || iz == 122) {
            //
            //           ASCII character set
            //
	    if (ic >= 97 && ic <= 122) {
		*subnam = (char) (ic - 32);
		for (i = 2; i <= 6; i++) {
		    ic = *&subnam[i - 1];
		    if (ic >= 97 && ic <= 122) {
			*subnam[i - 1] = (char) (ic - 32);
		    }
		}
	    }
        }
        //
	if (strncmp(subnam + 1, "GGHRD", 5) == 0 || strncmp(
		subnam + 1, "GGHD3", 5) == 0) {
	    return_value = 1;
	    if (nh >= 14) {
		return_value = 2;
	    }
	} else if (strncmp(subnam + 3, "EXC", 3) == 0) {
	    if (nh >= 14) {
		return_value = 1;
	    }
	    if (nh >= 14) {
		return_value = 2;
	    }
	} else if (strncmp(subnam + 1, "HSEQR", 5) == 0 || 
		strncmp(subnam + 1, "LAQR", 4) == 0) {
	    if (ns >= 14) {
		return_value = 1;
	    }
	    if (ns >= 14) {
		return_value = 2;
	    }
	}
        //
    } else {
        //        ===== invalid value of ispec =====
        return_value = -1;
        //
    }
    return return_value;
    //
    //     ==== End of IPARMQ ====
    //
}
