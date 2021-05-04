/*
 * Copyright (c) 2021
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

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_eig.h>

#include <mplapack_debug.h>

INTEGER iMlaenv(INTEGER const ispec, const char *name, const char *opts, INTEGER const n1, INTEGER const n2, INTEGER const n3, INTEGER const n4) {
    INTEGER return_value = 0;
    // COMMON claenv
    //
    //
    //  -- LAPACK test routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Arrays in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Save statement ..
    //     ..
    //     .. Executable Statements ..
    //
    if (ispec >= 1 && ispec <= 5) {
        //
        //        Return a value from the common block.
        //
        return_value = iparms[ispec - 1];
        //
    } else if (ispec == 6) {
        //
        //        Compute SVD crossover point.
        //
        return_value = int(real[(min(n1 - 1) + (n2)-1) * ldreal] * 1.6e0f);
        //
    } else if (ispec >= 7 && ispec <= 9) {
        //
        //        Return a value from the common block.
        //
        return_value = iparms[ispec - 1];
        //
    } else if (ispec == 10) {
        //
        //        IEEE NaN arithmetic can be trusted not to trap
        //
        //        iMlaenv = 0
        return_value = 1;
        if (return_value == 1) {
            return_value = ieeeck1, 0.0f, 1.0f;
        }
        //
    } else if (ispec == 11) {
        //
        //        Infinity arithmetic can be trusted not to trap
        //
        //        iMlaenv = 0
        return_value = 1;
        if (return_value == 1) {
            return_value = ieeeck0, 0.0f, 1.0f;
        }
        //
    } else if ((ispec >= 12) && (ispec <= 16)) {
        //
        //     12 <= ISPEC <= 16: xHSEQR or one of its subroutines.
        //
        return_value = iparms[ispec - 1];
        //         WRITE(*,*) 'ISPEC = ',ISPEC,' iMlaenv =',iMlaenv
        //         iMlaenv = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
        //
    } else if ((ispec >= 17) && (ispec <= 21)) {
        //
        //     17 <= ISPEC <= 21: 2stage eigenvalues SVD routines.
        //
        if (ispec == 17) {
            return_value = iparms[1 - 1];
        } else {
            return_value = iparam2stage(ispec, name, opts, n1, n2, n3, n4);
        }
        //
    } else {
        //
        //        Invalid value for ISPEC
        //
        return_value = -1;
    }
    //
    return return_value;
    //
    //     End of iMlaenv
    //
}

INTEGER iMlaenv2stage(INTEGER const ispec, const char *name, const char *opts, INTEGER const n1, INTEGER const n2, INTEGER const n3, INTEGER const n4) {
    INTEGER return_value = 0;
    // COMMON claenv
    //
    //     .. Scalar Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Local variables ..
    //     .. External Functions ..
    //     ..
    //     .. Arrays in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Save statement ..
    //     ..
    //     .. Executable Statements ..
    //
    INTEGER iispec = 0;
    if ((ispec >= 1) && (ispec <= 5)) {
        //
        //     1 <= ISPEC <= 5: 2stage eigenvalues SVD routines.
        //
        if (ispec == 1) {
            return_value = iparms[1 - 1];
        } else {
            iispec = 16 + ispec;
            return_value = iparam2stage(iispec, name, opts, n1, n2, n3, n4);
        }
        //
    } else {
        //
        //        Invalid value for ISPEC
        //
        return_value = -1;
    }
    //
    return return_value;
}

INTEGER iparmq(INTEGER const ispec, const char * /* name */, const char * /* opts */, INTEGER const  /* n */, INTEGER const ilo, INTEGER const ihi, INTEGER const  /* lwork */) {
    INTEGER return_value = 0;
    //
    //     ..
    //     .. Scalar Arguments ..
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
            ns = max((INTEGER)10, nh / nint(log(real[nh - 1]) / log(two)));
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
        ns = max({(INTEGER)2, ns - mod(ns, 2)});
    }
    //
    const INTEGER inmin = 12;
    const INTEGER nmin = 11;
    const INTEGER inibl = 14;
    const INTEGER nibble = 14;
    const INTEGER knwswp = 500;
    const INTEGER kacmin = 14;
    const INTEGER k22min = 14;
    if (ispec == inmin) {
        //
        //        ===== Matrices of order smaller than NMIN get sent
        //        .     to LAHQR, the classic REAL shift algorithm.
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
        return_value = 0;
        if (ns >= kacmin) {
            return_value = 1;
        }
        if (ns >= k22min) {
            return_value = 2;
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
