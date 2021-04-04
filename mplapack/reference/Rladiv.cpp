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

//
//     .. Parameters ..
//
//     .. Local Scalars ..
//     ..
//     .. Executable Statements ..
//
const REAL zero = 0.0;
REAL br = 0.0;
if (r != zero) {
    br = b * r;
    if (br != zero) {
        return_value = (a + br) * t;
    } else {
        return_value = a * t + (b * t) * r;
    }
} else {
    return_value = (a + d * (b / c)) * t;
}
//
return return_value;
//
//     End of Rladiv12
//
}

//
//> \ingroup REALOTHERauxiliary
//
void Rladiv1(REAL &a, REAL const &b, REAL const &c, REAL const &d, REAL &p, REAL &q) {
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Parameters ..
    //
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    REAL r = d / c;
    const REAL one = 1.0;
    REAL t = one / (c + d * r);
    p = Rladiv2[(a - 1) + (b - 1) * ldRladiv2];
    a = -a;
    q = Rladiv2[(b - 1) + (a - 1) * ldRladiv2];
    //
    //     End of Rladiv1
    //
}

//> \brief \b Rladiv performs complex division in real arithmetic, avoiding unnecessary overflow.
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//> \htmlonly
//> Download Rladiv + dependencies
//> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/Rladiv.f">
//> [TGZ]</a>
//> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/Rladiv.f">
//> [ZIP]</a>
//> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/Rladiv.f">
//> [TXT]</a>
//> \endhtmlonly
//
//  Definition:
//  ===========
//
//       SUBROUTINE Rladiv( A, B, C, D, P, Q )
//
//       .. Scalar Arguments ..
//       DOUBLE PRECISION   A, B, C, D, P, Q
//       ..
//
//> \par Purpose:
//  =============
//>
//> \verbatim
//>
//> Rladiv performs complex division in  real arithmetic
//>
//>                       a + i*b
//>            p + i*q = ---------
//>                       c + i*d
//>
//> The algorithm is due to Michael Baudin and Robert L. Smith
//> and can be found in the paper
//> "A Robust Complex Division in Scilab"
//> \endverbatim
//
//  Arguments:
//  ==========
//
//> \param[in] A
//> \verbatim
//>          A is DOUBLE PRECISION
//> \endverbatim
//>
//> \param[in] B
//> \verbatim
//>          B is DOUBLE PRECISION
//> \endverbatim
//>
//> \param[in] C
//> \verbatim
//>          C is DOUBLE PRECISION
//> \endverbatim
//>
//> \param[in] D
//> \verbatim
//>          D is DOUBLE PRECISION
//>          The scalars a, b, c, and d in the above expression.
//> \endverbatim
//>
//> \param[out] P
//> \verbatim
//>          P is DOUBLE PRECISION
//> \endverbatim
//>
//> \param[out] Q
//> \verbatim
//>          Q is DOUBLE PRECISION
//>          The scalars p and q in the above expression.
//> \endverbatim
//
//  Authors:
//  ========
//
//> \author Univ. of Tennessee
//> \author Univ. of California Berkeley
//> \author Univ. of Colorado Denver
//> \author NAG Ltd.
//
//> \ingroup REALOTHERauxiliary
//
//  =====================================================================
void Rladiv(REAL const &a, REAL const &b, REAL const &c, REAL const &d, REAL &p, REAL &q) {
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Parameters ..
    //
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    REAL aa = a;
    REAL bb = b;
    REAL cc = c;
    REAL dd = d;
    REAL ab = max(abs(a), abs(b));
    REAL cd = max(abs(c), abs(d));
    REAL s = 1.0;
    //
    REAL ov = dlamch("Overflow threshold");
    REAL un = dlamch("Safe minimum");
    REAL eps = dlamch("Epsilon");
    const REAL bs = 2.0;
    REAL be = bs / (eps * eps);
    //
    const REAL half = 0.5e0;
    const REAL two = 2.0;
    if (ab >= half * ov) {
        aa = half * aa;
        bb = half * bb;
        s = two * s;
    }
    if (cd >= half * ov) {
        cc = half * cc;
        dd = half * dd;
        s = half * s;
    }
    if (ab <= un * bs / eps) {
        aa = aa * be;
        bb = bb * be;
        s = s / be;
    }
    if (cd <= un * bs / eps) {
        cc = cc * be;
        dd = dd * be;
        s = s * be;
    }
    if (abs(d) <= abs(c)) {
        Rladiv1(aa, bb, cc, dd, p, q);
    } else {
        Rladiv1(bb, aa, dd, cc, p, q);
        q = -q;
    }
    p = p * s;
    q = q * s;
    //
    //     End of Rladiv
    //
}
