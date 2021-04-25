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

COMPLEX
Clarnd(INTEGER const idist, INTEGER *iseed) {
    COMPLEX return_value = (0.0, 0.0);
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Generate a pair of real random numbers from a uniform (0,1)
    //     distribution
    //
    REAL t1 = Rlaran[iseed - 1];
    REAL t2 = Rlaran[iseed - 1];
    //
    const REAL two = 2.0e+0;
    const REAL one = 1.0;
    const REAL zero = 0.0;
    const REAL twopi = 6.28318530717958647692528676655900576839e+0;
    if (idist == 1) {
        //
        //        real and imaginary parts each uniform (0,1)
        //
        return_value = COMPLEX(t1, t2);
    } else if (idist == 2) {
        //
        //        real and imaginary parts each uniform (-1,1)
        //
        return_value = COMPLEX(two * t1 - one, two * t2 - one);
    } else if (idist == 3) {
        //
        //        real and imaginary parts each normal (0,1)
        //
        return_value = sqrt(-two * log(t1)) * exp(COMPLEX(zero, twopi * t2));
    } else if (idist == 4) {
        //
        //        uniform distribution on the unit disc abs(z) <= 1
        //
        return_value = sqrt(t1) * exp(COMPLEX(zero, twopi * t2));
    } else if (idist == 5) {
        //
        //        uniform distribution on the unit circle abs(z) = 1
        //
        return_value = exp(COMPLEX(zero, twopi * t2));
    }
    return return_value;
    //
    //     End of Clarnd
    //
}
