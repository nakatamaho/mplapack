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
#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;
#include <mplapack_lin.h>
#include <mplapack.h>

void Clatb5(const char *path, INTEGER const imat, INTEGER const n, char *type, INTEGER &kl, INTEGER &ku, REAL &anorm, INTEGER &mode, REAL &cndnum, char *dist) {
    FEM_CMN_SVE(Clatb5);
    // SAVE
    REAL &badc1 = sve.badc1;
    REAL &badc2 = sve.badc2;
    REAL &eps = sve.eps;
    bool &first = sve.first;
    REAL &large = sve.large;
    REAL &small = sve.small;
    //
    if (is_called_first_time) {
        first = true;
    }
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
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Save statement ..
    //     ..
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Set some constants for use in the subroutine.
    //
    const REAL tenth = 0.1e+0;
    const REAL one = 1.0;
    const REAL shrink = 0.25e0;
    if (first) {
        first = false;
        eps = Rlamch("Precision");
        badc2 = tenth / eps;
        badc1 = sqrt(badc2);
        small = Rlamch("Safe minimum");
        large = one / small;
        //
        //        If it looks like we're on a Cray, take the square root of
        //        SMALL and LARGE to avoid overflow and underflow problems.
        //
        Rlabad(small, large);
        small = shrink * (small / eps);
        large = one / small;
    }
    //
    char[2] c2 = path[(2 - 1) + (3 - 1) * ldpath];
    //
    //     Set some parameters
    //
    dist = "S";
    mode = 3;
    //
    //     Set TYPE, the type of matrix to be generated.
    //
    type = c2[(1 - 1)];
    //
    //     Set the lower and upper bandwidths.
    //
    if (imat == 1) {
        kl = 0;
    } else {
        kl = max(n - 1, 0);
    }
    ku = kl;
    //
    //     Set the condition number and norm.etc
    //
    const REAL two = 2.0e+0;
    if (imat == 3) {
        cndnum = 1.0e12;
        mode = 2;
    } else if (imat == 4) {
        cndnum = 1.0e12;
        mode = 1;
    } else if (imat == 5) {
        cndnum = 1.0e12;
        mode = 3;
    } else if (imat == 6) {
        cndnum = badc1;
    } else if (imat == 7) {
        cndnum = badc2;
    } else {
        cndnum = two;
    }
    //
    if (imat == 8) {
        anorm = small;
    } else if (imat == 9) {
        anorm = large;
    } else {
        anorm = one;
    }
    //
    if (n <= 1) {
        cndnum = one;
    }
    //
    //     End of Clatb5
    //
}
