/*
 * Copyright (c) 2008-2025
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

struct common : fem::common {
    fem::cmn_sve drotm_sve;

    common(int argc, char const *argv[]) : fem::common(argc, argv) {}
};

struct drotm_save {
    double two;
    double zero;

    drotm_save() : two(0.0), zero(0.0) {}
};

// > \brief \b DROTM
// >
// > \verbatim
// >
// >    APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
// >
// >    (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN
// >    (DY**T)
// >
// >    DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE
// >    LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.
// >    WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
// >
// >    DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
// >
// >      (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
// >    H=(          )    (          )    (          )    (          )
// >      (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
// >    SEE DROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.
// > \endverbatim
//
// Arguments:
// > \param[in] N
// > \verbatim
// >          N is INTEGER
// >         number of elements in input vector(s)
// > \endverbatim
// >
// > \param[in,out] DX
// > \verbatim
// >          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
// > \endverbatim
// >
// > \param[in] INCX
// > \verbatim
// >          INCX is INTEGER
// >         storage spacing between elements of DX
// > \endverbatim
// >
// > \param[in,out] DY
// > \verbatim
// >          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
// > \endverbatim
// >
// > \param[in] INCY
// > \verbatim
// >          INCY is INTEGER
// >         storage spacing between elements of DY
// > \endverbatim
// >
// > \param[in] DPARAM
// > \verbatim
// >          DPARAM is DOUBLE PRECISION array, dimension (5)
// >     DPARAM(1)=DFLAG
// >     DPARAM(2)=DH11
// >     DPARAM(3)=DH21
// >     DPARAM(4)=DH12
// >     DPARAM(5)=DH22
// > \endverbatim
//
// Authors:
void Rrotm(common &cmn, INTEGER const &n, REAL *dx, INTEGER const &incx, REAL *dy, INTEGER const &incy, REAL *dparam) {
    FEM_CMN_SVE(drotm);
    // SAVE
    double &two = sve.two;
    double &zero = sve.zero;
    //
    if (is_called_first_time) {
        zero = 0.0;
        two = 2.0;
    }
    //
    REAL dflag = dparam[1 - 1];
    if (n <= 0 || (dflag + two == zero)) {
        return;
    }
    INTEGER nsteps = 0;
    REAL dh11 = 0.0;
    REAL dh12 = 0.0;
    REAL dh21 = 0.0;
    REAL dh22 = 0.0;
    INTEGER i = 0;
    REAL w = 0.0;
    REAL z = 0.0;
    INTEGER kx = 0;
    INTEGER ky = 0;
    if (incx == incy && incx > 0) {
        //
        nsteps = n * incx;
        if (dflag < zero) {
            dh11 = dparam[2 - 1];
            dh12 = dparam[4 - 1];
            dh21 = dparam[3 - 1];
            dh22 = dparam[5 - 1];
            for (i = 1; i <= nsteps; i = i + incx) {
                w = dx[i - 1];
                z = dy[i - 1];
                dx[i - 1] = w * dh11 + z * dh12;
                dy[i - 1] = w * dh21 + z * dh22;
            }
        } else if (dflag == zero) {
            dh12 = dparam[4 - 1];
            dh21 = dparam[3 - 1];
            for (i = 1; i <= nsteps; i = i + incx) {
                w = dx[i - 1];
                z = dy[i - 1];
                dx[i - 1] = w + z * dh12;
                dy[i - 1] = w * dh21 + z;
            }
        } else {
            dh11 = dparam[2 - 1];
            dh22 = dparam[5 - 1];
            for (i = 1; i <= nsteps; i = i + incx) {
                w = dx[i - 1];
                z = dy[i - 1];
                dx[i - 1] = w * dh11 + z;
                dy[i - 1] = -w + dh22 * z;
            }
        }
    } else {
        kx = 1;
        ky = 1;
        if (incx < 0) {
            kx = 1 + (1 - n) * incx;
        }
        if (incy < 0) {
            ky = 1 + (1 - n) * incy;
        }
        //
        if (dflag < zero) {
            dh11 = dparam[2 - 1];
            dh12 = dparam[4 - 1];
            dh21 = dparam[3 - 1];
            dh22 = dparam[5 - 1];
            for (i = 1; i <= n; i = i + 1) {
                w = dx[kx - 1];
                z = dy[ky - 1];
                dx[kx - 1] = w * dh11 + z * dh12;
                dy[ky - 1] = w * dh21 + z * dh22;
                kx += incx;
                ky += incy;
            }
        } else if (dflag == zero) {
            dh12 = dparam[4 - 1];
            dh21 = dparam[3 - 1];
            for (i = 1; i <= n; i = i + 1) {
                w = dx[kx - 1];
                z = dy[ky - 1];
                dx[kx - 1] = w + z * dh12;
                dy[ky - 1] = w * dh21 + z;
                kx += incx;
                ky += incy;
            }
        } else {
            dh11 = dparam[2 - 1];
            dh22 = dparam[5 - 1];
            for (i = 1; i <= n; i = i + 1) {
                w = dx[kx - 1];
                z = dy[ky - 1];
                dx[kx - 1] = w * dh11 + z;
                dy[ky - 1] = -w + dh22 * z;
                kx += incx;
                ky += incy;
            }
        }
    }
}
