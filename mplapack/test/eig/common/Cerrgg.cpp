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
#include <mplapack_lin.h>
#include <mplapack_eig.h>

#include <mplapack_debug.h>

bool _Clctsx(COMPLEX const /* alpha */, COMPLEX const /* beta */) {
    bool return_value = false;
    return return_value;
}
bool _Clctes(COMPLEX const z, COMPLEX const d) {
    bool return_value = false;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    const REAL zero = 0.0;
    const REAL one = 1.0;
    REAL zmax = 0.0;
    if (d == czero) {
        return_value = (z.real() < zero);
    } else {
        if (z.real() == zero || d.real() == zero) {
            return_value = (sign(one, z.imag()) != sign(one, d.imag()));
        } else if (z.imag() == zero || d.imag() == zero) {
            return_value = (sign(one, z.real()) != sign(one, d.real()));
        } else {
            zmax = max(abs(z.real()), abs(z.imag()));
            return_value = ((z.real() / zmax) * d.real() + (z.imag() / zmax) * d.imag() < zero);
        }
    }
    //
    return return_value;
    //
    //     End of Clctes
    //
}

void Cerrgg(const char *path, INTEGER const nunit) {
    common cmn;
    common_write write(cmn);
    INTEGER infot;
    INTEGER nout;
    bool ok;
    bool lerr;
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
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Scalars in Common ..
    //     ..
    //     .. Common blocks ..
    //     ..
    //     .. Executable Statements ..
    //
    nout = nunit;
    char c2[2];
    c2[0] = path[1];
    c2[1] = path[2];
    //
    //     Set the variables to innocuous values.
    //
    INTEGER j = 0;
    const INTEGER nmax = 3;
    bool sel[nmax];
    INTEGER i = 0;
    const REAL zero = 0.0;
    COMPLEX a[nmax * nmax];
    COMPLEX b[nmax * nmax];
    INTEGER lda = nmax;
    INTEGER ldb = nmax;
    for (j = 1; j <= nmax; j = j + 1) {
        sel[j - 1] = true;
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = zero;
            b[(i - 1) + (j - 1) * ldb] = zero;
        }
    }
    const REAL one = 1.0;
    for (i = 1; i <= nmax; i = i + 1) {
        a[(i - 1) + (i - 1) * lda] = one;
        b[(i - 1) + (i - 1) * ldb] = one;
    }
    ok = true;
    REAL tola = 1.0;
    REAL tolb = 1.0;
    INTEGER ifst = 1;
    INTEGER ilst = 1;
    INTEGER nt = 0;
    INTEGER lwork = 1;
    //
    //     Test error exits for the GG path.
    //
    COMPLEX q[nmax * nmax];
    COMPLEX z[nmax * nmax];
    INTEGER info = 0;
    const INTEGER lw = 6 * nmax;
    COMPLEX w[lw];
    COMPLEX alpha[nmax];
    COMPLEX beta[nmax];
    REAL rw[lw];
    INTEGER m = 0;
    INTEGER dummyk = 0;
    INTEGER dummyl = 0;
    REAL r1[nmax];
    REAL r2[nmax];
    COMPLEX u[nmax * nmax];
    COMPLEX v[nmax * nmax];
    INTEGER idum[nmax];
    INTEGER iw[lw];
    COMPLEX tau[nmax];
    INTEGER ncycle = 0;
    REAL rs[nmax];
    INTEGER sdim = 0;
    bool bw[nmax];
    REAL rce[nmax];
    REAL rcv[nmax];
    INTEGER ilo = 0;
    INTEGER ihi = 0;
    REAL ls[nmax];
    REAL anrm = 0.0;
    REAL bnrm = 0.0;
    REAL scale = 0.0;
    REAL dif = 0.0;
    if (Mlsamen(2, c2, "GG")) {
        //
        //        Cgghrd
        //
        infot = 1;
        Cgghrd("/", "N", 0, 1, 0, a, 1, b, 1, q, 1, z, 1, info);
        chkxer("Cgghrd", infot, nout, lerr, ok);
        infot = 2;
        Cgghrd("N", "/", 0, 1, 0, a, 1, b, 1, q, 1, z, 1, info);
        chkxer("Cgghrd", infot, nout, lerr, ok);
        infot = 3;
        Cgghrd("N", "N", -1, 0, 0, a, 1, b, 1, q, 1, z, 1, info);
        chkxer("Cgghrd", infot, nout, lerr, ok);
        infot = 4;
        Cgghrd("N", "N", 0, 0, 0, a, 1, b, 1, q, 1, z, 1, info);
        chkxer("Cgghrd", infot, nout, lerr, ok);
        infot = 5;
        Cgghrd("N", "N", 0, 1, 1, a, 1, b, 1, q, 1, z, 1, info);
        chkxer("Cgghrd", infot, nout, lerr, ok);
        infot = 7;
        Cgghrd("N", "N", 2, 1, 1, a, 1, b, 2, q, 1, z, 1, info);
        chkxer("Cgghrd", infot, nout, lerr, ok);
        infot = 9;
        Cgghrd("N", "N", 2, 1, 1, a, 2, b, 1, q, 1, z, 1, info);
        chkxer("Cgghrd", infot, nout, lerr, ok);
        infot = 11;
        Cgghrd("V", "N", 2, 1, 1, a, 2, b, 2, q, 1, z, 1, info);
        chkxer("Cgghrd", infot, nout, lerr, ok);
        infot = 13;
        Cgghrd("N", "V", 2, 1, 1, a, 2, b, 2, q, 1, z, 1, info);
        chkxer("Cgghrd", infot, nout, lerr, ok);
        nt += 9;
        //
        //        Cgghd3
        //
        infot = 1;
        Cgghd3("/", "N", 0, 1, 0, a, 1, b, 1, q, 1, z, 1, w, lw, info);
        chkxer("Cgghd3", infot, nout, lerr, ok);
        infot = 2;
        Cgghd3("N", "/", 0, 1, 0, a, 1, b, 1, q, 1, z, 1, w, lw, info);
        chkxer("Cgghd3", infot, nout, lerr, ok);
        infot = 3;
        Cgghd3("N", "N", -1, 0, 0, a, 1, b, 1, q, 1, z, 1, w, lw, info);
        chkxer("Cgghd3", infot, nout, lerr, ok);
        infot = 4;
        Cgghd3("N", "N", 0, 0, 0, a, 1, b, 1, q, 1, z, 1, w, lw, info);
        chkxer("Cgghd3", infot, nout, lerr, ok);
        infot = 5;
        Cgghd3("N", "N", 0, 1, 1, a, 1, b, 1, q, 1, z, 1, w, lw, info);
        chkxer("Cgghd3", infot, nout, lerr, ok);
        infot = 7;
        Cgghd3("N", "N", 2, 1, 1, a, 1, b, 2, q, 1, z, 1, w, lw, info);
        chkxer("Cgghd3", infot, nout, lerr, ok);
        infot = 9;
        Cgghd3("N", "N", 2, 1, 1, a, 2, b, 1, q, 1, z, 1, w, lw, info);
        chkxer("Cgghd3", infot, nout, lerr, ok);
        infot = 11;
        Cgghd3("V", "N", 2, 1, 1, a, 2, b, 2, q, 1, z, 1, w, lw, info);
        chkxer("Cgghd3", infot, nout, lerr, ok);
        infot = 13;
        Cgghd3("N", "V", 2, 1, 1, a, 2, b, 2, q, 1, z, 1, w, lw, info);
        chkxer("Cgghd3", infot, nout, lerr, ok);
        nt += 9;
        //
        //        Chgeqz
        //
        infot = 1;
        Chgeqz("/", "N", "N", 0, 1, 0, a, 1, b, 1, alpha, beta, q, 1, z, 1, w, 1, rw, info);
        chkxer("Chgeqz", infot, nout, lerr, ok);
        infot = 2;
        Chgeqz("E", "/", "N", 0, 1, 0, a, 1, b, 1, alpha, beta, q, 1, z, 1, w, 1, rw, info);
        chkxer("Chgeqz", infot, nout, lerr, ok);
        infot = 3;
        Chgeqz("E", "N", "/", 0, 1, 0, a, 1, b, 1, alpha, beta, q, 1, z, 1, w, 1, rw, info);
        chkxer("Chgeqz", infot, nout, lerr, ok);
        infot = 4;
        Chgeqz("E", "N", "N", -1, 0, 0, a, 1, b, 1, alpha, beta, q, 1, z, 1, w, 1, rw, info);
        chkxer("Chgeqz", infot, nout, lerr, ok);
        infot = 5;
        Chgeqz("E", "N", "N", 0, 0, 0, a, 1, b, 1, alpha, beta, q, 1, z, 1, w, 1, rw, info);
        chkxer("Chgeqz", infot, nout, lerr, ok);
        infot = 6;
        Chgeqz("E", "N", "N", 0, 1, 1, a, 1, b, 1, alpha, beta, q, 1, z, 1, w, 1, rw, info);
        chkxer("Chgeqz", infot, nout, lerr, ok);
        infot = 8;
        Chgeqz("E", "N", "N", 2, 1, 1, a, 1, b, 2, alpha, beta, q, 1, z, 1, w, 1, rw, info);
        chkxer("Chgeqz", infot, nout, lerr, ok);
        infot = 10;
        Chgeqz("E", "N", "N", 2, 1, 1, a, 2, b, 1, alpha, beta, q, 1, z, 1, w, 1, rw, info);
        chkxer("Chgeqz", infot, nout, lerr, ok);
        infot = 14;
        Chgeqz("E", "V", "N", 2, 1, 1, a, 2, b, 2, alpha, beta, q, 1, z, 1, w, 1, rw, info);
        chkxer("Chgeqz", infot, nout, lerr, ok);
        infot = 16;
        Chgeqz("E", "N", "V", 2, 1, 1, a, 2, b, 2, alpha, beta, q, 1, z, 1, w, 1, rw, info);
        chkxer("Chgeqz", infot, nout, lerr, ok);
        nt += 10;
        //
        //        Ctgevc
        //
        infot = 1;
        Ctgevc("/", "A", sel, 0, a, 1, b, 1, q, 1, z, 1, 0, m, w, rw, info);
        chkxer("Ctgevc", infot, nout, lerr, ok);
        infot = 2;
        Ctgevc("R", "/", sel, 0, a, 1, b, 1, q, 1, z, 1, 0, m, w, rw, info);
        chkxer("Ctgevc", infot, nout, lerr, ok);
        infot = 4;
        Ctgevc("R", "A", sel, -1, a, 1, b, 1, q, 1, z, 1, 0, m, w, rw, info);
        chkxer("Ctgevc", infot, nout, lerr, ok);
        infot = 6;
        Ctgevc("R", "A", sel, 2, a, 1, b, 2, q, 1, z, 2, 0, m, w, rw, info);
        chkxer("Ctgevc", infot, nout, lerr, ok);
        infot = 8;
        Ctgevc("R", "A", sel, 2, a, 2, b, 1, q, 1, z, 2, 0, m, w, rw, info);
        chkxer("Ctgevc", infot, nout, lerr, ok);
        infot = 10;
        Ctgevc("L", "A", sel, 2, a, 2, b, 2, q, 1, z, 1, 0, m, w, rw, info);
        chkxer("Ctgevc", infot, nout, lerr, ok);
        infot = 12;
        Ctgevc("R", "A", sel, 2, a, 2, b, 2, q, 1, z, 1, 0, m, w, rw, info);
        chkxer("Ctgevc", infot, nout, lerr, ok);
        infot = 13;
        Ctgevc("R", "A", sel, 2, a, 2, b, 2, q, 1, z, 2, 1, m, w, rw, info);
        chkxer("Ctgevc", infot, nout, lerr, ok);
        nt += 8;
        //
        //     Test error exits for the GSV path.
        //
    } else if (Mlsamen(3, path, "GSV")) {
        //
        //        Cggsvd3
        //
        infot = 1;
        Cggsvd3("/", "N", "N", 0, 0, 0, dummyk, dummyl, a, 1, b, 1, r1, r2, u, 1, v, 1, q, 1, w, lwork, rw, idum, info);
        chkxer("Cggsvd3", infot, nout, lerr, ok);
        infot = 2;
        Cggsvd3("N", "/", "N", 0, 0, 0, dummyk, dummyl, a, 1, b, 1, r1, r2, u, 1, v, 1, q, 1, w, lwork, rw, idum, info);
        chkxer("Cggsvd3", infot, nout, lerr, ok);
        infot = 3;
        Cggsvd3("N", "N", "/", 0, 0, 0, dummyk, dummyl, a, 1, b, 1, r1, r2, u, 1, v, 1, q, 1, w, lwork, rw, idum, info);
        chkxer("Cggsvd3", infot, nout, lerr, ok);
        infot = 4;
        Cggsvd3("N", "N", "N", -1, 0, 0, dummyk, dummyl, a, 1, b, 1, r1, r2, u, 1, v, 1, q, 1, w, lwork, rw, idum, info);
        chkxer("Cggsvd3", infot, nout, lerr, ok);
        infot = 5;
        Cggsvd3("N", "N", "N", 0, -1, 0, dummyk, dummyl, a, 1, b, 1, r1, r2, u, 1, v, 1, q, 1, w, lwork, rw, idum, info);
        chkxer("Cggsvd3", infot, nout, lerr, ok);
        infot = 6;
        Cggsvd3("N", "N", "N", 0, 0, -1, dummyk, dummyl, a, 1, b, 1, r1, r2, u, 1, v, 1, q, 1, w, lwork, rw, idum, info);
        chkxer("Cggsvd3", infot, nout, lerr, ok);
        infot = 10;
        Cggsvd3("N", "N", "N", 2, 1, 1, dummyk, dummyl, a, 1, b, 1, r1, r2, u, 1, v, 1, q, 1, w, lwork, rw, idum, info);
        chkxer("Cggsvd3", infot, nout, lerr, ok);
        infot = 12;
        Cggsvd3("N", "N", "N", 1, 1, 2, dummyk, dummyl, a, 1, b, 1, r1, r2, u, 1, v, 1, q, 1, w, lwork, rw, idum, info);
        chkxer("Cggsvd3", infot, nout, lerr, ok);
        infot = 16;
        Cggsvd3("U", "N", "N", 2, 2, 2, dummyk, dummyl, a, 2, b, 2, r1, r2, u, 1, v, 1, q, 1, w, lwork, rw, idum, info);
        chkxer("Cggsvd3", infot, nout, lerr, ok);
        infot = 18;
        Cggsvd3("N", "V", "N", 2, 2, 2, dummyk, dummyl, a, 2, b, 2, r1, r2, u, 2, v, 1, q, 1, w, lwork, rw, idum, info);
        chkxer("Cggsvd3", infot, nout, lerr, ok);
        infot = 20;
        Cggsvd3("N", "N", "Q", 2, 2, 2, dummyk, dummyl, a, 2, b, 2, r1, r2, u, 2, v, 2, q, 1, w, lwork, rw, idum, info);
        chkxer("Cggsvd3", infot, nout, lerr, ok);
        nt += 11;
        //
        //        Cggsvp3
        //
        infot = 1;
        Cggsvp3("/", "N", "N", 0, 0, 0, a, 1, b, 1, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, rw, tau, w, lwork, info);
        chkxer("Cggsvp3", infot, nout, lerr, ok);
        infot = 2;
        Cggsvp3("N", "/", "N", 0, 0, 0, a, 1, b, 1, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, rw, tau, w, lwork, info);
        chkxer("Cggsvp3", infot, nout, lerr, ok);
        infot = 3;
        Cggsvp3("N", "N", "/", 0, 0, 0, a, 1, b, 1, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, rw, tau, w, lwork, info);
        chkxer("Cggsvp3", infot, nout, lerr, ok);
        infot = 4;
        Cggsvp3("N", "N", "N", -1, 0, 0, a, 1, b, 1, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, rw, tau, w, lwork, info);
        chkxer("Cggsvp3", infot, nout, lerr, ok);
        infot = 5;
        Cggsvp3("N", "N", "N", 0, -1, 0, a, 1, b, 1, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, rw, tau, w, lwork, info);
        chkxer("Cggsvp3", infot, nout, lerr, ok);
        infot = 6;
        Cggsvp3("N", "N", "N", 0, 0, -1, a, 1, b, 1, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, rw, tau, w, lwork, info);
        chkxer("Cggsvp3", infot, nout, lerr, ok);
        infot = 8;
        Cggsvp3("N", "N", "N", 2, 1, 1, a, 1, b, 1, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, rw, tau, w, lwork, info);
        chkxer("Cggsvp3", infot, nout, lerr, ok);
        infot = 10;
        Cggsvp3("N", "N", "N", 1, 2, 1, a, 1, b, 1, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, rw, tau, w, lwork, info);
        chkxer("Cggsvp3", infot, nout, lerr, ok);
        infot = 16;
        Cggsvp3("U", "N", "N", 2, 2, 2, a, 2, b, 2, tola, tolb, dummyk, dummyl, u, 1, v, 1, q, 1, iw, rw, tau, w, lwork, info);
        chkxer("Cggsvp3", infot, nout, lerr, ok);
        infot = 18;
        Cggsvp3("N", "V", "N", 2, 2, 2, a, 2, b, 2, tola, tolb, dummyk, dummyl, u, 2, v, 1, q, 1, iw, rw, tau, w, lwork, info);
        chkxer("Cggsvp3", infot, nout, lerr, ok);
        infot = 20;
        Cggsvp3("N", "N", "Q", 2, 2, 2, a, 2, b, 2, tola, tolb, dummyk, dummyl, u, 2, v, 2, q, 1, iw, rw, tau, w, lwork, info);
        chkxer("Cggsvp3", infot, nout, lerr, ok);
        nt += 11;
        //
        //        Ctgsja
        //
        infot = 1;
        Ctgsja("/", "N", "N", 0, 0, 0, dummyk, dummyl, a, 1, b, 1, tola, tolb, r1, r2, u, 1, v, 1, q, 1, w, ncycle, info);
        chkxer("Ctgsja", infot, nout, lerr, ok);
        infot = 2;
        Ctgsja("N", "/", "N", 0, 0, 0, dummyk, dummyl, a, 1, b, 1, tola, tolb, r1, r2, u, 1, v, 1, q, 1, w, ncycle, info);
        chkxer("Ctgsja", infot, nout, lerr, ok);
        infot = 3;
        Ctgsja("N", "N", "/", 0, 0, 0, dummyk, dummyl, a, 1, b, 1, tola, tolb, r1, r2, u, 1, v, 1, q, 1, w, ncycle, info);
        chkxer("Ctgsja", infot, nout, lerr, ok);
        infot = 4;
        Ctgsja("N", "N", "N", -1, 0, 0, dummyk, dummyl, a, 1, b, 1, tola, tolb, r1, r2, u, 1, v, 1, q, 1, w, ncycle, info);
        chkxer("Ctgsja", infot, nout, lerr, ok);
        infot = 5;
        Ctgsja("N", "N", "N", 0, -1, 0, dummyk, dummyl, a, 1, b, 1, tola, tolb, r1, r2, u, 1, v, 1, q, 1, w, ncycle, info);
        chkxer("Ctgsja", infot, nout, lerr, ok);
        infot = 6;
        Ctgsja("N", "N", "N", 0, 0, -1, dummyk, dummyl, a, 1, b, 1, tola, tolb, r1, r2, u, 1, v, 1, q, 1, w, ncycle, info);
        chkxer("Ctgsja", infot, nout, lerr, ok);
        infot = 10;
        Ctgsja("N", "N", "N", 0, 0, 0, dummyk, dummyl, a, 0, b, 1, tola, tolb, r1, r2, u, 1, v, 1, q, 1, w, ncycle, info);
        chkxer("Ctgsja", infot, nout, lerr, ok);
        infot = 12;
        Ctgsja("N", "N", "N", 0, 0, 0, dummyk, dummyl, a, 1, b, 0, tola, tolb, r1, r2, u, 1, v, 1, q, 1, w, ncycle, info);
        chkxer("Ctgsja", infot, nout, lerr, ok);
        infot = 18;
        Ctgsja("U", "N", "N", 0, 0, 0, dummyk, dummyl, a, 1, b, 1, tola, tolb, r1, r2, u, 0, v, 1, q, 1, w, ncycle, info);
        chkxer("Ctgsja", infot, nout, lerr, ok);
        infot = 20;
        Ctgsja("N", "V", "N", 0, 0, 0, dummyk, dummyl, a, 1, b, 1, tola, tolb, r1, r2, u, 1, v, 0, q, 1, w, ncycle, info);
        chkxer("Ctgsja", infot, nout, lerr, ok);
        infot = 22;
        Ctgsja("N", "N", "Q", 0, 0, 0, dummyk, dummyl, a, 1, b, 1, tola, tolb, r1, r2, u, 1, v, 1, q, 0, w, ncycle, info);
        chkxer("Ctgsja", infot, nout, lerr, ok);
        nt += 11;
        //
        //     Test error exits for the GLM path.
        //
    } else if (Mlsamen(3, path, "GLM")) {
        //
        //        Cggglm
        //
        infot = 1;
        Cggglm(-1, 0, 0, a, 1, b, 1, tau, alpha, beta, w, lw, info);
        chkxer("Cggglm", infot, nout, lerr, ok);
        infot = 2;
        Cggglm(0, -1, 0, a, 1, b, 1, tau, alpha, beta, w, lw, info);
        chkxer("Cggglm", infot, nout, lerr, ok);
        infot = 2;
        Cggglm(0, 1, 0, a, 1, b, 1, tau, alpha, beta, w, lw, info);
        chkxer("Cggglm", infot, nout, lerr, ok);
        infot = 3;
        Cggglm(0, 0, -1, a, 1, b, 1, tau, alpha, beta, w, lw, info);
        chkxer("Cggglm", infot, nout, lerr, ok);
        infot = 3;
        Cggglm(1, 0, 0, a, 1, b, 1, tau, alpha, beta, w, lw, info);
        chkxer("Cggglm", infot, nout, lerr, ok);
        infot = 5;
        Cggglm(0, 0, 0, a, 0, b, 1, tau, alpha, beta, w, lw, info);
        chkxer("Cggglm", infot, nout, lerr, ok);
        infot = 7;
        Cggglm(0, 0, 0, a, 1, b, 0, tau, alpha, beta, w, lw, info);
        chkxer("Cggglm", infot, nout, lerr, ok);
        infot = 12;
        Cggglm(1, 1, 1, a, 1, b, 1, tau, alpha, beta, w, 1, info);
        chkxer("Cggglm", infot, nout, lerr, ok);
        nt += 8;
        //
        //     Test error exits for the LSE path.
        //
    } else if (Mlsamen(3, path, "LSE")) {
        //
        //        Cgglse
        //
        infot = 1;
        Cgglse(-1, 0, 0, a, 1, b, 1, tau, alpha, beta, w, lw, info);
        chkxer("Cgglse", infot, nout, lerr, ok);
        infot = 2;
        Cgglse(0, -1, 0, a, 1, b, 1, tau, alpha, beta, w, lw, info);
        chkxer("Cgglse", infot, nout, lerr, ok);
        infot = 3;
        Cgglse(0, 0, -1, a, 1, b, 1, tau, alpha, beta, w, lw, info);
        chkxer("Cgglse", infot, nout, lerr, ok);
        infot = 3;
        Cgglse(0, 0, 1, a, 1, b, 1, tau, alpha, beta, w, lw, info);
        chkxer("Cgglse", infot, nout, lerr, ok);
        infot = 3;
        Cgglse(0, 1, 0, a, 1, b, 1, tau, alpha, beta, w, lw, info);
        chkxer("Cgglse", infot, nout, lerr, ok);
        infot = 5;
        Cgglse(0, 0, 0, a, 0, b, 1, tau, alpha, beta, w, lw, info);
        chkxer("Cgglse", infot, nout, lerr, ok);
        infot = 7;
        Cgglse(0, 0, 0, a, 1, b, 0, tau, alpha, beta, w, lw, info);
        chkxer("Cgglse", infot, nout, lerr, ok);
        infot = 12;
        Cgglse(1, 1, 1, a, 1, b, 1, tau, alpha, beta, w, 1, info);
        chkxer("Cgglse", infot, nout, lerr, ok);
        nt += 8;
        //
        //     Test error exits for the CSD path.
        //
    } else if (Mlsamen(3, path, "CSD")) {
        //
        //        Cuncsd
        //
        infot = 7;
        Cuncsd("Y", "Y", "Y", "Y", "N", "N", -1, 0, 0, a, 1, a, 1, a, 1, a, 1, rs, a, 1, a, 1, a, 1, a, 1, w, lw, rw, lw, iw, info);
        chkxer("Cuncsd", infot, nout, lerr, ok);
        infot = 8;
        Cuncsd("Y", "Y", "Y", "Y", "N", "N", 1, -1, 0, a, 1, a, 1, a, 1, a, 1, rs, a, 1, a, 1, a, 1, a, 1, w, lw, rw, lw, iw, info);
        chkxer("Cuncsd", infot, nout, lerr, ok);
        infot = 9;
        Cuncsd("Y", "Y", "Y", "Y", "N", "N", 1, 1, -1, a, 1, a, 1, a, 1, a, 1, rs, a, 1, a, 1, a, 1, a, 1, w, lw, rw, lw, iw, info);
        chkxer("Cuncsd", infot, nout, lerr, ok);
        infot = 11;
        Cuncsd("Y", "Y", "Y", "Y", "N", "N", 1, 1, 1, a, -1, a, 1, a, 1, a, 1, rs, a, 1, a, 1, a, 1, a, 1, w, lw, rw, lw, iw, info);
        chkxer("Cuncsd", infot, nout, lerr, ok);
        infot = 20;
        Cuncsd("Y", "Y", "Y", "Y", "N", "N", 1, 1, 1, a, 1, a, 1, a, 1, a, 1, rs, a, -1, a, 1, a, 1, a, 1, w, lw, rw, lw, iw, info);
        chkxer("Cuncsd", infot, nout, lerr, ok);
        infot = 22;
        Cuncsd("Y", "Y", "Y", "Y", "N", "N", 1, 1, 1, a, 1, a, 1, a, 1, a, 1, rs, a, 1, a, -1, a, 1, a, 1, w, lw, rw, lw, iw, info);
        chkxer("Cuncsd", infot, nout, lerr, ok);
        infot = 24;
        Cuncsd("Y", "Y", "Y", "Y", "N", "N", 1, 1, 1, a, 1, a, 1, a, 1, a, 1, rs, a, 1, a, 1, a, -1, a, 1, w, lw, rw, lw, iw, info);
        chkxer("Cuncsd", infot, nout, lerr, ok);
        infot = 26;
        Cuncsd("Y", "Y", "Y", "Y", "N", "N", 1, 1, 1, a, 1, a, 1, a, 1, a, 1, rs, a, 1, a, 1, a, 1, a, -1, w, lw, rw, lw, iw, info);
        chkxer("Cuncsd", infot, nout, lerr, ok);
        nt += 8;
        //
        //     Test error exits for the GQR path.
        //
    } else if (Mlsamen(3, path, "GQR")) {
        //
        //        Cggqrf
        //
        infot = 1;
        Cggqrf(-1, 0, 0, a, 1, alpha, b, 1, beta, w, lw, info);
        chkxer("Cggqrf", infot, nout, lerr, ok);
        infot = 2;
        Cggqrf(0, -1, 0, a, 1, alpha, b, 1, beta, w, lw, info);
        chkxer("Cggqrf", infot, nout, lerr, ok);
        infot = 3;
        Cggqrf(0, 0, -1, a, 1, alpha, b, 1, beta, w, lw, info);
        chkxer("Cggqrf", infot, nout, lerr, ok);
        infot = 5;
        Cggqrf(0, 0, 0, a, 0, alpha, b, 1, beta, w, lw, info);
        chkxer("Cggqrf", infot, nout, lerr, ok);
        infot = 8;
        Cggqrf(0, 0, 0, a, 1, alpha, b, 0, beta, w, lw, info);
        chkxer("Cggqrf", infot, nout, lerr, ok);
        infot = 11;
        Cggqrf(1, 1, 2, a, 1, alpha, b, 1, beta, w, 1, info);
        chkxer("Cggqrf", infot, nout, lerr, ok);
        nt += 6;
        //
        //        Cggrqf
        //
        infot = 1;
        Cggrqf(-1, 0, 0, a, 1, alpha, b, 1, beta, w, lw, info);
        chkxer("Cggrqf", infot, nout, lerr, ok);
        infot = 2;
        Cggrqf(0, -1, 0, a, 1, alpha, b, 1, beta, w, lw, info);
        chkxer("Cggrqf", infot, nout, lerr, ok);
        infot = 3;
        Cggrqf(0, 0, -1, a, 1, alpha, b, 1, beta, w, lw, info);
        chkxer("Cggrqf", infot, nout, lerr, ok);
        infot = 5;
        Cggrqf(0, 0, 0, a, 0, alpha, b, 1, beta, w, lw, info);
        chkxer("Cggrqf", infot, nout, lerr, ok);
        infot = 8;
        Cggrqf(0, 0, 0, a, 1, alpha, b, 0, beta, w, lw, info);
        chkxer("Cggrqf", infot, nout, lerr, ok);
        infot = 11;
        Cggrqf(1, 1, 2, a, 1, alpha, b, 1, beta, w, 1, info);
        chkxer("Cggrqf", infot, nout, lerr, ok);
        nt += 6;
        //
        //     Test error exits for the ZGS, ZGV, ZGX, and ZXV paths.
        //
    } else if (Mlsamen(3, path, "ZGS") || Mlsamen(3, path, "ZGV") || Mlsamen(3, path, "ZGX") || Mlsamen(3, path, "ZXV")) {
        //
        //        Cgges
        //
        infot = 1;
        Cgges("/", "N", "S", _Clctes, 1, a, 1, b, 1, sdim, alpha, beta, q, 1, u, 1, w, 1, rw, bw, info);
        chkxer("Cgges ", infot, nout, lerr, ok);
        infot = 2;
        Cgges("N", "/", "S", _Clctes, 1, a, 1, b, 1, sdim, alpha, beta, q, 1, u, 1, w, 1, rw, bw, info);
        chkxer("Cgges ", infot, nout, lerr, ok);
        infot = 3;
        Cgges("N", "V", "/", _Clctes, 1, a, 1, b, 1, sdim, alpha, beta, q, 1, u, 1, w, 1, rw, bw, info);
        chkxer("Cgges ", infot, nout, lerr, ok);
        infot = 5;
        Cgges("N", "V", "S", _Clctes, -1, a, 1, b, 1, sdim, alpha, beta, q, 1, u, 1, w, 1, rw, bw, info);
        chkxer("Cgges ", infot, nout, lerr, ok);
        infot = 7;
        Cgges("N", "V", "S", _Clctes, 1, a, 0, b, 1, sdim, alpha, beta, q, 1, u, 1, w, 1, rw, bw, info);
        chkxer("Cgges ", infot, nout, lerr, ok);
        infot = 9;
        Cgges("N", "V", "S", _Clctes, 1, a, 1, b, 0, sdim, alpha, beta, q, 1, u, 1, w, 1, rw, bw, info);
        chkxer("Cgges ", infot, nout, lerr, ok);
        infot = 14;
        Cgges("N", "V", "S", _Clctes, 1, a, 1, b, 1, sdim, alpha, beta, q, 0, u, 1, w, 1, rw, bw, info);
        chkxer("Cgges ", infot, nout, lerr, ok);
        infot = 14;
        Cgges("V", "V", "S", _Clctes, 2, a, 2, b, 2, sdim, alpha, beta, q, 1, u, 2, w, 1, rw, bw, info);
        chkxer("Cgges ", infot, nout, lerr, ok);
        infot = 16;
        Cgges("N", "V", "S", _Clctes, 1, a, 1, b, 1, sdim, alpha, beta, q, 1, u, 0, w, 1, rw, bw, info);
        chkxer("Cgges ", infot, nout, lerr, ok);
        infot = 16;
        Cgges("V", "V", "S", _Clctes, 2, a, 2, b, 2, sdim, alpha, beta, q, 2, u, 1, w, 1, rw, bw, info);
        chkxer("Cgges ", infot, nout, lerr, ok);
        infot = 18;
        Cgges("V", "V", "S", _Clctes, 2, a, 2, b, 2, sdim, alpha, beta, q, 2, u, 2, w, 1, rw, bw, info);
        chkxer("Cgges ", infot, nout, lerr, ok);
        nt += 11;
        //
        //        Cgges3
        //
        infot = 1;
        Cgges3("/", "N", "S", _Clctes, 1, a, 1, b, 1, sdim, alpha, beta, q, 1, u, 1, w, 1, rw, bw, info);
        chkxer("Cgges3", infot, nout, lerr, ok);
        infot = 2;
        Cgges3("N", "/", "S", _Clctes, 1, a, 1, b, 1, sdim, alpha, beta, q, 1, u, 1, w, 1, rw, bw, info);
        chkxer("Cgges3", infot, nout, lerr, ok);
        infot = 3;
        Cgges3("N", "V", "/", _Clctes, 1, a, 1, b, 1, sdim, alpha, beta, q, 1, u, 1, w, 1, rw, bw, info);
        chkxer("Cgges3", infot, nout, lerr, ok);
        infot = 5;
        Cgges3("N", "V", "S", _Clctes, -1, a, 1, b, 1, sdim, alpha, beta, q, 1, u, 1, w, 1, rw, bw, info);
        chkxer("Cgges3", infot, nout, lerr, ok);
        infot = 7;
        Cgges3("N", "V", "S", _Clctes, 1, a, 0, b, 1, sdim, alpha, beta, q, 1, u, 1, w, 1, rw, bw, info);
        chkxer("Cgges3", infot, nout, lerr, ok);
        infot = 9;
        Cgges3("N", "V", "S", _Clctes, 1, a, 1, b, 0, sdim, alpha, beta, q, 1, u, 1, w, 1, rw, bw, info);
        chkxer("Cgges3", infot, nout, lerr, ok);
        infot = 14;
        Cgges3("N", "V", "S", _Clctes, 1, a, 1, b, 1, sdim, alpha, beta, q, 0, u, 1, w, 1, rw, bw, info);
        chkxer("Cgges3", infot, nout, lerr, ok);
        infot = 14;
        Cgges3("V", "V", "S", _Clctes, 2, a, 2, b, 2, sdim, alpha, beta, q, 1, u, 2, w, 1, rw, bw, info);
        chkxer("Cgges3", infot, nout, lerr, ok);
        infot = 16;
        Cgges3("N", "V", "S", _Clctes, 1, a, 1, b, 1, sdim, alpha, beta, q, 1, u, 0, w, 1, rw, bw, info);
        chkxer("Cgges3", infot, nout, lerr, ok);
        infot = 16;
        Cgges3("V", "V", "S", _Clctes, 2, a, 2, b, 2, sdim, alpha, beta, q, 2, u, 1, w, 1, rw, bw, info);
        chkxer("Cgges3", infot, nout, lerr, ok);
        infot = 18;
        Cgges3("V", "V", "S", _Clctes, 2, a, 2, b, 2, sdim, alpha, beta, q, 2, u, 2, w, 1, rw, bw, info);
        chkxer("Cgges3", infot, nout, lerr, ok);
        nt += 11;
        //
        //        Cggesx
        //
        infot = 1;
        Cggesx("/", "N", "S", _Clctsx, "N", 1, a, 1, b, 1, sdim, alpha, beta, q, 1, u, 1, rce, rcv, w, 1, rw, iw, 1, bw, info);
        chkxer("Cggesx", infot, nout, lerr, ok);
        infot = 2;
        Cggesx("N", "/", "S", _Clctsx, "N", 1, a, 1, b, 1, sdim, alpha, beta, q, 1, u, 1, rce, rcv, w, 1, rw, iw, 1, bw, info);
        chkxer("Cggesx", infot, nout, lerr, ok);
        infot = 3;
        Cggesx("V", "V", "/", _Clctsx, "N", 1, a, 1, b, 1, sdim, alpha, beta, q, 1, u, 1, rce, rcv, w, 1, rw, iw, 1, bw, info);
        chkxer("Cggesx", infot, nout, lerr, ok);
        infot = 5;
        Cggesx("V", "V", "S", _Clctsx, "/", 1, a, 1, b, 1, sdim, alpha, beta, q, 1, u, 1, rce, rcv, w, 1, rw, iw, 1, bw, info);
        chkxer("Cggesx", infot, nout, lerr, ok);
        infot = 6;
        Cggesx("V", "V", "S", _Clctsx, "B", -1, a, 1, b, 1, sdim, alpha, beta, q, 1, u, 1, rce, rcv, w, 1, rw, iw, 1, bw, info);
        chkxer("Cggesx", infot, nout, lerr, ok);
        infot = 8;
        Cggesx("V", "V", "S", _Clctsx, "B", 1, a, 0, b, 1, sdim, alpha, beta, q, 1, u, 1, rce, rcv, w, 1, rw, iw, 1, bw, info);
        chkxer("Cggesx", infot, nout, lerr, ok);
        infot = 10;
        Cggesx("V", "V", "S", _Clctsx, "B", 1, a, 1, b, 0, sdim, alpha, beta, q, 1, u, 1, rce, rcv, w, 1, rw, iw, 1, bw, info);
        chkxer("Cggesx", infot, nout, lerr, ok);
        infot = 15;
        Cggesx("V", "V", "S", _Clctsx, "B", 1, a, 1, b, 1, sdim, alpha, beta, q, 0, u, 1, rce, rcv, w, 1, rw, iw, 1, bw, info);
        chkxer("Cggesx", infot, nout, lerr, ok);
        infot = 15;
        Cggesx("V", "V", "S", _Clctsx, "B", 2, a, 2, b, 2, sdim, alpha, beta, q, 1, u, 1, rce, rcv, w, 1, rw, iw, 1, bw, info);
        chkxer("Cggesx", infot, nout, lerr, ok);
        infot = 17;
        Cggesx("V", "V", "S", _Clctsx, "B", 1, a, 1, b, 1, sdim, alpha, beta, q, 1, u, 0, rce, rcv, w, 1, rw, iw, 1, bw, info);
        chkxer("Cggesx", infot, nout, lerr, ok);
        infot = 17;
        Cggesx("V", "V", "S", _Clctsx, "B", 2, a, 2, b, 2, sdim, alpha, beta, q, 2, u, 1, rce, rcv, w, 1, rw, iw, 1, bw, info);
        chkxer("Cggesx", infot, nout, lerr, ok);
        infot = 21;
        Cggesx("V", "V", "S", _Clctsx, "B", 2, a, 2, b, 2, sdim, alpha, beta, q, 2, u, 2, rce, rcv, w, 1, rw, iw, 1, bw, info);
        chkxer("Cggesx", infot, nout, lerr, ok);
        infot = 24;
        Cggesx("V", "V", "S", _Clctsx, "V", 1, a, 1, b, 1, sdim, alpha, beta, q, 1, u, 1, rce, rcv, w, 32, rw, iw, 0, bw, info);
        chkxer("Cggesx", infot, nout, lerr, ok);
        nt += 13;
        //
        //        Cggev
        //
        infot = 1;
        Cggev("/", "N", 1, a, 1, b, 1, alpha, beta, q, 1, u, 1, w, 1, rw, info);
        chkxer("Cggev ", infot, nout, lerr, ok);
        infot = 2;
        Cggev("N", "/", 1, a, 1, b, 1, alpha, beta, q, 1, u, 1, w, 1, rw, info);
        chkxer("Cggev ", infot, nout, lerr, ok);
        infot = 3;
        Cggev("V", "V", -1, a, 1, b, 1, alpha, beta, q, 1, u, 1, w, 1, rw, info);
        chkxer("Cggev ", infot, nout, lerr, ok);
        infot = 5;
        Cggev("V", "V", 1, a, 0, b, 1, alpha, beta, q, 1, u, 1, w, 1, rw, info);
        chkxer("Cggev ", infot, nout, lerr, ok);
        infot = 7;
        Cggev("V", "V", 1, a, 1, b, 0, alpha, beta, q, 1, u, 1, w, 1, rw, info);
        chkxer("Cggev ", infot, nout, lerr, ok);
        infot = 11;
        Cggev("N", "V", 1, a, 1, b, 1, alpha, beta, q, 0, u, 1, w, 1, rw, info);
        chkxer("Cggev ", infot, nout, lerr, ok);
        infot = 11;
        Cggev("V", "V", 2, a, 2, b, 2, alpha, beta, q, 1, u, 2, w, 1, rw, info);
        chkxer("Cggev ", infot, nout, lerr, ok);
        infot = 13;
        Cggev("V", "N", 2, a, 2, b, 2, alpha, beta, q, 2, u, 0, w, 1, rw, info);
        chkxer("Cggev ", infot, nout, lerr, ok);
        infot = 13;
        Cggev("V", "V", 2, a, 2, b, 2, alpha, beta, q, 2, u, 1, w, 1, rw, info);
        chkxer("Cggev ", infot, nout, lerr, ok);
        infot = 15;
        Cggev("V", "V", 1, a, 1, b, 1, alpha, beta, q, 1, u, 1, w, 1, rw, info);
        chkxer("Cggev ", infot, nout, lerr, ok);
        nt += 10;
        //
        //        Cggev3
        //
        infot = 1;
        Cggev3("/", "N", 1, a, 1, b, 1, alpha, beta, q, 1, u, 1, w, 1, rw, info);
        chkxer("Cggev3", infot, nout, lerr, ok);
        infot = 2;
        Cggev3("N", "/", 1, a, 1, b, 1, alpha, beta, q, 1, u, 1, w, 1, rw, info);
        chkxer("Cggev3", infot, nout, lerr, ok);
        infot = 3;
        Cggev3("V", "V", -1, a, 1, b, 1, alpha, beta, q, 1, u, 1, w, 1, rw, info);
        chkxer("Cggev3", infot, nout, lerr, ok);
        infot = 5;
        Cggev3("V", "V", 1, a, 0, b, 1, alpha, beta, q, 1, u, 1, w, 1, rw, info);
        chkxer("Cggev3", infot, nout, lerr, ok);
        infot = 7;
        Cggev3("V", "V", 1, a, 1, b, 0, alpha, beta, q, 1, u, 1, w, 1, rw, info);
        chkxer("Cggev3", infot, nout, lerr, ok);
        infot = 11;
        Cggev3("N", "V", 1, a, 1, b, 1, alpha, beta, q, 0, u, 1, w, 1, rw, info);
        chkxer("Cggev3", infot, nout, lerr, ok);
        infot = 11;
        Cggev3("V", "V", 2, a, 2, b, 2, alpha, beta, q, 1, u, 2, w, 1, rw, info);
        chkxer("Cggev3", infot, nout, lerr, ok);
        infot = 13;
        Cggev3("V", "N", 2, a, 2, b, 2, alpha, beta, q, 2, u, 0, w, 1, rw, info);
        chkxer("Cggev3", infot, nout, lerr, ok);
        infot = 13;
        Cggev3("V", "V", 2, a, 2, b, 2, alpha, beta, q, 2, u, 1, w, 1, rw, info);
        chkxer("Cggev3", infot, nout, lerr, ok);
        infot = 15;
        Cggev3("V", "V", 1, a, 1, b, 1, alpha, beta, q, 1, u, 1, w, 1, rw, info);
        chkxer("Cggev3", infot, nout, lerr, ok);
        nt += 10;
        //
        //        Cggevx
        //
        infot = 1;
        Cggevx("/", "N", "N", "N", 1, a, 1, b, 1, alpha, beta, q, 1, u, 1, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, rw, iw, bw, info);
        chkxer("Cggevx", infot, nout, lerr, ok);
        infot = 2;
        Cggevx("N", "/", "N", "N", 1, a, 1, b, 1, alpha, beta, q, 1, u, 1, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, rw, iw, bw, info);
        chkxer("Cggevx", infot, nout, lerr, ok);
        infot = 3;
        Cggevx("N", "N", "/", "N", 1, a, 1, b, 1, alpha, beta, q, 1, u, 1, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, rw, iw, bw, info);
        chkxer("Cggevx", infot, nout, lerr, ok);
        infot = 4;
        Cggevx("N", "N", "N", "/", 1, a, 1, b, 1, alpha, beta, q, 1, u, 1, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, rw, iw, bw, info);
        chkxer("Cggevx", infot, nout, lerr, ok);
        infot = 5;
        Cggevx("N", "N", "N", "N", -1, a, 1, b, 1, alpha, beta, q, 1, u, 1, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, rw, iw, bw, info);
        chkxer("Cggevx", infot, nout, lerr, ok);
        infot = 7;
        Cggevx("N", "N", "N", "N", 1, a, 0, b, 1, alpha, beta, q, 1, u, 1, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, rw, iw, bw, info);
        chkxer("Cggevx", infot, nout, lerr, ok);
        infot = 9;
        Cggevx("N", "N", "N", "N", 1, a, 1, b, 0, alpha, beta, q, 1, u, 1, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, rw, iw, bw, info);
        chkxer("Cggevx", infot, nout, lerr, ok);
        infot = 13;
        Cggevx("N", "N", "N", "N", 1, a, 1, b, 1, alpha, beta, q, 0, u, 1, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, rw, iw, bw, info);
        chkxer("Cggevx", infot, nout, lerr, ok);
        infot = 13;
        Cggevx("N", "V", "N", "N", 2, a, 2, b, 2, alpha, beta, q, 1, u, 2, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, rw, iw, bw, info);
        chkxer("Cggevx", infot, nout, lerr, ok);
        infot = 15;
        Cggevx("N", "N", "N", "N", 1, a, 1, b, 1, alpha, beta, q, 1, u, 0, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, rw, iw, bw, info);
        chkxer("Cggevx", infot, nout, lerr, ok);
        infot = 15;
        Cggevx("N", "N", "V", "N", 2, a, 2, b, 2, alpha, beta, q, 2, u, 1, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 1, rw, iw, bw, info);
        chkxer("Cggevx", infot, nout, lerr, ok);
        infot = 25;
        Cggevx("N", "N", "V", "N", 2, a, 2, b, 2, alpha, beta, q, 2, u, 2, ilo, ihi, ls, rs, anrm, bnrm, rce, rcv, w, 0, rw, iw, bw, info);
        chkxer("Cggevx", infot, nout, lerr, ok);
        nt += 12;
        //
        //        Ctgexc
        //
        infot = 3;
        Ctgexc(true, true, -1, a, 1, b, 1, q, 1, z, 1, ifst, ilst, info);
        chkxer("Ctgexc", infot, nout, lerr, ok);
        infot = 5;
        Ctgexc(true, true, 1, a, 0, b, 1, q, 1, z, 1, ifst, ilst, info);
        chkxer("Ctgexc", infot, nout, lerr, ok);
        infot = 7;
        Ctgexc(true, true, 1, a, 1, b, 0, q, 1, z, 1, ifst, ilst, info);
        chkxer("Ctgexc", infot, nout, lerr, ok);
        infot = 9;
        Ctgexc(false, true, 1, a, 1, b, 1, q, 0, z, 1, ifst, ilst, info);
        chkxer("Ctgexc", infot, nout, lerr, ok);
        infot = 9;
        Ctgexc(true, true, 1, a, 1, b, 1, q, 0, z, 1, ifst, ilst, info);
        chkxer("Ctgexc", infot, nout, lerr, ok);
        infot = 11;
        Ctgexc(true, false, 1, a, 1, b, 1, q, 1, z, 0, ifst, ilst, info);
        chkxer("Ctgexc", infot, nout, lerr, ok);
        infot = 11;
        Ctgexc(true, true, 1, a, 1, b, 1, q, 1, z, 0, ifst, ilst, info);
        chkxer("Ctgexc", infot, nout, lerr, ok);
        nt += 7;
        //
        //        Ctgsen
        //
        infot = 1;
        Ctgsen(-1, true, true, sel, 1, a, 1, b, 1, alpha, beta, q, 1, z, 1, m, tola, tolb, rcv, w, 1, iw, 1, info);
        chkxer("Ctgsen", infot, nout, lerr, ok);
        infot = 5;
        Ctgsen(1, true, true, sel, -1, a, 1, b, 1, alpha, beta, q, 1, z, 1, m, tola, tolb, rcv, w, 1, iw, 1, info);
        chkxer("Ctgsen", infot, nout, lerr, ok);
        infot = 7;
        Ctgsen(1, true, true, sel, 1, a, 0, b, 1, alpha, beta, q, 1, z, 1, m, tola, tolb, rcv, w, 1, iw, 1, info);
        chkxer("Ctgsen", infot, nout, lerr, ok);
        infot = 9;
        Ctgsen(1, true, true, sel, 1, a, 1, b, 0, alpha, beta, q, 1, z, 1, m, tola, tolb, rcv, w, 1, iw, 1, info);
        chkxer("Ctgsen", infot, nout, lerr, ok);
        infot = 13;
        Ctgsen(1, true, true, sel, 1, a, 1, b, 1, alpha, beta, q, 0, z, 1, m, tola, tolb, rcv, w, 1, iw, 1, info);
        chkxer("Ctgsen", infot, nout, lerr, ok);
        infot = 15;
        Ctgsen(1, true, true, sel, 1, a, 1, b, 1, alpha, beta, q, 1, z, 0, m, tola, tolb, rcv, w, 1, iw, 1, info);
        chkxer("Ctgsen", infot, nout, lerr, ok);
        infot = 21;
        Ctgsen(3, true, true, sel, 1, a, 1, b, 1, alpha, beta, q, 1, z, 1, m, tola, tolb, rcv, w, -5, iw, 1, info);
        chkxer("Ctgsen", infot, nout, lerr, ok);
        infot = 23;
        Ctgsen(0, true, true, sel, 1, a, 1, b, 1, alpha, beta, q, 1, z, 1, m, tola, tolb, rcv, w, 20, iw, 0, info);
        chkxer("Ctgsen", infot, nout, lerr, ok);
        infot = 23;
        Ctgsen(1, true, true, sel, 1, a, 1, b, 1, alpha, beta, q, 1, z, 1, m, tola, tolb, rcv, w, 20, iw, 0, info);
        chkxer("Ctgsen", infot, nout, lerr, ok);
        infot = 23;
        Ctgsen(5, true, true, sel, 1, a, 1, b, 1, alpha, beta, q, 1, z, 1, m, tola, tolb, rcv, w, 20, iw, 1, info);
        chkxer("Ctgsen", infot, nout, lerr, ok);
        nt += 11;
        //
        //        Ctgsna
        //
        infot = 1;
        Ctgsna("/", "A", sel, 1, a, 1, b, 1, q, 1, u, 1, r1, r2, 1, m, w, 1, iw, info);
        chkxer("Ctgsna", infot, nout, lerr, ok);
        infot = 2;
        Ctgsna("B", "/", sel, 1, a, 1, b, 1, q, 1, u, 1, r1, r2, 1, m, w, 1, iw, info);
        chkxer("Ctgsna", infot, nout, lerr, ok);
        infot = 4;
        Ctgsna("B", "A", sel, -1, a, 1, b, 1, q, 1, u, 1, r1, r2, 1, m, w, 1, iw, info);
        chkxer("Ctgsna", infot, nout, lerr, ok);
        infot = 6;
        Ctgsna("B", "A", sel, 1, a, 0, b, 1, q, 1, u, 1, r1, r2, 1, m, w, 1, iw, info);
        chkxer("Ctgsna", infot, nout, lerr, ok);
        infot = 8;
        Ctgsna("B", "A", sel, 1, a, 1, b, 0, q, 1, u, 1, r1, r2, 1, m, w, 1, iw, info);
        chkxer("Ctgsna", infot, nout, lerr, ok);
        infot = 10;
        Ctgsna("E", "A", sel, 1, a, 1, b, 1, q, 0, u, 1, r1, r2, 1, m, w, 1, iw, info);
        chkxer("Ctgsna", infot, nout, lerr, ok);
        infot = 12;
        Ctgsna("E", "A", sel, 1, a, 1, b, 1, q, 1, u, 0, r1, r2, 1, m, w, 1, iw, info);
        chkxer("Ctgsna", infot, nout, lerr, ok);
        infot = 15;
        Ctgsna("E", "A", sel, 1, a, 1, b, 1, q, 1, u, 1, r1, r2, 0, m, w, 1, iw, info);
        chkxer("Ctgsna", infot, nout, lerr, ok);
        infot = 18;
        Ctgsna("E", "A", sel, 1, a, 1, b, 1, q, 1, u, 1, r1, r2, 1, m, w, 0, iw, info);
        chkxer("Ctgsna", infot, nout, lerr, ok);
        nt += 9;
        //
        //        Ctgsyl
        //
        infot = 1;
        Ctgsyl("/", 0, 1, 1, a, 1, b, 1, q, 1, u, 1, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Ctgsyl", infot, nout, lerr, ok);
        infot = 2;
        Ctgsyl("N", -1, 1, 1, a, 1, b, 1, q, 1, u, 1, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Ctgsyl", infot, nout, lerr, ok);
        infot = 3;
        Ctgsyl("N", 0, 0, 1, a, 1, b, 1, q, 1, u, 1, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Ctgsyl", infot, nout, lerr, ok);
        infot = 4;
        Ctgsyl("N", 0, 1, 0, a, 1, b, 1, q, 1, u, 1, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Ctgsyl", infot, nout, lerr, ok);
        infot = 6;
        Ctgsyl("N", 0, 1, 1, a, 0, b, 1, q, 1, u, 1, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Ctgsyl", infot, nout, lerr, ok);
        infot = 8;
        Ctgsyl("N", 0, 1, 1, a, 1, b, 0, q, 1, u, 1, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Ctgsyl", infot, nout, lerr, ok);
        infot = 10;
        Ctgsyl("N", 0, 1, 1, a, 1, b, 1, q, 0, u, 1, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Ctgsyl", infot, nout, lerr, ok);
        infot = 12;
        Ctgsyl("N", 0, 1, 1, a, 1, b, 1, q, 1, u, 0, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Ctgsyl", infot, nout, lerr, ok);
        infot = 14;
        Ctgsyl("N", 0, 1, 1, a, 1, b, 1, q, 1, u, 1, v, 0, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Ctgsyl", infot, nout, lerr, ok);
        infot = 16;
        Ctgsyl("N", 0, 1, 1, a, 1, b, 1, q, 1, u, 1, v, 1, z, 0, scale, dif, w, 1, iw, info);
        chkxer("Ctgsyl", infot, nout, lerr, ok);
        infot = 20;
        Ctgsyl("N", 1, 1, 1, a, 1, b, 1, q, 1, u, 1, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Ctgsyl", infot, nout, lerr, ok);
        infot = 20;
        Ctgsyl("N", 2, 1, 1, a, 1, b, 1, q, 1, u, 1, v, 1, z, 1, scale, dif, w, 1, iw, info);
        chkxer("Ctgsyl", infot, nout, lerr, ok);
        nt += 12;
    }
    //
    //     Print a summary line.
    //
    if (ok) {
        write(nout, "(1x,a3,' routines passed the tests of the error exits (',i3,"
                    "' tests done)')"),
            path, nt;
    } else {
        write(nout, "(' *** ',a3,' routines failed the tests of the error ','exits ***')"), path;
    }
    //
    //     End of Cerrgg
    //
}
