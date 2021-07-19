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
#include <mplapack_debug.h>

void Cerrvx(const char *path, INTEGER const nunit) {
    common cmn;
    common_write write(cmn);
    //
    nout = nunit;
    char c2[2];
    c2[0] = path[1];
    c2[1] = path[2];
    //
    //     Set the variables to innocuous values.
    //
    INTEGER j = 0;
    const INTEGER nmax = 4;
    INTEGER i = 0;
    COMPLEX a[nmax * nmax];
    COMPLEX af[nmax * nmax];
    INTEGER lda = nmax;
    INTEGER ldaf = nmax;
    COMPLEX b[nmax];
    COMPLEX e[nmax];
    REAL r1[nmax];
    REAL r2[nmax];
    COMPLEX w[2 * nmax];
    COMPLEX x[nmax];
    REAL c[nmax];
    REAL r[nmax];
    INTEGER ip[nmax];
    for (j = 1; j <= nmax; j = j + 1) {
        for (i = 1; i <= nmax; i = i + 1) {
            a[(i - 1) + (j - 1) * lda] = COMPLEX(1.0 / castREAL(i + j), -1.0 / castREAL(i + j));
            af[(i - 1) + (j - 1) * ldaf] = COMPLEX(1.0 / castREAL(i + j), -1.0 / castREAL(i + j));
        }
        b[j - 1] = 0.0;
        e[j - 1] = 0.0;
        r1[j - 1] = 0.0;
        r2[j - 1] = 0.0;
        w[j - 1] = 0.0;
        x[j - 1] = 0.0;
        c[j - 1] = 0.0;
        r[j - 1] = 0.0;
        ip[j - 1] = j;
    }
    char eq=' ';
    ok = true;
    //
    INTEGER info = 0;
    REAL rcond = 0.0;
    REAL rw[nmax];
    REAL rf[nmax];
    if (Mlsamen(2, c2, "GE")) {
        //
        //        Cgesv
        //
        strncpy(srnamt, "Cgesv", srnamt_len);
        infot = 1;
        Cgesv(-1, 0, a, 1, ip, b, 1, info);
        chkxer("Cgesv", infot, nout, lerr, ok);
        infot = 2;
        Cgesv(0, -1, a, 1, ip, b, 1, info);
        chkxer("Cgesv", infot, nout, lerr, ok);
        infot = 4;
        Cgesv(2, 1, a, 1, ip, b, 2, info);
        chkxer("Cgesv", infot, nout, lerr, ok);
        infot = 7;
        Cgesv(2, 1, a, 2, ip, b, 1, info);
        chkxer("Cgesv", infot, nout, lerr, ok);
        //
        //        Cgesvx
        //
        strncpy(srnamt, "Cgesvx", srnamt_len);
        infot = 1;
        Cgesvx("/", "N", 0, 0, a, 1, af, 1, ip, &eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgesvx", infot, nout, lerr, ok);
        infot = 2;
        Cgesvx("N", "/", 0, 0, a, 1, af, 1, ip, &eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgesvx", infot, nout, lerr, ok);
        infot = 3;
        Cgesvx("N", "N", -1, 0, a, 1, af, 1, ip, &eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgesvx", infot, nout, lerr, ok);
        infot = 4;
        Cgesvx("N", "N", 0, -1, a, 1, af, 1, ip, &eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgesvx", infot, nout, lerr, ok);
        infot = 6;
        Cgesvx("N", "N", 2, 1, a, 1, af, 2, ip, &eq, r, c, b, 2, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("Cgesvx", infot, nout, lerr, ok);
        infot = 8;
        Cgesvx("N", "N", 2, 1, a, 2, af, 1, ip, &eq, r, c, b, 2, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("Cgesvx", infot, nout, lerr, ok);
        infot = 10;
        eq = '/';
        Cgesvx("F", "N", 0, 0, a, 1, af, 1, ip, &eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgesvx", infot, nout, lerr, ok);
        infot = 11;
        eq = 'R';
        Cgesvx("F", "N", 1, 0, a, 1, af, 1, ip, &eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgesvx", infot, nout, lerr, ok);
        infot = 12;
        eq = 'C';
        Cgesvx("F", "N", 1, 0, a, 1, af, 1, ip, &eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgesvx", infot, nout, lerr, ok);
        infot = 14;
        Cgesvx("N", "N", 2, 1, a, 2, af, 2, ip, &eq, r, c, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("Cgesvx", infot, nout, lerr, ok);
        infot = 16;
        Cgesvx("N", "N", 2, 1, a, 2, af, 2, ip, &eq, r, c, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgesvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "GB")) {
        //
        //        Cgbsv
        //
        strncpy(srnamt, "Cgbsv", srnamt_len);
        infot = 1;
        Cgbsv(-1, 0, 0, 0, a, 1, ip, b, 1, info);
        chkxer("Cgbsv ", infot, nout, lerr, ok);
        infot = 2;
        Cgbsv(1, -1, 0, 0, a, 1, ip, b, 1, info);
        chkxer("Cgbsv ", infot, nout, lerr, ok);
        infot = 3;
        Cgbsv(1, 0, -1, 0, a, 1, ip, b, 1, info);
        chkxer("Cgbsv ", infot, nout, lerr, ok);
        infot = 4;
        Cgbsv(0, 0, 0, -1, a, 1, ip, b, 1, info);
        chkxer("Cgbsv ", infot, nout, lerr, ok);
        infot = 6;
        Cgbsv(1, 1, 1, 0, a, 3, ip, b, 1, info);
        chkxer("Cgbsv ", infot, nout, lerr, ok);
        infot = 9;
        Cgbsv(2, 0, 0, 0, a, 1, ip, b, 1, info);
        chkxer("Cgbsv ", infot, nout, lerr, ok);
        //
        //        Cgbsvx
        //
        strncpy(srnamt, "Cgbsvx", srnamt_len);
        infot = 1;
        Cgbsvx("/", "N", 0, 0, 0, 0, a, 1, af, 1, ip, &eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgbsvx", infot, nout, lerr, ok);
        infot = 2;
        Cgbsvx("N", "/", 0, 0, 0, 0, a, 1, af, 1, ip, &eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgbsvx", infot, nout, lerr, ok);
        infot = 3;
        Cgbsvx("N", "N", -1, 0, 0, 0, a, 1, af, 1, ip, &eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgbsvx", infot, nout, lerr, ok);
        infot = 4;
        Cgbsvx("N", "N", 1, -1, 0, 0, a, 1, af, 1, ip, &eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgbsvx", infot, nout, lerr, ok);
        infot = 5;
        Cgbsvx("N", "N", 1, 0, -1, 0, a, 1, af, 1, ip, &eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgbsvx", infot, nout, lerr, ok);
        infot = 6;
        Cgbsvx("N", "N", 0, 0, 0, -1, a, 1, af, 1, ip, &eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgbsvx", infot, nout, lerr, ok);
        infot = 8;
        Cgbsvx("N", "N", 1, 1, 1, 0, a, 2, af, 4, ip, &eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgbsvx", infot, nout, lerr, ok);
        infot = 10;
        Cgbsvx("N", "N", 1, 1, 1, 0, a, 3, af, 3, ip, &eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgbsvx", infot, nout, lerr, ok);
        infot = 12;
        eq = '/';
        Cgbsvx("F", "N", 0, 0, 0, 0, a, 1, af, 1, ip, &eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgbsvx", infot, nout, lerr, ok);
        infot = 13;
        eq = 'R';
        Cgbsvx("F", "N", 1, 0, 0, 0, a, 1, af, 1, ip, &eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgbsvx", infot, nout, lerr, ok);
        infot = 14;
        eq = 'C';
        Cgbsvx("F", "N", 1, 0, 0, 0, a, 1, af, 1, ip, &eq, r, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgbsvx", infot, nout, lerr, ok);
        infot = 16;
        Cgbsvx("N", "N", 2, 0, 0, 0, a, 1, af, 1, ip, &eq, r, c, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("Cgbsvx", infot, nout, lerr, ok);
        infot = 18;
        Cgbsvx("N", "N", 2, 0, 0, 0, a, 1, af, 1, ip, &eq, r, c, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgbsvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "GT")) {
        //
        //        Cgtsv
        //
        strncpy(srnamt, "Cgtsv", srnamt_len);
        infot = 1;
        Cgtsv(-1, 0, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], b, 1, info);
        chkxer("Cgtsv ", infot, nout, lerr, ok);
        infot = 2;
        Cgtsv(0, -1, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], b, 1, info);
        chkxer("Cgtsv ", infot, nout, lerr, ok);
        infot = 7;
        Cgtsv(2, 0, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], b, 1, info);
        chkxer("Cgtsv ", infot, nout, lerr, ok);
        //
        //        Cgtsvx
        //
        strncpy(srnamt, "Cgtsvx", srnamt_len);
        infot = 1;
        Cgtsvx("/", "N", 0, 0, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &af[(1 - 1)], &af[(2 - 1) * ldaf], &af[(3 - 1) * ldaf], &af[(4 - 1) * ldaf], ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgtsvx", infot, nout, lerr, ok);
        infot = 2;
        Cgtsvx("N", "/", 0, 0, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &af[(1 - 1)], &af[(2 - 1) * ldaf], &af[(3 - 1) * ldaf], &af[(4 - 1) * ldaf], ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgtsvx", infot, nout, lerr, ok);
        infot = 3;
        Cgtsvx("N", "N", -1, 0, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &af[(1 - 1)], &af[(2 - 1) * ldaf], &af[(3 - 1) * ldaf], &af[(4 - 1) * ldaf], ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgtsvx", infot, nout, lerr, ok);
        infot = 4;
        Cgtsvx("N", "N", 0, -1, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &af[(1 - 1)], &af[(2 - 1) * ldaf], &af[(3 - 1) * ldaf], &af[(4 - 1) * ldaf], ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgtsvx", infot, nout, lerr, ok);
        infot = 14;
        Cgtsvx("N", "N", 2, 0, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &af[(1 - 1)], &af[(2 - 1) * ldaf], &af[(3 - 1) * ldaf], &af[(4 - 1) * ldaf], ip, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("Cgtsvx", infot, nout, lerr, ok);
        infot = 16;
        Cgtsvx("N", "N", 2, 0, &a[(1 - 1) + (1 - 1) * lda], &a[(1 - 1) + (2 - 1) * lda], &a[(1 - 1) + (3 - 1) * lda], &af[(1 - 1)], &af[(2 - 1) * ldaf], &af[(3 - 1) * ldaf], &af[(4 - 1) * ldaf], ip, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cgtsvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "PO")) {
        //
        //        Cposv
        //
        strncpy(srnamt, "Cposv", srnamt_len);
        infot = 1;
        Cposv("/", 0, 0, a, 1, b, 1, info);
        chkxer("Cposv", infot, nout, lerr, ok);
        infot = 2;
        Cposv("U", -1, 0, a, 1, b, 1, info);
        chkxer("Cposv", infot, nout, lerr, ok);
        infot = 3;
        Cposv("U", 0, -1, a, 1, b, 1, info);
        chkxer("Cposv", infot, nout, lerr, ok);
        infot = 5;
        Cposv("U", 2, 0, a, 1, b, 2, info);
        chkxer("Cposv", infot, nout, lerr, ok);
        infot = 7;
        Cposv("U", 2, 0, a, 2, b, 1, info);
        chkxer("Cposv", infot, nout, lerr, ok);
        //
        //        Cposvx
        //
        strncpy(srnamt, "Cposvx", srnamt_len);
        infot = 1;
        Cposvx("/", "U", 0, 0, a, 1, af, 1, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cposvx", infot, nout, lerr, ok);
        infot = 2;
        Cposvx("N", "/", 0, 0, a, 1, af, 1, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cposvx", infot, nout, lerr, ok);
        infot = 3;
        Cposvx("N", "U", -1, 0, a, 1, af, 1, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cposvx", infot, nout, lerr, ok);
        infot = 4;
        Cposvx("N", "U", 0, -1, a, 1, af, 1, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cposvx", infot, nout, lerr, ok);
        infot = 6;
        Cposvx("N", "U", 2, 0, a, 1, af, 2, &eq, c, b, 2, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("Cposvx", infot, nout, lerr, ok);
        infot = 8;
        Cposvx("N", "U", 2, 0, a, 2, af, 1, &eq, c, b, 2, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("Cposvx", infot, nout, lerr, ok);
        infot = 9;
        eq = '/';
        Cposvx("F", "U", 0, 0, a, 1, af, 1, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cposvx", infot, nout, lerr, ok);
        infot = 10;
        eq = 'Y';
        Cposvx("F", "U", 1, 0, a, 1, af, 1, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cposvx", infot, nout, lerr, ok);
        infot = 12;
        Cposvx("N", "U", 2, 0, a, 2, af, 2, &eq, c, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("Cposvx", infot, nout, lerr, ok);
        infot = 14;
        Cposvx("N", "U", 2, 0, a, 2, af, 2, &eq, c, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cposvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "PP")) {
        //
        //        Cppsv
        //
        strncpy(srnamt, "Cppsv", srnamt_len);
        infot = 1;
        Cppsv("/", 0, 0, a, b, 1, info);
        chkxer("Cppsv ", infot, nout, lerr, ok);
        infot = 2;
        Cppsv("U", -1, 0, a, b, 1, info);
        chkxer("Cppsv ", infot, nout, lerr, ok);
        infot = 3;
        Cppsv("U", 0, -1, a, b, 1, info);
        chkxer("Cppsv ", infot, nout, lerr, ok);
        infot = 6;
        Cppsv("U", 2, 0, a, b, 1, info);
        chkxer("Cppsv ", infot, nout, lerr, ok);
        //
        //        Cppsvx
        //
        strncpy(srnamt, "Cppsvx", srnamt_len);
        infot = 1;
        Cppsvx("/", "U", 0, 0, a, af, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cppsvx", infot, nout, lerr, ok);
        infot = 2;
        Cppsvx("N", "/", 0, 0, a, af, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cppsvx", infot, nout, lerr, ok);
        infot = 3;
        Cppsvx("N", "U", -1, 0, a, af, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cppsvx", infot, nout, lerr, ok);
        infot = 4;
        Cppsvx("N", "U", 0, -1, a, af, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cppsvx", infot, nout, lerr, ok);
        infot = 7;
        eq = '/';
        Cppsvx("F", "U", 0, 0, a, af, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cppsvx", infot, nout, lerr, ok);
        infot = 8;
        eq = 'Y';
        Cppsvx("F", "U", 1, 0, a, af, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cppsvx", infot, nout, lerr, ok);
        infot = 10;
        Cppsvx("N", "U", 2, 0, a, af, &eq, c, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("Cppsvx", infot, nout, lerr, ok);
        infot = 12;
        Cppsvx("N", "U", 2, 0, a, af, &eq, c, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cppsvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "PB")) {
        //
        //        Cpbsv
        //
        strncpy(srnamt, "Cpbsv", srnamt_len);
        infot = 1;
        Cpbsv("/", 0, 0, 0, a, 1, b, 1, info);
        chkxer("Cpbsv ", infot, nout, lerr, ok);
        infot = 2;
        Cpbsv("U", -1, 0, 0, a, 1, b, 1, info);
        chkxer("Cpbsv ", infot, nout, lerr, ok);
        infot = 3;
        Cpbsv("U", 1, -1, 0, a, 1, b, 1, info);
        chkxer("Cpbsv ", infot, nout, lerr, ok);
        infot = 4;
        Cpbsv("U", 0, 0, -1, a, 1, b, 1, info);
        chkxer("Cpbsv ", infot, nout, lerr, ok);
        infot = 6;
        Cpbsv("U", 1, 1, 0, a, 1, b, 2, info);
        chkxer("Cpbsv ", infot, nout, lerr, ok);
        infot = 8;
        Cpbsv("U", 2, 0, 0, a, 1, b, 1, info);
        chkxer("Cpbsv ", infot, nout, lerr, ok);
        //
        //        Cpbsvx
        //
        strncpy(srnamt, "Cpbsvx", srnamt_len);
        infot = 1;
        Cpbsvx("/", "U", 0, 0, 0, a, 1, af, 1, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cpbsvx", infot, nout, lerr, ok);
        infot = 2;
        Cpbsvx("N", "/", 0, 0, 0, a, 1, af, 1, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cpbsvx", infot, nout, lerr, ok);
        infot = 3;
        Cpbsvx("N", "U", -1, 0, 0, a, 1, af, 1, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cpbsvx", infot, nout, lerr, ok);
        infot = 4;
        Cpbsvx("N", "U", 1, -1, 0, a, 1, af, 1, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cpbsvx", infot, nout, lerr, ok);
        infot = 5;
        Cpbsvx("N", "U", 0, 0, -1, a, 1, af, 1, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cpbsvx", infot, nout, lerr, ok);
        infot = 7;
        Cpbsvx("N", "U", 1, 1, 0, a, 1, af, 2, &eq, c, b, 2, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("Cpbsvx", infot, nout, lerr, ok);
        infot = 9;
        Cpbsvx("N", "U", 1, 1, 0, a, 2, af, 1, &eq, c, b, 2, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("Cpbsvx", infot, nout, lerr, ok);
        infot = 10;
        eq = '/';
        Cpbsvx("F", "U", 0, 0, 0, a, 1, af, 1, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cpbsvx", infot, nout, lerr, ok);
        infot = 11;
        eq = 'Y';
        Cpbsvx("F", "U", 1, 0, 0, a, 1, af, 1, &eq, c, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cpbsvx", infot, nout, lerr, ok);
        infot = 13;
        Cpbsvx("N", "U", 2, 0, 0, a, 1, af, 1, &eq, c, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("Cpbsvx", infot, nout, lerr, ok);
        infot = 15;
        Cpbsvx("N", "U", 2, 0, 0, a, 1, af, 1, &eq, c, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cpbsvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "PT")) {
        //
        //        Cptsv
        //
        strncpy(srnamt, "Cptsv", srnamt_len);
        infot = 1;
        Cptsv(-1, 0, r, &a[(1 - 1) + (1 - 1) * lda], b, 1, info);
        chkxer("Cptsv ", infot, nout, lerr, ok);
        infot = 2;
        Cptsv(0, -1, r, &a[(1 - 1) + (1 - 1) * lda], b, 1, info);
        chkxer("Cptsv ", infot, nout, lerr, ok);
        infot = 6;
        Cptsv(2, 0, r, &a[(1 - 1) + (1 - 1) * lda], b, 1, info);
        chkxer("Cptsv ", infot, nout, lerr, ok);
        //
        //        Cptsvx
        //
        strncpy(srnamt, "Cptsvx", srnamt_len);
        infot = 1;
        Cptsvx("/", 0, 0, r, &a[(1 - 1) + (1 - 1) * lda], rf, &af[(1 - 1)], b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cptsvx", infot, nout, lerr, ok);
        infot = 2;
        Cptsvx("N", -1, 0, r, &a[(1 - 1) + (1 - 1) * lda], rf, &af[(1 - 1)], b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cptsvx", infot, nout, lerr, ok);
        infot = 3;
        Cptsvx("N", 0, -1, r, &a[(1 - 1) + (1 - 1) * lda], rf, &af[(1 - 1)], b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cptsvx", infot, nout, lerr, ok);
        infot = 9;
        Cptsvx("N", 2, 0, r, &a[(1 - 1) + (1 - 1) * lda], rf, &af[(1 - 1)], b, 1, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("Cptsvx", infot, nout, lerr, ok);
        infot = 11;
        Cptsvx("N", 2, 0, r, &a[(1 - 1) + (1 - 1) * lda], rf, &af[(1 - 1)], b, 2, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cptsvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "HE")) {
        //
        //        Chesv
        //
        strncpy(srnamt, "Chesv", srnamt_len);
        infot = 1;
        Chesv("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Chesv ", infot, nout, lerr, ok);
        infot = 2;
        Chesv("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Chesv ", infot, nout, lerr, ok);
        infot = 3;
        Chesv("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        chkxer("Chesv ", infot, nout, lerr, ok);
        infot = 5;
        Chesv("U", 2, 0, a, 1, ip, b, 2, w, 1, info);
        chkxer("Chesv ", infot, nout, lerr, ok);
        infot = 8;
        Chesv("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        chkxer("Chesv ", infot, nout, lerr, ok);
        infot = 10;
        Chesv("U", 0, 0, a, 1, ip, b, 1, w, 0, info);
        chkxer("Chesv ", infot, nout, lerr, ok);
        infot = 10;
        Chesv("U", 0, 0, a, 1, ip, b, 1, w, -2, info);
        chkxer("Chesv ", infot, nout, lerr, ok);
        //
        //        Chesvx
        //
        strncpy(srnamt, "Chesvx", srnamt_len);
        infot = 1;
        Chesvx("/", "U", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        chkxer("Chesvx", infot, nout, lerr, ok);
        infot = 2;
        Chesvx("N", "/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        chkxer("Chesvx", infot, nout, lerr, ok);
        infot = 3;
        Chesvx("N", "U", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        chkxer("Chesvx", infot, nout, lerr, ok);
        infot = 4;
        Chesvx("N", "U", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        chkxer("Chesvx", infot, nout, lerr, ok);
        infot = 6;
        Chesvx("N", "U", 2, 0, a, 1, af, 2, ip, b, 2, x, 2, rcond, r1, r2, w, 4, rw, info);
        chkxer("Chesvx", infot, nout, lerr, ok);
        infot = 8;
        Chesvx("N", "U", 2, 0, a, 2, af, 1, ip, b, 2, x, 2, rcond, r1, r2, w, 4, rw, info);
        chkxer("Chesvx", infot, nout, lerr, ok);
        infot = 11;
        Chesvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 1, x, 2, rcond, r1, r2, w, 4, rw, info);
        chkxer("Chesvx", infot, nout, lerr, ok);
        infot = 13;
        Chesvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 2, x, 1, rcond, r1, r2, w, 4, rw, info);
        chkxer("Chesvx", infot, nout, lerr, ok);
        infot = 18;
        Chesvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 2, x, 2, rcond, r1, r2, w, 3, rw, info);
        chkxer("Chesvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "HR")) {
        //
        //        Chesv_rook
        //
        strncpy(srnamt, "Chesv_rook", srnamt_len);
        infot = 1;
        Chesv_rook("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Chesv_rook", infot, nout, lerr, ok);
        infot = 2;
        Chesv_rook("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Chesv_rook", infot, nout, lerr, ok);
        infot = 3;
        Chesv_rook("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        chkxer("Chesv_rook", infot, nout, lerr, ok);
        infot = 8;
        Chesv_rook("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        chkxer("Chesv_rook", infot, nout, lerr, ok);
        infot = 10;
        Chesv_rook("U", 0, 0, a, 1, ip, b, 1, w, 0, info);
        chkxer("Chesv_rook", infot, nout, lerr, ok);
        infot = 10;
        Chesv_rook("U", 0, 0, a, 1, ip, b, 1, w, -2, info);
        chkxer("Chesv_rook", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "HK")) {
        //
        //        Csysv_rk
        //
        //        Test error exits of the driver that uses factorization
        //        of a Hermitian indefinite matrix with rook
        //        (bounded Bunch-Kaufman) pivoting with the new storage
        //        format for factors L ( or U) and D.
        //
        //        L (or U) is stored in A, diagonal of D is stored on the
        //        diagonal of A, subdiagonal of D is stored in a separate array E.
        //
        strncpy(srnamt, "Chesv_rk", srnamt_len);
        infot = 1;
        Chesv_rk("/", 0, 0, a, 1, e, ip, b, 1, w, 1, info);
        chkxer("Chesv_rk", infot, nout, lerr, ok);
        infot = 2;
        Chesv_rk("U", -1, 0, a, 1, e, ip, b, 1, w, 1, info);
        chkxer("Chesv_rk", infot, nout, lerr, ok);
        infot = 3;
        Chesv_rk("U", 0, -1, a, 1, e, ip, b, 1, w, 1, info);
        chkxer("Chesv_rk", infot, nout, lerr, ok);
        infot = 5;
        Chesv_rk("U", 2, 0, a, 1, e, ip, b, 2, w, 1, info);
        chkxer("Chesv_rk", infot, nout, lerr, ok);
        infot = 9;
        Chesv_rk("U", 2, 0, a, 2, e, ip, b, 1, w, 1, info);
        chkxer("Chesv_rk", infot, nout, lerr, ok);
        infot = 11;
        Chesv_rk("U", 0, 0, a, 1, e, ip, b, 1, w, 0, info);
        chkxer("Chesv_rk", infot, nout, lerr, ok);
        infot = 11;
        Chesv_rk("U", 0, 0, a, 1, e, ip, b, 1, w, -2, info);
        chkxer("Chesv_rk", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "HA")) {
        //
        //        Chesv_aa
        //
        strncpy(srnamt, "Chesv_aa", srnamt_len);
        infot = 1;
        Chesv_aa("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Chesv_aa", infot, nout, lerr, ok);
        infot = 2;
        Chesv_aa("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Chesv_aa", infot, nout, lerr, ok);
        infot = 3;
        Chesv_aa("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        chkxer("Chesv_aa", infot, nout, lerr, ok);
        infot = 8;
        Chesv_aa("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        chkxer("Chesv_aa", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "H2")) {
        //
        //        Chesv_aa_2stage
        //
        strncpy(srnamt, "Chesv_aa_2stage", srnamt_len);
        infot = 1;
        Chesv_aa_2stage("/", 0, 0, a, 1, a, 1, ip, ip, b, 1, w, 1, info);
        chkxer("Chesv_aa_2stage", infot, nout, lerr, ok);
        infot = 2;
        Chesv_aa_2stage("U", -1, 0, a, 1, a, 1, ip, ip, b, 1, w, 1, info);
        chkxer("Chesv_aa_2stage", infot, nout, lerr, ok);
        infot = 3;
        Chesv_aa_2stage("U", 0, -1, a, 1, a, 1, ip, ip, b, 1, w, 1, info);
        chkxer("Chesv_aa_2stage", infot, nout, lerr, ok);
        infot = 5;
        Chesv_aa_2stage("U", 2, 1, a, 1, a, 1, ip, ip, b, 1, w, 1, info);
        chkxer("Chesv_aa_2stage", infot, nout, lerr, ok);
        infot = 11;
        Chesv_aa_2stage("U", 2, 1, a, 2, a, 8, ip, ip, b, 1, w, 1, info);
        chkxer("Chesv_aa_2stage", infot, nout, lerr, ok);
        infot = 7;
        Chesv_aa_2stage("U", 2, 1, a, 2, a, 1, ip, ip, b, 2, w, 1, info);
        chkxer("Chesv_aa_2stage", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "S2")) {
        //
        //        Csysv_aa_2stage
        //
        strncpy(srnamt, "Csysv_aa_2stage", srnamt_len);
        infot = 1;
        Csysv_aa_2stage("/", 0, 0, a, 1, a, 1, ip, ip, b, 1, w, 1, info);
        chkxer("Csysv_aa_2stage", infot, nout, lerr, ok);
        infot = 2;
        Csysv_aa_2stage("U", -1, 0, a, 1, a, 1, ip, ip, b, 1, w, 1, info);
        chkxer("Csysv_aa_2stage", infot, nout, lerr, ok);
        infot = 3;
        Csysv_aa_2stage("U", 0, -1, a, 1, a, 1, ip, ip, b, 1, w, 1, info);
        chkxer("Csysv_aa_2stage", infot, nout, lerr, ok);
        infot = 5;
        Csysv_aa_2stage("U", 2, 1, a, 1, a, 1, ip, ip, b, 1, w, 1, info);
        chkxer("Csysv_aa_2stage", infot, nout, lerr, ok);
        infot = 11;
        Csysv_aa_2stage("U", 2, 1, a, 2, a, 8, ip, ip, b, 1, w, 1, info);
        chkxer("Csysv_aa_2stage", infot, nout, lerr, ok);
        infot = 7;
        Csysv_aa_2stage("U", 2, 1, a, 2, a, 1, ip, ip, b, 2, w, 1, info);
        chkxer("Csysv_aa_2stage", infot, nout, lerr, ok);
        //*
    } else if (Mlsamen(2, c2, "HP")) {
        //
        //        Chpsv
        //
        strncpy(srnamt, "Chpsv", srnamt_len);
        infot = 1;
        Chpsv("/", 0, 0, a, ip, b, 1, info);
        chkxer("Chpsv", infot, nout, lerr, ok);
        infot = 2;
        Chpsv("U", -1, 0, a, ip, b, 1, info);
        chkxer("Chpsv", infot, nout, lerr, ok);
        infot = 3;
        Chpsv("U", 0, -1, a, ip, b, 1, info);
        chkxer("Chpsv", infot, nout, lerr, ok);
        infot = 7;
        Chpsv("U", 2, 0, a, ip, b, 1, info);
        chkxer("Chpsv", infot, nout, lerr, ok);
        //
        //        Chpsvx
        //
        strncpy(srnamt, "Chpsvx", srnamt_len);
        infot = 1;
        Chpsvx("/", "U", 0, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Chpsvx", infot, nout, lerr, ok);
        infot = 2;
        Chpsvx("N", "/", 0, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Chpsvx", infot, nout, lerr, ok);
        infot = 3;
        Chpsvx("N", "U", -1, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Chpsvx", infot, nout, lerr, ok);
        infot = 4;
        Chpsvx("N", "U", 0, -1, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Chpsvx", infot, nout, lerr, ok);
        infot = 9;
        Chpsvx("N", "U", 2, 0, a, af, ip, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("Chpsvx", infot, nout, lerr, ok);
        infot = 11;
        Chpsvx("N", "U", 2, 0, a, af, ip, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Chpsvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "SY")) {
        //
        //        Csysv
        //
        strncpy(srnamt, "Csysv", srnamt_len);
        infot = 1;
        Csysv("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Csysv", infot, nout, lerr, ok);
        infot = 2;
        Csysv("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Csysv", infot, nout, lerr, ok);
        infot = 3;
        Csysv("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        chkxer("Csysv", infot, nout, lerr, ok);
        infot = 8;
        Csysv("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        chkxer("Csysv", infot, nout, lerr, ok);
        infot = 10;
        Csysv("U", 0, 0, a, 1, ip, b, 1, w, 0, info);
        chkxer("Csysv", infot, nout, lerr, ok);
        infot = 10;
        Csysv("U", 0, 0, a, 1, ip, b, 1, w, -2, info);
        chkxer("Csysv", infot, nout, lerr, ok);
        //
        //        Csysvx
        //
        strncpy(srnamt, "Csysvx", srnamt_len);
        infot = 1;
        Csysvx("/", "U", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        chkxer("Csysvx", infot, nout, lerr, ok);
        infot = 2;
        Csysvx("N", "/", 0, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        chkxer("Csysvx", infot, nout, lerr, ok);
        infot = 3;
        Csysvx("N", "U", -1, 0, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        chkxer("Csysvx", infot, nout, lerr, ok);
        infot = 4;
        Csysvx("N", "U", 0, -1, a, 1, af, 1, ip, b, 1, x, 1, rcond, r1, r2, w, 1, rw, info);
        chkxer("Csysvx", infot, nout, lerr, ok);
        infot = 6;
        Csysvx("N", "U", 2, 0, a, 1, af, 2, ip, b, 2, x, 2, rcond, r1, r2, w, 4, rw, info);
        chkxer("Csysvx", infot, nout, lerr, ok);
        infot = 8;
        Csysvx("N", "U", 2, 0, a, 2, af, 1, ip, b, 2, x, 2, rcond, r1, r2, w, 4, rw, info);
        chkxer("Csysvx", infot, nout, lerr, ok);
        infot = 11;
        Csysvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 1, x, 2, rcond, r1, r2, w, 4, rw, info);
        chkxer("Csysvx", infot, nout, lerr, ok);
        infot = 13;
        Csysvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 2, x, 1, rcond, r1, r2, w, 4, rw, info);
        chkxer("Csysvx", infot, nout, lerr, ok);
        infot = 18;
        Csysvx("N", "U", 2, 0, a, 2, af, 2, ip, b, 2, x, 2, rcond, r1, r2, w, 3, rw, info);
        chkxer("Csysvx", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "SR")) {
        //
        //        Csysv_rook
        //
        strncpy(srnamt, "Csysv_rook", srnamt_len);
        infot = 1;
        Csysv_rook("/", 0, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Csysv_rook", infot, nout, lerr, ok);
        infot = 2;
        Csysv_rook("U", -1, 0, a, 1, ip, b, 1, w, 1, info);
        chkxer("Csysv_rook", infot, nout, lerr, ok);
        infot = 3;
        Csysv_rook("U", 0, -1, a, 1, ip, b, 1, w, 1, info);
        chkxer("Csysv_rook", infot, nout, lerr, ok);
        infot = 8;
        Csysv_rook("U", 2, 0, a, 2, ip, b, 1, w, 1, info);
        chkxer("Csysv_rook", infot, nout, lerr, ok);
        infot = 10;
        Csysv_rook("U", 0, 0, a, 1, ip, b, 1, w, 0, info);
        chkxer("Csysv_rook", infot, nout, lerr, ok);
        infot = 10;
        Csysv_rook("U", 0, 0, a, 1, ip, b, 1, w, -2, info);
        //
    } else if (Mlsamen(2, c2, "SK")) {
        //
        //        Csysv_rk
        //
        //        Test error exits of the driver that uses factorization
        //        of a symmetric indefinite matrix with rook
        //        (bounded Bunch-Kaufman) pivoting with the new storage
        //        format for factors L ( or U) and D.
        //
        //        L (or U) is stored in A, diagonal of D is stored on the
        //        diagonal of A, subdiagonal of D is stored in a separate array E.
        //
        strncpy(srnamt, "Csysv_rk", srnamt_len);
        infot = 1;
        Csysv_rk("/", 0, 0, a, 1, e, ip, b, 1, w, 1, info);
        chkxer("Csysv_rk", infot, nout, lerr, ok);
        infot = 2;
        Csysv_rk("U", -1, 0, a, 1, e, ip, b, 1, w, 1, info);
        chkxer("Csysv_rk", infot, nout, lerr, ok);
        infot = 3;
        Csysv_rk("U", 0, -1, a, 1, e, ip, b, 1, w, 1, info);
        chkxer("Csysv_rk", infot, nout, lerr, ok);
        infot = 5;
        Csysv_rk("U", 2, 0, a, 1, e, ip, b, 2, w, 1, info);
        chkxer("Csysv_rk", infot, nout, lerr, ok);
        infot = 9;
        Csysv_rk("U", 2, 0, a, 2, e, ip, b, 1, w, 1, info);
        chkxer("Csysv_rk", infot, nout, lerr, ok);
        infot = 11;
        Csysv_rk("U", 0, 0, a, 1, e, ip, b, 1, w, 0, info);
        chkxer("Csysv_rk", infot, nout, lerr, ok);
        infot = 11;
        Csysv_rk("U", 0, 0, a, 1, e, ip, b, 1, w, -2, info);
        chkxer("Csysv_rk", infot, nout, lerr, ok);
        //
    } else if (Mlsamen(2, c2, "SP")) {
        //
        //        Cspsv
        //
        strncpy(srnamt, "Cspsv", srnamt_len);
        infot = 1;
        Cspsv("/", 0, 0, a, ip, b, 1, info);
        chkxer("Cspsv", infot, nout, lerr, ok);
        infot = 2;
        Cspsv("U", -1, 0, a, ip, b, 1, info);
        chkxer("Cspsv", infot, nout, lerr, ok);
        infot = 3;
        Cspsv("U", 0, -1, a, ip, b, 1, info);
        chkxer("Cspsv", infot, nout, lerr, ok);
        infot = 7;
        Cspsv("U", 2, 0, a, ip, b, 1, info);
        chkxer("Cspsv", infot, nout, lerr, ok);
        //
        //        Cspsvx
        //
        strncpy(srnamt, "Cspsvx", srnamt_len);
        infot = 1;
        Cspsvx("/", "U", 0, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cspsvx", infot, nout, lerr, ok);
        infot = 2;
        Cspsvx("N", "/", 0, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cspsvx", infot, nout, lerr, ok);
        infot = 3;
        Cspsvx("N", "U", -1, 0, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cspsvx", infot, nout, lerr, ok);
        infot = 4;
        Cspsvx("N", "U", 0, -1, a, af, ip, b, 1, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cspsvx", infot, nout, lerr, ok);
        infot = 9;
        Cspsvx("N", "U", 2, 0, a, af, ip, b, 1, x, 2, rcond, r1, r2, w, rw, info);
        chkxer("Cspsvx", infot, nout, lerr, ok);
        infot = 11;
        Cspsvx("N", "U", 2, 0, a, af, ip, b, 2, x, 1, rcond, r1, r2, w, rw, info);
        chkxer("Cspsvx", infot, nout, lerr, ok);
    }
    //
    //     Print a summary line.
    //
    if (ok) {
        write(nout, "(1x,a3,' drivers passed the tests of the error exits')"), path;
    } else {
        write(nout, "(' *** ',a3,' drivers failed the tests of the error ','exits ***')"), path;
    }
    //
    //     End of Cerrvx
    //
}
