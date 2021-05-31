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

void Clarnv(INTEGER const idist, INTEGER *iseed, INTEGER const n, COMPLEX *x) {
    //
    INTEGER iv = 0;
    const INTEGER lv = 128;
    INTEGER il = 0;
    REAL u[lv];
    INTEGER i = 0;
    const REAL two = 2.0e+0;
    const REAL one = 1.0;
    const REAL zero = 0.0;
    const REAL twopi = two * pi(zero);
    for (iv = 1; iv <= n; iv = iv + lv / 2) {
        il = min(lv / 2, n - iv + 1);
        //
        //        Call DLARUV to generate 2*IL real numbers from a uniform (0,1)
        //        distribution (2*IL <= LV)
        //
        Rlaruv(iseed, 2 * il, u);
        //
        if (idist == 1) {
            //
            //           Copy generated numbers
            //
            for (i = 1; i <= il; i = i + 1) {
                x[(iv + i - 1) - 1] = COMPLEX(u[(2 * i - 1) - 1], u[(2 * i) - 1]);
            }
        } else if (idist == 2) {
            //
            //           Convert generated numbers to uniform (-1,1) distribution
            //
            for (i = 1; i <= il; i = i + 1) {
                x[(iv + i - 1) - 1] = COMPLEX(two * u[(2 * i - 1) - 1] - one, two * u[(2 * i) - 1] - one);
            }
        } else if (idist == 3) {
            //
            //           Convert generated numbers to normal (0,1) distribution
            //
            for (i = 1; i <= il; i = i + 1) {
                x[(iv + i - 1) - 1] = sqrt(-two * log(u[(2 * i - 1) - 1])) * exp(COMPLEX(zero, twopi * u[(2 * i) - 1]));
            }
        } else if (idist == 4) {
            //
            //           Convert generated numbers to complex numbers uniformly
            //           distributed on the unit disk
            //
            for (i = 1; i <= il; i = i + 1) {
                x[(iv + i - 1) - 1] = sqrt(u[(2 * i - 1) - 1]) * exp(COMPLEX(zero, twopi * u[(2 * i) - 1]));
            }
        } else if (idist == 5) {
            //
            //           Convert generated numbers to complex numbers uniformly
            //           distributed on the unit circle
            //
            for (i = 1; i <= il; i = i + 1) {
                x[(iv + i - 1) - 1] = exp(COMPLEX(zero, twopi * u[(2 * i) - 1]));
            }
        }
    }
    //
    //     End of Clarnv
    //
}
