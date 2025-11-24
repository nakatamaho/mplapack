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

void Crotg(COMPLEX const &ca, COMPLEX const &cb, REAL const &c, COMPLEX const &s) {
    //
    REAL scale = 0.0;
    REAL norm = 0.0;
    COMPLEX alpha = 0.0;
    if (fem::cdabs(ca) == 0.0) {
        c = 0.0;
        s = COMPLEX(1.0, 0.0);
        ca = cb;
    } else {
        scale = fem::cdabs(ca) + fem::cdabs(cb);
        norm = scale * sqrt(pow2((fem::cdabs(ca / COMPLEX(scale, 0.0)))) + pow2((fem::cdabs(cb / COMPLEX(scale, 0.0)))));
        alpha = ca / fem::cdabs(ca);
        c = fem::cdabs(ca) / norm;
        s = alpha * conj(cb) / norm;
        ca = alpha * norm;
    }
}
