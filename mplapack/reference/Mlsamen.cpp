/*
 * Copyright (c) 2008-2025
 *	Nakata, Maho
 * 	All rights reserved.
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

/*
Based on http://www.netlib.org/blas/lsame.f
Mlsame returns 1 if CA is the same letter as CB regardless of case.
*/

#include <mpblas.h>
#include <mplapack.h>

bool Mlsamen(INTEGER n, const char *a, const char *b) {
    // Compare first n characters case-insensitively (ASCII only).
    // This matches BLAS/LAPACK intent better than locale-dependent toupper().
    if (n <= 0)
        return true;
    if (a == nullptr || b == nullptr)
        return false;

    for (INTEGER i = 0; i < n; ++i) {
        unsigned char ca = static_cast<unsigned char>(a[i]);
        unsigned char cb = static_cast<unsigned char>(b[i]);

        // Inline ASCII upper: convert only 'a'..'z' to 'A'..'Z'
        if (ca >= static_cast<unsigned char>('a') && ca <= static_cast<unsigned char>('z')) {
            ca = static_cast<unsigned char>(ca - (static_cast<unsigned char>('a') - static_cast<unsigned char>('A')));
        }
        if (cb >= static_cast<unsigned char>('a') && cb <= static_cast<unsigned char>('z')) {
            cb = static_cast<unsigned char>(cb - (static_cast<unsigned char>('a') - static_cast<unsigned char>('A')));
        }

        if (ca != cb)
            return false;
    }
    return true;
}
