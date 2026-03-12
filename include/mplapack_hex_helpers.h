/*
 * Copyright (c) 2026
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

#pragma once

inline void format_hex_double_fixedexp(char *buf, size_t n, double x) {
    if (n == 0) {
        return;
    }

    if (std::isnan(x)) {
        std::snprintf(buf, n, "@NaN@");
        return;
    }
    if (std::isinf(x)) {
        std::snprintf(buf, n, "%s@Inf@", std::signbit(x) ? "-" : "+");
        return;
    }

    // %+.13a gives normalized hexadecimal floating-point for normal doubles:
    //   +0x1.fffffffffffffp+1023
    // with exactly 13 hex digits after the point.
    char tmp[128];
    std::snprintf(tmp, sizeof(tmp), "%+.13a", x);

    // Find exponent marker.
    char *p = std::strchr(tmp, 'p');
    if (p == nullptr) {
        // Fallback: should not happen for valid %a output.
        std::snprintf(buf, n, "%s", tmp);
        return;
    }

    // Parse exponent and rewrite with fixed width like p+0000 / p-0105 / p+1023.
    long exp_val = std::strtol(p + 1, nullptr, 10);
    *p = '\0';

    std::snprintf(buf, n, "%sp%+05ld", tmp, exp_val);
}
