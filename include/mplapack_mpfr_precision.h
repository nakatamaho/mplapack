/*
 * Copyright (c) 2008-2026
 *	Nakata, Maho
 *	All rights reserved.
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
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef _MPLAPACK_MPFR_PRECISION_H_
#define _MPLAPACK_MPFR_PRECISION_H_

#include <mpfrxx_mkII.h>
#include <stdexcept>

// Propagate caller-selected precision into an MPFR arithmetic worker.
// The MPLAPACK MPFR contract requires all operands of one invocation to
// have the same precision.  This helper therefore accepts one precision
// value and performs no operand inspection or consistency checking.
class MplapackMpfrPrecisionScope {
public:
    explicit MplapackMpfrPrecisionScope(mpfr_prec_t precision)
        : previous_precision_(mpfrxx::default_precision_bits()) {
        if (precision < MPFR_PREC_MIN || precision > MPFR_PREC_MAX) {
            throw std::invalid_argument("invalid MPFR precision for MPLAPACK scope");
        }
        mpfrxx::set_default_precision_bits(precision);
    }

    ~MplapackMpfrPrecisionScope() noexcept {
        mpfrxx::set_default_precision_bits(previous_precision_);
    }

    MplapackMpfrPrecisionScope(const MplapackMpfrPrecisionScope &) = delete;
    MplapackMpfrPrecisionScope &operator=(const MplapackMpfrPrecisionScope &) = delete;

private:
    mpfr_prec_t previous_precision_;
};

#endif
