/*
 * Copyright (c) 2026
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

// Device code of the MPFR complex BLAS with the runtime-precision types
// (MPFR/MPC functions; Cblas_mpfr_dev_types.h) for 2, 4, 8, 16, 32 limbs (up to 2048 bits).

#include "Cblas_mpfr_device_cuda.h"

namespace mplapack_blas_dev {
MPLAPACK_MPFR_CBLAS_JOBS(MPLAPACK_MPFR_CBLAS_RT_JOB, 2)
MPLAPACK_MPFR_CBLAS_JOBS(MPLAPACK_MPFR_CBLAS_RT_JOB, 4)
MPLAPACK_MPFR_CBLAS_JOBS(MPLAPACK_MPFR_CBLAS_RT_JOB, 8)
MPLAPACK_MPFR_CBLAS_JOBS(MPLAPACK_MPFR_CBLAS_RT_JOB, 16)
MPLAPACK_MPFR_CBLAS_JOBS(MPLAPACK_MPFR_CBLAS_RT_JOB, 32)
} // namespace mplapack_blas_dev
