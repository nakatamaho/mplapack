/*
 * Copyright (c) 2010-2012
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

#include <mpreal.h>
#include <mpcomplex.h>
#include <complex>
#include <mpc.h>

namespace mpfr{
#if MPC_VERSION_MAJOR > 0
  mpc_rnd_t  mpcomplex::default_rnd = MPC_RND(mpfr_get_default_rounding_mode(), mpfr_get_default_rounding_mode());
#else
  mpc_rnd_t  mpcomplex::default_rnd = RNDC(mpfr_get_default_rounding_mode(), mpfr_get_default_rounding_mode());
#endif
  mp_prec_t  mpcomplex::default_real_prec = mpfr_get_default_prec();	
  mp_prec_t  mpcomplex::default_imag_prec = mpfr_get_default_prec();	
  int        mpcomplex::default_base = 10;
  int        mpcomplex::double_bits = -1;
}
