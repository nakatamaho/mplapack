/*
 * Copyright (c) 2008-2021
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: mplapack_config.h,v 1.15 2010/08/07 03:15:46 nakatamaho Exp $
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

/* work in progress */
/* put some definitons on mplapack */

/* should depend on C compiler and environment
   our intention is that use 64bit int when USE64BITINT is set.
   This should be the default on 64bit environment.
*/

#ifndef _MPLAPACK_CONFIG_H_
#define _MPLAPACK_CONFIG_H_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <complex>
#include <inttypes.h>
#include <stdlib.h>

#undef F77_FUNC
#undef F77_FUNC_
#undef F77_FUNC_EQUIV

#define USE64BITINT

#ifdef USE64BITINT
    #if defined _WIN32  //workaround for Windows. long int is 32bit and int64_t is long long. Supporting GMP version is not straightfoward.
        typedef long int mplapackint;
    #elif defined __APPLE__ //on apple int64_t doesn't work, but it long works and it is 8 bytes.
        typedef long mplapackint;
    #else
        typedef int64_t mplapackint;
    #endif
#endif

typedef mplapackint mplapacklogical;

#ifdef __cplusplus
typedef mplapacklogical (*LFP)(...);
#else
typedef mplapacklogical(*LFP);
#endif

#endif
