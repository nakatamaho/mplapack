/*
 * Copyright (c) 2008-2010
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Mxerbla.cpp,v 1.7 2010/08/07 05:50:10 nakatamaho Exp $
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
Based on http://www.netlib.org/blas/xerbla.f
Mxerbla is an error handler for the Mplapack routines.
*/

#if !defined  __MPLAPACK_ERRNO__
#define __MPLAPACK_ERRNO__
#endif

#include <stdio.h>

int mplapack_errno;

void Mxerbla___float128(const char *srname, int info)
{
    mplapack_errno = info;
    return;
}

void Mxerbla_mpfr(const char *srname, int info)
{
    mplapack_errno = info;
    return;
}

void Mxerbla_double(const char *srname, int info)
{
    mplapack_errno = info;
    return;
}

void Mxerbla_dd(const char *srname, int info)
{
    mplapack_errno = info;
    return;
}

void Mxerbla_qd(const char *srname, int info)
{
    mplapack_errno = info;
    return;
}

void Mxerbla_gmp(const char *srname, int info)
{
    mplapack_errno = info;
    return;
}
