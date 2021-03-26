/*
 * Copyright (c) 2010, 2012
 *	Nakata, Maho
 * 	All rights reserved.
 *
 * $Id: Raxpy.cpp,v 1.11 2010/08/07 05:50:09 nakatamaho Exp $
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

#include <mpblas__Float64x.h>

void Raxpy_omp(mplapackint n, _Float64x da, _Float64x * dx, mplapackint incx, _Float64x * dy, mplapackint incy);
void Raxpy_ref(mplapackint n, _Float64x da, _Float64x * dx, mplapackint incx, _Float64x * dy, mplapackint incy);

#define SINGLEOROMP 1000

void Raxpy(mplapackint n, _Float64x da, _Float64x * dx, mplapackint incx, _Float64x * dy, mplapackint incy)
{
    _Float64x Zero = 0.0;
    mplapackint i;

    if (n <= 0)
	return;
    if (da == Zero)
	return;

    if (1) {
	Raxpy_omp(n, da, dx, incx, dy, incy);
    } else {
	Raxpy_ref(n, da, dx, incx, dy, incy);
    }

    return;
}
