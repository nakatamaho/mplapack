/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasrt.cpp,v 1.9 2010/08/07 04:48:33 nakatamaho Exp $ 
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

#include <mblas.h>
#include <mlapack.h>
#include <stdlib.h>

int compare_mpf_gt(const REAL * a, const REAL * b)
{
    if (*a > *b)
	return 1;
    if (*a == *b)
	return 0;
    if (*a < *b)
	return -1;
    return 0;			//never occurs
}

int compare_mpf_lt(const REAL * a, const REAL * b)
{
    if (*a > *b)
	return -1;
    if (*a == *b)
	return 0;
    if (*a < *b)
	return 1;
    return 0;			//never occurs
}

void Rlasrt(const char *id, INTEGER n, REAL * d, INTEGER * info)
{
    //Error check
    if (!Mlsame(id, "I") && !Mlsame(id, "D")) {
	*info = -1;
	Mxerbla("Rlasrt", -(*info));
	return;
    }
    if (n < 0) {
	*info = -2;
	Mxerbla("Rlasrt", -(*info));
	return;
    }
    if (Mlsame(id, "I")) {
	qsort(d, n, sizeof(REAL), (int (*)(const void *, const void *)) compare_mpf_gt);
    }
    if (Mlsame(id, "d")) {
	qsort(d, n, sizeof(REAL), (int (*)(const void *, const void *)) compare_mpf_lt);
    }
    *info = 0;
}
