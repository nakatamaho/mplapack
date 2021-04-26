/*
 * Copyright (c) 2008-2010
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
#include <mpblas.h>
#include <mplapack.h>
#include <mplapack_debug.h>

#define VERBOSE_TEST

#if defined VERBOSE_TEST
#include <iostream>
#endif

void Mlsame_test()
{
//  char a="A";
//  char b="A";
    int errorflag = FALSE;

    if (Mlsamen(2, "AB", "AB")) {
#if defined VERBOSE_TEST
	printf("same letter ok\n");
#endif
    } else {
#if defined VERBOSE_TEST
	printf("same letter NG\n");
#endif
	errorflag = TRUE;
    }
    if (Mlsamen(2, "AB", "ab")) {
#if defined VERBOSE_TEST
	printf("Uppercase/lowercase ok\n");
#endif
    } else {
#if defined VERBOSE_TEST
	printf("Uppercase/lowercase NG\n");
#endif
	errorflag = TRUE;
    }

    if (Mlsamen(3, "abc", "def")) {
#if defined VERBOSE_TEST
	printf("same! but should not be the same ng\n");
	errorflag = TRUE;
#endif
    } else {
#if defined VERBOSE_TEST
	printf("different ok\n");
#endif
    }
    if (errorflag == TRUE) {
        printf("*** Testing Mlsamen failed ***\n");
	exit(1);
    }
}

int main(int argc, char *argv[])
{
    printf("*** Testing Mlsamen start ***\n");
    Mlsame_test();
    printf("*** Testing Mlsamen successful ***\n");
    return (0);
}
