/*
 * Copyright (c) 2008-2021
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
#include <mplapack_compare_debug.h>

#include <blas.h>
#include <lapack.h>

#if defined VERBOSE_TEST
#include <iostream>
#endif

void iMlaenv_test(void) {
    char name[16], name_ref[16], opts[16];
    int ispec_ref, n1_ref, n2_ref, n3_ref, n4_ref, ret_ref;
    INTEGER ispec, n1, n2, n3, n4, ret;

    ispec_ref = 1;
    n1_ref = 16;
    n2_ref = -1;
    n3_ref = -1;
    n4_ref = -1;
    char teststr[] = "RTRTRI";
    char testopts[] = "UN";
    strncpy(name, teststr, strlen(teststr));
    strncpy(name_ref, teststr, strlen(teststr));
    /* we set ilaenv (d,z) and iMlaenv (R, C) are equivalent */
    if (name_ref[0] == 'R') {
        name_ref[0] = 'd';
    }
    if (name_ref[0] == 'C') {
        name_ref[0] = 'z';
    }
    strncpy(opts, testopts, strlen(testopts));
    ispec = (INTEGER)ispec_ref;
    n1 = (INTEGER)n1_ref;
    n2 = (INTEGER)n2_ref;
    n3 = (INTEGER)n3_ref;
    n4 = (INTEGER)n4_ref;

    ret_ref = ilaenv_f77(&ispec_ref, name_ref, opts, &n1_ref, &n2_ref, &n3_ref, &n4_ref);
    ret = iMlaenv(ispec, name, opts, n1, n2, n3, n4);
    printf("%s %s, LAPACK:%d  MPLAPACK:%d\n", name, opts, (int)ret_ref, (int)ret);
    if ((int)ret_ref != (int)ret) {
        printf("*** Testing iMlaenv failed ***\n");
        printf("%s %s, LAPACK:%d  MPLAPACK:%d\n", name, opts, (int)ret_ref, (int)ret);
    }
}

int main(int argc, char *argv[]) {
    printf("*** Testing iMlaenv start ***\n");
    iMlaenv_test();
    printf("*** Testing iMlaenv successful ***\n");
    return (0);
}
