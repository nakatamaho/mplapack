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
#include <mplapack_compare_debug.h>
#include <complex>

#if defined ___MPLAPACK_BUILD_WITH_QD___ || defined ___MPLAPACK_BUILD_WITH_DD___
#include <qd/fpu.h>
static unsigned int g_oldcw;
static void __attribute__((constructor)) mplapack_test_fpu_init(void) { fpu_fix_start(&g_oldcw); }
static void __attribute__((destructor))  mplapack_test_fpu_fini(void) { fpu_fix_end(&g_oldcw); }
#endif

#define VERBOSE_TEST
#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MAX_ITER 10

void mpc_ptr_conversion_test() {
    printf("mpcomplex <=> mpc_ptr conversion test \n");

    int flag = 0;
    mpcomplex z(mpreal(1.25), mpreal(-2.5));
    mpc_t raw;
    mpc_init2(raw, mpcomplex::default_real_prec);

    mpc_set(raw, z, mpcomplex::default_rnd);
    if (abs(mpcomplex(raw) - z) > mpreal(EPSILON)) {
        flag = 1;
    }

    mpc_add(raw, raw, z, mpcomplex::default_rnd);
    mpc_set(z, raw, mpcomplex::default_rnd);
    if (abs(z - mpcomplex(mpreal(2.5), mpreal(-5.0))) > mpreal(EPSILON)) {
        flag = 1;
    }

    mpc_clear(raw);
    if (flag) {
        printf("mpcomplex <=> mpc_ptr conversion test failed\n");
        exit(1);
    } else {
        printf("mpcomplex <=> mpc_ptr conversion test passed\n");
    }
}

static int mpcomplex_precision_mismatch(const mpcomplex &z, mp_prec_t expected_re, mp_prec_t expected_im) {
    return z.get_prec_re() != expected_re || z.get_prec_im() != expected_im;
}

void mpc_constructor_precision_test() {
    printf("mpcomplex constructor precision test \n");

    int flag = 0;
    const mpc_rnd_t mode = mpcomplex::default_rnd;

    const mpcomplex from_pair(1.0, 2.0, 101, 137, mode);
    if (mpcomplex_precision_mismatch(from_pair, 101, 137)) {
        flag = 1;
    }

    const mpcomplex from_double(1.0, 103, 139, mode);
    if (mpcomplex_precision_mismatch(from_double, 103, 139)) {
        flag = 1;
    }

    const mpcomplex from_std(std::complex<double>(1.0, 2.0), 107, 149, mode);
    if (mpcomplex_precision_mismatch(from_std, 107, 149)) {
        flag = 1;
    }

    const mpcomplex from_ld(std::complex<long double>(1.0L, 2.0L), 109, 151, mode);
    if (mpcomplex_precision_mismatch(from_ld, 109, 151)) {
        flag = 1;
    }

    const mpcomplex from_strs("1.0", "2.0", 113, 157, mode);
    if (mpcomplex_precision_mismatch(from_strs, 113, 157)) {
        flag = 1;
    }

    if (flag) {
        printf("mpcomplex constructor precision test failed\n");
        exit(1);
    } else {
        printf("mpcomplex constructor precision test passed\n");
    }
}

void mpc_subst_test1() {
    printf("Complex <= Real test \n");
    REAL_REF diff;
    COMPLEX Ctemp1, Ctemp2, cdummy;
    REAL Ftemp;
    REAL rdummy = 0.0;

    int flag = 0;
    for (int i = 0; i < MAX_ITER; i++) {
        Ctemp1 = mpc_randomnumber(cdummy);
#if defined VERBOSE_TEST
        cout << "c1       ";
        printnum(Ctemp1);
        cout << endl;
#endif
        Ftemp = mpf_randomnumber(rdummy);
        Ctemp1 = Ftemp; // this is what we want to test
#if defined VERBOSE_TEST
        cout << "r1=      ";
        printnum(Ftemp);
        cout << endl;
        cout << "c2=      ";
        printnum(Ctemp1);
        cout << endl;
        cout << "residue=r1-c2 " << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
#endif
        Ctemp2 = Ctemp1 - Ftemp;
        diff = abs(Ctemp2);
#if defined VERBOSE_TEST
        cout << "DIFF ";
        printnum(diff);
        cout << endl;
#endif
        if (abs(diff) > EPSILON) {
            flag = 1;
        }
    }
    if (flag) {
        printf("Complex <= Real subst test failed\n");
        exit(1);
    } else {
        printf("Complex <= Real subst test passed\n");
    }
}

void mpc_abs_test() {
    printf("Complex abs test \n");
    COMPLEX_REF ctemp1;
    REAL_REF ftemp1;
    REAL_REF diff;
    COMPLEX Ctemp1, Ctemp2, cdummy;
    REAL Ftemp1;

    int flag = 0;

    for (int i = 0; i < MAX_ITER; i++) {
        set_random_number(ctemp1, Ctemp1);
        Ftemp1 = abs(Ctemp1);
        ftemp1 = abs(ctemp1);
        diff = abs(Ftemp1 - ftemp1);
#if defined VERBOSE_TEST
        cout << "C1 = ";
        printnum(Ctemp1);
        cout << endl;
        cout << "F1 = ";
        printnum(Ftemp1);
        cout << endl;
        cout << endl;
        cout << "c1 = ";
        printnum(ctemp1);
        cout << endl;
        cout << "f1 = ";
        printnum(ftemp1);
        cout << endl;
        cout << "residue=sqrt(c1*conj(c1))-f1 " << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "residue=F1-f1 " << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "DIFF = ";
        printnum(diff);
        cout << endl;
        cout << endl;
#endif
        if (diff > EPSILON) {
            flag = 1;
        }
    }
    if (flag) {
        printf("Complex ABS test failed\n");
        exit(1);
    } else {
        printf("Complex ABS test passed\n");
    }
}

void mpc_add_test1() {
    printf("Complex + Complex test \n");
    COMPLEX_REF ctemp1;
    COMPLEX_REF ctemp2;
    COMPLEX_REF ctemp3;

    COMPLEX Ctemp1;
    COMPLEX Ctemp2;
    COMPLEX Ctemp3;

    REAL_REF diff;

    int flag = 0;

    for (int i = 0; i < MAX_ITER; i++) {
        set_random_number(ctemp1, Ctemp1);
        set_random_number(ctemp2, Ctemp2);

        Ctemp3 = Ctemp1 + Ctemp2;
        ctemp3 = ctemp1 + ctemp2;
        diff = abs(Ctemp3 - ctemp3);

#if defined VERBOSE_TEST
        cout << "C1 = ";
        printnum(Ctemp1);
        cout << endl;
        cout << "C2 = ";
        printnum(Ctemp2);
        cout << endl;
        cout << "C3 = ";
        printnum(Ctemp3);
        cout << endl;
        cout << endl;
        cout << "c1 = ";
        printnum(ctemp1);
        cout << endl;
        cout << "c2 = ";
        printnum(ctemp2);
        cout << endl;
        cout << "c3 = ";
        printnum(ctemp3);
        cout << endl;

        cout << "residue=C1+C2-C3" << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "residue=C3-c3 " << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "DIFF = ";
        printnum(diff);
        cout << endl;
        cout << endl;
#endif
        if (abs(diff) > EPSILON) {
            flag = 1;
        }
    }
    if (flag) {
        printf("Complex + Complex ADD test failed\n");
        exit(1);
    } else {
        printf("Complex + Complex ADD test passed\n");
    }
}

void mpc_add_test2() {
    printf("Complex + Real test \n");
    COMPLEX_REF ctemp1;
    REAL_REF ftemp2;
    COMPLEX_REF ctemp3;
    REAL_REF diff;

    COMPLEX Ctemp1;
    REAL Ftemp2;
    COMPLEX Ctemp3;

    int flag = 0;

    for (int i = 0; i < MAX_ITER; i++) {
        set_random_number(ctemp1, Ctemp1);
        set_random_number(ftemp2, Ftemp2);

        Ctemp3 = Ctemp1 + Ftemp2;
        ctemp3 = ctemp1 + ftemp2;
        diff = abs(Ctemp3 - ctemp3);
#if defined VERBOSE_TEST
        cout << "C1 = ";
        printnum(Ctemp1);
        cout << endl;
        cout << "F2 = ";
        printnum(Ftemp2);
        cout << endl;
        cout << "C3 = ";
        printnum(Ctemp3);
        cout << endl;
        cout << endl;
        cout << "c1 = ";
        printnum(ctemp1);
        cout << endl;
        cout << "f2 = ";
        printnum(ftemp2);
        cout << endl;
        cout << "c3 = ";
        printnum(ctemp3);
        cout << endl;

        cout << "residue=C1+F2-C3" << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "residue=C3-c3 " << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "DIFF = ";
        printnum(diff);
        cout << endl;
        cout << endl;
#endif
        if (diff > EPSILON) {
            flag = 1;
        }
    }
    if (flag) {
        printf("Complex + Real ADD test failed\n");
        exit(1);
    } else {
        printf("Complex + Real ADD test passed\n");
    }
}

void mpc_add_test3() {
    printf("Real + Complex test \n");
    COMPLEX_REF ctemp2;
    REAL_REF ftemp1;
    COMPLEX_REF ctemp3;
    REAL_REF diff;

    COMPLEX Ctemp2;
    REAL Ftemp1;
    COMPLEX Ctemp3;

    int flag = 0;

    for (int i = 0; i < MAX_ITER; i++) {
        set_random_number(ftemp1, Ftemp1);
        set_random_number(ctemp2, Ctemp2);

        Ctemp3 = Ftemp1 + Ctemp2;
        ctemp3 = ftemp1 + ctemp2;
#if defined VERBOSE_TEST
        cout << "F1 = ";
        printnum(Ftemp1);
        cout << endl;
        cout << "C2 = ";
        printnum(Ctemp2);
        cout << endl;
        cout << "C3 = ";
        printnum(Ctemp3);
        cout << endl;
        cout << endl;
        cout << "f1 = ";
        printnum(ftemp1);
        cout << endl;
        cout << "c2 = ";
        printnum(ctemp2);
        cout << endl;
        cout << "c3 = ";
        printnum(ctemp3);
        cout << endl;
#endif
        diff = abs(Ctemp3 - ctemp3);
#if defined VERBOSE_TEST
        cout << "residue=F1+C2-C3" << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "residue=C3-c3 " << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "DIFF = ";
        printnum(diff);
        cout << endl;
        cout << endl;
#endif
        if (abs(diff) > EPSILON) {
            flag = 1;
        }
    }
    if (flag) {
        printf("Real + Complex test failed\n");
        exit(1);
    } else {
        printf("Real + Complex test passed\n");
    }
}

void mpc_mul_test1() {
    printf("Complex * Complex test \n");
    COMPLEX_REF ctemp1;
    COMPLEX_REF ctemp2;
    COMPLEX_REF ctemp3;
    REAL_REF diff;

    COMPLEX Ctemp1;
    COMPLEX Ctemp2;
    COMPLEX Ctemp3;

    int flag = 0;

    for (int i = 0; i < MAX_ITER; i++) {
        set_random_number(ctemp1, Ctemp1);
        set_random_number(ctemp2, Ctemp2);

        Ctemp3 = Ctemp1 * Ctemp2;
        ctemp3 = ctemp1 * ctemp2;

#if defined VERBOSE_TEST
        cout << "C1 = ";
        printnum(Ctemp1);
        cout << endl;
        cout << "C2 = ";
        printnum(Ctemp2);
        cout << endl;
        cout << "C3 = ";
        printnum(Ctemp3);
        cout << endl;
        cout << endl;
        cout << "c1 = ";
        printnum(ctemp1);
        cout << endl;
        cout << "c2 = ";
        printnum(ctemp2);
        cout << endl;
        cout << "c3 = ";
        printnum(ctemp3);
        cout << endl;
#endif
        diff = abs(Ctemp3 - ctemp3);
#if defined VERBOSE_TEST
        cout << "residue=C1*C2-C3" << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "residue=C3-c3 " << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "DIFF = ";
        printnum(diff);
        cout << endl;
        cout << endl;
#endif
        if (abs(diff) > EPSILON) {
            flag = 1;
        }
    }
    if (flag) {
        printf("Complex * Complex test failed\n");
        exit(1);
    } else {
        printf("Complex * Complex test passed\n");
    }
}

void mpc_mul_test2() {
    printf("Complex * Real test \n");
    COMPLEX_REF ctemp1;
    REAL_REF ftemp2;
    COMPLEX_REF ctemp3;
    REAL_REF diff;

    COMPLEX Ctemp1;
    REAL Ftemp2;
    COMPLEX Ctemp3;

    int flag = 0;

    for (int i = 0; i < MAX_ITER; i++) {
        set_random_number(ctemp1, Ctemp1);
        set_random_number(ftemp2, Ftemp2);

        ctemp1 = Ctemp1;
        ftemp2 = Ftemp2;
        Ctemp3 = Ctemp1 * Ftemp2;
        ctemp3 = ctemp1 * ftemp2;

#if defined VERBOSE_TEST
        cout << "C1 = ";
        printnum(Ctemp1);
        cout << endl;
        cout << "F2 = ";
        printnum(Ftemp2);
        cout << endl;
        cout << "C3 = ";
        printnum(Ctemp3);
        cout << endl;
        cout << endl;
        cout << "c1 = ";
        printnum(ctemp1);
        cout << endl;
        cout << "f2 = ";
        printnum(ftemp2);
        cout << endl;
        cout << "c3 = ";
        printnum(ctemp3);
        cout << endl;
#endif
        diff = abs(Ctemp3 - ctemp3);
#if defined VERBOSE_TEST
        cout << "residue=C1*F2-C3" << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "residue=C3-c3 " << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "DIFF = ";
        printnum(diff);
        cout << endl;
        cout << endl;
#endif
        if (abs(diff) > EPSILON) {
            flag = 1;
        }
    }
    if (flag) {
        printf("Complex * Real test failed\n");
        exit(1);
    } else {
        printf("Complex * Real test passed\n");
    }
}

void mpc_mul_test3() {
    printf("Real * Complex test \n");
    COMPLEX_REF ctemp2;
    REAL_REF ftemp1;
    COMPLEX_REF ctemp3;
    REAL_REF diff;

    COMPLEX Ctemp2;
    REAL Ftemp1;
    COMPLEX Ctemp3;

    int flag = 0;

    for (int i = 0; i < MAX_ITER; i++) {
        set_random_number(ftemp1, Ftemp1);
        set_random_number(ctemp2, Ctemp2);

        ctemp2 = Ctemp2;
        ftemp1 = Ftemp1;
        Ctemp3 = Ftemp1 * Ctemp2;
        ctemp3 = ftemp1 * ctemp2;
#if defined VERBOSE_TEST
        cout << "F1 = ";
        printnum(Ftemp1);
        cout << endl;
        cout << "C2 = ";
        printnum(Ctemp2);
        cout << endl;
        cout << "C3 = ";
        printnum(Ctemp3);
        cout << endl;
        cout << endl;
        cout << "f1 = ";
        printnum(ftemp1);
        cout << endl;
        cout << "c2 = ";
        printnum(ctemp2);
        cout << endl;
        cout << "c3 = ";
        printnum(ctemp3);
        cout << endl;
#endif
        diff = abs(Ctemp3 - ctemp3);
#if defined VERBOSE_TEST
        cout << "residue=C1*F2-C3" << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "residue=C3-c3 " << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "DIFF = ";
        printnum(diff);
        cout << endl;
        cout << endl;
#endif
        if (abs(diff) > EPSILON) {
            flag = 1;
        }
    }
    if (flag) {
        printf("Complex * Real test failed\n");
        exit(1);
    } else {
        printf("Complex * Real test passed\n");
    }
}

void mpc_div_test1() {
    printf("Complex / Complex test \n");
    REAL_REF diff;
    COMPLEX_REF ctemp1;
    COMPLEX_REF ctemp2;
    COMPLEX_REF ctemp3;

    COMPLEX Ctemp1;
    COMPLEX Ctemp2;
    COMPLEX Ctemp3;

    int flag = 0;

    for (int i = 0; i < MAX_ITER; i++) {
        set_random_number(ctemp1, Ctemp1);
        set_random_number(ctemp2, Ctemp2);

        Ctemp3 = Ctemp1 / Ctemp2;
        ctemp3 = ctemp1 / ctemp2;

#if defined VERBOSE_TEST
        cout << "C1 = ";
        printnum(Ctemp1);
        cout << endl;
        cout << "C2 = ";
        printnum(Ctemp2);
        cout << endl;
        cout << "C3 = ";
        printnum(Ctemp3);
        cout << endl;
        cout << endl;
        cout << "c1 = ";
        printnum(ctemp1);
        cout << endl;
        cout << "c2 = ";
        printnum(ctemp2);
        cout << endl;
        cout << "c3 = ";
        printnum(ctemp3);
        cout << endl;
#endif
        diff = abs(Ctemp3 - ctemp3);
#if defined VERBOSE_TEST
        cout << "residue=C1/F2-C3" << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "residue=C3-c3 " << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "DIFF = ";
        printnum(diff);
        cout << endl;
        cout << endl;
#endif
        if (abs(diff) > EPSILON) {
            flag = 1;
        }
    }
    if (flag) {
        printf("Complex / Complex test failed\n");
        exit(1);
    } else {
        printf("Complex / Complex test passed\n");
    }
}

void mpc_div_test2() {
    printf("Complex / Real test \n");

    COMPLEX_REF ctemp1;
    REAL_REF ftemp2;
    COMPLEX_REF ctemp3;
    COMPLEX Ctemp1;
    REAL Ftemp2;
    COMPLEX Ctemp3;

    REAL_REF diff;

    int flag = 0;

    for (int i = 0; i < MAX_ITER; i++) {
        set_random_number(ctemp1, Ctemp1);
        set_random_number(ftemp2, Ftemp2);

        Ctemp3 = Ctemp1 / Ftemp2;
        ctemp3 = ctemp1 / ftemp2;
#if defined VERBOSE_TEST
        cout << "C1 = ";
        printnum(Ctemp1);
        cout << endl;
        cout << "f2 = ";
        printnum(Ftemp2);
        cout << endl;
        cout << "C3 = ";
        printnum(Ctemp3);
        cout << endl;
        cout << endl;
        cout << "c1 = ";
        printnum(ctemp1);
        cout << endl;
        cout << "f2 = ";
        printnum(ftemp2);
        cout << endl;
        cout << "c3 = ";
        printnum(ctemp3);
        cout << endl;
#endif
        diff = abs(Ctemp3 - ctemp3);
#if defined VERBOSE_TEST
        cout << "residue=C1/F2-C3" << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "residue=C3-c3 " << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "DIFF = ";
        printnum(diff);
        cout << endl;
        cout << endl;
#endif
        if (abs(diff) > EPSILON) {
            flag = 1;
        }
    }
    if (flag) {
        printf("Complex / Real test failed\n");
    } else {
        printf("Complex / Real test passed\n");
    }
}

void mpc_div_test3() {
    printf("Real / Complex test \n");
    COMPLEX_REF ctemp2;
    REAL_REF ftemp1;
    COMPLEX_REF ctemp3;
    REAL_REF diff;

    COMPLEX Ctemp2;
    REAL Ftemp1;
    COMPLEX Ctemp3;

    int flag = 0;

    for (int i = 0; i < MAX_ITER; i++) {
        set_random_number(ftemp1, Ftemp1);
        set_random_number(ctemp2, Ctemp2);

        Ctemp3 = Ftemp1 / Ctemp2;
        ctemp3 = ftemp1 / ctemp2;
#if defined VERBOSE_TEST
        cout << "F1 = ";
        printnum(Ftemp1);
        cout << endl;
        cout << "C2 = ";
        printnum(Ctemp2);
        cout << endl;
        cout << "C3 = ";
        printnum(Ctemp3);
        cout << endl;
        cout << endl;
        cout << "f1 = ";
        printnum(ftemp1);
        cout << endl;
        cout << "c2 = ";
        printnum(ctemp2);
        cout << endl;
        cout << "c3 = ";
        printnum(ctemp3);
        cout << endl;
#endif
        diff = abs(Ctemp3 - ctemp3);
#if defined VERBOSE_TEST
        cout << "residue=F1/C2-C3" << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "residue=C3-c3 " << endl;
        cout << "if(abs(residue)< " << EPSILON << ") printf \"ok\\n\"; else printf \"ng\\n\"; endif" << endl;
        cout << "DIFF = ";
        printnum(diff);
        cout << endl;
        cout << endl;
#endif
        if (abs(diff) > EPSILON) {
            flag = 1;
        }
    }
    if (flag) {
        printf("Real / Complex test failed\n");
        exit(1);
    } else {
        printf("Real / Complex test passed\n");
    }
}

#if defined ___MPLAPACK_BUILD_WITH_GMP___
void mpc_algebraic_test() {
    printf("GMP complex algebraic helper test \n");

    int flag = 0;
    const REAL tolerance = EPSILON;

    mpc_class z(mpf_class(6.0), mpf_class(-8.0));
    mpc_class z_div_assign(z);
    z_div_assign /= mpf_class(2.0);
    if (abs(z_div_assign - mpc_class(mpf_class(3.0), mpf_class(-4.0))) > tolerance) {
        flag = 1;
    }

    mpc_class real_one(mpf_class(1.0), mpf_class(0.0));
    mpc_class complex_one(mpf_class(1.0), mpf_class(1.0));
    if (real_one != 1.0) {
        flag = 1;
    }
    if (1.0 != complex_one) {
        // expected
    } else {
        flag = 1;
    }

    mpc_class c(mpf_class(3.0), mpf_class(4.0));
    if (abs(real(c) - mpf_class(3.0)) > tolerance) {
        flag = 1;
    }
    if (abs(imag(c) - mpf_class(4.0)) > tolerance) {
        flag = 1;
    }
    if (abs(norm(c) - mpf_class(25.0)) > tolerance) {
        flag = 1;
    }
    if (abs(abs(c) - mpf_class(5.0)) > tolerance) {
        flag = 1;
    }
    if (abs(conj(c) - mpc_class(mpf_class(3.0), mpf_class(-4.0))) > tolerance) {
        flag = 1;
    }

    mpc_class swap_a(mpf_class(1.0), mpf_class(2.0));
    mpc_class swap_b(mpf_class(3.0), mpf_class(4.0));
    swap(swap_a, swap_b);
    if (abs(swap_a - mpc_class(mpf_class(3.0), mpf_class(4.0))) > tolerance) {
        flag = 1;
    }
    if (abs(swap_b - mpc_class(mpf_class(1.0), mpf_class(2.0))) > tolerance) {
        flag = 1;
    }

    mpc_class base(mpf_class(2.0), mpf_class(3.0));
    std::complex<double> scalar(1.5, -0.5);
    if (abs((base + scalar) - mpc_class(mpf_class(3.5), mpf_class(2.5))) > tolerance) {
        flag = 1;
    }
    if (abs((scalar + base) - mpc_class(mpf_class(3.5), mpf_class(2.5))) > tolerance) {
        flag = 1;
    }
    if (abs((base - scalar) - mpc_class(mpf_class(0.5), mpf_class(3.5))) > tolerance) {
        flag = 1;
    }
    if (abs((scalar - base) - mpc_class(mpf_class(-0.5), mpf_class(-3.5))) > tolerance) {
        flag = 1;
    }
    if (abs((base * scalar) - mpc_class(mpf_class(4.5), mpf_class(3.5))) > tolerance) {
        flag = 1;
    }

    const mpc_class i(mpf_class(0.0), mpf_class(1.0));
    const mpc_class w(mpf_class(0.25), mpf_class(0.125));
    const REAL transcendental_tolerance = EPSILON * 1024.0;

    if (abs(exp(log(w)) - w) > transcendental_tolerance) {
        flag = 1;
    }
    if (abs(sin(i * mpf_class(0.25)) - i * sinh(mpf_class(0.25))) > transcendental_tolerance) {
        flag = 1;
    }
    if (abs(cos(i * mpf_class(0.25)) - mpc_class(cosh(mpf_class(0.25)), mpf_class(0.0))) > transcendental_tolerance) {
        flag = 1;
    }
    if (abs(tan(w) - sin(w) / cos(w)) > transcendental_tolerance) {
        flag = 1;
    }
    if (abs(sinh(i * mpf_class(0.25)) - i * sin(mpf_class(0.25))) > transcendental_tolerance) {
        flag = 1;
    }
    if (abs(cosh(i * mpf_class(0.25)) - mpc_class(cos(mpf_class(0.25)), mpf_class(0.0))) > transcendental_tolerance) {
        flag = 1;
    }
    if (abs(tanh(w) - sinh(w) / cosh(w)) > transcendental_tolerance) {
        flag = 1;
    }
    if (abs(polar(abs(w), arg(w)) - w) > transcendental_tolerance) {
        flag = 1;
    }
    if (abs(pow(w, mpf_class(2.0)) - w * w) > transcendental_tolerance) {
        flag = 1;
    }
    if (abs(pow(w, (mplapackint)2) - w * w) > transcendental_tolerance) {
        flag = 1;
    }
    if (abs(asin(sin(w)) - w) > transcendental_tolerance) {
        flag = 1;
    }
    if (abs(acos(cos(w)) - w) > transcendental_tolerance) {
        flag = 1;
    }
    if (abs(atan(tan(w)) - w) > transcendental_tolerance) {
        flag = 1;
    }
    if (abs(asinh(sinh(w)) - w) > transcendental_tolerance) {
        flag = 1;
    }
    if (abs(atanh(tanh(w)) - w) > transcendental_tolerance) {
        flag = 1;
    }
    if (abs(acosh(cosh(w)) - w) > transcendental_tolerance) {
        flag = 1;
    }

    mpcomplex mp_a(mpc_class(mpf_class(2.0), mpf_class(3.0)));
    mpc_class gmp_b(mpf_class(1.0), mpf_class(0.5));
    if (abs((mp_a - gmp_b) - mpcomplex(mpc_class(mpf_class(1.0), mpf_class(2.5)))) > mpreal(EPSILON)) {
        flag = 1;
    }
    if (abs((gmp_b - mp_a) - mpcomplex(mpc_class(mpf_class(-1.0), mpf_class(-2.5)))) > mpreal(EPSILON)) {
        flag = 1;
    }

    if (flag) {
        printf("GMP complex algebraic helper test failed\n");
        exit(1);
    } else {
        printf("GMP complex algebraic helper test passed\n");
    }
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_DD___ || defined ___MPLAPACK_BUILD_WITH_QD___
void qd_dd_complex_helper_test() {
    printf("QD/DD complex helper test \n");

    int flag = 0;
    const REAL tolerance = EPSILON;

    COMPLEX z(REAL(6.0), REAL(-8.0));
    z /= REAL(2.0);
    if (abs(z - COMPLEX(REAL(3.0), REAL(-4.0))) > tolerance) {
        flag = 1;
    }

    const COMPLEX base(REAL(2.0), REAL(3.0));
    const std::complex<double> scalar(1.5, -0.5);
    if (abs((base + scalar) - COMPLEX(REAL(3.5), REAL(2.5))) > tolerance) {
        flag = 1;
    }
    if (abs((scalar + base) - COMPLEX(REAL(3.5), REAL(2.5))) > tolerance) {
        flag = 1;
    }
    if (abs((base - scalar) - COMPLEX(REAL(0.5), REAL(3.5))) > tolerance) {
        flag = 1;
    }
    if (abs((scalar - base) - COMPLEX(REAL(-0.5), REAL(-3.5))) > tolerance) {
        flag = 1;
    }
    if (abs((base * scalar) - COMPLEX(REAL(4.5), REAL(3.5))) > tolerance) {
        flag = 1;
    }
    if (abs((scalar * base) - COMPLEX(REAL(4.5), REAL(3.5))) > tolerance) {
        flag = 1;
    }
    if (abs((base / scalar) * scalar - base) > tolerance) {
        flag = 1;
    }
    if (abs((scalar / base) * base - COMPLEX(scalar)) > tolerance) {
        flag = 1;
    }

    const mpcomplex mp_a(base);
    const COMPLEX b(REAL(1.0), REAL(0.5));
    if (abs((mp_a - b) - mpcomplex(COMPLEX(REAL(1.0), REAL(2.5)))) > mpreal(EPSILON)) {
        flag = 1;
    }
    if (abs((b - mp_a) - mpcomplex(COMPLEX(REAL(-1.0), REAL(-2.5)))) > mpreal(EPSILON)) {
        flag = 1;
    }

    if (flag) {
        printf("QD/DD complex helper test failed\n");
        exit(1);
    } else {
        printf("QD/DD complex helper test passed\n");
    }
}
#endif

int main(int argc, char *argv[]) {
    printf("*** Testing complex start ***\n");

#if defined ___MPLAPACK_BUILD_WITH_GMP___
    mpf_set_default_prec(___MPLAPACK_GMP_DEFAULT_PRECISION___);
#endif

    // we need to specify explicitly.
    mpreal::default_prec = ___MPLAPACK_MPFR_DEFAULT_PRECISION___;
    mpcomplex::default_real_prec = ___MPLAPACK_MPFR_DEFAULT_PRECISION___;
    mpcomplex::default_imag_prec = ___MPLAPACK_MPFR_DEFAULT_PRECISION___;

    mpc_ptr_conversion_test();
    mpc_constructor_precision_test();
    mpc_subst_test1();
    mpc_abs_test();
    mpc_add_test1();
    mpc_add_test2();
    mpc_add_test3();
    mpc_mul_test1();
    mpc_mul_test2();
    mpc_mul_test3();
    mpc_div_test1();
    mpc_div_test2();
    mpc_div_test3();
#if defined ___MPLAPACK_BUILD_WITH_GMP___
    mpc_algebraic_test();
#endif
#if defined ___MPLAPACK_BUILD_WITH_DD___ || defined ___MPLAPACK_BUILD_WITH_QD___
    qd_dd_complex_helper_test();
#endif
    printf("*** Testing complex successful ***\n");
    return (0);
}
