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
#include <mplapack_debug.h>
#include <complex>

#if defined VERBOSE_TEST
#include <iostream>
#endif

#define MAX_ITER 10

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

int main(int argc, char *argv[]) {
    printf("*** Testing complex start ***\n");

#if defined ___MPLAPACK_BUILD_WITH_GMP___
    mpf_set_default_prec(___MPLAPACK_DEFAULT_PRECISION___);
#endif

    // we need to specify explicitly.
    mpreal::set_default_prec(___MPLAPACK_DEFAULT_PRECISION___);
    mpcomplex::set_default_prec(___MPLAPACK_DEFAULT_PRECISION___);

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
    printf("*** Testing complex successful ***\n");
    return (0);
}
