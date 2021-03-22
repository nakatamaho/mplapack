/*
 * Copyright (c) 2008-2012
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

#define __MPLAPACK_BUILD_DEBUG_CPP__

#include <mpblas.h>
#include <mplapack_debug.h>

void __attribute__ ((constructor)) mplapack_debug_initialize(void);
void mplapack_debug_initialize(void)
{
    gmp_randinit_default(uniformrandomstate_mpfr);
}

#if defined ___MPLAPACK_BUILD_WITH_GMP___
void __attribute__ ((constructor)) mplapack_debug_initialize_gmp(void);
void mplapack_debug_initialize_gmp(void)
{
    uniformrandomstate_gmp = new gmp_randclass(gmp_randinit_default);
}
#endif

mpreal mpf_randomnumber(mpreal dummy)
{
    mpreal mtmp;

    mtmp = urandomb(uniformrandomstate_mpfr);
    mtmp = 2.0 * mtmp - 1.0;

    return mtmp;
}

mpcomplex mpc_randomnumber(mpcomplex dummy)
{
    mpcomplex ctmp = urandom_c(uniformrandomstate_mpfr);
    return ctmp;
}

void printnum(mpreal rtmp)
{
    mpfr_printf(MPFR_P_FORMAT, mpfr_ptr(rtmp));
    return;
}

void printnum(mpcomplex ctmp)
{
    mpreal cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_printf(MPFR_P_FORMAT MPFR_P_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
    return;
}

void sprintnum(char *buf, mpreal rtmp)
{
    mpfr_snprintf(buf, BUFLEN, MPFR_P_FORMAT, mpfr_ptr(rtmp));
    return;
}

void sprintnum(char *buf, mpcomplex ctmp)
{
    mpreal cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_snprintf(buf, BUFLEN, MPFR_P_FORMAT MPFR_P_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
    return;
}

#if defined ___MPLAPACK_BUILD_WITH_MPFR___
double mpf_randomnumber(double dummy)
{
    double mtmp = drand48();
    return mtmp;
}

complex < double >mpc_randomnumber(complex < double >)
{
    std::complex < double >ctmp;
    ctmp.real(drand48());
    ctmp.imag(drand48());
    return ctmp;
}

void set_random_number(double &a, mpreal & b)
{
    double dummy = 0.0;
    a = mpf_randomnumber(dummy);
    b = a;
}

void set_random_number(complex<double> &a, mpcomplex & b)
{
    complex<double> dummy, p;
    a = mpc_randomnumber(dummy);
    b = a;
}

void set_random_number1to2(double &a, mpreal & b)
{
    double dummy = 0.0;
    a = mpf_randomnumber(dummy);
    if (a > 0.0) a = a + 1.0; else a = a - 1.0;
    b = a;
}

void set_random_number1to2(complex<double> &a, mpcomplex & b)
{
    complex<double> dummy;
    double p, q;
    a = mpc_randomnumber(dummy);
    if (a.real() > 0.0) p = 1.0; else p = -1.0;
    if (a.imag() > 0.0) q = 1.0; else q = -1.0;
    a = a + complex<double>(p, q);
    b = a;
}

#endif

#if defined ___MPLAPACK_BUILD_WITH_GMP___
void printnum(mpf_class rtmp)
{
    gmp_printf(GMP_P_FORMAT, rtmp.get_mpf_t());
    return;
}

void printnum(mpc_class ctmp)
{
    gmp_printf(GMP_P_FORMAT GMP_P_FORMAT "i", ctmp.real().get_mpf_t(), ctmp.imag().get_mpf_t());
    return;
}

void sprintnum(char *buf, mpf_class rtmp)
{
    gmp_snprintf(buf, BUFLEN, GMP_P_FORMAT, rtmp.get_mpf_t());
    return;
}

void sprintnum(char *buf, mpc_class ctmp)
{
    gmp_snprintf(buf, BUFLEN, GMP_P_FORMAT GMP_P_FORMAT "i", ctmp.real().get_mpf_t(), ctmp.imag().get_mpf_t());
    return;
}

mpf_class mpf_randomnumber(mpf_class dummy)
{
    mpf_class mtmp;

    mtmp = uniformrandomstate_gmp->get_f();
    mtmp = 2.0 * mtmp - 1.0;
    return mtmp;
}

mpc_class mpc_randomnumber(mpc_class dummy)
{
    mpc_class ctmp;
    mpf_class mtmp1;
    mpf_class mtmp2;

    mtmp1 = uniformrandomstate_gmp->get_f();
    mtmp1 = 2.0 * mtmp1 - 1.0;

    mtmp2 = uniformrandomstate_gmp->get_f();
    mtmp2 = 2.0 * mtmp2 - 1.0;

    ctmp.real(mtmp1);
    ctmp.imag(mtmp2);
    return ctmp;
}

void set_random_number(mpreal & a, mpf_class & b)
{
    mpreal dummy;
    a = mpf_randomnumber(dummy);
    b = cast2mpf_class(a);
}

void set_random_number(mpcomplex & a, mpc_class & b)
{
    mpcomplex dummy;
    a = mpc_randomnumber(dummy);
    b = cast2mpc_class(a);
}

void set_random_number1to2(mpreal & a, mpf_class & b)
{
    mpreal dummy;
    a = mpf_randomnumber(dummy);
    if (a > 0.0) a = a + 1.0; else a = a - 1.0;
    b = cast2mpf_class(a);
}

void set_random_number1to2(mpcomplex & a, mpc_class & b)
{
    mpcomplex dummy;
    double p, q;
    a = mpc_randomnumber(dummy);
    if (a.real() > 0.0) p = 1.0; else p = -1.0;
    if (a.imag() > 0.0) q = 1.0; else q = -1.0;
    a = a + complex<double>(p, q);
    b = cast2mpc_class(a);
}

#endif

#if defined ___MPLAPACK_BUILD_WITH_QD___
void printnum(qd_real rtmp)
{
    cout.precision(QD_PRECISION);
    if (rtmp >= 0.0) {
	cout << "+" << rtmp;
    } else {
	cout << rtmp;
    }
    return;
}

void printnum(qd_complex rtmp)
{
    cout.precision(QD_PRECISION);
    if (rtmp.real() >= 0.0) {
	cout << "+" << rtmp.real();
    } else {
	cout << rtmp.real();
    }
    if (rtmp.imag() >= 0.0) {
	cout << "+" << rtmp.imag() << "i";
    } else {
	cout << rtmp.imag() << "i";
    }
    return;
}

void sprintnum(char *buf, qd_real rtmp)
{
    rtmp.write(buf, BUFLEN, 64);
    return;
}

void sprintnum(char *buf, qd_complex rtmp)
{
printf("XXX");
    if (rtmp.real() >= 0.0) {
        rtmp.real().write(buf, BUFLEN, QD_PRECISION);
    } else {
        rtmp.real().write(buf, BUFLEN, QD_PRECISION);
    }
    if (rtmp.imag() >= 0.0) {
        rtmp.imag().write(buf, BUFLEN, QD_PRECISION);
    } else {
        rtmp.imag().write(buf, BUFLEN, QD_PRECISION);
    }
    return;
}

qd_real mpf_randomnumber(qd_real dummy)
{
    qd_real mtmp;

    mtmp = qdrand();		//uniform random between [0,1] via lrand48
    mtmp = 2.0 * mtmp - 1.0;
    return mtmp;
}

qd_complex mpc_randomnumber(qd_complex dummy)
{
    qd_complex ctmp;
    qd_real mtmp1;
    qd_real mtmp2;

    mtmp1 = qdrand();
    mtmp1 = 2.0 * mtmp1 - 1.0;

    mtmp2 = qdrand();
    mtmp2 = 2.0 * mtmp2 - 1.0;

    ctmp.real() = mtmp1;
    ctmp.imag() = mtmp2;
    return ctmp;
}

void set_random_number(mpreal &a, qd_real & b)
{
    mpreal dummy;
    a = mpf_randomnumber(dummy);
    b = cast2qd_real(a);
}

void set_random_number(mpcomplex & a, qd_complex &b)
{
    mpcomplex dummy;
    a = mpc_randomnumber(dummy);
    b = cast2qd_complex(a);
}

void set_random_number1to2(mpreal &a, qd_real & b)
{
    mpreal dummy;
    a = mpf_randomnumber(dummy);
    if (a > 0.0) a = a + 1.0; else a = a - 1.0;
    b = cast2qd_real(a);
}

void set_random_number1to2(mpcomplex & a, qd_complex &b)
{
    mpcomplex dummy;
    double p, q;
    a = mpc_randomnumber(dummy);
    if (a.real() > 0.0) p = 1.0; else p = -1.0;
    if (a.imag() > 0.0) q = 1.0; else q = -1.0;
    a = a + complex<double>(p, q);
    b = cast2qd_complex(a);
}

#endif

#if defined ___MPLAPACK_BUILD_WITH_DD___
void printnum(dd_real rtmp)
{
    cout.precision(DD_PRECISION);
    if (rtmp >= 0.0) {
	cout << "+" << rtmp;
    } else {
	cout << rtmp;
    }
    return;
}

void printnum(dd_complex rtmp)
{
    cout.precision(DD_PRECISION);
    if (rtmp.real() >= 0.0) {
	cout << "+" << rtmp.real();
    } else {
	cout << rtmp.real();
    }
    if (rtmp.imag() >= 0.0) {
	cout << "+" << rtmp.imag() << "i";
    } else {
	cout << rtmp.imag() << "i";
    }
    return;
}

void sprintnum(char *buf, dd_real rtmp)
{
    rtmp.write(buf, BUFLEN, DD_PRECISION);
    return;
}

void sprintnum(char *buf, dd_complex rtmp)
{
printf("XXX");
    if (rtmp.real() >= 0.0) {
        rtmp.real().write(buf, BUFLEN, DD_PRECISION);
    } else {
        rtmp.real().write(buf, BUFLEN, DD_PRECISION);
    }
    if (rtmp.imag() >= 0.0) {
        rtmp.imag().write(buf, BUFLEN, DD_PRECISION);
    } else {
        rtmp.imag().write(buf, BUFLEN, DD_PRECISION);
    }
    return;
}

dd_real mpf_randomnumber(dd_real dummy)
{
    dd_real mtmp;

    mtmp = ddrand();		//uniform random between [0,1] via lrand48
    mtmp = 2.0 * mtmp - 1.0;
    return mtmp;
}

dd_complex mpc_randomnumber(dd_complex dummy)
{
    dd_complex ctmp;
    dd_real mtmp1;
    dd_real mtmp2;

    mtmp1 = ddrand();
    mtmp1 = 2.0 * mtmp1 - 1.0;

    mtmp2 = ddrand();
    mtmp2 = 2.0 * mtmp2 - 1.0;

    ctmp.real() = mtmp1;
    ctmp.imag() = mtmp2;
    return ctmp;
}

void set_random_number(mpreal &a, dd_real & b)
{
    mpreal dummy;
    a = mpf_randomnumber(dummy);
    b = cast2dd_real(a);
}

void set_random_number(mpcomplex & a, dd_complex &b)
{
    mpcomplex dummy;
    a = mpc_randomnumber(dummy);
    b = cast2dd_complex(a);
}

void set_random_number1to2(mpreal &a, dd_real & b)
{
    mpreal dummy;
    a = mpf_randomnumber(dummy);
    if (a > 0.0) a = a + 1.0; else a = a - 1.0;
    b = cast2dd_real(a);
}

void set_random_number1to2(mpcomplex & a, dd_complex &b)
{
    mpcomplex dummy;
    double p, q;
    a = mpc_randomnumber(dummy);
    if (a.real() > 0.0) p = 1.0; else p = -1.0;
    if (a.imag() > 0.0) q = 1.0; else q = -1.0;
    a = a + complex<double>(p, q);
    b = cast2dd_complex(a);
}
#endif

#if defined ___MPLAPACK_BUILD_WITH_DOUBLE___
double mpf_randomnumber(double dummy)
{
    double mtmp = drand48();
    return mtmp;
}

complex < double >mpc_randomnumber(complex < double >)
{
    std::complex < double >ctmp;
    ctmp.real(drand48());
    ctmp.imag(drand48());
    return ctmp;
}

void set_random_number(mpreal &a, double & b)
{
    mpreal dummy;
    a = mpf_randomnumber(dummy);
    b = a;
}

void set_random_number(mpcomplex & a, complex<double> &b)
{
    mpcomplex dummy;
    a = mpc_randomnumber(dummy);
    b = a;
}

void set_random_number1to2(mpreal &a, double & b)
{
    mpreal dummy;
    a = mpf_randomnumber(dummy);
    if (a > 0.0) a = a + 1.0; else a = a - 1.0;
    b = a;
}

void set_random_number1to2(mpcomplex & a, complex<double> &b)
{
    mpcomplex dummy;
    double p, q;
    a = mpc_randomnumber(dummy);
    if (a.real() > 0.0) p = 1.0; else p = -1.0;
    if (a.imag() > 0.0) q = 1.0; else q = -1.0;
    a = a + complex<double>(p, q);
    b = a;
}

#endif

#if defined ___MPLAPACK_BUILD_WITH__FLOAT64X___
_Float64x mpf_randomnumber(_Float64x dummy)
{
    _Float64x mtmp = drand48();
    return mtmp;
}

complex < _Float64x >mpc_randomnumber(complex < _Float64x >)
{
    std::complex < _Float64x >ctmp;
    ctmp.real(drand48());
    ctmp.imag(drand48());
    return ctmp;
}

void set_random_number(mpreal &a, _Float64x & b)
{
    mpreal dummy;
    a = mpf_randomnumber(dummy);
    b = a;
}

void set_random_number(mpcomplex & a, complex<_Float64x> &b)
{
    mpcomplex dummy;
    a = mpc_randomnumber(dummy);
    b = a;
}

void set_random_number1to2(mpreal &a, _Float64x & b)
{
    mpreal dummy;
    a = mpf_randomnumber(dummy);
    if (a > 0.0) a = a + 1.0; else a = a - 1.0;
    b = a;
}

void set_random_number1to2(mpcomplex & a, complex<_Float64x> &b)
{
    mpcomplex dummy;
    double p, q;
    a = mpc_randomnumber(dummy);
    if (a.real() > 0.0) p = 1.0; else p = -1.0;
    if (a.imag() > 0.0) q = 1.0; else q = -1.0;
    a = a + complex<double>(p, q);
    b = a;
}

#endif

#if defined ___MPLAPACK_BUILD_WITH__FLOAT128___ && defined _MPLAPACK_WANT_LIBQUADMATH_
void printnum(_Float128 rtmp)
{
    int width = 42;
    char buf[BUFLEN];
    int n = quadmath_snprintf (buf, sizeof buf, "%+-#*.35Qe", width, rtmp);
    if ((size_t) n < sizeof buf)
    printf ("%s", buf);
    return;
}

void printnum(std::complex<_Float128> rtmp)
{
    int width = 42, n;
    char buf[BUFLEN], buf2[BUFLEN];
    n = quadmath_snprintf (buf, sizeof buf, "%+-#*.35Qe", width, rtmp.real());
    if ((size_t) n < sizeof buf)
    printf ("%s", buf);
    n = quadmath_snprintf (buf2, sizeof buf2, "%+-#*.35Qe", width, rtmp.imag());
    if ((size_t) n < sizeof buf2)
    printf ("%s", buf2);
    printf ("i");
    return;
}

void sprintnum(char *buf, _Float128 rtmp)
{
    int width = 42;
    quadmath_snprintf (buf, BUFLEN, "%+-#*.35Qe", width, rtmp);
    return;
}

void sprintnum(char *buf, std::complex<_Float128> rtmp)
{
    int width = 42;
    quadmath_snprintf (buf, BUFLEN, "%+-#*.35Qe" "%+-#*.35Qe", width, rtmp.real(), rtmp.imag());
    return;
}
#endif

//native support of _Float128.
#if defined ___MPLAPACK_BUILD_WITH__FLOAT128___ && !defined _MPLAPACK__FLOAT128_IS_LONGDOUBLE_
void printnum(_Float128 rtmp)
{
    int width = 42;
    char buf[BUFLEN];
    strfromf128 (buf, sizeof buf, "%+-#*.35Qe", rtmp);
    printf ("%s", buf);
    return;
}

void printnum(std::complex<_Float128> rtmp)
{
    int width = 42, n;
    char buf[BUFLEN], buf2[BUFLEN];
    strfromf128 (buf, sizeof buf, "%+-#*.35Qe", rtmp.real());
    printf ("%s", buf);
    strfromf128 (buf2, sizeof buf2, "%+-#*.35Qe", rtmp.imag());
    printf ("%s", buf2);
    printf ("i");
    return;
}

void sprintnum(char *buf, _Float128 rtmp)
{
    int width = 42;
    strfromf128 (buf, BUFLEN, "%+-#*.35Qe", rtmp);
    return;
}

void sprintnum(char *buf, std::complex<_Float128> rtmp)
{
    int width = 42, n;
    char buf2[BUFLEN];
    strfromf128 (buf, sizeof buf, "%+-#*.35Qe", rtmp.real());
    strfromf128 (buf2, sizeof buf2, "%+-#*.35Qe", rtmp.imag());
    n = strlen(buf);  
    strncat(buf, buf2, n);
    return;
}

_Float128 mpf_randomnumber(_Float128 dummy)
{
    _Float128 mtmp;
    mtmp = lrand48();		//uniform random between [0,1] via lrand48
    mtmp += lrand48()*1e-16;	//uniform random between [0,1] via lrand48
    mtmp += lrand48()*1e-32;	//uniform random between [0,1] via lrand48
    mtmp = 2.0 * mtmp - 1.0;
    return mtmp;
}

std::complex<_Float128> mpc_randomnumber(std::complex<_Float128> dummy)
{
    complex<_Float128> ctmp;
    _Float128 mtmp1;
    _Float128 mtmp2;

    mtmp1 = lrand48();		//uniform random between [0,1] via lrand48
    mtmp1 += lrand48()*1e-16;	//uniform random between [0,1] via lrand48
    mtmp1 += lrand48()*1e-32;	//uniform random between [0,1] via lrand48
    mtmp1 = 2.0 * mtmp1 - 1.0;

    mtmp2 = lrand48();		//uniform random between [0,1] via lrand48
    mtmp2 += lrand48()*1e-16;	//uniform random between [0,1] via lrand48
    mtmp2 += lrand48()*1e-32;	//uniform random between [0,1] via lrand48
    mtmp2 = 2.0 * mtmp2 - 1.0;

    ctmp.real(mtmp1);
    ctmp.imag(mtmp2);
    return ctmp;
}

void set_random_number(mpreal &a, _Float128 & b)
{
    mpreal dummy;
    a = mpf_randomnumber(dummy);
    b = cast2_Float128(a);
}

void set_random_number(mpcomplex & a, std::complex<_Float128> &b)
{
    mpcomplex dummy;
    a = mpc_randomnumber(dummy);
    b.real(cast2_Float128(a.real()));
    b.imag(cast2_Float128(a.imag()));
}

void set_random_number1to2(mpreal &a, _Float128 & b)
{
    mpreal dummy;
    a = mpf_randomnumber(dummy);
    if (a > 0.0) a = a + 1.0; else a = a - 1.0;
    b = cast2_Float128(a);
}

void set_random_number1to2(mpcomplex & a, std::complex<_Float128> &b)
{
    mpcomplex dummy;
    _Float128 p, q;
    a = mpc_randomnumber(dummy);
    if (a.real() > 0.0) p = 1.0; else p = -1.0;
    if (a.imag() > 0.0) q = 1.0; else q = -1.0;
    a = a + std::complex<_Float128>(p, q);
    b.real(cast2_Float128(a.real()));
    b.imag(cast2_Float128(a.imag()));
}
#endif

//not depend on any mp libs.
void printnum(complex < double >ctmp)
{
    printf(P_FORMAT P_FORMAT "i", ctmp.real(), ctmp.imag());
    return;
}

void printnum(double rtmp)
{
    printf(P_FORMAT, rtmp);
    return;
}

void sprintnum(char *buf, double rtmp)
{
    snprintf(buf, BUFLEN, P_FORMAT, rtmp);
    return;
}

void sprintnum(char *buf, complex < double >ctmp)
{
    snprintf(buf, BUFLEN, P_FORMAT P_FORMAT "i", ctmp.real(), ctmp.imag());
    return;
}

void printnum(complex < long double >ctmp)
{
    printf(LP_FORMAT LP_FORMAT "i", ctmp.real(), ctmp.imag());
    return;
}

void printnum(long double rtmp)
{
    printf(LP_FORMAT, rtmp);
    return;
}

void sprintnum(char *buf, long double rtmp)
{
    snprintf(buf, BUFLEN, LP_FORMAT, rtmp);
    return;
}

void sprintnum(char *buf, complex < long double >ctmp)
{
    snprintf(buf, BUFLEN, LP_FORMAT LP_FORMAT "i", ctmp.real(), ctmp.imag());
    return;
}

REAL_REF infnorm(COMPLEX_REF * vec_ref, COMPLEX * vec, int len, int inc)
{
    COMPLEX_REF ctemp;
    REAL_REF inorm = 0.0;

    for (int i = 0; i < len * inc; i = i + inc) {
	ctemp = vec[i] - vec_ref[i];
	if (inorm < abs(ctemp)) {
	    inorm = abs(ctemp);
	}
    }
    return inorm;
}

REAL_REF infnorm(REAL_REF * vec_ref, REAL * vec, int len, int inc)
{
    REAL_REF ctemp;
    REAL_REF inorm = 0.0;

    for (int i = 0; i < len * inc; i = i + inc) {
	ctemp = vec[i] - vec_ref[i];
	if (inorm < abs(ctemp)) {
	    inorm = abs(ctemp);
	}
    }
    return inorm;
}

REAL_REF infnorm_mat(REAL_REF * vec_ref, REAL * vec, int n, int m, int lda)
{
    REAL_REF ctemp;
    REAL_REF inorm = 0.0;

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
	ctemp = vec[i + j * lda] - vec_ref[i + j * lda];
	if (inorm < abs(ctemp)) {
	    inorm = abs(ctemp);
	}
      }
    }
    return inorm;
}

REAL_REF reldiff(REAL_REF * vec_ref, REAL * vec, int len, int inc)
{
    REAL_REF ctemp1, ctemp2, ctemp3;
    REAL_REF inorm = 0.0;

    for (int i = 0; i < len * inc; i = i + inc) {
	ctemp1 = (REAL_REF)abs(vec[i]);
        ctemp2 = abs(vec_ref[i]);
	ctemp3 = max (ctemp1, ctemp2);
        ctemp1 = vec[i] - vec_ref[i];
	if (inorm < abs(ctemp1) / ctemp3 ) {
	    inorm = abs(ctemp1) / ctemp3;
	}
    }
    return inorm;
}

INTEGER_REF infnorm(INTEGER_REF * vec_ref, INTEGER * vec, int len, int inc)
{
    INTEGER_REF inorm = 0;
    INTEGER_REF temp = 0;
    for (int i = 0; i < len * inc; i = i + inc) {
	temp = vec[i] - vec_ref[i];
	if (inorm < abs(temp)) {
	    inorm = abs(temp);
	}
    }
    return inorm;
}
