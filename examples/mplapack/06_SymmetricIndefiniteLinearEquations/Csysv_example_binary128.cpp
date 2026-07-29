//public domain
#include <mpblas_binary128.h>
#include <mplapack_binary128.h>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <algorithm>

#define BUFLEN 1024

void printnum(mplapack_binary128_t rtmp)
{
    int width = 42;
    char buf[BUFLEN];
#if MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_STRFROMF128
    strfromf128(buf, sizeof(buf), "%.35e", rtmp);
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_QUADMATH_SNPRINTF
    int n = quadmath_snprintf (buf, sizeof buf, "%*.35Qe", width, rtmp);
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_SNPRINTF_LDBL
    snprintf (buf, sizeof buf, "%.35Le", rtmp);
#else
    #error "unsupported binary128 type"
#endif
    if (rtmp >= 0.0)
        printf ("+%s", buf);
    else
        printf ("%s", buf);
    return;
}

void printnum(std::complex<mplapack_binary128_t> rtmp)
{
    int width = 42;
    char buf[BUFLEN];
#if MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_STRFROMF128
    strfromf128(buf, sizeof(buf), "%.35e", rtmp.real());
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_QUADMATH_SNPRINTF
    int n = quadmath_snprintf (buf, sizeof buf, "%*.35Qe", width, rtmp.real());
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_SNPRINTF_LDBL
    snprintf (buf, sizeof buf, "%.35Le", rtmp.real());
#else
    #error "unsupported binary128 type"
#endif
    if (rtmp.real() >= 0.0)
        printf ("+%s", buf);
    else
        printf ("%s", buf);
#if MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_STRFROMF128
    strfromf128(buf, sizeof(buf), "%.35e", rtmp.imag());
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_QUADMATH_SNPRINTF
    n = quadmath_snprintf (buf, sizeof buf, "%*.35Qe", width, rtmp.imag());
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_SNPRINTF_LDBL
    snprintf (buf, sizeof buf, "%.35Le", rtmp.imag());
#else
    #error "unsupported binary128 type"
#endif
    if (rtmp.imag() >= 0.0)
        printf ("+%si", buf);
    else
        printf ("%si", buf);
    return;
}

//Matlab/Octave format
template <class X> void printvec(X *a, int len) {
    X tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

template <class X> void printmat(int n, int m, X *a, int lda)
{
    X mtmp;

    printf("[ ");
    for (int i = 0; i < n; i++) {
        printf("[ ");
        for (int j = 0; j < m; j++) {
            mtmp = a[i + j * lda];
            printnum(mtmp);
            if (j < m - 1)
                printf(", ");
        }
        if (i < n - 1)
            printf("]; ");
        else
            printf("] ");
    }
    printf("]");
}
mplapack_binary128_t max_solution_error(mplapackint n, mplapackint nrhs, std::complex<mplapack_binary128_t> *x, mplapackint ldx, std::complex<mplapack_binary128_t> *xexact, mplapackint ldxexact) {
    mplapack_binary128_t err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            mplapack_binary128_t d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

mplapack_binary128_t max_residual(mplapackint m, mplapackint n, mplapackint nrhs, std::complex<mplapack_binary128_t> *a, mplapackint lda, std::complex<mplapack_binary128_t> *x, mplapackint ldx, std::complex<mplapack_binary128_t> *b, mplapackint ldb) {
    mplapack_binary128_t err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            std::complex<mplapack_binary128_t> s = std::complex<mplapack_binary128_t>(0.0, 0.0);
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            mplapack_binary128_t d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

int main() {
    mplapackint n = 2, nrhs = 1, lda = n, ldb = n, info, lwork = -1;
    std::complex<mplapack_binary128_t> *a = new std::complex<mplapack_binary128_t>[n * n];
    std::complex<mplapack_binary128_t> *aorg = new std::complex<mplapack_binary128_t>[n * n];
    std::complex<mplapack_binary128_t> *b = new std::complex<mplapack_binary128_t>[n];
    std::complex<mplapack_binary128_t> *borg = new std::complex<mplapack_binary128_t>[n];
    std::complex<mplapack_binary128_t> *xexact = new std::complex<mplapack_binary128_t>[n];
    mplapackint *ipiv = new mplapackint[n];
    a[0] = std::complex<mplapack_binary128_t>(2.0, 0.0);
    a[1] = std::complex<mplapack_binary128_t>(1.0, 1.0);
    a[2] = std::complex<mplapack_binary128_t>(1.0, 1.0);
    a[3] = std::complex<mplapack_binary128_t>(-3.0, 0.0);
    xexact[0] = std::complex<mplapack_binary128_t>(1.0, -1.0);
    xexact[1] = std::complex<mplapack_binary128_t>(2.0, 1.0);
    for (mplapackint i = 0; i < n * n; i++)
        aorg[i] = a[i];
    for (mplapackint i = 0; i < n; i++) {
        b[i] = std::complex<mplapack_binary128_t>(0.0, 0.0);
        for (mplapackint k = 0; k < n; k++)
            b[i] = b[i] + aorg[i + k * lda] * xexact[k];
        borg[i] = b[i];
    }
    std::complex<mplapack_binary128_t> wk;
    Csysv("U", n, nrhs, a, lda, ipiv, b, ldb, &wk, lwork, info);
    lwork = castINTEGER_binary128(wk.real());
    std::complex<mplapack_binary128_t> *work = new std::complex<mplapack_binary128_t>[lwork];
    Csysv("U", n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
    printf("A = "); printmat(n, n, aorg, lda); printf("\n");
    printf("x = "); printvec(b, n); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, b, ldb, xexact, n)); printf("\n");
    printf("max residual = "); printnum(max_residual(n, n, nrhs, aorg, lda, b, ldb, borg, ldb)); printf("\n");
    delete[] work;
    delete[] ipiv;
    delete[] xexact;
    delete[] borg;
    delete[] b;
    delete[] aorg;
    delete[] a;
    return info != 0 ? 1 : 0;
}
