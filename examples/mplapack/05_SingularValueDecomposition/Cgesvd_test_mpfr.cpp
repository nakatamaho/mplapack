//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_mpfr.h>
#include <mplapack_mpfr.h>

#define MPFR_FORMAT "%+68.64Re"
#define MPFR_SHORT_FORMAT "%+20.16Re"

inline void printnum(mpreal rtmp) { mpfr_printf(MPFR_FORMAT, mpfr_ptr(rtmp)); }
inline void printnum_short(mpreal rtmp) { mpfr_printf(MPFR_SHORT_FORMAT, mpfr_ptr(rtmp)); }
inline void printnum(mpcomplex ctmp) {
    mpreal cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_printf(MPFR_SHORT_FORMAT MPFR_SHORT_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
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
int main() {
    mplapackint n = 4;
    mplapackint m = 4;

    mpcomplex *a = new mpcomplex[m * n];
    mpreal *s = new mpreal[std::min(m, n)];
    mpcomplex *u = new mpcomplex[m * m];
    mpcomplex *vt = new mpcomplex[n * n];
    mplapackint lwork = std::max((mplapackint)1, 2 * std::min(m, n) + std::max(m, n));
    mpcomplex *work = new mpcomplex[lwork];
    mpreal *rwork = new mpreal[5 * std::min(m, n)];
    mplapackint info;

    // setting A matrix
    a[0 + 0 * n] = mpcomplex(0.9, -1.0); a[0 + 1 * n] = mpcomplex(20.0, -2.25);  a[0 + 2 * n] = mpcomplex(1.75, -0.5);  a[0 + 3 * n] = mpcomplex(0.0, 0.5);
    a[1 + 0 * n] = mpcomplex(8.0,-2.25); a[1 + 1 * n] = mpcomplex(-0.25, 0.0);   a[1 + 2 * n] = mpcomplex(1.25, -0.25); a[1 + 3 * n] = mpcomplex(-3.75, 0.0);
    a[2 + 0 * n] = mpcomplex(-1.75,0.0); a[2 + 1 * n] = mpcomplex(-80.0,  1.25); a[2 + 2 * n] = mpcomplex(1.5, 0.0);    a[2 + 3 * n] = mpcomplex(30.0, 2.25);
    a[3 + 0 * n] = mpcomplex(3.0, 0.25); a[3 + 1 * n] = mpcomplex(1.75, 0.0);    a[3 + 2 * n] = mpcomplex(0.0, 2.25);   a[3 + 3 * n] = mpcomplex(-0.25, -80.0);
    
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(m, n, a, m); printf("\n");
    Cgesvd("A", "A", m, n, a, m, s, u, m, vt, n, work, lwork, rwork, info);
    printf("s="); printvec(s, std::min(m, n)); printf("\n");
    if (m < n)
        printf("padding=zeros(%d, %d-%d)\n", (int)m, (int)n, (int)m);
    if (n < m)
        printf("padding=zeros(%d-%d,%d)\n", (int)m, (int)n, (int)n);
    printf("u ="); printmat(m, m, u, m); printf("\n");
    printf("vt ="); printmat(n, n, vt, n); printf("\n");
    printf("svd(a)\n");
    if (m < n)
        printf("sigma=[diag(s) padding] \n");
    if (n < m)
        printf("sigma=[diag(s); padding] \n");
    if (n == m)
        printf("sigma=[diag(s)] \n");
    printf("sigma \n");
    printf("u * sigma  * vt\n");
    delete[] rwork;
    delete[] work;
    delete[] vt;
    delete[] u;
    delete[] s;
    delete[] a;
}
