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
//taking from Collection of Matrices for Testing Computational Algorithms 1969 Robert T. Gregory, David L. Karney pp.30
int main()
{
    mplapackint n = 3;
    mplapackint lwork, info;

    mpcomplex *a = new mpcomplex[n * n];
    mplapackint *ipiv = new mplapackint[n];

//setting a matrix


    a[0 + 0 * n] = mpcomplex(1.0, 0.0);   a[0 + 1 * n] = mpcomplex(1.0, 2.0);    a[0 + 2 * n] = mpcomplex(2.0, 10.0);
    a[1 + 0 * n] = mpcomplex(1.0, 1.0);   a[1 + 1 * n] = mpcomplex(0.0, 3.0);    a[1 + 2 * n] = mpcomplex(-5.0, 14.0);
    a[2 + 0 * n] = mpcomplex(1.0, 1.0);   a[2 + 1 * n] = mpcomplex(0.0, 5.0);    a[2 + 2 * n] = mpcomplex(-8.0, 20.0);

    printf("a ="); printmat(n, n, a, n); printf("\n");

//work space query
    lwork = -1;
    mpcomplex *work = new mpcomplex[1];

    Cgetri(n, a, n, ipiv, work, lwork, info);
    lwork = castINTEGER_mpfr (work[0].real());
    delete[]work;
    work = new mpcomplex[std::max(1, (int) lwork)];

//inverse matrix
    Cgetrf(n, n, a, n, ipiv, info);
    Cgetri(n, a, n, ipiv, work, lwork, info);

    printf("ainv ="); printmat(n, n, a, n); printf("\n");
    delete[]work;
    delete[]ipiv;
    delete[]a;
}
