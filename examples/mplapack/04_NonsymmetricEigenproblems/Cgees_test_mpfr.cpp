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
bool cselect(mpcomplex a) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 4;

    mpcomplex *a = new mpcomplex[n * n];
    mplapackint sdim = 0;
    mplapackint lwork = 2 * n;
    mpcomplex *w = new mpcomplex[n];
    mpcomplex *vs = new mpcomplex[n * n];
    mpcomplex *work = new mpcomplex[lwork];
    mpreal *rwork = new mpreal[n];
    bool bwork[n];
    mplapackint info;

    // setting A matrix
    a[0 + 0 * n] = mpcomplex(5.0,  9.0); a[0 + 1 * n] = mpcomplex(5.0, 5.0);   a[0 + 2 * n] = mpcomplex(-6.0, -6.0); a[0 + 3 * n] = mpcomplex(-7.0,-7.0);
    a[1 + 0 * n] = mpcomplex(3.0,  3.0); a[1 + 1 * n] = mpcomplex(6.0,10.0);   a[1 + 2 * n] = mpcomplex(-5.0, -5.0); a[1 + 3 * n] = mpcomplex(-6.0,-6.0);
    a[2 + 0 * n] = mpcomplex(2.0,  2.0); a[2 + 1 * n] = mpcomplex(3.0, 3.0);   a[2 + 2 * n] = mpcomplex(-1.0, 3.0);  a[2 + 3 * n] = mpcomplex(-5.0,-5.0);
    a[3 + 0 * n] = mpcomplex(1.0,  1.0); a[3 + 1 * n] = mpcomplex(2.0, 2.0);   a[3 + 2 * n] = mpcomplex(-3.0,-3.0);  a[3 + 3 * n] = mpcomplex(0.0, 4.0);

    printf("# Ex. 6.5 p. 116, Collection of Matrices for Testing Computational Algorithms, Robert T. Gregory, David L. Karney\n");
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Cgees("V", "S", cselect, n, a, n, sdim, w, vs, n, work, lwork, rwork, bwork, info);
    printf("w ="); printvec(w, n); printf("\n");
    printf("vs ="); printmat(n, n, vs, n); printf("\n");
    printf("t ="); printmat(n, n, a, n); printf("\n");
    printf("vs*t*vs'\n");
    printf("eig(a)\n");

    delete[] rwork;
    delete[] work;
    delete[] vs;
    delete[] w;
    delete[] a;
}
