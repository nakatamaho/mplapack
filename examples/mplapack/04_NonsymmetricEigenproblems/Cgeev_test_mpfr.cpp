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
bool rselect(mpreal ar, mpreal ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 4;
    mpcomplex *a = new mpcomplex[n * n];
    mpcomplex *w = new mpcomplex[n];
    mpcomplex *vl = new mpcomplex[n * n];
    mpcomplex *vr = new mpcomplex[n * n];
    mplapackint lwork = 4 * n;
    mpcomplex *work = new mpcomplex[lwork];    
    mpreal *rwork = new mpreal[lwork];
    mplapackint info;
    // setting A matrix
    //# Example 6.5 "Collection of Matrices for Testing Computational Algorithms", Robert T. Gregory, David L. Karney    
    a[0 + 0 * n] = mpcomplex(5.0, 9.0); a[0 + 1 * n] = mpcomplex(5.0, 5.0);   a[0 + 2 * n] = mpcomplex(-6.0, -6.0); a[0 + 3 * n] = mpcomplex(-7.0, -7.0);
    a[1 + 0 * n] = mpcomplex(3.0, 3.0); a[1 + 1 * n] = mpcomplex(6.0, 10.0);  a[1 + 2 * n] = mpcomplex(-5.0, -5.0); a[1 + 3 * n] = mpcomplex(-6.0, -6.0);
    a[2 + 0 * n] = mpcomplex(2.0, 2.0); a[2 + 1 * n] = mpcomplex(3.0, 3.0);   a[2 + 2 * n] = mpcomplex(-1.0,  3.0); a[2 + 3 * n] = mpcomplex(-5.0, -5.0);
    a[3 + 0 * n] = mpcomplex(1.0, 1.0); a[3 + 1 * n] = mpcomplex(2.0, 2.0);   a[3 + 2 * n] = mpcomplex(-3.0, -3.0); a[3 + 3 * n] = mpcomplex(0.0, 4.0); 

    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Cgeev("V", "V", n, a, n, w, vl, n, vr, n, work, lwork, rwork, info);
    printf("lambda ="); printvec(w,n); printf("\n");    
    printf("vr ="); printmat(n,n,vr,n); printf("\n");    

    delete[] rwork;
    delete[] work;
    delete[] vr;
    delete[] vl;
    delete[] w;
    delete[] a;
}
