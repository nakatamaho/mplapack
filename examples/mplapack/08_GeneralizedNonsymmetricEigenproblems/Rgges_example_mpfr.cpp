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

// Matlab/Octave format
void printvec(mpreal *a, int len) {
    mpreal tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mpreal *a, int lda) {
    mpreal mtmp;
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

bool select_none(mpreal ar, mpreal ai, mpreal beta) {
    return false;
}
int main() {
    mplapackint n = 3, lda = n, ldb = n, ldv = 1, sdim, info, lwork = -1;
    mpreal *a = new mpreal[n * n];
    mpreal *b = new mpreal[n * n];
    mpreal *alphar = new mpreal[n];
    mpreal *alphai = new mpreal[n];
    mpreal *beta = new mpreal[n];
    mpreal *vsl = new mpreal[1];
    mpreal *vsr = new mpreal[1];
    bool *bwork = new bool[n];
    for (mplapackint i = 0; i < n * n; i++) {
        a[i] = 0;
        b[i] = 0;
    }
    a[0] = 1;
    a[4] = 2;
    a[8] = 3;
    b[0] = 1;
    b[4] = 1;
    b[8] = 0;
    mpreal wk;
    Rgges("N", "N", "N", select_none, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldv, vsr, ldv, &wk, lwork, bwork, info);
    lwork = castINTEGER_mpfr(wk);
    mpreal *work = new mpreal[lwork];
    Rgges("N", "N", "N", select_none, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldv, vsr, ldv, work, lwork, bwork, info);
    printf("S = "); printmat(n, n, a, lda); printf("\n");
    printf("T = "); printmat(n, n, b, ldb); printf("\n");
    for (mplapackint i = 0; i < n; i++) {
        printf("lambda[%ld] = ", (long)i);
        if (abs(beta[i]) <= Rlamch_mpfr("E"))
            printf("Inf\n");
        else {
            printnum(alphar[i] / beta[i]);
            printf(" + "); printnum(alphai[i] / beta[i]); printf("i\n");
        }
    }
    delete[] work;
    delete[] bwork;
    delete[] vsr;
    delete[] vsl;
    delete[] beta;
    delete[] alphai;
    delete[] alphar;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
