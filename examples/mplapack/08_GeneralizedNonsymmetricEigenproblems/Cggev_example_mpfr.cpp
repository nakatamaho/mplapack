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

inline void printnum(mpfr_class rtmp) { mpfr_printf(MPFR_FORMAT, mpfr_ptr(rtmp)); }
inline void printnum_short(mpfr_class rtmp) { mpfr_printf(MPFR_SHORT_FORMAT, mpfr_ptr(rtmp)); }
inline void printnum(mpc_class ctmp) {
    mpfr_class cre, cim;
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
    mplapackint n = 3, lda = n, ldb = n, ldv = 1, info, lwork = -1;
    mpc_class *a = new mpc_class[n * n];
    mpc_class *b = new mpc_class[n * n];
    mpc_class *alpha = new mpc_class[n];
    mpc_class *beta = new mpc_class[n];
    mpc_class *vl = new mpc_class[1];
    mpc_class *vr = new mpc_class[1];
    mpfr_class *rwork = new mpfr_class[8 * n];
    for (mplapackint i = 0; i < n * n; i++) {
        a[i] = mpc_class(0.0, 0.0);
        b[i] = mpc_class(0.0, 0.0);
    }
    a[0] = mpc_class(1.0, 1.0);
    a[4] = mpc_class(2.0, 0.0);
    a[8] = mpc_class(3.0, -1.0);
    b[0] = mpc_class(1.0, 0.0);
    b[4] = mpc_class(1.0, 0.0);
    b[8] = mpc_class(0.0, 0.0);
    mpc_class wk;
    Cggev("N", "N", n, a, lda, b, ldb, alpha, beta, vl, ldv, vr, ldv, &wk, lwork, rwork, info);
    lwork = castINTEGER_mpfr(wk.real());
    mpc_class *work = new mpc_class[lwork];
    Cggev("N", "N", n, a, lda, b, ldb, alpha, beta, vl, ldv, vr, ldv, work, lwork, rwork, info);
    for (mplapackint i = 0; i < n; i++) {
        printf("alpha[%ld] = ", (long)i); printnum(alpha[i]); printf(", beta = "); printnum(beta[i]);
        if (abs(beta[i]) <= Rlamch_mpfr("E"))
            printf(", lambda = Inf\n");
        else {
            printf(", lambda = "); printnum(alpha[i] / beta[i]); printf("\n");
        }
    }
    delete[] work;
    delete[] rwork;
    delete[] vr;
    delete[] vl;
    delete[] beta;
    delete[] alpha;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
