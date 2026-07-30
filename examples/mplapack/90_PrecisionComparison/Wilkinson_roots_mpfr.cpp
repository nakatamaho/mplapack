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

// Matlab/Octave format
void printvec(mpfr_class *a, int len) {
    mpfr_class tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mpfr_class *a, int lda) {
    mpfr_class mtmp;
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

void sort_real(mplapackint n, mpfr_class *x) {
    for (mplapackint i = 0; i < n; i++)
        for (mplapackint j = i + 1; j < n; j++)
            if (x[j] < x[i]) {
                mpfr_class t = x[i];
                x[i] = x[j];
                x[j] = t;
            }
}
int main() {
    mplapackint n = 20, lda = n, ldv = 1, info, lwork = -1;
    mpfr_class *coef = new mpfr_class[n + 1];
    for (mplapackint i = 0; i <= n; i++)
        coef[i] = 0.0;
    coef[0] = 1.0;
    for (mplapackint k = 1; k <= n; k++) {
        for (mplapackint j = k; j >= 1; j--)
            coef[j] = coef[j] - mpfr_class((double)k) * coef[j - 1];
    }
    mpfr_class *a = new mpfr_class[n * n];
    for (mplapackint i = 0; i < n * n; i++)
        a[i] = 0.0;
    for (mplapackint i = 1; i < n; i++)
        a[i + (i - 1) * lda] = 1.0;
    for (mplapackint j = 0; j < n; j++)
        a[j + (n - 1) * lda] = -coef[n - j] / coef[0];
    mpfr_class *wr = new mpfr_class[n];
    mpfr_class *wi = new mpfr_class[n];
    mpfr_class *vl = new mpfr_class[1];
    mpfr_class *vr = new mpfr_class[1];
    mpfr_class wk;
    Rgeev("N", "N", n, a, lda, wr, wi, vl, ldv, vr, ldv, &wk, lwork, info);
    lwork = castINTEGER_mpfr(wk);
    mpfr_class *work = new mpfr_class[lwork];
    Rgeev("N", "N", n, a, lda, wr, wi, vl, ldv, vr, ldv, work, lwork, info);
    sort_real(n, wr);
    mpfr_class maxerr = 0.0;
    for (mplapackint i = 0; i < n; i++) {
        mpfr_class err = abs(wr[i] - mpfr_class((double)(i + 1)));
        if (maxerr < err)
            maxerr = err;
        printf("root[%ld] = ", (long)i); printnum(wr[i]); printf(", error = "); printnum(err); printf("\n");
    }
    printf("max root error = "); printnum(maxerr); printf("\n");
    delete[] work;
    delete[] vr;
    delete[] vl;
    delete[] wi;
    delete[] wr;
    delete[] a;
    delete[] coef;
    return info != 0 ? 1 : 0;
}
