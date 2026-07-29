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

mpfr_class max_solution_error(mplapackint n, mpfr_class *x, mpfr_class *xexact) {
    mpfr_class err = 0;
    for (mplapackint i = 0; i < n; i++) {
        mpfr_class d = abs(x[i] - xexact[i]);
        if (err < d)
            err = d;
    }
    return err;
}
int main() {
    mplapackint n = 15, lda = n, ldb = n, info;
    mpfr_class *a = new mpfr_class[n * n];
    mpfr_class *b = new mpfr_class[n];
    mpfr_class *xexact = new mpfr_class[n];
    mplapackint *ipiv = new mplapackint[n];
    for (mplapackint j = 0; j < n; j++)
        xexact[j] = mpfr_class(j % 3 - 1);
    for (mplapackint i = 0; i < n; i++) {
        mpfr_class node = mpfr_class(i + 1);
        mpfr_class p = 1;
        for (mplapackint j = 0; j < n; j++) {
            a[i + j * lda] = p;
            p = p * node;
        }
    }
    for (mplapackint i = 0; i < n; i++) {
        b[i] = 0;
        for (mplapackint j = 0; j < n; j++)
            b[i] = b[i] + a[i + j * lda] * xexact[j];
    }
    Rgesv(n, (mplapackint)1, a, lda, ipiv, b, ldb, info);
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, b, xexact)); printf("\n");
    delete[] ipiv;
    delete[] xexact;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
