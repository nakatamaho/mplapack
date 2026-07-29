//public domain
#include <mpblas_double.h>
#include <mplapack_double.h>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <algorithm>

#define DOUBLE_FORMAT "%+20.16e"
#define DOUBLE_SHORT_FORMAT "%+20.16e"

inline void printnum(double rtmp) { printf(DOUBLE_FORMAT, rtmp); }

// Matlab/Octave format
void printvec(double *a, int len) {
    double tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, double *a, int lda)
{
    double mtmp;

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
    double *a = new double[n * n];
    double *b = new double[n * n];
    double *alphar = new double[n];
    double *alphai = new double[n];
    double *beta = new double[n];
    double *vl = new double[1];
    double *vr = new double[1];
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
    double wk;
    Rggev("N", "N", n, a, lda, b, ldb, alphar, alphai, beta, vl, ldv, vr, ldv, &wk, lwork, info);
    lwork = castINTEGER_double(wk);
    double *work = new double[lwork];
    Rggev("N", "N", n, a, lda, b, ldb, alphar, alphai, beta, vl, ldv, vr, ldv, work, lwork, info);
    for (mplapackint i = 0; i < n; i++) {
        printf("alpha[%ld] = ", (long)i); printnum(alphar[i]); printf(" + "); printnum(alphai[i]); printf("i, beta = "); printnum(beta[i]);
        if (abs(beta[i]) <= Rlamch_double("E")) {
            printf(", lambda = Inf\n");
        } else {
            printf(", lambda = "); printnum(alphar[i] / beta[i]); printf(" + "); printnum(alphai[i] / beta[i]); printf("i\n");
        }
    }
    delete[] work;
    delete[] vr;
    delete[] vl;
    delete[] beta;
    delete[] alphai;
    delete[] alphar;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
