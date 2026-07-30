//public domain
#include <mpblas_double.h>
#include <mplapack_double.h>
#include <iostream>
#include <sstream>
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
bool select_none(double ar, double ai, double beta) {
    return false;
}
int main() {
    mplapackint n = 3, lda = n, ldb = n, ldv = 1, sdim, info, lwork = -1;
    double *a = new double[n * n];
    double *b = new double[n * n];
    double *alphar = new double[n];
    double *alphai = new double[n];
    double *beta = new double[n];
    double *vsl = new double[1];
    double *vsr = new double[1];
    bool *bwork = new bool[n];
    for (mplapackint i = 0; i < n * n; i++) {
        a[i] = 0.0;
        b[i] = 0.0;
    }
    a[0] = 1.0;
    a[4] = 2.0;
    a[8] = 3.0;
    b[0] = 1.0;
    b[4] = 1.0;
    b[8] = 0.0;
    double wk;
    Rgges("N", "N", "N", select_none, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldv, vsr, ldv, &wk, lwork, bwork, info);
    lwork = castINTEGER_double(wk);
    double *work = new double[lwork];
    Rgges("N", "N", "N", select_none, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldv, vsr, ldv, work, lwork, bwork, info);
    printf("S = "); printmat(n, n, a, lda); printf("\n");
    printf("T = "); printmat(n, n, b, ldb); printf("\n");
    for (mplapackint i = 0; i < n; i++) {
        printf("lambda[%ld] = ", (long)i);
        if (abs(beta[i]) <= Rlamch_double("E"))
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
