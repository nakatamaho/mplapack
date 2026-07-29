//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_dd.h>
#include <mplapack_dd.h>

#define DD_PRECISION_SHORT 16

inline void printnum(dd_real rtmp) {
    std::cout.precision(DD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

//Matlab/Octave format
void printvec(dd_real *a, int len) {
    dd_real tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, dd_real * a, int lda)
{
    dd_real mtmp;
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
    dd_real *a = new dd_real[n * n];
    dd_real *b = new dd_real[n * n];
    dd_real *alphar = new dd_real[n];
    dd_real *alphai = new dd_real[n];
    dd_real *beta = new dd_real[n];
    dd_real *vl = new dd_real[1];
    dd_real *vr = new dd_real[1];
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
    dd_real wk;
    Rggev("N", "N", n, a, lda, b, ldb, alphar, alphai, beta, vl, ldv, vr, ldv, &wk, lwork, info);
    lwork = castINTEGER_dd(wk);
    dd_real *work = new dd_real[lwork];
    Rggev("N", "N", n, a, lda, b, ldb, alphar, alphai, beta, vl, ldv, vr, ldv, work, lwork, info);
    for (mplapackint i = 0; i < n; i++) {
        printf("alpha[%ld] = ", (long)i); printnum(alphar[i]); printf(" + "); printnum(alphai[i]); printf("i, beta = "); printnum(beta[i]);
        if (abs(beta[i]) <= Rlamch_dd("E")) {
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
