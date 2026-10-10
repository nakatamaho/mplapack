//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_td.h>
#include <mplapack_td.h>

#define TD_PRECISION_SHORT 16

inline void printnum(td_real rtmp) {
    std::cout.precision(TD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

//Matlab/Octave format
void printvec(td_real *a, int len) {
    td_real tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, td_real * a, int lda)
{
    td_real mtmp;
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
    td_real *a = new td_real[n * n];
    td_real *b = new td_real[n * n];
    td_real *alphar = new td_real[n];
    td_real *alphai = new td_real[n];
    td_real *beta = new td_real[n];
    td_real *vl = new td_real[1];
    td_real *vr = new td_real[1];
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
    td_real wk;
    Rggev("N", "N", n, a, lda, b, ldb, alphar, alphai, beta, vl, ldv, vr, ldv, &wk, lwork, info);
    lwork = castINTEGER_td(wk);
    td_real *work = new td_real[lwork];
    Rggev("N", "N", n, a, lda, b, ldb, alphar, alphai, beta, vl, ldv, vr, ldv, work, lwork, info);
    for (mplapackint i = 0; i < n; i++) {
        printf("alpha[%ld] = ", (long)i); printnum(alphar[i]); printf(" + "); printnum(alphai[i]); printf("i, beta = "); printnum(beta[i]);
        if (abs(beta[i]) <= Rlamch_td("E")) {
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
