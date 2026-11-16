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
bool select_none(td_real ar, td_real ai, td_real beta) {
    return false;
}
int main() {
    mplapackint n = 3, lda = n, ldb = n, ldv = 1, sdim, info, lwork = -1;
    td_real *a = new td_real[n * n];
    td_real *b = new td_real[n * n];
    td_real *alphar = new td_real[n];
    td_real *alphai = new td_real[n];
    td_real *beta = new td_real[n];
    td_real *vsl = new td_real[1];
    td_real *vsr = new td_real[1];
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
    td_real wk;
    Rgges("N", "N", "N", select_none, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldv, vsr, ldv, &wk, lwork, bwork, info);
    lwork = castINTEGER_td(wk);
    td_real *work = new td_real[lwork];
    Rgges("N", "N", "N", select_none, n, a, lda, b, ldb, sdim, alphar, alphai, beta, vsl, ldv, vsr, ldv, work, lwork, bwork, info);
    printf("S = "); printmat(n, n, a, lda); printf("\n");
    printf("T = "); printmat(n, n, b, ldb); printf("\n");
    for (mplapackint i = 0; i < n; i++) {
        printf("lambda[%ld] = ", (long)i);
        if (abs(beta[i]) <= Rlamch_td("E"))
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
