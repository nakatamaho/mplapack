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
td_real max_solution_error(mplapackint n, td_real *x, td_real *xexact) {
    td_real err = 0.0;
    for (mplapackint i = 0; i < n; i++) {
        td_real d = abs(x[i] - xexact[i]);
        if (err < d)
            err = d;
    }
    return err;
}
int main() {
    mplapackint n = 15, lda = n, ldb = n, info;
    td_real *a = new td_real[n * n];
    td_real *b = new td_real[n];
    td_real *xexact = new td_real[n];
    mplapackint *ipiv = new mplapackint[n];
    for (mplapackint j = 0; j < n; j++)
        xexact[j] = td_real((double)(j % 3 - 1));
    for (mplapackint i = 0; i < n; i++) {
        td_real node = td_real((double)(i + 1));
        td_real p = 1.0;
        for (mplapackint j = 0; j < n; j++) {
            a[i + j * lda] = p;
            p = p * node;
        }
    }
    for (mplapackint i = 0; i < n; i++) {
        b[i] = 0.0;
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
