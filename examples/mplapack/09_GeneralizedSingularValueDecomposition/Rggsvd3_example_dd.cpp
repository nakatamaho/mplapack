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
    mplapackint m = 2, n = 3, p = 2, k, l, lda = m, ldb = p, ldu = m, ldv = p, ldq = n, info, lwork = -1;
    dd_real *a = new dd_real[lda * n];
    dd_real *b = new dd_real[ldb * n];
    dd_real *alpha = new dd_real[n];
    dd_real *beta = new dd_real[n];
    dd_real *u = new dd_real[ldu * m];
    dd_real *v = new dd_real[ldv * p];
    dd_real *q = new dd_real[ldq * n];
    mplapackint *iwork = new mplapackint[n];
    for (mplapackint i = 0; i < lda * n; i++)
        a[i] = 0.0;
    for (mplapackint i = 0; i < ldb * n; i++)
        b[i] = 0.0;
    a[0] = 1.0;
    a[1 + lda] = 2.0;
    a[0 + 2 * lda] = 1.0;
    b[0] = 1.0;
    b[1 + ldb] = 3.0;
    b[0 + 2 * ldb] = 1.0;
    dd_real wk;
    Rggsvd3("U", "V", "Q", m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, &wk, lwork, iwork, info);
    lwork = castINTEGER_dd(wk);
    dd_real *work = new dd_real[lwork];
    Rggsvd3("U", "V", "Q", m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, iwork, info);
    printf("k = %ld, l = %ld\n", (long)k, (long)l);
    printf("alpha = "); printvec(alpha, n); printf("\n");
    printf("beta = "); printvec(beta, n); printf("\n");
    for (mplapackint i = k; i < k + l; i++) {
        printf("gsv[%ld] = ", (long)i);
        if (abs(beta[i]) <= Rlamch_dd("E"))
            printf("Inf\n");
        else {
            printnum(alpha[i] / beta[i]);
            printf("\n");
        }
    }
    delete[] work;
    delete[] iwork;
    delete[] q;
    delete[] v;
    delete[] u;
    delete[] beta;
    delete[] alpha;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
