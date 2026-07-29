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
int main() {
    mplapackint m = 2, n = 3, p = 2, k, l, lda = m, ldb = p, ldu = m, ldv = p, ldq = n, info, lwork = -1;
    double *a = new double[lda * n];
    double *b = new double[ldb * n];
    double *alpha = new double[n];
    double *beta = new double[n];
    double *u = new double[ldu * m];
    double *v = new double[ldv * p];
    double *q = new double[ldq * n];
    mplapackint *iwork = new mplapackint[n];
    for (mplapackint i = 0; i < lda * n; i++)
        a[i] = 0;
    for (mplapackint i = 0; i < ldb * n; i++)
        b[i] = 0;
    a[0] = 1;
    a[1 + lda] = 2;
    a[0 + 2 * lda] = 1;
    b[0] = 1;
    b[1 + ldb] = 3;
    b[0 + 2 * ldb] = 1;
    double wk;
    Rggsvd3("U", "V", "Q", m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, &wk, lwork, iwork, info);
    lwork = castINTEGER_double(wk);
    double *work = new double[lwork];
    Rggsvd3("U", "V", "Q", m, n, p, k, l, a, lda, b, ldb, alpha, beta, u, ldu, v, ldv, q, ldq, work, lwork, iwork, info);
    printf("k = %ld, l = %ld\n", (long)k, (long)l);
    printf("alpha = "); printvec(alpha, n); printf("\n");
    printf("beta = "); printvec(beta, n); printf("\n");
    for (mplapackint i = k; i < k + l; i++) {
        printf("gsv[%ld] = ", (long)i);
        if (abs(beta[i]) <= Rlamch_double("E"))
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
