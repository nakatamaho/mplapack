//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_qd.h>
#include <mplapack_qd.h>

#define QD_PRECISION_SHORT 16

inline void printnum(qd_real rtmp) {
    std::cout.precision(QD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

//Matlab/Octave format
void printvec(qd_real *a, int len) {
    qd_real tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, qd_real * a, int lda)
{
    qd_real mtmp;
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
qd_real binom(mplapackint n, mplapackint k) {
    qd_real r = 1.0;
    for (mplapackint i = 1; i <= k; i++)
        r = r * qd_real((double)(n - k + i)) / qd_real((double)i);
    return r;
}
qd_real nearest_integer_error(qd_real x) {
    mplapackint nearest = castINTEGER_qd(x >= qd_real(0.0) ? x + qd_real(0.5) : x - qd_real(0.5));
    return abs(x - qd_real((double)nearest));
}
int main() {
    mplapackint n = 8, lda = n, info, lwork = -1;
    qd_real *a = new qd_real[n * n];
    mplapackint *ipiv = new mplapackint[n];
    for (mplapackint j = 0; j < n; j++)
        for (mplapackint i = 0; i < n; i++)
            a[i + j * lda] = binom(i + j, i);
    Rgetrf(n, n, a, lda, ipiv, info);
    qd_real wk;
    if (info == 0)
        Rgetri(n, a, lda, ipiv, &wk, lwork, info);
    lwork = castINTEGER_qd(wk);
    qd_real *work = new qd_real[lwork];
    if (info == 0)
        Rgetri(n, a, lda, ipiv, work, lwork, info);
    qd_real err = 0.0;
    for (mplapackint i = 0; i < n * n; i++) {
        qd_real d = nearest_integer_error(a[i]);
        if (err < d)
            err = d;
    }
    printf("P inverse = "); printmat(n, n, a, lda); printf("\n");
    printf("max distance to integer = "); printnum(err); printf("\n");
    delete[] work;
    delete[] ipiv;
    delete[] a;
    return info != 0 ? 1 : 0;
}
