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
dd_real binom(mplapackint n, mplapackint k) {
    dd_real r = 1;
    for (mplapackint i = 1; i <= k; i++)
        r = r * dd_real(n - k + i) / dd_real(i);
    return r;
}
dd_real nearest_integer_error(dd_real x) {
    mplapackint nearest = castINTEGER_dd(x >= dd_real(0.0) ? x + dd_real(0.5) : x - dd_real(0.5));
    return abs(x - dd_real(nearest));
}
int main() {
    mplapackint n = 8, lda = n, info, lwork = -1;
    dd_real *a = new dd_real[n * n];
    mplapackint *ipiv = new mplapackint[n];
    for (mplapackint j = 0; j < n; j++)
        for (mplapackint i = 0; i < n; i++)
            a[i + j * lda] = binom(i + j, i);
    Rgetrf(n, n, a, lda, ipiv, info);
    dd_real wk;
    if (info == 0)
        Rgetri(n, a, lda, ipiv, &wk, lwork, info);
    lwork = castINTEGER_dd(wk);
    dd_real *work = new dd_real[lwork];
    if (info == 0)
        Rgetri(n, a, lda, ipiv, work, lwork, info);
    dd_real err = 0;
    for (mplapackint i = 0; i < n * n; i++) {
        dd_real d = nearest_integer_error(a[i]);
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
