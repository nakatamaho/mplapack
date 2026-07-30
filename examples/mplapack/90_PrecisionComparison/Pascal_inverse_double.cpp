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
double binom(mplapackint n, mplapackint k) {
    double r = 1.0;
    for (mplapackint i = 1; i <= k; i++)
        r = r * double(n - k + i) / double(i);
    return r;
}
double nearest_integer_error(double x) {
    mplapackint nearest = castINTEGER_double(x >= double(0.0) ? x + double(0.5) : x - double(0.5));
    return abs(x - double(nearest));
}
int main() {
    mplapackint n = 8, lda = n, info, lwork = -1;
    double *a = new double[n * n];
    mplapackint *ipiv = new mplapackint[n];
    for (mplapackint j = 0; j < n; j++)
        for (mplapackint i = 0; i < n; i++)
            a[i + j * lda] = binom(i + j, i);
    Rgetrf(n, n, a, lda, ipiv, info);
    double wk;
    if (info == 0)
        Rgetri(n, a, lda, ipiv, &wk, lwork, info);
    lwork = castINTEGER_double(wk);
    double *work = new double[lwork];
    if (info == 0)
        Rgetri(n, a, lda, ipiv, work, lwork, info);
    double err = 0.0;
    for (mplapackint i = 0; i < n * n; i++) {
        double d = nearest_integer_error(a[i]);
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
