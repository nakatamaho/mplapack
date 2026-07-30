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
void sort_real(mplapackint n, double *x) {
    for (mplapackint i = 0; i < n; i++)
        for (mplapackint j = i + 1; j < n; j++)
            if (x[j] < x[i]) {
                double t = x[i];
                x[i] = x[j];
                x[j] = t;
            }
}
int main() {
    mplapackint n = 20, lda = n, ldv = 1, info, lwork = -1;
    double *coef = new double[n + 1];
    for (mplapackint i = 0; i <= n; i++)
        coef[i] = 0.0;
    coef[0] = 1.0;
    for (mplapackint k = 1; k <= n; k++) {
        for (mplapackint j = k; j >= 1; j--)
            coef[j] = coef[j] - double(k) * coef[j - 1];
    }
    double *a = new double[n * n];
    for (mplapackint i = 0; i < n * n; i++)
        a[i] = 0.0;
    for (mplapackint i = 1; i < n; i++)
        a[i + (i - 1) * lda] = 1.0;
    for (mplapackint j = 0; j < n; j++)
        a[j + (n - 1) * lda] = -coef[n - j] / coef[0];
    double *wr = new double[n];
    double *wi = new double[n];
    double *vl = new double[1];
    double *vr = new double[1];
    double wk;
    Rgeev("N", "N", n, a, lda, wr, wi, vl, ldv, vr, ldv, &wk, lwork, info);
    lwork = castINTEGER_double(wk);
    double *work = new double[lwork];
    Rgeev("N", "N", n, a, lda, wr, wi, vl, ldv, vr, ldv, work, lwork, info);
    sort_real(n, wr);
    double maxerr = 0.0;
    for (mplapackint i = 0; i < n; i++) {
        double err = abs(wr[i] - double(i + 1));
        if (maxerr < err)
            maxerr = err;
        printf("root[%ld] = ", (long)i); printnum(wr[i]); printf(", error = "); printnum(err); printf("\n");
    }
    printf("max root error = "); printnum(maxerr); printf("\n");
    delete[] work;
    delete[] vr;
    delete[] vl;
    delete[] wi;
    delete[] wr;
    delete[] a;
    delete[] coef;
    return info != 0 ? 1 : 0;
}
