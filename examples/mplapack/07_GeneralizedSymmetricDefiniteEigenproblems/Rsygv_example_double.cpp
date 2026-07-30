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
double max_eigen_residual(mplapackint n, double *a, double *b, double lambda, double *z) {
    double err = 0.0;
    for (mplapackint i = 0; i < n; i++) {
        double s = 0.0, t = 0.0;
        for (mplapackint j = 0; j < n; j++) {
            s = s + a[i + j * n] * z[j];
            t = t + b[i + j * n] * z[j];
        }
        double d = abs(s - lambda * t);
        if (err < d)
            err = d;
    }
    return err;
}

int main() {
    mplapackint n = 2, lda = n, ldb = n, info, lwork = -1;
    double *a = new double[n * n];
    double *b = new double[n * n];
    double *aorg = new double[n * n];
    double *borg = new double[n * n];
    double *w = new double[n];
    for (mplapackint i = 0; i < n * n; i++) {
        a[i] = 0.0;
        b[i] = 0.0;
    }
    a[0] = 2.0;
    a[3] = 6.0;
    b[0] = 1.0;
    b[3] = 2.0;
    for (mplapackint i = 0; i < n * n; i++) {
        aorg[i] = a[i];
        borg[i] = b[i];
    }
    double wk;
    Rsygv((mplapackint)1, "V", "U", n, a, lda, b, ldb, w, &wk, lwork, info);
    lwork = castINTEGER_double(wk);
    double *work = new double[lwork];
    Rsygv((mplapackint)1, "V", "U", n, a, lda, b, ldb, w, work, lwork, info);
    printf("eigenvalues = "); printvec(w, n); printf("\n");
    printf("eigenvectors = "); printmat(n, n, a, lda); printf("\n");
    for (mplapackint j = 0; j < n; j++) {
        printf("residual[%ld] = ", (long)j); printnum(max_eigen_residual(n, aorg, borg, w[j], &a[j * lda])); printf("\n");
    }
    delete[] work;
    delete[] w;
    delete[] borg;
    delete[] aorg;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
