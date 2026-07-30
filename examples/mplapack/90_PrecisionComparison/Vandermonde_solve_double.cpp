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
double max_solution_error(mplapackint n, double *x, double *xexact) {
    double err = 0.0;
    for (mplapackint i = 0; i < n; i++) {
        double d = abs(x[i] - xexact[i]);
        if (err < d)
            err = d;
    }
    return err;
}
int main() {
    mplapackint n = 15, lda = n, ldb = n, info;
    double *a = new double[n * n];
    double *b = new double[n];
    double *xexact = new double[n];
    mplapackint *ipiv = new mplapackint[n];
    for (mplapackint j = 0; j < n; j++)
        xexact[j] = double(j % 3 - 1);
    for (mplapackint i = 0; i < n; i++) {
        double node = double(i + 1);
        double p = 1.0;
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
