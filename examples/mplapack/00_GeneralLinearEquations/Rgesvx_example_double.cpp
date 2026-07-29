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
double maxabs(double a, double b) {
    double d = abs(a - b);
    return d;
}

double max_solution_error(mplapackint n, mplapackint nrhs, double *x, mplapackint ldx, double *xexact, mplapackint ldxexact) {
    double err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            double d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

double max_residual(mplapackint m, mplapackint n, mplapackint nrhs, double *a, mplapackint lda, double *x, mplapackint ldx, double *b, mplapackint ldb) {
    double err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            double s = 0.0;
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            double d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

int main() {
    mplapackint n = 5, nrhs = 1, lda = n, ldb = n, info;
    double *a = new double[n * n];
    double *af = new double[n * n];
    double *b = new double[n];
    double *x = new double[n];
    double *xexact = new double[n];
    double *r = new double[n];
    double *c = new double[n];
    double *ferr = new double[nrhs];
    double *berr = new double[nrhs];
    double *work = new double[4 * n];
    mplapackint *iwork = new mplapackint[n];
    mplapackint *ipiv = new mplapackint[n];
    char equed[2];
    equed[0] = 'N';
    equed[1] = '\0';
    double rcond;
    for (mplapackint j = 0; j < n; j++) for (mplapackint i = 0; i < n; i++) a[i + j * lda] = double(1.0) / double(i + j + 1);
    for (mplapackint i = 0; i < n; i++) xexact[i] = (i % 2 == 0) ? double(1.0) : double(-1.0);
    for (mplapackint i = 0; i < n; i++) {
        b[i] = 0;
        for (mplapackint k = 0; k < n; k++)
            b[i] = b[i] + a[i + k * lda] * xexact[k];
    }
    printf("A = "); printmat(n, n, a, lda); printf("\n");
    printf("b = "); printvec(b, n); printf("\n");
    Rgesvx("N", "N", n, nrhs, a, lda, af, lda, ipiv, equed, r, c, b, ldb, x, ldb, rcond, ferr, berr, work, iwork, info);
    printf("x = "); printvec(x, n); printf("\n");
    printf("rcond = "); printnum(rcond); printf("\n");
    printf("ferr = "); printvec(ferr, nrhs); printf("\n");
    printf("berr = "); printvec(berr, nrhs); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, x, ldb, xexact, n)); printf("\n");
    delete[] ipiv;
    delete[] iwork;
    delete[] work;
    delete[] berr;
    delete[] ferr;
    delete[] c;
    delete[] r;
    delete[] xexact;
    delete[] x;
    delete[] b;
    delete[] af;
    delete[] a;
    return info != 0 ? 1 : 0;
}
