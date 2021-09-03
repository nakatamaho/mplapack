//public domain
#include <mpblas_double.h>
#include <mplapack_double.h>
#include <iostream>
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
void inv_hilbert_matrix(int n) {
    mplapackint lwork, info;
    double *ainv = new double[n * n];
    double *aorg = new double[n * n];
    double *c = new double[n * n];
    mplapackint *ipiv = new mplapackint[n];
    double one = 1.0, zero = 0.0, mtmp;

    // setting A matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mtmp = (i + 1) + (j + 1) - 1;
            ainv[i + j * n] = one / mtmp;
            aorg[i + j * n] = one / mtmp;
        }
    }

    printf("a = "); printmat(n, n, ainv, n); printf("\n");
    // work space query
    lwork = -1;
    double *work = new double[1];
    Rgetri(n, ainv, n, ipiv, work, lwork, info);
    lwork = castINTEGER_double(work[0]);
    delete[] work;
    work = new double[std::max(1, (int)lwork)];

    // inverse matrix
    Rgetrf(n, n, ainv, n, ipiv, info);
    Rgetri(n, ainv, n, ipiv, work, lwork, info);
    printf("ainv = "); printmat(n, n, ainv, n); printf("\n");

    // Left residual |ainv * a -I|/|ainv||a|
    // is usually accurate than Right residual |a * ainv -I|/|ainv||a|  See chap.14 of
    // Accuracy and Stability of Numerical Algorithms by Nicholas J. Higham
    // https://doi.org/10.1137/1.9780898718027
    one = 1.0, zero = 0.0;
    Rgemm("N", "N", n, n, n, one, aorg, n, ainv, n, zero, c, n);
    printf("a * ainv ="); printmat(n, n, c, n); printf("\n");
    printf("InfnormR:(a * ainv - I)=");
    mtmp = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (mtmp < abs(c[i + j * n] - ((i == j) ? 1.0 : 0.0)))
                mtmp = abs(c[i + j * n] - ((i == j) ? 1.0 : 0.0));
        }
    }
    printnum(mtmp); printf("\n");

    Rgemm("N", "N", n, n, n, one, ainv, n, aorg, n, zero, c, n);
    printf("ainv * a ="); printmat(n, n, c, n); printf("\n");
    printf("InfnormL:(ainv * a - I)=");
    mtmp = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (mtmp < abs(c[i + j * n] - ((i == j) ? 1.0 : 0.0)))
                mtmp = abs(c[i + j * n] - ((i == j) ? 1.0 : 0.0));
        }
    }
    printnum(mtmp); printf("\n");

    delete[] ainv;
    delete[] aorg;
    delete[] c;
    delete[] ipiv;
    delete[] work;
}

int main()
{
    for (int n = 1; n < 50; n++) {
	printf("# inversion of Hilbert matrix of order n=%d\n", n);
	inv_hilbert_matrix(n);
    }
}
