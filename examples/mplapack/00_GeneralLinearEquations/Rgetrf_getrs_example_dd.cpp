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
dd_real maxabs(dd_real a, dd_real b) {
    dd_real d = abs(a - b);
    return d;
}

dd_real max_solution_error(mplapackint n, mplapackint nrhs, dd_real *x, mplapackint ldx, dd_real *xexact, mplapackint ldxexact) {
    dd_real err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            dd_real d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

dd_real max_residual(mplapackint m, mplapackint n, mplapackint nrhs, dd_real *a, mplapackint lda, dd_real *x, mplapackint ldx, dd_real *b, mplapackint ldb) {
    dd_real err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            dd_real s = 0.0;
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            dd_real d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

int main() {
    mplapackint n = 3, nrhs = 2, lda = n, ldb = n, info;
    dd_real *a = new dd_real[n * n];
    dd_real *aorg = new dd_real[n * n];
    dd_real *b = new dd_real[nrhs * ldb];
    dd_real *borg = new dd_real[nrhs * ldb];
    dd_real *xexact = new dd_real[nrhs * n];
    mplapackint *ipiv = new mplapackint[n];

    a[0 + 0 * lda] = 4; a[0 + 1 * lda] = 1; a[0 + 2 * lda] = 2;
    a[1 + 0 * lda] = 0; a[1 + 1 * lda] = 3; a[1 + 2 * lda] = -1;
    a[2 + 0 * lda] = 2; a[2 + 1 * lda] = 1; a[2 + 2 * lda] = 5;
    xexact[0] = 1; xexact[1] = -1; xexact[2] = 2;
    xexact[0 + 1 * n] = 2; xexact[1 + 1 * n] = 0; xexact[2 + 1 * n] = -1;
    for (mplapackint i = 0; i < n * n; i++) aorg[i] = a[i];
    for (mplapackint j = 0; j < nrhs; j++) for (mplapackint i = 0; i < n; i++) {
        b[i + j * ldb] = 0;
        for (mplapackint k = 0; k < n; k++) b[i + j * ldb] = b[i + j * ldb] + aorg[i + k * lda] * xexact[k + j * n];
        borg[i + j * ldb] = b[i + j * ldb];
    }

    printf("A = "); printmat(n, n, aorg, lda); printf("\n");
    printf("B = "); printmat(n, nrhs, b, ldb); printf("\n");
    Rgetrf(n, n, a, lda, ipiv, info);
    if (info == 0) Rgetrs("N", n, nrhs, a, lda, ipiv, b, ldb, info);
    printf("X = "); printmat(n, nrhs, b, ldb); printf("\n");
    printf("max |X-X_exact| = "); printnum(max_solution_error(n, nrhs, b, ldb, xexact, n)); printf("\n");
    printf("max residual = "); printnum(max_residual(n, n, nrhs, aorg, lda, b, ldb, borg, ldb)); printf("\n");

    delete[] ipiv;
    delete[] xexact;
    delete[] borg;
    delete[] b;
    delete[] aorg;
    delete[] a;
    return info != 0 ? 1 : 0;
}
