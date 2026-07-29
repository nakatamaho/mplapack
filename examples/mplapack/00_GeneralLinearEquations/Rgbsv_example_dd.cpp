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
    mplapackint n = 6, kl = 2, ku = 1, nrhs = 1, ldab = 2 * kl + ku + 1, ldb = n, info;
    dd_real *a = new dd_real[n * n];
    dd_real *ab = new dd_real[ldab * n];
    dd_real *b = new dd_real[n];
    dd_real *borg = new dd_real[n];
    dd_real *xexact = new dd_real[n];
    mplapackint *ipiv = new mplapackint[n];
    for (mplapackint i = 0; i < n * n; i++) a[i] = 0;
    for (mplapackint i = 0; i < ldab * n; i++) ab[i] = 0;
    for (mplapackint i = 0; i < n; i++) {
        xexact[i] = i + 1;
        for (mplapackint j = 0; j < n; j++) {
            if (i == j) a[i + j * n] = 6;
            else if (i > j && i - j <= kl) a[i + j * n] = -1;
            else if (j > i && j - i <= ku) a[i + j * n] = 2;
        }
    }
    for (mplapackint j = 0; j < n; j++) for (mplapackint i = 0; i < n; i++) if (a[i + j * n] != 0) {
        mplapackint row = kl + ku + i - j;
        ab[row + j * ldab] = a[i + j * n];
    }
    for (mplapackint i = 0; i < n; i++) {
        b[i] = 0;
        for (mplapackint k = 0; k < n; k++) b[i] = b[i] + a[i + k * n] * xexact[k];
        borg[i] = b[i];
    }
    printf("A = "); printmat(n, n, a, n); printf("\n");
    printf("AB = "); printmat(ldab, n, ab, ldab); printf("\n");
    Rgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
    printf("x = "); printvec(b, n); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, b, ldb, xexact, n)); printf("\n");
    printf("max residual = "); printnum(max_residual(n, n, nrhs, a, n, b, ldb, borg, ldb)); printf("\n");
    delete[] ipiv;
    delete[] xexact;
    delete[] borg;
    delete[] b;
    delete[] ab;
    delete[] a;
    return info != 0 ? 1 : 0;
}
