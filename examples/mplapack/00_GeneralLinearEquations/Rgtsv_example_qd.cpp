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
qd_real maxabs(qd_real a, qd_real b) {
    qd_real d = abs(a - b);
    return d;
}

qd_real max_solution_error(mplapackint n, mplapackint nrhs, qd_real *x, mplapackint ldx, qd_real *xexact, mplapackint ldxexact) {
    qd_real err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            qd_real d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

qd_real max_residual(mplapackint m, mplapackint n, mplapackint nrhs, qd_real *a, mplapackint lda, qd_real *x, mplapackint ldx, qd_real *b, mplapackint ldb) {
    qd_real err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            qd_real s = 0.0;
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            qd_real d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

int main() {
    mplapackint n = 6, nrhs = 1, ldb = n, info;
    qd_real *dl = new qd_real[n - 1];
    qd_real *d = new qd_real[n];
    qd_real *du = new qd_real[n - 1];
    qd_real *a = new qd_real[n * n];
    qd_real *b = new qd_real[n];
    qd_real *borg = new qd_real[n];
    qd_real *xexact = new qd_real[n];
    for (mplapackint i = 0; i < n * n; i++) a[i] = 0;
    for (mplapackint i = 0; i < n; i++) { d[i] = 2; xexact[i] = i + 1; a[i + i * n] = 2; }
    for (mplapackint i = 0; i < n - 1; i++) { dl[i] = -1; du[i] = -1; a[i + 1 + i * n] = -1; a[i + (i + 1) * n] = -1; }
    for (mplapackint i = 0; i < n; i++) {
        b[i] = 0;
        for (mplapackint k = 0; k < n; k++) b[i] = b[i] + a[i + k * n] * xexact[k];
        borg[i] = b[i];
    }
    printf("A = "); printmat(n, n, a, n); printf("\n");
    Rgtsv(n, nrhs, dl, d, du, b, ldb, info);
    printf("x = "); printvec(b, n); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, b, ldb, xexact, n)); printf("\n");
    printf("max residual = "); printnum(max_residual(n, n, nrhs, a, n, b, ldb, borg, ldb)); printf("\n");
    delete[] xexact; delete[] borg; delete[] b; delete[] a; delete[] du; delete[] d; delete[] dl;
    return info != 0 ? 1 : 0;
}
