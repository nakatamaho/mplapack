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
    mplapackint n = 5, nrhs = 1, lda = n, ldb = n, info;
    qd_real *a = new qd_real[n * n];
    qd_real *af = new qd_real[n * n];
    qd_real *b = new qd_real[n];
    qd_real *x = new qd_real[n];
    qd_real *xexact = new qd_real[n];
    qd_real *r = new qd_real[n];
    qd_real *c = new qd_real[n];
    qd_real *ferr = new qd_real[nrhs];
    qd_real *berr = new qd_real[nrhs];
    qd_real *work = new qd_real[4 * n];
    mplapackint *iwork = new mplapackint[n];
    mplapackint *ipiv = new mplapackint[n];
    char equed[2]; equed[0] = 'N'; equed[1] = '\0';
    qd_real rcond;
    for (mplapackint j = 0; j < n; j++) for (mplapackint i = 0; i < n; i++) a[i + j * lda] = qd_real(1.0) / qd_real(i + j + 1);
    for (mplapackint i = 0; i < n; i++) xexact[i] = (i % 2 == 0) ? qd_real(1.0) : qd_real(-1.0);
    for (mplapackint i = 0; i < n; i++) { b[i] = 0; for (mplapackint k = 0; k < n; k++) b[i] = b[i] + a[i + k * lda] * xexact[k]; }
    printf("A = "); printmat(n, n, a, lda); printf("\n");
    printf("b = "); printvec(b, n); printf("\n");
    Rgesvx("N", "N", n, nrhs, a, lda, af, lda, ipiv, equed, r, c, b, ldb, x, ldb, rcond, ferr, berr, work, iwork, info);
    printf("x = "); printvec(x, n); printf("\n");
    printf("rcond = "); printnum(rcond); printf("\n");
    printf("ferr = "); printvec(ferr, nrhs); printf("\n");
    printf("berr = "); printvec(berr, nrhs); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, x, ldb, xexact, n)); printf("\n");
    delete[] ipiv; delete[] iwork; delete[] work; delete[] berr; delete[] ferr; delete[] c; delete[] r; delete[] xexact; delete[] x; delete[] b; delete[] af; delete[] a;
    return info != 0 ? 1 : 0;
}
