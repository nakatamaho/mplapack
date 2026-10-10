//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_td.h>
#include <mplapack_td.h>

#define TD_PRECISION_SHORT 16

inline void printnum(td_real rtmp) {
    std::cout.precision(TD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

//Matlab/Octave format
void printvec(td_real *a, int len) {
    td_real tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, td_real * a, int lda)
{
    td_real mtmp;
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
td_real maxabs(td_real a, td_real b) {
    td_real d = abs(a - b);
    return d;
}

td_real max_solution_error(mplapackint n, mplapackint nrhs, td_real *x, mplapackint ldx, td_real *xexact, mplapackint ldxexact) {
    td_real err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            td_real d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

td_real max_residual(mplapackint m, mplapackint n, mplapackint nrhs, td_real *a, mplapackint lda, td_real *x, mplapackint ldx, td_real *b, mplapackint ldb) {
    td_real err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            td_real s = 0.0;
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            td_real d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

int main() {
    mplapackint n = 5, nrhs = 1, lda = n, ldb = n, info;
    td_real *a = new td_real[n * n];
    td_real *af = new td_real[n * n];
    td_real *b = new td_real[n];
    td_real *x = new td_real[n];
    td_real *xexact = new td_real[n];
    td_real *r = new td_real[n];
    td_real *c = new td_real[n];
    td_real *ferr = new td_real[nrhs];
    td_real *berr = new td_real[nrhs];
    td_real *work = new td_real[4 * n];
    mplapackint *iwork = new mplapackint[n];
    mplapackint *ipiv = new mplapackint[n];
    char equed[2];
    equed[0] = 'N';
    equed[1] = '\0';
    td_real rcond;
    for (mplapackint j = 0; j < n; j++) for (mplapackint i = 0; i < n; i++) a[i + j * lda] = td_real(1.0) / td_real((double)(i + j + 1));
    for (mplapackint i = 0; i < n; i++) xexact[i] = (i % 2 == 0) ? td_real(1.0) : td_real(-1.0);
    for (mplapackint i = 0; i < n; i++) {
        b[i] = 0.0;
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
