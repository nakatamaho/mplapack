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
    mplapackint n = 3, m = 2, p = 3, lda = n, ldb = n, info, lwork = -1;
    qd_real *a = new qd_real[lda * m];
    qd_real *b = new qd_real[ldb * p];
    qd_real *d = new qd_real[n];
    qd_real *x = new qd_real[m];
    qd_real *y = new qd_real[p];
    qd_real *xexact = new qd_real[m];
    qd_real *yexact = new qd_real[p];
    for (mplapackint i = 0; i < lda * m; i++)
        a[i] = 0;
    a[0] = 1;
    a[1] = 0;
    a[2] = 1;
    a[0 + lda] = 0;
    a[1 + lda] = 1;
    a[2 + lda] = 1;
    for (mplapackint i = 0; i < ldb * p; i++)
        b[i] = 0;
    for (mplapackint i = 0; i < n; i++)
        b[i + i * ldb] = 1;
    xexact[0] = 1;
    xexact[1] = 2;
    yexact[0] = qd_real(0.5);
    yexact[1] = qd_real(-0.5);
    yexact[2] = 1;
    for (mplapackint i = 0; i < n; i++)
        d[i] = a[i] * xexact[0] + a[i + lda] * xexact[1] + yexact[i];
    qd_real wk;
    Rggglm(n, m, p, a, lda, b, ldb, d, x, y, &wk, lwork, info);
    lwork = castINTEGER_qd(wk);
    qd_real *work = new qd_real[lwork];
    Rggglm(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info);
    printf("x = "); printvec(x, m); printf("\n");
    printf("y = "); printvec(y, p); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(m, (mplapackint)1, x, m, xexact, m)); printf("\n");
    printf("max |y-y_exact| = "); printnum(max_solution_error(p, (mplapackint)1, y, p, yexact, p)); printf("\n");
    delete[] work;
    delete[] yexact;
    delete[] xexact;
    delete[] y;
    delete[] x;
    delete[] d;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
