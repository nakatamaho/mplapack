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

qd_real one_norm(mplapackint n, qd_real *a, mplapackint lda) {
    qd_real anorm = 0.0;
    for (mplapackint j = 0; j < n; j++) {
        qd_real s = 0.0;
        for (mplapackint i = 0; i < n; i++)
            s = s + abs(a[i + j * lda]);
        if (anorm < s)
            anorm = s;
    }
    return anorm;
}
qd_real inf_norm(mplapackint n, qd_real *a, mplapackint lda) {
    qd_real anorm = 0.0;
    for (mplapackint i = 0; i < n; i++) {
        qd_real s = 0.0;
        for (mplapackint j = 0; j < n; j++)
            s = s + abs(a[i + j * lda]);
        if (anorm < s)
            anorm = s;
    }
    return anorm;
}
int main() {
    mplapackint n = 3, lda = n, info;
    qd_real *a = new qd_real[n * n];
    qd_real *lu = new qd_real[n * n];
    qd_real *work = new qd_real[4 * n];
    mplapackint *iwork = new mplapackint[n];
    mplapackint *ipiv = new mplapackint[n];
    for (mplapackint i = 0; i < n * n; i++) a[i] = 0.0;
    a[0 + 0 * lda] = 1.0; a[1 + 1 * lda] = qd_real(1.0e-3); a[2 + 2 * lda] = qd_real(1.0e-6);
    for (mplapackint i = 0; i < n * n; i++) lu[i] = a[i];
    Rgetrf(n, n, lu, lda, ipiv, info);
    qd_real rcond1 = 0.0, rcondi = 0.0;
    if (info == 0) Rgecon("1", n, lu, lda, one_norm(n, a, lda), rcond1, work, iwork, info);
    if (info == 0) Rgecon("I", n, lu, lda, inf_norm(n, a, lda), rcondi, work, iwork, info);
    printf("A = "); printmat(n, n, a, lda); printf("\n");
    printf("true cond_1 = "); printnum(qd_real(1.0e6)); printf("\n");
    printf("rcond_1 = "); printnum(rcond1); printf("\n");
    printf("rcond_inf = "); printnum(rcondi); printf("\n");
    delete[] ipiv;
    delete[] iwork;
    delete[] work;
    delete[] lu;
    delete[] a;
    return info != 0 ? 1 : 0;
}
