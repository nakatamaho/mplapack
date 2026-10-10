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

td_real one_norm(mplapackint n, td_real *a, mplapackint lda) {
    td_real anorm = 0.0;
    for (mplapackint j = 0; j < n; j++) {
        td_real s = 0.0;
        for (mplapackint i = 0; i < n; i++)
            s = s + abs(a[i + j * lda]);
        if (anorm < s)
            anorm = s;
    }
    return anorm;
}
td_real inf_norm(mplapackint n, td_real *a, mplapackint lda) {
    td_real anorm = 0.0;
    for (mplapackint i = 0; i < n; i++) {
        td_real s = 0.0;
        for (mplapackint j = 0; j < n; j++)
            s = s + abs(a[i + j * lda]);
        if (anorm < s)
            anorm = s;
    }
    return anorm;
}
int main() {
    mplapackint n = 3, lda = n, info;
    td_real *a = new td_real[n * n];
    td_real *lu = new td_real[n * n];
    td_real *work = new td_real[4 * n];
    mplapackint *iwork = new mplapackint[n];
    mplapackint *ipiv = new mplapackint[n];
    for (mplapackint i = 0; i < n * n; i++) a[i] = 0.0;
    a[0 + 0 * lda] = 1.0; a[1 + 1 * lda] = td_real(1.0e-3); a[2 + 2 * lda] = td_real(1.0e-6);
    for (mplapackint i = 0; i < n * n; i++) lu[i] = a[i];
    Rgetrf(n, n, lu, lda, ipiv, info);
    td_real rcond1 = 0.0, rcondi = 0.0;
    if (info == 0) Rgecon("1", n, lu, lda, one_norm(n, a, lda), rcond1, work, iwork, info);
    if (info == 0) Rgecon("I", n, lu, lda, inf_norm(n, a, lda), rcondi, work, iwork, info);
    printf("A = "); printmat(n, n, a, lda); printf("\n");
    printf("true cond_1 = "); printnum(td_real(1.0e6)); printf("\n");
    printf("rcond_1 = "); printnum(rcond1); printf("\n");
    printf("rcond_inf = "); printnum(rcondi); printf("\n");
    delete[] ipiv;
    delete[] iwork;
    delete[] work;
    delete[] lu;
    delete[] a;
    return info != 0 ? 1 : 0;
}
