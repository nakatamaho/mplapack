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
dd_real max_eigen_residual(mplapackint n, dd_real *a, dd_real *b, dd_real lambda, dd_real *z) {
    dd_real err = 0.0;
    for (mplapackint i = 0; i < n; i++) {
        dd_real s = 0.0, t = 0.0;
        for (mplapackint j = 0; j < n; j++) {
            s = s + a[i + j * n] * z[j];
            t = t + b[i + j * n] * z[j];
        }
        dd_real d = abs(s - lambda * t);
        if (err < d)
            err = d;
    }
    return err;
}

int main() {
    mplapackint n = 2, lda = n, ldb = n, info, lwork = -1;
    dd_real *a = new dd_real[n * n];
    dd_real *b = new dd_real[n * n];
    dd_real *aorg = new dd_real[n * n];
    dd_real *borg = new dd_real[n * n];
    dd_real *w = new dd_real[n];
    for (mplapackint i = 0; i < n * n; i++) {
        a[i] = 0;
        b[i] = 0;
    }
    a[0] = 2;
    a[3] = 6;
    b[0] = 1;
    b[3] = 2;
    for (mplapackint i = 0; i < n * n; i++) {
        aorg[i] = a[i];
        borg[i] = b[i];
    }
    dd_real wk;
    Rsygv((mplapackint)1, "V", "U", n, a, lda, b, ldb, w, &wk, lwork, info);
    lwork = castINTEGER_dd(wk);
    dd_real *work = new dd_real[lwork];
    Rsygv((mplapackint)1, "V", "U", n, a, lda, b, ldb, w, work, lwork, info);
    printf("eigenvalues = "); printvec(w, n); printf("\n");
    printf("eigenvectors = "); printmat(n, n, a, lda); printf("\n");
    for (mplapackint j = 0; j < n; j++) {
        printf("residual[%ld] = ", (long)j); printnum(max_eigen_residual(n, aorg, borg, w[j], &a[j * lda])); printf("\n");
    }
    delete[] work;
    delete[] w;
    delete[] borg;
    delete[] aorg;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
