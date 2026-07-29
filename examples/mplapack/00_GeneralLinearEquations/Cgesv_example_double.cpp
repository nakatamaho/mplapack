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
inline void printnum(std::complex<double> ctmp) { printf(DOUBLE_FORMAT DOUBLE_FORMAT "i", ctmp.real(), ctmp.imag()); }

//Matlab/Octave format
template <class X> void printvec(X *a, int len) {
    X tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

template <class X> void printmat(int n, int m, X *a, int lda)
{
    X mtmp;

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
double max_solution_error(mplapackint n, mplapackint nrhs, std::complex<double> *x, mplapackint ldx, std::complex<double> *xexact, mplapackint ldxexact) {
    double err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            double d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

double max_residual(mplapackint m, mplapackint n, mplapackint nrhs, std::complex<double> *a, mplapackint lda, std::complex<double> *x, mplapackint ldx, std::complex<double> *b, mplapackint ldb) {
    double err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            std::complex<double> s = std::complex<double>(0.0, 0.0);
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            double d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

int main() {
    mplapackint n = 2, nrhs = 1, lda = n, ldb = n, info;
    std::complex<double> *a = new std::complex<double>[n * n];
    std::complex<double> *aorg = new std::complex<double>[n * n];
    std::complex<double> *b = new std::complex<double>[nrhs * ldb];
    std::complex<double> *borg = new std::complex<double>[nrhs * ldb];
    std::complex<double> *xexact = new std::complex<double>[nrhs * n];
    mplapackint *ipiv = new mplapackint[n];

    a[0 + 0 * lda] = std::complex<double>(3.0, 0.0);  a[0 + 1 * lda] = std::complex<double>(1.0, -1.0);
    a[1 + 0 * lda] = std::complex<double>(1.0, 1.0);  a[1 + 1 * lda] = std::complex<double>(4.0, 0.0);
    xexact[0] = std::complex<double>(1.0, 2.0);
    xexact[1] = std::complex<double>(-1.0, 1.0);
    for (mplapackint i = 0; i < n * n; i++) aorg[i] = a[i];
    for (mplapackint i = 0; i < n; i++) {
        b[i] = std::complex<double>(0.0, 0.0);
        for (mplapackint k = 0; k < n; k++) b[i] = b[i] + aorg[i + k * lda] * xexact[k];
        borg[i] = b[i];
    }

    printf("A = "); printmat(n, n, aorg, lda); printf("\n");
    printf("b = "); printvec(b, n); printf("\n");
    Cgesv(n, (mplapackint)1, a, lda, ipiv, b, ldb, info);
    printf("x = "); printvec(b, n); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, b, ldb, xexact, n)); printf("\n");
    printf("max residual = "); printnum(max_residual(n, n, nrhs, aorg, lda, b, ldb, borg, ldb)); printf("\n");

    delete[] ipiv; delete[] xexact; delete[] borg; delete[] b; delete[] aorg; delete[] a;
    return info != 0 ? 1 : 0;
}
