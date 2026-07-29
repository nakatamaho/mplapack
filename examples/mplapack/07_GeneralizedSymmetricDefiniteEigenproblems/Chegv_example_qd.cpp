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

inline void printnum(qd_complex rtmp) {
    std::cout.precision(QD_PRECISION_SHORT);
    if (rtmp.real() >= 0.0) {
        std::cout << "+" << rtmp.real();
    } else {
        std::cout << rtmp.real();
    }
    if (rtmp.imag() >= 0.0) {
        std::cout << "+" << rtmp.imag() << "i";
    } else {
        std::cout << rtmp.imag() << "i";
    }
    return;
}

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
qd_real max_eigen_residual(mplapackint n, qd_complex *a, qd_complex *b, qd_real lambda, qd_complex *z) {
    qd_real err = 0.0;
    for (mplapackint i = 0; i < n; i++) {
        qd_complex s = qd_complex(0.0, 0.0), t = qd_complex(0.0, 0.0);
        for (mplapackint j = 0; j < n; j++) {
            s = s + a[i + j * n] * z[j];
            t = t + b[i + j * n] * z[j];
        }
        qd_real d = abs(s - lambda * t);
        if (err < d)
            err = d;
    }
    return err;
}

int main() {
    mplapackint n = 2, lda = n, ldb = n, info, lwork = -1;
    qd_complex *a = new qd_complex[n * n];
    qd_complex *b = new qd_complex[n * n];
    qd_complex *aorg = new qd_complex[n * n];
    qd_complex *borg = new qd_complex[n * n];
    qd_real *w = new qd_real[n];
    qd_real *rwork = new qd_real[3 * n];
    for (mplapackint i = 0; i < n * n; i++) {
        a[i] = qd_complex(0.0, 0.0);
        b[i] = qd_complex(0.0, 0.0);
    }
    a[0] = qd_complex(2.0, 0.0);
    a[3] = qd_complex(6.0, 0.0);
    b[0] = qd_complex(1.0, 0.0);
    b[3] = qd_complex(2.0, 0.0);
    for (mplapackint i = 0; i < n * n; i++) {
        aorg[i] = a[i];
        borg[i] = b[i];
    }
    qd_complex wk;
    Chegv((mplapackint)1, "V", "U", n, a, lda, b, ldb, w, &wk, lwork, rwork, info);
    lwork = castINTEGER_qd(wk.real());
    qd_complex *work = new qd_complex[lwork];
    Chegv((mplapackint)1, "V", "U", n, a, lda, b, ldb, w, work, lwork, rwork, info);
    printf("eigenvalues = "); printvec(w, n); printf("\n");
    printf("eigenvectors = "); printmat(n, n, a, lda); printf("\n");
    for (mplapackint j = 0; j < n; j++) {
        printf("residual[%ld] = ", (long)j); printnum(max_eigen_residual(n, aorg, borg, w[j], &a[j * lda])); printf("\n");
    }
    delete[] work;
    delete[] rwork;
    delete[] w;
    delete[] borg;
    delete[] aorg;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
