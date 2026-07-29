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
qd_real max_solution_error(mplapackint n, mplapackint nrhs, qd_complex *x, mplapackint ldx, qd_complex *xexact, mplapackint ldxexact) {
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

qd_real max_residual(mplapackint m, mplapackint n, mplapackint nrhs, qd_complex *a, mplapackint lda, qd_complex *x, mplapackint ldx, qd_complex *b, mplapackint ldb) {
    qd_real err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            qd_complex s = qd_complex(0.0, 0.0);
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            qd_real d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

void set_problem(mplapackint m, mplapackint n, qd_complex *a, mplapackint lda, qd_complex *b, mplapackint ldb, qd_complex *xexact) {
    for (mplapackint i = 0; i < m; i++) {
        a[i + 0 * lda] = qd_complex(1.0, 0.0);
        a[i + 1 * lda] = qd_complex(i, 1.0);
    }
    xexact[0] = qd_complex(1.0, -1.0);
    xexact[1] = qd_complex(2.0, 1.0);
    for (mplapackint i = 0; i < m; i++)
        b[i] = a[i + 0 * lda] * xexact[0] + a[i + 1 * lda] * xexact[1];
}

int main() {
    mplapackint m = 4, n = 2, nrhs = 1, lda = m, ldb = m, info, lwork = -1;
    qd_complex *a = new qd_complex[lda * n];
    qd_complex *aorg = new qd_complex[lda * n];
    qd_complex *b = new qd_complex[ldb];
    qd_complex *borg = new qd_complex[ldb];
    qd_complex *xexact = new qd_complex[n];
    set_problem(m, n, a, lda, b, ldb, xexact);
    for (mplapackint i = 0; i < lda * n; i++)
        aorg[i] = a[i];
    for (mplapackint i = 0; i < ldb; i++)
        borg[i] = b[i];
    qd_complex wk;
    Cgels("N", m, n, nrhs, a, lda, b, ldb, &wk, lwork, info);
    lwork = castINTEGER_qd(wk.real());
    qd_complex *work = new qd_complex[lwork];
    Cgels("N", m, n, nrhs, a, lda, b, ldb, work, lwork, info);
    printf("A = "); printmat(m, n, aorg, lda); printf("\n");
    printf("x = "); printvec(b, n); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, b, ldb, xexact, n)); printf("\n");
    printf("max residual = "); printnum(max_residual(m, n, nrhs, aorg, lda, b, ldb, borg, ldb)); printf("\n");
    delete[] work;
    delete[] xexact;
    delete[] borg;
    delete[] b;
    delete[] aorg;
    delete[] a;
    return info != 0 ? 1 : 0;
}
