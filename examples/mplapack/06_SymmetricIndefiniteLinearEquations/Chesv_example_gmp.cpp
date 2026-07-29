//public domain
#include <mpblas_gmp.h>
#include <mplapack_gmp.h>
#include <iostream>
#include <cstring>
#include <algorithm>

#define GMP_FORMAT "%+68.64Fe"
#define GMP_SHORT_FORMAT "%+20.16Fe"

inline void printnum(mpf_class rtmp) { gmp_printf(GMP_FORMAT, rtmp.get_mpf_t()); }
inline void printnum_short(mpf_class rtmp) { gmp_printf(GMP_SHORT_FORMAT, rtmp.get_mpf_t()); }
inline void printnum(mpc_class ctmp) { gmp_printf(GMP_FORMAT GMP_FORMAT "i", ctmp.real().get_mpf_t(), ctmp.imag().get_mpf_t()); }

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
mpf_class max_solution_error(mplapackint n, mplapackint nrhs, mpc_class *x, mplapackint ldx, mpc_class *xexact, mplapackint ldxexact) {
    mpf_class err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            mpf_class d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

mpf_class max_residual(mplapackint m, mplapackint n, mplapackint nrhs, mpc_class *a, mplapackint lda, mpc_class *x, mplapackint ldx, mpc_class *b, mplapackint ldb) {
    mpf_class err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            mpc_class s = mpc_class(0.0, 0.0);
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            mpf_class d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

int main() {
    mplapackint n = 2, nrhs = 1, lda = n, ldb = n, info, lwork = -1;
    mpc_class *a = new mpc_class[n * n];
    mpc_class *aorg = new mpc_class[n * n];
    mpc_class *b = new mpc_class[n];
    mpc_class *borg = new mpc_class[n];
    mpc_class *xexact = new mpc_class[n];
    mplapackint *ipiv = new mplapackint[n];
    a[0] = mpc_class(2.0, 0.0);
    a[1] = mpc_class(1.0, -1.0);
    a[2] = mpc_class(1.0, 1.0);
    a[3] = mpc_class(-3.0, 0.0);
    xexact[0] = mpc_class(1.0, -1.0);
    xexact[1] = mpc_class(2.0, 1.0);
    for (mplapackint i = 0; i < n * n; i++)
        aorg[i] = a[i];
    for (mplapackint i = 0; i < n; i++) {
        b[i] = mpc_class(0.0, 0.0);
        for (mplapackint k = 0; k < n; k++)
            b[i] = b[i] + aorg[i + k * lda] * xexact[k];
        borg[i] = b[i];
    }
    mpc_class wk;
    Chesv("U", n, nrhs, a, lda, ipiv, b, ldb, &wk, lwork, info);
    lwork = castINTEGER_gmp(wk.real());
    mpc_class *work = new mpc_class[lwork];
    Chesv("U", n, nrhs, a, lda, ipiv, b, ldb, work, lwork, info);
    printf("A = "); printmat(n, n, aorg, lda); printf("\n");
    printf("x = "); printvec(b, n); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, b, ldb, xexact, n)); printf("\n");
    printf("max residual = "); printnum(max_residual(n, n, nrhs, aorg, lda, b, ldb, borg, ldb)); printf("\n");
    delete[] work;
    delete[] ipiv;
    delete[] xexact;
    delete[] borg;
    delete[] b;
    delete[] aorg;
    delete[] a;
    return info != 0 ? 1 : 0;
}
