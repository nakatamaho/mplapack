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
inline void printnum(mpfc_class ctmp) { gmp_printf(GMP_FORMAT GMP_FORMAT "i", ctmp.real().get_mpf_t(), ctmp.imag().get_mpf_t()); }

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
mpf_class max_solution_error(mplapackint n, mplapackint nrhs, mpfc_class *x, mplapackint ldx, mpfc_class *xexact, mplapackint ldxexact) {
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

mpf_class max_residual(mplapackint m, mplapackint n, mplapackint nrhs, mpfc_class *a, mplapackint lda, mpfc_class *x, mplapackint ldx, mpfc_class *b, mplapackint ldb) {
    mpf_class err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            mpfc_class s = mpfc_class(0.0, 0.0);
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            mpf_class d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

void set_problem(mplapackint m, mplapackint n, mpfc_class *a, mplapackint lda, mpfc_class *b, mplapackint ldb, mpfc_class *xexact) {
    for (mplapackint i = 0; i < m; i++) {
        a[i + 0 * lda] = mpfc_class(1.0, 0.0);
        a[i + 1 * lda] = mpfc_class(i, 1.0);
    }
    xexact[0] = mpfc_class(1.0, -1.0);
    xexact[1] = mpfc_class(2.0, 1.0);
    for (mplapackint i = 0; i < m; i++)
        b[i] = a[i + 0 * lda] * xexact[0] + a[i + 1 * lda] * xexact[1];
}

int main() {
    mplapackint m = 4, n = 2, nrhs = 1, lda = m, ldb = m, info, lwork = -1, rank;
    mpfc_class *a = new mpfc_class[lda * n];
    mpfc_class *aorg = new mpfc_class[lda * n];
    mpfc_class *b = new mpfc_class[ldb];
    mpfc_class *borg = new mpfc_class[ldb];
    mpfc_class *xexact = new mpfc_class[n];
    mpf_class *s = new mpf_class[n];
    mpf_class *rwork = new mpf_class[5 * n];
    set_problem(m, n, a, lda, b, ldb, xexact);
    for (mplapackint i = 0; i < lda * n; i++)
        aorg[i] = a[i];
    for (mplapackint i = 0; i < ldb; i++)
        borg[i] = b[i];
    mpfc_class wk;
    Cgelss(m, n, nrhs, a, lda, b, ldb, s, mpf_class(-1.0), rank, &wk, lwork, rwork, info);
    lwork = castINTEGER_gmp(wk.real());
    mpfc_class *work = new mpfc_class[lwork];
    Cgelss(m, n, nrhs, a, lda, b, ldb, s, mpf_class(-1.0), rank, work, lwork, rwork, info);
    printf("singular values = "); printvec(s, n); printf("\n");
    printf("rank = %ld\n", (long)rank);
    printf("x = "); printvec(b, n); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, b, ldb, xexact, n)); printf("\n");
    printf("max residual = "); printnum(max_residual(m, n, nrhs, aorg, lda, b, ldb, borg, ldb)); printf("\n");
    delete[] work;
    delete[] rwork;
    delete[] s;
    delete[] xexact;
    delete[] borg;
    delete[] b;
    delete[] aorg;
    delete[] a;
    return info != 0 ? 1 : 0;
}
