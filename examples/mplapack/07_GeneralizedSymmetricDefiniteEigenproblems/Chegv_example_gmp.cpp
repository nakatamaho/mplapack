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
mpf_class max_eigen_residual(mplapackint n, mpfc_class *a, mpfc_class *b, mpf_class lambda, mpfc_class *z) {
    mpf_class err = 0.0;
    for (mplapackint i = 0; i < n; i++) {
        mpfc_class s = mpfc_class(0.0, 0.0), t = mpfc_class(0.0, 0.0);
        for (mplapackint j = 0; j < n; j++) {
            s = s + a[i + j * n] * z[j];
            t = t + b[i + j * n] * z[j];
        }
        mpf_class d = abs(s - lambda * t);
        if (err < d)
            err = d;
    }
    return err;
}

int main() {
    mplapackint n = 2, lda = n, ldb = n, info, lwork = -1;
    mpfc_class *a = new mpfc_class[n * n];
    mpfc_class *b = new mpfc_class[n * n];
    mpfc_class *aorg = new mpfc_class[n * n];
    mpfc_class *borg = new mpfc_class[n * n];
    mpf_class *w = new mpf_class[n];
    mpf_class *rwork = new mpf_class[3 * n];
    for (mplapackint i = 0; i < n * n; i++) {
        a[i] = mpfc_class(0.0, 0.0);
        b[i] = mpfc_class(0.0, 0.0);
    }
    a[0] = mpfc_class(2.0, 0.0);
    a[3] = mpfc_class(6.0, 0.0);
    b[0] = mpfc_class(1.0, 0.0);
    b[3] = mpfc_class(2.0, 0.0);
    for (mplapackint i = 0; i < n * n; i++) {
        aorg[i] = a[i];
        borg[i] = b[i];
    }
    mpfc_class wk;
    Chegv((mplapackint)1, "V", "U", n, a, lda, b, ldb, w, &wk, lwork, rwork, info);
    lwork = castINTEGER_gmp(wk.real());
    mpfc_class *work = new mpfc_class[lwork];
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
