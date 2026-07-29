//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_mpfr.h>
#include <mplapack_mpfr.h>

#define MPFR_FORMAT "%+68.64Re"
#define MPFR_SHORT_FORMAT "%+20.16Re"

inline void printnum(mpfr_class rtmp) { mpfr_printf(MPFR_FORMAT, mpfr_ptr(rtmp)); }
inline void printnum_short(mpfr_class rtmp) { mpfr_printf(MPFR_SHORT_FORMAT, mpfr_ptr(rtmp)); }

// Matlab/Octave format
void printvec(mpfr_class *a, int len) {
    mpfr_class tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mpfr_class *a, int lda) {
    mpfr_class mtmp;
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

mpfr_class maxabs(mpfr_class a, mpfr_class b) {
    mpfr_class d = abs(a - b);
    return d;
}

mpfr_class max_solution_error(mplapackint n, mplapackint nrhs, mpfr_class *x, mplapackint ldx, mpfr_class *xexact, mplapackint ldxexact) {
    mpfr_class err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            mpfr_class d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

mpfr_class max_residual(mplapackint m, mplapackint n, mplapackint nrhs, mpfr_class *a, mplapackint lda, mpfr_class *x, mplapackint ldx, mpfr_class *b, mplapackint ldb) {
    mpfr_class err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            mpfr_class s = 0.0;
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            mpfr_class d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

mpfr_class one_norm(mplapackint n, mpfr_class *a, mplapackint lda) {
    mpfr_class anorm = 0.0;
    for (mplapackint j = 0; j < n; j++) {
        mpfr_class s = 0.0;
        for (mplapackint i = 0; i < n; i++)
            s = s + abs(a[i + j * lda]);
        if (anorm < s)
            anorm = s;
    }
    return anorm;
}
mpfr_class inf_norm(mplapackint n, mpfr_class *a, mplapackint lda) {
    mpfr_class anorm = 0.0;
    for (mplapackint i = 0; i < n; i++) {
        mpfr_class s = 0.0;
        for (mplapackint j = 0; j < n; j++)
            s = s + abs(a[i + j * lda]);
        if (anorm < s)
            anorm = s;
    }
    return anorm;
}
int main() {
    mplapackint n = 3, lda = n, info;
    mpfr_class *a = new mpfr_class[n * n];
    mpfr_class *lu = new mpfr_class[n * n];
    mpfr_class *work = new mpfr_class[4 * n];
    mplapackint *iwork = new mplapackint[n];
    mplapackint *ipiv = new mplapackint[n];
    for (mplapackint i = 0; i < n * n; i++) a[i] = 0;
    a[0 + 0 * lda] = 1; a[1 + 1 * lda] = mpfr_class(1.0e-3); a[2 + 2 * lda] = mpfr_class(1.0e-6);
    for (mplapackint i = 0; i < n * n; i++) lu[i] = a[i];
    Rgetrf(n, n, lu, lda, ipiv, info);
    mpfr_class rcond1 = 0.0, rcondi = 0.0;
    if (info == 0) Rgecon("1", n, lu, lda, one_norm(n, a, lda), rcond1, work, iwork, info);
    if (info == 0) Rgecon("I", n, lu, lda, inf_norm(n, a, lda), rcondi, work, iwork, info);
    printf("A = "); printmat(n, n, a, lda); printf("\n");
    printf("true cond_1 = "); printnum(mpfr_class(1.0e6)); printf("\n");
    printf("rcond_1 = "); printnum(rcond1); printf("\n");
    printf("rcond_inf = "); printnum(rcondi); printf("\n");
    delete[] ipiv;
    delete[] iwork;
    delete[] work;
    delete[] lu;
    delete[] a;
    return info != 0 ? 1 : 0;
}
