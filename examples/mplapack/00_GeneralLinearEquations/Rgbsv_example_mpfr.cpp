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

int main() {
    mplapackint n = 6, kl = 2, ku = 1, nrhs = 1, ldab = 2 * kl + ku + 1, ldb = n, info;
    mpfr_class *a = new mpfr_class[n * n];
    mpfr_class *ab = new mpfr_class[ldab * n];
    mpfr_class *b = new mpfr_class[n];
    mpfr_class *borg = new mpfr_class[n];
    mpfr_class *xexact = new mpfr_class[n];
    mplapackint *ipiv = new mplapackint[n];
    for (mplapackint i = 0; i < n * n; i++) a[i] = 0.0;
    for (mplapackint i = 0; i < ldab * n; i++) ab[i] = 0.0;
    for (mplapackint i = 0; i < n; i++) {
        xexact[i] = i + 1.0;
        for (mplapackint j = 0; j < n; j++) {
            if (i == j) a[i + j * n] = 6.0;
            else if (i > j && i - j <= kl) a[i + j * n] = -1.0;
            else if (j > i && j - i <= ku) a[i + j * n] = 2.0;
        }
    }
    for (mplapackint j = 0; j < n; j++) for (mplapackint i = 0; i < n; i++) if (a[i + j * n] != 0.0) {
        mplapackint row = kl + ku + i - j;
        ab[row + j * ldab] = a[i + j * n];
    }
    for (mplapackint i = 0; i < n; i++) {
        b[i] = 0.0;
        for (mplapackint k = 0; k < n; k++) b[i] = b[i] + a[i + k * n] * xexact[k];
        borg[i] = b[i];
    }
    printf("A = "); printmat(n, n, a, n); printf("\n");
    printf("AB = "); printmat(ldab, n, ab, ldab); printf("\n");
    Rgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
    printf("x = "); printvec(b, n); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, b, ldb, xexact, n)); printf("\n");
    printf("max residual = "); printnum(max_residual(n, n, nrhs, a, n, b, ldb, borg, ldb)); printf("\n");
    delete[] ipiv;
    delete[] xexact;
    delete[] borg;
    delete[] b;
    delete[] ab;
    delete[] a;
    return info != 0 ? 1 : 0;
}
