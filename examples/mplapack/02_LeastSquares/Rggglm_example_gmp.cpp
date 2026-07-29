//public domain
#include <sstream>
#include <mpblas_gmp.h>
#include <mplapack_gmp.h>

#define GMP_FORMAT "%+68.64Fe"
#define GMP_SHORT_FORMAT "%+20.16Fe"

inline void printnum(mpf_class rtmp) { gmp_printf(GMP_FORMAT, rtmp.get_mpf_t()); }
inline void printnum_short(mpf_class rtmp) { gmp_printf(GMP_SHORT_FORMAT, rtmp.get_mpf_t()); }

//Matlab/Octave format
void printvec(mpf_class *a, int len) {
    mpf_class tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mpf_class * a, int lda)
{
    mpf_class mtmp;

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
mpf_class maxabs(mpf_class a, mpf_class b) {
    mpf_class d = abs(a - b);
    return d;
}

mpf_class max_solution_error(mplapackint n, mplapackint nrhs, mpf_class *x, mplapackint ldx, mpf_class *xexact, mplapackint ldxexact) {
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

mpf_class max_residual(mplapackint m, mplapackint n, mplapackint nrhs, mpf_class *a, mplapackint lda, mpf_class *x, mplapackint ldx, mpf_class *b, mplapackint ldb) {
    mpf_class err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            mpf_class s = 0.0;
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
    mplapackint n = 3, m = 2, p = 3, lda = n, ldb = n, info, lwork = -1;
    mpf_class *a = new mpf_class[lda * m];
    mpf_class *b = new mpf_class[ldb * p];
    mpf_class *d = new mpf_class[n];
    mpf_class *x = new mpf_class[m];
    mpf_class *y = new mpf_class[p];
    mpf_class *xexact = new mpf_class[m];
    mpf_class *yexact = new mpf_class[p];
    for (mplapackint i = 0; i < lda * m; i++)
        a[i] = 0;
    a[0] = 1;
    a[1] = 0;
    a[2] = 1;
    a[0 + lda] = 0;
    a[1 + lda] = 1;
    a[2 + lda] = 1;
    for (mplapackint i = 0; i < ldb * p; i++)
        b[i] = 0;
    for (mplapackint i = 0; i < n; i++)
        b[i + i * ldb] = 1;
    xexact[0] = 1;
    xexact[1] = 2;
    yexact[0] = mpf_class(0.5);
    yexact[1] = mpf_class(-0.5);
    yexact[2] = 1;
    for (mplapackint i = 0; i < n; i++)
        d[i] = a[i] * xexact[0] + a[i + lda] * xexact[1] + yexact[i];
    mpf_class wk;
    Rggglm(n, m, p, a, lda, b, ldb, d, x, y, &wk, lwork, info);
    lwork = castINTEGER_gmp(wk);
    mpf_class *work = new mpf_class[lwork];
    Rggglm(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info);
    printf("x = "); printvec(x, m); printf("\n");
    printf("y = "); printvec(y, p); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(m, (mplapackint)1, x, m, xexact, m)); printf("\n");
    printf("max |y-y_exact| = "); printnum(max_solution_error(p, (mplapackint)1, y, p, yexact, p)); printf("\n");
    delete[] work;
    delete[] yexact;
    delete[] xexact;
    delete[] y;
    delete[] x;
    delete[] d;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
