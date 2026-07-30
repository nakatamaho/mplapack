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
    mplapackint n = 5, nrhs = 1, lda = n, ldb = n, info;
    mpf_class *a = new mpf_class[n * n];
    mpf_class *af = new mpf_class[n * n];
    mpf_class *b = new mpf_class[n];
    mpf_class *x = new mpf_class[n];
    mpf_class *xexact = new mpf_class[n];
    mpf_class *r = new mpf_class[n];
    mpf_class *c = new mpf_class[n];
    mpf_class *ferr = new mpf_class[nrhs];
    mpf_class *berr = new mpf_class[nrhs];
    mpf_class *work = new mpf_class[4 * n];
    mplapackint *iwork = new mplapackint[n];
    mplapackint *ipiv = new mplapackint[n];
    char equed[2];
    equed[0] = 'N';
    equed[1] = '\0';
    mpf_class rcond;
    for (mplapackint j = 0; j < n; j++) for (mplapackint i = 0; i < n; i++) a[i + j * lda] = mpf_class(1.0) / mpf_class(i + j + 1);
    for (mplapackint i = 0; i < n; i++) xexact[i] = (i % 2 == 0) ? mpf_class(1.0) : mpf_class(-1.0);
    for (mplapackint i = 0; i < n; i++) {
        b[i] = 0.0;
        for (mplapackint k = 0; k < n; k++)
            b[i] = b[i] + a[i + k * lda] * xexact[k];
    }
    printf("A = "); printmat(n, n, a, lda); printf("\n");
    printf("b = "); printvec(b, n); printf("\n");
    Rgesvx("N", "N", n, nrhs, a, lda, af, lda, ipiv, equed, r, c, b, ldb, x, ldb, rcond, ferr, berr, work, iwork, info);
    printf("x = "); printvec(x, n); printf("\n");
    printf("rcond = "); printnum(rcond); printf("\n");
    printf("ferr = "); printvec(ferr, nrhs); printf("\n");
    printf("berr = "); printvec(berr, nrhs); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n, nrhs, x, ldb, xexact, n)); printf("\n");
    delete[] ipiv;
    delete[] iwork;
    delete[] work;
    delete[] berr;
    delete[] ferr;
    delete[] c;
    delete[] r;
    delete[] xexact;
    delete[] x;
    delete[] b;
    delete[] af;
    delete[] a;
    return info != 0 ? 1 : 0;
}
