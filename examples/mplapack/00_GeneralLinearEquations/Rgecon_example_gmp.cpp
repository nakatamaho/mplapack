//public domain
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

mpf_class one_norm(mplapackint n, mpf_class *a, mplapackint lda) {
    mpf_class anorm = 0.0;
    for (mplapackint j = 0; j < n; j++) { mpf_class s = 0.0; for (mplapackint i = 0; i < n; i++) s = s + abs(a[i + j * lda]); if (anorm < s) anorm = s; }
    return anorm;
}
mpf_class inf_norm(mplapackint n, mpf_class *a, mplapackint lda) {
    mpf_class anorm = 0.0;
    for (mplapackint i = 0; i < n; i++) { mpf_class s = 0.0; for (mplapackint j = 0; j < n; j++) s = s + abs(a[i + j * lda]); if (anorm < s) anorm = s; }
    return anorm;
}
int main() {
    mplapackint n = 3, lda = n, info;
    mpf_class *a = new mpf_class[n * n];
    mpf_class *lu = new mpf_class[n * n];
    mpf_class *work = new mpf_class[4 * n];
    mplapackint *iwork = new mplapackint[n];
    mplapackint *ipiv = new mplapackint[n];
    for (mplapackint i = 0; i < n * n; i++) a[i] = 0;
    a[0 + 0 * lda] = 1; a[1 + 1 * lda] = mpf_class(1.0e-3); a[2 + 2 * lda] = mpf_class(1.0e-6);
    for (mplapackint i = 0; i < n * n; i++) lu[i] = a[i];
    Rgetrf(n, n, lu, lda, ipiv, info);
    mpf_class rcond1 = 0.0, rcondi = 0.0;
    if (info == 0) Rgecon("1", n, lu, lda, one_norm(n, a, lda), rcond1, work, iwork, info);
    if (info == 0) Rgecon("I", n, lu, lda, inf_norm(n, a, lda), rcondi, work, iwork, info);
    printf("A = "); printmat(n, n, a, lda); printf("\n");
    printf("true cond_1 = "); printnum(mpf_class(1.0e6)); printf("\n");
    printf("rcond_1 = "); printnum(rcond1); printf("\n");
    printf("rcond_inf = "); printnum(rcondi); printf("\n");
    delete[] ipiv; delete[] iwork; delete[] work; delete[] lu; delete[] a;
    return info != 0 ? 1 : 0;
}
