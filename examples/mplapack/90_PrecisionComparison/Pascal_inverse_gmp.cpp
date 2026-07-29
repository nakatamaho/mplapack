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
mpf_class binom(mplapackint n, mplapackint k) {
    mpf_class r = 1;
    for (mplapackint i = 1; i <= k; i++)
        r = r * mpf_class(n - k + i) / mpf_class(i);
    return r;
}
mpf_class nearest_integer_error(mpf_class x) {
    mplapackint nearest = castINTEGER_gmp(x >= mpf_class(0.0) ? x + mpf_class(0.5) : x - mpf_class(0.5));
    return abs(x - mpf_class(nearest));
}
int main() {
    mplapackint n = 8, lda = n, info, lwork = -1;
    mpf_class *a = new mpf_class[n * n];
    mplapackint *ipiv = new mplapackint[n];
    for (mplapackint j = 0; j < n; j++)
        for (mplapackint i = 0; i < n; i++)
            a[i + j * lda] = binom(i + j, i);
    Rgetrf(n, n, a, lda, ipiv, info);
    mpf_class wk;
    if (info == 0)
        Rgetri(n, a, lda, ipiv, &wk, lwork, info);
    lwork = castINTEGER_gmp(wk);
    mpf_class *work = new mpf_class[lwork];
    if (info == 0)
        Rgetri(n, a, lda, ipiv, work, lwork, info);
    mpf_class err = 0;
    for (mplapackint i = 0; i < n * n; i++) {
        mpf_class d = nearest_integer_error(a[i]);
        if (err < d)
            err = d;
    }
    printf("P inverse = "); printmat(n, n, a, lda); printf("\n");
    printf("max distance to integer = "); printnum(err); printf("\n");
    delete[] work;
    delete[] ipiv;
    delete[] a;
    return info != 0 ? 1 : 0;
}
