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
void sort_real(mplapackint n, mpf_class *x) {
    for (mplapackint i = 0; i < n; i++)
        for (mplapackint j = i + 1; j < n; j++)
            if (x[j] < x[i]) {
                mpf_class t = x[i];
                x[i] = x[j];
                x[j] = t;
            }
}
int main() {
    mplapackint n = 20, lda = n, ldv = 1, info, lwork = -1;
    mpf_class *coef = new mpf_class[n + 1];
    for (mplapackint i = 0; i <= n; i++)
        coef[i] = 0;
    coef[0] = 1;
    for (mplapackint k = 1; k <= n; k++) {
        for (mplapackint j = k; j >= 1; j--)
            coef[j] = coef[j] - mpf_class(k) * coef[j - 1];
    }
    mpf_class *a = new mpf_class[n * n];
    for (mplapackint i = 0; i < n * n; i++)
        a[i] = 0;
    for (mplapackint i = 1; i < n; i++)
        a[i + (i - 1) * lda] = 1;
    for (mplapackint j = 0; j < n; j++)
        a[j + (n - 1) * lda] = -coef[n - j] / coef[0];
    mpf_class *wr = new mpf_class[n];
    mpf_class *wi = new mpf_class[n];
    mpf_class *vl = new mpf_class[1];
    mpf_class *vr = new mpf_class[1];
    mpf_class wk;
    Rgeev("N", "N", n, a, lda, wr, wi, vl, ldv, vr, ldv, &wk, lwork, info);
    lwork = castINTEGER_gmp(wk);
    mpf_class *work = new mpf_class[lwork];
    Rgeev("N", "N", n, a, lda, wr, wi, vl, ldv, vr, ldv, work, lwork, info);
    sort_real(n, wr);
    mpf_class maxerr = 0;
    for (mplapackint i = 0; i < n; i++) {
        mpf_class err = abs(wr[i] - mpf_class(i + 1));
        if (maxerr < err)
            maxerr = err;
        printf("root[%ld] = ", (long)i); printnum(wr[i]); printf(", error = "); printnum(err); printf("\n");
    }
    printf("max root error = "); printnum(maxerr); printf("\n");
    delete[] work;
    delete[] vr;
    delete[] vl;
    delete[] wi;
    delete[] wr;
    delete[] a;
    delete[] coef;
    return info != 0 ? 1 : 0;
}
