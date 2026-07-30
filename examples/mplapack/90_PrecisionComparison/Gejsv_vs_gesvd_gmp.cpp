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
int main() {
    mplapackint m = 8, n = 4, lda = m, info, lwork = -1;
    mpf_class *a = new mpf_class[m * n];
    mpf_class *b = new mpf_class[m * n];
    for (mplapackint j = 0; j < n; j++)
        for (mplapackint i = 0; i < m; i++) {
            mpf_class scale = pow(mpf_class(10.0), mpf_class((double)(i - 4)));
            a[i + j * lda] = scale / (mpf_class((double)(i + j + 1)));
            b[i + j * lda] = a[i + j * lda];
        }
    mpf_class *s1 = new mpf_class[n];
    mpf_class *s2 = new mpf_class[n];
    mpf_class *u = new mpf_class[1];
    mpf_class *vt = new mpf_class[1];
    mpf_class wk;
    Rgesvd("N", "N", m, n, a, lda, s1, u, (mplapackint)1, vt, (mplapackint)1, &wk, lwork, info);
    lwork = castINTEGER_gmp(wk);
    mpf_class *work = new mpf_class[lwork];
    Rgesvd("N", "N", m, n, a, lda, s1, u, (mplapackint)1, vt, (mplapackint)1, work, lwork, info);
    delete[] work;
    mplapackint *iwork = new mplapackint[m + 3 * n + 10];
    lwork = 5 * n * n + 9 * n + m;
    if (lwork < 2 * m + n)
        lwork = 2 * m + n;
    if (lwork < 4 * n + n * n)
        lwork = 4 * n + n * n;
    if (lwork < 7)
        lwork = 7;
    work = new mpf_class[lwork];
    Rgejsv("G", "N", "N", "N", "N", "N", m, n, b, lda, s2, u, (mplapackint)1, vt, (mplapackint)1, work, lwork, iwork, info);
    printf("Rgesvd singular values = "); printvec(s1, n); printf("\n");
    printf("Rgejsv singular values = "); printvec(s2, n); printf("\n");
    for (mplapackint i = 0; i < n; i++) {
        mpf_class rel = abs(s1[i] - s2[i]) / (abs(s2[i]) + Rlamch_gmp("S"));
        printf("relative difference[%ld] = ", (long)i); printnum(rel); printf("\n");
    }
    delete[] work;
    delete[] iwork;
    delete[] vt;
    delete[] u;
    delete[] s2;
    delete[] s1;
    delete[] b;
    delete[] a;
    return info != 0 ? 1 : 0;
}
