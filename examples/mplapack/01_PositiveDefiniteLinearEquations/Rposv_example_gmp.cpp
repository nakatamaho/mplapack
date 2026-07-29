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

int main() {
    mplapackint n = 2, nrhs = 1, lda = n, ldb = n, info;
    mpf_class *a = new mpf_class[n * n]; mpf_class *aorg = new mpf_class[n * n]; mpf_class *b = new mpf_class[n]; mpf_class *borg = new mpf_class[n]; mpf_class *xexact = new mpf_class[n];
    a[0]=4; a[1]=2; a[2]=2; a[3]=5;
    xexact[0]=1; xexact[1]=2;
    for (mplapackint i=0;i<n*n;i++) aorg[i]=a[i];
    for (mplapackint i=0;i<n;i++){ b[i]=0; for(mplapackint k=0;k<n;k++) b[i]+=aorg[i+k*lda]*xexact[k]; borg[i]=b[i]; }
    printf("A = "); printmat(n,n,aorg,lda); printf("\n");
    Rposv("L", n, nrhs, a, lda, b, ldb, info);
    printf("x = "); printvec(b,n); printf("\n");
    printf("max |x-x_exact| = "); printnum(max_solution_error(n,nrhs,b,ldb,xexact,n)); printf("\n");
    printf("max residual = "); printnum(max_residual(n,n,nrhs,aorg,lda,b,ldb,borg,ldb)); printf("\n");
    delete[] xexact; delete[] borg; delete[] b; delete[] aorg; delete[] a; return info != 0 ? 1 : 0;
}
