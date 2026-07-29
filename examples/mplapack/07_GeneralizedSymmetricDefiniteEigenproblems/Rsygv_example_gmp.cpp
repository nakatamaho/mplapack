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
mpf_class max_eigen_residual(mplapackint n, mpf_class *a, mpf_class *b, mpf_class lambda, mpf_class *z) {
    mpf_class err = 0.0;
    for (mplapackint i = 0; i < n; i++) {
        mpf_class s = 0.0, t = 0.0;
        for (mplapackint j = 0; j < n; j++) { s = s + a[i + j * n] * z[j]; t = t + b[i + j * n] * z[j]; }
        mpf_class d = abs(s - lambda * t);
        if (err < d) err = d;
    }
    return err;
}

int main(){ mplapackint n=2,lda=n,ldb=n,info,lwork=-1; mpf_class *a=new mpf_class[n*n]; mpf_class *b=new mpf_class[n*n]; mpf_class *aorg=new mpf_class[n*n]; mpf_class *borg=new mpf_class[n*n]; mpf_class *w=new mpf_class[n]; for(mplapackint i=0;i<n*n;i++){a[i]=0;b[i]=0;} a[0]=2; a[3]=6; b[0]=1; b[3]=2; for(mplapackint i=0;i<n*n;i++){aorg[i]=a[i]; borg[i]=b[i];} mpf_class wk; Rsygv((mplapackint)1,"V","U",n,a,lda,b,ldb,w,&wk,lwork,info); lwork=castINTEGER_gmp(wk); mpf_class *work=new mpf_class[lwork]; Rsygv((mplapackint)1,"V","U",n,a,lda,b,ldb,w,work,lwork,info); printf("eigenvalues = "); printvec(w,n); printf("\n"); printf("eigenvectors = "); printmat(n,n,a,lda); printf("\n"); for(mplapackint j=0;j<n;j++){ printf("residual[%ld] = ",(long)j); printnum(max_eigen_residual(n,aorg,borg,w[j],&a[j*lda])); printf("\n"); } delete[] work; delete[] w; delete[] borg; delete[] aorg; delete[] b; delete[] a; return info!=0?1:0; }
