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
mpf_class max_solution_error(mplapackint n, mpf_class *x, mpf_class *xexact){ mpf_class err=0; for(mplapackint i=0;i<n;i++){ mpf_class d=abs(x[i]-xexact[i]); if(err<d) err=d; } return err; }
int main(){ mplapackint n=15,lda=n,ldb=n,info; mpf_class *a=new mpf_class[n*n]; mpf_class *b=new mpf_class[n]; mpf_class *xexact=new mpf_class[n]; mplapackint *ipiv=new mplapackint[n]; for(mplapackint j=0;j<n;j++) xexact[j]=mpf_class(j%3-1); for(mplapackint i=0;i<n;i++){ mpf_class node=mpf_class(i+1); mpf_class p=1; for(mplapackint j=0;j<n;j++){ a[i+j*lda]=p; p=p*node; }} for(mplapackint i=0;i<n;i++){ b[i]=0; for(mplapackint j=0;j<n;j++) b[i]=b[i]+a[i+j*lda]*xexact[j]; } Rgesv(n,(mplapackint)1,a,lda,ipiv,b,ldb,info); printf("max |x-x_exact| = "); printnum(max_solution_error(n,b,xexact)); printf("\n"); delete[] ipiv; delete[] xexact; delete[] b; delete[] a; return info!=0?1:0; }
