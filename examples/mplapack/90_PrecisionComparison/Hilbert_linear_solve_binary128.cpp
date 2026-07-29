//public domain
#include <mpblas_binary128.h>
#include <mplapack_binary128.h>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <algorithm>

#define BUFLEN 1024

void printnum(mplapack_binary128_t rtmp)
{
    int width = 42;
    char buf[BUFLEN];
#if MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_STRFROMF128
    strfromf128(buf, sizeof(buf), "%.35e", rtmp);
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_QUADMATH_SNPRINTF
    int n = quadmath_snprintf (buf, sizeof buf, "%*.35Qe", width, rtmp);
#elif MPLAPACK_BINARY128_IO == MPLAPACK_BINARY128_IO_SNPRINTF_LDBL
    snprintf (buf, sizeof buf, "%.35Le", rtmp);
#else
    #error "unsupported binary128 type"
#endif
    if (rtmp >= 0.0)
        printf ("+%s", buf);
    else
        printf ("%s", buf);
    return;
}

//Matlab/Octave format
void printvec(mplapack_binary128_t *a, int len) {
    mplapack_binary128_t tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mplapack_binary128_t *a, int lda)
{
    mplapack_binary128_t mtmp;

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
mplapack_binary128_t max_solution_error(mplapackint n, mplapack_binary128_t *x, mplapack_binary128_t *xexact){ mplapack_binary128_t err=0; for(mplapackint i=0;i<n;i++){ mplapack_binary128_t d=abs(x[i]-xexact[i]); if(err<d) err=d; } return err; }
int main(){ mplapackint n=12,lda=n,ldb=n,info; mplapack_binary128_t *a=new mplapack_binary128_t[n*n]; mplapack_binary128_t *aorg=new mplapack_binary128_t[n*n]; mplapack_binary128_t *b=new mplapack_binary128_t[n]; mplapack_binary128_t *xexact=new mplapack_binary128_t[n]; mplapackint *ipiv=new mplapackint[n]; for(mplapackint j=0;j<n;j++) for(mplapackint i=0;i<n;i++){ a[i+j*lda]=mplapack_binary128_t(1.0)/mplapack_binary128_t(i+j+1); aorg[i+j*lda]=a[i+j*lda]; } for(mplapackint i=0;i<n;i++) xexact[i]=(i%2==0)?mplapack_binary128_t(1.0):mplapack_binary128_t(-1.0); for(mplapackint i=0;i<n;i++){ b[i]=0; for(mplapackint k=0;k<n;k++) b[i]=b[i]+aorg[i+k*lda]*xexact[k]; } printf("Hilbert n = %ld\n",(long)n); Rgesv(n,(mplapackint)1,a,lda,ipiv,b,ldb,info); printf("max |x-x_exact| = "); printnum(max_solution_error(n,b,xexact)); printf("\n"); delete[] ipiv; delete[] xexact; delete[] b; delete[] aorg; delete[] a; return info!=0?1:0; }
