//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_dd.h>
#include <mplapack_dd.h>

#define DD_PRECISION_SHORT 16

inline void printnum(dd_real rtmp) {
    std::cout.precision(DD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

//Matlab/Octave format
void printvec(dd_real *a, int len) {
    dd_real tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, dd_real * a, int lda)
{
    dd_real mtmp;
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
dd_real max_solution_error(mplapackint n, dd_real *x, dd_real *xexact){ dd_real err=0; for(mplapackint i=0;i<n;i++){ dd_real d=abs(x[i]-xexact[i]); if(err<d) err=d; } return err; }
int main(){ mplapackint n=15,lda=n,ldb=n,info; dd_real *a=new dd_real[n*n]; dd_real *b=new dd_real[n]; dd_real *xexact=new dd_real[n]; mplapackint *ipiv=new mplapackint[n]; for(mplapackint j=0;j<n;j++) xexact[j]=dd_real(j%3-1); for(mplapackint i=0;i<n;i++){ dd_real node=dd_real(i+1); dd_real p=1; for(mplapackint j=0;j<n;j++){ a[i+j*lda]=p; p=p*node; }} for(mplapackint i=0;i<n;i++){ b[i]=0; for(mplapackint j=0;j<n;j++) b[i]=b[i]+a[i+j*lda]*xexact[j]; } Rgesv(n,(mplapackint)1,a,lda,ipiv,b,ldb,info); printf("max |x-x_exact| = "); printnum(max_solution_error(n,b,xexact)); printf("\n"); delete[] ipiv; delete[] xexact; delete[] b; delete[] a; return info!=0?1:0; }
