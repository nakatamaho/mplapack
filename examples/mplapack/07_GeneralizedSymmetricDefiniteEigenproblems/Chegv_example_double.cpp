//public domain
#include <mpblas_double.h>
#include <mplapack_double.h>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <algorithm>

#define DOUBLE_FORMAT "%+20.16e"
#define DOUBLE_SHORT_FORMAT "%+20.16e"

inline void printnum(double rtmp) { printf(DOUBLE_FORMAT, rtmp); }
inline void printnum(std::complex<double> ctmp) { printf(DOUBLE_FORMAT DOUBLE_FORMAT "i", ctmp.real(), ctmp.imag()); }

//Matlab/Octave format
template <class X> void printvec(X *a, int len) {
    X tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

template <class X> void printmat(int n, int m, X *a, int lda)
{
    X mtmp;

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
double max_eigen_residual(mplapackint n, std::complex<double> *a, std::complex<double> *b, double lambda, std::complex<double> *z) {
    double err = 0.0;
    for (mplapackint i = 0; i < n; i++) {
        std::complex<double> s = std::complex<double>(0.0, 0.0), t = std::complex<double>(0.0, 0.0);
        for (mplapackint j = 0; j < n; j++) { s = s + a[i + j * n] * z[j]; t = t + b[i + j * n] * z[j]; }
        double d = abs(s - lambda * t);
        if (err < d) err = d;
    }
    return err;
}

int main(){ mplapackint n=2,lda=n,ldb=n,info,lwork=-1; std::complex<double> *a=new std::complex<double>[n*n]; std::complex<double> *b=new std::complex<double>[n*n]; std::complex<double> *aorg=new std::complex<double>[n*n]; std::complex<double> *borg=new std::complex<double>[n*n]; double *w=new double[n]; double *rwork=new double[3*n]; for(mplapackint i=0;i<n*n;i++){a[i]=std::complex<double>(0.0,0.0);b[i]=std::complex<double>(0.0,0.0);} a[0]=std::complex<double>(2.0,0.0); a[3]=std::complex<double>(6.0,0.0); b[0]=std::complex<double>(1.0,0.0); b[3]=std::complex<double>(2.0,0.0); for(mplapackint i=0;i<n*n;i++){aorg[i]=a[i]; borg[i]=b[i];} std::complex<double> wk; Chegv((mplapackint)1,"V","U",n,a,lda,b,ldb,w,&wk,lwork,rwork,info); lwork=castINTEGER_double(wk.real()); std::complex<double> *work=new std::complex<double>[lwork]; Chegv((mplapackint)1,"V","U",n,a,lda,b,ldb,w,work,lwork,rwork,info); printf("eigenvalues = "); printvec(w,n); printf("\n"); printf("eigenvectors = "); printmat(n,n,a,lda); printf("\n"); for(mplapackint j=0;j<n;j++){ printf("residual[%ld] = ",(long)j); printnum(max_eigen_residual(n,aorg,borg,w[j],&a[j*lda])); printf("\n"); } delete[] work; delete[] rwork; delete[] w; delete[] borg; delete[] aorg; delete[] b; delete[] a; return info!=0?1:0; }
