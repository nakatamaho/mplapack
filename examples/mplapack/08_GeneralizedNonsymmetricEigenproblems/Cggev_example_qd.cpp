//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_qd.h>
#include <mplapack_qd.h>

#define QD_PRECISION_SHORT 16

inline void printnum(qd_real rtmp) {
    std::cout.precision(QD_PRECISION_SHORT);
    if (rtmp >= 0.0) {
        std::cout << "+" << rtmp;
    } else {
        std::cout << rtmp;
    }
    return;
}

inline void printnum(qd_complex rtmp) {
    std::cout.precision(QD_PRECISION_SHORT);
    if (rtmp.real() >= 0.0) {
        std::cout << "+" << rtmp.real();
    } else {
        std::cout << rtmp.real();
    }
    if (rtmp.imag() >= 0.0) {
        std::cout << "+" << rtmp.imag() << "i";
    } else {
        std::cout << rtmp.imag() << "i";
    }
    return;
}

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
int main(){ mplapackint n=3,lda=n,ldb=n,ldv=1,info,lwork=-1; qd_complex *a=new qd_complex[n*n]; qd_complex *b=new qd_complex[n*n]; qd_complex *alpha=new qd_complex[n]; qd_complex *beta=new qd_complex[n]; qd_complex *vl=new qd_complex[1]; qd_complex *vr=new qd_complex[1]; qd_real *rwork=new qd_real[8*n]; for(mplapackint i=0;i<n*n;i++){a[i]=qd_complex(0.0,0.0);b[i]=qd_complex(0.0,0.0);} a[0]=qd_complex(1.0,1.0); a[4]=qd_complex(2.0,0.0); a[8]=qd_complex(3.0,-1.0); b[0]=qd_complex(1.0,0.0); b[4]=qd_complex(1.0,0.0); b[8]=qd_complex(0.0,0.0); qd_complex wk; Cggev("N","N",n,a,lda,b,ldb,alpha,beta,vl,ldv,vr,ldv,&wk,lwork,rwork,info); lwork=castINTEGER_qd(wk.real()); qd_complex *work=new qd_complex[lwork]; Cggev("N","N",n,a,lda,b,ldb,alpha,beta,vl,ldv,vr,ldv,work,lwork,rwork,info); for(mplapackint i=0;i<n;i++){ printf("alpha[%ld] = ",(long)i); printnum(alpha[i]); printf(", beta = "); printnum(beta[i]); if(abs(beta[i])<=Rlamch_qd("E")) printf(", lambda = Inf\n"); else { printf(", lambda = "); printnum(alpha[i]/beta[i]); printf("\n"); }} delete[] work; delete[] rwork; delete[] vr; delete[] vl; delete[] beta; delete[] alpha; delete[] b; delete[] a; return info!=0?1:0; }
