//public domain
#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <mpblas_mpfr.h>
#include <mplapack_mpfr.h>

#define MPFR_FORMAT "%+68.64Re"
#define MPFR_SHORT_FORMAT "%+20.16Re"

inline void printnum(mpreal rtmp) { mpfr_printf(MPFR_FORMAT, mpfr_ptr(rtmp)); }
inline void printnum_short(mpreal rtmp) { mpfr_printf(MPFR_SHORT_FORMAT, mpfr_ptr(rtmp)); }
inline void printnum(mpcomplex ctmp) {
    mpreal cre, cim;
    cre = ctmp.real();
    cim = ctmp.imag();
    mpfr_printf(MPFR_SHORT_FORMAT MPFR_SHORT_FORMAT "i", mpfr_ptr(cre), mpfr_ptr(cim));
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
mpreal max_solution_error(mplapackint n, mplapackint nrhs, mpcomplex *x, mplapackint ldx, mpcomplex *xexact, mplapackint ldxexact) {
    mpreal err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < n; i++) {
            mpreal d = abs(x[i + j * ldx] - xexact[i + j * ldxexact]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

mpreal max_residual(mplapackint m, mplapackint n, mplapackint nrhs, mpcomplex *a, mplapackint lda, mpcomplex *x, mplapackint ldx, mpcomplex *b, mplapackint ldb) {
    mpreal err = 0.0;
    for (mplapackint j = 0; j < nrhs; j++) {
        for (mplapackint i = 0; i < m; i++) {
            mpcomplex s = mpcomplex(0.0, 0.0);
            for (mplapackint k = 0; k < n; k++)
                s = s + a[i + k * lda] * x[k + j * ldx];
            mpreal d = abs(s - b[i + j * ldb]);
            if (err < d)
                err = d;
        }
    }
    return err;
}

void set_problem(mplapackint m, mplapackint n, mpcomplex *a, mplapackint lda, mpcomplex *b, mplapackint ldb, mpcomplex *xexact) {
    for (mplapackint i = 0; i < m; i++) { a[i + 0 * lda] = mpcomplex(1.0, 0.0); a[i + 1 * lda] = mpcomplex(i, 1.0); }
    xexact[0] = mpcomplex(1.0, -1.0); xexact[1] = mpcomplex(2.0, 1.0);
    for (mplapackint i = 0; i < m; i++) b[i] = a[i + 0 * lda] * xexact[0] + a[i + 1 * lda] * xexact[1];
}

int main(){ mplapackint m=4,n=2,nrhs=1,lda=m,ldb=m,info,lwork=-1,rank; mpcomplex *a=new mpcomplex[lda*n]; mpcomplex *aorg=new mpcomplex[lda*n]; mpcomplex *b=new mpcomplex[ldb]; mpcomplex *borg=new mpcomplex[ldb]; mpcomplex *xexact=new mpcomplex[n]; mpreal *s=new mpreal[n]; mpreal *rwork=new mpreal[5*n]; set_problem(m,n,a,lda,b,ldb,xexact); for(mplapackint i=0;i<lda*n;i++) aorg[i]=a[i]; for(mplapackint i=0;i<ldb;i++) borg[i]=b[i]; mpcomplex wk; Cgelss(m,n,nrhs,a,lda,b,ldb,s,mpreal(-1.0),rank,&wk,lwork,rwork,info); lwork=castINTEGER_mpfr(wk.real()); mpcomplex *work=new mpcomplex[lwork]; Cgelss(m,n,nrhs,a,lda,b,ldb,s,mpreal(-1.0),rank,work,lwork,rwork,info); printf("singular values = "); printvec(s,n); printf("\n"); printf("rank = %ld\n", (long)rank); printf("x = "); printvec(b,n); printf("\n"); printf("max |x-x_exact| = "); printnum(max_solution_error(n,nrhs,b,ldb,xexact,n)); printf("\n"); printf("max residual = "); printnum(max_residual(m,n,nrhs,aorg,lda,b,ldb,borg,ldb)); printf("\n"); delete[] work; delete[] rwork; delete[] s; delete[] xexact; delete[] borg; delete[] b; delete[] aorg; delete[] a; return info!=0?1:0; }
