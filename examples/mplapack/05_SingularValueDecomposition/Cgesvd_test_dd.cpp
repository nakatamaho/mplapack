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

inline void printnum(dd_complex rtmp) {
    std::cout.precision(DD_PRECISION_SHORT);
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
int main() {
    mplapackint n = 4;
    mplapackint m = 4;

    dd_complex *a = new dd_complex[m * n];
    dd_real *s = new dd_real[std::min(m, n)];
    dd_complex *u = new dd_complex[m * m];
    dd_complex *vt = new dd_complex[n * n];
    mplapackint lwork = std::max((mplapackint)1, 2 * std::min(m, n) + std::max(m, n));
    dd_complex *work = new dd_complex[lwork];
    dd_real *rwork = new dd_real[5 * std::min(m, n)];
    mplapackint info;

    // setting A matrix
    a[0 + 0 * n] = dd_complex(0.9, -1.0); a[0 + 1 * n] = dd_complex(20.0, -2.25);  a[0 + 2 * n] = dd_complex(1.75, -0.5);  a[0 + 3 * n] = dd_complex(0.0, 0.5);
    a[1 + 0 * n] = dd_complex(8.0,-2.25); a[1 + 1 * n] = dd_complex(-0.25, 0.0);   a[1 + 2 * n] = dd_complex(1.25, -0.25); a[1 + 3 * n] = dd_complex(-3.75, 0.0);
    a[2 + 0 * n] = dd_complex(-1.75,0.0); a[2 + 1 * n] = dd_complex(-80.0,  1.25); a[2 + 2 * n] = dd_complex(1.5, 0.0);    a[2 + 3 * n] = dd_complex(30.0, 2.25);
    a[3 + 0 * n] = dd_complex(3.0, 0.25); a[3 + 1 * n] = dd_complex(1.75, 0.0);    a[3 + 2 * n] = dd_complex(0.0, 2.25);   a[3 + 3 * n] = dd_complex(-0.25, -80.0);
    
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(m, n, a, m); printf("\n");
    Cgesvd("A", "A", m, n, a, m, s, u, m, vt, n, work, lwork, rwork, info);
    printf("s="); printvec(s, std::min(m, n)); printf("\n");
    if (m < n)
        printf("padding=zeros(%d, %d-%d)\n", (int)m, (int)n, (int)m);
    if (n < m)
        printf("padding=zeros(%d-%d,%d)\n", (int)m, (int)n, (int)n);
    printf("u ="); printmat(m, m, u, m); printf("\n");
    printf("vt ="); printmat(n, n, vt, n); printf("\n");
    printf("svd(a)\n");
    if (m < n)
        printf("sigma=[diag(s) padding] \n");
    if (n < m)
        printf("sigma=[diag(s); padding] \n");
    if (n == m)
        printf("sigma=[diag(s)] \n");
    printf("sigma \n");
    printf("u * sigma  * vt\n");
    delete[] rwork;
    delete[] work;
    delete[] vt;
    delete[] u;
    delete[] s;
    delete[] a;
}
