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
bool cselect(dd_complex a) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 4;

    dd_complex *a = new dd_complex[n * n];
    mplapackint sdim = 0;
    mplapackint lwork = 2 * n;
    dd_complex *w = new dd_complex[n];
    dd_complex *vs = new dd_complex[n * n];
    dd_complex *work = new dd_complex[lwork];
    dd_real *rwork = new dd_real[n];
    bool bwork[n];
    mplapackint info;

    // setting A matrix
    a[0 + 0 * n] = dd_complex(3.0,  0.0); a[0 + 1 * n] = dd_complex(1.0, 0.0);   a[0 + 2 * n] = dd_complex(0.0, 0.0);  a[0 + 3 * n] = dd_complex(0.0, 2.0);
    a[1 + 0 * n] = dd_complex(1.0,  0.0); a[1 + 1 * n] = dd_complex(3.0, 0.0);   a[1 + 2 * n] = dd_complex(0.0, -2.0); a[1 + 3 * n] = dd_complex(0.0, 0.0);
    a[2 + 0 * n] = dd_complex(0.0,  0.0); a[2 + 1 * n] = dd_complex(0.0, 2.0);   a[2 + 2 * n] = dd_complex(1.0, 0.0);  a[2 + 3 * n] = dd_complex(1.0, 0.0);
    a[3 + 0 * n] = dd_complex(0.0, -2.0); a[3 + 1 * n] = dd_complex(0.0, 0.0);   a[3 + 2 * n] = dd_complex(1.0, 0.0);  a[3 + 3 * n] = dd_complex(1.0, 0.0); 

    printf("# Ex. 6.6 p. 116, Collection of Matrices for Testing Computational Algorithms, Robert T. Gregory, David L. Karney\n");
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Cgees("V", "S", cselect, n, a, n, sdim, w, vs, n, work, lwork, rwork, bwork, info);
    printf("w ="); printvec(w, n); printf("\n");
    printf("vs ="); printmat(n, n, vs, n); printf("\n");
    printf("t ="); printmat(n, n, a, n); printf("\n");
    printf("vs*t*vs'\n");
    printf("eig(a)\n");

    delete[] rwork;
    delete[] work;
    delete[] vs;
    delete[] w;
    delete[] a;
}
