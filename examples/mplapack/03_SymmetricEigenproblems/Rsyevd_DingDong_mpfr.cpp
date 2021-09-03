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

// Matlab/Octave format
void printvec(mpreal *a, int len) {
    mpreal tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mpreal *a, int lda) {
    mpreal mtmp;
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

//https://math.nist.gov/MatrixMarket/deli/DingDong/
//J.C. Nash, Compact Numerical Methods for Computers: Linear Algebra and Function Minimisation, second edition, Adam Hilger, Bristol, 1990 (Appendix 1). 

void DingDong(mplapackint n) {
    mplapackint lwork, liwork, info, m;
    mpreal *a = new mpreal[n * n];
    mpreal *w = new mpreal[n];
    mpreal PI;
    PI = pi(PI);

    // setting A matrix
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            a[(i - 1) + (j - 1) * n] = 1.0 / mpreal( 2.0 * ( n - i - j + 3.0 / 2.0 ));
        }
    }
    printf("a ="); printmat(n, n, a, n); printf("\n");

    // work space query
    lwork = -1;
    mpreal *work = new mpreal[1];
    liwork = -1;
    mplapackint *iwork = new mplapackint[1];

    Rsyevd("N", "U", n, a, n, w, work, lwork, iwork, liwork, info);
    lwork = (int)cast2double(work[0]);
    delete[] work;
    work = new mpreal[std::max((mplapackint)1, lwork)];
    liwork = iwork[0];
    delete[] iwork;
    iwork = new mplapackint[std::max((mplapackint)1, liwork)];

    // diagonalize matrix
    Rsyevd("N", "U", n, a, n, w, work, lwork, iwork, liwork, info);

    // print out
    printf("#eigenvalues \n");
    printf("w ="); printvec(w, n); printf("\n");
    printf("w_smallest ="); printnum(w[0]); printf("\n");
    printf("w_largest  ="); printnum(w[n-1]); printf("\n");

    printf("w_relerror_to_halfPI ="); printnum( (w[n-1] - PI / 2.0) /  (PI / 2.0) ); printf("\n");

    delete[] iwork;
    delete[] work;
    delete[] w;
    delete[] a;
}

int main(int argc, char *argv[]) {
    int STARTN = 5;
    int ENDN = 1000;
    int STEPN = 1;
    if (argc != 1) {
        for (int i = 1; i < argc; i++) {
            if (strcmp("-STEPN", argv[i]) == 0) {
                STEPN = atoi(argv[++i]);
            } else if (strcmp("-STARTN", argv[i]) == 0) {
                STARTN = atoi(argv[++i]);
            } else if (strcmp("-ENDN", argv[i]) == 0) {
                ENDN = atoi(argv[++i]);
            }
        }
    }
    for (int n = STARTN; n <= ENDN; n = n + STEPN) {
        printf("# Eigenvalues of DingDong matrix of order n=%d\n", n);
        DingDong((mplapackint)n);
    }
}
