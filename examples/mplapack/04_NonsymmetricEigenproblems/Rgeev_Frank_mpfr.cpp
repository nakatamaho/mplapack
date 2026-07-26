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

inline void printnum(mpfr_class rtmp) { mpfr_printf(MPFR_FORMAT, mpfr_ptr(rtmp)); }
inline void printnum_short(mpfr_class rtmp) { mpfr_printf(MPFR_SHORT_FORMAT, mpfr_ptr(rtmp)); }

// Matlab/Octave format
void printvec(mpfr_class *a, int len) {
    mpfr_class tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, mpfr_class *a, int lda) {
    mpfr_class mtmp;
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

void Frank(mplapackint n) {
    mplapackint lwork, liwork, info, m;
    mpfr_class *a = new mpfr_class[n * n];
    mpfr_class *vl = new mpfr_class[n * n]; //not used
    mpfr_class *vr = new mpfr_class[n * n]; //not used
    mpfr_class *wr = new mpfr_class[n];
    mpfr_class *wi = new mpfr_class[n];

    // setting A matrix
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            a[(i - 1) + (j - 1) * n] = 0.0;
	}
    }
    for (int i = 1; i <= n; i++) {
        for (int j = i; j <= n; j++) {
            a[(i - 1) + (j - 1) * n] = n - std::max(i, j) + 1;
        }
    }
    for (int i = 1; i <= n - 1 ; i++) {
        a[i + (i - 1) * n] = n - i;
    }
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");

    // work space query
    lwork = -1;
    mpfr_class *work = new mpfr_class[1];
    Rgeev("N", "N", n, a, n, wr, wi, vl, n, vr, n, work, lwork, info);
    lwork = (int)cast2double(work[0]);
    delete[] work;
    work = new mpfr_class[std::max((mplapackint)1, lwork)];

    // diagonalize matrix
    Rgeev("N", "N", n, a, n, wr, wi, vl, n, vr, n, work, lwork, info);

    // print out
    printf("#eigenvalues \n");
    for (int i = 1; i <= n; i = i + 1) {
        printf("w_%d = ", (int)i); printnum(wr[i - 1]); printf(" "); printnum(wi[i - 1]); printf("i\n");
    }

    delete[] work;
    delete[] wr;
    delete[] wi;
    delete[] vr;
    delete[] vl;
    delete[] a;
}

bool rselect(mpfr_class ar, mpfr_class ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main(int argc, char *argv[]) {
    int STARTN = 5;
    int ENDN = 25;
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
        printf("# Eigenvalues of Frank matrix of order n=%d\n", n);
        Frank((mplapackint)n);
    }
}
