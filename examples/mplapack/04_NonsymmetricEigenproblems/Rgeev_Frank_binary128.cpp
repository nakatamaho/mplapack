//public domain
#include <mpblas_binary128.h>
#include <mplapack_binary128.h>
#include <iostream>
#include <sstream>
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
void Frank(mplapackint n) {
    mplapackint lwork, liwork, info, m;
    mplapack_binary128_t *a = new mplapack_binary128_t[n * n];
    mplapack_binary128_t *vl = new mplapack_binary128_t[n * n]; //not used
    mplapack_binary128_t *vr = new mplapack_binary128_t[n * n]; //not used
    mplapack_binary128_t *wr = new mplapack_binary128_t[n];
    mplapack_binary128_t *wi = new mplapack_binary128_t[n];

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
    mplapack_binary128_t *work = new mplapack_binary128_t[1];
    Rgeev("N", "N", n, a, n, wr, wi, vl, n, vr, n, work, lwork, info);
    lwork = (int)cast2double(work[0]);
    delete[] work;
    work = new mplapack_binary128_t[std::max((mplapackint)1, lwork)];

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

bool rselect(mplapack_binary128_t ar, mplapack_binary128_t ai) {
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
