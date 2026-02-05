//public domain
#include <mpblas__Float128.h>
#include <mplapack__Float128.h>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <algorithm>

#define BUFLEN 1024

void printnum(_Float128 rtmp)
{
    int width = 42;
    char buf[BUFLEN];
#if defined ___MPLAPACK_WANT_LIBQUADMATH___
    int n = quadmath_snprintf (buf, sizeof buf, "%*.35Qe", width, rtmp);
#elif defined ___MPLAPACK_LONGDOUBLE_IS_BINARY128___
    snprintf (buf, sizeof buf, "%.35Le", rtmp);
#else
    strfromf128(buf, sizeof(buf), "%.35e", rtmp);
#endif
    if (rtmp >= 0.0)
        printf ("+%s", buf);
    else
        printf ("%s", buf);
    return;
}

//Matlab/Octave format
void printvec(_Float128 *a, int len) {
    _Float128 tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = a[i];
        printnum(tmp);
        if (i < len - 1)
            printf(", ");
    }
    printf("]");
}

void printmat(int n, int m, _Float128 *a, int lda)
{
    _Float128 mtmp;

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
bool rselect(_Float128 ar, _Float128 ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 100;

    _Float128 *a = new _Float128[n * n];
    _Float128 *vs = new _Float128[n * n];
    mplapackint sdim = 0;
    mplapackint lwork = 3 * n;
    _Float128 *wr = new _Float128[n];
    _Float128 *wi = new _Float128[n];
    _Float128 *work = new _Float128[lwork];
    bool bwork[n];
    mplapackint info;
    // setting A matrix
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            a[(i - 1) + (j - 1) * n] = 0.0;
            if ( i <= j && j <= i + 3 ) a[(i - 1) + (j - 1) * n] = 1.0;
            if ( i - 1 == j ) a[(i - 1) + (j - 1) * n] = -1.0;
	}
    }
    printf("# octave check\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Rgees("V", "N", rselect, n, a, n, sdim, wr, wi, vs, n, work, lwork, bwork, info);
    printf("vs ="); printmat(n, n, vs, n); printf("\n");
    printf("vs*vs'\n");
    printf("eig(a)\n");
    for (int i = 1; i <= n; i = i + 1) {
        printf("w_%d = ", (int)i); printnum(wr[i - 1]); printf(" "); printnum(wi[i - 1]); printf("i\n");
    }
    delete[] work;
    delete[] wr;
    delete[] wi;
    delete[] vs;
    delete[] a;
}
