//public domain
#include <mpblas_binary128.h>
#include <mplapack_binary128.h>
#include <iostream>
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
bool rselect(mplapack_binary128_t ar, mplapack_binary128_t ai) {
    // sorting rule for eigenvalues.
    return false;
}

#include <iostream>
using namespace std;

int main() {
    mplapackint n;

    string str;
    getline(cin, str);
    cout << str << endl;
    getline(cin, str);
    stringstream ss(str);
    ss >> n;
    printf("# n %d\n", (int)n);

    mplapack_binary128_t *a = new mplapack_binary128_t[n * n];
    mplapack_binary128_t *vl = new mplapack_binary128_t[n * n];
    mplapack_binary128_t *vr = new mplapack_binary128_t[n * n];
    mplapackint sdim = 0;
    mplapackint lwork = 4 * n;
    mplapack_binary128_t *wr = new mplapack_binary128_t[n];
    mplapack_binary128_t *wi = new mplapack_binary128_t[n];
    mplapack_binary128_t *work = new mplapack_binary128_t[lwork];
    mplapackint info;
    double dtmp;
    for (int i = 0; i < n; i++) {
        getline(cin, str);
        stringstream ss(str);
        for (int j = 0; j < n; j++) {
            ss >> dtmp;
            a[i + j * n] = dtmp;
        }
    }
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Rgeev("V", "V", n, a, n, wr, wi, vl, n, vr, n, work, lwork, info);
    printf("# left vectors\n");
    printf("vl ="); printmat(n, n, vl, n); printf("\n");
    printf("# left vectors\n");
    printf("vr ="); printmat(n, n, vr, n); printf("\n");
    printf("[vl, d, w] = eig(a)\n");
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
