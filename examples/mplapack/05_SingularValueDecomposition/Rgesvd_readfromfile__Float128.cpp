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
    int n = quadmath_snprintf (buf, sizeof buf, "%+-#*.35Qe", width, rtmp);
#elif defined ___MPLAPACK_LONGDOUBLE_IS_BINARY128___
    snprintf (buf, sizeof buf, "%.35Le", rtmp);
#else
    strfromf128(buf, sizeof(buf), "%.35e", rtmp);
#endif
    if (rtmp > 0.0)
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
using namespace std;

int main() {
    mplapackint m, n;

    string str;
    getline(cin, str);
    cout << str << endl;
    getline(cin, str);
    stringstream ss(str);
    ss >> m;
    ss >> n;
    printf("# m n %d %d \n", (int)m, (int)m);

    _Float128 *a = new _Float128[m * n];
    _Float128 *s = new _Float128[std::min(m, n)];
    _Float128 *u = new _Float128[m * m];
    _Float128 *vt = new _Float128[n * n];
    mplapackint lwork = std::max({(mplapackint)1, 3 * std::min(m, n) + std::max(m, n), 5 * std::min(m, n)});
    _Float128 *work = new _Float128[lwork];
    mplapackint info;
    double dtmp;
    for (int i = 0; i < m; i++) {
        getline(cin, str);
        stringstream ss(str);
        for (int j = 0; j < n; j++) {
            ss >> dtmp;
            a[i + j * m] = dtmp;
        }
    }
    printf("# octave check\n");
    printf("a ="); printmat(m, n, a, m); printf("\n");
    Rgesvd("A", "A", m, n, a, m, s, u, m, vt, n, work, lwork, info);
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
    delete[] work;
    delete[] vt;
    delete[] u;
    delete[] s;
    delete[] a;
}
