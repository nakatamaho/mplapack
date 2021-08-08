// Public domain

#include <iostream>
#include <string>
#include <sstream>

#include <mpblas_mpfr.h>
#include <mplapack_mpfr.h>

#define MPFR_FORMAT "%+68.64Re"
#define MPFR_SHORT_FORMAT "%+20.16Re"

inline void printnum(mpreal rtmp) { mpfr_printf(MPFR_FORMAT, mpfr_ptr(rtmp)); }
inline void printnum_short(mpreal rtmp) { mpfr_printf(MPFR_SHORT_FORMAT, mpfr_ptr(rtmp)); }

// Matlab/Octave format
void printmat(int N, int M, mpreal *A, int LDA) {
    mpreal mtmp;
    printf("[ ");
    for (int i = 0; i < N; i++) {
        printf("[ ");
        for (int j = 0; j < M; j++) {
            mtmp = A[i + j * LDA];
            printnum_short(mtmp);
            if (j < M - 1)
                printf(", ");
        }
        if (i < N - 1)
            printf("]; ");
        else
            printf("] ");
    }
    printf("]");
}

void printvec(mpreal *A, int len) {
    mpreal tmp;
    printf("[ ");
    for (int i = 0; i < len; i++) {
        tmp = A[i];
        printnum_short(tmp);
        if (i < len - 1)
            printf(", ");
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

    mpreal *a = new mpreal[m * n];
    mpreal *s = new mpreal[std::min(m, n)];
    mpreal *u = new mpreal[m * m];
    mpreal *vt = new mpreal[n * n];
    mplapackint lwork = std::max({(mplapackint)1, 3 * std::min(m, n) + std::max(m, n), 5 * std::min(m, n)});
    mpreal *work = new mpreal[lwork];
    mplapackint info;
    for (int i = 0; i < m; i++) {
        getline(cin, str);
        stringstream ss(str);
        for (int j = 0; j < n; j++) {
            ss >> a[i + j * m]; // string to mpreal (mpfr)
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
