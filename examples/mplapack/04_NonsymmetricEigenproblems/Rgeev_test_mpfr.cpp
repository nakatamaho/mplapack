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

inline void printnum(mpreal rtmp) { mpfr_printf(MPFR_SHORT_FORMAT, mpfr_ptr(rtmp)); }
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

bool rselect(mpreal ar, mpreal ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    mplapackint n = 4;
    mpreal *a = new mpreal[n * n];
    mpreal *vl = new mpreal[n * n];
    mpreal *vr = new mpreal[n * n];
    mplapackint lwork = 4 * n;
    mpreal *wr = new mpreal[n];
    mpreal *wi = new mpreal[n];
    mpreal *work = new mpreal[lwork];
    mplapackint info;
    // setting A matrix
    a[0 + 0 * n] = 4.0; a[0 + 1 * n] = -5.0; a[0 + 2 * n] =  0.0;  a[0 + 3 * n] =  3.0;
    a[1 + 0 * n] = 0.0; a[1 + 1 * n] =  4.0; a[1 + 2 * n] = -3.0;  a[1 + 3 * n] = -5.0;
    a[2 + 0 * n] = 5.0; a[2 + 1 * n] = -3.0; a[2 + 2 * n] =  4.0;  a[2 + 3 * n] =  0.0;
    a[3 + 0 * n] = 3.0; a[3 + 1 * n] =  0.0; a[3 + 2 * n] =  5.0;  a[3 + 3 * n] =  4.0;

    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a =");
    printmat(n, n, a, n);
    printf("\n");
    Rgeev("V", "V", n, a, n, wr, wi, vl, n, vr, n, work, lwork, info);
    printf("# right vectors\n");
    for (int j = 1; j <= n; j = j + 1) {
        if (abs(wi[j - 1]) < 1e-15) {
	    printf("vr_%d =[ ", j);
            for (int i = 1; i <= n - 1; i = i + 1) {
                printnum(vr[(i - 1) + (j - 1) * n]); printf(", ");
            }
            printnum(vr[(n - 1) + (j - 1) * n]); printf("];\n");
        } else {
	    printf("vr_%d =[ ", j);
            for (int i = 1; i <= n - 1; i = i + 1) {
                printnum(vr[(i - 1) + (j - 1) * n]); printnum(-vr[(i - 1) + j * n]); printf("i, ");
            }
            printnum(vr[(n - 1) + (j - 1) * n]); printnum(-vr[(n - 1) + j * n]); printf("i ];\n");
            printf("vr_%d =[ ", j + 1);
            for (int i = 1; i <= n - 1; i = i + 1) {
                printnum(vr[(i - 1) + (j - 1) * n]); printnum(vr[(i - 1) + j * n]); printf("i, ");
            }
            printnum(vr[(n - 1) + (j - 1) * n]); printnum(vr[(n - 1) + j * n]); printf("i ];\n");
            j++; 
	}
    }
    printf("# left vectors\n");
    for (int j = 1; j <= n; j = j + 1) {
        if (abs(wi[j - 1]) < 1e-15) {
	    printf("vl_%d =[ ", j);
            for (int i = 1; i <= n - 1; i = i + 1) {
                printnum(vl[(i - 1) + (j - 1) * n]); printf(", ");
            }
            printnum(vl[(n - 1) + (j - 1) * n]); printf("];\n");
        } else {
	    printf("vl_%d =[ ", j);
            for (int i = 1; i <= n - 1; i = i + 1) {
                printnum(vl[(i - 1) + (j - 1) * n]); printnum(-vl[(i - 1) + j * n]); printf("i, ");
            }
            printnum(vl[(n - 1) + (j - 1) * n]); printnum(-vl[(n - 1) + j * n]); printf("i ];\n");
            printf("vl_%d =[ ", j + 1);
            for (int i = 1; i <= n - 1; i = i + 1) {
                printnum(vl[(i - 1) + (j - 1) * n]); printnum(vl[(i - 1) + j * n]); printf("i, ");
            }
            printnum(vl[(n - 1) + (j - 1) * n]); printnum(vl[(n - 1) + j * n]); printf("i ];\n");
            j++; 
	}
    }

    for (int i = 1; i <= n; i = i + 1) {
        printf("w_%d = ", (int)i); printnum(wr[i - 1]); printf(" "); printnum(wi[i - 1]); printf("i\n");
    }

    for (int i = 1; i <= n; i = i + 1) {
        printf("disp (\"a * vr_%d\")\n", i);
        printf("a * vr_%d'\n", i);
        printf("disp (\"w_%d * vr_%d\")\n", i, i);
        printf("w_%d * vr_%d'\n", i, i);
        printf("disp (\"vr_%d\")\n", i);
        printf("vr_%d'\n", i);
    }

    for (int i = 1; i <= n; i = i + 1) {
        printf("disp (\"vl_%d * a \")\n", i);
        printf("vl_%d * a\n", i);
        printf("disp (\"w_%d * vl_%d \")\n", i, i);
        printf("w_%d * vl_%d\n", i, i);
        printf("disp (\"vl_%d\")\n", i);
        printf("vl_%d\n", i);
    }

    delete[] work;
    delete[] wr;
    delete[] wi;
    delete[] vr;
    delete[] vl;
    delete[] a;
}
