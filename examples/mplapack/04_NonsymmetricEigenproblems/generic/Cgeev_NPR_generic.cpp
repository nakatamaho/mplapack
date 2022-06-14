#include <mplapack_utils_%%MPLIB%%.h>

bool rselect(REAL ar, REAL ai) {
    // sorting rule for eigenvalues.
    return false;
}

int main() {
    INTEGER n = 10;
    COMPLEX *a = new COMPLEX[n * n];
    COMPLEX *w = new COMPLEX[n];
    COMPLEX *vl = new COMPLEX[n * n];
    COMPLEX *vr = new COMPLEX[n * n];
    INTEGER lwork = 4 * n;
    COMPLEX *work = new COMPLEX[lwork];    
    REAL *rwork = new REAL[lwork];
    INTEGER info;
    // setting A matrix
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            a [ (i - 1) + (j - 1) * n ] = 0.0;
        }
    }
    //Tridiagonal Toeplitz matrices: properties and novel applications 
    //https://doi.org/10.1002/nla.1811
    //http://www.math.kent.edu/~reichel/publications/toep3.pdf

    COMPLEX sigma = COMPLEX(4.0, 3.0) / REAL(8.0);
    COMPLEX delta = COMPLEX(16.0, -3.0);
    COMPLEX tau   = COMPLEX(0.0, -5.0);

    for (int i = 1; i <= n; i++) {
        a [ (i - 1) + (i - 1) * n ] = delta;
    }

    for (int i = 1; i <= n - 1; i++) {
        a [ (i - 1) + i * n ] = sigma;
        a [ i + (i - 1) * n ] = tau;
    }

    printf("# Tridiagonal Toeplitz matrices: properties and novel applications, https://doi.org/10.1002/nla.1811 http://www.math.kent.edu/~reichel/publications/toep3.pdf\n");
    printf("# octave check\n");
    printf("split_long_rows(0)\n");
    printf("a ="); printmat(n, n, a, n); printf("\n");
    Cgeev("V", "V", n, a, n, w, vl, n, vr, n, work, lwork, rwork, info);
    printf("lambda ="); printvec(w,n); printf("\n");

    COMPLEX _pi = pi(REAL(0.0));
    COMPLEX *lambda = new COMPLEX[n];
    for (int h = 1; h <= n; h++) {
        lambda [h - 1] = delta + COMPLEX(2.0, 0.0) * sqrt (sigma * tau) * cos( (REAL(h) * _pi) / REAL((int)n + 1) );
    }
    printf("lambda_true = "); printvec(lambda, n); printf("\n");
    printf("vr ="); printmat(n,n,vr,n); printf("\n");    

    delete[] lambda;
    delete[] rwork;
    delete[] work;
    delete[] vr;
    delete[] vl;
    delete[] w;
    delete[] a;
}
