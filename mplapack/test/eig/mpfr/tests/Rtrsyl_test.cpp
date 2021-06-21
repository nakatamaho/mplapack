// Get eigenvalues and eigenvecs
// of Matrix A via Rsyev, using MPFR
// This file is freely usable.
// written by Nakata Maho, 2010/5/17.

#include <mpblas_mpfr.h>
#include <mplapack_mpfr.h>

#include <lapacke.h>

void printnum(mpreal a) { mpfr_printf("%10.8Re", mpfr_ptr(a)); }
void printnum(double a) { printf("%10.8e", a); }

void printnum_short(mpreal a) { mpfr_printf("%10.8Re", mpfr_ptr(a)); }
void printnum_short(double a) { printf("%10.8e", a); }

// Matlab/Octave format
void printmat(int N, int M, mpreal *A, int LDA) {
    mpreal mtmp;
    printf("[ ");
    for (int i = 0; i < N; i++) {
        printf("[ ");
        for (int j = 0; j < M; j++) {
            mtmp = A[i + j * LDA];
            mpfr_printf("%3.16Re", mpfr_ptr(mtmp));
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

void printmat(int N, int M, double *A, int LDA) {
    double mtmp;
    printf("[ ");
    for (int i = 0; i < N; i++) {
        printf("[ ");
        for (int j = 0; j < M; j++) {
            mtmp = A[i + j * LDA];
            printf("%3.16e", mtmp);
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

int main() {
    mplapackint m = 6;
    mplapackint n = 1;
    mpreal *a = new mpreal[m * m];
    mpreal *b = new mpreal[n * n];
    mpreal *c = new mpreal[n * m];
    mpreal *corg = new mpreal[n * m];
    mpreal *x = new mpreal[m * n];
    mpreal *ax = new mpreal[m * n];
    mpreal *xb = new mpreal[m * n];
    mplapackint lda = m;
    mplapackint ldb = n;
    mplapackint ldc = m;
    mplapackint ldx = m;
    mplapackint ldax = n;
    mplapackint ldxb = n;
    mpreal scale;
    mplapackint info;
    mpreal rtmp;

    double *a_d = new double[m * m];
    double *b_d = new double[n * n];
    double *c_d = new double[n * m];
    double *corg_d = new double[n * m];
    double *x_d = new double[m * n];
    double *ax_d = new double[m * n];
    double *xb_d = new double[m * n];
    double scale_d;
    double dtmp;

    // setting A matrix
a[0+0*m]=+4.1021326741398649e+00; a[0+1*m]=+5.0996766117784975e-01; a[0+2*m]=-6.7302159454745936e-01; a[0+3*m]=+5.9082894514562767e-01; a[0+4*m]=+7.8022372922007899e-01; a[0+5*m]=-2.2113379881473688e-01;
a[1+0*m]=+0.0000000000000000e+00; a[1+1*m]=+4.7997488713909542e-01; a[1+2*m]=+5.1640603807363783e-01; a[1+3*m]=+1.4272590149420855e-01; a[1+4*m]=-5.4159446710538156e-02; a[1+5*m]=-9.2743230929985238e-02;
a[2+0*m]=+0.0000000000000000e+00; a[2+1*m]=-2.7943872570386069e-01; a[2+2*m]=+4.7997488713909542e-01; a[2+3*m]=-9.7136361689257267e-02; a[2+4*m]=-3.1940791285376230e-01; a[2+5*m]=+1.9264373531514606e-01;
a[3+0*m]=+0.0000000000000000e+00; a[3+1*m]=+0.0000000000000000e+00; a[3+2*m]=+0.0000000000000000e+00; a[3+3*m]=-1.0335122289616158e-01; a[3+4*m]=+4.9449318918641155e-01; a[3+5*m]=-3.5325744510308614e-02;
a[4+0*m]=+0.0000000000000000e+00; a[4+1*m]=+0.0000000000000000e+00; a[4+2*m]=+0.0000000000000000e+00; a[4+3*m]=-3.7199971231898915e-02; a[4+4*m]=-1.0335122289616158e-01; a[4+5*m]=-4.6906480611652918e-01;
a[5+0*m]=+0.0000000000000000e+00; a[5+1*m]=+0.0000000000000000e+00; a[5+2*m]=+0.0000000000000000e+00; a[5+3*m]=+0.0000000000000000e+00; a[5+4*m]=+0.0000000000000000e+00; a[5+5*m]=-5.8784000131286651e-01;

    b[0] = -5.8784000131286651e-01;

    c[0 + 0 * m] = -8.7150060559905582e-02;
    c[1 + 0 * m] = -8.8900291107698124e-03;
    c[2 + 0 * m] = -1.4901746432198207e-01;
    c[3 + 0 * m] = +1.6869563062999970e-01;
    c[4 + 0 * m] = +4.6123530182483677e-02;
    c[5 + 0 * m] = -3.8913537081260602e-01;

    for (int pp = 0; pp < m * m; pp++) a_d[pp] = a[pp];
    for (int pp = 0; pp < n * n; pp++) b_d[pp] = b[pp];
    for (int pp = 0; pp < n * m; pp++) c_d[pp] = c[pp];
    for (int pp = 0; pp < n * m; pp++) corg[pp] = c[pp];
    for (int pp = 0; pp < n * m; pp++) corg_d[pp] = corg[pp];

    printf("a ="); printmat(m, m, a, lda); printf("\n");
    printf("b ="); printmat(n, n, b, ldb); printf("\n");
    printf("c ="); printmat(m, n, c, ldc); printf("\n");

    Rtrsyl("N", "N", -1, m, n, a, lda, b, ldb, c, ldc, scale, info);
    LAPACKE_dtrsyl(LAPACK_COL_MAJOR, 'N', 'N', -1, (int)m, (int)n, a_d, (int)lda, b_d, (int)ldb, c_d, (int)ldc, &scale_d);

    printf("scale = "); printnum(scale); printf("\n");
    printf("scale_d = "); printnum(scale_d); printf("\n");

    for (int pp = 0; pp < n * m; pp++) x[pp] = c[pp];
    for (int pp = 0; pp < n * m; pp++) x_d[pp] = c_d[pp];
    printf("x ="); printmat(m, n, c, ldc); printf("\n");
    printf("x_d ="); printmat(m, n, c_d, ldc); printf("\n");

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            rtmp = 0.0;
            dtmp = 0.0;
            for (int k = 0; k < m; k++) {
                rtmp = rtmp + a[i + k * lda] * x[k + j * ldx];
                dtmp = dtmp + a_d[i + k * lda] * x_d[k + j * ldx];
            }
            ax[i + j * ldax] = rtmp;
            ax_d[i + j * ldax] = dtmp;
        }
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            rtmp = 0.0;
            dtmp = 0.0;
            for (int k = 0; k < n; k++) {
                rtmp = rtmp + x[i + k * ldx] * b[k + j * ldb];
                dtmp = dtmp + x_d[i + k * ldx] * b_d[k + j * ldb];
            }
            xb[i + j * ldxb] = rtmp;
            xb_d[i + j * ldxb] = dtmp;
        }
    }
    printf("ax ="); printmat(m, n, ax, ldax); printf("\n");
    printf("xb ="); printmat(m, n, xb, ldxb); printf("\n");

    printf("ax_d ="); printmat(m, n, ax_d, ldax); printf("\n");
    printf("xb_d ="); printmat(m, n, xb_d, ldxb); printf("\n");

    for (int pp = 0; pp < n * m; pp++) c[pp] = ax[pp] - xb[pp];
    for (int pp = 0; pp < n * m; pp++) c_d[pp] = ax_d[pp] - xb_d[pp];
    printf("ax - bx="); printmat(m, n, c, ldc); printf("\n");
    printf("ax_d - bx_d="); printmat(m, n, c_d, ldc); printf("\n");
    printf("c="); printmat(m, n, corg, ldc); printf("\n");
    printf("c_d="); printmat(m, n, corg_d, ldc); printf("\n");

    delete[] c_d;
    delete[] b_d;
    delete[] a_d;
    delete[] xb;
    delete[] ax;
    delete[] c;
    delete[] b;
    delete[] a;
}
