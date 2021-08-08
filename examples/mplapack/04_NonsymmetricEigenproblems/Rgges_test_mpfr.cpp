//Public domain
#include <mpblas_mpfr.h>
#include <mplapack_mpfr.h>

#include <lapacke.h>

void printnum_short(mpreal a)
{
  mpfr_printf("%10.8Re", mpfr_ptr(a));
}

void printnum(mpreal a) { mpfr_printf("%10.8Re", mpfr_ptr(a)); }

// Matlab/Octave format
void printmat(int N, int M, mpreal *A, int LDA) {
    mpreal mtmp;
    printf("[ ");
    for (int i = 0; i < N; i++) {
        printf("[ ");
        for (int j = 0; j < M; j++) {
            mtmp = A[i + j * LDA];
            mpfr_printf("%8.6Re", mpfr_ptr(mtmp));
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

int main()
{
    mplapackint n = 4;

    mpreal *t = new mpreal[n * n];
    mpreal *q = new mpreal[n * n];    
    mpreal *wi = new mpreal[n];
    mpreal *wr = new mpreal[n];
    mpreal *work = new mpreal[n];
    mplapackint ldt = n;
    mplapackint ldq = n;
    mplapackint lwork = n;
    mplapackint info;
//setting A matrix
    t[0 + 0 * n] =  0.35;   t[0 + 1 * n] =-0.1160;  t[0 + 2 * n] = -0.3886; t[0 + 3 * n] = -0.2942;
    t[1 + 0 * n] = -0.5140; t[1 + 1 * n] = 0.1225;  t[1 + 2 * n] =  0.1004; t[1 + 3 * n] = 0.1126;    
    t[2 + 0 * n] = 0.0;     t[2 + 1 * n] = 0.6443;  t[2 + 2 * n] = -0.1357; t[2 + 3 * n] = -0.0977;    
    t[3 + 0 * n] = 0.0;     t[3 + 1 * n] = 0.0;     t[3 + 2 * n] =  0.4262; t[3 + 3 * n] = 0.1632; 
    
    q[0 + 0 * n] = 1.0;    q[0 + 1 * n] = 0.0;     q[0 + 2 * n] = 0.0;     q[0 + 3 * n] = 0.0;
    q[1 + 0 * n] = 0.0;    q[1 + 1 * n] = 1.0;     q[1 + 2 * n] = 0.0;     q[1 + 3 * n] = 0.0;
    q[2 + 0 * n] = 0.0;    q[2 + 1 * n] = 0.0;     q[2 + 2 * n] = 1.0;     q[2 + 3 * n] = 0.0;
    q[3 + 0 * n] = 0.0;    q[3 + 1 * n] = 0.0;     q[3 + 2 * n] = 0.0;     q[3 + 3 * n] = 1.0;

    printf("t =");  printmat(n, n, t, n); printf("\n");
    printf("q =");  printmat(n, n, q, n); printf("\n");    
    Rgges("S", "V", n, 1, n, t, ldt, wr, wi, q, ldt, work, lwork, info);
    for (int i = 1; i <= n; i = i + 1) {
      printf("wr, wi %d = ", (int)i); printnum(wr[i-1]), printf(" "); printnum(wi[i-1]); printf("\n");
    }
	
    delete[]work;
    delete[]wr;
    delete[]wi;
    delete[]q;
    delete[]t;
}
