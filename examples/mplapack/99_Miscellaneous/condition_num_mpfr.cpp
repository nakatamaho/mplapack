// Get condition number of Hibert matrix via Rgecon
// written by Nakata Maho, 2010/8/19.

#include <mpblas_mpfr.h>
#include <mplapack_mpfr.h>

#define MPFR_FORMAT "%+68.64Re"
#define MPFR_SHORT_FORMAT "%+20.16Re"

inline void printnum(mpreal rtmp) { mpfr_printf(MPFR_FORMAT, mpfr_ptr(rtmp)); }
inline void printnum_short(mpreal rtmp) { mpfr_printf(MPFR_SHORT_FORMAT, mpfr_ptr(rtmp)); }

//Matlab/Octave format
void printmat(int N, int M, mpreal * A, int LDA)
{
    mpreal mtmp;

    printf("[ ");
    for (int i = 0; i < N; i++) {
	printf("[ ");
	for (int j = 0; j < M; j++) {
	    mtmp = A[i + j * LDA];
	    mpfr_printf("%5.2Re", mpfr_ptr(mtmp));
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

void getAinv(mplapackint n, mplapackint lda, mpreal * A)
{
    mplapackint info;
    mplapackint lwork;
    /* pivot vector allocation */
    mplapackint *ipiv = new mplapackint[n];
    lwork = -1;
    mpreal *work = new mpreal[1];
    /* query work space */
    Rgetri(n, A, lda, ipiv, work, lwork, info);
    lwork = (int) cast2double(work[0]);
    delete[]work;
    work = new mpreal[std::max((mplapackint) 1, lwork)];
    /* do inversion */
    Rgetrf(n, n, A, lda, ipiv, info);
    Rgetri(n, A, lda, ipiv, work, lwork, info);
    delete[]ipiv;

    if (info == 0)
	return;
    if (info > 0)
	printf("matrix is singular\n");
    exit(1);
    if (info < 0)
	printf("%d th argument had an illegal value\n", (int) info);
    exit(1);
}

mpreal get_estimated_condition_num(const char *norm, mplapackint n, mplapackint lda, mpreal * A)
{
    mpreal anorm, cond, rcond;
    mpreal *work = new mpreal[std::max((mplapackint) 1, n * 4)];
    mplapackint *iwork = new mplapackint[std::max((mplapackint) 1, n)];
    mplapackint info;

    /* First, calculate norm */
    anorm = Rlange(norm, n, n, A, lda, work);
    /* Second, do LU factorization */
    Rgetrf(n, n, A, lda, iwork, info);
    /* Third, calculate estimated condition number */
    Rgecon(norm, n, A, lda, anorm, rcond, work, iwork, info);

    cond = 1.0 / rcond;

    delete[]work;
    delete[]iwork;
    return cond;
}

mpreal get_exact_condition_num(const char *norm, mplapackint n, mplapackint lda, mpreal * A)
{
    mpreal *Ainv = new mpreal[n * n];
    mpreal *work = new mpreal[std::max((mplapackint) 1, n)];
    mpreal anorm, ainvnorm, cond;
//save A matrix
    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    Ainv[i + j * n] = A[i + j * n];
	}
    }
//get inversion of A matrix
    getAinv(n, lda, Ainv);
//get norm of A and Ainv
    anorm = Rlange(norm, n, n, A, lda, work);
    ainvnorm = Rlange(norm, n, n, Ainv, lda, work);
//condition number
    cond = anorm * ainvnorm;

    delete[]Ainv;
    delete[]work;
    return cond;
}

void condition_number_demo(mplapackint n)
{
    mplapackint lwork, info;
    mpreal *A = new mpreal[n * n];
    mpreal mtmp;
    mpreal cond_est, cond;

//setting Hilbert matrix in A.
    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    mtmp = (i + 1) + (j + 1) - 1;
	    A[i + j * n] = 1.0 / mtmp;
	}
    }
    cond = get_exact_condition_num("1", n, n, A);
    cond_est = get_estimated_condition_num("1", n, n, A);

    printf("Hilbert matrix of order %d\n", (int) n);
    mpfr_printf("Estimated condition number: %10.6Re\n", mpfr_ptr(cond_est));
    mpfr_printf("    Exact Condition number: %10.6Re\n", mpfr_ptr(cond));
    delete[]A;
}

int main()
{
//initialization of MPFR
    int default_prec = 512;
    mpfr_set_default_prec(default_prec);

    for (int n = 2; n < 100; n++) {
	condition_number_demo(n);
    }
}
