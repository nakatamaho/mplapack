// Get condition number of Hibert matrix via Rgecon
// written by Nakata Maho, 2010/8/19.

#include <mpblas_gmp.h>
#include <mplapack_gmp.h>

//Matlab/Octave format
void printmat(int N, int M, mpf_class * A, int LDA)
{
    mpf_class mtmp;

    printf("[ ");
    for (int i = 0; i < N; i++) {
	printf("[ ");
	for (int j = 0; j < M; j++) {
	    mtmp = A[i + j * LDA];
	    gmp_printf("%5.2Fe", mtmp.get_mpf_t());
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

void getAinv(mplapackint n, mplapackint lda, mpf_class * A)
{
    mplapackint info;
    mplapackint lwork;
    /* pivot vector allocation */
    mplapackint *ipiv = new mplapackint[n];
    lwork = -1;
    mpf_class *work = new mpf_class[1];
    /* query work space */
    Rgetri(n, A, lda, ipiv, work, lwork, &info);
    lwork = (int) cast2double(work[0]);
    delete[]work;
    work = new mpf_class[std::max((mplapackint) 1, lwork)];
    /* do inversion */
    Rgetrf(n, n, A, lda, ipiv, &info);
    Rgetri(n, A, lda, ipiv, work, lwork, &info);
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

mpf_class get_estimated_condition_num_psd(mplapackint n, mplapackint lda, mpf_class * A)
{
    mpf_class anorm, cond, rcond;
    mpf_class *work = new mpf_class[std::max((mplapackint) 1, n * 3)];
    mplapackint *iwork = new mplapackint[std::max((mplapackint) 1, n)];
    mplapackint info;

    /* First, calculate norm */
    anorm = Rlange("1", n, n, A, lda, work);
    /* Second, do Cholesky factorization */
    Rpotrf("U", n, A, lda, &info);
    /* Third, calculate estimated condition number */
    Rpocon("U", n, A, lda, anorm, &rcond, work, iwork, &info);

    cond = 1.0 / rcond;

    delete[]work;
    delete[]iwork;
    return cond;
}

mpf_class get_exact_condition_num(const char *norm, mplapackint n, mplapackint lda, mpf_class * A)
{
    mpf_class *Ainv = new mpf_class[n * n];
    mpf_class *work = new mpf_class[std::max((mplapackint) 1, n)];
    mpf_class anorm, ainvnorm, cond;
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
    mpf_class *A = new mpf_class[n * n];
    mpf_class mtmp;
    mpf_class cond_est, cond;

//setting Hilbert matrix in A.
    for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	    mtmp = (i + 1) + (j + 1) - 1;
	    A[i + j * n] = 1.0 / mtmp;
	}
    }
    cond = get_exact_condition_num("1", n, n, A);
    cond_est = get_estimated_condition_num_psd(n, n, A);

    printf("Hilbert matrix of order %d\n", (int) n);
    gmp_printf("Estimated condition number: %10.6Fe\n", cond_est.get_mpf_t());
    gmp_printf("    Exact Condition number: %10.6Fe\n", cond.get_mpf_t());
    delete[]A;
}

int main()
{
//initialization of GMP
    int default_prec = 1024;
    mpf_set_default_prec(default_prec);
    for (int n = 2; n < 30; n++) {
	condition_number_demo(n);
    }
}

