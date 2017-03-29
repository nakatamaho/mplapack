// Get condition number of Hibert matrix via Rgecon
// written by Nakata Maho, 2010/8/19.

#include <mblas_dd.h>
#include <mlapack_dd.h>
#include <stdio.h>

//Matlab/Octave format
void printmat(int N, int M, dd_real * A, int LDA)
{
    dd_real mtmp;

    printf("[ ");
    for (int i = 0; i < N; i++) {
	printf("[ ");
	for (int j = 0; j < M; j++) {
	    mtmp = A[i + j * LDA];
	    printf("%5.2e", mtmp.x[0]);
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

void getAinv(mpackint n, mpackint lda, dd_real * A)
{
    mpackint info;
    mpackint lwork;
    /* pivot vector allocation */
    mpackint *ipiv = new mpackint[n];
    lwork = -1;
    dd_real *work = new dd_real[1];
    /* query work space */
    Rgetri(n, A, lda, ipiv, work, lwork, &info);
    lwork = (int) cast2double(work[0]);
    delete[]work;
    work = new dd_real[std::max((mpackint) 1, lwork)];
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

dd_real get_estimated_condition_num(const char *norm, mpackint n, mpackint lda, dd_real * A)
{
    dd_real anorm, cond, rcond;
    dd_real *work = new dd_real[std::max((mpackint) 1, n * 4)];
    mpackint *iwork = new mpackint[std::max((mpackint) 1, n)];
    mpackint info;

    /* First, calculate norm */
    anorm = Rlange(norm, n, n, A, lda, work);
    /* Second, do LU factorization */
    Rgetrf(n, n, A, lda, iwork, &info);
    /* Third, calculate estimated condition number */
    Rgecon(norm, n, A, lda, anorm, &rcond, work, iwork, &info);

    cond = 1.0 / rcond;

    delete[]work;
    delete[]iwork;
    return cond;
}

dd_real get_exact_condition_num(const char *norm, mpackint n, mpackint lda, dd_real * A)
{
    dd_real *Ainv = new dd_real[n * n];
    dd_real *work = new dd_real[std::max((mpackint) 1, n)];
    dd_real anorm, ainvnorm, cond;
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

void condition_number_demo(mpackint n)
{
    mpackint lwork, info;
    dd_real *A = new dd_real[n * n];
    dd_real mtmp;
    dd_real cond_est, cond;

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
    printf("Estimated condition number: %10.6e\n", cond_est.x[0]);
    printf("    Exact Condition number: %10.6e\n", cond.x[0]);
    delete[]A;
}

int main()
{
    for (int n = 2; n < 20; n++) {
	condition_number_demo(n);
    }
}
