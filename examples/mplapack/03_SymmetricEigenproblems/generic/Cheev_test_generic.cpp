int main()
{
    INTEGER n = 3;
    INTEGER lwork, info;

    COMPLEX *A = new COMPLEX[n * n];
    REAL *w = new REAL[n];
    REAL *rwork = new REAL[3 * n - 1];

//setting A matrix
    A[0 + 0 * n] = 2.0;               A[0 + 1 * n] = COMPLEX(0.0, -1.0);    A[0 + 2 * n] = 0.0;
    A[1 + 0 * n] = COMPLEX(0.0, 1.0); A[1 + 1 * n] = 2.0;                   A[1 + 2 * n] = 0.0;
    A[2 + 0 * n] = 0.0;               A[2 + 1 * n] = 0.0;                   A[2 + 2 * n] = 3.0;

    printf("A ="); printmat(n, n, A, n); printf("\n");
//work space query
    lwork = -1;
    COMPLEX *work = new COMPLEX[1];

    Cheev("V", "U", n, A, n, w, work, lwork, rwork, info);
    lwork = (int) cast2double (work[0].real());
    delete[]work;
    work = new COMPLEX[std::max((INTEGER) 1, lwork)];
//inverse matrix
    Cheev("V", "U", n, A, n, w, work, lwork, rwork, info);
//print out some results.
    printf("#eigenvalues \n");
    printf("w ="); printmat(n, 1, w, 1); printf("\n");

    printf("#eigenvecs \n");
    printf("U ="); printmat(n, n, A, n); printf("\n");
    printf("#you can check eigenvalues using octave/Matlab by:\n");
    printf("eig(A)\n");
    printf("#you can check eigenvectors using octave/Matlab by:\n");
    printf("U'*A*U\n");

    delete[]work;
    delete[]w;
    delete[]A;
}
