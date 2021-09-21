//taking from Collection of Matrices for Testing Computational Algorithms 1969 Robert T. Gregory, David L. Karney pp.30
int main()
{
    INTEGER n = 3;
    INTEGER lwork, info;

    COMPLEX *a = new COMPLEX[n * n];
    INTEGER *ipiv = new INTEGER[n];

//setting a matrix


    a[0 + 0 * n] = COMPLEX(1.0, 0.0);   a[0 + 1 * n] = COMPLEX(1.0, 2.0);    a[0 + 2 * n] = COMPLEX(2.0, 10.0);
    a[1 + 0 * n] = COMPLEX(1.0, 1.0);   a[1 + 1 * n] = COMPLEX(0.0, 3.0);    a[1 + 2 * n] = COMPLEX(-5.0, 14.0);
    a[2 + 0 * n] = COMPLEX(1.0, 1.0);   a[2 + 1 * n] = COMPLEX(0.0, 5.0);    a[2 + 2 * n] = COMPLEX(-8.0, 20.0);

    printf("a ="); printmat(n, n, a, n); printf("\n");

//work space query
    lwork = -1;
    COMPLEX *work = new COMPLEX[1];

    Cgetri(n, a, n, ipiv, work, lwork, info);
    lwork = castInTEGER (work[0].real());
    delete[]work;
    work = new COMPLEX[std::max(1, (int) lwork)];

//inverse matrix
    Cgetrf(n, n, a, n, ipiv, info);
    Cgetri(n, a, n, ipiv, work, lwork, info);

    printf("ainv ="); printmat(n, n, a, n); printf("\n");
    delete[]work;
    delete[]ipiv;
    delete[]a;
}
