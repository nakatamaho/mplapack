#include <fem.hpp> // Fortran EMulation library of fable module

namespace placeholder_please_replace {

using namespace major_types;

void xerbla(...) { throw std::runtime_error("Missing function implementation: xerbla"); }

using common;

//> \brief <b> DGTSV computes the solution to system of linear equations A * X = B for GT matrices </b>
//
//  =========== DOCUMENTATION ===========
//
// Online html documentation available at
//            http://www.netlib.org/lapack/explore-html/
//
//> \htmlonly
//> Download DGTSV + dependencies
//> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgtsv.f">
//> [TGZ]</a>
//> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgtsv.f">
//> [ZIP]</a>
//> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgtsv.f">
//> [TXT]</a>
//> \endhtmlonly
//
//  Definition:
//  ===========
//
//       SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
//
//       .. Scalar Arguments ..
//       INTEGER            INFO, LDB, N, NRHS
//       ..
//       .. Array Arguments ..
//       DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * )
//       ..
//
//> \par Purpose:
//  =============
//>
//> \verbatim
//>
//> DGTSV  solves the equation
//>
//>    A*X = B,
//>
//> where A is an n by n tridiagonal matrix, by Gaussian elimination with
//> partial pivoting.
//>
//> Note that the equation  A**T*X = B  may be solved by interchanging the
//> order of the arguments DU and DL.
//> \endverbatim
//
//  Arguments:
//  ==========
//
//> \param[in] N
//> \verbatim
//>          N is INTEGER
//>          The order of the matrix A.  N >= 0.
//> \endverbatim
//>
//> \param[in] NRHS
//> \verbatim
//>          NRHS is INTEGER
//>          The number of right hand sides, i.e., the number of columns
//>          of the matrix B.  NRHS >= 0.
//> \endverbatim
//>
//> \param[in,out] DL
//> \verbatim
//>          On entry, DL must contain the (n-1) sub-diagonal elements of
//>          A.
//>
//>          On exit, DL is overwritten by the (n-2) elements of the
//>          second super-diagonal of the upper triangular matrix U from
//>          the LU factorization of A, in DL(1), ..., DL(n-2).
//> \endverbatim
//>
//> \param[in,out] D
//> \verbatim
//>          On entry, D must contain the diagonal elements of A.
//>
//>          On exit, D is overwritten by the n diagonal elements of U.
//> \endverbatim
//>
//> \param[in,out] DU
//> \verbatim
//>          On entry, DU must contain the (n-1) super-diagonal elements
//>          of A.
//>
//>          On exit, DU is overwritten by the (n-1) elements of the first
//>          super-diagonal of U.
//> \endverbatim
//>
//> \param[in,out] B
//> \verbatim
//>          On entry, the N by NRHS matrix of right hand side matrix B.
//>          On exit, if INFO = 0, the N by NRHS solution matrix X.
//> \endverbatim
//>
//> \param[in] LDB
//> \verbatim
//>          LDB is INTEGER
//> \endverbatim
//>
//> \param[out] INFO
//> \verbatim
//>          INFO is INTEGER
//>          = 0: successful exit
//>          < 0: if INFO = -i, the i-th argument had an illegal value
//>          > 0: if INFO = i, U(i,i) is exactly zero, and the solution
//>               has not been computed.  The factorization has not been
//>               completed unless i = N.
//> \endverbatim
//
//  Authors:
//  ========
//
//> \author Univ. of Tennessee
//> \author Univ. of California Berkeley
//> \author Univ. of Colorado Denver
//> \author NAG Ltd.
//
//> \ingroup REALGTsolve
//
//  =====================================================================
void dgtsv(INTEGER const &n, INTEGER const &nrhs, REAL *dl, REAL *d, REAL *du, REAL *b, INTEGER const &ldb, int &info) {
    INTEGER i = 0;
    const REAL zero = 0.0;
    REAL fact = 0.0;
    REAL temp = 0.0;
    INTEGER j = 0;
    //
    //  -- LAPACK driver routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Executable Statements ..
    //
    info = 0;
    if (n < 0) {
        info = -1;
    } else if (nrhs < 0) {
        info = -2;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -7;
    }
    if (info != 0) {
        xerbla("DGTSV ", -info);
        return;
    }
    //
    if (n == 0) {
        return;
    }
    //
    if (nrhs == 1) {
        for (i = 1; i <= n - 2; i = i + 1) {
            if (abs(d[i - 1]) >= abs(dl[i - 1])) {
                //
                //              No row interchange required
                //
                if (d[i - 1] != zero) {
                    fact = dl[i - 1] / d[i - 1];
                    d[(i + 1) - 1] = d[(i + 1) - 1] - fact * du[i - 1];
                    b[((i + 1) - 1)] = b[((i + 1) - 1)] - fact * b[(i - 1)];
                } else {
                    info = i;
                    return;
                }
                dl[i - 1] = zero;
            } else {
                //
                //              Interchange rows I and I+1
                //
                fact = d[i - 1] / dl[i - 1];
                d[i - 1] = dl[i - 1];
                temp = d[(i + 1) - 1];
                d[(i + 1) - 1] = du[i - 1] - fact * temp;
                dl[i - 1] = du[(i + 1) - 1];
                du[(i + 1) - 1] = -fact * dl[i - 1];
                du[i - 1] = temp;
                temp = b[(i - 1)];
                b[(i - 1)] = b[((i + 1) - 1)];
                b[((i + 1) - 1)] = temp - fact * b[((i + 1) - 1)];
            }
        }
        if (n > 1) {
            i = n - 1;
            if (abs(d[i - 1]) >= abs(dl[i - 1])) {
                if (d[i - 1] != zero) {
                    fact = dl[i - 1] / d[i - 1];
                    d[(i + 1) - 1] = d[(i + 1) - 1] - fact * du[i - 1];
                    b[((i + 1) - 1)] = b[((i + 1) - 1)] - fact * b[(i - 1)];
                } else {
                    info = i;
                    return;
                }
            } else {
                fact = d[i - 1] / dl[i - 1];
                d[i - 1] = dl[i - 1];
                temp = d[(i + 1) - 1];
                d[(i + 1) - 1] = du[i - 1] - fact * temp;
                du[i - 1] = temp;
                temp = b[(i - 1)];
                b[(i - 1)] = b[((i + 1) - 1)];
                b[((i + 1) - 1)] = temp - fact * b[((i + 1) - 1)];
            }
        }
        if (d[n - 1] == zero) {
            info = n;
            return;
        }
    } else {
        for (i = 1; i <= n - 2; i = i + 1) {
            if (abs(d[i - 1]) >= abs(dl[i - 1])) {
                //
                //              No row interchange required
                //
                if (d[i - 1] != zero) {
                    fact = dl[i - 1] / d[i - 1];
                    d[(i + 1) - 1] = d[(i + 1) - 1] - fact * du[i - 1];
                    for (j = 1; j <= nrhs; j = j + 1) {
                        b[((i + 1) - 1) + (j - 1) * ldb] = b[((i + 1) - 1) + (j - 1) * ldb] - fact * b[(i - 1) + (j - 1) * ldb];
                    }
                } else {
                    info = i;
                    return;
                }
                dl[i - 1] = zero;
            } else {
                //
                //              Interchange rows I and I+1
                //
                fact = d[i - 1] / dl[i - 1];
                d[i - 1] = dl[i - 1];
                temp = d[(i + 1) - 1];
                d[(i + 1) - 1] = du[i - 1] - fact * temp;
                dl[i - 1] = du[(i + 1) - 1];
                du[(i + 1) - 1] = -fact * dl[i - 1];
                du[i - 1] = temp;
                for (j = 1; j <= nrhs; j = j + 1) {
                    temp = b[(i - 1) + (j - 1) * ldb];
                    b[(i - 1) + (j - 1) * ldb] = b[((i + 1) - 1) + (j - 1) * ldb];
                    b[((i + 1) - 1) + (j - 1) * ldb] = temp - fact * b[((i + 1) - 1) + (j - 1) * ldb];
                }
            }
        }
        if (n > 1) {
            i = n - 1;
            if (abs(d[i - 1]) >= abs(dl[i - 1])) {
                if (d[i - 1] != zero) {
                    fact = dl[i - 1] / d[i - 1];
                    d[(i + 1) - 1] = d[(i + 1) - 1] - fact * du[i - 1];
                    for (j = 1; j <= nrhs; j = j + 1) {
                        b[((i + 1) - 1) + (j - 1) * ldb] = b[((i + 1) - 1) + (j - 1) * ldb] - fact * b[(i - 1) + (j - 1) * ldb];
                    }
                } else {
                    info = i;
                    return;
                }
            } else {
                fact = d[i - 1] / dl[i - 1];
                d[i - 1] = dl[i - 1];
                temp = d[(i + 1) - 1];
                d[(i + 1) - 1] = du[i - 1] - fact * temp;
                du[i - 1] = temp;
                for (j = 1; j <= nrhs; j = j + 1) {
                    temp = b[(i - 1) + (j - 1) * ldb];
                    b[(i - 1) + (j - 1) * ldb] = b[((i + 1) - 1) + (j - 1) * ldb];
                    b[((i + 1) - 1) + (j - 1) * ldb] = temp - fact * b[((i + 1) - 1) + (j - 1) * ldb];
                }
            }
        }
        if (d[n - 1] == zero) {
            info = n;
            return;
        }
    }
    //
    //     Back solve with the matrix U from the factorization.
    //
    if (nrhs <= 2) {
        j = 1;
    statement_70:
        b[(n - 1) + (j - 1) * ldb] = b[(n - 1) + (j - 1) * ldb] / d[n - 1];
        if (n > 1) {
            b[((n - 1) - 1) + (j - 1) * ldb] = (b[((n - 1) - 1) + (j - 1) * ldb] - du[(n - 1) - 1] * b[(n - 1) + (j - 1) * ldb]) / d[(n - 1) - 1];
        }
        for (i = n - 2; i >= 1; i = i - 1) {
            b[(i - 1) + (j - 1) * ldb] = (b[(i - 1) + (j - 1) * ldb] - du[i - 1] * b[((i + 1) - 1) + (j - 1) * ldb] - dl[i - 1] * b[((i + 2) - 1) + (j - 1) * ldb]) / d[i - 1];
        }
        if (j < nrhs) {
            j++;
            goto statement_70;
        }
    } else {
        for (j = 1; j <= nrhs; j = j + 1) {
            b[(n - 1) + (j - 1) * ldb] = b[(n - 1) + (j - 1) * ldb] / d[n - 1];
            if (n > 1) {
                b[((n - 1) - 1) + (j - 1) * ldb] = (b[((n - 1) - 1) + (j - 1) * ldb] - du[(n - 1) - 1] * b[(n - 1) + (j - 1) * ldb]) / d[(n - 1) - 1];
            }
            for (i = n - 2; i >= 1; i = i - 1) {
                b[(i - 1) + (j - 1) * ldb] = (b[(i - 1) + (j - 1) * ldb] - du[i - 1] * b[((i + 1) - 1) + (j - 1) * ldb] - dl[i - 1] * b[((i + 2) - 1) + (j - 1) * ldb]) / d[i - 1];
            }
        }
    }
    //
    //     End of DGTSV
    //
}

} // namespace placeholder_please_replace
