/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rggglm.cpp,v 1.4 2010/08/07 04:48:32 nakatamaho Exp $ 
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */
/*
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

    $COPYRIGHT$

    Additional copyrights may follow

    $HEADER$

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer. 

    - Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer listed
    in this license in the documentation and/or other materials
    provided with the distribution.

    - Neither the name of the copyright holders nor the names of its
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
  */

#include <mblas.h>
#include <mlapack.h>

void Rggglm(INTEGER n, INTEGER m, INTEGER p, REAL * A, INTEGER lda, REAL * B, INTEGER ldb, REAL * d, REAL * x, REAL * y, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i, nb, np, nb1, nb2, nb3, nb4, lopt;
    INTEGER lwkmin;
    INTEGER lwkopt;
    INTEGER lquery;
    REAL Zero = 0.0, One = 1.0;

//Test the input parameters
    *info = 0;
    np = min(n, p);
    lquery = lwork == -1;
    if (n < 0) {
	*info = -1;
    } else if (m < 0 || m > n) {
	*info = -2;
    } else if (p < 0 || p < n - m) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -5;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -7;
    }
//Calculate workspace
    if (*info == 0) {
	if (n == 0) {
	    lwkmin = 1;
	    lwkopt = 1;
	} else {
	    nb1 = iMlaenv(1, "Rgeqrf", " ", n, m, -1, -1);
	    nb2 = iMlaenv(1, "Rgerqf", " ", n, m, -1, -1);
	    nb3 = iMlaenv(1, "Rormqr", " ", n, m, p, -1);
	    nb4 = iMlaenv(1, "Rormrq", " ", n, m, p, -1);
	    nb = max(max(max(nb1, nb2), nb3), nb4);
	    lwkmin = m + n + p;
	    lwkopt = m + np + max(n, p) * nb;
	}
	work[1] = lwkopt;
	if (lwork < lwkmin && !lquery) {
	    *info = -12;
	}
    }
    if (*info != 0) {
	Mxerbla("Rggglm", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//Compute the GQR factorization of matrices A and B:
//       Q'*A = ( R11 ) M,    Q'*B*Z' = ( T11   T12 ) M
//              (  0  ) N-M             (  0    T22 ) N-M
//                 M                     M+P-N  N-M

//where R11 and T22 are upper triangular, and Q and Z are
//orthogonal.
    Rggqrf(n, m, p, &A[0], lda, &work[0], &B[0], ldb, &work[m + 1], &work[m + np + 1], lwork - m - np, info);
    lopt = (INTEGER) cast2double(work[m + np + 1]);
//Update left-hand-side vector d = Q'*d = ( d1 ) M
//                                        ( d2 ) N-M
    Rormqr("Left", "Transpose", n, 1, m, &A[0], lda, &work[0], &d[0], max((INTEGER) 1, n), &work[m + np + 1], lwork - m - np, info);
    lopt = max(lopt, (INTEGER) cast2double(work[m + np + 1]));
//Solve T22*y2 = d2 for y2
    if (n > m) {
	Rtrtrs("Upper", "No transpose", "Non unit", n - m, 1, &B[m + 1 + (m + p - n + 1) * ldb], ldb, &d[m + 1], n - m, info);
	if (*info > 0) {
	    *info = 1;
	    return;
	}
	Rcopy(n - m, &d[m + 1], 1, &y[m + p - n + 1], 1);
    }
//Set y1 = 0
    for (i = 0; i < m + p - n; i++) {
	y[i] = Zero;
    }
//Update d1 = d1 - T12*y2
    Rgemv("No transpose", m, n - m, -One, &B[(m + p - n + 1) * ldb + 1], ldb, &y[m + p - n + 1], 1, One, &d[0], 1);
//Solve triangular system: R11*x = d1
    if (m > 0) {
	Rtrtrs("Upper", "No Transpose", "Non unit", m, 1, &A[0], lda, &d[0], m, info);
	if (*info > 0) {
	    *info = 2;
	    return;
	}
//Copy D to X
	Rcopy(m, &d[0], 1, &x[0], 1);
    }
//Backward transformation y = Z'*y
    Rormrq("Left", "Transpose", p, 1, np, &B[max((INTEGER) 1, n - p + 1) + ldb], ldb, &work[m + 1], &y[0], max((INTEGER) 1, n - p + 1), &work[m + np + 1], lwork - m - np, info);
    work[1] = m + np + max(lopt, (INTEGER) cast2double(work[m + np + 1]));
    return;
}
