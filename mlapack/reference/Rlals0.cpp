/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlals0.cpp,v 1.8 2010/08/07 04:48:32 nakatamaho Exp $ 
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

void
Rlals0(INTEGER icompq, INTEGER nl, INTEGER nr, INTEGER sqre, INTEGER nrhs, REAL * B, INTEGER ldb,
       REAL * bx, INTEGER ldbx, INTEGER * perm, INTEGER givptr, INTEGER * givcol, INTEGER ldgcol,
       REAL * givnum, INTEGER ldgnum, REAL * poles, REAL * difl, REAL * difr, REAL * z, INTEGER k, REAL c, REAL s, REAL * work, INTEGER * info)
{
    INTEGER i, j, m, n;
    REAL dj;
    INTEGER nlp1;
    REAL temp;
    REAL diflj, difrj = 0.0, dsigj;
    REAL dsigjp = 0.0;
    REAL One = 1.0, Zero = 0.0;

    *info = 0;
    if (icompq < 0 || icompq > 1) {
	*info = -1;
    } else if (nl < 1) {
	*info = -2;
    } else if (nr < 1) {
	*info = -3;
    } else if (sqre < 0 || sqre > 1) {
	*info = -4;
    }

    n = nl + nr + 1;

    if (nrhs < 1) {
	*info = -5;
    } else if (ldb < n) {
	*info = -7;
    } else if (ldbx < n) {
	*info = -9;
    } else if (givptr < 0) {
	*info = -11;
    } else if (ldgcol < n) {
	*info = -13;
    } else if (ldgnum < n) {
	*info = -15;
    } else if (k < 1) {
	*info = -20;
    }
    if (*info != 0) {
	Mxerbla("Rlals0", -(*info));
	return;
    }

    m = n + sqre;
    nlp1 = nl + 1;

    if (icompq == 0) {
//Apply back orthogonal transformations from the left.
//Step (1L): apply back the Givens rotations performed.
	for (i = 0; i < givptr; i++) {
	    Rrot(nrhs, &B[givcol[i + (ldgcol * 2)] + ldb], ldb, &B[givcol[i + ldgcol] + ldb], ldb, givnum[i + (ldgnum * 2)], givnum[i + ldgnum]);
	}
//Step (2L): permute rows of B.
	Rcopy(nrhs, &B[nlp1 + ldb], ldb, &bx[ldbx + 1], ldbx);
	for (i = 1; i < n; i++) {
	    Rcopy(nrhs, &B[perm[i] + ldb], ldb, &bx[i + ldbx], ldbx);
	}
//Step (3L): apply the inverse of the left singular vector
//matrix to BX.
	if (k == 1) {
	    Rcopy(nrhs, &bx[0], ldbx, &B[0], ldb);
	    if (z[1] < Zero) {
		Rscal(nrhs, -One, &B[0], ldb);
	    }
	} else {
	    for (j = 0; j < k; j++) {
		diflj = difl[j];
		dj = poles[j + ldgnum];
		dsigj = -poles[j + (ldgnum * 2)];
		if (j < k) {
		    difrj = -difr[j + ldgnum];
		    dsigjp = -poles[j + 1 + (ldgnum * 2)];
		}
		if (z[j] == Zero || poles[j + (ldgnum * 2)] == Zero) {
		    work[j] = Zero;
		} else {
		    work[j] = -poles[j + (ldgnum * 2)] * z[j] / diflj / (poles[j + (ldgnum * 2)] + dj);
		}
		for (i = 0; i < j - 1; i++) {
		    if (z[i] == Zero || poles[i + (ldgnum * 2)] == Zero) {
			work[i] = Zero;
		    } else {
			work[i] = poles[i + (ldgnum * 2)] * z[i]
			    / (Rlamc3(poles[i + (ldgnum * 2)], dsigj) - diflj) / (poles[i + (ldgnum * 2)] + dj);
		    }
		}
		for (i = j + 1; i <= k; i++) {
		    if (z[i] == Zero || poles[i + (ldgnum * 2)] == Zero) {
			work[i] = Zero;
		    } else {
			work[i] = poles[i + (ldgnum * 2)] * z[i]
			    / (Rlamc3(poles[i + (ldgnum * 2)], dsigjp) + difrj) / (poles[i + (ldgnum * 2)] + dj);
		    }
		}
		work[1] = -One;
		temp = Rnrm2(k, &work[0], 1);
		Rgemv("T", k, nrhs, One, &bx[0], ldbx, &work[0], 1, Zero, &B[j + ldb], ldb);
		Rlascl("G", 0, 0, temp, One, 1, nrhs, &B[j + ldb], ldb, info);
	    }
	}
//Move the deflated rows of BX to B also.
	if (k < max(m, n)) {
	    Rlacpy("A", n - k, nrhs, &bx[k + 1 + ldbx], ldbx, &B[k + 1 + ldb], ldb);
	}
    } else {
//Apply back the right orthogonal transformations.
//Step (1R): apply back the new right singular vector matrix
//to B.
	if (k == 1) {
	    Rcopy(nrhs, &B[0], ldb, &bx[0], ldbx);
	} else {
	    for (j = 0; j < k; j++) {
		dsigj = poles[j + (ldgnum * 2)];
		if (z[j] == Zero) {
		    work[j] = Zero;
		} else {
		    work[j] = -z[j] / difl[j] / (dsigj + poles[j + ldgnum]) / difr[j + (ldgnum * 2)];
		}
		for (i = 0; i < j - 1; i++) {
		    if (z[j] == Zero) {
			work[i] = Zero;
		    } else {
			work[i] = z[j] / (Rlamc3(dsigj, -poles[i + 1 + (ldgnum * 2)]) - difr[i + ldgnum]) / (dsigj + poles[i + ldgnum]) / difr[i + (ldgnum * 2)];
		    }

		}
		for (i = j + 1; i <= k; i++) {
		    if (z[j] == Zero) {
			work[i] = Zero;
		    } else {
			work[i] = z[j] / (Rlamc3(dsigj, -poles[i + (ldgnum * 2)]) - difl[i]) / (dsigj + poles[i + ldgnum]) / difr[i + (ldgnum * 2)];
		    }
		}
		Rgemv("T", k, nrhs, One, &B[0], ldb, &work[0], 1, Zero, &bx[j + ldbx], ldbx);
	    }
	}
//Step (2R): if SQRE = 1, apply back the rotation that is
//related to the right null space of the subproblem.
	if (sqre == 1) {
	    Rcopy(nrhs, &B[m + ldb], ldb, &bx[m + ldbx], ldbx);
	    Rrot(nrhs, &bx[ldbx + 1], ldbx, &bx[m + ldbx], ldbx, c, s);
	}
	if (k < max(m, n)) {
	    Rlacpy("A", n - k, nrhs, &B[k + 1 + ldb], ldb, &bx[k + 1 + ldbx], ldbx);
	}
//Step (3R): permute rows of B.
	Rcopy(nrhs, &bx[ldbx + 1], ldbx, &B[nlp1 + ldb], ldb);
	if (sqre == 1) {
	    Rcopy(nrhs, &bx[m + ldbx], ldbx, &B[m + ldb], ldb);
	}
	for (i = 1; i < n; i++) {
	    Rcopy(nrhs, &bx[i + ldbx], ldbx, &B[perm[i] + ldb], ldb);
	}
//Step (4R): apply back the Givens rotations performed.
	for (i = givptr; i >= 1; i--) {
	    Rrot(nrhs, &B[givcol[i + (ldgcol * 2)] + ldb], ldb, &B[givcol[i + ldgcol] + ldb], ldb, givnum[i + (ldgnum * 2)], -givnum[i + ldgnum]);
	}
    }
    return;
}
