/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rgesvd.cpp,v 1.6 2010/08/07 04:48:32 nakatamaho Exp $ 
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
  
- Redistributions in binary form * must reproduce the above copyright
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

#include <mpblas.h>
#include <mplapack.h>

void Rgesvd(const char *jobu, const char *jobvt, INTEGER m, INTEGER n,
	    REAL * A, INTEGER lda, REAL * s, REAL * u, INTEGER ldu, REAL * vt, INTEGER ldvt, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER i, ie = 0, ir, iu, blk, ncu = 0;
    REAL dum, eps;
    INTEGER nru = 0, iscl;
    REAL anrm;
    INTEGER ierr, itau, ncvt = 0, nrvt = 0;
    INTEGER chunk, minmn, wrkbl = 0, itaup, itauq, mnthr = 0, iwork;
    LOGICAL wntua, wntva, wntun, wntuo, wntvn, wntvo, wntus, wntvs;
    INTEGER bdspac = 0;
    REAL bignum;
    INTEGER ldwrkr, minwrk, ldwrku, maxwrk;
    REAL smlnum;
    LOGICAL lquery, wntuas, wntvas;
    REAL Zero = 0.0, One = 1.0;

    char jobu_jobvt[3];
//Test the input arguments
    *info = 0;
    minmn = min(m, n);
    wntua = Mlsame(jobu, "A");
    wntus = Mlsame(jobu, "S");
    wntuas = wntua || wntus;
    wntuo = Mlsame(jobu, "O");
    wntun = Mlsame(jobu, "N");
    wntva = Mlsame(jobvt, "A");
    wntvs = Mlsame(jobvt, "S");
    wntvas = wntva || wntvs;
    wntvo = Mlsame(jobvt, "O");
    wntvn = Mlsame(jobvt, "N");
    lquery = lwork == -1;

    if (!(wntua || wntus || wntuo || wntun)) {
	*info = -1;
    } else if (!(wntva || wntvs || wntvo || wntvn) || (wntvo && wntuo)) {
	*info = -2;
    } else if (m < 0) {
	*info = -3;
    } else if (n < 0) {
	*info = -4;
    } else if (lda < max((INTEGER) 1, m)) {
	*info = -6;
    } else if (ldu < 1 || (wntuas && ldu < m)) {
	*info = -9;
    } else if (ldvt < 1 || (wntva && ldvt < n) || (wntvs && ldvt < minmn)) {
	*info = -11;
    }
//Compute workspace
// (Note: Comments in the code beginning "Workspace:" describe the
//  minimal amount of workspace needed at that point in the code,
//  as well as the preferred amount for good performance.
//  NB refers to the optimal block size for the immediately
//  following subroutine, as returned by ILAENV.)
    if (*info == 0) {
	minwrk = 0;
	maxwrk = 0;
	if (m >= n && minmn > 0) {
//Compute space needed for DBDSQR
	    jobu_jobvt[0] = jobu[0];
	    jobu_jobvt[1] = jobvt[0];
	    jobu_jobvt[2] = '\0';
	    mnthr = iMlaenv(6, "RGESVD", jobu_jobvt, m, n, 0, 0);
	    bdspac = n * 5;
	    if (m >= mnthr) {
		if (wntun) {
//Path 1 (M much larger than N, JOBU='N')
		    maxwrk = n + n * iMlaenv(1, "RGEQRF", " ", m, n, -1, -1);
		    maxwrk = max(maxwrk, n * 3 + (n << 1) * iMlaenv(1, "RGEBRD", " ", n, n, -1, -1));
		    if (wntvo || wntvas) {
			maxwrk = max(maxwrk, n * 3 + (n - 1) * iMlaenv(1, "RORGBR", "P", n, n, n, -1));
		    }
		    maxwrk = max(maxwrk, bdspac);
		    minwrk = max(n << 2, bdspac);
		} else if (wntuo && wntvn) {
//Path 2 (M much larger than N, JOBU='O', JOBVT='N')
		    wrkbl = n + n * iMlaenv(1, "RGEQRF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, n + n * iMlaenv(1, "RORGQR", " ", m, n, n, -1));
		    wrkbl = max(wrkbl, n * 3 + (n << 1) * iMlaenv(1, "RGEBRD", " ", n, n, -1, -1));
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORGBR", "Q", n, n, n, -1));
		    wrkbl = max(wrkbl, bdspac);
		    maxwrk = max(n * n + wrkbl, n * n + m * n + n);
		    minwrk = max(n * 3 + m, bdspac);
		} else if (wntuo && wntvas) {
//Path 3 (M much larger than * n, JOBU='O', JOBVT='S' or
//'A')
		    wrkbl = n + n * iMlaenv(1, "RGEQRF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, n + n * iMlaenv(1, "RORGQR", " ", m, n, n, -1));
		    wrkbl = max(wrkbl, n * 3 + (n << 1) * iMlaenv(1, "RGEBRD", " ", n, n, -1, -1));
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORGBR", "Q", n, n, n, -1));
		    wrkbl = max(wrkbl, n * 3 + (n - 1) * iMlaenv(1, "RORGBR", "P", n, n, n, -1));
		    wrkbl = max(wrkbl, bdspac);
		    maxwrk = max(n * n + wrkbl, n * n + m * n + n);
		    minwrk = max(n * 3 + m, bdspac);
		} else if (wntus && wntvn) {
//Path 4 (M much larger than * n, JOBU='S', JOBVT='N')
		    wrkbl = n + n * iMlaenv(1, "RGEQRF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, n + n * iMlaenv(1, "RORGQR", " ", m, n, n, -1));
		    wrkbl = max(wrkbl, n * 3 + (n << 1) * iMlaenv(1, "RGEBRD", " ", n, n, -1, -1));
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORGBR", "Q", n, n, n, -1));
		    wrkbl = max(wrkbl, bdspac);
		    maxwrk = n * n + wrkbl;
		    minwrk = max(n * 3 + m, bdspac);
		} else if (wntus && wntvo) {
//Path 5 (M much larger than * n, JOBU='S', JOBVT='O')
		    wrkbl = n + n * iMlaenv(1, "RGEQRF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, n + n * iMlaenv(1, "RORGQR", " ", m, n, n, -1));
		    wrkbl = max(wrkbl, n * 3 + (n << 1) * iMlaenv(1, "RGEBRD", " ", n, n, -1, -1));
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORGBR", "Q", n, n, n, -1));
		    wrkbl = max(wrkbl, n * 3 + (n - 1) * iMlaenv(1, "RORGBR", "P", n, n, n, -1));
		    wrkbl = max(wrkbl, bdspac);
		    maxwrk = (n << 1) * n + wrkbl;
		    minwrk = max(n * 3 + m, bdspac);
		} else if (wntus && wntvas) {
//Path 6 (M much larger than * n, JOBU='S', JOBVT='S' or 'A')
		    wrkbl = n + n * iMlaenv(1, "RGEQRF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, n + n * iMlaenv(1, "RORGQR", " ", m, n, n, -1));
		    wrkbl = max(wrkbl, n * 3 + (n << 1) * iMlaenv(1, "RGEBRD", " ", n, n, -1, -1));
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORGBR", "Q", n, n, n, -1));
		    wrkbl = max(wrkbl, n * 3 + (n - 1) * iMlaenv(1, "RORGBR", "P", n, n, n, -1));
		    wrkbl = max(wrkbl, bdspac);
		    maxwrk = n * n + wrkbl;
		    minwrk = max(n * 3 + m, bdspac);
		} else if (wntua && wntvn) {
//Path 7 (M much larger than * n, JOBU='A', JOBVT='N')
		    wrkbl = n + n * iMlaenv(1, "RGEQRF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, n + m * iMlaenv(1, "RORGQR", " ", m, m, n, -1));
		    wrkbl = max(wrkbl, n * 3 + (n << 1) * iMlaenv(1, "RGEBRD", " ", n, n, -1, -1));
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORGBR", "Q", n, n, n, -1));
		    wrkbl = max(wrkbl, bdspac);
		    maxwrk = n * n + wrkbl;
		    minwrk = max(n * 3 + m, bdspac);
		} else if (wntua && wntvo) {
//Path 8 (M much larger than * n, JOBU='A', JOBVT='O')
		    wrkbl = n + n * iMlaenv(1, "RGEQRF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, n + m * iMlaenv(1, "RORGQR", " ", m, m, n, -1));
		    wrkbl = max(wrkbl, n * 3 + (n << 1) * iMlaenv(1, "RGEBRD", " ", n, n, -1, -1));
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORGBR", "Q", n, n, n, -1));
		    wrkbl = max(wrkbl, n * 3 + (n - 1) * iMlaenv(1, "RORGBR", "P", n, n, n, -1));
		    wrkbl = max(wrkbl, bdspac);
		    maxwrk = (n << 1) * n + wrkbl;
		    minwrk = max(n * 3 + m, bdspac);
		} else if (wntua && wntvas) {
//Path 9 (M much larger than * n, JOBU='A', JOBVT='S' or 'A')
		    wrkbl = n + n * iMlaenv(1, "RGEQRF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, n + m * iMlaenv(1, "RORGQR", " ", m, m, n, -1));
		    wrkbl = max(wrkbl, n * 3 + (n << 1) * iMlaenv(1, "RGEBRD", " ", n, n, -1, -1));
		    wrkbl = max(wrkbl, n * 3 + n * iMlaenv(1, "RORGBR", "Q", n, n, n, -1));
		    wrkbl = max(wrkbl, n * 3 + (n - 1) * iMlaenv(1, "RORGBR", "P", n, n, n, -1));
		    wrkbl = max(wrkbl, bdspac);
		    maxwrk = n * n + wrkbl;
		    minwrk = max(n * 3 + m, bdspac);
		}
	    } else {
//Path 10 (M at least N, but not much larger)
		maxwrk = n * 3 + (m + n) * iMlaenv(1, "RGEBRD", " ", m, n, -1, -1);
		if (wntus || wntuo) {
		    maxwrk = max(maxwrk, n * 3 + n * iMlaenv(1, "RORG" "BR", "Q", m, n, n, -1));
		}
		if (wntua) {
		    maxwrk = max(maxwrk, n * 3 + m * iMlaenv(1, "RORG" "BR", "Q", m, m, n, -1));
		}
		if (!wntvn) {
		    maxwrk = max(maxwrk, n * 3 + (n - 1) * iMlaenv(1, "RORGBR", "P", n, n, n, -1));
		}
		maxwrk = max(maxwrk, bdspac);
		minwrk = max(n * 3 + m, bdspac);
	    }
	} else if (minmn > 0) {
//Compute space needed for DBDSQR
	    jobu_jobvt[0] = jobu[0];
	    jobu_jobvt[1] = jobvt[0];
	    jobu_jobvt[2] = '\0';
	    mnthr = iMlaenv(6, "RGESVD", jobu_jobvt, m, n, 0, 0);
	    bdspac = m * 5;
	    if (n >= mnthr) {
		if (wntvn) {
//Path 1t(N much larger than M, JOBVT='N')
		    maxwrk = m + m * iMlaenv(1, "RGELQF", " ", m, n, -1, -1);
		    maxwrk = max(maxwrk, m * 3 + (m << 1) * iMlaenv(1, "RGEBRD", " ", m, m, -1, -1));
		    if (wntuo || wntuas) {
			maxwrk = max(maxwrk, m * 3 + m * iMlaenv(1, "RORGBR", "Q", m, m, m, -1));
		    }
		    maxwrk = max(maxwrk, bdspac);
		    minwrk = max(m << 2, bdspac);
		} else if (wntvo && wntun) {
//Path 2t(N much larger than M, JOBU='N', JOBVT='O')
		    wrkbl = m + m * iMlaenv(1, "RGELQF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, m + m * iMlaenv(1, "RORGLQ", " ", m, n, m, -1));
		    wrkbl = max(wrkbl, m * 3 + (m << 1) * iMlaenv(1, "RGEBRD", " ", m, m, -1, -1));
		    wrkbl = max(wrkbl, m * 3 + (m - 1) * iMlaenv(1, "RORGBR", "P", m, m, m, -1));
		    wrkbl = max(wrkbl, bdspac);
		    maxwrk = max(m * m + wrkbl, m * m + m * n + m);
		    minwrk = max(m * 3 + n, bdspac);
		} else if (wntvo && wntuas) {
//Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O')
		    wrkbl = m + m * iMlaenv(1, "RGELQF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, m + m * iMlaenv(1, "RORGLQ", " ", m, n, m, -1));
		    wrkbl = max(wrkbl, m * 3 + (m << 1) * iMlaenv(1, "RGEBRD", " ", m, m, -1, -1));
		    wrkbl = max(wrkbl, m * 3 + (m - 1) * iMlaenv(1, "RORGBR", "P", m, m, m, -1));
		    wrkbl = max(wrkbl, m * 3 + m * iMlaenv(1, "RORGBR", "Q", m, m, m, -1));
		    wrkbl = max(wrkbl, bdspac);
		    maxwrk = max(m * m + wrkbl, m * m + m * n + m);
		    minwrk = max(m * 3 + n, bdspac);
		} else if (wntvs && wntun) {
//Path 4t(N much larger than M, JOBU='N', JOBVT='S')
		    wrkbl = m + m * iMlaenv(1, "RGELQF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, m + m * iMlaenv(1, "RORGLQ", " ", m, n, m, -1));
		    wrkbl = max(wrkbl, m * 3 + (m << 1) * iMlaenv(1, "RGEBRD", " ", m, m, -1, -1));
		    wrkbl = max(wrkbl, m * 3 + (m - 1) * iMlaenv(1, "RORGBR", "P", m, m, m, -1));
		    wrkbl = max(wrkbl, bdspac);
		    maxwrk = m * m + wrkbl;
		    minwrk = max(m * 3 + n, bdspac);
		} else if (wntvs && wntuo) {
//Path 5t(N much larger than M, JOBU='O', JOBVT='S')
		    wrkbl = m + m * iMlaenv(1, "RGELQF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, m + m * iMlaenv(1, "RORGLQ", " ", m, n, m, -1));
		    wrkbl = max(wrkbl, m * 3 + (m << 1) * iMlaenv(1, "RGEBRD", " ", m, m, -1, -1));
		    wrkbl = max(wrkbl, m * 3 + (m - 1) * iMlaenv(1, "RORGBR", "P", m, m, m, -1));
		    wrkbl = max(wrkbl, m * 3 + m * iMlaenv(1, "RORGBR", "Q", m, m, m, -1));
		    wrkbl = max(wrkbl, bdspac);
		    maxwrk = (m << 1) * m + wrkbl;
		    minwrk = max(m * 3 + n, bdspac);
		} else if (wntvs && wntuas) {
//Path 6t(N much larger than M, JOBU='S' or 'A', JOBVT='S')
		    wrkbl = m + m * iMlaenv(1, "RGELQF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, m + m * iMlaenv(1, "RORGLQ", " ", m, n, m, -1));
		    wrkbl = max(wrkbl, m * 3 + (m << 1) * iMlaenv(1, "RGEBRD", " ", m, m, -1, -1));
		    wrkbl = max(wrkbl, m * 3 + (m - 1) * iMlaenv(1, "RORGBR", "P", m, m, m, -1));
		    wrkbl = max(wrkbl, m * 3 + m * iMlaenv(1, "RORGBR", "Q", m, m, m, -1));
		    wrkbl = max(wrkbl, bdspac);
		    maxwrk = m * m + wrkbl;
		    minwrk = max(m * 3 + n, bdspac);
		} else if (wntva && wntun) {
//Path 7t(N much larger than M, JOBU='N', JOBVT='A')
		    wrkbl = m + m * iMlaenv(1, "RGELQF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, m + n * iMlaenv(1, "RORGLQ", " ", n, n, m, -1));
		    wrkbl = max(wrkbl, m * 3 + (m << 1) * iMlaenv(1, "RGEBRD", " ", m, m, -1, -1));
		    wrkbl = max(wrkbl, m * 3 + (m - 1) * iMlaenv(1, "RORGBR", "P", m, m, m, -1));
		    wrkbl = max(wrkbl, bdspac);
		    maxwrk = m * m + wrkbl;
		    minwrk = max(m * 3 + n, bdspac);
		} else if (wntva && wntuo) {
//Path 8t(N much larger than M, JOBU='O', JOBVT='A')
		    wrkbl = m + m * iMlaenv(1, "RGELQF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, m + n * iMlaenv(1, "RORGLQ", " ", n, n, m, -1));
		    wrkbl = max(wrkbl, m * 3 + (m << 1) * iMlaenv(1, "RGEBRD", " ", m, m, -1, -1));
		    wrkbl = max(wrkbl, m * 3 + (m - 1) * iMlaenv(1, "RORGBR", "P", m, m, m, -1));
		    wrkbl = max(wrkbl, m * 3 + m * iMlaenv(1, "RORGBR", "Q", m, m, m, -1));
		    wrkbl = max(wrkbl, bdspac);
		    maxwrk = (m << 1) * m + wrkbl;
		    minwrk = max(m * 3 + n, bdspac);
		} else if (wntva && wntuas) {
//Path 9t(N much larger than M, JOBU='S' or 'A', JOBVT='A')
		    wrkbl = m + m * iMlaenv(1, "RGELQF", " ", m, n, -1, -1);
		    wrkbl = max(wrkbl, m + n * iMlaenv(1, "RORGLQ", " ", n, n, m, -1));
		    wrkbl = max(wrkbl, m * 3 + (m << 1) * iMlaenv(1, "RGEBRD", " ", m, m, -1, -1));
		    wrkbl = max(wrkbl, m * 3 + (m - 1) * iMlaenv(1, "RORGBR", "P", m, m, m, -1));
		    wrkbl = max(wrkbl, m * 3 + m * iMlaenv(1, "RORGBR", "Q", m, m, m, -1));
		    wrkbl = max(wrkbl, bdspac);
		    maxwrk = m * m + wrkbl;
		    minwrk = max(m * 3 + n, bdspac);
		}
	    } else {
//Path 10t(N greater than M, but not much larger)
		maxwrk = m * 3 + (m + n) * iMlaenv(1, "RGEBRD", " ", m, n, -1, -1);
		if (wntvs || wntvo) {
		    maxwrk = max(maxwrk, m * 3 + m * iMlaenv(1, "RORG" "BR", "P", m, n, m, -1));
		}
		if (wntva) {
		    maxwrk = max(maxwrk, m * 3 + n * iMlaenv(1, "RORG" "BR", "P", n, n, m, -1));
		}
		if (!wntun) {
		    maxwrk = max(maxwrk, m * 3 + (m - 1) * iMlaenv(1, "RORGBR", "Q", m, m, m, -1));
		}
		maxwrk = max(maxwrk, bdspac);
		minwrk = max(m * 3 + n, bdspac);
	    }
	}
	maxwrk = max(maxwrk, minwrk);
	work[1] = maxwrk;
	if (lwork < minwrk && !lquery) {
	    *info = -13;
	}
    }
    if (*info != 0) {
	Mxerbla("Rgesvd", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (m == 0 || n == 0) {
	return;
    }
//Get machine constants
    eps = Rlamch("P");
    smlnum = sqrt(Rlamch("S")) / eps;
    bignum = One / smlnum;
//Scale A if max element outside range [SMLNUM,BIGNUM]
    anrm = Rlange("M", m, n, &A[0], lda, &dum);
    iscl = 0;
    if (anrm > Zero && anrm < smlnum) {
	iscl = 0;
	Rlascl("G", 0, 0, anrm, smlnum, m, n, &A[0], lda, &ierr);
    } else if (anrm > bignum) {
	iscl = 0;
	Rlascl("G", 0, 0, anrm, bignum, m, n, &A[0], lda, &ierr);
    }
    if (m >= n) {
//A has at least as many rows as columns. If A has sufficiently
//more rows than columns, first reduce using the QR
//decomposition (if sufficient workspace available)
	if (m >= mnthr) {
	    if (wntun) {
//Path 1 (m * much larger than * n, JOBU='N')
//No left singular vectors to be computed
		itau = 1;
		iwork = itau + n;
//Compute A=Q*R
//(Workspace: need 2*N, prefer N+N*NB)
		Rgeqrf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Zero out below R
		Rlaset("L", n - 1, n - 1, Zero, Zero, &A[lda + 2], lda);
		ie = 1;
		itauq = ie + n;
		itaup = itauq + n;
		iwork = itaup + n;
//Bidiagonalize R in A
//(Workspace: need 4*N, prefer 3*N+2*N*NB)
		Rgebrd(n, n, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
		ncvt = 0;
		if (wntvo || wntvas) {
//If right singular vectors desired, generate P'.
//(Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
		    Rorgbr("P", n, n, n, &A[0], lda, &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
		    ncvt = n;
		}
		iwork = ie + n;
//Perform bidiagonal QR iteration, computing right
//singular vectors of A in A if desired
//(Workspace: need BDSPAC)
		Rbdsqr("U", n, ncvt, 0, 0, &s[1], &work[ie], &A[0], lda, &dum, 1, &dum, 1, &work[iwork], info);
//If right singular vectors desired in VT, copy them there
		if (wntvas) {
		    Rlacpy("F", n, n, &A[0], lda, &vt[0], ldvt);
		}
	    } else if (wntuo && wntvn) {
//Path 2 (m * much larger than * n, JOBU='O', JOBVT='N')
//N left singular vectors to be overwritten on A and
//no right singular vectors to be computed
		if (lwork >= n * n + max(n << 2, bdspac)) {
//Sufficient workspace for a fast algorithm
		    ir = 1;
		    if (lwork >= max(wrkbl, lda * n + n) + lda * n) {
//WORK(IU) is LDA by N, WORK(IR) is LDA by N
			ldwrku = lda;
			ldwrkr = lda;
		    } else {
			if (lwork >= max(wrkbl, lda * n + n) + n * n) {
//WORK(IU) is LDA by N, WORK(IR) is N by N
			    ldwrku = lda;
			    ldwrkr = n;
			} else {
//WORK(IU) is LDWRKU by N, WORK(IR) is N by N
			    ldwrku = (lwork - n * n - n) / n;
			    ldwrkr = n;
			}
		    }
		    itau = ir + ldwrkr * n;
		    iwork = itau + n;
//Compute A=Q*R
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
		    Rgeqrf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Copy R to WORK(IR) and zero out below it
		    Rlacpy("U", n, n, &A[0], lda, &work[ir], ldwrkr);
		    Rlaset("L", n - 1, n - 1, Zero, Zero, &work[ir + 1], ldwrkr);
//Generate Q in A
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
		    Rorgqr(m, n, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
		    ie = itau;
		    itauq = ie + n;
		    itaup = itauq + n;
		    iwork = itaup + n;
//Bidiagonalize R in WORK(IR)
//(Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
		    Rgebrd(n, n, &work[ir], ldwrkr, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Generate left vectors bidiagonalizing R
//(Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)
		    Rorgbr("Q", n, n, n, &work[ir], ldwrkr, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
		    iwork = ie + n;
//Perform bidiagonal QR iteration, computing left
//singular vectors of R in WORK(IR)
//(Workspace: need N*N+BDSPAC)
		    Rbdsqr("U", n, 0, n, 0, &s[1], &work[ie], &dum, 1, &work[ir], ldwrkr, &dum, 1, &work[iwork], info);
		    iu = ie + n;
//Multiply Q in A by left singular vectors of R in
//WORK(IR), storing result in WORK(IU) and copying to A
//(Workspace: need N*N+2*N, prefer N*N+M*N+N)
		    for (i = 1; i <= m; i += ldwrku) {
			chunk = min(m - i + 1, ldwrku);
			Rgemm("N", "N", chunk, n, n, One, &A[i + lda], lda, &work[ir], ldwrkr, Zero, &work[iu], ldwrku);
			Rlacpy("F", chunk, n, &work[iu], ldwrku, &A[i + lda], lda);
		    }
		} else {
//Insufficient workspace for a fast algorithm
		    ie = 1;
		    itauq = ie + n;
		    itaup = itauq + n;
		    iwork = itaup + n;
//Bidiagonalize A
//(Workspace: need 3*N+M, prefer 3*N+(M+N)*NB)
		    Rgebrd(m, n, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Generate left vectors bidiagonalizing A
//(Workspace: need 4*N, prefer 3*N+N*NB)
		    Rorgbr("Q", m, n, n, &A[0], lda, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
		    iwork = ie + n;
//Perform bidiagonal QR iteration, computing left
//singular vectors of A in A
//(Workspace: need BDSPAC)
		    Rbdsqr("U", n, 0, m, 0, &s[1], &work[ie], &dum, 1, &A[0], lda, &dum, 1, &work[iwork], info);
		}
	    } else if (wntuo && wntvas) {
//Path 3 (m * much larger than * n, JOBU='O', JOBVT='S' or 'A')
//N left singular vectors to be overwritten on A and
//N right singular vectors to be computed in VT
		if (lwork >= n * n + max(n << 2, bdspac)) {
//Sufficient workspace for a fast algorithm
		    ir = 1;
		    if (lwork >= max(wrkbl, lda * n + n) + lda * n) {
//WORK(IU) is LDA by N and WORK(IR) is LDA by N
			ldwrku = lda;
			ldwrkr = lda;
		    } else {
			if (lwork >= max(wrkbl, lda * n + n) + n * n) {
//WORK(IU) is LDA by N and WORK(IR) is N by N
			    ldwrku = lda;
			    ldwrkr = n;
			} else {
//WORK(IU) is LDWRKU by N and WORK(IR) is N by N
			    ldwrku = (lwork - n * n - n) / n;
			    ldwrkr = n;
			}
		    }
		    itau = ir + ldwrkr * n;
		    iwork = itau + n;
//Compute A=Q*R
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
		    Rgeqrf(m, n, &A[0], lda, &work[itau], &work[iwork]
			   , lwork - iwork + 1, &ierr);
//Copy R to VT, zeroing out below it
		    Rlacpy("U", n, n, &A[0], lda, &vt[0], ldvt);
		    if (n > 1) {
			Rlaset("L", n - 1, n - 1, Zero, Zero, &vt[ldvt + 2], ldvt);
		    }
//Generate Q in A
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
		    Rorgqr(m, n, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
		    ie = itau;
		    itauq = ie + n;
		    itaup = itauq + n;
		    iwork = itaup + n;
//Bidiagonalize R in VT, copying result to WORK(IR)
//(Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
		    Rgebrd(n, n, &vt[0], ldvt, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
		    Rlacpy("L", n, n, &vt[0], ldvt, &work[ir], ldwrkr);
//Generate left vectors bidiagonalizing R in WORK(IR)
//(Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)
		    Rorgbr("Q", n, n, n, &work[ir], ldwrkr, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
//Generate right vectors bidiagonalizing R in VT
//(Workspace: need N*N+4*N-1, prefer N*N+3*N+(N-1)*NB)
		    Rorgbr("P", n, n, n, &vt[0], ldvt, &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
		    iwork = ie + n;
//Perform bidiagonal QR iteration, computing left
//singular vectors of R in WORK(IR) and computing right
//singular vectors of R in VT
//(Workspace: need N*N+BDSPAC)
		    Rbdsqr("U", n, n, n, 0, &s[1], &work[ie], &vt[0], ldvt, &work[ir], ldwrkr, &dum, 1, &work[iwork], info);
		    iu = ie + n;
//Multiply Q in A by left singular vectors of R in
//WORK(IR), storing result in WORK(IU) and copying to A
//(Workspace: need N*N+2*N, prefer N*N+M*N+N)
		    for (i = 1; i <= m; i += ldwrku) {
			chunk = min(m - i + 1, ldwrku);
			Rgemm("N", "N", chunk, n, n, One, &A[i + lda], lda, &work[ir], ldwrkr, Zero, &work[iu], ldwrku);
			Rlacpy("F", chunk, n, &work[iu], ldwrku, &A[i + lda], lda);
		    }
		} else {
//Insufficient workspace for a fast algorithm
		    itau = 1;
		    iwork = itau + n;
//Compute A=Q*R
//(Workspace: need 2*N, prefer N+N*NB)
		    Rgeqrf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Copy R to VT, zeroing out below it
		    Rlacpy("U", n, n, &A[0], lda, &vt[0], ldvt);
		    if (n > 1) {
			Rlaset("L", n - 1, n - 1, Zero, Zero, &vt[ldvt + 2], ldvt);
		    }
//Generate Q in A
//(Workspace: need 2*N, prefer N+N*NB)
		    Rorgqr(m, n, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
		    ie = itau;
		    itauq = ie + n;
		    itaup = itauq + n;
		    iwork = itaup + n;
//Bidiagonalize R in VT
//(Workspace: need 4*N, prefer 3*N+2*N*NB)
		    Rgebrd(n, n, &vt[0], ldvt, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Multiply Q in A by left vectors bidiagonalizing R
//(Workspace: need 3*N+M, prefer 3*N+M*NB)
		    Rormbr("Q", "R", "N", m, n, n, &vt[0], ldvt, &work[itauq], &A[0], lda, &work[iwork], lwork - iwork + 1, &ierr);
//Generate right vectors bidiagonalizing R in VT
//(Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
		    Rorgbr("P", n, n, n, &vt[0], ldvt, &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
		    iwork = ie + n;
//Perform bidiagonal QR iteration, computing left
//singular vectors of A in A and computing right
//singular vectors of A in VT
//(Workspace: need BDSPAC)
		    Rbdsqr("U", n, n, m, 0, &s[1], &work[ie], &vt[0], ldvt, &A[0], lda, &dum, 1, &work[iwork], info);
		}
	    } else if (wntus) {
		if (wntvn) {
//Path 4 (m * much larger than * n, JOBU='S', JOBVT='N')
//N left singular vectors to be computed in U and
//no right singular vectors to be computed
		    if (lwork >= n * n + max(n << 2, bdspac)) {
//Sufficient workspace for a fast algorithm
			ir = 1;
			if (lwork >= wrkbl + lda * n) {
//WORK(IR) is LDA by N
			    ldwrkr = lda;
			} else {
//WORK(IR) is N by N
			    ldwrkr = n;
			}
			itau = ir + ldwrkr * n;
			iwork = itau + n;
//Compute A=Q*R
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
			Rgeqrf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Copy R to WORK(IR), zeroing out below it
			Rlacpy("U", n, n, &A[0], lda, &work[ir], ldwrkr);
			Rlaset("L", n - 1, n - 1, Zero, Zero, &work[ir + 1], ldwrkr);
//Generate Q in A
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
			Rorgqr(m, n, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			ie = itau;
			itauq = ie + n;
			itaup = itauq + n;
			iwork = itaup + n;
//Bidiagonalize R in WORK(IR)
//(Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
			Rgebrd(n, n, &work[ir], ldwrkr, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Generate left vectors bidiagonalizing R in WORK(IR)
//(Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)
			Rorgbr("Q", n, n, n, &work[ir], ldwrkr, &work[itauq]
			       , &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + n;
//Perform bidiagonal QR iteration, computing left
//singular vectors of R in WORK(IR)
//(Workspace: need N*N+BDSPAC)
			Rbdsqr("U", n, 0, n, 0, &s[1], &work[ie], &dum, 1, &work[ir], ldwrkr, &dum, 1, &work[iwork], info);
//Multiply Q in A by left singular vectors of R in
//WORK(IR), storing result in U
//(Workspace: need N*N)
			Rgemm("N", "N", m, n, n, One, &A[0], lda, &work[ir], ldwrkr, Zero, &u[0], ldu);
		    } else {
//Insufficient workspace for a fast algorithm
			itau = 1;
			iwork = itau + n;
//Compute A=Q*R, copying result to U
//(Workspace: need 2*N, prefer N+N*NB)
			Rgeqrf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("L", m, n, &A[0], lda, &u[0], ldu);
//Generate Q in U
//(Workspace: need 2*N, prefer N+N*NB)
			Rorgqr(m, n, n, &u[0], ldu, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			ie = itau;
			itauq = ie + n;
			itaup = itauq + n;
			iwork = itaup + n;
//Zero out below R in A
			Rlaset("L", n - 1, n - 1, Zero, Zero, &A[lda + 2], lda);
//Bidiagonalize R in A
//(Workspace: need 4*N, prefer 3*N+2*N*NB)
			Rgebrd(n, n, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Multiply Q in U by left vectors bidiagonalizing R
//(Workspace: need 3*N+M, prefer 3*N+M*NB)
			Rormbr("Q", "R", "N", m, n, n, &A[0], lda, &work[itauq], &u[0], ldu, &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + n;
//Perform bidiagonal QR iteration, computing left
//singular vectors of A in U
//(Workspace: need BDSPAC)
			Rbdsqr("U", n, 0, m, 0, &s[1], &work[ie], &dum, 1, &u[0], ldu, &dum, 1, &work[iwork], info);
		    }
		} else if (wntvo) {
//Path 5 (m * much larger than * n, JOBU='S', JOBVT='O')
//N left singular vectors to be computed in U and
//N right singular vectors to be overwritten on A
		    if (lwork >= (n << 1) * n + max(n << 2, bdspac)) {
//Sufficient workspace for a fast algorithm
			iu = 1;
			if (lwork >= wrkbl + (lda << 1) * n) {
//WORK(IU) is LDA by N and WORK(IR) is LDA by N
			    ldwrku = lda;
			    ir = iu + ldwrku * n;
			    ldwrkr = lda;
			} else if (lwork >= wrkbl + (lda + n) * n) {
//WORK(IU) is LDA by N and WORK(IR) is N by N
			    ldwrku = lda;
			    ir = iu + ldwrku * n;
			    ldwrkr = n;
			} else {
//WORK(IU) is N by N and WORK(IR) is N by N
			    ldwrku = n;
			    ir = iu + ldwrku * n;
			    ldwrkr = n;
			}
			itau = ir + ldwrkr * n;
			iwork = itau + n;

//Compute A=Q*R
//(Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
			Rgeqrf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Copy R to WORK(IU), zeroing out below it
			Rlacpy("U", n, n, &A[0], lda, &work[iu], ldwrku);
			Rlaset("L", n - 1, n - 1, Zero, Zero, &work[iu + 1], ldwrku);
//Generate Q in A
//(Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
			Rorgqr(m, n, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			ie = itau;
			itauq = ie + n;
			itaup = itauq + n;
			iwork = itaup + n;
//Bidiagonalize R in WORK(IU), copying result to
//WORK(IR)
//(Workspace: need 2*N*N+4*N,
//            prefer 2*N*N+3*N+2*N*NB)
			Rgebrd(n, n, &work[iu], ldwrku, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("U", n, n, &work[iu], ldwrku, &work[ir], ldwrkr);
//Generate left bidiagonalizing vectors in WORK(IU)
//(Workspace: need 2*N*N+4*N, prefer 2*N*N+3*N+N*NB)
			Rorgbr("Q", n, n, n, &work[iu], ldwrku, &work[itauq]
			       , &work[iwork], lwork - iwork + 1, &ierr);
//Generate right bidiagonalizing vectors in WORK(IR)
//(Workspace: need 2*N*N+4*N-1,
//            prefer 2*N*N+3*N+(N-1)*NB)
			Rorgbr("P", n, n, n, &work[ir], ldwrkr, &work[itaup]
			       , &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + n;
//Perform bidiagonal QR iteration, computing left
//singular vectors of R in WORK(IU) and computing
//right singular vectors of R in WORK(IR)
//(Workspace: need 2*N*N+BDSPAC)
			Rbdsqr("U", n, n, n, 0, &s[1], &work[ie], &work[ir], ldwrkr, &work[iu], ldwrku, &dum, 1, &work[iwork], info);
//Multiply Q in A by left singular vectors of R in
//WORK(IU), storing result in U
//(Workspace: need N*N)
			Rgemm("N", "N", m, n, n, One, &A[0], lda, &work[iu], ldwrku, Zero, &u[0], ldu);
//Copy right singular vectors of R to A
//(Workspace: need N*N)
			Rlacpy("F", n, n, &work[ir], ldwrkr, &A[0], lda);
		    } else {
//Insufficient workspace for a fast algorithm
			itau = 1;
			iwork = itau + n;
//Compute A=Q*R, copying result to U
//(Workspace: need 2*N, prefer N+N*NB)
			Rgeqrf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("L", m, n, &A[0], lda, &u[0], ldu);
//Generate Q in U
//(Workspace: need 2*N, prefer N+N*NB)
			Rorgqr(m, n, n, &u[0], ldu, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			ie = itau;
			itauq = ie + n;
			itaup = itauq + n;
			iwork = itaup + n;
//Zero out below R in A
			Rlaset("L", n - 1, n - 1, Zero, Zero, &A[lda + 2], lda);
//Bidiagonalize R in A
//(Workspace: need 4*N, prefer 3*N+2*N*NB)
			Rgebrd(n, n, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Multiply Q in U by left vectors bidiagonalizing R
//(Workspace: need 3*N+M, prefer 3*N+M*NB)
			Rormbr("Q", "R", "N", m, n, n, &A[0], lda, &work[itauq], &u[0], ldu, &work[iwork], lwork - iwork + 1, &ierr);
//Generate right vectors bidiagonalizing R in A
//(Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
			Rorgbr("P", n, n, n, &A[0], lda, &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + n;
//Perform bidiagonal QR iteration, computing left
//singular vectors of A in U and computing right
//singular vectors of A in A
//(Workspace: need BDSPAC)
			Rbdsqr("U", n, n, m, 0, &s[1], &work[ie], &A[0], lda, &u[0], ldu, &dum, 1, &work[iwork], info);
		    }
		} else if (wntvas) {
//Path 6 (m * much larger than * n, JOBU='S', JOBVT='S'
//     or 'A')
//N left singular vectors to be computed in U and
//N right singular vectors to be computed in VT
		    if (lwork >= n * n + max(n << 2, bdspac)) {
//Sufficient workspace for a fast algorithm
			iu = 1;
			if (lwork >= wrkbl + lda * n) {
//WORK(IU) is LDA by N
			    ldwrku = lda;
			} else {
//WORK(IU) is N by N
			    ldwrku = n;
			}
			itau = iu + ldwrku * n;
			iwork = itau + n;
//Compute A=Q*R
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
			Rgeqrf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Copy R to WORK(IU), zeroing out below it
			Rlacpy("U", n, n, &A[0], lda, &work[iu], ldwrku);
			Rlaset("L", n - 1, n - 1, Zero, Zero, &work[iu + 1], ldwrku);
//Generate Q in A
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
			Rorgqr(m, n, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			ie = itau;
			itauq = ie + n;
			itaup = itauq + n;
			iwork = itaup + n;
//Bidiagonalize R in WORK(IU), copying result to VT
//(Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
			Rgebrd(n, n, &work[iu], ldwrku, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("U", n, n, &work[iu], ldwrku, &vt[0], ldvt);
//Generate left bidiagonalizing vectors in WORK(IU)
//(Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)
			Rorgbr("Q", n, n, n, &work[iu], ldwrku, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
//Generate right bidiagonalizing vectors in VT
//(Workspace: need N*N+4*N-1,
//            prefer N*N+3*N+(N-1)*NB)
			Rorgbr("P", n, n, n, &vt[0], ldvt, &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + n;
//Perform bidiagonal QR iteration, computing left
//singular vectors of R in WORK(IU) and computing
//right singular vectors of R in VT
//(Workspace: need N*N+BDSPAC)
			Rbdsqr("U", n, n, n, 0, &s[1], &work[ie], &vt[0], ldvt, &work[iu], ldwrku, &dum, 1, &work[iwork], info);
//Multiply Q in A by left singular vectors of R in
//WORK(IU), storing result in U
//(Workspace: need N*N)
			Rgemm("N", "N", m, n, n, One, &A[0], lda, &work[iu], ldwrku, Zero, &u[0], ldu);
		    } else {
//Insufficient workspace for a fast algorithm
			itau = 1;
			iwork = itau + n;
//Compute A=Q*R, copying result to U
//(Workspace: need 2*N, prefer N+N*NB)
			Rgeqrf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("L", m, n, &A[0], lda, &u[0], ldu);

//Generate Q in U
//(Workspace: need 2*N, prefer N+N*NB)
			Rorgqr(m, n, n, &u[0], ldu, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Copy R to VT, zeroing out below it
			Rlacpy("U", n, n, &A[0], lda, &vt[0], ldvt);
			if (n > 1) {
			    Rlaset("L", n - 1, n - 1, Zero, Zero, &vt[ldvt + 2], ldvt);
			}
			ie = itau;
			itauq = ie + n;
			itaup = itauq + n;
			iwork = itaup + n;
//Bidiagonalize R in VT
//(Workspace: need 4*N, prefer 3*N+2*N*NB)
			Rgebrd(n, n, &vt[0], ldvt, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Multiply Q in U by left bidiagonalizing vectors
//in VT
//(Workspace: need 3*N+M, prefer 3*N+M*NB)
			Rormbr("Q", "R", "N", m, n, n, &vt[0], ldvt, &work[itauq], &u[0], ldu, &work[iwork], lwork - iwork + 1, &ierr);
//Generate right bidiagonalizing vectors in VT
//(Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
			Rorgbr("P", n, n, n, &vt[0], ldvt, &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + n;
//Perform bidiagonal QR iteration, computing left
//singular vectors of A in U and computing right
//singular vectors of A in VT
//(Workspace: need BDSPAC)
			Rbdsqr("U", n, n, m, 0, &s[1], &work[ie], &vt[0], ldvt, &u[0], ldu, &dum, 1, &work[iwork], info);
		    }
		}
	    } else if (wntua) {
		if (wntvn) {
//Path 7 (m * much larger than * n, JOBU='A', JOBVT='N')
//M left singular vectors to be computed in U and
//no right singular vectors to be computed
		    if (lwork >= n * n + max(max(n + m, n << 2), bdspac)) {
//Sufficient workspace for a fast algorithm
			ir = 1;
			if (lwork >= wrkbl + lda * n) {
//WORK(IR) is LDA by N
			    ldwrkr = lda;
			} else {
//WORK(IR) is N by N
			    ldwrkr = n;
			}
			itau = ir + ldwrkr * n;
			iwork = itau + n;
//Compute A=Q*R, copying result to U
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
			Rgeqrf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("L", m, n, &A[0], lda, &u[0], ldu);
//Copy R to WORK(IR), zeroing out below it
			Rlacpy("U", n, n, &A[0], lda, &work[ir], ldwrkr);
			Rlaset("L", n - 1, n - 1, Zero, Zero, &work[ir + 1], ldwrkr);
//Generate Q in U
//(Workspace: need N*N+N+M, prefer N*N+N+M*NB)
			Rorgqr(m, m, n, &u[0], ldu, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			ie = itau;
			itauq = ie + n;
			itaup = itauq + n;
			iwork = itaup + n;
//Bidiagonalize R in WORK(IR)
//(Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
			Rgebrd(n, n, &work[ir], ldwrkr, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Generate left bidiagonalizing vectors in WORK(IR)
//(Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)
			Rorgbr("Q", n, n, n, &work[ir], ldwrkr, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + n;
//Perform bidiagonal QR iteration, computing left
//singular vectors of R in WORK(IR)
//(Workspace: need N*N+BDSPAC)
			Rbdsqr("U", n, 0, n, 0, &s[1], &work[ie], &dum, 1, &work[ir], ldwrkr, &dum, 1, &work[iwork], info);
//Multiply Q in U by left singular vectors of R in
//WORK(IR), storing result in A
//(Workspace: need N*N)
			Rgemm("N", "N", m, n, n, One, &u[0], ldu, &work[ir], ldwrkr, Zero, &A[0], lda);
//Copy left singular vectors of A from A to U
			Rlacpy("F", m, n, &A[0], lda, &u[0], ldu);
		    } else {
//Insufficient workspace for a fast algorithm
			itau = 1;
			iwork = itau + n;
//Compute A=Q*R, copying result to U
//(Workspace: need 2*N, prefer N+N*NB)
			Rgeqrf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("L", m, n, &A[0], lda, &u[0], ldu);
//Generate Q in U
//(Workspace: need N+M, prefer N+M*NB)
			Rorgqr(m, m, n, &u[0], ldu, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			ie = itau;
			itauq = ie + n;
			itaup = itauq + n;
			iwork = itaup + n;
//Zero out below R in A
			Rlaset("L", n - 1, n - 1, Zero, Zero, &A[lda + 2], lda);
//Bidiagonalize R in A
//(Workspace: need 4*N, prefer 3*N+2*N*NB)
			Rgebrd(n, n, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Multiply Q in U by left bidiagonalizing vectors
//in A
//(Workspace: need 3*N+M, prefer 3*N+M*NB)
			Rormbr("Q", "R", "N", m, n, n, &A[0], lda, &work[itauq], &u[0], ldu, &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + n;
//Perform bidiagonal QR iteration, computing left
//singular vectors of A in U
//(Workspace: need BDSPAC)
			Rbdsqr("U", n, 0, m, 0, &s[1], &work[ie], &dum, 1, &u[0], ldu, &dum, 1, &work[iwork], info);
		    }
		} else if (wntvo) {
//Path 8 (m * much larger than * n, JOBU='A', JOBVT='O')
//M left singular vectors to be computed in U and
//N right singular vectors to be overwritten on A
		    if (lwork >= (n << 1) * n + max(max(n + m, n << 2), bdspac)) {
//Sufficient workspace for a fast algorithm
			iu = 1;
			if (lwork >= wrkbl + (lda << 1) * n) {
//WORK(IU) is LDA by N and WORK(IR) is LDA by N
			    ldwrku = lda;
			    ir = iu + ldwrku * n;
			    ldwrkr = lda;
			} else if (lwork >= wrkbl + (lda + n) * n) {
//WORK(IU) is LDA by N and WORK(IR) is N by N
			    ldwrku = lda;
			    ir = iu + ldwrku * n;
			    ldwrkr = n;
			} else {
//WORK(IU) is N by N and WORK(IR) is N by N
			    ldwrku = n;
			    ir = iu + ldwrku * n;
			    ldwrkr = n;
			}
			itau = ir + ldwrkr * n;
			iwork = itau + n;
//Compute A=Q*R, copying result to U
//(Workspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
			Rgeqrf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("L", m, n, &A[0], lda, &u[0], ldu);
//Generate Q in U
//(Workspace: need 2*N*N+N+M, prefer 2*N*N+N+M*NB)
			Rorgqr(m, m, n, &u[0], ldu, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Copy R to WORK(IU), zeroing out below it
			Rlacpy("U", n, n, &A[0], lda, &work[iu], ldwrku);
			Rlaset("L", n - 1, n - 1, Zero, Zero, &work[iu + 1], ldwrku);
			ie = itau;
			itauq = ie + n;
			itaup = itauq + n;
			iwork = itaup + n;
//Bidiagonalize R in WORK(IU), copying result to
//WORK(IR)
//(Workspace: need 2*N*N+4*N,
//            prefer 2*N*N+3*N+2*N*NB)
			Rgebrd(n, n, &work[iu], ldwrku, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("U", n, n, &work[iu], ldwrku, &work[ir], ldwrkr);
//Generate left bidiagonalizing vectors in WORK(IU)
//(Workspace: need 2*N*N+4*N, prefer 2*N*N+3*N+N*NB)
			Rorgbr("Q", n, n, n, &work[iu], ldwrku, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
//Generate right bidiagonalizing vectors in WORK(IR)
//(Workspace: need 2*N*N+4*N-1,
//            prefer 2*N*N+3*N+(N-1)*NB)
			Rorgbr("P", n, n, n, &work[ir], ldwrkr, &work[itaup]
			       , &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + n;
//Perform bidiagonal QR iteration, computing left
//singular vectors of R in WORK(IU) and computing
//right singular vectors of R in WORK(IR)
//(Workspace: need 2*N*N+BDSPAC)
			Rbdsqr("U", n, n, n, 0, &s[1], &work[ie], &work[ir], ldwrkr, &work[iu], ldwrku, &dum, 1, &work[iwork], info);
//Multiply Q in U by left singular vectors of R in
//WORK(IU), storing result in A
//(Workspace: need N*N)
			Rgemm("N", "N", m, n, n, One, &u[0], ldu, &work[iu], ldwrku, Zero, &A[0], lda);
//Copy left singular vectors of A from A to U
			Rlacpy("F", m, n, &A[0], lda, &u[0], ldu);
//Copy right singular vectors of R from WORK(IR) to A
			Rlacpy("F", n, n, &work[ir], ldwrkr, &A[0], lda);
		    } else {
//Insufficient workspace for a fast algorithm
			itau = 1;
			iwork = itau + n;
//Compute A=Q*R, copying result to U
//(Workspace: need 2*N, prefer N+N*NB)
			Rgeqrf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("L", m, n, &A[0], lda, &u[0], ldu);
//Generate Q in U
//(Workspace: need N+M, prefer N+M*NB)
			Rorgqr(m, m, n, &u[0], ldu, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			ie = itau;
			itauq = ie + n;
			itaup = itauq + n;
			iwork = itaup + n;
//Zero out below R in A
			Rlaset("L", n - 1, n - 1, Zero, Zero, &A[lda + 2], lda);
//Bidiagonalize R in A
//(Workspace: need 4*N, prefer 3*N+2*N*NB)
			Rgebrd(n, n, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Multiply Q in U by left bidiagonalizing vectors
//in A
//(Workspace: need 3*N+M, prefer 3*N+M*NB)
			Rormbr("Q", "R", "N", m, n, n, &A[0], lda, &work[itauq], &u[0], ldu, &work[iwork], lwork - iwork + 1, &ierr);
//Generate right bidiagonalizing vectors in A
//(Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
			Rorgbr("P", n, n, n, &A[0], lda, &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + n;
//Perform bidiagonal QR iteration, computing left
//singular vectors of A in U and computing right
//singular vectors of A in A
//(Workspace: need BDSPAC)
			Rbdsqr("U", n, n, m, 0, &s[1], &work[ie], &A[0], lda, &u[0], ldu, &dum, 1, &work[iwork], info);
		    }
		} else if (wntvas) {
//Path 9 (m * much larger than * n, JOBU='A', JOBVT='S'
//or 'A')
//M left singular vectors to be computed in U and
//N right singular vectors to be computed in VT
		    if (lwork >= n * n + max(max(n + m, n << 2), bdspac)) {
//Sufficient workspace for a fast algorithm
			iu = 1;
			if (lwork >= wrkbl + lda * n) {
//WORK(IU) is LDA by N
			    ldwrku = lda;
			} else {
//WORK(IU) is N by N
			    ldwrku = n;
			}
			itau = iu + ldwrku * n;
			iwork = itau + n;
//Compute A=Q*R, copying result to U
//(Workspace: need N*N+2*N, prefer N*N+N+N*NB)
			Rgeqrf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("L", m, n, &A[0], lda, &u[0], ldu);
//Generate Q in U
//(Workspace: need N*N+N+M, prefer N*N+N+M*NB)
			Rorgqr(m, m, n, &u[0], ldu, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Copy R to WORK(IU), zeroing out below it
			Rlacpy("U", n, n, &A[0], lda, &work[iu], ldwrku);
			Rlaset("L", n - 1, n - 1, Zero, Zero, &work[iu + 1], ldwrku);
			ie = itau;
			itauq = ie + n;
			itaup = itauq + n;
			iwork = itaup + n;
//Bidiagonalize R in WORK(IU), copying result to VT
//(Workspace: need N*N+4*N, prefer N*N+3*N+2*N*NB)
			Rgebrd(n, n, &work[iu], ldwrku, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("U", n, n, &work[iu], ldwrku, &vt[0], ldvt);
//Generate left bidiagonalizing vectors in WORK(IU)
//(Workspace: need N*N+4*N, prefer N*N+3*N+N*NB)
			Rorgbr("Q", n, n, n, &work[iu], ldwrku, &work[itauq]
			       , &work[iwork], lwork - iwork + 1, &ierr);
//Generate right bidiagonalizing vectors in VT
//(Workspace: need N*N+4*N-1,
//            prefer N*N+3*N+(N-1)*NB)
			Rorgbr("P", n, n, n, &vt[0], ldvt, &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + n;
//Perform bidiagonal QR iteration, computing left
//singular vectors of R in WORK(IU) and computing
//right singular vectors of R in VT
//(Workspace: need N*N+BDSPAC)
			Rbdsqr("U", n, n, n, 0, &s[1], &work[ie], &vt[0], ldvt, &work[iu], ldwrku, &dum, 1, &work[iwork], info);
//Multiply Q in U by left singular vectors of R in
//WORK(IU), storing result in A
//(Workspace: need N*N)
			Rgemm("N", "N", m, n, n, One, &u[0], ldu, &work[iu], ldwrku, Zero, &A[0], lda);
//Copy left singular vectors of A from A to U
			Rlacpy("F", m, n, &A[0], lda, &u[0], ldu);
		    } else {
//Insufficient workspace for a fast algorithm
			itau = 1;
			iwork = itau + n;
//Compute A=Q*R, copying result to U
//(Workspace: need 2*N, prefer N+N*NB)
			Rgeqrf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("L", m, n, &A[0], lda, &u[0], ldu);
//Generate Q in U
//(Workspace: need N+M, prefer N+M*NB)
			Rorgqr(m, m, n, &u[0], ldu, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Copy R from A to VT, zeroing out below it
			Rlacpy("U", n, n, &A[0], lda, &vt[0], ldvt);
			if (n > 1) {
			    Rlaset("L", n - 1, n - 1, Zero, Zero, &vt[ldvt + 2], ldvt);
			}
			ie = itau;
			itauq = ie + n;
			itaup = itauq + n;
			iwork = itaup + n;
//Bidiagonalize R in VT
//(Workspace: need 4*N, prefer 3*N+2*N*NB)
			Rgebrd(n, n, &vt[0], ldvt, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Multiply Q in U by left bidiagonalizing vectors
//in VT
//(Workspace: need 3*N+M, prefer 3*N+M*NB)
			Rormbr("Q", "R", "N", m, n, n, &vt[0], ldvt, &work[itauq], &u[0], ldu, &work[iwork], lwork - iwork + 1, &ierr);
//Generate right bidiagonalizing vectors in VT
//(Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
			Rorgbr("P", n, n, n, &vt[0], ldvt, &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + n;
//Perform bidiagonal QR iteration, computing left
//singular vectors of A in U and computing right
//singular vectors of A in VT
//(Workspace: need BDSPAC)
			Rbdsqr("U", n, n, m, 0, &s[1], &work[ie], &vt[0], ldvt, &u[0], ldu, &dum, 1, &work[iwork], info);
		    }
		}
	    }
	} else {
//M .LT. MNTHR
//Path 10 (M at least N, but not much larger)
//Reduce to bidiagonal form without QR decomposition
	    ie = 1;
	    itauq = ie + n;
	    itaup = itauq + n;
	    iwork = itaup + n;
//Bidiagonalize A
//(Workspace: need 3*N+M, prefer 3*N+(M+N)*NB)
	    Rgebrd(m, n, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
	    if (wntuas) {
//If left singular vectors desired in U, copy result to U
//and generate left bidiagonalizing vectors in U
//(Workspace: need 3*N+NCU, prefer 3*N+NCU*NB)
		Rlacpy("L", m, n, &A[0], lda, &u[0], ldu);
		if (wntus) {
		    ncu = n;
		}
		if (wntua) {
		    ncu = m;
		}
		Rorgbr("Q", m, ncu, n, &u[0], ldu, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
	    }
	    if (wntvas) {
//If right singular vectors desired in VT, copy result to
//VT and generate right bidiagonalizing vectors in VT
//(Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
		Rlacpy("U", n, n, &A[0], lda, &vt[0], ldvt);
		Rorgbr("P", n, n, n, &vt[0], ldvt, &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
	    }
	    if (wntuo) {
//If left singular vectors desired in A, generate left
//bidiagonalizing vectors in A
//(Workspace: need 4*N, prefer 3*N+N*NB)
		Rorgbr("Q", m, n, n, &A[0], lda, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
	    }
	    if (wntvo) {
//If right singular vectors desired in A, generate right
//bidiagonalizing vectors in A
//(Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
		Rorgbr("P", n, n, n, &A[0], lda, &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
	    }
	    iwork = ie + n;
	    if (wntuas || wntuo) {
		nru = m;
	    }
	    if (wntun) {
		nru = 0;
	    }
	    if (wntvas || wntvo) {
		ncvt = n;
	    }
	    if (wntvn) {
		ncvt = 0;
	    }
	    if (!wntuo && !wntvo) {
//Perform bidiagonal QR iteration, if desired, computing
//left singular vectors in U and computing right singular
//vectors in VT
//(Workspace: need BDSPAC)
		Rbdsqr("U", n, ncvt, nru, 0, &s[1], &work[ie], &vt[0], ldvt, &u[0], ldu, &dum, 1, &work[iwork], info);
	    } else if (!wntuo && wntvo) {
//Perform bidiagonal QR iteration, if desired, computing
//left singular vectors in U and computing right singular
//vectors in A
//(Workspace: need BDSPAC)
		Rbdsqr("U", n, ncvt, nru, 0, &s[1], &work[ie], &A[0], lda, &u[0], ldu, &dum, 1, &work[iwork], info);
	    } else {
//Perform bidiagonal QR iteration, if desired, computing
//left singular vectors in A and computing right singular
//vectors in VT
//(Workspace: need BDSPAC)
		Rbdsqr("U", n, ncvt, nru, 0, &s[1], &work[ie], &vt[0], ldvt, &A[0], lda, &dum, 1, &work[iwork], info);
	    }
	}
    } else {
//A has more columns than rows. If A has sufficiently more
//columns than rows, first reduce using the LQ decomposition (if
//sufficient workspace available)
	if (n >= mnthr) {
	    if (wntvn) {
//Path 1t(N much larger than M, JOBVT='N')
//No right singular vectors to be computed
		itau = 1;
		iwork = itau + m;
//Compute A=L*Q
//(Workspace: need 2*M, prefer M+M*NB)
		Rgelqf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Zero out above L
		Rlaset("U", m - 1, m - 1, Zero, Zero, &A[(lda << 1) + 1], lda);
		ie = 1;
		itauq = ie + m;
		itaup = itauq + m;
		iwork = itaup + m;
//Bidiagonalize L in A
//(Workspace: need 4*M, prefer 3*M+2*M*NB)
		Rgebrd(m, m, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
		if (wntuo || wntuas) {
//   If left singular vectors desired, generate Q
//   (Workspace: need 4*M, prefer 3*M+M*NB)
		    Rorgbr("Q", m, m, m, &A[0], lda, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
		}
		iwork = ie + m;
		nru = 0;
		if (wntuo || wntuas) {
		    nru = m;
		}
//Perform bidiagonal QR iteration, computing left singular
//vectors of A in A if desired
//(Workspace: need BDSPAC)
		Rbdsqr("U", m, 0, nru, 0, &s[1], &work[ie], &dum, 1, &A[0], lda, &dum, 1, &work[iwork], info);
//If left singular vectors desired in U, copy them there
		if (wntuas) {
		    Rlacpy("F", m, m, &A[0], lda, &u[0], ldu);
		}
	    } else if (wntvo && wntun) {
//Path 2t(N much larger than M, JOBU='N', JOBVT='O')
//M right singular vectors to be overwritten on A and
//no left singular vectors to be computed
		if (lwork >= m * m + max(m << 2, bdspac)) {
//Sufficient workspace for a fast algorithm
		    ir = 1;
		    if (lwork >= max(wrkbl, lda * n + m) + lda * m) {
//WORK(IU) is LDA by N and WORK(IR) is LDA by M
			ldwrku = lda;
			chunk = n;
			ldwrkr = lda;
		    } else {
			if (lwork >= max(wrkbl, lda * n + m) + m * m) {
//WORK(IU) is LDA by N and WORK(IR) is M by M
			    ldwrku = lda;
			    chunk = n;
			    ldwrkr = m;
			} else {
//WORK(IU) is M by CHUNK and WORK(IR) is M by M
			    ldwrku = m;
			    chunk = (lwork - m * m - m) / m;
			    ldwrkr = m;
			}
		    }
		    itau = ir + ldwrkr * m;
		    iwork = itau + m;
//Compute A=L*Q
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
		    Rgelqf(m, n, &A[0], lda, &work[itau], &work[iwork]
			   , lwork - iwork + 1, &ierr);
//Copy L to WORK(IR) and zero out above it
		    Rlacpy("L", m, m, &A[0], lda, &work[ir], ldwrkr);
		    Rlaset("U", m - 1, m - 1, Zero, Zero, &work[ir + ldwrkr], ldwrkr);
//Generate Q in A
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
		    Rorglq(m, n, m, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
		    ie = itau;
		    itauq = ie + m;
		    itaup = itauq + m;
		    iwork = itaup + m;
//Bidiagonalize L in WORK(IR)
//(Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
		    Rgebrd(m, m, &work[ir], ldwrkr, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Generate right vectors bidiagonalizing L
//(Workspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB)
		    Rorgbr("P", m, m, m, &work[ir], ldwrkr, &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
		    iwork = ie + m;
//Perform bidiagonal QR iteration, computing right
//singular vectors of L in WORK(IR)
//(Workspace: need M*M+BDSPAC)
		    Rbdsqr("U", m, m, 0, 0, &s[1], &work[ie], &work[ir], ldwrkr, &dum, 1, &dum, 1, &work[iwork]
			   , info);
		    iu = ie + m;
//Multiply right singular vectors of L in WORK(IR) by Q
//in A, storing result in WORK(IU) and copying to A
//(Workspace: need M*M+2*M, prefer M*M+M*N+M)
		    for (i = 1; i <= n; i += chunk) {
			blk = min(n - i + 1, chunk);
			Rgemm("N", "N", m, blk, m, One, &work[ir], ldwrkr, &A[i * lda + 1], lda, Zero, &work[iu], ldwrku);
			Rlacpy("F", m, blk, &work[iu], ldwrku, &A[i * lda + 1], lda);
		    }
		} else {
//Insufficient workspace for a fast algorithm
		    ie = 1;
		    itauq = ie + m;
		    itaup = itauq + m;
		    iwork = itaup + m;
//Bidiagonalize A
//(Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
		    Rgebrd(m, n, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Generate right vectors bidiagonalizing A
//(Workspace: need 4*M, prefer 3*M+M*NB)
		    Rorgbr("P", m, n, m, &A[0], lda, &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
		    iwork = ie + m;
//Perform bidiagonal QR iteration, computing right
//singular vectors of A in A
//(Workspace: need BDSPAC)
		    Rbdsqr("L", m, n, 0, 0, &s[1], &work[ie], &A[0], lda, &dum, 1, &dum, 1, &work[iwork], info);
		}
	    } else if (wntvo && wntuas) {
//Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O')
//M right singular vectors to be overwritten on A and
//M left singular vectors to be computed in U
		if (lwork >= m * m + max(m << 2, bdspac)) {
//Sufficient workspace for a fast algorithm
		    ir = 1;
		    if (lwork >= max(wrkbl, lda * n + m) + lda * m) {
//WORK(IU) is LDA by N and WORK(IR) is LDA by M
			ldwrku = lda;
			chunk = n;
			ldwrkr = lda;
		    } else {
			if (lwork >= max(wrkbl, lda * n + m) + m * m) {
//WORK(IU) is LDA by N and WORK(IR) is M by M
			    ldwrku = lda;
			    chunk = n;
			    ldwrkr = m;
			} else {
//WORK(IU) is M by CHUNK and WORK(IR) is M by M
			    ldwrku = m;
			    chunk = (lwork - m * m - m) / m;
			    ldwrkr = m;
			}
		    }
		    itau = ir + ldwrkr * m;
		    iwork = itau + m;
//Compute A=L*Q
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
		    Rgelqf(m, n, &A[0], lda, &work[itau], &work[iwork]
			   , lwork - iwork + 1, &ierr);
//Copy L to U, zeroing about above it
		    Rlacpy("L", m, m, &A[0], lda, &u[0], ldu);
		    Rlaset("U", m - 1, m - 1, Zero, Zero, &u[(ldu << 1) + 1], ldu);
//Generate Q in A
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
		    Rorglq(m, n, m, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
		    ie = itau;
		    itauq = ie + m;
		    itaup = itauq + m;
		    iwork = itaup + m;
//Bidiagonalize L in U, copying result to WORK(IR)
//(Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
		    Rgebrd(m, m, &u[0], ldu, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
		    Rlacpy("U", m, m, &u[0], ldu, &work[ir], ldwrkr);
//Generate right vectors bidiagonalizing L in WORK(IR)
//(Workspace: need M*M+4*M-1, prefer M*M+3*M+(M-1)*NB)
		    Rorgbr("P", m, m, m, &work[ir], ldwrkr, &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Generate left vectors bidiagonalizing L in U
//(Workspace: need M*M+4*M, prefer M*M+3*M+M*NB)
		    Rorgbr("Q", m, m, m, &u[0], ldu, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
		    iwork = ie + m;
//Perform bidiagonal QR iteration, computing left
//singular vectors of L in U, and computing right
//singular vectors of L in WORK(IR)
//(Workspace: need M*M+BDSPAC)
		    Rbdsqr("U", m, m, m, 0, &s[1], &work[ie], &work[ir], ldwrkr, &u[0], ldu, &dum, 1, &work[iwork], info);
		    iu = ie + m;
//Multiply right singular vectors of L in WORK(IR) by Q
//in A, storing result in WORK(IU) and copying to A
//(Workspace: need M*M+2*M, prefer M*M+M*N+M))
		    for (i = 1; i <= n; i += chunk) {
			blk = min(n - i + 1, chunk);
			Rgemm("N", "N", m, blk, m, One, &work[ir], ldwrkr, &A[i * lda + 1], lda, Zero, &work[iu], ldwrku);
			Rlacpy("F", m, blk, &work[iu], ldwrku, &A[i * lda + 1], lda);
		    }
		} else {
//Insufficient workspace for a fast algorithm
		    itau = 1;
		    iwork = itau + m;
//Compute A=L*Q
//(Workspace: need 2*M, prefer M+M*NB)
		    Rgelqf(m, n, &A[0], lda, &work[itau], &work[iwork]
			   , lwork - iwork + 1, &ierr);
//Copy L to U, zeroing out above it
		    Rlacpy("L", m, m, &A[0], lda, &u[0], ldu);
		    Rlaset("U", m - 1, m - 1, Zero, Zero, &u[(ldu << 1) + 1], ldu);
//Generate Q in A
//(Workspace: need 2*M, prefer M+M*NB)
		    Rorglq(m, n, m, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
		    ie = itau;
		    itauq = ie + m;
		    itaup = itauq + m;
		    iwork = itaup + m;
//Bidiagonalize L in U
//(Workspace: need 4*M, prefer 3*M+2*M*NB)
		    Rgebrd(m, m, &u[0], ldu, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Multiply right vectors bidiagonalizing L by Q in A
//(Workspace: need 3*M+N, prefer 3*M+N*NB)
		    Rormbr("P", "L", "T", m, n, m, &u[0], ldu, &work[itaup], &A[0], lda, &work[iwork], lwork - iwork + 1, &ierr);
//Generate left vectors bidiagonalizing L in U
//(Workspace: need 4*M, prefer 3*M+M*NB)
		    Rorgbr("Q", m, m, m, &u[0], ldu, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
		    iwork = ie + m;
//Perform bidiagonal QR iteration, computing left
//singular vectors of A in U and computing right
//singular vectors of A in A
//(Workspace: need BDSPAC)
		    Rbdsqr("U", m, n, m, 0, &s[1], &work[ie], &A[0], lda, &u[0], ldu, &dum, 1, &work[iwork], info);
		}
	    } else if (wntvs) {
		if (wntun) {
//Path 4t(N much larger than M, JOBU='N', JOBVT='S')
//M right singular vectors to be computed in VT and
//no left singular vectors to be computed
		    if (lwork >= m * m + max(m << 2, bdspac)) {
//Sufficient workspace for a fast algorithm
			ir = 1;
			if (lwork >= wrkbl + lda * m) {
//WORK(IR) is LDA by M
			    ldwrkr = lda;
			} else {
//WORK(IR) is M by M
			    ldwrkr = m;
			}
			itau = ir + ldwrkr * m;
			iwork = itau + m;
//Compute A=L*Q
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
			Rgelqf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Copy L to WORK(IR), zeroing out above it
			Rlacpy("L", m, m, &A[0], lda, &work[ir], ldwrkr);
			Rlaset("U", m - 1, m - 1, Zero, Zero, &work[ir + ldwrkr], ldwrkr);
//Generate Q in A
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
			Rorglq(m, n, m, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			ie = itau;
			itauq = ie + m;
			itaup = itauq + m;
			iwork = itaup + m;
//Bidiagonalize L in WORK(IR)
//(Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
			Rgebrd(m, m, &work[ir], ldwrkr, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Generate right vectors bidiagonalizing L in
//WORK(IR)
//(Workspace: need M*M+4*M, prefer M*M+3*M+(M-1)*NB)
			Rorgbr("P", m, m, m, &work[ir], ldwrkr, &work[itaup]
			       , &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + m;
//Perform bidiagonal QR iteration, computing right
//singular vectors of L in WORK(IR)
//(Workspace: need M*M+BDSPAC)
			Rbdsqr("U", m, m, 0, 0, &s[1], &work[ie], &work[ir], ldwrkr, &dum, 1, &dum, 1, &work[iwork], info);
//Multiply right singular vectors of L in WORK(IR) by
//Q in A, storing result in VT
//(Workspace: need M*M)
			Rgemm("N", "N", m, n, m, One, &work[ir], ldwrkr, &A[0], lda, Zero, &vt[0], ldvt);
		    } else {
//Insufficient workspace for a fast algorithm
			itau = 1;
			iwork = itau + m;
//Compute A=L*Q
//(Workspace: need 2*M, prefer M+M*NB)
			Rgelqf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Copy result to VT
			Rlacpy("U", m, n, &A[0], lda, &vt[0], ldvt);
//Generate Q in VT
//(Workspace: need 2*M, prefer M+M*NB)
			Rorglq(m, n, m, &vt[0], ldvt, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			ie = itau;
			itauq = ie + m;
			itaup = itauq + m;
			iwork = itaup + m;
//Zero out above L in A
			Rlaset("U", m - 1, m - 1, Zero, Zero, &A[(lda << 1) + 1], lda);
//Bidiagonalize L in A
//(Workspace: need 4*M, prefer 3*M+2*M*NB)
			Rgebrd(m, m, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Multiply right vectors bidiagonalizing L by Q in VT
//(Workspace: need 3*M+N, prefer 3*M+N*NB)
			Rormbr("P", "L", "T", m, n, m, &A[0], lda, &work[itaup], &vt[0], ldvt, &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + m;
//Perform bidiagonal QR iteration, computing right
//singular vectors of A in VT
//(Workspace: need BDSPAC)
			Rbdsqr("U", m, n, 0, 0, &s[1], &work[ie], &vt[0], ldvt, &dum, 1, &dum, 1, &work[iwork], info);
		    }
		} else if (wntuo) {
//Path 5t(N much larger than M, JOBU='O', JOBVT='S')
//M right singular vectors to be computed in VT and
//M left singular vectors to be overwritten on A
		    if (lwork >= (m << 1) * m + max(m << 2, bdspac)) {
//Sufficient workspace for a fast algorithm
			iu = 1;
			if (lwork >= wrkbl + (lda << 1) * m) {
//WORK(IU) is LDA by M and WORK(IR) is LDA by M
			    ldwrku = lda;
			    ir = iu + ldwrku * m;
			    ldwrkr = lda;
			} else if (lwork >= wrkbl + (lda + m) * m) {
//WORK(IU) is LDA by M and WORK(IR) is M by M
			    ldwrku = lda;
			    ir = iu + ldwrku * m;
			    ldwrkr = m;
			} else {
//WORK(IU) is M by M and WORK(IR) is M by M
			    ldwrku = m;
			    ir = iu + ldwrku * m;
			    ldwrkr = m;
			}
			itau = ir + ldwrkr * m;
			iwork = itau + m;
//Compute A=L*Q
//(Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
			Rgelqf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Copy L to WORK(IU), zeroing out below it
			Rlacpy("L", m, m, &A[0], lda, &work[iu], ldwrku);
			Rlaset("U", m - 1, m - 1, Zero, Zero, &work[iu + ldwrku], ldwrku);
//Generate Q in A
//(Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
			Rorglq(m, n, m, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			ie = itau;
			itauq = ie + m;
			itaup = itauq + m;
			iwork = itaup + m;
//Bidiagonalize L in WORK(IU), copying result to
//WORK(IR)
//(Workspace: need 2*M*M+4*M,
//         prefer 2*M*M+3*M+2*M*NB)
			Rgebrd(m, m, &work[iu], ldwrku, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("L", m, m, &work[iu], ldwrku, &work[ir], ldwrkr);
//Generate right bidiagonalizing vectors in WORK(IU)
//(Workspace: need 2*M*M+4*M-1,
//         prefer 2*M*M+3*M+(M-1)*NB)
			Rorgbr("P", m, m, m, &work[iu], ldwrku, &work[itaup]
			       , &work[iwork], lwork - iwork + 1, &ierr);
//Generate left bidiagonalizing vectors in WORK(IR)
//(Workspace: need 2*M*M+4*M, prefer 2*M*M+3*M+M*NB)
			Rorgbr("Q", m, m, m, &work[ir], ldwrkr, &work[itauq]
			       , &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + m;
//Perform bidiagonal QR iteration, computing left
//singular vectors of L in WORK(IR) and computing
//right singular vectors of L in WORK(IU)
//(Workspace: need 2*M*M+BDSPAC)
			Rbdsqr("U", m, m, m, 0, &s[1], &work[ie], &work[iu], ldwrku, &work[ir], ldwrkr, &dum, 1, &work[iwork], info);
//Multiply right singular vectors of L in WORK(IU) by
//Q in A, storing result in VT
//(Workspace: need M*M)
			Rgemm("N", "N", m, n, m, One, &work[iu], ldwrku, &A[0], lda, Zero, &vt[0], ldvt);
//Copy left singular vectors of L to A
//(Workspace: need M*M)
			Rlacpy("F", m, m, &work[ir], ldwrkr, &A[0], lda);
		    } else {
//Insufficient workspace for a fast algorithm
			itau = 1;
			iwork = itau + m;
//Compute A=L*Q, copying result to VT
//(Workspace: need 2*M, prefer M+M*NB)
			Rgelqf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("U", m, n, &A[0], lda, &vt[0], ldvt);
//Generate Q in VT
//(Workspace: need 2*M, prefer M+M*NB)
			Rorglq(m, n, m, &vt[0], ldvt, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			ie = itau;
			itauq = ie + m;
			itaup = itauq + m;
			iwork = itaup + m;
//Zero out above L in A
			Rlaset("U", m - 1, m - 1, Zero, Zero, &A[(lda << 1) + 1], lda);
//Bidiagonalize L in A
//(Workspace: need 4*M, prefer 3*M+2*M*NB)
			Rgebrd(m, m, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Multiply right vectors bidiagonalizing L by Q in VT
//(Workspace: need 3*M+N, prefer 3*M+N*NB)
			Rormbr("P", "L", "T", m, n, m, &A[0], lda, &work[itaup], &vt[0], ldvt, &work[iwork], lwork - iwork + 1, &ierr);
//Generate left bidiagonalizing vectors of L in A
//(Workspace: need 4*M, prefer 3*M+M*NB)
			Rorgbr("Q", m, m, m, &A[0], lda, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + m;
//Perform bidiagonal QR iteration, compute left
//singular vectors of A in A and compute right
//singular vectors of A in VT
//(Workspace: need BDSPAC)
			Rbdsqr("U", m, n, m, 0, &s[1], &work[ie], &vt[0], ldvt, &A[0], lda, &dum, 1, &work[iwork], info);
		    }
		} else if (wntuas) {
//Path 6t(N much larger than M, JOBU='S' or 'A',
//JOBVT='S')
//M right singular vectors to be computed in VT and
//M left singular vectors to be computed in U
		    if (lwork >= m * m + max(m << 2, bdspac)) {
//Sufficient workspace for a fast algorithm
			iu = 1;
			if (lwork >= wrkbl + lda * m) {
//WORK(IU) is LDA by N
			    ldwrku = lda;
			} else {
//WORK(IU) is LDA by M
			    ldwrku = m;
			}
			itau = iu + ldwrku * m;
			iwork = itau + m;
//Compute A=L*Q
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
			Rgelqf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Copy L to WORK(IU), zeroing out above it
			Rlacpy("L", m, m, &A[0], lda, &work[iu], ldwrku);
			Rlaset("U", m - 1, m - 1, Zero, Zero, &work[iu + ldwrku], ldwrku);
//Generate Q in A
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
			Rorglq(m, n, m, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			ie = itau;
			itauq = ie + m;
			itaup = itauq + m;
			iwork = itaup + m;
//Bidiagonalize L in WORK(IU), copying result to U
//(Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
			Rgebrd(m, m, &work[iu], ldwrku, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("L", m, m, &work[iu], ldwrku, &u[0], ldu);
//Generate right bidiagonalizing vectors in WORK(IU)
//(Workspace: need M*M+4*M-1,
//         prefer M*M+3*M+(M-1)*NB)
			Rorgbr("P", m, m, m, &work[iu], ldwrku, &work[itaup]
			       , &work[iwork], lwork - iwork + 1, &ierr);
//Generate left bidiagonalizing vectors in U
//(Workspace: need M*M+4*M, prefer M*M+3*M+M*NB)
			Rorgbr("Q", m, m, m, &u[0], ldu, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + m;
//Perform bidiagonal QR iteration, computing left
//singular vectors of L in U and computing right
//singular vectors of L in WORK(IU)
//(Workspace: need M*M+BDSPAC)
			Rbdsqr("U", m, m, m, 0, &s[1], &work[ie], &work[iu], ldwrku, &u[0], ldu, &dum, 1, &work[iwork], info);
//Multiply right singular vectors of L in WORK(IU) by
//Q in A, storing result in VT
//(Workspace: need M*M)
			Rgemm("N", "N", m, n, m, One, &work[iu], ldwrku, &A[0], lda, Zero, &vt[0], ldvt);
		    } else {
//Insufficient workspace for a fast algorithm
			itau = 1;
			iwork = itau + m;
//Compute A=L*Q, copying result to VT
//(Workspace: need 2*M, prefer M+M*NB)
			Rgelqf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("U", m, n, &A[0], lda, &vt[0], ldvt);
//Generate Q in VT
//(Workspace: need 2*M, prefer M+M*NB)
			Rorglq(m, n, m, &vt[0], ldvt, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Copy L to U, zeroing out above it
			Rlacpy("L", m, m, &A[0], lda, &u[0], ldu);
			Rlaset("U", m - 1, m - 1, Zero, Zero, &u[(ldu << 1) + 1], ldu);
			ie = itau;
			itauq = ie + m;
			itaup = itauq + m;
			iwork = itaup + m;
//Bidiagonalize L in U
//(Workspace: need 4*M, prefer 3*M+2*M*NB)
			Rgebrd(m, m, &u[0], ldu, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Multiply right bidiagonalizing vectors in U by Q
//in VT
//(Workspace: need 3*M+N, prefer 3*M+N*NB)
			Rormbr("P", "L", "T", m, n, m, &u[0], ldu, &work[itaup], &vt[0], ldvt, &work[iwork], lwork - iwork + 1, &ierr);
//Generate left bidiagonalizing vectors in U
//(Workspace: need 4*M, prefer 3*M+M*NB)
			Rorgbr("Q", m, m, m, &u[0], ldu, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + m;
//Perform bidiagonal QR iteration, computing left
//singular vectors of A in U and computing right
//singular vectors of A in VT
//(Workspace: need BDSPAC)
			Rbdsqr("U", m, n, m, 0, &s[1], &work[ie], &vt[0], ldvt, &u[0], ldu, &dum, 1, &work[iwork], info);
		    }
		}
	    } else if (wntva) {
		if (wntun) {
//Path 7t(N much larger than M, JOBU='N', JOBVT='A')
//N right singular vectors to be computed in VT and
//no left singular vectors to be computed
		    if (lwork >= m * m + max(max(n + m, m << 2), bdspac)) {
//Sufficient workspace for a fast algorithm
			ir = 1;
			if (lwork >= wrkbl + lda * m) {
//WORK(IR) is LDA by M
			    ldwrkr = lda;
			} else {
//WORK(IR) is M by M
			    ldwrkr = m;
			}
			itau = ir + ldwrkr * m;
			iwork = itau + m;
//Compute A=L*Q, copying result to VT
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
			Rgelqf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("U", m, n, &A[0], lda, &vt[0], ldvt);
//Copy L to WORK(IR), zeroing out above it
			Rlacpy("L", m, m, &A[0], lda, &work[ir], ldwrkr);
			Rlaset("U", m - 1, m - 1, Zero, Zero, &work[ir + ldwrkr], ldwrkr);
//Generate Q in VT
//(Workspace: need M*M+M+N, prefer M*M+M+N*NB)
			Rorglq(n, n, m, &vt[0], ldvt, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			ie = itau;
			itauq = ie + m;
			itaup = itauq + m;
			iwork = itaup + m;
//Bidiagonalize L in WORK(IR)
//(Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
			Rgebrd(m, m, &work[ir], ldwrkr, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Generate right bidiagonalizing vectors in WORK(IR)
//(Workspace: need M*M+4*M-1,
//         prefer M*M+3*M+(M-1)*NB)
			Rorgbr("P", m, m, m, &work[ir], ldwrkr, &work[itaup]
			       , &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + m;
//Perform bidiagonal QR iteration, computing right
//singular vectors of L in WORK(IR)
//(Workspace: need M*M+BDSPAC)
			Rbdsqr("U", m, m, 0, 0, &s[1], &work[ie], &work[ir], ldwrkr, &dum, 1, &dum, 1, &work[iwork], info);
//Multiply right singular vectors of L in WORK(IR) by
//Q in VT, storing result in A
//(Workspace: need M*M)
			Rgemm("N", "N", m, n, m, One, &work[ir], ldwrkr, &vt[0], ldvt, Zero, &A[0], lda);
//Copy right singular vectors of A from A to VT
			Rlacpy("F", m, n, &A[0], lda, &vt[0], ldvt);
		    } else {
//Insufficient workspace for a fast algorithm
			itau = 1;
			iwork = itau + m;
//Compute A=L*Q, copying result to VT
//(Workspace: need 2*M, prefer M+M*NB)
			Rgelqf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("U", m, n, &A[0], lda, &vt[0], ldvt);
//Generate Q in VT
//(Workspace: need M+N, prefer M+N*NB)
			Rorglq(n, n, m, &vt[0], ldvt, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			ie = itau;
			itauq = ie + m;
			itaup = itauq + m;
			iwork = itaup + m;
//Zero out above L in A
			Rlaset("U", m - 1, m - 1, Zero, Zero, &A[(lda << 1) + 1], lda);

//Bidiagonalize L in A
//(Workspace: need 4*M, prefer 3*M+2*M*NB)
			Rgebrd(m, m, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Multiply right bidiagonalizing vectors in A by Q
//in VT
//(Workspace: need 3*M+N, prefer 3*M+N*NB)
			Rormbr("P", "L", "T", m, n, m, &A[0], lda, &work[itaup], &vt[0], ldvt, &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + m;
//Perform bidiagonal QR iteration, computing right
//singular vectors of A in VT
//(Workspace: need BDSPAC)
			Rbdsqr("U", m, n, 0, 0, &s[1], &work[ie], &vt[0], ldvt, &dum, 1, &dum, 1, &work[iwork], info);
		    }
		} else if (wntuo) {
//Path 8t(N much larger than M, JOBU='O', JOBVT='A')
//N right singular vectors to be computed in VT and
//M left singular vectors to be overwritten on A
		    if (lwork >= (m << 1) * m + max(max(n + m, m << 2), bdspac)) {
//Sufficient workspace for a fast algorithm
			iu = 1;
			if (lwork >= wrkbl + (lda << 1) * m) {
//WORK(IU) is LDA by M and WORK(IR) is LDA by M
			    ldwrku = lda;
			    ir = iu + ldwrku * m;
			    ldwrkr = lda;
			} else if (lwork >= wrkbl + (lda + m) * m) {
//WORK(IU) is LDA by M and WORK(IR) is M by M
			    ldwrku = lda;
			    ir = iu + ldwrku * m;
			    ldwrkr = m;
			} else {
//WORK(IU) is M by M and WORK(IR) is M by M
			    ldwrku = m;
			    ir = iu + ldwrku * m;
			    ldwrkr = m;
			}
			itau = ir + ldwrkr * m;
			iwork = itau + m;
//Compute A=L*Q, copying result to VT
//(Workspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
			Rgelqf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("U", m, n, &A[0], lda, &vt[0], ldvt);
//Generate Q in VT
//(Workspace: need 2*M*M+M+N, prefer 2*M*M+M+N*NB)
			Rorglq(n, n, m, &vt[0], ldvt, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Copy L to WORK(IU), zeroing out above it
			Rlacpy("L", m, m, &A[0], lda, &work[iu], ldwrku);
			Rlaset("U", m - 1, m - 1, Zero, Zero, &work[iu + ldwrku], ldwrku);
			ie = itau;
			itauq = ie + m;
			itaup = itauq + m;
			iwork = itaup + m;
//Bidiagonalize L in WORK(IU), copying result to
//WORK(IR)
//(Workspace: need 2*M*M+4*M,
//         prefer 2*M*M+3*M+2*M*NB)
			Rgebrd(m, m, &work[iu], ldwrku, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("L", m, m, &work[iu], ldwrku, &work[ir], ldwrkr);
//Generate right bidiagonalizing vectors in WORK(IU)
//(Workspace: need 2*M*M+4*M-1,
//         prefer 2*M*M+3*M+(M-1)*NB)
			Rorgbr("P", m, m, m, &work[iu], ldwrku, &work[itaup]
			       , &work[iwork], lwork - iwork + 1, &ierr);
//Generate left bidiagonalizing vectors in WORK(IR)
//(Workspace: need 2*M*M+4*M, prefer 2*M*M+3*M+M*NB)
			Rorgbr("Q", m, m, m, &work[ir], ldwrkr, &work[itauq]
			       , &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + m;
//Perform bidiagonal QR iteration, computing left
//singular vectors of L in WORK(IR) and computing
//right singular vectors of L in WORK(IU)
//(Workspace: need 2*M*M+BDSPAC)
			Rbdsqr("U", m, m, m, 0, &s[1], &work[ie], &work[iu], ldwrku, &work[ir], ldwrkr, &dum, 1, &work[iwork], info);
//Multiply right singular vectors of L in WORK(IU) by
//Q in VT, storing result in A
//(Workspace: need M*M)
			Rgemm("N", "N", m, n, m, One, &work[iu], ldwrku, &vt[0], ldvt, Zero, &A[0], lda);
//Copy right singular vectors of A from A to VT
			Rlacpy("F", m, n, &A[0], lda, &vt[0], ldvt);
//Copy left singular vectors of A from WORK(IR) to A
			Rlacpy("F", m, m, &work[ir], ldwrkr, &A[0], lda);
		    } else {
//Insufficient workspace for a fast algorithm
			itau = 1;
			iwork = itau + m;
//Compute A=L*Q, copying result to VT
//(Workspace: need 2*M, prefer M+M*NB)
			Rgelqf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("U", m, n, &A[0], lda, &vt[0], ldvt);
//Generate Q in VT
//(Workspace: need M+N, prefer M+N*NB)
			Rorglq(n, n, m, &vt[0], ldvt, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			ie = itau;
			itauq = ie + m;
			itaup = itauq + m;
			iwork = itaup + m;
//Zero out above L in A
			Rlaset("U", m - 1, m - 1, Zero, Zero, &A[(lda << 1) + 1], lda);
//Bidiagonalize L in A
//(Workspace: need 4*M, prefer 3*M+2*M*NB)
			Rgebrd(m, m, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Multiply right bidiagonalizing vectors in A by Q
//in VT
//(Workspace: need 3*M+N, prefer 3*M+N*NB)
			Rormbr("P", "L", "T", m, n, m, &A[0], lda, &work[itaup], &vt[0], ldvt, &work[iwork], lwork - iwork + 1, &ierr);
//Generate left bidiagonalizing vectors in A
//(Workspace: need 4*M, prefer 3*M+M*NB)
			Rorgbr("Q", m, m, m, &A[0], lda, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + m;
//Perform bidiagonal QR iteration, computing left
//singular vectors of A in A and computing right
//singular vectors of A in VT
//(Workspace: need BDSPAC)
			Rbdsqr("U", m, n, m, 0, &s[1], &work[ie], &vt[0], ldvt, &A[0], lda, &dum, 1, &work[iwork], info);
		    }
		} else if (wntuas) {
//Path 9t(N much larger than M, JOBU='S' or 'A',
//  JOBVT='A')
//   N right singular vectors to be computed in VT and
//   M left singular vectors to be computed in U
		    if (lwork >= m * m + max(max(n + m, m << 2), bdspac)) {
//Sufficient workspace for a fast algorithm
			iu = 1;
			if (lwork >= wrkbl + lda * m) {
//WORK(IU) is LDA by M
			    ldwrku = lda;
			} else {
//WORK(IU) is M by M
			    ldwrku = m;
			}
			itau = iu + ldwrku * m;
			iwork = itau + m;
//Compute A=L*Q, copying result to VT
//(Workspace: need M*M+2*M, prefer M*M+M+M*NB)
			Rgelqf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("U", m, n, &A[0], lda, &vt[0], ldvt);
//Generate Q in VT
//(Workspace: need M*M+M+N, prefer M*M+M+N*NB)
			Rorglq(n, n, m, &vt[0], ldvt, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Copy L to WORK(IU), zeroing out above it
			Rlacpy("L", m, m, &A[0], lda, &work[iu], ldwrku);
			Rlaset("U", m - 1, m - 1, Zero, Zero, &work[iu + ldwrku], ldwrku);
			ie = itau;
			itauq = ie + m;
			itaup = itauq + m;
			iwork = itaup + m;
//Bidiagonalize L in WORK(IU), copying result to U
//(Workspace: need M*M+4*M, prefer M*M+3*M+2*M*NB)
			Rgebrd(m, m, &work[iu], ldwrku, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("L", m, m, &work[iu], ldwrku, &u[0], ldu);
//Generate right bidiagonalizing vectors in WORK(IU)
//(Workspace: need M*M+4*M, prefer M*M+3*M+(M-1)*NB)
			Rorgbr("P", m, m, m, &work[iu], ldwrku, &work[itaup]
			       , &work[iwork], lwork - iwork + 1, &ierr);
//Generate left bidiagonalizing vectors in U
//(Workspace: need M*M+4*M, prefer M*M+3*M+M*NB)
			Rorgbr("Q", m, m, m, &u[0], ldu, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + m;
//Perform bidiagonal QR iteration, computing left
//singular vectors of L in U and computing right
//singular vectors of L in WORK(IU)
//(Workspace: need M*M+BDSPAC)
			Rbdsqr("U", m, m, m, 0, &s[1], &work[ie], &work[iu], ldwrku, &u[0], ldu, &dum, 1, &work[iwork], info);
//Multiply right singular vectors of L in WORK(IU) by
//Q in VT, storing result in A
//(Workspace: need M*M)
			Rgemm("N", "N", m, n, m, One, &work[iu], ldwrku, &vt[0], ldvt, Zero, &A[0], lda);
//Copy right singular vectors of A from A to VT
			Rlacpy("F", m, n, &A[0], lda, &vt[0], ldvt);
		    } else {
//Insufficient workspace for a fast algorithm
			itau = 1;
			iwork = itau + m;
//Compute A=L*Q, copying result to VT
//(Workspace: need 2*M, prefer M+M*NB)
			Rgelqf(m, n, &A[0], lda, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
			Rlacpy("U", m, n, &A[0], lda, &vt[0], ldvt);
//Generate Q in VT
//(Workspace: need M+N, prefer M+N*NB)
			Rorglq(n, n, m, &vt[0], ldvt, &work[itau], &work[iwork], lwork - iwork + 1, &ierr);
//Copy L to U, zeroing out above it
			Rlacpy("L", m, m, &A[0], lda, &u[0], ldu);
			Rlaset("U", m - 1, m - 1, Zero, Zero, &u[(ldu << 1) + 1], ldu);
			ie = itau;
			itauq = ie + m;
			itaup = itauq + m;
			iwork = itaup + m;
//Bidiagonalize L in U
//(Workspace: need 4*M, prefer 3*M+2*M*NB)
			Rgebrd(m, m, &u[0], ldu, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
//Multiply right bidiagonalizing vectors in U by Q
//in VT
//(Workspace: need 3*M+N, prefer 3*M+N*NB)
			Rormbr("P", "L", "T", m, n, m, &u[0], ldu, &work[itaup], &vt[0], ldvt, &work[iwork], lwork - iwork + 1, &ierr);
//Generate left bidiagonalizing vectors in U
//(Workspace: need 4*M, prefer 3*M+M*NB)
			Rorgbr("Q", m, m, m, &u[0], ldu, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
			iwork = ie + m;

//Perform bidiagonal QR iteration, computing left
//singular vectors of A in U and computing right
//singular vectors of A in VT
//(Workspace: need BDSPAC)
			Rbdsqr("U", m, n, m, 0, &s[1], &work[ie], &vt[0], ldvt, &u[0], ldu, &dum, 1, &work[iwork], info);
		    }
		}
	    }
	} else {
//N .LT. MNTHR
//Path 10t(N greater than M, but not much larger)
//Reduce to bidiagonal form without LQ decomposition
	    ie = 1;
	    itauq = ie + m;
	    itaup = itauq + m;
	    iwork = itaup + m;
//Bidiagonalize A
//(Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
	    Rgebrd(m, n, &A[0], lda, &s[1], &work[ie], &work[itauq], &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
	    if (wntuas) {
//If left singular vectors desired in U, copy result to U
//and generate left bidiagonalizing vectors in U
//(Workspace: need 4*M-1, prefer 3*M+(M-1)*NB)
		Rlacpy("L", m, m, &A[0], lda, &u[0], ldu);
		Rorgbr("Q", m, m, n, &u[0], ldu, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
	    }
	    if (wntvas) {
//If right singular vectors desired in VT, copy result to
//VT and generate right bidiagonalizing vectors in VT
//(Workspace: need 3*M+NRVT, prefer 3*M+NRVT*NB)
		Rlacpy("U", m, n, &A[0], lda, &vt[0], ldvt);
		if (wntva) {
		    nrvt = n;
		}
		if (wntvs) {
		    nrvt = m;
		}
		Rorgbr("P", nrvt, n, m, &vt[0], ldvt, &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
	    }
	    if (wntuo) {
//If left singular vectors desired in A, generate left
//bidiagonalizing vectors in A
//(Workspace: need 4*M-1, prefer 3*M+(M-1)*NB)
		Rorgbr("Q", m, m, n, &A[0], lda, &work[itauq], &work[iwork], lwork - iwork + 1, &ierr);
	    }
	    if (wntvo) {
//If right singular vectors desired in A, generate right
//bidiagonalizing vectors in A
//(Workspace: need 4*M, prefer 3*M+M*NB)
		Rorgbr("P", m, n, m, &A[0], lda, &work[itaup], &work[iwork], lwork - iwork + 1, &ierr);
	    }
	    iwork = ie + m;
	    if (wntuas || wntuo) {
		nru = m;
	    }
	    if (wntun) {
		nru = 0;
	    }
	    if (wntvas || wntvo) {
		ncvt = n;
	    }
	    if (wntvn) {
		ncvt = 0;
	    }
	    if (!wntuo && !wntvo) {
//Perform bidiagonal QR iteration, if desired, computing
//left singular vectors in U and computing right singular
//vectors in VT
//(Workspace: need BDSPAC)
		Rbdsqr("L", m, ncvt, nru, 0, &s[1], &work[ie], &vt[0], ldvt, &u[0], ldu, &dum, 1, &work[iwork], info);
	    } else if (!wntuo && wntvo) {
//Perform bidiagonal QR iteration, if desired, computing
//left singular vectors in U and computing right singular
//vectors in A
//(Workspace: need BDSPAC)
		Rbdsqr("L", m, ncvt, nru, 0, &s[1], &work[ie], &A[0], lda, &u[0], ldu, &dum, 1, &work[iwork], info);
	    } else {

//Perform bidiagonal QR iteration, if desired, computing
//left singular vectors in A and computing right singular
//vectors in VT
//(Workspace: need BDSPAC)
		Rbdsqr("L", m, ncvt, nru, 0, &s[1], &work[ie], &vt[0], ldvt, &A[0], lda, &dum, 1, &work[iwork], info);
	    }
	}
    }
//If DBDSQR failed to converge, copy unconverged superdiagonals
//to WORK( 2:MINMN )!
    if (*info != 0) {
	if (ie > 2) {
	    for (i = 0; i < minmn - 1; i++) {
		work[i + 1] = work[i + ie - 1];
	    }
	}
	if (ie < 2) {
	    for (i = minmn - 1; i >= 1; i--) {
		work[i + 1] = work[i + ie - 1];
	    }
	}
    }
//Undo scaling if necessary
    if (iscl == 1) {
	if (anrm > bignum) {
	    Rlascl("G", 0, 0, bignum, anrm, minmn, 1, &s[1], minmn, &ierr);
	}
	if (*info != 0 && anrm > bignum) {
	    Rlascl("G", 0, 0, bignum, anrm, minmn - 1, 1, &work[2], minmn, &ierr);
	}
	if (anrm < smlnum) {
	    Rlascl("G", 0, 0, smlnum, anrm, minmn, 1, &s[1], minmn, &ierr);
	}
	if (*info != 0 && anrm < smlnum) {
	    Rlascl("G", 0, 0, smlnum, anrm, minmn - 1, 1, &work[2], minmn, &ierr);
	}
    }
//Return optimal workspace in WORK(1)
    work[1] = maxwrk;
    return;
}
