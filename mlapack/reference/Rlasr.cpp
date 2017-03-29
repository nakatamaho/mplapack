/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasr.cpp,v 1.8 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rlasr(const char *side, const char *pivot, const char *direct, INTEGER m, INTEGER n, REAL * c, REAL * s, REAL * A, INTEGER lda)
{
    INTEGER i, j, info;
    REAL ctemp, stemp, temp;
    REAL Zero = 0.0, One = 1.0;

//Test the input parameters
    info = 0;
    if (!(Mlsame(side, "L") || Mlsame(side, "R"))) {
	info = 1;
    } else if (!(Mlsame(pivot, "V") || Mlsame(pivot, "T")
		 || Mlsame(pivot, "B"))) {
	info = 2;
    } else if (!(Mlsame(direct, "F") || Mlsame(direct, "B"))) {
	info = 3;
    } else if (m < 0) {
	info = 4;
    } else if (n < 0) {
	info = 5;
    } else if (lda < max((INTEGER) 1, m)) {
	info = 9;
    }
    if (info != 0) {
	Mxerbla("Rlasr ", info);
	return;
    }
//Quick return if possible
    if (m == 0 || n == 0) {
	return;
    }
    if (Mlsame(side, "L")) {
//Form  P * A
	if (Mlsame(pivot, "V")) {
	    if (Mlsame(direct, "F")) {
		for (j = 0; j < m - 1; j++) {
		    ctemp = c[j];
		    stemp = s[j];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < n; i++) {
			    temp = A[(j + 1) + i * lda];
			    A[(j + 1) + i * lda] = ctemp * temp - stemp * A[j + i * lda];
			    A[j + i * lda] = stemp * temp + ctemp * A[j + i * lda];
			}
		    }
		}
	    } else if (Mlsame(direct, "B")) {
		for (j = m - 2; j >= 0; j--) {
		    ctemp = c[j];
		    stemp = s[j];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < n; i++) {
			    temp = A[(j + 1) + i * lda];
			    A[(j + 1) + i * lda] = ctemp * temp - stemp * A[j + i * lda];
			    A[j + i * lda] = stemp * temp + ctemp * A[j + i * lda];
			}
		    }
		}
	    }
	} else if (Mlsame(pivot, "T")) {
	    if (Mlsame(direct, "F")) {
		for (j = 1; j < m; j++) {
		    ctemp = c[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < n; i++) {
			    temp = A[j + i * lda];
			    A[j + i * lda] = ctemp * temp - stemp * A[i * lda];
			    A[i * lda] = stemp * temp + ctemp * A[i * lda];
			}
		    }
		}
	    } else if (Mlsame(direct, "B")) {
		for (j = m - 1; j >= 1; j--) {
		    ctemp = c[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < n; i++) {
			    temp = A[j + i * lda];
			    A[j + i * lda] = ctemp * temp - stemp * A[i * lda];
			    A[i * lda] = stemp * temp + ctemp * A[i * lda];
			}
		    }
		}
	    }
	} else if (Mlsame(pivot, "B")) {
	    if (Mlsame(direct, "F")) {
		for (j = 0; j < m - 1; j++) {
		    ctemp = c[j];
		    stemp = s[j];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < n; i++) {
			    temp = A[j + i * lda];
			    A[j + i * lda] = stemp * A[(m - 1) + i * lda]
				+ ctemp * temp;
			    A[(m - 1) + i * lda] = ctemp * A[(m - 1) + i * lda] - stemp * temp;
			}
		    }
		}
	    } else if (Mlsame(direct, "B")) {
		for (j = m - 2; j >= 0; j--) {
		    ctemp = c[j];
		    stemp = s[j];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < n; i++) {
			    temp = A[j + i * lda];
			    A[j + i * lda] = stemp * A[(m - 1) + i * lda]
				+ ctemp * temp;
			    A[(m - 1) + i * lda] = ctemp * A[(m - 1) + i * lda] - stemp * temp;
			}
		    }
		}
	    }
	}
    } else if (Mlsame(side, "R")) {
//Form A * P'
	if (Mlsame(pivot, "V")) {
	    if (Mlsame(direct, "F")) {
		for (j = 0; j < n - 1; j++) {
		    ctemp = c[j];
		    stemp = s[j];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < m; i++) {
			    temp = A[i + (j + 1) * lda];
			    A[i + (j + 1) * lda] = ctemp * temp - stemp * A[i + j * lda];
			    A[i + j * lda] = stemp * temp + ctemp * A[i + j * lda];
			}
		    }
		}
	    } else if (Mlsame(direct, "B")) {
		for (j = n - 2; j >= 0; j--) {
		    ctemp = c[j];
		    stemp = s[j];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < m; i++) {
			    temp = A[i + (j + 1) * lda];
			    A[i + (j + 1) * lda] = ctemp * temp - stemp * A[i + j * lda];
			    A[i + j * lda] = stemp * temp + ctemp * A[i + j * lda];
			}
		    }
		}
	    }
	} else if (Mlsame(pivot, "T")) {
	    if (Mlsame(direct, "F")) {
		for (j = 1; j < n; j++) {
		    ctemp = c[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < m; i++) {
			    temp = A[i + j * lda];
			    A[i + j * lda] = ctemp * temp - stemp * A[i];
			    A[i] = stemp * temp + ctemp * A[i];
			}
		    }
		}
	    } else if (Mlsame(direct, "B")) {
		for (j = n - 1; j >= 1; j--) {
		    ctemp = c[j - 1];
		    stemp = s[j - 1];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < m; i++) {
			    temp = A[i + j * lda];
			    A[i + j * lda] = ctemp * temp - stemp * A[i];
			    A[i] = stemp * temp + ctemp * A[i];
			}
		    }
		}
	    }
	} else if (Mlsame(pivot, "B")) {
	    if (Mlsame(direct, "F")) {
		for (j = 0; j < n - 1; j++) {
		    ctemp = c[j];
		    stemp = s[j];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < m; i++) {
			    temp = A[i + j * lda];
			    A[i + j * lda] = stemp * A[i + (n - 1) * lda]
				+ ctemp * temp;
			    A[i + (n - 1) * lda] = ctemp * A[i + (n - 1) * lda] - stemp * temp;
			}
		    }
		}
	    } else if (Mlsame(direct, "B")) {
		for (j = n - 2; j >= 0; j--) {
		    ctemp = c[j];
		    stemp = s[j];
		    if (ctemp != One || stemp != Zero) {
			for (i = 0; i < m; i++) {
			    temp = A[i + j * lda];
			    A[i + j * lda] = stemp * A[i + (n - 1) * lda]
				+ ctemp * temp;
			    A[i + (n - 1) * lda] = ctemp * A[i + (n - 1) * lda] - stemp * temp;
			}
		    }
		}
	    }
	}
    }
    return;
}
