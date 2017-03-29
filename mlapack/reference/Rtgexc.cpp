/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rtgexc.cpp,v 1.4 2010/08/07 04:48:33 nakatamaho Exp $ 
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

void Rtgexc(LOGICAL wantq, LOGICAL wantz, INTEGER n,
	    REAL * A, INTEGER lda, REAL * B, INTEGER ldb, REAL * q, INTEGER ldq, REAL * z, INTEGER ldz, INTEGER * ifst, INTEGER * ilst, REAL * work, INTEGER lwork, INTEGER * info)
{
    INTEGER nbf, nbl, here, lwmin;
    INTEGER nbnext;
    LOGICAL lquery;
    REAL Zero = 0.0;

//Decode and test input arguments.
    *info = 0;
    lquery = lwork == -1;
    if (n < 0) {
	*info = -3;
    } else if (lda < max((INTEGER) 1, n)) {
	*info = -5;
    } else if (ldb < max((INTEGER) 1, n)) {
	*info = -7;
    } else if (ldq < 1 || (wantq && ldq < max((INTEGER) 1, n))) {
	*info = -9;
    } else if (ldz < 1 || (wantz && ldz < max((INTEGER) 1, n))) {
	*info = -11;
    } else if (*ifst < 1 || *ifst > n) {
	*info = -12;
    } else if (*ilst < 1 || *ilst > n) {
	*info = -13;
    }
    if (*info == 0) {
	if (n <= 1) {
	    lwmin = 1;
	} else {
	    lwmin = (n << 2) + 16;
	}
	work[1] = lwmin;
	if (lwork < lwmin && !lquery) {
	    *info = -15;
	}
    }
    if (*info != 0) {
	Mxerbla("Rtgexc", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n <= 1) {
	return;
    }
//Determine the first row of the specified block and find out
//if it is 1-by-1 or 2-by-2
    if (*ifst > 1) {
	if (A[*ifst + (*ifst - 1) * lda] != Zero) {
	    --(*ifst);
	}
    }
    nbf = 1;
    if (*ifst < n) {
	if (A[*ifst + 1 + *ifst * lda] != Zero) {
	    nbf = 2;
	}
    }
//Determine the first row of the final block
//and find out if it is 1-by-1 or 2-by-2
    if (*ilst > 1) {
	if (A[*ilst + (*ilst - 1) * lda] != Zero) {
	    --(*ilst);
	}
    }
    nbl = 0;
    if (*ilst < n) {
	if (A[*ilst + 1 + *ilst * lda] != Zero) {
	    nbl = 2;
	}
    }
    if (*ifst == *ilst) {
	return;
    }
    if (*ifst < *ilst) {
//Update ILST.
	if (nbf == 2 && nbl == 1) {
	    --(*ilst);
	}
	if (nbf == 1 && nbl == 2) {
	    ++(*ilst);
	}
	here = *ifst;
      L10:
//Swap with next one below.
	if (nbf == 1 || nbf == 2) {
//Current block either 1-by-1 or 2-by-2
	    nbnext = 1;
	    if (here + nbf + 1 <= n) {
		if (A[here + nbf + 1 + (here + nbf) * lda] != Zero) {
		    nbnext = 2;
		}
	    }
	    Rtgex2(wantq, wantz, n, &A[0], lda, &B[0], ldb, &q[0], ldq, &z[0], ldz, here, nbf, nbnext, &work[0], lwork, info);
	    if (*info != 0) {
		*ilst = here;
		return;
	    }
	    here = here + nbnext;
//Test if 2-by-2 block breaks into two 1-by-1 blocks.
	    if (nbf == 2) {
		if (A[here + 1 + here * lda] == Zero) {
		    nbf = 3;
		}
	    }
	} else {
//Current block consists of two 1-by-1 blocks, each of which
//must be swapped individually.
	    nbnext = 1;
	    if (here + 3 <= n) {
		if (A[here + 3 + (here + 2) * lda] != Zero) {
		    nbnext = 2;
		}
	    }
	    Rtgex2(wantq, wantz, n, &A[0], lda, &B[0], ldb, &q[0], ldq, &z[0], ldz, here + 1, 1, nbnext, &work[0], lwork, info);
	    if (*info != 0) {
		*ilst = here;
		return;
	    }
	    if (nbnext == 1) {
//Swap two 1-by-1 blocks.
		Rtgex2(wantq, wantz, n, &A[0], lda, &B[0], ldb, &q[0], ldq, &z[0], ldz, here, 1, 1, &work[0], lwork, info);
		if (*info != 0) {
		    *ilst = here;
		    return;
		}
		++here;
	    } else {
//Recompute NBNEXT in case of 2-by-2 split.
		if (A[here + 2 + (here + 1) * lda] == Zero) {
		    nbnext = 1;
		}
		if (nbnext == 2) {
//2-by-2 block did not split.
		    Rtgex2(wantq, wantz, n, &A[0], lda, &B[0], ldb, &q[0], ldq, &z[0], ldz, here, 1, nbnext, &work[0], lwork, info);
		    if (*info != 0) {
			*ilst = here;
			return;
		    }
		    here += 2;
		} else {
//2-by-2 block did split.
		    Rtgex2(wantq, wantz, n, &A[0], lda, &B[0], ldb, &q[0], ldq, &z[0], ldz, here, 1, 1, &work[0], lwork, info);
		    if (*info != 0) {
			*ilst = here;
			return;
		    }
		    ++here;
		    Rtgex2(wantq, wantz, n, &A[0], lda, &B[0], ldb, &q[0], ldq, &z[0], ldz, here, 1, 1, &work[0], lwork, info);
		    if (*info != 0) {
			*ilst = here;
			return;
		    }
		    ++here;
		}
	    }
	}
	if (here < *ilst) {
	    goto L10;
	}
    } else {
	here = *ifst;
      L20:
//Swap with next one below.
	if (nbf == 1 || nbf == 2) {
//Current block either 1-by-1 or 2-by-2
	    nbnext = 1;
	    if (here >= 3) {
		if (A[here - 1 + (here - 2) * lda] != Zero) {
		    nbnext = 2;
		}
	    }
	    Rtgex2(wantq, wantz, n, &A[0], lda, &B[0], ldb, &q[0], ldq, &z[0], ldz, here - nbnext, nbnext, nbf, &work[0], lwork, info);
	    if (*info != 0) {
		*ilst = here;
		return;
	    }
	    here -= nbnext;
//Test if 2-by-2 block breaks into two 1-by-1 blocks.
	    if (nbf == 2) {
		if (A[here + 1 + here * lda] == Zero) {
		    nbf = 3;
		}
	    }
	} else {
//Current block consists of two 1-by-1 blocks, each of which
//must be swapped individually.
	    nbnext = 1;
	    if (here >= 3) {
		if (A[here - 1 + (here - 2) * lda] != Zero) {
		    nbnext = 2;
		}
	    }
	    Rtgex2(wantq, wantz, n, &A[0], lda, &B[0], ldb, &q[0], ldq, &z[0], ldz, here - nbnext, nbnext, 1, &work[0], lwork, info);
	    if (*info != 0) {
		*ilst = here;
		return;
	    }
	    if (nbnext == 1) {
//Swap two 1-by-1 blocks.
		Rtgex2(wantq, wantz, n, &A[0], lda, &B[0], ldb, &q[0], ldq, &z[0], ldz, here, nbnext, 1, &work[0], lwork, info);
		if (*info != 0) {
		    *ilst = here;
		    return;
		}
		--here;
	    } else {
//Recompute NBNEXT in case of 2-by-2 split.
		if (A[here + (here - 1) * lda] == Zero) {
		    nbnext = 1;
		}
		if (nbnext == 2) {
//2-by-2 block did not split.
		    Rtgex2(wantq, wantz, n, &A[0], lda, &B[0], ldb, &q[0], ldq, &z[0], ldz, here - 1, 2, 1, &work[0], lwork, info);
		    if (*info != 0) {
			*ilst = here;
			return;
		    }
		    here = here + -2;
		} else {
//2-by-2 block did split.
		    Rtgex2(wantq, wantz, n, &A[0], lda, &B[0], ldb, &q[0], ldq, &z[0], ldz, here, 1, 1, &work[0], lwork, info);
		    if (*info != 0) {
			*ilst = here;
			return;
		    }
		    --here;
		    Rtgex2(wantq, wantz, n, &A[0], lda, &B[0], ldb, &q[0], ldq, &z[0], ldz, here, 1, 1, &work[0], lwork, info);
		    if (*info != 0) {
			*ilst = here;
			return;
		    }
		    --here;
		}
	    }
	}
	if (here > *ilst) {
	    goto L20;
	}
    }
    *ilst = here;
    work[1] = lwmin;
    return;
}
