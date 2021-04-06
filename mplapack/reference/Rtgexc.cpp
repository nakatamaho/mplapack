/*
 * Copyright (c) 2008-2021
 *      Nakata, Maho
 *      All rights reserved.
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

#include <mpblas.h>
#include <mplapack.h>

void Rtgexc(bool const wantq, bool const wantz, INTEGER const n, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *q, INTEGER const ldq, REAL *z, INTEGER const ldz, INTEGER &ifst, INTEGER &ilst, REAL *work, INTEGER const lwork, INTEGER &info) {
    bool lquery = false;
    INTEGER lwmin = 0;
    const REAL zero = 0.0;
    INTEGER nbf = 0;
    INTEGER nbl = 0;
    INTEGER here = 0;
    INTEGER nbnext = 0;
    //
    //  -- LAPACK computational routine --
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
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Decode and test input arguments.
    //
    info = 0;
    lquery = (lwork == -1);
    if (n < 0) {
        info = -3;
    } else if (lda < max((INTEGER)1, n)) {
        info = -5;
    } else if (ldb < max((INTEGER)1, n)) {
        info = -7;
    } else if (ldq < 1 || wantq && (ldq < max((INTEGER)1, n))) {
        info = -9;
    } else if (ldz < 1 || wantz && (ldz < max((INTEGER)1, n))) {
        info = -11;
    } else if (ifst < 1 || ifst > n) {
        info = -12;
    } else if (ilst < 1 || ilst > n) {
        info = -13;
    }
    //
    if (info == 0) {
        if (n <= 1) {
            lwmin = 1;
        } else {
            lwmin = 4 * n + 16;
        }
        work[1 - 1] = lwmin;
        //
        if (lwork < lwmin && !lquery) {
            info = -15;
        }
    }
    //
    if (info != 0) {
        Mxerbla("Rtgexc", -info);
        return;
    } else if (lquery) {
        return;
    }
    //
    //     Quick return if possible
    //
    if (n <= 1) {
        return;
    }
    //
    //     Determine the first row of the specified block and find out
    //     if it is 1-by-1 or 2-by-2.
    //
    if (ifst > 1) {
        if (a[(ifst - 1) + ((ifst - 1) - 1) * lda] != zero) {
            ifst = ifst - 1;
        }
    }
    nbf = 1;
    if (ifst < n) {
        if (a[((ifst + 1) - 1) + (ifst - 1) * lda] != zero) {
            nbf = 2;
        }
    }
    //
    //     Determine the first row of the final block
    //     and find out if it is 1-by-1 or 2-by-2.
    //
    if (ilst > 1) {
        if (a[(ilst - 1) + ((ilst - 1) - 1) * lda] != zero) {
            ilst = ilst - 1;
        }
    }
    nbl = 1;
    if (ilst < n) {
        if (a[((ilst + 1) - 1) + (ilst - 1) * lda] != zero) {
            nbl = 2;
        }
    }
    if (ifst == ilst) {
        return;
    }
    //
    if (ifst < ilst) {
        //
        //        Update ILST.
        //
        if (nbf == 2 && nbl == 1) {
            ilst = ilst - 1;
        }
        if (nbf == 1 && nbl == 2) {
            ilst++;
        }
        //
        here = ifst;
    //
    statement_10:
        //
        //        Swap with next one below.
        //
        if (nbf == 1 || nbf == 2) {
            //
            //           Current block either 1-by-1 or 2-by-2.
            //
            nbnext = 1;
            if (here + nbf + 1 <= n) {
                if (a[((here + nbf + 1) - 1) + ((here + nbf) - 1) * lda] != zero) {
                    nbnext = 2;
                }
            }
            Rtgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here, nbf, nbnext, work, lwork, info);
            if (info != 0) {
                ilst = here;
                return;
            }
            here += nbnext;
            //
            //           Test if 2-by-2 block breaks INTEGERo two 1-by-1 blocks.
            //
            if (nbf == 2) {
                if (a[((here + 1) - 1) + (here - 1) * lda] == zero) {
                    nbf = 3;
                }
            }
            //
        } else {
            //
            //           Current block consists of two 1-by-1 blocks, each of which
            //           must be swapped individually.
            //
            nbnext = 1;
            if (here + 3 <= n) {
                if (a[((here + 3) - 1) + ((here + 2) - 1) * lda] != zero) {
                    nbnext = 2;
                }
            }
            Rtgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here + 1, 1, nbnext, work, lwork, info);
            if (info != 0) {
                ilst = here;
                return;
            }
            if (nbnext == 1) {
                //
                //              Swap two 1-by-1 blocks.
                //
                Rtgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here, 1, 1, work, lwork, info);
                if (info != 0) {
                    ilst = here;
                    return;
                }
                here++;
                //
            } else {
                //
                //              Recompute NBNEXT in case of 2-by-2 split.
                //
                if (a[((here + 2) - 1) + ((here + 1) - 1) * lda] == zero) {
                    nbnext = 1;
                }
                if (nbnext == 2) {
                    //
                    //                 2-by-2 block did not split.
                    //
                    Rtgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here, 1, nbnext, work, lwork, info);
                    if (info != 0) {
                        ilst = here;
                        return;
                    }
                    here += 2;
                } else {
                    //
                    //                 2-by-2 block did split.
                    //
                    Rtgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here, 1, 1, work, lwork, info);
                    if (info != 0) {
                        ilst = here;
                        return;
                    }
                    here++;
                    Rtgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here, 1, 1, work, lwork, info);
                    if (info != 0) {
                        ilst = here;
                        return;
                    }
                    here++;
                }
                //
            }
        }
        if (here < ilst) {
            goto statement_10;
        }
    } else {
        here = ifst;
    //
    statement_20:
        //
        //        Swap with next one below.
        //
        if (nbf == 1 || nbf == 2) {
            //
            //           Current block either 1-by-1 or 2-by-2.
            //
            nbnext = 1;
            if (here >= 3) {
                if (a[((here - 1) - 1) + ((here - 2) - 1) * lda] != zero) {
                    nbnext = 2;
                }
            }
            Rtgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here - nbnext, nbnext, nbf, work, lwork, info);
            if (info != 0) {
                ilst = here;
                return;
            }
            here = here - nbnext;
            //
            //           Test if 2-by-2 block breaks INTEGERo two 1-by-1 blocks.
            //
            if (nbf == 2) {
                if (a[((here + 1) - 1) + (here - 1) * lda] == zero) {
                    nbf = 3;
                }
            }
            //
        } else {
            //
            //           Current block consists of two 1-by-1 blocks, each of which
            //           must be swapped individually.
            //
            nbnext = 1;
            if (here >= 3) {
                if (a[((here - 1) - 1) + ((here - 2) - 1) * lda] != zero) {
                    nbnext = 2;
                }
            }
            Rtgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here - nbnext, nbnext, 1, work, lwork, info);
            if (info != 0) {
                ilst = here;
                return;
            }
            if (nbnext == 1) {
                //
                //              Swap two 1-by-1 blocks.
                //
                Rtgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here, nbnext, 1, work, lwork, info);
                if (info != 0) {
                    ilst = here;
                    return;
                }
                here = here - 1;
            } else {
                //
                //             Recompute NBNEXT in case of 2-by-2 split.
                //
                if (a[(here - 1) + ((here - 1) - 1) * lda] == zero) {
                    nbnext = 1;
                }
                if (nbnext == 2) {
                    //
                    //                 2-by-2 block did not split.
                    //
                    Rtgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here - 1, 2, 1, work, lwork, info);
                    if (info != 0) {
                        ilst = here;
                        return;
                    }
                    here = here - 2;
                } else {
                    //
                    //                 2-by-2 block did split.
                    //
                    Rtgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here, 1, 1, work, lwork, info);
                    if (info != 0) {
                        ilst = here;
                        return;
                    }
                    here = here - 1;
                    Rtgex2(wantq, wantz, n, a, lda, b, ldb, q, ldq, z, ldz, here, 1, 1, work, lwork, info);
                    if (info != 0) {
                        ilst = here;
                        return;
                    }
                    here = here - 1;
                }
            }
        }
        if (here > ilst) {
            goto statement_20;
        }
    }
    ilst = here;
    work[1 - 1] = lwmin;
    //
    //     End of Rtgexc
    //
}
