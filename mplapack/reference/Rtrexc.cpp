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

void Rtrexc(const char *compq, INTEGER const n, REAL *t, INTEGER const ldt, REAL *q, INTEGER const ldq, INTEGER &ifst, INTEGER &ilst, REAL *work, INTEGER &info) {
    bool wantq = false;
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
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    //     Decode and test the input arguments.
    //
    info = 0;
    wantq = Mlsame(compq, "V");
    if (!wantq && !Mlsame(compq, "N")) {
        info = -1;
    } else if (n < 0) {
        info = -2;
    } else if (ldt < max((INTEGER)1, n)) {
        info = -4;
    } else if (ldq < 1 || (wantq && ldq < max((INTEGER)1, n))) {
        info = -6;
    } else if ((ifst < 1 || ifst > n) && (n > 0)) {
        info = -7;
    } else if ((ilst < 1 || ilst > n) && (n > 0)) {
        info = -8;
    }
    if (info != 0) {
        Mxerbla("Rtrexc", -info);
        return;
    }
    //
    //     Quick return if possible
    //
    if (n <= 1) {
        return;
    }
    //
    //     Determine the first row of specified block
    //     and find out it is 1 by 1 or 2 by 2.
    //
    if (ifst > 1) {
        if (t[(ifst - 1) + ((ifst - 1) - 1) * ldt] != zero) {
            ifst = ifst - 1;
        }
    }
    nbf = 1;
    if (ifst < n) {
        if (t[((ifst + 1) - 1) + (ifst - 1) * ldt] != zero) {
            nbf = 2;
        }
    }
    //
    //     Determine the first row of the final block
    //     and find out it is 1 by 1 or 2 by 2.
    //
    if (ilst > 1) {
        if (t[(ilst - 1) + ((ilst - 1) - 1) * ldt] != zero) {
            ilst = ilst - 1;
        }
    }
    nbl = 1;
    if (ilst < n) {
        if (t[((ilst + 1) - 1) + (ilst - 1) * ldt] != zero) {
            nbl = 2;
        }
    }
    //
    if (ifst == ilst) {
        return;
    }
    //
    if (ifst < ilst) {
        //
        //        Update ILST
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
        //        Swap block with next one below
        //
        if (nbf == 1 || nbf == 2) {
            //
            //           Current block either 1 by 1 or 2 by 2
            //
            nbnext = 1;
            if (here + nbf + 1 <= n) {
                if (t[((here + nbf + 1) - 1) + ((here + nbf) - 1) * ldt] != zero) {
                    nbnext = 2;
                }
            }
            Rlaexc(wantq, n, t, ldt, q, ldq, here, nbf, nbnext, work, info);
            if (info != 0) {
                ilst = here;
                return;
            }
            here += nbnext;
            //
            //           Test if 2 by 2 block breaks into two 1 by 1 blocks
            //
            if (nbf == 2) {
                if (t[((here + 1) - 1) + (here - 1) * ldt] == zero) {
                    nbf = 3;
                }
            }
            //
        } else {
            //
            //           Current block consists of two 1 by 1 blocks each of which
            //           must be swapped individually
            //
            nbnext = 1;
            if (here + 3 <= n) {
                if (t[((here + 3) - 1) + ((here + 2) - 1) * ldt] != zero) {
                    nbnext = 2;
                }
            }
            Rlaexc(wantq, n, t, ldt, q, ldq, here + 1, 1, nbnext, work, info);
            if (info != 0) {
                ilst = here;
                return;
            }
            if (nbnext == 1) {
                //
                //              Swap two 1 by 1 blocks, no problems possible
                //
                Rlaexc(wantq, n, t, ldt, q, ldq, here, 1, nbnext, work, info);
                here++;
            } else {
                //
                //              Recompute NBNEXT in case 2 by 2 split
                //
                if (t[((here + 2) - 1) + ((here + 1) - 1) * ldt] == zero) {
                    nbnext = 1;
                }
                if (nbnext == 2) {
                    //
                    //                 2 by 2 Block did not split
                    //
                    Rlaexc(wantq, n, t, ldt, q, ldq, here, 1, nbnext, work, info);
                    if (info != 0) {
                        ilst = here;
                        return;
                    }
                    here += 2;
                } else {
                    //
                    //                 2 by 2 Block did split
                    //
                    Rlaexc(wantq, n, t, ldt, q, ldq, here, 1, 1, work, info);
                    Rlaexc(wantq, n, t, ldt, q, ldq, here + 1, 1, 1, work, info);
                    here += 2;
                }
            }
        }
        if (here < ilst) {
            goto statement_10;
        }
        //
    } else {
        //
        here = ifst;
    statement_20:
        //
        //        Swap block with next one above
        //
        if (nbf == 1 || nbf == 2) {
            //
            //           Current block either 1 by 1 or 2 by 2
            //
            nbnext = 1;
            if (here >= 3) {
                if (t[((here - 1) - 1) + ((here - 2) - 1) * ldt] != zero) {
                    nbnext = 2;
                }
            }
            Rlaexc(wantq, n, t, ldt, q, ldq, here - nbnext, nbnext, nbf, work, info);
            if (info != 0) {
                ilst = here;
                return;
            }
            here = here - nbnext;
            //
            //           Test if 2 by 2 block breaks into two 1 by 1 blocks
            //
            if (nbf == 2) {
                if (t[((here + 1) - 1) + (here - 1) * ldt] == zero) {
                    nbf = 3;
                }
            }
            //
        } else {
            //
            //           Current block consists of two 1 by 1 blocks each of which
            //           must be swapped individually
            //
            nbnext = 1;
            if (here >= 3) {
                if (t[((here - 1) - 1) + ((here - 2) - 1) * ldt] != zero) {
                    nbnext = 2;
                }
            }
            Rlaexc(wantq, n, t, ldt, q, ldq, here - nbnext, nbnext, 1, work, info);
            if (info != 0) {
                ilst = here;
                return;
            }
            if (nbnext == 1) {
                //
                //              Swap two 1 by 1 blocks, no problems possible
                //
                Rlaexc(wantq, n, t, ldt, q, ldq, here, nbnext, 1, work, info);
                here = here - 1;
            } else {
                //
                //              Recompute NBNEXT in case 2 by 2 split
                //
                if (t[(here - 1) + ((here - 1) - 1) * ldt] == zero) {
                    nbnext = 1;
                }
                if (nbnext == 2) {
                    //
                    //                 2 by 2 Block did not split
                    //
                    Rlaexc(wantq, n, t, ldt, q, ldq, here - 1, 2, 1, work, info);
                    if (info != 0) {
                        ilst = here;
                        return;
                    }
                    here = here - 2;
                } else {
                    //
                    //                 2 by 2 Block did split
                    //
                    Rlaexc(wantq, n, t, ldt, q, ldq, here, 1, 1, work, info);
                    Rlaexc(wantq, n, t, ldt, q, ldq, here - 1, 1, 1, work, info);
                    here = here - 2;
                }
            }
        }
        if (here > ilst) {
            goto statement_20;
        }
    }
    ilst = here;
    //
    //     End of Rtrexc
    //
}
