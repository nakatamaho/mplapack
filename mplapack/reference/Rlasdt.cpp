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

void Rlasdt(INTEGER const n, INTEGER &lvl, INTEGER &nd, INTEGER *inode, INTEGER *ndiml, INTEGER *ndimr, INTEGER const msub) {
    //
    //     Find the number of levels on the tree.
    //
    INTEGER maxn = max((INTEGER)1, n);
    const REAL two = 2.0e+0;
    REAL temp = log(castREAL(maxn) / castREAL(msub + (INTEGER)1)) / log(two);
    lvl = castINTEGER(temp) + 1;
    //
    INTEGER i = n / 2;
    inode[1 - 1] = i + 1;
    ndiml[1 - 1] = i;
    ndimr[1 - 1] = n - i - 1;
    INTEGER il = 0;
    INTEGER ir = 1;
    INTEGER llst = 1;
    INTEGER nlvl = 0;
    INTEGER ncrnt = 0;
    for (nlvl = 1; nlvl <= lvl - 1; nlvl = nlvl + 1) {
        //
        //        Constructing the tree at (NLVL+1)-st level. The number of
        //        nodes created on this level is LLST * 2.
        //
        for (i = 0; i <= llst - 1; i = i + 1) {
            il += 2;
            ir += 2;
            ncrnt = llst + i;
            ndiml[il - 1] = ndiml[ncrnt - 1] / 2;
            ndimr[il - 1] = ndiml[ncrnt - 1] - ndiml[il - 1] - 1;
            inode[il - 1] = inode[ncrnt - 1] - ndimr[il - 1] - 1;
            ndiml[ir - 1] = ndimr[ncrnt - 1] / 2;
            ndimr[ir - 1] = ndimr[ncrnt - 1] - ndiml[ir - 1] - 1;
            inode[ir - 1] = inode[ncrnt - 1] + ndiml[ir - 1] + 1;
        }
        llst = llst * 2;
    }
    nd = llst * 2 - 1;
    //
    //     End of Rlasdt
    //
}
