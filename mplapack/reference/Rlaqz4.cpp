/*
 * Copyright (c) 2008-2025
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

// Derived from LAPACK routine DLAQZ4.
// Original LAPACK authors:
//   Thijs Steel, KU Leuven

#include <mpblas.h>
#include <mplapack.h>

void Rlaqz4(bool const ilschur, bool const ilq, bool const ilz, INTEGER const n, INTEGER const ilo, INTEGER const ihi, INTEGER const nshifts, INTEGER const nblock_desired, REAL *sr, REAL *si, REAL *ss, REAL *a, INTEGER const lda, REAL *b, INTEGER const ldb, REAL *q, INTEGER const ldq, REAL *z, INTEGER const ldz, REAL *qc, INTEGER const ldqc, REAL *zc, INTEGER const ldzc, REAL *work, INTEGER const lwork, INTEGER &info) {
    //
    // Function arguments
    //
    // Parameters
    //
    // Local scalars
    //
    // External functions
    //
    info = 0;
    if (nblock_desired < nshifts + 1) {
        info = -8;
    }
    if (lwork == -1) {
        // workspace query, quick return
        work[1 - 1] = n * nblock_desired;
        return;
    } else if (lwork < n * nblock_desired) {
        info = -25;
    }
    //
    if (info != 0) {
        Mxerbla("Rlaqz4", -info);
        return;
    }
    //
    // Executable statements
    //
    if (nshifts < 2) {
        return;
    }
    //
    if (ilo >= ihi) {
        return;
    }
    //
    INTEGER istartm = 0;
    INTEGER istopm = 0;
    if (ilschur) {
        istartm = 1;
        istopm = n;
    } else {
        istartm = ilo;
        istopm = ihi;
    }
    //
    // Shuffle shifts into pairs of real shifts and pairs
    // of complex conjugate shifts assuming complex
    // conjugate shifts are already adjacent to one
    // another
    //
    INTEGER i = 0;
    REAL swap = 0.0;
    for (i = 1; i <= nshifts - 2; i = i + 2) {
        if (si[i - 1] != -si[(i + 1) - 1]) {
            //
            swap = sr[i - 1];
            sr[i - 1] = sr[(i + 1) - 1];
            sr[(i + 1) - 1] = sr[(i + 2) - 1];
            sr[(i + 2) - 1] = swap;
            //
            swap = si[i - 1];
            si[i - 1] = si[(i + 1) - 1];
            si[(i + 1) - 1] = si[(i + 2) - 1];
            si[(i + 2) - 1] = swap;
            //
            swap = ss[i - 1];
            ss[i - 1] = ss[(i + 1) - 1];
            ss[(i + 1) - 1] = ss[(i + 2) - 1];
            ss[(i + 2) - 1] = swap;
        }
    }
    //
    // NSHFTS is supposed to be even, but if it is odd,
    // then simply reduce it by one.  The shuffle above
    // ensures that the dropped shift is real and that
    // the remaining shifts are paired.
    //
    INTEGER ns = nshifts - mod(nshifts, 2);
    INTEGER npos = max(nblock_desired - ns, (INTEGER)1);
    //
    // The following block introduces the shifts and chases
    // them down one by one just enough to make space for
    // the other shifts. The near-the-diagonal block is
    // of size (ns+1) x ns.
    //
    const REAL zero = 0.0;
    const REAL one = 1.0;
    Rlaset("FULL", ns + 1, ns + 1, zero, one, qc, ldqc);
    Rlaset("FULL", ns, ns, zero, one, zc, ldzc);
    //
    REAL v[3];
    REAL temp = 0.0;
    REAL c1 = 0.0;
    REAL s1 = 0.0;
    REAL c2 = 0.0;
    REAL s2 = 0.0;
    INTEGER j = 0;
    for (i = 1; i <= ns; i = i + 2) {
        // Introduce the shift
        Rlaqz1(&a[(ilo - 1) + (ilo - 1) * lda], lda, &b[(ilo - 1) + (ilo - 1) * ldb], ldb, sr[i - 1], sr[(i + 1) - 1], si[i - 1], ss[i - 1], ss[(i + 1) - 1], v);
        //
        temp = v[2 - 1];
        Rlartg(temp, v[3 - 1], c1, s1, v[2 - 1]);
        Rlartg(v[1 - 1], v[2 - 1], c2, s2, temp);
        //
        Rrot(ns, &a[((ilo + 1) - 1) + (ilo - 1) * lda], lda, &a[((ilo + 2) - 1) + (ilo - 1) * lda], lda, c1, s1);
        Rrot(ns, &a[(ilo - 1) + (ilo - 1) * lda], lda, &a[((ilo + 1) - 1) + (ilo - 1) * lda], lda, c2, s2);
        Rrot(ns, &b[((ilo + 1) - 1) + (ilo - 1) * ldb], ldb, &b[((ilo + 2) - 1) + (ilo - 1) * ldb], ldb, c1, s1);
        Rrot(ns, &b[(ilo - 1) + (ilo - 1) * ldb], ldb, &b[((ilo + 1) - 1) + (ilo - 1) * ldb], ldb, c2, s2);
        Rrot(ns + 1, &qc[(2 - 1) * ldqc], 1, &qc[(3 - 1) * ldqc], 1, c1, s1);
        Rrot(ns + 1, &qc[0], 1, &qc[(2 - 1) * ldqc], 1, c2, s2);
        //
        // Chase the shift down
        for (j = 1; j <= ns - 1 - i; j = j + 1) {
            //
            Rlaqz2(true, true, j, 1, ns, ihi - ilo + 1, &a[(ilo - 1) + (ilo - 1) * lda], lda, &b[(ilo - 1) + (ilo - 1) * ldb], ldb, ns + 1, 1, qc, ldqc, ns, 1, zc, ldzc);
            //
        }
        //
    }
    //
    // Update the rest of the pencil
    //
    // Update A(ilo:ilo+ns,ilo+ns:istopm) and B(ilo:ilo+ns,ilo+ns:istopm)
    // from the left with Qc(1:ns+1,1:ns+1)'
    INTEGER sheight = ns + 1;
    INTEGER swidth = istopm - (ilo + ns) + 1;
    if (swidth > 0) {
        Rgemm("T", "N", sheight, swidth, sheight, one, qc, ldqc, &a[(ilo - 1) + ((ilo + ns) - 1) * lda], lda, zero, work, sheight);
        Rlacpy("ALL", sheight, swidth, work, sheight, &a[(ilo - 1) + ((ilo + ns) - 1) * lda], lda);
        Rgemm("T", "N", sheight, swidth, sheight, one, qc, ldqc, &b[(ilo - 1) + ((ilo + ns) - 1) * ldb], ldb, zero, work, sheight);
        Rlacpy("ALL", sheight, swidth, work, sheight, &b[(ilo - 1) + ((ilo + ns) - 1) * ldb], ldb);
    }
    if (ilq) {
        Rgemm("N", "N", n, sheight, sheight, one, &q[(ilo - 1) * ldq], ldq, qc, ldqc, zero, work, n);
        Rlacpy("ALL", n, sheight, work, n, &q[(ilo - 1) * ldq], ldq);
    }
    //
    // Update A(istartm:ilo-1,ilo:ilo+ns-1) and B(istartm:ilo-1,ilo:ilo+ns-1)
    // from the right with Zc(1:ns,1:ns)
    sheight = ilo - 1 - istartm + 1;
    swidth = ns;
    if (sheight > 0) {
        Rgemm("N", "N", sheight, swidth, swidth, one, &a[(istartm - 1) + (ilo - 1) * lda], lda, zc, ldzc, zero, work, sheight);
        Rlacpy("ALL", sheight, swidth, work, sheight, &a[(istartm - 1) + (ilo - 1) * lda], lda);
        Rgemm("N", "N", sheight, swidth, swidth, one, &b[(istartm - 1) + (ilo - 1) * ldb], ldb, zc, ldzc, zero, work, sheight);
        Rlacpy("ALL", sheight, swidth, work, sheight, &b[(istartm - 1) + (ilo - 1) * ldb], ldb);
    }
    if (ilz) {
        Rgemm("N", "N", n, swidth, swidth, one, &z[(ilo - 1) * ldz], ldz, zc, ldzc, zero, work, n);
        Rlacpy("ALL", n, swidth, work, n, &z[(ilo - 1) * ldz], ldz);
    }
    //
    // The following block chases the shifts down to the bottom
    // right block. If possible, a shift is moved down npos
    // positions at a time
    //
    INTEGER k = ilo;
    INTEGER np = 0;
    INTEGER nblock = 0;
    INTEGER istartb = 0;
    INTEGER istopb = 0;
    while (k < ihi - ns) {
        np = min(ihi - ns - k, npos);
        // Size of the near-the-diagonal block
        nblock = ns + np;
        // istartb points to the first row we will be updating
        istartb = k + 1;
        // istopb points to the last column we will be updating
        istopb = k + nblock - 1;
        //
        Rlaset("FULL", ns + np, ns + np, zero, one, qc, ldqc);
        Rlaset("FULL", ns + np, ns + np, zero, one, zc, ldzc);
        //
        // Near the diagonal shift chase
        for (i = ns - 1; i >= 0; i = i - 2) {
            for (j = 0; j <= np - 1; j = j + 1) {
                // Move down the block with index k+i+j-1, updating
                // the (ns+np x ns+np) block:
                // (k:k+ns+np,k:k+ns+np-1)
                Rlaqz2(true, true, k + i + j - 1, istartb, istopb, ihi, a, lda, b, ldb, nblock, k + 1, qc, ldqc, nblock, k, zc, ldzc);
            }
        }
        //
        // Update rest of the pencil
        //
        // Update A(k+1:k+ns+np, k+ns+np:istopm) and
        // B(k+1:k+ns+np, k+ns+np:istopm)
        // from the left with Qc(1:ns+np,1:ns+np)'
        sheight = ns + np;
        swidth = istopm - (k + ns + np) + 1;
        if (swidth > 0) {
            Rgemm("T", "N", sheight, swidth, sheight, one, qc, ldqc, &a[((k + 1) - 1) + ((k + ns + np) - 1) * lda], lda, zero, work, sheight);
            Rlacpy("ALL", sheight, swidth, work, sheight, &a[((k + 1) - 1) + ((k + ns + np) - 1) * lda], lda);
            Rgemm("T", "N", sheight, swidth, sheight, one, qc, ldqc, &b[((k + 1) - 1) + ((k + ns + np) - 1) * ldb], ldb, zero, work, sheight);
            Rlacpy("ALL", sheight, swidth, work, sheight, &b[((k + 1) - 1) + ((k + ns + np) - 1) * ldb], ldb);
        }
        if (ilq) {
            Rgemm("N", "N", n, nblock, nblock, one, &q[((k + 1) - 1) * ldq], ldq, qc, ldqc, zero, work, n);
            Rlacpy("ALL", n, nblock, work, n, &q[((k + 1) - 1) * ldq], ldq);
        }
        //
        // Update A(istartm:k,k:k+ns+npos-1) and B(istartm:k,k:k+ns+npos-1)
        // from the right with Zc(1:ns+np,1:ns+np)
        sheight = k - istartm + 1;
        swidth = nblock;
        if (sheight > 0) {
            Rgemm("N", "N", sheight, swidth, swidth, one, &a[(istartm - 1) + (k - 1) * lda], lda, zc, ldzc, zero, work, sheight);
            Rlacpy("ALL", sheight, swidth, work, sheight, &a[(istartm - 1) + (k - 1) * lda], lda);
            Rgemm("N", "N", sheight, swidth, swidth, one, &b[(istartm - 1) + (k - 1) * ldb], ldb, zc, ldzc, zero, work, sheight);
            Rlacpy("ALL", sheight, swidth, work, sheight, &b[(istartm - 1) + (k - 1) * ldb], ldb);
        }
        if (ilz) {
            Rgemm("N", "N", n, nblock, nblock, one, &z[(k - 1) * ldz], ldz, zc, ldzc, zero, work, n);
            Rlacpy("ALL", n, nblock, work, n, &z[(k - 1) * ldz], ldz);
        }
        //
        k += np;
        //
    }
    //
    // The following block removes the shifts from the bottom right corner
    // one by one. Updates are initially applied to A(ihi-ns+1:ihi,ihi-ns:ihi).
    //
    Rlaset("FULL", ns, ns, zero, one, qc, ldqc);
    Rlaset("FULL", ns + 1, ns + 1, zero, one, zc, ldzc);
    //
    // istartb points to the first row we will be updating
    istartb = ihi - ns + 1;
    // istopb points to the last column we will be updating
    istopb = ihi;
    //
    INTEGER ishift = 0;
    for (i = 1; i <= ns; i = i + 2) {
        // Chase the shift down to the bottom right corner
        for (ishift = ihi - i - 1; ishift <= ihi - 2; ishift = ishift + 1) {
            Rlaqz2(true, true, ishift, istartb, istopb, ihi, a, lda, b, ldb, ns, ihi - ns + 1, qc, ldqc, ns + 1, ihi - ns, zc, ldzc);
        }
        //
    }
    //
    // Update rest of the pencil
    //
    // Update A(ihi-ns+1:ihi, ihi+1:istopm)
    // from the left with Qc(1:ns,1:ns)'
    sheight = ns;
    swidth = istopm - (ihi + 1) + 1;
    if (swidth > 0) {
        Rgemm("T", "N", sheight, swidth, sheight, one, qc, ldqc, &a[((ihi - ns + 1) - 1) + ((ihi + 1) - 1) * lda], lda, zero, work, sheight);
        Rlacpy("ALL", sheight, swidth, work, sheight, &a[((ihi - ns + 1) - 1) + ((ihi + 1) - 1) * lda], lda);
        Rgemm("T", "N", sheight, swidth, sheight, one, qc, ldqc, &b[((ihi - ns + 1) - 1) + ((ihi + 1) - 1) * ldb], ldb, zero, work, sheight);
        Rlacpy("ALL", sheight, swidth, work, sheight, &b[((ihi - ns + 1) - 1) + ((ihi + 1) - 1) * ldb], ldb);
    }
    if (ilq) {
        Rgemm("N", "N", n, ns, ns, one, &q[((ihi - ns + 1) - 1) * ldq], ldq, qc, ldqc, zero, work, n);
        Rlacpy("ALL", n, ns, work, n, &q[((ihi - ns + 1) - 1) * ldq], ldq);
    }
    //
    // Update A(istartm:ihi-ns,ihi-ns:ihi)
    // from the right with Zc(1:ns+1,1:ns+1)
    sheight = ihi - ns - istartm + 1;
    swidth = ns + 1;
    if (sheight > 0) {
        Rgemm("N", "N", sheight, swidth, swidth, one, &a[(istartm - 1) + ((ihi - ns) - 1) * lda], lda, zc, ldzc, zero, work, sheight);
        Rlacpy("ALL", sheight, swidth, work, sheight, &a[(istartm - 1) + ((ihi - ns) - 1) * lda], lda);
        Rgemm("N", "N", sheight, swidth, swidth, one, &b[(istartm - 1) + ((ihi - ns) - 1) * ldb], ldb, zc, ldzc, zero, work, sheight);
        Rlacpy("ALL", sheight, swidth, work, sheight, &b[(istartm - 1) + ((ihi - ns) - 1) * ldb], ldb);
    }
    if (ilz) {
        Rgemm("N", "N", n, ns + 1, ns + 1, one, &z[((ihi - ns) - 1) * ldz], ldz, zc, ldzc, zero, work, n);
        Rlacpy("ALL", n, ns + 1, work, n, &z[((ihi - ns) - 1) * ldz], ldz);
    }
    //
}
