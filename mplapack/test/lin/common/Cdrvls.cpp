/*
 * Copyright (c) 2021
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

#include <fem.hpp> // Fortran EMulation library of fable module
using namespace fem::major_types;
using fem::common;

#include <mplapack_matgen.h>
#include <mplapack_lin.h>

#include <mplapack_debug.h>

void Cdrvls(bool *dotype, INTEGER const nm, INTEGER *mval, INTEGER const nn, INTEGER *nval, INTEGER const nns, INTEGER *nsval, INTEGER const nnb, INTEGER *nbval, INTEGER *nxval, REAL const thresh, bool const tsterr, COMPLEX *a, COMPLEX *copya, COMPLEX *b, COMPLEX *copyb, COMPLEX *c, REAL *s, REAL *copys, INTEGER const nout) {
    common cmn;
    common_write write(cmn);
    //
    INTEGER iseedy[] = {1988, 1989, 1990, 1991};
    char path[3];
    char buf[1024];
    INTEGER nrun = 0;
    INTEGER nfail = 0;
    INTEGER nerrs = 0;
    INTEGER i = 0;
    INTEGER iseed[4];
    REAL eps = 0.0;
    REAL rcond = 0.0;
    const INTEGER smlsiz = 25;
    const REAL zero = 0.0;
    INTEGER nmax = 0;
    INTEGER mmax = 0;
    INTEGER nsmax = 0;
    INTEGER m = 0;
    INTEGER n = 0;
    INTEGER nrhs = 0;
    INTEGER mnmin = 0;
    INTEGER lwork = 0;
    INTEGER lrwork = 0;
    INTEGER liwork = 0;
    INTEGER im = 0;
    INTEGER lda = 0;
    INTEGER in = 0;
    INTEGER ldb = 0;
    INTEGER ins = 0;
    INTEGER irank = 0;
    INTEGER iscale = 0;
    INTEGER itype = 0;
    INTEGER itran = 0;
    char trans[1];
    COMPLEX wq[1];
    INTEGER info = 0;
    INTEGER lwork_Cgels = 0;
    INTEGER lwork_Cgetsls = 0;
    INTEGER iwq[1];
    INTEGER crank = 0;
    REAL rwq[1];
    INTEGER lwork_Cgelsy = 0;
    INTEGER lrwork_Cgelsy = 0;
    INTEGER lwork_Cgelss = 0;
    INTEGER lrwork_Cgelss = 0;
    INTEGER lwork_Cgelsd = 0;
    INTEGER lrwork_Cgelsd = 0;
    INTEGER lwlsy = 0;
    INTEGER mb = 0;
    REAL norma = 0.0;
    INTEGER inb = 0;
    INTEGER nb = 0;
    INTEGER nrows = 0;
    INTEGER ncols = 0;
    INTEGER ldwork = 0;
    COMPLEX work[lwork];
    const REAL one = 1.0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    REAL rwork[lrwork];
    const INTEGER ntests = 16;
    REAL result[ntests];
    INTEGER k = 0;
    INTEGER imb = 0;
    REAL work2[2 * lwork];
    INTEGER rank = 0;
    REAL normb = 0.0;
    INTEGER j = 0;
    INTEGER iwork[liwork];
    //
    //     Initialize constants and the random number seed.
    //
    path[0] = 'C';
    path[1] = 'L';
    path[2] = 'S';
    nrun = 0;
    nfail = 0;
    nerrs = 0;
    for (i = 1; i <= 4; i = i + 1) {
        iseed[i - 1] = iseedy[i - 1];
    }
    eps = Rlamch("Epsilon");
    //
    //     Threshold for rank estimation
    //
    rcond = sqrt(eps) - (sqrt(eps) - eps) / 2.0;
    //
    //     Test the error exits
    //
    xlaenv(9, smlsiz);
    if (tsterr) {
        Cerrls(path, nout);
    }
    //
    //     Print the header if NM = 0 or NN = 0 and THRESH = 0.
    //
    if ((nm == 0 || nn == 0) && thresh == zero) {
        Alahd(nout, path);
    }
    infot = 0;
    //
    //     Compute maximal workspace needed for all routines
    //
    nmax = 0;
    mmax = 0;
    nsmax = 0;
    for (i = 1; i <= nm; i = i + 1) {
        if (mval[i - 1] > mmax) {
            mmax = mval[i - 1];
        }
    }
    for (i = 1; i <= nn; i = i + 1) {
        if (nval[i - 1] > nmax) {
            nmax = nval[i - 1];
        }
    }
    for (i = 1; i <= nns; i = i + 1) {
        if (nsval[i - 1] > nsmax) {
            nsmax = nsval[i - 1];
        }
    }
    m = mmax;
    n = nmax;
    nrhs = nsmax;
    mnmin = max(min(m, n), (INTEGER)1);
    //
    //     Compute workspace needed for routines
    //     Cqrt14, Cqrt17 (two side cases), Cqrt15 and Cqrt12
    //
    lwork = max({(INTEGER)1, (m + n) * nrhs, (n + nrhs) * (m + 2), (m + nrhs) * (n + 2), max({m + mnmin, nrhs * mnmin, 2 * n + m}), max({m * n + 4 * mnmin + max(m, n), m * n + 2 * mnmin + 4 * n})});
    lrwork = 1;
    liwork = 1;
    //
    //     Iterate through all test cases and compute necessary workspace
    //     sizes for ?GELS, ?GETSLS, ?GELSY, ?GELSS and ?GELSD routines.
    //
    for (im = 1; im <= nm; im = im + 1) {
        m = mval[im - 1];
        lda = max((INTEGER)1, m);
        for (in = 1; in <= nn; in = in + 1) {
            n = nval[in - 1];
            mnmin = max(min(m, n), (INTEGER)1);
            ldb = max({(INTEGER)1, m, n});
            for (ins = 1; ins <= nns; ins = ins + 1) {
                nrhs = nsval[ins - 1];
                for (irank = 1; irank <= 2; irank = irank + 1) {
                    for (iscale = 1; iscale <= 3; iscale = iscale + 1) {
                        itype = (irank - 1) * 3 + iscale;
                        if (dotype[itype - 1]) {
                            if (irank == 1) {
                                for (itran = 1; itran <= 2; itran = itran + 1) {
                                    if (itran == 1) {
                                        trans[0] = 'N';
                                    } else {
                                        trans[0] = 'C';
                                    }
                                    //
                                    //                             Compute workspace needed for Cgels
                                    Cgels(trans, m, n, nrhs, a, lda, b, ldb, wq, -1, info);
                                    lwork_Cgels = castINTEGER(wq[1 - 1].real());
                                    //                             Compute workspace needed for Cgetsls
                                    Cgetsls(trans, m, n, nrhs, a, lda, b, ldb, wq, -1, info);
                                    lwork_Cgetsls = castINTEGER(wq[1 - 1].real());
                                }
                            }
                            //                       Compute workspace needed for Cgelsy
                            Cgelsy(m, n, nrhs, a, lda, b, ldb, iwq, rcond, crank, wq, -1, rwq, info);
                            lwork_Cgelsy = castINTEGER(wq[1 - 1].real());
                            lrwork_Cgelsy = 2 * n;
                            //                       Compute workspace needed for Cgelss
                            Cgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, crank, wq, -1, rwq, info);
                            lwork_Cgelss = castINTEGER(wq[1 - 1].real());
                            lrwork_Cgelss = 5 * mnmin;
                            //                       Compute workspace needed for Cgelsd
                            Cgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, crank, wq, -1, rwq, iwq, info);
                            lwork_Cgelsd = castINTEGER(wq[1 - 1].real());
                            lrwork_Cgelsd = castINTEGER(rwq[1 - 1]);
                            //                       Compute LIWORK workspace needed for Cgelsy and Cgelsd
                            liwork = max({liwork, n, iwq[1 - 1]});
                            //                       Compute LRWORK workspace needed for Cgelsy, Cgelss and Cgelsd
                            lrwork = max({lrwork, lrwork_Cgelsy, lrwork_Cgelss, lrwork_Cgelsd});
                            //                       Compute LWORK workspace needed for all functions
                            lwork = max({lwork, lwork_Cgels, lwork_Cgetsls, lwork_Cgelsy, lwork_Cgelss, lwork_Cgelsd});
                        }
                    }
                }
            }
        }
    }
    //
    lwlsy = lwork;
    //
    for (im = 1; im <= nm; im = im + 1) {
        m = mval[im - 1];
        lda = max((INTEGER)1, m);
        //
        for (in = 1; in <= nn; in = in + 1) {
            n = nval[in - 1];
            mnmin = max(min(m, n), (INTEGER)1);
            ldb = max({(INTEGER)1, m, n});
            mb = (mnmin + 1);
            //
            for (ins = 1; ins <= nns; ins = ins + 1) {
                nrhs = nsval[ins - 1];
                //
                for (irank = 1; irank <= 2; irank = irank + 1) {
                    for (iscale = 1; iscale <= 3; iscale = iscale + 1) {
                        itype = (irank - 1) * 3 + iscale;
                        if (!dotype[itype - 1]) {
                            goto statement_100;
                        }
                        //
                        if (irank == 1) {
                            //
                            //                       Test Cgels
                            //
                            //                       Generate a matrix of scaling type ISCALE
                            //
                            Cqrt13(iscale, m, n, copya, lda, norma, iseed);
                            for (inb = 1; inb <= nnb; inb = inb + 1) {
                                nb = nbval[inb - 1];
                                //
                                for (itran = 1; itran <= 2; itran = itran + 1) {
                                    if (itran == 1) {
                                        trans[0] = 'N';
                                        nrows = m;
                                        ncols = n;
                                    } else {
                                        trans[0] = 'C';
                                        nrows = n;
                                        ncols = m;
                                    }
                                    ldwork = max((INTEGER)1, ncols);
                                    //
                                    //                             Set up a consistent rhs
                                    //
                                    if (ncols > 0) {
                                        Clarnv(2, iseed, ncols * nrhs, work);
                                        CRscal(ncols * nrhs, one / castREAL(ncols), work, 1);
                                    }
                                    Cgemm(trans, "No transpose", nrows, nrhs, ncols, cone, copya, lda, work, ldwork, czero, b, ldb);
                                    Clacpy("Full", nrows, nrhs, b, ldb, copyb, ldb);
                                    //
                                    //                             Solve LS or overdetermined system
                                    //
                                    if (m > 0 && n > 0) {
                                        Clacpy("Full", m, n, copya, lda, a, lda);
                                        Clacpy("Full", nrows, nrhs, copyb, ldb, b, ldb);
                                    }
                                    Cgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
                                    //
                                    if (info != 0) {
                                        Alaerh(path, "Cgels ", info, 0, trans, m, n, nrhs, -1, nb, itype, nfail, nerrs, nout);
                                    }
                                    //
                                    //                             Check correctness of results
                                    //
                                    ldwork = max((INTEGER)1, nrows);
                                    if (nrows > 0 && nrhs > 0) {
                                        Clacpy("Full", nrows, nrhs, copyb, ldb, c, ldb);
                                    }
                                    Cqrt16(trans, m, n, nrhs, copya, lda, b, ldb, c, ldb, rwork, result[1 - 1]);
                                    //
                                    if ((itran == 1 && m >= n) || (itran == 2 && m < n)) {
                                        //
                                        //                                Solving LS system
                                        //
                                        result[2 - 1] = Cqrt17(trans, 1, m, n, nrhs, copya, lda, b, ldb, copyb, ldb, c, work, lwork);
                                    } else {
                                        //
                                        //                                Solving overdetermined system
                                        //
                                        result[2 - 1] = Cqrt14(trans, m, n, nrhs, copya, lda, b, ldb, work, lwork);
                                    }
                                    //
                                    //                             Print information about the tests that
                                    //                             did not pass the threshold.
                                    //
                                    for (k = 1; k <= 2; k = k + 1) {
                                        if (result[k - 1] >= thresh) {
                                            if (nfail == 0 && nerrs == 0) {
                                                Alahd(nout, path);
                                            }
                                            sprintnum_short(buf, result[k - 1]);
                                            write(nout, "(' TRANS=''',a1,''', M=',i5,', N=',i5,', NRHS=',i4,"
                                                        "', NB=',i4,', type',i2,', test(',i2,')=',a)"),
                                                trans, m, n, nrhs, nb, itype, k, buf;
                                            nfail++;
                                        }
                                    }
                                    nrun += 2;
                                }
                            }
                            //
                            //                       Test Cgetsls
                            //
                            //                       Generate a matrix of scaling type ISCALE
                            //
                            Cqrt13(iscale, m, n, copya, lda, norma, iseed);
                            for (inb = 1; inb <= nnb; inb = inb + 1) {
                                mb = nbval[inb - 1];
                                for (imb = 1; imb <= nnb; imb = imb + 1) {
                                    nb = nbval[imb - 1];
                                    //
                                    for (itran = 1; itran <= 2; itran = itran + 1) {
                                        if (itran == 1) {
                                            trans[0] = 'N';
                                            nrows = m;
                                            ncols = n;
                                        } else {
                                            trans[0] = 'C';
                                            nrows = n;
                                            ncols = m;
                                        }
                                        ldwork = max((INTEGER)1, ncols);
                                        //
                                        //                             Set up a consistent rhs
                                        //
                                        if (ncols > 0) {
                                            Clarnv(2, iseed, ncols * nrhs, work);
                                            Cscal(ncols * nrhs, cone / castREAL(ncols), work, 1);
                                        }
                                        Cgemm(trans, "No transpose", nrows, nrhs, ncols, cone, copya, lda, work, ldwork, czero, b, ldb);
                                        Clacpy("Full", nrows, nrhs, b, ldb, copyb, ldb);
                                        //
                                        //                             Solve LS or overdetermined system
                                        //
                                        if (m > 0 && n > 0) {
                                            Clacpy("Full", m, n, copya, lda, a, lda);
                                            Clacpy("Full", nrows, nrhs, copyb, ldb, b, ldb);
                                        }
                                        Cgetsls(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
                                        if (info != 0) {
                                            Alaerh(path, "Cgetsls ", info, 0, trans, m, n, nrhs, -1, nb, itype, nfail, nerrs, nout);
                                        }
                                        //
                                        //                             Check correctness of results
                                        //
                                        ldwork = max((INTEGER)1, nrows);
                                        if (nrows > 0 && nrhs > 0) {
                                            Clacpy("Full", nrows, nrhs, copyb, ldb, c, ldb);
                                        }
                                        Cqrt16(trans, m, n, nrhs, copya, lda, b, ldb, c, ldb, work2, result[15 - 1]);
                                        //
                                        if ((itran == 1 && m >= n) || (itran == 2 && m < n)) {
                                            //
                                            //                                Solving LS system
                                            //
                                            result[16 - 1] = Cqrt17(trans, 1, m, n, nrhs, copya, lda, b, ldb, copyb, ldb, c, work, lwork);
                                        } else {
                                            //
                                            //                                Solving overdetermined system
                                            //
                                            result[16 - 1] = Cqrt14(trans, m, n, nrhs, copya, lda, b, ldb, work, lwork);
                                        }
                                        //
                                        //                             Print information about the tests that
                                        //                             did not pass the threshold.
                                        //
                                        for (k = 15; k <= 16; k = k + 1) {
                                            if (result[k - 1] >= thresh) {
                                                if (nfail == 0 && nerrs == 0) {
                                                    Alahd(nout, path);
                                                }
                                                sprintnum_short(buf, result[k - 1]);
                                                write(nout, "(' TRANS=''',a1,' M=',i5,', N=',i5,', NRHS=',i4,"
                                                            "', MB=',i4,', NB=',i4,', type',i2,', test(',i2,"
                                                            "')=',a)"),
                                                    trans, m, n, nrhs, mb, nb, itype, k, buf;
                                                nfail++;
                                            }
                                        }
                                        nrun += 2;
                                    }
                                }
                            }
                        }
                        //
                        //                    Generate a matrix of scaling type ISCALE and rank
                        //                    type IRANK.
                        //
                        Cqrt15(iscale, irank, m, n, nrhs, copya, lda, copyb, ldb, copys, rank, norma, normb, iseed, work, lwork);
                        //
                        //                    workspace used: MAX(M+MIN(M,N),NRHS*MIN(M,N),2*N+M)
                        //
                        ldwork = max((INTEGER)1, m);
                        //
                        //                    Loop for testing different block sizes.
                        //
                        for (inb = 1; inb <= nnb; inb = inb + 1) {
                            nb = nbval[inb - 1];
                            //
                            //                       Test Cgelsy
                            //
                            //                       Cgelsy:  Compute the minimum-norm solution
                            //                       X to min( norm( A * X - B ) )
                            //                       using the rank-revealing orthogonal
                            //                       factorization.
                            //
                            Clacpy("Full", m, n, copya, lda, a, lda);
                            Clacpy("Full", m, nrhs, copyb, ldb, b, ldb);
                            //
                            //                       Initialize vector IWORK.
                            //
                            for (j = 1; j <= n; j = j + 1) {
                                iwork[j - 1] = 0;
                            }
                            //
                            Cgelsy(m, n, nrhs, a, lda, b, ldb, iwork, rcond, crank, work, lwlsy, rwork, info);
                            if (info != 0) {
                                Alaerh(path, "Cgelsy", info, 0, " ", m, n, nrhs, -1, nb, itype, nfail, nerrs, nout);
                            }
                            //
                            //                       workspace used: 2*MNMIN+NB*NB+NB*MAX(N,NRHS)
                            //
                            //                       Test 3:  Compute relative error in svd
                            //                                workspace: M*N + 4*MIN(M,N) + MAX(M,N)
                            //
                            result[3 - 1] = Cqrt12(crank, crank, a, lda, copys, work, lwork, rwork);
                            //
                            //                       Test 4:  Compute error in solution
                            //                                workspace:  M*NRHS + M
                            //
                            Clacpy("Full", m, nrhs, copyb, ldb, work, ldwork);
                            Cqrt16("No transpose", m, n, nrhs, copya, lda, b, ldb, work, ldwork, rwork, result[4 - 1]);
                            //
                            //                       Test 5:  Check norm of r'*A
                            //                                workspace: NRHS*(M+N)
                            //
                            result[5 - 1] = zero;
                            if (m > crank) {
                                result[5 - 1] = Cqrt17("No transpose", 1, m, n, nrhs, copya, lda, b, ldb, copyb, ldb, c, work, lwork);
                            }
                            //
                            //                       Test 6:  Check if x is in the rowspace of A
                            //                                workspace: (M+NRHS)*(N+2)
                            //
                            result[6 - 1] = zero;
                            //
                            if (n > crank) {
                                result[6 - 1] = Cqrt14("No transpose", m, n, nrhs, copya, lda, b, ldb, work, lwork);
                            }
                            //
                            //                       Test Cgelss
                            //
                            //                       Cgelss:  Compute the minimum-norm solution
                            //                       X to min( norm( A * X - B ) )
                            //                       using the SVD.
                            //
                            Clacpy("Full", m, n, copya, lda, a, lda);
                            Clacpy("Full", m, nrhs, copyb, ldb, b, ldb);
                            Cgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, crank, work, lwork, rwork, info);
                            //
                            if (info != 0) {
                                Alaerh(path, "Cgelss", info, 0, " ", m, n, nrhs, -1, nb, itype, nfail, nerrs, nout);
                            }
                            //
                            //                       workspace used: 3*min(m,n) +
                            //                                       max((INTEGER)2*min(m,n),nrhs,max(m,n))
                            //
                            //                       Test 7:  Compute relative error in svd
                            //
                            if (rank > 0) {
                                Raxpy(mnmin, -one, copys, 1, s, 1);
                                result[7 - 1] = Rasum(mnmin, s, 1) / Rasum(mnmin, copys, 1) / (eps * castREAL(mnmin));
                            } else {
                                result[7 - 1] = zero;
                            }
                            //
                            //                       Test 8:  Compute error in solution
                            //
                            Clacpy("Full", m, nrhs, copyb, ldb, work, ldwork);
                            Cqrt16("No transpose", m, n, nrhs, copya, lda, b, ldb, work, ldwork, rwork, result[8 - 1]);
                            //
                            //                       Test 9:  Check norm of r'*A
                            //
                            result[9 - 1] = zero;
                            if (m > crank) {
                                result[9 - 1] = Cqrt17("No transpose", 1, m, n, nrhs, copya, lda, b, ldb, copyb, ldb, c, work, lwork);
                            }
                            //
                            //                       Test 10:  Check if x is in the rowspace of A
                            //
                            result[10 - 1] = zero;
                            if (n > crank) {
                                result[10 - 1] = Cqrt14("No transpose", m, n, nrhs, copya, lda, b, ldb, work, lwork);
                            }
                            //
                            //                       Test Cgelsd
                            //
                            //                       Cgelsd:  Compute the minimum-norm solution X
                            //                       to min( norm( A * X - B ) ) using a
                            //                       divide and conquer SVD.
                            //
                            //
                            Clacpy("Full", m, n, copya, lda, a, lda);
                            Clacpy("Full", m, nrhs, copyb, ldb, b, ldb);
                            //
                            Cgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, crank, work, lwork, rwork, iwork, info);
                            if (info != 0) {
                                Alaerh(path, "Cgelsd", info, 0, " ", m, n, nrhs, -1, nb, itype, nfail, nerrs, nout);
                            }
                            //
                            //                       Test 11:  Compute relative error in svd
                            //
                            if (rank > 0) {
                                Raxpy(mnmin, -one, copys, 1, s, 1);
                                result[11 - 1] = Rasum(mnmin, s, 1) / Rasum(mnmin, copys, 1) / (eps * castREAL(mnmin));
                            } else {
                                result[11 - 1] = zero;
                            }
                            //
                            //                       Test 12:  Compute error in solution
                            //
                            Clacpy("Full", m, nrhs, copyb, ldb, work, ldwork);
                            Cqrt16("No transpose", m, n, nrhs, copya, lda, b, ldb, work, ldwork, rwork, result[12 - 1]);
                            //
                            //                       Test 13:  Check norm of r'*A
                            //
                            result[13 - 1] = zero;
                            if (m > crank) {
                                result[13 - 1] = Cqrt17("No transpose", 1, m, n, nrhs, copya, lda, b, ldb, copyb, ldb, c, work, lwork);
                            }
                            //
                            //                       Test 14:  Check if x is in the rowspace of A
                            //
                            result[14 - 1] = zero;
                            if (n > crank) {
                                result[14 - 1] = Cqrt14("No transpose", m, n, nrhs, copya, lda, b, ldb, work, lwork);
                            }
                            //
                            //                       Print information about the tests that did not
                            //                       pass the threshold.
                            //
                            for (k = 3; k <= 14; k = k + 1) {
                                if (result[k - 1] >= thresh) {
                                    if (nfail == 0 && nerrs == 0) {
                                        Alahd(nout, path);
                                    }
                                    sprintnum_short(buf, result[k - 1]);
                                    write(nout, "(' M=',i5,', N=',i5,', NRHS=',i4,', NB=',i4,', type',i2,"
                                                "', test(',i2,')=',a)"),
                                        m, n, nrhs, nb, itype, k, buf;
                                    nfail++;
                                }
                            }
                            nrun += 12;
                            //
                        }
                    statement_100:;
                    }
                }
            }
        }
    }
    //
    //     Print a summary of the results.
    //
    Alasvm(path, nout, nfail, nrun, nerrs);
    //
    FEM_THROW_UNHANDLED("executable deallocate: deallocate(work)");
    FEM_THROW_UNHANDLED("executable deallocate: deallocate(iwork)");
    FEM_THROW_UNHANDLED("executable deallocate: deallocate(rwork)");
    //
    //     End of Cdrvls
    //
}
