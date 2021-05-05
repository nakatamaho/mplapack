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
#include <mplapack_eig.h>

#include <mplapack_debug.h>

void Cdrvev(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *dotype, INTEGER *iseed, REAL const thresh, INTEGER const nounit, COMPLEX *a, INTEGER const lda, COMPLEX *h, COMPLEX *w, COMPLEX *w1, COMPLEX *vl, INTEGER const ldvl, COMPLEX *vr, INTEGER const ldvr, COMPLEX *lre, INTEGER const ldlre, REAL *result, COMPLEX *work, INTEGER const nwork, REAL *rwork, INTEGER *iwork, INTEGER &info) {
    FEM_CMN_SVE(Cdrvev);
    iseed([4]);
    a([lda * star]);
    h([lda * star]);
    vl([ldvl * star]);
    vr([ldvr * star]);
    lre([ldlre * star]);
    result([7]);
    common_write write(cmn);
    const INTEGER maxtyp = 21;
    INTEGER *kconds(sve.kconds, [maxtyp]);
    INTEGER *kmagn(sve.kmagn, [maxtyp]);
    INTEGER *kmode(sve.kmode, [maxtyp]);
    INTEGER *ktype(sve.ktype, [maxtyp]);
    if (is_called_first_time) {
        data((values, 1, 2, 3, 5 * datum(4), 4 * datum(6), 6 * datum(6), 3 * datum(9))), ktype;
        {
            data_values data;
            data.values, 3 * datum(1), 1, 1, 1, 2, 3, 4 * datum(1), 1;
            data.values, 1, 1, 1, 2, 3, 1, 2, 3;
            data, kmagn;
        }
        {
            data_values data;
            data.values, 3 * datum(0), 4, 3, 1, 4, 4, 4, 3;
            data.values, 1, 5, 4, 3, 1, 5, 5, 5;
            data.values, 4, 3, 1;
            data, kmode;
        }
        data((values, 3 * datum(0), 5 * datum(0), 4 * datum(1), 6 * datum(2), 3 * datum(0))), kconds;
    }
    char path[3];
    INTEGER ntestt = 0;
    INTEGER ntestf = 0;
    bool badnn = false;
    INTEGER nmax = 0;
    INTEGER j = 0;
    const REAL zero = 0.0;
    REAL unfl = 0.0;
    const REAL one = 1.0;
    REAL ovfl = 0.0;
    REAL ulp = 0.0;
    REAL ulpinv = 0.0;
    REAL rtulp = 0.0;
    REAL rtulpi = 0.0;
    INTEGER nerrs = 0;
    INTEGER jsize = 0;
    INTEGER n = 0;
    INTEGER mtypes = 0;
    INTEGER jtype = 0;
    INTEGER ioldsd[4];
    INTEGER itype = 0;
    INTEGER imode = 0;
    REAL anorm = 0.0;
    const COMPLEX czero = COMPLEX(0.0, 0.0);
    INTEGER iinfo = 0;
    REAL cond = 0.0;
    INTEGER jcol = 0;
    const COMPLEX cone = COMPLEX(1.0, 0.0);
    REAL conds = 0.0;
    INTEGER idumma[1];
    INTEGER iwk = 0;
    INTEGER nnwork = 0;
    REAL res[2];
    REAL tnrm = 0.0;
    REAL vmx = 0.0;
    REAL vrmx = 0.0;
    INTEGER jj = 0;
    REAL vtst = 0.0;
    const REAL two = 2.0e+0;
    COMPLEX dum[1];
    INTEGER ntest = 0;
    INTEGER nfail = 0;
    static const char *format_9993 = "(' Cdrvev: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,', JTYPE=',i6,"
                                     "', ISEED=(',3(i5,','),i5,')')";
    //
    //  -- LAPACK test routine --
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
    //     .. Local Arrays ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. External Subroutines ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..
    //
    path[(1 - 1)] = "Zomplex precision";
    path[(2 - 1) + (3 - 1) * ldpath] = "EV";
    //
    //     Check for errors
    //
    ntestt = 0;
    ntestf = 0;
    info = 0;
    //
    //     Important constants
    //
    badnn = false;
    nmax = 0;
    for (j = 1; j <= nsizes; j = j + 1) {
        nmax = max(nmax, nn[j - 1]);
        if (nn[j - 1] < 0) {
            badnn = true;
        }
    }
    //
    //     Check for errors
    //
    if (nsizes < 0) {
        info = -1;
    } else if (badnn) {
        info = -2;
    } else if (ntypes < 0) {
        info = -3;
    } else if (thresh < zero) {
        info = -6;
    } else if (nounit <= 0) {
        info = -7;
    } else if (lda < 1 || lda < nmax) {
        info = -9;
    } else if (ldvl < 1 || ldvl < nmax) {
        info = -14;
    } else if (ldvr < 1 || ldvr < nmax) {
        info = -16;
    } else if (ldlre < 1 || ldlre < nmax) {
        info = -28;
    } else if (5 * nmax + 2 * pow2(nmax) > nwork) {
        info = -21;
    }
    //
    if (info != 0) {
        Mxerbla("Cdrvev", -info);
        return;
    }
    //
    //     Quick return if nothing to do
    //
    if (nsizes == 0 || ntypes == 0) {
        return;
    }
    //
    //     More Important constants
    //
    unfl = Rlamch("Safe minimum");
    ovfl = one / unfl;
    Rlabad(unfl, ovfl);
    ulp = Rlamch("Precision");
    ulpinv = one / ulp;
    rtulp = sqrt(ulp);
    rtulpi = one / rtulp;
    //
    //     Loop over sizes, types
    //
    nerrs = 0;
    //
    for (jsize = 1; jsize <= nsizes; jsize = jsize + 1) {
        n = nn[jsize - 1];
        if (nsizes != 1) {
            mtypes = min(maxtyp, ntypes);
        } else {
            mtypes = min(maxtyp + 1, ntypes);
        }
        //
        for (jtype = 1; jtype <= mtypes; jtype = jtype + 1) {
            if (!dotype[jtype - 1]) {
                goto statement_260;
            }
            //
            //           Save ISEED in case of an error.
            //
            for (j = 1; j <= 4; j = j + 1) {
                ioldsd[j - 1] = iseed[j - 1];
            }
            //
            //           Compute "A"
            //
            //           Control parameters:
            //
            //           KMAGN  KCONDS  KMODE        KTYPE
            //       =1  O(1)   1       clustered 1  zero
            //       =2  large  large   clustered 2  identity
            //       =3  small          exponential  Jordan
            //       =4                 arithmetic   diagonal, (w/ eigenvalues)
            //       =5                 random log   symmetric, w/ eigenvalues
            //       =6                 random       general, w/ eigenvalues
            //       =7                              random diagonal
            //       =8                              random symmetric
            //       =9                              random general
            //       =10                             random triangular
            //
            if (mtypes > maxtyp) {
                goto statement_90;
            }
            //
            itype = ktype[jtype - 1];
            imode = kmode[jtype - 1];
            //
            //           Compute norm
            //
            switch (kmagn[jtype - 1]) {
            case 1:
                goto statement_30;
            case 2:
                goto statement_40;
            case 3:
                goto statement_50;
            default:
                break;
            }
        //
        statement_30:
            anorm = one;
            goto statement_60;
        //
        statement_40:
            anorm = ovfl * ulp;
            goto statement_60;
        //
        statement_50:
            anorm = unfl * ulpinv;
            goto statement_60;
        //
        statement_60:
            //
            Claset("Full", lda, n, czero, czero, a, lda);
            iinfo = 0;
            cond = ulpinv;
            //
            //           Special Matrices -- Identity & Jordan block
            //
            //              Zero
            //
            if (itype == 1) {
                iinfo = 0;
                //
            } else if (itype == 2) {
                //
                //              Identity
                //
                for (jcol = 1; jcol <= n; jcol = jcol + 1) {
                    a[(jcol - 1) + (jcol - 1) * lda] = COMPLEX(anorm);
                }
                //
            } else if (itype == 3) {
                //
                //              Jordan Block
                //
                for (jcol = 1; jcol <= n; jcol = jcol + 1) {
                    a[(jcol - 1) + (jcol - 1) * lda] = COMPLEX(anorm);
                    if (jcol > 1) {
                        a[(jcol - 1) + ((jcol - 1) - 1) * lda] = cone;
                    }
                }
                //
            } else if (itype == 4) {
                //
                //              Diagonal Matrix, [Eigen]values Specified
                //
                zlatms(n, n, "S", iseed, "H", rwork, imode, cond, anorm, 0, 0, "N", a, lda, &work[(n + 1) - 1], iinfo);
                //
            } else if (itype == 5) {
                //
                //              Hermitian, eigenvalues specified
                //
                zlatms(n, n, "S", iseed, "H", rwork, imode, cond, anorm, n, n, "N", a, lda, &work[(n + 1) - 1], iinfo);
                //
            } else if (itype == 6) {
                //
                //              General, eigenvalues specified
                //
                if (kconds[jtype - 1] == 1) {
                    conds = one;
                } else if (kconds[jtype - 1] == 2) {
                    conds = rtulpi;
                } else {
                    conds = zero;
                }
                //
                zlatme(n, "D", iseed, work, imode, cond, cone, "T", "T", "T", rwork, 4, conds, n, n, anorm, a, lda, &work[(2 * n + 1) - 1], iinfo);
                //
            } else if (itype == 7) {
                //
                //              Diagonal, random eigenvalues
                //
                zlatmr(n, n, "D", iseed, "N", work, 6, one, cone, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, 0, 0, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 8) {
                //
                //              Symmetric, random eigenvalues
                //
                zlatmr(n, n, "D", iseed, "H", work, 6, one, cone, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, n, n, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 9) {
                //
                //              General, random eigenvalues
                //
                zlatmr(n, n, "D", iseed, "N", work, 6, one, cone, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, n, n, zero, anorm, "NO", a, lda, iwork, iinfo);
                if (n >= 4) {
                    Claset("Full", 2, n, czero, czero, a, lda);
                    Claset("Full", n - 3, 1, czero, czero, &a[(3 - 1)], lda);
                    Claset("Full", n - 3, 2, czero, czero, &a[(3 - 1) + ((n - 1) - 1) * lda], lda);
                    Claset("Full", 1, n, czero, czero, &a[(n - 1)], lda);
                }
                //
            } else if (itype == 10) {
                //
                //              Triangular, random eigenvalues
                //
                zlatmr(n, n, "D", iseed, "N", work, 6, one, cone, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, n, 0, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else {
                //
                iinfo = 1;
            }
            //
            if (iinfo != 0) {
                write(nounit, format_9993), "Generator", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                return;
            }
        //
        statement_90:
            //
            //           Test for minimal and generous workspace
            //
            for (iwk = 1; iwk <= 2; iwk = iwk + 1) {
                if (iwk == 1) {
                    nnwork = 2 * n;
                } else {
                    nnwork = 5 * n + 2 * pow2(n);
                }
                nnwork = max(nnwork, 1);
                //
                //              Initialize RESULT
                //
                for (j = 1; j <= 7; j = j + 1) {
                    result[j - 1] = -one;
                }
                //
                //              Compute eigenvalues and eigenvectors, and test them
                //
                Clacpy("F", n, n, a, lda, h, lda);
                Cgeev("V", "V", n, h, lda, w, vl, ldvl, vr, ldvr, work, nnwork, rwork, iinfo);
                if (iinfo != 0) {
                    result[1 - 1] = ulpinv;
                    write(nounit, format_9993), "Cgeev1", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    goto statement_220;
                }
                //
                //              Do Test (1)
                //
                Cget22("N", "N", "N", n, a, lda, vr, ldvr, w, work, rwork, res);
                result[1 - 1] = res[1 - 1];
                //
                //              Do Test (2)
                //
                Cget22("C", "N", "C", n, a, lda, vl, ldvl, w, work, rwork, res);
                result[2 - 1] = res[1 - 1];
                //
                //              Do Test (3)
                //
                for (j = 1; j <= n; j = j + 1) {
                    tnrm = RCnrm2(n, &vr[(j - 1) * ldvr], 1);
                    result[3 - 1] = max({result[3 - 1], min(ulpinv, abs(tnrm - one) / ulp)});
                    vmx = zero;
                    vrmx = zero;
                    for (jj = 1; jj <= n; jj = jj + 1) {
                        vtst = abs(vr[(jj - 1) + (j - 1) * ldvr]);
                        if (vtst > vmx) {
                            vmx = vtst;
                        }
                        if (vr[(jj - 1) + (j - 1) * ldvr].imag() == zero && abs(vr[(jj - 1) + (j - 1) * ldvr].real()) > vrmx) {
                            vrmx = abs(vr[(jj - 1) + (j - 1) * ldvr].real());
                        }
                    }
                    if (vrmx / vmx < one - two * ulp) {
                        result[3 - 1] = ulpinv;
                    }
                }
                //
                //              Do Test (4)
                //
                for (j = 1; j <= n; j = j + 1) {
                    tnrm = RCnrm2(n, &vl[(j - 1) * ldvl], 1);
                    result[4 - 1] = max({result[4 - 1], min(ulpinv, abs(tnrm - one) / ulp)});
                    vmx = zero;
                    vrmx = zero;
                    for (jj = 1; jj <= n; jj = jj + 1) {
                        vtst = abs(vl[(jj - 1) + (j - 1) * ldvl]);
                        if (vtst > vmx) {
                            vmx = vtst;
                        }
                        if (vl[(jj - 1) + (j - 1) * ldvl].imag() == zero && abs(vl[(jj - 1) + (j - 1) * ldvl].real()) > vrmx) {
                            vrmx = abs(vl[(jj - 1) + (j - 1) * ldvl].real());
                        }
                    }
                    if (vrmx / vmx < one - two * ulp) {
                        result[4 - 1] = ulpinv;
                    }
                }
                //
                //              Compute eigenvalues only, and test them
                //
                Clacpy("F", n, n, a, lda, h, lda);
                Cgeev("N", "N", n, h, lda, w1, dum, 1, dum, 1, work, nnwork, rwork, iinfo);
                if (iinfo != 0) {
                    result[1 - 1] = ulpinv;
                    write(nounit, format_9993), "Cgeev2", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    goto statement_220;
                }
                //
                //              Do Test (5)
                //
                for (j = 1; j <= n; j = j + 1) {
                    if (w[j - 1] != w1[j - 1]) {
                        result[5 - 1] = ulpinv;
                    }
                }
                //
                //              Compute eigenvalues and right eigenvectors, and test them
                //
                Clacpy("F", n, n, a, lda, h, lda);
                Cgeev("N", "V", n, h, lda, w1, dum, 1, lre, ldlre, work, nnwork, rwork, iinfo);
                if (iinfo != 0) {
                    result[1 - 1] = ulpinv;
                    write(nounit, format_9993), "Cgeev3", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    goto statement_220;
                }
                //
                //              Do Test (5) again
                //
                for (j = 1; j <= n; j = j + 1) {
                    if (w[j - 1] != w1[j - 1]) {
                        result[5 - 1] = ulpinv;
                    }
                }
                //
                //              Do Test (6)
                //
                for (j = 1; j <= n; j = j + 1) {
                    for (jj = 1; jj <= n; jj = jj + 1) {
                        if (vr[(j - 1) + (jj - 1) * ldvr] != lre[(j - 1) + (jj - 1) * ldlre]) {
                            result[6 - 1] = ulpinv;
                        }
                    }
                }
                //
                //              Compute eigenvalues and left eigenvectors, and test them
                //
                Clacpy("F", n, n, a, lda, h, lda);
                Cgeev("V", "N", n, h, lda, w1, lre, ldlre, dum, 1, work, nnwork, rwork, iinfo);
                if (iinfo != 0) {
                    result[1 - 1] = ulpinv;
                    write(nounit, format_9993), "Cgeev4", iinfo, n, jtype, ioldsd;
                    info = abs(iinfo);
                    goto statement_220;
                }
                //
                //              Do Test (5) again
                //
                for (j = 1; j <= n; j = j + 1) {
                    if (w[j - 1] != w1[j - 1]) {
                        result[5 - 1] = ulpinv;
                    }
                }
                //
                //              Do Test (7)
                //
                for (j = 1; j <= n; j = j + 1) {
                    for (jj = 1; jj <= n; jj = jj + 1) {
                        if (vl[(j - 1) + (jj - 1) * ldvl] != lre[(j - 1) + (jj - 1) * ldlre]) {
                            result[7 - 1] = ulpinv;
                        }
                    }
                }
            //
            //              End of Loop -- Check for RESULT(j) > THRESH
            //
            statement_220:
                //
                ntest = 0;
                nfail = 0;
                for (j = 1; j <= 7; j = j + 1) {
                    if (result[j - 1] >= zero) {
                        ntest++;
                    }
                    if (result[j - 1] >= thresh) {
                        nfail++;
                    }
                }
                //
                if (nfail > 0) {
                    ntestf++;
                }
                if (ntestf == 1) {
                    write(nounit, "(/,1x,a3,' -- Complex Eigenvalue-Eigenvector ',"
                                  "'Decomposition Driver',/,"
                                  "' Matrix types (see Cdrvev for details): ')"),
                        path;
                    write(nounit, "(/,' Special Matrices:',/,'  1=Zero matrix.             ',"
                                  "'           ','  5=Diagonal: geometr. spaced entries.',/,"
                                  "'  2=Identity matrix.                    ','  6=Diagona',"
                                  "'l: clustered entries.',/,'  3=Transposed Jordan block.  ',"
                                  "'          ','  7=Diagonal: large, evenly spaced.',/,'  ',"
                                  "'4=Diagonal: evenly spaced entries.    ','  8=Diagonal: s',"
                                  "'mall, evenly spaced.')");
                    write(nounit, "(' Dense, Non-Symmetric Matrices:',/,'  9=Well-cond., ev',"
                                  "'enly spaced eigenvals.',' 14=Ill-cond., geomet. spaced e',"
                                  "'igenals.',/,' 10=Well-cond., geom. spaced eigenvals. ',"
                                  "' 15=Ill-conditioned, clustered e.vals.',/,' 11=Well-cond',"
                                  "'itioned, clustered e.vals. ',' 16=Ill-cond., random comp',"
                                  "'lex ',a6,/,' 12=Well-cond., random complex ',a6,'   ',"
                                  "' 17=Ill-cond., large rand. complx ',a4,/,' 13=Ill-condi',"
                                  "'tioned, evenly spaced.     ',' 18=Ill-cond., small rand.',"
                                  "' complx ',a4)");
                    write(nounit, "(' 19=Matrix with random O(1) entries.    ',' 21=Matrix ',"
                                  "'with small random entries.',/,' 20=Matrix with large ran',"
                                  "'dom entries.   ',/)");
                    write(nounit, "(' Tests performed with test threshold =',f8.2,/,/,"
                                  "' 1 = | A VR - VR W | / ( n |A| ulp ) ',/,"
                                  "' 2 = | conj-trans(A) VL - VL conj-trans(W) | /',"
                                  "' ( n |A| ulp ) ',/,' 3 = | |VR(i)| - 1 | / ulp ',/,"
                                  "' 4 = | |VL(i)| - 1 | / ulp ',/,"
                                  "' 5 = 0 if W same no matter if VR or VL computed,',"
                                  "' 1/ulp otherwise',/,"
                                  "' 6 = 0 if VR same no matter if VL computed,',"
                                  "'  1/ulp otherwise',/,"
                                  "' 7 = 0 if VL same no matter if VR computed,',"
                                  "'  1/ulp otherwise',/)"),
                        thresh;
                    ntestf = 2;
                }
                //
                for (j = 1; j <= 7; j = j + 1) {
                    if (result[j - 1] >= thresh) {
                        write(nounit, "(' N=',i5,', IWK=',i2,', seed=',4(i4,','),' type ',i2,"
                                      "', test(',i2,')=',g10.3)"),
                            n, iwk, ioldsd, jtype, j, result(j);
                    }
                }
                //
                nerrs += nfail;
                ntestt += ntest;
                //
            }
        statement_260:;
        }
    }
    //
    //     Summary
    //
    Rlasum(path, nounit, nerrs, ntestt);
    //
    //     End of Cdrvev
    //
}
