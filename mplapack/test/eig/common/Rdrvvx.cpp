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

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>

using namespace std;
using std::regex;
using std::regex_replace;

void Rdrvvx(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *dotype, INTEGER *iseed, REAL const thresh, INTEGER const niunit, INTEGER const nounit, REAL *a, INTEGER const lda, REAL *h, REAL *wr, REAL *wi, REAL *wr1, REAL *wi1, REAL *vl, INTEGER const ldvl, REAL *vr, INTEGER const ldvr, REAL *lre, INTEGER const ldlre, REAL *rcondv, REAL *rcndv1, REAL *rcdvin, REAL *rconde, REAL *rcnde1, REAL *rcdein, REAL *scale, REAL *scale1, REAL *result, REAL *work, INTEGER const nwork, INTEGER *iwork, INTEGER &info) {

    INTEGER ldh = lda;
    common cmn;
    common_read read(cmn);
    common_write write(cmn);
    const INTEGER maxtyp = 21;
    INTEGER ktype[21] = {1, 2, 3, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9};
    INTEGER kmagn[21] = {1, 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 1, 2, 3};
    INTEGER kmode[21] = {0, 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1};
    INTEGER kconds[21] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0};
    char bal[1 * 4] = {'N', 'P', 'S', 'B'};
    char buf[1024];
    double dtmp;
    double dtmp1, dtmp2, dtmp3, dtmp4;
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
    INTEGER iinfo = 0;
    REAL cond = 0.0;
    INTEGER jcol = 0;
    REAL conds = 0.0;
    char adumma[1];
    INTEGER idumma[1];
    INTEGER iwk = 0;
    INTEGER nnwork = 0;
    INTEGER ibal = 0;
    char balanc;
    INTEGER ntest = 0;
    INTEGER nfail = 0;
    INTEGER i = 0;
    static const char *format_9995 = "(' Tests performed with test threshold =',a,/,/,"
                                     "' 1 = | A VR - VR W | / ( n |A| ulp ) ',/,"
                                     "' 2 = | transpose(A) VL - VL W | / ( n |A| ulp ) ',/,"
                                     "' 3 = | |VR(i)| - 1 | / ulp ',/,' 4 = | |VL(i)| - 1 | / ulp ',/,"
                                     "' 5 = 0 if W same no matter if VR or VL computed,',' 1/ulp otherwise',/,"
                                     "' 6 = 0 if VR same no matter what else computed,','  1/ulp otherwise',/,"
                                     "' 7 = 0 if VL same no matter what else computed,','  1/ulp otherwise',/,"
                                     "' 8 = 0 if RCONDV same no matter what else computed,',"
                                     "'  1/ulp otherwise',/,"
                                     "' 9 = 0 if SCALE, ILO, IHI, ABNRM same no matter what else',"
                                     "' computed,  1/ulp otherwise',/,"
                                     "' 10 = | RCONDV - RCONDV(precomputed) | / cond(RCONDV),',/,"
                                     "' 11 = | RCONDE - RCONDE(precomputed) | / cond(RCONDE),')";
    static const char *format_9996 = "(' 19=Matrix with random O(1) entries.    ',' 21=Matrix ',"
                                     "'with small random entries.',/,' 20=Matrix with large ran',"
                                     "'dom entries.   ',' 22=Matrix read from input file',/)";
    static const char *format_9997 = "(' Dense, Non-Symmetric Matrices:',/,'  9=Well-cond., ev',"
                                     "'enly spaced eigenvals.',' 14=Ill-cond., geomet. spaced e','igenals.',/,"
                                     "' 10=Well-cond., geom. spaced eigenvals. ',"
                                     "' 15=Ill-conditioned, clustered e.vals.',/,' 11=Well-cond',"
                                     "'itioned, clustered e.vals. ',' 16=Ill-cond., random comp','lex ',/,"
                                     "' 12=Well-cond., random complex ','         ',"
                                     "' 17=Ill-cond., large rand. complx ',/,' 13=Ill-condi',"
                                     "'tioned, evenly spaced.     ',' 18=Ill-cond., small rand.',' complx ')";
    static const char *format_9998 = "(/,' Special Matrices:',/,'  1=Zero matrix.             ','           ',"
                                     "'  5=Diagonal: geometr. spaced entries.',/,"
                                     "'  2=Identity matrix.                    ','  6=Diagona',"
                                     "'l: clustered entries.',/,'  3=Transposed Jordan block.  ','          ',"
                                     "'  7=Diagonal: large, evenly spaced.',/,'  ',"
                                     "'4=Diagonal: evenly spaced entries.    ','  8=Diagonal: s',"
                                     "'mall, evenly spaced.')";
    static const char *format_9999 = "(/,1x,a3,' -- Real Eigenvalue-Eigenvector Decomposition',"
                                     "' Expert Driver',/,' Matrix types (see Rdrvvx for details): ')";
    //
    //
    //     .. Executable Statements ..
    //
    path[0] = 'R';
    path[1] = 'V';
    path[2] = 'X';
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
    //
    //     problems
    //
    nmax = 12;
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
    } else if (lda < 1 || lda < nmax) {
        info = -10;
    } else if (ldvl < 1 || ldvl < nmax) {
        info = -17;
    } else if (ldvr < 1 || ldvr < nmax) {
        info = -19;
    } else if (ldlre < 1 || ldlre < nmax) {
        info = -21;
    } else if (6 * nmax + 2 * pow2(nmax) > nwork) {
        info = -32;
    }
    //
    if (info != 0) {
        Mxerbla("Rdrvvx", -info);
        return;
    }
    //
    //     If nothing to do check on NIUNIT
    //
    if (nsizes == 0 || ntypes == 0) {
        goto statement_160;
    }
    //
    //     More Important constants
    //
    unfl = Rlamch("Safe minimum");
    ovfl = one / unfl;
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
                goto statement_140;
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
            Rlaset("Full", lda, n, zero, zero, a, lda);
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
                    a[(jcol - 1) + (jcol - 1) * lda] = anorm;
                }
                //
            } else if (itype == 3) {
                //
                //              Jordan Block
                //
                for (jcol = 1; jcol <= n; jcol = jcol + 1) {
                    a[(jcol - 1) + (jcol - 1) * lda] = anorm;
                    if (jcol > 1) {
                        a[(jcol - 1) + ((jcol - 1) - 1) * lda] = one;
                    }
                }
                //
            } else if (itype == 4) {
                //
                //              Diagonal Matrix, [Eigen]values Specified
                //
                Rlatms(n, n, "S", iseed, "S", work, imode, cond, anorm, 0, 0, "N", a, lda, &work[(n + 1) - 1], iinfo);
                //
            } else if (itype == 5) {
                //
                //              Symmetric, eigenvalues specified
                //
                Rlatms(n, n, "S", iseed, "S", work, imode, cond, anorm, n, n, "N", a, lda, &work[(n + 1) - 1], iinfo);
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
                adumma[1 - 1] = ' ';
                Rlatme(n, "S", iseed, work, imode, cond, one, adumma, "T", "T", "T", &work[(n + 1) - 1], 4, conds, n, n, anorm, a, lda, &work[(2 * n + 1) - 1], iinfo);
                //
            } else if (itype == 7) {
                //
                //              Diagonal, random eigenvalues
                //
                Rlatmr(n, n, "S", iseed, "S", work, 6, one, one, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, 0, 0, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 8) {
                //
                //              Symmetric, random eigenvalues
                //
                Rlatmr(n, n, "S", iseed, "S", work, 6, one, one, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, n, n, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else if (itype == 9) {
                //
                //              General, random eigenvalues
                //
                Rlatmr(n, n, "S", iseed, "N", work, 6, one, one, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, n, n, zero, anorm, "NO", a, lda, iwork, iinfo);
                if (n >= 4) {
                    Rlaset("Full", 2, n, zero, zero, a, lda);
                    Rlaset("Full", n - 3, 1, zero, zero, &a[(3 - 1)], lda);
                    Rlaset("Full", n - 3, 2, zero, zero, &a[(3 - 1) + ((n - 1) - 1) * lda], lda);
                    Rlaset("Full", 1, n, zero, zero, &a[(n - 1)], lda);
                }
                //
            } else if (itype == 10) {
                //
                //              Triangular, random eigenvalues
                //
                Rlatmr(n, n, "S", iseed, "N", work, 6, one, one, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, n, 0, zero, anorm, "NO", a, lda, iwork, iinfo);
                //
            } else {
                //
                iinfo = 1;
            }
            //
            if (iinfo != 0) {
                write(nounit, "(' Rdrvvx: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,', JTYPE=',i6,"
                              "', ISEED=(',3(i5,','),i5,')')"),
                    "Generator", iinfo, n, jtype, ioldsd;
                info = abs(iinfo);
                return;
            }
        //
        statement_90:
            //
            //           Test for minimal and generous workspace
            //
            for (iwk = 1; iwk <= 3; iwk = iwk + 1) {
                if (iwk == 1) {
                    nnwork = 3 * n;
                } else if (iwk == 2) {
                    nnwork = 6 * n + pow2(n);
                } else {
                    nnwork = 6 * n + 2 * pow2(n);
                }
                nnwork = max(nnwork, 1);
                //
                //              Test for all balancing options
                //
                for (ibal = 1; ibal <= 4; ibal = ibal + 1) {
                    balanc = bal[ibal - 1];
                    //
                    //                 Perform tests
                    //
                    Rget23(false, &balanc, jtype, thresh, ioldsd, nounit, n, a, lda, h, wr, wi, wr1, wi1, vl, ldvl, vr, ldvr, lre, ldlre, rcondv, rcndv1, rcdvin, rconde, rcnde1, rcdein, scale, scale1, result, work, nnwork, iwork, info);
                    //
                    //                 Check for RESULT(j) > THRESH
                    //
                    ntest = 0;
                    nfail = 0;
                    for (j = 1; j <= 9; j = j + 1) {
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
                        sprintnum_short(buf, thresh);
                        write(nounit, format_9999), path;
                        write(nounit, format_9998);
                        write(nounit, format_9997);
                        write(nounit, format_9996);
                        write(nounit, format_9995), buf;
                        ntestf = 2;
                    }
                    //
                    for (j = 1; j <= 9; j = j + 1) {
                        if (result[j - 1] >= thresh) {
                            sprintnum_short(buf, result[j - 1]);
                            write(nounit, "(' BALANC=''',a1,''',N=',i4,',IWK=',i1,', seed=',4(i4,','),"
                                          "' type ',i2,', test(',i2,')=',a)"),
                                balanc, n, iwk, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3], jtype, j, buf;
                        }
                    }
                    //
                    nerrs += nfail;
                    ntestt += ntest;
                    //
                }
            }
        statement_140:;
        }
    }
//
statement_160:
    //
    //     Read in data from file to check accuracy of condition estimation.
    //     Assume input eigenvalues are sorted lexicographically (increasing
    //     by real part, then decreasing by imaginary part)
    //
    jtype = 0;
    string str;
    istringstream iss;
    while (1) {
        getline(cin, str);
        iss.clear();
        iss.str(str);
        iss >> n;
        //
        //     Read input data until N=0
        //
        if (n == 0)
            break;
        jtype++;
        iseed[1 - 1] = jtype;
        for (i = 1; i <= n; i = i + 1) {
            getline(cin, str);
            string _r = regex_replace(str, regex("D\\+"), "e+");
            str = regex_replace(_r, regex("D\\-"), "e-");
            iss.clear();
            iss.str(str);
            for (j = 1; j <= n; j = j + 1) {
                iss >> dtmp;
                a[(i - 1) + (j - 1) * lda] = dtmp;
            }
        }
        for (i = 1; i <= n; i = i + 1) {
            getline(cin, str);
            string _r = regex_replace(str, regex("D\\+"), "e+");
            str = regex_replace(_r, regex("D\\-"), "e-");
            iss.clear();
            iss.str(str);
            iss >> dtmp1;
            iss >> dtmp2;
            iss >> dtmp3;
            iss >> dtmp4;
            wr1[i - 1] = dtmp1;
            wi1[i - 1] = dtmp2;
            rcdein[i - 1] = dtmp3;
            rcdvin[i - 1] = dtmp4;
        }
        Rget23(true, "N", 22, thresh, iseed, nounit, n, a, lda, h, wr, wi, wr1, wi1, vl, ldvl, vr, ldvr, lre, ldlre, rcondv, rcndv1, rcdvin, rconde, rcnde1, rcdein, scale, scale1, result, work, 6 * n + 2 * pow2(n), iwork, info);
        //
        //     Check for RESULT(j) > THRESH
        //
        ntest = 0;
        nfail = 0;
        for (j = 1; j <= 11; j = j + 1) {
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
            sprintnum_short(buf, thresh);
            write(nounit, format_9999), path;
            write(nounit, format_9998);
            write(nounit, format_9997);
            write(nounit, format_9996);
            write(nounit, format_9995), buf;
            ntestf = 2;
        }
        //
        for (j = 1; j <= 11; j = j + 1) {
            if (result[j - 1] >= thresh) {
                sprintnum_short(buf, result[j - 1]);
                write(nounit, "(' N=',i5,', input example =',i3,',  test(',i2,')=',a)"), n, jtype, j, buf;
            }
        }
        //
        nerrs += nfail;
        ntestt += ntest;
    }
statement_220:
    //
    //     Summary
    //
    Rlasum(path, nounit, nerrs, ntestt);
    //
    //     End of Rdrvvx
    //
}
