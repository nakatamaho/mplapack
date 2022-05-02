/*
 * Copyright (c) 2021-2022
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

void Cdrvsx(INTEGER const nsizes, INTEGER *nn, INTEGER const ntypes, bool *dotype, INTEGER *iseed, REAL const thresh, INTEGER const niunit, INTEGER const nounit, COMPLEX *a, INTEGER const lda, COMPLEX *h, COMPLEX *ht, COMPLEX *w, COMPLEX *wt, COMPLEX *wtmp, COMPLEX *vs, INTEGER const ldvs, COMPLEX *vs1, REAL *result, COMPLEX *work, INTEGER const lwork, REAL *rwork, bool *bwork, INTEGER &info) {
    INTEGER ldh = lda;
    INTEGER ldht = lda;
    INTEGER ldvs1 = ldvs;
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
    double dtmp1, dtmp2;
    double dtmp_r, dtmp_i;
    std::complex<double> ctmp;
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
    REAL rcdein = 0.0;
    REAL rcdvin = 0.0;
    INTEGER nslct = 0;
    INTEGER islct[20];
    INTEGER ntest = 0;
    INTEGER nfail = 0;
    INTEGER isrt = 0;
    INTEGER i = 0;
    static const char *format_9994 = "(' 7 = 0 if T in Schur form (sort), ','  1/ulp otherwise',/,"
                                     "' 8 = | A - VS T transpose(VS) | / ( n |A| ulp ) (sort)',/,"
                                     "' 9 = | I - VS transpose(VS) | / ( n ulp ) (sort) ',/,"
                                     "' 10 = 0 if W are eigenvalues of T (sort),','  1/ulp otherwise',/,"
                                     "' 11 = 0 if T same no matter what else computed (sort),',"
                                     "'  1/ulp otherwise',/,' 12 = 0 if W same no matter what else computed ',"
                                     "'(sort), 1/ulp otherwise',/,"
                                     "' 13 = 0 if sorting successful, 1/ulp otherwise',/,"
                                     "' 14 = 0 if RCONDE same no matter what else computed,',"
                                     "' 1/ulp otherwise',/,"
                                     "' 15 = 0 if RCONDv same no matter what else computed,',"
                                     "' 1/ulp otherwise',/,"
                                     "' 16 = | RCONDE - RCONDE(precomputed) | / cond(RCONDE),',/,"
                                     "' 17 = | RCONDV - RCONDV(precomputed) | / cond(RCONDV),')";
    static const char *format_9995 = "(' Tests performed with test threshold =',a,/,"
                                     "' ( A denotes A on input and T denotes A on output)',/,/,"
                                     "' 1 = 0 if T in Schur form (no sort), ','  1/ulp otherwise',/,"
                                     "' 2 = | A - VS T transpose(VS) | / ( n |A| ulp ) (no sort)',/,"
                                     "' 3 = | I - VS transpose(VS) | / ( n ulp ) (no sort) ',/,"
                                     "' 4 = 0 if W are eigenvalues of T (no sort),','  1/ulp otherwise',/,"
                                     "' 5 = 0 if T same no matter if VS computed (no sort),',"
                                     "'  1/ulp otherwise',/,"
                                     "' 6 = 0 if W same no matter if VS computed (no sort)',"
                                     "',  1/ulp otherwise')";
    static const char *format_9996 = "(' 19=Matrix with random O(1) entries.    ',' 21=Matrix ',"
                                     "'with small random entries.',/,' 20=Matrix with large ran',"
                                     "'dom entries.   ',/)";
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
    static const char *format_9999 = "(/,1x,a3,' -- Complex Schur Form Decomposition Expert ','Driver',/,"
                                     "' Matrix types (see Cdrvsx for details): ')";
    //
    //
    path[0] = 'C';
    path[1] = 'S';
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
    nmax = 8;
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
    } else if (niunit <= 0) {
        info = -7;
    } else if (nounit <= 0) {
        info = -8;
    } else if (lda < 1 || lda < nmax) {
        info = -10;
    } else if (ldvs < 1 || ldvs < nmax) {
        info = -20;
    } else if (max(3 * nmax, 2 * nmax * nmax) > lwork) {
        info = -24;
    }
    //
    if (info != 0) {
        Mxerbla("Cdrvsx", -info);
        return;
    }
    //
    //     If nothing to do check on NIUNIT
    //
    if (nsizes == 0 || ntypes == 0) {
        goto statement_150;
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
                goto statement_130;
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
            if (itype == 1) {
                //
                //              Zero
                //
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
                        a[(jcol - 1) + ((jcol - 1) - 1) * lda] = cone;
                    }
                }
                //
            } else if (itype == 4) {
                //
                //              Diagonal Matrix, [Eigen]values Specified
                //
                Clatms(n, n, "S", iseed, "H", rwork, imode, cond, anorm, 0, 0, "N", a, lda, &work[(n + 1) - 1], iinfo);
                //
            } else if (itype == 5) {
                //
                //              Symmetric, eigenvalues specified
                //
                Clatms(n, n, "S", iseed, "H", rwork, imode, cond, anorm, n, n, "N", a, lda, &work[(n + 1) - 1], iinfo);
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
                Clatme(n, "D", iseed, work, imode, cond, cone, "T", "T", "T", rwork, 4, conds, n, n, anorm, a, lda, &work[(2 * n + 1) - 1], iinfo);
                //
            } else if (itype == 7) {
                //
                //              Diagonal, random eigenvalues
                //
                Clatmr(n, n, "D", iseed, "N", work, 6, one, cone, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, 0, 0, zero, anorm, "NO", a, lda, idumma, iinfo);
                //
            } else if (itype == 8) {
                //
                //              Symmetric, random eigenvalues
                //
                Clatmr(n, n, "D", iseed, "H", work, 6, one, cone, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, n, n, zero, anorm, "NO", a, lda, idumma, iinfo);
                //
            } else if (itype == 9) {
                //
                //              General, random eigenvalues
                //
                Clatmr(n, n, "D", iseed, "N", work, 6, one, cone, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, n, n, zero, anorm, "NO", a, lda, idumma, iinfo);
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
                Clatmr(n, n, "D", iseed, "N", work, 6, one, cone, "T", "N", &work[(n + 1) - 1], 1, one, &work[(2 * n + 1) - 1], 1, one, "N", idumma, n, 0, zero, anorm, "NO", a, lda, idumma, iinfo);
                //
            } else {
                //
                iinfo = 1;
            }
            //
            if (iinfo != 0) {
                write(nounit, "(' Cdrvsx: ',a,' returned INFO=',i6,'.',/,9x,'N=',i6,', JTYPE=',i6,"
                              "', ISEED=(',3(i5,','),i5,')')"),
                    "Generator", iinfo, n, jtype, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3];
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
                    nnwork = max((INTEGER)2 * n, n * (n + 1) / 2);
                }
                nnwork = max(nnwork, (INTEGER)1);
                //
                Cget24(false, jtype, thresh, ioldsd, nounit, n, a, lda, h, ht, w, wt, wtmp, vs, ldvs, vs1, rcdein, rcdvin, nslct, islct, 0, result, work, nnwork, rwork, bwork, info);
                //
                //              Check for RESULT(j) > THRESH
                //
                ntest = 0;
                nfail = 0;
                for (j = 1; j <= 15; j = j + 1) {
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
                    write(nounit, format_9994);
                    ntestf = 2;
                }
                //
                for (j = 1; j <= 15; j = j + 1) {
                    if (result[j - 1] >= thresh) {
                        sprintnum_short(buf, result[j - 1]);
                        write(nounit, "(' N=',i5,', IWK=',i2,', seed=',4(i4,','),' type ',i2,"
                                      "', test(',i2,')=',a)"),
                            n, iwk, ioldsd[0], ioldsd[1], ioldsd[2], ioldsd[3], jtype, j, buf;
                    }
                }
                //
                nerrs += nfail;
                ntestt += ntest;
                //
            }
        statement_130:;
        }
    }
//
statement_150:
    //
    //     Read in data from file to check accuracy of condition estimation
    //     Read input data until N=0
    //
    jtype = 0;
    string str;
    istringstream iss;
    while (1) {
        getline(cin, str);
        iss.clear();
        iss.str(str);
        iss >> n;
        iss >> nslct;
        iss >> isrt;
//        printf("# n = %d, nslct = %d, isrt = %d\n", (int)n, int(nslct), int(isrt));
        if (n == 0)
            break;
        jtype++;
        iseed[1 - 1] = jtype;
        getline(cin, str);
        string __r = regex_replace(str, regex(","), " ");
        string _r = regex_replace(__r, regex("\\)"), " ");
        string r = regex_replace(_r, regex("\\("), " ");
        str = regex_replace(r, regex("D"), "e");
        iss.clear();
        iss.str(str);
        for (i = 1; i <= nslct; i = i + 1) {
            iss >> islct[i - 1];
        }

//        for (i = 1; i <= nslct; i = i + 1) {
//            cout << "islct: " << islct[i - 1] << "\n";
//        }

        for (i = 1; i <= n; i = i + 1) {
            for (j = 1; j <= n; j = j + 1) {
                getline(cin, str);
                string __r = regex_replace(str, regex(","), " ");
                string _r = regex_replace(__r, regex("\\)"), " ");
                string r = regex_replace(_r, regex("\\("), " ");
                str = regex_replace(r, regex("D"), "e");
                iss.clear();
                iss.str(str);
                iss >> dtmp_r;
                iss >> dtmp_i;
                a[(i - 1) + (j - 1) * lda] = COMPLEX(dtmp_r, dtmp_i);
            }
        }
        // printf("a="); printmat(n, n, a, lda); printf("\n");
        getline(cin, str);
        iss.clear();
        iss.str(str);
        iss >> dtmp1;
        iss >> dtmp2;
        rcdein = dtmp1;
        rcdvin = dtmp2;
        //
        Cget24(true, 22, thresh, iseed, nounit, n, a, lda, h, ht, w, wt, wtmp, vs, ldvs, vs1, rcdein, rcdvin, nslct, islct, isrt, result, work, lwork, rwork, bwork, info);
        //
        //     Check for RESULT(j) > THRESH
        //
        ntest = 0;
        nfail = 0;
        for (j = 1; j <= 17; j = j + 1) {
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
            write(nounit, format_9994);
            ntestf = 2;
        }
        for (j = 1; j <= 17; j = j + 1) {
            if (result[j - 1] >= thresh) {
                sprintnum_short(buf, result[j - 1]);
                write(nounit, "(' N=',i5,', input example =',i3,',  test(',i2,')=',a)"), n, jtype, j, buf;
            }
        }
        //
        nerrs += nfail;
        ntestt += ntest;
    }
    //
    //     Summary
    //
    Rlasum(path, nounit, nerrs, ntestt);
    //
    //     End of Cdrvsx
    //
}
