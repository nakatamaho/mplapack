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

INTEGER iMlaenv(INTEGER const ispec, const char *name, const char *opts, INTEGER const n1, INTEGER const n2, INTEGER const n3, INTEGER const n4) {
    INTEGER return_value = 0;
    char subnam[16];
    INTEGER ic = 0;
    INTEGER iz = 0;
    INTEGER i = 0;
    char c1;
    bool sname = false;
    bool cname = false;
    char c2;
    char c3;
    char c4;
    bool twostage = false;
    INTEGER nb = 0;
    INTEGER nbmin = 0;
    INTEGER nx = 0;
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Local Scalars ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. External Functions ..
    //     ..
    //     .. Executable Statements ..
    //
    switch (ispec) {
    case 1:
        goto statement_10;
    case 2:
        goto statement_10;
    case 3:
        goto statement_10;
    case 4:
        goto statement_80;
    case 5:
        goto statement_90;
    case 6:
        goto statement_100;
    case 7:
        goto statement_110;
    case 8:
        goto statement_120;
    case 9:
        goto statement_130;
    case 10:
        goto statement_140;
    case 11:
        goto statement_150;
    case 12:
        goto statement_160;
    case 13:
        goto statement_160;
    case 14:
        goto statement_160;
    case 15:
        goto statement_160;
    case 16:
        goto statement_160;
    default:
        break;
    }
    //
    //     Invalid value for ISPEC
    //
    return_value = -1;
    return return_value;
//
statement_10:
    //
    //     Convert NAME to upper case if the first character is lower case.
    //
    return_value = 1;
    subnam = name;
    ic = ichar[subnam[(1 - 1)] - 1];
    iz = ichar["Z" - 1];
    if (iz == 90 || iz == 122) {
        //
        //        ASCII character set
        //
        if (ic >= 97 && ic <= 122) {
            subnam[(1 - 1)] = fchar[(ic - 32) - 1];
            for (i = 2; i <= 6; i = i + 1) {
                ic = ichar[subnam[(i - 1) + (i - 1) * ldsubnam] - 1];
                if (ic >= 97 && ic <= 122) {
                    subnam[(i - 1) + (i - 1) * ldsubnam] = fchar[(ic - 32) - 1];
                }
            }
        }
        //
    } else if (iz == 233 || iz == 169) {
        //
        //        EBCDIC character set
        //
        if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 && ic <= 169)) {
            subnam[(1 - 1)] = fchar[(ic + 64) - 1];
            for (i = 2; i <= 6; i = i + 1) {
                ic = ichar[subnam[(i - 1) + (i - 1) * ldsubnam] - 1];
                if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 && ic <= 169)) {
                    subnam[(i - 1) + (i - 1) * ldsubnam] = fchar[(ic + 64) - 1];
                }
            }
        }
        //
    } else if (iz == 218 || iz == 250) {
        //
        //        Prime machines:  ASCII+128
        //
        if (ic >= 225 && ic <= 250) {
            subnam[(1 - 1)] = fchar[(ic - 32) - 1];
            for (i = 2; i <= 6; i = i + 1) {
                ic = ichar[subnam[(i - 1) + (i - 1) * ldsubnam] - 1];
                if (ic >= 225 && ic <= 250) {
                    subnam[(i - 1) + (i - 1) * ldsubnam] = fchar[(ic - 32) - 1];
                }
            }
        }
    }
    //
    c1 = subnam[(1 - 1)];
    sname = c1 == "S" || c1 == "D";
    cname = c1 == "C" || c1 == "Z";
    if (!(cname || sname)) {
        return return_value;
    }
    c2 = subnam[(2 - 1) + (3 - 1) * ldsubnam];
    c3 = subnam[(4 - 1) + (6 - 1) * ldsubnam];
    c4 = c3[(2 - 1) + (3 - 1) * ldc3];
    twostage = len[subnam - 1] >= 11 && subnam[(11 - 1) + (11 - 1) * ldsubnam] == "2";
    //
    switch (ispec) {
    case 1:
        goto statement_50;
    case 2:
        goto statement_60;
    case 3:
        goto statement_70;
    default:
        break;
    }
//
statement_50:
    //
    //     ISPEC = 1:  block size
    //
    //     In these examples, separate code is provided for setting NB for
    //     real and complex.  We assume that NB will take the same value in
    //     single or REAL precision.
    //
    nb = 1;
    //
    if (subnam[(2 - 1) + (6 - 1) * ldsubnam] == "LAORH") {
        //
        //        This is for *LAORHR_GETRFNP routine
        //
        if (sname) {
            nb = 32;
        } else {
            nb = 32;
        }
    } else if (c2 == "GE") {
        if (c3 == "TRF") {
            if (sname) {
                nb = 64;
            } else {
                nb = 64;
            }
        } else if (c3 == "QRF" || c3 == "RQF" || c3 == "LQF" || c3 == "QLF") {
            if (sname) {
                nb = 32;
            } else {
                nb = 32;
            }
        } else if (c3 == "QR ") {
            if (n3 == 1) {
                if (sname) {
                    //     M*N
                    if ((n1 * n2 <= 131072) || (n1 <= 8192)) {
                        nb = n1;
                    } else {
                        nb = 32768 / n2;
                    }
                } else {
                    if ((n1 * n2 <= 131072) || (n1 <= 8192)) {
                        nb = n1;
                    } else {
                        nb = 32768 / n2;
                    }
                }
            } else {
                if (sname) {
                    nb = 1;
                } else {
                    nb = 1;
                }
            }
        } else if (c3 == "LQ ") {
            if (n3 == 2) {
                if (sname) {
                    //     M*N
                    if ((n1 * n2 <= 131072) || (n1 <= 8192)) {
                        nb = n1;
                    } else {
                        nb = 32768 / n2;
                    }
                } else {
                    if ((n1 * n2 <= 131072) || (n1 <= 8192)) {
                        nb = n1;
                    } else {
                        nb = 32768 / n2;
                    }
                }
            } else {
                if (sname) {
                    nb = 1;
                } else {
                    nb = 1;
                }
            }
        } else if (c3 == "HRD") {
            if (sname) {
                nb = 32;
            } else {
                nb = 32;
            }
        } else if (c3 == "BRD") {
            if (sname) {
                nb = 32;
            } else {
                nb = 32;
            }
        } else if (c3 == "TRI") {
            if (sname) {
                nb = 64;
            } else {
                nb = 64;
            }
        }
    } else if (c2 == "PO") {
        if (c3 == "TRF") {
            if (sname) {
                nb = 64;
            } else {
                nb = 64;
            }
        }
    } else if (c2 == "SY") {
        if (c3 == "TRF") {
            if (sname) {
                if (twostage) {
                    nb = 192;
                } else {
                    nb = 64;
                }
            } else {
                if (twostage) {
                    nb = 192;
                } else {
                    nb = 64;
                }
            }
        } else if (sname && c3 == "TRD") {
            nb = 32;
        } else if (sname && c3 == "GST") {
            nb = 64;
        }
    } else if (cname && c2 == "HE") {
        if (c3 == "TRF") {
            if (twostage) {
                nb = 192;
            } else {
                nb = 64;
            }
        } else if (c3 == "TRD") {
            nb = 32;
        } else if (c3 == "GST") {
            nb = 64;
        }
    } else if (sname && c2 == "OR") {
        if (c3[(1 - 1)] == "G") {
            if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") {
                nb = 32;
            }
        } else if (c3[(1 - 1)] == "M") {
            if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") {
                nb = 32;
            }
        }
    } else if (cname && c2 == "UN") {
        if (c3[(1 - 1)] == "G") {
            if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") {
                nb = 32;
            }
        } else if (c3[(1 - 1)] == "M") {
            if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") {
                nb = 32;
            }
        }
    } else if (c2 == "GB") {
        if (c3 == "TRF") {
            if (sname) {
                if (n4 <= 64) {
                    nb = 1;
                } else {
                    nb = 32;
                }
            } else {
                if (n4 <= 64) {
                    nb = 1;
                } else {
                    nb = 32;
                }
            }
        }
    } else if (c2 == "PB") {
        if (c3 == "TRF") {
            if (sname) {
                if (n2 <= 64) {
                    nb = 1;
                } else {
                    nb = 32;
                }
            } else {
                if (n2 <= 64) {
                    nb = 1;
                } else {
                    nb = 32;
                }
            }
        }
    } else if (c2 == "TR") {
        if (c3 == "TRI") {
            if (sname) {
                nb = 64;
            } else {
                nb = 64;
            }
        } else if (c3 == "EVC") {
            if (sname) {
                nb = 64;
            } else {
                nb = 64;
            }
        }
    } else if (c2 == "LA") {
        if (c3 == "UUM") {
            if (sname) {
                nb = 64;
            } else {
                nb = 64;
            }
        }
    } else if (sname && c2 == "ST") {
        if (c3 == "EBZ") {
            nb = 1;
        }
    } else if (c2 == "GG") {
        nb = 32;
        if (c3 == "HD3") {
            if (sname) {
                nb = 32;
            } else {
                nb = 32;
            }
        }
    }
    return_value = nb;
    return return_value;
//
statement_60:
    //
    //     ISPEC = 2:  minimum block size
    //
    nbmin = 2;
    if (c2 == "GE") {
        if (c3 == "QRF" || c3 == "RQF" || c3 == "LQF" || c3 == "QLF") {
            if (sname) {
                nbmin = 2;
            } else {
                nbmin = 2;
            }
        } else if (c3 == "HRD") {
            if (sname) {
                nbmin = 2;
            } else {
                nbmin = 2;
            }
        } else if (c3 == "BRD") {
            if (sname) {
                nbmin = 2;
            } else {
                nbmin = 2;
            }
        } else if (c3 == "TRI") {
            if (sname) {
                nbmin = 2;
            } else {
                nbmin = 2;
            }
        }
    } else if (c2 == "SY") {
        if (c3 == "TRF") {
            if (sname) {
                nbmin = 8;
            } else {
                nbmin = 8;
            }
        } else if (sname && c3 == "TRD") {
            nbmin = 2;
        }
    } else if (cname && c2 == "HE") {
        if (c3 == "TRD") {
            nbmin = 2;
        }
    } else if (sname && c2 == "OR") {
        if (c3[(1 - 1)] == "G") {
            if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") {
                nbmin = 2;
            }
        } else if (c3[(1 - 1)] == "M") {
            if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") {
                nbmin = 2;
            }
        }
    } else if (cname && c2 == "UN") {
        if (c3[(1 - 1)] == "G") {
            if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") {
                nbmin = 2;
            }
        } else if (c3[(1 - 1)] == "M") {
            if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") {
                nbmin = 2;
            }
        }
    } else if (c2 == "GG") {
        nbmin = 2;
        if (c3 == "HD3") {
            nbmin = 2;
        }
    }
    return_value = nbmin;
    return return_value;
//
statement_70:
    //
    //     ISPEC = 3:  crossover point
    //
    nx = 0;
    if (c2 == "GE") {
        if (c3 == "QRF" || c3 == "RQF" || c3 == "LQF" || c3 == "QLF") {
            if (sname) {
                nx = 128;
            } else {
                nx = 128;
            }
        } else if (c3 == "HRD") {
            if (sname) {
                nx = 128;
            } else {
                nx = 128;
            }
        } else if (c3 == "BRD") {
            if (sname) {
                nx = 128;
            } else {
                nx = 128;
            }
        }
    } else if (c2 == "SY") {
        if (sname && c3 == "TRD") {
            nx = 32;
        }
    } else if (cname && c2 == "HE") {
        if (c3 == "TRD") {
            nx = 32;
        }
    } else if (sname && c2 == "OR") {
        if (c3[(1 - 1)] == "G") {
            if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") {
                nx = 128;
            }
        }
    } else if (cname && c2 == "UN") {
        if (c3[(1 - 1)] == "G") {
            if (c4 == "QR" || c4 == "RQ" || c4 == "LQ" || c4 == "QL" || c4 == "HR" || c4 == "TR" || c4 == "BR") {
                nx = 128;
            }
        }
    } else if (c2 == "GG") {
        nx = 128;
        if (c3 == "HD3") {
            nx = 128;
        }
    }
    return_value = nx;
    return return_value;
//
statement_80:
    //
    //     ISPEC = 4:  number of shifts (used by xHSEQR)
    //
    return_value = 6;
    return return_value;
//
statement_90:
    //
    //
    return_value = 2;
    return return_value;
//
statement_100:
    //
    //     ISPEC = 6:  crossover poINTEGER for SVD (used by xGELSS and xGESVD)
    //
    return_value = int(real[(min(n1 - 1) + (n2)-1) * ldreal] * 1.6e0f);
    return return_value;
//
statement_110:
    //
    //     ISPEC = 7:  number of processors (not used)
    //
    return_value = 1;
    return return_value;
//
statement_120:
    //
    //     ISPEC = 8:  crossover poINTEGER for multishift (used by xHSEQR)
    //
    return_value = 50;
    return return_value;
//
statement_130:
    //
    //     ISPEC = 9:  maximum size of the subproblems at the bottom of the
    //                 computation tree in the divide-and-conquer algorithm
    //                 (used by xGELSD and xGESDD)
    //
    return_value = 25;
    return return_value;
//
statement_140:
    //
    //     ISPEC = 10: ieee and infinity NaN arithmetic can be trusted not to trap
    //
    //     iMlaenv = 0
    return_value = 1;
    if (return_value == 1) {
        return_value = ieeeck[(0.0f - 1) * ldieeeck];
    }
    return return_value;
//
statement_150:
    //
    //     ISPEC = 11: ieee infinity arithmetic can be trusted not to trap
    //
    //     iMlaenv = 0
    return_value = 1;
    if (return_value == 1) {
        return_value = ieeeck[(0 - 1) + (0.0f - 1) * ldieeeck];
    }
    return return_value;
//
statement_160:
    //
    //     12 <= ISPEC <= 16: xHSEQR or related subroutines.
    //
    return_value = iparmq[(ispec - 1) + (name - 1) * ldiparmq];
    return return_value;
    //
    //     End of iMlaenv
    //
}
