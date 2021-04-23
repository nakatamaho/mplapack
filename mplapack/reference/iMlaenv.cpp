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
#include <string.h>

#define subnamlen 16

INTEGER iMlaenv(INTEGER const ispec, const char *name, const char *opts, INTEGER const n1, INTEGER const n2, INTEGER const n3, INTEGER const n4) {
    INTEGER return_value = 0;
    char subnam[subnamlen];
    memset(subnam, '\0', sizeof(subnam));
    INTEGER ic = 0;
    INTEGER i = 0;
    bool sname = false;
    bool cname = false;
    char c1[1];
    char c2[2];
    char c3[3];
    char c4[2];
    bool twostage = false;
    INTEGER nb = 0;
    INTEGER nbmin = 0;
    INTEGER nx = 0;
    INTEGER name_len;
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
        goto L10;
    case 2:
        goto L10;
    case 3:
        goto L10;
    case 4:
        goto L80;
    case 5:
        goto L90;
    case 6:
        goto L100;
    case 7:
        goto L110;
    case 8:
        goto L120;
    case 9:
        goto L130;
    case 10:
        goto L140;
    case 11:
        goto L150;
    case 12:
        goto L160;
    case 13:
        goto L160;
    case 14:
        goto L160;
    case 15:
        goto L160;
    case 16:
        goto L160;
    }
    //
    //     Invalid value for ISPEC
    //
    return_value = -1;
    return return_value;
    //
L10:
    //
    //     Convert NAME to upper case if the first character is lower case.
    //
    return_value = 1;
    name_len = min((int)strlen(name), subnamlen);
    strncpy(subnam, name, name_len);
    ic = *subnam;

    for (int i=0 ; i < strlen(subnam) ; i++)
    {
	subnam[i] = toupper(subnam[i]);
    }
    *c1 = *subnam;
    sname = *c1 == 'R';
    cname = *c1 == 'C';
    if (!(cname || sname)) {
        return return_value;
    }
    strncpy(c2, subnam + 1, 2);
    strncpy(c3, subnam + 3, 3);
    strncpy(c4, c3 + 1, 2);
    twostage = strlen(subnam) >= 11 && subnam[10] == '2';

    switch (ispec) {
    case 1:
        goto L50;
    case 2:
        goto L60;
    case 3:
        goto L70;
    }
L50:
    //
    //     ISPEC = 1:  block size
    //
    //     In these examples, separate code is provided for setting NB for
    //     real and complex.  We assume that NB will take the same value in
    //     single or REAL precision.
    //
    nb = 1;
    //
    if (strncmp(subnam + 1, "LAORH", 5) == 0) {
        //
        //        This is for *LAORHR_GETRFNP routine
        //
        if (sname) {
            nb = 32;
        } else {
            nb = 32;
        }
    } else if (strncmp(c2, "GE", 2) == 0) {
        if (strncmp(c3, "TRF", 3) == 0) {
            if (sname) {
                nb = 64;
            } else {
                nb = 64;
            }
        } else if (strncmp(c3, "QRF", 3) == 0 || strncmp(c3, "RQF", 3) == 0 || strncmp(c3, "LQF", 3) == 0 || strncmp(c3, "QLF", 3) == 0) {
            if (sname) {
                nb = 32;
            } else {
                nb = 32;
            }
        } else if (strncmp(c3, "QR ", 3) == 0) {
            if (n3 == 1) {
                if (sname) {
                    //    M*N
                    if (n1 * n2 <= 131072 || n1 <= 8192) {
                        nb = n1;
                    } else {
                        nb = 32768 / n2;
                    }
                } else {
                    if (n1 * n2 <= 131072 || n1 <= 8192) {
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
        } else if (strncmp(c3, "LQ ", 3) == 0) {
            if (n3 == 2) {
                if (sname) {
                    //    M*N
                    if (n1 * n2 <= 131072 || n1 <= 8192) {
                        nb = n1;
                    } else {
                        nb = 32768 / n2;
                    }
                } else {
                    if (n1 * n2 <= 131072 || n1 <= 8192) {
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
        } else if (strncmp(c3, "HRD", 3) == 0) {
            if (sname) {
                nb = 32;
            } else {
                nb = 32;
            }
        } else if (strncmp(c3, "BRD", 3) == 0) {
            if (sname) {
                nb = 32;
            } else {
                nb = 32;
            }
        } else if (strncmp(c3, "TRI", 3) == 0) {
            if (sname) {
                nb = 64;
            } else {
                nb = 64;
            }
        }
    } else if (strncmp(c2, "PO", 2) == 0) {
        if (strncmp(c3, "TRF", 3) == 0) {
            if (sname) {
                nb = 64;
            } else {
                nb = 64;
            }
        }
    } else if (strncmp(c2, "SY", 2) == 0) {
        if (strncmp(c3, "TRF", 3) == 0) {
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
        } else if (sname && strncmp(c3, "TRD", 3) == 0) {
            nb = 32;
        } else if (sname && strncmp(c3, "GST", 3) == 0) {
            nb = 64;
        }
    } else if (cname && strncmp(c2, "HE", 2) == 0) {
        if (strncmp(c3, "TRF", 3) == 0) {
            if (twostage) {
                nb = 192;
            } else {
                nb = 64;
            }
        } else if (strncmp(c3, "TRD", 3) == 0) {
            nb = 32;
        } else if (strncmp(c3, "GST", 3) == 0) {
            nb = 64;
        }
    } else if (sname && strncmp(c2, "OR", 2) == 0) {
        if (*c3 == 'G') {
            if (strncmp(c4, "QR", 2) == 0 || strncmp(c4, "RQ", 2) == 0 || strncmp(c4, "LQ", 2) == 0 || strncmp(c4, "QL", 2) == 0 || strncmp(c4, "HR", 2) == 0 || strncmp(c4, "TR", 2) == 0 || strncmp(c4, "BR", 2) == 0) {
                nb = 32;
            }
        } else if (*c3 == 'M') {
            if (strncmp(c4, "QR", 2) == 0 || strncmp(c4, "RQ", 2) == 0 || strncmp(c4, "LQ", 2) == 0 || strncmp(c4, "QL", 2) == 0 || strncmp(c4, "HR", 2) == 0 || strncmp(c4, "TR", 2) == 0 || strncmp(c4, "BR", 2) == 0) {
                nb = 32;
            }
        }
    } else if (cname && strncmp(c2, "UN", 2) == 0) {
        if (*c3 == 'G') {
            if (strncmp(c4, "QR", 2) == 0 || strncmp(c4, "RQ", 2) == 0 || strncmp(c4, "LQ", 2) == 0 || strncmp(c4, "QL", 2) == 0 || strncmp(c4, "HR", 2) == 0 || strncmp(c4, "TR", 2) == 0 || strncmp(c4, "BR", 2) == 0) {
                nb = 32;
            }
        } else if (*c3 == 'M') {
            if (strncmp(c4, "QR", 2) == 0 || strncmp(c4, "RQ", 2) == 0 || strncmp(c4, "LQ", 2) == 0 || strncmp(c4, "QL", 2) == 0 || strncmp(c4, "HR", 2) == 0 || strncmp(c4, "TR", 2) == 0 || strncmp(c4, "BR", 2) == 0) {
                nb = 32;
            }
        }
    } else if (strncmp(c2, "GB", 2) == 0) {
        if (strncmp(c3, "TRF", 3) == 0) {
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
    } else if (strncmp(c2, "PB", 2) == 0) {
        if (strncmp(c3, "TRF", 3) == 0) {
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
    } else if (strncmp(c2, "TR", 2) == 0) {
        if (strncmp(c3, "TRI", 3) == 0) {
            if (sname) {
                nb = 64;
            } else {
                nb = 64;
            }
        } else if (strncmp(c3, "EVC", 3) == 0) {
            if (sname) {
                nb = 64;
            } else {
                nb = 64;
            }
        }
    } else if (strncmp(c2, "LA", 2) == 0) {
        if (strncmp(c3, "UUM", 3) == 0) {
            if (sname) {
                nb = 64;
            } else {
                nb = 64;
            }
        }
    } else if (sname && strncmp(c2, "ST", 2) == 0) {
        if (strncmp(c3, "EBZ", 3) == 0) {
            nb = 1;
        }
    } else if (strncmp(c2, "GG", 2) == 0) {
        nb = 32;
        if (strncmp(c3, "HD3", 3) == 0) {
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
L60:
    //
    //     ISPEC = 2:  minimum block size
    //
    nbmin = 2;
    if (strncmp(c2, "GE", 2) == 0) {
        if (strncmp(c3, "QRF", 3) == 0 || strncmp(c3, "RQF", 3) == 0 || strncmp(c3, "LQF", 3) == 0 || strncmp(c3, "QLF", 3) == 0) {
            if (sname) {
                nbmin = 2;
            } else {
                nbmin = 2;
            }
        } else if (strncmp(c3, "HRD", 3) == 0) {
            if (sname) {
                nbmin = 2;
            } else {
                nbmin = 2;
            }
        } else if (strncmp(c3, "BRD", 3) == 0) {
            if (sname) {
                nbmin = 2;
            } else {
                nbmin = 2;
            }
        } else if (strncmp(c3, "TRI", 3) == 0) {
            if (sname) {
                nbmin = 2;
            } else {
                nbmin = 2;
            }
        }
    } else if (strncmp(c2, "SY", 2) == 0) {
        if (strncmp(c3, "TRF", 3) == 0) {
            if (sname) {
                nbmin = 8;
            } else {
                nbmin = 8;
            }
        } else if (sname && strncmp(c3, "TRD", 3) == 0) {
            nbmin = 2;
        }
    } else if (cname && strncmp(c2, "HE", 2) == 0) {
        if (strncmp(c3, "TRD", 3) == 0) {
            nbmin = 2;
        }
    } else if (sname && strncmp(c2, "OR", 2) == 0) {
        if (*c3 == 'G') {
            if (strncmp(c4, "QR", 2) == 0 || strncmp(c4, "RQ", 2) == 0 || strncmp(c4, "LQ", 2) == 0 || strncmp(c4, "QL", 2) == 0 || strncmp(c4, "HR", 2) == 0 || strncmp(c4, "TR", 2) == 0 || strncmp(c4, "BR", 2) == 0) {
                nbmin = 2;
            }
        } else if (*c3 == 'M') {
            if (strncmp(c4, "QR", 2) == 0 || strncmp(c4, "RQ", 2) == 0 || strncmp(c4, "LQ", 2) == 0 || strncmp(c4, "QL", 2) == 0 || strncmp(c4, "HR", 2) == 0 || strncmp(c4, "TR", 2) == 0 || strncmp(c4, "BR", 2) == 0) {
                nbmin = 2;
            }
        }
    } else if (cname && strncmp(c2, "UN", 2) == 0) {
        if (*c3 == 'G') {
            if (strncmp(c4, "QR", 2) == 0 || strncmp(c4, "RQ", 2) == 0 || strncmp(c4, "LQ", 2) == 0 || strncmp(c4, "QL", 2) == 0 || strncmp(c4, "HR", 2) == 0 || strncmp(c4, "TR", 2) == 0 || strncmp(c4, "BR", 2) == 0) {
                nbmin = 2;
            }
        } else if (*c3 == 'M') {
            if (strncmp(c4, "QR", 2) == 0 || strncmp(c4, "RQ", 2) == 0 || strncmp(c4, "LQ", 2) == 0 || strncmp(c4, "QL", 2) == 0 || strncmp(c4, "HR", 2) == 0 || strncmp(c4, "TR", 2) == 0 || strncmp(c4, "BR", 2) == 0) {
                nbmin = 2;
            }
        }
    } else if (strncmp(c2, "GG", 2) == 0) {
        nbmin = 2;
        if (strncmp(c3, "HD3", 3) == 0) {
            nbmin = 2;
        }
    }
    return_value = nbmin;
    return return_value;
//
L70:
    //
    //     ISPEC = 3:  crossover point
    //
    nx = 0;
    if (strncmp(c2, "GE", 2) == 0) {
        if (strncmp(c3, "QRF", 3) == 0 || strncmp(c3, "RQF", 3) == 0 || strncmp(c3, "LQF", 3) == 0 || strncmp(c3, "QLF", 3) == 0) {
            if (sname) {
                nx = 128;
            } else {
                nx = 128;
            }
        } else if (strncmp(c3, "HRD", 3) == 0) {
            if (sname) {
                nx = 128;
            } else {
                nx = 128;
            }
        } else if (strncmp(c3, "BRD", 3) == 0) {
            if (sname) {
                nx = 128;
            } else {
                nx = 128;
            }
        }
    } else if (strncmp(c2, "SY", 2) == 0) {
        if (sname && strncmp(c3, "TRD", 3) == 0) {
            nx = 32;
        }
    } else if (cname && strncmp(c2, "HE", 2) == 0) {
        if (strncmp(c3, "TRD", 3) == 0) {
            nx = 32;
        }
    } else if (sname && strncmp(c2, "OR", 2) == 0) {
        if (*c3 == 'G') {
            if (strncmp(c4, "QR", 2) == 0 || strncmp(c4, "RQ", 2) == 0 || strncmp(c4, "LQ", 2) == 0 || strncmp(c4, "QL", 2) == 0 || strncmp(c4, "HR", 2) == 0 || strncmp(c4, "TR", 2) == 0 || strncmp(c4, "BR", 2) == 0) {
                nx = 128;
            }
        }
    } else if (cname && strncmp(c2, "UN", 2) == 0) {
        if (*c3 == 'G') {
            if (strncmp(c4, "QR", 2) == 0 || strncmp(c4, "RQ", 2) == 0 || strncmp(c4, "LQ", 2) == 0 || strncmp(c4, "QL", 2) == 0 || strncmp(c4, "HR", 2) == 0 || strncmp(c4, "TR", 2) == 0 || strncmp(c4, "BR", 2) == 0) {
                nx = 128;
            }
        }
    } else if (strncmp(c2, "GG", 2) == 0) {
        nx = 128;
        if (strncmp(c3, "HD3", 3) == 0) {
            nx = 128;
        }
    }
    return_value = nx;
    return return_value;
//
L80:
    //
    //     ISPEC = 4:  number of shifts (used by xHSEQR)
    //
    return_value = 6;
    return return_value;
//
L90:
    return_value = 2;
    return return_value;

L100:
    //
    //     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
    //
    return_value = castINTEGER(castREAL(min(n1, n2)) * 1.6);
    return return_value;
    //
L110:
    //
    //     ISPEC = 7:  number of processors (not used)
    //
    return_value = 1;
    return return_value;
//
L120:
    //
    //     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
    //
    return_value = 50;
    return return_value;
//
L130:
    //
    //     ISPEC = 9:  maximum size of the subproblems at the bottom of the
    //                 computation tree in the divide-and-conquer algorithm
    //                 (used by xGELSD and xGESDD)
    //
    return_value = 25;
    return return_value;
//
L140:
    //
    //     ISPEC = 10: ieee and infinity NaN arithmetic can be trusted not to trap
    //
    //     iMlaenv = 0
    return_value = 1;
    if (return_value == 1) {
        return_value = iMieeeck(1, 0.0, 1.0);
    }
    return return_value;
    //
L150:
    //
    //     ISPEC = 11: ieee infinity arithmetic can be trusted not to trap
    //
    //     iMlaenv = 0
    return_value = 1;
    if (return_value == 1) {
        return_value = iMieeeck(0, 0.0, 1.0);
    }
    return return_value;
//
L160:
    //
    //     12 <= ISPEC <= 16: xHSEQR or related subroutines.
    //
    return_value = iMparmq(ispec, name, opts, n1, n2, n3, n4);
    return return_value;
    //
    //     End of iMlaenv
    //
}
