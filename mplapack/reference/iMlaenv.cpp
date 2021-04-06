/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: iMlaenv.cpp,v 1.12 2010/08/19 01:20:10 nakatamaho Exp $
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
/*
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

$COPYRIGHT$

Additional copyrights may follow

$HEADER$

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.

- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <ctype.h>
#include <mpblas.h>
#include <mplapack.h>
#include <stdio.h>
#include <string.h>

#define MLANAMESIZE 6

// ISPEC = 1:  block size
// In these examples, separate code is provided for setting NB for
// real and complex.  We assume that NB will take the same value in
// single or double precision.

INTEGER iMlaenv1(const char *Mlaname, const char *opts, INTEGER n1, INTEGER n2, INTEGER n3, INTEGER n4) {
    INTEGER nb = 1;
#if !defined(IMLAENV_DEBUG)
    if (strcmp(&Mlaname[1], "orgqr") == 0) {
        nb = 32;
        return nb;
    }
    if (strcmp(&Mlaname[1], "ungqr") == 0) {
        nb = 32;
        return nb;
    }
    if (strcmp(&Mlaname[1], "orgql") == 0) {
        nb = 32;
        return nb;
    }
    if (strcmp(&Mlaname[1], "ungql") == 0) {
        nb = 32;
        return nb;
    }
    if (strcmp(&Mlaname[1], "potrf") == 0) {
        nb = 64;
        return nb;
    }
    if (strcmp(&Mlaname[1], "trtri") == 0) {
        nb = 64;
        return nb;
    }
    if (strcmp(&Mlaname[0], "rsytrd") == 0) {
        nb = 32;
        return nb;
    }
    if (strcmp(&Mlaname[0], "chetrd") == 0) {
        nb = 32;
        return nb;
    }
    if (strcmp(&Mlaname[1], "getrf") == 0) {
        nb = 64;
        return nb;
    }
    if (strcmp(&Mlaname[1], "getri") == 0) {
        nb = 64;
        return nb;
    }
    if (strcmp(&Mlaname[1], "lauum") == 0) {
        nb = 64;
        return nb;
    }
#else
    if (strcmp(&Mlaname[1], "orgqr") == 0) {
        nb = 8;
        return nb;
    }
    if (strcmp(&Mlaname[1], "ungqr") == 0) {
        nb = 8;
        return nb;
    }
    if (strcmp(&Mlaname[1], "orgql") == 0) {
        nb = 8;
        return nb;
    }
    if (strcmp(&Mlaname[1], "ungql") == 0) {
        nb = 8;
        return nb;
    }
    if (strcmp(&Mlaname[1], "potrf") == 0) {
        nb = 8;
        return nb;
    }
    if (strcmp(&Mlaname[1], "trtri") == 0) {
        nb = 8;
        return nb;
    }
    if (strcmp(&Mlaname[0], "rsytrd") == 0) {
        nb = 8;
        return nb;
    }
    if (strcmp(&Mlaname[0], "chetrd") == 0) {
        nb = 8;
        return nb;
    }
    if (strcmp(&Mlaname[1], "getrf") == 0) {
        nb = 8;
        return nb;
    }
    if (strcmp(&Mlaname[1], "getri") == 0) {
        nb = 8;
        return nb;
    }
    if (strcmp(&Mlaname[1], "lauum") == 0) {
        nb = 8;
        return nb;
    }
#endif
    return nb;
}

// ISPEC = 2:  minimum block size
INTEGER iMlaenv2(const char *Mlaname, const char *opts, INTEGER n1, INTEGER n2, INTEGER n3, INTEGER n4) {
    INTEGER nbmin = 1;
    if (strcmp(&Mlaname[1], "orgqr") == 0) {
        nbmin = 2;
        return nbmin;
    }
    if (strcmp(&Mlaname[1], "ungqr") == 0) {
        nbmin = 2;
        return nbmin;
    }
    if (strcmp(&Mlaname[1], "orgql") == 0) {
        nbmin = 2;
        return nbmin;
    }
    if (strcmp(&Mlaname[1], "ungql") == 0) {
        nbmin = 2;
        return nbmin;
    }
    if (strcmp(&Mlaname[1], "trtri") == 0) {
        nbmin = 2;
        return nbmin;
    }
    if (strcmp(&Mlaname[0], "rsytrd") == 0) {
        nbmin = 2;
        return nbmin;
    }
    if (strcmp(&Mlaname[0], "chetrd") == 0) {
        nbmin = 2;
        return nbmin;
    }
    if (strcmp(&Mlaname[1], "getri") == 0) {
        nbmin = 2;
        return nbmin;
    }
    return nbmin;
}

// ISPEC = 3:  crossover point
INTEGER iMlaenv3(const char *Mlaname, const char *opts, INTEGER n1, INTEGER n2, INTEGER n3, INTEGER n4) {
    INTEGER nx = 1;
#if !defined(IMLAENV_DEBUG)
    if (strcmp(&Mlaname[1], "orgqr") == 0) {
        nx = 128;
        return nx;
    }
    if (strcmp(&Mlaname[1], "ungqr") == 0) {
        nx = 128;
        return nx;
    }
    if (strcmp(&Mlaname[1], "orgql") == 0) {
        nx = 128;
        return nx;
    }
    if (strcmp(&Mlaname[1], "ungql") == 0) {
        nx = 128;
        return nx;
    }
    if (strcmp(&Mlaname[0], "rsytrd") == 0) {
        nx = 32;
        return nx;
    }
    if (strcmp(&Mlaname[0], "chetrd") == 0) {
        nx = 32;
        return nx;
    }
#else
    if (strcmp(&Mlaname[1], "orgqr") == 0) {
        nx = 6;
        return nx;
    }
    if (strcmp(&Mlaname[1], "ungqr") == 0) {
        nx = 6;
        return nx;
    }
    if (strcmp(&Mlaname[1], "orgql") == 0) {
        nx = 6;
        return nx;
    }
    if (strcmp(&Mlaname[1], "ungql") == 0) {
        nx = 6;
        return nx;
    }
    if (strcmp(&Mlaname[0], "rsytrd") == 0) {
        nx = 6;
        return nx;
    }
    if (strcmp(&Mlaname[0], "chetrd") == 0) {
        nx = 6;
        return nx;
    }
#endif
    return nx;
}

INTEGER iMlaenv4(const char *Mlaname, const char *opts, INTEGER n1, INTEGER n2, INTEGER n3, INTEGER n4) { return 1; }

INTEGER iMlaenv5(const char *Mlaname, const char *opts, INTEGER n1, INTEGER n2, INTEGER n3, INTEGER n4) { return 1; }

INTEGER iMlaenv6(const char *Mlaname, const char *opts, INTEGER n1, INTEGER n2, INTEGER n3, INTEGER n4) { return 1; }

INTEGER iMlaenv7(const char *Mlaname, const char *opts, INTEGER n1, INTEGER n2, INTEGER n3, INTEGER n4) { return 1; }

INTEGER iMlaenv8(const char *Mlaname, const char *opts, INTEGER n1, INTEGER n2, INTEGER n3, INTEGER n4) { return 1; }

INTEGER iMlaenv9(const char *Mlaname, const char *opts, INTEGER n1, INTEGER n2, INTEGER n3, INTEGER n4) { return 1; }

INTEGER iMlaenv10(const char *Mlaname, const char *opts, INTEGER n1, INTEGER n2, INTEGER n3, INTEGER n4) { return 1; }

INTEGER iMlaenv11(const char *Mlaname, const char *opts, INTEGER n1, INTEGER n2, INTEGER n3, INTEGER n4) { return 1; }

INTEGER iMlaenv12(const char *Mlaname, const char *opts, INTEGER n1, INTEGER n2, INTEGER n3, INTEGER n4) { return 1; }

INTEGER iMlaenv13(const char *Mlaname, const char *opts, INTEGER n1, INTEGER n2, INTEGER n3, INTEGER n4) { return 1; }

INTEGER iMlaenv14(const char *Mlaname, const char *opts, INTEGER n1, INTEGER n2, INTEGER n3, INTEGER n4) { return 1; }

INTEGER iMlaenv15(const char *Mlaname, const char *opts, INTEGER n1, INTEGER n2, INTEGER n3, INTEGER n4) { return 1; }

INTEGER iMlaenv16(const char *Mlaname, const char *opts, INTEGER n1, INTEGER n2, INTEGER n3, INTEGER n4) { return 1; }

INTEGER iMlaenv(INTEGER ispec, const char *name, const char *opts, INTEGER n1, INTEGER n2, INTEGER n3, INTEGER n4) {
    INTEGER iret, i, up, len;
    iret = -1;
    char Mlaname[MLANAMESIZE + 1] = "000000";
    // buggy
    len = strlen(name);
    strncpy(Mlaname, name, (len > MLANAMESIZE) ? MLANAMESIZE : len);
    for (i = 0; i < MLANAMESIZE; i++) {
        up = tolower(Mlaname[i]);
        Mlaname[i] = up;
    }
    Mlaname[MLANAMESIZE] = '\0';

    if (!Mlsame(Mlaname, "r") && !Mlsame(Mlaname, "c"))
        return iret;
    switch (ispec) {
    case 1:
        iret = iMlaenv1(Mlaname, opts, n1, n2, n3, n4);
        break;
    case 2:
        iret = iMlaenv2(Mlaname, opts, n1, n2, n3, n4);
        break;
    case 3:
        iret = iMlaenv3(Mlaname, opts, n1, n2, n3, n4);
        break;
    case 4:
        iret = iMlaenv4(Mlaname, opts, n1, n2, n3, n4);
        break;
    case 5:
        iret = iMlaenv5(Mlaname, opts, n1, n2, n3, n4);
        break;
    case 6:
        iret = iMlaenv6(Mlaname, opts, n1, n2, n3, n4);
        break;
    case 7:
        iret = iMlaenv7(Mlaname, opts, n1, n2, n3, n4);
        break;
    case 8:
        iret = iMlaenv8(Mlaname, opts, n1, n2, n3, n4);
        break;
    case 9:
        iret = iMlaenv9(Mlaname, opts, n1, n2, n3, n4);
        break;
    case 10:
        iret = iMlaenv10(Mlaname, opts, n1, n2, n3, n4);
        break;
    case 11:
        iret = iMlaenv11(Mlaname, opts, n1, n2, n3, n4);
        break;
    case 12:
        iret = iMlaenv12(Mlaname, opts, n1, n2, n3, n4);
        break;
    case 13:
        iret = iMlaenv13(Mlaname, opts, n1, n2, n3, n4);
        break;
    case 14:
        iret = iMlaenv14(Mlaname, opts, n1, n2, n3, n4);
        break;
    case 15:
        iret = iMlaenv15(Mlaname, opts, n1, n2, n3, n4);
        break;
    case 16:
        iret = iMlaenv16(Mlaname, opts, n1, n2, n3, n4);
        break;
    default:
        iret = -1;
    }
    return iret;
}
