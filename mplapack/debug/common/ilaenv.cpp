/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 * $Id: ilaenv.cpp,v 1.4 2010/08/19 01:17:55 nakatamaho Exp $
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
 *
 * $Id: ilaenv.cpp,v 1.4 2010/08/19 01:17:55 nakatamaho Exp $

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

#include <string.h>
#include <ctype.h>

#define MLANAMESIZE 6

//ISPEC = 1:  block size
//In these examples, separate code is provided for setting NB for
//real and complex.  We assume that NB will take the same value in
//single or double precision.

int
ilaenv1(const char *Mlaname, const char *opts, int n1, int n2, int n3, int n4)
{

    int nb = 1;
#if !defined (IMLAENV_DEBUG)
    if (strcmp(&Mlaname[1],"orgqr") == 0) { nb = 32; return nb; }
    if (strcmp(&Mlaname[1],"orgql") == 0) { nb = 32; return nb; }
    if (strcmp(&Mlaname[1],"potrf") == 0) { nb = 64; return nb; }
    if (strcmp(&Mlaname[1],"trtri") == 0) { nb = 64; return nb; }
    if (strcmp(&Mlaname[1],"dsytrd") == 0) { nb = 32;return nb;  }
    if (strcmp(&Mlaname[1],"getrf") == 0)  { nb = 64;return nb;  }
    if (strcmp(&Mlaname[1],"getri") == 0)  { nb = 64;return nb;  }
    if (strcmp(&Mlaname[1],"lauum") == 0)  { nb = 64;return nb;  }
#else
    if (strcmp(&Mlaname[1],"potrf") == 0)  { nb = 8;return nb; }
    if (strcmp(&Mlaname[1],"orgqr") == 0)  { nb = 8;return nb; }
    if (strcmp(&Mlaname[1],"orgql") == 0)  { nb = 8;return nb; }
    if (strcmp(&Mlaname[1],"trtri") == 0)  { nb = 8;return nb; }
    if (strcmp(&Mlaname[0],"dsytrd") == 0) { nb = 8;return nb; }
    if (strcmp(&Mlaname[1],"getrf") == 0)  { nb = 8;return nb; }
    if (strcmp(&Mlaname[1],"getri") == 0)  { nb = 8;return nb; }
    if (strcmp(&Mlaname[1],"lauum") == 0)  { nb = 8;return nb;  }
#endif
    return nb;
}

//*     ISPEC = 2:  minimum block size
int
ilaenv2(const char *Mlaname, const char *opts, int n1, int n2, int n3, int n4)
{
    int nbmin = 1;
    if (strcmp(&Mlaname[1], "orgqr") == 0)  { nbmin = 2; return nbmin; }
    if (strcmp(&Mlaname[1], "orgql") == 0)  { nbmin = 2; return nbmin; }
    if (strcmp(&Mlaname[1], "trtri") == 0)  { nbmin = 2; return nbmin; }
    if (strcmp(&Mlaname[0], "dsytrd") == 0) { nbmin = 2; return nbmin; }
    if (strcmp(&Mlaname[0], "getri") == 0)  { nbmin = 2; return nbmin; }

    return nbmin;
}

//*     ISPEC = 3:  crossover point
int
ilaenv3(const char *Mlaname, const char *opts, int n1, int n2, int n3, int n4)
{
    int nx = 1;
#if !defined (IMLAENV_DEBUG)
    if (strcmp(&Mlaname[1],"orgqr")==0) { nx = 128; return nx; }
    if (strcmp(&Mlaname[1],"orgql")==0) { nx = 128; return nx; }
    if (strcmp(&Mlaname[0],"dsytrd")==0){ nx = 32; return nx; }
#else
    if (strcmp(&Mlaname[1],"orgqr") == 0) { nx = 6; return nx; }
    if (strcmp(&Mlaname[1],"orgql") == 0) { nx = 6; return nx; }
    if (strcmp(&Mlaname[0], "dsytrd")== 0){ nx = 6; return nx; }
#endif
    return nx;
}

int
ilaenv4(const char *Mlaname, const char *opts, int n1, int n2, int n3, int n4)
{
    return 1;
}

int
ilaenv5(const char *Mlaname, const char *opts, int n1, int n2, int n3, int n4)
{
    return 1;
}

int
ilaenv6(const char *Mlaname, const char *opts, int n1, int n2, int n3, int n4)
{
    return 1;
}

int
ilaenv7(const char *Mlaname, const char *opts, int n1, int n2, int n3, int n4)
{
    return 1;
}

int
ilaenv8(const char *Mlaname, const char *opts, int n1, int n2, int n3, int n4)
{
    return 1;
}

int
ilaenv9(const char *Mlaname, const char *opts, int n1, int n2, int n3, int n4)
{
    return 1;
}

int
ilaenv10(const char *Mlaname, const char *opts, int n1, int n2, int n3, int n4)
{
    return 1;
}

int
ilaenv11(const char *Mlaname, const char *opts, int n1, int n2, int n3, int n4)
{
    return 1;
}

int
ilaenv12(const char *Mlaname, const char *opts, int n1, int n2, int n3, int n4)
{
    return 1;
}

int
ilaenv13(const char *Mlaname, const char *opts, int n1, int n2, int n3, int n4)
{
    return 1;
}

int
ilaenv14(const char *Mlaname, const char *opts, int n1, int n2, int n3, int n4)
{
    return 1;
}

int
ilaenv15(const char *Mlaname, const char *opts, int n1, int n2, int n3, int n4)
{
    return 1;
}

int
ilaenv16(const char *Mlaname, const char *opts, int n1, int n2, int n3, int n4)
{
    return 1;
}

int
ilaenv_f77(int *ispec, const char *name, const char *opts, int *n1, int *n2,
    int *n3, int *n4)
{
    int iret, i, up;

    iret = -1;
    char Mlaname[MLANAMESIZE + 1];

    strncpy(Mlaname, name, MLANAMESIZE);
    for (i = 0; i < MLANAMESIZE; i++) {
	up = tolower(Mlaname[i]);
	Mlaname[i] = up;
    }
/*
    if (!Mlsame(Mlaname, "r") && !Mlsame(Mlaname, "c"))
	return iret;
*/
    switch (*ispec) {
    case 1:
	iret = ilaenv1(Mlaname, opts, *n1, *n2, *n3, *n4);
	break;
    case 2:
	iret = ilaenv2(Mlaname, opts, *n1, *n2, *n3, *n4);
	break;
    case 3:
	iret = ilaenv3(Mlaname, opts, *n1, *n2, *n3, *n4);
	break;
    case 4:
	iret = ilaenv4(Mlaname, opts, *n1, *n2, *n3, *n4);
	break;
    case 5:
	iret = ilaenv5(Mlaname, opts, *n1, *n2, *n3, *n4);
	break;
    case 6:
	iret = ilaenv6(Mlaname, opts, *n1, *n2, *n3, *n4);
	break;
    case 7:
	iret = ilaenv7(Mlaname, opts, *n1, *n2, *n3, *n4);
	break;
    case 8:
	iret = ilaenv8(Mlaname, opts, *n1, *n2, *n3, *n4);
	break;
    case 9:
	iret = ilaenv9(Mlaname, opts, *n1, *n2, *n3, *n4);
	break;
    case 10:
	iret = ilaenv10(Mlaname, opts, *n1, *n2, *n3, *n4);
	break;
    case 11:
	iret = ilaenv11(Mlaname, opts, *n1, *n2, *n3, *n4);
	break;
    case 12:
	iret = ilaenv12(Mlaname, opts, *n1, *n2, *n3, *n4);
	break;
    case 13:
	iret = ilaenv13(Mlaname, opts, *n1, *n2, *n3, *n4);
	break;
    case 14:
	iret = ilaenv14(Mlaname, opts, *n1, *n2, *n3, *n4);
	break;
    case 15:
	iret = ilaenv15(Mlaname, opts, *n1, *n2, *n3, *n4);
	break;
    case 16:
	iret = ilaenv16(Mlaname, opts, *n1, *n2, *n3, *n4);
	break;
    default:
	iret = -1;
    }
    return iret;
}
