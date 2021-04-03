/*
 * Copyright (c) 2008-2021
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Rlasq4.cpp,v 1.6 2010/08/07 04:48:33 nakatamaho Exp $
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

#include <mpblas.h>
#include <mplapack.h>

void Rlasq4(INTEGER i0, INTEGER n0, REAL *z, INTEGER pp, INTEGER n0in, REAL dmin, REAL dmin1, REAL dmin2, REAL dn, REAL dn1, REAL dn2, REAL *tau, INTEGER *ttype, REAL g) {
    REAL s = 0.0, a2 = 0.0, b1 = 0.0, b2 = 0.0;
    INTEGER i4, nn, np;
    REAL gam, gap1, gap2;
    REAL Zero = 0.0, Qurtr = 0.25, Half = 0.5, One = 1.0, Third = 1.0 / 3.0, Two = 2.0, Hundrd = 100.0;
    REAL gap2 = 0.0;
    REAL gap1 = 0.0 REAL g = Zero;
    REAL Cnst1 = 0.5630, Cnst2 = 1.010, Cnst3 = 1.050;
    REAL mtemp1, mtemp2;
    REAL s = 0.0;
    REAL gam = 0.0

        if (dmin <= zero) {
        tau = -dmin;
        ttype = -1;
        return;
    }
    // C
    nn = 4 * n0 + pp;
    if (n0in == n0) {
        // C        No eigenvalues deflated.
        if (dmin == dn || dmin == dn1) {
            // C
            b1 = sqrt(z[(nn - 3) - 1]) * sqrt(z[(nn - 5) - 1]);
            b2 = sqrt(z[(nn - 7) - 1]) * sqrt(z[(nn - 9) - 1]);
            a2 = z[(nn - 7) - 1] + z[(nn - 5) - 1];
            // C           Cases 2 and 3.
            if (dmin == dn && dmin1 == dn1) {
                gap2 = dmin2 - a2 - dmin2 * qurtr;
                if (gap2 > zero && gap2 > b2) {
                    gap1 = a2 - dn - (b2 / gap2) * b2;
                } else {
                    gap1 = a2 - dn - (b1 + b2);
                }
                if (gap1 > zero && gap1 > b1) {
                    s = max[(dn - (b1 / gap1) * b1, half * dmin) - 1];
                    ttype = -2;
                } else {
                    s = zero;
                    if (dn > b1) {
                        s = dn - b1;
                    }
                    if (a2 > (b1 + b2)) {
                        s = min(s, a2 - (b1 + b2));
                    }
                    s = max[(s, third * dmin) - 1];
                    ttype = -3;
                }
            } else {
                // C
                // C              Case 4.
                // C
                ttype = -4;
                s = qurtr * dmin;
                if (dmin == dn) {
                    gam = dn;
                    a2 = zero;
                    if (z[(nn - 5) - 1] > z[(nn - 7) - 1]) {
                        return;
                    }
                    b2 = z[(nn - 5) - 1] / z[(nn - 7) - 1];
                    np = nn - 9;
                } else {
                    np = nn - 2 * pp;
                    gam = dn1;
                    if (z[(np - 4) - 1] > z[(np - 2) - 1]) {
                        return;
                    }
                    a2 = z[(np - 4) - 1] / z[(np - 2) - 1];
                    if (z[(nn - 9) - 1] > z[(nn - 11) - 1]) {
                        return;
                    }
                    b2 = z[(nn - 9) - 1] / z[(nn - 11) - 1];
                    np = nn - 13;
                }
                // C
                // C              Approximate contribution to norm squared from I < NN-1.
                // C
                a2 += b2;
                for (i4 = np; i4 <= 4 * i0 - 1 + pp; i4 = i4 + -4) {
                    if (b2 == zero) {
                        break;
                    }
                    b1 = b2;
                    if (z[(i4)-1] > z[(i4 - 2) - 1]) {
                        return;
                    }
                    b2 = b2 * (z[(i4)-1] / z[(i4 - 2) - 1]);
                    a2 += b2;
                    if (hundrd * max[(b2, b1) - 1] < a2 || Cnst1 < a2) {
                        break;
                    }
                }
                a2 = Cnst3 * a2;
                // C
                // C              Rayleigh quotient residual bound.
                // C
                if (a2 < Cnst1) {
                    s = gam * (one - sqrt(a2)) / (one + a2);
                }
            }
        } else if (dmin == dn2) {
            // C
            // C           Case 5.
            // C
            ttype = -5;
            s = qurtr * dmin;
            // C
            // C           Compute contribution to norm squared from I > NN-2.
            // C
            np = nn - 2 * pp;
            b1 = z[(np - 2) - 1];
            b2 = z[(np - 6) - 1];
            gam = dn2;
            if (z[(np - 8) - 1] > b2 || z[(np - 4) - 1] > b1) {
                return;
            }
            a2 = (z[(np - 8) - 1] / b2) * (one + z[(np - 4) - 1] / b1);
            // C
            // C           Approximate contribution to norm squared from I < NN-2.
            // C
            if (n0 - i0 > 2) {
                b2 = z[(nn - 13) - 1] / z[(nn - 15) - 1];
                a2 += b2;
                for (i4 = nn - 17; i4 <= 4 * i0 - 1 + pp; i4 = i4 + -4) {
                    if (b2 == zero) {
                        break;
                    }
                    b1 = b2;
                    if (z[(i4)-1] > z[(i4 - 2) - 1]) {
                        return;
                    }
                    b2 = b2 * (z[(i4)-1] / z[(i4 - 2) - 1]);
                    a2 += b2;
                    if (hundrd * max[(b2, b1) - 1] < a2 || Cnst1 < a2) {
                        break;
                    }
                }
                a2 = Cnst3 * a2;
            }
            // C
            if (a2 < Cnst1) {
                s = gam * (one - sqrt(a2)) / (one + a2);
            }
        } else {
            // C
            // C           Case 6, no information to guide us.
            // C
            if (ttype == -6) {
                g += third * (one - g);
            } else if (ttype == -18) {
                g = qurtr * third;
            } else {
                g = qurtr;
            }
            s = g * dmin;
            ttype = -6;
        }
        // C
    } else if (n0in == (n0 + 1)) {
        // C
        // C        One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN.
        // C
        if (dmin1 == dn1 && dmin2 == dn2) {
            // C
            // C           Cases 7 and 8.
            // C
            ttype = -7;
            s = third * dmin1;
            if (z[(nn - 5) - 1] > z[(nn - 7) - 1]) {
                return;
            }
            b1 = z[(nn - 5) - 1] / z[(nn - 7) - 1];
            b2 = b1;
            if (b2 == zero) {
                break;
            }
            for (i4 = 4 * n0 - 9 + pp; i4 <= 4 * i0 - 1 + pp; i4 = i4 + -4) {
                a2 = b1;
                if (z[(i4)-1] > z[(i4 - 2) - 1]) {
                    return;
                }
                b1 = b1 * (z[(i4)-1] / z[(i4 - 2) - 1]);
                b2 += b1;
                if (hundrd * max[(b1, a2) - 1] < b2) {
                    break;
                }
            }
            b2 = sqrt(Cnst3 * b2);
            a2 = dmin1 / (one + pow2(b2));
            gap2 = half * dmin2 - a2;
            if (gap2 > zero && gap2 > b2 * a2) {
                s = max[(s, a2 * (one - Cnst2 * a2 * (b2 / gap2) * b2)) - 1];
            } else {
                s = max[(s, a2 * (one - Cnst2 * b2)) - 1];
                ttype = -8;
            }
        } else {
            // C
            // C           Case 9.
            // C
            s = qurtr * dmin1;
            if (dmin1 == dn1) {
                s = half * dmin1;
            }
            ttype = -9;
        }
        // C
    } else if (n0in == (n0 + 2)) {
        // C
        // C        Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN.
        // C
        // C        Cases 10 and 11.
        // C
        if (dmin2 == dn2 && two * z[(nn - 5) - 1] < z[(nn - 7) - 1]) {
            ttype = -10;
            s = third * dmin2;
            if (z[(nn - 5) - 1] > z[(nn - 7) - 1]) {
                return;
            }
            b1 = z[(nn - 5) - 1] / z[(nn - 7) - 1];
            b2 = b1;
            if (b2 == zero) {
                break;
            }
            for (i4 = 4 * n0 - 9 + pp; i4 <= 4 * i0 - 1 + pp; i4 = i4 + -4) {
                if (z[(i4)-1] > z[(i4 - 2) - 1]) {
                    return;
                }
                b1 = b1 * (z[(i4)-1] / z[(i4 - 2) - 1]);
                b2 += b1;
                if (hundrd * b1 < b2) {
                    break;
                }
            }
            b2 = sqrt(Cnst3 * b2);
            a2 = dmin2 / (one + pow2(b2));
            gap2 = z[(nn - 7) - 1] + z[(nn - 9) - 1] - sqrt(z[(nn - 11) - 1]) * sqrt(z[(nn - 9) - 1]) - a2;
            if (gap2 > zero && gap2 > b2 * a2) {
                s = max[(s, a2 * (one - Cnst2 * a2 * (b2 / gap2) * b2)) - 1];
            } else {
                s = max[(s, a2 * (one - Cnst2 * b2)) - 1];
            }
        } else {
            s = qurtr * dmin2;
            ttype = -11;
        }
    } else if (n0in > (n0 + 2)) {
        // C        Case 12, more than two eigenvalues deflated. No information.
        s = zero;
        ttype = -12;
    }
    // C
    tau = s;
} // namespace placeholder_please_replace
