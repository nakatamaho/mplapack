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

INTEGER
iMparmq(INTEGER const ispec, const char *name, const char * /* opts */, INTEGER const /* n */, INTEGER const ilo, INTEGER const ihi, INTEGER const /* lwork */) {
    INTEGER return_value = 0;
    const INTEGER ishfts = 15;
    const INTEGER inwin = 13;
    const INTEGER iacc22 = 16;
    INTEGER nh = 0;
    INTEGER ns = 0;
    const REAL two = 2.0;
    if ((ispec == ishfts) || (ispec == inwin) || (ispec == iacc22)) {
        //
        // ==== Set the number simultaneous shifts ====
        //
        nh = ihi - ilo + 1;
        ns = 2;
        if (nh >= 30) {
            ns = 4;
        }
        if (nh >= 60) {
            ns = 10;
        }
        if (nh >= 150) {
            ns = max((INTEGER)10, nh / nint(log(castREAL(nh)) / log(two)));
        }
        if (nh >= 590) {
            ns = 64;
        }
        if (nh >= 3000) {
            ns = 128;
        }
        if (nh >= 6000) {
            ns = 256;
        }
        ns = max((INTEGER)2, ns - mod(ns, 2));
    }
    //
    const INTEGER inmin = 12;
    const INTEGER nmin = 75;
    const INTEGER inibl = 14;
    const INTEGER nibble = 14;
    const INTEGER knwswp = 500;
    char subnam[6];
    INTEGER ic = 0;
    INTEGER iz = 0;
    INTEGER i = 0;
    const INTEGER k22min = 14;
    const INTEGER kacmin = 14;
    if (ispec == inmin) {
        //
        // ===== Matrices of order smaller than NMIN get sent
        // .     to xLAHQR, the classic double shift algorithm.
        // .     This must be at least 11. ====
        //
        return_value = nmin;
        //
    } else if (ispec == inibl) {
        //
        // ==== INIBL: skip a multi-shift qr iteration and
        // .    whenever aggressive early deflation finds
        // .    at least (NIBBLE*(window size)/100) deflations. ====
        //
        return_value = nibble;
        //
    } else if (ispec == ishfts) {
        //
        // ==== NSHFTS: The number of simultaneous shifts =====
        //
        return_value = ns;
        //
    } else if (ispec == inwin) {
        //
        // ==== NW: deflation window size.  ====
        //
        if (nh <= knwswp) {
            return_value = ns;
        } else {
            return_value = 3 * ns / 2;
        }
        //
    } else if (ispec == iacc22) {
        //
        // ==== IACC22: Whether to accumulate reflections
        // .     before updating the far-from-diagonal elements
        // .     and whether to use 2-by-2 block structure while
        // .     doing it.  A small amount of work could be saved
        // .     by making this choice dependent also upon the
        // .     NH=IHI-ILO+1.
        //
        // Convert NAME to upper case if the first character is lower case.
        //
        return_value = 0;
        subnam = name;
        ic = fem::ichar(subnam[0]);
        iz = fem::ichar("Z");
        if (iz == 90 || iz == 122) {
            //
            // ASCII character set
            //
            if (ic >= 97 && ic <= 122) {
                subnam[0] = fem::fchar(ic - 32);
                for (i = 2; i <= 6; i = i + 1) {
                    ic = fem::ichar(subnam(i, i));
                    if (ic >= 97 && ic <= 122) {
                        subnam(i, i) = fem::fchar(ic - 32);
                    }
                }
            }
            //
        } else if (iz == 233 || iz == 169) {
            //
            // EBCDIC character set
            //
            if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 && ic <= 169)) {
                subnam[0] = fem::fchar(ic + 64);
                for (i = 2; i <= 6; i = i + 1) {
                    ic = fem::ichar(subnam(i, i));
                    if ((ic >= 129 && ic <= 137) || (ic >= 145 && ic <= 153) || (ic >= 162 && ic <= 169)) {
                        subnam(i, i) = fem::fchar(ic + 64);
                    }
                }
            }
            //
        } else if (iz == 218 || iz == 250) {
            //
            // Prime machines:  ASCII+128
            //
            if (ic >= 225 && ic <= 250) {
                subnam[0] = fem::fchar(ic - 32);
                for (i = 2; i <= 6; i = i + 1) {
                    ic = fem::ichar(subnam(i, i));
                    if (ic >= 225 && ic <= 250) {
                        subnam(i, i) = fem::fchar(ic - 32);
                    }
                }
            }
        }
        //
        if (subnam(2, 6) == "GGHRD" || subnam(2, 6) == "GGHD3") {
            return_value = 1;
            if (nh >= k22min) {
                return_value = 2;
            }
        } else if (subnam(4, 6) == "EXC") {
            if (nh >= kacmin) {
                return_value = 1;
            }
            if (nh >= k22min) {
                return_value = 2;
            }
        } else if (subnam(2, 6) == "HSEQR" || subnam(2, 5) == "LAQR") {
            if (ns >= kacmin) {
                return_value = 1;
            }
            if (ns >= k22min) {
                return_value = 2;
            }
        }
        //
    } else {
        // ===== invalid value of ispec =====
        return_value = -1;
        //
    }
    return return_value;
    //
    // ==== End of iMparmq ====
    //
}
