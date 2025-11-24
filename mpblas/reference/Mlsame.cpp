/*
 * Copyright (c) 2008-2025
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

LOGICAL
Mlsame(const char *ca, const char *cb) {
    LOGICAL return_value = false;
    //
    // Test if the characters are equal
    //
    return_value = ca == cb;
    if (return_value) {
        return return_value;
    }
    //
    // Now test for equivalence if both characters are alphabetic.
    //
    INTEGER zcode = fem::ichar("Z");
    //
    // Use 'Z' rather than 'A' so that ASCII can be detected on Prime
    // machines, on which ICHAR returns a value with bit 8 set.
    // ICHAR('A') on Prime machines returns 193 which is the same as
    // ICHAR('A') on an EBCDIC machine.
    //
    INTEGER inta = fem::ichar(ca);
    INTEGER intb = fem::ichar(cb);
    //
    if (zcode == 90 || zcode == 122) {
        //
        // ASCII is assumed - ZCODE is the ASCII code of either lower or
        // upper case 'Z'.
        //
        if (inta >= 97 && inta <= 122) {
            inta = inta - 32;
        }
        if (intb >= 97 && intb <= 122) {
            intb = intb - 32;
        }
        //
    } else if (zcode == 233 || zcode == 169) {
        //
        // EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
        // upper case 'Z'.
        //
        if (inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta >= 162 && inta <= 169) {
            inta += 64;
        }
        if (intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb >= 162 && intb <= 169) {
            intb += 64;
        }
        //
    } else if (zcode == 218 || zcode == 250) {
        //
        // ASCII is assumed, on Prime machines - ZCODE is the ASCII code
        // plus 128 of either lower or upper case 'Z'.
        //
        if (inta >= 225 && inta <= 250) {
            inta = inta - 32;
        }
        if (intb >= 225 && intb <= 250) {
            intb = intb - 32;
        }
    }
    return_value = inta == intb;
    return return_value;
    //
    // RETURN
    //
    // End of Mlsame
    //
}
