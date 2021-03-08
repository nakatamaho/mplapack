//https://www.johndcook.com/blog/2019/11/12/rump-floating-point/
/*
 * Copyright (c) 2010-2021
 *	Nakata, Maho
 * 	All rights reserved.
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

#include "mpreal.h"
#include <iostream>

using namespace mpfr;
using namespace std;

void rump_example(int prec) {
    mpreal::set_default_prec(prec);
    mpreal x = "77617";
    mpreal y = "33096";
    mpreal z;
    mpreal x2 = x * x;
    mpreal y2 = y * y;
    mpreal y4 = y2 * y2;
    mpreal y6 = y2 * y4;
    mpreal y8 = y4 * y4;
    z = 333.75 * y6 + x2 * (11 * x2 * y2 - y6 - 121 * y4 - 2.0) + 5.5 * y8 + x / (2.0 * y);
    mpfr_printf("Precision %4ld, result: %2.8RNe\n", prec, (mpfr_ptr)z);
}

int main(int argc, char *argv[]) {
    cout << "Rump's example\n";
    rump_example(24);
    rump_example(53);
    rump_example(116);
    rump_example(128);
    rump_example(256);
    rump_example(1024);
    return 0;
}
