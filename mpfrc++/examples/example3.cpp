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

// https://www.johndcook.com/blog/2019/11/01/arbitrary-precision/
// https://oeis.org/A249836
// tan(n) > n
// 1, 260515, 37362253, 122925461, 534483448, 3083975227, 902209779836, 74357078147863, 214112296674652, 642336890023956, 18190586279576483, 248319196091979065, 1108341089274117551, 118554299812338354516058, 1428599129020608582548671, 4285797387061825747646013

template <typename TYPE, std::size_t SIZE> std::size_t array_length(const TYPE (&)[SIZE]) { return SIZE; }

void tan_n_larger_n(int prec, mpreal n) {
    mpreal::set_default_prec(prec);
    mpreal _n;
    _n = mpreal(n);
    mpfr_printf("prec %d:  tan(%.0RF)= %2.8Re tan(n)-n = %2.8Re ", prec, _n, tan(_n), tan(_n) - (_n));
    if (tan(_n) > _n)
        printf("ok\n");
    else
        printf("ng\n");
}

int main(int argc, char *argv[]) {

    mpreal known_list[] = {"1", "260515", "37362253", "122925461", "534483448", "3083975227", "902209779836", "74357078147863", "214112296674652", "642336890023956", "18190586279576483", "248319196091979065", "1108341089274117551", "118554299812338354516058", "1428599129020608582548671", "4285797387061825747646013"};

    int prec[] = {23, 53, 128, 256, 512, 1024, 100000, 100000000};

    for (int j = 0; j < array_length(prec); j++) {
        for (int i = 0; i < array_length(known_list); i++) {
            tan_n_larger_n(prec[j], known_list[i]);
        }
    }
    return 0;
}
