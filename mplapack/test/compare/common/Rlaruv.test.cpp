
/*
 * Copyright (c) 2008-2026
 *	Nakata, Maho
 * 	All rights reserved.
 *
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
#include <mplapack_compare_debug.h>

#include <fstream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <climits>
#include <limits>

#define MAX_N 128
#define outname "Rlaruv.txt"

static bool parse_int(const char *s, long &out) {
    errno = 0;
    char *end = nullptr;
    long v = std::strtol(s, &end, 10);
    if (errno != 0 || end == s || *end != '\0')
        return false;
    out = v;
    return true;
}

static void usage(const char *prog) {
    std::fprintf(stderr,
                 "Usage:\n"
                 "  %s [iseed1 iseed2 iseed3 iseed4]\n\n"
                 "If omitted, defaults to: 4 3 2 1\n"
                 "Seed rules (MPLAPACK extension):\n"
                 "  - 0 0 0 0 : non-deterministic (seeded from current time / random_device); NOT reproducible\n"
                 "  - otherwise: deterministic (same iseed => same sequence); reproducible\n"
                 "LAPACK convention (deterministic mode): 0<=iseed1..3<=4095 and iseed4 must be odd.\n",
                 prog);
}

int main(int argc, char *argv[]) {

    printf("*** Testing Rlaruv start ***\n");

    // Default seed
    INTEGER seed[4] = {4, 3, 2, 1};

    // Parse optional command-line seeds
    if (argc != 1 && argc != 5) {
        usage(argv[0]);
        return 2;
    }
    if (argc == 5) {
        for (int i = 0; i < 4; ++i) {
            long v = 0;
            if (!parse_int(argv[i + 1], v)) {
                std::fprintf(stderr, "Error: invalid integer for iseed%d: '%s'\n", i + 1, argv[i + 1]);
                return 2;
            }
            if (v < static_cast<long>(std::numeric_limits<INTEGER>::min()) || v > static_cast<long>(std::numeric_limits<INTEGER>::max())) {
                std::fprintf(stderr, "Error: iseed%d out of range for INTEGER: %ld\n", i + 1, v);
                return 2;
            }
            seed[i] = static_cast<INTEGER>(v);
        }
    }

    // Validate seed according to LAPACK convention
    for (int i = 0; i < 3; ++i) {
        if (seed[i] < 0 || seed[i] > 4095) {
            std::fprintf(stderr, "Error: iseed%d must be in [0,4095], got %ld\n", i + 1, static_cast<long>(seed[i]));
            return 2;
        }
    }
    if ((seed[3] % 2) == 0) {
        std::fprintf(stderr, "Error: iseed4 must be odd, got %ld\n", static_cast<long>(seed[3]));
        return 2;
    }

    std::ofstream outputfile(outname);
    if (!outputfile) {
        std::fprintf(stderr, "Failed to open output file: %s\n", outname);
        return 1;
    }

    // Write a small header to the output file (no stdout)
    outputfile << "# Rlaruv test output\n";
    outputfile << "# iseed = " << static_cast<long>(seed[0]) << " " << static_cast<long>(seed[1]) << " " << static_cast<long>(seed[2]) << " " << static_cast<long>(seed[3]) << "\n";
    outputfile << "# Format: one random number per line\n";

    // Optional: generate small vectors (len=1..MAX_N-1) and dump them
    for (int len = 1; len < MAX_N; ++len) {
        std::vector<REAL> x(static_cast<size_t>(len));
        Rlaruv(seed, len, x.data());
        for (int i = 0; i < len; ++i) {
            char buf[4096];
            sprintnum(buf, x[i]);
            outputfile << buf << "\n";
        }
    }

    // Large output block
    const int big_n = 100000;
    std::vector<REAL> x(big_n);
    Rlaruv(seed, big_n, x.data());
    for (int i = 0; i < big_n; ++i) {
        char buf[4096];
        sprintnum(buf, x[i]);
        outputfile << buf << "\n";
    }

    outputfile.close();

    printf("*** Testing Rlaruv successful ***\n");

    return 0;
}
