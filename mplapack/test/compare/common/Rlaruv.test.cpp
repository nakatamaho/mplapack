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
#include <cstring>
#include <array>
#include <iomanip>

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
    std::fprintf(stderr, "Usage: %s [iseed1 iseed2 iseed3 iseed4]\n", prog);
    std::fprintf(stderr, "  If no seeds are given, defaults to 4 3 2 1.\n");
}

int main(int argc, char *argv[]) {

    std::fprintf(stderr, "*** Testing Rlaruv start ***\n");

    // Default seed
    INTEGER seed[4] = {4, 3, 2, 1};

    // Parse optional command-line seeds
    if (argc != 1 && argc != 5) {
        usage(argv[0]);
        return 1;
    }
    if (argc == 5) {
        long tmp = 0;
        for (int i = 0; i < 4; ++i) {
            if (!parse_int(argv[i + 1], tmp) || tmp < 0 || tmp > INT_MAX) {
                std::fprintf(stderr, "Error: invalid seed argument '%s'\n", argv[i + 1]);
                return 1;
            }
            seed[i] = static_cast<INTEGER>(tmp);
        }
    }

    std::ofstream outputfile(outname, std::ios::binary);
    if (!outputfile.is_open()) {
        std::fprintf(stderr, "Error: cannot open output file '%s'\n", outname);
        return 1;
    }

    // Write a small header to the output file (no stdout)
    outputfile << "# Rlaruv test output\n";
    outputfile << "# iseed = " << static_cast<long>(seed[0]) << " " << static_cast<long>(seed[1]) << " " << static_cast<long>(seed[2]) << " " << static_cast<long>(seed[3]) << "\n";
    outputfile << "# Format: one random number per line\n";

    // -------------------------------------------------------------------------
    // Simple statistical sanity checks for the generated U(0,1) sequence.
    // Requirements:
    //   - Mean/variance via Welford's online algorithm
    //   - Track min/max
    //   - Build a 100-bin histogram on [0,1) and compute a simple chi-square value
    //
    // Note: We rely on conversion to double for accumulation. This is intended as a
    // lightweight smoke test, not a high-precision distribution proof.
    // -------------------------------------------------------------------------
    long long stats_n = 0;
    double stats_mean = 0.0;
    double stats_M2 = 0.0;
    double stats_min = std::numeric_limits<double>::infinity();
    double stats_max = -std::numeric_limits<double>::infinity();
    std::array<long long, 100> stats_hist{};
    stats_hist.fill(0);

    auto stats_update = [&](const REAL &v) {
        const double xd = cast2double(v);
        const double x = xd;

        // Welford update (population moments)
        ++stats_n;
        const double delta = x - stats_mean;
        stats_mean += delta / static_cast<double>(stats_n);
        const double delta2 = x - stats_mean;
        stats_M2 += delta * delta2;

        if (x < stats_min)
            stats_min = x;
        if (x > stats_max)
            stats_max = x;

        // Histogram binning for U(0,1). Clamp out-of-range values defensively.
        int bin = 0;
        if (x <= 0.0) {
            bin = 0;
        } else if (x >= 1.0) {
            bin = 99;
        } else {
            bin = static_cast<int>(x * 100.0); // x in (0,1) -> [0,99]
            if (bin < 0)
                bin = 0;
            if (bin > 99)
                bin = 99;
        }
        ++stats_hist[static_cast<size_t>(bin)];
    };

    // Optional: generate small vectors (len=1..MAX_N-1) and dump them
    for (int len = 1; len < MAX_N; ++len) {
        std::vector<REAL> x(static_cast<size_t>(len));
        Rlaruv(seed, len, x.data());
        for (int i = 0; i < len; ++i) {
            stats_update(x[static_cast<size_t>(i)]);
            char buf[4096];
            sprintnum(buf, x[i]);
            // Trim a leading '+' to keep output stable and minimal.
            if (buf[0] == '+') {
                std::memmove(buf, buf + 1, std::strlen(buf));
            }
            outputfile << buf << "\n";
        }
    }

    // Large output block
    const int big_n = 100000;
    std::vector<REAL> x(big_n);
    Rlaruv(seed, big_n, x.data());
    for (int i = 0; i < big_n; ++i) {
        stats_update(x[static_cast<size_t>(i)]);
        char buf[4096];
        sprintnum(buf, x[i]);
        // Trim a leading '+' to keep output stable and minimal.
        if (buf[0] == '+') {
            std::memmove(buf, buf + 1, std::strlen(buf));
        }
        outputfile << buf << "\n";
    }

    // -------------------------------------------------------------------------
    // Append stats to the end of the same output file (stdout is forbidden).
    // Theory for U(0,1): mean = 0.5, variance(population) = 1/12.
    // -------------------------------------------------------------------------
    const double stats_theory_mean = 0.5;
    const double stats_theory_var_pop = 1.0 / 12.0;

    const double stats_var_pop = (stats_n > 0) ? (stats_M2 / static_cast<double>(stats_n)) : std::numeric_limits<double>::quiet_NaN();
    const double stats_var_sample = (stats_n > 1) ? (stats_M2 / static_cast<double>(stats_n - 1)) : std::numeric_limits<double>::quiet_NaN();

    double stats_chi2 = std::numeric_limits<double>::quiet_NaN();
    if (stats_n > 0) {
        const double expected = static_cast<double>(stats_n) / 100.0;
        if (expected > 0.0) {
            stats_chi2 = 0.0;
            for (size_t b = 0; b < stats_hist.size(); ++b) {
                const double obs = static_cast<double>(stats_hist[b]);
                const double diff = obs - expected;
                stats_chi2 += (diff * diff) / expected;
            }
        }
    }

    outputfile << std::setprecision(10);
    outputfile << "# stats_n " << stats_n << "\n";
    outputfile << "# stats_mean " << stats_mean << " theory " << stats_theory_mean << " diff " << (stats_mean - stats_theory_mean) << "\n";
    outputfile << "# stats_var_pop " << stats_var_pop << " theory " << stats_theory_var_pop << " diff " << (stats_var_pop - stats_theory_var_pop) << "\n";
    outputfile << "# stats_var_sample " << stats_var_sample << "\n";
    outputfile << "# stats_min " << stats_min << "\n";
    outputfile << "# stats_max " << stats_max << "\n";
    outputfile << "# stats_chi2_100bins " << stats_chi2 << "\n";

    outputfile << "# stats_hist_100bins";
    for (size_t b = 0; b < stats_hist.size(); ++b) {
        outputfile << " " << stats_hist[b];
    }
    outputfile << "\n";

    outputfile.close();

    std::fprintf(stderr, "*** Testing Rlaruv successful ***\n");

    return 0;
}
