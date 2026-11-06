#!/usr/bin/env python3
"""MPFR (cu_freal) arithmetic peaks for benchmark/roofline.py, in two levels.

    mpfr_level_peaks.py PEAKS.csv PATHS.csv MIX.txt LEVEL2_OUT.csv LEVEL1_OUT.csv

PEAKS.csv  benchmark/peak_bench.cu output (dd/td/qd, memory, MPFR div are kept)
PATHS.csv  benchmark/mpfr_path_bench.cu output (limb primitives, cu_freal paths)
MIX.txt    benchmark/mpfr_gemm_paths.py output (share of each addition path)

Level 2 (cu_freal as implemented): mul = cu_fmul rate with full-width random
significands; add = the path rates weighted by MIX (harmonic: time per add =
sum share/rate); madd = 1 / (1/mul + 1/add) (cu_freal has no fused operation).
Level 1 (hardware): mul = mad64 / N^2 (full schoolbook product, N = PB/64 limbs),
add = 1 / (N/shf64 + N/addc64), madd likewise.  div is not modelled; both levels
keep the peak_bench div (an upper bound on the rate, so the time bound stays a
lower bound only approximately for division-heavy kernels: trsm, potrf2, inverses).
The other rows of PEAKS.csv are copied unchanged.
"""
import ast, csv, sys

peaks_in, paths_in, mix_in, out2, out1 = sys.argv[1:6]
rate = {}
for r in csv.DictReader(open(paths_in)):
    rate[(int(r["pb"]), r["op"])] = float(r["ops_per_second"])
mix = {}
for line in open(mix_in):
    pb, d = line.split(" ", 1)
    mix[int(pb)] = ast.literal_eval(d.strip())
rows = list(csv.DictReader(open(peaks_in)))


def level(pb, which):
    n = pb // 64
    if which == 2:
        mul = rate[(pb, "mul")]
        w = {k: v for k, v in mix[pb].items() if k.startswith("add_")}
        tot = sum(w.values())
        add = tot / sum(v / rate[(pb, k)] for k, v in w.items())
    else:
        mul = rate[(0, "mad64")] / n ** 2
        add = 1 / (n / rate[(0, "shf64")] + n / rate[(0, "addc64")])
    return {"mul": mul, "add": add, "madd": 1 / (1 / mul + 1 / add)}


for out, which in ((out2, 2), (out1, 1)):
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["type", "op", "ops_per_second", "note"])
        for r in rows:
            t = r["type"]
            if t.startswith("mpfr") and r["op"] in ("madd", "add", "mul"):
                v = level(int(t[4:]), which)[r["op"]]
                note = ("level 2: cu_freal, full-width significands, GEMM sign/gap mix" if which == 2
                        else "level 1: hardware limb operations (full product)")
                w.writerow([t, r["op"], "%.6e" % v, note])
            else:
                w.writerow([t, r["op"], r["ops_per_second"], r["note"]])
