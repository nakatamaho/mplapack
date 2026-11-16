#!/usr/bin/env python3
"""Roofline analysis of the GPU kernels: how close each kernel is to the A100's limits.

    roofline.py PEAKS.csv KERNELS.csv GEMM_SWEEP_DIR OUTDIR

PEAKS.csv     output of benchmark/peak_bench.cu (type,op,ops_per_second,note):
              arithmetic throughput of madd (c += a*b), add, mul, div per type,
              and the device-memory bandwidth.
KERNELS.csv   benchmark/results/.../kernel-bench/kernels.csv (GPU kernel-only times); several
              files separated by commas (kernel-bench-2/kernels.csv: only phase = sweep).
GEMM_SWEEP    benchmark/results/.../gemm-sweep/ (<type>_<FORM>.csv).

For every GPU point, the lower bound of the time is

    t_bound = max(t_arith, t_mem)
    t_arith = pairs / madd + (muls - pairs) / mul + (adds - pairs) / add + divs / div,
              pairs = min(muls, adds)            (a multiply-add pair runs at the madd rate)
    t_mem   = bytes / bandwidth                  (each operand read once, each output written once)

with the multiplication and addition counts of LAPACK Working Note 41
(https://netlib.org/lapack/lawnspdf/lawn41.pdf, the operation counts of the
BLAS and LAPACK routines, multiplications and additions separately).  A complex
multiplication counts 4 real multiplications and 2 additions, a complex
addition 2 additions.  Efficiency = t_bound / t_measured (1 = at the limit).
A second number relates every kernel to what GEMM of the same type actually
reaches (practical peak = the GEMM NN plateau of the sweep).

Writes OUTDIR/roofline.csv, OUTDIR/roofline_summary.md and the graphs
OUTDIR/roofline_<type>.png (efficiency against n) and OUTDIR/roofline_summary.png.
Needs matplotlib for the graphs.
"""

import csv
import math
import os
import sys
from collections import defaultdict

NB = 64
TYPES = ["dd", "td", "qd", "mpfr256", "mpfr512", "mpfr1024", "mpfr200rt"]
SIZE = {"dd": 16, "td": 24, "qd": 32}
VECTOR = {"Raxpy", "Rcopy", "iRamax", "iCamax", "Rdot", "Cdotc", "Cdotu", "Rasum", "RCasum", "Rnrm2", "RCnrm2"}  # bytes of one real element on the device


def elem_bytes(t):
    if t in SIZE:
        return SIZE[t]
    if t.startswith("mpfr"):
        bits = 200 if t.endswith("rt") else int(t[4:])
        return 8 * math.ceil(bits / 64) + 16  # limbs + exponent + sign (cu_freal / packed cu_mpfr)
    raise ValueError(t)


# ---- LAWN 41 counts (real arithmetic): (muls, adds, divs, elements read, elements written) ----
def counts(kernel, form, n):
    x = float(n)
    k = kernel[1:] if kernel[0] == "c" else kernel
    if k == "gemm":
        return x ** 3, x ** 3, 0, 2 * x * x + x * x, x * x
    if k == "gemm_update":
        return x * x * NB, x * x * NB, 0, 2 * x * NB + x * x, x * x
    if k in ("syrk", "herk"):
        kk = NB if form.endswith("k64") else x
        tri = x * (x + 1) / 2
        return tri * kk, tri * kk, 0, x * kk + tri, tri
    if k == "trsm":
        m = NB if form.endswith("p64") else x  # the triangle
        other = x
        return other * m * (m + 1) / 2, other * m * (m - 1) / 2, m, m * (m + 1) / 2 + m * other, m * other
    if k == "getrf2":  # LU of n x 64 (LAWN 41 getrf, m >= n)
        mm, nn = x, float(NB)
        muls = 0.5 * nn * (nn * (mm - nn / 3 - 1) + mm) + 2 * nn / 3
        adds = 0.5 * nn * (nn * (mm - nn / 3) - mm) + nn / 6
        return muls, adds, nn, mm * nn, mm * nn
    if k == "potrf2":
        return x ** 3 / 6 + x * x / 2 + x / 3, x ** 3 / 6 - x / 6, x, x * (x + 1) / 2, x * (x + 1) / 2
    if k == "lauum":
        return x ** 3 / 6 + x * x / 2 + x / 3, x ** 3 / 6 - x / 6, 0, x * (x + 1) / 2, x * (x + 1) / 2
    if k == "trtri":
        return x ** 3 / 6 + x * x / 2 + x / 3, x ** 3 / 6 - x * x / 2 + x / 3, x, x * (x + 1) / 2, x * (x + 1) / 2
    if k == "laswp":
        return 0, 0, 0, 2 * NB * x, 2 * NB * x
    return None


# ---- kernel-bench-2 names (goal 11): BLAS1, BLAS2 solves, inverses ----
def counts2(kernel, form, n):
    x = float(n)
    k = kernel.lower()
    band = NB
    if k in ("raxpy",):
        c = (x, x, 0, 2 * x, x)
    elif k in ("rcopy",):
        c = (0, 0, 0, x, x)
    elif k == "iramax":
        c = (0, 0, 0, x, 0)
    elif k == "icamax":
        return 0, x, 0, 2 * x, 0  # |re| + |im|
    elif k in ("rdot", "cdotc", "cdotu"):
        c = (x, x, 0, 2 * x, 0)
    elif k == "rasum":
        c = (0, x, 0, x, 0)
    elif k == "rcasum":
        return 0, 2 * x, 0, 2 * x, 0
    elif k == "rnrm2":
        c = (2 * x, x, 0, x, 0)  # scaled square and sum
    elif k == "rcnrm2":
        return 4 * x, 2 * x, 0, 2 * x, 0
    elif k in ("rtrsv", "ctrsv", "rtpsv", "ctpsv"):
        c = (x * (x - 1) / 2, x * (x - 1) / 2, x, x * (x + 1) / 2 + x, x)
    elif k in ("rtbsv", "ctbsv"):
        c = (x * band, x * band, x, x * (band + 1) + x, x)
    elif k in ("rtrtri", "ctrtri"):
        c = counts("trtri", "", n)
    elif k in ("rlauum", "clauum"):
        c = counts("lauum", "", n)
    elif k in ("rpotri", "cpotri"):
        a, b = counts("trtri", "", n), counts("lauum", "", n)
        c = (a[0] + b[0], a[1] + b[1], a[2], a[3], a[4])
    elif k in ("rgetri", "cgetri"):  # LAWN 41 GETRI
        c = (2 * x ** 3 / 3 + x * x / 2 + 5 * x / 6, 2 * x ** 3 / 3 - 3 * x * x / 2 + 5 * x / 6, x, x * x, x * x)
    else:
        return None
    return complexify(c) if kernel[0] == "C" else c


def complexify(c):
    muls, adds, divs, rd, wr = c
    # complex mul = 4 mul + 2 add; complex add = 2 add; complex division ~ 6 mul + 3 add + 2 div
    return 4 * muls, 2 * muls + 2 * adds, 2 * divs, 2 * rd, 2 * wr


def load_peaks(path):
    p = defaultdict(dict)
    for r in csv.DictReader(open(path)):
        p[r["type"]][r["op"]] = float(r["ops_per_second"])
    return p


def gemm_plateau(sweep, t):
    path = os.path.join(sweep, "%s_NN.csv" % t)
    if not os.path.exists(path):
        return None
    rows = [(int(r["n"]), float(r["gflops"])) for r in csv.DictReader(open(path))]
    return max(g for _, g in rows) * 1e9 if rows else None


def main():
    peaks_csv, kernels_csv, sweep, out = sys.argv[1:5]
    os.makedirs(out, exist_ok=True)
    peaks = load_peaks(peaks_csv)
    bw = peaks["memory"]["copy_bytes"]
    rows = []
    for path in kernels_csv.split(","):
      for r in csv.DictReader(open(path)):
        if r["method"] != "gpu" or r["type"] not in TYPES or r.get("phase", "sweep") != "sweep":
            continue
        t, kern, form, n, sec = r["type"], r["kernel"], r["form"], int(r["n"]), float(r["seconds"])
        if kern[0].isupper() or kern.startswith("i"):
            c = counts2(kern, form, n)
        else:
            c = counts(kern, form, n)
            if c is not None and kern[0] == "c":
                c = complexify(c)
        if c is None or sec <= 0:
            continue
        muls, adds, divs, rd, wr = c
        cplx = 2 if kern[0] == "c" else 1
        t_mem = (rd + wr) * elem_bytes(t) / bw
        pk = peaks.get(t, {})
        if pk:
            pairs = min(muls, adds)
            t_ar = pairs / pk["madd"] + (muls - pairs) / pk["mul"] + (adds - pairs) / pk["add"] + divs / pk["div"]
        else:
            t_ar = float("nan")  # no arithmetic peak (runtime-precision MPFR)
        t_bound = max(t_ar, t_mem) if not math.isnan(t_ar) else float("nan")  # no roofline without an arithmetic peak
        bound = "memory" if t_mem >= (t_ar if not math.isnan(t_ar) else 0) else "arith"
        gp = gemm_plateau(sweep, t)
        practical = ((muls + adds) / sec) / gp if gp and muls + adds > 0 else float("nan")
        rows.append(dict(type=t, kernel=kern, form=form, n=n, seconds=sec, t_arith=t_ar, t_mem=t_mem,
                         bound=bound, efficiency=t_bound / sec, vs_gemm=practical))
    # GEMM sweep itself (NN): efficiency of the GEMM plateau against the madd peak
    for t in TYPES:
        path = os.path.join(sweep, "%s_NN.csv" % t)
        if not os.path.exists(path) or "madd" not in peaks.get(t, {}):
            continue
        for r in csv.DictReader(open(path)):
            n, sec = int(r["n"]), float(r["seconds"])
            if n < 32 or sec <= 0:
                continue
            x = float(n)
            t_ar = x ** 3 / peaks[t]["madd"]
            t_mem = 4 * x * x * elem_bytes(t) / bw
            rows.append(dict(type=t, kernel="gemm_sweep", form="NN", n=n, seconds=sec, t_arith=t_ar, t_mem=t_mem,
                             bound="arith" if t_ar >= t_mem else "memory", efficiency=max(t_ar, t_mem) / sec,
                             vs_gemm=float("nan")))
    with open(os.path.join(out, "roofline.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # summary at a representative size: the largest n <= 2048 (MPFR 1024 / runtime: <= 1024)
    best = {}
    for r in rows:
        lim = 1024 if r["type"] in ("mpfr1024", "mpfr200rt") else 2048
        if r["kernel"] in VECTOR:
            lim = 1 << 21  # vectors: n = 2^21
        if r["n"] > lim:
            continue
        key = (r["type"], r["kernel"], r["form"])
        if key not in best or r["n"] > best[key]["n"]:
            best[key] = r
    with open(os.path.join(out, "roofline_summary.md"), "w") as f:
        f.write("| type | kernel | form | n | bound | % of roofline | % of GEMM plateau |\n|---|---|---|---:|---|---:|---:|\n")
        for t in TYPES:
            for key in sorted(k for k in best if k[0] == t):
                r = best[key]
                f.write("| %s | %s | %s | %d | %s | %s | %s |\n" % (t, r["kernel"], r["form"], r["n"], r["bound"],
                        "-" if math.isnan(r["efficiency"]) else "%.1f" % (100 * r["efficiency"]),
                        "" if math.isnan(r["vs_gemm"]) else "%.1f" % (100 * r["vs_gemm"])))
    try:
        plot(rows, best, out)
    except ImportError:
        print("matplotlib not found: no graphs", file=sys.stderr)


def plot(rows, best, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 11})
    page, grid = "#f9f9f7", "#e1e0d9"
    for t in TYPES:
        series = defaultdict(list)
        for r in rows:
            if r["type"] == t:
                if not math.isnan(r["efficiency"]):
                    series[(r["kernel"], r["form"])].append((r["n"], 100 * r["efficiency"]))
        if not series:
            continue
        fig, ax = plt.subplots(figsize=(13.8, 7.2), dpi=100)
        fig.patch.set_facecolor(page)
        ax.set_facecolor(page)
        for (k, fm), pts in sorted(series.items()):
            pts.sort()
            ax.plot([p[0] for p in pts], [p[1] for p in pts], linewidth=1.6, label="%s %s" % (k, fm))
        ax.axhline(100, color="#52514e", linewidth=1)
        ax.set_yscale("log")
        ax.set_xscale("log", base=2)
        ax.grid(True, color=grid, which="both")
        ax.set_xlabel("n")
        ax.set_ylabel("ルーフラインに対する割合 (%)")
        ax.set_title("%s: GPU カーネル単体の性能 / ルーフライン(演算と帯域の上限、A100)" % t, loc="left", fontsize=15, fontweight="bold")
        ax.legend(loc="lower right", ncol=3, fontsize=8, frameon=False)
        fig.tight_layout()
        fig.savefig(os.path.join(out, "roofline_%s.png" % t), facecolor=page)
        plt.close(fig)
    # summary heat table
    kernels = sorted({(k[1], k[2]) for k in best})
    types = [t for t in TYPES if any(k[0] == t and not math.isnan(best[k]["efficiency"]) for k in best)]
    fig, ax = plt.subplots(figsize=(1.0 + 0.55 * len(kernels), 0.6 + 0.5 * len(types)), dpi=100)
    data = [[100 * best[(t,) + k]["efficiency"] if (t,) + k in best else float("nan") for k in kernels] for t in types]
    im = ax.imshow(data, cmap="viridis", vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(kernels)))
    ax.set_xticklabels(["%s %s" % k for k in kernels], rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(len(types)))
    ax.set_yticklabels(types)
    for i, row in enumerate(data):
        for j, v in enumerate(row):
            if not math.isnan(v):
                ax.text(j, i, "%.0f" % v, ha="center", va="center", fontsize=7, color="white" if v < 60 else "black")
    ax.set_title("ルーフラインに対する割合 (%)、n = 2048(MPFR 1024 と実行時精度は 1024)", loc="left", fontsize=11)
    fig.colorbar(im, ax=ax, fraction=0.02)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "roofline_summary.png"))
    plt.close(fig)


if __name__ == "__main__":
    main()
