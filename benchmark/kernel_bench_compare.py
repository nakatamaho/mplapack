#!/usr/bin/env python3
"""Bit comparison of the result dumps of Rkernel_bench (GPU against CPU).

    kernel_bench_compare.py DUMPDIR [CSV]

Rkernel_bench --dump DIR writes <label>_<kernel>_<form>_<n>_<gpu|cpu>.txt:
"info", "ipiv" lines, then one element per line in hexadecimal (dd/td/qd:
the limbs with %a, comma separated; MPFR: %Ra; complex: real;imag).  For
every pair this prints label,kernel,form,n,elements,mismatches,zero_sign,
status, where zero_sign counts the elements that differ only in the sign
of a zero (limb), and mismatches every other difference (info and ipiv
included).  status: identical, zero-sign, DIFFERENT or missing.
"""

import csv
import glob
import os
import sys


def is_zero(token):
    t = token.lstrip("+-")
    return t in ("0x0p+0", "0x0.0p+0", "0") or (t.startswith("0x0p") or t.startswith("0x0.p"))


def compare_line(a, b):
    """0: equal, 1: only signed zeros differ, 2: different."""
    if a == b:
        return 0
    ta = a.replace(";", ",").split(",")
    tb = b.replace(";", ",").split(",")
    if len(ta) != len(tb):
        return 2
    for x, y in zip(ta, tb):
        if x != y and not (is_zero(x) and is_zero(y)):
            return 2
    return 1


def compare(gpu, cpu):
    with open(gpu) as f:
        g = f.read().splitlines()
    with open(cpu) as f:
        c = f.read().splitlines()
    if len(g) != len(c):
        return len(c), max(len(g), len(c)), 0
    elements = mismatches = zero_sign = 0
    for a, b in zip(g, c):
        header = a.startswith("info") or a.startswith("ipiv")
        if not header:
            elements += 1
        r = compare_line(a, b) if not header else (0 if a == b else 2)
        if r == 1:
            zero_sign += 1
        elif r == 2:
            mismatches += 1
    return elements, mismatches, zero_sign


def main():
    d = sys.argv[1]
    rows = []
    for gpu in sorted(glob.glob(os.path.join(d, "*_gpu.txt"))):
        cpu = gpu[: -len("_gpu.txt")] + "_cpu.txt"
        stem = os.path.basename(gpu)[: -len("_gpu.txt")]
        # label_kernel_form_n; the kernel may contain "_" (gemm_update), the form too (_p64, _k64)
        parts = stem.split("_")
        label, n = parts[0], parts[-1]
        rest = parts[1:-1]
        if rest[0] == "gemm" and len(rest) > 1 and rest[1] == "update":
            kernel, form = "gemm_update", "_".join(rest[2:])
        else:
            kernel, form = rest[0], "_".join(rest[1:])
        if not os.path.exists(cpu):
            rows.append([label, kernel, form, n, "", "", "", "missing"])
            continue
        e, m, z = compare(gpu, cpu)
        status = "DIFFERENT" if m else ("zero-sign" if z else "identical")
        rows.append([label, kernel, form, n, e, m, z, status])
    out = open(sys.argv[2], "w", newline="") if len(sys.argv) > 2 else sys.stdout
    w = csv.writer(out)
    w.writerow(["label", "kernel", "form", "n", "elements", "mismatches", "zero_sign", "status"])
    w.writerows(rows)


if __name__ == "__main__":
    main()
