#!/usr/bin/env python3
"""Run the kernel benchmark (Rkernel_bench.*): CPU against GPU, per type.

    kernel_bench.py probe  --bindir B --out DIR [--types ...]
    kernel_bench.py run    --bindir B --out DIR [--types ...] [--plan DIR/plan.json] [--only gpu|cpu]
    kernel_bench.py check  --bindir B --out DIR [--types ...]
    kernel_bench.py merge  --out DIR

B is the directory of the programs (benchmark/ of an autotools build, or the
CMake build's benchmark/).  DIR receives:

  raw/<type>_<gpu|cpu>.csv   one line per point (label,kernel,form,n,method,
                             seconds,gflops,runs), appended as they finish
  raw/<type>_<gpu|cpu>.log   stderr of the programs and the commands run
  gpu-processes.log          nvidia-smi before and after every GPU series
  probe.csv, plan.json       "probe": two sizes per series, the time model and
                             the size lists chosen ("run" uses them)
  dumps/, checks.csv         "check": results at n = 64 and 200 compared bit
                             for bit (kernel_bench_compare.py)
  kernels.csv                "merge": every point, with the bit comparison and
                             the GPU-busy flag of its series

GPU: one warm-up run and the best of 3 (one run when the warm-up exceeds
10 s).  CPU: best of up to 3 (one run when a run exceeds 10 s), OpenMP with
OMP_NUM_THREADS (default 32).  "run" executes the GPU series and the CPU
series in two concurrent streams (--only selects one): the CPU programs do
not use the GPU, the GPU programs use one host thread.
"""

import argparse
import collections
import concurrent.futures
import csv
import json
import math
import os
import subprocess
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))

# type: CPU program, GPU program, environment, GPU size limit, CPU size limit
TYPES = {
    "dd": ("Rkernel_bench.dd_opt", "Rkernel_bench.dd_cuda_kernel", {}, 4096, 4096),
    "td": ("Rkernel_bench.td_opt", "Rkernel_bench.td_cuda_kernel", {}, 4096, 4096),
    "qd": ("Rkernel_bench.qd_opt", "Rkernel_bench.qd_cuda_kernel", {}, 4096, 4096),
    "mpfr256": ("Rkernel_bench.mpfr_opt", "Rkernel_bench.mpfr_cuda_kernel", {"MPFRXX_DEFAULT_PRECISION_BITS": "256"}, 4096, 2048),
    "mpfr512": ("Rkernel_bench.mpfr_opt", "Rkernel_bench.mpfr_cuda_kernel", {"MPFRXX_DEFAULT_PRECISION_BITS": "512"}, 4096, 2048),
    "mpfr1024": ("Rkernel_bench.mpfr_opt", "Rkernel_bench.mpfr_cuda_kernel", {"MPFRXX_DEFAULT_PRECISION_BITS": "1024"}, 2048, 2048),
    "mpfr200rt": ("Rkernel_bench.mpfr_opt", "Rkernel_bench.mpfr_cuda_kernel", {"MPFRXX_DEFAULT_PRECISION_BITS": "200"}, 2048, 2048),
}
ALL_TYPES = list(TYPES)

# (kernel, form): the shapes of the factorizations (see Rkernel_bench.cpp)
REAL_SERIES = [
    ("gemm", "NN"), ("gemm", "NT"), ("gemm", "TN"), ("gemm_update", "NN"),
    ("syrk", "LN"), ("syrk", "UT"), ("syrk", "LN_k64"), ("syrk", "UT_k64"),
    ("trsm", "LLNU"), ("trsm", "RLTN"), ("trsm", "LUTN"), ("trsm", "LUNN"),
    ("trsm", "LLNU_p64"), ("trsm", "RLTN_p64"), ("trsm", "LUTN_p64"), ("trsm", "LUNN_p64"),
    ("getrf2", "P64"), ("potrf2", "L"), ("laswp", "R64"), ("trtri", "UN"), ("lauum", "L"),
]
# complex: GEMM, HERK, TRSM in the shapes of the complex LU / Cholesky (time
# budget; the other forms of Rkernel_bench work too, see --series)
COMPLEX_SERIES = [
    ("cgemm", "NN"), ("cgemm", "NC"), ("cherk", "LN"), ("cherk", "LN_k64"),
    ("ctrsm", "LLNU"), ("ctrsm", "RLCN_p64"),
]
SERIES = REAL_SERIES + COMPLEX_SERIES
# MPFR: no GPU laswp (not reachable alone), no GPU inverses yet
MPFR_NO_GPU = {("laswp", "R64"), ("trtri", "UN"), ("lauum", "L")}
POTRF2_SIZES = [64, 256]
CHECK_SIZES = [64, 200]


def is_mpfr(t):
    return t.startswith("mpfr")


def big_sweep(kernel, form):
    """dd/td/qd GEMM and SYRK/HERK (k = n) go to 8192."""
    return kernel in ("gemm", "cgemm", "syrk", "cherk") and not form.endswith("_k64")


# series measured on the fine grid of the specification; the others on a
# coarser grid (time budget, see README.kernel_bench.md)
FINE = {("gemm", "NN"), ("gemm_update", "NN"), ("syrk", "LN"), ("syrk", "LN_k64"), ("trsm", "LLNU"),
        ("trsm", "RLTN_p64"), ("getrf2", "P64")}


def gpu_grid(t, kernel, form, top=None, coarse=False):
    if kernel == "potrf2":
        return list(POTRF2_SIZES)
    limit = top if top else TYPES[t][3]
    s1, s2, s3 = (64, 512 if is_mpfr(t) else 256, 512) if coarse else (32, 128, 256)
    sizes = list(range(32, min(1024, limit) + 1, s1)) + list(range(1024 + s2, min(4096, limit) + 1, s2))
    if not is_mpfr(t) and big_sweep(kernel, form) and top is None:
        sizes += list(range(4096 + s3, 8192 + 1, s3))
    top_n = max(sizes)
    sizes += [1 << e for e in range(5, 14) if (1 << e) <= top_n]  # every power of two
    return sorted(set(sizes))


def cpu_grid(t, kernel, form, top=None, coarse=False):
    if kernel == "potrf2":
        return list(POTRF2_SIZES)
    limit = top if top else TYPES[t][4]
    return [n for n in (128, 256, 512, 1024, 2048, 4096) if n <= limit]


def flops(kernel, form, n):
    """the counts of Rkernel_bench.cpp"""
    k = kernel[1:] if kernel[0] == "c" else kernel
    f = 4.0 if kernel[0] == "c" else 1.0
    x = float(n)
    if k == "gemm":
        return 2 * f * x ** 3
    if k == "gemm_update":
        return 2 * f * x * x * 64
    if k in ("syrk", "herk"):
        kk = 64 if form.endswith("_k64") else x
        return f * kk * x * (x + 1)
    if k == "trsm":
        tri = 64 if form.endswith("_p64") else x
        return f * x * tri * tri
    if k == "getrf2":
        m, nn = x, 64.0
        if m < nn:
            return 0.5 * m * (m * (nn - m / 3 - 1) + nn) + 2 * m / 3 + 0.5 * m * (m * (nn - m / 3) - nn) + m / 6
        return 0.5 * nn * (nn * (m - nn / 3 - 1) + m) + 2 * nn / 3 + 0.5 * nn * (nn * (m - nn / 3) - m) + nn / 6
    if k in ("potrf2", "lauum"):
        return x ** 3 / 3 + x * x / 2 + x / 6
    if k == "trtri":
        return x ** 3 / 3 + 2 * x / 3
    return 0.0


def series_for(t, method):
    out = []
    for kernel, form in SERIES:
        if method == "gpu" and is_mpfr(t) and (kernel, form) in MPFR_NO_GPU:
            continue
        out.append((kernel, form))
    return out


class Runner:
    def __init__(self, a):
        self.bindir = os.path.abspath(a.bindir)
        self.out = os.path.abspath(a.out)
        self.threads = str(a.threads)
        self.pin = getattr(a, "pin", False)
        os.makedirs(os.path.join(self.out, "raw"), exist_ok=True)
        self.lock = threading.Lock()

    def env(self, t, method):
        e = dict(os.environ)
        for k in list(e):
            if k.startswith("MPLAPACK_") or k.startswith("MPFRXX_"):
                del e[k]
        e.update(TYPES[t][2])
        # libtool wrappers add only the MPLAPACK libraries: the bundled QD/MPFR/GMP/MPC too
        top = os.path.dirname(self.bindir)
        extra = [os.path.join(top, "external", "i", x, "lib") for x in ("QD", "MPFR", "GMP", "MPC")] + ["/usr/local/cuda/lib64"]
        extra = [d for d in extra if os.path.isdir(d)]
        e["LD_LIBRARY_PATH"] = ":".join(extra + ([e["LD_LIBRARY_PATH"]] if e.get("LD_LIBRARY_PATH") else []))
        # the GPU programs run one host thread (packing in the MPFR entries may use OpenMP)
        e["OMP_NUM_THREADS"] = self.threads if method == "cpu" else "1"
        return e

    def program(self, t, method):
        return os.path.join(self.bindir, TYPES[t][1 if method == "gpu" else 0])

    def log(self, name, text):
        with self.lock:
            with open(os.path.join(self.out, name), "a") as f:
                f.write(text)

    def gpu_processes(self, tag):
        q = subprocess.run(["nvidia-smi", "--query-compute-apps=pid,process_name,used_memory", "--format=csv,noheader"],
                           capture_output=True, text=True).stdout.strip()
        g = subprocess.run(["nvidia-smi", "--query-gpu=utilization.gpu,memory.used", "--format=csv,noheader"],
                           capture_output=True, text=True).stdout.strip()
        others = [l for l in q.splitlines() if l.strip() and "Rkernel_bench" not in l]
        self.log("gpu-processes.log", "%s %s | gpu %s | apps: %s\n" % (time.strftime("%Y-%m-%dT%H:%M:%S"), tag, g, q.replace("\n", "; ") or "none"))
        return others

    def invoke(self, t, method, kernel, form, sizes, extra=(), csvname=None):
        cmd = [self.program(t, method), kernel, form, ",".join(str(n) for n in sizes), "--label", t] + list(extra)
        if method == "cpu" and self.pin:
            # cpus 0..threads-1 (not used for the results: on a shared host the
            # barrier-heavy panels then wait for descheduled threads)
            cmd = ["taskset", "-c", "0-%d" % (int(self.threads) - 1)] + cmd
        start = time.time()
        p = subprocess.run(cmd, env=self.env(t, method), capture_output=True, text=True)
        elapsed = time.time() - start
        name = csvname or "raw/%s_%s.csv" % (t, method)
        self.log(name, p.stdout)
        self.log("raw/%s_%s.log" % (t, method), "$ %s\n%s# exit %d, %.1f s\n" % (" ".join(cmd), p.stderr, p.returncode, elapsed))
        return p.returncode, p.stdout, elapsed


def parse_rows(text):
    rows = []
    for line in text.splitlines():
        f = line.split(",")
        if len(f) == 8 and f[5]:
            rows.append({"label": f[0], "kernel": f[1], "form": f[2], "n": int(f[3]), "method": f[4],
                         "seconds": float(f[5]), "gflops": f[6], "runs": int(f[7])})
    return rows


# ---- probe and plan ----

def probe_sizes(t, method, kernel=""):
    if t == "mpfr200rt" and kernel.startswith("c"):
        return (64, 128) if kernel == "cgemm" else (32, 64)  # the complex cu_mpfr kernels are very slow
    if method == "gpu":
        return (128, 512) if is_mpfr(t) else (256, 1024)
    return (128, 256) if is_mpfr(t) else (256, 512)


def model(points):
    """time(n) = a n^b through two points; b clamped to [1, 3.3]"""
    (n1, t1), (n2, t2) = points
    b = math.log(max(t2, 1e-9) / max(t1, 1e-9)) / math.log(n2 / n1)
    b = min(max(b, 1.0), 3.3)
    return t2, n2, b


def predict(m, n):
    """the fit up to the larger probe size, at least cubic above it (calibrated
    against the first measured series: within about 20 % on the GPU)"""
    t2, n2, b = m
    return t2 * (n / n2) ** (b if n <= n2 else max(b, 3.0))


def point_cost(t_run, method, mpfr_overhead):
    """wall time of one point: GPU warm-up + 3 (or warm-up + 1), CPU up to 3"""
    if method == "gpu":
        runs = 4 if t_run <= 10 else 2
    else:
        runs = 3 if t_run <= 10 else 1
    return runs * (t_run + mpfr_overhead)


def mpfr_host_overhead(t, method, n):
    """packing / restoring n x n MPFR operands per run (measured order of magnitude)"""
    if not is_mpfr(t):
        return 0.0
    bits = int(TYPES[t][2]["MPFRXX_DEFAULT_PRECISION_BITS"])
    per = 6e-8 * (1 + bits / 512.0) * (2.0 if method == "gpu" else 1.0)
    return per * n * n


def estimate(probe, t, method, kernel, form, sizes):
    key = "%s|%s|%s|%s" % (t, method, kernel, form)
    if key not in probe:
        return float("nan")
    m = model(probe[key])
    total = 0.0
    for n in sizes:
        total += point_cost(predict(m, n), method, mpfr_host_overhead(t, method, n)) + 0.1
    return total


def cmd_probe(a):
    r = Runner(a)
    out = {}
    path = os.path.join(r.out, "probe.csv")
    if os.path.exists(path):
        for row in csv.DictReader(open(path)):
            out.setdefault("%s|%s|%s|%s" % (row["type"], row["method"], row["kernel"], row["form"]), []).append((int(row["n"]), float(row["seconds"])))
    else:
        with open(path, "w") as f:
            f.write("type,method,kernel,form,n,seconds\n")
    for t in a.types:
        for method in a.methods:
            for kernel, form in series_for(t, method):
                key = "%s|%s|%s|%s" % (t, method, kernel, form)
                if key in out:
                    continue
                sizes = POTRF2_SIZES if kernel == "potrf2" else probe_sizes(t, method, kernel)
                extra = ["--reps", "1"] + (["--warmup", "1"] if method == "gpu" else [])
                if method == "gpu":
                    r.gpu_processes("probe %s %s %s start" % (t, kernel, form))
                rc, text, el = r.invoke(t, method, kernel, form, sizes, extra, csvname="raw/probe_%s_%s.csv" % (t, method))
                rows = parse_rows(text)
                pts = [(x["n"], x["seconds"]) for x in rows]
                with open(path, "a") as f:
                    for n, s in pts:
                        f.write("%s,%s,%s,%s,%d,%.6e\n" % (t, method, kernel, form, n, s))
                out[key] = pts
                print("probe %-9s %s %-12s %-9s %s  (%.1f s)" % (t, method, kernel, form, " ".join("%d:%.3g" % p for p in pts), el), flush=True)


def load_probe(outdir):
    out = {}
    for row in csv.DictReader(open(os.path.join(outdir, "probe.csv"))):
        out.setdefault("%s|%s|%s|%s" % (row["type"], row["method"], row["kernel"], row["form"]), []).append((int(row["n"]), float(row["seconds"])))
    return {k: sorted(v) for k, v in out.items() if len(v) >= 2}


def cmd_plan(a):
    """size lists of every series and the estimated time; --top caps the sizes"""
    probe = load_probe(a.out)
    plan = {"series": [], "note": a.note}
    totals = {}
    for t in a.types:
        for method in ("gpu", "cpu"):
            for kernel, form in series_for(t, method):
                top = None
                for spec in a.top or []:
                    # TYPE:METHOD:KERNEL[:FORM]=N   ("*" matches all)
                    lhs, n = spec.split("=")
                    f = (lhs.split(":") + ["*"] * 4)[:4]
                    if all(x == "*" or x == y for x, y in zip(f, (t, method, kernel, form))):
                        top = int(n)
                coarse = a.coarse and (kernel, form) not in FINE
                sizes = (gpu_grid if method == "gpu" else cpu_grid)(t, kernel, form, top, coarse)
                key = "%s|%s|%s|%s" % (t, method, kernel, form)
                cap = a.max_point_fine if (kernel, form) in FINE else a.max_point
                if (kernel, form) == ("gemm", "NN") or (not is_mpfr(t) and (kernel, form) == ("syrk", "LN")):
                    cap = 0  # the full sweeps of the specification (dd/td/qd: 8192)
                if cap and method == "gpu" and key in probe and kernel != "potrf2":
                    # drop the sizes whose predicted run exceeds the cap (seconds)
                    m = model(probe[key])
                    kept = [n for n in sizes if predict(m, n) <= cap]
                    sizes = kept if len(kept) >= 2 else sizes[:2]
                est = estimate(probe, t, method, kernel, form, sizes)
                plan["series"].append({"type": t, "method": method, "kernel": kernel, "form": form, "sizes": sizes, "estimate": est})
                totals[(t, method)] = totals.get((t, method), 0) + (0 if math.isnan(est) else est)
    plan["totals"] = {"%s %s" % k: v for k, v in totals.items()}
    gpu = sum(v for k, v in totals.items() if k[1] == "gpu")
    cpu = sum(v for k, v in totals.items() if k[1] == "cpu")
    plan["gpu_hours"], plan["cpu_hours"] = gpu / 3600, cpu / 3600
    json.dump(plan, open(os.path.join(a.out, "plan.json"), "w"), indent=1)
    for k, v in sorted(totals.items()):
        print("%-10s %s %7.1f min" % (k[0], k[1], v / 60))
    print("GPU stream %.2f h, CPU stream %.2f h" % (gpu / 3600, cpu / 3600))
    if a.verbose:
        for s in plan["series"]:
            print("%-9s %s %-12s %-9s %5.1f min  n<=%d" % (s["type"], s["method"], s["kernel"], s["form"], s["estimate"] / 60, max(s["sizes"])))


# ---- run ----

def cmd_run(a):
    r = Runner(a)
    plan_path = a.plan or os.path.join(r.out, "plan.json")
    plan = json.load(open(plan_path))
    done = set()
    for t in TYPES:
        for method in ("gpu", "cpu"):
            p = os.path.join(r.out, "raw", "%s_%s.csv" % (t, method))
            if os.path.exists(p):
                for row in parse_rows(open(p).read()):
                    done.add((t, method, row["kernel"], row["form"], row["n"]))
    jobs = {"gpu": [], "cpu": []}
    deadline = time.time() + a.hours * 3600 if a.hours else None

    def priority(s):
        k = (s["kernel"], s["form"])
        return (0 if k in FINE else 1 if not s["kernel"].startswith("c") else 2, list(TYPES).index(s["type"]))
    for s in sorted(plan["series"], key=priority):
        if s["type"] not in a.types or (a.only and s["method"] != a.only):
            continue
        sizes = [n for n in s["sizes"] if (s["type"], s["method"], s["kernel"], s["form"], n) not in done]
        if sizes:
            jobs[s["method"]].append((s, sizes))

    def stream(method):
        for s, sizes in jobs[method]:
            t, kernel, form = s["type"], s["kernel"], s["form"]
            tag = "%s %s %s" % (t, kernel, form)
            if deadline and time.time() > deadline:
                r.log("omitted.txt", "%s %s (deadline)\n" % (method, tag))
                continue
            # the plan may be revised while running (sizes of the series not started yet)
            current = {(x["type"], x["method"], x["kernel"], x["form"]): x for x in json.load(open(plan_path))["series"]}
            if (t, method, kernel, form) not in current:
                r.log("omitted.txt", "%s %s (removed from the plan)\n" % (method, tag))
                continue
            sizes = [n for n in current[(t, method, kernel, form)]["sizes"] if n in sizes or n not in s["sizes"]]
            if not sizes:
                continue
            busy = False
            if method == "gpu":
                busy = bool(r.gpu_processes("before " + tag))
            rc, text, el = r.invoke(t, method, kernel, form, sizes)
            if method == "gpu":
                busy = bool(r.gpu_processes("after " + tag)) or busy
                if busy:
                    r.log("gpu-busy.txt", "%s\n" % tag)
            print("%s %-4s %-9s %-12s %-9s %3d points %8.1f s (est %.1f)%s" % (time.strftime("%H:%M:%S"), method, t, kernel, form, len(sizes), el, s["estimate"],
                                                                            "  OTHER GPU PROCESS" if busy else ""), flush=True)

    with concurrent.futures.ThreadPoolExecutor(2) as ex:
        futures = [ex.submit(stream, m) for m in ("gpu", "cpu") if jobs[m]]
        for f in futures:
            f.result()


# ---- bit comparison ----

def cmd_check(a):
    r = Runner(a)
    d = os.path.join(r.out, "dumps")
    os.makedirs(d, exist_ok=True)
    for t in a.types:
        for method in ("gpu", "cpu"):
            for kernel, form in series_for(t, method):
                if is_mpfr(t) and (kernel, form) in MPFR_NO_GPU:
                    continue  # nothing to compare
                extra = ["--reps", "1", "--warmup", "0", "--dump", d, "--dump-sizes", ",".join(map(str, CHECK_SIZES))]
                if method == "gpu":
                    r.gpu_processes("check %s %s %s" % (t, kernel, form))
                rc, text, el = r.invoke(t, method, kernel, form, CHECK_SIZES, extra, csvname="raw/check_%s_%s.csv" % (t, method))
                print("check %-9s %s %-12s %-9s rc=%d %.1f s" % (t, method, kernel, form, rc, el), flush=True)
    subprocess.run([sys.executable, os.path.join(HERE, "kernel_bench_compare.py"), d, os.path.join(r.out, "checks.csv")], check=True)


# ---- merge ----

def cmd_merge(a):
    """kernels.csv: the GPU points of raw/, the CPU points of the first source
    of --cpu-sources (NAME=DIR,...) that has the series (default raw/)"""
    out = os.path.abspath(a.out)
    checks = {}
    p = os.path.join(out, "checks.csv")
    if os.path.exists(p):
        for row in csv.DictReader(open(p)):
            key = (row["label"], row["kernel"], row["form"])
            checks.setdefault(key, []).append("n=%s:%s" % (row["n"], row["status"] if row["status"] != "zero-sign" else "zero-sign(%s)" % row["zero_sign"]))
    busy = set()
    p = os.path.join(out, "gpu-busy.txt")
    if os.path.exists(p):
        busy = {tuple(l.split()) for l in open(p) if l.strip()}
    sources = [x.split("=", 1) for x in a.cpu_sources.split(",")] if a.cpu_sources else [("", "raw")]

    def points(path):
        series = collections.OrderedDict()
        if os.path.exists(path):
            for line in open(path):
                f = line.strip().split(",")
                if len(f) == 8:
                    series.setdefault((f[0], f[1], f[2]), []).append(f)
        return series

    rows = []
    for t in TYPES:
        cpu, taken = [], set()
        for name, d in sources:
            for key, pts in points(os.path.join(out, d, "%s_cpu.csv" % t)).items():
                if key not in taken:
                    taken.add(key)
                    cpu += [(f, name) for f in pts]
        gpu = [(f, "") for pts in points(os.path.join(out, "raw", "%s_gpu.csv" % t)).values() for f in pts]
        for f, host in gpu + cpu:
            label, kernel, form, n, m, sec, gf, runs = f
            key = (label, kernel, form)
            rows.append([label, kernel, form, n, m, sec, gf, runs, " ".join(checks.get(key, [])) if sec else "FAILED",
                         "1" if m == "gpu" and key in busy else "", host])
    with open(os.path.join(out, "kernels.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["type", "kernel", "form", "n", "method", "seconds", "gflops", "runs", "bitcheck", "gpu_busy", "cpu_host"])
        w.writerows(rows)
    print("%d rows -> %s" % (len(rows), os.path.join(out, "kernels.csv")))


# ---- tables for the README ----

def cmd_table(a):
    """markdown: per type, GFLOPS of GPU and CPU and their ratio at n0 (2048,
    MPFR 1024; potrf2 256; the largest common size below n0 otherwise)"""
    out = os.path.abspath(a.out)
    rows = list(csv.DictReader(open(os.path.join(out, "kernels.csv"))))
    data = {}
    for r in rows:
        if r["seconds"]:
            data.setdefault((r["type"], r["kernel"], r["form"], r["method"]), {})[int(r["n"])] = r
    for t in a.types:
        n0 = 1024 if is_mpfr(t) else 2048
        print("\n#### %s (n = %d unless noted)\n" % (t, n0))
        print("| kernel | form | n | GPU GFLOPS | CPU GFLOPS | GPU/CPU | GPU max GFLOPS (n) | CPU host | bit check (n = 64, 200) |")
        print("|---|---|---:|---:|---:|---:|---:|---|---|")
        for kernel, form in SERIES:
            g = data.get((t, kernel, form, "gpu"), {})
            c = data.get((t, kernel, form, "cpu"), {})
            if not g and not c:
                continue
            want = 256 if kernel == "potrf2" else n0
            common = sorted(n for n in g if n in c and n <= want)
            n = want if want in g and want in c else (common[-1] if common else want)
            gg, cc = g.get(n), c.get(n)
            ratio = float(cc["seconds"]) / float(gg["seconds"]) if gg and cc else None
            def gf(x):
                if not x:
                    return "—"
                if kernel == "laswp":
                    return "%.3g s" % float(x["seconds"])
                return "%.3g" % float(x["gflops"])
            best = max(g.items(), key=lambda kv: float(kv[1]["gflops"] or 0)) if g and kernel != "laswp" else None
            check = (gg or cc or {}).get("bitcheck", "")
            check = "identical" if check and all(x.endswith(":identical") for x in check.split()) else (check or "— (no GPU)")
            print("| %s | %s | %d | %s | %s | %s | %s | %s | %s |" % (kernel, form, n, gf(gg), gf(cc), "%.3g" % ratio if ratio else "—",
                                                            "%.3g (%d)" % (float(best[1]["gflops"]), best[0]) if best else "—",
                                                            (cc or {}).get("cpu_host", "") or "—", check))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("command", choices=["probe", "plan", "run", "check", "merge", "table"])
    p.add_argument("--bindir", default=HERE)
    p.add_argument("--out", required=True)
    p.add_argument("--types", default=",".join(ALL_TYPES))
    p.add_argument("--methods", default="gpu,cpu")
    p.add_argument("--threads", type=int, default=32)
    p.add_argument("--pin", action="store_true", help="run the CPU programs on cpus 0..threads-1 (taskset)")
    p.add_argument("--plan")
    p.add_argument("--only", choices=["gpu", "cpu"])
    p.add_argument("--cpu-sources", default="", help="merge: NAME=DIR,... CPU results in order of preference (default raw)")
    p.add_argument("--hours", type=float, default=0, help="run: start no series after this many hours")
    p.add_argument("--top", action="append", help="plan: TYPE:METHOD:KERNEL[:FORM]=N caps the sizes (* matches all)")
    p.add_argument("--note", default="")
    p.add_argument("--max-point", type=float, default=0, help="plan: drop GPU sizes predicted to run longer (s)")
    p.add_argument("--max-point-fine", type=float, default=0, help="plan: the same for the FINE series")
    p.add_argument("--coarse", action="store_true", help="plan: coarser GPU grid for the series not in FINE")
    p.add_argument("--verbose", action="store_true")
    a = p.parse_args()
    a.types = a.types.split(",")
    a.methods = a.methods.split(",")
    {"probe": cmd_probe, "plan": cmd_plan, "run": cmd_run, "check": cmd_check, "merge": cmd_merge, "table": cmd_table}[a.command](a)


if __name__ == "__main__":
    main()
