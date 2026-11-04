#!/usr/bin/env python3
"""Graphs of the kernel benchmark (kernel_bench.py merge -> kernels.csv).

    plot_kernel_bench.py kernels.csv OUTDIR [--gpu-name "NVIDIA A100 80GB PCIe"]

1. OUTDIR/<type>/<kernel>_<form>.png: GFLOPS against n, drawn by
   plot_sweep.draw() as benchmark/results/2026/.../rgemm/mpfr512_sweep10k.png:
   the GPU points joined (solid), the powers of two marked and annotated, the
   CPU powers of two dashed (laswp: G elements swapped per second).
2. OUTDIR/ratio_<type>.png: GPU time / CPU time inverted (speed ratio), log
   scale, at the powers of two measured on both, one colour per kernel.
3. OUTDIR/summary.png: GFLOPS at n = 2048 (MPFR: 1024), CPU and GPU side by
   side, every type and kernel.
The PNG files are stored with a 256-colour palette (PIL), a third of the size.
Needs matplotlib (PYTHONPATH=/home/docker/plotlib on the benchmark host) and
the Noto Sans CJK JP font.
"""

import argparse
import collections
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plot_sweep  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

TYPE_NAMES = {
    "dd": "dd (double-double)", "td": "td (triple-double)", "qd": "qd (quad-double)",
    "mpfr256": "MPFR 256ビット (cu_freal)", "mpfr512": "MPFR 512ビット (cu_freal)",
    "mpfr1024": "MPFR 1024ビット (cu_freal)", "mpfr200rt": "MPFR 200ビット (cu_mpfr、実行時精度)",
}
ROUTINE = {"gemm": "Rgemm", "gemm_update": "Rgemm", "syrk": "Rsyrk", "trsm": "Rtrsm", "getrf2": "Rgetrf2",
           "potrf2": "Rpotrf2", "laswp": "Rlaswp", "trtri": "Rtrtri", "lauum": "Rlauum",
           "cgemm": "Cgemm", "cherk": "Cherk", "ctrsm": "Ctrsm"}
FLOPS_TEXT = {"gemm": "2n³", "gemm_update": "2·64·n²", "syrk": "k·n(n+1)", "trsm": "n³ (Left: n·m², Right: m·n²)",
              "getrf2": "LU m×64 (LAWN 41)", "potrf2": "n³/3+n²/2+n/6", "trtri": "n³/3+2n/3", "lauum": "n³/3+n²/2+n/6",
              "cgemm": "8n³", "cherk": "4·k·n(n+1)", "ctrsm": "4·(TRSM)", "laswp": "64·n 交換 / 時間"}


def compact(png):
    """256-colour palette PNG (line graphs: no visible change, a third of the size)"""
    try:
        from PIL import Image
    except ImportError:
        return
    im = Image.open(png).convert("RGB")
    im.quantize(colors=256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE).save(png, optimize=True)


def describe(kernel, form):
    """routine call of the series, in Japanese"""
    base = form.split("_")[0]
    p64 = form.endswith("_p64")
    k64 = form.endswith("_k64")
    if kernel in ("gemm", "cgemm"):
        return '%s("%s","%s",n,n,n)' % (ROUTINE[kernel], base[0], base[1]), "n (m = n = k)"
    if kernel == "gemm_update":
        return 'Rgemm("N","N",n,n,64)(分解の更新)', "n (m = n、k = 64)"
    if kernel in ("syrk", "cherk"):
        return '%s("%s","%s",n,%s)' % (ROUTINE[kernel], base[0], base[1], "64" if k64 else "n"), "n (k = %s)" % ("64" if k64 else "n")
    if kernel in ("trsm", "ctrsm"):
        side = base[0]
        call = '%s("%s","%s","%s","%s",' % (ROUTINE[kernel], base[0], base[1], base[2], base[3])
        if p64:
            call += ("64,n" if side == "L" else "n,64") + ")、三角 64×64"
            return call, "n (B は %s)" % ("64×n" if side == "L" else "n×64")
        return call + "n,n)", "n (三角 n×n、B n×n)"
    if kernel == "getrf2":
        return "LU パネル Rgetrf2(n,64)", "n (パネルの行数 m = n、64 列)"
    if kernel == "potrf2":
        return "Cholesky 対角ブロック Rpotrf2(\"L\",n)", "n"
    if kernel == "laswp":
        return "Rlaswp(n,A,n,1,64,ipiv,1)(64 行の入れ替え)", "n (n×n)"
    if kernel == "trtri":
        return 'Rtrtri("U","N",n)', "n"
    if kernel == "lauum":
        return 'Rlauum("L",n)', "n"
    return kernel + " " + form, "n"


def load(path):
    data = collections.defaultdict(dict)  # (type, kernel, form, method) -> n -> row
    for r in csv.DictReader(open(path)):
        if not r["seconds"]:
            continue
        data[(r["type"], r["kernel"], r["form"], r["method"])][int(r["n"])] = r
    return data


LOADED_NOTE = "負荷のあるホストで測った参考値、静かなホストより約1.5倍遅い"


def cpu_label(cpu):
    """legend of the CPU line: quiet host, or the loaded-host reference values"""
    hosts = {r.get("cpu_host", "") for r in cpu.values()}
    if "loaded" in hosts:
        return "CPU 32コア (%s)" % LOADED_NOTE
    return "CPU 32コア (2 の冪、静かなホスト)" if "quiet" in hosts else "CPU 32コア (2 の冪)"


def rate(kernel, n, row):
    if kernel == "laswp":
        return 64.0 * n / float(row["seconds"]) * 1e-9
    return float(row["gflops"])


def series_graphs(data, out, gpu_name, meta):
    made = []
    keys = sorted({(t, k, f) for (t, k, f, m) in data})
    for t, kernel, form in keys:
        gpu = data.get((t, kernel, form, "gpu"), {})
        cpu = data.get((t, kernel, form, "cpu"), {})
        if kernel == "potrf2" or not (gpu or cpu):
            continue  # two sizes only: in the summary
        call, xlabel = describe(kernel, form)
        rows = sorted((n, rate(kernel, n, r)) for n, r in gpu.items())
        cpu_pts = sorted((n, rate(kernel, n, r)) for n, r in cpu.items() if n & (n - 1) == 0)
        cpu_opt = ",".join("%d:%.6g" % p for p in cpu_pts) if len(cpu_pts) > 1 else ("%.6g" % cpu_pts[0][1] if cpu_pts else "")
        runs = sorted({r["runs"] for r in gpu.values()})
        minutes = sum(float(r["seconds"]) * (int(r["runs"]) + 1) for r in gpu.values()) / 60
        if kernel == "laswp":
            ylabel = "G要素/秒 (64·n 交換 / 時間)"
        else:
            ylabel = "GFLOPS (%s / 時間)" % FLOPS_TEXT[kernel]
        gpu_label = "GPU" if not t.startswith("mpfr") else "GPU (%s)" % ("cu_mpfr" if t == "mpfr200rt" else "cu_freal")
        how = "CUPTI、転送・詰め直しを含まない" if t.startswith("mpfr") else "GPU 上のデータ、CUDA イベント、転送を含まない"
        subtitle = ("%s\nGPU: %s、予熱1回+最良3回(10秒超は1回)、%d点、計測約%.1f分。\nCPU: libmplapack_%s_opt 32スレッド、2 の冪の点%s。"
                    % (call, how, len(rows), minutes, "mpfr" if t.startswith("mpfr") else t,
                       "、" + LOADED_NOTE if "負荷" in cpu_label(cpu) else "、静かなホスト"))
        if meta.get((t, kernel, form)):
            subtitle += " " + meta[(t, kernel, form)]
        title = "%s %s、GPU と CPU (%s)" % (TYPE_NAMES.get(t, t), ROUTINE[kernel] + ("" if kernel in ("getrf2", "laswp", "trtri", "lauum") else " " + form), gpu_name)
        d = os.path.join(out, t)
        os.makedirs(d, exist_ok=True)
        png = os.path.join(d, "%s_%s.png" % (kernel, form))
        if not rows:  # CPU only (MPFR laswp, trtri, lauum: no GPU version)
            rows = cpu_pts
            cpu_opt = ""
            title = "%s %s、CPU のみ (GPU 版なし)" % (TYPE_NAMES.get(t, t), ROUTINE[kernel])
            subtitle = "%s\nCPU: libmplapack_mpfr_opt 32スレッド、2 の冪の点、3回中の最良値(10秒を超える点は1回)\n%s。" % (
                call, LOADED_NOTE if "負荷" in cpu_label(cpu) else "静かなホスト")
        if not gpu:
            gpu_label = cpu_label(cpu).replace("CPU 32コア", "CPU 32コア、GPU 版なし")
        plot_sweep.draw(rows, png, title=title, subtitle=subtitle, xlabel=xlabel, ylabel=ylabel,
                        series=gpu_label, cpu=cpu_opt, cpu_label=cpu_label(cpu),
                        annotate_min=256, ymax=None, subtitle_lines=3)
        compact(png)
        made.append(png)
    return made


COLORS = {"gemm": "#1f5fbf", "cgemm": "#1f5fbf", "gemm_update": "#4aa3df", "syrk": "#2e9e4f", "cherk": "#2e9e4f",
          "trsm": "#e0731c", "ctrsm": "#e0731c", "getrf2": "#8b3fbf", "laswp": "#7a7a7a", "trtri": "#c23b5a", "lauum": "#b08a1a"}
DASHES = ["-", "--", ":", "-."]


def ratio_graphs(data, out, gpu_name):
    made = []
    types = sorted({k[0] for k in data}, key=lambda t: list(TYPE_NAMES).index(t) if t in TYPE_NAMES else 99)
    for t in types:
        plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 10})
        fig, axes = plt.subplots(1, 2, figsize=(16, 9.0), dpi=100, sharey=True)
        for ax, cplx in zip(axes, (False, True)):
            used = collections.Counter()
            for (tt, kernel, form, m), gpu in sorted(data.items()):
                if tt != t or m != "gpu" or kernel == "potrf2" or kernel.startswith("c") != cplx:
                    continue
                cpu = data.get((t, kernel, form, "cpu"), {})
                pts = sorted((n, float(cpu[n]["seconds"]) / float(gpu[n]["seconds"])) for n in gpu if n in cpu and n & (n - 1) == 0)
                if not pts:
                    continue
                style = DASHES[used[kernel] % len(DASHES)]
                used[kernel] += 1
                loaded = any(r.get("cpu_host") == "loaded" for r in cpu.values())
                ax.plot([p[0] for p in pts], [p[1] for p in pts], style, marker="^" if loaded else "o", ms=4, color=COLORS.get(kernel, "#444"),
                        label="%s%s" % ("Rgemm 更新 k=64" if kernel == "gemm_update" else ROUTINE[kernel] + " " + form, " (CPU 負荷下)" if loaded else ""))
            ax.axhline(1.0, color="#0b0b0b", linewidth=1.2)
            ax.set_xscale("log", base=2)
            ax.set_yscale("log")
            ax.grid(True, which="both", color="#e1e0d9")
            ax.set_axisbelow(True)
            ax.set_xlabel("n")
            ax.set_title("複素" if cplx else "実数", fontsize=13)
            if ax.get_legend_handles_labels()[0]:
                ax.legend(fontsize=8, ncol=3, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.1))
        axes[0].set_ylabel("GPU / CPU の速度比 (CPU 時間 / GPU 時間、1 より上で GPU が速い)")
        fig.suptitle("%s: GPU カーネル (%s) と CPU 32コアの速度比、2 の冪の点" % (TYPE_NAMES.get(t, t), gpu_name), x=0.02, y=0.99, ha="left", va="top", fontsize=15, fontweight="bold")
        fig.text(0.02, 0.945, "「CPU 負荷下」(▲) の系列の CPU は%s(その分、比は大きめに出る)。GPU は転送を含まない。" % LOADED_NOTE, fontsize=10, color="#52514e")
        fig.tight_layout(rect=(0, 0, 1, 0.93))
        png = os.path.join(out, "ratio_%s.png" % t)
        fig.savefig(png)
        plt.close(fig)
        compact(png)
        made.append(png)
    return made


SUMMARY = [("gemm", "NN"), ("gemm_update", "NN"), ("syrk", "LN"), ("syrk", "LN_k64"), ("trsm", "LLNU"), ("trsm", "RLTN"),
           ("trsm", "RLTN_p64"), ("getrf2", "P64"), ("potrf2", "L"), ("trtri", "UN"), ("lauum", "L"),
           ("cgemm", "NN"), ("cherk", "LN"), ("ctrsm", "LLNU")]


def summary_graph(data, out, gpu_name):
    types = [t for t in TYPE_NAMES if any(k[0] == t for k in data)]
    plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 10})
    fig, axes = plt.subplots(len(types), 1, figsize=(16, 3.0 * len(types)), dpi=100)
    for ax, t in zip(axes if len(types) > 1 else [axes], types):
        n0 = 1024 if t.startswith("mpfr") else 2048
        labels, cpu_v, gpu_v, cpu_host = [], [], [], []
        for kernel, form in SUMMARY:
            want = 256 if kernel == "potrf2" else n0
            gs, cs = data.get((t, kernel, form, "gpu"), {}), data.get((t, kernel, form, "cpu"), {})
            # the largest size <= want measured on both (GPU sizes capped by the time budget)
            common = sorted(m for m in gs if m in cs and m <= want)
            n = want if (want in cs and (want in gs or not gs)) or not common else common[-1]
            g, c = gs.get(n), cs.get(n)
            name = "Rgemm 更新 k=64" if kernel == "gemm_update" else "%s %s" % (ROUTINE[kernel], form)
            labels.append("%s\n(n=%d)" % (name, n))
            gpu_v.append(float(g["gflops"]) if g else float("nan"))
            cpu_v.append(float(c["gflops"]) if c else float("nan"))
            cpu_host.append(c.get("cpu_host", "") if c else "")
        x = range(len(labels))
        bars = ax.bar([i - 0.2 for i in x], cpu_v, width=0.4, color="#9a9890", label="CPU 32コア")
        for b, h in zip(bars, cpu_host):
            if h == "loaded":
                b.set_hatch("//")
        ax.bar([i + 0.2 for i in x], gpu_v, width=0.4, color="#2a78d6", label="GPU カーネル")
        for i, (c, g) in enumerate(zip(cpu_v, gpu_v)):
            for v, dx in ((c, -0.2), (g, 0.2)):
                if not math.isnan(v):
                    ax.annotate("%.3g" % v, (i + dx, v), ha="center", va="bottom", fontsize=7)
        ax.set_yscale("log")
        ax.set_axisbelow(True)
        ax.set_xticks(list(x))
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_ylabel("GFLOPS")
        ax.set_title(TYPE_NAMES.get(t, t), loc="left", fontsize=12)
        ax.grid(True, axis="y", which="both", color="#e1e0d9")
        ax.legend(fontsize=8, frameon=False, loc="upper right")
    fig.suptitle("カーネル単体の GFLOPS、n = 2048 (MPFR は 1024、Rpotrf2 は 256)、CPU と GPU (%s)。GPU は転送を含まない" % gpu_name,
                 x=0.02, y=1 - 0.2 / fig.get_figheight(), ha="left", va="top", fontsize=15, fontweight="bold")
    fig.text(0.02, 1 - 0.75 / fig.get_figheight(), "斜線の CPU 棒は%s。GPU の棒が無いものは GPU 版が無い (MPFR の Rtrtri、Rlauum) か、"
             "その大きさの GPU を測っていない (時間の上限、括弧内の n は両方を測った最大の大きさ)。" % LOADED_NOTE, fontsize=10, color="#52514e")
    fig.tight_layout(rect=(0, 0, 1, 1 - 1.0 / fig.get_figheight()))
    png = os.path.join(out, "summary.png")
    fig.savefig(png)
    plt.close(fig)
    compact(png)
    return [png]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv")
    p.add_argument("out")
    p.add_argument("--gpu-name", default="NVIDIA A100")
    p.add_argument("--only", choices=["series", "ratio", "summary"])
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    data = load(a.csv)
    meta = {}
    made = []
    if a.only in (None, "series"):
        made += series_graphs(data, a.out, a.gpu_name, meta)
    if a.only in (None, "ratio"):
        made += ratio_graphs(data, a.out, a.gpu_name)
    if a.only in (None, "summary"):
        made += summary_graph(data, a.out, a.gpu_name)
    print("%d graphs in %s" % (len(made), a.out))


if __name__ == "__main__":
    main()
