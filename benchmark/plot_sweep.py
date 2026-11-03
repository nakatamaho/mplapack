#!/usr/bin/env python3
"""Plot a GFLOPS-versus-n sweep (CSV with columns n,seconds,gflops).

    plot_sweep.py SWEEP.csv OUT.png [options]

Draws the GPU sweep as a line, marks the powers of two with dots and their
values, and draws the CPU reference as a dashed line (one value: horizontal;
several "n:gflops" values: a dashed polyline through them).  Used for
benchmark/results/rgemm-a100-2026-10/mpfr512_sweep10k.png and as the model of
the kernel benchmark graphs.  Needs matplotlib and a Japanese font (Noto Sans
CJK JP) for Japanese labels.
"""

import argparse
import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv")
    p.add_argument("png")
    p.add_argument("--title", default="MPFR 512ビット Rgemm、GPU 固定精度 (cu_freal) の速度とサイズ (NVIDIA A100)")
    p.add_argument("--subtitle", default="")
    p.add_argument("--xlabel", default="n (m = n = k)")
    p.add_argument("--ylabel", default="GFLOPS (2n³ / 時間)")
    p.add_argument("--series", default="GPU", help="legend label of the sweep line")
    p.add_argument("--cpu", default="", help='CPU reference: "GFLOPS" or "n:GFLOPS,n:GFLOPS,..."')
    p.add_argument("--cpu-label", default="CPU 32コア")
    p.add_argument("--min-n", type=int, default=1)
    p.add_argument("--dark", action="store_true")
    a = p.parse_args()

    rows = [(int(r["n"]), float(r["gflops"])) for r in csv.DictReader(open(a.csv))]
    rows = sorted(r for r in rows if r[0] >= a.min_n)
    pow2 = [r for r in rows if r[0] & (r[0] - 1) == 0]

    if a.dark:
        page, text, sub, grid, s1, s2 = "#0d0d0d", "#f0efec", "#c3c2b7", "#2c2c2a", "#3987e5", "#d95926"
    else:
        page, text, sub, grid, s1, s2 = "#f9f9f7", "#0b0b0b", "#52514e", "#e1e0d9", "#2a78d6", "#eb6834"
    plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 11,
                         "text.color": text, "axes.labelcolor": sub, "xtick.color": sub, "ytick.color": sub})
    fig, ax = plt.subplots(figsize=(13.8, 7.2), dpi=100)
    fig.patch.set_facecolor(page)
    ax.set_facecolor(page)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(grid)
    ax.grid(True, color=grid, linewidth=1)
    ax.set_axisbelow(True)

    ax.plot([r[0] for r in rows], [r[1] for r in rows], color=s1, linewidth=2.6,
            label="%s (%d点)" % (a.series, len(rows)))
    ax.scatter([r[0] for r in pow2], [r[1] for r in pow2], color=s2, s=60, zorder=3,
               label="n = 2 の冪 (%d点)" % len(pow2))
    for n, g in pow2:
        if n >= 1024:
            ax.annotate("%d: %.2f" % (n, g), (n, g), xytext=(8, -16), textcoords="offset points", fontsize=11)
    if rows:
        n, g = rows[-1]
        ax.annotate("%d: %.2f" % (n, g), (n, g), xytext=(-60, 12), textcoords="offset points", fontsize=11)

    ax.set_xlim(0, rows[-1][0] * 1.03 if rows else 1)
    if a.cpu:
        if ":" in a.cpu:
            pts = sorted((int(x.split(":")[0]), float(x.split(":")[1])) for x in a.cpu.split(","))
            ax.plot([q[0] for q in pts], [q[1] for q in pts], "--", color=sub, linewidth=1.4, label=a.cpu_label)
        else:
            ax.axhline(float(a.cpu), linestyle="--", color=sub, linewidth=1.4)
            ax.annotate("%s %.3g GFLOPS" % (a.cpu_label, float(a.cpu)), (ax.get_xlim()[1], float(a.cpu)),
                        xytext=(-4, 4), textcoords="offset points", ha="right", fontsize=11, color=sub)

    ax.set_xlabel(a.xlabel)
    ax.set_ylabel(a.ylabel)
    ax.set_ylim(0, None)
    leg = ax.legend(loc="upper left", frameon=False)
    for t in leg.get_texts():
        t.set_color(text)
    fig.suptitle(a.title, x=0.02, ha="left", fontsize=17, fontweight="bold", color=text)
    if a.subtitle:
        fig.text(0.05, 0.915, a.subtitle, fontsize=12, color=sub)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(a.png, facecolor=page)


if __name__ == "__main__":
    main()
