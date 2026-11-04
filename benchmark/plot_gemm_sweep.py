#!/usr/bin/env python3
"""Overlay GPU kernel-only Rgemm sweeps (CSV n,seconds,gflops), in the style of
benchmark/plot_sweep.py (one line per series, powers of two marked).

    plot_gemm_sweep.py OUT.png --title T --subtitle S LABEL=FILE.csv [LABEL=FILE.csv ...]

Used for benchmark/results/2026/<host>/gemm-sweep/: one graph per type with
NN/NT/TN/TT, and one graph of NN for every type (--logy recommended there).
Needs matplotlib and the Noto Sans CJK JP font.
"""

import argparse
import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

COLORS = ["#2a78d6", "#eb6834", "#1b9e77", "#8e44ad", "#c0392b", "#7f8c8d"]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("png")
    p.add_argument("series", nargs="+", help="LABEL=FILE.csv")
    p.add_argument("--title", default="")
    p.add_argument("--subtitle", default="")
    p.add_argument("--xlabel", default="n (m = n = k)")
    p.add_argument("--ylabel", default="GFLOPS (2n³ / 時間)")
    p.add_argument("--logy", action="store_true")
    p.add_argument("--annotate", default="", help="label of the series whose powers of two are annotated (default: first)")
    a = p.parse_args()

    page, text, sub, grid = "#f9f9f7", "#0b0b0b", "#52514e", "#e1e0d9"
    plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 11,
                         "text.color": text, "axes.labelcolor": sub, "xtick.color": sub, "ytick.color": sub})
    fig, ax = plt.subplots(figsize=(13.8, 7.2), dpi=100)
    fig.patch.set_facecolor(page)
    ax.set_facecolor(page)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(grid)
    ax.grid(True, color=grid, linewidth=1, which="both" if a.logy else "major")
    ax.set_axisbelow(True)

    nmax = 1
    first = None
    for i, spec in enumerate(a.series):
        label, path = spec.split("=", 1)
        rows = sorted((int(r["n"]), float(r["gflops"])) for r in csv.DictReader(open(path)))
        rows = [r for r in rows if r[1] > 0]
        if not rows:
            continue
        nmax = max(nmax, rows[-1][0])
        color = COLORS[i % len(COLORS)]
        ax.plot([r[0] for r in rows], [r[1] for r in rows], color=color, linewidth=2.2,
                label="%s (%d点、最大 %.3g)" % (label, len(rows), max(r[1] for r in rows)))
        pow2 = [r for r in rows if r[0] & (r[0] - 1) == 0 and r[0] >= 64]
        ax.scatter([r[0] for r in pow2], [r[1] for r in pow2], color=color, s=36, zorder=3)
        if (a.annotate and label == a.annotate) or (not a.annotate and first is None):
            first = label
            for n, g in pow2:
                if n >= 512:
                    ax.annotate("%d: %.3g" % (n, g), (n, g), xytext=(6, -15), textcoords="offset points", fontsize=10)

    ax.set_xlim(0, nmax * 1.03)
    if a.logy:
        ax.set_yscale("log")
    else:
        ax.set_ylim(0, None)
    ax.set_xlabel(a.xlabel)
    ax.set_ylabel(a.ylabel)
    leg = ax.legend(loc="lower right", frameon=False)
    for t in leg.get_texts():
        t.set_color(text)
    fig.suptitle(a.title, x=0.02, ha="left", fontsize=17, fontweight="bold", color=text)
    if a.subtitle:
        fig.text(0.05, 0.915, a.subtitle, fontsize=12, color=sub)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(a.png, facecolor=page)


if __name__ == "__main__":
    main()
