"""Figure 2: (A) primer maps with median mismatches per group and binding site (fig_primer_heatmap.py), (B, C) limit of
detection of the uidA and papGII LAMP assays (median and interquartile range of three replicates; fluorescence relative
to the first cycle). Layout as in the original figure (A left, divider, B above C); panel A replaces the colorimetric panel.

Output: results/figures/fig2.{pdf,png,svg}
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

import fig_primer_heatmap as H
from fig_primer_lollipop import draw_zoom_lines
from lamp_common import DATA, RESULTS

OUT = RESULTS / "figures"
CONDS = ["10^6", "10^5", "10^4", "10^3", "10^2", "NTC"]   # default colour cycle, as in the notebook


def lod(target):
    d = pd.concat([pd.read_csv(DATA / f"{'uida' if target == 'uidA' else 'papg'}_lod_{i}.csv") for i in (1, 2, 3)], ignore_index=True)
    d.columns = ["sample", "time_min", "fluorescence"]
    d["condition"] = d["sample"].str.replace(r"_\d+$", "", regex=True)
    d["rel"] = d["fluorescence"] / d.groupby("sample")["fluorescence"].transform("first")
    g = d.groupby(["condition", "time_min"]).rel
    return g.median().reset_index(name="median").merge(g.quantile(0.25).reset_index(name="q25")).merge(g.quantile(0.75).reset_index(name="q75"))


def lod_panel(ax, target, letter):
    m = lod(target)
    for c in CONDS:
        d = m[m.condition == c]
        ax.plot(d.time_min, d["median"], label=c, lw=1.0)
        ax.fill_between(d.time_min, d.q25, d.q75, alpha=0.2)
    ax.set_xlim(0, 60); ax.set_ylim(0, None); ax.grid(True)
    ax.set_xlabel("Time [min]", fontsize=7); ax.set_ylabel("Relative Fluorescence (fold-change)", fontsize=7)
    ax.tick_params(labelsize=6.5)
    ax.set_title(f"LAMP limit of detection for $\\it{{{target}}}$ (median ± IQR)", fontsize=7.5)
    ax.legend(title="CFU/ml", title_fontsize=6.5, fontsize=6.5, loc="upper left")
    ax.text(-0.1, 1.12, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
    return ax


def main():
    dims = (4.6, 17.0, 12.0, 6.0, 0.5, 20.0, 20.0)
    ROW, TRK, GENE, ZOOM, TIGHT, BETWEEN, CB = dims
    hA = (GENE + ZOOM + TRK + TIGHT) * 2 + ROW * 12 + BETWEEN + CB
    W = 295.0
    fig = plt.figure(figsize=(W / 25.4, (hA + 10) / 25.4))
    top = fig.add_gridspec(1, 2, width_ratios=[100, 118], wspace=0.3, left=40 / W, right=0.99, top=0.975, bottom=0.03)
    H.SUB_SEP, H.LETTER = None, "(A)"
    H.draw_panel(fig, top[0], dims, title_x=-0.3)
    right = top[1].subgridspec(2, 1, hspace=0.42)
    axb = lod_panel(fig.add_subplot(right[0]), "uidA", "(B)")
    axc = lod_panel(fig.add_subplot(right[1]), "papGII", "(C)")
    draw_zoom_lines(fig)
    x = axb.get_position().x0 - 20 / W                      # divider between A and B/C
    fig.add_artist(Line2D([x, x], [0.0, 1.0], transform=fig.transFigure, color="black", lw=1.0))
    for ext in ("pdf", "png", "svg"):
        meta = {"pdf": {"CreationDate": None, "Creator": "fig2_composite.py"}, "svg": {"Date": None}}.get(ext)
        fig.savefig(OUT / f"fig2.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight", metadata=meta)


if __name__ == "__main__":
    main()
