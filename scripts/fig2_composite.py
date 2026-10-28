"""Figure 2: (A) primer maps with median mismatches per group and binding site (fig_primer_heatmap.py), (B, C) limit of
detection of the uidA and papGII LAMP assays (mean +- SD of three replicates; fluorescence relative to the first cycle).

Output: results/figures/fig2.{pdf,png,svg}
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

import fig_primer_heatmap as H
from fig_primer_lollipop import INK, INK2, draw_zoom_lines
from lamp_common import DATA, RESULTS

OUT = RESULTS / "figures"
CONDS = ["10^6", "10^5", "10^4", "10^3", "10^2", "NTC"]   # as in the notebook (default colour cycle, mean ± SD)


def lod(target):
    d = pd.concat([pd.read_csv(DATA / f"{'uida' if target == 'uidA' else 'papg'}_lod_{i}.csv") for i in (1, 2, 3)], ignore_index=True)
    d.columns = ["sample", "time_min", "fluorescence"]
    d["condition"] = d["sample"].str.replace(r"_\d+$", "", regex=True)
    d["rel"] = d["fluorescence"] / d.groupby("sample")["fluorescence"].transform("first")
    g = d.groupby(["condition", "time_min"]).rel
    return g.mean().reset_index(name="mean").merge(g.std().reset_index(name="sd"))


def lod_panel(ax, target, letter):
    m = lod(target)
    for c in CONDS:
        d = m[m.condition == c]
        ax.plot(d.time_min, d["mean"], label=c)
        ax.fill_between(d.time_min, d["mean"] - d.sd, d["mean"] + d.sd, alpha=0.2)
    ax.set_xlim(0, 60); ax.set_ylim(0, None); ax.grid(True)
    ax.set_xlabel("Time [min]", fontsize=8.5); ax.set_ylabel("Relative Fluorescence (fold-change)", fontsize=8.5)
    ax.tick_params(labelsize=8)
    ax.set_title(f"LAMP limit of detection\nfor $\\it{{{target}}}$ (mean ± SD)", fontsize=8.5, pad=6)
    ax.legend(title="CFU/ml", title_fontsize=7.5, fontsize=7.5, loc="upper left")
    ax.text(-0.27, 1.22, letter, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom", color=INK)
    return ax


def main():
    dims = (6.5, 22.0, 14.0, 8.0, 0.5, 30.0, 22.0)
    ROW, TRK, GENE, ZOOM, TIGHT, BETWEEN, CB = dims
    h_u = GENE + ZOOM + TRK + TIGHT + ROW * 8
    h_p = GENE + ZOOM + TRK + TIGHT + ROW * 4
    hA = h_u + BETWEEN + h_p + CB
    W = 215.0
    fig = plt.figure(figsize=(W / 25.4, (hA + 14) / 25.4))
    top = fig.add_gridspec(1, 2, width_ratios=[118, 62], wspace=0.2, left=38 / W, right=0.985, top=0.985, bottom=0.04)
    H.SUB_SEP, H.LETTER = "\n", "A"
    H.draw_panel(fig, top[0], dims, title_x=-0.3)
    right = top[1].subgridspec(4, 1, height_ratios=[h_u, BETWEEN, h_p, CB], hspace=0)
    for spec, target, letter, h in ((right[0], "uidA", "B", h_u), (right[2], "papGII", "C", h_p)):
        inner = spec.subgridspec(2, 1, height_ratios=[h - (TRK + TIGHT + (ROW * (8 if target == "uidA" else 4))) - 0.0, TRK + TIGHT + ROW * (8 if target == "uidA" else 4)], hspace=0)
        lod_panel(fig.add_subplot(inner[1]), target, letter)
    draw_zoom_lines(fig)
    for ext in ("pdf", "png", "svg"):
        meta = {"pdf": {"CreationDate": None, "Creator": "fig2_composite.py"}, "svg": {"Date": None}}.get(ext)
        fig.savefig(OUT / f"fig2.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight", metadata=meta)


if __name__ == "__main__":
    main()
