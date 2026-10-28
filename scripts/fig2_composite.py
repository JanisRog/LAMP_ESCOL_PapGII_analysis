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
CONDS = ["10^6", "10^5", "10^4", "10^3", "10^2", "NTC"]
LABELS = ["10$^6$", "10$^5$", "10$^4$", "10$^3$", "10$^2$", "NTC"]
RAMP = ["#08306b", "#2171b5", "#4292c6", "#6baed6", "#9ecae1"]        # one hue, dense -> dilute
COLORS = RAMP + ["#8a8984"]


def lod(target):
    d = pd.concat([pd.read_csv(DATA / f"{'uida' if target == 'uidA' else 'papg'}_lod_{i}.csv") for i in (1, 2, 3)], ignore_index=True)
    d.columns = ["sample", "time_min", "fluorescence"]
    d["condition"] = d["sample"].str.replace(r"_\d+$", "", regex=True)
    d["rel"] = d["fluorescence"] / d.groupby("sample")["fluorescence"].transform("first")
    g = d.groupby(["condition", "time_min"]).rel
    return g.mean().reset_index(name="mean").merge(g.std().reset_index(name="sd"))


def lod_panel(ax, target, letter):
    m = lod(target)
    for c, lab, col in zip(CONDS, LABELS, COLORS):
        s = m[m.condition == c]
        ax.plot(s.time_min, s["mean"], color=col, lw=1.4, ls="--" if c == "NTC" else "-", label=lab)
        ax.fill_between(s.time_min, s["mean"] - s.sd, s["mean"] + s.sd, color=col, alpha=0.13, lw=0)
    ax.set_xlim(0, 60); ax.set_ylim(0, None)
    ax.set_xlabel("Time (min)", fontsize=8.5, color=INK2); ax.set_ylabel("Relative fluorescence (fold change)", fontsize=8.5, color=INK2)
    ax.tick_params(labelsize=8, colors=INK2, length=2)
    ax.grid(color="#e1e0d9", lw=0.5); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color("#c3c2b7")
    ax.set_title(target, loc="left", fontsize=12, fontstyle="italic", color=INK, pad=8)
    ax.text(-0.2 if letter == "B" else -0.23, 1.2, letter, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom", color=INK)
    return ax


def main():
    dims = (6.5, 22.0, 14.0, 8.0, 0.5, 24.0, 26.0)
    ROW, TRK, GENE, ZOOM, TIGHT, BETWEEN, CB = dims
    hA = (GENE + ZOOM + TRK + TIGHT) * 2 + ROW * 12 + BETWEEN + CB
    LOD, GAP = 58.0, 20.0
    fig = plt.figure(figsize=(180 / 25.4, (hA + GAP + LOD + 30) / 25.4))
    top = fig.add_gridspec(3, 1, height_ratios=[hA, GAP, LOD], hspace=0)
    H.draw_panel(fig, top[0], dims, title_x=-0.17)
    fig.text(0.0, 0.985, "A", fontsize=14, fontweight="bold", va="top", color=INK)
    bot = top[2].subgridspec(1, 2, wspace=0.32)
    axb = lod_panel(fig.add_subplot(bot[0]), "uidA", "B")
    axc = lod_panel(fig.add_subplot(bot[1]), "papGII", "C")
    h, l = axc.get_legend_handles_labels()
    fig.legend(h, l, title="CFU/ml", title_fontsize=8.5, fontsize=8.5, frameon=False, ncol=6, loc="lower center", bbox_to_anchor=(0.575, 0.0))
    fig.subplots_adjust(left=0.17, right=0.98, top=0.96, bottom=0.085)
    draw_zoom_lines(fig)
    for ext in ("pdf", "png", "svg"):
        meta = {"pdf": {"CreationDate": None, "Creator": "fig2_composite.py"}, "svg": {"Date": None}}.get(ext)
        fig.savefig(OUT / f"fig2.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight", metadata=meta)


if __name__ == "__main__":
    main()
