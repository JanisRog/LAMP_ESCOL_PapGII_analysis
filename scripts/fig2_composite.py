"""Figure 2: (A) primer maps with median mismatches per group and binding site (fig_primer_heatmap.py), (B, C) limit of
detection of the uidA and papGII LAMP assays (median and interquartile range of three replicates; fluorescence relative
to the first cycle). Layout as in the original figure (A left, divider, B above C); panel A replaces the colorimetric panel.

Output: results/figures/fig2.{pdf,png,svg}
"""
import io

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

import fig_primer_heatmap as H
import fig_primer_map as FM
from fig_primer_lollipop import draw_zoom_lines
from lamp_common import DATA, RESULTS

OUT = RESULTS / "figures"
CONDS = ["10^6", "10^5", "10^4", "10^3", "10^2", "NTC"]   # default colour cycle, as in the notebook


def original_lod_plot(target):
    """The notebook's median + IQR plot (notebooks/LAMP_analysis.ipynb, cells 1 and 3), code unchanged apart from
    file paths, drawn with matplotlib's default style. Returns the rendered image (RGBA array)."""
    with mpl.rc_context(mpl.rcParamsDefault):
        mpl.rcParams["figure.dpi"] = 300
        files = [DATA / f"{'uida' if target == 'uidA' else 'papg'}_lod_{i}.csv" for i in (1, 2, 3)]
        data = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
        data.columns = ["sample", "time_min", "fluorescence"]
        data["replicate"] = data["sample"].str.extract(r"_(\d+)$")
        data["condition"] = data["sample"].str.replace(r"_\d+$", "", regex=True)
        baseline = data.groupby("sample")["fluorescence"].transform("first")
        data["fluorescence_norm"] = data["fluorescence"] / baseline
        grouped = data.groupby(["condition", "time_min"])["fluorescence_norm"]
        q25 = grouped.quantile(0.25).reset_index(name="q25")
        q75 = grouped.quantile(0.75).reset_index(name="q75")
        median = grouped.median().reset_index(name="median")
        order = ["10^6", "10^5", "10^4", "10^3", "10^2", "NTC"]
        f = plt.figure(figsize=(10, 6))
        for cond in order:
            m = median[median["condition"] == cond]
            q1 = q25[q25["condition"] == cond]
            q3 = q75[q75["condition"] == cond]
            plt.plot(m["time_min"], m["median"], label=cond)
            plt.fill_between(m["time_min"], q1["q25"], q3["q75"], alpha=0.2)
        plt.xlabel("Time [min]")
        plt.ylabel("Relative Fluorescence (fold-change)")
        plt.title(f"LAMP limit of detection for $\\it{{{target}}}$ (median ± IQR)")
        plt.legend(title="CFU/ml", loc="upper left")
        plt.xlim(0, 60)
        plt.ylim(0, None)
        plt.grid(True)
        plt.tight_layout()
        buf = io.BytesIO()
        f.savefig(buf, format="png", dpi=600)
        plt.close(f)
    buf.seek(0)
    return plt.imread(buf)


def lod_panel(ax, target, letter):
    ax.imshow(original_lod_plot(target), interpolation="lanczos")
    ax.set_anchor("NW"); ax.axis("off")
    if letter != "(B)":   # (B) is placed at the height of (A) in main()
        ax.text(-0.03, 1.02, letter, transform=ax.transAxes, va="bottom", **H.LETTER_KW)
    return ax


def main():
    dims = (4.6, 10.0, 12.0, 6.0, 1.5, 30.0, 0.0)
    ROW, TRK, GENE, ZOOM, TIGHT, BETWEEN, CB = dims
    hA = (GENE + ZOOM + TRK + TIGHT) * 2 + ROW * 12 + BETWEEN + CB
    W = 295.0
    fig = plt.figure(figsize=(W / 25.4, (hA + 10) / 25.4))
    top = fig.add_gridspec(1, 2, width_ratios=[92, 130], wspace=0.4, left=40 / W, right=0.99, top=0.975, bottom=0.03)
    H.SUB_SEP, H.LETTER = None, "(A)"
    FM.LABEL_INSIDE, FM.PRIMER_FS = True, 7.5
    H.draw_panel(fig, top[0], dims, title_x=-0.3)
    right = top[1].subgridspec(2, 1, hspace=0.12)
    axb = lod_panel(fig.add_subplot(right[0]), "uidA", "(B)")
    axc = lod_panel(fig.add_subplot(right[1]), "papGII", "(C)")
    draw_zoom_lines(fig)
    # vertical colour bar right of panel A, centred between the uidA and papGII heatmaps
    norm, vmax = H.SCALE
    pa = top[0].get_position(fig)
    yc, hgt = (pa.y0 + pa.y1) / 2, 0.36 * (pa.y1 - pa.y0)
    cax = fig.add_axes([pa.x1 + 5 / W, yc - hgt / 2, 3 / W, hgt])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=H.CMAP), cax=cax, ticks=range(0, vmax + 1))
    cb.outline.set_visible(False); cb.ax.tick_params(labelsize=7, colors=H.INK2, length=2)
    cb.set_label("Median mismatches to the primer", fontsize=7.5, color=H.INK2)
    a = H.LETTER_ARTIST                                     # (B) at the same height as (A)
    ya = fig.transFigure.inverted().transform(a.get_transform().transform(a.get_position()))[1]
    pb = axb.get_position()
    fig.text(pb.x0 - 0.03 * pb.width, ya, "(B)", va="bottom", **H.LETTER_KW)
    x = axb.get_position().x0 - 8 / W                      # divider between A and B/C
    fig.add_artist(Line2D([x, x], [0.0, 1.0], transform=fig.transFigure, color="black", lw=1.0))
    for ext in ("pdf", "png", "svg"):
        meta = {"pdf": {"CreationDate": None, "Creator": "fig2_composite.py"}, "svg": {"Date": None}}.get(ext)
        fig.savefig(OUT / f"fig2.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight", metadata=meta)


if __name__ == "__main__":
    main()
