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
EXTRA = (0.55, 0.40, 0.35)   # extra canvas (inch) left / bottom / top around the unchanged 10 x 6 in plot box, for the larger text
LOD_FONTS = {"axes.titlesize": 21, "axes.labelsize": 19, "xtick.labelsize": 17, "ytick.labelsize": 17,
            "axes.titlepad": 14, "axes.labelpad": 10}
HSPACE_BC = 0.12   # vertical gap between B and C
C_LIFT = -0.062      # figure fraction by which panel C (with its legend) is shifted (negative = down)
LEG_DROP = 0.016       # figure fraction by which C and the legend are lowered (more room between B and C)
LEG_ROW = 0.001     # relative height of the empty row kept free below C for the shared legend
CONDS = ["10^6", "10^5", "10^4", "10^3", "10^2", "NTC"]   # default colour cycle, as in the notebook


def original_box():
    """Position (inch: left, bottom, width, height) of the plot box in the notebook's 10 x 6 in figure (default fonts, tight_layout)."""
    with mpl.rc_context(mpl.rcParamsDefault):
        f = plt.figure(figsize=(10, 6))
        plt.plot([0, 1], [0, 1]); plt.xlabel("Time [min]"); plt.ylabel("Relative Fluorescence (fold-change)")
        plt.title("LAMP limit of detection for $\\it{uidA}$ (median ± IQR)"); plt.tight_layout()
        b = f.axes[0].get_position(); plt.close(f)
    return b.x0 * 10, b.y0 * 6, b.width * 10, b.height * 6


def original_lod_plot(target):
    """The notebook's median + IQR plot (notebooks/LAMP_analysis.ipynb, cells 1 and 3), code unchanged apart from
    file paths and without the legend (shared, below panel C), drawn with matplotlib's default style. Returns the rendered image (RGBA array)."""
    with mpl.rc_context(mpl.rcParamsDefault):
        mpl.rcParams["figure.dpi"] = 300
        mpl.rcParams.update(LOD_FONTS)   # larger title / axis / tick fonts for readability at panel size
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
        ex, eyb, eyt = EXTRA
        l0, b0, w0, h0 = original_box()
        cw, ch = 10 + ex, 6 + eyb + eyt
        f = plt.figure(figsize=(cw, ch))
        plt.axes([(l0 + ex) / cw, (b0 + eyb) / ch, w0 / cw, h0 / ch])   # plot box: same size and place as the notebook figure
        for cond in order:
            m = median[median["condition"] == cond]
            q1 = q25[q25["condition"] == cond]
            q3 = q75[q75["condition"] == cond]
            plt.plot(m["time_min"], m["median"], label=cond)
            plt.fill_between(m["time_min"], q1["q25"], q3["q75"], alpha=0.2)
        plt.xlabel("Time [min]")
        plt.ylabel("Relative Fluorescence (fold-change)")
        plt.title(f"LAMP limit of detection for $\\it{{{target}}}$ (median ± IQR)")
        # legend omitted here: one shared legend for B and C is drawn below (C) in main()
        plt.xlim(0, 60)
        plt.ylim(0, None)
        plt.grid(True)
        buf = io.BytesIO()
        f.savefig(buf, format="png", dpi=600)
        plt.close(f)
    buf.seek(0)
    return plt.imread(buf)


def lod_panel(ax, target, letter):
    img = original_lod_plot(target)
    ax.imshow(img, interpolation="lanczos", aspect="auto")
    ax.left_frac = float((img[..., :3] < 0.98).any(axis=2).any(axis=0).argmax()) / img.shape[1]   # blank margin left of the plot
    ax.axis("off")
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
    right = top[1].subgridspec(3, 1, height_ratios=[1, 1, LEG_ROW], hspace=HSPACE_BC)
    axb = lod_panel(fig.add_subplot(right[0]), "uidA", "(B)")
    axc = lod_panel(fig.add_subplot(right[1]), "papGII", "(C)")
    draw_zoom_lines(fig)
    win, hin = fig.get_figwidth(), fig.get_figheight()
    ex, eyb, eyt = EXTRA

    def base_box(spec):    # where the notebook's 10 x 6 in figure sits in this cell (fitted, anchored top-left); s = fig. inch per plot inch
        pos = spec.get_position(fig)
        sc = min(pos.width * win / 10, pos.height * hin / 6)
        return [pos.x0, pos.y1 - 6 * sc / hin, 10 * sc / win, 6 * sc / hin], sc

    def extended(box, sc):  # the enlarged image: same plot box, extra canvas left / bottom / top
        return [box[0] - ex * sc / win, box[1] - eyb * sc / hin, (10 + ex) * sc / win, (6 + eyb + eyt) * sc / hin]

    bB, sc = base_box(right[0])
    bC, _ = base_box(right[1])
    a = H.LETTER_ARTIST                                     # (B) at the same height as (A)
    ya = fig.transFigure.inverted().transform(a.get_transform().transform(a.get_position()))[1]
    nh = (6 + eyb + eyt) * sc / hin
    bB[1] = ya - 0.02 * nh - nh + eyb * sc / hin            # B up: gap to (B) equals C's gap to (C) (letter bottom = image top + 2 %)
    bC[1] += C_LIFT
    nB, nC = extended(bB, sc), extended(bC, sc)
    axb.set_position(nB); axc.set_position(nC)
    cyc = mpl.rcParamsDefault["axes.prop_cycle"].by_key()["color"]
    handles = [Line2D([], [], color=cyc[i], lw=3) for i in range(len(CONDS))]
    labels = [c.replace("10^", "10$^{") + "}$" if c.startswith("10^") else c for c in CONDS]
    leg = fig.legend(handles, labels, title="CFU/ml", ncol=len(CONDS), loc="upper center", frameon=False,
               bbox_to_anchor=(bC[0] + bC[2] / 2, nC[1] + 0.004 - LEG_DROP), bbox_transform=fig.transFigure,
               prop={"family": "DejaVu Sans", "size": 11}, title_fontproperties={"family": "DejaVu Sans", "size": 12},
               handlelength=1.6, handletextpad=0.4, columnspacing=1.3, borderaxespad=0)
    # vertical colour bar right of panel A, centred between the uidA and papGII heatmaps
    norm, vmax = H.SCALE
    pa = top[0].get_position(fig)
    yc, hgt = (pa.y0 + pa.y1) / 2, 0.36 * (pa.y1 - pa.y0)
    cax = fig.add_axes([pa.x1 + 5 / W, yc - hgt / 2, 3 / W, hgt])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=H.CMAP), cax=cax, ticks=range(0, vmax + 1))
    cb.outline.set_visible(False); cb.ax.tick_params(labelsize=7, colors=H.INK2, length=2)
    cb.set_label("Median mismatches to the primer", fontsize=7.5, color=H.INK2)
    fig.text(nB[0] - 0.03 * nB[2], ya, "(B)", va="bottom", **H.LETTER_KW)
    fig.canvas.draw()                                                             # lay out the colour-bar label before measuring it
    inv = fig.transFigure.inverted()
    a_right = inv.transform((cb.ax.yaxis.label.get_window_extent().x1, 0))[0]     # right edge of panel A (colour-bar label)
    bc_left = nB[0] + axb.left_frac * nB[2]                                    # left edge of the B/C plots (axis label)
    x = (a_right + bc_left) / 2                                                   # divider midway between A and B/C
    y_low = inv.transform(leg.get_window_extent())[0, 1]                          # divider reaches down to the legend
    fig.add_artist(Line2D([x, x], [min(0.0, y_low), 1.0], transform=fig.transFigure, color="black", lw=1.0))
    for ext in ("pdf", "png", "svg"):
        meta = {"pdf": {"CreationDate": None, "Creator": "fig2_composite.py"}, "svg": {"Date": None}}.get(ext)
        fig.savefig(OUT / f"fig2.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight", metadata=meta)


if __name__ == "__main__":
    main()
