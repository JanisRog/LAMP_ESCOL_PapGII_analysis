"""Figure 4: diagnostic performance of the uidA and papGII LAMP assays as a function of the thresholds.

Left column uidA (A-C), right column papGII (D-F). Per target: (top) Youden's index, sensitivity, specificity and PPV vs time
threshold at a fixed fluorescence threshold (>= 3-fold); (middle) the same metrics vs fluorescence threshold at a fixed time
threshold (uidA 30 min, papGII 23 min); (bottom) ROC curve of the maximum fold-change up to the selected time threshold.
Same definitions as the notebook (notebooks/LAMP_analysis.ipynb); only the 191 samples are used (the no-template control, which the
notebook counted as a negative sample, is excluded). One shared legend below the figure.
Style: Arial, sizes as in Figure 2 B/C (scale 0.465).

Output: results/figures/fig4.{pdf,png,svg}
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from sklearn.metrics import auc, roc_curve

from lamp_common import DATA, MEASUREMENTS, RESULTS

OUT = RESULTS / "figures"
SC = 0.465   # Figure 2 B/C scale: fonts 21 / 19 / 17 pt on the 10 x 6 in notebook canvas; line weights scaled the same way
TITLE, LABEL, TICK, LEGEND = 21 * SC, 19 * SC, 17 * SC, 21 * SC
mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
                     "pdf.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": "lamp-fig4",
                     "lines.linewidth": 1.5 * SC, "axes.linewidth": 0.8 * SC, "grid.linewidth": 0.8 * SC,
                     "xtick.major.size": 3.5 * SC, "ytick.major.size": 3.5 * SC, "xtick.major.width": 0.8 * SC, "ytick.major.width": 0.8 * SC})
# colour-blind-safe palette (Okabe-Ito; the four metric colours pass the CVD check of the dataviz validator, the old blue/orange/green/red did not:
# orange and green were indistinguishable for protan viewers). Line styles (solid / dotted / dash-dot) add a second cue.
BLUE, ORANGE, GREEN, VERM = "#0072B2", "#E69F00", "#009E73", "#D55E00"
SEL, UIDA = "black", "#CC79A7"          # selected threshold / uidA time threshold (reference marks): black and reddish purple
LETTER = {"fontsize": 16, "fontweight": "bold", "family": "DejaVu Sans"}
FIXED_F = 3.0
TARGETS = [("uidA", "culture_results.csv", 30), ("papGII", "papg_results.csv", 23)]
UIDA_T = 30


def load(target, ref_file):
    m = pd.read_csv(DATA / MEASUREMENTS[target], index_col="Cycle")
    m = m[[c for c in m.columns if c.strip().isdigit()]]                 # samples only (no NTC)
    fc = m / m.iloc[0]
    ref = pd.read_csv(DATA / ref_file); ref["sample"] = ref["sample"].astype(str).str.strip()
    y = ref[ref["sample"] != "NTC"].set_index("sample").iloc[:, 0].reindex(fc.columns.str.strip()).values.astype(int)
    return fc, y


def metrics(score, y, f):
    pred = score >= f
    tp, fn = int((pred & (y == 1)).sum()), int((~pred & (y == 1)).sum())
    tn, fp = int((~pred & (y == 0)).sum()), int((pred & (y == 0)).sum())
    se = tp / (tp + fn) if tp + fn else 0.0
    sp = tn / (tn + fp) if tn + fp else 0.0
    ppv = tp / (tp + fp) if tp + fp else 0.0
    return se + sp - 1, se, sp, ppv


def style(ax, xlabel, ylabel, title):
    ax.grid(True); ax.tick_params(labelsize=TICK)
    ax.set_xlabel(xlabel, fontsize=LABEL); ax.set_ylabel(ylabel, fontsize=LABEL)
    ax.set_title(title, fontsize=TITLE, pad=7)


LW = 2.6 * SC   # metric curves in panels A, B, D, E (thicker than the default for readability)


def metric_lines(ax, x, vals):
    ax.plot(x, vals[0], color=BLUE, lw=LW * 1.25)
    ax.plot(x, vals[1], color=ORANGE, ls=":", lw=LW * 1.2)
    ax.plot(x, vals[2], color=GREEN, ls="-.", lw=LW)
    ax.plot(x, vals[3], color=VERM, lw=LW)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    left, box_w, gap, right = 17.0, 105.0, 30.0, 10.0                     # mm
    row_h, pitch, top, legend_h = 63.0, 90.0, 12.0, 22.0                  # box height, row pitch (box + title + labels), top margin, legend room
    W = left + 2 * box_w + gap + right
    roc = 86.0                                                            # ROC plot box (square), larger than the metric plots are high
    H = top + 2 * pitch + (roc + 21) + legend_h - 8 - 2                  # whitespace above the shared legend
    fig = plt.figure(figsize=(W / 25.4, H / 25.4))
    mpl.rcParams["axes.titlepad"] = 7

    def axes_at(col, row, w, h, dx=0.0):                                   # row 0 = top
        x0 = left + col * (box_w + gap) + dx
        y_top = H - top - row * pitch
        return fig.add_axes([x0 / W, (y_top - h) / H, w / W, h / H]), x0

    for col, (target, ref_file, T) in enumerate(TARGETS):
        fc, y = load(target, ref_file)
        gene = f"$\\it{{{target}}}$"
        # --- time sweep at f = 3
        times = fc.index.values
        tv = np.array([metrics(fc.loc[:t].max().values, y, FIXED_F) for t in times]).T
        ax, x0 = axes_at(col, 0, box_w, row_h)
        metric_lines(ax, times, tv)
        ax.axvline(T, color=SEL, ls="--", lw=LW)
        if target == "papGII": ax.axvline(UIDA_T, color=UIDA, ls="--", lw=LW)
        ax.set_xlim(-1, 41); style(ax, "Time [min]", "Score", f"Diagnostic Metrics for {gene} vs Time (Fluorescence ≥ 3-fold change)")
        fig.text((x0 - 14) / W, (H - 1.5 - 0 * pitch) / H, f"({'ABCDEF'[col * 3]})", va="top", **LETTER)
        # --- fluorescence sweep at the selected T
        mx = fc.loc[:T].max().values
        fluos = np.linspace(fc.values.min(), fc.values.max(), 200)
        fv = np.array([metrics(mx, y, f) for f in fluos]).T
        ax, x0 = axes_at(col, 1, box_w, row_h)
        metric_lines(ax, fluos, fv)
        ax.axvline(FIXED_F, color=SEL, ls="--", lw=LW)
        ax.set_xlim(-0.5, fluos.max() + 0.5)
        style(ax, "Fluorescence Threshold (fold-change)", "Score", f"Diagnostic Metrics for {gene} vs Fluorescence (Time = {T} min)")
        fig.text((x0 - 14) / W, (H - 1.5 - pitch) / H, f"({'ABCDEF'[col * 3 + 1]})", va="top", **LETTER)
        # --- ROC at T
        fpr, tpr, _ = roc_curve(y, mx); a = auc(fpr, tpr)
        ax, x0 = axes_at(col, 2, roc, roc, dx=(box_w - roc) / 2)
        ax.plot(fpr, tpr, color=BLUE, lw=LW * 1.25)
        ax.plot([0, 1], [0, 1], ls="--", color="gray", lw=LW * 0.6)
        se, sp = metrics(mx, y, FIXED_F)[1:3]
        ax.scatter(1 - sp, se, s=80, facecolors="none", edgecolors=SEL, linewidths=LW, zorder=3)
        if target == "papGII":
            m30 = fc.loc[:UIDA_T].max().values
            se30, sp30 = metrics(m30, y, FIXED_F)[1:3]
            ax.scatter(1 - sp30, se30, s=80, facecolors="none", edgecolors=UIDA, linewidths=LW, zorder=3)
        roc_bottom = (H - top - 2 * pitch - roc) / H                          # bottom of the ROC axes (figure fraction)
        ax.text(0.97, 0.06, f"AUC = {a:.3f}", ha="right", va="bottom", transform=ax.transAxes, fontsize=LABEL)
        ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.05, 1.05)
        style(ax, "False Positive Rate (1 - Specificity)", "True Positive Rate (Sensitivity)", f"ROC Curve for {gene} at {T} min")
        ax.grid(True, ls="--", alpha=0.6)
        fig.text((x0 - 14 - (box_w - roc) / 2) / W, (H - 1.5 - 2 * pitch) / H, f"({'ABCDEF'[col * 3 + 2]})", va="top", **LETTER)

    # shared legend below the figure (row 1: metrics; row 2: thresholds and cutoffs)
    h = [Line2D([], [], color=BLUE, lw=3), Line2D([], [], color=ORANGE, lw=3, ls=":"), Line2D([], [], color=GREEN, lw=3, ls="-."),
         Line2D([], [], color=VERM, lw=3),
         Line2D([], [], color=SEL, lw=2.4, ls="--"), Line2D([], [], color=SEL, marker="o", mfc="none", mew=1.5, ms=8, ls=""),
         Line2D([], [], color=UIDA, lw=2.4, ls="--"), Line2D([], [], color=UIDA, marker="o", mfc="none", mew=1.5, ms=8, ls="")]
    lab = ["Youden's index", "Sensitivity", "Specificity", "PPV",
           "Selected threshold", "Selected cutoff (ROC)", "uidA time threshold (30 min)", "uidA cutoff (ROC)"]
    order = [0, 4, 1, 5, 2, 6, 3, 7]                                       # column-major fill -> two tidy rows
    leg = fig.legend([h[i] for i in order], [lab[i] for i in order], ncol=4, loc="lower center", frameon=False, fontsize=LEGEND,
               bbox_to_anchor=(0.5, 0.0), handlelength=2.4, columnspacing=1.8, borderaxespad=0.2)
    fig.canvas.draw()
    y_leg = roc_bottom - 12 / H - 6 / H                                                    # divider ends just below the ROC axis labels
    x_div = (left + box_w + (left + box_w + gap - 14)) / 2 / W                              # midway between the left column and the (D) letter
    fig.add_artist(Line2D([x_div, x_div], [y_leg + 6 / H, 1.0], transform=fig.transFigure, color="black", lw=1.0))
    for ext in ("pdf", "png", "svg"):
        meta = {"pdf": {"CreationDate": None, "Creator": "fig4_metrics.py"}, "svg": {"Date": None}}.get(ext)
        fig.savefig(OUT / f"fig4.{ext}", dpi=600 if ext == "png" else None, metadata=meta)


if __name__ == "__main__":
    main()
