"""Figure 5: confusion matrices of the LAMP assays against the reference methods.

(A) E. coli (uidA; time threshold 30 min, fluorescence threshold >= 3-fold) vs culture with MALDI-TOF;
(B) papGII (23 min, >= 3-fold) vs whole genome sequencing. Cells: numbers of samples; margins: sensitivity, specificity, PPV
and NPV in % (whole numbers) with exact (Clopper-Pearson) 95% confidence intervals. Same layout as the original figure; Arial. TP/FP/TN/FN follow lamp_common.time_to_threshold.

Output: results/figures/fig5.{pdf,png,svg}
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle

from diagnostic_metrics_ci import clopper_pearson
from lamp_common import RESULTS, TIME_THRESHOLD, load_samples

OUT = RESULTS / "figures"
mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
                     "pdf.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": "lamp-fig5"})
LETTER = {"fontsize": 16, "fontweight": "bold", "family": "DejaVu Sans"}
PANELS = [("uidA", "ecoli", "$\\it{E. coli}$", "culture", "(A)"), ("papGII", "wgs_papGII", "$\\it{papGII}$", "WGS", "(B)")]


def counts(meta, target, ref_col):
    pos = (meta[f"ttt_{target}"] <= TIME_THRESHOLD[target]).values          # never crossed (NaN) -> negative
    ref = meta[ref_col].values == 1
    return int((pos & ref).sum()), int((pos & ~ref).sum()), int((~pos & ref).sum()), int((~pos & ~ref).sum())   # TP, FP, FN, TN


def pct(a, b):
    """Proportion a/b in % (whole number) with the exact Clopper-Pearson 95% CI, e.g. '87%\n(72–96%)'."""
    lo, hi = clopper_pearson(a, b)
    return f"{100 * a / b:.0f}%\n({100 * lo:.0f}–{100 * hi:.0f}%)"


def panel(ax, tp, fp, fn, tn, gene, ref_name):
    cm = np.array([[tp, fp], [fn, tn]])
    ax.imshow(cm, cmap="Blues", vmin=0, vmax=cm.max(), extent=(0, 2, 2, 0), interpolation="none", aspect="equal", zorder=0)
    lw = 1.8
    b = Rectangle((0, 0), 2, 2, fill=False, linewidth=lw, edgecolor="black", antialiased=False, zorder=3); b.set_clip_on(False); ax.add_patch(b)
    ax.add_collection(LineCollection([[(0, 1), (2, 1)], [(1, 0), (1, 2)]], colors="black", linewidths=lw, antialiased=False, zorder=3))
    for i in range(2):
        for j in range(2):
            ax.text(j + 0.5, i + 0.5, f"{cm[i, j]}", ha="center", va="center", fontsize=14, zorder=4,
                    color="white" if cm[i, j] >= cm.max() * 0.5 else "black")
    ax.set_xlim(0, 2); ax.set_ylim(2, 0)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_xticks([0.5, 1.5]); ax.set_xticklabels(["positive", "negative"], fontsize=12)
    ax.xaxis.tick_top(); ax.tick_params(length=0, pad=4)
    ax.set_yticks([0.5, 1.5]); ax.set_yticklabels(["positive", "negative"], fontsize=12)
    ax.set_title(f"{gene} {ref_name}", fontsize=12.5, pad=30)                   # reference method on top
    ax.set_ylabel(f"{gene} LAMP", fontsize=12.5, labelpad=8)
    for y, lab, val in ((0.5, "PPV", pct(tp, tp + fp)), (1.5, "NPV", pct(tn, tn + fn))):
        ax.text(2.12, y, f"{lab}\n{val}", ha="left", va="center", fontsize=12, fontweight="bold", linespacing=1.15, clip_on=False)
    for x, lab, val in ((0.5, "Sensitivity", pct(tp, tp + fn)), (1.5, "Specificity", pct(tn, tn + fp))):
        ax.text(x, 2.1, f"{lab}\n{val}", ha="center", va="top", fontsize=12, fontweight="bold", linespacing=1.15, clip_on=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    meta = load_samples()
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 3.9), gridspec_kw={"wspace": 1.05})
    for ax, (target, ref_col, gene, ref_name, letter) in zip(axes, PANELS):
        tp, fp, fn, tn = counts(meta, target, ref_col)
        panel(ax, tp, fp, fn, tn, gene, ref_name)
        ax.text(-0.55, 1.28, letter, transform=ax.transAxes, va="bottom", ha="left", **LETTER)
        print(target, dict(TP=tp, FP=fp, FN=fn, TN=tn))
    fig.subplots_adjust(left=0.11, right=0.90, top=0.78, bottom=0.2)
    for ext in ("pdf", "png", "svg"):
        meta_ = {"pdf": {"CreationDate": None, "Creator": "fig5_confusion.py"}, "svg": {"Date": None}}.get(ext)
        fig.savefig(OUT / f"fig5.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight", metadata=meta_)


if __name__ == "__main__":
    main()
