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


# Layout of the original Figure 5 (pixel positions in the 2000-px-wide rendering of the original, 9.8 in wide; converted to inches)
IN = 9.8 / 2000
FIG_W_PX, FIG_H_PX = 2095, 708            # a little wider than the original (room for the CI line next to PPV / NPV of panel B)
MATRIX = {"A": (205, 792), "B": (1290, 1877)}            # left / right edge of the matrix (px)
M_TOP, M_BOT = 142, 607                                  # top / bottom edge of the matrices (px)
X_LETTER, X_YLAB, X_RIGHT = {"A": 30, "B": 1115}, {"A": 35, "B": 1118}, {"A": 833, "B": 1918}
FS_NUM, FS_TICK, FS_TITLE, FS_MARGIN = 13, 11.5, 11.5, 10.6


def panel(fig, key, tp, fp, fn, tn, gene, ref_name, letter):
    x0, x1 = MATRIX[key]
    px = lambda x, y: (x / FIG_W_PX, 1 - y / FIG_H_PX)                          # pixel -> figure fraction
    ax = fig.add_axes([x0 / FIG_W_PX, 1 - M_BOT / FIG_H_PX, (x1 - x0) / FIG_W_PX, (M_BOT - M_TOP) / FIG_H_PX])
    cm = np.array([[tp, fp], [fn, tn]])
    ax.imshow(cm, cmap="Blues", vmin=0, vmax=cm.max(), extent=(0, 2, 2, 0), interpolation="none", aspect="auto", zorder=0)
    lw = 2.0
    bd = Rectangle((0, 0), 2, 2, fill=False, linewidth=lw, edgecolor="black", antialiased=False, zorder=3); bd.set_clip_on(False); ax.add_patch(bd)
    ax.add_collection(LineCollection([[(0, 1), (2, 1)], [(1, 0), (1, 2)]], colors="black", linewidths=lw, antialiased=False, zorder=3))
    for i in range(2):
        for j in range(2):
            ax.text(j + 0.5, i + 0.5, f"{cm[i, j]}", ha="center", va="center", fontsize=FS_NUM, zorder=4,
                    color="white" if cm[i, j] >= cm.max() * 0.5 else "black")
    ax.set_xlim(0, 2); ax.set_ylim(2, 0); ax.axis("off")
    cx = [(x0 + (x1 - x0) * f) for f in (0.25, 0.75)]                          # column centres (px)
    ry = [M_TOP + (M_BOT - M_TOP) * f for f in (0.25, 0.75)]                   # row centres (px)
    for x, lab in zip(cx, ("positive", "negative")):
        fig.text(*px(x, 117), lab, ha="center", va="center", fontsize=FS_TICK)                         # top: reference
    for y, lab in zip(ry, ("positive", "negative")):
        fig.text(*px(x0 - 10, y), lab, ha="right", va="center", fontsize=FS_TICK)                      # left: LAMP
    fig.text(*px((x0 + x1) / 2, 72), f"{gene} {ref_name}", ha="center", va="center", fontsize=FS_TITLE)
    fig.text(*px(X_YLAB[key], (M_TOP + M_BOT) / 2), f"{gene} LAMP", rotation=90, ha="center", va="center", fontsize=FS_TITLE)
    fig.text(*px(X_LETTER[key], 38), letter, ha="left", va="center", **LETTER)
    for y, lab, val in zip(ry, ("PPV", "NPV"), (pct(tp, tp + fp), pct(tn, tn + fn))):                 # right margin, left-aligned as in the original
        fig.text(*px(X_RIGHT[key], y), f"{lab}\n{val}", ha="left", va="center", fontsize=FS_MARGIN, fontweight="bold", linespacing=1.2)
    for x, lab, val in zip(cx, ("Sensitivity", "Specificity"), (pct(tp, tp + fn), pct(tn, tn + fp))):  # bottom, centred under the columns
        fig.text(*px(x, 628), f"{lab}\n{val}", ha="center", va="top", fontsize=FS_MARGIN, fontweight="bold", linespacing=1.2)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    meta = load_samples()
    fig = plt.figure(figsize=(FIG_W_PX * IN, (FIG_H_PX + 40) * IN))
    for target, ref_col, gene, ref_name, letter in PANELS:
        tp, fp, fn, tn = counts(meta, target, ref_col)
        panel(fig, "A" if target == "uidA" else "B", tp, fp, fn, tn, gene, ref_name, letter)
        print(target, dict(TP=tp, FP=fp, FN=fn, TN=tn))
    for ext in ("pdf", "png", "svg"):
        meta_ = {"pdf": {"CreationDate": None, "Creator": "fig5_confusion.py"}, "svg": {"Date": None}}.get(ext)
        fig.savefig(OUT / f"fig5.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight", metadata=meta_)


if __name__ == "__main__":
    main()
