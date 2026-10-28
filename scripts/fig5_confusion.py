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
PANELS = [("uidA", "ecoli", "$\\it{E. coli}$", "culture", "(A)"), ("papGII", "wgs_papGII", "$\\it{papGII}$", "WGS", "(B)")]


def counts(meta, target, ref_col):
    pos = (meta[f"ttt_{target}"] <= TIME_THRESHOLD[target]).values          # never crossed (NaN) -> negative
    ref = meta[ref_col].values == 1
    return int((pos & ref).sum()), int((pos & ~ref).sum()), int((~pos & ref).sum()), int((~pos & ~ref).sum())   # TP, FP, FN, TN


def pct(a, b):
    """Proportion a/b in % (whole number) with the exact Clopper-Pearson 95% CI, e.g. '87%\n(72–96%)'."""
    lo, hi = clopper_pearson(a, b)
    return f"{100 * a / b:.0f}%\n({100 * lo:.0f}–{100 * hi:.0f}%)"


# Geometry of the original Figure 5, measured pixel by pixel on the 2944 x 1042 px original (300 dpi); panel B = panel A shifted by DX.
DPI = 300
W_PX, H_PX = 3075, 1070                       # canvas: the original (2944 x 1042) plus room for the CI lines (right of PPV/NPV, below Sens./Spec.)
FRAME = {"x0": 301.5, "x1": 1164.5, "y0": 208.5, "y1": 892.5}          # centre lines of the outer frame (px)
DX = {"A": 0, "B": 1597}
LW = 1.50                                      # frame lines: 6 px at 300 dpi (non-antialiased rendering rounds 1.44 pt down to 5 px)
FS_NUM, FS_TICK, FS_TITLE, FS_YLAB, FS_RIGHT, FS_BOTTOM = 12.5, 11.45, 11.4, 10.1, 9.7, 10.7    # pt (Arial sizes that reproduce the widths of the original text)
LETTER = {"fontsize": 14.2, "fontweight": "bold", "family": "DejaVu Sans"}


def panel(fig, key, tp, fp, fn, tn, gene, ref_name, letter):
    dx = DX[key]
    px = lambda x, y: ((x + dx) / W_PX, 1 - y / H_PX)                       # original-pixel -> figure fraction
    f = FRAME
    ax = fig.add_axes([(f["x0"] + dx) / W_PX, 1 - f["y1"] / H_PX, (f["x1"] - f["x0"]) / W_PX, (f["y1"] - f["y0"]) / H_PX])
    cm = np.array([[tp, fp], [fn, tn]])
    ax.imshow(cm, cmap="Blues", vmin=0, vmax=cm.max(), extent=(0, 2, 2, 0), interpolation="none", aspect="auto", zorder=0)
    bd = Rectangle((0, 0), 2, 2, fill=False, linewidth=LW, edgecolor="black", antialiased=False, joinstyle="miter", zorder=3); bd.set_clip_on(False); ax.add_patch(bd)
    lc = LineCollection([[(0, 1), (2, 1)], [(1, 0), (1, 2)]], colors="black", linewidths=LW, antialiased=False, capstyle="butt", zorder=3)
    ax.add_collection(lc)
    for i in range(2):
        for j in range(2):
            ax.text(j + 0.5, i + 0.5 + 0.006, f"{cm[i, j]}", ha="center", va="center", fontsize=FS_NUM, zorder=4,
                    color="white" if cm[i, j] >= cm.max() * 0.5 else "black")
    ax.set_xlim(0, 2); ax.set_ylim(2, 0); ax.axis("off")
    cx = [517.5, 949.5]                                                     # column centres (px)
    ry = [379.5, 721.5]                                                     # row centres (px)
    for x, lab in zip(cx, ("positive", "negative")):
        fig.text(*px(x, 187), lab, ha="center", va="baseline", fontsize=FS_TICK)                         # top: reference result
    for y, lab in zip(ry, ("positive", "negative")):
        fig.text(*px(287, y + 2), lab, ha="right", va="center_baseline", fontsize=FS_TICK)                   # left: LAMP result
    fig.text(*px(733, 120), f"{gene} {ref_name}", ha="center", va="baseline", fontsize=FS_TITLE)
    fig.text(*px(52, 552), f"{gene} LAMP", rotation=90, ha="center", va="center", fontsize=FS_YLAB)
    fig.text(*px(39 if key == "A" else 44, 53), letter, ha="left", va="center", **LETTER)
    for y, lab, val in zip((400.5, 742.5), ("PPV", "NPV"), (pct(tp, tp + fp), pct(tn, tn + fn))):       # right margin, left-aligned
        fig.text(*px(1229, y), f"{lab}\n{val}", ha="left", va="center", fontsize=FS_RIGHT, fontweight="bold", linespacing=1.14)
    for x, lab, val in zip(cx, ("Sensitivity", "Specificity"), (pct(tp, tp + fn), pct(tn, tn + fp))):  # bottom, centred
        fig.text(*px(x, 935), f"{lab}\n{val}", ha="center", va="top", fontsize=FS_BOTTOM, fontweight="bold", linespacing=1.0)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    meta = load_samples()
    fig = plt.figure(figsize=(W_PX / DPI, H_PX / DPI), dpi=DPI)
    for target, ref_col, gene, ref_name, letter in PANELS:
        tp, fp, fn, tn = counts(meta, target, ref_col)
        panel(fig, "A" if target == "uidA" else "B", tp, fp, fn, tn, gene, ref_name, letter)
        print(target, dict(TP=tp, FP=fp, FN=fn, TN=tn))
    for ext in ("pdf", "png", "svg"):
        meta_ = {"pdf": {"CreationDate": None, "Creator": "fig5_confusion.py"}, "svg": {"Date": None}}.get(ext)
        fig.savefig(OUT / f"fig5.{ext}", dpi=2 * DPI if ext == "png" else None, metadata=meta_)       # PNG at 600 dpi (twice the 300-dpi geometry)


if __name__ == "__main__":
    main()
