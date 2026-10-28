"""Figure 6: distribution of the time to threshold (fluorescence >= 3-fold change) of the clinical samples.

(A) uidA, E. coli culture-positive vs culture-negative samples; (B) papGII, WGS-positive vs WGS-negative samples.
Same design as the notebook (notebooks/LAMP_analysis.ipynb: overlaid histograms, 30 bins between 0 and 40 min); the no-template control
is not a sample and is excluded. Samples that never reached the threshold within the 40-minute run have no time to threshold and are not
shown (their number is given in the legend of the figure). Style: Arial, sizes as in Figure 2 B/C (scale 0.465), legends below the plots.

Output: results/figures/fig6.{pdf,png,svg}
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

from lamp_common import RESULTS, load_samples

OUT = RESULTS / "figures"
SC = 0.465   # Figure 2 B/C scale: fonts 21 / 19 / 17 pt on the 10 x 6 in notebook canvas; line weights scaled the same way
TITLE, LABEL, TICK, LEGEND = 21 * SC, 19 * SC, 17 * SC, 21 * SC
mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
                     "pdf.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": "lamp-fig6",
                     "axes.linewidth": 0.8 * SC, "grid.linewidth": 0.8 * SC, "xtick.major.size": 3.5 * SC, "ytick.major.size": 3.5 * SC,
                     "xtick.major.width": 0.8 * SC, "ytick.major.width": 0.8 * SC})
LETTER = {"fontsize": 16, "fontweight": "bold", "family": "DejaVu Sans"}
LETTER_SHIFT_B = 4.0   # mm: panel letter (B) further left, as in Figure 3
PANELS = [("uidA", "ecoli", "$\\it{E. coli}$ Culture", "(A)"), ("papGII", "wgs_papGII", "$\\it{papGII}$ WGS", "(B)")]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    meta = load_samples()
    left, box_w, gap, right, bottom, box_h, top = 17.0, 105.0, 30.0, 14.0, 26.0, 63.0, 13.0      # mm (as in Figure 3)
    W, H = left + 2 * box_w + gap + right, bottom + box_h + top
    fig = plt.figure(figsize=(W / 25.4, H / 25.4))
    bins = np.linspace(0, 40, 30)
    for i, (target, ref_col, name, letter) in enumerate(PANELS):
        x0 = left + i * (box_w + gap)
        ax = fig.add_axes([x0 / W, bottom / H, box_w / W, box_h / H])
        ref = meta[ref_col].values == 1
        ttt = meta[f"ttt_{target}"].values
        pos, neg = ttt[ref], ttt[~ref]
        n_pos, n_neg = int(ref.sum()), int((~ref).sum())
        ax.hist(pos[~np.isnan(pos)], bins=bins, alpha=0.6, color="red")
        ax.hist(neg[~np.isnan(neg)], bins=bins, alpha=0.6, color="blue")
        ax.set_xlim(0, 40); ax.grid(True, ls="--", alpha=0.6)
        ax.set_xlabel("Time to Threshold [min]", fontsize=LABEL); ax.set_ylabel("Number of Samples", fontsize=LABEL)
        ax.tick_params(labelsize=TICK)
        ax.set_title(f"Distribution of Time-to-Threshold $\\it{{{target}}}$ (Fluorescence $\\geq$ 3-fold change)", fontsize=TITLE, pad=7)
        handles = [Patch(facecolor="red", alpha=0.6, edgecolor="none"), Patch(facecolor="blue", alpha=0.6, edgecolor="none")]
        labels = [f"{name} Positive (n = {n_pos})", f"{name} Negative (n = {n_neg})"]
        ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.23), ncol=2, frameon=False, fontsize=LEGEND,
                  handlelength=1.4, handleheight=1.0, columnspacing=1.6, borderaxespad=0)
        fig.text((x0 - 14 - (LETTER_SHIFT_B if i == 1 else 0)) / W, (H - 1.5) / H, letter, va="top", **LETTER)
        print(target, "never crossed: positives", int(np.isnan(pos).sum()), "negatives", int(np.isnan(neg).sum()), "of", n_neg)
    for ext in ("pdf", "png", "svg"):
        meta_ = {"pdf": {"CreationDate": None, "Creator": "fig6_ttt.py"}, "svg": {"Date": None}}.get(ext)
        fig.savefig(OUT / f"fig6.{ext}", dpi=600 if ext == "png" else None, metadata=meta_)


if __name__ == "__main__":
    main()
