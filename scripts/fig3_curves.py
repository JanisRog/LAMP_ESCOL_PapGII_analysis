"""Figure 3: fluorescence kinetics of the clinical urine samples (median and IQR of the fold-change F(t)/F(1)).

(A) uidA, E. coli culture-positive vs culture-negative samples; (B) papGII, WGS-positive vs WGS-negative samples.
Same data and layout as the notebook (notebooks/LAMP_analysis.ipynb); the no-template control is not a sample and is excluded.
Style: Arial, sizes as in Figure 2 (title 9.8 pt, axis labels 8.8 pt, tick labels 7.9 pt).

Output: results/figures/fig3.{pdf,png,svg}
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from lamp_common import DATA, MEASUREMENTS, RESULTS

OUT = RESULTS / "figures"
mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
                     "pdf.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": "lamp-fig3"})
TITLE, LABEL, TICK, LEGEND = 9.8, 8.8, 7.9, 9.8
LETTER = {"fontsize": 16, "fontweight": "bold", "family": "DejaVu Sans"}
PANELS = [("uidA", "culture_results.csv", "$\\it{E. coli}$ Culture", "(A)"), ("papGII", "papg_results.csv", "$\\it{papGII}$ WGS", "(B)")]


def curves(target, ref_file):
    m = pd.read_csv(DATA / MEASUREMENTS[target], index_col="Cycle")
    m = m[[c for c in m.columns if c.strip().isdigit()]]                  # samples only (drops the NTC column)
    fc = m / m.iloc[0]
    ref = pd.read_csv(DATA / ref_file); ref["sample"] = ref["sample"].astype(str).str.strip()
    ref = ref[ref["sample"] != "NTC"].set_index("sample").iloc[:, 0]
    out = {}
    for label in (1, 0):
        cols = [c for c in fc.columns if ref.get(c.strip()) == label]
        g = fc[cols]
        out[label] = (g.median(axis=1), g.quantile(0.25, axis=1), g.quantile(0.75, axis=1), len(cols))
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    left, box_w, gap, right, bottom, box_h, top = 17.0, 105.0, 30.0, 14.0, 26.0, 63.0, 13.0   # mm (bottom: room for the legends)
    W, H = left + 2 * box_w + gap + right, bottom + box_h + top
    fig = plt.figure(figsize=(W / 25.4, H / 25.4))
    for i, (target, ref_file, name, letter) in enumerate(PANELS):
        x0 = left + i * (box_w + gap)
        ax = fig.add_axes([x0 / W, bottom / H, box_w / W, box_h / H])
        cv = curves(target, ref_file)
        for label, color, word in ((1, "red", "Positive"), (0, "blue", "Negative")):
            med, q1, q3, n = cv[label]
            ax.plot(med.index, med.values, color=color, lw=1.4, label=f"{name} {word} (n = {n})")
            ax.fill_between(med.index, q1.values, q3.values, color=color, alpha=0.2, lw=0)
        ax.set_xlim(0, 40); ax.set_ylim(bottom=0); ax.grid(True)
        ax.set_xlabel("Time [min]", fontsize=LABEL); ax.set_ylabel("Relative Fluorescence (fold-change)", fontsize=LABEL)
        ax.tick_params(labelsize=TICK)
        ax.set_title(f"Fluorescence change over time $\\it{{{target}}}$ (Median ± IQR)", fontsize=TITLE, pad=7)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.23), ncol=2, frameon=False, fontsize=LEGEND,
                  handlelength=1.8, columnspacing=1.6, borderaxespad=0)      # legend below the plot, as in Figure 2
        fig.text((x0 - 14) / W, (H - 1.5) / H, letter, va="top", **LETTER)
    for ext in ("pdf", "png", "svg"):
        meta = {"pdf": {"CreationDate": None, "Creator": "fig3_curves.py"}, "svg": {"Date": None}}.get(ext)
        fig.savefig(OUT / f"fig3.{ext}", dpi=600 if ext == "png" else None, metadata=meta)


if __name__ == "__main__":
    main()
