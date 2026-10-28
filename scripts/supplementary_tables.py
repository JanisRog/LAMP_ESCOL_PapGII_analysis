"""Supplementary Tables 1-8 of the manuscript as one Excel workbook (results/supplement/Supplementary_Tables.xlsx).

1  all clinical samples: culture, density, WGS papGII and papG alleles, phylogroup, LAMP calls, times to threshold
2  uidA: highest-ranking threshold combinations (Youden's index) and the selected thresholds
3  papGII: metrics at different time thresholds (fluorescence threshold 3-fold)
4  binding-site mismatches of the primers (UPEC + NCBI genomes) per phylogroup / papG class (data of Figure 2A) and presence of the sites
5  plate counts of the limit-of-detection dilution series and the resulting densities
6  times to threshold of the E. coli-positive urines by density, papGII status and primer set
7  bootstrap of the threshold selection
8  clinical isolates: phylogroup, mismatches at the LF binding site and uidA time to threshold

Run with a Python that has pandas, numpy, scipy and openpyxl.
"""
import re

import numpy as np
import pandas as pd

from lamp_common import DATA, FOLD_CHANGE, MEASUREMENTS, RESULTS, TIME_THRESHOLD, load_samples
from threshold_bootstrap import score_matrix, counts, youden
from diagnostic_metrics_ci import clopper_pearson

OUT = RESULTS / "supplement"
REF = {"uidA": "ecoli", "papGII": "wgs_papGII"}


def metric_table(target, combos):
    """combos: list of (T, f) -> Se, Sp, PPV, NPV, J with counts, from the 191 samples."""
    M, y = score_matrix(target)
    rows = []
    for T, f in combos:
        tp, fp, fn, tn = (int(a[T - 1, 0]) for a in counts(M, y, np.array([f]), np.ones(len(y))))
        se, sp = tp / (tp + fn), tn / (tn + fp)
        ppv = tp / (tp + fp) if tp + fp else float("nan"); npv = tn / (tn + fn) if tn + fn else float("nan")
        rows.append({"Time threshold (min)": T, "Fluorescence threshold (fold-change)": round(float(f), 1), "TP": tp, "FP": fp, "TN": tn, "FN": fn,
                     "Sensitivity (%)": round(100 * se, 1), "Specificity (%)": round(100 * sp, 1), "PPV (%)": round(100 * ppv, 1),
                     "NPV (%)": round(100 * npv, 1), "Youden's index": round(se + sp - 1, 3)})
    return pd.DataFrame(rows)


def table1():
    meta = load_samples()
    meta = meta[meta["sample_id"].notna()].copy()
    ty = pd.read_csv(RESULTS / "papg_site_scan_lamp" / "locus_site_mismatches.tsv", sep="\t")
    alleles = ty.groupby("genome").allele.apply(lambda a: sorted(set(a)))
    ph = pd.read_csv(DATA / "lamp_phylogroups.tsv", sep="\t").set_index("sample_id").phylogroup

    def papg(sid):
        got = set()
        for g in (sid, f"{sid}_1"):                       # the _2 sequencing replicates are not used
            got |= set(alleles.get(g, []))
        return ", ".join(sorted(a.replace("PapG-", "papG") for a in got)) if got else ("none" if meta.loc[meta.sample_id == sid, "ecoli"].iloc[0] == 1 else "")
    meta["papG alleles (WGS)"] = meta.sample_id.map(papg)
    meta["Phylogroup"] = meta.sample_id.map(ph).fillna("")
    out = pd.DataFrame({
        "Sample": meta["sample"], "Sample ID": meta.sample_id,
        "Bacterial density (CFU/ml)": meta.cfu_ml.map(lambda v: "no growth" if v == 0 else ("10^4" if v == 10000 else "10^5 (includes >10^5)")),
        "Culture (MALDI-TOF)": meta.culture_maldi, "E. coli by culture (1 = positive)": meta.ecoli,
        "papGII by WGS (1 = positive)": meta.wgs_papGII, "papG alleles by WGS": meta["papG alleles (WGS)"], "Phylogroup": meta["Phylogroup"],
        "uidA LAMP result": meta.lamp_uidA, "uidA time to threshold (min)": meta.ttt_uidA,
        "papGII LAMP result": meta.lamp_papGII, "papGII time to threshold (min)": meta.ttt_papGII})
    out["uidA LAMP correct"] = (out["uidA LAMP result"] == out["E. coli by culture (1 = positive)"]).map({True: "yes", False: "no"})
    out["papGII LAMP correct"] = (out["papGII LAMP result"] == out["papGII by WGS (1 = positive)"]).map({True: "yes", False: "no"})
    return out


def table4():
    mm = pd.read_csv(RESULTS / "figures" / "tableS_binding_site_mismatches.tsv", sep="\t")
    lol = pd.read_csv(RESULTS / "figures" / "fig2a_lollipop_uidA.tsv", sep="\t"); pg = pd.read_csv(RESULTS / "figures" / "fig2a_lollipop_papGII.tsv", sep="\t")
    return lol, pg, mm


PLATES = {   # colonies per plate (columns 10^6, 10^5, 10^4, 10^3, 10^2; None = too many to count); rows 100 / 10 / 1 ul; source: LOS_counts.xlsx
    "Replicate 1": [[None, None, None, 128, 15], [None, None, 232, 17, 0], [None, 178, 34, 1, 0]],
    "Replicate 2": [[None, None, None, 98, 12], [None, None, 153, 11, 1], [None, 97, 12, 0, 0]],
    "Replicate 3": [[None, None, None, 104, 10], [None, None, 122, 12, 0], [None, 109, 14, 1, 0]]}
VOLUMES_ML = (0.1, 0.01, 0.001)


def plate_tables():
    raw_rows, dens_rows = [], []
    for rep, mat in PLATES.items():
        for vol_label, vals in zip(("100 ul", "10 ul", "1 ul"), mat):
            raw_rows.append({"Replicate": rep, "Plated volume": vol_label,
                             **{f"10^{6 - k}": ("too many to count" if v is None else v) for k, v in enumerate(vals)}})
        row = {"Replicate": rep}
        for k in range(1, 5):                                  # 10^5 ... 10^2
            conc = 6 - k
            usable = [(mat[r][k], VOLUMES_ML[r]) for r in range(3) if mat[r][k] is not None and 30 <= mat[r][k] <= 300]
            n, vol = usable[0] if usable else (mat[0][k], VOLUMES_ML[0])      # < 30 colonies: 100 ul plate as an estimate
            est = n / vol
            row[f"10^{conc}: CFU/ml"] = f"{est:.2e}" + ("" if usable else " (<30 colonies)")
            row[f"10^{conc}: log10 difference"] = round(float(np.log10(est / 10 ** conc)), 2)
        dens_rows.append(row)
    return pd.DataFrame(raw_rows), pd.DataFrame(dens_rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    f = FOLD_CHANGE
    # Table 2: distinct classification outcomes of all (time, fluorescence) combinations, ranked by Youden's index (many combinations give the same outcome)
    M, y = score_matrix("uidA"); fl = np.sort(np.concatenate([np.linspace(M.min(), M.max(), 200), [FOLD_CHANGE, 11.4]]))   # grid + the two thresholds discussed in the text
    tp, fp, fn, tn = counts(M, y, fl, np.ones(len(y))); j = youden(tp, fp, fn, tn)
    keys = np.stack([tp, fp, fn, tn], -1).reshape(-1, 4).astype(int)
    TT, FF = np.meshgrid(np.arange(1, 41), fl, indexing="ij")
    df = pd.DataFrame({"TP": keys[:, 0], "FP": keys[:, 1], "TN": keys[:, 3], "FN": keys[:, 2], "J": j.reshape(-1), "T": TT.reshape(-1), "f": FF.reshape(-1)})
    g = df.groupby(["TP", "FP", "TN", "FN"]).agg(J=("J", "max"), T_min=("T", "min"), T_max=("T", "max"), f_min=("f", "min"), f_max=("f", "max"), n=("T", "size")).reset_index()
    g = g.sort_values(["J", "n"], ascending=[False, False]).head(10)
    t2 = pd.DataFrame({"Rank": range(1, len(g) + 1), "Time thresholds giving this result (min)": [f"{a}-{b}" if a != b else str(a) for a, b in zip(g.T_min, g.T_max)],
                       "Fluorescence thresholds (fold-change)": [f"{a:.1f}-{b:.1f}" for a, b in zip(g.f_min, g.f_max)], "Number of combinations": g.n.values,
                       "TP": g.TP.values, "FP": g.FP.values, "TN": g.TN.values, "FN": g.FN.values,
                       "Sensitivity (%)": (100 * g.TP / (g.TP + g.FN)).round(1).values, "Specificity (%)": (100 * g.TN / (g.TN + g.FP)).round(1).values,
                       "PPV (%)": (100 * g.TP / (g.TP + g.FP)).round(1).values, "NPV (%)": (100 * g.TN / (g.TN + g.FN)).round(1).values, "Youden's index": g.J.round(3).values})
    sel = metric_table("uidA", [(30, f), (39, 11.4)]); sel.insert(0, "Combination", ["selected (30 min, 3-fold)", "39 min, 11.4-fold"])
    # Table 3: papGII at time thresholds, 3-fold
    t3 = metric_table("papGII", [(T, f) for T in list(range(18, 28)) + [30, 39]])
    # Table 6
    cf = pd.read_csv(RESULTS / "cfu_ttt" / "ecoli_time_to_threshold.tsv", sep="\t")
    cf["cfu_ml"] = cf.cfu_ml.map({10000: "10^4", 100000: "10^5 (includes >10^5)"})
    cf = cf.rename(columns={"sample_id": "Sample ID", "cfu_ml": "Bacterial density (CFU/ml)", "culture_maldi": "Culture (MALDI-TOF)", "mixed_culture": "Mixed culture (1 = yes)",
                            "wgs_papGII": "papGII by WGS", "ttt_uidA": "uidA time to threshold (min)", "ttt_papGII": "papGII time to threshold (min)"})
    summ = (RESULTS / "cfu_ttt" / "summary.txt").read_text().splitlines()
    # Table 7
    bs = {t: pd.read_csv(RESULTS / "threshold_bootstrap" / f"selected_thresholds_{t}.tsv", sep="\t") for t in ("uidA", "papGII")}
    rows = []
    for t, d in bs.items():
        T = TIME_THRESHOLD[t]
        rows.append({"Target": t, "Time threshold of the study (min)": T, "Selected time threshold (min), median": d.T_fixed_f.median(), "2.5-97.5% range": f"{d.T_fixed_f.quantile(.025):.0f}-{d.T_fixed_f.quantile(.975):.0f}",
                     "Resamples selecting the study threshold (%)": round(100 * (d.T_fixed_f == T).mean()), "Resamples within +-3 min of it (%)": round(100 * ((d.T_fixed_f - T).abs() <= 3).mean()),
                     "Youden's index in the resample (mean)": round(d.J_fixed_f_in_sample.mean(), 3), "Youden's index out-of-bag (mean)": round(d.J_fixed_f_out_of_bag.mean(), 3),
                     "Optimism": round(d.J_fixed_f_in_sample.mean() - d.J_fixed_f_out_of_bag.mean(), 3)})
    t7 = pd.DataFrame(rows)
    t7b = pd.DataFrame([{"Target": t, "Selected T, median (min)": d.T_selected.median(), "T 2.5-97.5%": f"{d.T_selected.quantile(.025):.0f}-{d.T_selected.quantile(.975):.0f}",
                         "Selected f, median": round(d.f_selected.median(), 1), "Youden's index in the resample": round(d.J_in_sample.mean(), 3),
                         "Youden's index out-of-bag": round(d.J_out_of_bag.mean(), 3)} for t, d in bs.items()])
    # Table 8
    iso = pd.read_csv(RESULTS / "clinical_phylogroups" / "isolates.tsv", sep="\t")
    pres = pd.read_csv(RESULTS / "figures" / "tableS_binding_site_presence.tsv", sep="\t")
    lol, pg, mm = table4()
    t1 = table1()
    path = OUT / "Supplementary_Tables.xlsx"
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        def put(df, name, title, notes=(), startrow=None, first=True):
            r0 = 0 if first else startrow
            df.to_excel(xw, sheet_name=name, index=False, startrow=(len(notes) + 2) if first else startrow)
            ws = xw.sheets[name]
            if first:
                ws.cell(1, 1, title)
                for i, n in enumerate(notes): ws.cell(2 + i, 1, n)
            for col in ws.columns:
                ws.column_dimensions[col[0].column_letter].width = min(max(len(str(c.value)) if c.value is not None else 0 for c in col[2:30]) + 2, 46)
        put(t1, "Supp Table 1", "Supplementary Table 1. Overview of all 191 clinical urine samples, reference results and LAMP results.",
            ["LAMP result: 1 = positive (uidA: <= 30 min and >= 3-fold; papGII: <= 23 min and >= 3-fold); time to threshold empty = fluorescence threshold not reached within 40 min. "
             "Bacterial density from semi-quantitative culture of 10 ul; the category 10^5 includes >10^5 CFU/ml. papG alleles: BLASTN vs papG classes I-IV (>= 90% identity, >= 80% coverage); phylogroups: EzClermont."])
        put(t2, "Supp Table 2", "Supplementary Table 2. uidA: the ten highest-ranking results of the threshold search by Youden's index (191 samples).",
            ["All combinations of 40 time thresholds (1-40 min) and 200 fluorescence thresholds (fold-change relative to the first measurement; equally spaced between the lowest and highest observed value) were evaluated; "
             "combinations with the same classification result are summarised in one row. The selected thresholds (30 min, 3-fold) and the combination of 39 min and 11.4-fold are given below the table."])
        put(sel, "Supp Table 2", "", first=False, startrow=len(t2) + 6)
        put(t3, "Supp Table 3", "Supplementary Table 3. papGII: diagnostic metrics at different time thresholds (fluorescence threshold >= 3-fold change; 191 samples).",
            ["The selected time threshold of 23 min is the midpoint of the range (20-25 min) in which the Youden's index was 1.00."])
        put(lol, "Supp Table 4", "Supplementary Table 4. Median number of mismatches between the primers and their binding sites (UPEC + NCBI RefSeq genomes pooled; data of Figure 2A).",
            ["uidA: per Clermont phylogroup (genomes with the target region present). papGII: per PapG class (best-matching papG locus of each genome), with minimum and maximum."])
        put(pg, "Supp Table 4", "", first=False, startrow=len(lol) + 6)
        put(mm, "Supp Table 4", "", first=False, startrow=len(lol) + len(pg) + 9)
        put(pres, "Supp Table 4", "", first=False, startrow=len(lol) + len(pg) + len(mm) + 12)
        raw, dens = plate_tables()
        put(raw, "Supp Table 5", "Supplementary Table 5. Plate counts of the dilution series used for the limit of detection and the resulting bacterial densities.",
            ["Colonies per plate after plating 100, 10 and 1 ul of each dilution on COS agar. Densities were calculated from plates with 30-300 colonies; for 10^2 CFU/ml (< 30 colonies) the count of the 100 ul plate is given as an estimate. "
             "The density of the 10^6 CFU/ml dilution was inferred from the 10^5 CFU/ml count (10-fold)."])
        put(dens, "Supp Table 5", "", first=False, startrow=len(raw) + 7)
        put(cf, "Supp Table 6", "Supplementary Table 6. Time to threshold of the 38 E. coli-positive urine samples by bacterial density, papGII status and primer set.",
            ["Time to threshold: first minute with a fold-change >= 3 (empty: not reached within 40 min)."])
        put(pd.DataFrame({"Summary (results/cfu_ttt/summary.txt)": summ}), "Supp Table 6", "", first=False, startrow=len(cf) + 7)
        put(t7, "Supp Table 7", "Supplementary Table 7. Bootstrap of the threshold selection (2,000 stratified resamples of the 191 samples).",
            ["(a) Time threshold selected at a fluorescence threshold of 3-fold change (maximum Youden's index; midpoint if reached over a range) and its Youden's index in the resample and in the samples not drawn (out-of-bag)."])
        put(pd.DataFrame({"": ["(b) Both thresholds selected simultaneously from the full grid (ties: higher PPV, then shorter time): the Youden's index is similar over a wide range of combinations, so that this selection is not stable."]}), "Supp Table 7", "", first=False, startrow=len(t7) + 6)
        put(t7b, "Supp Table 7", "", first=False, startrow=len(t7) + 9)
        put(iso, "Supp Table 8", "Supplementary Table 8. E. coli isolates of the clinical urine samples: phylogroup, mismatches at the LF binding site and uidA time to threshold.",
            ["LF_mismatches_3prime: mismatches in the 3'-terminal part of the LF primer binding site."])
    print("saved", path)


if __name__ == "__main__":
    main()
