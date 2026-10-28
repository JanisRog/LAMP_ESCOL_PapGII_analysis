"""Stability of the threshold selection (bootstrap) and out-of-sample performance of the selected thresholds.

The thresholds (time T = 1-40 min, fluorescence f = 200 equally spaced values between the lowest and highest observed fold-change) were chosen by the
maximum Youden's index in the same 191 samples. To show how stable this choice is, the complete selection is repeated in B bootstrap resamples
(stratified by reference status, so that each resample has the same number of positive and negative samples): the combination (T, f) with the highest
Youden's index in the resample is selected (ties: higher PPV, then shorter T), and the selected combination is evaluated on the samples that were not
drawn (out-of-bag). The same is done for the fixed rule used in the study (uidA: 30 min / 3-fold, papGII: 23 min / 3-fold).

Output: results/threshold_bootstrap/{summary.txt, selected_thresholds_<target>.tsv, README.md}
"""
import numpy as np
import pandas as pd

from lamp_common import DATA, MEASUREMENTS, RESULTS, TIME_THRESHOLD, FOLD_CHANGE, load_samples, write_readme

OUT = RESULTS / "threshold_bootstrap"
B, SEED = 2000, 20261028
REF = {"uidA": "ecoli", "papGII": "wgs_papGII"}


def score_matrix(target):
    """Per-sample maximum fold-change up to each minute (rows: 1..40 min) for the 191 samples, and the reference labels."""
    m = pd.read_csv(DATA / MEASUREMENTS[target], index_col="Cycle")
    m = m[[c for c in m.columns if c.strip().isdigit()]]
    fc = (m / m.iloc[0]).cummax()
    meta = load_samples().set_index(load_samples()["sample"].astype(str))
    cols = fc.columns.str.strip()
    y = meta.loc[cols, REF[target]].values.astype(int)
    return fc.values, y


def counts(M, y, f, w):
    """TP, FP, FN, TN for every (T, f) with sample weights w. M: (40, n); f: (nf,); returns arrays of shape (40, nf)."""
    pred = M[:, None, :] >= f[None, :, None]
    wp, wn = w * (y == 1), w * (y == 0)
    tp = (pred * wp).sum(-1); fn = wp.sum() - tp
    fp = (pred * wn).sum(-1); tn = wn.sum() - fp
    return tp, fp, fn, tn


def youden(tp, fp, fn, tn):
    return tp / np.maximum(tp + fn, 1e-9) + tn / np.maximum(tn + fp, 1e-9) - 1


def select(M, y, f, w):
    tp, fp, fn, tn = counts(M, y, f, w)
    j = youden(tp, fp, fn, tn); ppv = tp / np.maximum(tp + fp, 1e-9)
    best = j.max()
    cand = np.argwhere(j >= best - 1e-12)
    # ties: higher PPV, then shorter time
    order = sorted(cand, key=lambda ij: (-ppv[ij[0], ij[1]], ij[0]))
    return order[0][0], order[0][1], best


def select_T_fixed_f(M, y, w, f=FOLD_CHANGE):
    """Study rule at a fixed fluorescence threshold: the time threshold with the maximum Youden's index; where the maximum is reached over a range of
    times, the midpoint of the range (rounded up). Returns the index of T (0-based) and the Youden's index."""
    tp, fp, fn, tn = counts(M, y, np.array([f]), w)
    j = youden(tp, fp, fn, tn)[:, 0]
    top = np.where(j >= j.max() - 1e-12)[0]
    return int(np.ceil((top.min() + top.max()) / 2)), j.max()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    lines = [f"Bootstrap of the threshold selection: B = {B} stratified resamples of the 191 samples (seed {SEED}); selection by the maximum Youden's index "
             f"(ties: PPV, then shorter time); evaluation on the out-of-bag samples.\n"]
    for target in ("uidA", "papGII"):
        M, y = score_matrix(target)
        f = np.linspace(M.min(), M.max(), 200)
        fixed_T, fixed_f = TIME_THRESHOLD[target], FOLD_CHANGE
        iT = fixed_T - 1
        n = len(y); idx_pos, idx_neg = np.where(y == 1)[0], np.where(y == 0)[0]
        rows = []
        full = select(M, y, f, np.ones(n))
        for b in range(B):
            draw = np.concatenate([rng.choice(idx_pos, len(idx_pos)), rng.choice(idx_neg, len(idx_neg))])
            w = np.bincount(draw, minlength=n).astype(float)
            oob = (w == 0).astype(float)
            iT_b, if_b, j_in = select(M, y, f, w)
            tp, fp, fn, tn = counts(M[:, :], y, f[if_b:if_b + 1], oob)
            j_oob = float(youden(tp, fp, fn, tn)[iT_b, 0])
            tp, fp, fn, tn = counts(M, y, np.array([fixed_f]), oob); j_fix_oob = float(youden(tp, fp, fn, tn)[iT, 0])
            tp, fp, fn, tn = counts(M, y, np.array([fixed_f]), w); j_fix_in = float(youden(tp, fp, fn, tn)[iT, 0])
            iT_m, j_m = select_T_fixed_f(M, y, w)
            tp, fp, fn, tn = counts(M, y, np.array([fixed_f]), oob); j_m_oob = float(youden(tp, fp, fn, tn)[iT_m, 0])
            rows.append(dict(T_fixed_f=iT_m + 1, J_fixed_f_in_sample=j_m, J_fixed_f_out_of_bag=j_m_oob, T_selected=iT_b + 1, f_selected=f[if_b], J_in_sample=j_in, J_out_of_bag=j_oob,
                             J_fixed_rule_in_sample=j_fix_in, J_fixed_rule_out_of_bag=j_fix_oob))
        d = pd.DataFrame(rows)
        d.to_csv(OUT / f"selected_thresholds_{target}.tsv", sep="\t", index=False, lineterminator="\n", float_format="%.4f")
        q = lambda s, p: s.quantile(p)
        near = ((d.T_selected - fixed_T).abs() <= 3).mean()
        lines += [f"== {target} - (a) time threshold selected at the fixed fluorescence threshold of {fixed_f:g}-fold (study rule: maximum Youden's index, midpoint of a plateau)",
                  f"  selected T: median {d.T_fixed_f.median():.0f} min (IQR {q(d.T_fixed_f, .25):.0f}-{q(d.T_fixed_f, .75):.0f}; 2.5-97.5%: {q(d.T_fixed_f, .025):.0f}-{q(d.T_fixed_f, .975):.0f}); "
                  f"within +-3 min of {fixed_T} min: {100 * ((d.T_fixed_f - fixed_T).abs() <= 3).mean():.0f}%; exactly {fixed_T}: {100 * (d.T_fixed_f == fixed_T).mean():.0f}%",
                  f"  Youden's index of the selected T: in the resample {d.J_fixed_f_in_sample.mean():.3f}, out-of-bag {d.J_fixed_f_out_of_bag.mean():.3f} (optimism {d.J_fixed_f_in_sample.mean() - d.J_fixed_f_out_of_bag.mean():.3f})",
                  f"== {target} - (b) both thresholds selected on the full grid (maximum Youden's index; ties: PPV, then shorter T)",
                  f"  (study rule: {fixed_T} min, >= {fixed_f:g}-fold; full-data optimum: T = {full[0] + 1} min, f = {f[full[1]]:.1f}, J = {full[2]:.3f}; the Youden index is flat over a wide range of combinations)",
                  f"  selected time threshold T: median {d.T_selected.median():.0f} min (IQR {q(d.T_selected, .25):.0f}-{q(d.T_selected, .75):.0f}; 2.5-97.5%: {q(d.T_selected, .025):.0f}-{q(d.T_selected, .975):.0f})",
                  f"  resamples selecting T within +-3 min of the fixed rule: {100 * near:.0f}%; selecting exactly T = {fixed_T}: {100 * (d.T_selected == fixed_T).mean():.0f}%",
                  f"  selected fluorescence threshold f: median {d.f_selected.median():.1f} (IQR {q(d.f_selected, .25):.1f}-{q(d.f_selected, .75):.1f}); f <= 4: {100 * (d.f_selected <= 4).mean():.0f}%",
                  f"  Youden's index of the selected combination: in the resample {d.J_in_sample.mean():.3f}, out-of-bag {d.J_out_of_bag.mean():.3f} "
                  f"(optimism {d.J_in_sample.mean() - d.J_out_of_bag.mean():.3f})",
                  f"  Youden's index of the fixed rule: in the resample {d.J_fixed_rule_in_sample.mean():.3f}, out-of-bag {d.J_fixed_rule_out_of_bag.mean():.3f}; "
                  f"2.5-97.5% out-of-bag: {q(d.J_fixed_rule_out_of_bag, .025):.3f}-{q(d.J_fixed_rule_out_of_bag, .975):.3f}\n"]
    (OUT / "summary.txt").write_text("\n".join(lines))
    write_readme(OUT, f"""# Bootstrap of the threshold selection

Generated by `scripts/threshold_bootstrap.py`. {B} stratified bootstrap resamples of the 191 clinical samples; in each, the combination of time threshold (1-40 min) and
fluorescence threshold (200 values) with the highest Youden's index is selected and evaluated on the out-of-bag samples. See `summary.txt`.
- selected_thresholds_<target>.tsv - per resample: selected T and f, Youden's index in the resample and out-of-bag, and the same for the fixed rule used in the study.
""")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
