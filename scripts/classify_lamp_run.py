"""Classify the samples of a fluorescent LAMP run with the thresholds of the study.

Input: a CSV exported from the real-time instrument, one row per cycle (1 cycle = 1 min) and one column per sample,
plus a column 'Cycle' and, optionally, a no-template control column 'NTC' (same layout as data/all_measurements_*.csv).
Per sample: fluorescence is normalised to the first cycle, the time to threshold (TTT) is the first cycle with a
fold-change >= 3, and the sample is positive if TTT <= the time threshold of the target (uidA 30 min, papGII 23 min).
The run is valid if the NTC does not reach the fluorescence threshold within the time threshold.

Usage: python classify_lamp_run.py --target uidA|papGII --input run.csv [--output calls.tsv]
       [--fold 3] [--time-threshold 30]
Example (clinical data of the study): python classify_lamp_run.py --target papGII --input ../data/all_measurements_papg.csv
"""
import argparse

import numpy as np
import pandas as pd

from lamp_common import FOLD_CHANGE, TIME_THRESHOLD


def classify(csv, target, fold=FOLD_CHANGE, time_threshold=None):
    time_threshold = TIME_THRESHOLD[target] if time_threshold is None else time_threshold
    m = pd.read_csv(csv, index_col="Cycle")
    fc = m / m.iloc[0]

    def first_crossing(col):
        hit = col.index[col.values >= fold]
        return float(hit[0]) if len(hit) else np.nan

    ttt = fc.apply(first_crossing)
    ntc_ttt = ttt.pop("NTC") if "NTC" in ttt.index else np.nan
    calls = pd.DataFrame({"sample": ttt.index, "ttt_min": ttt.values})
    calls["call"] = np.where(calls.ttt_min <= time_threshold, "positive", "negative")      # NaN (never crossed) -> negative
    valid = bool(np.isnan(ntc_ttt) or ntc_ttt > time_threshold) if "NTC" in fc.columns else None
    return calls, valid, ntc_ttt, time_threshold


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", required=True, choices=sorted(TIME_THRESHOLD))
    ap.add_argument("--input", required=True)
    ap.add_argument("--output")
    ap.add_argument("--fold", type=float, default=FOLD_CHANGE)
    ap.add_argument("--time-threshold", type=float)
    a = ap.parse_args()
    calls, valid, ntc_ttt, tt = classify(a.input, a.target, a.fold, a.time_threshold)
    n_pos = int((calls.call == "positive").sum())
    ntc_txt = "no NTC column" if valid is None else ("NTC did not cross the threshold" if np.isnan(ntc_ttt) else f"NTC crossed at {ntc_ttt:.0f} min")
    print(f"{a.target}: {len(calls)} samples, {n_pos} positive (fold >= {a.fold:g} within {tt:g} min); "
          f"run {'VALID' if valid else 'INVALID' if valid is False else 'not checked'} ({ntc_txt})")
    if a.output:
        calls.to_csv(a.output, sep="\t", index=False, lineterminator="\n")


if __name__ == "__main__":
    main()
