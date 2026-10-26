"""Shared helpers: load clinical LAMP fluorescence data and compute time to threshold (TTT)."""
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"
RESULTS = REPO / "results"

FOLD_CHANGE = 3.0                         # fluorescence threshold (fold-change vs first cycle)
TIME_THRESHOLD = {"uidA": 30, "papGII": 23}  # minutes; a sample is positive if TTT <= threshold
MEASUREMENTS = {"uidA": "all_measurements_uida.csv", "papGII": "all_measurements_papg.csv"}


def time_to_threshold(target, fold=FOLD_CHANGE):
    """TTT (min) per clinical sample: first cycle (1 cycle = 1 min) with fluorescence >= fold x cycle 1.
    NaN = never crossed within the run. The water negative control (column 'NTC') is excluded."""
    m = pd.read_csv(DATA / MEASUREMENTS[target], index_col="Cycle")
    m = m.drop(columns=[c for c in m.columns if not c.strip().isdigit()])
    fc = m / m.iloc[0]
    ttt = fc.apply(lambda c: float(c.index[c.values >= fold][0]) if (c.values >= fold).any() else np.nan)
    ttt.index = ttt.index.astype(int)
    return ttt.rename(f"ttt_{target}")


def load_samples():
    """Sample metadata (Supplementary Table 1) joined with TTT for both targets."""
    meta = pd.read_csv(DATA / "sample_metadata.csv")
    meta["ecoli"] = meta.culture_maldi.str.contains("Escherichia coli").astype(int)
    meta["mixed_culture"] = meta.culture_maldi.str.contains(",").astype(int)
    for target in MEASUREMENTS:
        meta = meta.merge(time_to_threshold(target), left_on="sample", right_index=True, how="left")
    return meta


def write_readme(folder, text):
    (folder / "README.md").write_text(text.strip() + "\n")
