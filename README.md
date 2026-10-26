# LAMP assay for E. coli (uidA) and papGII – analysis code

Code and data to reproduce the analyses of the LAMP assay for rapid detection of *Escherichia coli* (uidA) and
papGII directly from urine (Rogenmoser et al.). Raw sequencing reads: ENA project PRJEB106156.

## Data (`data/`)

| File | Content |
|---|---|
| `all_measurements_uida.csv`, `all_measurements_papg.csv` | Clinical fluorescent LAMP: raw fluorescence per cycle (1 cycle = 1 min, 40 min) for 191 urines (columns = sample number) and the water negative control (`NTC`) |
| `culture_results.csv`, `papg_results.csv` | Reference results per sample: E. coli by culture + MALDI-TOF, papGII by WGS (1 = positive) |
| `sample_metadata.csv` | Supplementary Table 1: semi-quantitative culture load (CFU/ml), culture result (MALDI-TOF; comma = mixed culture), papGII by WGS, LAMP result |
| `uida_lod_*.csv`, `papg_lod_*.csv` | Limit of detection: fluorescence over 60 min for the dilution series (10⁶–10² CFU/ml, NTC), three replicates |

## Analyses

| Code | Output | Content |
|---|---|---|
| `notebooks/LAMP_analysis.ipynb` | figures | Limit of detection, amplification curves, threshold selection (Youden's index), confusion matrices |
| `scripts/diagnostic_metrics_ci.py` | `results/diagnostic_metrics/` | Sensitivity, specificity, PPV, NPV with exact (Clopper–Pearson) 95% CI |
| `scripts/cfu_time_to_threshold.py` | `results/cfu_ttt/` | Time to threshold vs bacterial load, papGII status, mixed culture and primer set |
| `scripts/primer_screen.py` (+ `primer_screen.sbatch`) | `results/primer_screen/` | In silico screen of both primer sets against 738 UPEC genomes and the 38 clinical isolates, with papG allele typing |

Shared settings (fluorescence threshold ≥ 3-fold of cycle 1; time thresholds uidA 30 min, papGII 23 min) are in
`scripts/lamp_common.py`. Each results folder contains a README describing its files and columns.

## Running

```bash
pip install -r requirements.txt
cd scripts
python diagnostic_metrics_ci.py
python cfu_time_to_threshold.py
```

The primer screen needs BLAST+ (`blastn`) and the genome assemblies; it was run on the UZH ScienceCluster via
`primer_screen.sbatch`. Inputs, versions and parameters of the run are recorded in `results/primer_screen/*/manifest.txt`.
