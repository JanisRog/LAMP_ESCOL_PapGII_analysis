# LAMP assay for E. coli (uidA) and papGII – analysis code

Code and data to reproduce the analyses of the LAMP assay for rapid detection of *Escherichia coli* (uidA) and
papGII directly from urine (Rogenmoser et al.). Raw sequencing reads: ENA project PRJEB106156.

## Data (`data/`)

| File | Content |
|---|---|
| `all_measurements_uida.csv`, `all_measurements_papg.csv` | Clinical fluorescent LAMP: raw fluorescence per cycle (1 cycle = 1 min, 40 min) for 191 urines (columns = sample number) and the water negative control (`NTC`) |
| `culture_results.csv`, `papg_results.csv` | Reference results per sample: E. coli by culture + MALDI-TOF, papGII by WGS (1 = positive) |
| `sample_metadata.csv` | Supplementary Table 1: semi-quantitative culture load (CFU/ml), culture result (MALDI-TOF; comma = mixed culture), papGII by WGS, LAMP result |
| `lamp_phylogroups.tsv` | Clermont phylogroups of the 38 clinical E. coli isolates (EzClermont 1.0.0) |
| `upec_phylogroups.tsv` | Clermont phylogroups of the 738 UPEC genomes (EzClermont 1.0.0; `scripts/phylogroups.sbatch`) |
| `papG_reference_alleles.fasta` | PapG-I, PapG-II and PapG-III reference alleles used for papG typing and papGII primer design |
| `uida_lod_*.csv`, `papg_lod_*.csv` | Limit of detection: fluorescence over 60 min for the dilution series (10⁶–10² CFU/ml, NTC), three replicates |

## Analyses

| Code | Output | Content |
|---|---|---|
| `notebooks/LAMP_analysis.ipynb` | figures | Limit of detection, amplification curves, threshold selection (Youden's index), confusion matrices |
| `scripts/diagnostic_metrics_ci.py` | `results/diagnostic_metrics/` | Sensitivity, specificity, PPV, NPV with exact (Clopper–Pearson) 95% CI |
| `scripts/cfu_time_to_threshold.py` | `results/cfu_ttt/` | Time to threshold vs bacterial load, papGII status, mixed culture and primer set |
| `scripts/primer_screen.py` (+ `primer_screen.sbatch`) | `results/primer_screen/` | In silico screen of both primer sets against 738 UPEC genomes and the 38 clinical isolates, with papG allele typing |
| `scripts/gene_scan.py` (+ `gene_scan.sbatch`) | `results/gene_scan/` | Per-position variability of uidA (738 genomes) and papGII (140 papGII genomes); amplicon and binding sites in gene coordinates |
| `scripts/uidA_candidate_sets.py` | `results/uidA_candidate_sets/` | Comparison of the five uidA primer sets proposed by the NEB LAMP Primer Design Tool |
| `scripts/clinical_phylogroups.py` | `results/clinical_phylogroups/` | Clinical isolates: phylogroup, uidA LF primer mismatches and time to threshold |
| `scripts/papg_site_scan.py` (+ `papg_site_scan.sbatch`) | `results/papg_site_scan/` | Mismatches of the papGII binding sites to every papG locus (papGII, papGIII) of the 738 UPEC genomes |
| `scripts/fig_primer_map.py` | `results/figures/` | Figure 2A: primer binding sites on both amplicons; median mismatches of the uidA sites per phylogroup and of the papGII sites per papG allele |
| `scripts/tables_s_primer_conservation.py` | `results/figures/` | Supplementary tables: per-site mismatches to the primers in 738 UPEC genomes (descriptive) and presence of the papGII binding sites per papG allele group |

Shared settings (fluorescence threshold ≥ 3-fold of cycle 1; time thresholds uidA 30 min, papGII 23 min) are in
`scripts/lamp_common.py`. Each results folder contains a README describing its files and columns.

## Running

```bash
pip install -r requirements.txt
cd scripts
python diagnostic_metrics_ci.py
python cfu_time_to_threshold.py
python uidA_candidate_sets.py
python clinical_phylogroups.py
python fig_primer_map.py
python tables_s_primer_conservation.py
```

The primer screen, gene scan and papG site scan need BLAST+ (`blastn`) and the genome assemblies; all were run on the UZH ScienceCluster via
`primer_screen.sbatch`, `gene_scan.sbatch` and `papg_site_scan.sbatch`. Inputs, versions and parameters of each run are recorded in the
`manifest.txt` of its results folder.
