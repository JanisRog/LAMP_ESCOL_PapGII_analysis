# LAMP assay for E. coli (uidA) and papGII – analysis code

Code and data to reproduce the analyses of the LAMP assay for rapid detection of *Escherichia coli* (uidA) and papGII
directly from urine (Rogenmoser et al.). Raw sequencing reads: ENA project PRJEB106156.

## Quick start

```bash
pip install -r requirements.txt          # Python 3.9.6 was used
bash scripts/run_local_analyses.sh       # regenerates results/ and the figures from data/ and the stored cluster outputs
```

The local analyses need only the files in `data/` and `results/` (cluster outputs are stored in the repository). The
cluster-side steps (genome analyses) are listed below under "Genome analyses" and need BLAST+ and the genome assemblies.

## Which code produces which manuscript item

| Manuscript item | Code | Output |
|---|---|---|
| Fig. 2B/C, 3–6; threshold selection (Youden's index); confusion matrices | `notebooks/LAMP_analysis.ipynb` | figures, `notebooks/FP_*`, `FN_*` |
| Sensitivity, specificity, PPV, NPV with exact 95% CI | `scripts/diagnostic_metrics_ci.py` | `results/diagnostic_metrics/` |
| Time to threshold vs bacterial density, papGII status, primer set (Supp. table) | `scripts/cfu_time_to_threshold.py` | `results/cfu_ttt/` |
| Fig. 2A (primer maps and binding-site conservation) | `scripts/fig_primer_map.py` | `results/figures/fig2a_*` |
| Supp. tables on binding-site mismatches and presence | `scripts/tables_s_primer_conservation.py` | `results/figures/tableS_*` |
| Clinical isolates: phylogroup, LF mismatches, time to threshold (Supp. table) | `scripts/clinical_phylogroups.py` | `results/clinical_phylogroups/` |
| uidA target region intact/absent/truncated; contig-end check | `scripts/uida_ncbi_screen.py`, `scripts/uida_nonintact_check.py` | `results/ncbi_ecoli_uida/`, `results/uida_nonintact/` |

Exploratory analyses not shown in the manuscript: `scripts/uidA_candidate_sets.py` (comparison of the five uidA primer
sets proposed by the NEB LAMP Primer Design Tool, `results/uidA_candidate_sets/`).

## Data (`data/`)

| File | Content |
|---|---|
| `all_measurements_uida.csv`, `all_measurements_papg.csv` | Clinical fluorescent LAMP: raw fluorescence per cycle (1 cycle = 1 min, 40 min) for 191 urines (columns = sample number) and the water negative control (`NTC`) |
| `culture_results.csv`, `papg_results.csv` | Reference results per sample: E. coli by culture + MALDI-TOF, papGII by WGS (1 = positive) |
| `sample_metadata.csv` | Supplementary Table 1: semi-quantitative culture density (CFU/ml; the category 10⁵ includes >10⁵), culture result (MALDI-TOF; comma = mixed culture), papGII by WGS, LAMP result |
| `uida_lod_*.csv`, `papg_lod_*.csv` | Limit of detection: fluorescence over 60 min for the dilution series (10⁶–10² CFU/ml, NTC), three replicates |
| `papG_reference_alleles.fasta` | PapG-I to PapG-IV reference alleles used for papG typing (PapG-IV: GenBank AF304159) |
| `upec_phylogroups.tsv`, `lamp_phylogroups.tsv` | Clermont phylogroups (EzClermont 1.0.0) of the 738 UPEC genomes and the clinical E. coli isolates |
| `lamp_ecoli_genome_ids.txt` | Genome IDs of the clinical E. coli isolates (input list for the contig-end check) |

Shared settings (fluorescence threshold ≥ 3-fold of cycle 1; time thresholds uidA 30 min, papGII 23 min) are in
`scripts/lamp_common.py`. Each results folder has a README describing its files and columns, and the cluster runs
record inputs, versions and thresholds in `manifest.txt`.

## Genome analyses (cluster)

The genome analyses were run on the UZH ScienceCluster with BLAST+ (2.17.0+) in a conda environment (`blast`: blast,
python 3, ncbi-datasets-cli, ezclermont). Paths at the top of each `.sbatch` file are defaults for our setup and can be
overridden with environment variables. Collections: 738 UPEC genomes (Cuénod et al., Genome Med 2023), the clinical
isolates (ENA PRJEB106156) and all 5,493 NCBI RefSeq complete *E. coli* assemblies available on 2026-10-27
(accession list: `results/ncbi_ecoli_uida/assemblies.tsv`).

| Step | Script | Results |
|---|---|---|
| Binding-site screen, papG typing | `primer_screen.py` / `.sbatch` | `results/primer_screen/` |
| papGII sites vs every papG locus | `papg_site_scan.py` / `.sbatch` | `results/papg_site_scan*/` |
| Gene-wide variability, amplicon position | `gene_scan.py` / `.sbatch` | `results/gene_scan/` |
| uidA in NCBI RefSeq genomes (download, screen, phylogroups) | `ncbi_ecoli_uida.sbatch` | `results/ncbi_ecoli_uida/` |
| Phylogroups | `phylogroups.sbatch` | `data/*_phylogroups.tsv` |
| Contig end vs genuine truncation | `uida_nonintact_check.py` | `results/uida_nonintact/` |
