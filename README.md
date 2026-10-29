# LAMP assay for E. coli (uidA) and papGII – analysis code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23740154.svg)](https://doi.org/10.5281/zenodo.23740154)

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
| Fig. 3 (fluorescence kinetics, clinical samples) | `scripts/fig3_curves.py` | `results/figures/fig3.*` |
| Fig. 4 (diagnostic metrics vs thresholds, ROC) | `scripts/fig4_metrics.py` | `results/figures/fig4.*` |
| Fig. 5 (confusion matrices) | `scripts/fig5_confusion.py` | `results/figures/fig5.*` |
| Fig. 6 (time-to-threshold distributions) | `scripts/fig6_ttt.py` | `results/figures/fig6.*` |
| Fig. 2B/C; threshold selection (Youden's index); confusion matrices | `notebooks/LAMP_analysis.ipynb` (runs top to bottom, figures displayed inline) | `notebooks/` holds the stored exports: confusion matrices (`confusion_matrix_*`) and lists of false-positive/false-negative samples (`FP_*`, `FN_*`) |
| Sensitivity, specificity, PPV, NPV with exact 95% CI | `scripts/diagnostic_metrics_ci.py` | `results/diagnostic_metrics/` |
| Classification of a new LAMP run with the study thresholds (TTT, positive/negative, NTC validity check) | `scripts/classify_lamp_run.py` (e.g. `python scripts/classify_lamp_run.py --target papGII --input data/all_measurements_papg.csv`) | table of calls on screen or `--output calls.tsv` |
| Supplementary Fig. S1 (colorimetric LAMP photos) | `scripts/figS1_colorimetric.py` (cut-outs in `data/colorimetric/`; `--photos DIR` rebuilds them from the original HEIC photos) | `results/figures/figS1.*` |
| Supplementary Tables 1-8 (one Excel workbook) | `scripts/supplementary_tables.py` (needs openpyxl) | `results/supplement/Supplementary_Tables.xlsx` |
| Stability of the threshold selection (bootstrap, out-of-bag performance) | `scripts/threshold_bootstrap.py` | `results/threshold_bootstrap/` |
| Time to threshold vs bacterial density, papGII status, primer set (Supp. table) | `scripts/cfu_time_to_threshold.py` | `results/cfu_ttt/` |
| Fig. 2 (A: primer maps and median mismatches per binding site by phylogroup and papG class; B, C: LoD amplification plots) | `scripts/fig2_composite.py` (panel A from `fig_primer_heatmap.py`, `fig_primer_map.py`; `fig_primer_lollipop.py` provides the shared helpers and the tables `fig2a_*.tsv`) | `results/figures/fig2.*`, `fig2a_*.tsv` |
| Agreement between the UPEC and NCBI genomes (Spearman) | `scripts/collection_agreement.py` | `results/collection_agreement/` |
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
| `colorimetric/crop_<min>_<E|P>.jpg` | Cut-outs of the photos of the colorimetric LAMP experiment (Supplementary Fig. S1): row E = uidA primers, row P = papGII primers |
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
| Gene-wide variability, amplicon position | `gene_scan.py` / `.sbatch` (UPEC); `ncbi_gene_scan.sbatch` (NCBI RefSeq genomes, by accession) | `results/gene_scan/`, `results/gene_scan_ncbi/` |
| uidA in NCBI RefSeq genomes (download, screen, phylogroups) | `ncbi_ecoli_uida.sbatch` | `results/ncbi_ecoli_uida/` |
| Phylogroups | `phylogroups.sbatch` | `data/*_phylogroups.tsv` |
| Contig end vs genuine truncation | `uida_nonintact_check.py` | `results/uida_nonintact/` |

## License

Code (`scripts/`, `notebooks/`): MIT License (`LICENSE`). Data and results (`data/`, `results/`): CC BY 4.0 (`LICENSE-DATA.md`). How to cite: `CITATION.cff`.

## Use of AI tools

Parts of the code were developed with the assistance of an AI tool (Claude, Anthropic). All code was reviewed and verified by the authors; the scripts reproduce the numbers reported in the manuscript.
