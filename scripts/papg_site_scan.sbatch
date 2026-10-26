#!/bin/bash
#SBATCH --job-name=papg_site_scan
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=papg_site_scan_%j.log
#
# Mismatches of the papGII LAMP binding sites to every papG locus of the UPEC genomes (conda env "blast").

set -euo pipefail
A=${A:-/home/jaroge/data/Aline_genomes}
SCRIPT_DIR=${SCRIPT_DIR:-$(pwd)}
OUT=${OUT:-$A/papg_site_scan_$(date +%Y%m%d)}
python3 "$SCRIPT_DIR/papg_site_scan.py" --genomes "${GENOMES:-$A/Aline_fastas}" --pattern "${PATTERN:-*.fasta}" \
  --papg-ref "${PAPG_REF:-/home/jaroge/data/papG_ref_seq/papG_db.fasta}" --out "$OUT" --threads "${SLURM_CPUS_PER_TASK:-16}"
cut -f3 "$OUT/locus_site_mismatches.tsv" | tail -n +2 | sort | uniq -c | awk '{print $2": "$1/8" loci"}'
