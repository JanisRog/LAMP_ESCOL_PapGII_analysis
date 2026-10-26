#!/bin/bash
#SBATCH --job-name=primer_screen
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=primer_screen_%j.log
#
# In silico screen of the uidA / papGII LAMP primers against the UPEC collection (ref 10).
# Needs blastn and python3 on PATH (conda env "blast").

set -euo pipefail
GENOMES=${GENOMES:-/home/jaroge/data/Aline_genomes/Aline_fastas}
PATTERN=${PATTERN:-*.fasta}
DESIGN_REF=${DESIGN_REF:-$(ls /home/jaroge/data/Aline_genomes/prokka_output_aline/700223-20/*.fna | head -1)}  # E. coli 700223-20: primer design strain
REF_UIDA=${REF_UIDA:-$DESIGN_REF}
REF_PAPGII=${REF_PAPGII:-$DESIGN_REF}
PAPG_REF=${PAPG_REF:-/home/jaroge/data/papG_ref_seq/papG_db.fasta}
OUT=${OUT:-/home/jaroge/data/Aline_genomes/primer_screen_$(date +%Y%m%d)}
SCRIPT_DIR=${SCRIPT_DIR:-$(pwd)}

python3 "$SCRIPT_DIR/primer_screen.py" --genomes "$GENOMES" --pattern "$PATTERN" \
  --ref-uidA "$REF_UIDA" --ref-papGII "$REF_PAPGII" --papg-ref "$PAPG_REF" --out "$OUT" --threads "${SLURM_CPUS_PER_TASK:-16}"
