#!/bin/bash
#SBATCH --job-name=gene_scan
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=gene_scan_%j.log
#
# Gene-wide variability of uidA (all UPEC genomes) and papGII (genomes typed papGII by primer_screen.py),
# with the LAMP amplicon located in each gene. Needs blastn and python3 (conda env "blast").

set -euo pipefail
A=${A:-/home/jaroge/data/Aline_genomes}
SCRIPT_DIR=${SCRIPT_DIR:-$(pwd)}
SCREEN=${SCREEN:-$(ls -d $A/primer_screen_20*/ | tail -1)}
DESIGN=${DESIGN:-$(ls $A/prokka_output_aline/700223-20/*.fna | head -1)}
UIDA_QUERY=${UIDA_QUERY:-$(ls $A/extracted_sequences/beta_glucuronidase_sequences/700223-20_*.fasta | head -1)}
PAPG_DB=${PAPG_DB:-/home/jaroge/data/papG_ref_seq/papG_db.fasta}
OUT=${OUT:-$A/gene_scan_$(date +%Y%m%d)}
T=${SLURM_CPUS_PER_TASK:-16}
mkdir -p "$OUT"

# PapG-II reference allele and the list of papGII genomes from the primer screen
awk '/^>/{p=($1==">PapG-II")} p' "$PAPG_DB" > "$OUT/PapG-II_query.fasta"
awk -F'\t' 'NR>1 && $12 ~ /PapG-II\(/ {print $1}' "$SCREEN/genome_summary.tsv" > "$OUT/papGII_genomes.txt"
echo "papGII genomes: $(wc -l < "$OUT/papGII_genomes.txt") (from $SCREEN)"

python3 "$SCRIPT_DIR/gene_scan.py" --target uidA --gene-query "$UIDA_QUERY" --design-genome "$DESIGN" \
  --genomes "$A/Aline_fastas" --out "$OUT/uidA" --threads "$T"
python3 "$SCRIPT_DIR/gene_scan.py" --target papGII --gene-query "$OUT/PapG-II_query.fasta" --design-genome "$DESIGN" \
  --genomes "$A/Aline_fastas" --genome-list "$OUT/papGII_genomes.txt" --out "$OUT/papGII" --threads "$T"
