#!/bin/bash
#SBATCH --job-name=ncbi_ecoli_uida
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=ncbi_ecoli_uida_%j.log
#
# uidA binding sites in all NCBI RefSeq complete E. coli genomes.
# 1. download (NCBI datasets, dehydrated + rehydrate)  2. uidA screen  3. Clermont phylogroups (EzClermont)
# 4. delete the genomes (KEEP_GENOMES=1 keeps them); results, accession list and logs stay in $OUT.
# Needs conda env "blast" with: blast, ncbi-datasets-cli, ezclermont (see the submit instructions).

set -euo pipefail
A=${A:-/home/jaroge/data/Aline_genomes}
SCRIPT_DIR=${SCRIPT_DIR:-$(pwd)}
DESIGN=${DESIGN:-$(ls $A/prokka_output_aline/700223-20/*.fna | head -1)}
D=$(date +%Y%m%d)
OUT=${OUT:-$A/ncbi_ecoli_uida_$D}
WORK=${WORK:-$A/ncbi_work_$D}   # ~30 GB; deleted at the end unless KEEP_GENOMES=1
T=${SLURM_CPUS_PER_TASK:-16}
mkdir -p "$OUT" "$WORK"

echo "== download $(date -Iseconds)"
cd "$WORK"
datasets summary genome taxon "Escherichia coli" --assembly-source RefSeq --assembly-level complete --as-json-lines \
  | dataformat tsv genome --fields accession,organism-name,assminfo-name,assmstats-total-sequence-len > "$OUT/assemblies.tsv"
echo "assemblies listed: $(($(wc -l < "$OUT/assemblies.tsv") - 1))"
datasets download genome taxon "Escherichia coli" --assembly-source RefSeq --assembly-level complete \
  --include genome --dehydrated --filename ecoli.zip
unzip -q -o ecoli.zip -d ecoli
datasets rehydrate --directory ecoli --max-workers 10
echo "genomes downloaded: $(find ecoli -name '*_genomic.fna' | wc -l); size: $(du -sh ecoli | cut -f1)"

echo "== uidA screen $(date -Iseconds)"
python3 "$SCRIPT_DIR/uida_ncbi_screen.py" --genomes ecoli --pattern "*_genomic.fna" --design-genome "$DESIGN" \
  --out "$OUT/uida_screen" --threads "$T"

echo "== phylogroups $(date -Iseconds)"
GENOMES="$WORK/ecoli" PATTERN="*_genomic.fna" OUT="$OUT/phylogroups" bash "$SCRIPT_DIR/phylogroups.sbatch"

if [ "${KEEP_GENOMES:-0}" != "1" ]; then rm -rf "$WORK"; echo "genomes deleted"; fi
echo "done: $OUT"
