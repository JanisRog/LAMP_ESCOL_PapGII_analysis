#!/bin/bash
#SBATCH --job-name=ncbi_gene_scan
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=ncbi_gene_scan_%j.log
#
# Gene-wide variability of uidA (all NCBI RefSeq complete E. coli genomes used for the uidA screen) and of papGII
# (the NCBI genomes in which papG typing found a PapG-II locus), relative to the gene of E. coli 700223-20.
# 1. re-download exactly the genomes listed in assemblies.tsv (by accession)  2. gene_scan.py for uidA and papGII
# 3. delete the genomes (KEEP_GENOMES=1 keeps them). Needs conda env "blast" (blast, ncbi-datasets-cli, python 3).

set -euo pipefail
A=${A:-/home/jaroge/data/Aline_genomes}
SCRIPT_DIR=${SCRIPT_DIR:-$(pwd)}
D=$(date +%Y%m%d)
OUT=${OUT:-$A/ncbi_gene_scan_$D}
WORK=${WORK:-$A/ncbi_work_$D}   # ~27 GB; deleted at the end unless KEEP_GENOMES=1
ASM=${ASM:-$(ls -d $A/ncbi_ecoli_uida_2026*/ | tail -1)/assemblies.tsv}
DESIGN=${DESIGN:-$(ls $A/prokka_output_aline/700223-20/*.fna | head -1)}
UIDA_QUERY=${UIDA_QUERY:-$(ls $A/extracted_sequences/beta_glucuronidase_sequences/700223-20_*.fasta | head -1)}
PAPG_DB=${PAPG_DB:-/home/jaroge/data/papG_ref_seq/papG_db_I-IV.fasta}
PAPG_SCAN=${PAPG_SCAN:-$A/papg_site_scan_ncbi/locus_site_mismatches.tsv}
T=${SLURM_CPUS_PER_TASK:-16}

for f in "$ASM" "$DESIGN" "$UIDA_QUERY" "$PAPG_DB" "$PAPG_SCAN" "$SCRIPT_DIR/gene_scan.py" "$SCRIPT_DIR/primer_screen.py"; do
  [ -f "$f" ] || { echo "ERROR: missing input $f"; exit 1; }
done
command -v datasets >/dev/null || { echo "ERROR: datasets not found (conda activate blast)"; exit 1; }
command -v blastn >/dev/null || { echo "ERROR: blastn not found (conda activate blast)"; exit 1; }
mkdir -p "$OUT" "$WORK"

# accession list (first column of assemblies.tsv, header removed) and the papGII genome list
tail -n +2 "$ASM" | cut -f1 | sort -u > "$OUT/accessions.txt"
N=$(wc -l < "$OUT/accessions.txt")
echo "accessions: $N (from $ASM)"
awk '/^>/{p=($1==">PapG-II")} p' "$PAPG_DB" > "$OUT/PapG-II_query.fasta"
[ -s "$OUT/PapG-II_query.fasta" ] || { echo "ERROR: PapG-II not found in $PAPG_DB"; exit 1; }
awk -F'\t' 'NR>1 && $3=="PapG-II"{print $1}' "$PAPG_SCAN" | sort -u > "$OUT/papGII_genomes.txt"
echo "papGII genomes: $(wc -l < "$OUT/papGII_genomes.txt")"

echo "== download $(date -Iseconds)"
cd "$WORK"
datasets download genome accession --inputfile "$OUT/accessions.txt" --include genome --dehydrated --filename ecoli.zip
unzip -q -o ecoli.zip -d ecoli
datasets rehydrate --directory ecoli --max-workers 10
G=$(find ecoli -name '*_genomic.fna' | wc -l)
echo "genomes downloaded: $G of $N; size: $(du -sh ecoli | cut -f1)"
[ "$G" -ge $((N * 95 / 100)) ] || { echo "ERROR: fewer than 95% of the genomes were downloaded"; exit 1; }

echo "== gene scan uidA $(date -Iseconds)"
python3 "$SCRIPT_DIR/gene_scan.py" --target uidA --gene-query "$UIDA_QUERY" --design-genome "$DESIGN" \
  --genomes "$WORK/ecoli" --pattern "*_genomic.fna" --out "$OUT/uidA" --threads "$T"
echo "== gene scan papGII $(date -Iseconds)"
python3 "$SCRIPT_DIR/gene_scan.py" --target papGII --gene-query "$OUT/PapG-II_query.fasta" --design-genome "$DESIGN" \
  --genomes "$WORK/ecoli" --pattern "*_genomic.fna" --genome-list "$OUT/papGII_genomes.txt" --out "$OUT/papGII" --threads "$T"

if [ "${KEEP_GENOMES:-0}" != "1" ]; then rm -rf "$WORK"; echo "genomes deleted"; fi
echo "done: $OUT"; ls "$OUT"
