#!/bin/bash
#SBATCH --job-name=phylogroups
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=phylogroups_%j.log
#
# In silico Clermont phylogroups of the UPEC genomes with EzClermont (conda env "blast" with ezclermont installed).
# Output: $OUT/phylogroups.tsv (genome, phylogroup) and $OUT/manifest.txt

set -euo pipefail
A=${A:-/home/jaroge/data/Aline_genomes}
GENOMES=${GENOMES:-$A/Aline_fastas}
OUT=${OUT:-$A/phylogroups_$(date +%Y%m%d)}
mkdir -p "$OUT/per_genome"

ls "$GENOMES"/*.fasta | xargs -P "${SLURM_CPUS_PER_TASK:-16}" -I{} sh -c \
  'g=$(basename {} .fasta); ezclermont {} -e "$g" > "'"$OUT"'/per_genome/$g.tsv" 2> "'"$OUT"'/per_genome/$g.log" || true'

printf "genome\tphylogroup\n" > "$OUT/phylogroups.tsv"
cat "$OUT"/per_genome/*.tsv | sort >> "$OUT/phylogroups.tsv"
cat > "$OUT/manifest.txt" <<EOF
date: $(date -Iseconds)
genomes: $GENOMES ($(ls "$GENOMES"/*.fasta | wc -l) *.fasta)
tool: $(ezclermont --version 2>&1 | head -1)
command: ezclermont <genome.fasta> -e <genome>
EOF
echo "typed: $(($(wc -l < "$OUT/phylogroups.tsv") - 1))"; cut -f2 "$OUT/phylogroups.tsv" | tail -n +2 | sort | uniq -c
