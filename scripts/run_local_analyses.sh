#!/bin/bash
# Regenerate all locally reproducible results (tables, statistics, Figure 2A) from data/ and the stored cluster outputs.
# Usage: bash scripts/run_local_analyses.sh   (from anywhere; needs the packages in requirements.txt)
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
for s in diagnostic_metrics_ci cfu_time_to_threshold clinical_phylogroups tables_s_primer_conservation uidA_candidate_sets fig_primer_map fig_primer_lollipop fig_primer_heatmap fig2_composite collection_agreement; do
  echo "== $s"; $PY $s.py > /dev/null
done
echo "done: results/ updated"
