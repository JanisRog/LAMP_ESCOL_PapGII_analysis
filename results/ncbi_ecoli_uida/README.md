# uidA binding sites in NCBI RefSeq complete E. coli genomes

Pipeline `scripts/ncbi_ecoli_uida.sbatch` (screen: `scripts/uida_ncbi_screen.py`), run 2026-10-27 on all 5,493 RefSeq
complete E. coli assemblies available at that date.

- assemblies.tsv – accession, organism name, assembly name and length of every genome analysed (NCBI datasets)
- genome_sites.tsv – per genome (file stem = accession_assemblyname_genomic): status (present = all eight binding sites covered
  by one alignment; incomplete = region partly present; absent = no alignment), best identity of the amplicon
  alignment, number of sites covered, and mismatches to the primer per site (empty = not covered)
- summary.txt – counts per status with exact 95% CI and the mismatch distribution per site
- phylogroups.tsv – Clermont phylogroup per genome (EzClermont 1.0.0; manifest: phylogroups_manifest.txt)
- uidA_amplicon.fasta – primer-derived query sequence
