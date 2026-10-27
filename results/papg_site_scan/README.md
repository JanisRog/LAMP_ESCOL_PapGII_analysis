# papG typing and papGII binding-site mismatches – UPEC collection (738 genomes)

`locus_site_mismatches.tsv` (`scripts/papg_site_scan.py`): one row per papG locus and papGII binding site – genome, locus
number, allele (PapG-I to PapG-IV, best match; >= 90% identity over >= 80% of the reference), position (contig:start-end),
site, and the fewest mismatches of the site sequence at any ungapped position on either strand within the locus
(+ 50 nt). Reference alleles: `data/papG_reference_alleles.fasta`. Inputs and versions: manifest.txt.
