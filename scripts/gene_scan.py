#!/usr/bin/env python3
"""Per-position variability of a target gene across a genome collection, with the LAMP amplicon located
in the gene. Answers whether the primer region is more or less variable than the rest of the gene.

1. Reference gene: best full-length hit of --gene-query in the primer design genome (E. coli 700223-20).
2. LAMP amplicon and binding sites located in the reference gene (same rules as primer_screen.py).
3. For every genome: blastn the reference gene against the assembly; along the best HSP, record per reference
   position whether the genome matches, mismatches (substitution) or has a gap; positions outside the HSP
   are 'not covered'. Optionally restricted to genomes listed in --genome-list.

Outputs (in --out): reference_gene.fasta, gene_sites.tsv (amplicon and binding sites in gene coordinates),
position_variability.tsv, genome_variants.tsv, genome_coverage.tsv, manifest.txt, README.md
Requires BLAST+ (blastn) and primer_screen.py in the same directory.
"""
import argparse, datetime, glob, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

from primer_screen import PRIMERS, ROLE, binding_sites, blast, rc, read_fasta, reference_amplicon


def reference_gene(query_fa, design_genome, out_fa):
    hits = sorted(blast(query_fa, design_genome), key=lambda h: -h["bitscore"])
    if not hits: sys.exit("ERROR: gene query not found in design genome")
    h = hits[0]
    contig = read_fasta(design_genome)[h["sseqid"]]
    lo, hi = min(h["sstart"], h["send"]), max(h["sstart"], h["send"])
    seq = contig[lo - 1:hi]
    if h["sstart"] > h["send"]: seq = rc(seq)
    with open(out_fa, "w") as f: f.write(f">reference_gene {h['sseqid']}:{lo}-{hi}\n{seq}\n")
    return seq, h


def write_gene_sites(target, gene, path):
    """Amplicon and binding sites in gene coordinates (1-based). The amplicon may lie on either strand of the
    gene; strand '-' means the primers' F3->B3 direction runs opposite to the gene."""
    L = len(gene)
    amp, _ = reference_amplicon(target, {"gene": gene})
    sites, _ = binding_sites(target, amp)
    j, strand = gene.find(amp), "+"
    if j < 0:
        j, strand = rc(gene).find(amp), "-"
    if j < 0: sys.exit(f"ERROR: {target} amplicon not found in reference gene")
    def to_gene(s, e):  # half-open interval [s, e) in amplicon -> 1-based closed interval in gene
        return (j + s + 1, j + e) if strand == "+" else (L - (j + e) + 1, L - (j + s))
    with open(path, "w") as fh:
        fh.write("feature\trole\tgene_start\tgene_end\tstrand\n")
        fh.write("amplicon\t-\t%d\t%d\t%s\n" % (*sorted(to_gene(sites["F3"][1], sites["B3"][2])), strand))
        for name, (_, s, e, o, _) in sites.items():
            fh.write(f"{name}\t{ROLE[name]}\t" + "\t".join(map(str, to_gene(s, e))) + f"\t{strand}\n")


def scan(args):
    gid, path, ref_fa, L = args
    hits = sorted(blast(ref_fa, path), key=lambda h: -h["bitscore"])
    if not hits: return gid, None, []
    h = hits[0]
    variants, qpos = [], h["qstart"] - 1
    for qc, sc in zip(h["qseq"], h["sseq"]):
        if qc == "-":                       # insertion in genome relative to reference
            variants.append((qpos, "ins", sc)); continue
        qpos += 1
        if sc == "-": variants.append((qpos, "del", "-"))
        elif sc != qc: variants.append((qpos, "sub", sc))
    return gid, (h["qstart"], h["qend"], h["pident"]), variants


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", required=True, choices=list(PRIMERS))
    ap.add_argument("--gene-query", help="FASTA of the gene/allele used to find the reference gene")
    ap.add_argument("--design-genome")
    ap.add_argument("--genomes")
    ap.add_argument("--sites-only", action="store_true",
                    help="only (re)write gene_sites.tsv from an existing OUT/reference_gene.fasta (no BLAST)")
    ap.add_argument("--pattern", default="*.fasta")
    ap.add_argument("--genome-list", help="optional file with genome IDs to include (one per line)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    ref_fa = os.path.join(a.out, "reference_gene.fasta")
    if a.sites_only:
        write_gene_sites(a.target, next(iter(read_fasta(ref_fa).values())), os.path.join(a.out, "gene_sites.tsv"))
        print(open(os.path.join(a.out, "gene_sites.tsv")).read()); return
    if not (a.gene_query and a.design_genome and a.genomes):
        sys.exit("ERROR: --gene-query, --design-genome and --genomes are required unless --sites-only")
    gene, ref_hit = reference_gene(a.gene_query, a.design_genome, ref_fa)
    L = len(gene)
    write_gene_sites(a.target, gene, os.path.join(a.out, "gene_sites.tsv"))

    genomes = sorted(glob.glob(os.path.join(a.genomes, "**", a.pattern), recursive=True))
    if a.genome_list:
        keep = {l.strip() for l in open(a.genome_list) if l.strip()}
        genomes = [g for g in genomes if os.path.splitext(os.path.basename(g))[0] in keep]
    if not genomes: sys.exit("ERROR: no genomes selected")
    jobs = [(os.path.splitext(os.path.basename(g))[0], g, ref_fa, L) for g in genomes]
    with ThreadPoolExecutor(a.threads) as ex:
        results = list(ex.map(scan, jobs))

    covered, sub, gap = [0] * (L + 1), [0] * (L + 1), [0] * (L + 1)
    with open(os.path.join(a.out, "genome_variants.tsv"), "w") as fv, \
         open(os.path.join(a.out, "genome_coverage.tsv"), "w") as fc:
        fv.write("genome\tposition\ttype\talt\n"); fc.write("genome\tqstart\tqend\tpident\n")
        for gid, hsp, variants in results:
            if hsp is None:
                fc.write(f"{gid}\t\t\t\n"); continue
            fc.write(f"{gid}\t{hsp[0]}\t{hsp[1]}\t{hsp[2]}\n")
            for p in range(hsp[0], hsp[1] + 1): covered[p] += 1
            seen_sub, seen_gap = set(), set()
            for pos, typ, alt in variants:
                fv.write(f"{gid}\t{pos}\t{typ}\t{alt}\n")
                if typ == "sub": seen_sub.add(pos)
                else: seen_gap.add(max(pos, 1))
            for p in seen_sub: sub[p] += 1
            for p in seen_gap: gap[p] += 1
    with open(os.path.join(a.out, "position_variability.tsv"), "w") as fh:
        fh.write("position\tref_base\tn_covered\tn_substitution\tn_gap\n")
        for p in range(1, L + 1):
            fh.write(f"{p}\t{gene[p - 1]}\t{covered[p]}\t{sub[p]}\t{gap[p]}\n")

    blast_version = subprocess.run(["blastn", "-version"], capture_output=True, text=True).stdout.splitlines()[0]
    with open(os.path.join(a.out, "manifest.txt"), "w") as fh:
        fh.write(f"date: {datetime.datetime.now().isoformat(timespec='seconds')}\nscript: {os.path.abspath(__file__)}\n"
                 f"target: {a.target}\ngene_query: {a.gene_query}\ndesign_genome: {a.design_genome}\n"
                 f"reference_gene: {ref_hit['sseqid']} ({L} bp, {ref_hit['pident']}% identity to query)\n"
                 f"genomes: {a.genomes} ({len(genomes)} files matching {a.pattern}"
                 f"{', restricted to ' + a.genome_list if a.genome_list else ''})\n{blast_version}\n"
                 f"per genome: best HSP of the reference gene (by bitscore)\n")
    with open(os.path.join(a.out, "README.md"), "w") as fh:
        fh.write(f"""# Gene-wide variability – {a.target}

Generated by gene_scan.py; inputs and versions in manifest.txt.

- reference_gene.fasta – {a.target} gene of the primer design strain (reference coordinates = positions 1..{L})
- gene_sites.tsv – LAMP amplicon (F3 to B3) and the eight binding sites in gene coordinates (1-based);
  strand '-' = primers run F3->B3 opposite to the gene orientation
- position_variability.tsv – per reference position: genomes whose best HSP covers it (n_covered), with a
  substitution (n_substitution) or an indel (n_gap; insertions are assigned to the preceding position)
- genome_variants.tsv – every substitution/indel per genome (position in reference gene, type, alternative base)
- genome_coverage.tsv – per genome: reference positions covered by the best HSP and its identity (empty = no hit)
""")
    print(f"done: {a.out} ({len(genomes)} genomes, gene {L} bp)")


if __name__ == "__main__":
    main()
