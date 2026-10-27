#!/usr/bin/env python3
"""Why is the uidA target region not intact? Contig end (assembly) vs genuine truncation vs absence.

For each genome listed in --ids: blastn of the full-length uidA of the design strain against the assembly.
The alignments of the gene are merged (any contig) and each missing end of the gene is checked:
if the alignment stops within --end-tol nt of a contig end, the gene may simply continue on another contig
(assembly artefact: contig_end); otherwise the gene is genuinely interrupted or truncated (internal).
Genomes without a nucleotide hit are searched at protein level (tblastx) to separate divergent from absent.

Classes
  gene_complete             >= 95% of the gene aligned (target region interrupted/diverged elsewhere)
  contig_end                gene truncated only where an alignment reaches a contig end (possible assembly break)
  internal_truncation       at least one missing end of the gene is not at a contig end (genuine truncation/deletion)
  internal_gap              both ends of the gene present but a part in between is missing (deletion/insertion)
  divergent_protein_hit     no nucleotide hit, but tblastx finds uidA at protein level
  absent                    no nucleotide and no protein hit

Output (in --out): nonintact_check.tsv, summary.txt. Requires BLAST+ and primer_screen.py in the same directory.
"""
import argparse, csv, glob, os, subprocess, sys

from primer_screen import blast, read_fasta

END_TOL = 50
MIN_FULL = 0.95


def classify(hsps, contig_len, gene_len, end_tol=END_TOL):
    """hsps: dicts with qstart, qend, sseqid, sstart, send (subject coordinates, send < sstart = minus strand).
    contig_len: {contig: length}. Returns (class, covered_fraction, five_prime, three_prime)."""
    if not hsps:
        return "no_hit", 0.0, "", ""
    cov = set()
    for h in hsps: cov.update(range(h["qstart"], h["qend"] + 1))
    frac = len(cov) / gene_len
    if frac >= MIN_FULL:
        return "gene_complete", frac, "complete", "complete"
    q_lo, q_hi = min(h["qstart"] for h in hsps), max(h["qend"] for h in hsps)

    def end_status(h, side):
        """side 'lo': the gene's 5' end (qstart); 'hi': the 3' end. Is the alignment at a contig end there?"""
        plus = h["sstart"] <= h["send"]
        clen = contig_len[h["sseqid"]]
        if (side == "lo") == plus:   # alignment's left subject edge continues the gene's missing part
            dist = min(h["sstart"], h["send"]) - 1
        else:
            dist = clen - max(h["sstart"], h["send"])
        return "contig_end" if dist < end_tol else "internal"

    five = "complete" if q_lo <= 20 else end_status(min(hsps, key=lambda h: h["qstart"]), "lo")
    three = "complete" if q_hi >= gene_len - 20 else end_status(max(hsps, key=lambda h: h["qend"]), "hi")
    if five == "complete" and three == "complete":      # both ends present, but part of the gene is missing in between
        return "internal_gap", frac, five, three
    internal = "internal" in (five, three)
    return ("internal_truncation" if internal else "contig_end"), frac, five, three


def protein_hit(gene_fa, path):
    out = subprocess.run(["tblastx", "-query", gene_fa, "-subject", path, "-outfmt", "6 qstart qend pident length evalue",
                          "-evalue", "1e-20"], capture_output=True, text=True).stdout
    return bool(out.strip())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ids", required=True, help="file with genome IDs (one per line)")
    ap.add_argument("--genomes", required=True, help="directory searched recursively for <id><ext>")
    ap.add_argument("--ext", default=".fna")
    ap.add_argument("--gene", required=True, help="full-length uidA of the design strain (FASTA)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-protein", action="store_true", help="skip tblastx for genomes without a nucleotide hit")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    gene_len = sum(len(s) for s in read_fasta(a.gene).values())
    rows = []
    for gid in (l.strip() for l in open(a.ids) if l.strip()):
        found = glob.glob(os.path.join(a.genomes, "**", gid + a.ext), recursive=True)
        if not found:
            rows.append(dict(genome=gid, cls="file_not_found")); continue
        path = found[0]
        contigs = read_fasta(path)
        hsps = [h for h in blast(a.gene, path) if h["pident"] >= 80]
        cls, frac, five, three = classify(hsps, {k: len(v) for k, v in contigs.items()}, gene_len)
        if cls == "no_hit":
            cls = "absent" if a.no_protein or not protein_hit(a.gene, path) else "divergent_protein_hit"
        rows.append(dict(genome=gid, n_contigs=len(contigs), gene_aligned_pct=round(100 * frac, 1),
                         best_pident=max((h["pident"] for h in hsps), default=""), five_prime=five, three_prime=three,
                         cls=cls))
    cols = ["genome", "n_contigs", "gene_aligned_pct", "best_pident", "five_prime", "three_prime", "cls"]
    with open(os.path.join(a.out, "nonintact_check.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    counts = {}
    for r in rows: counts[r["cls"]] = counts.get(r["cls"], 0) + 1
    lines = [f"genomes checked: {len(rows)}"] + [f"  {k}\t{v}" for k, v in sorted(counts.items())]
    open(os.path.join(a.out, "summary.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
