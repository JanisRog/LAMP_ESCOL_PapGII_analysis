#!/usr/bin/env python3
"""uidA LAMP binding sites in a large genome collection (e.g. NCBI RefSeq complete E. coli genomes).

For every genome: blastn of the uidA amplicon (primer-derived sequence, as in primer_screen.py) against the
assembly; for each of the eight binding sites the mismatches to the primer are counted if a single alignment
covers the site. Status per genome: present = all eight sites covered by one alignment; incomplete = amplicon
region only partly present (truncation, rearrangement or contig break); absent = no alignment.

Output (in --out): genome_sites.tsv (per genome: status, best identity, mismatches per site), summary.txt,
manifest.txt, README.md. Requires BLAST+ (blastn) and primer_screen.py in the same directory.
"""
import argparse, csv, datetime, glob, os, subprocess
from concurrent.futures import ThreadPoolExecutor

from primer_screen import blast, binding_sites, read_fasta, reference_amplicon, site_mismatches, THREE_PRIME_SITES

SITE_ORDER = ["F3", "F2", "LF", "F1c", "B1c", "LB", "B2", "B3"]


def screen(args):
    gid, path, amp_fa, sites = args
    hits = sorted(blast(amp_fa, path), key=lambda h: -h["bitscore"])
    row = {"genome": gid, "best_pident": hits[0]["pident"] if hits else ""}
    covered = 0
    for site in SITE_ORDER:
        _, s, e, o, _ = sites[site]
        res = next((r for r in (site_mismatches(h, s, e, o, site in THREE_PRIME_SITES) for h in hits) if r), None)
        row[site] = "" if res is None else res[0]
        covered += res is not None
    row["sites_covered"] = covered
    row["status"] = "absent" if not hits else "present" if covered == len(SITE_ORDER) else "incomplete"
    return row


def exact_ci(x, n, alpha=0.05):
    """Clopper-Pearson interval, binomial terms in log space (safe for large n)."""
    from math import lgamma, log, exp
    def term(i, p):
        if p <= 0: return 1.0 if i == 0 else 0.0
        if p >= 1: return 1.0 if i == n else 0.0
        return exp(lgamma(n + 1) - lgamma(i + 1) - lgamma(n - i + 1) + i * log(p) + (n - i) * log(1 - p))
    def cdf(k, p): return sum(term(i, p) for i in range(k + 1))
    def bis(f):
        lo, hi = 0.0, 1.0
        for _ in range(60):
            m = (lo + hi) / 2
            lo, hi = (m, hi) if f(m) else (lo, m)
        return lo
    lo = 0.0 if x == 0 else bis(lambda p: 1 - cdf(x - 1, p) < alpha / 2)
    hi = 1.0 if x == n else bis(lambda p: cdf(x, p) > alpha / 2)
    return 100 * lo, 100 * hi


def write_summary(rows, out):
    n = len(rows)
    lines = [f"uidA amplicon region in {n} genomes"]
    for st in ("present", "incomplete", "absent"):
        k = sum(r["status"] == st for r in rows)
        lo, hi = exact_ci(k, n)
        lines.append(f"  {st}\t{k}\t{100 * k / n:.1f}%\t(95% CI {lo:.1f}-{hi:.1f}%)")
    pres = [r for r in rows if r["status"] == "present"]
    lines.append(f"binding sites in genomes with the region present (n = {len(pres)}): genomes with 0 / 1 / 2 / >=3 mismatches")
    for site in SITE_ORDER:
        v = [int(r[site]) for r in pres]
        lines.append(f"  {site}\t" + "\t".join(str(sum(f(x) for x in v)) for f in
                     (lambda x: x == 0, lambda x: x == 1, lambda x: x == 2, lambda x: x >= 3)))
    open(os.path.join(out, "summary.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--pattern", default="*.fna")
    ap.add_argument("--design-genome", required=True, help="genome carrying the uidA used for primer design (700223-20)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--summary-only", action="store_true",
                    help="rebuild summary.txt from an existing genome_sites.tsv in --out (no BLAST)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.summary_only:
        with open(os.path.join(a.out, "genome_sites.tsv")) as fh:
            write_summary(list(csv.DictReader(fh, delimiter="\t")), a.out)
        return

    amp, offset = reference_amplicon("uidA", read_fasta(a.design_genome))
    sites, primer_amp = binding_sites("uidA", amp)
    amp_fa = os.path.join(a.out, "uidA_amplicon.fasta")
    with open(amp_fa, "w") as f: f.write(f">uidA_amplicon\n{primer_amp}\n")

    genomes = sorted(glob.glob(os.path.join(a.genomes, "**", a.pattern), recursive=True))
    jobs = [(os.path.splitext(os.path.basename(g))[0], g, amp_fa, sites) for g in genomes]
    with ThreadPoolExecutor(a.threads) as ex:
        rows = list(ex.map(screen, jobs))
    cols = ["genome", "status", "best_pident", "sites_covered"] + SITE_ORDER
    with open(os.path.join(a.out, "genome_sites.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    n = len(rows)
    write_summary(rows, a.out)

    blast_version = subprocess.run(["blastn", "-version"], capture_output=True, text=True).stdout.splitlines()[0]
    with open(os.path.join(a.out, "manifest.txt"), "w") as fh:
        fh.write(f"date: {datetime.datetime.now().isoformat(timespec='seconds')}\nscript: {os.path.abspath(__file__)}\n"
                 f"genomes: {a.genomes} ({n} files matching {a.pattern})\ndesign_genome: {a.design_genome}\n{blast_version}\n"
                 "per genome: best alignments of the primer-derived uidA amplicon (with 100 nt flanks); site covered if one\n"
                 "alignment spans it; mismatches counted against the primer sequence\n")
    with open(os.path.join(a.out, "README.md"), "w") as fh:
        fh.write("""# uidA binding sites in a large genome collection

Generated by uida_ncbi_screen.py (inputs and versions in manifest.txt).

- genome_sites.tsv - per genome: status (present = all eight binding sites covered by one alignment; incomplete =
  region only partly present; absent = no alignment), best identity of the amplicon alignment, number of sites
  covered, and the mismatches to the primer at each site (empty = site not covered)
- summary.txt - counts per status with exact 95% CI, and the mismatch distribution per site
""")


if __name__ == "__main__":
    main()
