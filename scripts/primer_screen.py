#!/usr/bin/env python3
"""In silico screen of the uidA and papGII LAMP primer sets against a genome collection.

Steps
1. Reference amplicons: locate F3 and B3 of each primer set in a reference gene sequence (<= 2 mismatches)
   and extract the F3..B3 region in F3-forward orientation. The primer-derived sequence is then written
   into every binding site, so genomes are compared with the primers rather than with the reference isolate.
2. Binding sites: split FIP into F1c/F2 and BIP into B1c/B2 by matching against the amplicon (fewest mismatches);
   all eight sites (F3, F2, LF, F1c, B1c, LB, B2, B3) get amplicon coordinates and an orientation.
3. For every genome: blastn the amplicon against the assembly; for each binding site, take the
   best HSP that covers it fully and count mismatches/gaps, separately in the 3'-terminal 5 nt
   of primers with a free 3' end (F3, F2, LF, LB, B2, B3). Genomes are classified by primer role:
   inner sites (F2/F1c/B1c/B2) are essential, outer (F3/B3) and loop (LF/LB) primers are tolerant.
4. papG typing of the same genomes, per locus: HSPs of the reference alleles are grouped by genome position
   (so a gene split by an indel counts once), the best-scoring allele is called if its HSPs (>= 90% identity)
   cover >= 80% of the allele. Used to report papGII coverage and papGI/papGIII cross-reactivity.

Only the Python standard library and BLAST+ (blastn on PATH) are required.
"""
import argparse, csv, datetime, glob, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

PRIMERS = {  # Table 1 of the manuscript
    "uidA": {
        "F3": "CCTGATCAACCACAAACCGT", "B3": "CGACAGCTGCAGTTTCATCA",
        "FIP": "TCGTGCACCATCAGCACGTTATACTTTACTGGCTTTGGTCGT",
        "BIP": "AATGGACTGGATTGGGGCCAATGTTCATCTGCCCAGTCGA",
        "LF": "CTTTGCCGCGCAAATCCG", "LB": "CGCATTACCCTTACGCTGAAGAG",
    },
    "papGII": {
        "F3": "ACACAGAAAGAGTCTGAGC", "B3": "TGTTACAATTCCATACACTTCC",
        "FIP": "TTAAGAATATCGGCGGATGCCGGCTATTAATAGACAGATCACCAT",
        "BIP": "AACATTTCATTCTCTCTTGGGAGAGCATTTCGCGAGTTACTTGG",
        "LF": "TCTGCACAGTCTCTGGAAATAAAGC", "LB": "ATTGTATGGGATTTTAAAACGGGCC",
    },
}
# Length of the 5' part of FIP (F1c) and BIP (B1c). uidA: NEB LAMP Primer Design Tool output (the junction
# base fits both halves, so it cannot be inferred); papGII: the only split matching the PapG-II reference.
SPLIT = {"uidA": {"FIP": 22, "BIP": 21}, "papGII": {"FIP": 22, "BIP": 25}}
THREE_PRIME_SITES = {"F3", "F2", "LF", "LB", "B2", "B3"}  # sites carrying a primer's free 3' end
ROLE = {"F3": "outer", "B3": "outer", "LF": "loop", "LB": "loop",
        "F2": "inner", "F1c": "inner", "B1c": "inner", "B2": "inner"}
INITIATING = {"F2", "B2"}  # 3' ends of FIP/BIP: extension from these sites is essential for amplification
END_NT = 5
FLANK = 100  # reference flank added on both sides of the amplicon so that F3/B3 lie inside the BLAST alignment
BLAST_FMT = "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qseq sseq"
LOCUS_GAP = 200  # HSPs on the same contig closer than this belong to one papG locus (indel-split gene),
                 # while neighbouring papG copies (>= ~600 nt apart in tandem pap operons) stay separate

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(s): return s.translate(COMP)[::-1]

def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name: seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name: seqs[name] = "".join(buf).upper()
    return seqs

MAX_REF_MM = 2  # mismatches tolerated when placing primers on the reference sequence

def hamming_find(seq, text, max_mm):
    """Best (mismatches, position) of seq in text with <= max_mm mismatches, or None."""
    best, n = None, len(seq)
    for i in range(len(text) - n + 1):
        mm = 0
        for a, b in zip(seq, text[i:i + n]):
            if a != b:
                mm += 1
                if mm > max_mm: break
        if mm <= max_mm and (best is None or mm < best[0]):
            best = (mm, i)
            if mm == 0: break
    return best

def reference_amplicon(target, ref):
    """F3..B3 region of the reference (F3-forward orientation) plus FLANK nt on both sides, allowing
    MAX_REF_MM mismatches per primer. Returns (sequence, offset of F3 within sequence)."""
    def pad(strand, i, j):
        lo = max(0, i - FLANK)
        return strand[lo:j + FLANK], i - lo
    p = PRIMERS[target]
    for contig in ref.values():  # fast path: exact F3 and B3
        for strand in (contig, rc(contig)):
            i = strand.find(p["F3"])
            j = strand.find(rc(p["B3"]), i, i + 600) if i >= 0 else -1
            if j >= 0: return pad(strand, i, j + len(p["B3"]))
    for contig in ref.values():
        for strand in (contig, rc(contig)):
            f = hamming_find(p["F3"], strand, MAX_REF_MM)
            if not f: continue
            b = hamming_find(rc(p["B3"]), strand[f[1]:f[1] + 600], MAX_REF_MM)
            if b: return pad(strand, f[1], f[1] + b[1] + len(p["B3"]))
    sys.exit(f"ERROR: F3/B3 of {target} not found (<= {MAX_REF_MM} mismatches) in the reference sequence")

def locate(seq, amp):
    """(start, end, orientation, mismatches) of seq in amp; orientation '+' = primer has sense sequence."""
    hits = [(h[0], h[1], o) for o, q in (("+", seq), ("-", rc(seq))) for h in [hamming_find(q, amp, MAX_REF_MM)] if h]
    if not hits: return None
    mm, i, o = min(hits)
    return i, i + len(seq), o, mm

def binding_sites(target, amp):
    """Eight binding sites with amplicon coordinates; returns (sites, primer_amplicon).
    primer_amplicon = reference amplicon with the primer-derived sequence written into every site,
    so that genomes are compared with the primers, not with the reference isolate."""
    p, sites = PRIMERS[target], {}
    for name in ("F3", "B3", "LF", "LB"):
        loc = locate(p[name], amp)
        if not loc: sys.exit(f"ERROR: {target} {name} not found in reference amplicon")
        sites[name] = (p[name],) + loc
    for comp, (five, three) in (("FIP", ("F1c", "F2")), ("BIP", ("B1c", "B2"))):
        seq, best = p[comp], None
        splits = [SPLIT[target][comp]] if comp in SPLIT.get(target, {}) else range(15, len(seq) - 14)
        for k in splits:
            a, b = locate(seq[:k], amp), locate(seq[k:], amp)
            if a and b and (best is None or a[3] + b[3] < best[0]):
                best = (a[3] + b[3], k, a, b)
        if not best: sys.exit(f"ERROR: could not split {target} {comp} against the reference amplicon")
        _, k, a, b = best
        sites[five] = (seq[:k],) + a
        sites[three] = (seq[k:],) + b
    sites = dict(sorted(sites.items(), key=lambda kv: kv[1][1]))
    ivals = sorted((v[1], v[2], k) for k, v in sites.items())
    for (s1, e1, k1), (s2, e2, k2) in zip(ivals, ivals[1:]):
        if s2 < e1: sys.exit(f"ERROR: {target} sites {k1} and {k2} overlap")
    seq = list(amp)
    for name, (prim, s, e, o, _) in sites.items():
        seq[s:e] = prim if o == "+" else rc(prim)
    return sites, "".join(seq)

def blast(query, subject):
    out = subprocess.run(["blastn", "-query", query, "-subject", subject, "-outfmt", BLAST_FMT,
                          "-evalue", "1e-10", "-dust", "no"], capture_output=True, text=True, check=True).stdout
    hits = []
    for line in out.splitlines():
        f = line.split("\t")
        hits.append(dict(qseqid=f[0], sseqid=f[1], pident=float(f[2]), length=int(f[3]), qlen=int(f[4]),
                         qstart=int(f[5]), qend=int(f[6]), sstart=int(f[7]), send=int(f[8]),
                         bitscore=float(f[10]), qseq=f[11], sseq=f[12]))
    return hits

def site_mismatches(hit, start, end, orient, three_prime):
    """Mismatches/gaps of one HSP inside amplicon interval [start, end) (0-based). None if not fully covered."""
    if hit["qstart"] - 1 > start or hit["qend"] < end:
        return None
    qpos, mm, mm3 = hit["qstart"] - 2, 0, 0
    end_zone = range(end - END_NT, end) if orient == "+" else range(start, start + END_NT)
    for qc, sc in zip(hit["qseq"], hit["sseq"]):
        if qc != "-": qpos += 1
        if start <= qpos < end and qc != sc:
            mm += 1
            if three_prime and qpos in end_zone: mm3 += 1
    return mm, mm3

def type_papg(hits, min_id, min_cov):
    """papG allele per locus. HSPs (>= min_id identity) are grouped into loci by contig position, so a gene
    split into several HSPs by an indel is assessed as one; at each locus the allele with the highest summed
    bitscore is called if its HSPs together cover >= min_cov % of the reference allele."""
    calls, loci_out, by_contig = [], [], {}
    for h in hits:
        if h["pident"] >= min_id:
            by_contig.setdefault(h["sseqid"], []).append(
                dict(h, lo=min(h["sstart"], h["send"]), hi=max(h["sstart"], h["send"])))
    for hs in by_contig.values():
        hs.sort(key=lambda h: h["lo"])
        loci, cur = [], [hs[0]]
        for h in hs[1:]:
            if h["lo"] > max(x["hi"] for x in cur) + LOCUS_GAP:
                loci.append(cur); cur = [h]
            else:
                cur.append(h)
        loci.append(cur)
        for locus in loci:
            score, qlen, ivals = {}, {}, {}
            for h in locus:
                score[h["qseqid"]] = score.get(h["qseqid"], 0) + h["bitscore"]
                qlen[h["qseqid"]] = h["qlen"]
                ivals.setdefault(h["qseqid"], []).append((h["qstart"], h["qend"]))
            best = max(score, key=score.get)
            covered = len(set().union(*(range(a, b + 1) for a, b in ivals[best])))
            cov_pct = 100 * covered / qlen[best]
            if cov_pct >= min_cov:
                calls.append(best)
            loci_out.append(dict(contig=locus[0]["sseqid"], start=min(h["lo"] for h in locus),
                                 end=max(h["hi"] for h in locus), best_allele=best, coverage=round(cov_pct, 1),
                                 n_hsps=len(locus), called=int(cov_pct >= min_cov)))
    return calls, loci_out

def screen_genome(args):
    gid, path, amp_files, sites, papg_ref, min_id, min_cov = args
    rows, summary = [], {"genome": gid}
    for target, amp_fa in amp_files.items():
        hits = sorted(blast(amp_fa, path), key=lambda h: -h["bitscore"])
        n_cov, mm_inner, mm_other, critical, uncovered_inner, uncovered_other = 0, 0, 0, False, False, False
        for site, (seq, s, e, o, _) in sites[target].items():
            res = next((r for r in (site_mismatches(h, s, e, o, site in THREE_PRIME_SITES) for h in hits) if r), None)
            inner = ROLE[site] == "inner"
            if res is None:
                uncovered_inner |= inner; uncovered_other |= not inner
            else:
                n_cov += 1
                if inner: mm_inner += res[0]
                else: mm_other += res[0]
                critical |= (site in INITIATING and res[1] > 0) or (inner and res[0] >= 2)
            rows.append(dict(genome=gid, target=target, site=site, role=ROLE[site], covered=int(res is not None),
                             mismatches="" if res is None else res[0], mismatches_3prime="" if res is None else res[1]))
        status = ("absent" if not hits else "incomplete" if uncovered_inner else "critical" if critical
                  else "minor" if (mm_inner or mm_other or uncovered_other) else "intact")
        summary.update({f"{target}_status": status, f"{target}_sites_covered": n_cov,
                        f"{target}_mismatches_inner": mm_inner, f"{target}_mismatches_outer_loop": mm_other,
                        f"{target}_best_pident": hits[0]["pident"] if hits else ""})
    alleles = {}
    calls, loci = type_papg(blast(papg_ref, path), min_id, min_cov)
    for allele in calls:
        alleles[allele] = alleles.get(allele, 0) + 1
    summary["papG_alleles"] = ";".join(f"{a}({n})" for a, n in sorted(alleles.items())) or "none"
    return summary, rows, [dict(genome=gid, **l) for l in loci]

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", required=True, help="directory searched recursively for assemblies")
    ap.add_argument("--pattern", default="*.fasta", help="assembly file pattern (default *.fasta)")
    ap.add_argument("--ref-uidA", required=True, help="reference uidA sequence (FASTA) for placing the primers")
    ap.add_argument("--ref-papGII", required=True, help="reference papGII sequence (FASTA) for placing the primers")
    ap.add_argument("--papg-ref", required=True, help="PapG reference alleles (FASTA)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-id", type=float, default=90)
    ap.add_argument("--min-cov", type=float, default=80)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    refs = {"uidA": a.ref_uidA, "papGII": a.ref_papGII}
    amp_files, sites = {}, {}
    with open(os.path.join(a.out, "binding_sites.tsv"), "w") as fh:
        fh.write("target\tsite\trole\tstart\tend\torientation\tsequence\tthree_prime_end_scored\treference_mismatches\n")
        for target in PRIMERS:
            amp, offset = reference_amplicon(target, read_fasta(refs[target]))
            sites[target], primer_amp = binding_sites(target, amp)
            amp_files[target] = os.path.join(a.out, f"{target}_amplicon.fasta")
            for site, (seq, s, e, o, mm) in sites[target].items():
                fh.write(f"{target}\t{site}\t{ROLE[site]}\t{s - offset + 1}\t{e - offset}\t{o}\t{seq}"
                         f"\t{int(site in THREE_PRIME_SITES)}\t{mm}\n")
            with open(amp_files[target], "w") as f:  # query file for blastn: primer amplicon only
                f.write(f">{target}_amplicon\n{primer_amp}\n")
            with open(os.path.join(a.out, f"{target}_reference_isolate_amplicon.fasta"), "w") as f:
                f.write(f">{target}_reference_isolate\n{amp}\n")

    genomes = sorted(glob.glob(os.path.join(a.genomes, "**", a.pattern), recursive=True))
    if not genomes: sys.exit(f"ERROR: no {a.pattern} files under {a.genomes}")
    jobs = [(os.path.splitext(os.path.basename(g))[0], g, amp_files, sites, a.papg_ref, a.min_id, a.min_cov) for g in genomes]
    with ThreadPoolExecutor(a.threads) as ex:
        results = list(ex.map(screen_genome, jobs))

    summaries = [s for s, _, _ in results]
    with open(os.path.join(a.out, "papG_loci.tsv"), "w", newline="") as fh:
        cols = ["genome", "contig", "start", "end", "best_allele", "coverage", "n_hsps", "called"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader()
        for _, _, loci in results: w.writerows(loci)
    with open(os.path.join(a.out, "genome_summary.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summaries[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(summaries)
    with open(os.path.join(a.out, "site_mismatches.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0][1][0]), delimiter="\t", lineterminator="\n"); w.writeheader()
        for _, rows, _ in results: w.writerows(rows)

    # counts: target status overall, and papGII status per papG allele group
    lines, n = [], len(summaries)
    for target in PRIMERS:
        lines.append(f"{target} (all {n} genomes)")
        for st in ("intact", "minor", "critical", "incomplete", "absent"):
            k = sum(s[f"{target}_status"] == st for s in summaries)
            lines.append(f"  {st}\t{k}\t{100 * k / n:.1f}%")
    groups = {}
    for s in summaries:
        g = "+".join(sorted({x.split("(")[0] for x in s["papG_alleles"].split(";")}))
        groups.setdefault(g, []).append(s)
    lines.append("papGII primer set by papG allele group (status counts)")
    for g, members in sorted(groups.items()):
        counts = {}
        for s in members: counts[s["papGII_status"]] = counts.get(s["papGII_status"], 0) + 1
        lines.append(f"  {g} (n={len(members)}): " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    lines.append("per-site (genomes covering the site): number with 1 / 2 / >=3 mismatches, and with a 3'-end mismatch")
    for target in PRIMERS:
        for site in sites[target]:
            cov = [r for _, rows, _ in results for r in rows if r["target"] == target and r["site"] == site and r["covered"]]
            lines.append(f"  {target}\t{site}\t{ROLE[site]}\tcovered={len(cov)}"
                         f"\tmm1={sum(r['mismatches'] == 1 for r in cov)}\tmm2={sum(r['mismatches'] == 2 for r in cov)}"
                         f"\tmm3plus={sum(r['mismatches'] >= 3 for r in cov)}\tend3={sum(r['mismatches_3prime'] > 0 for r in cov)}")
    open(os.path.join(a.out, "summary.txt"), "w").write("\n".join(lines) + "\n")

    blast_version = subprocess.run(["blastn", "-version"], capture_output=True, text=True).stdout.splitlines()[0]
    with open(os.path.join(a.out, "manifest.txt"), "w") as fh:
        fh.write(f"date: {datetime.datetime.now().isoformat(timespec='seconds')}\n"
                 f"script: {os.path.abspath(__file__)}\ngenomes: {a.genomes} ({n} files matching {a.pattern})\n"
                 f"reference_uidA: {a.ref_uidA}\nreference_papGII: {a.ref_papGII}\npapG_reference: {a.papg_ref}\n{blast_version}\n"
                 f"blast_options: -outfmt '{BLAST_FMT}' -evalue 1e-10 -dust no\n"
                 f"papG_call: per locus (HSPs within {LOCUS_GAP} nt merged), best allele by summed bitscore, pident >= {a.min_id}, combined coverage >= {a.min_cov} %\n"
                 f"three_prime_zone: last {END_NT} nt of {sorted(THREE_PRIME_SITES)}\nflank: {FLANK} nt reference flank on both sides of the amplicon query\n"
                 f"status rules: incomplete = inner site (F2/F1c/B1c/B2) not covered; critical = 3'-end mismatch in F2/B2 or >= 2 mismatches in an inner site; "
                 f"minor = any other mismatch or uncovered outer/loop site\n")
    with open(os.path.join(a.out, "README.md"), "w") as fh:
        fh.write(f"""# In silico LAMP primer screen

Generated by primer_screen.py; paths, versions and thresholds in manifest.txt.

- uidA_amplicon.fasta / papGII_amplicon.fasta – F3..B3 region of the reference sequence (+ {FLANK} nt flanks) with the
  primer-derived sequence written into every binding site (BLAST query; genomes are compared with the primers)
- *_reference_isolate_amplicon.fasta – the unmodified region of the reference isolate
- binding_sites.tsv – the eight binding sites per primer set (FIP = F1c+F2, BIP = B1c+B2) with role
  (outer F3/B3, loop LF/LB, inner F2/F1c/B1c/B2), coordinates (1-based, F3 start = 1); orientation '+' = primer has the sense sequence of the amplicon;
  three_prime_end_scored = 1 for sites carrying a primer's free 3' end;
  reference_mismatches = mismatches between primer and reference isolate at that site
- site_mismatches.tsv – per genome, target and site (with role): covered (1 = a single HSP spans the site),
  mismatches (mismatches + gap positions in the site), mismatches_3prime (within the last {END_NT} nt at the 3' end)
- genome_summary.tsv – per genome: status per target, mismatches in inner vs outer/loop sites, best identity,
  papG alleles with copies (per locus; see manifest). Status: intact = all sites covered, 0 mismatches; minor = mismatches only in
  outer/loop sites or single internal mismatches in inner sites (or an outer/loop site not covered);
  critical = 3'-end mismatch in F2 or B2, or >= 2 mismatches in one inner site; incomplete = an inner site not
  covered (truncation or contig break); absent = no hit. The status rules are a heuristic, not a measured
  amplification outcome.
- papG_loci.tsv – every papG locus (HSPs >= 90% identity within {LOCUS_GAP} nt merged): position, best allele,
  combined coverage of that allele, number of HSPs, called (1 = coverage >= threshold)
- summary.txt – counts and percentages per status, papGII status per papG allele group, per-site mismatch distribution
""")
    print(open(os.path.join(a.out, "summary.txt")).read())

if __name__ == "__main__":
    main()
