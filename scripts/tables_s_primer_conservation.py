"""Supplementary tables: descriptive in silico conservation of the LAMP binding sites in 738 UPEC genomes.

Per target and binding site: genomes with 0, 1, 2 or >= 3 mismatches to the primer and with a mismatch in the
3'-terminal 5 nt of primers with a free 3' end. Descriptive only: the effect of mismatches on amplification was
not tested experimentally.
Genomes: uidA - all genomes in which the amplicon region is present; papGII - genomes typed papGII.
Presence of binding sites per papG group: coverage from results/primer_screen/upec, groups from the per-locus
typing against PapG-I to PapG-IV (results/papg_site_scan).

Output: results/figures/tableS_binding_site_mismatches.tsv, tableS_binding_site_presence.tsv
"""
import pandas as pd

from lamp_common import RESULTS
from primer_screen import ROLE, THREE_PRIME_SITES, binding_sites, rc, read_fasta, reference_amplicon

SCAN, SCREEN, OUT = RESULTS / "gene_scan", RESULTS / "primer_screen" / "upec", RESULTS / "figures"
END_NT = 5


def load(target):
    """Binding sites (amplicon coordinates) and per-genome mismatches to the primer at each site."""
    d = SCAN / target
    gene = next(iter(read_fasta(d / "reference_gene.fasta").values()))
    L = len(gene)
    amp, offset = reference_amplicon(target, {"gene": gene})
    sites, primer_amp = binding_sites(target, amp)
    j, strand = gene.find(amp), "+"
    if j < 0: j, strand = rc(gene).find(amp), "-"
    a0 = j + offset
    to_amp = (lambda g: g - a0) if strand == "+" else (lambda g: (L - g + 1) - a0)
    st = pd.DataFrame([dict(site=k, role=ROLE[k], start=s - offset + 1, end=e - offset, orientation=o)
                       for k, (seq, s, e, o, _) in sites.items()]).sort_values("start")
    amp_len = int(st.end.max())
    gpos = [a0 + 1, a0 + amp_len] if strand == "+" else [L - a0 - amp_len + 1, L - a0]
    gene_lo, gene_hi = min(gpos), max(gpos)
    cov = pd.read_csv(d / "genome_coverage.tsv", sep="\t")
    var = pd.read_csv(d / "genome_variants.tsv", sep="\t", keep_default_na=False)
    var["apos"] = var.position.map(to_amp)
    if strand == "-":
        var["alt"] = var.alt.map(lambda b: rc(b) if b not in ("", "-") else b)
    # variants are called against 700223-20; where the primer differs from it, a genome mismatches the primer
    # unless it carries the primer base
    diff = {k - offset + 1: (amp[k], primer_amp[k]) for k in range(len(amp)) if amp[k] != primer_amp[k]}
    vpos = {g: [] for g in cov.genome}; vsub = {g: {} for g in cov.genome}
    for gid, grp in var.groupby("genome"):
        vpos[gid] = [p for p, typ in zip(grp.apos, grp.type) if p not in diff or typ != "sub"]
        vsub[gid] = {p: a for p, typ, a in zip(grp.apos, grp.type, grp.alt) if typ == "sub" and p in diff}
    for gid in cov.genome:
        vpos[gid] += [p for p, (ref_b, prim_b) in diff.items() if vsub[gid].get(p, ref_b) != prim_b]
    present = cov[(cov.qstart <= gene_lo) & (cov.qend >= gene_hi)].genome
    rows = []
    for gid in present:
        for _, r in st.iterrows():
            n = sum(r.start <= p <= r.end for p in vpos[gid])
            zone = range(r.end - END_NT + 1, r.end + 1) if r.orientation == "+" else range(r.start, r.start + END_NT)
            n3 = sum(p in zone for p in vpos[gid]) if r.site in THREE_PRIME_SITES else 0
            rows.append(dict(genome=gid, site=r.site, mismatches=n, mismatches_3prime=n3))
    mm = pd.DataFrame(rows)
    g = mm.groupby("site")
    agg = pd.DataFrame({"genomes": g.size(), "mm0": g.mismatches.apply(lambda x: (x == 0).sum()),
                        "mm1": g.mismatches.apply(lambda x: (x == 1).sum()),
                        "mm2": g.mismatches.apply(lambda x: (x == 2).sum()),
                        "mm3plus": g.mismatches.apply(lambda x: (x >= 3).sum()),
                        "end3": g.mismatches_3prime.apply(lambda x: (x > 0).sum())}).reset_index()
    st = st.merge(agg, on="site")
    for c in ("mm1", "mm2", "mm3plus", "end3"): st[f"pct_{c}"] = 100 * st[c] / st.genomes
    st.loc[~st.site.isin(THREE_PRIME_SITES), "end3"] = pd.NA
    pv = pd.read_csv(d / "position_variability.tsv", sep="\t")
    pv["pct"] = 100 * (pv.n_substitution + pv.n_gap) / pv.n_covered.where(pv.n_covered > 0)
    return dict(L=L, strand=strand, sites=st, amp_len=amp_len, gene_lo=gene_lo, gene_hi=gene_hi, pv=pv,
                n_total=len(cov), n_present=len(present), per_genome=mm)


def presence_table():
    """Binding sites of the papGII primer set present (covered by an alignment, primer_screen) per papG group.
    Groups from the per-locus typing against PapG-I to PapG-IV (results/papg_site_scan); genomes carrying papGII
    together with another allele are counted as papGII."""
    s = pd.read_csv(SCREEN / "genome_summary.tsv", sep="\t")
    m = pd.read_csv(SCREEN / "site_mismatches.tsv", sep="\t")
    m = m[m.target == "papGII"]
    loci = pd.read_csv(RESULTS / "papg_site_scan" / "locus_site_mismatches.tsv", sep="\t")
    alleles = loci.groupby("genome").allele.apply(set).reindex(s.genome).apply(
        lambda a: a if isinstance(a, set) else set())
    def group(a):
        if "PapG-II" in a: return "papGII (alone or with papGIII/papGIV)"
        if a == {"PapG-III"}: return "papGIII only"
        if a == {"PapG-IV"}: return "papGIV only"
        if not a: return "no papG allele"
        return "+".join(sorted(a))
    groups = alleles.map(group)
    rows = []
    for name in ["papGII (alone or with papGIII/papGIV)", "papGIII only", "papGIV only", "no papG allele"] + \
            sorted(set(groups) - {"papGII (alone or with papGIII/papGIV)", "papGIII only", "papGIV only", "no papG allele"}):
        gids = groups[groups == name].index
        if not len(gids): continue
        x = m[m.genome.isin(gids)].pivot(index="genome", columns="site", values="covered").reindex(gids).fillna(0)
        rows.append(dict(group=name, genomes=len(gids), **{site: int(x[site].sum()) for site in
                                                          ["F3", "F2", "LF", "F1c", "B1c", "LB", "B2", "B3"]}))
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = {t: load(t) for t in ("uidA", "papGII")}
    cols = ["target", "site", "role", "start", "end", "genomes", "mm0", "mm1", "mm2", "mm3plus", "end3"]
    tab = pd.concat([data[t]["sites"].assign(target=t) for t in data])[cols]
    tab.to_csv(OUT / "tableS_binding_site_mismatches.tsv", sep="\t", index=False)
    pres = presence_table()
    pres.to_csv(OUT / "tableS_binding_site_presence.tsv", sep="\t", index=False)
    print(tab.to_string(index=False)); print(pres.to_string(index=False))


if __name__ == "__main__":
    main()
