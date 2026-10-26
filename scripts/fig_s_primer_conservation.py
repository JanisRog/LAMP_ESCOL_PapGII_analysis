"""Supplementary figure and table: descriptive in silico conservation of the LAMP binding sites in 738 UPEC genomes.

Figure, per target: (1) % of genomes with a variant at each position of the reference gene (E. coli 700223-20),
amplicon shaded; (2) amplicon map with binding sites; (3) per site, % of genomes with 1, 2 or >= 3 mismatches to
the primer and with a mismatch in the 3'-terminal 5 nt of primers with a free 3' end. Descriptive only: the effect
of mismatches on amplification was not tested experimentally.
Genomes: uidA - all genomes in which the amplicon region is present; papGII - genomes typed papGII.
Presence of binding sites per papG allele group comes from results/primer_screen/upec.

Output: results/figures/figS_primer_conservation.{pdf,png,svg}, tableS_binding_site_mismatches.tsv,
        tableS_binding_site_presence.tsv
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, FancyArrow, Patch

from lamp_common import RESULTS
from primer_screen import ROLE, THREE_PRIME_SITES, binding_sites, rc, read_fasta, reference_amplicon

SCAN, SCREEN, OUT = RESULTS / "gene_scan", RESULTS / "primer_screen" / "upec", RESULTS / "figures"
ROLE_COLOR = {"outer": "#2a78d6", "inner": "#eb6834", "loop": "#1baf7a"}   # categorical slots 1-3
ROLE_LABEL = {"outer": "Outer (F3, B3)", "inner": "Inner (FIP = F1c+F2, BIP = B1c+B2)", "loop": "Loop (LF, LB)"}
MM_COLOR = {"mm1": "#c3c2b7", "mm2": "#898781", "mm3plus": "#52514e"}     # neutral ordinal ramp
INK, INK2, GRID, SHADE = "#0b0b0b", "#52514e", "#e1e0d9", "#f0efec"
END_NT, XPAD = 5, 4

mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 8, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "pdf.fonttype": 42, "svg.fonttype": "none"})


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
                n_total=len(cov), n_present=len(present))


def draw_gene(ax, t):
    pv = t["pv"]
    ax.axvspan(t["gene_lo"], t["gene_hi"], color=SHADE, zorder=0, lw=0)
    ax.bar(pv.position, pv.pct, width=1.0, color=INK2, lw=0, zorder=2)
    ax.set_xlim((t["L"] + 1, 0) if t["strand"] == "-" else (0, t["L"] + 1))
    ax.set_ylim(0, 100); ax.set_yticks([0, 50, 100])
    ax.set_ylabel("Variant vs\n700223-20 (%)", fontsize=7)
    ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=6.5, pad=1)


def draw_map(ax, st, xmax):
    ax.plot([1, st.end.max()], [0, 0], color=INK2, lw=1.2, solid_capstyle="butt", zorder=1)
    for _, r in st.iterrows():
        y = 0.55 if r.orientation == "+" else -0.55
        x0, x1 = (r.start - 1, r.end) if r.orientation == "+" else (r.end, r.start - 1)
        ax.add_patch(FancyArrow(x0, y, x1 - x0, 0, width=0.32, head_width=0.32, length_includes_head=True,
                                head_length=min(6, abs(x1 - x0) * 0.35), color=ROLE_COLOR[r.role], lw=0, zorder=2))
        up = y > 0
        ax.text((r.start + r.end) / 2, y + (0.4 if up else -0.4), r.site, ha="center", va="bottom" if up else "top",
                fontsize=7.5, color=INK, fontweight="bold")
    ax.set_xlim(-XPAD, xmax); ax.set_ylim(-1.6, 1.6); ax.axis("off")
    ax.text(-XPAD, 0, "5′", ha="right", va="center", fontsize=7, color=INK2)
    ax.text(st.end.max() + 2, 0, "3′", ha="left", va="center", fontsize=7, color=INK2)


def draw_bars(ax, st, xmax, xlabel):
    for _, r in st.iterrows():
        x, w, bottom = (r.start - 1 + r.end) / 2, r.end - r.start + 1, 0
        for c in ("mm1", "mm2", "mm3plus"):
            h = r[f"pct_{c}"]
            if h > 0: ax.bar(x, h, width=w, bottom=bottom, color=MM_COLOR[c], edgecolor="white", linewidth=0.6, zorder=2)
            bottom += h
        if pd.notna(r.end3) and r.pct_end3 > 0:
            ax.plot(x, r.pct_end3, marker="D", ms=3.5, color=INK, markeredgecolor="white", markeredgewidth=0.6,
                    zorder=3, lw=0)
    ax.set_xlim(-XPAD, xmax); ax.set_ylim(0, 60); ax.set_yticks([0, 20, 40, 60])
    ax.set_ylabel("Genomes with\nmismatches (%)", fontsize=7)
    if xlabel: ax.set_xlabel("Position in amplicon (nt)", fontsize=7)
    ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)
    ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=6.5)


def presence_table():
    """Binding sites present (covered by an alignment) per papG allele group, papGII primer set."""
    s = pd.read_csv(SCREEN / "genome_summary.tsv", sep="\t")
    m = pd.read_csv(SCREEN / "site_mismatches.tsv", sep="\t")
    m = m[m.target == "papGII"]
    groups = {"papGII": s[s.papG_alleles.str.contains(r"PapG-II\(", regex=True)].genome,
              "papGIII only": s[s.papG_alleles.str.match(r"^PapG-III\(\d+\)$")].genome,
              "no papG allele called": s[s.papG_alleles == "none"].genome}
    rows = []
    for name, gids in groups.items():
        x = m[m.genome.isin(gids)].pivot(index="genome", columns="site", values="covered").reindex(gids).fillna(0)
        rows.append(dict(group=name, genomes=len(gids), **{site: int(x[site].sum()) for site in
                                                          ["F3", "F2", "LF", "F1c", "B1c", "LB", "B2", "B3"]}))
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = {t: load(t) for t in ("uidA", "papGII")}
    xmax = max(t["amp_len"] for t in data.values()) + 8
    fig = plt.figure(figsize=(180 / 25.4, 150 / 25.4))
    gs = fig.add_gridspec(7, 1, height_ratios=[0.6, 1.0, 0.9, 0.3, 0.6, 1.0, 0.9], hspace=0.45)
    tables = []
    for i, target in enumerate(("uidA", "papGII")):
        t, st = data[target], data[target]["sites"]
        ax_g, ax_m, ax_b = (fig.add_subplot(gs[4 * i + k]) for k in range(3))
        draw_gene(ax_g, t); draw_map(ax_m, st, xmax); draw_bars(ax_b, st, xmax, xlabel=(target == "papGII"))
        ends = ((t["gene_lo"], 1), (t["gene_hi"], t["amp_len"])) if t["strand"] == "+" else \
               ((t["gene_hi"], 1), (t["gene_lo"], t["amp_len"]))
        for xg, xa in ends:
            fig.add_artist(ConnectionPatch((xg, -0.38), (xa, 1.5), ax_g.get_xaxis_transform(), "data",
                                           axesA=ax_g, axesB=ax_m, color="#c3c2b7", lw=0.6, ls="--"))
        ax_g.set_title(target, loc="left", fontsize=9, fontstyle="italic", color=INK, pad=3)
        note = (f"amplicon nt {t['gene_lo']}–{t['gene_hi']} · region present in {t['n_present']}/{t['n_total']} UPEC genomes"
                if target == "uidA" else
                f"amplicon nt {t['gene_lo']}–{t['gene_hi']} · all binding sites present in {t['n_present']}/{t['n_total']} papGII genomes")
        ax_g.text(1.0, 1.06, note, transform=ax_g.transAxes, ha="right", va="bottom", fontsize=6.5, color=INK2)
        if t["strand"] == "-":
            ax_g.text(0.075, 1.06, "gene axis reversed (primers on the reverse strand)", transform=ax_g.transAxes,
                      ha="left", va="bottom", fontsize=6.5, color=INK2)
        tables.append(st.assign(target=target))
    handles = [Patch(color=ROLE_COLOR[k], label=ROLE_LABEL[k]) for k in ROLE_COLOR]
    handles += [Patch(color=MM_COLOR["mm1"], label="1 mismatch"), Patch(color=MM_COLOR["mm2"], label="2 mismatches"),
                Patch(color=MM_COLOR["mm3plus"], label="≥3 mismatches"),
                Line2D([], [], marker="D", ms=4, color=INK, lw=0, label="Mismatch in 3′-terminal 5 nt")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=7, bbox_to_anchor=(0.5, -0.03))
    fig.subplots_adjust(left=0.09, right=0.98, top=0.96, bottom=0.11)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(OUT / f"figS_primer_conservation.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight")
    cols = ["target", "site", "role", "start", "end", "genomes", "mm0", "mm1", "mm2", "mm3plus", "end3"]
    tab = pd.concat(tables)[cols]
    tab.to_csv(OUT / "tableS_binding_site_mismatches.tsv", sep="\t", index=False)
    pres = presence_table()
    pres.to_csv(OUT / "tableS_binding_site_presence.tsv", sep="\t", index=False)
    print(tab.to_string(index=False)); print(pres.to_string(index=False))


if __name__ == "__main__":
    main()
