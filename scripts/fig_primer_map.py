"""Figure 2A: LAMP primer binding sites on the uidA and papGII amplicons.
uidA: median number of mismatches to the primer per binding site within each Clermont phylogroup (EzClermont,
data/upec_phylogroups.tsv; Clermont phylogroups with >= MIN_GROUP genomes in which the amplicon region is present; unassigned 'U'/cryptic
not shown).
papGII: median mismatches of each binding site per papG allele (papGII, papGIII, papGIV) across the papG loci of
the same genomes (results/papg_site_scan: loci typed against the PapG-I to PapG-IV reference alleles; best-matching
locus per genome and allele). No papGI occurs in the collection, so its row shows the PapG-I reference allele.

Amplicons and binding-site coordinates come from the reference genes of the design strain 700223-20
(results/gene_scan/*/reference_gene.fasta) with the FIP/BIP split lengths of the design (primer_screen.SPLIT).
Mismatches: best ungapped match of each binding-site sequence (either strand) within the full allele.

Output: results/figures/fig2a_primer_map.{pdf,png,svg} and fig2a_binding_sites.tsv
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrow, Patch, Rectangle

from lamp_common import DATA, RESULTS
from primer_screen import ROLE, binding_sites, rc, read_fasta, reference_amplicon
from tables_s_primer_conservation import load as load_site_mismatches

OUT = RESULTS / "figures"
ROLE_COLOR = {"outer": "#2a78d6", "inner": "#eb6834", "loop": "#1baf7a"}   # categorical slots 1-3
ROLE_LABEL = {"outer": "Outer primers (F3, B3)", "inner": "Inner primers (FIP = F1c+F2, BIP = B1c+B2)",
              "loop": "Loop primers (LF, LB)"}
# mismatch share of the site -> ordinal neutral ramp (light = few, dark = many); text switches to white on dark
MM_STEPS = [(0.0, "#f0efec", "#0b0b0b"), (0.10, "#c3c2b7", "#0b0b0b"), (0.30, "#898781", "#ffffff"),
            (1.01, "#52514e", "#ffffff")]
INK, INK2 = "#0b0b0b", "#52514e"
XPAD = 4
MIN_GROUP = 10

mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none"})


def best_mismatches(q, text):
    return min(sum(a != b for a, b in zip(qq, text[i:i + len(qq)]))
               for qq in (q, rc(q)) for i in range(len(text) - len(qq) + 1))


def sites_for(target):
    gene = next(iter(read_fasta(RESULTS / "gene_scan" / target / "reference_gene.fasta").values()))
    amp, offset = reference_amplicon(target, {"gene": gene})
    sites, _ = binding_sites(target, amp)
    gs = pd.read_csv(RESULTS / "gene_scan" / target / "gene_sites.tsv", sep="\t").set_index("feature")
    st = pd.DataFrame([dict(target=target, site=k, role=ROLE[k], start=s - offset + 1, end=e - offset, orientation=o,
                            length=e - s, sequence=seq) for k, (seq, s, e, o, _) in sites.items()]).sort_values("start")
    return st, (int(gs.loc["amplicon", "gene_start"]), int(gs.loc["amplicon", "gene_end"])), len(gene)


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
    ax.set_xlim(-XPAD, xmax); ax.set_ylim(-1.5, 1.5); ax.axis("off")
    ax.text(-XPAD, 0, "5′", ha="right", va="center", fontsize=7, color=INK2)
    ax.text(st.end.max() + 2, 0, "3′", ha="left", va="center", fontsize=7, color=INK2)


def fmt(v):
    return str(int(v)) if float(v).is_integer() else f"{v:.1f}"


def draw_mismatch_rows(ax, st, rows, xmax, caption):
    """One row per (label, {site: mismatches}); a cell under each binding site, shaded by mismatches / site length."""
    for row, (label, values) in enumerate(rows):
        y = len(rows) - row - 1
        ax.text(-XPAD - 2, y + 0.5, label, ha="right", va="center", fontsize=7, color=INK2)
        for _, r in st.iterrows():
            mm = values[r.site]
            frac = mm / r.length
            fill, ink = next((f, k) for lim, f, k in MM_STEPS if frac < lim) if mm else MM_STEPS[0][1:]
            ax.add_patch(Rectangle((r.start - 1 + 1.2, y + 0.08), r.end - r.start + 1 - 2.4, 0.84, color=fill, lw=0))
            ax.text((r.start - 1 + r.end) / 2, y + 0.5, fmt(mm), ha="center", va="center", fontsize=7, color=ink)
    ax.set_xlim(-XPAD, xmax); ax.set_ylim(0, len(rows)); ax.axis("off")
    ax.text(1.0, -0.06 * 6 / max(len(rows), 1), caption, transform=ax.transAxes, ha="right", va="top", fontsize=6.5,
            color=INK2)


def papg_rows(p, alleles):
    """Median mismatches per papG allele group (one locus per genome and allele: the best-matching one)."""
    d = pd.read_csv(RESULTS / "papg_site_scan" / "locus_site_mismatches.tsv", sep="\t")
    tot = d.groupby(["genome", "locus", "allele"]).mismatches.sum().reset_index()
    best = tot.loc[tot.groupby(["genome", "allele"]).mismatches.idxmin(), ["genome", "locus", "allele"]]
    d = d.merge(best, on=["genome", "locus", "allele"])
    med = d.groupby(["allele", "site"]).mismatches.median().unstack()
    n = d.groupby("allele").genome.nunique()
    label = {"PapG-I": "papGI", "PapG-II": "papGII", "PapG-III": "papGIII", "PapG-IV": "papGIV"}
    rows = [(f"{label[a]} (n={n[a]})", med.loc[a].to_dict()) for a in ("PapG-II", "PapG-III", "PapG-IV", "PapG-I")
            if a in n]
    if "PapG-I" not in n:
        rows.append(("papGI (reference)", dict(zip(p.site, [best_mismatches(s, alleles["PapG-I"]) for s in p.sequence]))))
    return rows, med.assign(genomes=n)


def phylogroup_rows():
    """Median mismatches to the uidA primers per binding site within each phylogroup (>= MIN_GROUP genomes)."""
    mm = load_site_mismatches("uidA")["per_genome"]
    pg = pd.read_csv(DATA / "upec_phylogroups.tsv", sep="\t")
    mm = mm.merge(pg, on="genome")
    n = mm.groupby("phylogroup").genome.nunique().sort_values(ascending=False)
    named = [g for g in n.index if g in {"A", "B1", "B2", "C", "D", "E", "F", "G"}]  # Clermont phylogroups only
    keep = [g for g in named if n[g] >= MIN_GROUP]
    med = mm[mm.phylogroup.isin(keep)].groupby(["phylogroup", "site"]).mismatches.median().unstack()
    rows = [(f"{g} (n={n[g]})", med.loc[g].to_dict()) for g in keep]
    excluded = int(n.drop(keep).sum())
    table = med.loc[keep].assign(genomes=n[keep])
    return rows, excluded, table, n


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    alleles = read_fasta(DATA / "papG_reference_alleles.fasta")
    u, u_pos, u_len = sites_for("uidA")
    p, p_pos, p_len = sites_for("papGII")
    for name in ("PapG-III", "PapG-I", "PapG-II"):
        p[f"mm_{name}_reference"] = [best_mismatches(s, alleles[name]) for s in p.sequence]
    grows, gtable = papg_rows(p, alleles)
    xmax = max(u.end.max(), p.end.max()) + 8
    prows, excluded, ptable, pn = phylogroup_rows()

    fig = plt.figure(figsize=(180 / 25.4, 125 / 25.4))
    gs = fig.add_gridspec(5, 1, height_ratios=[1.0, 1.35, 0.32, 1.0, 0.95], hspace=0.12)
    ax_u, ax_ut, ax_p, ax_t = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[3]), fig.add_subplot(gs[4])
    draw_map(ax_u, u, xmax)
    draw_mismatch_rows(ax_ut, u, prows, xmax,
                       f"Median mismatches to the primer per Clermont phylogroup ({excluded} genomes in smaller "
                       f"or unassigned groups not shown)")
    draw_map(ax_p, p, xmax)
    draw_mismatch_rows(ax_t, p, grows, xmax, "Median mismatches to the primer per papG allele in the same genomes "
                       "(papGI: reference allele, absent from the collection)")
    for ax, name, (lo, hi), L, st in ((ax_u, "uidA", u_pos, u_len, u), (ax_p, "papGII", p_pos, p_len, p)):
        ax.text(-XPAD, 1.45, name, ha="left", va="bottom", fontsize=9, fontstyle="italic", color=INK)
        ax.text(xmax, 1.45, f"{st.end.max()}-bp amplicon · nt {lo}–{hi} of {L}-bp gene", ha="right", va="bottom",
                fontsize=6.5, color=INK2)
    handles = [Patch(color=ROLE_COLOR[k], label=ROLE_LABEL[k]) for k in ROLE_COLOR]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=7, bbox_to_anchor=(0.5, -0.02))
    fig.subplots_adjust(left=0.1, right=0.98, top=0.95, bottom=0.1)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(OUT / f"fig2a_primer_map.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight")
    tab = pd.concat([u, p])
    tab.to_csv(OUT / "fig2a_binding_sites.tsv", sep="\t", index=False)
    ptable.to_csv(OUT / "fig2a_uidA_phylogroup_median_mismatches.tsv", sep="\t")
    gtable.to_csv(OUT / "fig2a_papG_allele_median_mismatches.tsv", sep="\t")
    print(gtable.to_string())
    print(ptable.to_string()); print("excluded genomes:", excluded)
    print(tab.drop(columns="sequence").to_string(index=False))


if __name__ == "__main__":
    main()
