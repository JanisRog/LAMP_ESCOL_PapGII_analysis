"""Figure 2A (lollipop version): LAMP primer binding sites and median mismatches to the primers.

Per target: the primer map (top; forward primers above, reverse primers below the amplicon) and one lollipop row per
group underneath, aligned to it:
- uidA: median number of mismatches to the primer at each binding site, per Clermont phylogroup (UPEC collection +
  NCBI RefSeq complete genomes pooled; genomes with the amplicon region present).
- papGII: median number of mismatches at the best-matching papG locus of each genome, per papG class (both
  collections pooled).
All rows share one y-axis from 0 to the length of the longest primer binding site (25 nt), so the height of a lollipop
is the part of the primer that does not match. Groups differ by colour and marker shape. Both collections are also
reported separately in the output tables.

Output: results/figures/fig2a_primer_lollipop.{pdf,png,svg}; fig2a_lollipop_{uidA,papGII}.tsv (pooled, with n, and
range for papGII); fig2a_lollipop_{uidA,papGII}_by_collection.tsv
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from fig_primer_map import MIN_GROUP, NAMED, ROLE_COLOR, ROLE_LABEL, XPAD, draw_map, fmt, load_site_mismatches, sites_for
from lamp_common import DATA, RESULTS

OUT = RESULTS / "figures"
ORDER = ["F3", "F2", "LF", "F1c", "B1c", "LB", "B2", "B3"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e1e0d9"
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]   # fixed order
MARKERS = ["o", "s", "^", "D", "v", "P", "h", "X"]
PHYLO = ["A", "B1", "B2", "C", "D", "E", "F", "G"]
PAPG = [("PapG-II", "papGII"), ("PapG-III", "papGIII"), ("PapG-IV", "papGIV"), ("PapG-I", "papGI")]

mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": "lamp-fig2a"})


def uida_long():
    u = load_site_mismatches("uidA")["per_genome"].merge(pd.read_csv(DATA / "upec_phylogroups.tsv", sep="\t"), on="genome")
    u = u.assign(collection="UPEC")[["genome", "site", "mismatches", "phylogroup", "collection"]]
    g = pd.read_csv(RESULTS / "ncbi_ecoli_uida" / "genome_sites.tsv", sep="\t")
    g = g[g.status == "present"]
    n = g.melt(id_vars="genome", value_vars=ORDER, var_name="site", value_name="mismatches")
    n = n.merge(pd.read_csv(RESULTS / "ncbi_ecoli_uida" / "phylogroups.tsv", sep="\t"), on="genome").assign(collection="NCBI")
    return pd.concat([u, n], ignore_index=True)


def papg_long():
    parts = []
    for folder, coll in (("papg_site_scan", "UPEC"), ("papg_site_scan_ncbi", "NCBI")):
        d = pd.read_csv(RESULTS / folder / "locus_site_mismatches.tsv", sep="\t")
        tot = d.groupby(["genome", "locus", "allele"]).mismatches.sum().reset_index()
        best = tot.loc[tot.groupby(["genome", "allele"]).mismatches.idxmin(), ["genome", "locus", "allele"]]
        parts.append(d.merge(best, on=["genome", "locus", "allele"]).assign(collection=coll))
    return pd.concat(parts, ignore_index=True)


CONNECT = []   # (gene axis, primer-map axis, [(gene position, amplicon position), ...]) -> zoom lines drawn after layout


def draw_zoom_lines(fig):
    """Dashed lines from the amplicon edges in the gene strip to the ends of the primer map (figure coordinates)."""
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    for axg, axt, ends in CONNECT:
        for xg, xa in ends:
            x0, y0 = inv.transform(axg.transData.transform((xg, 0)))
            x1, y1 = inv.transform(axt.transData.transform((xa, 1.58)))
            fig.add_artist(Line2D([x0, x1], [y0, y1], transform=fig.transFigure, color="#898781", lw=0.7,
                                  ls=(0, (3, 2)), zorder=5))


def pooled_variability(target):
    """Per-position counts of the UPEC collection (results/gene_scan) and the NCBI genomes (results/gene_scan_ncbi), added
    up (both scans use the same reference gene of E. coli 700223-20)."""
    parts = [pd.read_csv(RESULTS / folder / target / "position_variability.tsv", sep="\t")
             for folder in ("gene_scan", "gene_scan_ncbi")]
    assert (parts[0].position == parts[1].position).all() and (parts[0].ref_base == parts[1].ref_base).all()
    pv = parts[0].copy()
    for c in ("n_covered", "n_substitution", "n_gap"):
        pv[c] = parts[0][c] + parts[1][c]
    return pv


def gene_strip(ax, target, st):
    """Whole gene of the design strain: % of genomes (UPEC + NCBI) with a variant at each position; amplicon shaded.
    For papGII the primers lie on the reverse strand, so the gene axis is reversed (F3 -> B3 runs left to right)."""
    pv = pooled_variability(target)
    gs_ = pd.read_csv(RESULTS / "gene_scan" / target / "gene_sites.tsv", sep="\t").set_index("feature")
    lo, hi, strand = int(gs_.loc["amplicon", "gene_start"]), int(gs_.loc["amplicon", "gene_end"]), gs_.loc["amplicon", "strand"]
    pct = 100 * (pv.n_substitution + pv.n_gap) / pv.n_covered.where(pv.n_covered > 0)
    ax.axvspan(lo, hi, color="#e8e6df", lw=0, zorder=0)
    ax.bar(pv.position, pct, width=1.0, color="#52514e", lw=0, zorder=2)
    L = len(pv)
    ax.set_xlim((L + 1, 0) if strand == "-" else (0, L + 1)); ax.set_ylim(0, 100); ax.set_yticks([0, 50, 100])
    ax.set_xticks([1, lo, hi, L]); ax.xaxis.tick_top()
    ax.tick_params(axis="x", labelsize=6, colors=INK2, length=2, top=True, bottom=False, labeltop=True, labelbottom=False)
    ax.tick_params(axis="y", labelsize=6, colors=INK2, length=2); ax.set_ylabel("Variant (%)", fontsize=6.5, color=INK2)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color("#c3c2b7")
    return lo, hi, strand


def block(fig, gs, st, xmax, groups, values, ymax, title, sub, target=None, xlabel=True):
    """gs: subgridspec with 1 + len(groups) rows (track, group rows)."""
    axg = fig.add_subplot(gs[0])
    lo, hi, strand = gene_strip(axg, target, st)
    axg.set_title(title, loc="left", fontsize=9, fontstyle="italic", color=INK, pad=17)
    axg.text(1.0, 1.36, sub, transform=axg.transAxes, ha="right", va="bottom", fontsize=6.2, color=INK2)
    axt = fig.add_subplot(gs[1]); draw_map(axt, st, xmax)
    amp_len = st.end.max()
    ends = ((lo, 1), (hi, amp_len)) if strand == "+" else ((hi, 1), (lo, amp_len))
    CONNECT.append((axg, axt, ends))
    axes = []
    for i, (key, label, c, m) in enumerate(groups):
        ax = fig.add_subplot(gs[2 + i], sharex=axt)
        for _, r in st.iterrows():
            ax.axvspan(r.start - 1, r.end, color="#f4f3ef", lw=0, zorder=0)
            y, x = values[key][r.site], (r.start - 1 + r.end) / 2
            if y > 0:
                ax.plot([x, x], [0, y], color=c, lw=1.3, zorder=2)
                ax.text(x + 2.6, y, fmt(y), ha="left", va="center", fontsize=6.3, color=INK)
            ax.plot(x, y, marker=m, color=c, ms=4.4 if y > 0 else 2.6, alpha=1 if y > 0 else 0.55,
                    mec="white", mew=0.4 if y > 0 else 0, zorder=3, ls="none")
        ax.set_ylim(-ymax * 0.1, ymax * 1.02); ax.set_xlim(-XPAD, xmax)
        ax.set_yticks([0, ymax]); ax.tick_params(axis="y", labelsize=6, colors=INK2, length=2)
        ax.axhline(0, color="#c3c2b7", lw=0.6, zorder=1)
        for sp in ("top", "right", "bottom"): ax.spines[sp].set_visible(False)
        ax.spines["left"].set_color("#c3c2b7"); ax.tick_params(axis="x", bottom=False, labelbottom=False)
        ax.text(-0.105, 0.5, label, transform=ax.transAxes, ha="right", va="center", fontsize=7, color=INK2)
        ax.plot(-0.075, 0.5, marker=m, color=c, ms=5, transform=ax.transAxes, clip_on=False, ls="none")
        axes.append(ax)
    last = axes[-1]
    last.spines["bottom"].set_visible(True); last.spines["bottom"].set_color("#c3c2b7")
    last.tick_params(axis="x", bottom=True, labelbottom=True, labelsize=7, colors=INK2)
    last.set_xlabel("Position in amplicon (nt)", fontsize=7, color=INK2)
    return axt, axes


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    u, u_pos, u_len = sites_for("uidA")
    p, p_pos, p_len = sites_for("papGII")
    xmax = max(u.end.max(), p.end.max()) + 6
    ymax = int(max((u.end - u.start + 1).max(), (p.end - p.start + 1).max()))          # longest binding site = 25 nt

    long = uida_long()
    n_u = long.groupby("phylogroup").genome.nunique()
    keep = [g for g in PHYLO if g in NAMED and n_u.get(g, 0) >= MIN_GROUP]
    med = long.groupby(["phylogroup", "site"]).mismatches.median().unstack()[ORDER]
    med.loc[keep].assign(genomes=n_u[keep]).to_csv(OUT / "fig2a_lollipop_uidA.tsv", sep="\t")
    long.groupby(["collection", "phylogroup", "site"]).mismatches.median().unstack()[ORDER].to_csv(
        OUT / "fig2a_lollipop_uidA_by_collection.tsv", sep="\t")
    excl = int(n_u.drop(keep).sum())

    pl = papg_long()
    n_g = pl.groupby("allele").genome.nunique()
    stat = pl.groupby(["allele", "site"]).mismatches.agg(["median", "min", "max"]).reset_index()
    stat.assign(genomes=stat.allele.map(n_g)).to_csv(OUT / "fig2a_lollipop_papGII.tsv", sep="\t", index=False)
    pl.groupby(["collection", "allele", "site"]).mismatches.median().unstack()[ORDER].to_csv(
        OUT / "fig2a_lollipop_papGII_by_collection.tsv", sep="\t")
    classes = [(a, l) for a, l in PAPG if a in n_g]

    gu = [(g, f"{g} (n={n_u[g]:,})", COLORS[i], MARKERS[i]) for i, g in enumerate(keep)]
    vu = {g: med.loc[g].to_dict() for g in keep}
    gp = [(a, f"{l} (n={n_g[a]:,})", COLORS[i], MARKERS[i]) for i, (a, l) in enumerate(classes)]
    vp = {a: stat[stat.allele == a].set_index("site")["median"].to_dict() for a, _ in classes}

    ROW, TRK, GENE = 8.5, 24.0, 17.0                       # mm per group row / primer map / gene-wide strip
    h_u, h_p = GENE + TRK + ROW * len(gu), GENE + TRK + ROW * len(gp)
    total = h_u + h_p + 84
    fig = plt.figure(figsize=(180 / 25.4, total / 25.4))
    outer = fig.add_gridspec(2, 1, height_ratios=[h_u, h_p], hspace=0.2)
    gs_u = outer[0].subgridspec(2 + len(gu), 1, height_ratios=[GENE, TRK] + [ROW] * len(gu), hspace=0.3)
    gs_p = outer[1].subgridspec(2 + len(gp), 1, height_ratios=[GENE, TRK] + [ROW] * len(gp), hspace=0.3)
    block(fig, gs_u, u, xmax, gu, vu, ymax, "uidA",
          f"{u.end.max()}-bp amplicon (nt {u_pos[0]}–{u_pos[1]}) of the {u_len}-bp gene · variants: UPEC + NCBI genomes (n = 6,231)", target="uidA")
    block(fig, gs_p, p, xmax, gp, vp, ymax, "papGII",
          f"{p.end.max()}-bp amplicon (nt {p_pos[0]}–{p_pos[1]}) of the {p_len}-bp gene · variants: papGII genomes of UPEC + NCBI (n = 645); axis reversed", target="papGII")
    fig.text(0.5, 0.043, f"Rows: {'Clermont phylogroups' if False else 'groups'}, UPEC + NCBI pooled. Lollipop height = median number of mismatches to the primer (y-axis: 0–{ymax} nt, the length of the longest primer binding site)",
             ha="center", va="center", fontsize=6.8, color=INK2)
    role_handles = [Patch(color=ROLE_COLOR[k], label=ROLE_LABEL[k]) for k in ROLE_COLOR]
    fig.legend(handles=role_handles, loc="lower center", ncol=3, frameon=False, fontsize=7, bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(left=0.17, right=0.98, top=0.96, bottom=0.105)
    draw_zoom_lines(fig)
    for ext in ("pdf", "png", "svg"):
        meta = {"pdf": {"CreationDate": None, "Creator": "fig_primer_lollipop.py"}, "svg": {"Date": None}}.get(ext)
        fig.savefig(OUT / f"fig2a_primer_lollipop.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight", metadata=meta)
    print("ymax", ymax, "| groups", keep, [c for c, _ in classes], "| excluded phylogroup genomes", excl)


if __name__ == "__main__":
    main()
