"""Figure 2A (lollipop version): LAMP primer binding sites and their conservation / off-target divergence.

uidA: per Clermont phylogroup (UPEC collection + NCBI RefSeq complete genomes pooled), one lollipop per binding site,
aligned under the primer map: % of genomes with >= 1 mismatch to the primer.
papGII: per papG class (both collections pooled), one lollipop per binding site: median number of mismatches to the
primer (thin line = range) at the best-matching papG locus of each genome.
Both collections are also reported separately (results/figures/fig2a_lollipop_*_by_collection.tsv).

Output: results/figures/fig2a_primer_lollipop.{pdf,png,svg}, fig2a_lollipop_uidA.tsv, fig2a_lollipop_papGII.tsv
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from fig_primer_map import (INK, INK2, MIN_GROUP, NAMED, ROLE_COLOR, ROLE_LABEL, XPAD, draw_map, fmt, sites_for,
                            load_site_mismatches)
from lamp_common import DATA, RESULTS

OUT = RESULTS / "figures"
ORDER = ["F3", "F2", "LF", "F1c", "B1c", "LB", "B2", "B3"]
GRID = "#e1e0d9"
mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": "lamp-fig2a"})


def uida_long():
    """One row per genome and binding site: collection, phylogroup, mismatches (genomes with the region present)."""
    u = load_site_mismatches("uidA")["per_genome"].merge(pd.read_csv(DATA / "upec_phylogroups.tsv", sep="\t"), on="genome")
    u = u.assign(collection="UPEC")[["genome", "site", "mismatches", "phylogroup", "collection"]]
    g = pd.read_csv(RESULTS / "ncbi_ecoli_uida" / "genome_sites.tsv", sep="\t")
    g = g[g.status == "present"]
    n = g.melt(id_vars="genome", value_vars=ORDER, var_name="site", value_name="mismatches")
    n = n.merge(pd.read_csv(RESULTS / "ncbi_ecoli_uida" / "phylogroups.tsv", sep="\t"), on="genome").assign(collection="NCBI")
    return pd.concat([u, n], ignore_index=True)


def uida_table(long):
    share = long.assign(any=long.mismatches > 0).groupby(["phylogroup", "site"])["any"].mean().unstack()[ORDER] * 100
    n = long.groupby("phylogroup").genome.nunique()
    keep = sorted(g for g in n.index if g in NAMED and n[g] >= MIN_GROUP)
    return share.loc[keep], n[keep], int(n.drop(keep).sum())


def papg_long():
    parts = []
    for folder, coll in (("papg_site_scan", "UPEC"), ("papg_site_scan_ncbi", "NCBI")):
        d = pd.read_csv(RESULTS / folder / "locus_site_mismatches.tsv", sep="\t")
        tot = d.groupby(["genome", "locus", "allele"]).mismatches.sum().reset_index()
        best = tot.loc[tot.groupby(["genome", "allele"]).mismatches.idxmin(), ["genome", "locus", "allele"]]
        parts.append(d.merge(best, on=["genome", "locus", "allele"]).assign(collection=coll))
    return pd.concat(parts, ignore_index=True)


def lollipops(ax, st, values, ymax, lo=None, hi=None, labels=True, thresh=0):
    """values: {site: y}; lo/hi: optional whiskers."""
    for _, r in st.iterrows():
        x, y = (r.start - 1 + r.end) / 2, values[r.site]
        c = ROLE_COLOR[r.role]
        small = y <= thresh
        if not small: ax.plot([x, x], [0, y], color=c, lw=1.4, zorder=2)
        ax.plot(x, y, "o", color=c, ms=2.6 if small else 5.4, alpha=0.5 if small else 1, mec="white",
                mew=0.0 if small else 0.6, zorder=3)
        if labels and y > thresh:
            ax.text(x + 4.5, y, fmt(round(y)) if ymax > 20 else fmt(y), ha="left", va="center", fontsize=6.5, color=INK)


def block(fig, gs_rows, st, xmax, rows, ymax, ylabel, **kw):
    axes = []
    for i, (label, values, extra) in enumerate(rows):
        ax = fig.add_subplot(gs_rows[i])
        lollipops(ax, st, values, ymax, **extra, **kw)
        ax.set_xlim(-XPAD, xmax); ax.set_ylim(-ymax * 0.12, ymax * 1.35)
        ax.axhline(0, color="#c3c2b7", lw=0.6)
        ax.set_yticks([]); ax.set_xticks([])
        for sp in ax.spines.values(): sp.set_visible(False)
        ax.text(-XPAD - 2, ymax * 0.15, label, ha="right", va="center", fontsize=7, color=INK2)
        axes.append(ax)
    return axes


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    u, u_pos, u_len = sites_for("uidA")
    p, p_pos, p_len = sites_for("papGII")
    xmax = max(u.end.max(), p.end.max()) + 8

    long = uida_long()
    share, n_u, excl = uida_table(long)
    by_coll = long.assign(any=long.mismatches > 0).groupby(["collection", "phylogroup", "site"])["any"].mean().unstack()[ORDER] * 100
    share.assign(genomes=n_u).to_csv(OUT / "fig2a_lollipop_uidA.tsv", sep="\t")
    by_coll.to_csv(OUT / "fig2a_lollipop_uidA_by_collection.tsv", sep="\t")

    pl = papg_long()
    n_g = pl.groupby("allele").genome.nunique()
    stats = pl.groupby(["allele", "site"]).mismatches.agg(["median", "min", "max"]).unstack()
    pl.groupby(["collection", "allele", "site"]).mismatches.median().unstack()[ORDER].to_csv(
        OUT / "fig2a_lollipop_papGII_by_collection.tsv", sep="\t")
    classes = [a for a in ("PapG-II", "PapG-III", "PapG-IV", "PapG-I") if a in n_g]
    label = {"PapG-I": "papGI", "PapG-II": "papGII", "PapG-III": "papGIII", "PapG-IV": "papGIV"}
    stats_out = pd.concat({a: stats.loc[a].unstack(0)[ORDER] if False else stats.loc[a] for a in classes}).to_csv(
        OUT / "fig2a_lollipop_papGII.tsv", sep="\t")

    rows_u = [(f"{g} (n={n_u[g]:,})", {s: share.loc[g, s] for s in ORDER}, {}) for g in share.index]
    rows_p = [(f"{label[a]} (n={n_g[a]:,})", {s: stats.loc[a, ("median", s)] for s in ORDER}, {}) for a in classes]

    RU, RP = 0.5, 0.58
    ratios = [1.0] + [RU] * len(rows_u) + [0.45, 1.0] + [RP] * len(rows_p)
    fig = plt.figure(figsize=(180 / 25.4, (sum(ratios) * 25 + 16) / 25.4))
    gs = fig.add_gridspec(len(ratios), 1, height_ratios=ratios, hspace=0.12)
    ax_u = fig.add_subplot(gs[0]); draw_map(ax_u, u, xmax)
    ax_p = fig.add_subplot(gs[len(rows_u) + 2]); draw_map(ax_p, p, xmax)
    axes_u = block(fig, [gs[1 + i] for i in range(len(rows_u))], u, xmax, rows_u, 100, "", thresh=9)
    axes_p = block(fig, [gs[len(rows_u) + 3 + i] for i in range(len(rows_p))], p, xmax, rows_p, 12, "", thresh=0.4)
    for ax, name, (lo, hi), L, st in ((ax_u, "uidA", u_pos, u_len, u), (ax_p, "papGII", p_pos, p_len, p)):
        ax.text(-XPAD, 1.45, name, ha="left", va="bottom", fontsize=9, fontstyle="italic", color=INK)
        ax.text(xmax, 1.45, f"{st.end.max()}-bp amplicon · nt {lo}–{hi} of {L}-bp gene", ha="right", va="bottom",
                fontsize=6.5, color=INK2)
    handles = [Patch(color=ROLE_COLOR[k], label=ROLE_LABEL[k]) for k in ROLE_COLOR]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=7, bbox_to_anchor=(0.5, -0.01))
    fig.subplots_adjust(left=0.1, right=0.98, top=0.97, bottom=0.07)
    fig.canvas.draw()
    caps = ((axes_u[-1], "% of genomes with ≥1 mismatch to the primer at each binding site; UPEC and NCBI genomes pooled; phylogroups with "
                         f"<{MIN_GROUP} genomes or unassigned not shown ({excl:,} genomes)"),
            (axes_p[-1], "median number of mismatches to the primer at the best-matching papG locus of each genome; UPEC and NCBI genomes pooled "
                         "(ranges: Supplementary Table)"))
    for ax_last, text in caps:
        pos = ax_last.get_position()
        fig.text(0.98, pos.y0 - 0.004, text, ha="right", va="top", fontsize=6.3, color=INK2)
    for ext in ("pdf", "png", "svg"):
        meta = {"pdf": {"CreationDate": None, "Creator": "fig_primer_lollipop.py"}, "svg": {"Date": None}}.get(ext)
        fig.savefig(OUT / f"fig2a_primer_lollipop.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight", metadata=meta)
    print(share.round(0).assign(n=n_u).to_string()); print(stats.round(1).to_string())


if __name__ == "__main__":
    main()
