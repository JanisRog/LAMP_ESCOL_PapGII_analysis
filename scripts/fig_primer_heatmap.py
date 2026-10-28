"""Figure 2A (heatmap version): LAMP primer maps with a group x binding-site heatmap of median mismatches.

Same data and layout as fig_primer_lollipop.py (gene-wide variant strip, primer map, one row per group aligned to the
binding sites), but each binding site is a cell coloured and labelled by the median number of mismatches to the primer
(UPEC + NCBI genomes pooled). One colour scale (0 to the highest median) is shared by both targets.

Output: results/figures/fig2a_primer_heatmap.{pdf,png,svg}
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle

import fig_primer_lollipop as L
from fig_primer_lollipop import CONNECT, INK, INK2, ORDER, PAPG, PHYLO, draw_zoom_lines, gene_strip, papg_long, uida_long
from fig_primer_map import MIN_GROUP, NAMED, XPAD, draw_map, fmt, sites_for
from lamp_common import RESULTS

OUT = RESULTS / "figures"
CMAP = mpl.colormaps["YlOrRd"]


def block(fig, gs, st, xmax, groups, values, norm, title, sub, target):
    axg = fig.add_subplot(gs[0])
    lo, hi, strand = gene_strip(axg, target, st)
    axg.set_title(title, loc="left", fontsize=9, fontstyle="italic", color=INK, pad=17)
    axg.text(1.0, 1.36, sub, transform=axg.transAxes, ha="right", va="bottom", fontsize=6.2, color=INK2)
    axt = fig.add_subplot(gs[2]); draw_map(axt, st, xmax)
    amp_len = st.end.max()
    CONNECT.append((axg, axt, ((lo, 1), (hi, amp_len)) if strand == "+" else ((hi, 1), (lo, amp_len))))
    ax = fig.add_subplot(gs[4], sharex=axt)
    n = len(groups)
    for i, (key, label) in enumerate(groups):
        for _, r in st.iterrows():
            v = values[key][r.site]
            ax.add_patch(Rectangle((r.start - 1, n - 1 - i - 0.5), r.end - r.start + 1, 1, facecolor=CMAP(norm(v)),
                                   edgecolor="white", lw=1.4))
            ax.text((r.start - 1 + r.end) / 2, n - 1 - i, fmt(v), ha="center", va="center", fontsize=6.8,
                    color="white" if norm(v) > 0.55 else INK)
    ax.set_xlim(-XPAD, xmax); ax.set_ylim(-0.5, n - 0.5)
    ax.set_yticks(range(n)); ax.set_yticklabels([g[1] for g in groups][::-1], fontsize=7, color=INK2)
    ax.tick_params(axis="y", length=0); ax.tick_params(axis="x", labelsize=7, colors=INK2)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.spines["bottom"].set_visible(True); ax.spines["bottom"].set_color("#c3c2b7")
    ax.set_xlabel("Position in amplicon (nt)", fontsize=7, color=INK2)
    return axt, ax


def main():
    u, u_pos, u_len = sites_for("uidA")
    p, p_pos, p_len = sites_for("papGII")
    xmax = max(u.end.max(), p.end.max()) + 6
    long = uida_long()
    n_u = long.groupby("phylogroup").genome.nunique()
    keep = [g for g in PHYLO if g in NAMED and n_u.get(g, 0) >= MIN_GROUP]
    med = long.groupby(["phylogroup", "site"]).mismatches.median().unstack()[ORDER]
    pl = papg_long()
    n_g = pl.groupby("allele").genome.nunique()
    stat = pl.groupby(["allele", "site"]).mismatches.median().unstack()[ORDER]
    classes = [(a, l) for a, l in PAPG if a in n_g]
    vmax = int(max(med.loc[keep].max().max(), stat.max().max()))
    norm = Normalize(0, vmax)
    gu = [(g, f"{g} (n={n_u[g]:,})") for g in keep]
    gp = [(a, f"{l} (n={n_g[a]:,})") for a, l in classes]
    vu = {g: med.loc[g].to_dict() for g in keep}
    vp = {a: stat.loc[a].to_dict() for a, _ in classes}

    ROW, TRK, GENE, ZOOM, TIGHT, BETWEEN = 7.5, 24.0, 17.0, 10.0, 0.5, 26.0   # mm; ZOOM = room for the zoom lines
    h_u, h_p = GENE + ZOOM + TRK + TIGHT + ROW * len(gu), GENE + ZOOM + TRK + TIGHT + ROW * len(gp)
    fig = plt.figure(figsize=(180 / 25.4, (h_u + h_p + BETWEEN + 62) / 25.4))
    outer = fig.add_gridspec(3, 1, height_ratios=[h_u, BETWEEN, h_p], hspace=0)
    gs_u = outer[0].subgridspec(5, 1, height_ratios=[GENE, ZOOM, TRK, TIGHT, ROW * len(gu)], hspace=0)
    gs_p = outer[2].subgridspec(5, 1, height_ratios=[GENE, ZOOM, TRK, TIGHT, ROW * len(gp)], hspace=0)
    block(fig, gs_u, u, xmax, gu, vu, norm, "uidA",
          f"{u.end.max()}-bp amplicon (nt {u_pos[0]}–{u_pos[1]}) of the {u_len}-bp gene · variants: UPEC + NCBI genomes (n = 6,231)", "uidA")
    block(fig, gs_p, p, xmax, gp, vp, norm, "papGII",
          f"{p.end.max()}-bp amplicon (nt {p_pos[0]}–{p_pos[1]}) of the {p_len}-bp gene · variants: papGII genomes of UPEC + NCBI (n = 645); axis reversed", "papGII")
    fig.subplots_adjust(left=0.17, right=0.98, top=0.96, bottom=0.085)
    cax = fig.add_axes([0.37, 0.035, 0.26, 0.007])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cax, orientation="horizontal", ticks=range(0, vmax + 1, 1))
    cb.outline.set_visible(False); cb.ax.tick_params(labelsize=6, colors=INK2, length=2)
    cb.set_label("Median number of mismatches to the primer (UPEC + NCBI genomes pooled)", fontsize=6.8, color=INK2)
    draw_zoom_lines(fig)
    for ext in ("pdf", "png", "svg"):
        meta = {"pdf": {"CreationDate": None, "Creator": "fig_primer_heatmap.py"}, "svg": {"Date": None}}.get(ext)
        fig.savefig(OUT / f"fig2a_primer_heatmap.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight", metadata=meta)


if __name__ == "__main__":
    main()
