"""Figure 2A (lollipop version): LAMP primer binding sites and median mismatches to the primers.

One lollipop plot per target, with the binding sites as a track below (as in mutation-lollipop plots):
- uidA: median number of mismatches to the primer at each binding site, one lollipop per Clermont phylogroup
  (UPEC collection + NCBI RefSeq complete genomes pooled; genomes with the amplicon region present).
- papGII: median number of mismatches at the best-matching papG locus of each genome, one lollipop per papG class
  (both collections pooled).
Groups differ by colour and marker shape. Both collections are also reported separately in the output tables.

Output: results/figures/fig2a_primer_lollipop.{pdf,png,svg}; fig2a_lollipop_{uidA,papGII}.tsv (pooled, with n, and
range for papGII); fig2a_lollipop_{uidA,papGII}_by_collection.tsv
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon

from fig_primer_map import MIN_GROUP, NAMED, XPAD, load_site_mismatches, sites_for
from lamp_common import DATA, RESULTS

OUT = RESULTS / "figures"
ORDER = ["F3", "F2", "LF", "F1c", "B1c", "LB", "B2", "B3"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e1e0d9"
ROLE_TINT = {"outer": "#cfe0f7", "inner": "#fbd9cb", "loop": "#cdeee1"}          # pale tints of the role colours
ROLE_LABEL = {"outer": "Outer (F3, B3)", "inner": "Inner (F2, F1c, B1c, B2)", "loop": "Loop (LF, LB)"}
# categorical palette (fixed order) + marker shapes
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
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


def track(ax, st, y0=0.0, h=1.0):
    """Binding sites as pointed boxes (arrow tip = primer direction on the amplicon), labelled, tinted by primer role."""
    for _, r in st.iterrows():
        tip = min(3.0, (r.end - r.start + 1) * 0.25)
        if r.orientation == "+":
            pts = [(r.start - 1, y0), (r.end - tip, y0), (r.end, y0 + h / 2), (r.end - tip, y0 + h), (r.start - 1, y0 + h)]
        else:
            pts = [(r.end, y0), (r.start - 1 + tip, y0), (r.start - 1, y0 + h / 2), (r.start - 1 + tip, y0 + h), (r.end, y0 + h)]
        ax.add_patch(Polygon(pts, closed=True, fc=ROLE_TINT[r.role], ec="#898781", lw=0.5))
        ax.text((r.start - 1 + r.end) / 2, y0 + h / 2, r.site, ha="center", va="center", fontsize=6.8, color=INK, fontweight="bold")


def lollipop_panel(fig, gs, st, xmax, groups, values, ymax, yticks, title, sub, ylabel):
    """groups: [(key, label, color, marker)], values: {key: {site: y}}"""
    ax = fig.add_subplot(gs[0]); axt = fig.add_subplot(gs[1], sharex=ax)
    k = len(groups)
    for _, r in st.iterrows():
        ax.axvspan(r.start - 1, r.end, color="#f4f3ef", lw=0, zorder=0)
        w = r.end - r.start + 1
        offs = np.linspace(-0.36 * w, 0.36 * w, k) if k > 1 else [0]
        for (key, label, c, m), off in zip(groups, offs):
            y, x = values[key][r.site], (r.start - 1 + r.end) / 2 + off
            if y > 0:
                ax.plot([x, x], [0, y], color=c, lw=1.0, zorder=2)
            ax.plot(x, y, marker=m, color=c, ms=4.0 if y > 0 else 2.4, alpha=1 if y > 0 else 0.55, mec="white",
                    mew=0.4 if y > 0 else 0, zorder=3, ls="none")
    ax.set_ylim(-ymax * 0.05, ymax * 1.06); ax.set_xlim(-XPAD, xmax)
    ax.set_yticks(yticks); ax.tick_params(labelsize=7, colors=INK2)
    ax.set_ylabel(ylabel, fontsize=7, color=INK2)
    ax.grid(axis="y", color=GRID, lw=0.5, zorder=0); ax.set_axisbelow(True)
    for sp in ("top", "right", "bottom"): ax.spines[sp].set_visible(False)
    ax.spines["left"].set_color("#c3c2b7"); ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.set_title(title, loc="left", fontsize=9, fontstyle="italic", color=INK, pad=37)
    ax.text(1.0, 1.215, sub, transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5, color=INK2)
    track(axt, st); axt.set_ylim(-0.1, 1.1)
    axt.set_yticks([]); axt.set_xlim(-XPAD, xmax)
    for sp in ("top", "right", "left"): axt.spines[sp].set_visible(False)
    axt.spines["bottom"].set_color("#c3c2b7"); axt.tick_params(labelsize=7, colors=INK2)
    return ax, axt


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    u, u_pos, u_len = sites_for("uidA")
    p, p_pos, p_len = sites_for("papGII")
    xmax = max(u.end.max(), p.end.max()) + 6

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

    fig = plt.figure(figsize=(180 / 25.4, 168 / 25.4))
    outer = fig.add_gridspec(2, 1, hspace=0.62)
    gs_u = outer[0].subgridspec(2, 1, height_ratios=[1, 0.16], hspace=0.04)
    gs_p = outer[1].subgridspec(2, 1, height_ratios=[1, 0.16], hspace=0.04)
    ax_u, axt_u = lollipop_panel(fig, gs_u, u, xmax, gu, vu, 3, [0, 1, 2, 3], "uidA",
                                 f"{u.end.max()}-bp amplicon · nt {u_pos[0]}–{u_pos[1]} of {u_len}-bp gene", "Median mismatches\nto the primer")
    ax_p, axt_p = lollipop_panel(fig, gs_p, p, xmax, gp, vp, 12, [0, 3, 6, 9, 12], "papGII",
                                 f"{p.end.max()}-bp amplicon · nt {p_pos[0]}–{p_pos[1]} of {p_len}-bp gene", "Median mismatches\nto the primer")
    axt_u.set_xlabel("Position in amplicon (nt)", fontsize=7, color=INK2)
    axt_p.set_xlabel("Position in amplicon (nt)", fontsize=7, color=INK2)
    for ax, groups, ncol, name in ((ax_u, gu, 4, "Clermont phylogroup"), (ax_p, gp, 4, "papG class")):
        handles = [Line2D([], [], marker=m, color=c, ls="none", ms=4.6, label=l) for _, l, c, m in groups]
        leg = ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=ncol, frameon=False, fontsize=6.8,
                        handletextpad=0.2, columnspacing=1.0, borderaxespad=0.2, title=name, title_fontsize=6.5)
        leg._legend_box.align = "left"
    role_handles = [Line2D([], [], marker="s", color=ROLE_TINT[k], mec="#898781", mew=0.5, ls="none", ms=6, label=ROLE_LABEL[k])
                    for k in ROLE_TINT]
    fig.legend(handles=role_handles, loc="lower center", ncol=3, frameon=False, fontsize=7, bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(left=0.1, right=0.98, top=0.94, bottom=0.095)
    for ext in ("pdf", "png", "svg"):
        meta = {"pdf": {"CreationDate": None, "Creator": "fig_primer_lollipop.py"}, "svg": {"Date": None}}.get(ext)
        fig.savefig(OUT / f"fig2a_primer_lollipop.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight", metadata=meta)
    print(med.loc[keep].assign(n=n_u[keep]).to_string()); print(stat.pivot(index="allele", columns="site", values="median")[ORDER].to_string())


if __name__ == "__main__":
    main()
