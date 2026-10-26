"""Figure 2A: LAMP primer binding sites on uidA and papGII with in silico conservation (results/gene_scan).

Per target, three aligned rows:
1. Gene track: % of genomes with a variant at each position of the reference gene (700223-20), amplicon shaded.
   For papGII the amplicon lies on the reverse strand, so the gene axis runs right-to-left (F3 -> B3 direction).
2. Amplicon map: binding sites as arrows (forward above, reverse below), coloured by primer role.
   papGII: mismatches of each site against the PapG-III reference allele next to the site label.
3. Per site: % of genomes with 1, 2 or >= 3 mismatches (grey, light -> dark) and with a mismatch in the
   3'-terminal 5 nt (diamond). Genomes: uidA - all with the amplicon covered; papGII - genomes typed papGII.
Binding-site coordinates follow the NEB design (FIP/BIP split lengths in primer_screen.SPLIT).

Output: results/figures/fig2a_primer_map.{pdf,png,svg}, fig2a_site_mismatches.tsv, fig2a_genome_status_*.tsv
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, FancyArrow, Patch

from lamp_common import DATA, RESULTS
from primer_screen import INITIATING, ROLE, THREE_PRIME_SITES, binding_sites, rc, read_fasta, reference_amplicon

SCAN = RESULTS / "gene_scan"
OUT = RESULTS / "figures"
ROLE_COLOR = {"outer": "#2a78d6", "inner": "#eb6834", "loop": "#1baf7a"}   # categorical slots 1-3
ROLE_LABEL = {"outer": "Outer (F3, B3)", "inner": "Inner (FIP = F1c+F2, BIP = B1c+B2)", "loop": "Loop (LF, LB)"}
MM_COLOR = {"mm1": "#c3c2b7", "mm2": "#898781", "mm3plus": "#52514e"}     # neutral ordinal ramp
INK, INK2, GRID, SHADE = "#0b0b0b", "#52514e", "#e1e0d9", "#f0efec"
END_NT, XPAD = 5, 4

mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 8, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "pdf.fonttype": 42, "svg.fonttype": "none"})


def best_mismatches(q, text):
    return min(sum(a != b for a, b in zip(qq, text[i:i + len(qq)]))
               for qq in (q, rc(q)) for i in range(len(text) - len(qq) + 1))


def load(target):
    """Sites in amplicon coordinates, per-genome site mismatches and statuses, per-position variability."""
    d = SCAN / target
    gene = next(iter(read_fasta(d / "reference_gene.fasta").values()))
    L = len(gene)
    amp, offset = reference_amplicon(target, {"gene": gene})
    sites, primer_amp = binding_sites(target, amp)
    j, strand = gene.find(amp), "+"
    if j < 0: j, strand = rc(gene).find(amp), "-"
    a0 = j + offset  # 0-based start of F3 on the strand carrying the amplicon
    to_amp = (lambda g: g - a0) if strand == "+" else (lambda g: (L - g + 1) - a0)
    st = pd.DataFrame([dict(site=k, role=ROLE[k], start=s - offset + 1, end=e - offset, orientation=o, seq=seq)
                       for k, (seq, s, e, o, _) in sites.items()]).sort_values("start")
    amp_len = int(st.end.max())
    gpos = [a0 + 1, a0 + amp_len] if strand == "+" else [L - a0 - amp_len + 1, L - a0]
    gene_lo, gene_hi = min(gpos), max(gpos)
    cov = pd.read_csv(d / "genome_coverage.tsv", sep="\t")
    var = pd.read_csv(d / "genome_variants.tsv", sep="\t", keep_default_na=False)
    var["apos"] = var.position.map(to_amp)
    if strand == "-":  # alternative bases are reported on the gene strand
        var["alt"] = var.alt.map(lambda b: rc(b) if isinstance(b, str) and b not in "-" else b)
    # positions where the primer differs from the reference isolate: a genome mismatches the primer there unless it
    # carries the primer base (variants are called against the reference isolate, not against the primer)
    diff = {k - offset + 1: (amp[k], primer_amp[k]) for k in range(len(amp)) if amp[k] != primer_amp[k]}
    vpos, vsub = {}, {}
    for gid, grp in var.groupby("genome"):
        vpos[gid] = [p for p, typ in zip(grp.apos, grp.type) if p not in diff or typ != "sub"]
        vsub[gid] = {p: a for p, typ, a in zip(grp.apos, grp.type, grp.alt) if typ == "sub" and p in diff}
    for gid in cov.genome:
        vpos.setdefault(gid, []); vsub.setdefault(gid, {})
        vpos[gid] += [p for p, (ref_b, prim_b) in diff.items() if vsub[gid].get(p, ref_b) != prim_b]
    rows, status = [], []
    for _, c in cov.iterrows():
        if pd.isna(c.qstart):
            status.append(dict(genome=c.genome, status="absent")); continue
        if c.qstart > gene_lo or c.qend < gene_hi:
            status.append(dict(genome=c.genome, status="incomplete")); continue
        P, crit, anymm = vpos.get(c.genome, []), False, False
        for _, r in st.iterrows():
            n = sum(r.start <= p <= r.end for p in P)
            zone = range(r.end - END_NT + 1, r.end + 1) if r.orientation == "+" else range(r.start, r.start + END_NT)
            n3 = sum(p in zone for p in P) if r.site in THREE_PRIME_SITES else 0
            rows.append(dict(genome=c.genome, site=r.site, mismatches=n, mismatches_3prime=n3))
            anymm |= n > 0
            crit |= (r.site in INITIATING and n3 > 0) or (r.role == "inner" and n >= 2)
        status.append(dict(genome=c.genome, status="critical" if crit else "minor" if anymm else "intact"))
    mm = pd.DataFrame(rows)
    g = mm.groupby("site")
    agg = pd.DataFrame({"n": g.size(), "mm1": g.mismatches.apply(lambda x: (x == 1).sum()),
                        "mm2": g.mismatches.apply(lambda x: (x == 2).sum()),
                        "mm3plus": g.mismatches.apply(lambda x: (x >= 3).sum()),
                        "end3": g.mismatches_3prime.apply(lambda x: (x > 0).sum())}).reset_index()
    st = st.merge(agg, on="site")
    for c in ("mm1", "mm2", "mm3plus", "end3"): st[f"pct_{c}"] = 100 * st[c] / st.n
    pv = pd.read_csv(d / "position_variability.tsv", sep="\t")
    pv["pct"] = 100 * (pv.n_substitution + pv.n_gap) / pv.n_covered.where(pv.n_covered > 0)
    return dict(L=L, strand=strand, sites=st, amp_len=amp_len, gene_lo=gene_lo, gene_hi=gene_hi,
                status=pd.DataFrame(status), pv=pv)


def draw_gene(ax, t, name):
    pv = t["pv"]
    ax.axvspan(t["gene_lo"], t["gene_hi"], color=SHADE, zorder=0, lw=0)
    ax.bar(pv.position, pv.pct, width=1.0, color=INK2, lw=0, zorder=2)
    ax.set_xlim((t["L"] + 1, 0) if t["strand"] == "-" else (0, t["L"] + 1))
    ax.set_ylim(0, 100); ax.set_yticks([0, 50, 100])
    ax.set_ylabel("Variant vs\n700223-20 (%)", fontsize=7)
    if t["strand"] == "-":
        ax.text(0.075, 1.06, "gene axis reversed (primers on the reverse strand)", transform=ax.transAxes,
                ha="left", va="bottom", fontsize=6.5, color=INK2)
    ax.set_xlabel(None)
    ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=6.5, pad=1)


def draw_map(ax, st, xmax, papg3=None):
    ax.plot([1, st.end.max()], [0, 0], color=INK2, lw=1.2, solid_capstyle="butt", zorder=1)
    for _, r in st.iterrows():
        y = 0.55 if r.orientation == "+" else -0.55
        x0, x1 = (r.start - 1, r.end) if r.orientation == "+" else (r.end, r.start - 1)
        ax.add_patch(FancyArrow(x0, y, x1 - x0, 0, width=0.32, head_width=0.32, length_includes_head=True,
                                head_length=min(6, abs(x1 - x0) * 0.35), color=ROLE_COLOR[r.role], lw=0, zorder=2))
        up = y > 0
        label = r.site if papg3 is None else f"{r.site}  ({papg3[r.site]})"
        ax.text((r.start + r.end) / 2, y + (0.4 if up else -0.4), label, ha="center", va="bottom" if up else "top",
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
        if r.pct_end3 > 0:
            ax.plot(x, r.pct_end3, marker="D", ms=3.5, color=INK, markeredgecolor="white", markeredgewidth=0.6,
                    zorder=3, lw=0)
    ax.set_xlim(-XPAD, xmax); ax.set_ylim(0, 60); ax.set_yticks([0, 20, 40, 60])
    ax.set_ylabel("Genomes with\nmismatches (%)", fontsize=7)
    if xlabel: ax.set_xlabel("Position in amplicon (nt)", fontsize=7)
    ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)
    ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=6.5)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    alleles = read_fasta(DATA / "papG_reference_alleles.fasta")
    data = {t: load(t) for t in ("uidA", "papGII")}
    xmax = max(t["amp_len"] for t in data.values()) + 8
    fig = plt.figure(figsize=(180 / 25.4, 150 / 25.4))
    gs = fig.add_gridspec(7, 1, height_ratios=[0.6, 1.0, 0.9, 0.3, 0.6, 1.0, 0.9], hspace=0.45)
    tables = []
    for i, target in enumerate(("uidA", "papGII")):
        t, st = data[target], data[target]["sites"]
        ax_g, ax_m, ax_b = (fig.add_subplot(gs[4 * i + k]) for k in range(3))
        draw_gene(ax_g, t, target)
        papg3 = {r.site: best_mismatches(r.seq, alleles["PapG-III"]) for _, r in st.iterrows()} if target == "papGII" else None
        draw_map(ax_m, st, xmax, papg3)
        draw_bars(ax_b, st, xmax, xlabel=(target == "papGII"))
        # zoom lines: amplicon ends in the gene track -> amplicon map
        ends = ((t["gene_lo"], 1), (t["gene_hi"], t["amp_len"])) if t["strand"] == "+" else \
               ((t["gene_hi"], 1), (t["gene_lo"], t["amp_len"]))
        for xg, xa in ends:
            fig.add_artist(ConnectionPatch((xg, -0.38), (xa, 1.5), ax_g.get_xaxis_transform(), "data",
                                           axesA=ax_g, axesB=ax_m, color="#c3c2b7", lw=0.6, ls="--"))
        s = t["status"].status.value_counts()
        n = len(t["status"]); ok = s.get("intact", 0) + s.get("minor", 0)
        who = "UPEC genomes" if target == "uidA" else "papGII genomes"
        ax_g.set_title(target, loc="left", fontsize=9, fontstyle="italic", color=INK, pad=3)
        ax_g.text(1.0, 1.06, f"{t['amp_len']}-bp amplicon (nt {t['gene_lo']}–{t['gene_hi']}) · no critical mismatch in "
                  f"{ok}/{n} {who} ({100 * ok / n:.1f}%)", transform=ax_g.transAxes, ha="right", va="bottom",
                  fontsize=6.5, color=INK2)
        if papg3:
            ax_m.text(1.0, -0.02, "(n) = mismatches to the PapG-III allele", transform=ax_m.transAxes, ha="right",
                      va="top", fontsize=6, color=INK2)
        tables.append(st.assign(target=target, papGIII_mismatches=st.site.map(papg3) if papg3 else None))
        t["status"].assign(target=target).to_csv(OUT / f"fig2a_genome_status_{target}.tsv", sep="\t", index=False)

    handles = [Patch(color=ROLE_COLOR[k], label=ROLE_LABEL[k]) for k in ROLE_COLOR]
    handles += [Patch(color=MM_COLOR["mm1"], label="1 mismatch"), Patch(color=MM_COLOR["mm2"], label="2 mismatches"),
                Patch(color=MM_COLOR["mm3plus"], label="≥3 mismatches"),
                Line2D([], [], marker="D", ms=4, color=INK, lw=0, label="Mismatch in 3′-terminal 5 nt")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=7, bbox_to_anchor=(0.5, -0.03))
    fig.subplots_adjust(left=0.09, right=0.98, top=0.96, bottom=0.11)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(OUT / f"fig2a_primer_map.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight")
    tab = pd.concat(tables)
    tab.to_csv(OUT / "fig2a_site_mismatches.tsv", sep="\t", index=False)
    for target in data:
        print(target, data[target]["status"].status.value_counts().to_dict(), "amplicon nt",
              data[target]["gene_lo"], "-", data[target]["gene_hi"])
    print(tab[["target", "site", "start", "end", "orientation", "n", "pct_mm1", "pct_mm2", "pct_mm3plus", "pct_end3",
               "papGIII_mismatches"]].round(1).to_string(index=False))


if __name__ == "__main__":
    main()
