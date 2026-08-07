#!/usr/bin/env python3
"""
06_figures.py — main display items for the manuscript.

Reads the tables produced by 04 (source stats) and 05 (inequity), plus feasibility
tables, and renders 5 colourblind-safe (Okabe-Ito), 300-dpi figures to figures/
as both .png and .pdf.

  Fig 1  Global visibility inequity (records/million + country coverage by income)
  Fig 2  One Health source composition of deposited genomes
  Fig 3  AMR/virulence gene prevalence heatmap (gene x source)
  Fig 4  Host-differentiation forest plot (livestock vs human odds ratios)
  Fig 5  Mixed-source SNP cluster overlap at the One Health interface
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(__file__).resolve().parent.parent
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"
FIGS.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "Arial", "font.size": 9, "axes.titlesize": 10,
    "axes.titleweight": "bold", "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
})

# Okabe-Ito
OI = {"blue": "#0072B2", "sky": "#56B4E9", "green": "#009E73", "orange": "#E69F00",
      "vermillion": "#D55E00", "purple": "#CC79A7", "yellow": "#F0E442", "grey": "#999999"}
INCOME = ["High income", "Upper middle income", "Lower middle income", "Low income"]
INCOME_SHORT = ["High", "Upper-mid", "Lower-mid", "Low"]
INCOME_COL = [OI["blue"], OI["sky"], OI["orange"], OI["vermillion"]]
SRC_COL = {"human_or_clinical": OI["blue"], "livestock": OI["vermillion"],
           "companion_animal": OI["green"], "food": OI["orange"],
           "environment": OI["sky"], "unknown": OI["grey"]}


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(FIGS / f"{name}.{ext}")
    plt.close(fig)
    print(f"  wrote figures/{name}.png / .pdf")


def fig1_visibility():
    vis = pd.read_csv(TABLES / "30_visibility_by_income.csv")
    cov = pd.read_csv(TABLES / "36_country_coverage_by_income.csv")
    vis = vis.set_index("income_group").reindex(INCOME)
    cov = cov.set_index("income_group").reindex(INCOME)

    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.2))
    x = np.arange(len(INCOME))
    a.bar(x, vis["records_per_million"], color=INCOME_COL, width=0.7)
    a.set_yscale("log")
    a.set_ylabel("Genomes per million population\n(log scale)")
    a.set_xticks(x); a.set_xticklabels(INCOME_SHORT)
    a.set_title("A  Public genome visibility")
    for xi, v in zip(x, vis["records_per_million"]):
        a.text(xi, v * 1.15, f"{v:.2f}", ha="center", va="bottom", fontsize=8)
    a.set_ylim(0.2, vis["records_per_million"].max() * 3)

    b.bar(x, cov["coverage_pct"], color=INCOME_COL, width=0.7)
    b.set_ylabel("Countries with ≥1 deposited genome (%)")
    b.set_xticks(x); b.set_xticklabels(INCOME_SHORT)
    b.set_title("B  Country-level coverage")
    b.set_ylim(0, 100)
    for xi, v, w, t in zip(x, cov["coverage_pct"], cov["countries_with_genome"], cov["countries_total"]):
        b.text(xi, v + 2, f"{v:.0f}%\n({w}/{t})", ha="center", va="bottom", fontsize=7.5)
    fig.suptitle("Global inequity in public MRSA/S. aureus genomic sequence visibility, 2015–2025",
                 fontsize=10, fontweight="bold", y=1.03)
    fig.tight_layout()
    save(fig, "fig1_visibility_inequity")


def fig2_sources():
    vis = pd.read_csv(TABLES / "30_visibility_by_income.csv").set_index("income_group").reindex(INCOME)
    # panel A: overall source-category counts (from feasibility table 03/ processed)
    src = pd.read_csv(TABLES / "03_source_mrsa_summary.csv")
    src = src[["source_category", "total_records"]].rename(
        columns={"source_category": "src", "total_records": "n"})
    order = ["human_or_clinical", "unknown", "livestock", "food", "companion_animal", "environment"]
    src = src.set_index("src").reindex(order).dropna()

    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.2))
    ypos = np.arange(len(src))[::-1]
    a.barh(ypos, src["n"], color=[SRC_COL[s] for s in src.index])
    a.set_xscale("log")
    a.set_yticks(ypos); a.set_yticklabels([s.replace("_", " ") for s in src.index])
    a.set_xlabel("Deposited genomes (log scale)")
    a.set_title("A  Source composition")
    for yi, v in zip(ypos, src["n"]):
        a.text(v * 1.1, yi, f"{int(v):,}", va="center", fontsize=7.5)

    # panel B: non-human genomes per million by income (One Health visibility)
    x = np.arange(len(INCOME))
    b.bar(x, vis["non_human_per_million"], color=INCOME_COL, width=0.7)
    b.set_yscale("log")
    b.set_ylabel("Non-human genomes per million\n(log scale)")
    b.set_xticks(x); b.set_xticklabels(INCOME_SHORT)
    b.set_title("B  One Health interface visibility")
    for xi, v in zip(x, vis["non_human_per_million"]):
        b.text(xi, v * 1.15, f"{v:.3f}", ha="center", va="bottom", fontsize=7.5)
    b.set_ylim(0.02, vis["non_human_per_million"].max() * 4)
    fig.tight_layout()
    save(fig, "fig2_source_composition")


def fig3_heatmap():
    prev = pd.read_csv(TABLES / "26_gene_prevalence_by_source.csv").set_index("gene")
    genes = ["mecA (methicillin)", "blaZ (penicillinase)", "any acquired tetracycline [tet(K/L/M/O)]",
             "tet(M)", "any MLSb [erm(*)]", "msr(A)", "mph(C)", "any aminoglycoside [aac/ant/aph]",
             "dfrG (trimethoprim)", "gyrA/parC QRDR mutation (fluoroquinolone)", "fusC (fusidic acid)",
             "PVL [lukS-PV/lukF-PV]", "scn (IEC / staph. complement inhibitor)",
             "sak (IEC / staphylokinase)", "tst (TSST-1)", "sea (enterotoxin A)", "cna (collagen adhesin)"]
    srcs = ["human_or_clinical", "livestock", "companion_animal", "food", "environment"]
    M = np.array([[prev.loc[g, f"{s}_pct"] for s in srcs] for g in genes], dtype=float)
    short = [g.split(" (")[0].split(" [")[0] for g in genes]

    cmap = LinearSegmentedColormap.from_list("v", ["#f7fbff", OI["sky"], OI["blue"], "#08306b"])
    fig, ax = plt.subplots(figsize=(5.2, 6.2))
    im = ax.imshow(M, cmap=cmap, aspect="auto", vmin=0, vmax=100)
    ax.set_xticks(range(len(srcs)))
    ax.set_xticklabels(["Human", "Livestock", "Companion", "Food", "Environment"], rotation=30, ha="right")
    ax.set_yticks(range(len(genes))); ax.set_yticklabels(short, fontsize=8)
    ax.axhline(10.5, color="k", lw=1.2)  # AMR | virulence divider
    ax.text(-1.9, 5, "AMR", rotation=90, va="center", fontweight="bold", fontsize=9)
    ax.text(-1.9, 14, "Virulence", rotation=90, va="center", fontweight="bold", fontsize=9)
    for i in range(len(genes)):
        for j in range(len(srcs)):
            ax.text(j, i, f"{M[i,j]:.0f}", ha="center", va="center", fontsize=6.5,
                    color="white" if M[i, j] > 55 else "black")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Prevalence in deposited genomes (%)")
    ax.set_title("Gene prevalence across One Health source categories", fontsize=9.5)
    fig.tight_layout()
    save(fig, "fig3_gene_heatmap")


def fig4_forest():
    orr = pd.read_csv(TABLES / "28_gene_or_vs_human.csv")
    liv = orr[orr["source"] == "livestock"].set_index("gene")
    genes = ["scn (IEC / staph. complement inhibitor)", "sak (IEC / staphylokinase)",
             "PVL [lukS-PV/lukF-PV]", "msr(A)", "mph(C)", "fusC (fusidic acid)",
             "sea (enterotoxin A)", "any MLSb [erm(*)]", "dfrG (trimethoprim)",
             "cna (collagen adhesin)", "tet(M)", "any acquired tetracycline [tet(K/L/M/O)]"]
    genes = [g for g in genes if g in liv.index]
    d = liv.loc[genes]
    short = [g.split(" (")[0].split(" [")[0] for g in genes]
    y = np.arange(len(genes))[::-1]

    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    for yi, (_, r) in zip(y, d.iterrows()):
        col = OI["vermillion"] if r["odds_ratio_vs_human"] > 1 else OI["blue"]
        ax.plot([r["ci95_low"], r["ci95_high"]], [yi, yi], color=col, lw=1.6, zorder=1)
        ax.scatter(r["odds_ratio_vs_human"], yi, color=col, s=26, zorder=2)
    ax.axvline(1, color="k", ls="--", lw=0.9)
    ax.set_xscale("log")
    ax.set_yticks(y); ax.set_yticklabels(short, fontsize=8)
    ax.set_xlabel("Odds ratio, livestock vs human (log scale)")
    ax.set_title("Host differentiation of gene carriage\n(livestock vs human reservoirs)", fontsize=9.5)
    ax.text(0.02, 0.02, "lower in livestock", transform=ax.transAxes, color=OI["blue"], fontsize=7.5)
    ax.text(0.98, 0.02, "higher in livestock", transform=ax.transAxes, color=OI["vermillion"],
            fontsize=7.5, ha="right")
    fig.tight_layout()
    save(fig, "fig4_host_forest")


def fig5_clusters():
    c = pd.read_csv(TABLES / "14_cross_interface_cluster_patterns.csv")
    labels = {"human_livestock": "Human–Livestock", "human_companion": "Human–Companion",
              "human_food": "Human–Food", "human_environment": "Human–Environment",
              "livestock_food": "Livestock–Food", "livestock_environment": "Livestock–Environment",
              "three_or_more_sources": "≥3 source categories"}
    c["label"] = c["pattern"].map(labels)
    c = c.sort_values("records_in_clusters")
    y = np.arange(len(c))
    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    ax.barh(y, c["records_in_clusters"], color=OI["green"])
    ax.set_yticks(y); ax.set_yticklabels(c["label"])
    ax.set_xlabel("Genomes in mixed-source SNP clusters")
    ax.set_title("Genomic overlap across the One Health interface", fontsize=9.5)
    for yi, rec, cl in zip(y, c["records_in_clusters"], c["clusters"]):
        ax.text(rec + 20, yi, f"{int(rec):,} ({int(cl)} clusters)", va="center", fontsize=7.5)
    ax.set_xlim(0, c["records_in_clusters"].max() * 1.25)
    ax.text(0.98, 0.04, "Genomic overlap only — not proof of direct transmission",
            transform=ax.transAxes, ha="right", fontsize=7, style="italic", color=OI["grey"])
    fig.tight_layout()
    save(fig, "fig5_cluster_overlap")


def main():
    print("Rendering figures ->", FIGS)
    fig1_visibility()
    fig2_sources()
    fig3_heatmap()
    fig4_forest()
    fig5_clusters()
    print("Done.")


if __name__ == "__main__":
    main()
