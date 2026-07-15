# zmRBA predicted vs experimental protein abundance scatterplot

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

prot = pd.read_excel('RBA_Protein_Syn_Analysis.xlsx')
subsystem = pd.read_excel('GPR_Subsytem.xlsx')    
ribo = pd.read_excel('RIBOSOME_ZM.xlsx')
MU  = 0.36                              # growth rate (h-1)

  
#protein data
prot = prot[(prot.Exp_PROSYN != 0) & (prot.Model_PROSYN != 0)].copy()
prot["Exp_abund"]   = prot.Exp_PROSYN   / MU    # synthesis rate -> abundance (mmol gDCW^-1)
prot["Model_abund"] = prot.Model_PROSYN / MU

#mapping genes/proteins to subsystems (as in iZM489, inherited from iZM487)
gene = subsystem.dropna().copy()
gene["GPR"] = gene["GPR"].astype(str)           # so multiple gene entries can be split
for sep in ("or", "and"):
    gene["GPR"] = gene["GPR"].str.split(sep); gene = gene.explode("GPR")
gene["GPR"] = (gene["GPR"].str.replace("(", "", regex=False)
                    .str.replace(")", "", regex=False).str.strip())
gene = gene[gene["GPR"] != ""].drop_duplicates()
gene2sub = gene.drop_duplicates("GPR").set_index("GPR")["Subsystem"]
prot["Subsystem"] = prot["Gene"].map(gene2sub).fillna("Unassigned")

#ribosomal proteins (not in the GSM subsystem assignments)
ribo_genes = set()
for col in ribo.columns:
    v = ribo[col].astype(str)
    ribo_genes |= set(x.strip() for x in v[v.str.contains("ZMO")])

#collapsing the BiGG subsystems into broader categories
BROAD_PRO_CATEGORIES = {
    "Central carbon metabolism": [
        "Glycolysis/Gluconeogenesis", "Pentose Phosphate Pathway", "Citric Acid Cycle",
        "Alternate Carbon Metabolism", "Pyruvate Metabolism", "Anaplerotic Reactions",
    ],
    "Energy / oxidative phosphorylation": [
        "Oxidative Phosphorylation",
    ],
    "Amino acid metabolism": [
        "Tyrosine, Tryptophan, and Phenylalanine Metabolism", "Arginine and Proline Metabolism",
        "Histidine Metabolism", "Threonine and Lysine Metabolism",
        "Valine, Leucine, and Isoleucine Metabolism", "Cysteine Metabolism",
        "Methionine Metabolism", "Glutamate Metabolism", "Alanine and Aspartate Metabolism",
        "Glycine and Serine Metabolism",
    ],
    "Nucleotide metabolism": [
        "Purine and Pyrimidine Biosynthesis", "Nucleotide Salvage Pathway",
    ],
    "Cofactor & vitamin biosynthesis": [
        "Cofactor and Prosthetic Group Biosynthesis", "Folate Metabolism",
    ],
    "Cell envelope & lipid": [
        "Cell Envelope Biosynthesis", "Glycerophospholipid Metabolism", "Membrane Lipid Metabolism",
    ],
    "Transport": [
        "Inorganic Ion Transport and Metabolism", "Transport, Inner Membrane",
        "Transport, Outer Membrane Porin", "Transport, Outer Membrane",
    ],
}

TRANSLATION = "Translation & ribosome"
OTHER = "Other / unassigned"


# order of categories on the legend
ORDER = [
    "Central carbon metabolism", "Energy / oxidative phosphorylation", "Amino acid metabolism",
    "Nucleotide metabolism", "Cofactor & vitamin biosynthesis", "Cell envelope & lipid",
    "Translation & ribosome", "Transport", "Other / unassigned",
]


sub2broad = {s: b for b, subs in BROAD_PRO_CATEGORIES.items() for s in subs}
def categorize(row):
    if row["Gene"] in ribo_genes or row["Subsystem"] == "tRNA Charging":
        return TRANSLATION
    return sub2broad.get(row["Subsystem"], OTHER)
prot["Category"] = prot.apply(categorize, axis=1)

#stats
x, y = np.log10(prot.Exp_abund.values), np.log10(prot.Model_abund.values)
r2 = np.corrcoef(x, y)[0, 1] ** 2
within2 = np.mean(np.abs(x - y) <= np.log10(2)) * 100   # predictions within 2-fold of experiment

# plotting
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 13, "axes.linewidth": 1.1,
    "xtick.direction": "out", "ytick.direction": "out",
})
CB = {
    "Central carbon metabolism": "#0072B2",
    "Energy / oxidative phosphorylation": "#E69F00",
    "Amino acid metabolism": "#009E73",
    "Nucleotide metabolism": "#D55E00",
    "Cofactor & vitamin biosynthesis": "#CC79A7",
    "Cell envelope & lipid": "#8C564B",
    "Translation & ribosome": "#56B4E9",
    "Transport": "#000000",
    OTHER: "#BEBEBE",
}
fig, ax = plt.subplots(figsize=(8.4, 6))
present = [c for c in ORDER if (prot.Category == c).any()]
for c in present:
    sub = prot[prot.Category == c]
    ax.scatter(sub.Exp_abund, sub.Model_abund, s=46, alpha=0.85,
               edgecolors="white", linewidths=0.5, color=CB.get(c, "#BEBEBE"),
               label=f"{c} (n={len(sub)})", zorder=2 if c == OTHER else 3)

allv = np.concatenate([prot.Exp_abund.values, prot.Model_abund.values])
lo, hi = allv.min() * 0.6, allv.max() * 1.6
ax.plot([lo, hi], [lo, hi], "--", color="0.35", lw=1.3, zorder=1)

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
ax.spines[["top", "right"]].set_visible(False)
ax.tick_params(which="both", length=4, width=1.1)
ax.set_xlabel("Experimental protein abundance (mmol gDCW$^{-1}$)")
ax.set_ylabel("zmRBA predicted protein abundance (mmol gDCW$^{-1}$)")
ax.text(0.03, 0.97, f"$R^2$ = {r2:.2f}",
        transform=ax.transAxes, va="top", ha="left", fontsize=11)
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False,
          fontsize=9.5, title="Functional category", title_fontsize=10.5,
          handletextpad=0.3, labelspacing=0.35)
fig.tight_layout()

OUT = "Figure5_prosyn_scatter"  
for ext in ("pdf", "png", "svg"):
    fig.savefig(f"{OUT}.{ext}", dpi=300, bbox_inches="tight")

print(prot.Category.value_counts().to_string())
print(f"n={len(prot)}  R^2(log10)={r2:.3f}  within-2-fold={within2:.1f}%")

plt.show()  