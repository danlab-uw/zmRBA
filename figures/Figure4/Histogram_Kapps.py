import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#Apply Seaborn theme
sns.set_theme(style="ticks", context="talk")

#Load data
df = pd.read_excel('Kapp comparison and compilation.xlsx', sheet_name='Compiled Kapps')
df = df.dropna(subset=['kapp (1/s)', 'rxnid']).copy()

#Extract gene (enzyme) ID from rxnid
#rxnid format: "RXN-<rxn_name>_<direction>-<gene_id>"  e.g. "RXN-AHCi_FWD-ZMO0182"
#gene ID is whatever comes after the final '-'
df['gene_id'] = df['rxnid'].str.rsplit('-', n=1).str[-1]

#remove by gene_id
enzyme_df = df.drop_duplicates(subset='gene_id', keep='first').copy()

kapps_per_rxn     = df['kapp (1/s)']
kapps_per_enzyme  = enzyme_df['kapp (1/s)']

print(f"Total reaction entries:   {len(kapps_per_rxn)}")
print(f"Unique enzymes (genes):   {len(kapps_per_enzyme)}")
print(f"Per-reaction median kapp: {np.median(kapps_per_rxn):.3f} s^-1")
print(f"Per-enzyme median kapp:   {np.median(kapps_per_enzyme):.3f} s^-1")

#log-spaced bins 
bins = np.logspace(np.log10(kapps_per_rxn.min()), np.log10(kapps_per_rxn.max()), 25)


# per-enzyme distribution 
plt.figure(figsize=(10, 6))
plt.hist(kapps_per_enzyme, bins=bins, color="#388fd1",
         edgecolor='white', linewidth=1.2, alpha=0.9)

median_kapp = np.median(kapps_per_enzyme)
plt.axvline(median_kapp, color="#000000", linestyle='--', linewidth=3,
            label=f'Median: {median_kapp:.2f} s$^{{-1}}$ (n={len(kapps_per_enzyme)} enzymes)')

plt.xscale('log')
plt.xlim(1e-2, 1e4)
plt.tick_params(axis='x', which='minor', bottom=False)

plt.xlabel('Kapp (s$^{-1}$)', fontsize=14, labelpad=10)
plt.ylabel('Frequency', fontsize=14, labelpad=10)
plt.title('Distribution of Estimated Kapp Values (per enzyme)',
          fontsize=18, weight='bold', pad=15)

plt.grid(axis='y', linestyle='--', alpha=0.4)
sns.despine()
plt.legend(frameon=False, fontsize=12)

plt.tight_layout()
plt.savefig("Kapp_histogram_per_enzyme.pdf", bbox_inches='tight')
plt.show()


#side-by-side comparison of per-reaction vs per-enzyme
fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=False)

axes[0].hist(kapps_per_rxn, bins=bins, color="#c94f4f",
             edgecolor='white', linewidth=1.2, alpha=0.9)
axes[0].axvline(np.median(kapps_per_rxn), color='k', linestyle='--', linewidth=2.5,
                label=f'Median: {np.median(kapps_per_rxn):.2f} s$^{{-1}}$')
axes[0].set_xscale('log')
axes[0].set_xlim(1e-2, 1e4)
axes[0].tick_params(axis='x', which='minor', bottom=False)
axes[0].set_xlabel('Kapp (s$^{-1}$)', fontsize=14, labelpad=10)
axes[0].set_ylabel('Frequency', fontsize=14, labelpad=10)
axes[0].set_title(f'Per reaction (n={len(kapps_per_rxn)})', fontsize=15, weight='bold')
axes[0].grid(axis='y', linestyle='--', alpha=0.4)
axes[0].legend(frameon=False, fontsize=11)

axes[1].hist(kapps_per_enzyme, bins=bins, color="#388fd1",
             edgecolor='white', linewidth=1.2, alpha=0.9)
axes[1].axvline(np.median(kapps_per_enzyme), color='k', linestyle='--', linewidth=2.5,
                label=f'Median: {np.median(kapps_per_enzyme):.2f} s$^{{-1}}$')
axes[1].set_xscale('log')
axes[1].set_xlim(1e-2, 1e4)
axes[1].tick_params(axis='x', which='minor', bottom=False)
axes[1].set_xlabel('Kapp (s$^{-1}$)', fontsize=14, labelpad=10)
axes[1].set_ylabel('Frequency', fontsize=14, labelpad=10)
axes[1].set_title(f'Per unique enzyme (n={len(kapps_per_enzyme)})', fontsize=15, weight='bold')
axes[1].grid(axis='y', linestyle='--', alpha=0.4)
axes[1].legend(frameon=False, fontsize=11)

sns.despine(fig=fig)
plt.tight_layout()
plt.savefig("Kapp_histogram_comparison.pdf", bbox_inches='tight')
plt.show()



#extreme enzymes at the per-enzyme level
slow_threshold = 0.1    # Everything below 0.1 s^-1
fast_threshold = 150    # Everything above 150 s^-1

#number reactions does each enzyme catalyze
rxn_count = df.groupby('gene_id').size().rename('n_reactions')
enzyme_df = enzyme_df.merge(rxn_count, left_on='gene_id', right_index=True)

slow_enzymes = (enzyme_df[enzyme_df['kapp (1/s)'] < slow_threshold]
                [['rxnid', 'gene_id', 'kapp (1/s)', 'n_reactions']]
                .sort_values(by='kapp (1/s)', ascending=True))

fast_enzymes = (enzyme_df[enzyme_df['kapp (1/s)'] > fast_threshold]
                [['rxnid', 'gene_id', 'kapp (1/s)', 'n_reactions']]
                .sort_values(by='kapp (1/s)', ascending=False))

print(f"\n--- SLOW ENZYMES (< {slow_threshold} s^-1) ---")
print(slow_enzymes.head(15).to_string(index=False))

print(f"\n--- FAST ENZYMES (> {fast_threshold} s^-1) ---")
print(fast_enzymes.head(40).to_string(index=False))

with pd.ExcelWriter('Extreme_Kapps_Analysis_per_enzyme.xlsx') as writer:
    slow_enzymes.to_excel(writer, sheet_name='Slow Enzymes', index=False)
    fast_enzymes.to_excel(writer, sheet_name='Fast Enzymes', index=False)
