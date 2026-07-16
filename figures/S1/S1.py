#Isoprene molar yield vs growth rate (mu) at fixed glucose (55.1 mmol/gDCW/h)
import pandas as pd
import matplotlib.pyplot as plt

d = pd.read_csv("mu_sweep_results.csv").sort_values("mu")

plt.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Arial","Helvetica","Liberation Sans","DejaVu Sans"],
    "font.size":13,"axes.linewidth":1.1,
    "xtick.direction":"out","ytick.direction":"out"})

BLUE="#2E6CA4"

fig, ax = plt.subplots(figsize=(6.6,4.6))

ax.plot(d.mu, d.molar_yield_pct, "-o", color=BLUE, ms=6, lw=1.8,
        markeredgecolor="white", markeredgewidth=0.6, zorder=3)
ax.axvline(0.359, ls="--", color="0.55", lw=1.1, zorder=1)

#Case 3
ax.annotate("$\\mu$ = 0.359 h$^{-1}$\nmolar yield: 0.13%",
            xy=(0.3585,0.131), xytext=(0.335,0.75), fontsize=10, color="0.35",
            ha="left", arrowprops=dict(arrowstyle="->", color="0.55", lw=1.5))

#5% decrease mu 15-fold increase
ax.annotate("$\\mu$ = 0.34 h$^{-1}$\nmolar yield: 1.95%",
            xy=(0.3405,1.948), xytext=(0.315,2.62), fontsize=10, color=BLUE,
            ha="left", arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.5))

ax.set_xlabel("Growth rate (h$^{-1}$)")
ax.set_ylabel("Isoprene molar yield (%)")
ax.set_xlim(0.372, 0.03); ax.set_ylim(0, 3.05)
ax.spines[["top","right"]].set_visible(False)
ax.tick_params(which="both", length=4, width=1.1)
fig.tight_layout()
for ext in ("pdf","png","svg", "eps"):
    fig.savefig(f"FigureS1.{ext}", dpi=300, bbox_inches="tight")
print("saved")

plt.show()
