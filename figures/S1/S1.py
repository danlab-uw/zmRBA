#Ribosome abundance vs growth rate (mu) 
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtik

#To get these values vglc was set to 30, 35, 40, 45, 50, 55.1, 60, and 65 mmol/gDCW/hr and zmRBA binary search run 
rib_prots = [1.42e-3, 1.68e-3, 1.95e-3, 2.22e-3, 2.49e-3, 2.77e-3, 3.02e-3, 3.31e-3 ]
rib_cap = [19.606790021671716, 23.334455894967474, 27.066228018870786, 30.79971806857517, 34.53681495322886,38.352491977762976, 42.02120201249414, 45.767485099259975]
mu = [0.18500, 0.21995, 0.25488, 0.28975, 0.32458, 0.36008, 0.39413, 0.42885]

plt.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":["Arial","Helvetica","Liberation Sans","DejaVu Sans"],
    "font.size":13,"axes.linewidth":1.1,
    "xtick.direction":"out","ytick.direction":"out"})

BLUE="#2E6CA4"

fig, ax = plt.subplots(figsize=(6.8,4.7))

ax.plot(mu, rib_prots, "-o", color=BLUE, ms=6, lw=1.8,
        markeredgecolor="white", markeredgewidth=0.6, zorder=3)


ax.set_xlabel("Growth rate (h$^{-1}$)")
ax.set_ylabel("Ribosomal protein abundance \n (mmol gDCW$^{-1}$)")
ax.set_xlim(0.15, 0.45)
ax.set_ylim(0, 3.6e-3)
ax.spines[["top","right"]].set_visible(False)
ax.tick_params(which="both", length=4, width=1.1)
ax.yaxis.set_major_formatter(mtik.FuncFormatter(
    lambda v, pos: "0" if v == 0 else rf"${v/1e-3:.1f}{{\!\times\!}}10^{{-3}}$"))
fig.tight_layout()
for ext in ("pdf","png","svg", "eps"):
    fig.savefig(f"FigureS1.{ext}", dpi=300, bbox_inches="tight")
print("saved")

plt.show()
