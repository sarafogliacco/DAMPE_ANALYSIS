import uproot 
import awkward as ak
import numpy as np
from hist import Hist
from hist.axis import Variable
import matplotlib.pyplot as plt

# =========================
# Silicon
# =========================

file_Si = "Out_DATA_120months_SILICON_6bin.root"
weights_Si = np.loadtxt(
    "weights_data_Si.txt",
    usecols=(0,)
)

h_Si = uproot.open(file_Si)["h1Nobs"]

values_Si = h_Si.values()
edges_Si = h_Si.axes[0].edges()

weights_Si_safe = np.where(weights_Si > 0, weights_Si, np.nan)

values_norm_Si = 4.0 / weights_Si_safe

# =========================
# Sulfur
# =========================

file_S = "Out_DATA_120months_SULFUR_6bin_v2.root"
weights_S = np.loadtxt(
    "weights_data_S.txt",
    usecols=(0,)
)

h_S = uproot.open(file_S)["h1Nobs"]

values_S = h_S.values()
edges_S = h_S.axes[0].edges()

weights_S_safe  = np.where(weights_S  > 0, weights_S,  np.nan)

values_norm_S = 1 / weights_S_safe

# =========================
# Plot
# =========================

centers_Si = 0.5 * (edges_Si[1:] + edges_Si[:-1])
centers_S = 0.5 * (edges_S[1:] + edges_S[:-1])

plt.figure(figsize=(8,6))

plt.step(
    centers_Si,
    values_norm_Si,
    where="mid",
    label="Silicon"
)

plt.step(
    centers_S,
    values_norm_S,
    where="mid",
    label="Sulfur"
)

plt.xscale("log")
plt.yscale("log")

plt.xlabel("Et")
plt.ylabel("Normalized entries")

plt.legend()
plt.tight_layout()

#plt.savefig("Si_vs_Sulfur.png", dpi=300)
print("\nBin-by-bin comparison:")
print(f"{'Et center':>12} {'Si':>15} {'S':>15} {'Diff [%]':>15}")

for x, si, s in zip(centers_Si, values_norm_Si, values_norm_S):

    if (not np.isfinite(si)) or (not np.isfinite(s)):
        continue

    if si <= 0 or s <= 0:
        continue

    diff = 100.0 * (s - si) / ((s + si) / 2.0)

    print(
        f"{x:12.2e} "
        f"{si:15.6e} "
        f"{s:15.6e} "
        f"{diff:15.2f}"
    )

#plt.savefig("Si_vs_Sulfur_4.png", dpi=300)
plt.show()