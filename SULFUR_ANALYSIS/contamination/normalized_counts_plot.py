import uproot
import numpy as np
import matplotlib.pyplot as plt


def load_and_normalize(filename, histname, weights):

    with uproot.open(filename) as f:
        h = f[histname]

        counts = h.values().astype(float)
        variances = h.variances()
        edges = h.axes[0].edges()

    errors = np.sqrt(variances)
    centers = 0.5 * (edges[:-1] + edges[1:])
    
    for i, x in enumerate(centers):

        for (xmin, xmax), w in weights.items():

            if xmin <= x < xmax:
                counts[i] *= w
                errors[i] *= w
                break
    
    return edges, centers, counts, errors 

weights_Si = {
    (1e2, 1e3): 4.0 / 685625,
    (1e3, 1e4): 4.0 / 586530,
    (1e4, 1e5): 4.0 / 108468,
    (1e5, 5e6): 4.0 / 34810,
}

weights_S = {
    (1e2, 1e3): 1.0 / 232080,
    (1e3, 1e4): 1.0 / 592975,
    (1e4, 1e5): 1.0 / 122556,
    (1e5, 5e6): 1.0 / 18610,
}

# file 1
edges1, centers1, counts1, errors1 = load_and_normalize(
    "Out_MC_FTFP_S_E2e7_6_bin_500TeV_v2.root",
    "h1Ngen",
    weights_S
    #"Out_DATA_120months_SULFUR_6bin_v2.root",
    #"h1Nobs",
)

# file 2
edges2, centers2, counts2, errors2 = load_and_normalize(
    "Out_MC_FTFP_Si_E2e7_6bin_500TeV.root",
    "h1Ngen",
    weights_Si
    #" Out_DATA_120months_SILICON_6bin.root",
    #"h1Nobs",
)

# plot
fig, ax = plt.subplots(figsize=(8,6))

ax.stairs(counts1, edges1, label="SULFUR")
#ax.errorbar(centers1, counts1, yerr=errors1, fmt="o")

ax.stairs(counts2, edges2, label="SILICON")
#ax.errorbar(centers2, counts2, yerr=errors2, fmt="s")

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel("Energy [GeV]")
ax.set_ylabel("Weighted events MC")
#ax.set_ylabel("Entries DATA")
ax.legend()

plt.tight_layout()
plt.savefig("comparisonMC4.png", dpi=300)
plt.show()