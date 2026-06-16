import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt


files = [ 
    "unfold_results_argon_3bin_smooth.root",
    #"unfold_results_argon_4bin_smooth.root ",
    #"unfold_results_argon_5-4bin_smooth.root",
    #"unfold_results_argon_6-4bin_smooth_v3.root",
    #"unfold_results_argon_6bin_E27_smooth.root",
    #"unfold_results_argon_8bin_smooth.root"
    ]

hist_name = "hacc1"

scales = [1, 1 ,1, 1] # [1, 1, 1, 1, 1]

plt.figure(figsize=(8, 6))

for i, filename in enumerate(files):
    with uproot.open(filename) as f:
        h = f[hist_name]

        values = h.values()
        edges = h.axis().edges()
        variances = h.variances()

    x = 0.5 * (edges[:-1] + edges[1:])
    xerr = 0.5 * np.diff(edges)

    scale = scales[i]

    y = values * scale
    #y = values
    #yerr = np.sqrt(variances) if variances is not None else np.sqrt(np.clip(y, 0, None))

    # Use filename (without path) as legend label
    label = ["3 bins", "4 bins", "6 bins", "8 bins"][i] # "5-4 bins", "6-4 bins",

    # Draw as step histogram
    #plt.step(centers, values, where="mid", linewidth=1.8, label=label)
    plt.errorbar(x, y, xerr=xerr, fmt='.', label=f"{label} x {scale}", markersize=5, capsize=0)

# Plot formatting
plt.xlabel("Energy (GeV)")
plt.ylabel(r"Acceptance (m$^2$ sr)")

plt.xscale("log")
#plt.ylim(0.01, 0.11)

#plt.title(f"Acceptance")
plt.legend(loc='lower center',
    frameon=False,
    prop={'size': 12}
    )
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

#plt.savefig("../plots/acceptance_argon_6bin_2e7.png", dpi=300)
plt.show()