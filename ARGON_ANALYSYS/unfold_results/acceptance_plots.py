import uproot
import numpy as np
import matplotlib.pyplot as plt

# List of ROOT files to process
root_files = [
    "unfold_results_argon_6bin_E27_het.root",
    "unfold_results_argon_6bin_E27_stk.root",
    "unfold_results_argon_6bin_E27_psd.root",
    "unfold_results_argon_6bin_E27_charge.root"
]

hist_name = "hacc1"

plt.figure(figsize=(8, 6))

for filename in root_files:
    # Open file and read histogram
    with uproot.open(filename) as f:
        h = f[hist_name]

        # Get bin contents and edges
        values = h.values()
        edges = h.axis().edges()
        variances = h.variances()

    # Compute bin centers
    x = 0.5 * (edges[:-1] + edges[1:])
    xerr = 0.5 * np.diff(edges)

    y = values
    yerr = np.sqrt(variances) if variances is not None else np.sqrt(np.clip(y, 0, None))

    # Use filename (without path) as legend label
    label = ["HET", "STK cut", "PSD cut", "charge cut"][root_files.index(filename)]

    # Draw as step histogram
    #plt.step(centers, values, where="mid", linewidth=1.8, label=label)
    plt.errorbar(x, y, xerr=xerr, fmt='.', label=label, markersize=5, capsize=0)

# Plot formatting
plt.xlabel("Energy (GeV)")
plt.ylabel(r"Acceptance (m$^2$ sr)")

plt.xscale("log")
plt.ylim(0.0, 0.25)

#plt.title(f"Acceptance")
plt.legend(loc='lower center',
    frameon=False,
    prop={'size': 12}
    )
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

plt.savefig("../plots/acceptance_argon_6bin_2e7.png", dpi=300)
plt.show()