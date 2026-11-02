import uproot
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# List of ROOT files to process
#root_files = [
#    "unfold_results/unfold_results_argon_6bin_smooth_cut01.root",
#    "unfold_results/unfold_results_argon_6bin_smooth_cutstk.root",
#    "unfold_results/unfold_results_argon_6bin_smooth_cutpsd.root",
#    "unfold_results/unfold_results_argon_6bin_smooth_bkg_v2.root"
#]

root_files = [
    "/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/ARGON_ANALYSIS/unfold_results/unfold_results_argon_6bin_low1_up2_v2.root",
    "/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/SULFUR_ANALYSIS/unfold_results/unfold_results_flux_spectrum_sulfur_120months_6bin_smooth_2e5up_1e5low_bkg_v2.root",
    "/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/GROUP_ANALYSIS/unfold_results/unfold_results_flux_spectrum_SArCa_120months_6bin_smooth_mix_v2.root"
]

hist_name = "hacc1"

plt.figure(figsize=(8, 6))
colors = sns.color_palette("icefire", n_colors=3)  #viridis, icefire, mako

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
    #label = ["HET + fiducial volume", "STK cut", "PSD cut", "charge cut"][root_files.index(filename)]
    #label = ["ARGON", "SULFUR"][root_files.index(filename)]
    label = ["ARGON", "SULFUR", "GROUP"][root_files.index(filename)]

    # Draw as step histogram
    #plt.step(centers, values, where="mid", linewidth=1.8, label=label)
    plt.errorbar(x, y, xerr=xerr, fmt='.', label=label, markersize=5, capsize=0, color=colors[root_files.index(filename)])

# Plot formatting
plt.xlabel("Energy (GeV)")
plt.ylabel(r"Acceptance (m$^2$ sr)")

plt.xscale("log")
plt.ylim(0.001, 0.12)

#plt.title(f"Acceptance")
plt.legend(loc='lower center',
    frameon=False,
    prop={'size': 12}
    )
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

#plt.savefig("plots/acceptance_argon_oct26.png", dpi=300)
plt.show()