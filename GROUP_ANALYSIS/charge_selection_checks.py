####check cut on charge selcetion 
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import proton_mass

gamma = 2.7 #spectral index
Z = 16. 
A = 32.06
m_p = proton_mass
m_S = m_p *A

def load_spectrum(filename, gamma=2.7):
    E_med, flux, flux_err = np.loadtxt(
        filename,
        skiprows=1,
        unpack=True,
        usecols=(0, 3, 4)
    )

    return {
        "E_med": E_med,
        "flux": flux,
        "flux_err": flux_err,
        "flux_pow": flux * E_med**gamma,
        "flux_pow_err": flux_err * E_med**gamma
    }


S = load_spectrum('../SULFUR_ANALYSIS/data/flux_spectrum_sulfur_120months_6-4bin_smooth_v4.dat')
Ar = load_spectrum("../ARGON_ANALYSYS/data/flux_spectrum_argon_6-4bin_smooth_v2.dat")
sarca = load_spectrum("data/flux_spectrum_SArCa_120months_6-4bin_smooth.dat")


def plot_error(dataset, label, color,fmt):

    plt.errorbar(
        dataset["E_med"][:-1], #[:-1]
        dataset["flux_pow"][:-1],
        yerr=dataset["flux_pow_err"][:-1],
        fmt=fmt,
        color=color,
        ecolor=color,
        elinewidth=1,
        capsize=3,
        markersize=6,
        label=label
    )





plt.figure(figsize=(8,5))

plot_error(S, "SULFUR DAMPE", "blue", "x")
plot_error(Ar, "ARGON DAMPE", "#FF6A00", "+")
plot_error(sarca, "SArCa DAMPE", "#00C130", "*")
#plot_error(S_2up, "S +2sigma", "green", "x")
#plot_error(S_1, "S +1sigma", "red", "o")


plt.xlabel("Energy [GeV]", fontsize=14)
plt.ylabel(r"$E^{2.7}$ Flux [m$^{-2}$ s$^{-1}$ sr$^{-1}$ (GeV)$^{1.7}$]", fontsize=14)
plt.xscale("log")
#plt.yscale("log")
#plt.ylim(-10, 2000)
plt.xlim(1e2,1e6)

plt.legend(
    loc='upper left',
    frameon=False,
    prop={'size': 14}
    )
'''
plt.text(0.05, 0.95, 'Z=16-20', #'Sulfur Z=16'
         transform=plt.gca().transAxes,
         fontsize=20,
         color='k',
         verticalalignment='top')


plt.text(0.05, 0.8, 'Preliminary',
         transform=plt.gca().transAxes,
         fontsize=20,
         color='gray',
         verticalalignment='top')
'''

#plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.tight_layout()

#plt.savefig("plots/flux_SULFUR_6-4bin_2e7.png", dpi=300)
plt.show()