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


S = load_spectrum('data/flux_spectrum_sulfur_120months_6-4bin_smooth_v4.dat')
S_2e5up = load_spectrum("data/flux_spectrum_sulfur_120months_6-4bin_+2e5sigma_smooth.dat")
S_2up = load_spectrum("data/flux_spectrum_sulfur_120months_6-4bin_+2sigma_smooth.dat")
S_1 = load_spectrum("data/flux_spectrum_sulfur_120months_6-4bin_1sigma_smooth.dat")

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


#AMS 
#data_ams = "/mnt/c/Users/saraf/Desktop/dampe/sulfur/AMS/sulfur_ams.txt"
data_ams = "data/AMS_data_2023.txt"
#Quantity R-GV bin-range Flux-(GV m^2 s sr)^-1 stat-err syst-err
R_low = np.loadtxt(data_ams, skiprows=1, usecols=(0,), unpack=True)
R_high = np.loadtxt(data_ams, skiprows=1, usecols=(1,), unpack=True)
R = (R_low + R_high) / 2.

Flux_R = np.loadtxt(data_ams, skiprows=1, usecols=(2,), unpack=True)

stat_err = np.loadtxt(data_ams, skiprows=1, usecols=(3,), unpack=True)
acc_err = np.loadtxt(data_ams, skiprows=1, usecols=(4,), unpack=True)
unf_err = np.loadtxt(data_ams, skiprows=1, usecols=(5,), unpack=True)
scale_err = np.loadtxt(data_ams, skiprows=1, usecols=(6,), unpack=True)
syst_err = np.loadtxt(data_ams, skiprows=1, usecols=(7,), unpack=True)
tot_err = np.sqrt(stat_err**2 + acc_err**2 + unf_err**2 + scale_err**2 + syst_err**2)


E_ams = np.sqrt((Z * R)**2 + m_S**2) - m_S
dR_dE = np.sqrt((Z * R)**2 + m_S**2) / (Z**2 * R)
Flux_AMS = Flux_R * dR_dE
tot_err_E = tot_err * dR_dE

Flux_AMS_E_27 = Flux_AMS * E_ams**2.7
tot_err_E_27 = tot_err_E * E_ams**2.7



plt.figure(figsize=(8,5))

plot_error(S, "SULFUR DAMPE", "blue", "v")
#plot_error(S_2e5up, "S +2e5sigma", "orange", "+")
#plot_error(S_2up, "S +2sigma", "green", "x")
#plot_error(S_1, "S +1sigma", "red", "o")

plt.errorbar(
    E_ams, Flux_AMS_E_27,
    yerr=[tot_err_E_27, tot_err_E_27],
    fmt="d",
    color="grey",
    ecolor="grey",
    elinewidth=1,
    capsize=3,
    markersize=4,
    label="AMS (2023)"
)

plt.xlabel("Energy [GeV]", fontsize=14)
plt.ylabel(r"$E^{2.7}$ Flux [m$^{-2}$ s$^{-1}$ sr$^{-1}$ (GeV)$^{1.7}$]", fontsize=14)
plt.xscale("log")
#plt.yscale("log")
plt.ylim(-10, 2000)
#plt.xlim(1e1,1e6)

plt.legend(
    loc='center left',
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