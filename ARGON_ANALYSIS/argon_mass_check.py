import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import proton_mass

Z = 18. 
m_p = proton_mass

gamma = 2.7


def load_flux_data(filepath):
    # legge tutte le colonne insieme
    Emed,  flux,  flux_err = np.loadtxt(
        filepath,
        skiprows=1,
        unpack=True,
        usecols=(0, 1, 2)
    )
    E_nucleon = Emed / A
    flux_pow = flux * Emed**2.7  #* A**(-1.7)
    flux_pow_err = flux_err * Emed**2.7 #* A**(-1.7)
    return {
        "Emed": Emed,
        #"E_nucleon": E_nucleon,
        "flux_pow": flux_pow,
        "flux_pow_err": flux_pow_err
    }


plt.figure(figsize=(10,8))
colors = plt.cm.cividis(np.linspace(0.2, 0.8, 6))


def plot_error(dataset, label, color):

    plt.errorbar(
        dataset["Emed"][:-2], #[:-1]
        dataset["flux_pow"][:-2], #[:-1]
        yerr=dataset["flux_pow_err"][:-2], #[:-1]
        fmt="o",
        color=color,
        ecolor=color,
        elinewidth=1,
        capsize=3,
        markersize=6,
        label=label
    )



###### AMS
##### AMS published data 2026
AMS_26 = "/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/ARGON_ANALYSIS/data/other_exp/Ar_AMS_2026.csv"
R_min, R_max, flux_R, stat_err_R, acc_err_R, unf_err_R, rigidity_err_R, syst_err_R = np.loadtxt(AMS_26, skiprows=1, unpack=True, delimiter=",")

def convert_rigidity_to_energy(R_min, R_max, flux_R, acc_err_R, unf_err_R, rigidity_err_R, syst_err_R, A, Z, gamma, m_p):

    m_Ar = m_p * A
    R_med = np.sqrt(R_min * R_max)
    E_AMS = np.sqrt((Z * R_med)**2 + m_Ar**2) - m_Ar
    dR_dE = np.sqrt((Z * R_med)**2 + m_Ar**2) / (Z**2 * R_med)
    flux_AMS_E = flux_R * dR_dE
    tot_err_R = np.sqrt( acc_err_R**2 + unf_err_R**2 + rigidity_err_R**2 + syst_err_R**2)
    tot_err_E = tot_err_R * dR_dE
    flux_AMS_pow = flux_AMS_E * E_AMS**gamma
    tot_err_pow = tot_err_E * E_AMS**gamma

    return E_AMS, flux_AMS_pow, tot_err_pow

E_AMS, flux_AMS_pow, tot_err_pow = convert_rigidity_to_energy(R_min, R_max, flux_R, acc_err_R, unf_err_R, rigidity_err_R, syst_err_R, A=36., Z=18., gamma=2.7, m_p=m_p)
E_AMS_2, flux_AMS_pow_2, tot_err_pow_2 = convert_rigidity_to_energy(R_min, R_max, flux_R, acc_err_R, unf_err_R, rigidity_err_R, syst_err_R, A=39.948, Z=18., gamma=2.7, m_p=m_p)

plt.errorbar(E_AMS, flux_AMS_pow, yerr=tot_err_pow, markerfacecolor='none', markeredgecolor="r", ecolor="r", fmt='s', markersize=4, label="A = 36 ")
plt.errorbar(E_AMS_2, flux_AMS_pow_2, yerr=tot_err_pow_2  , markerfacecolor='none', markeredgecolor="b", ecolor="b", fmt='s', markersize=4, label="A = 40")


plt.xlabel("Energy [GeV]", fontsize=14)
plt.ylabel(r"$E^{2.7}$ Flux [m$^{-2}$ s$^{-1}$ sr$^{-1}$ (GeV)$^{1.7}$]", fontsize=14)
plt.xscale("log")
#plt.yscale("log")
#plt.ylim(-10, 1500)
#plt.xlim(3*1e1,1e6)

plt.legend(
    loc='center left',
    frameon=False,
    prop={'size': 14}
    )

plt.text(0.05, 0.9, 'ARGON', #'Sulfur Z=16'
         transform=plt.gca().transAxes,
         fontsize=30,
         color='k',
         verticalalignment='top')


plt.text(0.5, 0.9, 'Preliminary',
         transform=plt.gca().transAxes,
         fontsize=30,
         color='gray',
         verticalalignment='top')


plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

#plt.savefig("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/ARGON_ANALYSIS/plots/flux_ARGON_oct26.png", dpi=300)
#plt.show()


print("m_p:", m_p)
print("Massa A=36:", m_p * 36.)
print("Massa A=39.948:", m_p * 39.948)

print("Prime energie A=36:", E_AMS[:5])
print("Prime energie A=40:", E_AMS_2[:5])

print("Differenza assoluta energia:",
      np.max(np.abs(E_AMS_2 - E_AMS)))

print("Differenza relativa energia [%]:",
      np.max(np.abs((E_AMS_2 - E_AMS) / E_AMS)) * 100)

