import numpy as np
import matplotlib.pyplot as plt

Z = 16. #16. 
A = 32.
m_p = 0.938
m_S = m_p *A
Ebr = 15300 * Z
err = 1600 * Z



def load_flux_data(filepath):
    # legge tutte le colonne insieme
    Emin,  Emax,  Emed,  flux_pow,  flux_pow_err = np.loadtxt(
        filepath,
        skiprows=2,
        unpack=True
    )
    return {
        "Emin": Emin,
        "Emax": Emax,
        "Emed": Emed,
        "flux_pow": flux_pow,
        "flux_pow_err": flux_pow_err
    }


S2 = load_flux_data('data/flux_pow_SULFUR_6bin_2e7_2sigma.txt')
S2e5 = load_flux_data('data/flux_pow_SULFUR_6bin_2e7_2e5sigma.txt')
S3 = load_flux_data('data/flux_pow_SULFUR_6bin_2e7_3sigma.txt')

S = load_flux_data('data/flux_pow_SULFUR_6-4bin_2e7.txt')

plt.figure(figsize=(8,5))
'''
data_S = "/mnt/c/Users/saraf/Desktop/dampe/range-16-20/z16-20_noreg_6bin.txt"
Emean_S = np.loadtxt(data_S, skiprows=1, usecols=(0,), unpack=True)
Flux_S = np.loadtxt(data_S,skiprows=1,usecols=(3,),unpack=True)
Flux_low_S = np.loadtxt(data_S,skiprows=1,usecols=(4,),unpack=True)
Flux_up_S = np.loadtxt(data_S,skiprows=1,usecols=(4,),unpack=True)

Flux_27_S = Flux_S * Emean_S**2.7
Flux_low_27_S = Flux_low_S * Emean_S**2.7
Flux_up_27_S = Flux_up_S * Emean_S**2.7
'''

'''
plt.errorbar(
    Emean_S, 
    Flux_27_S, 
    yerr=[Flux_low_27_S,Flux_up_27_S],
    fmt="d",
    color="grey",
    ecolor="grey",
    elinewidth=1,
    capsize=3,
    markersize=4,
    label="S-Ar-Ca no smoothing")
'''

def plot_error(dataset, label, color):

    plt.errorbar(
        dataset["Emed"][:-1], #[:-1]
        dataset["flux_pow"][:-1],
        yerr=dataset["flux_pow_err"][:-1],
        fmt="o",
        color=color,
        ecolor=color,
        elinewidth=1,
        capsize=3,
        markersize=6,
        label=label
    )

plot_error(S, "SULFUR 6-4bin 2e7", "#1900FFC3")

'''
plot_error(S2, "SULFUR 2sigma", "#FF6F00C3")
plot_error(S2e5, "SULFUR 2.5sigma", "#0062FFC7")
plot_error(S3, "SULFUR 3sigma", "#51FF00C7")
'''

data_ams = "/mnt/c/Users/saraf/Desktop/dampe/sulfur/AMS/sulfur_ams.txt"
#Quantity R-GV bin-range Flux-(GV m^2 s sr)^-1 stat-err syst-err
R = np.loadtxt(data_ams, skiprows=1, usecols=(1,), unpack=True)
Flux_ams_R = np.loadtxt(data_ams, skiprows=1, usecols=(4,), unpack=True)
stat_err = np.loadtxt(data_ams, skiprows=1, usecols=(5,), unpack=True)
sys_err = np.loadtxt(data_ams, skiprows=1, usecols=(6,), unpack=True)

E_ams = np.sqrt((Z * R)**2 + m_S**2) - m_S
dR_dE = np.sqrt((Z * R)**2 + m_S**2) / (Z**2 * R)
Flux_ams_E = Flux_ams_R * dR_dE
stat_err_E = stat_err * dR_dE
sys_err_E  = sys_err  * dR_dE

Flux_ams_E_27 = Flux_ams_E * E_ams**2.7
stat_err_E_27 = stat_err_E * E_ams**2.7
sys_err_E_27  = sys_err_E  * E_ams**2.7
err_tot_E_27 = np.sqrt(stat_err_E_27**2 + sys_err_E_27**2)


plt.errorbar(
    E_ams, Flux_ams_E_27,
    yerr=[err_tot_E_27, err_tot_E_27],
    fmt="d",
    color="grey",
    ecolor="grey",
    elinewidth=1,
    capsize=3,
    markersize=4,
    label="AMS (2023)"
)


#softening position
#plt.vlines(Ebr, 0, 2000, linestyles='--', linewidth= 2, alpha=0.5, colors="b", label="Energy Break Softening")
#plt.vlines(Ebr-err, 0, 1400, linestyles='--', colors="skyblue")
#plt.vlines(Ebr+err, 0, 1400, linestyles='--', colors="skyblue")
#plt.axvspan(Ebr-err,Ebr+err, color='skyblue', alpha=0.5)

plt.xlabel("Energy [GeV]", fontsize=14)
plt.ylabel(r"$E^{2.7}$ Flux [m$^{-2}$ s$^{-1}$ sr$^{-1}$ (GeV)$^{1.7}$]", fontsize=14)
plt.xscale("log")
#plt.yscale("log")
#plt.ylim(-10, 2000)
plt.xlim(1e1,1e6)

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
'''

plt.text(0.05, 0.8, 'Preliminary',
         transform=plt.gca().transAxes,
         fontsize=20,
         color='gray',
         verticalalignment='top')


plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.tight_layout()

plt.savefig("plots/flux_SULFUR_6-4bin_2e7.png", dpi=300)
plt.show()

