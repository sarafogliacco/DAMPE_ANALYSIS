import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import proton_mass

Z = 18. 
A = 39.948
m_p = proton_mass
m_Ar = m_p *A
Ebr = 15400 * Z
err = 4800 * Z
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

Ar = load_flux_data('data/flux_spectrum_argon_6bin_smooth.dat')
#Ar_b = load_flux_data('data/flux_spectrum_argon_5-4bin_smooth_v2.dat')
#Ar_c = load_flux_data('data/flux_spectrum_argon_3bin_smooth.dat')
Ar_old = "data/v0/flux_spectrum_argon_6-4bin_smooth_v2.dat"
E_old, flux_old, flux_err_old = np.loadtxt(Ar_old, skiprows=1, usecols=(0,3,4), unpack=True)

'''
decade_scale = {
    2: -(0.0001 + 0.007),
    3: -(0.0004 + 0.006),
    4: -(0.08 + 0.01),
    5: -(0.2 + 0.11)
}

logE = np.floor(np.log10(Ar["Emed"])).astype(int)

scale_factors = np.array([
    1.0 + decade_scale.get(dec, 0.0)
    for dec in logE
])
'''
Ar["flux_scaled"] = Ar["flux_pow"] * 0.9
Ar["flux_err_scaled"] = Ar["flux_pow_err"] * 0.9



#S = load_flux_data('../SULFUR_ANALYSIS/data/flux_spectrum_sulfur_120months_6-4bin_smooth_v5.dat')

plt.figure(figsize=(10,8))
colors = plt.cm.cividis(np.linspace(0.2, 0.8, 6))


def plot_error(dataset, label, color):

    plt.errorbar(
        dataset["Emed"][:-1], #[:-1]
        dataset["flux_pow"][:-1], #[:-1]
        yerr=dataset["flux_pow_err"][:-1], #[:-1]
        fmt="o",
        color=color,
        ecolor=color,
        elinewidth=1,
        capsize=3,
        markersize=6,
        label=label
    )

plot_error(Ar, "DAMPE (ML_new)", "r")
#plot_error(Ar_b, "ARGON (v2)", "#FF0000")
#plot_error(Ar_c, "ARGON (3bin)", "#006800")
plt.errorbar(E_old, flux_old*E_old**2.7, yerr=flux_err_old*E_old**2.7, color="b", label="DAMPE (ML_old)", fmt="x")
#plt.errorbar(Ar["Emed"],Ar["flux_pow"],yerr=Ar["flux_pow_err"],color="g",fmt="+", label="DAMPE ARGON")

#plot_error(S, "SULFUR", "#FF0000")

###### AMS
##### AMS published data 2026
AMS_26 = "data/other_exp/Ar_AMS_2026.csv"
R_min, R_max, flux_R, stat_err_R, acc_err_R, unf_err_R, rigidity_err_R, syst_err_R = np.loadtxt(AMS_26, skiprows=1, unpack=True, delimiter=",")

R_med = np.sqrt(R_min * R_max)

E_AMS = np.sqrt((Z * R_med)**2 + m_Ar**2) - m_Ar
dR_dE = np.sqrt((Z * R_med)**2 + m_Ar**2) / (Z**2 * R_med)

flux_AMS_E = flux_R * dR_dE

tot_err_R = np.sqrt(acc_err_R**2 + unf_err_R**2 + rigidity_err_R**2 + syst_err_R**2)
tot_err_E = tot_err_R * dR_dE


flux_AMS_pow = flux_AMS_E * E_AMS**gamma
tot_err_pow = tot_err_E * E_AMS**gamma 


#plt.errorbar(E_ams, flux_ams, yerr=[flux_ams_err_d, flux_ams_err_u], fmt=".", color="grey", ecolor="grey", elinewidth=1, capsize=3, markersize=6, label="AMS (ICRC 2025)")
#plt.scatter(E_ams, flux_ams, color="grey", label="AMS (ICRC 2025)", marker = ".")
plt.errorbar(E_AMS, flux_AMS_pow, yerr=tot_err_pow, markerfacecolor='none', markeredgecolor=colors[2], ecolor=colors[2], fmt='s', markersize=4, label="AMS (2026)")

#########TRACER

##### TRACER
tracer_s = "/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/ARGON_ANALYSYS/data/other_exp/Ar_TRACER_2008.txt"
E_tracer_n, flux_tracer, err_up_tracer, err_dw_tracer = np.loadtxt(tracer_s, skiprows=1, usecols=(3,5,6,7), delimiter=",", unpack=True)

E_tracer = E_tracer_n * A
flux_pow_tracer =   flux_tracer    * E_tracer_n**gamma * A**(gamma-1)
err_up_pow_tracer = err_up_tracer  * E_tracer_n**gamma * A**(gamma-1)
err_dw_pow_tracer = err_dw_tracer  * E_tracer_n**gamma * A**(gamma-1)

plt.errorbar(E_tracer, flux_pow_tracer, yerr=[err_dw_pow_tracer,err_up_pow_tracer], markerfacecolor='none', markeredgecolor=colors[3], ecolor=colors[3], fmt='*', markersize=6, label="TRACER (2008)")

tracer_s_99 = "/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/ARGON_ANALYSYS/data/other_exp/Ar_TRACER_1999.txt"
E_tracer_n_99, flux_tracer_99, err_tracer_99 = np.loadtxt(tracer_s_99, skiprows=1, usecols=(1,3,4), unpack=True)

E_tracer_99 = E_tracer_n_99 * A
flux_pow_tracer_99 =   flux_tracer_99    * E_tracer_n_99**gamma * A**(gamma-1)
err_pow_tracer_99 = err_tracer_99  * E_tracer_n_99**gamma * A**(gamma-1)

#plt.errorbar(E_tracer_99, flux_pow_tracer_99, yerr=err_pow_tracer_99, markerfacecolor='none', markeredgecolor=colors[4], ecolor=colors[4], fmt='+', markersize=6, label="TRACER (1999)")



#softening position
#plt.vlines(Ebr, -10, 2500, linestyles='--', linewidth= 1, alpha=0.5, colors="b", label="Energy Break Softening")
#plt.vlines(Ebr-err, -10, 2000, linestyles='--', colors="skyblue")
#plt.vlines(Ebr+err, -10, 2000, linestyles='--', colors="skyblue")
plt.axvspan(Ebr-err,Ebr+err, color='skyblue', alpha=0.5)

#plt.axvspan(2*1e5,1e6, color='gray', alpha=0.3)

plt.xlabel("Energy [GeV]", fontsize=14)
plt.ylabel(r"$E^{2.7}$ Flux [m$^{-2}$ s$^{-1}$ sr$^{-1}$ (GeV)$^{1.7}$]", fontsize=14)
plt.xscale("log")
#plt.yscale("log")
#plt.ylim(-10, 1500)
plt.xlim(1e1,1e6)

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

#plt.savefig("plots/flux_ARGON_band_240626.png", dpi=300)
plt.show()

