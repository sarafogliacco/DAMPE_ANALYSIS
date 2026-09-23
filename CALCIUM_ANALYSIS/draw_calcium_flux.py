import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import proton_mass


def load_spectrum(filename, gamma=2.7):
    E_med, flux, flux_err = np.loadtxt(
        filename,
        skiprows=1,
        unpack=True,
        usecols=(0, 1, 2)
    )

    return {
        "E_med": E_med,
        "flux": flux,
        "flux_err": flux_err,
        "flux_pow": flux * E_med**gamma *0.9,
        "flux_pow_err": flux_err * E_med**gamma
    }

def plot_error(dataset, label, color,fmt):

    plt.errorbar(
        dataset["E_med"], #[:-4]
        dataset["flux_pow"],
        yerr=dataset["flux_pow_err"],
        fmt=fmt,
        color=color,
        ecolor=color,
        elinewidth=1,
        capsize=3,
        markersize=6,
        label=label
    )

#####AMS
# rigidity_min GV,rigidity_max GV,S_flux m^-2sr^-1s^-1GV^-1,S_flux_error_statistical m^-2sr^-1s^-1GV^-1,S_flux_error_acceptance m^-2sr^-1s^-1GV^-1,S_flux_error_unfolding m^-2sr^-1s^-1GV^-1,S_flux_error_rigidity_scale m^-2sr^-1s^-1GV^-1,S_flux_error_systematic_total m^-2sr^-1s^-1GV^-1

def load_ams_spectrum(filename, Z, A, gamma=2.7):
    R_min, R_max, Flux, Flux_err_stat, Flux_err_acc, Flux_err_unf, Flux_err_rig, Flux_err_sys = np.loadtxt(
        filename,
        skiprows=1,
        unpack=True,
        delimiter=",",
        usecols=(0, 1, 2, 3, 4, 5, 6, 7)
    )

    R_med = (R_min + R_max) / 2
    Flux_err_r = np.sqrt(Flux_err_stat**2 + Flux_err_acc**2 + Flux_err_unf**2 + Flux_err_rig**2 + Flux_err_sys**2)
    m=proton_mass*A
    E_ams = np.sqrt((Z * R_med)**2 + m**2) - m
    dR_dE = np.sqrt((Z * R_med)**2 + m**2) / (Z**2 * R_med)
    Flux_ams_E = Flux * dR_dE
    stat_err_E = Flux_err_r * dR_dE

    return {
        "E_med": E_ams,
        "flux_pow": Flux_ams_E * E_ams**gamma,
        "flux_pow_err": stat_err_E * E_ams**gamma
    }
Ca_ams = load_ams_spectrum("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/GROUP_ANALYSIS/data/ams/ca_ams.csv", Z=20, A=40)

flux_ams =  Ca_ams["flux_pow"]
err_ams = Ca_ams["flux_pow_err"]
E_ams = Ca_ams["E_med"]

plt.figure(figsize=(10,8))

#Ca_4bin = load_spectrum("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/CALCIUM_ANALYSIS/data/flux_spectrum_calcium_4bin_smooth.dat")
#Ca_5bin = load_spectrum("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/CALCIUM_ANALYSIS/data/flux_spectrum_calcium_5bin_smooth.dat")
Ca_6bin = load_spectrum("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/CALCIUM_ANALYSIS/data/flux_spectrum_calcium_6bin_smooth.dat")
#Ca_7bin = load_spectrum("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/CALCIUM_ANALYSIS/data/flux_spectrum_calcium_7bin_smooth.dat")
#plot_error(DAMPE_Ca, label="DAMPE", color="blue", fmt="o")

#plot_error(Ca_4bin, label="DAMPE 4 bins", color="blue", fmt="o")
#plot_error(Ca_5bin, label="DAMPE 5 bins", color="orange", fmt="s")
plot_error(Ca_6bin, label="DAMPE 6 bins", color="green", fmt="^")
#plot_error(Ca_7bin, label="DAMPE 7 bins", color="red", fmt="d")


plt.errorbar( E_ams, flux_ams, yerr=err_ams, fmt="o", color="grey", ecolor="grey", elinewidth=1, capsize=3, markersize=4, label="AMS")


plt.xlabel("Energy [GeV]", fontsize=14)
plt.ylabel(r"$E^{2.7}$ Flux [m$^{-2}$ s$^{-1}$ sr$^{-1}$ (GeV)$^{1.7}$]", fontsize=14)
plt.xscale("log")
#plt.xlim(10, 1e6)

plt.legend()

plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.show()