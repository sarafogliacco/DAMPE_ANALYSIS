####check cut on charge selcetion 
import numpy as np
import matplotlib.pyplot as plt


gamma = 2.7 #spectral index
Z = (16.+18.+20.)/3. 
Ebr = 15400. * Z
err = 4800. * Z

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


S = load_spectrum('/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS//SULFUR_ANALYSIS/data/flux_spectrum_sulfur_120months_6-4bin_smooth_v4.dat')
Ar = load_spectrum("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS//ARGON_ANALYSYS/data/flux_spectrum_argon_6-4bin_smooth_v2.dat")

sarca_a = load_spectrum("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/GROUP_ANALYSIS/data/flux_spectrum_SArCa_120months_6-4bin_smooth.dat")
sarca_b = load_spectrum("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/GROUP_ANALYSIS/data/flux_spectrum_SArCa_120months_6-4bin_smooth_v2.dat")

sarca_6bin = load_spectrum("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/GROUP_ANALYSIS/data/flux_spectrum_SArCa_120months_6bin_smooth.dat")

def plot_error(dataset, label, color,fmt):

    plt.errorbar(
        dataset["E_med"][:-2], #
        dataset["flux_pow"][:-2],
        yerr=dataset["flux_pow_err"][:-2],
        fmt=fmt,
        color=color,
        ecolor=color,
        elinewidth=1,
        capsize=3,
        markersize=6,
        label=label
    )

################## old 
def load_flux_data(filepath):
    # Emin  Emax  Emed  flux_pow  flux_pow_err
    Emed, Emin, Emax, flux_pow, flux_pow_err = np.loadtxt(
        filepath,
        skiprows=2,
        unpack=True
    )
    return {
        "Emin": Emin,
        "Emax": Emax,
        "Emean": Emed,
        "Flux": flux_pow,
        "Flux_stat": flux_pow_err
    }

data_sarca = "/mnt/c/Users/saraf/Desktop/dampe/unfolding/data/z16_20/flux_pow_SArCa_6bin.txt"
E_min, E_max, E_old, flux_old, flux_err_old = np.loadtxt(data_sarca, skiprows=1, usecols=(0,1,2,3,4), unpack=True, delimiter=" ")

#DAMPE data plotting


plt.figure(figsize=(10,8))

#plot_error(sarca_a, "SArCa from Ar", "orange", "v")
#plot_error(sarca_b, "SArCa from S", "blue", "^")

#plot_error(sarca_b, "6-4 bins", "orange", "v")
#plot_error(sarca_6bin, "6bins", "red", "^")
plot_error(sarca_b, "DAMPE (this work)", "green", "o")

#plt.errorbar(E_old, flux_old, yerr=flux_err_old, color="g", fmt=".", label= "SArCa from Si")




#softening position
plt.vlines(Ebr, 0, 6000, linestyles='--', linewidth= 2, alpha=0.5, colors="k", label="Energy Break Softening")
plt.axvspan(Ebr-err,Ebr+err, color='skyblue', alpha=0.5)

#plt.axvspan(3.5*1e5,1e6, color='gray', alpha=0.3)

plt.xlabel("Energy [GeV]", fontsize=14)
plt.ylabel(r"$E^{2.7}$ Flux [m$^{-2}$ s$^{-1}$ sr$^{-1}$ (GeV)$^{1.7}$]", fontsize=14)
plt.xscale("log")
#plt.yscale("log")
plt.ylim(5e2, 3000)
plt.xlim(1e2,1e6)

plt.legend(
    loc='center left',
    frameon=False,
    prop={'size': 14}
    )

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


plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()

plt.savefig("plots/flux_sarca_Ebr_240626.png", dpi=300)
plt.show()