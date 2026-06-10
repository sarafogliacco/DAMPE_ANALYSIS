import numpy as np
import matplotlib.pyplot as plt

langaus_data = np.loadtxt('results_langaus_fit_MC_SULFUR_CORRECTED.txt')
mpv = langaus_data[:, 3]
mpv_err = langaus_data[:, 4]
energy = [133.3521653, 237.1374976, 421.6965375, 749.8939925, 1333.521653, 2371.374976, 4216.965375, 7498.939925, 17782.80068, 56234.15332]

plt.errorbar(energy, mpv, yerr=mpv_err, fmt='o', label='MPV with error bars')

plt.xscale('log')
plt.ylim(15, 18)
#plt.show()


print("MPV:", mpv)
print("MPV Error:", mpv_err)