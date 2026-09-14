import numpy as np
import matplotlib.pyplot as plt


def load_bkg(filename):
        E_min, E_max, C, C_err, O, O_err, Ne, Ne_err, Mg, Mg_err, Si, Si_err, Ar, Ar_err, Ca, Ca_err, Tot, Tot_err = np.genfromtxt(
        filename,
        skip_header=1,
        delimiter=',',
        usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),
        unpack=True
    )
        return {
        "E_min": E_min,
        "E_max": E_max,
        "C": C,
        "C_err": C_err,
        "O": O,
        "O_err": O_err,
        "Ne": Ne,
        "Ne_err": Ne_err,
        "Mg": Mg,
        "Mg_err": Mg_err,
        "Si": Si,
        "Si_err": Si_err,
        "Ar": Ar,
        "Ar_err": Ar_err,
        "Ca": Ca,
        "Ca_err": Ca_err,
        "Tot": Tot,
        "Tot_err": Tot_err
    }

#data = load_bkg("txt/sulfur/bkg_est.csv")
data = load_bkg("txt/argon/bkg_est.csv")

sum = data["C"] + data["O"] + data["Ne"] + data["Mg"] + data["Si"] + data["Ar"] + data["Ca"]
sum_err = np.sqrt(data["C_err"]**2 + data["O_err"]**2 + data["Ne_err"]**2 + data["Mg_err"]**2 + data["Si_err"]**2 + data["Ar_err"]**2 + data["Ca_err"]**2)

plt.figure(figsize=(10, 6))

energy_bins = (data["E_min"] + data["E_max"]) / 2
#E = ([121.1527961, 177.8277918, 261.015534 , 383.1188748, 562.3413673, 825.4041434, 1467.79767 , 3162.274615, 6812.921547]) 

def plot_with_errorbars(x, y, yerr, label, color):
    plt.errorbar(x, y*100, yerr=yerr*100, marker='.', label=label, color=color)

plot_with_errorbars(energy_bins, data["C"],  data["C_err"],  label='C', color='blue')
plot_with_errorbars(energy_bins, data["O"],  data["O_err"],  label='O', color='orange')
plot_with_errorbars(energy_bins, data["Ne"], data["Ne_err"], label='Ne', color='green')
plot_with_errorbars(energy_bins, data["Mg"], data["Mg_err"], label='Mg', color='red')
plot_with_errorbars(energy_bins, data["Si"], data["Si_err"], label='Si', color='purple')
plot_with_errorbars(energy_bins, data["Ar"], data["Ar_err"], label='S', color='brown')
plot_with_errorbars(energy_bins, data["Ca"], data["Ca_err"], label='Ca', color='magenta')

plt.errorbar(energy_bins, data["Tot"]*100, yerr=data["Tot_err"]*100, marker='.', linestyle='--', label='Total', color='black')
#plt.errorbar(E, sum*100, yerr=sum_err*100, marker='.', linestyle=':', label='Sum', color='grey')

'''
plt.figure(figsize=(10, 6))
def plot(x, y, label, color):
    plt.plot(x, y*100, label=label, color=color, marker='.')

plot(energy_bins, data["C"],  label='C', color='blue')
plot(energy_bins, data["O"],  label='O', color='orange')
plot(energy_bins, data["Ne"], label='Ne', color='green')
plot(energy_bins, data["Mg"], label='Mg', color='red')
plot(energy_bins, data["Si"], label='Si', color='purple')
plot(energy_bins, data["Ar"], label='Ar', color='brown')
plot(energy_bins, data["Ca"], label='Ca', color='magenta')

plt.plot(energy_bins, data["Tot"]*100, label='Total', color='black', linestyle='--')
'''

plt.xscale('log')
#plt.yscale('log')

plt.xlabel('DepositedEnergy (GeV)')
plt.ylabel('Background Rate %')

plt.xlim(90, 3*1e4)

plt.text(0.05, 0.9, 'ARGON',
         transform=plt.gca().transAxes,
         fontsize=20,
         color='grey',
         verticalalignment='top')

plt.legend()
plt.grid(True,  linestyle="--", alpha=0.5)

plt.savefig("PLOTS/bkg_est_argon.png", dpi=300, bbox_inches='tight')
plt.show()
