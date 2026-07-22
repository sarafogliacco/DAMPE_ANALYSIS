import numpy as np
import matplotlib.pyplot as plt 
import uproot as up
import awkward as ak


file = up.open("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/ARGON_ANALYSYS/out_root/PSDcharge_perCUT_DATA_120months_SAr_v0.root")

#print(file.keys())
#print(file.classnames())

all_hist = file["hqpsdPathWeight_all"].to_numpy()
cut00_hist = file["hqpsdPathWeight_cut00"].to_numpy()   
cut01_hist = file["hqpsdPathWeight_cut01"].to_numpy()
cut05_hist = file["hqpsdPathWeight_cut05"].to_numpy()
cut06_hist = file["hqpsdPathWeight_cut06"].to_numpy()
stkcut_hist = file["hqpsdPathWeight_stkcut"].to_numpy()
psdcut_hist = file["hqpsdPathWeight_psdcut"].to_numpy()
cutS_hist = file["hqpsdPathWeight_cutS"].to_numpy()
cutAr_hist = file["hqpsdPathWeight_cutAr"].to_numpy()
#cutCa_hist = file["hqpsdPathWeight_cutCa"].to_numpy()

values, edges = file["hqpsdPathWeight_cut00"].to_numpy()

plt.figure(figsize=(10,6))

def plot_hist(name, label, color):
    values, edges = file[name].to_numpy()
    centers = 0.5 * (edges[:-1] + edges[1:])
    plt.step(centers, values, where="mid", label=label, color=color)

plot_hist("hqpsdPathWeight_cut00", "Pre-selection + HET", color='k')
plot_hist("hqpsdPathWeight_cut01", "PSD hits on YZ and XZ views", color="#34B800")
plot_hist("hqpsdPathWeight_stkcut", "STK 1st layer > 1200 ADC", color="#26a1ff")
plot_hist("hqpsdPathWeight_psdcut", r"|PSD$_Y$-PSD$_X$| < 2", color="#ff751f")
#plot_hist("hqpsdPathWeight_cutS", "S cut", color="#626262")
#plot_hist("hqpsdPathWeight_cutAr", "Ar cut", color="#626262")

plt.axvspan(15.5,21, color="#ffb81f", alpha=0.3)

plt.text(0.7, 700000, "H", fontsize=10)
plt.text(1.7, 500000, "He", fontsize=10)
plt.text(3.0, 50000, "Li", fontsize=10)
plt.text(3.7, 40000, "Be", fontsize=10)
plt.text(4.9, 40000, "B", fontsize=10)
plt.text(5.7, 200000, "C", fontsize=10)
plt.text(6.7, 60000, "N", fontsize=10)
plt.text(7.7, 300000, "O", fontsize=10)

plt.text(9.8, 50000, "Ne", fontsize=14)
plt.text(11.7, 60000, "Mg", fontsize=14)
plt.text(14.0, 60000, "Si", fontsize=14)
plt.text(16.0, 20000, "S", fontsize=14)
plt.text(17.8, 10000, "Ar", fontsize=14)
plt.text(19.8, 20000, "Ca", fontsize=14)

plt.text(26.0, 70000, "Fe", fontsize=10)
plt.text(28.0, 6000, "Ni", fontsize=10)

plt.ylim(200, 1.5e6)
plt.xlabel("PSD charge", fontsize=14)
plt.ylabel("events", fontsize=14)
plt.yscale("log")
plt.legend( loc="upper right", fontsize=12)
#plt.grid(True, which="both", ls="--", lw=0.5)
#plt.savefig("PSD_charge_perCUT_DATA_120months_220626.png", dpi=300)
plt.show()