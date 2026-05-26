import uproot
import numpy as np
import matplotlib.pyplot as plt
'''
file = uproot.open("/mnt/c/Users/saraf/Desktop/dampe/MC/MC_SULFUR/skim_MC_allS32_100GeV_1TeV.root")
tree = file["newtree"]   

PSD_charge = tree["PSD_PathWeighted_Charge"].array(library="np")

#print("PSD charge:", PSD_charge)

plt.hist(PSD_charge, bins=100, histtype="step")
plt.xlabel("PSD Charge")
plt.ylabel("Events")
#plt.title("Distribuzione PSD Charge")
plt.show()
'''

files = ["/mnt/c/Users/saraf/Desktop/dampe/MC/MC_SULFUR/skim_MC_allS32_100GeV_1TeV.root", 
         "/mnt/c/Users/saraf/Desktop/dampe/MC/MC_SILICON/skim_new/skim_MC_allSi28-v6r0p15_100GeV_1TeV_FTFP_merged.root", 
         "/mnt/c/Users/saraf/Desktop/dampe/MC/MC_ARGON/skim_MC_allAr40_100GeV_1TeV.root"
         ]

all_PSD_charge = []

for f in files:
    file = uproot.open(f)
    tree = file["newtree"]   # nome del TTree

    PSD_charge = tree["PSD_PathWeighted_Charge"].array(library="np")
    all_PSD_charge.append(PSD_charge)


plt.hist(all_PSD_charge[0], bins=1000, histtype="step", label="Sulfur")
plt.hist(all_PSD_charge[1], bins=1000, histtype="step", label="Silicon")
plt.hist(all_PSD_charge[2], bins=1000, histtype="step", label="Argon")

plt.xlim(12, 20)

plt.legend()
plt.xlabel("PSD Charge")
plt.ylabel("Events")

plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()


plt.show()