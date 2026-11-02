#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <iostream>
#include <iomanip>

// Z e A effettivi del gruppo S-Ar-Ca pesati con i flussi, CON propagazione degli errori.
//
// caOverAr  = F_Ca / F_Ar (da letteratura). Con 0 il Ca viene ignorato.
// dCaOverAr = incertezza su caOverAr (0 = trascurata).
// isDifferential = true se "flux" e' un flusso differenziale (per unita' di energia):
//                  in tal caso la media globale pesa con la larghezza del bin.
//                  Se "flux" e' gia' integrato nel bin, mettere false.
//
// Ipotesi: errori di S e Ar indipendenti, bin indipendenti (le correlazioni dell'unfolding
// sono ignorate), GetBinError = incertezza totale che si vuole propagare.
void group_ZA_err(const char* fAr = "/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/ARGON_ANALYSIS/unfold_results/unfold_results_argon_6bin_low1_up2_v2.root",
                  const char* fS  = "/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/SULFUR_ANALYSIS/unfold_results/unfold_results_flux_spectrum_sulfur_120months_6bin_smooth_2e5up_1e5low_bkg_v2.root",
                  double caOverAr = 0.0, double dCaOverAr = 0.0, bool isDifferential = true)
{
    const double ZS = 16, ZAr = 18, ZCa = 20;
    const double AS = 32.06, AAr = 39.948, ACa = 40.5;
    const double r = caOverAr;

    TFile* s = TFile::Open(fS);  TFile* a = TFile::Open(fAr);
    if (!s || s->IsZombie() || !a || a->IsZombie()) { std::cerr << "Errore apertura file\n"; return; }
    TH1D* hS  = (TH1D*)s->Get("flux");
    TH1D* hAr = (TH1D*)a->Get("flux");
    if (!hS || !hAr || hS->GetNbinsX() != hAr->GetNbinsX()) { std::cerr << "Istogrammi non validi\n"; return; }

    std::cout << std::setw(10) << "Elow" << std::setw(10) << "Eup"
              << std::setw(9) << "wS" << std::setw(9) << "wAr" << std::setw(9) << "wCa"
              << std::setw(9) << "Z_eff" << std::setw(9) << "err_Z" << std::setw(9) << "sig_Z"
              << std::setw(9) << "A_eff" << std::setw(9) << "err_A" << std::setw(9) << "sig_A" << "\n";

    // Accumulatori per la media globale
    double Ptot = 0;                 // sum_i T_i * w_i
    double NZ = 0, NA = 0;           // sum_i (sum_k Z_k F_k) * w_i   (e analogo per A)
    double Z2N = 0, A2N = 0;         // per la dispersione globale
    double SfaW = 0;                 // sum_i w_i * fa_i  (per la derivata rispetto a r)
    const int nb = hS->GetNbinsX();
    std::vector<double> vw(nb + 1, 0.), vfs(nb + 1, 0.), vfa(nb + 1, 0.);
    std::vector<double> vefs(nb + 1, 0.), vefa(nb + 1, 0.);
    std::vector<bool> ok(nb + 1, false);

    for (int i = 1; i <= nb; i++) {
        double fs = hS->GetBinContent(i),  efs = hS->GetBinError(i);
        double fa = hAr->GetBinContent(i), efa = hAr->GetBinError(i);
        if (fs <= 0 || fa <= 0) continue;
        double fc  = r * fa;
        double w   = isDifferential ? hS->GetBinWidth(i) : 1.0;
        double tot = fs + fa + fc;
        double wS = fs / tot, wA = fa / tot, wC = fc / tot;

        double Z  = wS*ZS + wA*ZAr + wC*ZCa;
        double Z2 = wS*ZS*ZS + wA*ZAr*ZAr + wC*ZCa*ZCa;
        double A  = wS*AS + wA*AAr + wC*ACa;
        double A2 = wS*AS*AS + wA*AAr*AAr + wC*ACa*ACa;

        // Derivate parziali rispetto alle variabili indipendenti fs, fa (con F_Ca = r*fa):
        //   dZ/dfs = (ZS - Z)/tot ;  dZ/dfa = ((ZAr + r*ZCa) - (1+r)*Z)/tot ;  dZ/dr = fa*(ZCa - Z)/tot
        double dZ_s = (ZS - Z) / tot;
        double dZ_a = ((ZAr + r*ZCa) - (1 + r)*Z) / tot;
        double dZ_r = fa * (ZCa - Z) / tot;
        double dA_s = (AS - A) / tot;
        double dA_a = ((AAr + r*ACa) - (1 + r)*A) / tot;
        double dA_r = fa * (ACa - A) / tot;

        double errZ = TMath::Sqrt(dZ_s*dZ_s*efs*efs + dZ_a*dZ_a*efa*efa + dZ_r*dZ_r*dCaOverAr*dCaOverAr);
        double errA = TMath::Sqrt(dA_s*dA_s*efs*efs + dA_a*dA_a*efa*efa + dA_r*dA_r*dCaOverAr*dCaOverAr);

        std::cout << std::fixed << std::setprecision(1)
                  << std::setw(10) << hS->GetBinLowEdge(i) << std::setw(10) << hS->GetBinLowEdge(i) + hS->GetBinWidth(i)
                  << std::setprecision(3)
                  << std::setw(9) << wS << std::setw(9) << wA << std::setw(9) << wC
                  << std::setw(9) << Z << std::setw(9) << errZ << std::setw(9) << TMath::Sqrt(TMath::Max(0., Z2 - Z*Z))
                  << std::setw(9) << A << std::setw(9) << errA << std::setw(9) << TMath::Sqrt(TMath::Max(0., A2 - A*A)) << "\n";

        Ptot += tot * w;
        NZ   += Z * tot * w;   Z2N += Z2 * tot * w;
        NA   += A * tot * w;   A2N += A2 * tot * w;
        SfaW += fa * w;
        vw[i] = w; vfs[i] = fs; vfa[i] = fa; vefs[i] = efs; vefa[i] = efa; ok[i] = true;
    }
    if (Ptot <= 0) return;

    double Zg = NZ / Ptot, Ag = NA / Ptot;

    // Errore globale: somma in quadratura sui bin (indipendenti) delle derivate
    //   dZg/dfs_i = w_i (ZS - Zg)/Ptot ,  dZg/dfa_i = w_i ((ZAr + r ZCa) - (1+r) Zg)/Ptot
    double varZ = 0, varA = 0;
    for (int i = 1; i <= nb; i++) {
        if (!ok[i]) continue;
        double gZs = vw[i] * (ZS - Zg) / Ptot;
        double gZa = vw[i] * ((ZAr + r*ZCa) - (1 + r)*Zg) / Ptot;
        double gAs = vw[i] * (AS - Ag) / Ptot;
        double gAa = vw[i] * ((AAr + r*ACa) - (1 + r)*Ag) / Ptot;
        varZ += gZs*gZs*vefs[i]*vefs[i] + gZa*gZa*vefa[i]*vefa[i];
        varA += gAs*gAs*vefs[i]*vefs[i] + gAa*gAa*vefa[i]*vefa[i];
    }
    // Contributo di r: completamente correlato tra i bin -> si somma una sola volta
    double gZr = (ZCa - Zg) * SfaW / Ptot;
    double gAr = (ACa - Ag) * SfaW / Ptot;
    varZ += gZr*gZr*dCaOverAr*dCaOverAr;
    varA += gAr*gAr*dCaOverAr*dCaOverAr;

    std::cout << "\nMedia pesata sul flusso" << (isDifferential ? " integrato (flusso differenziale x larghezza bin)" : " (flusso gia' integrato nel bin)") << ":\n"
              << "  Z_eff = " << Zg << " +- " << TMath::Sqrt(varZ)
              << "   (dispersione tra specie sigma = " << TMath::Sqrt(TMath::Max(0., Z2N/Ptot - Zg*Zg)) << ")\n"
              << "  A_eff = " << Ag << " +- " << TMath::Sqrt(varA)
              << "   (dispersione tra specie sigma = " << TMath::Sqrt(TMath::Max(0., A2N/Ptot - Ag*Ag)) << ")\n";

    s->Close(); a->Close();
}