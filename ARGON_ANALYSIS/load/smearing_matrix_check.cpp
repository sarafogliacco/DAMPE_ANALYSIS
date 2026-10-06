#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TAxis.h"
#include "TMath.h"
#include "TString.h"
#include <iostream>
#include <vector>
#include <cmath>

// =====================================================================
//  Macro: Ntrig(Et,Eo) con binnatura fine in un range ristretto.
//  Ngen viene ottenuto rebinnando h_energy_truth_120bins_weight_E2e7_cut_0
// =====================================================================

void load_MC_rebin()
{
    // ----------------------------------------------------------------- CONFIGURAZIONE
    const int    BINS_PER_DECADE = 20;      // 10 (oppure 20 se l'istogramma fine lo permette)
    const double EMIN_REQ        = 1e2;     // GeV  
    const double EMAX_REQ        = 1e4;     // GeV  (500 TeV)
    const double index           = 1.7;     // indice spettrale
    const TString histName       = "h_energy_truth_120bins_weight_E2e7_cut_0";
    const TString outName        = "/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/ARGON_ANALYSIS/load/out_load/Out_MC_FTFP_Ar_E2e7_20rebin.root";

    // fattore di correzione sul bin che contiene 500 TeV (file 3 e 4 si sovrappongono)
    //const double E_OVERLAP = 5e5;
    //const double SCALE_H03 = 3027360. / 1669300.;   // <-- DA RICONTROLLARE con la nuova binnatura

    const int NUM_SET = 4;
    std::vector<TString> files = {
        "/mnt/c/Users/saraf/Desktop/dampe/MC/MC_ARGON/skim_MC_allAr40_100GeV_1TeV_merged_v2.root",
        "/mnt/c/Users/saraf/Desktop/dampe/MC/MC_ARGON/skim_MC_allAr40_1TeV_10TeV_merged_v2.root",
        "/mnt/c/Users/saraf/Desktop/dampe/MC/MC_ARGON/skim_MC_allAr40_10TeV_100TeV_merged_v2.root",
        "/mnt/c/Users/saraf/Desktop/dampe/MC/MC_ARGON/skim_MC_allAr40_100TeV_1PeV_merged_v2.root",
    };
    // range in energia coperto da ciascun file (GeV)
    const double fileEmin[NUM_SET] = {1e2, 1e3, 1e4, 1e5};
    const double fileEmax[NUM_SET] = {1e3, 1e4, 1e5, 1e6};

    // ----------------------------------------------------------------- BIN DI ENERGIA
    // Griglia logaritmica: Ebin[i] = EMIN * 10^(i/BPD). Con 10 bin/decade il bordo
    // superiore piu' vicino a 500 TeV e' 10^5.7 GeV = 501.2 TeV (17 bin da 10 TeV).
    const double EMIN = EMIN_REQ;
    const int    noe  = TMath::Nint(BINS_PER_DECADE * TMath::Log10(EMAX_REQ / EMIN));
    std::vector<double> Ebin(noe + 1);
    for (int i = 0; i <= noe; i++) Ebin[i] = EMIN * TMath::Power(10, double(i) / BINS_PER_DECADE);
    const double EMAX = Ebin[noe];

    std::cout << "Number of energy bins: " << noe << "  (" << EMIN << " - " << EMAX << " GeV)" << std::endl;
    for (int i = 0; i <= noe; i++) std::cout << "  Ebin[" << i << "] = " << Ebin[i] << std::endl;

    // ----------------------------------------------------------------- FILE USATI
    std::vector<bool> used(NUM_SET, false);
    for (int f = 0; f < NUM_SET; f++)
        used[f] = !(fileEmax[f] <= EMIN || fileEmin[f] >= EMAX);

    TChain *skim[NUM_SET];
    TFile  *pf[NUM_SET];
    TH1F   *hfine[NUM_SET];
    for (int f = 0; f < NUM_SET; f++) { skim[f] = nullptr; pf[f] = nullptr; hfine[f] = nullptr; }

    std::cout << "\n-----------------------------------" << std::endl;
    for (int f = 0; f < NUM_SET; f++) {
        if (!used[f]) { std::cout << "file " << f << ": fuori range, saltato" << std::endl; continue; }

        skim[f] = new TChain("newtree");
        skim[f]->Add(files[f]);
        std::cout << "file " << f << " (" << fileEmin[f] << " - " << fileEmax[f]
                  << " GeV): " << skim[f]->GetEntries() << " entries" << std::endl;

        pf[f] = TFile::Open(files[f], "READ");
        if (!pf[f] || pf[f]->IsZombie()) {
            std::cerr << "Errore: non riesco ad aprire il file " << f << std::endl;
            return;
        }
        hfine[f] = (TH1F*) pf[f]->Get(histName);
        if (!hfine[f]) {
            std::cerr << "Errore: istogramma " << histName << " mancante nel file " << f << std::endl;
            return;
        }
    }
    std::cout << "-----------------------------------\n" << std::endl;

    // ----------------------------------------------------------------- GEOMETRIA REBIN
    int ref = -1;
    for (int f = 0; f < NUM_SET; f++) if (used[f]) { ref = f; break; }

    TAxis *ax       = hfine[ref]->GetXaxis();
    double decadi   = TMath::Log10(ax->GetXmax() / ax->GetXmin());
    double fineBPD  = hfine[ref]->GetNbinsX() / decadi;   // bin fini per decade
    double ratio    = fineBPD / BINS_PER_DECADE;          // bin fini per bin nuovo
    int    nMerge   = TMath::Nint(ratio);
    double offset   = TMath::Log10(EMIN / ax->GetXmin()) * fineBPD;

    std::cout << "Fine hist: " << hfine[ref]->GetNbinsX() << " bins, " << fineBPD
              << " bins/decade, range " << ax->GetXmin() << " - " << ax->GetXmax()
              << " GeV, merge factor = " << ratio << std::endl;

    if (fabs(ratio - nMerge) > 1e-6 || nMerge < 1) {
        std::cerr << "Errore: " << fineBPD << " bin/decade non rebinnabili a "
                  << BINS_PER_DECADE << " bin/decade (prova un altro valore)" << std::endl;
        return;
    }
    if (fabs(offset - TMath::Nint(offset)) > 1e-6 || EMIN < ax->GetXmin() || EMAX > ax->GetXmax() * (1. + 1e-6)) {
        std::cerr << "Errore: EMIN/EMAX non allineati ai bin dell'istogramma fine" << std::endl;
        return;
    }
    // controllo che tutti gli istogrammi abbiano la stessa binnatura
    for (int f = 0; f < NUM_SET; f++) {
        if (!used[f]) continue;
        if (hfine[f]->GetNbinsX() != hfine[ref]->GetNbinsX() ||
            fabs(hfine[f]->GetXaxis()->GetXmin() - ax->GetXmin()) > 1e-6 * ax->GetXmin() ||
            fabs(hfine[f]->GetXaxis()->GetXmax() - ax->GetXmax()) > 1e-6 * ax->GetXmax()) {
            std::cerr << "Errore: binnatura diversa tra gli istogrammi fini (file " << f << ")" << std::endl;
            return;
        }
    }

    // ----------------------------------------------------------------- OUTPUT E Ngen
    TFile *fout = new TFile(outName, "RECREATE");
    fout->cd();

    TH1F *h1Ngen = new TH1F("h1Ngen", "Ngen(Et)", noe, Ebin.data());

    const double tol = 1e-2;   // tolleranza per il test "il bin sta dentro il file"
    for (int j = 0; j < noe; j++) {
        int b0 = ax->FindBin(Ebin[j] * (1. + 1e-9));   // primo bin fine del bin j
        //bool overlap = (Ebin[j] < E_OVERLAP && E_OVERLAP < Ebin[j + 1]);

        double ngen = 0.;
        for (int f = 0; f < NUM_SET; f++) {
            if (!used[f]) continue;

            bool covers = (Ebin[j]     >= fileEmin[f] * (1. - tol) &&
                           Ebin[j + 1] <= fileEmax[f] * (1. + tol));
            //bool inOverlap = overlap && (f == 3 || f == 4);
            //if (!covers && !inOverlap) continue;

            double s = 0.;
            for (int b = b0; b < b0 + nMerge; b++) s += hfine[f]->GetBinContent(b);
            //if (f == 3 && overlap) s *= SCALE_H03;
            ngen += s;
        }
        h1Ngen->SetBinContent(j + 1, ngen);
    }

    std::vector<double> wNgen(noe);
    for (int i = 0; i < noe; i++) {
        double ngen = h1Ngen->GetBinContent(i + 1);
        wNgen[i] = (ngen > 0.0) ? 1.0 / ngen : 0.0;
        std::cout << "Bin " << i << " [" << Ebin[i] << ", " << Ebin[i + 1] << "]  Ngen = " << ngen << std::endl;
    }

    // ----------------------------------------------------------------- PESI
    // Peso integrale del bin i
    auto WeightIntegral = [&](int i) {
        double E0 = Ebin[i];
        double E1 = Ebin[i + 1];
        return index * log(E1 / E0) / (pow(E0, -index) - pow(E1, -index));
    };
    // Peso per evento: E^-index
    auto WeightSpectral = [&](double E) { return pow(E, -index); };

    // ----------------------------------------------------------------- FUNZIONI DI CARICA
    auto MPV_MC = [](double E) {
        double x = log10(E);
        return 18.6749 + 0.182178*x;
    };
    auto sigma_MC = [](double E) {
        double x = log10(E);
        return 0.131044 + 0.0168788*x;
    };
    auto MPV_DATA = [](double E) {
        double x = log10(E);
        return 17.905 + 0.0985843*x;
    };
    auto sigma_DATA = [](double E) {
        double x = log10(E);
        return 0.300046 + 0.035818*x;
    };
    auto MPV_MC_CORR = [](double E) {
        double x = log10(E);
        return 17.9623 + 0.0942476*x;
    };
    auto sigma_MC_CORR = [](double E) {
        double x = log10(E);
        return 0.303748 + 0.0348213*x;
    };

    // ----------------------------------------------------------------- ISTOGRAMMI 2D
    fout->cd();
    TH2F *h2Ntrig_cut00  = new TH2F("h2Ntrig_cut00",  "Ntrig(Eo,Et) cut00",  noe, Ebin.data(), noe, Ebin.data());
    TH2F *h2Ntrig_cut01  = new TH2F("h2Ntrig_cut01",  "Ntrig(Eo,Et) cut01",  noe, Ebin.data(), noe, Ebin.data());
    TH2F *h2Ntrig_cut05  = new TH2F("h2Ntrig_cut05",  "Ntrig(Eo,Et) cut05",  noe, Ebin.data(), noe, Ebin.data());
    TH2F *h2Ntrig_cut06  = new TH2F("h2Ntrig_cut06",  "Ntrig(Eo,Et) cut06",  noe, Ebin.data(), noe, Ebin.data());
    TH2F *h2Ntrig_stk    = new TH2F("h2Ntrig_stk",    "Ntrig(Eo,Et) stk",    noe, Ebin.data(), noe, Ebin.data());
    TH2F *h2Ntrig_psd    = new TH2F("h2Ntrig_psd",    "Ntrig(Eo,Et) psd",    noe, Ebin.data(), noe, Ebin.data());
    TH2F *h2Ntrig_charge = new TH2F("h2Ntrig_charge", "Ntrig(Eo,Et) charge", noe, Ebin.data(), noe, Ebin.data());
    TH2F *h2Ntrig_wgt    = new TH2F("h2Ntrig_wgt",    "Ntrig(Eo,Et) -S- (skim)", noe, Ebin.data(), noe, Ebin.data());

    // ----------------------------------------------------------------- LOOP SUGLI EVENTI
    for (int i = 0; i < NUM_SET; i++) {
        if (!used[i]) continue;

        std::cout << "... looping on events, set " << i << " ..." << std::endl;

        double MC_EnergyT;
        double BGO_Energy;
        int    BGO_HET;

        double PSD_ChargeY0, PSD_ChargeY1;
        double PSD_ChargeX0, PSD_ChargeX1;

        Double_t BGO_cbgomax[14];
        Double_t BGO_cbgostk[14];

        double STKtrack_to_PSD_topY;
        double STKtrack_to_PSD_topX;

        Double_t STK_chargeY_etaCorr[6];
        Double_t STK_chargeX_etaCorr[6];

        double PSD_PW_Charge;

        skim[i]->SetBranchAddress("MC_EnergyT", &MC_EnergyT);
        skim[i]->SetBranchAddress("BGO_EnergyG_QuenchSatCorr_ML_ions2", &BGO_Energy);
        skim[i]->SetBranchAddress("BGO_HET", &BGO_HET);

        skim[i]->SetBranchAddress("PSD_ChargeY0", &PSD_ChargeY0);
        skim[i]->SetBranchAddress("PSD_ChargeY1", &PSD_ChargeY1);
        skim[i]->SetBranchAddress("PSD_ChargeX0", &PSD_ChargeX0);
        skim[i]->SetBranchAddress("PSD_ChargeX1", &PSD_ChargeX1);

        skim[i]->SetBranchAddress("BGO_cbgomax", BGO_cbgomax);
        skim[i]->SetBranchAddress("BGO_cbgostk", BGO_cbgostk);

        skim[i]->SetBranchAddress("STKtrack_to_PSD_topY", &STKtrack_to_PSD_topY);
        skim[i]->SetBranchAddress("STKtrack_to_PSD_topX", &STKtrack_to_PSD_topX);

        skim[i]->SetBranchAddress("STK_chargeY_etaCorr", STK_chargeY_etaCorr);
        skim[i]->SetBranchAddress("STK_chargeX_etaCorr", STK_chargeX_etaCorr);

        skim[i]->SetBranchAddress("PSD_PathWeighted_Charge", &PSD_PW_Charge);

        Long64_t nentries = skim[i]->GetEntries();

        for (Long64_t ev = 0; ev < nentries; ev++) {

            skim[i]->GetEntry(ev);
            if (ev % 100000 == 0) std::cout << "event " << ev << " out of " << nentries << std::endl;

            // ---- range in energia vera
            if (MC_EnergyT < EMIN || MC_EnergyT >= EMAX) continue;

            // ---- cut 00
            if (BGO_Energy < 100.) continue;
            if (BGO_HET <= 0) continue;

            // ---- bin di energia vera e peso totale
            int j = h1Ngen->FindBin(MC_EnergyT) - 1;
            if (j < 0 || j >= noe) continue;

            double w_tot = wNgen[j] * WeightIntegral(j) * WeightSpectral(MC_EnergyT);

            h2Ntrig_cut00->Fill(MC_EnergyT, BGO_Energy, w_tot);

            // ---- cut 01
            bool hasY = (PSD_ChargeY0 > 0.0 || PSD_ChargeY1 > 0.0);
            bool hasX = (PSD_ChargeX0 > 0.0 || PSD_ChargeX1 > 0.0);
            if (!(hasY && hasX)) continue;
            h2Ntrig_cut01->Fill(MC_EnergyT, BGO_Energy, w_tot);

            // ---- cut 05
            if (fabs(BGO_cbgomax[0] - BGO_cbgostk[0]) >= 30.0) continue;
            if (fabs(BGO_cbgomax[1] - BGO_cbgostk[1]) >= 30.0) continue;
            h2Ntrig_cut05->Fill(MC_EnergyT, BGO_Energy, w_tot);

            // ---- cut 06
            if (fabs(STKtrack_to_PSD_topY) >= 400.) continue;
            if (fabs(STKtrack_to_PSD_topX) >= 400.) continue;
            h2Ntrig_cut06->Fill(MC_EnergyT, BGO_Energy, w_tot);

            // ---- stk cut
            double wY = (TMath::Sign(1., STK_chargeY_etaCorr[0]) + 1.) / 2.;
            double wX = (TMath::Sign(1., STK_chargeX_etaCorr[0]) + 1.) / 2.;
            if (wY + wX <= 0.) continue;
            double stk_charge = (wY * STK_chargeY_etaCorr[0] + wX * STK_chargeX_etaCorr[0]) / (wY + wX);
            if (stk_charge <= 1200.) continue;
            h2Ntrig_stk->Fill(MC_EnergyT, BGO_Energy, w_tot);

            // ---- psd cut
            double wY0 = (TMath::Sign(1., PSD_ChargeY0) + 1.) / 2.;
            double wY1 = (TMath::Sign(1., PSD_ChargeY1) + 1.) / 2.;
            double wX0 = (TMath::Sign(1., PSD_ChargeX0) + 1.) / 2.;
            double wX1 = (TMath::Sign(1., PSD_ChargeX1) + 1.) / 2.;
            double psdY = (wY0 * PSD_ChargeY0 + wY1 * PSD_ChargeY1) / (wY0 + wY1);
            double psdX = (wX0 * PSD_ChargeX0 + wX1 * PSD_ChargeX1) / (wX0 + wX1);
            if (fabs(psdY - psdX) >= 2.0) continue;
            h2Ntrig_psd->Fill(MC_EnergyT, BGO_Energy, w_tot);

            // ---- charge cut
            double QcorrMC = ((PSD_PW_Charge - MPV_MC(BGO_Energy)) * (sigma_DATA(BGO_Energy) / sigma_MC(BGO_Energy)))
                             + MPV_DATA(BGO_Energy);

            double cutLow     = 1.0;
            double cutUp      = 2.0;
            double mpv_corr   = MPV_MC_CORR(BGO_Energy);
            double sigma_corr = sigma_MC_CORR(BGO_Energy);

            if (QcorrMC <= mpv_corr - cutLow * sigma_corr) continue;
            if (QcorrMC >= mpv_corr + cutUp  * sigma_corr) continue;
            h2Ntrig_charge->Fill(MC_EnergyT, BGO_Energy, w_tot);

            // ---- 2D finale
            h2Ntrig_wgt->Fill(MC_EnergyT, BGO_Energy, w_tot);
        }
    }

    // ----------------------------------------------------------------- SALVATAGGIO
    std::cout << " -------------------------------------------" << std::endl;
    std::cout << "       Done!  THIS IS THE END!  MIAO        " << std::endl;
    std::cout << " -------------------------------------------" << std::endl;

    fout->cd();
    fout->Write();
    fout->Close();
}