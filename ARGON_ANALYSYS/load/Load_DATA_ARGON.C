{
        #include "TChain.h"
        #include "TCut.h"
        #include "TString.h"
        #include "TMath.h"
        #include "TH1F.h"
        #include "TFile.h"
        #include "math.h"
        #include <string>
        #include <cstring>
        #include <fstream>
        #include <vector>

        // --------------------------------------------------------------       DEFINIZIONE BIN DI ENERGIA bin uguali
       
        Int_t nbd = 6;
        Int_t ndec = 4;
        Double_t arg = 1./double(nbd);
        int noe = nbd*ndec;
        //cout << "Number of energy bin " << noe << endl;
        Double_t e0 = 100.;
        //Double_t Ebin[noe+1];
        std::vector<Double_t> Ebin(noe + 1);
        Ebin[0] = e0;
        for (int i=1; i<noe+1; i++) { Ebin[i] = Ebin[i-1]*TMath::Power(10., arg);
        }
        for (int i=1; i<noe+1; i++) {
            std::cout << "Bin " << i << ": "
              << Ebin[i-1] << " - " << Ebin[i] << std::endl;

        }
        
         /*
        Int_t n1 = 6;
        Int_t n2 = 6;
        Int_t n3 = 4;
        Int_t n4 = 4;

        Double_t E0 = 100.;
        Double_t E1 = 1e3;    // fine primo blocco
        Double_t E2 = 1e4;   // fine secondo blocco
        Double_t E3 = 1e5;  // fine terzo blocco
        Double_t E4 = 1e6;   // fine quarto blocco

        Int_t noe = n1 + n2 + n3 + n4;

        std::vector<Double_t> Ebin(noe + 1);
        Ebin[0] = E0;

        // ---------- blocco 1 ----------
        Double_t arg1 = TMath::Log10(E1/E0) / n1;

        for (int i = 1; i <= n1; i++)
             Ebin[i] = Ebin[i-1] * TMath::Power(10., arg1);

        // ---------- blocco 2 ----------
        Double_t arg2 = TMath::Log10(E2/E1) / n2;

        for (int i = n1+1; i <= n1+n2; i++)
             Ebin[i] = Ebin[i-1] * TMath::Power(10., arg2);

        // ---------- blocco 3 ----------
        Double_t arg3 = TMath::Log10(E3/E2) / n3;

        for (int i = n1+n2+1; i <= n1+n2+n3; i++)
             Ebin[i] = Ebin[i-1] * TMath::Power(10., arg3);

        // ---------- blocco 4  ----------
        //for (int i = n1+n2+n3+1; i <= noe; i++)
        //     Ebin[i] = Ebin[i-1] * TMath::Power(10., arg3);

        // ---------- blocco 4  ----------
        Double_t arg4 = TMath::Log10(E4/E3) / n4;

        for (int i = n1+n2+n3+1; i <= noe; i++)
             Ebin[i] = Ebin[i-1] * TMath::Power(10., arg4);

        // ---------- stampa ----------
        for (int i = 1; i <= noe; i++) {
             std::cout << "Bin " << i << ": "
                       << Ebin[i-1] << " - " << Ebin[i] << std::endl;
        }
        */
        
        
        // --------------------------------------------------------------       DATA LOADING...
        TChain *skim = new TChain("newtree");
        // ------- 2016 (Jan-Dec)
        skim->Add("/mnt/c/Users/saraf/Desktop/dampe/FD_SKIM/skim_FD_2016_merged.root");
        // ------- 2017 (Jan-Dec)
        skim->Add("/mnt/c/Users/saraf/Desktop/dampe/FD_SKIM/skim_FD_2017_merged.root");
        // ------- 2018 (Jan-Dec)
        skim->Add("/mnt/c/Users/saraf/Desktop/dampe/FD_SKIM/skim_FD_2018_merged.root");
        // ------- 2019 (Jan-Dec)
        skim->Add("/mnt/c/Users/saraf/Desktop/dampe/FD_SKIM/skim_FD_2019_merged.root");
        // ------- 2020 (Jan-Dec)
        skim->Add("/mnt/c/Users/saraf/Desktop/dampe/FD_SKIM/skim_FD_2020_merged.root");
        // ------- 2021 (Jan-Dec)
        skim->Add("/mnt/c/Users/saraf/Desktop/dampe/FD_SKIM/skim_FD_2021_merged.root");
        // ------- 2022 (Jan-Dec)
        skim->Add("/mnt/c/Users/saraf/Desktop/dampe/FD_SKIM/skim_FD_2022_merged.root");
        // ------- 2023 (Jan-Dec)
        skim->Add("/mnt/c/Users/saraf/Desktop/dampe/FD_SKIM/skim_FD_2023_merged.root");
        // ------- 2024 (Jan-Dec)
        skim->Add("/mnt/c/Users/saraf/Desktop/dampe/FD_SKIM/skim_FD_2024_merged.root");
        // ------- 2025 (Jan-Dec)
        skim->Add("/mnt/c/Users/saraf/Desktop/dampe/FD_SKIM/skim_FD_2025_merged.root");
        // ----------------------
        cout << " " << endl;
        cout << "------------------------------------------------------" << endl;
        cout << "  Number of entries (10 years / 120 months): " << skim->GetEntries() << endl;
        cout << "------------------------------------------------------" << endl;
        cout << " " << endl;


        // --------------------------------------------------------------       CUT DEFINITION...

        TCut cTrig_HET = "BGO_HET>0.";
        TCut cEne = "BGO_EnergyG_SatCorr_ML_ions>100."; // Min deposited energy [GeV]
        TCut cut00  =  cTrig_HET*cEne;
        TCut cut01  = "(PSD_ChargeY0>0.0 || PSD_ChargeY1>0.0) && (PSD_ChargeX0>0.0 || PSD_ChargeX1>0.0)";  //(Elisabetta 02-12-2021)
        TCut cut05  = "fabs(BGO_cbgomax[0]-BGO_cbgostk[0])<30.0 && fabs(BGO_cbgomax[1]-BGO_cbgostk[1])<30.0";
        TCut cut06 = "fabs(STKtrack_to_PSD_topY)< 400. && fabs(STKtrack_to_PSD_topX) < 400";
        // STK cut (Elisabetta 22-02-2024)
        TCut cutSTK1200 = "((((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2.*STK_chargeY_etaCorr[0]+(TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.*STK_chargeX_etaCorr[0])/((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2.+(TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.))>1200.)";

        // 06-11-2023 --> Provo ad implementare un taglio sul PSD scegliendo solo gli eventi che hanno una consistenza tra le due viste del PSD!
        //                                In particolare: impongo che la differenza tra la carica in una vista (singolo hit o media) abbia una differenza minore di 2 con l'altra vista.
        TCut cutPSD = "fabs((((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_ChargeY0+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_ChargeY1)/((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.))-(((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.*PSD_ChargeX0+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.*PSD_ChargeX1)/((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.)))<2.";

        // 27/05/2026

        TString MPV_DATA = "(16.762+(1.14769*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.264239*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))+(0.00423957*pow(log10(BGO_EnergyG_SatCorr_ML_ions),4)))";
        TString sigma_DATA = "(0.242923+(0.0369508*log10(BGO_EnergyG_SatCorr_ML_ions)))";

        // --------------------

        TString cutLow = "1.5";
        TString cutUp = "2.0";
        TString PSDcharge = "PSD_PathWeighted_Charge";
        //TCut chargeUP = PSDcharge+"<("+MPV_DATA+"+"+cutUp +"*"+sigma_DATA+")";
        //TCut chargeDW = PSDcharge+">("+MPV_DATA+"-"+cutLow+"*"+sigma_DATA+")";
        TCut chargeUP = (PSDcharge + "<(" + MPV_DATA + "+" + cutUp + "*" + sigma_DATA + ")").Data();
        TCut chargeDW = (PSDcharge + ">(" + MPV_DATA + "-" + cutLow + "*" + sigma_DATA + ")").Data();
        TCut chargeCUT = chargeUP*chargeDW;

        // -------   TAGLIO COMPLESSIVO   -------
        //TCut c0   = cut00;
        //TCut c01   = cut00*cut01;
        //TCut c05   = cut00*cut01*cut05;
        //TCut c06   = cut00*cut01*cut05*cut06;
        //TCut cstk   = cut00*cut01*cut05*cut06*cutSTK1200*;
        //TCut cpsd   = cut00*cut01*cut05*cut06*cutSTK1200*cutPSD;
        TCut ctot   = cut00*cut01*cut05*cut06*cutSTK1200*cutPSD*chargeCUT;
        //TCut ctot = cut00*cut01*cut05*cut06*cutSTK1200*cutPSD; 

        // --------------------------------------------------------------       FILL HISTOGRAMS...

        TFile *foutp = new TFile("out_load/Out_DATA_120months_ARGON_6bin.root","RECREATE");

        TH1D *h1Nobs = new TH1D("h1Nobs", "Nobs(Eo) -S- (skim)", noe, Ebin.data());
        //skim->Draw("BGO_EnergyG_SatCorr_ML_ions>>h1Nobs",c0,"");
        //skim->Draw("BGO_EnergyG_SatCorr_ML_ions>>h1Nobs",c01,"");
        //skim->Draw("BGO_EnergyG_SatCorr_ML_ions>>h1Nobs",c05,"");
        //skim->Draw("BGO_EnergyG_SatCorr_ML_ions>>h1Nobs",c06,"");
        //skim->Draw("BGO_EnergyG_SatCorr_ML_ions>>h1Nobs",cstk,"");
        //skim->Draw("BGO_EnergyG_SatCorr_ML_ions>>h1Nobs",cpsd,"");
        skim->Draw("BGO_EnergyG_SatCorr_ML_ions>>h1Nobs",ctot,"");
        //skim->Draw("BGO_EnergyG_SatCorr_ML[13]>>h1Nobs",ctot,"");

        /*
        double bkg[noe] = {
                0.0273069, 0.0280229, 0.0277414, 0.0264806,
                0.0252237, 0.0261614, 0.0228536, 0.0228536,
                0.0228659, 0.0228659, 0.0252042, 0.0252042,
                0.0222354, 0.0222354, 0.0222354, 0.0222354
        };

        double bkg_err[noe] = {
                0.000326105, 0.000437967, 0.000491302, 0.000602619,
                0.000784365, 0.00118562, 0.00126267,  0.00126267,
                0.00315496,  0.00315496,  0.00630662,  0.00630662,
                0.00760391,  0.00760391,  0.00760391,   0.00760391
        };

        TH1F *h1Nobs_nobkg = new TH1F("h1Nobs_nobkg", "Nobs(Eo) -S- (skim)", noe, Ebin);
        for (int i=1; i<noe+1; i++) {
                double counts = h1Nobs->GetBinContent(i);
                double bkg_counts = counts * bkg[i-1];
                h1Nobs_nobkg->SetBinContent(i,(counts-bkg_counts));
        }
        */

        cout << " -------------------------------------------" << endl;
        cout << "           Done!  THIS IS THE END!          " << endl;
        cout << " -------------------------------------------" << endl;

        // --------------------------------------------------------------       DATA SAVING...

        foutp->cd();
        foutp->Write();

}