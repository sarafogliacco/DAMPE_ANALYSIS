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
       /*
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
        */

     std::vector<int> binsPerDecade = {6, 6, 4, 4};
     
     double e0 = 100.0;
     std::vector<double> Ebin;
     Ebin.push_back(e0);
     
     double currentE = e0;
     
     for (size_t idec = 0; idec < binsPerDecade.size(); ++idec) {
     
         int nbin = binsPerDecade[idec];
         double factor = TMath::Power(10.0, 1.0 / nbin);
     
         for (int ib = 0; ib < nbin; ++ib) {
             currentE *= factor;
             Ebin.push_back(currentE);
         }
     }
     
     for (size_t i = 1; i < Ebin.size(); ++i) {
         std::cout << "Bin " << i
                   << ": " << Ebin[i-1]
                   << " - " << Ebin[i]
                   << std::endl;
     }

     int noe = 0;
     for (auto n : binsPerDecade)
     noe += n;
         std::cout << "Total bins = " << noe << std::endl;

        
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

        // 08/06/2026
        // SULFUR LOWER CUT
        TString MPV_DATA_S =  "(14.5821+(1.51344*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.48318*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))+(0.0520877*pow(log10(BGO_EnergyG_SatCorr_ML_ions),3)))";

        TString sigma_DATA_S = "(0.231952+(0.0318629*log10(BGO_EnergyG_SatCorr_ML_ions)))";

        // CALCIUM UPPER CUT
        TString MPV_DATA_Ca =  "(17.1818+(2.97939*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.976515*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))+(0.10604*pow(log10(BGO_EnergyG_SatCorr_ML_ions),3)))";

        TString sigma_DATA_Ca = "(0.206041+(0.0693299*log10(BGO_EnergyG_SatCorr_ML_ions)))";

        // --------------------

        TString cutLow = "2.0";
        TString cutUp = "2.0";
        TString PSDcharge = "PSD_PathWeighted_Charge";
        //TCut chargeUP = PSDcharge+"<("+MPV_DATA+"+"+cutUp +"*"+sigma_DATA+")";
        //TCut chargeDW = PSDcharge+">("+MPV_DATA+"-"+cutLow+"*"+sigma_DATA+")";
        TCut chargeUP = (PSDcharge + "<(" + MPV_DATA_Ca + "+" + cutUp + "*" + sigma_DATA_Ca + ")").Data();
        TCut chargeDW = (PSDcharge + ">(" + MPV_DATA_S + "-" + cutLow + "*" + sigma_DATA_S + ")").Data();
        TCut chargeCUT = chargeUP*chargeDW;

        // -------   TAGLIO COMPLESSIVO   -------
        TCut ctot   = cut00*cut01*cut05*cut06*cutSTK1200*cutPSD*chargeCUT;
        //TCut ctot = cut00*cut01*cut05*cut06*cutSTK1200*cutPSD; 

        // --------------------------------------------------------------       FILL HISTOGRAMS...

        TFile *foutp = new TFile("../out_root/Out_DATA_120months_SArCa_6-4bin.root","RECREATE");

        TH1D *h1Nobs = new TH1D("h1Nobs", "Nobs(Eo) -S- (skim)", noe, Ebin.data());
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
	cout << "       Done!  THIS IS THE END!  MIAO        " << endl;
	cout << " -------------------------------------------" << endl;

        // --------------------------------------------------------------       DATA SAVING...

        foutp->cd();
        foutp->Write();

}