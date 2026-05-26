
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

	// -------------------------------------------------------------- 	DEFINIZIONE BIN DI ENERGIA
void HET_efficiency_DATA_SULFUR() {

	const int nbd = 8;
	const int ndec = 4;
	Double_t arg = 1./double(nbd);
	const int noe = nbd*ndec;
	//cout << "Number of energy bin " << noe << endl;
	Double_t e0 = 100.;
	Double_t Ebin[noe+1];
	Ebin[0] = e0;
	for (int i=1; i<noe+1; i++) { Ebin[i] = Ebin[i-1]*TMath::Power(10., arg); }

	for (int i=1; i<noe+1; i++) {
            std::cout << "Bin " << i << ": "
              << Ebin[i-1] << " - " << Ebin[i] << std::endl;

        }

	// -------------------------------------------------------------- 	DATA LOADING...
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

	cout << " " << endl;
	cout << "--------------------------------------------------------" << endl;
	cout << "  Number of entries (10 years / 120 months): " << skim->GetEntries() << endl;
	cout << "--------------------------------------------------------" << endl;
	cout << " " << endl;

	TFile *foutp = new TFile("out_eff/HET_eff_DATA_SULFUR_120months_8bin.root","RECREATE");
	
	// -------------------------------------------------------------- 	CUT DEFINITION...

	TCut cTrig_HET = "BGO_HET>0.";
	
	TCut cTrig_LET="BGO_LET>0.";
	TCut cEne = "BGO_EnergyG_SatCorr_ML_ions>100."; // Min deposited energy [GeV]
	TCut cut00  =  cTrig_HET*cEne;
	TCut cut01  = "(PSD_ChargeY0>0.0 || PSD_ChargeY1>0.0) && (PSD_ChargeX0>0.0 || PSD_ChargeX1>0.0)";  //(Elisabetta 02-12-2021)   
	
	TCut cut05  = "fabs(BGO_cbgomax[0]-BGO_cbgostk[0])<30.0 && fabs(BGO_cbgomax[1]-BGO_cbgostk[1])<30.0";
	TCut cut06 = "fabs(STKtrack_to_PSD_topY)<400. && fabs(STKtrack_to_PSD_topX)<400."; 

	// STK cut (Elisabetta 22-02-2024)
	TCut cutSTK1200 = "((((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2.*STK_chargeY_etaCorr[0] + (TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.*STK_chargeX_etaCorr[0]) / ((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2. + (TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.) ) >1200.)";

	// 06-11-2023 --> Provo ad implementare un taglio sul PSD scegliendo solo gli eventi che hanno una consistenza tra le due viste del PSD!
	// 				  In particolare: impongo che la differenza tra la carica in una vista (singolo hit o media) abbia una differenza minore di 2 con l'altra vista.
	TCut cutPSD = "fabs(( ( (TMath::Sign(1.,PSD_ChargeY0) + 1.)/2. * PSD_ChargeY0 + (TMath::Sign(1.,PSD_ChargeY1) + 1.)/2. * PSD_ChargeY1 ) / ((TMath::Sign(1.,PSD_ChargeY0) + 1.)/2. + (TMath::Sign(1.,PSD_ChargeY1) + 1.)/2.) ) - ( ( (TMath::Sign(1.,PSD_ChargeX0) + 1.)/2. * PSD_ChargeX0 + (TMath::Sign(1.,PSD_ChargeX1) + 1.)/2. * PSD_ChargeX1 ) / ((TMath::Sign(1.,PSD_ChargeX0) + 1.)/2. + (TMath::Sign(1.,PSD_ChargeX1) + 1.)/2.) ) ) < 2.";
	
	///////////////////////////////// PSD charge - pathlength average
	
	// SULFUR 01/2026
	TString MPV_DATA = "(14.5821+(1.51344*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.48318*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))+(0.0520877*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2)))";
	//"(15.9434+(0.0730006*log10(BGO_EnergyG_SatCorr_ML_ions))+(0.00533945*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))+(-0.000744557*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2)))";
	//"(14.9723+(1.09723*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.342737*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))+(0.0375106*pow(log10(BGO_EnergyG_SatCorr_ML_ions),3)))";
	TString sigma_DATA =  "(0.20977+(0.0478004*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.00269758*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2)))";
	//"(0.355433+(-0.0688415*log10(BGO_EnergyG_SatCorr_ML_ions))+(0.019678*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2)))";
	//"(0.138753+(0.158506*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.0559774*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))+(0.00789775*pow(log10(BGO_EnergyG_SatCorr_ML_ions),3)))";

	//CALCIUM 01/2026
	//TString MPV_DATA_ca = "(20.0184+(0.0481385*log10(BGO_EnergyG_SatCorr_ML_ions)))";
    //TString sigma_DATA_ca = "(0.331059+(0.0108657*log10(BGO_EnergyG_SatCorr_ML_ions)))";
	
	TString cutLow = "2.0";
	TString cutUp = "2.5";
	TString PSDcharge = "PSD_PathWeighted_Charge";
	TCut chargeUP = TCut( (PSDcharge+"<("+MPV_DATA+"+"+cutUp +"*"+sigma_DATA+")").Data() );
	TCut chargeDW = TCut( (PSDcharge+">("+MPV_DATA+"-"+cutLow+"*"+sigma_DATA+")").Data() );
	TCut chargeCUT = chargeUP*chargeDW;

	// -------   TAGLI COMPLESSIVO   -------
	TCut ctot = cut00*cut01*cut05*cut06*cutSTK1200*cutPSD*chargeCUT;

	TCut cut_00s  = cEne*cTrig_LET;
	TCut cut00_ANDs  = cEne*cTrig_HET*cTrig_LET;

	TCut ctoto_0s      = cut_00s*cut01*cut05*cut06*cutSTK1200*cutPSD*chargeCUT;
	TCut ctoto_ANDs = cut00_ANDs*cut01*cut05*cut06*cutSTK1200*cutPSD*chargeCUT;

	// -------------------------------------------------------------- 	FILL HISTOGRAMS...

	TH1F *h1_Unb = new TH1F("h1_Unb", "Nobs(Eo) - LET-Trig", noe, Ebin); h1_Unb->Sumw2();
	skim->Draw("BGO_EnergyG_SatCorr_ML_ions>>h1_Unb",ctoto_0s,""); 

	TH1F *h1Nobs_And = new TH1F("h1Nobs_And", "Nobs(Eo) - LET*HET", noe, Ebin); h1Nobs_And->Sumw2();
	skim->Draw("BGO_EnergyG_SatCorr_ML_ions>>h1Nobs_And",ctoto_ANDs,""); 

	// ------- EFFICIENCY

	TH1F *h1_HET_Eff = new TH1F("h1_HET_Eff", "HET-Eff vs Ebgo", noe, Ebin); h1_HET_Eff->Sumw2();
	h1_HET_Eff->Divide(h1Nobs_And,h1_Unb);

	cout << " -------------------------------------------" << endl;
	cout << "           Done!  THIS IS THE END!          " << endl;
	cout << " -------------------------------------------" << endl;

	// -------------------------------------------------------------- 	DATA SAVING...

	foutp->cd();
	foutp->Write();
	foutp->ls();
	delete foutp; 

}
