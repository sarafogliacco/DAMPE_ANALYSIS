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

	// -------------------------------------------------------------- 	DEFINIZIONE BIN DI ENERGIA
	const int nbd = 8;
	const int ndec = 4;
	Double_t arg = 1./double(nbd);
	const int noe = nbd*ndec;
	//cout << "Number of energy bin " << noe << endl;
	Double_t e0 = 100.;
	Double_t Ebin[noe+1];
	Ebin[0] = e0;
	for (int i=1; i<noe+1; i++) { Ebin[i] = Ebin[i-1]*TMath::Power(10., arg); }

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
	// ----------------------
	cout << " " << endl;
	cout << "--------------------------------------------------------" << endl;
	cout << "  Number of entries (10 years / 120 months): " << skim->GetEntries() << endl;
	cout << "--------------------------------------------------------" << endl;
	cout << " " << endl;
	
	// -------------------------------------------------------------- 	CUT DEFINITION...

	TCut cTrig_HET = "BGO_HET>0.";
	TCut cEne = "BGO_EnergyG_SatCorr_ML_ions>100."; // Min deposited energy [GeV]
	TCut cut00  =  cTrig_HET*cEne;
	TCut cut01  = "(PSD_ChargeY0>0.0 || PSD_ChargeY1>0.0) && (PSD_ChargeX0>0.0 || PSD_ChargeX1>0.0)";  //(Elisabetta 02-12-2021)   
	TCut cut05  = "fabs(BGO_cbgomax[0]-BGO_cbgostk[0])<30.0 && fabs(BGO_cbgomax[1]-BGO_cbgostk[1])<30.0";
	TCut cut06 = "fabs(STKtrack_to_PSD_topY)<400. && fabs(STKtrack_to_PSD_topX)<400."; 

	// STK cut (Elisabetta 22-02-2024)
	TCut cutSTK1st = "((((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2.*STK_chargeY_etaCorr[0] + (TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.*STK_chargeX_etaCorr[0]) / ((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2. + (TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.) ) >1200.)";

	TCut cutSTK2nd = "((((TMath::Sign(1.,STK_chargeY_etaCorr[1])+1.)/2.*STK_chargeY_etaCorr[1] + (TMath::Sign(1.,STK_chargeX_etaCorr[1])+1.)/2.*STK_chargeX_etaCorr[1]) / ((TMath::Sign(1.,STK_chargeY_etaCorr[1])+1.)/2. + (TMath::Sign(1.,STK_chargeX_etaCorr[1])+1.)/2.) ) >1200.)";

	// 06-11-2023 --> Provo ad implementare un taglio sul PSD scegliendo solo gli eventi che hanno una consistenza tra le due viste del PSD!
	// 				  In particolare: impongo che la differenza tra la carica in una vista (singolo hit o media) abbia una differenza minore di 2 con l'altra vista.
	TCut cutPSD = "fabs(( ( (TMath::Sign(1.,PSD_ChargeY0) + 1.)/2. * PSD_ChargeY0 + (TMath::Sign(1.,PSD_ChargeY1) + 1.)/2. * PSD_ChargeY1 ) / ((TMath::Sign(1.,PSD_ChargeY0) + 1.)/2. + (TMath::Sign(1.,PSD_ChargeY1) + 1.)/2.) ) - ( ( (TMath::Sign(1.,PSD_ChargeX0) + 1.)/2. * PSD_ChargeX0 + (TMath::Sign(1.,PSD_ChargeX1) + 1.)/2. * PSD_ChargeX1 ) / ((TMath::Sign(1.,PSD_ChargeX0) + 1.)/2. + (TMath::Sign(1.,PSD_ChargeX1) + 1.)/2.) ) ) < 2.";
	
	TString PSDcharge = "PSD_PathWeighted_Charge";

	/*
	TString Qlow= "(15.232534 + (0.2106836)*BGO_EnergyG_SatCorr_ML_ions + (-0.03401655)*pow(BGO_EnergyG_SatCorr_ML_ions,2) + (-0.000744557)*pow(BGO_EnergyG_SatCorr_ML_ions,3))";
	//"(15.472350 + (0.024609*log10(BGO_EnergyG_SatCorr_ML_ions)))";

	TString Qhigh = "(17.009699 + (-0.1335239)*BGO_EnergyG_SatCorr_ML_ions + (0.06437345)*pow(BGO_EnergyG_SatCorr_ML_ions,2) + (-0.000744557)*pow(BGO_EnergyG_SatCorr_ML_ions,3))";
	//"(16.626725 + (0.090761*log10(BGO_EnergyG_SatCorr_ML_ions)))";
	
	TCut chargeUP = TCut( (PSDcharge+"<("+Qhigh+")").Data() );
	TCut chargeDW = TCut( (PSDcharge+">("+Qlow+")").Data() );
	TCut chargeCUT = chargeUP*chargeDW;
	*/

	TString MPV_DATA = "(14.5821+(1.51344*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.48318*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))+(0.0520877*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2)))";
	TString sigma_DATA =  "(0.20977+(0.0478004*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.00269758*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2)))";

	TString cutLow = "2.0";
    TString cutUp = "2.5";
	TCut chargeUP = (PSDcharge + "<(" + MPV_DATA + "+" + cutUp + "*" + sigma_DATA + ")").Data();
    TCut chargeDW = (PSDcharge + ">(" + MPV_DATA + "-" + cutLow + "*" + sigma_DATA + ")").Data();
    TCut chargeCUT = chargeUP*chargeDW;

	// -------   TAGLI COMPLESSIVO   -------
	TCut ctot = cut00*cut01*cut05*cut06*cutSTK1st*cutPSD*chargeCUT;

	TCut ctoto_0s   = cut00*cut01*cut05*cut06*cutSTK2nd*cutPSD*chargeCUT;
	TCut ctoto_ANDs = cut00*cut01*cut05*cut06*cutSTK1st*cutSTK2nd*cutPSD*chargeCUT;

	// -------------------------------------------------------------- 	FILL HISTOGRAMS...

	TFile *foutp = new TFile("out_eff/eff_DATA_120months_SULFUR_STK2nd_layer_8bin.root","RECREATE");

	TH1F *h1_Unb = new TH1F("h1_Unb", "Nobs(Eo) - STK 2nd layer (S-peak)", noe, Ebin); h1_Unb->Sumw2();
	skim->Draw("BGO_EnergyG_SatCorr_ML_ions>>h1_Unb",ctoto_0s,""); 

	TH1F *h1Nobs_And = new TH1F("h1Nobs_And", "Nobs(Eo) - STK 2nd*1st (S-peak)", noe, Ebin); h1Nobs_And->Sumw2();
	skim->Draw("BGO_EnergyG_SatCorr_ML_ions>>h1Nobs_And",ctoto_ANDs,""); 

	// ------- EFFICIENCY

	TH1F *h1_HET_Eff = new TH1F("h1_HET_Eff", "Eff vs Ebgo (S-peak)", noe, Ebin); h1_HET_Eff->Sumw2();
	h1_HET_Eff->Divide(h1Nobs_And,h1_Unb);

	cout << " -------------------------------------------" << endl;
	cout << "           Done!  THIS IS THE END!          " << endl;
	cout << " -------------------------------------------" << endl;
	/*
	cout << "cut00: " << skim->GetEntries(cut00) << endl;
	cout << "cut01: " << skim->GetEntries(cut00*cut01) << endl;
	cout << "cut05: " << skim->GetEntries(cut00*cut01*cut05) << endl;
	cout << "cut06: " << skim->GetEntries(cut00*cut01*cut05*cut06) << endl;
	cout << "STK2nd: " << skim->GetEntries(cut00*cut01*cut05*cut06*cutSTK2nd) << endl;
	cout << "PSD: " << skim->GetEntries(cut00*cut01*cut05*cut06*cutSTK2nd*cutPSD) << endl;
	cout << "charge: " << skim->GetEntries(ctoto_0s) << endl;
	*/
	// -------------------------------------------------------------- 	DATA SAVING...

	foutp->cd();
	foutp->Write();
	foutp->ls();
	delete foutp; 

}