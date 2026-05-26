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

		
	const int NUM_SET = 4;

	// -------------------------------------------------------------- 	DATA LOADING...
	TChain *skim[NUM_SET];
	for(int i=0; i < NUM_SET; i++){ skim[i] = new TChain("newtree"); }
	// ------- 100 GeV - 1 TeV
	skim[0]->Add("/mnt/c/Users/saraf/Desktop/dampe/MC/MC_SULFUR/skim_MC_allS32_100GeV_1TeV.root");
	// ------- 1 - 10 TeV 
	skim[1]->Add("/mnt/c/Users/saraf/Desktop/dampe/MC/MC_SULFUR/skim_MC_allS32_1TeV_10TeV.root");
	// ------- 10 - 100 TeV 
	skim[2]->Add("/mnt/c/Users/saraf/Desktop/dampe/MC/MC_SULFUR/skim_MC_allS32_10TeV_100TeV.root");
	// ------- 100 - 500 TeV 
	skim[3]->Add("/mnt/c/Users/saraf/Desktop/dampe/MC/MC_SULFUR/skim_MC_allS32_100TeV_500TeV.root");
	// ------- 500 TeV - 1 PeV
	//skim[4]->Add("/mnt/c/Users/saraf/Desktop/dampe/MC/MC_SULFUR/skim_MC_allS32_500TeV_1PeV_EPOSLHC_FTFP_merged.root");

	// Number of entries ...
	cout << " " << endl;
	cout << "-----------------------------------" << endl;
	cout << "Number of entries: " << endl;
	cout << "100GeV - 1TeV  : " << skim[0]->GetEntries() << endl;
	cout << "  1TeV - 10TeV : " << skim[1]->GetEntries() << endl;
	cout << " 10TeV - 100TeV: " << skim[2]->GetEntries() << endl;
	cout << "100TeV - 500TeV:  " <<skim[3]->GetEntries() << endl;
	//cout << "500TeV - 1 PeV :  " <<skim[4]->GetEntries() << endl;
	cout << "-----------------------------------" << endl;
	cout << " " << endl;

	// -------------------------------------------------------------- 	WEIGHTS DEFINITION...

	TCut wS[NUM_SET]; 
	wS[0] = "(1./37980000.)*log(10.)";
	wS[1] = "(1./10778800.)*log(10.)";
	wS[2] = "(1./5110000.)*log(10.)";
	wS[3] = "(1./1669300.)*log(5.)";
	//wS[4] = "(1./3050580.)*log(2.)";

	TCut wEnergy = "(MC_EnergyT)**(-1.7)";

	TCut wSN[NUM_SET];
	for (int i=0; i<NUM_SET; i++) { wSN[i] = wS[i]*wEnergy; };

	// -------------------------------------------------------------- 	CUTS DEFINITION...

	TCut Emin = "MC_EnergyT>100.";

	TCut cTrig_HET = "BGO_HET>0.";
	TCut cEne = "BGO_EnergyG_QuenchSatCorr_ML_ions>100."; // Min deposited energy [GeV]
	TCut cut00  = cTrig_HET*cEne;
	TCut cut01  = "(PSD_ChargeY0>0.0 || PSD_ChargeY1>0.0) && (PSD_ChargeX0>0.0 || PSD_ChargeX1>0.0)";    //(Elisabetta 02-12-2021)   
	TCut cut05  = "fabs(BGO_cbgomax[0]-BGO_cbgostk[0])<30.0 && fabs(BGO_cbgomax[1]-BGO_cbgostk[1])<30.0";
	TCut cut06 = "fabs(STKtrack_to_PSD_topY)<400. && fabs(STKtrack_to_PSD_topX)<400."; 

	// STK cut (Elisabetta 22-02-2024)
	TCut cutSTK1st = "((((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2.*STK_chargeY_etaCorr[0] + (TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.*STK_chargeX_etaCorr[0]) / ((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2. + (TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.) ) >1200.)";

	TCut cutSTK2nd = "((((TMath::Sign(1.,STK_chargeY_etaCorr[1])+1.)/2.*STK_chargeY_etaCorr[1] + (TMath::Sign(1.,STK_chargeX_etaCorr[1])+1.)/2.*STK_chargeX_etaCorr[1]) / ((TMath::Sign(1.,STK_chargeY_etaCorr[1])+1.)/2. + (TMath::Sign(1.,STK_chargeX_etaCorr[1])+1.)/2.) ) >1200.)";

	// 06-11-2023 --> Provo ad implementare un taglio sul PSD scegliendo solo gli eventi che hanno una consistenza tra le due viste del PSD!
	// 				  In particolare: impongo che la differenza tra la carica in una vista (singolo hit o media) abbia una differenza minore di 2 con l'altra vista.
	TCut cutPSD = "fabs(( ( (TMath::Sign(1.,PSD_ChargeY0) + 1.)/2. * PSD_ChargeY0 + (TMath::Sign(1.,PSD_ChargeY1) + 1.)/2. * PSD_ChargeY1 ) / ((TMath::Sign(1.,PSD_ChargeY0) + 1.)/2. + (TMath::Sign(1.,PSD_ChargeY1) + 1.)/2.) ) - ( ( (TMath::Sign(1.,PSD_ChargeX0) + 1.)/2. * PSD_ChargeX0 + (TMath::Sign(1.,PSD_ChargeX1) + 1.)/2. * PSD_ChargeX1 ) / ((TMath::Sign(1.,PSD_ChargeX0) + 1.)/2. + (TMath::Sign(1.,PSD_ChargeX1) + 1.)/2.) ) ) < 2.";
	
	// SULFUR
	TString MPV_MC = "(14.976+(1.78884*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(-0.528537*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))+(0.0554816*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3)))";
	TString sigma_MC = "(0.0454068+(0.0415077*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(-0.00228883*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2)))";
	

	TString MPV_DATA = "(14.5821+(1.51344*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.48318*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))+(0.0520877*pow(log10(BGO_EnergyG_SatCorr_ML_ions),3)))";
	TString sigma_DATA = "(0.20977+(0.0478004*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.00269758*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2)))";
	

	TString MPV_MC_CORR = "(15.0971+(1.00316*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(-0.323067*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))+(0.0374537*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3)))";
	TString sigma_MC_CORR = "(0.305322+(-0.0188628*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(0.00897858*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2)))";
	
	TString PSDcharge = "PSD_PathWeighted_Charge";
	TString QcorrMC = "(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")";

	/*
	TString Qhigh= "(17.009699 + (-0.1335239)*BGO_EnergyG_SatCorr_ML_ions + (0.06437345)*pow(BGO_EnergyG_SatCorr_ML_ions,2) + (-0.000744557)*pow(BGO_EnergyG_SatCorr_ML_ions,3))";
	TString Qlow = "(15.232534 + (0.2106836)*BGO_EnergyG_SatCorr_ML_ions + (-0.03401655)*pow(BGO_EnergyG_SatCorr_ML_ions,2) + (-0.000744557)*pow(BGO_EnergyG_SatCorr_ML_ions,3))";
	TCut chargeUP = TCut( (PSDcharge+"<("+Qhigh+")").Data() );
	TCut chargeDW = TCut( (PSDcharge+">("+Qlow+")").Data() );
	TCut chargeCUT = chargeUP*chargeDW;
	*/

	TString cutLow = "2.0";
    TString cutUp = "2.5";
	TCut chargeUP = (PSDcharge + "<(" + MPV_DATA + "+" + cutUp + "*" + sigma_DATA + ")").Data();
    TCut chargeDW = (PSDcharge + ">(" + MPV_DATA + "-" + cutLow + "*" + sigma_DATA + ")").Data();
    TCut chargeCUT = chargeUP*chargeDW;
	
	// -------   TAGLIO COMPLESSIVO   -------
	TCut ctot = cut00*cut01*cut05*cut06*cutSTK1st*cutPSD*chargeCUT;

	TCut ctoto_0s   = cut00*cut01*cut05*cut06*cutSTK2nd*cutPSD*chargeCUT;
	TCut ctoto_ANDs = cut00*cut01*cut05*cut06*cutSTK1st*cutSTK2nd*cutPSD*chargeCUT;
	
	// -------------------------------------------------------------- 	HERE, THE MAGIC! (FILLING HISTOGRAMS)

	TFile *foutp = new TFile("out_eff/eff_MC_FTFP_SULFUR_STK2nd_layer_8bin.root","RECREATE");

	TH1F *h1_Unb = new TH1F("h1_Unb", "Nobs(Eo) - STK 2nd layer (SULFUR-peak)", noe, Ebin); h1_Unb->Sumw2();
	for (int i=0; i< NUM_SET; i++) { skim[i]->Draw("BGO_EnergyG_QuenchSatCorr_ML_ions>>+h1_Unb",ctoto_0s*wSN[i],""); }

	TH1F *h1Nobs_And = new TH1F("h1Nobs_And", "Nobs(Eo) - STK 2nd*1st (SULFUR-peak)", noe, Ebin); h1Nobs_And->Sumw2();
	for (int i=0; i< NUM_SET; i++) { skim[i]->Draw("BGO_EnergyG_QuenchSatCorr_ML_ions>>+h1Nobs_And",ctoto_ANDs*wSN[i],"");  }


	TH1F *h1_HET_Eff = new TH1F("h1_HET_Eff", "Eff vs Ebgo (SULFUR-peak)", noe, Ebin); h1_HET_Eff->Sumw2();
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