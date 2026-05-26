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
	Int_t nbd = 8;
	Int_t ndec = 4;
	Double_t arg = 1./double(nbd);
	int noe = nbd*ndec;
	//cout << "Number of energy bin " << noe << endl;
	Double_t e0 = 100.;
	std::vector<Double_t> Ebin(noe + 1);
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
	cout << "------------------------------------------------------" << endl;
	cout << "  Number of entries (10 years / 120 months): " << skim->GetEntries() << endl;
	cout << "------------------------------------------------------" << endl;
	cout << " " << endl;

	TFile *foutp = new TFile("out_eff/CHARGE_eff_DATA_120months_SULFUR_SATcorr_STKcut1200_PSDcutXY_8bin.root","RECREATE");
	
	// -------------------------------------------------------------- 	CUT DEFINITION...

	TCut cTrig_HET = "BGO_HET>0.";
    TCut cEne = "BGO_EnergyG_SatCorr_ML_ions>100."; // Min deposited energy [GeV]
    TCut cut00  =  cTrig_HET*cEne;
    TCut cut01  = "(PSD_ChargeY0>0.0 || PSD_ChargeY1>0.0) && (PSD_ChargeX0>0.0 || PSD_ChargeX1>0.0)";  //(Elisabetta 02-12-2021)   
	//TCut cut02  = "BGO_EneLay[0]+BGO_EneLay[1]<BGO_EneLay[2]+BGO_EneLay[3]";  // sciame corretto nel bgo SECONDO LA CORREZIONE DI CHUAN
	//TCut cut03  = "(PSD_EnergyX_corr[0]+PSD_EnergyX_corr[1])*10./(PSD_psdX_pathlength[0]+PSD_psdX_pathlength[1])>0.0 && (PSD_EnergyY_corr[0]+PSD_EnergyY_corr[1])*10/(PSD_psdY_pathlength[0]+PSD_psdY_pathlength[1])>0.0";  // un segnale su entrambe i due piani di PSD + traccia scelta [QUESTO CUT È INUTILE QUANDO SI HA A CHE FARE CON I FILE SKIMMATI -> LA RICHIESTA DI SEGNALE SU I DUE PSD + LA TRACCIA È AUTOMATICA GIÀ DAL CUT 01 IN QUESTO CASO]
	TCut cut05  = "fabs(BGO_cbgomax[0]-BGO_cbgostk[0])<30.0 && fabs(BGO_cbgomax[1]-BGO_cbgostk[1])<30.0";
    TCut cut06 = "fabs(STKtrack_to_PSD_topY)< 400. && fabs(STKtrack_to_PSD_topX) < 400";

	// STK cut (Elisabetta 22-02-2024)
    TCut cutSTK1200 = "((((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2.*STK_chargeY_etaCorr[0]+(TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.*STK_chargeX_etaCorr[0])/((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2.+(TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.))>1200.)";

	// 06-11-2023 --> Provo ad implementare un taglio sul PSD scegliendo solo gli eventi che hanno una consistenza tra le due viste del PSD!
	// 				  In particolare: impongo che la differenza tra la carica in una vista (singolo hit o media) abbia una differenza minore di 2 con l'altra vista.
    TCut cutPSD = "fabs((((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_ChargeY0+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_ChargeY1)/((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.))-(((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.*PSD_ChargeX0+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.*PSD_ChargeX1)/((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.)))<2.";
	
	TString PSDcharge = "PSD_PathWeighted_Charge";
	
	TString avX = "(((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.*PSD_ChargeX0*PSD_psdX_pathlength[0]+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.*PSD_ChargeX1*PSD_psdX_pathlength[1])/((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.*PSD_psdX_pathlength[0]+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.*PSD_psdX_pathlength[1]))";
	TString avY = "(((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_ChargeY0*PSD_psdY_pathlength[0]+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_ChargeY1*PSD_psdY_pathlength[1])/((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_psdY_pathlength[0]+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_psdY_pathlength[1]))";

	// 29/08/2024 - PSD Y
	TString cutLow_PSDY = "2.0";
	TString cutUp_PSDY  = "2.5";
	TString MPV_DATA_PSDY = "(14.9817+(1.0822*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.37147*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))+(0.042988*pow(log10(BGO_EnergyG_SatCorr_ML_ions),3)))";
	TString sigma_DATA_PSDY = "(0.55017+(-0.134185*log10(BGO_EnergyG_SatCorr_ML_ions))+(0.0260233*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2)))";


	// Correct TCut construction using .Data()
	TCut chargeUP_PSDY = TCut( (avY+"<("+MPV_DATA_PSDY+"+"+cutUp_PSDY +"*"+sigma_DATA_PSDY+")").Data() );
	TCut chargeDW_PSDY = TCut( (avY+">("+MPV_DATA_PSDY+"-"+cutLow_PSDY+"*"+sigma_DATA_PSDY+")").Data() );
	TCut chargeCUT_PSDY = chargeUP_PSDY*chargeDW_PSDY;

	// 29/08/2024 - PSD X
	TString cutLow_PSDX = "2.0";
	TString cutUp_PSDX  = "2.5";
	TString MPV_DATA_PSDX = "(14.9041+(1.03958*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.275602*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))+(0.0244722*pow(log10(BGO_EnergyG_SatCorr_ML_ions),3)))";
	TString sigma_DATA_PSDX = "(-0.0661661+(0.308426*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.0437239*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2)))";
	TCut chargeUP_PSDX = TCut( (avX+"<("+MPV_DATA_PSDX+"+"+cutUp_PSDX +"*"+sigma_DATA_PSDX+")").Data() );
	TCut chargeDW_PSDX = TCut( (avX+">("+MPV_DATA_PSDX+"-"+cutLow_PSDX+"*"+sigma_DATA_PSDX+")").Data() );
	TCut chargeCUT_PSDX = chargeUP_PSDX*chargeDW_PSDX;

	// PSD charge
	TString cutLow = "2.0";
	TString cutUp = "2.5";
	TString MPV_DATA = "(14.5821+(1.51344*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2)))+(-0.48318*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2)))+(0.0520877*pow(log10(BGO_EnergyG_SatCorr_ML_ions),3))))";
	TString sigma_DATA = "(0.20977+(0.0478004*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2)))+(-0.00269758*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))))";

	TCut chargeUP = TCut( (PSDcharge+"<("+MPV_DATA+"+"+cutUp +"*"+sigma_DATA+")").Data() );
	TCut chargeDW = TCut( (PSDcharge+">("+MPV_DATA+"-"+cutLow+"*"+sigma_DATA+")").Data() );
	TCut chargeCUT = chargeUP*chargeDW;

	// -------   TAGLIO COMPLESSIVO   -------
	TCut ctot = cut00*cut01*cut05*cut06*cutSTK1200*cutPSD;
		
	TCut ctoto_0s   = ctot*chargeCUT_PSDY*chargeCUT_PSDX;
	TCut ctoto_ANDs = ctot*chargeCUT_PSDY*chargeCUT_PSDX*chargeCUT;

	// -------------------------------------------------------------- 	FILL HISTOGRAMS...
	TH1F *hCharge_XY = new TH1F("hCharge_XY", "Nobs(Eo) - PSD YX (S-peak)", noe, Ebin.data()); hCharge_XY->Sumw2();
	skim->Draw("BGO_EnergyG_SatCorr_ML_ions>>hCharge_XY",ctoto_0s,""); cout<<" ...hCharge_XY..."<<endl;

	TH1F *hCharge_Glob = new TH1F("hCharge_Glob", "Nobs(Eo) - PSD YX * PSD charge (S-peak)", noe, Ebin.data()); hCharge_Glob->Sumw2();
	skim->Draw("BGO_EnergyG_SatCorr_ML_ions>>hCharge_Glob",ctoto_ANDs,""); cout<<" ...hCharge_Glob..."<<endl;

	// ------- EFFICIENCY

	TH1F *h1_PSD_Eff = new TH1F("h1_PSD_Eff", "CHARGE-Eff vs Ebgo (S-peak)", noe, Ebin.data()); h1_PSD_Eff->Sumw2();
	h1_PSD_Eff->Divide(hCharge_Glob,hCharge_XY); cout<<" ...h1_PSD_Eff..."<<endl;
	

	cout << skim->GetEntries() << endl;

	cout << skim->GetEntries(cut00) << endl;
	cout << skim->GetEntries(cut00*cut01) << endl;
	cout << skim->GetEntries(cut00*cut01*cut05) << endl;

	// -------------------------------------------------------------- 	DATA SAVING...
	foutp->cd();
	foutp->Write();
	foutp->ls();
	delete foutp;

}
