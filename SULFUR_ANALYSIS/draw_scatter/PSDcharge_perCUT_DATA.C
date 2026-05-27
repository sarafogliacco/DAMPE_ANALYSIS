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
	// --------------------------------------------------------------	
	cout << " " << endl;
	cout << "------------------------------------------------------" << endl;
	cout << "  Number of entries (10 years / 120 months): " << skim->GetEntries() << endl;
	cout << "------------------------------------------------------" << endl;
	cout << " " << endl;

	TFile *foutp = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/ARGON_ANALYSYS/out_root/PSDcharge_perCUT_DATA_120months_SAr_v0.root","RECREATE");
	
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
	TCut cutSTK1200 = "((((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2.*STK_chargeY_etaCorr[0] + (TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.*STK_chargeX_etaCorr[0]) / ((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2. + (TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.) ) >1200.)";

	// 06-11-2023 --> Provo ad implementare un taglio sul PSD scegliendo solo gli eventi che hanno una consistenza tra le due viste del PSD!
	// 				  In particolare: impongo che la differenza tra la carica in una vista (singolo hit o media) abbia una differenza minore di 2 con l'altra vista.
	TCut cutPSD = "fabs(( ( (TMath::Sign(1.,PSD_ChargeY0) + 1.)/2. * PSD_ChargeY0 + (TMath::Sign(1.,PSD_ChargeY1) + 1.)/2. * PSD_ChargeY1 ) / ((TMath::Sign(1.,PSD_ChargeY0) + 1.)/2. + (TMath::Sign(1.,PSD_ChargeY1) + 1.)/2.) ) - ( ( (TMath::Sign(1.,PSD_ChargeX0) + 1.)/2. * PSD_ChargeX0 + (TMath::Sign(1.,PSD_ChargeX1) + 1.)/2. * PSD_ChargeX1 ) / ((TMath::Sign(1.,PSD_ChargeX0) + 1.)/2. + (TMath::Sign(1.,PSD_ChargeX1) + 1.)/2.) ) ) < 2.";
	
        // 23/04/2026 --- SULFUR
	//TString MPV_DATA_S = "(15.4568+(0.56732*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.162207*log10(BGO_EnergyG_SatCorr_ML_ions)**2)+(0.0187945*log10(BGO_EnergyG_SatCorr_ML_ions)**3))";
    //TString sigma_DATA_S = "(0.171732+(0.0556781*log10(BGO_EnergyG_SatCorr_ML_ions))+(0*log10(BGO_EnergyG_SatCorr_ML_ions)**2)+(0*log10(BGO_EnergyG_SatCorr_ML_ions)**3))";
	TString MPV_DATA_S =  "(15.9434+(0.0730006*log10(BGO_EnergyG_SatCorr_ML_ions))+(0.00533945*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))+(-0.000744557*pow(log10(BGO_EnergyG_SatCorr_ML_ions),3)))";
	//"(14.9723+(1.09723*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.342737*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))+(0.0375106*pow(log10(BGO_EnergyG_SatCorr_ML_ions),3)))";
    TString sigma_DATA_S = "(0.355433+(-0.0688415*log10(BGO_EnergyG_SatCorr_ML_ions))+(0.019678*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2)))";
	//"(0.138753+(0.158506*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.0559774*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))+(0.00789775*pow(log10(BGO_EnergyG_SatCorr_ML_ions),3)))";


        // 13/11/2025 --- ARGON
	//TString MPV_DATA_Ar = "(17.9443+(0.0876739*log10(BGO_EnergyG_SatCorr_ML_ions)))";
    //TString sigma_DATA_Ar = "(0.125492+(0.0779934*log10(BGO_EnergyG_SatCorr_ML_ions)))";

	// 27/05/2026 ARGOON 

	TString MPV_DATA_Ar = "(16.762+(1.14769*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.264239*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))+(0.00423957*pow(log10(BGO_EnergyG_SatCorr_ML_ions),4)))";
    TString sigma_DATA_Ar = "(0.242923+(0.0369508*log10(BGO_EnergyG_SatCorr_ML_ions)))";


	// 13/11/2025 ___ CALCIUM
	TString MPV_DATA_Ca = "(19.8698+(0.115209*log10(BGO_EnergyG_SatCorr_ML_ions)))";
    TString sigma_DATA_Ca = "(0.213149+(0.0331733*log10(BGO_EnergyG_SatCorr_ML_ions)))";

/*	// 28/08/2024  --- NEON
	TString MPV_DATA_Ne = "(9.52695+(0.449769*log10(BGO_EnergyG_SatCorrNeMgSi))+(-0.14382*log10(BGO_EnergyG_SatCorrNeMgSi)**2)+(0.0186715*log10(BGO_EnergyG_SatCorrNeMgSi)**3))";
	TString sigma_DATA_Ne = "(0.163224+(0.0132331*log10(BGO_EnergyG_SatCorrNeMgSi))+(0.000357618*log10(BGO_EnergyG_SatCorrNeMgSi)**2))";

	// 28/08/2024  --- MAGNESIUM
	TString MPV_DATA_Mg = "(11.203+(0.848194*log10(BGO_EnergyG_SatCorrNeMgSi))+(-0.274584*log10(BGO_EnergyG_SatCorrNeMgSi)**2)+(0.0325686*log10(BGO_EnergyG_SatCorrNeMgSi)**3))";
	TString sigma_DATA_Mg = "(0.266085+(-0.0435194*log10(BGO_EnergyG_SatCorrNeMgSi))+(0.00923917*log10(BGO_EnergyG_SatCorrNeMgSi)**2))";

	// 28/08/2024  --- SILICON
	TString MPV_DATA_Si = "(12.9701+(1.09124*log10(BGO_EnergyG_SatCorrNeMgSi))+(-0.351935*log10(BGO_EnergyG_SatCorrNeMgSi)**2)+(0.0405486*log10(BGO_EnergyG_SatCorrNeMgSi)**3))";
	TString sigma_DATA_Si = "(0.223611+(0.00107023*log10(BGO_EnergyG_SatCorrNeMgSi))+(0.000461643*log10(BGO_EnergyG_SatCorrNeMgSi)**2))";
*/	
	// --------------------
	
	TString cutLow_S = "2.0";
	TString cutUp_S = "2.5";
	TString PSDcharge = "PSD_PathWeighted_Charge";
	TCut chargeUP_S  = TCut(PSDcharge+"<("+MPV_DATA_S+"+"+cutUp_S +"*"+sigma_DATA_S+")");
	TCut chargeDW_S  = TCut(PSDcharge+">("+MPV_DATA_S+"-"+cutLow_S+"*"+sigma_DATA_S+")");
	TCut chargeCUT_S = chargeUP_S * chargeDW_S;

	TString cutLow_Ar = "1.5";
	TString cutUp_Ar = "2.0";
	TCut chargeUP_Ar  = TCut(PSDcharge+"<("+MPV_DATA_Ar+"+"+cutUp_Ar +"*"+sigma_DATA_Ar+")");
	TCut chargeDW_Ar  = TCut(PSDcharge+">("+MPV_DATA_Ar+"-"+cutLow_Ar+"*"+sigma_DATA_Ar+")");
	TCut chargeCUT_Ar = chargeUP_Ar * chargeDW_Ar;
/*
	TCut chargeUP_Ca  = TCut(PSDcharge+"<("+MPV_DATA_Ca+"+"+cutUp_Ca +"*"+sigma_DATA_Ca+")");
	TCut chargeDW_Ca  = TCut(PSDcharge+">("+MPV_DATA_Ca+"-"+cutLow_Ca+"*"+sigma_DATA_Ca+")");
	TCut chargeCUT_Ca = chargeUP_Ca * chargeDW_Ca;
*/
	// ******************************************************************************

	TCut ctot00 = cut00;
	TCut ctot01 = cut00*cut01;
	//TCut ctot02 = cut00*cut01;
	//TCut ctot03 = cut00*cut01*cut03;
	TCut ctot05 = cut00*cut01*cut05;
	TCut ctot06 = cut00*cut01*cut05*cut06;
	TCut ctotSTK= cut00*cut01*cut05*cut06*cutSTK1200;
	TCut ctotPSD= cut00*cut01*cut05*cut06*cutSTK1200*cutPSD;
	TCut ctotS = cut00*cut01*cut05*cut06*cutSTK1200*cutPSD*chargeCUT_S;
	TCut ctotAr = cut00*cut01*cut05*cut06*cutSTK1200*cutPSD*chargeCUT_Ar;
	//TCut ctotCa = cut00*cut01*cut05*cut06*cutSTK1200*cutPSD*chargeCUT_Ca;

	// -------------------------------------------------------------- 	FILL HISTOGRAMS...

	TH1F *hqpsdPathWeight_all = new TH1F("hqpsdPathWeight_all"," ", 1200, 0., 30.); hqpsdPathWeight_all->GetXaxis()->SetTitle("Q_{PSD}"); hqpsdPathWeight_all->GetYaxis()->SetTitle("events"); 
	skim->Draw("PSD_PathWeighted_Charge>>hqpsdPathWeight_all","","");

	TH1F *hqpsdPathWeight_cut00 = new TH1F("hqpsdPathWeight_cut00"," ", 1200, 0., 30.); hqpsdPathWeight_cut00->GetXaxis()->SetTitle("Q_{PSD}"); hqpsdPathWeight_cut00->GetYaxis()->SetTitle("events");
	skim->Draw("PSD_PathWeighted_Charge>>hqpsdPathWeight_cut00",ctot00,"");

	TH1F *hqpsdPathWeight_cut01 = new TH1F("hqpsdPathWeight_cut01"," ", 1200, 0., 30.); hqpsdPathWeight_cut01->GetXaxis()->SetTitle("Q_{PSD}"); hqpsdPathWeight_cut01->GetYaxis()->SetTitle("events");
	skim->Draw("PSD_PathWeighted_Charge>>hqpsdPathWeight_cut01",ctot01,"");

	TH1F *hqpsdPathWeight_cut05 = new TH1F("hqpsdPathWeight_cut05"," ", 1200, 0., 30.); hqpsdPathWeight_cut05->GetXaxis()->SetTitle("Q_{PSD}"); hqpsdPathWeight_cut05->GetYaxis()->SetTitle("events");
	skim->Draw("PSD_PathWeighted_Charge>>hqpsdPathWeight_cut05",ctot05,"");

	TH1F *hqpsdPathWeight_cut06 = new TH1F("hqpsdPathWeight_cut06"," ", 1200, 0., 30.); hqpsdPathWeight_cut06->GetXaxis()->SetTitle("Q_{PSD}"); hqpsdPathWeight_cut06->GetYaxis()->SetTitle("events");        
	skim->Draw("PSD_PathWeighted_Charge>>hqpsdPathWeight_cut06",ctot06,"");

	TH1F *hqpsdPathWeight_stkcut = new TH1F("hqpsdPathWeight_stkcut"," ", 1200, 0., 30.); hqpsdPathWeight_stkcut->GetXaxis()->SetTitle("Q_{PSD}"); hqpsdPathWeight_stkcut->GetYaxis()->SetTitle("events");
	skim->Draw("PSD_PathWeighted_Charge>>hqpsdPathWeight_stkcut",ctotSTK,"");

	TH1F *hqpsdPathWeight_psdcut = new TH1F("hqpsdPathWeight_psdcut"," ", 1200, 0., 30.); hqpsdPathWeight_psdcut->GetXaxis()->SetTitle("Q_{PSD}"); hqpsdPathWeight_psdcut->GetYaxis()->SetTitle("events");
	skim->Draw("PSD_PathWeighted_Charge>>hqpsdPathWeight_psdcut",ctotPSD,"");
	
	TH1F *hqpsdPathWeight_cutS = new TH1F("hqpsdPathWeight_cutS"," ", 1200, 0., 30.); hqpsdPathWeight_cutS->GetXaxis()->SetTitle("Q_{PSD}"); hqpsdPathWeight_cutS->GetYaxis()->SetTitle("events");
	skim->Draw("PSD_PathWeighted_Charge>>hqpsdPathWeight_cutS",ctotS,"");
	
	TH1F *hqpsdPathWeight_cutAr = new TH1F("hqpsdPathWeight_cutAr"," ", 1200, 0., 30.); hqpsdPathWeight_cutAr->GetXaxis()->SetTitle("Q_{PSD}"); hqpsdPathWeight_cutAr->GetYaxis()->SetTitle("events");
	skim->Draw("PSD_PathWeighted_Charge>>hqpsdPathWeight_cutAr",ctotAr,"");

	//TH1F *hqpsdPathWeight_cutCa = new TH1F("hqpsdPathWeight_cutCa"," ", 1200, 0., 30.); hqpsdPathWeight_cutCa->GetXaxis()->SetTitle("Q_{PSD}"); hqpsdPathWeight_cutCa->GetYaxis()->SetTitle("events");
	//skim->Draw("PSD_PathWeighted_Charge>>hqpsdPathWeight_cutCa",ctotCa,"");
	
	
	cout << " -------------------------------------------" << endl;
	cout << "           Done!  THIS IS THE END!          " << endl;
	cout << " -------------------------------------------" << endl;

	// -------------------------------------------------------------- 	DATA SAVING...

	foutp->cd();
	foutp->Write();  

}
