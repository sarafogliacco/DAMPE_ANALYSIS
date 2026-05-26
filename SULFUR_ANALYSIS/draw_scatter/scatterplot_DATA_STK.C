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
	cout << " " << endl;
	cout << "------------------------------------------------------" << endl;
	cout << "  Number of entries (10 years / 120 months): " << skim->GetEntries() << endl;
	cout << "------------------------------------------------------" << endl;
	cout << " " << endl;	

	TFile *foutp = new TFile("out_root/scatterplot_DATA_STK_120months_etaCorr.root","RECREATE");
	
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
	TCut cutSTK1200 = "( (((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2.*STK_chargeY_etaCorr[0] + (TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.*STK_chargeX_etaCorr[0]) / ((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2. + (TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.) ) >1200.)";

	// 06-11-2023 --> Provo ad implementare un taglio sul PSD scegliendo solo gli eventi che hanno una consistenza tra le due viste del PSD!
	// 				  In particolare: impongo che la differenza tra la carica in una vista (singolo hit o media) abbia una differenza minore di 2 con l'altra vista.
	TCut cutPSD = "fabs(( ( (TMath::Sign(1.,PSD_ChargeY0) + 1.)/2. * PSD_ChargeY0 + (TMath::Sign(1.,PSD_ChargeY1) + 1.)/2. * PSD_ChargeY1 ) / ((TMath::Sign(1.,PSD_ChargeY0) + 1.)/2. + (TMath::Sign(1.,PSD_ChargeY1) + 1.)/2.) ) - ( ( (TMath::Sign(1.,PSD_ChargeX0) + 1.)/2. * PSD_ChargeX0 + (TMath::Sign(1.,PSD_ChargeX1) + 1.)/2. * PSD_ChargeX1 ) / ((TMath::Sign(1.,PSD_ChargeX0) + 1.)/2. + (TMath::Sign(1.,PSD_ChargeX1) + 1.)/2.) ) ) < 2.";
	
	// -------   TAGLIO COMPLESSIVO   -------

	TCut precuts = cut00*cut01*cut05*cut06*cutPSD;
	TCut ctot    = cut00*cut01*cut05*cut06*cutSTK1200*cutPSD;

	// -------------------------------------------------------------- 	FILL HISTOGRAMS...

	TH2F *h2stk_precuts = new TH2F("h2stk_precuts","",300,0.,30.,600,0.,15000.); 
	h2stk_precuts->GetXaxis()->SetTitle("PSD charge"); h2stk_precuts->GetYaxis()->SetTitle("STK first layer average (ADC counts)");
	h2stk_precuts->GetYaxis()->SetTitleOffset(1.3); h2stk_precuts->GetXaxis()->SetTitleOffset(1.1); h2stk_precuts->GetZaxis()->SetLabelSize(0.03); 
	skim->Draw("((((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2.*STK_chargeY_etaCorr[0]+(TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.*STK_chargeX_etaCorr[0])/((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2.+(TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.))):PSD_PathWeighted_Charge>>h2stk_precuts",precuts,"");
	cout << " ...h2stk_precuts... " << endl;

	TH2F *h2stk = new TH2F("h2stk","",300,0.,30.,600,0.,15000.); 
	h2stk->GetXaxis()->SetTitle("PSD charge"); h2stk->GetYaxis()->SetTitle("STK first layer average (ADC counts)");
	h2stk->GetYaxis()->SetTitleOffset(1.3); h2stk->GetXaxis()->SetTitleOffset(1.1); h2stk->GetZaxis()->SetLabelSize(0.03); 
	skim->Draw("((((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2.*STK_chargeY_etaCorr[0]+(TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.*STK_chargeX_etaCorr[0])/((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2.+(TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.))):PSD_PathWeighted_Charge>>h2stk",ctot,"");
	cout << " ...h2stk... " << endl;


	cout << " -------------------------------------------" << endl;
	cout << "           Done!  THIS IS THE END!          " << endl;
	cout << " -------------------------------------------" << endl;

	// -------------------------------------------------------------- 	DATA SAVING...

	foutp->cd();
	foutp->Write();  

}
