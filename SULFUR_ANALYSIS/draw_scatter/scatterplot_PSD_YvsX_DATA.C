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


    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadLeftMargin(0.09);
    gStyle->SetPadRightMargin(0.11);
    gStyle->SetPadTopMargin(0.03);
    gStyle->SetPadBottomMargin(0.09);

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

	// ------------------------------------------------------------------------------

	TFile *fout = new TFile("out_root/scatterplot_DATA_PSD_120months_YvsX.root","RECREATE");

	// ------------------------------------------------------------------------------

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
	
	// -------   TAGLIO COMPLESSIVO   -------
	TCut ctot_stk = cut00*cut01*cut05*cut06*cutSTK1200;
	TCut ctot_TOT = cut00*cut01*cut05*cut06*cutSTK1200*cutPSD;

	// ---------------------

	TString meanX = "(((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.*PSD_ChargeX0+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.*PSD_ChargeX1)/((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.))";
	TString meanY = "(((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_ChargeY0+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_ChargeY1)/((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.))";
	
	TString avX = "(((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.*PSD_ChargeX0*PSD_psdX_pathlength[0]+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.*PSD_ChargeX1*PSD_psdX_pathlength[1])/((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.*PSD_psdX_pathlength[0]+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.*PSD_psdX_pathlength[1]))";
	TString avY = "(((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_ChargeY0*PSD_psdY_pathlength[0]+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_ChargeY1*PSD_psdY_pathlength[1])/((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_psdY_pathlength[0]+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_psdY_pathlength[1]))";

	TString PSDcharge = "PSD_PathWeighted_Charge";

	// --------------------

	TCanvas *c1 = new TCanvas("c1","c1",800,780);
	TH2F *h2psd_stk = new TH2F("h2psd_stk"," ", 600, 0., 30., 600, 0., 30.); 
	h2psd_stk->GetXaxis()->SetTitle("Q_{X}"); h2psd_stk->GetYaxis()->SetTitle("Q_{Y}"); 
	h2psd_stk->GetXaxis()->SetTitleOffset(1.2); h2psd_stk->GetYaxis()->SetTitleOffset(1.2); h2psd_stk->GetZaxis()->SetTitleOffset(1.2); 
	skim->Draw(meanY+":"+meanX+">>h2psd_stk",ctot_stk,"colz");
	h2psd_stk->SetStats(0);  c1->SetLogz(); c1->Update();


	TCanvas *c2 = new TCanvas("c2","c2",800,780);
	TH2F *h2psd_TOT = new TH2F("h2psd_TOT"," ", 600, 0., 30., 600, 0., 30.); 
	h2psd_TOT->GetXaxis()->SetTitle("Q_{X}"); h2psd_TOT->GetYaxis()->SetTitle("Q_{Y}"); 
	h2psd_TOT->GetXaxis()->SetTitleOffset(1.2); h2psd_TOT->GetYaxis()->SetTitleOffset(1.2); h2psd_TOT->GetZaxis()->SetTitleOffset(1.2); 
	skim->Draw(meanY+":"+meanX+">>h2psd_TOT",ctot_TOT,"colz");
	h2psd_TOT->SetStats(0);  c2->SetLogz(); c2->Update();

	// ---------------------------------------------------------------------------------------------------

	fout->cd();
	fout->Write();
	fout->Close();	


}
