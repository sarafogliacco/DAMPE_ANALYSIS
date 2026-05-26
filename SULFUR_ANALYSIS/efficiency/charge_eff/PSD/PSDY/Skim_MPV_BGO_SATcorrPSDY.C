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
	//-----------------------

	cout << " " << endl;
	cout << "------------------------------------------------------" << endl;
	cout << "  Number of entries (10 years / 120 months): " << skim->GetEntries() << endl;
	cout << "------------------------------------------------------" << endl;
	cout << " " << endl;

	// -------------------------------------------------------------- 	HISTOGRAMS DEFINITION...
	
	TH1F *h01=new TH1F("h01", "100 GeV < E_{BGO} < 178 GeV",350,15,22); h01->GetXaxis()->SetTitle("PSD Y charge"); h01->GetYaxis()->SetTitle("events"); h01->SetLineColor(kBlack); h01->SetMarkerColor(kBlack); h01->Sumw2();
	TH1F *h02=new TH1F("h02", "178 GeV < E_{BGO} < 316 GeV",350,15,22); h02->GetXaxis()->SetTitle("PSD Y charge"); h02->GetYaxis()->SetTitle("events"); h02->SetLineColor(kBlack); h02->SetMarkerColor(kBlack); h02->Sumw2();
	TH1F *h03=new TH1F("h03", "316 GeV < E_{BGO} < 562 GeV",350,15,22); h03->GetXaxis()->SetTitle("PSD Y charge"); h03->GetYaxis()->SetTitle("events"); h03->SetLineColor(kBlack); h03->SetMarkerColor(kBlack); h03->Sumw2();
	TH1F *h04=new TH1F("h04", "562 GeV < E_{BGO} < 1.0 TeV",350,15,22); h04->GetXaxis()->SetTitle("PSD Y charge"); h04->GetYaxis()->SetTitle("events"); h04->SetLineColor(kBlack); h04->SetMarkerColor(kBlack); h04->Sumw2();
	TH1F *h05=new TH1F("h05", "1.0 TeV < E_{BGO} < 1.8 TeV", 200,15,22); h05->GetXaxis()->SetTitle("PSD Y charge"); h05->GetYaxis()->SetTitle("events"); h05->SetLineColor(kBlack); h05->SetMarkerColor(kBlack); h05->Sumw2();
	TH1F *h06=new TH1F("h06", "1.8 TeV < E_{BGO} < 3.2 TeV", 200,15,22); h06->GetXaxis()->SetTitle("PSD Y charge"); h06->GetYaxis()->SetTitle("events"); h06->SetLineColor(kBlack); h06->SetMarkerColor(kBlack); h06->Sumw2();
	TH1F *h07=new TH1F("h07", "3.2 TeV < E_{BGO} < 10.0 TeV",200,15,22); h07->GetXaxis()->SetTitle("PSD Y charge"); h07->GetYaxis()->SetTitle("events"); h07->SetLineColor(kBlack); h07->SetMarkerColor(kBlack); h07->Sumw2();
    TH1F *h08=new TH1F("h08","10.0 TeV < E_{BGO} < 100.0 TeV",80,15,22);h08->GetXaxis()->SetTitle("PSD Y charge"); h08->GetYaxis()->SetTitle("events"); h08->SetLineColor(kBlack); h08->SetMarkerColor(kBlack); h08->Sumw2();
	
	// -------------------------------------------------------------- 	CUT DEFINITION...

	TCut cTrig_HET = "BGO_HET>0."; // trigger cut 
	TCut cEne = "BGO_EnergyG_SatCorr_ML_ions>100."; // Min deposited energy [GeV]
	TCut cut00  =  cTrig_HET*cEne; 
	TCut cut01  = "(PSD_ChargeY0>0.0 || PSD_ChargeY1>0.0) && (PSD_ChargeX0>0.0 || PSD_ChargeX1>0.0)";  //(Elisabetta 02-12-2021)   
	//TCut cut02  = "BGO_EneLay[0]+BGO_EneLay[1]<BGO_EneLay[2]+BGO_EneLay[3]";  // sciame corretto nel bgo SECONDO LA CORREZIONE DI CHUAN
	//TCut cut03  = "(PSD_EnergyX_corr[0]+PSD_EnergyX_corr[1])*10./(PSD_psdX_pathlength[0]+PSD_psdX_pathlength[1])>0.0 && (PSD_EnergyY_corr[0]+PSD_EnergyY_corr[1])*10/(PSD_psdY_pathlength[0]+PSD_psdY_pathlength[1])>0.0";  // un segnale su entrambe i due piani di PSD + traccia scelta [QUESTO CUT È INUTILE QUANDO SI HA A CHE FARE CON I FILE SKIMMATI -> LA RICHIESTA DI SEGNALE SU I DUE PSD + LA TRACCIA È AUTOMATICA GIÀ DAL CUT 01 IN QUESTO CASO]
	
	TCut cut05  = "fabs(BGO_cbgomax[0]-BGO_cbgostk[0])<30.0 && fabs(BGO_cbgomax[1]-BGO_cbgostk[1])<30.0"; // bgo fiducal volume

	TCut cut06 = "fabs(STKtrack_to_PSD_topY)< 400. && fabs(STKtrack_to_PSD_topX) < 400"; // consistency track 

	// STK cut (Elisabetta 22-02-2024) to remove lighter elements 
	TCut cutSTK1200 = "( (((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2.*STK_chargeY_etaCorr[0] + (TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.*STK_chargeX_etaCorr[0]) / ((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2. + (TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.) ) >1200.)";

	// 05-11-2025 --> Provo ad implementare un taglio sul PSD scegliendo solo gli eventi che hanno una consistenza tra le due viste del PSD!
	// 				  In particolare: impongo che la differenza tra la carica in una vista (singolo hit o media) abbia una differenza minore di 2 con l'altra vista.
	TCut cutPSD = "fabs(( ( (TMath::Sign(1.,PSD_ChargeY0) + 1.)/2. * PSD_ChargeY0 + (TMath::Sign(1.,PSD_ChargeY1) + 1.)/2. * PSD_ChargeY1 ) / ((TMath::Sign(1.,PSD_ChargeY0) + 1.)/2. + (TMath::Sign(1.,PSD_ChargeY1) + 1.)/2.) ) - ( ( (TMath::Sign(1.,PSD_ChargeX0) + 1.)/2. * PSD_ChargeX0 + (TMath::Sign(1.,PSD_ChargeX1) + 1.)/2. * PSD_ChargeX1 ) / ((TMath::Sign(1.,PSD_ChargeX0) + 1.)/2. + (TMath::Sign(1.,PSD_ChargeX1) + 1.)/2.) ) ) < 2.";
	
	
	// -------   TAGLIO COMPLESSIVO   -------
	TCut ctot = cut00*cut01*cut05*cut06*cutSTK1200*cutPSD;

	TCut bgo01 = "(BGO_EnergyG_SatCorr_ML_ions > 100.) && (BGO_EnergyG_SatCorr_ML_ions < 177.828)";
    TCut bgo02 = "(BGO_EnergyG_SatCorr_ML_ions > 177.828) && (BGO_EnergyG_SatCorr_ML_ions < 316.228)";
    TCut bgo03 = "(BGO_EnergyG_SatCorr_ML_ions > 316.228) && (BGO_EnergyG_SatCorr_ML_ions < 562.341)";
    TCut bgo04 = "(BGO_EnergyG_SatCorr_ML_ions > 562.341) && (BGO_EnergyG_SatCorr_ML_ions < 1000.0)";
    TCut bgo05 = "(BGO_EnergyG_SatCorr_ML_ions > 1000.0) && (BGO_EnergyG_SatCorr_ML_ions < 1778.28)";
    TCut bgo06 = "(BGO_EnergyG_SatCorr_ML_ions > 1778.28) && (BGO_EnergyG_SatCorr_ML_ions < 3162.28)";
    TCut bgo07 = "(BGO_EnergyG_SatCorr_ML_ions > 3162.28) && (BGO_EnergyG_SatCorr_ML_ions < 10000.0)";
    TCut bgo08 = "(BGO_EnergyG_SatCorr_ML_ions > 10000.0) && (BGO_EnergyG_SatCorr_ML_ions < 1000000.0)";

	TString avY = "(((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_ChargeY0*PSD_psdY_pathlength[0]+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_ChargeY1*PSD_psdY_pathlength[1])/((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_psdY_pathlength[0]+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_psdY_pathlength[1]))";


	// -------------------------------------------------------------- 	FILL HISTOGRAMS...
	
	TCanvas *c0=new TCanvas("c0","BGO-PSDglob",900,600); c0->Divide(2,4);
	c0_1->cd();  c0_1->SetTicks();  c0_2->cd();  c0_2->SetTicks();  c0_3->cd();  c0_3->SetTicks();  
	c0_4->cd();  c0_4->SetTicks();  c0_5->cd();  c0_5->SetTicks();  c0_6->cd();  c0_6->SetTicks();
	c0_7->cd();  c0_7->SetTicks();  c0_8->cd();  c0_8->SetTicks();  

	c0->cd(1);  skim->Draw(avY+">>+h01",ctot*bgo01,""); cout<<"... h01 ..."<<endl;
	c0->cd(2);  skim->Draw(avY+">>+h02",ctot*bgo02,""); cout<<"... h02 ..."<<endl;
	c0->cd(3);  skim->Draw(avY+">>+h03",ctot*bgo03,""); cout<<"... h03 ..."<<endl;
	c0->cd(4);  skim->Draw(avY+">>+h04",ctot*bgo04,""); cout<<"... h04 ..."<<endl;
	c0->cd(5);  skim->Draw(avY+">>+h05",ctot*bgo05,""); cout<<"... h05 ..."<<endl;
	c0->cd(6);  skim->Draw(avY+">>+h06",ctot*bgo06,""); cout<<"... h06 ..."<<endl;
	c0->cd(7);  skim->Draw(avY+">>+h07",ctot*bgo07,""); cout<<"... h07 ..."<<endl;
	c0->cd(8);  skim->Draw(avY+">>+h08",ctot*bgo08,""); cout<<"... h08 ..."<<endl;


	// -------------------------------------------------------------- 	DATA SAVING...
	TFile *f = TFile::Open("Skim_120months_SATcorr_STKcut1200_PSDY.root", "RECREATE");
	f->cd();
	h01->Write(); // save the histogram
	h02->Write(); // save the histogram
	h03->Write(); // save the histogram
	h04->Write(); // save the histogram
	h05->Write(); // save the histogram
	h06->Write(); // save the histogram
	h07->Write(); // save the histogram
	h08->Write(); // save the histogram

	c0->Write();
	f->ls();      // show the contents of the ROOT file
	delete f;     // close the ROOT file   

}
