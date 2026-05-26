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
	//skim[4]->Add("/mnt/c/Users/saraf/Desktop/dampe/MC/MC_SULFUR/skim_MC_allS32_500TeV_1PeV.root");

	// Number of entries ...
	cout << " " << endl;
	cout << "-----------------------------------" << endl;
	cout << "Number of entries: " << endl;
	cout << "100GeV - 1TeV  : " << skim[0]->GetEntries() << endl;
	cout << "  1TeV - 10TeV : " << skim[1]->GetEntries() << endl;
	cout << " 10TeV - 100TeV: " << skim[2]->GetEntries() << endl;
	cout << "100TeV - 500TeV:  " <<skim[3]->GetEntries() << endl;
	cout << "-----------------------------------" << endl;
	cout << " " << endl;

	TCut wS[NUM_SET]; 
	wS[0] = "(1./37980000.)*log(10.)";
	wS[1] = "(1./10778800.)*log(10.)";
	wS[2] = "(1./5110000.)*log(10.)";
	wS[3] = "(1./1669300.)*log(5.)";
	//wSi[4] = "(1./3050580.)*log(2.)";

	TCut wEnergy = "(MC_EnergyT)**(-1.7)";

	TCut wSN[NUM_SET];
	for (int i=0; i<NUM_SET; i++) { wSN[i] = wS[i]*wEnergy; };

	// -------------------------------------------------------------- 	CUTS DEFINITION...

	TCut cTrig_HEP = "BGO_HET>0.";
	TCut cEne = "BGO_EnergyG_QuenchSatCorr_ML_ions>100."; // Min deposited energy [GeV]
	TCut cut00  = cTrig_HEP*cEne;
	TCut cut01  = "(PSD_ChargeY0>0.0 || PSD_ChargeY1>0.0) && (PSD_ChargeX0>0.0 || PSD_ChargeX1>0.0)";    //(Elisabetta 02-12-2021)   
	TCut cut02  = "(BGO_EneLay[0]+BGO_EneLay[1]) < (BGO_EneLay[2]+BGO_EneLay[3])";         // sciame corretto nel bgo SECONDO LA CORREZIONE DI CHUAN
	TCut cut03  = "(PSD_EnergyX_corr[0]+PSD_EnergyX_corr[1])*10./(PSD_psdX_pathlength[0]+PSD_psdX_pathlength[1])>0.0 && (PSD_EnergyY_corr[0]+PSD_EnergyY_corr[1])*10/(PSD_psdY_pathlength[0]+PSD_psdY_pathlength[1])>0.0";  // un segnale su entrambe i due piani di PSD + traccia scelta [QUESTO CUT È INUTILE QUANDO SI HA A CHE FARE CON I FILE SKIMMATI]
	TCut cut05  = "fabs(BGO_cbgomax[0]-BGO_cbgostk[0])<30.0 && fabs(BGO_cbgomax[1]-BGO_cbgostk[1])<30.0";

	// STK cut (Elisabetta 22-02-2024)
	TCut cutSTK1200 = "((((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2.*STK_chargeY_etaCorr[0]+(TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.*STK_chargeX_etaCorr[0])/((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2.+(TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.))>1200.)";

	// 06-11-2023 --> Provo ad implementare un taglio sul PSD scegliendo solo gli eventi che hanno una consistenza tra le due viste del PSD!
	// 				  In particolare: impongo che la differenza tra la carica in una vista (singolo hit o media) abbia una differenza minore di 2 con l'altra vista.
	TCut cutPSD = "fabs((((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_ChargeY0+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_ChargeY1)/((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.))-(((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.*PSD_ChargeX0+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.*PSD_ChargeX1 )/((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.)))<2.";
	
	// -------   TAGLIO COMPLESSIVO   -------
	TCut ctot = cut00*cut01*cut02*cut03*cut05*cutSTK1200*cutPSD;

	TCut bgo01 =    "(BGO_EnergyG_QuenchSatCorr_ML_ions > 100.) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 177.828)";
	TCut bgo02 = "(BGO_EnergyG_QuenchSatCorr_ML_ions > 177.828) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 316.228)";
	TCut bgo03 = "(BGO_EnergyG_QuenchSatCorr_ML_ions > 316.228) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 562.341)";
	TCut bgo04 = "(BGO_EnergyG_QuenchSatCorr_ML_ions > 562.341) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 1000.0)";
	TCut bgo05 =  "(BGO_EnergyG_QuenchSatCorr_ML_ions > 1000.0) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 1778.28)";
	TCut bgo06 = "(BGO_EnergyG_QuenchSatCorr_ML_ions > 1778.28) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 3162.28)";
	TCut bgo07 = "(BGO_EnergyG_QuenchSatCorr_ML_ions > 3162.28) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 5623.41)";
	TCut bgo08 = "(BGO_EnergyG_QuenchSatCorr_ML_ions > 5623.41) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 10000.0)";
	TCut bgo09 = "(BGO_EnergyG_QuenchSatCorr_ML_ions > 10000.0) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 31622.8)";
	TCut bgo10 = "(BGO_EnergyG_QuenchSatCorr_ML_ions > 31622.8) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 100000.0)";

	TString avX = "(((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.*PSD_ChargeX0*PSD_psdX_pathlength[0]+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.*PSD_ChargeX1*PSD_psdX_pathlength[1])/((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.*PSD_psdX_pathlength[0]+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.*PSD_psdX_pathlength[1]))";
	TString avY = "(((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_ChargeY0*PSD_psdY_pathlength[0]+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_ChargeY1*PSD_psdY_pathlength[1])/((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_psdY_pathlength[0]+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_psdY_pathlength[1]))";

	//TString meanX = "(((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.*PSD_ChargeX0+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.*PSD_ChargeX1)/((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.))";
	//TString meanY = "(((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_ChargeY0+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_ChargeY1)/((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.))";

	// -------------------------------------------------------------- 	HISTOGRAMS DEFINITION...
	
	TH1F *h01=new TH1F("h01", "100 GeV < E_{BGO} < 178 GeV",1200,10.0,20.0); h01->GetXaxis()->SetTitle("PSD Y"); h01->GetYaxis()->SetTitle("normalized MC events"); h01->SetLineColor(kPink); h01->SetMarkerColor(kPink); h01->Sumw2();
	TH1F *h02=new TH1F("h02", "178 GeV < E_{BGO} < 316 GeV",1200,10.0,20.0); h02->GetXaxis()->SetTitle("PSD Y"); h02->GetYaxis()->SetTitle("normalized MC events"); h02->SetLineColor(kPink); h02->SetMarkerColor(kPink); h02->Sumw2();
	TH1F *h03=new TH1F("h03", "316 GeV < E_{BGO} < 562 GeV",1200,10.0,20.0); h03->GetXaxis()->SetTitle("PSD Y"); h03->GetYaxis()->SetTitle("normalized MC events"); h03->SetLineColor(kPink); h03->SetMarkerColor(kPink); h03->Sumw2();
	TH1F *h04=new TH1F("h04", "562 GeV < E_{BGO} < 1.0 TeV",1200,10.0,20.0); h04->GetXaxis()->SetTitle("PSD Y"); h04->GetYaxis()->SetTitle("normalized MC events"); h04->SetLineColor(kPink); h04->SetMarkerColor(kPink); h04->Sumw2();
	TH1F *h05=new TH1F("h05", "1.0 TeV < E_{BGO} < 1.8 TeV", 600,10.0,20.0); h05->GetXaxis()->SetTitle("PSD Y"); h05->GetYaxis()->SetTitle("normalized MC events"); h05->SetLineColor(kPink); h05->SetMarkerColor(kPink); h05->Sumw2();
	TH1F *h06=new TH1F("h06", "1.8 TeV < E_{BGO} < 3.2 TeV", 600,10.0,20.0); h06->GetXaxis()->SetTitle("PSD Y"); h06->GetYaxis()->SetTitle("normalized MC events"); h06->SetLineColor(kPink); h06->SetMarkerColor(kPink); h06->Sumw2();
	TH1F *h07=new TH1F("h07", "3.2 TeV < E_{BGO} < 5.6 TeV", 600,10.0,20.0); h07->GetXaxis()->SetTitle("PSD Y"); h07->GetYaxis()->SetTitle("normalized MC events"); h07->SetLineColor(kPink); h07->SetMarkerColor(kPink); h07->Sumw2();
	TH1F *h08=new TH1F("h08", "5.6 TeV < E_{BGO} < 10.0 TeV",600,10.0,20.0); h08->GetXaxis()->SetTitle("PSD Y"); h08->GetYaxis()->SetTitle("normalized MC events"); h08->SetLineColor(kPink); h08->SetMarkerColor(kPink); h08->Sumw2();
	TH1F *h09=new TH1F("h09","10.0 TeV < E_{BGO} < 31.6 TeV",600,10.0,20.0); h09->GetXaxis()->SetTitle("PSD Y"); h09->GetYaxis()->SetTitle("normalized MC events"); h09->SetLineColor(kPink); h09->SetMarkerColor(kPink); h09->Sumw2();
	TH1F *h10=new TH1F("h10","31.6 TeV < E_{BGO} < 100.0 TeV",600,10.0,22.0);h10->GetXaxis()->SetTitle("PSD Y"); h10->GetYaxis()->SetTitle("normalized MC events"); h10->SetLineColor(kPink); h10->SetMarkerColor(kPink); h10->Sumw2();
	
	// -------------------------------------------------------------- 	HERE, THE MAGIC! (FILLING HISTOGRAMS)
	TCanvas *c0=new TCanvas("c0","PSD",1400,800); c0->Divide(5,2);
	c0_1->cd();  c0_1->SetTicks();  c0_1->SetLogy();  c0_2->cd();  c0_2->SetTicks();  c0_2->SetLogy();
	c0_3->cd();  c0_3->SetTicks();  c0_3->SetLogy();  c0_4->cd();  c0_4->SetTicks();  c0_4->SetLogy();
	c0_5->cd();  c0_5->SetTicks();  c0_5->SetLogy();  c0_6->cd();  c0_6->SetTicks();  c0_6->SetLogy();
	c0_7->cd();  c0_7->SetTicks();  c0_7->SetLogy();  c0_8->cd();  c0_8->SetTicks();  c0_8->SetLogy();
	c0_9->cd();  c0_9->SetTicks();  c0_9->SetLogy();  c0_10->cd(); c0_10->SetTicks(); c0_10->SetLogy();
	
	// -------------- h01
	c0_1->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw(avY+">>+h01",ctot*wSN[i]*bgo01,"");} cout<<"... h01 ..."<<endl;
	// -------------- h02
	c0_2->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw(avY+">>+h02",ctot*wSN[i]*bgo02,"");} cout<<"... h02 ..."<<endl;
	// -------------- h03
	c0_3->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw(avY+">>+h03",ctot*wSN[i]*bgo03,"");} cout<<"... h03 ..."<<endl;
	// -------------- h04
	c0_4->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw(avY+">>+h04",ctot*wSN[i]*bgo04,"");} cout<<"... h04 ..."<<endl;
	// -------------- h05
	c0_5->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw(avY+">>+h05",ctot*wSN[i]*bgo05,"");} cout<<"... h05 ..."<<endl;
	// -------------- h06
	c0_6->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw(avY+">>+h06",ctot*wSN[i]*bgo06,"");} cout<<"... h06 ..."<<endl;
	// -------------- h07
	c0_7->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw(avY+">>+h07",ctot*wSN[i]*bgo07,"");} cout<<"... h07 ..."<<endl;
	// -------------- h08
	c0_8->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw(avY+">>+h08",ctot*wSN[i]*bgo08,"");} cout<<"... h08 ..."<<endl;
	// -------------- h09
	c0_9->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw(avY+">>+h09",ctot*wSN[i]*bgo09,"");} cout<<"... h09 ..."<<endl;
	// -------------- h10
	c0_10->cd(); for(int i=0; i<NUM_SET; i++){ skim[i]->Draw(avY+">>+h10",ctot*wSN[i]*bgo10,"");} cout<<"... h10 ..."<<endl;
	
	// -------------------------------------------------------------- 	SAVE ON A FILE...
	TFile *f = TFile::Open("MC_FTFP_SULFUR_PSDYpathweighted_SATcorr_STKcut1200_PSDY.root", "RECREATE");
	f->cd();
	h01->Write(); // save the histogram
	h02->Write(); // save the histogram
	h03->Write(); // save the histogram
	h04->Write(); // save the histogram
	h05->Write(); // save the histogram
	h06->Write(); // save the histogram
	h07->Write(); // save the histogram
	h08->Write(); // save the histogram
	h09->Write(); // save the histogram
	h10->Write(); // save the histogram

	c0->Write();
	f->ls();      // show the contents of the ROOT file
	delete f;     // close the ROOT file

}