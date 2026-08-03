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

	const int NUM_SET = 5;

	// -------------------------------------------------------------- 	DATA LOADING...
	TChain *skim[NUM_SET];
	for(int i=0; i < NUM_SET; i++){ skim[i] = new TChain("newtree"); }
	// ------- 100 GeV - 1 TeV
	skim[0]->Add("/mnt/c/Users/saraf/Desktop/dampe/MC/OXYGEN/O16_100GeV_1TeV_merged.root");
	// ------- 1 - 10 TeV 
	skim[1]->Add("/mnt/c/Users/saraf/Desktop/dampe/MC/OXYGEN/O16_1TeV_10TeV_merged.root");
	// ------- 10 - 100 TeV 
	skim[2]->Add("/mnt/c/Users/saraf/Desktop/dampe/MC/OXYGEN/O16_10TeV_100TeV_merged.root");
	// ------- 100 - 500 TeV 
	skim[3]->Add("/mnt/c/Users/saraf/Desktop/dampe/MC/OXYGEN/O16_100TeV_500TeV_merged.root");
	// ------- 500 TeV - 1 PeV
	skim[4]->Add("/mnt/c/Users/saraf/Desktop/dampe/MC/OXYGEN/O16_500TeV_1PeV_merged.root");


	// Number of entries ...
	cout << " " << endl;
	cout << "-----------------------------------" << endl;
	cout << "Number of entries: " << endl;
	cout << "100GeV - 1TeV  : " << skim[0]->GetEntries() << endl;
	cout << "  1TeV - 10TeV : " << skim[1]->GetEntries() << endl;
	cout << " 10TeV - 100TeV: " << skim[2]->GetEntries() << endl;
	cout << "100TeV - 500TeV:  " <<skim[3]->GetEntries() << endl;
	cout << "500TeV - 1 PeV :  " <<skim[4]->GetEntries() << endl;
	cout << "-----------------------------------" << endl;
	cout << " " << endl;

	// -------------------------------------------------------------- 	WEIGHTS DEFINITION...
	
	TCut wSi[NUM_SET]; 


	TCut wEnergy = "(MC_EnergyT)**(-1.7)";

	TCut wSiN[NUM_SET];
	for (int i=0; i<NUM_SET; i++) { wSiN[i] = wSi[i]*wEnergy; };

	// -------------------------------------------------------------- 	CUTS DEFINITION...

	TCut cTrig_HEP = "BGO_HET>0.";
	TCut cEne = "BGO_EnergyG_QuenchSatCorr_ML_ions2>100."; // Min deposited energy [GeV]
	TCut cut00  = cTrig_HEP*cEne;
	TCut cut01  = "(PSD_ChargeY0>0.0 || PSD_ChargeY1>0.0) && (PSD_ChargeX0>0.0 || PSD_ChargeX1>0.0)";    //(Elisabetta 02-12-2021)   

	TCut cut05  = "fabs(BGO_cbgomax[0]-BGO_cbgostk[0])<30.0 && fabs(BGO_cbgomax[1]-BGO_cbgostk[1])<30.0";

	TCut cut06 = "fabs(STKtrack_to_PSD_topY)< 400. && fabs(STKtrack_to_PSD_topX) < 400";

	// STK cut (Elisabetta 22-02-2024)
	TCut cutSTK1200 = "((((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2.*STK_chargeY_etaCorr[0] + (TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.*STK_chargeX_etaCorr[0]) / ((TMath::Sign(1.,STK_chargeY_etaCorr[0])+1.)/2. + (TMath::Sign(1.,STK_chargeX_etaCorr[0])+1.)/2.) ) >1200.)";

	// 06-11-2023 --> Provo ad implementare un taglio sul PSD scegliendo solo gli eventi che hanno una consistenza tra le due viste del PSD!
	// 				  In particolare: impongo che la differenza tra la carica in una vista (singolo hit o media) abbia una differenza minore di 2 con l'altra vista.
	TCut cutPSD = "fabs(( ( (TMath::Sign(1.,PSD_ChargeY0) + 1.)/2. * PSD_ChargeY0 + (TMath::Sign(1.,PSD_ChargeY1) + 1.)/2. * PSD_ChargeY1 ) / ((TMath::Sign(1.,PSD_ChargeY0) + 1.)/2. + (TMath::Sign(1.,PSD_ChargeY1) + 1.)/2.) ) - ( ( (TMath::Sign(1.,PSD_ChargeX0) + 1.)/2. * PSD_ChargeX0 + (TMath::Sign(1.,PSD_ChargeX1) + 1.)/2. * PSD_ChargeX1 ) / ((TMath::Sign(1.,PSD_ChargeX0) + 1.)/2. + (TMath::Sign(1.,PSD_ChargeX1) + 1.)/2.) ) ) < 2.";
		
	// -------   TAGLIO COMPLESSIVO   -------
	TCut ctot = cut00*cut01*cut05*cut06*cutSTK1200*cutPSD;

	TCut bgo01 = "(BGO_EnergyG_QuenchSatCorr_ML_ions2 > 100.)    && (BGO_EnergyG_QuenchSatCorr_ML_ions2 < 177.828)";              
	TCut bgo02 = "(BGO_EnergyG_QuenchSatCorr_ML_ions2 > 177.828) && (BGO_EnergyG_QuenchSatCorr_ML_ions2 < 316.228)";              
	TCut bgo03 = "(BGO_EnergyG_QuenchSatCorr_ML_ions2 > 316.228) && (BGO_EnergyG_QuenchSatCorr_ML_ions2 < 562.341)";              	
	TCut bgo04 = "(BGO_EnergyG_QuenchSatCorr_ML_ions2 > 562.341) && (BGO_EnergyG_QuenchSatCorr_ML_ions2 < 1000.0)";             
	TCut bgo05 = "(BGO_EnergyG_QuenchSatCorr_ML_ions2 > 1000.0)  && (BGO_EnergyG_QuenchSatCorr_ML_ions2 < 1778.28)";              
	TCut bgo06 = "(BGO_EnergyG_QuenchSatCorr_ML_ions2 > 1778.28) && (BGO_EnergyG_QuenchSatCorr_ML_ions2 < 3162.28)";              
	TCut bgo07 = "(BGO_EnergyG_QuenchSatCorr_ML_ions2 > 3162.28) && (BGO_EnergyG_QuenchSatCorr_ML_ions2 < 5623.41)";              
	TCut bgo08 = "(BGO_EnergyG_QuenchSatCorr_ML_ions2 > 5623.41) && (BGO_EnergyG_QuenchSatCorr_ML_ions2 < 10000.0)";              
	TCut bgo09 = "(BGO_EnergyG_QuenchSatCorr_ML_ions2 > 10000.0) && (BGO_EnergyG_QuenchSatCorr_ML_ions2 < 31622.8)";              
	TCut bgo10 = "(BGO_EnergyG_QuenchSatCorr_ML_ions2 > 31622.8) && (BGO_EnergyG_QuenchSatCorr_ML_ions2 < 100000.0)";
	TCut bgo11 = "(BGO_EnergyG_QuenchSatCorr_ML_ions2 > 100000.0) && (BGO_EnergyG_QuenchSatCorr_ML_ions2 < 316227.0)";
	TCut bgo12 = "(BGO_EnergyG_QuenchSatCorr_ML_ions2 > 316227.0) && (BGO_EnergyG_QuenchSatCorr_ML_ions2 < 1000000.0)";

	TString PSDcharge = "PSD_PathWeighted_Charge";

	// 21/10/2025
	TString MPV_MC =  
	TString sigma_MC =  

	TString MPV_DATA = 
	TString sigma_DATA =  

					 
	// -------------------------------------------------------------- 	HISTOGRAMS DEFINITION...
	
	TH1F *h01=new TH1F("h01", "100 GeV < E_{BGO} < 178 GeV",1200, 1.5,27.5); h01->GetXaxis()->SetTitle("PSD charge"); h01->GetYaxis()->SetTitle("normalized MC events"); h01->SetLineColor(kViolet); h01->SetMarkerColor(kViolet); h01->Sumw2();
	TH1F *h02=new TH1F("h02", "178 GeV < E_{BGO} < 316 GeV",1200, 1.5,27.5); h02->GetXaxis()->SetTitle("PSD charge"); h02->GetYaxis()->SetTitle("normalized MC events"); h02->SetLineColor(kViolet); h02->SetMarkerColor(kViolet); h02->Sumw2();
	TH1F *h03=new TH1F("h03", "316 GeV < E_{BGO} < 562 GeV",1200, 1.5,27.5); h03->GetXaxis()->SetTitle("PSD charge"); h03->GetYaxis()->SetTitle("normalized MC events"); h03->SetLineColor(kViolet); h03->SetMarkerColor(kViolet); h03->Sumw2();
	TH1F *h04=new TH1F("h04", "562 GeV < E_{BGO} < 1.0 TeV",1200, 1.5,27.5); h04->GetXaxis()->SetTitle("PSD charge"); h04->GetYaxis()->SetTitle("normalized MC events"); h04->SetLineColor(kViolet); h04->SetMarkerColor(kViolet); h04->Sumw2();
	TH1F *h05=new TH1F("h05", "1.0 TeV < E_{BGO} < 1.8 TeV", 1200,1.5,27.5); h05->GetXaxis()->SetTitle("PSD charge"); h05->GetYaxis()->SetTitle("normalized MC events"); h05->SetLineColor(kViolet); h05->SetMarkerColor(kViolet); h05->Sumw2();
	TH1F *h06=new TH1F("h06", "1.8 TeV < E_{BGO} < 3.2 TeV", 1200,1.5,27.5); h06->GetXaxis()->SetTitle("PSD charge"); h06->GetYaxis()->SetTitle("normalized MC events"); h06->SetLineColor(kViolet); h06->SetMarkerColor(kViolet); h06->Sumw2();
	TH1F *h07=new TH1F("h07", "3.2 TeV < E_{BGO} < 5.6 TeV", 1000,1.5,27.5); h07->GetXaxis()->SetTitle("PSD charge"); h07->GetYaxis()->SetTitle("normalized MC events"); h07->SetLineColor(kViolet); h07->SetMarkerColor(kViolet); h07->Sumw2();
	TH1F *h08=new TH1F("h08", "5.6 TeV < E_{BGO} < 10.0 TeV",1000,1.5,27.5); h08->GetXaxis()->SetTitle("PSD charge"); h08->GetYaxis()->SetTitle("normalized MC events"); h08->SetLineColor(kViolet); h08->SetMarkerColor(kViolet); h08->Sumw2();
	TH1F *h09=new TH1F("h09","10.0 TeV < E_{BGO} < 31.6 TeV",1000,1.5,27.5); h09->GetXaxis()->SetTitle("PSD charge"); h09->GetYaxis()->SetTitle("normalized MC events"); h09->SetLineColor(kViolet); h09->SetMarkerColor(kViolet); h09->Sumw2();
	TH1F *h10=new TH1F("h10","31.6 TeV < E_{BGO} < 100.0 TeV",500,1.5,27.5);h10->GetXaxis()->SetTitle("PSD charge"); h10->GetYaxis()->SetTitle("normalized MC events"); h10->SetLineColor(kViolet); h10->SetMarkerColor(kViolet); h10->Sumw2();
	//TH1F *h11=new TH1F("h11","100.0 TeV < E_{BGO} < 316.2 TeV",500,10.,18.);h11->GetXaxis()->SetTitle("PSD charge"); h11->GetYaxis()->SetTitle("normalized MC events"); h11->SetLineColor(kViolet); h11->SetMarkerColor(kViolet); h11->Sumw2();
	//TH1F *h12=new TH1F("h12","316.2 TeV < E_{BGO} < 1000.0 TeV",500,10.,20.);h12->GetXaxis()->SetTitle("PSD charge"); h12->GetYaxis()->SetTitle("normalized MC events"); h12->SetLineColor(kViolet); h12->SetMarkerColor(kViolet); h12->Sumw2();


	// -------------------------------------------------------------- 	HERE, THE MAGIC! (FILLING HISTOGRAMS)
	TCanvas *c0=new TCanvas("c0","PSD",1400,800); c0->Divide(5,2);
	c0_1->cd();  c0_1->SetTicks();  c0_1->SetLogy();  c0_2->cd();  c0_2->SetTicks();  c0_2->SetLogy();
	c0_3->cd();  c0_3->SetTicks();  c0_3->SetLogy();  c0_4->cd();  c0_4->SetTicks();  c0_4->SetLogy();
	c0_5->cd();  c0_5->SetTicks();  c0_5->SetLogy();  c0_6->cd();  c0_6->SetTicks();  c0_6->SetLogy();
	c0_7->cd();  c0_7->SetTicks();  c0_7->SetLogy();  c0_8->cd();  c0_8->SetTicks();  c0_8->SetLogy();
	c0_9->cd();  c0_9->SetTicks();  c0_9->SetLogy();  c0_10->cd(); c0_10->SetTicks(); c0_10->SetLogy();
	//c0_11->cd(); c0_11->SetTicks(); c0_11->SetLogy(); c0_12->cd(); c0_12->SetTicks(); c0_12->SetLogy();
	
	// -------------- h01
	c0_1->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h01",ctot*wSiN[i]*bgo01,"");} cout<<"... h01 ..."<<endl;
	// -------------- h02
	c0_2->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h02",ctot*wSiN[i]*bgo02,"");} cout<<"... h02 ..."<<endl;
	// -------------- h03
	c0_3->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h03",ctot*wSiN[i]*bgo03,"");} cout<<"... h03 ..."<<endl;
	// -------------- h04
	c0_4->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h04",ctot*wSiN[i]*bgo04,"");} cout<<"... h04 ..."<<endl;
	// -------------- h05
	c0_5->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h05",ctot*wSiN[i]*bgo05,"");} cout<<"... h05 ..."<<endl;
	// -------------- h06
	c0_6->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h06",ctot*wSiN[i]*bgo06,"");} cout<<"... h06 ..."<<endl;
	// -------------- h07
	c0_7->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h07",ctot*wSiN[i]*bgo07,"");} cout<<"... h07 ..."<<endl;
	// -------------- h08
	c0_8->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h08",ctot*wSiN[i]*bgo08,"");} cout<<"... h08 ..."<<endl;
	// -------------- h09
	c0_9->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h09",ctot*wSiN[i]*bgo09,"");} cout<<"... h09 ..."<<endl;
	// -------------- h10
	c0_10->cd(); for(int i=0; i<NUM_SET; i++){ skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h10",ctot*wSiN[i]*bgo10,"");} cout<<"... h10 ..."<<endl;
	// -------------- h11
	//c0_11->cd(); for(int i=0; i<NUM_SET; i++){ skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h11",ctot*wSiN[i]*bgo11,"");} cout<<"... h11 ..."<<endl;
	// -------------- h12
	//c0_12->cd(); for(int i=0; i<NUM_SET; i++){ skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h12",ctot*wSiN[i]*bgo12,"");} cout<<"... h12 ..."<<endl;


	// -------------------------------------------------------------- 	SAVE ON A FILE...
	TFile *f = TFile::Open("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/BACKGROUND/out_root/MC_FTFP_OXYGEN_SATcorr_STKcut1200_PSDcutXY_CORRECTED_mcfit.root", "RECREATE");
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
	//h11->Write(); // save the histogram
	//h12->Write(); // save the histogram

	c0->Write();
	f->ls();      // show the contents of the ROOT file
	delete f;     // close the ROOT file

}
