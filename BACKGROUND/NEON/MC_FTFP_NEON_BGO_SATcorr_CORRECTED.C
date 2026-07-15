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
	skim[2]->Add("/mnt/c/Users/saraf/Desktop/dampe/MC/MC_SULFUR/skim_MC_allS32_10TeV_100TeV_v2.root");
    // ------- 100 - 500 TeV
    skim[3]->Add("/mnt/c/Users/saraf/Desktop/dampe/MC/MC_SULFUR/skim_MC_allS32_100TeV_500TeV.root");
        // ------- 500 TeV - 1 PeV
        //skim[4]->Add("/storage/gpfs_data/dampe/users/fsara/simulazioni/zolfo/MC/Si28_500TeV_1PeV_EPOSLHC_FTFP_merged.root");

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
	wS[1] = "(1./51100000.)*log(10.)";
	wS[2] = "(1./10778800.)*log(10.)";
	wS[3] = "(1./1669300.)*log(5.)";
	//wS[4] = "(1./3050580.)*log(2.)";

	TCut wEnergy = "(MC_EnergyT)**(-1.7)";

	TCut wSN[NUM_SET];
	for (int i=0; i<NUM_SET; i++) { wSN[i] = wS[i]*wEnergy; };

	// -------------------------------------------------------------- 	CUTS DEFINITION...

	TCut cTrig_HEP = "BGO_HET>0.";
	TCut cEne = "BGO_EnergyG_QuenchSatCorr_ML_ions>100."; // Min deposited energy [GeV]
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

	TCut bgo01 = "(BGO_EnergyG_QuenchSatCorr_ML_ions > 100.) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 177.828)";              
	TCut bgo02 = "(BGO_EnergyG_QuenchSatCorr_ML_ions > 177.828) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 316.228)";              
	TCut bgo03 = "(BGO_EnergyG_QuenchSatCorr_ML_ions > 316.228) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 562.341)";              	
	TCut bgo04 = "(BGO_EnergyG_QuenchSatCorr_ML_ions > 562.341) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 1000.0)";             
	TCut bgo05 =  "(BGO_EnergyG_QuenchSatCorr_ML_ions > 1000.0) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 1778.28)";              
	TCut bgo06 = "(BGO_EnergyG_QuenchSatCorr_ML_ions > 1778.28) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 3162.28)";              
	TCut bgo07 = "(BGO_EnergyG_QuenchSatCorr_ML_ions > 3162.28) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 5623.41)";              
	TCut bgo08 = "(BGO_EnergyG_QuenchSatCorr_ML_ions > 5623.41) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 10000.0)";              
	TCut bgo09 = "(BGO_EnergyG_QuenchSatCorr_ML_ions > 10000.0) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 31622.8)";              
	TCut bgo10 = "(BGO_EnergyG_QuenchSatCorr_ML_ions > 31622.8) && (BGO_EnergyG_QuenchSatCorr_ML_ions < 100000.0)";

	TString PSDcharge = "PSD_PathWeighted_Charge";

	// 21/10/2025
	TString MPV_MC =  "(14.8272+(1.92572*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(-0.568508*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))+(0.0593389*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3)))";
	//"(14.976+(1.78884*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(-0.528537*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))+(0.0554816*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3)))";
	//"(14.8373+(1.91533*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(-0.566218*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))+(0.0591429*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3)))";
	//"(14.8877+(1.85879*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(-0.547347*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))+(0.0571742*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3)))";
	//"(15.0011+(1.73523*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(-0.502831*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))+(0.0519877*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3)))";
	//"(14.9923+(1.74485*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(-0.506289*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))+(0.0523946*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3)))";

	TString sigma_MC = "(0.248901+(-0.0820845*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(0.0177249*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2)))";
	//"(0.0454068+(0.0415077*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(-0.00228883*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2)))";
	//"(0.117412+(-0.00648231*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(0.00517629*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2)))";
	//"(0.0484749+(0.0695107*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(-0.0185577*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))+(0.00226874*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3)))";
	//"(-0.146673+(0.289607*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(-0.0989563*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))+(0.0115463*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3)))";
	//"(-0.146986+(0.289945*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(-0.0990751*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))+(0.0115598*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3)))";

	TString MPV_DATA = "(14.5821+(1.51344*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.48318*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))+(0.0520877*pow(log10(BGO_EnergyG_SatCorr_ML_ions),3)))";
	//"(14.5821+(1.51344*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.48318*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))+(0.0520877*pow(log10(BGO_EnergyG_SatCorr_ML_ions),3)))";
	//"(15.0798+(0.966632*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(-0.293964*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))+(0.0322678*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3)))";
	//"(15.9434+(0.0730006*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(0.00533945*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))+(-0.000744557*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3)))";
	//"(14.6412+(1.46648*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(-0.475708*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))+(0.0528252*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3)))";
	//"(14.9723+(1.09723*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(-0.342737*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))+(0.0375106*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3)))";

	TString sigma_DATA =  "(0.20977+(0.0478004*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.00269758*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2)))";
	//"(0.20977+(0.0478004*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.00269758*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2)))";
	//"(0.33239+(-0.0579834*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(0.0199986*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2)))";
	//"(0.355433+(-0.0688415*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(0.019678*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2)))";
	//"(-0.087395+(0.396037*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.136666*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))+(0.0168758*pow(log10(BGO_EnergyG_SatCorr_ML_ions),3)))";
	//"(0.138753+(0.158506*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))+(-0.0559774*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))+(0.00789775*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3)))";
					 
	// -------------------------------------------------------------- 	HISTOGRAMS DEFINITION...
	
	TH1F *h01=new TH1F("h01", "100 GeV < E_{BGO} < 178 GeV",1200,10.,20.); h01->GetXaxis()->SetTitle("PSD charge"); h01->GetYaxis()->SetTitle("normalized MC events"); h01->SetLineColor(kViolet); h01->SetMarkerColor(kViolet); h01->Sumw2();
	TH1F *h02=new TH1F("h02", "178 GeV < E_{BGO} < 316 GeV",1200,10.,20.); h02->GetXaxis()->SetTitle("PSD charge"); h02->GetYaxis()->SetTitle("normalized MC events"); h02->SetLineColor(kViolet); h02->SetMarkerColor(kViolet); h02->Sumw2();
	TH1F *h03=new TH1F("h03", "316 GeV < E_{BGO} < 562 GeV",1200,10.,20.); h03->GetXaxis()->SetTitle("PSD charge"); h03->GetYaxis()->SetTitle("normalized MC events"); h03->SetLineColor(kViolet); h03->SetMarkerColor(kViolet); h03->Sumw2();
	TH1F *h04=new TH1F("h04", "562 GeV < E_{BGO} < 1.0 TeV",1200,10.,20.); h04->GetXaxis()->SetTitle("PSD charge"); h04->GetYaxis()->SetTitle("normalized MC events"); h04->SetLineColor(kViolet); h04->SetMarkerColor(kViolet); h04->Sumw2();
	TH1F *h05=new TH1F("h05", "1.0 TeV < E_{BGO} < 1.8 TeV", 1200,10.,20.); h05->GetXaxis()->SetTitle("PSD charge"); h05->GetYaxis()->SetTitle("normalized MC events"); h05->SetLineColor(kViolet); h05->SetMarkerColor(kViolet); h05->Sumw2();
	TH1F *h06=new TH1F("h06", "1.8 TeV < E_{BGO} < 3.2 TeV", 1200,10.,20.); h06->GetXaxis()->SetTitle("PSD charge"); h06->GetYaxis()->SetTitle("normalized MC events"); h06->SetLineColor(kViolet); h06->SetMarkerColor(kViolet); h06->Sumw2();
	TH1F *h07=new TH1F("h07", "3.2 TeV < E_{BGO} < 5.6 TeV", 1000,10.,20.); h07->GetXaxis()->SetTitle("PSD charge"); h07->GetYaxis()->SetTitle("normalized MC events"); h07->SetLineColor(kViolet); h07->SetMarkerColor(kViolet); h07->Sumw2();
	TH1F *h08=new TH1F("h08", "5.6 TeV < E_{BGO} < 10.0 TeV",1000,10.,20.); h08->GetXaxis()->SetTitle("PSD charge"); h08->GetYaxis()->SetTitle("normalized MC events"); h08->SetLineColor(kViolet); h08->SetMarkerColor(kViolet); h08->Sumw2();
	TH1F *h09=new TH1F("h09","10.0 TeV < E_{BGO} < 31.6 TeV",1000,10.,20.); h09->GetXaxis()->SetTitle("PSD charge"); h09->GetYaxis()->SetTitle("normalized MC events"); h09->SetLineColor(kViolet); h09->SetMarkerColor(kViolet); h09->Sumw2();
	TH1F *h10=new TH1F("h10","31.6 TeV < E_{BGO} < 100.0 TeV",500,10.,22.);h10->GetXaxis()->SetTitle("PSD charge"); h10->GetYaxis()->SetTitle("normalized MC events"); h10->SetLineColor(kViolet); h10->SetMarkerColor(kViolet); h10->Sumw2();
	
	// -------------------------------------------------------------- 	HERE, THE MAGIC! (FILLING HISTOGRAMS)
	TCanvas *c0=new TCanvas("c0","PSD",1400,800); c0->Divide(5,2);
	c0_1->cd();  c0_1->SetTicks();  c0_1->SetLogy();  c0_2->cd();  c0_2->SetTicks();  c0_2->SetLogy();
	c0_3->cd();  c0_3->SetTicks();  c0_3->SetLogy();  c0_4->cd();  c0_4->SetTicks();  c0_4->SetLogy();
	c0_5->cd();  c0_5->SetTicks();  c0_5->SetLogy();  c0_6->cd();  c0_6->SetTicks();  c0_6->SetLogy();
	c0_7->cd();  c0_7->SetTicks();  c0_7->SetLogy();  c0_8->cd();  c0_8->SetTicks();  c0_8->SetLogy();
	c0_9->cd();  c0_9->SetTicks();  c0_9->SetLogy();  c0_10->cd(); c0_10->SetTicks(); c0_10->SetLogy();
	
	// -------------- h01
	c0_1->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h01",ctot*wSN[i]*bgo01,"");} cout<<"... h01 ..."<<endl;
	// -------------- h02
	c0_2->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h02",ctot*wSN[i]*bgo02,"");} cout<<"... h02 ..."<<endl;
	// -------------- h03
	c0_3->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h03",ctot*wSN[i]*bgo03,"");} cout<<"... h03 ..."<<endl;
	// -------------- h04
	c0_4->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h04",ctot*wSN[i]*bgo04,"");} cout<<"... h04 ..."<<endl;
	// -------------- h05
	c0_5->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h05",ctot*wSN[i]*bgo05,"");} cout<<"... h05 ..."<<endl;
	// -------------- h06
	c0_6->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h06",ctot*wSN[i]*bgo06,"");} cout<<"... h06 ..."<<endl;
	// -------------- h07
	c0_7->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h07",ctot*wSN[i]*bgo07,"");} cout<<"... h07 ..."<<endl;
	// -------------- h08
	c0_8->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h08",ctot*wSN[i]*bgo08,"");} cout<<"... h08 ..."<<endl;
	// -------------- h09
	c0_9->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h09",ctot*wSN[i]*bgo09,"");} cout<<"... h09 ..."<<endl;
	// -------------- h10
	c0_10->cd(); for(int i=0; i<NUM_SET; i++){ skim[i]->Draw("(("+PSDcharge+"-"+MPV_MC+")*("+sigma_DATA+"/"+sigma_MC+")+"+MPV_DATA+")>>+h10",ctot*wSN[i]*bgo10,"");} cout<<"... h10 ..."<<endl;
	
	// -------------------------------------------------------------- 	SAVE ON A FILE...
	TFile *f = TFile::Open("out_root/MC_FTFP_SULFUR_SATcorr_STKcut1200_PSDcutXY_CORRECTED_100GeV-500TeV_v4.root", "RECREATE");
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
