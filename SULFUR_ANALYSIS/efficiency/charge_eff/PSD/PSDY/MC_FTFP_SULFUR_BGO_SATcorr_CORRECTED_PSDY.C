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
	wS[1] = "(1./10778800.)*log(10.)";
	wS[2] = "(1./5110000.)*log(10.)";
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
	TString avX = "(((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.*PSD_ChargeX0*PSD_psdX_pathlength[0]+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.*PSD_ChargeX1*PSD_psdX_pathlength[1])/((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.*PSD_psdX_pathlength[0]+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.*PSD_psdX_pathlength[1]))";
	//TString avY = "(((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_ChargeY0*PSD_psdY_pathlength[0]+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_ChargeY1*PSD_psdY_pathlength[1])/((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_psdY_pathlength[0]+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_psdY_pathlength[1]))";

	//TString expr = "((" + avY + "-" + MPV_MC + ")*(" + sigma_DATA + "/" + sigma_MC + ")+" + MPV_DATA + ")";
	
	TString avY_expr = "(((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_ChargeY0*PSD_psdY_pathlength[0]+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_ChargeY1*PSD_psdY_pathlength[1])/((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_psdY_pathlength[0]+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_psdY_pathlength[1]))";

	

auto MPV_MC = [](double E){
    if (E <= 0) return 0.0;
    double x = log10(E);
    return 16.4732+ 0.140747*x + 0.0272845*x*x  - 0.0062576*x*x*x;
};

auto sigma_MC = [](double E){
    if (E <= 0) return 1.0;
    double x = log10(E);
    return 0.266589 - 0.0668461*x + 0.0147683*x*x;
};

auto MPV_DATA = [](double E){
    if (E <= 0) return 0.0;
    double x = log10(E);
    return 14.9041 + 1.03958*x - 0.275602*x*x + 0.0244722*x*x*x;
};

auto sigma_DATA = [](double E){
    if (E <= 0) return 1.0;
    double x = log10(E);
    return - 0.0661661 + 0.308426*x - 0.0437239*x*x;
};

					 
	// -------------------------------------------------------------- 	HISTOGRAMS DEFINITION...
	
	TH1F *h01=new TH1F("h01", "100 GeV < E_{BGO} < 178 GeV",1200,10.,20.); h01->GetXaxis()->SetTitle("PSD Y"); h01->GetYaxis()->SetTitle("normalized MC events"); h01->SetLineColor(kViolet); h01->SetMarkerColor(kViolet); h01->Sumw2();
	TH1F *h02=new TH1F("h02", "178 GeV < E_{BGO} < 316 GeV",1200,10.,20.); h02->GetXaxis()->SetTitle("PSD Y"); h02->GetYaxis()->SetTitle("normalized MC events"); h02->SetLineColor(kViolet); h02->SetMarkerColor(kViolet); h02->Sumw2();
	TH1F *h03=new TH1F("h03", "316 GeV < E_{BGO} < 562 GeV",1200,10.,20.); h03->GetXaxis()->SetTitle("PSD Y"); h03->GetYaxis()->SetTitle("normalized MC events"); h03->SetLineColor(kViolet); h03->SetMarkerColor(kViolet); h03->Sumw2();
	TH1F *h04=new TH1F("h04", "562 GeV < E_{BGO} < 1.0 TeV",1200,10.,20.); h04->GetXaxis()->SetTitle("PSD Y"); h04->GetYaxis()->SetTitle("normalized MC events"); h04->SetLineColor(kViolet); h04->SetMarkerColor(kViolet); h04->Sumw2();
	TH1F *h05=new TH1F("h05", "1.0 TeV < E_{BGO} < 1.8 TeV", 1200,10.,20.); h05->GetXaxis()->SetTitle("PSD Y"); h05->GetYaxis()->SetTitle("normalized MC events"); h05->SetLineColor(kViolet); h05->SetMarkerColor(kViolet); h05->Sumw2();
	TH1F *h06=new TH1F("h06", "1.8 TeV < E_{BGO} < 3.2 TeV", 1200,10.,20.); h06->GetXaxis()->SetTitle("PSD Y"); h06->GetYaxis()->SetTitle("normalized MC events"); h06->SetLineColor(kViolet); h06->SetMarkerColor(kViolet); h06->Sumw2();
	TH1F *h07=new TH1F("h07", "3.2 TeV < E_{BGO} < 5.6 TeV", 1000,10.,20.); h07->GetXaxis()->SetTitle("PSD Y"); h07->GetYaxis()->SetTitle("normalized MC events"); h07->SetLineColor(kViolet); h07->SetMarkerColor(kViolet); h07->Sumw2();
	TH1F *h08=new TH1F("h08", "5.6 TeV < E_{BGO} < 10.0 TeV",1000,10.,20.); h08->GetXaxis()->SetTitle("PSD Y"); h08->GetYaxis()->SetTitle("normalized MC events"); h08->SetLineColor(kViolet); h08->SetMarkerColor(kViolet); h08->Sumw2();
	TH1F *h09=new TH1F("h09","10.0 TeV < E_{BGO} < 31.6 TeV",1000,10.,20.); h09->GetXaxis()->SetTitle("PSD Y"); h09->GetYaxis()->SetTitle("normalized MC events"); h09->SetLineColor(kViolet); h09->SetMarkerColor(kViolet); h09->Sumw2();
	TH1F *h10=new TH1F("h10","31.6 TeV < E_{BGO} < 100.0 TeV",500,10.,22.);h10->GetXaxis()->SetTitle("PSD Y"); h10->GetYaxis()->SetTitle("normalized MC events"); h10->SetLineColor(kViolet); h10->SetMarkerColor(kViolet); h10->Sumw2();
	
	// -------------------------------------------------------------- 	HERE, THE MAGIC! (FILLING HISTOGRAMS)
TCanvas* c0 = new TCanvas("c0","PSD",1400,800);
c0->Divide(5,2);


TH1F* h[10] = {h01,h02,h03,h04,h05,h06,h07,h08,h09,h10};

double BGO_EnergyG_QuenchSatCorr_ML_ions;
double PSD_ChargeY0, PSD_ChargeY1;
double PSD_psdY_pathlength[2];

double Emin[10] = {100., 177.828, 316.228, 562.341, 1000., 1778.28, 3162.28, 5623.41, 10000., 31622.8};
double Emax[10] = {177.828, 316.228, 562.341, 1000., 1778.28, 3162.28, 5623.41, 10000., 31622.8, 100000.};

// ================= LOOP SU FILE =================
for (int i = 0; i < NUM_SET; i++) {
    skim[i]->SetBranchAddress("BGO_EnergyG_QuenchSatCorr_ML_ions", &BGO_EnergyG_QuenchSatCorr_ML_ions);
    skim[i]->SetBranchAddress("PSD_ChargeY0", &PSD_ChargeY0);
    skim[i]->SetBranchAddress("PSD_ChargeY1", &PSD_ChargeY1);
    skim[i]->SetBranchAddress("PSD_psdY_pathlength", PSD_psdY_pathlength);
}

for (int i = 0; i < NUM_SET; i++) {

    Long64_t n = skim[i]->GetEntries();

    for (Long64_t j = 0; j < n; j++) {

        skim[i]->GetEntry(j);

        double E = BGO_EnergyG_QuenchSatCorr_ML_ions;
        if (E <= 0) continue;

        // ================= PSD Y medio =================
        double sY0 = (TMath::Sign(1., PSD_ChargeY0) + 1.) / 2.;
        double sY1 = (TMath::Sign(1., PSD_ChargeY1) + 1.) / 2.;

        double den =
            sY0 * PSD_psdY_pathlength[0] +
            sY1 * PSD_psdY_pathlength[1];

        if (den == 0) continue;

        double num =
            sY0 * PSD_ChargeY0 * PSD_psdY_pathlength[0] +
            sY1 * PSD_ChargeY1 * PSD_psdY_pathlength[1];

        double avY = num / den;

        // ================= CORREZIONE =================
        double value =
            (avY - MPV_MC(E)) *
            (sigma_DATA(E) / sigma_MC(E)) +
            MPV_DATA(E);

        // ================= BIN ENERGIA =================
        for (int ic = 0; ic < 10; ic++) {

            if (E >= Emin[ic] && E < Emax[ic]) {
                h[ic]->Fill(value);
				c0->cd(ic);
				h[ic]->Draw();
                break;
            }
        }
    }

    cout << "file " << i << " done" << endl;
}

	/*
	// -------------- h01
	c0_1->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw(expr + ">>+h01", ctot*wSN[i]*bgo01,"");} cout<<"... h01 ..."<<endl;
	// -------------- h02
	c0_2->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw(expr + ">>+h02", ctot*wSN[i]*bgo02,"");} cout<<"... h02 ..."<<endl;
	// -------------- h03
	c0_3->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw(expr + ">>+h03", ctot*wSN[i]*bgo03,"");} cout<<"... h03 ..."<<endl;
	// -------------- h04
	c0_4->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw(expr + ">>+h04", ctot*wSN[i]*bgo04,"");} cout<<"... h04 ..."<<endl;
	// -------------- h05
	c0_5->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw(expr + ">>+h05", ctot*wSN[i]*bgo05,"");} cout<<"... h05 ..."<<endl;
	// -------------- h06
	c0_6->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw(expr + ">>+h06", ctot*wSN[i]*bgo06,"");} cout<<"... h06 ..."<<endl;
	// -------------- h07
	c0_7->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw(expr + ">>+h07", ctot*wSN[i]*bgo07,"");} cout<<"... h07 ..."<<endl;
	// -------------- h08
	c0_8->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw(expr + ">>+h08", ctot*wSN[i]*bgo08,"");} cout<<"... h08 ..."<<endl;
	// -------------- h09
	c0_9->cd(); for(int i=0; i<NUM_SET; i++){  skim[i]->Draw(expr + ">>+h09", ctot*wSN[i]*bgo09,"");} cout<<"... h09 ..."<<endl;
	// -------------- h10
	c0_10->cd(); for(int i=0; i<NUM_SET; i++){ skim[i]->Draw(expr + ">>+h10", ctot*wSN[i]*bgo10,"");} cout<<"... h10 ..."<<endl;
	*/
	// -------------------------------------------------------------- 	SAVE ON A FILE...
	TFile *f = TFile::Open("MC_FTFP_SULFUR_SATcorr_STKcut1200_PSDcutXY_CORRECTED_PSDY.root", "RECREATE");
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
