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
	Int_t nbd = 32;
	Int_t ndec = 8;
	Double_t arg = 1./double(nbd);
	int noe = nbd*ndec;
	//cout << "Number of energy bin " << noe << endl;
	Double_t e0 = 10.;
	std::vector<Double_t> Ebin(noe + 1);
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
	skim[2]->Add("/mnt/c/Users/saraf/Desktop/dampe/MC/MC_SULFUR/skim_MC_allS32_10TeV_100TeV_v2.root");
	// ------- 100 - 500 TeV 
	skim[3]->Add("/mnt/c/Users/saraf/Desktop/dampe/MC/MC_SULFUR/skim_MC_allS32_100TeV_500TeV.root");
	// ----------------------
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
	TFile *foutp = new TFile("../out_root/scatterplot_MC_500TeV_PSD_BGO.root","RECREATE");
	
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
	
	// -------   TAGLIO COMPLESSIVO   -------
	TCut ctot   = cut00*cut01*cut05*cut06*cutSTK1200*cutPSD;

	// -------------------------------------------------------------- 	FILL HISTOGRAMS...

	
/*	TCanvas *c1 = new TCanvas("c1","Canvas",800,600);
	c1->SetLogx();
	c1->SetLogx();
    	gStyle->SetPalette(55);

    // Histogram
    	TH2F *h2psdbgo = new TH2F("h2psdbgo","", noe, Ebin, 600, 10., 18.);
    	h2psdbgo->GetXaxis()->SetTitle("E_{BGO} (GeV)");
    	h2psdbgo->GetYaxis()->SetTitle("PSD charge");
    	h2psdbgo->GetYaxis()->SetTitleOffset(1.3);
    	h2psdbgo->GetXaxis()->SetTitleOffset(1.1);
    	h2psdbgo->GetZaxis()->SetLabelSize(0.03);

    	skim->Draw("PSD_PathWeighted_Charge:BGO_EnergyG_SatCorr_ML_ions >> h2psdbgo", ctot, "COLZ");

    	TFile *foutp = new TFile("scatterplot_DATA_108months_PSD_BGO.root","RECREATE");
    	foutp->cd();
    	foutp->Write();
    	foutp->Close();

	TString MPV_DATA_si = "(13.0366+(1.01399*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.321937*log10(BGO_EnergyG_SatCorr_ML_ions)**2)+(0.036754*log10(BGO_EnergyG_SatCorr_ML_ions)**3))";
   	TString sigma_DATA_si = "(-0.52482+(0.801866*log10(BGO_EnergyG_SatCorr_ML_ions))+(-0.278872*log10(BGO_EnergyG_SatCorr_ML_ions)**2)+(0.0318674*log10(BGO_EnergyG_SatCorr_ML_ions)**3))";

   	TF1 *mpv_si =  new TF1("mpv_si",MPV_DATA_si,100.,60000);
    	mpv_si->SetLineColor(1); mpv_si->SetLineWidth(4); mpv_si->Draw("same");
    	TF1 *slow_si = new TF1("slow_si",MPV_DATA_si+"-"+cut_dw+"*"+sigma_DATA_si,100.,60000);
    	slow_si->SetLineColor(1); slow_si->SetLineWidth(4); slow_si->SetLineStyle(9); slow_si->Draw("same");
    	TF1 *sup_si  = new TF1("sup_si", MPV_DATA_si+"+"+cut_up+"*"+sigma_DATA_si,100.,60000);
    	sup_si->SetLineColor(1); sup_si->SetLineWidth(4); sup_si->SetLineStyle(9); sup_si->Draw("same");

	TLatex *lat_Si_low = new TLatex(120.,13.4,"Si MPV - 2#sigma"); lat_Si_low->SetTextSize(0.027); lat_Si_low->Draw();
    	TLatex *lat_Si     = new TLatex(220.,14.2,"Si MPV"); lat_Si->SetTextSize(0.027); lat_Si->Draw();
    	TLatex *lat_Si_up  = new TLatex(120.,14.9,"Si MPV + 3#sigma"); lat_Si_up->SetTextSize(0.027);  lat_Si_up->Draw();


    	c1->SaveAs("plots/scatterplot_DATA_108months_PSD_BGO.png");
*/ 
    //	c1->SaveAs("scatterplot.png");

	TH2F *h2psdbgo = new TH2F("h2psdbgo","",noe,Ebin.data(),600,12.,25.); h2psdbgo->GetXaxis()->SetTitle("E_{BGO} (GeV)"); h2psdbgo->GetYaxis()->SetTitle("PSD charge ");
	h2psdbgo->GetYaxis()->SetTitleOffset(1.3); h2psdbgo->GetXaxis()->SetTitleOffset(1.1); h2psdbgo->GetZaxis()->SetLabelSize(0.03); 
	/*TH2F *h2psdbgo = new TH2F("h2psdbgo","", noe, Ebin, 600, 10., 18.);
	h2psdbgo->GetXaxis()->SetTitle("E_{BGO} (GeV)");
	h2psdbgo->GetYaxis()->SetTitle("PSD charge");
	h2psdbgo->GetYaxis()->SetTitleOffset(1.3);
	h2psdbgo->GetXaxis()->SetTitleOffset(1.1);
	h2psdbgo->GetZaxis()->SetLabelSize(0.03);
	gPad->SetLogx();
	gStyle->SetPalette(55);
	


	//skim->Draw("PSD_PathWeighted_Charge:BGO_EnergyG_SatCorr_ML_ions>>h2psdbgo",ctot,"");
	//h2psdbgo->Draw("COLZ SAME");
	skim->Draw("PSD_PathWeighted_Charge:BGO_EnergyG_SatCorr_ML_ions >> h2psdbgo",ctot,"COLZ");

*/

	for (int i = 0; i < NUM_SET; i++) {
    skim[i]->Draw("PSD_PathWeighted_Charge:BGO_EnergyG_SatCorr_ML_ions>>h2psdbgo",ctot,"");
}



	cout << " -------------------------------------------" << endl;
	cout << "        MIAO Done!  THIS IS THE END!        " << endl;
	cout << " -------------------------------------------" << endl;

	// -------------------------------------------------------------- 	DATA SAVING...

	foutp->cd();
	foutp->Write();  

}

