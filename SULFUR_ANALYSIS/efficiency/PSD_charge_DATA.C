#include "TChain.h"
#include "TCut.h"
#include "TString.h"
#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "math.h"
#include <string>
#include <iostream>
#include <cstring>
#include <fstream>
#include <vector>


using namespace std;

void PSD_charge_DATA(){

	TChain *skim = new TChain("newtree");
	skim->Add("/mnt/c/Users/saraf/Desktop/dampe/FD_SKIM/skim_FD_20*.root");

	cout<<"Number of entries: "<<skim->GetEntries()<<"\n"<<endl;

	cout<<"Setting branches addresses..."<<endl;
	double psd_y0, psd_y1, psd_x0, psd_x1;
	double psd_max_x, psd_max_y;
	double old_psd_y[2], old_psd_x[2];
	double psd_pl_y[2], psd_pl_x[2];
	double stk_theta;
	double stk_x[6], stk_y[6], stk_etacorr_x[6], stk_etacorr_y[6];
	double e_bgo;
	double e_bgo_lay[14], coo_bgo_max[14], coo_bgo_stk[14];
	int bgo_het;
	skim->SetBranchAddress("BGO_HET",&bgo_het);
	skim->SetBranchAddress("PSD_EnergyY_corr",old_psd_y);
	skim->SetBranchAddress("PSD_EnergyX_corr",old_psd_x);
	skim->SetBranchAddress("PSD_psdY_pathlength",psd_pl_y);
	skim->SetBranchAddress("PSD_psdX_pathlength",psd_pl_x);
	skim->SetBranchAddress("PSD_psdY_Emax_layer",&psd_max_y);
	skim->SetBranchAddress("PSD_psdX_Emax_layer",&psd_max_x);
	skim->SetBranchAddress("PSD_ChargeY0",&psd_y0);
	skim->SetBranchAddress("PSD_ChargeY1",&psd_y1);
	skim->SetBranchAddress("PSD_ChargeX0",&psd_x0);
	skim->SetBranchAddress("PSD_ChargeX1",&psd_x1);
	skim->SetBranchAddress("STK_theta_correction",&stk_theta);
	skim->SetBranchAddress("STK_chargeX",stk_x);
	skim->SetBranchAddress("STK_chargeY",stk_y);
	skim->SetBranchAddress("STK_chargeX_etaCorr",stk_etacorr_x);
	skim->SetBranchAddress("STK_chargeY_etaCorr",stk_etacorr_y);
	skim->SetBranchAddress("BGO_EnergyG",&e_bgo);
	skim->SetBranchAddress("BGO_EneLay",e_bgo_lay);
	skim->SetBranchAddress("BGO_cbgomax",coo_bgo_max);
	skim->SetBranchAddress("BGO_cbgostk",coo_bgo_stk);
	cout<<"done. Now defining histos.."<<endl;

	// Create output file
	string output_path = "/mnt/c/Users/saraf/Desktop/SULFUR_ANALYSIS/efficiency/";
	string output_file_name = "PSD_charge_DATA_STKcut_PSDcut.root";
	TFile *fout1 = new TFile((output_path + output_file_name).c_str(), "RECREATE");
	
	TH1F *hY0=new TH1F("hY0","",1200,0.0,30.);  hY0->GetXaxis()->SetTitle("PSD charge"); hY0->GetYaxis()->SetTitle("events"); hY0->SetLineColor(kAzure+1); hY0->SetMarkerColor(kAzure+1);
	hY0->Sumw2();

	TH1F *hY1=new TH1F("hY1","",1200,0.0,30.);  hY1->GetXaxis()->SetTitle("PSD charge"); hY1->GetYaxis()->SetTitle("events"); hY1->SetLineColor(kOrange+1); hY1->SetMarkerColor(kOrange+1);
	hY1->Sumw2();

	TH1F *hX0=new TH1F("hX0","",1200,0.0,30.);  hX0->GetXaxis()->SetTitle("PSD charge"); hX0->GetYaxis()->SetTitle("events"); hX0->SetLineColor(kViolet-1); hX0->SetMarkerColor(kViolet-1);
	hX0->Sumw2();

	TH1F *hX1=new TH1F("hX1","",1200,0.0,30.);  hX1->GetXaxis()->SetTitle("PSD charge"); hX1->GetYaxis()->SetTitle("events"); hX1->SetLineColor(kTeal-1); hX1->SetMarkerColor(kTeal-1);
	hX1->Sumw2();

	TH1F *hYmean=new TH1F("hYmean","",1200,0.0,30.);  hYmean->GetXaxis()->SetTitle("PSD charge"); hYmean->GetYaxis()->SetTitle("events"); hYmean->SetLineColor(kBlue); hYmean->SetMarkerColor(kBlue);
	hYmean->Sumw2();

	TH1F *hXmean=new TH1F("hXmean","",1200,0.0,30.);  hXmean->GetXaxis()->SetTitle("PSD charge"); hXmean->GetYaxis()->SetTitle("events"); hXmean->SetLineColor(kRed); hXmean->SetMarkerColor(kRed);
	hXmean->Sumw2();
	
	TH1F *hMean=new TH1F("hMean","",1200,0.0,30.);  hMean->GetXaxis()->SetTitle("PSD charge"); hMean->GetYaxis()->SetTitle("events"); 
	hMean->Sumw2();

	cout<<"Done. Reading tree and filling histos"<<endl;

	for(int iEv = 0; iEv<skim->GetEntries(); iEv++){

		bool verbose = 0;
		if(iEv%100000==0){
			verbose = 1;
			cout<<endl;
			cout<<"Processing event: "<<iEv<<endl;
		}

		skim->GetEntry(iEv);

		// List of CUTS
		bool cTrig_HEP = bgo_het>0.;
		// Min deposited energy [GeV]
		bool cEne = e_bgo>100.;
		bool cut00  =  cTrig_HEP*cEne;
		bool cut01  = (psd_y0>0.0 || psd_y1>0.0) && (psd_x0>0.0 || psd_x1>0.0);    //(Elisabetta 02-12-2021)   
		//bool cut02  = e_bgo_lay[0]+e_bgo_lay[1]<e_bgo_lay[2]+e_bgo_lay[3];         // sciame corretto nel bgo SECONDO LA CORREZIONE DI CHUAN
		//bool cut03  = (old_psd_x[0]+old_psd_x[1])*10./(psd_pl_x[0]+psd_pl_x[1])>0.0 && (old_psd_y[0]+old_psd_y[1])*10/(psd_pl_y[0]+psd_pl_y[1])>0.0;  // un segnale su entrambe i due piani di PSD + traccia scelta [QUESTO CUT È INUTILE QUANDO SI HA A CHE FARE CON I FILE SKIMMATI -> LA RICHIESTA DI SEGNALE SU I DUE PSD + LA TRACCIA È AUTOMATICA GIÀ DAL CUT 01 IN QUESTO CASO]
		bool cut05  = fabs(coo_bgo_max[0]-coo_bgo_stk[0])<30.0 && fabs(coo_bgo_max[1]-coo_bgo_stk[1])<30.0;
		bool cut06 = fabs(stk_y[0])< 400. && fabs(stk_x[0]) < 400; // consistency track
		//bool PreCut = LoadGlobalCut();
		bool PreCut = cut00*cut01*cut05*cut06;

		//TCut STKcut = "STK_chargeY[0] > 0. && STK_chargeX[0]>0. && (0.5*(STK_chargeX[0]+STK_chargeY[0])) > 500.";
		bool STKcut = (stk_y[0] > 0.) && (stk_x[0] > 0.) && ((0.5*(stk_x[0]+stk_y[0])) > 500.);

		bool finalCut = PreCut*STKcut;

		if(finalCut){

			// Define PSD X and Y charge variables
			double psd_mean_y = ( ( (TMath::Sign(1.,psd_y0) + 1.)/2. * psd_y0 + (TMath::Sign(1.,psd_y1) + 1.)/2. * psd_y1 ) / ((TMath::Sign(1.,psd_y0) + 1.)/2. + (TMath::Sign(1.,psd_y1) + 1.)/2.) );
			double psd_mean_x = ( ( (TMath::Sign(1.,psd_x0) + 1.)/2. * psd_x0 + (TMath::Sign(1.,psd_x1) + 1.)/2. * psd_x1 ) / ((TMath::Sign(1.,psd_x0) + 1.)/2. + (TMath::Sign(1.,psd_x1) + 1.)/2.) );
			hYmean->Fill(psd_mean_y);
			hXmean->Fill(psd_mean_x);

			bool isPSDY0 = psd_y0 > 0;
			bool isPSDY1 = psd_y1 > 0;
			bool isPSDX0 = psd_x0 > 0;
			bool isPSDX1 = psd_x1 > 0;

			if(isPSDY0) hY0->Fill(psd_y0);
			if(isPSDY1) hY1->Fill(psd_y1);
			if(isPSDX0) hX0->Fill(psd_x0);
			if(isPSDX1) hX1->Fill(psd_x1);

			bool isPSDlay[4] = {isPSDY0, isPSDY1, isPSDX0, isPSDX1};
			double PSDq[4] = {psd_y0, psd_y1, psd_x0, psd_x1};

			double psd_nlay = 0.;
			double psd_sum = 0.;

			std::vector <double> psdq_signal; // Fill with only the PSD charge from layers that have a signal

			for(int iLay = 0; iLay < 4; iLay++){
				if(isPSDlay[iLay]){
					psdq_signal.push_back(PSDq[iLay]);
				}
				if(verbose && isPSDlay[iLay]) cout<<"Signal in PSD layer "<<iLay<<" with value "<<PSDq[iLay]<<endl;
			}

			bool PSD_diff_global = 0;
			bool PSD_diff = 1;
			for(int iQ = 1; iQ < psdq_signal.size(); iQ++){
				PSD_diff *= (TMath::Abs( psdq_signal[iQ-1]- psdq_signal[iQ] ) < 1);
				if(verbose) cout<<"is |"<<psdq_signal[iQ-1]<<" - "<<psdq_signal[iQ]<<"| < 1? "<<PSD_diff<<endl;
				if(iQ == 1 && PSD_diff){
					PSD_diff_global = 1;
					psd_sum += psdq_signal[iQ-1];
					psd_sum += psdq_signal[iQ];
					psd_nlay += 2;
					if(verbose)cout<<"Adding "<<psdq_signal[iQ-1]<<" and "<<psdq_signal[iQ]<<" to PSD charge computation"<<endl;
				}
				else if(iQ != 0 && PSD_diff){
					psd_sum += psdq_signal[iQ];
					psd_nlay += 1;
					if(verbose)cout<<"Adding "<<psdq_signal[iQ]<<" to PSD charge computation"<<endl;
				}
				else break;
			} // for iQ

			double psd_global = psd_sum / psd_nlay;

			if(verbose){
				cout<<"Signal in "<<psd_nlay<<" PSD layers "<<endl;
				cout<<"PSD global charge value: "<<psd_sum<<" / "<<psd_nlay<<" = "<<psd_global<<endl;
				cout<<" Event will be kept? "<<PSD_diff_global<<endl;
			} // if(verbose)

			hMean->Fill(psd_global);

		} // if(finalCut)

	}

	cout<<"Done, now writing histos"<<endl;

	hY0->Write();
	hY1->Write();
	hX0->Write();
	hX1->Write();
	hYmean->Write();
	hXmean->Write();
	hMean->Write();
	cout << "Output written to file: " << output_file_name << endl;
	fout1->Close();

	cout<<"This is the END. Goodbye!"<<endl;

}
