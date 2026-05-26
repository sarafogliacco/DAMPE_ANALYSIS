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

#define bins 32 //
#define NUM_SET 4

void CHARGE_efficiency_MC_SULFUR(){

    // -------------------------------------------------------------- 	DEFINIZIONE BIN DI ENERGIA
	Int_t nbd = 8;
	Int_t ndec = 4;
	Double_t arg = 1./double(nbd);
	int noe = nbd*ndec;
	Double_t e0 = 100.;
	Double_t Ebin[bins];
	Ebin[0] = e0;
	for (int i=1; i<noe+1; i++) { Ebin[i] = Ebin[i-1]*TMath::Power(10., arg); }

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
	//cout << "500TeV - 1 PeV :  " <<skim[4]->GetEntries() << endl;
	cout << "-----------------------------------" << endl;
	cout << " " << endl;

	// -------------------------------------------------------------- 	WEIGHTS DEFINITION...
	double wS[NUM_SET]; 
	wS[0] = "(1./37980000.)*log(10.)";
	wS[1] = "(1./10778800.)*log(10.)";
	wS[2] = "(1./5110000.)*log(10.)";
	wS[3] = "(1./1669300.)*log(5.)";
	//wS[4] = (1./3050580.)*log(2.);

    // ------------------------------------------------------------   HISTOGRAMS DEFINITION
    
	TH1F *hCharge_XY = new TH1F("hCharge_XY", "Nobs(Eo) - PSD YX (S-peak)", noe, Ebin); hCharge_XY->Sumw2();
	TH1F *hCharge_Glob = new TH1F("hCharge_Glob", "Nobs(Eo) - PSD YX * PSD charge (S-peak)", noe, Ebin); hCharge_Glob->Sumw2();

    // -----------------------------------------------------------   HERE, THE MAGIC!
	for(int iRange=0; iRange < NUM_SET; iRange++){

		int nentries = skim[iRange]->GetEntries();
		double psd_y0, psd_y1, psd_x0, psd_x1, psd_glob_charge, psd_mean_y, psd_mean_x;
		bool ispsd_signal[4]; 
		double psd_charge_vec[4];
		double psd_max_x, psd_max_y;
		double old_psd_y[2], old_psd_x[2];
		double psd_pl_y[2], psd_pl_x[2];
		double stk_theta;
		double stk_x[6], stk_y[6], stk_etacorr_x[6], stk_etacorr_y[6];
		double e_bgo, e_bgo_quenchSATcor, e_truth, psd_charge;
		double e_bgo_lay[14], coo_bgo_max[14], coo_bgo_stk[14];
		int bgo_het, bgo_let, bgo_ubt;
		float STKtrack_to_PSD_topY;
		float STKtrack_to_PSD_topX;
		skim[iRange]->SetBranchAddress("BGO_HET",&bgo_het);
		skim[iRange]->SetBranchAddress("BGO_LET",&bgo_let);
		skim[iRange]->SetBranchAddress("BGO_UBT",&bgo_ubt);
		skim[iRange]->SetBranchAddress("PSD_EnergyY_corr",old_psd_y);
		skim[iRange]->SetBranchAddress("PSD_EnergyX_corr",old_psd_x);
		skim[iRange]->SetBranchAddress("PSD_psdY_pathlength",psd_pl_y);
		skim[iRange]->SetBranchAddress("PSD_psdX_pathlength",psd_pl_x);
		skim[iRange]->SetBranchAddress("PSD_psdY_Emax_layer",&psd_max_y);
		skim[iRange]->SetBranchAddress("PSD_psdX_Emax_layer",&psd_max_x);
		skim[iRange]->SetBranchAddress("PSD_ChargeY0",&psd_y0);
		skim[iRange]->SetBranchAddress("PSD_ChargeY1",&psd_y1);
		skim[iRange]->SetBranchAddress("PSD_ChargeX0",&psd_x0);
		skim[iRange]->SetBranchAddress("PSD_ChargeX1",&psd_x1);
		skim[iRange]->SetBranchAddress("STK_theta_correction",&stk_theta);
		skim[iRange]->SetBranchAddress("STKtrack_to_PSD_topY", &STKtrack_to_PSD_topY);
		skim[iRange]->SetBranchAddress("STKtrack_to_PSD_topX", &STKtrack_to_PSD_topX);
		skim[iRange]->SetBranchAddress("STK_chargeX_etaCorr",&stk_etacorr_x);
		skim[iRange]->SetBranchAddress("STK_chargeY_etaCorr",&stk_etacorr_y);
		skim[iRange]->SetBranchAddress("BGO_EnergyG_SatCorr_ML_ions",&e_bgo_quenchSATcor);
		skim[iRange]->SetBranchAddress("BGO_EneLay",e_bgo_lay);
		skim[iRange]->SetBranchAddress("BGO_cbgomax",coo_bgo_max);
		skim[iRange]->SetBranchAddress("BGO_cbgostk",coo_bgo_stk);
		skim[iRange]->SetBranchAddress("MC_EnergyT",&e_truth);
		skim[iRange]->SetBranchAddress("PSD_PathWeighted_Charge",&psd_charge);

		for(int iEv = 0; iEv<nentries; iEv++){

			bool verbose = 0;
			if(iEv%100000==0){
				verbose = 1;
				cout<<endl;
				cout<<"Processing event: "<<iEv<<endl;
			}
			skim[iRange]->GetEntry(iEv);

            double wEnergy = pow(e_truth, -1.7);
            // DEFINE TOTAL WEIGHT
			double wSiN = wS[iRange]*wEnergy;

			// List of CUTS
            bool Emin = e_truth > 100;
			bool cTrig_HET = bgo_het>0.;
			bool cEne = e_bgo_quenchSATcor>100.; // Min deposited energy [GeV]
			bool cut00  =  cTrig_HET*cEne;
			bool cut01  = (psd_y0>0.0 || psd_y1>0.0) && (psd_x0>0.0 || psd_x1>0.0);    //(Elisabetta 02-12-2021)   
			//bool cut02  = (e_bgo_lay[0]+e_bgo_lay[1]) < (e_bgo_lay[2]+e_bgo_lay[3]);         // sciame corretto nel bgo SECONDO LA CORREZIONE DI CHUAN
			//bool cut03  = (old_psd_x[0]+old_psd_x[1])*10./(psd_pl_x[0]+psd_pl_x[1])>0.0 && (old_psd_y[0]+old_psd_y[1])*10/(psd_pl_y[0]+psd_pl_y[1])>0.0;  // un segnale su entrambe i due piani di PSD + traccia scelta [QUESTO CUT È INUTILE QUANDO SI HA A CHE FARE CON I FILE SKIMMATI -> LA RICHIESTA DI SEGNALE SU I DUE PSD + LA TRACCIA È AUTOMATICA GIÀ DAL CUT 01 IN QUESTO CASO]
			bool cut05  = fabs(coo_bgo_max[0]-coo_bgo_stk[0])<30.0 && fabs(coo_bgo_max[1]-coo_bgo_stk[1])<30.0;
			bool cut06 = fabs(STKtrack_to_PSD_topY)< 400. && fabs(STKtrack_to_PSD_topX)< 400;

			double meanY = (((TMath::Sign(1.,psd_y0)+1.)/2.*psd_y0+(TMath::Sign(1.,psd_y1)+1.)/2.*psd_y1)/((TMath::Sign(1.,psd_y0)+1.)/2.+(TMath::Sign(1.,psd_y1)+1.)/2.));
			double meanX = (((TMath::Sign(1.,psd_x0)+1.)/2.*psd_x0+(TMath::Sign(1.,psd_x1)+1.)/2.*psd_x1)/((TMath::Sign(1.,psd_x0)+1.)/2.+(TMath::Sign(1.,psd_x1)+1.)/2.));
			
			double avX = (((TMath::Sign(1.,psd_x0)+1.)/2.*psd_x0*psd_pl_x[0]+(TMath::Sign(1.,psd_x1)+1.)/2.*psd_x1*psd_pl_x[1])/((TMath::Sign(1.,psd_x0)+1.)/2.*psd_pl_x[0]+(TMath::Sign(1.,psd_x1)+1.)/2.*psd_pl_x[1]));
			double avY = (((TMath::Sign(1.,psd_y0)+1.)/2.*psd_y0*psd_pl_y[0]+(TMath::Sign(1.,psd_y1)+1.)/2.*psd_y1*psd_pl_y[1])/((TMath::Sign(1.,psd_y0)+1.)/2.*psd_pl_y[0]+(TMath::Sign(1.,psd_y1)+1.)/2.*psd_pl_y[1]));

			double cutSTK1200 = ((((TMath::Sign(1.,stk_etacorr_y[0])+1.)/2.*stk_etacorr_y[0]+(TMath::Sign(1.,stk_etacorr_x[0])+1.)/2.*stk_etacorr_x[0])/((TMath::Sign(1.,stk_etacorr_y[0])+1.)/2.+(TMath::Sign(1.,stk_etacorr_x[0])+1.)/2.))>1200.);
			double cutPSD = fabs(meanY-meanX)<2.;

			bool PreCut = Emin*cut00*cut01*cut06*cut05*cutSTK1200*cutPSD;

			if(PreCut){	

				// 29/08/2024 - PSD Y
				double cutLow_PSDY = 2.0;
				double cutUp_PSDY  = 2.5;

				double MPV_MC_PSDY = (15.5725+(1.01584*log10(e_bgo_quenchSATcor))+(-0.248409*pow(log10(e_bgo_quenchSATcor),2))+(0.02147*pow(log10(e_bgo_quenchSATcor),3)));
				double sigma_MC_PSDY = (0.232387+(-0.05435*log10(e_bgo_quenchSATcor))+(0.0111726*pow(log10(e_bgo_quenchSATcor),2)));

				double MPV_DATA_PSDY = ;
				double sigma_DATA_PSDY = ;

				double QcorrMC_PSDY = ((avY-MPV_MC_PSDY)*(sigma_DATA_PSDY/sigma_MC_PSDY)+MPV_DATA_PSDY);
				double MPV_MC_CORR_PSDY = (13.4701+(0.557324*log10(e_bgo_quenchSATcor))+(-0.19535*pow(log10(e_bgo_quenchSATcor),2))+(0.0254733*pow(log10(e_bgo_quenchSATcor),3)));
				double sigma_MC_CORR_PSDY = (0.23415+(0.0159814*log10(e_bgo_quenchSATcor)));

				bool chargeUP_PSDY = QcorrMC_PSDY<(MPV_MC_CORR_PSDY+cutUp_PSDY *sigma_MC_CORR_PSDY);
				bool chargeDW_PSDY = QcorrMC_PSDY>(MPV_MC_CORR_PSDY-cutLow_PSDY*sigma_MC_CORR_PSDY);
				bool chargeCUT_PSDY = chargeUP_PSDY*chargeDW_PSDY;


				//29/08/2024 - PSD X
				double cutLow_PSDX = 2.0;
				double cutUp_PSDX  = 2.5;

				double MPV_MC_PSDX = (14.5692+(2.21805*log10(e_bgo_quenchSATcor))+(-0.65476*pow(log10(e_bgo_quenchSATcor),2))+(0.0680206*pow(log10(e_bgo_quenchSATcor),3)));
				double sigma_MC_PSDX = (0.0810492+(0.0505525*log10(e_bgo_quenchSATcor))+(-0.00352818*pow(log10(e_bgo_quenchSATcor),2)));

				double MPV_DATA_PSDX = (12.4454+(1.61755*log10(e_bgo_quenchSATcor))+(-0.518264*pow(log10(e_bgo_quenchSATcor),2))+(0.058086*pow(log10(e_bgo_quenchSATcor),3)));
				double sigma_DATA_PSDX = (0.295111+(0.00769167*log10(e_bgo_quenchSATcor)));

				double QcorrMC_PSDX = ((avX-MPV_MC_PSDX)*(sigma_DATA_PSDX/sigma_MC_PSDX)+MPV_DATA_PSDX);
				double MPV_MC_CORR_PSDX = (12.4454+(1.60961*log10(e_bgo_quenchSATcor))+(-0.514021*pow(log10(e_bgo_quenchSATcor),2))+(0.0576378*pow(log10(e_bgo_quenchSATcor),3)));
				double sigma_MC_CORR_PSDX = (0.292812+(0.00905601*log10(e_bgo_quenchSATcor)));

				bool chargeUP_PSDX = QcorrMC_PSDX<(MPV_MC_CORR_PSDX+cutUp_PSDX *sigma_MC_CORR_PSDX);

				bool chargeDW_PSDX = QcorrMC_PSDX>(MPV_MC_CORR_PSDX-cutLow_PSDX*sigma_MC_CORR_PSDX);
				bool chargeCUT_PSDX = chargeUP_PSDX*chargeDW_PSDX;


				// 28/08/2024
				double cutLow = 2.0;
				double cutUp = 2.5;
				double MPV_MC = (13.047+(1.6516*log10(e_bgo_quenchSATcor))+(-0.504739*pow(log10(e_bgo_quenchSATcor),2))+(0.0549877*pow(log10(e_bgo_quenchSATcor),3)));
				double sigma_MC = (-0.0760008+(0.233645*log10(e_bgo_quenchSATcor))+(-0.0875578*pow(log10(e_bgo_quenchSATcor),2))+(0.0114994*pow(log10(e_bgo_quenchSATcor),3)));
				double MPV_DATA = (13.0366+(1.01399*log10(e_bgo_quenchSATcor))+(-0.321937*pow(log10(e_bgo_quenchSATcor),2))+(0.036754*pow(log10(e_bgo_quenchSATcor),3)));
				double sigma_DATA = (0.205642+(0.00990178*log10(e_bgo_quenchSATcor)));
				double QcorrMC = ((psd_charge-MPV_MC)*(sigma_DATA/sigma_MC)+MPV_DATA);
				double MPV_MC_CORR = (13.0384+(1.03416*log10(e_bgo_quenchSATcor))+(-0.335197*pow(log10(e_bgo_quenchSATcor),2))+(0.0389124*pow(log10(e_bgo_quenchSATcor),3)));
				double sigma_MC_CORR = (-0.618361+(0.892525*log10(e_bgo_quenchSATcor))+(-0.307297*pow(log10(e_bgo_quenchSATcor),2))+(0.0346831*pow(log10(e_bgo_quenchSATcor),3)));

				bool chargeUP = QcorrMC<(MPV_MC_CORR+cutUp*sigma_MC_CORR);
				bool chargeDW = QcorrMC>(MPV_MC_CORR-cutLow*sigma_MC_CORR);
				bool chargeCUT = chargeUP*chargeDW;

				// --------------------------------------------------------------------------------------
				
				bool ctoto_0s   = chargeCUT_PSDY*chargeCUT_PSDX;
				bool ctoto_ANDs = chargeCUT_PSDY*chargeCUT_PSDX*chargeCUT;
				
				if(ctoto_0s) hCharge_XY->Fill(e_bgo_quenchSATcor,wSiN);
				if(ctoto_ANDs) hCharge_Glob->Fill(e_bgo_quenchSATcor,wSiN);
				
			} else continue; // if PreCuts && STK cut
		} // loop on tree events
	} // iRange 

	// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	TFile *fout = new TFile("out/CHARGE_eff_MC_SULFUR_SATcorr_STKcut1200_PSDcutXY_8bin.root","RECREATE"); 
    hCharge_XY->Write(); 
    hCharge_Glob->Write(); 

    TH1F *h1_PSD_Eff = new TH1F("h1_PSD_Eff", "CHARGE-Eff vs Ebgo (S-peak)", noe, Ebin); h1_PSD_Eff->Sumw2();
    h1_PSD_Eff->Divide(hCharge_Glob,hCharge_XY);
    h1_PSD_Eff->Write();
	
	fout->ls();
	fout->Close();

}
