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
	#include <cmath>

	const int NUM_SET = 5;

	std::vector<TString> files = {
        "/mnt/c/Users/saraf/Desktop/dampe/MC/MC_SILICON/skim_new/skim_MC_allSi28-v6r0p15_100GeV_1TeV_FTFP_merged.root",
        "/mnt/c/Users/saraf/Desktop/dampe/MC/MC_SILICON/skim_new/skim_MC_allSi28-v6r0p15_1TeV_10TeV_FTFP_merged.root",
        "/mnt/c/Users/saraf/Desktop/dampe/MC/MC_SILICON/skim_new/skim_MC_allSi28-v6r0p15_10TeV_100TeV_EPOSLHC_FTFP_merged.root",
        "/mnt/c/Users/saraf/Desktop/dampe/MC/MC_SILICON/skim_new/skim_MC_allSi28-v6r0p15_100TeV_500TeV_EPOSLHC_FTFP_merged.root",
        //"/mnt/c/Users/saraf/Desktop/dampe/MC/MC_SILICON/skim_new/skim_MC_allSi28-v6r0p15_500TeV_1PeV_EPOSLHC_FTFP_merged.root"
    };

	// -------------------------------------------------------------- 	DEFINIZIONE BIN DI ENERGIA
	
    const int n1=6, n2=6, n3=6, n4=6;
    const int noe = n1+n2+n3+n4;
    Double_t Ebin[noe+1];
    int k=0;
    // 100 GeV - 1 TeV 
    for (int i=0; i<=n1; i++) { Ebin[k++] = TMath::Power(10, 2 + i*(1.0/double(n1))); }
    // 1 TeV - 10 TeV
    for (int i=1; i<=n2; i++) { Ebin[k++] = TMath::Power(10, 3 + i*(1.0/double(n2))); }
    // 10 TeV - 100 TeV
    for (int i=1; i<=n3; i++) { Ebin[k++] = TMath::Power(10, 4 + i*(1.0/double(n3))); }
	// 100 TeV - 500 TeV
	/*
	double Emin = 1e5;   // 100 TeV
	double Emax = 5e5;   // 500 TeV
	double logEmin = log10(Emin);
	double logEmax = log10(Emax);
	for (int i=1; i<=n4; i++) {
    	double x = logEmin + i * (logEmax - logEmin) / double(n4);
    	Ebin[k++] = pow(10, x);
		}
	*/
    // 100 TeV - 1 PeV
    for (int i=1; i<=n4; i++) { Ebin[k++] = TMath::Power(10, 5 + i*(1.0/double(n4))); }
    cout << k << endl;
    cout << "Number of energy bins: " << noe << endl;
    for (int i=0; i<noe+1; i++) cout << Ebin[i] << endl; 
    

	// -------------------------------------------------------------  DATA LOADING...
	TChain *skim[NUM_SET];
	for(int i=0; i < NUM_SET; i++){ skim[i] = new TChain("newtree"); }
	// ------- 100 GeV - 1 TeV
	skim[0]->Add(files[0]);
	// ------- 1 - 10 TeV 
	skim[1]->Add(files[1]);
	// ------- 10 - 100 TeV 
	skim[2]->Add(files[2]);
	// ------- 100 - 500 TeV 
	skim[3]->Add(files[3]);
	// ------- 500 TeV - 1 PeV
	skim[4]->Add(files[4]);
	// ------- 500 TeV - 1 PeV
	//skim[5]->Add("../../DATA/Ne20_1PeV_5PeV_test1.root");

	// Number of entries ...
	cout << " " << endl;
	cout << "-----------------------------------" << endl;
	cout << "Number of entries: " << endl;
	cout << "100GeV - 1TeV  : " << skim[0]->GetEntries() << endl;
	cout << "  1TeV - 10TeV : " << skim[1]->GetEntries() << endl;
	cout << " 10TeV - 100TeV: " << skim[2]->GetEntries() << endl;
	cout << "100TeV - 500TeV:  " <<skim[3]->GetEntries() << endl;
	//cout << "500TeV - 1 PeV :  " <<skim[4]->GetEntries() << endl;
	//cout << " 1 PeV - 5 PeV :  " <<skim[5]->GetEntries() << endl;
	cout << "-----------------------------------" << endl;
	cout << " " << endl;

	// -------------------------------------------------------------- 	CUTS DEFINITION...

	    auto MPV_MC = [](double E){
             double x = log10(E);
             return 13.047 + 1.6516*x - 0.504739*x*x + 0.0549877*x*x*x;
        };

        auto sigma_MC = [](double E){
             double x = log10(E);
             return -0.0760008 + 0.233645*x - 0.0875578*x*x + 0.0114994*x*x*x;
        };

        auto MPV_DATA = [](double E){
             double x = log10(E);
             return 13.0366 + 1.01399*x - 0.321937*x*x + 0.036754*x*x*x;
        };

        auto sigma_DATA = [](double E){
             double x = log10(E);
             return 0.205642 + 0.00990178*x;
        };

        auto MPV_MC_CORR = [](double E){
             double x = log10(E);
             return 13.0384 + 1.03416*x - 0.335197*x*x + 0.0389124*x*x*x;
        };

        auto sigma_MC_CORR = [](double E){
             double x = log10(E);
             return -0.618361 + 0.892525*x - 0.307297*x*x + 0.0346831*x*x*x;
        };


	// -------------------------------------------------------------- 	HERE, THE MAGIC! (FILLING HISTOGRAMS)

	TFile *fout = new TFile("Out_MC_FTFP_Si_E2e7_6bin_500TeV.root","RECREATE");

	TFile *p0 = TFile::Open(files[0], "READ");
    if (!p0 || p0->IsZombie()) {
        std::cerr << "Errore: non riesco ad aprire il file\n";
        return;
    }
    TFile *p1 = TFile::Open(files[1], "READ");
    TFile *p2 = TFile::Open(files[2], "READ");
    TFile *p3 = TFile::Open(files[3], "READ");
    //TFile *p4 = TFile::Open(files[4], "READ");
    
    // 
    TH1F *h00 = (TH1F*) p0->Get("h_energy_truth_6bins_weight_E2e7_cut_0"); // 100 GeV - 1 TeV
    TH1F *h01 = (TH1F*) p1->Get("h_energy_truth_6bins_weight_E2e7_cut_0"); // 1 TeV - 10 TeV
    TH1F *h02 = (TH1F*) p2->Get("h_energy_truth_6bins_weight_E2e7_cut_0");  // 10 TeV - 100 TeV
    TH1F *h03 = (TH1F*) p3->Get("h_energy_truth_6bins_weight_E2e7_cut_0");  // 100 TeV - 500 TeV
    //TH1F *h04 = (TH1F*) p4->Get("h_energy_truth_6bins_weight_E2e7_cut_0");  // 500 TeV - 1 PeV


    if (!h00 || !h01 || !h02 || !h03 /*|| !h04*/) {
        std::cerr << "Errore: almeno un istogramma è NULL\n";
        if (!h00) std::cerr << "  manca h00\n";
        if (!h01) std::cerr << "  manca h01\n";
        if (!h02) std::cerr << "  manca h02\n";
        if (!h03) std::cerr << "  manca h03\n";
        //if (!h04) std::cerr << "  manca h04\n";
        return;
    }

	std::cout << "Entries h00 = " << h00->GetEntries() << std::endl;
	std::cout << "Entries h01 = " << h01->GetEntries() << std::endl;
	std::cout << "Entries h02 = " << h02->GetEntries() << std::endl;
	std::cout << "Entries h03 = " << h03->GetEntries() << std::endl;
	//std::cout << "Entries h04 = " << h04->GetEntries() << std::endl;

    fout->cd();
    TH1F *h1Ngen = new TH1F("h1Ngen", "Ngen(Et)", noe, Ebin);

    h1Ngen->SetBinContent(1,  h00->GetBinContent(7)); // 100 GeV - 1 TeV
    h1Ngen->SetBinContent(2,  h00->GetBinContent(8));
    h1Ngen->SetBinContent(3,  h00->GetBinContent(9));
    h1Ngen->SetBinContent(4,  h00->GetBinContent(10));
    h1Ngen->SetBinContent(5,  h00->GetBinContent(11)); 
    h1Ngen->SetBinContent(6,  h00->GetBinContent(12));

	h1Ngen->SetBinContent(7,  h01->GetBinContent(13)); 
    h1Ngen->SetBinContent(8,  h01->GetBinContent(14)); 
    h1Ngen->SetBinContent(9,  h01->GetBinContent(15)); // 1 TeV - 10 TeV
    h1Ngen->SetBinContent(10, h01->GetBinContent(16));
    h1Ngen->SetBinContent(11, h01->GetBinContent(17)); 
    h1Ngen->SetBinContent(12, h01->GetBinContent(18)); 

    h1Ngen->SetBinContent(13, h02->GetBinContent(19)); 
	h1Ngen->SetBinContent(14, h02->GetBinContent(20)); 
	h1Ngen->SetBinContent(15, h02->GetBinContent(21)); // 10 TeV - 100 TeV
	h1Ngen->SetBinContent(16, h02->GetBinContent(22));
	h1Ngen->SetBinContent(17, h02->GetBinContent(23));
	h1Ngen->SetBinContent(18, h02->GetBinContent(24));   
    //h1Ngen->SetBinContent(17, h02->GetBinContent(23));    
    //h1Ngen->SetBinContent(18, h02->GetBinContent(24));

    h1Ngen->SetBinContent(19, h03->GetBinContent(25)); //100 TeV - 500 TeV
    h1Ngen->SetBinContent(20, h03->GetBinContent(26));
    h1Ngen->SetBinContent(21, h03->GetBinContent(27)); 
    h1Ngen->SetBinContent(22, h03->GetBinContent(28));
    h1Ngen->SetBinContent(23, h03->GetBinContent(29)); 
    h1Ngen->SetBinContent(24, h03->GetBinContent(30));


	//h1Ngen->SetBinContent(25, h03->GetBinContent(27));
	//h1Ngen->SetBinContent(26, h03->GetBinContent(28)); 
	//h1Ngen->SetBinContent(27, h03->GetBinContent(29)); //+ h03->GetBinContent(29)*(3050580./3097200.)));
	//h1Ngen->SetBinContent(28, h03->GetBinContent(30));
	//h1Ngen->SetBinContent(29, h03->GetBinContent(37));
	//h1Ngen->SetBinContent(30, h03->GetBinContent(38));
	//h1Ngen->SetBinContent(31, h03->GetBinContent(39));
	//h1Ngen->SetBinContent(32, h03->GetBinContent(40));

	//for(int i=1;i<=24;i++)
    	//std::cout<<"Bin "<<i<<" = "<<h1Ngen->GetBinContent(i)<<std::endl;

 /*  
    std::vector<TH1*> hists = {h00, h01, h02, h03, h04};

   int nBins = h1Ngen->GetNbinsX();
   int blockSize = nBins / hists.size();  // oppure definisci tu
 
   for (int i = 1; i <= nBins; ++i)
   {
        int histIndex = std::min((i - 1) / blockSize, (int)hists.size() - 1);
        h1Ngen->SetBinContent(i, hists[histIndex]->GetBinContent(i));
   }
  */

    double index = 1.7; // indice spettrale

	// Peso integrale bin i (come nel TCut weights[i])
	auto WeightIntegral = [&](int i){
   		double E0 = Ebin[i];
   		double E1 = Ebin[i+1];
   		return (index) * log(E1/E0) / (pow(E0,-index) - pow(E1,-index));
	};

	// Peso per evento: E^-index
	auto WeightSpectral = [&](double E){
    	return pow(E, -index);
	};
    
	// -------

    // -------		
	//TH2F *h2Ntrig_wgt = new TH2F("h2Ntrig_wgt", "Ntrig(Eo,Et) -S- (skim)", noe, Ebin, noe, Ebin);
	TH2F *h2Ntrig_cut00 =  new TH2F("h2Ntrig_cut00", "Ntrig(Eo,Et) cut00", noe, Ebin, noe, Ebin);
    TH2F *h2Ntrig_cut01 =  new TH2F("h2Ntrig_cut01", "Ntrig(Eo,Et) cut01", noe, Ebin, noe, Ebin);
    TH2F *h2Ntrig_cut05 =  new TH2F("h2Ntrig_cut05", "Ntrig(Eo,Et) cut05", noe, Ebin, noe, Ebin);
    TH2F *h2Ntrig_cut06 =  new TH2F("h2Ntrig_cut06", "Ntrig(Eo,Et) cut06", noe, Ebin, noe, Ebin);
    TH2F *h2Ntrig_stk =    new TH2F("h2Ntrig_stk", "Ntrig(Eo,Et) stk", noe, Ebin, noe, Ebin);
    TH2F *h2Ntrig_psd =    new TH2F("h2Ntrig_psd", "Ntrig(Eo,Et) psd", noe, Ebin, noe, Ebin);
    TH2F *h2Ntrig_charge = new TH2F("h2Ntrig_charge","Ntrig(Eo,Et) charge",noe, Ebin,noe, Ebin);

	TH2F *h2Ntrig_wgt = new TH2F("h2Ntrig_wgt", "Ntrig(Eo,Et) -S- (skim)", noe, Ebin, noe, Ebin);

	std::vector<double> wNgen(noe);

	for (int i = 0; i < noe; i++) {
    	double ngen = h1Ngen->GetBinContent(i+1);
    	wNgen[i] = (ngen > 0.0) ? 1.0/ngen : 0.0;
    	std::cout << "wNgen[" << i << "] = " << wNgen[i] << std::endl;
	}

	for (int i=0; i<NUM_SET; i++) {

		std::cout << "... looping on events, set " << i << " ..." << std::endl;
		
		double MC_EnergyT;
		double BGO_Energy;
		int    BGO_HET;

		double PSD_ChargeY0, PSD_ChargeY1;
		double PSD_ChargeX0, PSD_ChargeX1;

		Double_t BGO_cbgomax[14];
		Double_t BGO_cbgostk[14];

		double STKtrack_to_PSD_topY;
		double STKtrack_to_PSD_topX;

		Double_t STK_chargeY_etaCorr[6];
		Double_t STK_chargeX_etaCorr[6];

		double PSD_PW_Charge;


		skim[i]->SetBranchAddress("MC_EnergyT", &MC_EnergyT);
		skim[i]->SetBranchAddress("BGO_EnergyG_QuenchSatCorr_ML_ions", &BGO_Energy);
		skim[i]->SetBranchAddress("BGO_HET", &BGO_HET);

		skim[i]->SetBranchAddress("PSD_ChargeY0", &PSD_ChargeY0);
		skim[i]->SetBranchAddress("PSD_ChargeY1", &PSD_ChargeY1);
		skim[i]->SetBranchAddress("PSD_ChargeX0", &PSD_ChargeX0);
		skim[i]->SetBranchAddress("PSD_ChargeX1", &PSD_ChargeX1);

		skim[i]->SetBranchAddress("BGO_cbgomax", BGO_cbgomax);
		skim[i]->SetBranchAddress("BGO_cbgostk", BGO_cbgostk);

		skim[i]->SetBranchAddress("STKtrack_to_PSD_topY", &STKtrack_to_PSD_topY);
		skim[i]->SetBranchAddress("STKtrack_to_PSD_topX", &STKtrack_to_PSD_topX);

		skim[i]->SetBranchAddress("STK_chargeY_etaCorr", &STK_chargeY_etaCorr);
		skim[i]->SetBranchAddress("STK_chargeX_etaCorr", &STK_chargeX_etaCorr);

		skim[i]->SetBranchAddress("PSD_PathWeighted_Charge", &PSD_PW_Charge);


		Long64_t nentries = skim[i]->GetEntries();

		for (Long64_t ev = 0; ev < nentries; ev++) {

        skim[i]->GetEntry(ev);
        if (ev%10000==0) cout << "event " << ev << " out of " << nentries << endl;

      // -----------------------------
        // cut 00
        if (MC_EnergyT < 100.) continue;
        if (BGO_Energy < 100.) continue;
        if (BGO_HET <= 0) continue;

		// --- true-energy bin
		int j = h1Ngen->FindBin(MC_EnergyT) - 1;
		if (j < 0 || j >= noe) continue;
				
		// --- total event weight
		double w_tot = wNgen[j] * WeightIntegral(j) *WeightSpectral(MC_EnergyT);

        h2Ntrig_cut00->Fill(MC_EnergyT, BGO_Energy, w_tot);

        // cut 01
        bool hasY = (PSD_ChargeY0 > 0.0 || PSD_ChargeY1 > 0.0);
		bool hasX = (PSD_ChargeX0 > 0.0 || PSD_ChargeX1 > 0.0);
		if (!(hasY && hasX)) continue;
		h2Ntrig_cut01->Fill(MC_EnergyT, BGO_Energy, w_tot);

		// cut 05
		if (fabs(BGO_cbgomax[0] - BGO_cbgostk[0]) >= 30.0) continue;
		if (fabs(BGO_cbgomax[1] - BGO_cbgostk[1]) >= 30.0) continue;
		h2Ntrig_cut05->Fill(MC_EnergyT, BGO_Energy, w_tot);

		// cut 06
		if (fabs(STKtrack_to_PSD_topY) >= 400.) continue;
		if (fabs(STKtrack_to_PSD_topX) >= 400.) continue;
		h2Ntrig_cut06->Fill(MC_EnergyT, BGO_Energy, w_tot);

		// stk cut
		double wY = (TMath::Sign(1., STK_chargeY_etaCorr[0]) + 1.) / 2.;
		double wX = (TMath::Sign(1., STK_chargeX_etaCorr[0]) + 1.) / 2.;
		if (wY + wX <= 0.) continue;
		double stk_charge =	(wY * STK_chargeY_etaCorr[0] + wX * STK_chargeX_etaCorr[0]) / (wY + wX);
		if (stk_charge <= 1200.) continue;
		h2Ntrig_stk->Fill(MC_EnergyT, BGO_Energy, w_tot);

		// psd cut
		double psdY =
		    ((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_ChargeY0 + (TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_ChargeY1) /
    		((TMath::Sign(1.,PSD_ChargeY0)+1.)/2. + (TMath::Sign(1.,PSD_ChargeY1)+1.)/2.);
		double psdX =
    		((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.*PSD_ChargeX0 + (TMath::Sign(1.,PSD_ChargeX1)+1.)/2.*PSD_ChargeX1) /
    		((TMath::Sign(1.,PSD_ChargeX0)+1.)/2. + (TMath::Sign(1.,PSD_ChargeX1)+1.)/2.);
		if (fabs(psdY - psdX) >= 2.0) continue;
		h2Ntrig_psd->Fill(MC_EnergyT, BGO_Energy, w_tot);

		// charge cut
		double QcorrMC =((PSD_PW_Charge-MPV_MC(BGO_Energy))*(sigma_DATA(BGO_Energy)/sigma_MC(BGO_Energy)))+MPV_DATA(BGO_Energy);

		double cutLow = 2.0;
		double cutUp  = 3.0;
		double mpv_corr = MPV_MC_CORR(BGO_Energy);
		double sigma_corr = sigma_MC_CORR(BGO_Energy);

		if (QcorrMC <= mpv_corr - cutLow * sigma_corr) continue;
		if (QcorrMC >= mpv_corr + cutUp  * sigma_corr) continue;
		h2Ntrig_charge->Fill(MC_EnergyT, BGO_Energy, w_tot);

		

        // --- Fill 2D
        h2Ntrig_wgt->Fill(MC_EnergyT, BGO_Energy, w_tot);
    }
		
	}

	// --------------------------------------------------------------

	cout << " -------------------------------------------" << endl;
	cout << "       Done!  THIS IS THE END!  MIAO        " << endl;
	cout << " -------------------------------------------" << endl;

	// -------------------------------------------------------------- 	SAVE ON A FILE...

	fout->cd();
	fout->Write();
	fout->Close();

}