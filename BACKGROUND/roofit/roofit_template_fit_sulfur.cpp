#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLine.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooFormulaVar.h"
#include "RooFitResult.h"
#include "RooPlot.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

using namespace RooFit;
using namespace std;

void roofit_template_fit_sulfur(int iBin = 4, int iRebin = 2, bool writefile = false){

	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPadLeftMargin(0.09);
	gStyle->SetPadRightMargin(0.19);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadBottomMargin(0.09);

	double BGO_EnergyG_QuenchSatCorr_ML_ions2 = 0.0;
	TString title;

	// 10 bins (6 - 3 - 1)
	if(iBin == 1) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 121.1527961; title = "100 GeV < E_{BGO} < 178 GeV"  ; }
	if(iBin == 2) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 177.8277918; title = "178 GeV < E_{BGO} < 316 GeV"  ; }
	if(iBin == 3) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 261.015534 ; title = "316 GeV < E_{BGO} < 562 GeV"  ; }
	if(iBin == 4) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 383.1188748; title = "562 GeV < E_{BGO} < 1.0 TeV"  ; }
	if(iBin == 5) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 562.3413673; title = "1.0 TeV < E_{BGO} < 1.8 TeV"  ; }
	if(iBin == 6) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 825.4041434; title = "1.8 TeV < E_{BGO} < 3.2 TeV"  ; }
	if(iBin == 7) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 1467.79767 ; title = "3.2 TeV < E_{BGO} < 5.6 TeV"  ; }
	if(iBin == 8) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 3162.274615; title = "5.6 TeV < E_{BGO} < 10.0 TeV" ; }
	if(iBin == 9) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 6812.921547; title = "10.0 TeV < E_{BGO} < 31.6 TeV"; }
	if(iBin ==10) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 31622.7766 ; title = "31.6 TeV < E_{BGO} < 100.0 TeV";}

	// SULFUR 09/2026
	Double_t S_high = (15.915774 + (0.6569678*log10(BGO_EnergyG_QuenchSatCorr_ML_ions2)) - (0.167878*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions2),2)) + (0.0179829*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions2),3)));
	//(16.0082 + (0.642601*log10(BGO_EnergyG_QuenchSatCorr_ML_ions2)) - (0.181615*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions2),2)) + (0.0204998*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions2),3)));
	//(16.2307 + (0.6795615*log10(BGO_EnergyG_QuenchSatCorr_ML_ions2)) - (0.181615*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions2),2)) + (0.0204998*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions2),3)));
	Double_t S_low  = (15.1162445 + (0.53318365*log10(BGO_EnergyG_QuenchSatCorr_ML_ions2)) - (0.167878*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions2),2)) + (0.0179829*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions2),3)));
	//(15.22945 + 0.51323925*log10(BGO_EnergyG_QuenchSatCorr_ML_ions2) - 0.181615*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions2),2) + 0.0204998*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions2),3));
	//(15.00695 + (0.47627875*log10(BGO_EnergyG_QuenchSatCorr_ML_ions2)) - (0.181615*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions2),2)) + (0.0204998*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions2),3)));

	cout << "High limit: " << S_high << endl;
	cout << "Low limit: "  << S_low  << endl;

	double range_low  = 5.5;
	double range_high = 21.0;

	double range_low_temp  = 15.;
	double range_high_temp = 18.;

	double Fe_range_low = 25.2;
	double He_range_high = 3.8;

	// retrieve FD and MC histograms

	// On-orbit DATA
	TFile *f_data = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/BACKGROUND/out_root/Skim_120months_SATcorr_STKcut1200_PSDcutXY_all_12histo.root");
	std::string h_data_name;
	if (iBin < 11){ h_data_name = Form("h0%d",iBin); }
	else h_data_name = Form("h%d",iBin);
	TH1D *h_data = (TH1D *)f_data->Get( (h_data_name).c_str() );
	h_data->Rebin(iRebin);

	// MC CARBON
	TFile *f_mc_C = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/BACKGROUND/out_root/MC_FTFP_CARBON_SATcorr_STKcut1200_PSDcutXY_CORRECTED_mcfit.root");
	std::string h_mc_C_name;
	if (iBin < 11){ h_mc_C_name = Form("h0%d",iBin);}
	else h_mc_C_name = Form("h0%d",iBin);
	TH1D *h_mc_C = (TH1D *)f_mc_C->Get( (h_mc_C_name).c_str());
	h_mc_C->Rebin(iRebin);

	// MC Oxygen
	TFile *f_mc_O = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/BACKGROUND/out_root/MC_FTFP_OXYGEN_SATcorr_STKcut1200_PSDcutXY_CORRECTED_mcfit.root");
	std::string h_mc_O_name;
	if (iBin < 11){ h_mc_O_name = Form("h0%d",iBin); }
	else h_mc_O_name = Form("h%d",iBin);
	TH1D *h_mc_O = (TH1D *)f_mc_O->Get( (h_mc_O_name).c_str() );
	h_mc_O->Rebin(iRebin);

	// MC Neon
	TFile *f_mc_Ne = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/BACKGROUND/out_root/MC_FTFP_NEON_SATcorr_STKcut1200_PSDcutXY_CORRECTED_mcfit.root");
	std::string h_mc_Ne_name;
	if (iBin < 11){ h_mc_Ne_name = Form("h0%d",iBin); }
	else h_mc_Ne_name = Form("h%d",iBin);
	TH1D *h_mc_Ne = (TH1D *)f_mc_Ne->Get( (h_mc_Ne_name).c_str() );
	h_mc_Ne->Rebin(iRebin);

	// MC Sodium
	TFile *f_mc_Na = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/BACKGROUND/out_root/MC_FTFP_SODIUM_SATcorr_STKcut1200_PSDcutXY_CORRECTED_mcfit.root");
	std::string h_mc_Na_name;
	if (iBin < 11){ h_mc_Na_name = Form("h0%d",iBin); }
	else h_mc_Na_name = Form("h%d",iBin);
	TH1D *h_mc_Na = (TH1D *)f_mc_Na->Get( (h_mc_Na_name).c_str() );
	h_mc_Na->Rebin(iRebin);

	// MC Magnesium
	TFile *f_mc_Mg = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/BACKGROUND/out_root/MC_FTFP_MAGNESIUM_SATcorr_STKcut1200_PSDcutXY_CORRECTED_mcfit.root");
	std::string h_mc_Mg_name;
	if (iBin < 11){ h_mc_Mg_name = Form("h0%d",iBin); }
	else h_mc_Mg_name = Form("h%d",iBin);
	TH1D *h_mc_Mg = (TH1D *)f_mc_Mg->Get( (h_mc_Mg_name).c_str() );
	h_mc_Mg->Rebin(iRebin);

	// MC Silicon
	TFile *f_mc_Si = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/BACKGROUND/out_root/MC_FTFP_SILICON_SATcorr_STKcut1200_PSDcutXY_CORRECTED_mcfit.root");
	std::string h_mc_Si_name;
	if (iBin < 11){ h_mc_Si_name = Form("h0%d",iBin); }
	else h_mc_Si_name = Form("h%d",iBin);
	TH1D *h_mc_Si = (TH1D *)f_mc_Si->Get( (h_mc_Si_name).c_str() );
	h_mc_Si->Rebin(iRebin);

	// MC Sulfur
	TFile *f_mc_S = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/SULFUR_ANALYSIS/out_root/MC_FTFP_SULFUR_SATcorr_STKcut1200_PSDcutXY_CORRECTED_1PeV_mcfit.root");
	std::string h_mc_S_name;
	if (iBin < 11){ h_mc_S_name = Form("h0%d",iBin); }
	else h_mc_S_name = Form("h%d",iBin);
	TH1D *h_mc_S = (TH1D *)f_mc_S->Get( (h_mc_S_name).c_str() );
	h_mc_S->Rebin(iRebin);

	// MC Chlorine
	TFile *f_mc_Cl = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/secondaries/out_root/MC_FTFP_CLORO_SATcorr_STKcut1200_PSDcutXY_CORRECTED.root");
	std::string h_mc_Cl_name;
	if (iBin < 11){ h_mc_Cl_name = Form("h0%d",iBin); }
	else h_mc_Cl_name = Form("h%d",iBin);
	TH1D *h_mc_Cl = (TH1D *)f_mc_Cl->Get( (h_mc_Cl_name).c_str() );
	h_mc_Cl->Rebin(iRebin);

	// MC Argon
	TFile *f_mc_Ar = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/ARGON_ANALYSIS/out_root/MC_FTFP_ARGON_SATcorr_STKcut1200_PSDcutXY_CORRECTED_mcfit.root");
	std::string h_mc_Ar_name;
	if (iBin < 11){ h_mc_Ar_name = Form("h0%d",iBin); }
	else h_mc_Ar_name = Form("h%d",iBin);
	TH1D *h_mc_Ar = (TH1D *)f_mc_Ar->Get( (h_mc_Ar_name).c_str() );
	h_mc_Ar->Rebin(iRebin);

	// MC Potassium
	TFile *f_mc_K = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/secondaries/out_root/MC_FTFP_POTASSIUM_SATcorr_STKcut1200_PSDcutXY_CORRECTED.root");
	std::string h_mc_K_name;
	if (iBin < 11){ h_mc_K_name = Form("h0%d",iBin); }
	else h_mc_K_name = Form("h%d",iBin);
	TH1D *h_mc_K = (TH1D *)f_mc_K->Get( (h_mc_K_name).c_str() );
	h_mc_K->Rebin(iRebin);

	// MC Calcium
	TFile *f_mc_Ca = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/CALCIUM_ANALYSIS/out_root/MC_FTFP_CALCIUM_SATcorr_STKcut1200_PSDcutXY_CORRECTED_100GeV-100TeV_mcfit.root");
	std::string h_mc_Ca_name;
	if (iBin < 11){ h_mc_Ca_name = Form("h0%d",iBin); }
	else h_mc_Ca_name = Form("h%d",iBin);
	TH1D *h_mc_Ca = (TH1D *)f_mc_Ca->Get( (h_mc_Ca_name).c_str() );
	h_mc_Ca->Rebin(iRebin);

	// Create RooFit objects
	RooRealVar x("x","PSD Charge ",range_low,range_high);

	RooDataHist data("data","data",x,h_data);
	RooDataHist mcC("mcC","scaled mc C",x,h_mc_C);
	RooDataHist mcO("mcO","scaled mc O",x,h_mc_O);
	RooDataHist mcNe("mcNe","scaled mc Ne",x,h_mc_Ne);
	RooDataHist mcNa("mcNa","scaled mc Na",x,h_mc_Na);
	RooDataHist mcMg("mcMg","scaled mc Mg",x,h_mc_Mg);
	RooDataHist mcSi("mcSi","scaled mc Si",x,h_mc_Si);
	RooDataHist mcS ("mcS", "scaled mc S", x,h_mc_S );
	RooDataHist mcCl("mcCl","scaled mc Cl",x,h_mc_Cl);
	RooDataHist mcAr("mcAr","scaled mc Ar",x,h_mc_Ar);
	RooDataHist mcK("mcK","scaled mc K",x,h_mc_K);
	RooDataHist mcCa("mcCa","scaled mc Ca",x,h_mc_Ca);

	//cout << "Range h_mc_C: "  << h_mc_C->GetXaxis()->GetXmin()  << " - " << h_mc_C->GetXaxis()->GetXmax()  << endl;
	//cout << "Range h_mc_O: "  << h_mc_O->GetXaxis()->GetXmin()  << " - " << h_mc_O->GetXaxis()->GetXmax()  << endl;
	//cout << "Range h_mc_Ne: " << h_mc_Ne->GetXaxis()->GetXmin() << " - " << h_mc_Ne->GetXaxis()->GetXmax() << endl;
	//cout << "Range h_mc_Na: " << h_mc_Na->GetXaxis()->GetXmin() << " - " << h_mc_Na->GetXaxis()->GetXmax() << endl;
	//cout << "Range h_mc_Mg: " << h_mc_Mg->GetXaxis()->GetXmin() << " - " << h_mc_Mg->GetXaxis()->GetXmax() << endl;
	//cout << "Range h_mc_Si: " << h_mc_Si->GetXaxis()->GetXmin() << " - " << h_mc_Si->GetXaxis()->GetXmax() << endl;
	//cout << "Range h_mc_S:  " << h_mc_S->GetXaxis()->GetXmin()  << " - " << h_mc_S->GetXaxis()->GetXmax()  << endl;
	//cout << "Range h_mc_Ar: " << h_mc_Ar->GetXaxis()->GetXmin() << " - " << h_mc_Ar->GetXaxis()->GetXmax() << endl;
	//cout << "Range h_mc_Ca: " << h_mc_Ca->GetXaxis()->GetXmin() << " - " << h_mc_Ca->GetXaxis()->GetXmax() << endl;

	// ------------------------------------------------------------------
	
	
	double seed_C  = h_mc_C->Integral();
	double seed_O  = h_mc_O->Integral();
	double seed_Ne = h_mc_Ne->Integral();
	double seed_Na = h_mc_Na->Integral();
	double seed_Mg = h_mc_Mg->Integral();
	double seed_Si = h_mc_Si->Integral();
	double seed_S  = h_mc_S->Integral();
	double seed_Cl = h_mc_Cl->Integral();
	double seed_Ar = h_mc_Ar->Integral();
	double seed_Ca = h_mc_Ca->Integral();
	double seed_K = h_mc_K->Integral();

	double seed_tot = seed_C + seed_O + seed_Ne + seed_Na + seed_Mg + seed_Si + seed_S + seed_Cl + seed_Ar + seed_Ca + seed_K;

	double N_data = h_data->GetEntries();

	seed_C = (seed_C / seed_tot) * N_data;
	seed_O = (seed_O / seed_tot) * N_data;
	seed_Ne = (seed_Ne / seed_tot) * N_data;
	seed_Na = (seed_Na / seed_tot) * N_data;
	seed_Mg = (seed_Mg / seed_tot) * N_data;
	seed_Si = (seed_Si / seed_tot) * N_data;
	seed_S  = (seed_S  / seed_tot) * N_data;
	seed_Cl = (seed_Cl / seed_tot) * N_data;
	seed_Ar = (seed_Ar / seed_tot) * N_data;
	seed_Ca = (seed_Ca / seed_tot) * N_data;
	seed_K = (seed_K / seed_tot) * N_data;

	//cout << "seed_C=" << seed_C << " seed_O=" << seed_O << " seed_Ne=" << seed_Ne << /*" seed_Na=" << seed_Na
	//     <<*/ " seed_Mg=" << seed_Mg << " seed_Si=" << seed_Si << " seed_S=" << seed_S << " seed_Cl=" << seed_Cl << " seed_Ar=" << seed_Ar << " seed_Ca=" << seed_Ca << endl;

	RooRealVar mcC_frac("mcC_frac","Fraction of mcC", seed_C, 0.0, 10.0 * h_data->GetEntries() );
	RooRealVar mcO_frac("mcO_frac","Fraction of mcO", seed_O, 0.0, 10.0 * h_data->GetEntries() );
	RooRealVar mcNe_frac("mcNe_frac","Fraction of mcNe", seed_Ne, 0.0, 10.0 * h_data->GetEntries() );
	RooRealVar mcNa_frac("mcNa_frac","Fraction of mcNa", seed_Na, 0.0, 10.0 * h_data->GetEntries() );
	RooRealVar mcMg_frac("mcMg_frac","Fraction of mcMg", seed_Mg, 0.0, 10.0 * h_data->GetEntries() );
	RooRealVar mcSi_frac("mcSi_frac","Fraction of mcSi", seed_Si, 0.0, 10.0 * h_data->GetEntries() );
	RooRealVar mcS_frac( "mcS_frac", "Fraction of mcS",  seed_S,  0.0, 10.0 * h_data->GetEntries() );
	RooRealVar mcCl_frac( "mcCl_frac", "Fraction of mcCl", seed_Cl,  0.0, 10.0 * h_data->GetEntries() );
	RooRealVar mcAr_frac( "mcAr_frac", "Fraction of mcAr", seed_Ar,  0.0, 10.0 * h_data->GetEntries() );
	RooRealVar mcCa_frac( "mcCa_frac", "Fraction of mcCa", seed_Ca,  0.0, 10.0 * h_data->GetEntries() );
	RooRealVar mcK_frac( "mcK_frac", "Fraction of mcK", seed_K,  0.0, 10.0 * h_data->GetEntries() );

	//cout << "Integral C = " << h_mc_C->Integral() << endl;
	//cout << "Integral O = " << h_mc_O->Integral() << endl;
	//cout << "Integral Ne = " << h_mc_Ne->Integral() << endl;
	//cout << "Integral Na = " << h_mc_Na->Integral() << endl;
	//cout << "Integral Mg = " << h_mc_Mg->Integral() << endl;
	//cout << "Integral Si = " << h_mc_Si->Integral() << endl;
	//cout << "Integral S  = " << h_mc_S->Integral() << endl;
	//cout << "Integral Ar = " << h_mc_Ar->Integral() << endl;
	//cout << "Integral Ca = " << h_mc_Ca->Integral() << endl;
	
	//cout << "Data integral = " << h_data->Integral() << endl;
	//cout << "Data entries  = " << h_data->GetEntries() << endl;
    
 
	// Create MC template model
	RooHistPdf modelmcC("modelmcC","modelmcC",x,mcC);
	RooHistPdf modelmcO("modelmcO","modelmcO",x, mcO);
	RooHistPdf modelmcNe("modelmcNe","modelmcNe",x, mcNe);
	RooHistPdf modelmcNa("modelmcNa","modelmcNa",x, mcNa);
	RooHistPdf modelmcMg("modelmcMg","modelmcMg",x, mcMg);
	RooHistPdf modelmcSi("modelmcSi","modelmcSi",x, mcSi);
	RooHistPdf modelmcS( "modelmcS", "modelmcS", x, mcS);
	RooHistPdf modelmcCl( "modelmcCl", "modelmcCl", x, mcCl);
	RooHistPdf modelmcAr( "modelmcAr", "modelmcAr", x, mcAr);
	RooHistPdf modelmcCa( "modelmcCa", "modelmcCa", x, mcCa);
	RooHistPdf modelmcK( "modelmcK", "modelmcK", x, mcK);
	RooAddPdf model("model","model",
		RooArgList(modelmcC, modelmcO, modelmcNe,modelmcNa,modelmcMg,modelmcSi,modelmcS,modelmcCl,modelmcAr,modelmcCa,modelmcK),
		RooArgList(mcC_frac, mcO_frac, mcNe_frac,mcNa_frac,mcMg_frac,mcSi_frac,mcS_frac,mcCl_frac,mcAr_frac,mcCa_frac,mcK_frac));

	// Define side band regions and full range
	x.setRange("Sulfur",  range_low_temp, range_high_temp);
	x.setRange("FULL",  range_low, range_high);
	x.setRange("Sulfur_2", S_low, S_high);

	// Fit template to data
	//model.fitTo(data, RooFit::Range("Iron,Neon"), RooFit::PrintLevel(0), RooFit::SumW2Error(false));
	//RooFitResult* res = model.fitTo(data, RooFit::Range("Helium,Iron,Sulfur"), RooFit::PrintLevel(0), RooFit::SumW2Error(false));
    RooFitResult* fitres = model.fitTo(data, RooFit::Range("FULL"), RooFit::Save(), RooFit::PrintLevel(0));

	// Integrals of each component in the purity region "Sulfur_2"
	RooAbsReal* I_C = modelmcC.createIntegral(x, NormSet(x), Range("Sulfur_2"));
	RooAbsReal* I_O = modelmcO.createIntegral(x, NormSet(x), Range("Sulfur_2"));
	RooAbsReal* I_Ne = modelmcNe.createIntegral(x, NormSet(x), Range("Sulfur_2"));
	RooAbsReal* I_Na = modelmcNa.createIntegral(x, NormSet(x), Range("Sulfur_2"));
	RooAbsReal* I_Mg = modelmcMg.createIntegral(x, NormSet(x), Range("Sulfur_2"));
	RooAbsReal* I_Si = modelmcSi.createIntegral(x, NormSet(x), Range("Sulfur_2"));
	RooAbsReal* I_S  = modelmcS.createIntegral(x, NormSet(x), Range("Sulfur_2"));
	RooAbsReal* I_Cl = modelmcCl.createIntegral(x, NormSet(x), Range("Sulfur_2"));
	RooAbsReal* I_Ar = modelmcAr.createIntegral(x, NormSet(x), Range("Sulfur_2"));
	RooAbsReal* I_Ca = modelmcCa.createIntegral(x, NormSet(x), Range("Sulfur_2"));
	RooAbsReal* I_K = modelmcK.createIntegral(x, NormSet(x), Range("Sulfur_2"));
/*
	std::vector<RooRealVar*> yields    = { &mcNe_frac, &mcMg_frac, &mcSi_frac, &mcS_frac, &mcAr_frac };
	std::vector<RooAbsReal*> integrals = { I_Ne, I_Mg, I_Si, I_S, I_Ar };
	std::vector<std::string> names     = { "Ne", "Mg", "Si", "S", "Ar" };
*/

	std::vector<RooRealVar*> yields = {
		&mcC_frac, &mcO_frac, &mcNe_frac, &mcNa_frac, &mcMg_frac, &mcSi_frac, &mcS_frac, &mcCl_frac, &mcAr_frac, &mcCa_frac, &mcK_frac
	};

	std::vector<RooAbsReal*> integrals = { I_C, I_O, I_Ne, I_Na, I_Mg, I_Si, I_S, I_Cl, I_Ar, I_Ca, I_K };
	std::vector<std::string> names = { "C", "O", "Ne", "Na", "Mg", "Si", "S", "Cl", "Ar", "Ca", "K" };

	RooArgList denomArgs;
	std::string denomExpr;

	for (size_t i = 0; i < yields.size(); ++i) {
		denomExpr += Form("@%zu*@%zu", 2*i, 2*i+1);
		if (i < yields.size()-1) denomExpr += " + ";
		denomArgs.add(*yields[i]);
		denomArgs.add(*integrals[i]);
	}

	RooFormulaVar Ntot_S("Ntot_S", denomExpr.c_str(), denomArgs);

	for (size_t i = 0; i < yields.size(); ++i) {

		if (names[i] == "S") continue;

		RooFormulaVar cont_i(
			Form("cont_%s", names[i].c_str()),
			"(@0*@1)/@2",
			RooArgList(*yields[i], *integrals[i], Ntot_S)
		);

		double val = cont_i.getVal();
		double err = cont_i.getPropagatedError(*fitres);

		std::cout << "Contamination of " << names[i]
				  << " in S = "
				  << val << " ± " << err << std::endl;
	}

	RooFormulaVar purity_S(
		"purity_S",
		"(@0*@1)/@2",
		RooArgList(mcS_frac, *I_S, Ntot_S)
	);

	double purity = purity_S.getVal();
	double purity_err = purity_S.getPropagatedError(*fitres);

	std::cout << "Sum of contaminations = " << (1.0 - purity)
			  << " ± " << purity_err << std::endl;


	//cout << "\n=== POST FIT ===" << endl;
	//cout << "C = " << mcC_frac.getVal() << endl;
	//cout << "O = " << mcO_frac.getVal() << endl;
	//cout << "Ne = " << mcNe_frac.getVal() << endl;
	//cout << "Na = " << mcNa_frac.getVal() << endl;
	//cout << "Mg = " << mcMg_frac.getVal() << endl;
	//cout << "Si = " << mcSi_frac.getVal() << endl;
	//cout << "S  = " << mcS_frac.getVal() << endl;
	//cout << "Ar = " << mcAr_frac.getVal() << endl;
	//cout << "Ca = " << mcCa_frac.getVal() << endl;



	fitres->Print("v");

	cout << "status  = " << fitres->status() << endl;
	cout << "covQual = " << fitres->covQual() << endl;
	cout << "EDM     = " << fitres->edm() << endl;
	
	// -------------------------- PLOT ----------------------------------
	TCanvas* c1 = new TCanvas("c1", "c1", 1300, 980);
	RooPlot* dframe = x.frame(range_low,range_high,100);
	dframe->SetTitle(" ");

	TString colCode[11] = {"#6a4c93", "#006627", "#0bb4ff", "#50e991", "#e6d800", "#81008d", "#ff6b6b", "#a775f8", "#ffa500", "#ff0000", "#c0c0c0" };
	Int_t col[11];
	for(int ii=0; ii<11; ii++){ col[ii]=TColor::GetColor(colCode[ii]); }

	//Color_t col[9];
	//for (int i = 0; i < 9; i++) {
	//    double x = i / 8.0;
	//    int gray = 30 + x * (210 - 30);
	//    col[i] = TColor::GetColor(gray, gray, gray);
	//}
	

	c1->SetLogy();
	data.plotOn(dframe);
	model.plotOn(dframe, Name("modelmc_postfit"), LineColor(kRed), LineWidth(4));
	model.plotOn(dframe, Components(modelmcC),  Name("modelmcC_postfit"),  LineStyle(1), LineColor(col[0]), FillColor(col[0]), FillStyle(3003), DrawOption("FL"), RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcO),  Name("modelmcO_postfit"),  LineStyle(1), LineColor(col[1]), FillColor(col[1]), FillStyle(3003), DrawOption("FL"), RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcNe), Name("modelmcNe_postfit"), LineStyle(1), LineColor(col[2]), FillColor(col[2]), FillStyle(3003), DrawOption("FL"), RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcNa), Name("modelmcNa_postfit"), LineStyle(1), LineColor(col[3]), FillColor(col[3]), FillStyle(3002), DrawOption("FL"), RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcMg), Name("modelmcMg_postfit"), LineStyle(1), LineColor(col[4]), FillColor(col[4]), FillStyle(3003), DrawOption("FL"), RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcSi), Name("modelmcSi_postfit"), LineStyle(1), LineColor(col[5]), FillColor(col[5]), FillStyle(3003), DrawOption("FL"), RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcS),  Name("modelmcS_postfit"),  LineStyle(1), LineColor(col[6]), FillColor(col[6]), FillStyle(3001), DrawOption("FL"), RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcCl), Name("modelmcCl_postfit"), LineStyle(1), LineColor(col[9]), FillColor(col[9]), FillStyle(3003), DrawOption("FL"), RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcAr), Name("modelmcAr_postfit"), LineStyle(1), LineColor(col[7]), FillColor(col[7]), FillStyle(3003), DrawOption("FL"), RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcCa), Name("modelmcCa_postfit"), LineStyle(1), LineColor(col[8]), FillColor(col[8]), FillStyle(3003), DrawOption("FL"), RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcK), Name("modelmcK_postfit"), LineStyle(1), LineColor(col[10]), FillColor(col[10]), FillStyle(3003), DrawOption("FL"), RooFit::VLines(), MoveToBack());
	dframe->Draw();

	TLine *lmax = new TLine(S_high,dframe->GetMinimum(),S_high,dframe->GetMaximum());
	lmax->SetLineColor(kGray+4);
	lmax->SetLineStyle(9);
	lmax->SetLineWidth(2);
	lmax->Draw();

	TLine *lmin = new TLine(S_low,dframe->GetMinimum(),S_low,dframe->GetMaximum());
	lmin->SetLineColor(kGray+4);
	lmin->SetLineStyle(9);
	lmin->SetLineWidth(2);
	lmin->Draw();

	// Draw Legend
	TLegend *l = new TLegend(0.82,0.33,0.99,0.9);
	l->AddEntry("modelmcC_postfit","Carbon","f");
	l->AddEntry("modelmcO_postfit","Oxygen","f");
	l->AddEntry("modelmcNe_postfit","Neon","f");
	l->AddEntry("modelmcNa_postfit","Sodium","f");
	l->AddEntry("modelmcMg_postfit","Magnesium","f");
	l->AddEntry("modelmcSi_postfit","Silicon","f");
	l->AddEntry("modelmcS_postfit","Sulfur","f");
	l->AddEntry("modelmcCl_postfit","Chlorine","f");
	l->AddEntry("modelmcAr_postfit","Argon","f");
	l->AddEntry("modelmcK_postfit","Potassium","f");
	l->AddEntry("modelmcCa_postfit","Calcium","f");
	l->AddEntry("modelmc_postfit","MC total ","l");
	l->AddEntry("data","Flight data","epl");
	l->SetLineColor(0);
	l->SetTextSize(0.035);
	l->Draw("same");

	TLatex latex; latex.SetTextSize(0.04); latex.SetNDC(true); latex.SetTextFont(42);
	latex.DrawLatex(0.5, 0.88, title);

	if (writefile) {

    // =========================
    //  FILE TXT
    // =========================

    TString outfileName = "contaminazioni_sulfur_10bins_ago26.txt";

    std::ofstream outFile(outfileName.Data(), std::ios::app);

    if (!outFile.is_open()) {
        std::cerr << "ERROR: unable to open "
                  << outfileName << std::endl;
    }
    else {

        outFile << std::fixed << std::setprecision(6);

        outFile << iBin << " "
                << BGO_EnergyG_QuenchSatCorr_ML_ions2 << " ";

        for (size_t i = 0; i < yields.size(); ++i) {

            if (names[i] == "S")
                continue;

            RooFormulaVar cont_i(
                Form("cont_%s", names[i].c_str()),
                "(@0*@1)/@2",
                RooArgList(*yields[i], *integrals[i], Ntot_S)
            );

            double val = cont_i.getVal();
            double err = cont_i.getPropagatedError(*fitres);

            outFile << val << " " << err << " ";
        }

        outFile << std::endl;

        outFile.close();

        std::cout << "TXT saved: "
                  << outfileName << std::endl;
    }


    // =========================
    //  FILE PNG
    // =========================

    TString plot_name1 =
        "PLOTS/template_sulfur_10bins_ago26_" +
        h_data_name + ".png";

    std::cout << "Trying to save: "
              << plot_name1 << std::endl;

    c1->SaveAs(plot_name1.Data());

    std::cout << "PNG saved." << std::endl;
}
	}
