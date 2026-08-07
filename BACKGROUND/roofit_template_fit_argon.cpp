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

void roofit_template_fit_argon(int iBin = 4, int iRebin = 2, bool writefile = false){

	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPadLeftMargin(0.09);
	gStyle->SetPadRightMargin(0.19);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadBottomMargin(0.09);

	double BGO_EnergyG_QuenchSatCorr_ML_ions2 = 0.0;
	TString title;

	// 10 bins (6 - 3 - 1)
	if(iBin == 1) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 121.1527961; title = "100 < E_{dep}/GeV < 147"; }
	if(iBin == 2) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 177.8277918; title = "147 < E_{dep}/GeV < 215"; }
	if(iBin == 3) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 261.015534 ; title = "215 < E_{dep}/GeV < 316"; }
	if(iBin == 4) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 383.1188748; title = "316 < E_{dep}/GeV < 464"; }
	if(iBin == 5) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 562.3413673; title = "464 < E_{dep}/GeV < 681"; }
	if(iBin == 6) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 825.4041434; title = "0.7 < E_{dep}/TeV < 1.0"; }
	if(iBin == 7) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 1467.79767 ; title = "1.0 < E_{dep}/TeV < 2.1"; }
	if(iBin == 8) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 3162.274615; title = "2.1 < E_{dep}/TeV < 4.6"; }
	if(iBin == 9) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 6812.921547; title = "4.6 < E_{dep}/TeV < 10";  }
	if(iBin ==10) { BGO_EnergyG_QuenchSatCorr_ML_ions2 = 31622.7766 ; title = " 10 < E_{dep}/TeV < 100"; }

	// ARGON 08/2026
	Double_t Ar_high = (18.3396475+(0.16117555*log10(BGO_EnergyG_QuenchSatCorr_ML_ions2)));
	Double_t Ar_low  = (17.684235+(0.0208433*log10(BGO_EnergyG_QuenchSatCorr_ML_ions2)));

	cout << "High limit: " << Ar_high << endl;
	cout << "Low limit: "  << Ar_low  << endl;

	double range_low  = 7.5;
	double range_high = 18.9;

	double range_low_temp  = 15.;
	double range_high_temp = 18.;

	double Fe_range_low = 25.2;
	double He_range_high = 3.8;

	// retrieve FD and MC histograms

	// On-orbit DATA
	TFile *f_data = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/BACKGROUND/out_root/Skim_120months_SATcorr_STKcut1200_PSDcutXY_10bin_v2.root");
	std::string h_data_name;
	if (iBin < 10){ h_data_name = Form("h0%d",iBin); }
	else h_data_name = Form("h%d",iBin);
	TH1D *h_data = (TH1D *)f_data->Get( (h_data_name).c_str() );
	h_data->Rebin(iRebin);

	// MC Oxygen
	TFile *f_mc_O = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/BACKGROUND/out_root/MC_FTFP_OXYGEN_SATcorr_STKcut1200_PSDcutXY_CORRECTED_mcfit.root");
	std::string h_mc_O_name;
	if (iBin < 10){ h_mc_O_name = Form("h0%d",iBin); }
	else h_mc_O_name = Form("h%d",iBin);
	TH1D *h_mc_O = (TH1D *)f_mc_O->Get( (h_mc_O_name).c_str() );
	h_mc_O->Rebin(iRebin);

	// MC Neon
	TFile *f_mc_Ne = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/BACKGROUND/out_root/MC_FTFP_NEON_SATcorr_STKcut1200_PSDcutXY_CORRECTED_mcfit.root");
	std::string h_mc_Ne_name;
	if (iBin < 10){ h_mc_Ne_name = Form("h0%d",iBin); }
	else h_mc_Ne_name = Form("h%d",iBin);
	TH1D *h_mc_Ne = (TH1D *)f_mc_Ne->Get( (h_mc_Ne_name).c_str() );
	h_mc_Ne->Rebin(iRebin);

	// MC Magnesium
	TFile *f_mc_Mg = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/BACKGROUND/out_root/MC_FTFP_MAGNESIUM_SATcorr_STKcut1200_PSDcutXY_CORRECTED_mcfit.root");
	std::string h_mc_Mg_name;
	if (iBin < 10){ h_mc_Mg_name = Form("h0%d",iBin); }
	else h_mc_Mg_name = Form("h%d",iBin);
	TH1D *h_mc_Mg = (TH1D *)f_mc_Mg->Get( (h_mc_Mg_name).c_str() );
	h_mc_Mg->Rebin(iRebin);

	// MC Silicon
	TFile *f_mc_Si = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/BACKGROUND/out_root/MC_FTFP_SILICON_SATcorr_STKcut1200_PSDcutXY_CORRECTED_mcfit.root");
	std::string h_mc_Si_name;
	if (iBin < 10){ h_mc_Si_name = Form("h0%d",iBin); }
	else h_mc_Si_name = Form("h%d",iBin);
	TH1D *h_mc_Si = (TH1D *)f_mc_Si->Get( (h_mc_Si_name).c_str() );
	h_mc_Si->Rebin(iRebin);

	// MC Sulfur
	TFile *f_mc_S = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/SULFUR_ANALYSIS/out_root/MC_FTFP_SULFUR_SATcorr_STKcut1200_PSDcutXY_CORRECTED_100GeV-500TeV_mcfit.root");
	std::string h_mc_S_name;
	if (iBin < 10){ h_mc_S_name = Form("h0%d",iBin); }
	else h_mc_S_name = Form("h%d",iBin);
	TH1D *h_mc_S = (TH1D *)f_mc_S->Get( (h_mc_S_name).c_str() );
	h_mc_S->Rebin(iRebin);

	// MC Argon
		TFile *f_mc_Ar = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/ARGON_ANALYSIS/out_root/MC_FTFP_ARGON_SATcorr_STKcut1200_PSDcutXY_CORRECTED_mcfit.root");
	std::string h_mc_Ar_name;
	if (iBin < 10){ h_mc_Ar_name = Form("h0%d",iBin); }
	else h_mc_Ar_name = Form("h%d",iBin);
	TH1D *h_mc_Ar = (TH1D *)f_mc_Ar->Get( (h_mc_Ar_name).c_str() );
	h_mc_Ar->Rebin(iRebin);

	// Create RooFit objects
	RooRealVar x("x","PSD Charge ",range_low,range_high);

	RooDataHist data("data","data",x,h_data);
	RooDataHist mcO("mcO","scaled mc O",x,h_mc_O);
	RooDataHist mcNe("mcNe","scaled mc Ne",x,h_mc_Ne);
	RooDataHist mcMg("mcMg","scaled mc Mg",x,h_mc_Mg);
	RooDataHist mcSi("mcSi","scaled mc Si",x,h_mc_Si);
	RooDataHist mcS ("mcS", "scaled mc S", x,h_mc_S );
	RooDataHist mcAr("mcAr","scaled mc Ar",x,h_mc_Ar);

	cout << "Range h_mc_O: "  << h_mc_O->GetXaxis()->GetXmin()  << " - " << h_mc_O->GetXaxis()->GetXmax()  << endl;
	cout << "Range h_mc_Ne: " << h_mc_Ne->GetXaxis()->GetXmin() << " - " << h_mc_Ne->GetXaxis()->GetXmax() << endl;
	cout << "Range h_mc_Mg: " << h_mc_Mg->GetXaxis()->GetXmin() << " - " << h_mc_Mg->GetXaxis()->GetXmax() << endl;
	cout << "Range h_mc_Si: " << h_mc_Si->GetXaxis()->GetXmin() << " - " << h_mc_Si->GetXaxis()->GetXmax() << endl;
	cout << "Range h_mc_S:  " << h_mc_S->GetXaxis()->GetXmin()  << " - " << h_mc_S->GetXaxis()->GetXmax()  << endl;
	cout << "Range h_mc_Ar: " << h_mc_Ar->GetXaxis()->GetXmin() << " - " << h_mc_Ar->GetXaxis()->GetXmax() << endl;

	// ------------------------------------------------------------------
	// Seed dinamico = integrale dell'istogramma MC. Sostituisce la vecchia
	// tabella hardcoded frMC[11] (pensata per un fit a 11 elementi He->Fe,
	// non coerente con questo fit a 4 componenti Ne/Mg/Si/S).
	// Richiede che gli MC siano già normalizzati alla stessa esposizione
	// dei dati di volo.
	// ------------------------------------------------------------------
	
	double seed_O  = h_mc_O->Integral();
	double seed_Ne = h_mc_Ne->Integral();
	double seed_Mg = h_mc_Mg->Integral();
	double seed_Si = h_mc_Si->Integral();
	double seed_S  = h_mc_S->Integral();
	double seed_Ar = h_mc_Ar->Integral();

	double seed_tot = seed_O + seed_Ne + seed_Mg + seed_Si + seed_S + seed_Ar;

	double N_data = h_data->GetEntries();

	seed_O = (seed_O / seed_tot) * N_data;
	seed_Ne = (seed_Ne / seed_tot) * N_data;
	seed_Mg = (seed_Mg / seed_tot) * N_data;
	seed_Si = (seed_Si / seed_tot) * N_data;
	seed_S  = (seed_S  / seed_tot) * N_data;
	seed_Ar = (seed_Ar / seed_tot) * N_data;

	cout << "seed_O=" << seed_O << " seed_Ne=" << seed_Ne << " seed_Mg=" << seed_Mg
	     << " seed_Si=" << seed_Si << " seed_S=" << seed_S << " seed_Ar=" << seed_Ar << endl;

	RooRealVar mcO_frac("mcO_frac","Fraction of mcO", seed_O, 0.0, 10.0 * h_data->GetEntries() );
	RooRealVar mcNe_frac("mcNe_frac","Fraction of mcNe", seed_Ne, 0.0, 10.0 * h_data->GetEntries() );
	RooRealVar mcMg_frac("mcMg_frac","Fraction of mcMg", seed_Mg, 0.0, 10.0 * h_data->GetEntries() );
	RooRealVar mcSi_frac("mcSi_frac","Fraction of mcSi", seed_Si, 0.0, 10.0 * h_data->GetEntries() );
	RooRealVar mcS_frac( "mcS_frac", "Fraction of mcS",  seed_S,  0.0, 10.0 * h_data->GetEntries() );
	RooRealVar mcAr_frac( "mcAr_frac", "Fraction of mcAr", seed_Ar,  0.0, 10.0 * h_data->GetEntries() );

	cout << "Integral O = " << h_mc_O->Integral() << endl;
	cout << "Integral Ne = " << h_mc_Ne->Integral() << endl;
	cout << "Integral Mg = " << h_mc_Mg->Integral() << endl;
	cout << "Integral Si = " << h_mc_Si->Integral() << endl;
	cout << "Integral S  = " << h_mc_S->Integral() << endl;
	cout << "Integral Ar = " << h_mc_Ar->Integral() << endl;
	
	cout << "Data integral = " << h_data->Integral() << endl;
	cout << "Data entries  = " << h_data->GetEntries() << endl;
    
 /*
	double frMC[5];
	if(iBin == 1) {
		double tmp[5] = { 0.24, 0.33, 0.32, 0.07, 0.01 };
		for(int ii=0; ii<5; ii++) frMC[ii] = tmp[ii];
	}
	if(iBin == 2) {
		double tmp[5] = { 0.23, 0.33, 0.32, 0.07, 0.03 };
		for(int ii=0; ii<5; ii++) frMC[ii] = tmp[ii];
	}
	if(iBin == 3) {
		double tmp[5] = { 0.23, 0.33, 0.32, 0.07, 0.03 }; 
		for(int ii=0; ii<5; ii++) frMC[ii] = tmp[ii];
	}
	if(iBin == 4) {
		double tmp[5] = { 0.23, 0.33, 0.32, 0.07, 0.03 }; 
		for(int ii=0; ii<5; ii++) frMC[ii] = tmp[ii];
	}
	if(iBin == 5) {
		double tmp[5] = { 0.23, 0.33, 0.32, 0.07, 0.03 }; 
		for(int ii=0; ii<5; ii++) frMC[ii] = tmp[ii];
	}
	if(iBin == 6) {
		double tmp[5] = { 0.23, 0.33, 0.32, 0.07, 0.03 }; 
		for(int ii=0; ii<5; ii++) frMC[ii] = tmp[ii];
	}
	if(iBin == 7) {
		double tmp[5] = {0.23, 0.33, 0.32, 0.07, 0.03 };  
		for(int ii=0; ii<5; ii++) frMC[ii] = tmp[ii];
	}
	if(iBin == 8) {
		double tmp[5] = {0.23, 0.33, 0.32, 0.07, 0.03 }; 
		for(int ii=0; ii<5; ii++) frMC[ii] = tmp[ii];
	}
	if(iBin == 9) {
		double tmp[5] = { 0.23, 0.33, 0.32, 0.07, 0.03 }; 
		for(int ii=0; ii<5; ii++) frMC[ii] = tmp[ii];
	}
	if(iBin == 10) {
		double tmp[5] = { 0.23, 0.33, 0.32, 0.07, 0.03 }; 
		for(int ii=0; ii<5; ii++) frMC[ii] = tmp[ii];
	}

	RooRealVar mcNe_frac("mcNe_frac","Fraction of mcNe",frMC[0]*( h_data->GetEntries()), 0.0, 10.0 * h_data->GetEntries() );	// These are variables for output
	RooRealVar mcMg_frac("mcMg_frac","Fraction of mcMg",frMC[1]*( h_data->GetEntries()), 0.0, 10.0 * h_data->GetEntries() );	// These are variables for output
	RooRealVar mcSi_frac("mcSi_frac","Fraction of mcSi",frMC[2]*( h_data->GetEntries()), 0.0, 10.0 * h_data->GetEntries() );	// These are variables for output
	RooRealVar mcS_frac("mcS_frac","Fraction of mcS",frMC[3]*(h_data->GetEntries()), 0.0, 10.0 * h_data->GetEntries() );	// These are variables for output
	RooRealVar mcAr_frac("mcAr_frac","Fraction of mcAr",frMC[4]*(h_data->GetEntries()), 0.0, 10.0 * h_data->GetEntries() );	// These are variables for output
*/
	// Create MC template model
	RooHistPdf modelmcO("modelmcO","modelmcO",x, mcO);
	RooHistPdf modelmcNe("modelmcNe","modelmcNe",x, mcNe);
	RooHistPdf modelmcMg("modelmcMg","modelmcMg",x, mcMg);
	RooHistPdf modelmcSi("modelmcSi","modelmcSi",x, mcSi);
	RooHistPdf modelmcS( "modelmcS", "modelmcS", x, mcS);
	RooHistPdf modelmcAr( "modelmcAr", "modelmcAr", x, mcAr);
	RooAddPdf model("model","model",
		RooArgList(modelmcO, modelmcNe,modelmcMg,modelmcSi,modelmcS,modelmcAr),
		RooArgList(mcO_frac, mcNe_frac,mcMg_frac,mcSi_frac,mcS_frac,mcAr_frac));

	// Define side band regions and full range
	x.setRange("Argon",  range_low_temp, range_high_temp);
	x.setRange("FULL",  range_low, range_high);
	x.setRange("Argon_2", Ar_low, Ar_high);

	// Fit template to data
	//model.fitTo(data, RooFit::Range("Iron,Neon"), RooFit::PrintLevel(0), RooFit::SumW2Error(false));
	//RooFitResult* res = model.fitTo(data, RooFit::Range("Helium,Iron,Sulfur"), RooFit::PrintLevel(0), RooFit::SumW2Error(false));
    RooFitResult* fitres = model.fitTo(data, RooFit::Range("FULL"), RooFit::Save(), RooFit::PrintLevel(0));

	// Integrals of each component in the purity region "Sulfur_2"
	RooAbsReal* I_O = modelmcO.createIntegral(x, NormSet(x), Range("Argon_2"));
	RooAbsReal* I_Ne = modelmcNe.createIntegral(x, NormSet(x), Range("Argon_2"));
	RooAbsReal* I_Mg = modelmcMg.createIntegral(x, NormSet(x), Range("Argon_2"));
	RooAbsReal* I_Si = modelmcSi.createIntegral(x, NormSet(x), Range("Argon_2"));
	RooAbsReal* I_S  = modelmcS.createIntegral(x, NormSet(x), Range( "Argon_2"));
	RooAbsReal* I_Ar = modelmcAr.createIntegral(x, NormSet(x), Range("Argon_2"));
/*
	std::vector<RooRealVar*> yields    = { &mcNe_frac, &mcMg_frac, &mcSi_frac, &mcS_frac, &mcAr_frac };
	std::vector<RooAbsReal*> integrals = { I_Ne, I_Mg, I_Si, I_S, I_Ar };
	std::vector<std::string> names     = { "Ne", "Mg", "Si", "S", "Ar" };
*/

	std::vector<RooRealVar*> yields = {
		&mcO_frac, &mcNe_frac,  &mcMg_frac,  &mcSi_frac, &mcS_frac, &mcAr_frac
	};

	std::vector<RooAbsReal*> integrals = { I_O, I_Ne,  I_Mg,  I_Si, I_S, I_Ar };
	std::vector<std::string> names = { "O", "Ne","Mg","Si","S", "Ar" };

	RooArgList denomArgs;
	std::string denomExpr;

	for (size_t i = 0; i < yields.size(); ++i) {
		denomExpr += Form("@%zu*@%zu", 2*i, 2*i+1);
		if (i < yields.size()-1) denomExpr += " + ";
		denomArgs.add(*yields[i]);
		denomArgs.add(*integrals[i]);
	}

	RooFormulaVar Ntot_Ar("Ntot_Ar", denomExpr.c_str(), denomArgs);

	for (size_t i = 0; i < yields.size(); ++i) {

		if (names[i] == "Ar") continue;

		RooFormulaVar cont_i(
			Form("cont_%s", names[i].c_str()),
			"(@0*@1)/@2",
			RooArgList(*yields[i], *integrals[i], Ntot_Ar)
		);

		double val = cont_i.getVal();
		double err = cont_i.getPropagatedError(*fitres);

		std::cout << "Contaminazione " << names[i]
				  << " nello Argon = "
				  << val << " ± " << err << std::endl;
	}

	RooFormulaVar purity_Ar(
		"purity_Ar",
		"(@0*@1)/@2",
		RooArgList(mcAr_frac, *I_Ar, Ntot_Ar)
	);

	double purity = purity_Ar.getVal();
	double purity_err = purity_Ar.getPropagatedError(*fitres);

	std::cout << "Somma contaminazioni = " << (1.0 - purity)
			  << " ± " << purity_err << std::endl;


	cout << "\n=== POST FIT ===" << endl;
	cout << "O = " << mcO_frac.getVal() << endl;
	cout << "Ne = " << mcNe_frac.getVal() << endl;
	cout << "Mg = " << mcMg_frac.getVal() << endl;
	cout << "Si = " << mcSi_frac.getVal() << endl;
	cout << "S  = " << mcS_frac.getVal() << endl;
	cout << "Ar = " << mcAr_frac.getVal() << endl;


	fitres->Print("v");

	cout << "status  = " << fitres->status() << endl;
	cout << "covQual = " << fitres->covQual() << endl;
	cout << "EDM     = " << fitres->edm() << endl;
	
	// -------------------------- PLOT ----------------------------------
	TCanvas* c1 = new TCanvas("c1", "c1", 1300, 980);
	RooPlot* dframe = x.frame(range_low,range_high,100);
	dframe->SetTitle(" ");

	TString colCode[6] = {"#6a4c93", "#b30000", "#0bb4ff", "#50e991", "#e6d800", "#81008d" };
	Int_t col[6];
	for(int ii=0; ii<6; ii++){ col[ii]=TColor::GetColor(colCode[ii]); }

	c1->SetLogy();
	data.plotOn(dframe);
	model.plotOn(dframe, Name("modelmc_postfit"), LineColor(kRed), LineWidth(4));
	model.plotOn(dframe, Components(modelmcO),  Name("modelmcO_postfit"),  LineStyle(1), LineColor(col[5]), FillColor(col[5]), FillStyle(3002), DrawOption("FL"), RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcNe), Name("modelmcNe_postfit"), LineStyle(1), LineColor(col[0]), FillColor(col[0]), FillStyle(3002), DrawOption("FL"), RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcMg), Name("modelmcMg_postfit"), LineStyle(1), LineColor(col[1]), FillColor(col[1]), FillStyle(3002), DrawOption("FL"), RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcSi), Name("modelmcSi_postfit"), LineStyle(1), LineColor(col[2]), FillColor(col[2]), FillStyle(3002), DrawOption("FL"), RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcS),  Name("modelmcS_postfit"),  LineStyle(1), LineColor(col[3]), FillColor(col[3]), FillStyle(3002), DrawOption("FL"), RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcAr), Name("modelmcAr_postfit"), LineStyle(1), LineColor(col[4]), FillColor(col[4]), FillStyle(3002), DrawOption("FL"), RooFit::VLines(), MoveToBack());
	dframe->GetXaxis()->SetLimits(7.,20.);
	dframe->Draw();

	TLine *lmax = new TLine(Ar_high,dframe->GetMinimum(),Ar_high,dframe->GetMaximum());
	lmax->SetLineColor(kGray+2);
	lmax->SetLineStyle(9);
	lmax->SetLineWidth(3);
	//lmax->Draw();

	TLine *lmin = new TLine(Ar_low,dframe->GetMinimum(),Ar_low,dframe->GetMaximum());
	lmin->SetLineColor(kGray+2);
	lmin->SetLineStyle(9);
	lmin->SetLineWidth(3);
	//lmin->Draw();

	// Draw Legend
	TLegend *l = new TLegend(0.82,0.33,0.99,0.9);
	l->AddEntry("modelmcO_postfit","Oxygen","f");
	l->AddEntry("modelmcNe_postfit","Neon","f");
	l->AddEntry("modelmcMg_postfit","Magnesium","f");
	l->AddEntry("modelmcSi_postfit","Silicon","f");
	l->AddEntry("modelmcS_postfit","Sulfur","f");
	l->AddEntry("modelmcAr_postfit","Argon","f");
	l->AddEntry("modelmc_postfit","MC total ","l");
	l->AddEntry("data","Flight data","epl");
	l->SetLineColor(0);
	l->SetTextSize(0.035);
	l->Draw("same");

	TLatex latex; latex.SetTextSize(0.04); latex.SetNDC(true); latex.SetTextFont(42);
	latex.DrawLatex(0.5, 0.88, title);

	if(writefile){

		TString outfileName = "contaminazioni_argon_10bins_Ago26.txt";
		std::ofstream outFile;
		outFile.open(outfileName.Data(), std::ios::app);

		outFile << std::fixed << std::setprecision(6);
		outFile << iBin << " " << BGO_EnergyG_QuenchSatCorr_ML_ions2<< " ";

		for (size_t i = 0; i < yields.size(); ++i) {
			if (names[i] == "Ar") continue;
			RooFormulaVar cont_i(
				Form("cont_%s", names[i].c_str()),
				"(@0*@1)/@2",
				RooArgList(*yields[i], *integrals[i], Ntot_Ar)
			);
			double val = cont_i.getVal();
			double err = cont_i.getPropagatedError(*fitres);

			outFile << val << " " << err << " ";
		}

		outFile << std::endl;
		outFile.close();

		//TString plot_name1 = "PLOTS/template_argon_10bins_Jul26_"+h_data_name+".png";
		//TString plot_name2 = "PLOTS/template_argon_10bins_Jul26_"+h_data_name+".eps";
		//c1->SaveAs(plot_name1);
		//c1->SaveAs(plot_name2);
	}
}