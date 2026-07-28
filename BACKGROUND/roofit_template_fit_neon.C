
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
#include "TCanvas.h"
#include "TLatex.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h" 
#include "RooLandau.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"

#include <fstream>
#include <iostream>

using namespace RooFit;
using namespace std;

void roofit_template_fit_neon(int iBin = 4, int iRebin = 2, bool writefile = false){

	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPadLeftMargin(0.09);
	gStyle->SetPadRightMargin(0.19);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadBottomMargin(0.09);

	double BGO_EnergyG = 0.0;
	TString title;
	/*
	// 9 bins (4 - 4 - 1) 
	if(iBin == 1) { BGO_EnergyG = 133.3521653; title = "100 GeV < E_{BGO} < 178 GeV"; }
	if(iBin == 2) { BGO_EnergyG = 237.1374976; title = "178 GeV < E_{BGO} < 316 GeV"; }
	if(iBin == 3) { BGO_EnergyG = 421.6965375; title = "316 GeV < E_{BGO} < 562 GeV"; }
	if(iBin == 4) { BGO_EnergyG = 749.8939925; title = "562 GeV < E_{BGO} < 1.0 TeV"; }
	if(iBin == 5) { BGO_EnergyG = 1333.521653; title = "1.0 TeV < E_{BGO} < 1.8 TeV"; }
	if(iBin == 6) { BGO_EnergyG = 2371.374976; title = "1.8 TeV < E_{BGO} < 3.1 TeV"; }
	if(iBin == 7) { BGO_EnergyG = 4216.965375; title = "3.1 TeV < E_{BGO} < 5.6 TeV"; }
	if(iBin == 8) { BGO_EnergyG = 7498.939925; title = "5.6 TeV < E_{BGO} < 10.0 TeV"; }
	if(iBin == 9) { BGO_EnergyG = 31622.7766; title = "10.0 TeV < E_{BGO} < 100.0 TeV"; }
	*/
	// 10 bins (6 - 3 - 1) 
	if(iBin == 1) { BGO_EnergyG = 121.1527961; title = "100 < E_{dep}/GeV < 147"; } // 100 GeV - 1 TeV
	if(iBin == 2) { BGO_EnergyG = 177.8277918; title = "147 < E_{dep}/GeV < 215"; }
	if(iBin == 3) { BGO_EnergyG = 261.015534 ; title = "215 < E_{dep}/GeV < 316"; }
	if(iBin == 4) { BGO_EnergyG = 383.1188748; title = "316 < E_{dep}/GeV < 464"; }
	if(iBin == 5) { BGO_EnergyG = 562.3413673; title = "464 < E_{dep}/GeV < 681"; }
	if(iBin == 6) { BGO_EnergyG = 825.4041434; title = "0.7 < E_{dep}/TeV < 1.0"; }
	if(iBin == 7) { BGO_EnergyG = 1467.79767 ; title = "1.0 < E_{dep}/TeV < 2.1"; } // 1 TeV - 10 TeV
	if(iBin == 8) { BGO_EnergyG = 3162.274615; title = "2.1 < E_{dep}/TeV < 4.6"; }
	if(iBin == 9) { BGO_EnergyG = 6812.921547; title = "4.6 < E_{dep}/TeV < 10";  }
	if(iBin ==10) { BGO_EnergyG = 31622.7766 ; title = " 10 < E_{dep}/TeV < 100"; } // 10 TeV - 100 TeV
	
	// NEON 02/2026
	Double_t Ne_high= (10.518 +(0.00615265*log10(BGO_EnergyG)*log10(BGO_EnergyG)*log10(BGO_EnergyG)) );
	Double_t Ne_low = (9.53331+ (0.0163751*log10(BGO_EnergyG)*log10(BGO_EnergyG)) );

	cout << "High limit: " << Ne_high << endl;
	cout << "Low limit: " << Ne_low << endl;

	double range_low  = 1.4;
	double range_high = 27.2;

	double range_low_temp  = 5.6;
    double range_high_temp = 15.;

	double Fe_range_low = 25.2;

	double He_range_high = 3.8;

	// retrieve FD and MC histograms 
	
	// On-orbit DATA
	TFile *f_data = new TFile("OUTPUTS/Skim_108months_MLSATcorr_STKcut1200_PSDcutXY_10bins.root");
	std::string h_data_name;
	if (iBin < 10){ h_data_name = Form("h0%d",iBin); }
	else h_data_name = Form("h%d",iBin);
	TH1D *h_data = (TH1D *)f_data->Get( (h_data_name).c_str() );
	h_data->Rebin(iRebin); 

	// MC Helium
	TFile *f_mc_He = new TFile("OUTPUTS/MC_FTFP_HELIUM_MLSATcorr_STKcut1200_PSDcutXY_CORRECTED_10bins.root");
	std::string h_mc_He_name;
	if (iBin < 10){ h_mc_He_name = Form("h0%d",iBin); }
	else h_mc_He_name = Form("h%d",iBin);
	TH1D *h_mc_He = (TH1D *)f_mc_He->Get( (h_mc_He_name).c_str() );
	h_mc_He->Rebin(iRebin); 

	// MC Carbon
	TFile *f_mc_C = new TFile("OUTPUTS/MC_FTFP_CARBON_MLSATcorr_STKcut1200_PSDcutXY_CORRECTED_10bins.root");
	std::string h_mc_C_name;
	if (iBin < 10){ h_mc_C_name = Form("h0%d",iBin); }
	else h_mc_C_name = Form("h%d",iBin);
	TH1D *h_mc_C = (TH1D *)f_mc_C->Get( (h_mc_C_name).c_str() );
	h_mc_C->Rebin(iRebin); 
	
	// MC Nitrogen 
	TFile *f_mc_N = new TFile("OUTPUTS/MC_FTFP_NITROGEN_MLSATcorr_STKcut1200_PSDcutXY_CORRECTED_10bins.root");
	std::string h_mc_N_name;
	if (iBin < 10){ h_mc_N_name = Form("h0%d",iBin); }
	else h_mc_N_name = Form("h%d",iBin);
	TH1D *h_mc_N = (TH1D *)f_mc_N->Get( (h_mc_N_name).c_str() );
	h_mc_N->Rebin(iRebin); 
	
	// MC Oxygen 
	TFile *f_mc_O = new TFile("OUTPUTS/MC_FTFP_OXYGEN_MLSATcorr_STKcut1200_PSDcutXY_CORRECTED_10bins.root");
	std::string h_mc_O_name;
	if (iBin < 10){ h_mc_O_name = Form("h0%d",iBin); }
	else h_mc_O_name = Form("h%d",iBin);
	TH1D *h_mc_O = (TH1D *)f_mc_O->Get( (h_mc_O_name).c_str() );
	h_mc_O->Rebin(iRebin); 
	
	// MC Fluorine
	TFile *f_mc_F = new TFile("OUTPUTS/MC_FTFP_FLUORINE_MLSATcorr_STKcut1200_PSDcutXY_CORRECTED_10bins.root");
	std::string h_mc_F_name;
	if (iBin < 10){ h_mc_F_name = Form("h0%d",iBin); }
	else h_mc_F_name = Form("h%d",iBin);
	TH1D *h_mc_F = (TH1D *)f_mc_F->Get( (h_mc_F_name).c_str() );
	h_mc_F->Rebin(iRebin); 

	// MC Neon 
	TFile *f_mc_Ne = new TFile("OUTPUTS/MC_FTFP_NEON_MLSATcorr_STKcut1200_PSDcutXY_CORRECTED_10bins.root");
	std::string h_mc_Ne_name;
	if (iBin < 10){ h_mc_Ne_name = Form("h0%d",iBin); }
	else h_mc_Ne_name = Form("h%d",iBin);
	TH1D *h_mc_Ne = (TH1D *)f_mc_Ne->Get( (h_mc_Ne_name).c_str() );
	h_mc_Ne->Rebin(iRebin); 

	// MC Sodium 
	TFile *f_mc_Na = new TFile("OUTPUTS/MC_FTFP_SODIUM_MLSATcorr_STKcut1200_PSDcutXY_CORRECTED_10bins_MPVfit.root");
	std::string h_mc_Na_name;
	if (iBin < 10){ h_mc_Na_name = Form("h0%d",iBin); }
	else h_mc_Na_name = Form("h%d",iBin);
	TH1D *h_mc_Na = (TH1D *)f_mc_Na->Get( (h_mc_Na_name).c_str() );
	h_mc_Na->Rebin(iRebin); 
	
	// MC Magnesium 
	TFile *f_mc_Mg = new TFile("OUTPUTS/MC_FTFP_MAGNESIUM_MLSATcorr_STKcut1200_PSDcutXY_CORRECTED_10bins.root");
	std::string h_mc_Mg_name;
	if (iBin < 10){ h_mc_Mg_name = Form("h0%d",iBin); }
	else h_mc_Mg_name = Form("h%d",iBin);
	TH1D *h_mc_Mg = (TH1D *)f_mc_Mg->Get( (h_mc_Mg_name).c_str() );
	h_mc_Mg->Rebin(iRebin); 

	// MC Aluminium 
	TFile *f_mc_Al = new TFile("OUTPUTS/MC_FTFP_ALUMINIUM_MLSATcorr_STKcut1200_PSDcutXY_CORRECTED_10bins.root");
	std::string h_mc_Al_name;
	if (iBin < 10){ h_mc_Al_name = Form("h0%d",iBin); }
	else h_mc_Al_name = Form("h%d",iBin);
	TH1D *h_mc_Al = (TH1D *)f_mc_Al->Get( (h_mc_Al_name).c_str() );
	h_mc_Al->Rebin(iRebin); 

	// MC Silicon 
	TFile *f_mc_Si = new TFile("OUTPUTS/MC_FTFP_SILICON_MLSATcorr_STKcut1200_PSDcutXY_CORRECTED_10bins.root");
	std::string h_mc_Si_name;
	if (iBin < 10){ h_mc_Si_name = Form("h0%d",iBin); }
	else h_mc_Si_name = Form("h%d",iBin);
	TH1D *h_mc_Si = (TH1D *)f_mc_Si->Get( (h_mc_Si_name).c_str() );
	h_mc_Si->Rebin(iRebin); 

	// MC Iron 
	TFile *f_mc_Fe = new TFile("OUTPUTS/MC_FTFP_IRON_MLSATcorr_STKcut1200_PSDcutXY_CORRECTED_10bins.root");
	std::string h_mc_Fe_name;
	if (iBin < 10){ h_mc_Fe_name = Form("h0%d",iBin); }
	else h_mc_Fe_name = Form("h%d",iBin);
	TH1D *h_mc_Fe = (TH1D *)f_mc_Fe->Get( (h_mc_Fe_name).c_str() );
	h_mc_Fe->Rebin(iRebin); 

	// Create RooFit objects 
	RooRealVar x("x","PSD Charge ",range_low,range_high);

	RooDataHist data("data","data",x,h_data);
	RooDataHist mcHe("mcHe","scaled mc He",x,h_mc_He);
	RooDataHist mcC( "mcC", "scaled mc C", x,h_mc_C);
	RooDataHist mcN( "mcN", "scaled mc N", x,h_mc_N);
	RooDataHist mcO( "mcO", "scaled mc O", x,h_mc_O);
	RooDataHist mcF( "mcF", "scaled mc F", x,h_mc_F);
	RooDataHist mcNe("mcNe","scaled mc Ne",x,h_mc_Ne);
	RooDataHist mcNa("mcNa","scaled mc Na",x,h_mc_Na);
	RooDataHist mcMg("mcMg","scaled mc Mg",x,h_mc_Mg);
	RooDataHist mcAl("mcAl","scaled mc Al",x,h_mc_Al);
	RooDataHist mcSi("mcSi","scaled mc Si",x,h_mc_Si);
	RooDataHist mcFe("mcFe","scaled mc Fe",x,h_mc_Fe);


	double frMC[11];
	if(iBin == 1) {
		double tmp[11] = { 0.00006, 0.002, 0.00042, 0.0035, 0.00007, 0.00071, 0.00014, 0.001, 0.00024, 0.0011, 0.0018 };
		for(int ii=0; ii<11; ii++) frMC[ii] = tmp[ii];
	}
	if(iBin == 2) {
		double tmp[11] = { 0.00008, 0.002, 0.00042, 0.0035, 0.00007, 0.00074, 0.00014, 0.0011, 0.00024, 0.0012, 0.002 };
		for(int ii=0; ii<11; ii++) frMC[ii] = tmp[ii];
	}
	if(iBin == 3) {
		double tmp[11] = { 0.00008, 0.002, 0.00042, 0.0035, 0.00007, 0.00074, 0.00014, 0.0011, 0.00024, 0.0012, 0.0022 }; 
		for(int ii=0; ii<11; ii++) frMC[ii] = tmp[ii];
	}
	if(iBin == 4) {
		double tmp[11] = { 0.0001, 0.002, 0.00042, 0.0035, 0.00007, 0.00078, 0.00014, 0.0012, 0.00024, 0.0013, 0.0026 }; 
		for(int ii=0; ii<11; ii++) frMC[ii] = tmp[ii];
	}
	if(iBin == 5) {
		double tmp[11] = { 0.0001, 0.002, 0.00042, 0.0035, 0.00007, 0.0008, 0.00014, 0.0012, 0.00024, 0.0014, 0.003 }; 
		for(int ii=0; ii<11; ii++) frMC[ii] = tmp[ii];
	}
	if(iBin == 6) {
		double tmp[11] = { 0.00012, 0.002, 0.00045, 0.0038, 0.00007, 0.0009, 0.00014, 0.0012, 0.00028, 0.0014, 0.0034 }; 
		for(int ii=0; ii<11; ii++) frMC[ii] = tmp[ii];
	}
	if(iBin == 7) {
		double tmp[11] = { 0.0003, 0.003, 0.00056, 0.0052, 0.0001, 0.0011, 0.00017, 0.0017, 0.00035, 0.002, 0.0045 };  
		for(int ii=0; ii<11; ii++) frMC[ii] = tmp[ii];
	}
	if(iBin == 8) {
		double tmp[11] = { 0.00083, 0.0034, 0.00064, 0.0057, 0.0001, 0.0014, 0.00017, 0.002, 0.00035, 0.0023, 0.0051 }; 
		for(int ii=0; ii<11; ii++) frMC[ii] = tmp[ii];
	}
	if(iBin == 9) {
		double tmp[11] = { 0.0018, 0.0043, 0.0007, 0.0074, 0.00020, 0.0018, 0.0002, 0.0022, 0.00044, 0.0023, 0.0047 }; 
		for(int ii=0; ii<11; ii++) frMC[ii] = tmp[ii];
	}
	if(iBin == 10) {
		double tmp[11] = { 0.0018, 0.0043, 0.001, 0.0074, 0.00020, 0.0018, 0.00025, 0.0026, 0.0005, 0.0026, 0.0041 }; 
		for(int ii=0; ii<11; ii++) frMC[ii] = tmp[ii];
	}
	

	//					He  	C 	  N     O     F     Ne    Na    Mg     Al    Si    Fe
	//double frMC[11] = { 0.01, 0.25, 0.06, 0.37, 0.01, 0.07, 0.02, 0.10, 0.025, 0.10, 0.03 }; 

	// Define MC fractions 
	RooRealVar mcHe_frac("mcHe_frac","Fraction of mcHe",frMC[0]*( h_data->GetEntries()), 0.0, 10.0 * h_data->GetEntries() );	// These are variables for output
	RooRealVar mcC_frac( "mcC_frac", "Fraction of mcC", frMC[1]*( h_data->GetEntries()), 0.0, 10.0 * h_data->GetEntries() );	// These are variables for output
	RooRealVar mcN_frac( "mcN_frac", "Fraction of mcN", frMC[2]*( h_data->GetEntries()), 0.0, 10.0 * h_data->GetEntries() );	// These are variables for output
	RooRealVar mcO_frac( "mcO_frac", "Fraction of mcO", frMC[3]*( h_data->GetEntries()), 0.0, 10.0 * h_data->GetEntries() );	// These are variables for output
	RooRealVar mcF_frac( "mcF_frac", "Fraction of mcF", frMC[4]*( h_data->GetEntries()), 0.0, 10.0 * h_data->GetEntries() );	// These are variables for output
	RooRealVar mcNe_frac("mcNe_frac","Fraction of mcNe",frMC[5]*( h_data->GetEntries()), 0.0, 10.0 * h_data->GetEntries() );	// These are variables for output
	RooRealVar mcNa_frac("mcNa_frac","Fraction of mcNa",frMC[6]*( h_data->GetEntries()), 0.0, 10.0 * h_data->GetEntries() );	// These are variables for output
	RooRealVar mcMg_frac("mcMg_frac","Fraction of mcMg",frMC[7]*( h_data->GetEntries()), 0.0, 10.0 * h_data->GetEntries() );	// These are variables for output
	RooRealVar mcAl_frac("mcAl_frac","Fraction of mcAl",frMC[8]*( h_data->GetEntries()), 0.0, 10.0 * h_data->GetEntries() );	// These are variables for output
	RooRealVar mcSi_frac("mcSi_frac","Fraction of mcSi",frMC[9]*( h_data->GetEntries()), 0.0, 10.0 * h_data->GetEntries() );	// These are variables for output
	RooRealVar mcFe_frac("mcFe_frac","Fraction of mcFe",frMC[10]*(h_data->GetEntries()), 0.0, 10.0 * h_data->GetEntries() );	// These are variables for output
	
	// Create MC template model 
	RooHistPdf modelmcHe("modelmcHe","modelmcHe",x, mcHe);
	RooHistPdf modelmcC( "modelmcC", "modelmcC", x, mcC);
	RooHistPdf modelmcN( "modelmcN", "modelmcN", x, mcN);
	RooHistPdf modelmcO( "modelmcO", "modelmcO", x, mcO);
	RooHistPdf modelmcF( "modelmcF", "modelmcF", x, mcF);
	RooHistPdf modelmcNe("modelmcNe","modelmcNe",x, mcNe);
	RooHistPdf modelmcNa("modelmcNa","modelmcNa",x, mcNa);
	RooHistPdf modelmcMg("modelmcMg","modelmcMg",x, mcMg);
	RooHistPdf modelmcAl("modelmcAl","modelmcAl",x, mcAl);
	RooHistPdf modelmcSi("modelmcSi","modelmcSi",x, mcSi);
	RooHistPdf modelmcFe("modelmcFe","modelmcFe",x, mcFe);

	RooAddPdf model("model","model",RooArgList(modelmcHe,modelmcC,modelmcN,modelmcO,modelmcF,modelmcNe,modelmcNa,modelmcMg,modelmcAl,modelmcSi,modelmcFe),RooArgList(mcHe_frac,mcC_frac,mcN_frac,mcO_frac,mcF_frac,mcNe_frac,mcNa_frac,mcMg_frac,mcAl_frac,mcSi_frac,mcFe_frac));	// Combines my MCs into one PDF model

	// Define side band regions and full range
	x.setRange("Neon",range_low_temp, range_high_temp);
 	x.setRange("Iron",Fe_range_low,range_high);
 	//x.setRange("FULL",range_low,range_high);
 	x.setRange("Helium",range_low,He_range_high);
 	x.setRange("Neon_2",Ne_low, Ne_high);
 	

	// Fit template to MC 
	//model.fitTo(data, RooFit::Range("Iron,Neon"), RooFit::PrintLevel(0), RooFit::SumW2Error(false));
	//RooFitResult* res = model.fitTo(data, RooFit::Range("Helium,Iron,Neon"), RooFit::PrintLevel(0), RooFit::SumW2Error(false));
 	RooFitResult* fitres = model.fitTo(data, RooFit::Range("Helium,Iron,Neon"), RooFit::Save(), RooFit::PrintLevel(0));

	

	RooAbsReal* I_He = modelmcHe.createIntegral(x, NormSet(x), Range("Neon_2"));
	RooAbsReal* I_C  = modelmcC .createIntegral(x, NormSet(x), Range("Neon_2"));
	RooAbsReal* I_N  = modelmcN .createIntegral(x, NormSet(x), Range("Neon_2"));
	RooAbsReal* I_O  = modelmcO .createIntegral(x, NormSet(x), Range("Neon_2"));
	RooAbsReal* I_F  = modelmcF .createIntegral(x, NormSet(x), Range("Neon_2"));
	RooAbsReal* I_Ne = modelmcNe.createIntegral(x, NormSet(x), Range("Neon_2"));
	RooAbsReal* I_Na = modelmcNa.createIntegral(x, NormSet(x), Range("Neon_2"));
	RooAbsReal* I_Mg = modelmcMg.createIntegral(x, NormSet(x), Range("Neon_2"));
	RooAbsReal* I_Al = modelmcAl.createIntegral(x, NormSet(x), Range("Neon_2"));
	RooAbsReal* I_Si = modelmcSi.createIntegral(x, NormSet(x), Range("Neon_2"));
	RooAbsReal* I_Fe = modelmcFe.createIntegral(x, NormSet(x), Range("Neon_2"));

	std::vector<RooRealVar*> yields = {
		&mcHe_frac, &mcC_frac, &mcN_frac, &mcO_frac, &mcF_frac,
		&mcNe_frac, &mcNa_frac, &mcMg_frac, &mcAl_frac, &mcSi_frac, &mcFe_frac
	};

	std::vector<RooAbsReal*> integrals = { I_He, I_C, I_N, I_O, I_F, I_Ne, I_Na, I_Mg, I_Al, I_Si, I_Fe };
	std::vector<std::string> names = { "He","C","N","O","F","Ne","Na","Mg","Al","Si","Fe" };


	RooArgList denomArgs;
	std::string denomExpr;

	for (size_t i = 0; i < yields.size(); ++i) {
		denomExpr += Form("@%zu*@%zu", 2*i, 2*i+1);
		if (i < yields.size()-1) denomExpr += " + ";
		denomArgs.add(*yields[i]);
		denomArgs.add(*integrals[i]);
	}

	RooFormulaVar Ntot_Ne("Ntot_Ne", denomExpr.c_str(), denomArgs);

	for (size_t i = 0; i < yields.size(); ++i) {

	if (names[i] == "Ne") continue;

	RooFormulaVar cont_i(
		Form("cont_%s", names[i].c_str()),
		"(@0*@1)/@2",
		RooArgList(*yields[i], *integrals[i], Ntot_Ne)
	);

	double val = cont_i.getVal();
	double err = cont_i.getPropagatedError(*fitres);

	std::cout << "Contaminazione " << names[i]
			  << " nel Neon = "
			  << val << " ± " << err << std::endl;
	}


	RooFormulaVar purity_Ne(
    "purity_Ne",
    "(@0*@1)/@2",
    RooArgList(mcNe_frac, *I_Ne, Ntot_Ne)
);

double purity = purity_Ne.getVal();
double purity_err = purity_Ne.getPropagatedError(*fitres);

double total_cont = 1.0 - purity;
double total_cont_err = purity_err;

std::cout << "Somma contaminazioni = "
          << (1.0 - purity) << std::endl;



	TCanvas* c1 = new TCanvas("c1", "c1", 1300, 980);
	RooPlot* dframe = x.frame(range_low,range_high,100);
	dframe->SetTitle(" ");
	//dframe->SetMinimum(4); //dframe->SetMaximum(iMax);

	TString colCode[11] = {"#6a4c93", "#b30000", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0", "#4C4C4C" };
	Int_t col[11];
	for(int ii=0; ii<11; ii++){ col[ii]=TColor::GetColor(colCode[ii]);} 

	c1->SetLogy();
	data.plotOn(dframe);
	model.plotOn(dframe, Name("modelmc_postfit"), LineColor(kRed), LineWidth(4));
	model.plotOn(dframe, Components(modelmcHe), Name("modelmcHe_postfit"), LineStyle(1), LineColor(col[0]), FillColor(col[0]), FillStyle(3002), DrawOption("FL"),RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcC),  Name("modelmcC_postfit"),  LineStyle(1), LineColor(col[1]), FillColor(col[1]), FillStyle(3002), DrawOption("FL"),RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcN),  Name("modelmcN_postfit"),  LineStyle(1), LineColor(col[2]), FillColor(col[2]), FillStyle(3002), DrawOption("FL"),RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcO),  Name("modelmcO_postfit"),  LineStyle(1), LineColor(col[3]), FillColor(col[3]), FillStyle(3002), DrawOption("FL"),RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcF),  Name("modelmcF_postfit"),  LineStyle(1), LineColor(col[4]), FillColor(col[4]), FillStyle(3002), DrawOption("FL"),RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcNe), Name("modelmcNe_postfit"), LineStyle(1), LineColor(col[5]), FillColor(col[5]), FillStyle(3002), DrawOption("FL"),RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcNa), Name("modelmcNa_postfit"), LineStyle(1), LineColor(col[6]), FillColor(col[6]), FillStyle(3002), DrawOption("FL"),RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcMg), Name("modelmcMg_postfit"), LineStyle(1), LineColor(col[7]), FillColor(col[7]), FillStyle(3002), DrawOption("FL"),RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcAl), Name("modelmcAl_postfit"), LineStyle(1), LineColor(col[8]), FillColor(col[8]), FillStyle(3002), DrawOption("FL"),RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcSi), Name("modelmcSi_postfit"), LineStyle(1), LineColor(col[9]), FillColor(col[9]), FillStyle(3002), DrawOption("FL"),RooFit::VLines(), MoveToBack());
	model.plotOn(dframe, Components(modelmcFe), Name("modelmcFe_postfit"), LineStyle(1), LineColor(col[10]), FillColor(col[10]), FillStyle(3002), DrawOption("FL"),RooFit::VLines(), MoveToBack());
	dframe->Draw(); 

	TLine *lmax = new TLine(Ne_high,dframe->GetMinimum(),Ne_high,dframe->GetMaximum());
	lmax->SetLineColor(kGray+2);
	lmax->SetLineStyle(9);
	lmax->SetLineWidth(3);
	//lmax->Draw();

	TLine *lmin = new TLine(Ne_low,dframe->GetMinimum(),Ne_low,dframe->GetMaximum());
	lmin->SetLineColor(kGray+2);
	lmin->SetLineStyle(9);
	lmin->SetLineWidth(3);
	//lmin->Draw();

	// Draw Legend 
	TLegend *l = new TLegend(0.82,0.33,0.99,0.9); 
	l->AddEntry("modelmcHe_postfit","Helium","f"); 
	l->AddEntry("modelmcC_postfit", "Carbon","f"); 
	l->AddEntry("modelmcN_postfit", "Nitrogen","f");
	l->AddEntry("modelmcO_postfit", "Oxygen","f");
	l->AddEntry("modelmcF_postfit", "Fluorine","f");
	l->AddEntry("modelmcNe_postfit","Neon","f");
	l->AddEntry("modelmcNa_postfit","Sodium","f");
	l->AddEntry("modelmcMg_postfit","Magnesium","f");
	l->AddEntry("modelmcAl_postfit","Aluminium","f");
	l->AddEntry("modelmcSi_postfit","Silicon","f"); 
	l->AddEntry("modelmcFe_postfit","Iron","f"); 
	l->AddEntry("modelmc_postfit","MC total ","l");
	l->AddEntry("data","Flight data","epl"); 
	l->SetLineColor(0);
	l->SetTextSize(0.035);
	l->Draw("same");
	TLatex latex; latex.SetTextSize(0.04); latex.SetNDC(true); latex.SetTextFont(42);
	latex.DrawLatex(0.5, 0.88, title);  

	if(writefile){

		// Nome del file di output (lo crea se non esiste, altrimenti aggiunge righe)
		TString outfileName = "contaminazioni_neon_10bins_Feb26.txt";
		std::ofstream outFile;
		outFile.open(outfileName.Data(), std::ios::app); // append mode

		// Imposta precisione decimale
		outFile << std::fixed << std::setprecision(6);

		// Scrive iBin e BGO_EnergyG
		outFile << iBin << " " << BGO_EnergyG << " ";

		// Loop per scrivere contaminazioni e errori
		for (size_t i = 0; i < yields.size(); ++i) {
    		if (names[i] == "Ne") continue; // salta il neon stesso
    		RooFormulaVar cont_i(
        		Form("cont_%s", names[i].c_str()),
        		"(@0*@1)/@2",
        		RooArgList(*yields[i], *integrals[i], Ntot_Ne)
    		);
    		double val = cont_i.getVal();
    		double err = cont_i.getPropagatedError(*fitres);

    		outFile << val << " " << err << " ";
		}

		// Aggiunge a fine riga un newline
		outFile << std::endl;
		// Chiude il file
		outFile.close();

		TString plot_name1 = "PLOTS/template_neon_10bins_Feb26_"+h_data_name+".png";
		TString plot_name2 = "PLOTS/template_neon_10bins_Feb26_"+h_data_name+".eps";
		c1->SaveAs(plot_name1);
		c1->SaveAs(plot_name2);

	}




}
