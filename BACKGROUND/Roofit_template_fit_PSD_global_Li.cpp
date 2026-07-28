// Template fit on PSD charge data for Li selection, in one BGO energy bin

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

// per i bin 4, 5 mettere un rebin di 2. Per gli altri va bene 1.
void Roofit_template_fit_PSD_global_Li(int iBin = 4., int iRebin = 2.)
{ 

	double BGO_EnergyG = 0.0;
	//Metto il valore massimo di contaminazione
	/*
	if (iBin == 0) BGO_EnergyG = 31.6; 
	if (iBin == 1) BGO_EnergyG = 100.0; 
	if (iBin == 2) BGO_EnergyG = 316.2; 
	if (iBin == 3) BGO_EnergyG = 1000.0; 
	if (iBin == 4) BGO_EnergyG = 3162.3; 
	if (iBin == 5) BGO_EnergyG = 10000.0; 
	*/
/*
//Metto il valore centrale di contaminazione (media aritmetica)
if (iBin == 0) BGO_EnergyG = 20.8; 
if (iBin == 1) BGO_EnergyG = 65.8; 
if (iBin == 2) BGO_EnergyG = 208.114; 
if (iBin == 3) BGO_EnergyG = 658.114; 
if (iBin == 4) BGO_EnergyG = 2081.139; 
if (iBin == 5) BGO_EnergyG = 6581.139; 
*/
//Metto il valore centrale di contaminazione (media logaritmica)
if (iBin == 0) BGO_EnergyG = 17.78; 
if (iBin == 1) BGO_EnergyG = 56.23; 
if (iBin == 2) BGO_EnergyG = 177.82; 
if (iBin == 3) BGO_EnergyG = 562.3; 
if (iBin == 4) BGO_EnergyG = 1778.2; 
if (iBin == 5) BGO_EnergyG = 5623.0;
if (iBin == 6) BGO_EnergyG = 17782.0;

Double_t HeMPVf=((7.346)+(1.866)*log10(BGO_EnergyG)+(-0.9664)*log10(BGO_EnergyG)**2+(0.2305)*log10(BGO_EnergyG)**3+(-0.01202)*log10(BGO_EnergyG)**4);
Double_t HeWidthf=((-0.4515)+(1.265)*log10(BGO_EnergyG)+(-0.5451)*log10(BGO_EnergyG)**2+(0.08947)*log10(BGO_EnergyG)**3+(5.695E-05)*log10(BGO_EnergyG)**4);
Double_t HeGSigmaf=(0.3339)+(-2.305e-22)*log10(BGO_EnergyG);	

Double_t HeFSig= sqrt(HeWidthf**2+HeGSigmaf**2);

Double_t he_high6= HeMPVf+(5.5)*HeFSig;

cout << "high limit: "<< he_high6 << endl;

  double range_low = 0.0;
  double range_high = 25.0; 

  // retrieve FD and MC histograms 
  TFile *f_data = new TFile("/afs/le.infn.it/user/f/falemann/LITHIUM_CONTAMINATION/FD_PSDHisto25_HE.root");
//  std::string h_data_name = Form("h_psd_charge_global_e_bgo_bin_%d",iBin);
  std::string h_data_name = Form("h_%d",iBin);
  TH1D *h_data = (TH1D *)f_data->Get( (h_data_name).c_str() );
  h_data->Rebin(iRebin); 

  TFile *f_mc_p = new TFile("/afs/le.infn.it/user/f/falemann/LITHIUM_CONTAMINATION/MCpHe_PSDHisto25_HE.root");
  std::string h_mc_p_name = Form("h_%d",iBin); 
  TH1D *h_mc_p = (TH1D *)f_mc_p->Get( (h_mc_p_name).c_str() );
//	h_mc_p->Scale(2*3.14**(2)*1.38*1.38*86400*1393.8606);
  h_mc_p->Rebin(iRebin); 

  TFile *f_mc_Li = new TFile("/afs/le.infn.it/user/f/falemann/LITHIUM_CONTAMINATION/MCLi_PSDHisto25_HE.root");
  std::string h_mc_Li_name = Form("h_%d",iBin); 
  TH1D *h_mc_Li = (TH1D *)f_mc_Li->Get( (h_mc_Li_name).c_str() );
//	h_mc_Li->Scale(2*3.14**(2)*86400*1393.8606);
  h_mc_Li->Rebin(iRebin); 
/*
  TFile *f_mc_Be = new TFile("/afs/le.infn.it/user/f/falemann/LITHIUM_CONTAMINATION/MC_Be7_Be9_10GeV_100TeV_PSD_Global_5bpd_E_BGOE_bins_wPreCut_SmearCorr_wetaCorr_w2_1Weight.root");
  std::string h_mc_Be_name = Form("h_psd_charge_global_e_bgo_bin_%d",iBin); 
  TH1D *h_mc_Be = (TH1D *)f_mc_Be->Get( (h_mc_Be_name).c_str() );
  h_mc_Be->Rebin(iRebin);
*/
  // Create RooFit objects 
  RooRealVar x("x","PSD Global Energy (MeV)",0.0,25.0);			

  RooDataHist data("data","data",x,h_data);
  RooDataHist mcp("mcp","scaled mc Proton",x,h_mc_p);
  RooDataHist mcLi("mcLi","scaled mc Li",x,h_mc_Li);
  //RooDataHist mcBe("mcBe","scaled mc Be",x,h_mc_Be);

  // Define MC fractions 
  RooRealVar mcp_frac("mcp_frac","Fraction of mcProton",0.7* (h_data->GetEntries()),0.7,0.7* (h_data->GetEntries()) );	// These are variables for output
  RooRealVar mcLi_frac("mcLi_frac","Fraction of mcLi",0.3*( h_data->GetEntries()),0.3,0.3*( h_data->GetEntries()) );	// These are variables for output
  //RooRealVar mcBe_frac("mcBe_frac","Fracrion of mcBe",0.1*( h_data->GetEntries()),0,0.8*( h_data->GetEntries()) ); 	// These are variables for output

  // Create MC template model 
  RooHistPdf modelmcp("modelmcp","modelmcp",x, mcp);
  RooHistPdf modelmcLi("modelmcLi","modelmcLi",x, mcLi);
  //RooHistPdf modelmcBe("modelmcBe","modelmcBe",x, mcBe);
 
  RooAddPdf model("model","model",RooArgList(modelmcp,modelmcLi),RooArgList(mcp_frac,mcLi_frac));	// Combines my MCs into one PDF model
  //RooAddPdf model("model","model",RooArgList(modelmcp,modelmcLi,modelmcBe),RooArgList(mcp_frac,mcLi_frac,mcBe_frac));	// Combines my MCs into one PDF model

  // Fit template to MC 
  model.fitTo(data, "L", RooFit::Range(range_low,range_high), PrintLevel(0), RooFit::SumW2Error(true));
  
  TCanvas* c1 = new TCanvas("c1", "c1", 1300, 1040);
  //RooPlot* dframe = x.frame(range_low,range_high,100);
  RooPlot* dframe = x.frame(0.,25.,100);
  dframe->SetTitle("");

	c1->SetLogy();
  data.plotOn(dframe);
  model.plotOn(dframe, LineColor(kRed), LineWidth(4));
//  model.plotOn(dframe);
  model.plotOn(dframe, Components(modelmcLi), Name("modelmcLi_postfit"), LineStyle(kDotted), LineColor(kOrange-2), FillColor(kOrange-2), FillStyle(3022));
  model.plotOn(dframe, Components(modelmcp), Name("modelmcp_postfit"), LineStyle(kDotted), LineColor(kCyan-2), FillColor(kCyan-2), FillStyle(3022));
//  model.plotOn(dframe, Components(modelmcBe), Name("modelmcBe_postfit"), LineStyle(kDashed), LineColor(kGreen-2), FillColor(kGreen-2), FillStyle(3022));

  dframe->Draw(); 

TLine *lmax = new TLine(he_high6,dframe->GetMinimum(),he_high6,dframe->GetMaximum());
lmax->SetLineColor(kRed+1);
lmax->SetLineStyle(2);
lmax->SetLineWidth(2);
lmax->Draw();


  // Draw Legend 
  TLegend *l = new TLegend(0.6,0.7,0.9,0.9); 
  l->AddEntry("data","Flight data","epl"); 
  l->AddEntry("modelmcp_postfit","MC Proton+Helium","l"); 
  l->AddEntry("modelmcLi_postfit","MC Lithium","l"); 
//  l->AddEntry("model","MC H+He+Li","l"); 
//  l->AddEntry("modelmcBe_postfit","Beryllium ","l"); 
  l->SetTextSize(0.027);
  l->Draw("same");

  // Compute number of contamination events 
  x.setRange("Li_range",0.0,he_high6); 

  RooAbsReal *integral[2]; 
  double counts[2], d_counts[2], percentage[2], d_percentage[2]; 
  double tot_counts = 0.;
  double d_tot_counts = 0.; 
  RooHistPdf models[2] = {modelmcp, modelmcLi}; 
  RooRealVar fracs[2] = {mcp_frac, mcLi_frac}; 
  std::string mc_names[2] = {"p+He","Li"}; 
  //RooHistPdf models[3] = {modelmcp, modelmcHe, modelmcLi}; 
  //RooRealVar fracs[3] = {mcp_frac, mcHe_frac, mcLi_frac}; 
  //std::string mc_names[3] = {"p","He","Li"}; 
 
  cout<<"\n*****************************************"<<endl;
  cout<<" _______________________________________________"<<endl;
  cout<<"|Bin: "<<iBin<<" - charge selection region: 0.0 - "<<he_high6<<"|"<<endl;
  cout<<" _______________________________________________"<<endl;

  for(int iMC = 0; iMC<2; iMC++){
    integral[iMC] = models[iMC].createIntegral(x,RooFit::NormSet(x),RooFit::Range("Li_range")); 
    counts[iMC] = integral[iMC]->getVal() * fracs[iMC].getVal(); 
    d_counts[iMC] = fracs[iMC].getError() / fracs[iMC].getVal() * counts[iMC]; 
    tot_counts += counts[iMC]; 
    d_tot_counts += d_counts[iMC]; 
    cout<<"\n"+mc_names[iMC]+" counts in charge selection region = " << integral[iMC]->getVal() <<" * "<<fracs[iMC].getVal()<<" = "<<counts[iMC]<<" +- "<<d_counts[iMC]<<endl;
  }

  auto model_int = model.createIntegral( x, RooFit::NormSet(x),RooFit::Range("Li_range")); 
  

  cout<<"count check with total model "<< model_int->getVal() <<endl;

  // Compute contamination percentage

  double tot_percentage = 0; 

  for(int jMC = 0; jMC<2; jMC++){
    percentage[jMC] = counts[jMC] / tot_counts; 
    d_percentage[jMC] = percentage[jMC] * (d_counts[jMC] / counts[jMC] + d_tot_counts / tot_counts); 
    cout<<"\n"+mc_names[jMC]+" percentage = " <<percentage[jMC] << " +- "<<d_percentage[jMC]<<endl;
    tot_percentage += percentage[jMC]; 
  }

  cout<<"\nCHECK! Tot percentage: "<<tot_percentage<<endl;

ofstream fout;
   fout.open("pollutionHE5_5sigma.txt",ios::app);

	fout << iBin << " " << BGO_EnergyG <<"\t" << percentage[1]*100 << " " << d_percentage[1]*100 << " " << endl;
	fout.close();


double err_rel, err;

err_rel = 1./sqrt(counts[1]);

err = percentage[1]*err_rel;

cout << err << endl;

}
