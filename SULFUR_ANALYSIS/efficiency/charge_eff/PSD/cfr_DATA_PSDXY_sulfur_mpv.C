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

void cfr_DATA_PSDXY_sulfur_mpv(){

    //gStyle->SetOptFit(111);
    gStyle->SetOptFit(0);

    int i=0; 
    double a,b,c,d,e,f,g,h,j,k; 
    Double_t En[30], mpv[30], errEn[30], errmpv[30];
    //Double_t vec[10] = { 133.3521653, 237.1374976, 421.6965375, 749.8939925, 1333.521653, 2371.374976, 4216.965375, 7498.939925, 17782.80068, 56234.15332};
    Double_t vec[8] = { 133.3521653, 237.1374976, 421.6965375, 749.8939925, 1333.521653, 2371.374976, 5623.415332, 31622.7766 };
    
    
    TCanvas *c1 = new TCanvas;
    c1->SetLogx(); //c1->SetGridx(); c1->SetGridy();
    c1->Update();

    TF1 *f1 = new TF1("f1","-0.5",50.,300000.); char tit[1000];
    sprintf(tit, "; BGO energy (GeV); SULFUR MPV ");
    f1->SetTitle(tit); 
    f1->GetXaxis()->SetTitleOffset(1.2); f1->GetYaxis()->SetTitleOffset(1.2);
    f1->SetMinimum(15.); 
    f1->SetMaximum(18.);
    f1->SetLineStyle(2); f1->SetLineWidth(2); f1->SetLineColor(2); 
    f1->Draw();

    //ifstream infileDATAX("PSDX/results_langaus_fit_DATA_PSDX.txt");
    ifstream infileDATAX("PSDX/results_langaus_fit_MC_PSDX.txt");
    //ifstream infileDATAX("PSDX/results_langaus_fit_MC_CORRECTED_PSDX.txt");
    
    while(infileDATAX >> a >> b >> c >> d >> e >> f >> g >> h >> k >> j){
        En[i] = vec[i]; 
        errEn[i] = 0;
        mpv[i] = d;
        errmpv[i] = e;
        i++;
    }

    TGraphErrors *gmpv_mc = new TGraphErrors(i,En,mpv,errEn,errmpv);
    gmpv_mc->SetMarkerStyle(20); gmpv_mc->SetMarkerColor(kViolet); gmpv_mc->SetLineColor(kViolet);
    gmpv_mc->SetMarkerSize(0.8); gmpv_mc->SetName("gmpv_mc"); 

    TF1 *fitDATAX = new TF1("fitDATAX","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3",100.,100000.);
    fitDATAX->SetLineColor(kViolet);

    int ii=0; 
    //ifstream infileDATAY("PSDY/results_langaus_fit_DATA_PSDY.txt");
    ifstream infileDATAY("PSDY/results_langaus_fit_MC_PSDY.txt");
    //ifstream infileDATAY("PSDY/results_langaus_fit_MC_CORRECTED_PSDY.txt");
    while(infileDATAY >> a >> b >> c >> d >> e >> f >> g >> h >> k >> j){
        En[ii] = vec[ii]; 
        errEn[ii] = 0;
        mpv[ii] = d;
        errmpv[ii] = e;
        ii++;
    }

    TGraphErrors *gmpv_data = new TGraphErrors(ii,En,mpv,errEn,errmpv);
    gmpv_data->SetMarkerStyle(20); gmpv_data->SetMarkerColor(39); gmpv_data->SetLineColor(39);
    gmpv_data->SetMarkerSize(0.8); gmpv_data->SetName("gmpv_data");

    TF1 *fitDATAY = new TF1("fitDATAY","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3",100.,100000.);
    fitDATAY->SetLineColor(39);

    ///////////////////////////////////////////////////////////////////////////////////
    cout << endl;
    cout << " *** DATAX fit *** " << endl;
    cout << endl;
    gmpv_mc->Draw("p");
    gmpv_mc->Fit("fitDATAX","r");
    cout << "TString MPV_DATAX = \"(" << fitDATAX->GetParameter(0) << "+(" << fitDATAX->GetParameter(1) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))" << "+(" << fitDATAX->GetParameter(2) << "*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))" << "+(" << fitDATAX->GetParameter(3) << "*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3))" << ")\";" << endl;
    cout << endl;

    cout << endl;
    cout << " *** DATAY fit *** " << endl;
    cout << endl;
    gmpv_data->Draw("p");
    gmpv_data->Fit("fitDATAY","r");
    cout << "TString MPV_DATAY = \"(" << fitDATAY->GetParameter(0) << "+(" << fitDATAY->GetParameter(1) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))" << "+(" << fitDATAY->GetParameter(2) << "*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))" << "+(" << fitDATAY->GetParameter(3) << "*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3))" << ")\";" << endl;
    cout << endl;

    auto legend = new TLegend(0.2,0.67,0.4,0.86);
    legend->AddEntry(gmpv_data," S DATAY","lep");
    legend->AddEntry(gmpv_mc," S DATAX","lep");
    legend->SetLineColor(0);
    legend->Draw();
    
    //c1->SaveAs("plots/cfr_DATA_MC_sulfur_mpv.png");
    //c1->SaveAs("plots/cfr_DATA_MC_sulfur_mpv_CORRECTED.png");
    
}
