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

void cfr_DATA_PSDXY_sulfur_sigma(){

    //gStyle->SetOptFit(111);
    gStyle->SetOptFit(0);

    int i=0; 
    double a,b,c,d,e,f,g,h,j,k; 
    Double_t En[20], errEn[20], width[20], errwidth[20], sigma[20], errsigma[20], totSigma[20], errtotSigma[20];
    Double_t vec[10] = { 133.3521653, 237.1374976, 421.6965375, 749.8939925, 1333.521653, 2371.374976, 4216.965375, 7498.939925, 17782.80068, 56234.15332 };
    //Double_t vec[8] = { 133.3521653, 237.1374976, 421.6965375, 749.8939925, 1333.521653, 2371.374976, 5623.415332, 31622.7766 };
    

    TCanvas *c1 = new TCanvas;
    c1->SetLogx(); //c1->SetGridx(); c1->SetGridy();
    c1->Update();

    TF1 *f1 = new TF1("f1","-2",50.,300000.); char tit[1000];
    sprintf(tit, "; BGO energy (GeV); SULFUR #sigma ");
    f1->SetTitle(tit); 
    f1->GetXaxis()->SetTitleOffset(1.2); f1->GetYaxis()->SetTitleOffset(1.2);
    f1->SetMinimum(-0.); 
    //f1->SetMaximum(0.6);
    f1->SetLineStyle(2); f1->SetLineWidth(2); f1->SetLineColor(2); 
    f1->Draw();

    ///////////////////////////////////// MONTE CARLO

    //ifstream infileMC("lanngaus_results/results_langaus_fit_MC_SULFUR.txt");
    //ifstream infileDATAX("PSDX/results_langaus_fit_DATA_PSDX.txt");
    ifstream infileDATAX("PSDX/results_langaus_fit_MC_PSDX.txt");
    //ifstream infileDATAX("PSDX/results_langaus_fit_MC_CORRECTED_PSDX.txt");
    
    while(infileDATAX >> a >> b >> c >> d >> e >> f >> g >> h >> k >> j){
        En[i] = vec[i]; 
        errEn[i] = 0;
        width[i] = b;
        errwidth[i] = c;
        sigma[i] = h;
        errsigma[i] = k;
        totSigma[i] = sqrt(b*b + h*h);
        errtotSigma[i] = sqrt( ( (b*c)*(b*c) + (h*k)*(h*k) )/(b*b + h*h) );
        i++;
    }

    TGraphErrors *gwidth_mc = new TGraphErrors(i,En,totSigma,errEn,errtotSigma);
    gwidth_mc->SetMarkerStyle(20); gwidth_mc->SetMarkerColor(kViolet); gwidth_mc->SetLineColor(kViolet);
    gwidth_mc->SetMarkerSize(0.8); gwidth_mc->SetName("gwidth_mc");

    //TF1 *fitMC = new TF1("fitMC","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3",100.,100000.);
    TF1 *fitDATAX = new TF1("fitDATAX","[0]+[1]*log10(x)+[2]*log10(x)**2",100.,100000.);
    //TF1 *fitMC = new TF1("fitMC","[0]+[1]*log10(x)",100.,100000.);
    fitDATAX->SetLineColor(kViolet);

    ///////////////////////////////////// FLIGHT DATA
    int ii=0; 
    //ifstream infileDATAY("PSDY/results_langaus_fit_DATA_PSDY.txt");
    ifstream infileDATAY("PSDY/results_langaus_fit_MC_PSDY.txt");
    //ifstream infileDATAY("PSDY/results_langaus_fit_MC_CORRECTED_PSDY.txt");

    while(infileDATAY >> a >> b >> c >> d >> e >> f >> g >> h >> k >> j){
        En[ii] = vec[ii]; 
        errEn[ii] = 0;
        width[ii] = b;
        errwidth[ii] = c;
        sigma[ii] = h;
        errsigma[ii] = k;
        totSigma[ii] = sqrt(b*b + h*h);
        errtotSigma[ii] = sqrt( ( (b*c)*(b*c) + (h*k)*(h*k) )/(b*b + h*h) );
        ii++;
    }

    TGraphErrors *gwidth_data = new TGraphErrors(ii,En,totSigma,errEn,errtotSigma);
    gwidth_data->SetMarkerStyle(20); gwidth_data->SetMarkerColor(39); gwidth_data->SetLineColor(39);
    gwidth_data->SetMarkerSize(0.8); gwidth_data->SetName("gwidth_data");

    //TF1 *fitDATA = new TF1("fitDATA","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3",100.,100000.);
    TF1 *fitDATAY = new TF1("fitDATAY","[0]+[1]*log10(x)+[2]*log10(x)**2",100.,100000.);
    //TF1 *fitDATA = new TF1("fitDATA","[0]+[1]*log10(x)",100.,100000.);
    fitDATAY->SetLineColor(39);

    ///////////////////////////////////////////////////////////////////////////////////
    cout << endl;
    cout << " *** MC fit *** " << endl;
    cout << endl;
    gwidth_mc->Draw("p");
    //gwidth_mc->Draw("lsame");
    gwidth_mc->Fit("fitDATAX","r");
    cout << "TString sigma_DATA_X = \"(" << fitDATAX->GetParameter(0) << "+(" << fitDATAX->GetParameter(1) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))" << "+(" << fitDATAX->GetParameter(2) << "*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))" << "+(" << fitDATAX->GetParameter(3) << "*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3))" << ")\";" << endl;
    cout << endl;

    cout << endl;
    cout << " *** DATA fit *** " << endl;
    cout << endl;
    gwidth_data->Draw("p");
    gwidth_data->Fit("fitDATAY","r");
    cout << "TString sigma_DATA_Y = \"(" << fitDATAY->GetParameter(0) << "+(" << fitDATAY->GetParameter(1) << "*log10(BGO_EnergyG_SatCorr_ML_ions))" << "+(" << fitDATAY->GetParameter(2) << "*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))" << "+(" << fitDATAY->GetParameter(3) << "*pow(log10(BGO_EnergyG_SatCorr_ML_ions),3))" << ")\";" << endl;
    cout << endl;

    auto legend = new TLegend(0.2,0.67,0.4,0.86);
    legend->AddEntry(gwidth_data," S DATA PSD_Y","lep");
    legend->AddEntry(gwidth_mc," S DATA PSD_X","lep");
    legend->SetLineColor(0);
    legend->Draw();
   
    //c1->SaveAs("plots/cfr_DATA_MC_sulfur_totSigma.png");
    //c1->SaveAs("plots/cfr_DATA_MC_silicon_totSigma_CORRECTED.png");
    
}
