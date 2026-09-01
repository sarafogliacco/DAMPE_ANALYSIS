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

void cfr_DATA_MC_argon_sigma(){

    //gStyle->SetOptFit(111);
    gStyle->SetOptFit(0);

    int i=0; 
    double a,b,c,d,e,f,g,h,j,k; 
    Double_t En[20], errEn[20], width[20], errwidth[20], sigma[20], errsigma[20], totSigma[20], errtotSigma[20];
    Double_t vecMC[10] = { 133.3521653, 237.1374976, 421.6965375, 749.8939925, 1333.521653, 2371.374976, 4216.965375, 7498.939925, 17782.80068, 56234.15332 };
    //Double_t vecMC[8] = { 133.3521653, 237.1374976, 421.6965375, 749.8939925, 1333.521653, 2371.374976, 5623.415332, 31622.7766 };

    TCanvas *c1 = new TCanvas;
    c1->SetLogx(); //c1->SetGridx(); c1->SetGridy();
    c1->Update();

    TF1 *f1 = new TF1("f1","-2",50.,300000.); char tit[1000];
    sprintf(tit, "; BGO energy (GeV); ARGON #sigma ");
    f1->SetTitle(tit); 
    f1->GetXaxis()->SetTitleOffset(1.2); f1->GetYaxis()->SetTitleOffset(1.2);
    f1->SetMinimum(-0.); 
    //f1->SetMaximum(0.6);
    f1->SetLineStyle(2); f1->SetLineWidth(2); f1->SetLineColor(2); 
    f1->Draw();

    ///////////////////////////////////// MONTE CARLO

   
    
    //ifstream infileMC("langaus_results/langaus_results_argon_MC.txt");
    ifstream infileMC("langaus_results/langaus_results_argon_MC_corrected_v2.txt");

    while(infileMC >> a >> b >> c >> d >> e >> f >> g >> h >> k >> j){
        En[i] = vecMC[i]; 
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
    TF1 *fitMC = new TF1("fitMC","[0]+[1]*log10(x)",100.,100000.);
    //TF1 *fitMC = new TF1("fitMC","[0]+[1]*log10(x)",100.,100000.);
    fitMC->SetLineColor(kViolet);

    ////////////////////////////////////
    for (int v=0; v<i; v++) En[v]=errEn[v]=width[v]=errwidth[v]=0.0;
    a=b=c=d=e=f=g=h=k=j=0.0;
    i=0;
    Double_t vec[8] = { 133.3521653, 237.1374976, 421.6965375, 749.8939925, 1333.521653, 2371.374976, 5623.415332, 31622.7766};

    ///////////////////////////////////// FLIGHT DATA
    ifstream infileDATA("langaus_results/langaus_results_argon_DATA_v2.txt");
    while(infileDATA >> a >> b >> c >> d >> e >> f >> g >> h >> k >> j){
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

    TGraphErrors *gwidth_data = new TGraphErrors(i,En,totSigma,errEn,errtotSigma);
    gwidth_data->SetMarkerStyle(20); gwidth_data->SetMarkerColor(39); gwidth_data->SetLineColor(39);
    gwidth_data->SetMarkerSize(0.8); gwidth_data->SetName("gwidth_data");

    //TF1 *fitDATA = new TF1("fitDATA","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3",100.,100000.);
    TF1 *fitDATA = new TF1("fitDATA","[0]+[1]*log10(x)",100.,100000.);
    //TF1 *fitDATA = new TF1("fitDATA","[0]+[1]*log10(x)",100.,100000.);
    fitDATA->SetLineColor(39);

    ///////////////////////////////////////////////////////////////////////////////////
    cout << endl;
    cout << " *** MC fit *** " << endl;
    cout << endl;
    gwidth_mc->Draw("p");
    //gwidth_mc->Draw("lsame");
    gwidth_mc->Fit("fitMC","r");
    cout << "TString sigma_MC = \"(" << fitMC->GetParameter(0) << "+(" << fitMC->GetParameter(1) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))" << "+(" << fitMC->GetParameter(2) << "*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))" << "+(" << fitMC->GetParameter(3) << "*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3))" << ")\";" << endl;
    cout << endl;

    cout << endl;
    cout << " *** DATA fit *** " << endl;
    cout << endl;
    gwidth_data->Draw("p");
    gwidth_data->Fit("fitDATA","r");
    cout << "TString sigma_DATA = \"(" << fitDATA->GetParameter(0) << "+(" << fitDATA->GetParameter(1) << "*log10(BGO_EnergyG_SatCorr_ML_ions))" << "+(" << fitDATA->GetParameter(2) << "*pow(log10(BGO_EnergyG_SatCorr_ML_ions),2))" << "+(" << fitDATA->GetParameter(3) << "*pow(log10(BGO_EnergyG_SatCorr_ML_ions),3))" << ")\";" << endl;
    cout << endl;

    auto legend = new TLegend(0.2,0.67,0.4,0.86);
    legend->AddEntry(gwidth_data," Ar DATA","lep");
    legend->AddEntry(gwidth_mc," Ar MC","lep");
    legend->SetLineColor(0);
    legend->Draw();
   
    //c1->SaveAs("plots/cfr_DATA_MC_argon_totSigma.png");
    //c1->SaveAs("plots/cfr_DATA_MC_argon_totSigma_CORRECTED.png");
    
}
