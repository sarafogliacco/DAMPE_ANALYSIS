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

void cfr_DATA_MC_sulfur_mpv(){

    //gStyle->SetOptFit(111);
    gStyle->SetOptFit(0);

    int i=0; 
    double a,b,c,d,e,f,g,h,j,k; 
    Double_t En[30], mpv[30], errEn[30], errmpv[30];
    Double_t vec[10] = { 133.3521653, 237.1374976, 421.6965375, 749.8939925, 1333.521653, 2371.374976, 4216.965375, 7498.939925, 17782.80068, 56234.15332};
    
    
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

    ///////////////////////////////////// MONTE CARLO

    //ifstream infileMC("langaus_results/results_langaus_fit_MC_SULFUR_v2.txt");
    ifstream infileMC("langaus_results/results_langaus_fit_MC_SULFUR_CORRECTED_v2.txt");
    while(infileMC >> a >> b >> c >> d >> e >> f >> g >> h >> k >> j){
        En[i] = vec[i]; 
        errEn[i] = 0;
        mpv[i] = d;
        errmpv[i] = e;
        i++;
    }

    TGraphErrors *gmpv_mc = new TGraphErrors(i,En,mpv,errEn,errmpv);
    gmpv_mc->SetMarkerStyle(20); gmpv_mc->SetMarkerColor(kViolet); gmpv_mc->SetLineColor(kViolet);
    gmpv_mc->SetMarkerSize(0.8); gmpv_mc->SetName("gmpv_mc"); 

    //TF1 *fitMC = new TF1("fitMC","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3",100.,100000.);
    TF1 *fitMC = new TF1("fitMC","[0]+[1]*log10(x)+[2]*log10(x)**2+[4]*log10(x)**4",100.,100000.);
    fitMC->SetLineColor(kViolet);

    ////////////////////////////////////
    for (int v=0; v<i; v++) En[v]=errEn[v]=mpv[v]=errmpv[v]=0.0;
    a=b=c=d=e=f=g=h=k=j=0.0;
    i=0;
    Double_t vecFD[8] = { 133.3521653, 237.1374976, 421.6965375, 749.8939925, 1333.521653, 2371.374976, 5623.415332, 31622.7766 };

    ///////////////////////////////////// FLIGHT DATA
    ifstream infileDATA("langaus_results/results_langaus_fit_DATA_SULFUR_v2.txt");
    while(infileDATA >> a >> b >> c >> d >> e >> f >> g >> h >> k >> j){
        En[i] = vecFD[i]; 
        errEn[i] = 0;
        mpv[i] = d;
        errmpv[i] = e;
        i++;
    }

    TGraphErrors *gmpv_data = new TGraphErrors(i,En,mpv,errEn,errmpv);
    gmpv_data->SetMarkerStyle(20); gmpv_data->SetMarkerColor(39); gmpv_data->SetLineColor(39);
    gmpv_data->SetMarkerSize(0.8); gmpv_data->SetName("gmpv_data");

    //TF1 *fitDATA = new TF1("fitDATA","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3",100.,100000.);
    TF1 *fitDATA = new TF1("fitDATA","[0]+[1]*log10(x)+[2]*log10(x)**2+[4]*log10(x)**4",100.,100000.);
    fitDATA->SetLineColor(39);

    ///////////////////////////////////////////////////////////////////////////////////
    cout << endl;
    cout << " *** MC fit *** " << endl;
    cout << endl;
    gmpv_mc->Draw("p");
    gmpv_mc->Fit("fitMC","r");
    //cout << "TString MPV_MC = \"(" << fitMC->GetParameter(0) << "+(" << fitMC->GetParameter(1) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))" << "+(" << fitMC->GetParameter(2) << "*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))" << "+(" << fitMC->GetParameter(3) << "*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3))" << ")\";" << endl;
    cout << "TString MPV_MC = \"(" << fitMC->GetParameter(0) << "+(" << fitMC->GetParameter(1) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))" << "+(" << fitMC->GetParameter(2) << "*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))" << "+(" << fitMC->GetParameter(4) << "*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),4))" << ")\";" << endl;
    cout << endl;

    cout << endl;
    cout << " *** DATA fit *** " << endl;
    cout << endl;
    gmpv_data->Draw("p");
    gmpv_data->Fit("fitDATA","r");
    //cout << "TString MPV_DATA = \"(" << fitDATA->GetParameter(0) << "+(" << fitDATA->GetParameter(1) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))" << "+(" << fitDATA->GetParameter(2) << "*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))" << "+(" << fitDATA->GetParameter(3) << "*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),3))" << ")\";" << endl;
    cout << "TString MPV_DATA = \"(" << fitDATA->GetParameter(0) << "+(" << fitDATA->GetParameter(1) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions))" << "+(" << fitDATA->GetParameter(2) << "*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),2))" << "+(" << fitDATA->GetParameter(4) << "*pow(log10(BGO_EnergyG_QuenchSatCorr_ML_ions),4))" << ")\";" << endl;
    cout << endl;

    auto legend = new TLegend(0.2,0.67,0.4,0.86);
    legend->AddEntry(gmpv_data," S DATA","lep");
    legend->AddEntry(gmpv_mc," S MC","lep");
    legend->SetLineColor(0);
    legend->Draw();
    
    //c1->SaveAs("plots/cfr_DATA_MC_sulfur_mpv.png");
    //c1->SaveAs("plots/cfr_DATA_MC_sulfur_mpv_CORRECTED.png");
    
}
