#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLatex.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <string>



//////////////////////////////////////////// LANGAUSS FUNCTION /////////////////////////////////////////////////

Double_t langaufun(Double_t *x, Double_t *par) {

   // DEFINIZIONE PARAMETRI:

   //   par[0]=Ampiezza della Landau
   //   par[1]=Most Probable (MP) della Landau
   //   par[2]=Area totale (costante di normalizzazione)
   //   par[3]=Sigma della Gaussiana

      // Costanti numeriche
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Costanti di controllo
      Double_t np = 100.0;      // numero di step di convoluzione
      Double_t sc =   5.0;      // convoluzione estesa a +-sc sigma della Gaussiana

      // Variabili
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;

      // correzione dello shift MP
      mpc = par[1] - mpshift * par[0];

      // range dell'integrale di convoluzione
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Integrale di convoluzione Landau+Gaussiana per somma

      for( i = 1.0; i <= np/2; i++) { xx = xlow + (i-.5) * step;
                                                  fland = TMath::Landau(xx,mpc,par[0]) / par[0];
                                                  sum += fland * TMath::Gaus(x[0],xx,par[3]);
                                                  xx = xupp - (i-.5) * step;
                                                  fland = TMath::Landau(xx,mpc,par[0]) / par[0];
                                                  sum += fland * TMath::Gaus(x[0],xx,par[3]);      }
      
      return (par[2] * step * sum * invsq2pi / par[3]); 
}

////////////////////////////////////////////////////LANGAUSS FIT///////////////////////////////////////////////////////////////////////////////////////

TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF){  
   // Variabili per il fit:
   //   his             istogramma da fittare
   //   fitrange[2]     range del fit
   //   startvalues[4]  valori iniziali di ciascun parametro
   //   parlimitslo[4]  limite inferiore di ciascun parametro
   //   parlimitshi[4]  limite superiore di ciascun parametro
   //   fitparams[4]    restituisce i valori finali dei parametri del fit
   //   fiterrors[4]    restituisce gli errori sui parametri del fit
   //   ChiSqr          restituisce il chi quadro
   //   NDF             restituisce ndf

   Int_t i;
   Char_t Datafit[100];

   sprintf(Datafit,"Fitfcn_%s",his->GetName());

   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(Datafit);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(Datafit,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MPV","Area","GSigma");

   for ( i=0; i < 4; i++) { ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);   }

   his->Fit(Datafit,"RB0");   

   ffit->GetParameters(fitparams);    
   for (i=0; i<4; i++) {   fiterrors[i] = ffit->GetParError(i);    }
   ChiSqr[0] = ffit->GetChisquare();  
   NDF[0] = ffit->GetNDF();         

   return (ffit);    
}

Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {

   // Cerca la x del massimo della convoluzione Landau-Gaussiana e la sua FWHM
   
   Double_t p,x,fy,fxr,fxl;
   Double_t step;
   Double_t l,lold;
   Int_t i = 0;
   Int_t MAXCALLS = 10000;

   // Ricerca del massimo

   p = params[1] - 0.1 * params[0];
   step = 0.05 * params[0];
   lold = -2.0;
   l    = -1.0;

   while ( (l != lold) && (i < MAXCALLS) ) { i++;
      lold = l;
      x = p + step;
      l = langaufun(&x,params);
      if (l < lold)         step = -step/10;
      p += step;  }

   if (i == MAXCALLS) return (-1);

   maxx = x;
   fy = l/2;

   // Search for right x location of fy

   p = maxx + params[0];
   step = params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;

   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;
      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
      if (l > lold)         step = -step/10;
      p += step;  }

   if (i == MAXCALLS)      return (-2);

   fxr = x;

   // Search for left x location of fy

   p = maxx - 0.5 * params[0];
   step = -params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;

   while ( (l != lold) && (i < MAXCALLS) ) {      i++;
      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
      if (l > lold) step = -step/10;
      p += step;   }

   if (i == MAXCALLS)      return (-3);

   fxl = x;

   FWHM = fxr - fxl;
   return (0); 
}


void langaus_fit_MC(){

   TFile *f = new TFile("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/BACKGROUND/out_root/MC_FTFP_NEON_SATcorr_STKcut1200_PSDcutXY.root");

   TH1F *h01 = (TH1F*)f->Get("h01");
   TH1F *h02 = (TH1F*)f->Get("h02");
   TH1F *h03 = (TH1F*)f->Get("h03");
   TH1F *h04 = (TH1F*)f->Get("h04");
   TH1F *h05 = (TH1F*)f->Get("h05");
   TH1F *h06 = (TH1F*)f->Get("h06");
   TH1F *h07 = (TH1F*)f->Get("h07");
   TH1F *h08 = (TH1F*)f->Get("h08");
   TH1F *h09 = (TH1F*)f->Get("h09");
   TH1F *h10 = (TH1F*)f->Get("h10");

   

	gStyle->SetOptStat(1111);
	gStyle->SetOptFit(111);

	int j;
	cout << " valore di j (=0 inizializzazione, =1-10 fit histo)" << endl;
	cin >> j ;

	TCanvas *c1=new TCanvas("c1","gaussian",800,800);// c1->Divide(4,3);

	Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
	Double_t fr[2];


	/////////////////// FIT RANGES DEFINITION
   
   if(j==1) {fr[0]=10.3; fr[1]=11.0;}
   if(j==2) {fr[0]=10.35; fr[1]=11.0;}
   if(j==3) {fr[0]=10.3; fr[1]=11.1;}
   if(j==4) {fr[0]=10.35; fr[1]=11.2;}
   if(j==5) {fr[0]=10.45; fr[1]=11.2;}
   if(j==6) {fr[0]=10.5; fr[1]=11.2;}
   if(j==7) {fr[0]=10.4; fr[1]=11.4;}
   if(j==8) {fr[0]=10.5; fr[1]=11.5;}
   if(j==9) {fr[0]=10.45; fr[1]=11.9;}
   if(j==10){fr[0]=10.5; fr[1]=12.2;}



   ////////////////// HISTOGRAM 
   TH1F *xxx = nullptr;

   if(j==1) xxx = (TH1F*)h01->Clone();
   if(j==2) xxx = (TH1F*)h02->Clone();
   if(j==3) xxx = (TH1F*)h03->Clone();
   if(j==4) xxx = (TH1F*)h04->Clone();
   if(j==5) xxx = (TH1F*)h05->Clone();
   if(j==6) xxx = (TH1F*)h06->Clone();
   if(j==7) xxx = (TH1F*)h07->Clone();
   if(j==8) xxx = (TH1F*)h08->Clone();
   if(j==9) xxx = (TH1F*)h09->Clone();
   if(j==10) xxx = (TH1F*)h10->Clone();


	xxx->Fit("landau","","",fr[0],fr[1]);

	cout << "Fitting...\n" << endl;

	// Qui la Langauss eredita i parametri di partenza direttamente dalla Landau 

	TF1 *fit =(TF1 *)xxx->GetFunction("landau");

	plhi[0] = fit->GetParameter(2);
	plhi[1] = fit->GetParameter(1);
	plhi[2] = 0.3*fit->GetParameter(0); //0.1 fino a h03, 0.3 da h04 
	plhi[3] = 1.5*fit->GetParameter(2);   // GSIGMA (sigma della gaussiana) = 0.9 sigma della landau

	for ( int i=0; i<4; i++) { pllo[i]=(0.2)*plhi[i]; sv[i]=plhi[i]; }
   /*
   sv[0] = fit->GetParameter(2); // Width
   sv[1] = fit->GetParameter(1); // MPV
   sv[2] = 0.3 * fit->GetParameter(0);
   sv[3] = 1.2 * fit->GetParameter(2);

   pllo[0] = 0.01;  plhi[0] = 0.06;   // Width
   pllo[1] = sv[1] - 0.5;  plhi[1] = sv[1] + 0.5; // MPV
   pllo[2] = 0.1 * sv[2];  plhi[2] = 10.0 * sv[2]; // Area
   pllo[3] = 0.1;   plhi[3] = 1.5;    // GSigma
   */
	/////////////////////////////////////////////////////////////////////////////////////////////////////

   Double_t chisqr; 
   Int_t    ndf;
   
   TF1 *fitsnr = langaufit(xxx,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);

   Double_t SNRPeak, SNRFWHM;
   langaupro(fp,SNRPeak,SNRFWHM);
   
   printf("Fitting done\nPlotting results...\n");

	/////////////////////// 2 BOX CON I PARAMETRI DEL FIT///////////////////////////////////

	// FIT

  	//xxx->Fit("landau"); xxx->GetFunction("landau")->SetLineColor(2);
    xxx->GetYaxis()->SetTitleOffset(1.3); xxx->GetXaxis()->SetRangeUser(10.,12.3);
  	TH1F *hxx = (TH1F*)(xxx->Clone());
   	
  	TF1* langaufun = (TF1*)(xxx->FindObject("langaufun"));
  	fitsnr->SetLineColor(kBlack);
   hxx->SetMarkerStyle(20); hxx->SetMarkerSize(0.6);
  	hxx->Draw("same");
   fitsnr->Draw("same"); fitsnr->SetLineWidth(3);
  	c1->Update();
 	 TPaveStats *stat = (TPaveStats*)(xxx->FindObject("stats"));
 	 TPaveStats *stat1 = (TPaveStats*)(hxx->FindObject("stats"));
  	if(stat && stat1) { 
  		stat->SetTextColor(8);  
		stat1->SetTextColor(8);
		float height = stat1->GetY2NDC() - stat1->GetY1NDC();
		stat1->SetY1NDC(stat->GetY1NDC() - height);
		stat1->SetY2NDC(stat->GetY1NDC() );
		stat1->Draw(); 
	}

   cout << "Chisqr/ndf: "<< chisqr/ndf << "\n" << endl;


   ///////////////// WRITE FIT PARAMETERS ON A FILE
/**/
   ofstream fout;
   fout.open("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/BACKGROUND/langaus_results/langaus_fit_MC_NEON_STKcut1200_PSDcutXY.txt",ios::app);
   fout << j<< " " << fp[0] << " " << fpe[0] << " " << fp[1] << " " << fpe[1] << " " << fp[2] << " " << fpe[2] << " " << fp[3] << " " << fpe[3] << " " << chisqr/ndf << endl;
   fout.close();
/**/	     

return; 

}
