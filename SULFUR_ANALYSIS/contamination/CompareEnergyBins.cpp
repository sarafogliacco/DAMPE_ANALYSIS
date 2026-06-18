#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>

#include <iostream>
#include <vector>

void CompareEnergyBins(const char* file1Name="MC_FTFP_SILICON_SATcorr_STKcut1200_PSDcutXY_CORRECTED.root",
                       const char* file2Name="../out_root/MC_FTFP_SULFUR_SATcorr_STKcut1200_PSDcutXY_CORRECTED_100GeV-500TeV_v4.root")
{
    TFile *f1 = TFile::Open(file1Name);
    TFile *f2 = TFile::Open(file2Name);

     if(!f1 || !f2) return;

    std::vector<TString> histNames =
    {
        "h01","h02","h03","h04","h05",
        "h06","h07","h08","h09","h10"
    };

    // UNA SOLA CANVAS
    TCanvas *c = new TCanvas("c","Comparison",900,700);

    while(true)
    {
        int idx;

        std::cout << "\nEnergia [1-10] (-1 esce): ";
        std::cin >> idx;

        if(idx==-1)
            break;

        if(idx<1 || idx>10)
            continue;

        TH1F *hA = (TH1F*)f1->Get(histNames[idx-1]);
        TH1F *hB = (TH1F*)f2->Get(histNames[idx-1]);

        if(!hA || !hB)
        {
            std::cout << "Istogramma mancante\n";
            continue;
        }

        c->Clear();

        hA->SetLineColor(kRed);
        hB->SetLineColor(kBlue);

        hA->Draw("hist");
        hB->Draw("hist same");

        auto leg = new TLegend(0.65,0.75,0.88,0.88);
        leg->AddEntry(hA,"Si ","l");
        leg->AddEntry(hB,"S","l");
        leg->Draw();

        c->Modified();
        c->Update();
    }
}