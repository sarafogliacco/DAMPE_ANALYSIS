#include <vector>

void efficienza_DATA_SULFUR(){
    /////  HIST ////
    TFile aa("out_eff/HET_eff_DATA_SULFUR_120months_8bin.root"); // sulfur/HET_eff_DATA_SULFUR_108months_6bin.root
    aa.ls();
    TH1F* h1_Unb = (TH1F*) aa.Get("h1_Unb");
    TH1F* h1Nobs_And = (TH1F*) aa.Get("h1Nobs_And");

    /////////////////////////////////////
    const int nbin_data=h1_Unb->GetNbinsX();
    std::vector<double> nnum(nbin_data+1, 0.0);
    std::vector<double> nden(nbin_data+1, 0.0);   
    // NNUM E NDEN SONO DOUBLE (NON INT COME ERA STATO SCRITTO ALL'INIZIO (PER VERIFICA SI OSSERVINO GLI ISTOGRAMMI DEI FILE IN INPUT

    for(int i=1; i<33; i++) {  // NEGLI ISTOGRAMMI IN INPUT h1Nobs_And & h1_Unb I BIN SONO 33 PER CUI IL CONTEGGIO CORRETTO VA DA i = 1 A i = 33 (*)
        nnum[i]=h1Nobs_And->GetBinContent(i);
        nden[i]=h1_Unb->GetBinContent(i);
        cout << i << " " << nnum[i] << " " << nden[i] << endl; 
    }

    //canvas only needed for the documentation
    TCanvas* c1 = new TCanvas("c1","",600,400);
    c1->Divide(2);
    c1->SetFillStyle(1001); c1->SetFillColor(kWhite);

    //create one-dimensional TEfficiency object with fixed bin size
    TEfficiency* pEff = new TEfficiency("eff"," ; ;HET Efficiency",32,2.,6.);
    bool bPassed;
    double x;

    for(int i=0; i<32; i++) { 
        x = 2.+ i*4/32.0;    
        for(int j=0; j<nden[i+1]; j++) { // SI INDICIZZANO CON  [i+1] A CAUSA DEL MOTIVO ESPOSTO IN (*)
            bPassed = j < nnum[i+1];
            pEff->Fill(bPassed,x); 
        }
    }

    double meanvalue, errlow;
    TH1F *Efficiency = new TH1F("Efficiency", "", 32, 2., 6.);

    for(int i=0; i<32; i++) {
        meanvalue = pEff->GetEfficiency(i);
        errlow = pEff->GetEfficiencyErrorLow(i);
        Efficiency->SetBinContent(i,meanvalue);
        Efficiency->SetBinError(i,errlow);

    }

    pEff->Draw();
    Efficiency->Draw("same");

    TFile *f = TFile::Open("out_eff/eff_HET_DATA_120months_SULFUR_8bin.root", "RECREATE");
    f->cd();
    pEff->Write(); Efficiency->Write(); // save the histogram

    c1->Write();
    f->ls();      // show the contents of the ROOT file
    delete f;     // close the ROOT file

}