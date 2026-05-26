{
    /////  HIST ////
    TFile aa("out/CHARGE_eff_DATA_108months_SULFUR_SATcorr_STKcut1200_PSDcutXY_4bin.root");
    aa.ls();

    /////////////////////////////////////
    int nbin_data = hCharge_Glob->GetNbinsX();

    std::vector<double> nnum(nbin_data);
    std::vector<double> nden(nbin_data);  
    // NNUM E NDEN SONO DOUBLE (NON INT COME ERA STATO SCRITTO ALL'INIZIO (PER VERIFICA SI OSSERVINO GLI ISTOGRAMMI DEI FILE IN INPUT

    for(int i=1; i<nbin_data+1; i++) {  // NEGLI ISTOGRAMMI IN INPUT hCharge_Glob & h1_Unb I BIN SONO 33 PER CUI IL CONTEGGIO CORRETTO VA DA i = 1 A i = 33 (*)
        nnum[i]=hCharge_Glob->GetBinContent(i);
        nden[i]=hCharge_XY->GetBinContent(i);
        cout << i << " " << nnum[i] << " " << nden[i] << endl; 
    }

    //canvas only needed for the documentation
    TCanvas* c1 = new TCanvas("c1","",600,400);
    c1->Divide(2);
    c1->SetFillStyle(1001); c1->SetFillColor(kWhite);

    //create one-dimensional TEfficiency object with fixed bin size
    TEfficiency* pEff = new TEfficiency("eff"," ; ; #varepsilon_{PSD} ",32,2.,6.);
    bool bPassed;
    double x;

    for(int i=0; i<16; i++) { 
        x = 1.+ i*4/16.0;    
        for(int j=0; j<nden[i+1]; j++) { // SI INDICIZZANO CON  [i+1] A CAUSA DEL MOTIVO ESPOSTO IN (*)
            bPassed = j < nnum[i+1];
            pEff->Fill(bPassed,x); 
        }
    }

    double meanvalue, errlow;
    TH1F *Efficiency = new TH1F("Efficiency", "", 16, 2., 6.);

    for(int i=0; i<16; i++) {
        meanvalue = pEff->GetEfficiency(i);
        errlow = pEff->GetEfficiencyErrorLow(i);
        Efficiency->SetBinContent(i,meanvalue);
        Efficiency->SetBinError(i,errlow);

    }

    pEff->Draw();
    Efficiency->Draw("same");

    TFile *f = TFile::Open("out/eff_CHARGE_DATA_SULFUR_4bin.root", "RECREATE");
    f->cd();
    pEff->Write(); Efficiency->Write(); // save the histogram

    c1->Write();
    f->ls();      // show the contents of the ROOT file
    delete f;     // close the ROOT file

}