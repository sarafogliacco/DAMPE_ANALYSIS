{
    /////  HIST ////
    TFile aa("out/CHARGE_eff_MC_SILICON_SATcorr_STKcut1200_PSDcutXY_4bin.root");
    aa.ls();

    /////////////////////////////////////
    int nbin_data = hCharge_Glob->GetNbinsX();

    std::vector<double> nnum(nbin_data);
    std::vector<double> nden(nbin_data);
    std::vector<double> W(nbin_data);
    std::vector<double> DEN(nbin_data);
    std::vector<double> epsilon(nbin_data);  
    // NNUM E NDEN SONO DOUBLE (NON INT COME ERA STATO SCRITTO ALL'INIZIO (PER VERIFICA SI OSSERVINO GLI ISTOGRAMMI DEI FILE IN INPUT

    for(int i=1; i<nbin_data+1; i++) {  // NEGLI ISTOGRAMMI IN INPUT hCharge_Glob & hCharge_XY I BIN SONO 33 PER CUI IL CONTEGGIO CORRETTO VA DA i = 1 A i = 33 (*)
        epsilon[i] = h1_PSD_Eff->GetBinContent(i);
        W[i] = hCharge_Glob->GetBinError(i);
        DEN[i] = hCharge_XY->GetBinContent(i);
        nnum[i]=(DEN[i]/W[i])*(DEN[i]/W[i])*(hCharge_XY->GetBinContent(i));
        nden[i]=(DEN[i]/W[i])*(DEN[i]/W[i])*(hCharge_Glob->GetBinContent(i));
        cout << i << " " << nnum[i] << " " << nden[i] << endl; 
    }

    //canvas only needed for the documentation
    TCanvas* c1 = new TCanvas("c1","",600,400);
    c1->Divide(2);
    c1->SetFillStyle(1001); c1->SetFillColor(kWhite);

    //create one-dimensional TEfficiency object with fixed bin size
    TEfficiency* pEff = new TEfficiency("eff"," ; ; PSD Charge Efficiency",32,2.,6.);
    bool bPassed;
    double x;
    double r = 0.;

    for(int i=0; i<16; i++) { 
        x = 1.+ i*4/16.0;    
        for(int j=0; j<(DEN[i+1]/W[i+1])*(DEN[i+1]/W[i+1]); j++) { // SI INDICIZZANO CON  [i+1] A CAUSA DEL MOTIVO ESPOSTO IN (*)
            bPassed = j < epsilon[i+1]*(DEN[i+1]/W[i+1])*(DEN[i+1]/W[i+1]);
            pEff->Fill(bPassed,x); 
            r = pEff->GetEfficiency(i);
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

    TFile *f = TFile::Open("out/eff_CHARGE_MC_SILICON_4bin.root", "RECREATE");
    f->cd();
    pEff->Write(); Efficiency->Write(); // save the histogram

    c1->Write();
    f->ls();      // show the contents of the ROOT file
    delete f;     // close the ROOT file

}