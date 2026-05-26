
void ratio_STK_charge_efficiency_SULFUR() {

    //==================== STYLE ====================
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadLeftMargin(0.08);
    gStyle->SetPadRightMargin(0.03);
    gStyle->SetPadTopMargin(0.03);
    gStyle->SetPadBottomMargin(0.0);

    //==================== FILE ====================
    TFile* fDATA = TFile::Open("out_eff/eff_DATA_SULFUR_STK2nd_layer_8bin.root");
    TFile* fMC   = TFile::Open("out_eff/eff_MC_SULFUR_STK2nd_layer_8bin.root");

    TEfficiency* effDATA = (TEfficiency*) fDATA->Get("eff");
    TEfficiency* effMC   = (TEfficiency*) fMC->Get("eff");

    TH1F* hDATA = (TH1F*) fDATA->Get("Efficiency");
    TH1F* hMC   = (TH1F*) fMC->Get("Efficiency");

    //==================== CANVAS ====================
    TCanvas* c = new TCanvas("c","HET efficiency SULFUR",1000,950);
    c->SetGridx();
    c->SetGridy();

    //==================== PAD1 (TEfficiency) ====================
    TPad* pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0.0); // niente etichette X
    pad1->SetTopMargin(0.04);
    pad1->SetGridx();
    pad1->SetGridy();
    pad1->Draw();
    

    // Disegna un frame per fissare gli assi
    auto frame1 = pad1->DrawFrame(2., 0.35, 5., 1.017); // xmin=2, xmax=5
    frame1->SetTitle("");
    //frame1->GetYaxis()->SetTitle("HET efficiency");
    frame1->GetYaxis()->CenterTitle();
    frame1->GetYaxis()->SetLabelSize(0.023);
    //frame1->GetXaxis()->SetLabelSize(0); // niente X in alto

    // Metti titolo Y con TLatex
    TLatex latexY;
    latexY.SetTextAngle(90); // ruota di 90 gradi
    latexY.SetTextFont(42);
    latexY.SetTextSize(0.03);
    latexY.DrawLatexNDC(0.027,0.5,"STK charge efficiency"); // 0.03 = offset a sinistra, 0.5 = centro verticale

    // Poi disegna le efficienze sopra il frame
    Int_t colDATA = TColor::GetColor("#b30000");
    effDATA->SetLineColor(kRed+1);
    effDATA->SetMarkerColor(kRed+1);
    effDATA->SetMarkerStyle(20);
    effDATA->SetMarkerSize(1.2);
    effDATA->SetLineWidth(2);
    effDATA->Draw("P same");

    Int_t colMC = TColor::GetColor("#0bb4ff");
    effMC->SetLineColor(kBlue+1);
    effMC->SetMarkerColor(kBlue+1);
    effMC->SetMarkerStyle(24);
    effMC->SetMarkerSize(1.2);
    effMC->SetLineWidth(2);
    effMC->Draw("P same");

    // Legend
    TLegend *leg = new TLegend(0.2,0.55,0.45,0.75);
    leg->SetLineColor(1);
    leg->SetFillColor(0);
    leg->SetTextSize(0.03);
    leg->SetHeader("#scale[1.3]{#bf{SULFUR}}","C");
    leg->AddEntry(effDATA,"FD","pel");
    leg->AddEntry(effMC,"MC","pel");
    leg->SetLineColor(1);
    leg->SetLineWidth(2);
    leg->Draw();


    //==================== PAD2 (TH1F ratio) ====================
    c->cd();
    TPad* pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.3);
    pad2->SetGridx();
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    TH1F* hRatio = (TH1F*) hDATA->Clone("hRatio");
    hRatio->Divide(hMC);
    hRatio->SetMarkerStyle(21);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetLineColor(kBlack);
    hRatio->SetMinimum(0.92);
    hRatio->SetMaximum(1.08);
    hRatio->SetStats(0);

    hRatio->Draw("EP");

    // 
    TF1* lUp = new TF1("lUp","1.032",hDATA->GetXaxis()->GetXmin(),hDATA->GetXaxis()->GetXmax());
    TF1* lDn = new TF1("lDn","0.968",hDATA->GetXaxis()->GetXmin(),hDATA->GetXaxis()->GetXmax());
    lUp->SetLineColor(kMagenta); lUp->SetLineStyle(2); lUp->SetLineWidth(2);
    lDn->SetLineColor(kMagenta); lDn->SetLineStyle(2); lDn->SetLineWidth(2);
    lUp->Draw("same"); lDn->Draw("same");

    // Assi ratio
    hRatio->GetYaxis()->SetLabelSize(0.076);
    hRatio->GetYaxis()->SetNdivisions(505);

    hRatio->GetXaxis()->SetTitle("Log(E_{dep} / GeV)");
    hRatio->GetXaxis()->CenterTitle();
    hRatio->GetXaxis()->SetTitleSize(0.10);
    hRatio->GetXaxis()->SetLabelSize(0.076);
    hRatio->GetXaxis()->SetRangeUser(2.,5.);

    // Metti titolo Y con TLatex
    TLatex latexY_2;
    latexY_2.SetTextAngle(90); // ruota di 90 gradi
    latexY_2.SetTextFont(42);
    latexY_2.SetTextSize(0.08);
    latexY_2.DrawLatexNDC(0.027,0.45," FD /  MC"); // 0.03 = offset a sinistra, 0.5 = centro verticale

    //==================== SAVE ====================
    c->SaveAs("plots/STKcharge_eff_S_DATA_MC_8bin.png");
    //c->SaveAs("STKcharge_eff_S_DATA_MC.pdf");
    //c->SaveAs("STKcharge_eff_S_DATA_MC.eps");


}
