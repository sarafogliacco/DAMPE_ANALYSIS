{

	gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadLeftMargin(0.08);
    gStyle->SetPadRightMargin(0.03);
    gStyle->SetPadTopMargin(0.03);
    gStyle->SetPadBottomMargin(0.08);
	//gStyle->SetTextFont(4);

	TCanvas *can = new TCanvas("can","can",900,780); can->SetLogy();

	TString colCode[9] = { "#b30000", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0" };
    Int_t col[9];
    for(int ii=0; ii<9; ii++){ col[ii]=TColor::GetColor(colCode[ii]);} 

	TFile aa("PSDcharge_perCUT_DATA_120months_SULFUR.root");
	TH1F *hcut00 = hqpsdPathWeight_cut00->Clone(); hcut00->SetStats(0);
	hcut00->SetLineWidth(3); hcut00->SetLineColor(kBlack); 
	hcut00->GetYaxis()->SetRangeUser(200.,1500000.);
	hcut00->GetXaxis()->SetTitle("PSD charge"); hcut00->GetYaxis()->SetTitle("events");
	hcut00->Draw("");

	TH1F *hcut03 = hqpsdPathWeight_cut03->Clone(); 
	hcut03->SetLineWidth(3); hcut03->SetLineColor(col[0]);
	hcut03->Draw("same");

	TH1F *hstkcut = hqpsdPathWeight_stkcut->Clone(); 
	hstkcut->SetLineWidth(3); hstkcut->SetLineColor(col[1]);
	hstkcut->Draw("same");

	TH1F *hpsdcut = hqpsdPathWeight_psdcut->Clone(); 
	hpsdcut->SetLineWidth(3); hpsdcut->SetLineColor(col[2]);
	hpsdcut->Draw("same");
	
	auto legend = new TLegend(0.51,0.70,0.84,0.94);
	legend->AddEntry(hcut00, "Pre-selection + HET","l");
	legend->AddEntry(hcut03,"PSD hits on YZ and XZ views","l");
	legend->AddEntry(hstkcut,"STK 1st layer > 1200 ADC","l");
	legend->AddEntry(hpsdcut,"|PSD_{Y} - PSD_{X}| < 2","l");
	legend->SetTextSize(0.035);
	legend->SetLineColor(0);
	legend->SetFillStyle(0);
	legend->Draw();

	/*
	tit = TLatex(); tit.SetTextAlign(13);
    tit.DrawLatex(5.5,1050000.,"Preliminary");
    tit.SetTextColor(1); //tit.SetTextSize(0.05);
    tit.Draw();
    //tit.SetTextFont(4);

    titH = TLatex(); titH.DrawLatex(0.7,600000.,"H");
    titH.SetTextColor(1); titH.SetTextSize(0.02); titH.Draw();

    titHe = TLatex(); titHe.DrawLatex(1.7,380000.,"He");
    titHe.SetTextColor(1); titHe.SetTextSize(0.02); titHe.Draw();

    titLi = TLatex(); titLi.DrawLatex(2.7,38000.,"Li");
    titLi.SetTextColor(1); titLi.SetTextSize(0.02); titLi.Draw();

    titBe = TLatex(); titBe.DrawLatex(3.7,24000.,"Be");
    titBe.SetTextColor(1); titBe.SetTextSize(0.02); titBe.Draw();

    titB = TLatex(); titB.DrawLatex(4.7,30000.,"B");
    titB.SetTextColor(1); titB.SetTextSize(0.02); titB.Draw();

    titC = TLatex(); titC.DrawLatex(5.7,122000.,"C");
    titC.SetTextColor(1); titC.SetTextSize(0.02); titC.Draw();

    titN = TLatex(); titN.DrawLatex(6.7,42000.,"N");
    titN.SetTextColor(1); titN.SetTextSize(0.02); titN.Draw();

    titO = TLatex(); titO.DrawLatex(7.7,168000.,"O");
    titO.SetTextColor(1); titO.SetTextSize(0.02); titO.Draw();

    titF = TLatex(); titF.DrawLatex(9.,30000.,"F");
    titF.SetTextColor(1); titF.SetTextSize(0.02); titF.Draw();

    titNe = TLatex(); titNe.DrawLatex(9.8,37300.,"Ne");
    titNe.SetTextColor(1); titNe.SetTextSize(0.02); titNe.Draw();

    titNa = TLatex(); titNa.DrawLatex(10.7,15000.,"Na");
    titNa.SetTextColor(1); titNa.SetTextSize(0.02); titNa.Draw();

    titMg = TLatex(); titMg.DrawLatex(11.7,48000.,"Mg");
    titMg.SetTextColor(1); titMg.SetTextSize(0.02); titMg.Draw();

    titAl = TLatex(); titAl.DrawLatex(12.7,17000.,"Al");
    titAl.SetTextColor(1); titAl.SetTextSize(0.02); titAl.Draw();

    titSi = TLatex(); titSi.DrawLatex(13.7,46500.,"Si");
    titSi.SetTextColor(1); titSi.SetTextSize(0.02); titSi.Draw();

    titS = TLatex(); titS.DrawLatex(15.,9500.,"S");
    titS.SetTextColor(1); titS.SetTextSize(0.02); titS.Draw();

    titFe = TLatex(); titFe.DrawLatex(25.7,60000.,"Fe");
    titFe.SetTextColor(1); titFe.SetTextSize(0.02); titFe.Draw();
*/


	//can.SaveAs("PLOTS/png_FILES/PSDcharge_perCUT_DATA.png");
	//can.SaveAs("PLOTS/PSDcharge_perCUT_DATA.pdf");
	//can.SaveAs("PLOTS/eps_FILES/PSDcharge_perCUT_DATA.eps");




}
