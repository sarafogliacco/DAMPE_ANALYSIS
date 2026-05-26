{
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   gStyle->SetPadLeftMargin(0.08);
   gStyle->SetPadRightMargin(0.03);
   gStyle->SetPadTopMargin(0.03);
   gStyle->SetPadBottomMargin(0.0);


/////HIST DATA MC////
 TFile aa("out/eff_CHARGE_DATA_SULFUR_4bin.root");
TFile bb("out/eff_CHARGE_MC_SILICON_4bin.root");

aa.ls();
bb.ls();

TEfficiency *hist1 = (TEfficiency*)aa.Get("eff");
TEfficiency *hist2 = (TEfficiency*)bb.Get("eff");

TH1F *hist11 = (TH1F*)aa.Get("Efficiency");
TH1F *hist22 = (TH1F*)bb.Get("Efficiency");

// Define the Canvas
   TCanvas *c = new TCanvas("c", "canvas", 700, 700);
	
// Upper plot will be in pad1
   TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);

   pad1->SetBottomMargin(0); // Upper and lower plot are joined
   pad1->SetBottomMargin(1.);
   pad1->SetGridx();         // Vertical grid
   pad1->SetGridy();         // horizontal grid
   //pad1->SetLogx();
   pad1->Draw();             // Draw the upper pad: pad1
   pad1->cd();               // pad1 becomes the current pad
  // hist1->SetStats(0);          // No statistics on upper plot
   hist1->Draw();               // Draw h1
   hist2->Draw("same");         // Draw h2 on top of h1

/////LEGEND///

   TLegend *legend = new TLegend(0.2,0.30,0.42,0.5);
   //legend->SetHeader("   Sulfur ");
   legend->AddEntry(hist1,"Flight data S","ple");
   legend->AddEntry(hist2,"MC Si ","lep");
   legend->SetLineWidth(2);
   legend->Draw();

   TLatex *title = new TLatex();
   title->SetTextAlign(13);  title->SetTextSize(0.06); 
   title->DrawLatex(2.25,0.85,"Sulfur");


   gPad->Update();
   auto graph = hist1->GetPaintedGraph();
   graph->SetMinimum(0.45);
   graph->SetMaximum(1.02);
   gPad->Update();

/////HIST DATA SETTINGS///////////////
   TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
   axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   axis->SetLabelSize(16);
   axis->Draw();

// lower plot will be in pad

   c->cd();          // Go back to the main canvas before defining pad2
   TPad *pad2 = new TPad("pad2", "pad2",  0, 0.0, 1, 0.3);
   pad2->SetTitle("DATA S/MC Si HET eff");
   pad2->SetTopMargin(0.0);
   pad2->SetBottomMargin(0.27);
   pad2->SetGridx(); // vertical grid
   pad2->SetGridy(); // horizontal grid
   //pad2->SetLogx();
   pad2->Draw();
   pad2->cd();       // pad2 becomes the current pad

/*
   TFile aa("out/eff_CHARGE_DATA_SULFUR.root");
	TH1F *hist11 = new TH1F("hist11", "", 32., 2., 6.);
	hist11 = Efficiency;
          hist11->SetTitle("High Energy Trigger efficiency");

/////HIST MC////

  TFile bb("out/eff_CHARGE_MC_SILICON.root");
	TH1F *hist22 = new TH1F("hist22", "", 32., 2., 6.);
	hist22 = Efficiency;
*/
/////HIST RATIO DATA/MC///////////////
   TH1F *h3 = (TH1F*)hist11->Clone("ratio");
   TF1 *cost = new TF1("cost", "1", 0., 10000.);
   cost->SetLineColor(kBlack);
   cost->SetLineStyle(2);
   cost->SetLineWidth(2);
   h3->SetLineColor(kBlack);
   h3->SetMinimum(0.76);  // Define Y ..
   h3->SetMaximum(1.24); // .. range
   h3->Sumw2();
   h3->SetStats(0);      // No statistics on lower plot
   h3->Divide(hist22);
   h3->SetMarkerStyle(21);
   h3->Draw("ep");       // Draw the ratio plot
   //cost->Draw("same");
   TF1 *cost1 = new TF1("cost1", "1.032", 0., 10000.);
   cost1->SetLineColor(TColor::GetColor("#dc0ab4"));
   cost1->SetLineStyle(2);
   cost1->SetLineWidth(2);
   cost1->Draw("same");
   TF1 *cost2 = new TF1("cost2", "0.968", 0., 10000.);
   cost2->SetLineColor(TColor::GetColor("#dc0ab4"));
   cost2->SetLineStyle(2);
   cost2->SetLineWidth(2);
   cost2->Draw("same");

// HIST DATA settings
   Int_t colDATA = TColor::GetColor("#b30000");
   hist1->SetLineColor(colDATA);
   hist1->SetMarkerColor(colDATA);
   hist1->SetMarkerStyle(20);
   hist1->SetMarkerSize(1.);
   hist1->SetLineWidth(2);

// HIST MC settings
   Int_t colMC = TColor::GetColor("#0bb4ff");
   hist2->SetLineColor(colMC);
   hist2->SetMarkerColor(colMC);
   hist2->SetMarkerStyle(20);
   hist2->SetMarkerSize(1.);
   hist2->SetLineWidth(2);

// HIST RATIO settings
   h3->SetTitle(""); // Remove the ratio title
   h3->GetYaxis()->SetTitleFont(43);
   h3->GetYaxis()->SetTitle("Ratio Data/MC");
   h3->GetXaxis()->SetTitle("Log(E_{BGO}/GeV)");
   h3->GetXaxis()->SetNdivisions(505);
   h3->GetYaxis()->SetNdivisions(505);
   h3->GetYaxis()->SetTitleSize(16);
   
   h3->GetYaxis()->SetTitleOffset(1.55);
   h3->GetYaxis()->SetLabelFont(43); 
   h3->GetYaxis()->SetLabelSize(16);

   h3->GetXaxis()->SetTitleSize(16);
   h3->GetXaxis()->SetTitleFont(43);
   h3->GetXaxis()->SetTitleOffset(4.);
   h3->GetXaxis()->SetLabelFont(43);
   h3->GetXaxis()->SetLabelSize(16);
   h3->GetXaxis()->SetRangeUser(2.,4.3);


   //c.SaveAs("../PLOTS/png_FILES/charge_eff_silicon_PhDthesis_v2.png");
   //c.SaveAs("../PLOTS/eps_FILES/charge_eff_silicon_PhDthesis_v2.eps");
}