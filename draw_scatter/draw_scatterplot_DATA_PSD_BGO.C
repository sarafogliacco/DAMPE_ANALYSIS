{
    gStyle->SetPalette(55);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadLeftMargin(0.09);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadTopMargin(0.03);
    gStyle->SetPadBottomMargin(0.09);

    TFile aa("/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/scatterplot_DATA_120months_PSD_BGO.root");

    TH2F *hDATA = (TH2F*)aa.Get("h2psdbgo"); 

    if (!hDATA) {
        std::cerr << "Missing histogram!" << std::endl;
        aa.ls();  // helpful debug
        return;
    }

    hDATA->SetStats(0);
    hDATA->GetXaxis()->SetTitleOffset(1.2);
    hDATA->GetYaxis()->SetTitle("Q^{PSD}");
    hDATA->GetYaxis()->SetTitleOffset(1.1);
    hDATA->GetYaxis()->SetRangeUser(15,22);
    hDATA->GetXaxis()->SetRangeUser(95.,67000.);

    TCanvas *cc = new TCanvas("cc","cc",780,780);
    cc->SetLogx();
    cc->SetLogz();

    hDATA->Draw("colz");

    //TString cut_dw = "2.0"; TString cut_up = "3.0"; 

    /////////////////////////////////////////////////////////////////////////////////
/*
    // NEON
    TString MPV_DATA_ne = "(10.0182+(-0.0839296*log10(x))+(0.0451273*pow(log10(x),2))+(-0.00272286*pow(log10(x),3)))";
    TString sigma_DATA_ne = "(0.233871+(-0.0557627*log10(x))+(0.0206737*pow(log10(x),2)))";

    TString cut_dw_ne = "2.0"; TString cut_up_ne = "3.0";

    TF1 *mpv_ne =  new TF1("mpv_ne",MPV_DATA_ne,100.,60000);
    mpv_ne->SetLineColor(1); mpv_ne->SetLineWidth(4); mpv_ne->Draw("same");
    TF1 *slow_ne = new TF1("slow_ne",MPV_DATA_ne+"-"+cut_dw_ne+"*"+sigma_DATA_ne,100.,60000);
    slow_ne->SetLineColor(1); slow_ne->SetLineWidth(4); slow_ne->SetLineStyle(9); slow_ne->Draw("same");
    TF1 *sup_ne  = new TF1("sup_ne", MPV_DATA_ne+"+"+cut_up_ne+"*"+sigma_DATA_ne,100.,60000);
    sup_ne->SetLineColor(1); sup_ne->SetLineWidth(4); sup_ne->SetLineStyle(9); sup_ne->Draw("same");

    TLatex *lat_Ne_low = new TLatex(120.,9.7,"Ne MPV - 2#sigma"); lat_Ne_low->SetTextSize(0.027); lat_Ne_low->Draw();
    TLatex *lat_Ne     = new TLatex(220.,10.2,"Ne MPV"); lat_Ne->SetTextSize(0.027); lat_Ne->Draw();
    TLatex *lat_Ne_up  = new TLatex(120.,10.6,"Ne MPV + 3#sigma"); lat_Ne_up->SetTextSize(0.027);  lat_Ne_up->Draw();
    
    // MAGNESIUM 
    TString MPV_DATA_mg = "(11.5043+(0.584417*log10(x))+(-0.203462*pow(log10(x),2))+(0.0273119*pow(log10(x),3)))";
    TString sigma_DATA_mg = "(0.378608+(-0.127283*log10(x))+(0.0285267*pow(log10(x),2)))";

    TString cut_dw_mg = "2.0"; TString cut_up_mg = "3.0";

    TF1 *mpv_mg =  new TF1("mpv_mg",MPV_DATA_mg,100.,60000);
    mpv_mg->SetLineColor(1); mpv_mg->SetLineWidth(4); mpv_mg->Draw("same");
    TF1 *slow_mg = new TF1("slow_mg",MPV_DATA_mg+"-"+cut_dw_mg+"*"+sigma_DATA_mg,100.,60000);
    slow_mg->SetLineColor(1); slow_mg->SetLineWidth(4); slow_mg->SetLineStyle(9); slow_mg->Draw("same");
    TF1 *sup_mg  = new TF1("sup_mg", MPV_DATA_mg+"+"+cut_up_mg+"*"+sigma_DATA_mg,100.,60000);
    sup_mg->SetLineColor(1); sup_mg->SetLineWidth(4); sup_mg->SetLineStyle(9); sup_mg->Draw("same");

    TLatex *lat_Mg_low = new TLatex(120.,11.7,"Mg MPV - 2#sigma"); lat_Mg_low->SetTextSize(0.027); lat_Mg_low->Draw();
    TLatex *lat_Mg     = new TLatex(220.,12.2,"Mg MPV"); lat_Mg->SetTextSize(0.027); lat_Mg->Draw();
    TLatex *lat_Mg_up  = new TLatex(120.,12.6,"Mg MPV + 3#sigma"); lat_Mg_up->SetTextSize(0.027);  lat_Mg_up->Draw();

    // SILICON *** 
	
    TString MPV_DATA_si = "(14.3409+(-0.313026*log10(x))+(0.122946*pow(log10(x),2))+(-0.0124761*pow(log10(x),3)))";
    TString sigma_DATA_si = "(0.315991+(-0.0518163*log10(x))+(0.0134847*pow(log10(x),2)))";

    TString cut_dw_si = "2.0"; TString cut_up_si = "3.0";

    TF1 *mpv_si =  new TF1("mpv_si",MPV_DATA_si,100.,60000);
    mpv_si->SetLineColor(1); mpv_si->SetLineWidth(4); mpv_si->Draw("same");
    TF1 *slow_si = new TF1("slow_si",MPV_DATA_si+"-"+cut_dw_si+"*"+sigma_DATA_si,100.,60000);
    slow_si->SetLineColor(1); slow_si->SetLineWidth(4); slow_si->SetLineStyle(9); slow_si->Draw("same");
    TF1 *sup_si  = new TF1("sup_si", MPV_DATA_si+"+"+cut_up_si+"*"+sigma_DATA_si,100.,60000);
    sup_si->SetLineColor(1); sup_si->SetLineWidth(4); sup_si->SetLineStyle(9); sup_si->Draw("same");

    TLatex *lat_Si_low = new TLatex(120.,13.7,"Si MPV - 2#sigma"); lat_Si_low->SetTextSize(0.027); lat_Si_low->Draw();
    TLatex *lat_Si     = new TLatex(220.,14.2,"Si MPV"); lat_Si->SetTextSize(0.027); lat_Si->Draw();
    TLatex *lat_Si_up  = new TLatex(120.,14.6,"Si MPV + 3#sigma"); lat_Si_up->SetTextSize(0.027);  lat_Si_up->Draw();
*/


    // SULFUR
    TString MPV_DATA_s =  "(15.0965+(0.926641*log10(x))+(-0.26595*pow(log10(x),2))+(0.0260414*pow(log10(x),3)))";
    TString sigma_DATA_s = "(0.194339+(0.0486462*log10(x)))";


    TString cut_dw_s = "2.0"; TString cut_up_s = "2.5";

    //TF1 *mpv_s =  new TF1("mpv_s",MPV_DATA_s,100.,60000); mpv_s->SetLineColor(1); mpv_s->SetLineWidth(4); mpv_s->Draw("same");
    TF1 *slow_s = new TF1("slow_s",MPV_DATA_s+"-"+cut_dw_s+"*"+sigma_DATA_s,100.,60000); slow_s->SetLineColor(1); slow_s->SetLineWidth(4); slow_s->SetLineStyle(9); slow_s->Draw("same");
    //TF1 *sup_s  = new TF1("sup_s", MPV_DATA_s+"+"+cut_up_s+"*"+sigma_DATA_s,100.,60000); sup_s->SetLineColor(1); sup_s->SetLineWidth(4); sup_s->SetLineStyle(9); sup_s->Draw("same");

    TLatex *lat_S_low = new TLatex(120.,15.6,"S MPV - 2#sigma"); lat_S_low->SetTextSize(0.027); lat_S_low->Draw();
    TLatex *lat_S     = new TLatex(220.,16.25,"S MPV"); lat_S->SetTextSize(0.027); lat_S->Draw();
    //TLatex *lat_S_up  = new TLatex(120.,16.99,"S MPV + 2.5#sigma"); lat_S_up->SetTextSize(0.027);  lat_S_up->Draw();

 
    // ARGON *** multigauss

    //TString MPV_DATA_ar = "(18.0376+(0.046027*log10(x)))";
    //TString sigma_DATA_ar = "(0.291262+(0.00339476*log10(x)))";


    
    // ARGON 27/05/2026

    TString MPV_DATA_ar = "(16.762+(1.14769*log10(x))+(-0.264239*pow(log10(x),2))+(0.00423957*pow(log10(x),4)))";
    TString sigma_DATA_ar = "(0.301538+(0.0364786*log10(x)))";
    
    TString cut_dw_ar = "1.0"; TString cut_up_ar = "1.5";

    //TF1 *mpv_ar =  new TF1("mpv_ar",MPV_DATA_ar,100.,60000); mpv_ar->SetLineColor(1); mpv_ar->SetLineWidth(4); mpv_ar->Draw("same");
    //TF1 *slow_ar = new TF1("slow_ar",MPV_DATA_ar+"-"+cut_dw_ar+"*"+sigma_DATA_ar,100.,60000); slow_ar->SetLineColor(1); slow_ar->SetLineWidth(4); slow_ar->SetLineStyle(9); slow_ar->Draw("same");
    //TF1 *sup_ar  = new TF1("sup_ar", MPV_DATA_ar+"+"+cut_up_ar+"*"+sigma_DATA_ar,100.,60000); sup_ar->SetLineColor(1); sup_ar->SetLineWidth(4); sup_ar->SetLineStyle(9); sup_ar->Draw("same");

    //TLatex *lat_Ar_low = new TLatex(120.,17.8,"Ar MPV - 1#sigma"); lat_Ar_low->SetTextSize(0.027); lat_Ar_low->Draw();
    TLatex *lat_Ar     = new TLatex(220.,18.2,"Ar MPV"); lat_Ar->SetTextSize(0.027); lat_Ar->Draw();
    //TLatex *lat_Ar_up  = new TLatex(120.,18.5,"Ar MPV + 1.5#sigma"); lat_Ar_up->SetTextSize(0.027);  lat_Ar_up->Draw();


    
    // CALCIUM *** multigauss
    TString MPV_DATA_ca = "(17.1818+(2.97939*log10(x))+(-0.976515*pow(log10(x),2))+(0.10604*pow(log10(x),3)))";
    TString sigma_DATA_ca = "(0.206041+(0.0693299*log10(x)))";
   

    TString cut_dw_ca = "1.0"; TString cut_up_ca = "2.0";
    
    //TF1 *mpv_ca =  new TF1("mpv_ca",MPV_DATA_ca,100.,60000);mpv_ca->SetLineColor(1); mpv_ca->SetLineWidth(4); mpv_ca->Draw("same");
    //TF1 *slow_ca = new TF1("slow_ca",MPV_DATA_ca+"-"+cut_dw_ca+"*"+sigma_DATA_ca,100.,60000);slow_ca->SetLineColor(1); slow_ca->SetLineWidth(4); slow_ca->SetLineStyle(9); slow_ca->Draw("same");
    TF1 *sup_ca  = new TF1("sup_ca", MPV_DATA_ca+"+"+cut_up_ca+"*"+sigma_DATA_ca,100.,60000);sup_ca->SetLineColor(1); sup_ca->SetLineWidth(4); sup_ca->SetLineStyle(9); sup_ca->Draw("same");

    //TLatex *lat_Ca_low = new TLatex(120.,19.6,"Ca MPV - 1#sigma"); lat_Ca_low->SetTextSize(0.027); lat_Ca_low->Draw();
    TLatex *lat_Ca     = new TLatex(220.,20.3,"Ca"); lat_Ca->SetTextSize(0.027); lat_Ca->Draw();
    TLatex *lat_Ca_up  = new TLatex(120.,21.,"Ca MPV + 2#sigma"); lat_Ca_up->SetTextSize(0.027);  lat_Ca_up->Draw();
 
/*
    TString low_cut = MPV_DATA_s + "-" + cut_dw_s + "*" + sigma_DATA_s;
    TString up_cut = MPV_DATA_ca + "-" + cut_up_cs + "*" + sigma_DATA_cs;

    
    cout << "low_cut = " << low_cut << endl;
    cout << "low_up = " << low_up << endl;


*/
    cc->SaveAs("plots/scatterplot_DATA_120months_PSD_BGO_group.png");
    //cc->SaveAs("PLOTS/scatterplot_DATA_96months_PSD_BGO_PhDthesis.pdf");
    //cc->SaveAs("PLOTS/eps_FILES/scatterplot_DATA_96months_PSD_BGO_PhDthesis.eps");


}
