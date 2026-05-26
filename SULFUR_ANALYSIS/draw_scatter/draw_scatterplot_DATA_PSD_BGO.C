{
    gStyle->SetPalette(55);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadLeftMargin(0.09);
    gStyle->SetPadRightMargin(0.1);
    gStyle->SetPadTopMargin(0.03);
    gStyle->SetPadBottomMargin(0.09);

    TFile aa("out_root/scatterplot_DATA_120months_PSD_BGO.root");

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
    hDATA->GetYaxis()->SetRangeUser(15,21.2);
    hDATA->GetXaxis()->SetRangeUser(95.,67000.);

    TCanvas *cc = new TCanvas("cc","cc",780,780);
    cc->SetLogx();
    cc->SetLogz();

    hDATA->Draw("colz");

    //TString cut_dw = "2.0"; TString cut_up = "3.0"; 

    /////////////////////////////////////////////////////////////////////////////////
    
    // SILICON *** multigauss
    //TString MPV_DATA_si = "(13.0366+(1.01399*log10(x))+(-0.321937*log10(x)**2)+(0.036754*log10(x)**3))";
    //TString sigma_DATA_si = "(-0.52482+(0.801866*log10(x))+(-0.278872*log10(x)**2)+(0.0318674*log10(x)**3))";
    //TString sigma_DATA_si = "(0.205642+(0.00990178*log10(x)))";	
   /* TString MPV_DATA_si = "(13.9252+(0.0535413*log10(x)))";
    TString sigma_DATA_si = "(0.262088+(0.00981284*log10(x)))";

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
    TString MPV_DATA_s =  "(14.5821+(1.51344*log10(x))+(-0.48318*pow(log10(x),2))+(0.0520877*pow(log10(x),3)))";
    //"(15.9434+(0.0730006*log10(x))+(0.00533945*pow(log10(x),2))+(-0.000744557*pow(log10(x),3)))";
    //"(14.9723+(1.09723*log10(x))+(-0.342737*pow(log10(x),2))+(0.0375106*pow(log10(x),3)))"; 
    TString sigma_DATA_s = "(0.20977+(0.0478004*log10(x))+(-0.00269758*pow(log10(x),2)))";
    //"(0.212645+(0.0405252*log10(x)))";
    //"(0.355433+(-0.0688415*log10(x))+(0.019678*pow(log10(x),2)))";
    //"(0.138753+(0.158506*log10(x))+(-0.0559774*pow(log10(x),2))+(0.00789775*pow(log10(x),3)))";;

    TString cut_dw_s = "2.0"; TString cut_up_s = "3.0";

    TF1 *mpv_s =  new TF1("mpv_s",MPV_DATA_s,100.,60000); mpv_s->SetLineColor(1); mpv_s->SetLineWidth(4); mpv_s->Draw("same");
    TF1 *slow_s = new TF1("slow_s",MPV_DATA_s+"-"+cut_dw_s+"*"+sigma_DATA_s,100.,60000); slow_s->SetLineColor(1); slow_s->SetLineWidth(4); slow_s->SetLineStyle(9); slow_s->Draw("same");
    TF1 *sup_s  = new TF1("sup_s", MPV_DATA_s+"+"+cut_up_s+"*"+sigma_DATA_s,100.,60000); sup_s->SetLineColor(1); sup_s->SetLineWidth(4); sup_s->SetLineStyle(9); sup_s->Draw("same");

    TLatex *lat_S_low = new TLatex(120.,15.6,"S MPV - 2#sigma"); lat_S_low->SetTextSize(0.027); lat_S_low->Draw();
    TLatex *lat_S     = new TLatex(220.,16.25,"S "); lat_S->SetTextSize(0.027); lat_S->Draw();
    TLatex *lat_S_up  = new TLatex(120.,16.85,"S MPV + 3#sigma"); lat_S_up->SetTextSize(0.027);  lat_S_up->Draw();

 
    // ARGON *** multigauss
    /*
    TString MPV_DATA_ar = "(18.0376+(0.046027*log10(x)))";
    TString sigma_DATA_ar = "(0.291262+(0.00339476*log10(x)))";

    TString cut_dw_ar = "1.0"; TString cut_up_ar = "1.0";

    TF1 *mpv_ar =  new TF1("mpv_ar",MPV_DATA_ar,100.,60000);
    mpv_ar->SetLineColor(1); mpv_ar->SetLineWidth(4); mpv_ar->Draw("same");
    TF1 *slow_ar = new TF1("slow_ar",MPV_DATA_ar+"-"+cut_dw_ar+"*"+sigma_DATA_ar,100.,60000);
    slow_ar->SetLineColor(1); slow_ar->SetLineWidth(4); slow_ar->SetLineStyle(9); slow_ar->Draw("same");
    TF1 *sup_ar  = new TF1("sup_ar", MPV_DATA_ar+"+"+cut_up_ar+"*"+sigma_DATA_ar,100.,60000);
    sup_ar->SetLineColor(1); sup_ar->SetLineWidth(4); sup_ar->SetLineStyle(9); sup_ar->Draw("same");

    TLatex *lat_Ar_low = new TLatex(120.,17.6,"Ar MPV - 1#sigma"); lat_Ar_low->SetTextSize(0.027); lat_Ar_low->Draw();*/
    TLatex *lat_Ar     = new TLatex(220.,18.2,"Ar"); lat_Ar->SetTextSize(0.027); lat_Ar->Draw();
/*    TLatex *lat_Ar_up  = new TLatex(120.,18.7,"Ar MPV + 1#sigma"); lat_Ar_up->SetTextSize(0.027);  lat_Ar_up->Draw();
    */


    // CALCIUM *** multigauss
    TString MPV_DATA_ca = "(20.0349+(0.0451529*log10(x)))";
    TString sigma_DATA_ca = "(0.400169+(-0.00286795*log10(x)))";

    TString cut_dw_ca = "1.0"; TString cut_up_ca = "2.0";

    //TF1 *mpv_ca =  new TF1("mpv_ca",MPV_DATA_ca,100.,60000);mpv_ca->SetLineColor(1); mpv_ca->SetLineWidth(4); mpv_ca->Draw("same");
    //TF1 *slow_ca = new TF1("slow_ca",MPV_DATA_ca+"-"+cut_dw_ca+"*"+sigma_DATA_ca,100.,60000);slow_ca->SetLineColor(1); slow_ca->SetLineWidth(4); slow_ca->SetLineStyle(9); slow_ca->Draw("same");
    //TF1 *sup_ca  = new TF1("sup_ca", MPV_DATA_ca+"+"+cut_up_ca+"*"+sigma_DATA_ca,100.,60000);sup_ca->SetLineColor(1); sup_ca->SetLineWidth(4); sup_ca->SetLineStyle(9); sup_ca->Draw("same");

    //TLatex *lat_Ca_low = new TLatex(120.,19.8,"Ca MPV - 1#sigma"); lat_Ca_low->SetTextSize(0.027); lat_Ca_low->Draw();
    TLatex *lat_Ca     = new TLatex(220.,20.3,"Ca"); lat_Ca->SetTextSize(0.027); lat_Ca->Draw();
    //TLatex *lat_Ca_up  = new TLatex(120.,20.7,"Ca MPV + 2#sigma"); lat_Ca_up->SetTextSize(0.027);  lat_Ca_up->Draw();



    //cc->SaveAs("plots/scatterplot_DATA_108months_PSD_BGO.png");
    //cc->SaveAs("PLOTS/scatterplot_DATA_96months_PSD_BGO_PhDthesis.pdf");
    //cc->SaveAs("PLOTS/eps_FILES/scatterplot_DATA_96months_PSD_BGO_PhDthesis.eps");


}
