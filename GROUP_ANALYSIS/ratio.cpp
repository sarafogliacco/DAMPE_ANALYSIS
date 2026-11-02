#include <TFile.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <algorithm>

double median_energy(double e_min, double e_max, double alpha = 2.7) {
    return TMath::Power((TMath::Power(e_max, 1.-alpha) - TMath::Power(e_min, 1.-alpha)) /
                        ((1.-alpha) * (e_max - e_min)), -1.0/alpha);
}

void fit_ratio_ArS(const char* fAr = "/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/ARGON_ANALYSIS/unfold_results/unfold_results_argon_6bin_low1_up2_v2.root",                                         
                   const char* fS  = "/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/SULFUR_ANALYSIS/unfold_results/unfold_results_flux_spectrum_sulfur_120months_6bin_smooth_2e5up_1e5low_bkg_v2.root",
                   double Efit_min = 1e2, double Efit_max = 1e5,
                   bool draw = true, const char* outname = "ratio_ArS")
{
    bool save = (outname && outname[0] != '\0');
    if (!draw) gROOT->SetBatch(kTRUE);

    TFile* a = TFile::Open(fAr); TFile* s = TFile::Open(fS);
    if (!a || a->IsZombie() || !s || s->IsZombie()) { std::cerr << "Error opening file\n"; return; }
    TH1D* hAr = (TH1D*)a->Get("flux");
    TH1D* hS  = (TH1D*)s->Get("flux");
    TH1D* hArP = (TH1D*)a->Get("flux_pow");
    TH1D* hSP  = (TH1D*)s->Get("flux_pow");
    if (!hAr || !hS) { std::cerr << "Histogram flux not found\n"; return; }


    if (hAr->GetNbinsX() != hS->GetNbinsX()) { std::cerr << "Binning not matching!\n"; return; }
    for (int i = 1; i <= hAr->GetNbinsX() + 1; i++)
        if (TMath::Abs(hAr->GetXaxis()->GetBinLowEdge(i) - hS->GetXaxis()->GetBinLowEdge(i)) > 1e-6) {
            std::cerr << "Bin edges different at bin " << i << "\n"; return; }

    std::vector<double> x, r, er, ex;
    for (int i = 1; i <= hAr->GetNbinsX(); i++) {
        double fa = hAr->GetBinContent(i), ea = hAr->GetBinError(i);
        double fs = hS->GetBinContent(i),  es = hS->GetBinError(i);
        double lo = hAr->GetXaxis()->GetBinLowEdge(i), up = hAr->GetXaxis()->GetBinUpEdge(i);
        double em = median_energy(lo, up);
        if (fa <= 0 || fs <= 0 || em < Efit_min || em > Efit_max) continue;
        double ri = fa / fs;
        double eri = ri * TMath::Sqrt(TMath::Power(ea/fa, 2) + TMath::Power(es/fs, 2));
        x.push_back(em); r.push_back(ri); er.push_back(eri);
        ex.push_back(0.5 * (up - lo));   
    }
    int n = x.size();
    if (n < 2) { std::cerr << "Too few points\n"; return; }

    // E0 = centro in log dell'intervallo, per decorrelare R0 e DGAMMA
    double E0 = TMath::Sqrt(x.front() * x.back());

    // Grafico per il fit (senza errori in x)
    TGraphErrors* g = new TGraphErrors(n);
    for (int i = 0; i < n; i++) { g->SetPoint(i, x[i], r[i]); g->SetPointError(i, 0, er[i]); }

    double xlo = x.front() * 0.9, xhi = x.back() * 1.1;

    // Modello 1: costante (DGAMMA = 0)
    TF1* f0 = new TF1("f0", "[0]*TMath::Power(x/[1], -[2])", xlo, xhi);
    f0->SetParameters(1., E0, 0.);
    f0->FixParameter(1, E0); f0->FixParameter(2, 0.);
    g->Fit(f0, "RN");

    // Modello 2: DGAMMA libero
    TF1* f1 = new TF1("f1", "[0]*TMath::Power(x/[1], -[2])", xlo, xhi);
    f1->SetParameters(f0->GetParameter(0), E0, 0.);
    f1->FixParameter(1, E0);
    g->Fit(f1, "RN");

    // Modello 3: broken power law  k*(E/Eb)^Delta per E<=Eb, k per E>Eb  (continua in Eb)
    // par: [0]=k  [1]=Eb  [2]=Delta
    bool doBreak = (n >= 5);
    TF1* f2 = new TF1("f2", "[0]*(x<=[1] ? TMath::Power(x/[1],[2]) : 1.)", xlo, xhi);
    if (doBreak) {
        // scansione in Eb (log) con Eb fisso: si fittano k e Delta, si tiene il chi2 minimo
        double bestChi2 = 1e300, bestEb = E0;
        const int NSCAN = 40;
        for (int k = 0; k <= NSCAN; k++) {
            double Eb = x[1] * TMath::Power(x[n-2] / x[1], double(k) / NSCAN);
            f2->SetParameters(f0->GetParameter(0), Eb, 0.);
            f2->FixParameter(1, Eb);
            g->Fit(f2, "RNQ");
            if (f2->GetChisquare() < bestChi2) { bestChi2 = f2->GetChisquare(); bestEb = Eb; }
        }
        // fit finale con Eb libero, partendo dal minimo della scansione
        f2->SetParameters(f0->GetParameter(0), bestEb, 0.);
        f2->ReleaseParameter(1);
        f2->SetParLimits(1, x[1], x[n-2]);
        g->Fit(f2, "RN");
        std::cout << "[broken PL]  k = " << f2->GetParameter(0) << " +- " << f2->GetParError(0)
                  << "  Eb = " << f2->GetParameter(1) << " +- " << f2->GetParError(1)
                  << "  Delta = " << f2->GetParameter(2) << " +- " << f2->GetParError(2)
                  << "  chi2/ndf = " << f2->GetChisquare() << "/" << f2->GetNDF() << "\n";
    }

    double dchi2 = f0->GetChisquare() - f1->GetChisquare();
    double pval  = TMath::Prob(dchi2, 1);
    std::cout << "E0 = " << E0 << "\n"
              << "[const]   R0 = " << f0->GetParameter(0) << " +- " << f0->GetParError(0)
              << "  chi2/ndf = " << f0->GetChisquare() << "/" << f0->GetNDF() << "\n"
              << "[power law]  R0 = " << f1->GetParameter(0) << " +- " << f1->GetParError(0)
              << "  DGAMMA = " << f1->GetParameter(2) << " +- " << f1->GetParError(2)
              << "  chi2/ndf = " << f1->GetChisquare() << "/" << f1->GetNDF() << "\n"
              << "Delta chi2 = " << dchi2 << "  p-value (1 dof) = " << pval << "\n";

    if (!draw && !save) return;

    // ---------------- Grafica ----------------
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // Grafico per il disegno (con barre di errore in x = semi-larghezza bin)
    TGraphErrors* gd = new TGraphErrors(n);
    TGraphErrors* p0 = new TGraphErrors(n);   // pull modello costante
    TGraphErrors* p1 = new TGraphErrors(n);   // pull modello power law
    TGraphErrors* p2 = new TGraphErrors(n);   // pull broken power law
    for (int i = 0; i < n; i++) {
        gd->SetPoint(i, x[i], r[i]); gd->SetPointError(i, ex[i], er[i]);
        double pu0 = (r[i] - f0->Eval(x[i])) / er[i];
        double pu1 = (r[i] - f1->Eval(x[i])) / er[i];
        p0->SetPoint(i, x[i] * 0.98, pu0); p0->SetPointError(i, 0, 1);
        p1->SetPoint(i, x[i] * 1.02, pu1); p1->SetPointError(i, 0, 1);
        if (doBreak) {
            p2->SetPoint(i, x[i] * 1.04, (r[i] - f2->Eval(x[i])) / er[i]); p2->SetPointError(i, 0, 1);
        }
    }
    gd->SetMarkerStyle(20); gd->SetMarkerSize(1.0); gd->SetLineWidth(2);
    f0->SetLineColor(kBlue + 1); f0->SetLineStyle(2); f0->SetLineWidth(2);
    f1->SetLineColor(kRed + 1);  f1->SetLineStyle(1); f1->SetLineWidth(2);
    p0->SetMarkerStyle(24); p0->SetMarkerColor(kBlue + 1); p0->SetLineColor(kBlue + 1);
    p1->SetMarkerStyle(20); p1->SetMarkerColor(kRed + 1);  p1->SetLineColor(kRed + 1);
    p2->SetMarkerStyle(21); p2->SetMarkerColor(kGreen + 2); p2->SetLineColor(kGreen + 2);
    f2->SetLineColor(kGreen + 2); f2->SetLineStyle(1); f2->SetLineWidth(3);

    TCanvas* c = new TCanvas("c_ratio", "Ratio Ar/S", 800, 1000);
    c->Divide(1, 3, 0, 0);

    // --- 1: flux_pow E^2.7  ---
    c->cd(1);
    gPad->SetLogx(); gPad->SetLogy(); gPad->SetGrid();
    gPad->SetLeftMargin(0.14); gPad->SetTopMargin(0.05);
    if (hArP && hSP) {
        TH1D* cAr = (TH1D*)hArP->Clone("cAr"); cAr->SetDirectory(0);
        TH1D* cS  = (TH1D*)hSP->Clone("cS");   cS->SetDirectory(0);
        double ymax = std::max(cAr->GetMaximum(), cS->GetMaximum());
        double ymin = 1e300;
        for (int i = 1; i <= cAr->GetNbinsX(); i++) {
            if (cAr->GetBinContent(i) > 0) ymin = std::min(ymin, cAr->GetBinContent(i));
            if (cS->GetBinContent(i)  > 0) ymin = std::min(ymin, cS->GetBinContent(i));
        }
        if (ymin > ymax) ymin = ymax * 1e-3;
        cS->SetMarkerStyle(21);  cS->SetMarkerColor(kGreen + 2);   cS->SetLineColor(kGreen + 2);
        cAr->SetMarkerStyle(22); cAr->SetMarkerColor(kMagenta + 1); cAr->SetLineColor(kMagenta + 1);
        cS->GetXaxis()->SetRangeUser(xlo, xhi);
        cS->SetMinimum(ymin * 0.5); cS->SetMaximum(ymax * 2.);
        cS->GetYaxis()->SetTitle("Flux #times E^{2.7}");
        cS->GetYaxis()->SetTitleSize(0.06); cS->GetYaxis()->SetTitleOffset(1.1);
        cS->GetYaxis()->SetLabelSize(0.05); cS->GetXaxis()->SetLabelSize(0);
        cS->Draw("E1");
        cAr->Draw("E1 SAME");
        TLegend* l1 = new TLegend(0.60, 0.72, 0.92, 0.92);
        l1->SetBorderSize(0); l1->SetFillStyle(0);
        l1->AddEntry(cS, "S", "lp"); l1->AddEntry(cAr, "Ar", "lp");
        l1->Draw();
    } else {
        gPad->SetLogy(0);
        TLegend* l1 = new TLegend(0.2, 0.4, 0.9, 0.6);
        l1->SetBorderSize(0); l1->AddEntry((TObject*)0, "flux_pow not found in files", "");
        l1->Draw();
    }

    // ---  2: ratio Ar/S + fit ---
    c->cd(2);
    gPad->SetLogx(); gPad->SetGrid();
    gPad->SetLeftMargin(0.14);
    double rmin = 1e300, rmax = -1e300;
    for (int i = 0; i < n; i++) {
        rmin = std::min(rmin, r[i] - er[i]); rmax = std::max(rmax, r[i] + er[i]);
    }
    double pad = 0.25 * (rmax - rmin);
    double ylo = std::max(0., rmin - pad), yhi = rmax + pad;
    TH1F* fr = gPad->DrawFrame(xlo, ylo, xhi, yhi);
    fr->GetYaxis()->SetTitle("F_{Ar} / F_{S}");
    fr->GetYaxis()->SetTitleSize(0.06); fr->GetYaxis()->SetTitleOffset(1.1);
    fr->GetYaxis()->SetLabelSize(0.05); fr->GetXaxis()->SetLabelSize(0);
    gd->Draw("P SAME");
    f0->Draw("SAME"); f1->Draw("SAME");
    if (doBreak) f2->Draw("SAME");
    TLine* lE0 = new TLine(E0, ylo, E0, yhi);
    lE0->SetLineColor(kGray + 1); lE0->SetLineStyle(3); lE0->Draw();

    TLegend* l2 = new TLegend(0.17, 0.52, 0.93, 0.93);
    l2->SetBorderSize(0); l2->SetFillStyle(0); l2->SetTextSize(0.042);
    l2->AddEntry(gd, "Ar/S", "pe");
    l2->AddEntry(f0, Form("costante: R_{0}=%.3g #pm %.2g,  #chi^{2}/ndf=%.2f/%d",
                          f0->GetParameter(0), f0->GetParError(0), f0->GetChisquare(), f0->GetNDF()), "l");
    l2->AddEntry(f1, Form("power law: R_{0}=%.3g #pm %.2g, #Delta#gamma=%.3f #pm %.3f, #chi^{2}/ndf=%.2f/%d",
                          f1->GetParameter(0), f1->GetParError(0), f1->GetParameter(2),
                          f1->GetParError(2), f1->GetChisquare(), f1->GetNDF()), "l");
    if (doBreak)
        l2->AddEntry(f2, Form("broken PL: k=%.3g #pm %.2g, E_{b}=%.3g GeV, #Delta=%.3f #pm %.3f, #chi^{2}/ndf=%.2f/%d",
                              f2->GetParameter(0), f2->GetParError(0), f2->GetParameter(1),
                              f2->GetParameter(2), f2->GetParError(2), f2->GetChisquare(), f2->GetNDF()), "l");
    l2->AddEntry((TObject*)0, Form("E_{0}=%.3g GeV,  #Delta#chi^{2}=%.2f (p=%.2g)", E0, dchi2, pval), "");
    l2->Draw();

    // --- 3: pull ---
    c->cd(3);
    gPad->SetLogx(); gPad->SetGrid();
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.22);
    TH1F* fp = gPad->DrawFrame(xlo, -3.5, xhi, 3.5);
    fp->GetXaxis()->SetTitle("Energy [GeV]");
    fp->GetXaxis()->SetTitleSize(0.07); fp->GetXaxis()->SetTitleOffset(1.2);
    fp->GetXaxis()->SetLabelSize(0.06);
    fp->GetYaxis()->SetTitle("(data - fit) / #sigma");
    fp->GetYaxis()->SetTitleSize(0.06); fp->GetYaxis()->SetTitleOffset(1.1);
    fp->GetYaxis()->SetLabelSize(0.05);
    TLine* l0 = new TLine(xlo, 0, xhi, 0); l0->SetLineStyle(2); l0->Draw();
    p0->Draw("P SAME"); p1->Draw("P SAME");
    if (doBreak) p2->Draw("P SAME");
    TLegend* l3 = new TLegend(0.60, 0.70, 0.92, 0.95);
    l3->SetBorderSize(0); l3->SetFillStyle(0);
    l3->AddEntry(p0, "constant", "p"); l3->AddEntry(p1, "power law", "p");
    if (doBreak) l3->AddEntry(p2, "broken power law", "p");
    l3->Draw();

    c->Update();
    /*
    if (save) {
        c->SaveAs(Form("%s.pdf", outname));
        c->SaveAs(Form("%s.png", outname));
        std::cout << "Grafici salvati in " << outname << ".pdf / .png\n";
    }
    */
}