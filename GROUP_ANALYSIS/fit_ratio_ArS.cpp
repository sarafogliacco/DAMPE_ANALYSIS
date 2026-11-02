// Ar/S flux ratio: constant vs power-law vs broken-power-law fits.
// Usage:  root -l 'fit_ratio_ArS.C+("ar.root","s.root")'
//         root -l -b -q 'fit_ratio_ArS.C+("ar.root","s.root",1e2,1e5,false,"ratio_ArS",1000)'

#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TString.h>
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>

namespace {

struct RatioData {
    std::vector<double> x, ex, r, er;   // energy, half bin width, ratio, ratio error
    size_t n() const { return x.size(); }
};

// Lafferty-Wyatt point: energy at which the true differential flux (E^-alpha)
// equals the bin-averaged flux. (This is NOT the median of the bin.)
double lafferty_wyatt(double lo, double hi, double alpha = 2.7) {
    const double a1 = 1. - alpha;
    return TMath::Power((TMath::Power(hi, a1) - TMath::Power(lo, a1)) / (a1 * (hi - lo)),
                        -1. / alpha);
}

// Read a histogram and detach it from the file so the file can be closed safely.
TH1D* load_hist(TFile* f, const char* name) {
    TH1D* h = f ? dynamic_cast<TH1D*>(f->Get(name)) : nullptr;
    if (h) h->SetDirectory(nullptr);
    return h;
}

bool same_binning(const TH1D* a, const TH1D* b) {
    if (a->GetNbinsX() != b->GetNbinsX()) {
        std::cerr << "Binning mismatch: " << a->GetNbinsX() << " vs " << b->GetNbinsX() << " bins\n";
        return false;
    }
    for (int i = 1; i <= a->GetNbinsX() + 1; i++) {
        double ea = a->GetXaxis()->GetBinLowEdge(i), eb = b->GetXaxis()->GetBinLowEdge(i);
        if (TMath::Abs(ea - eb) > 1e-6 * TMath::Max(1., TMath::Abs(ea))) {
            std::cerr << "Bin edges differ at bin " << i << ": " << ea << " vs " << eb << "\n";
            return false;
        }
    }
    return true;
}

// Ratio hA/hB with uncorrelated error propagation. Bins with non-positive flux
// or outside [emin, emax] are skipped (and reported).
RatioData build_ratio(const TH1D* hA, const TH1D* hB, double emin, double emax) {
    RatioData d;
    for (int i = 1; i <= hA->GetNbinsX(); i++) {
        double fa = hA->GetBinContent(i), ea = hA->GetBinError(i);
        double fb = hB->GetBinContent(i), eb = hB->GetBinError(i);
        double lo = hA->GetXaxis()->GetBinLowEdge(i), up = hA->GetXaxis()->GetBinUpEdge(i);
        double em = lafferty_wyatt(lo, up);
        if (em < emin || em > emax) continue;
        if (fa <= 0 || fb <= 0) {
            std::cerr << "Skipping bin " << i << " [" << lo << "," << up << "] (non-positive flux)\n";
            continue;
        }
        double ri  = fa / fb;
        double eri = ri * TMath::Hypot(ea / fa, eb / fb);
        d.x.push_back(em); d.r.push_back(ri); d.er.push_back(eri);
        d.ex.push_back(0.5 * (up - lo));
    }
    return d;
}

TGraphErrors* make_graph(const RatioData& d, bool xerr) {
    auto* g = new TGraphErrors((int)d.n());
    for (size_t i = 0; i < d.n(); i++) {
        g->SetPoint(i, d.x[i], d.r[i]);
        g->SetPointError(i, xerr ? d.ex[i] : 0., d.er[i]);
    }
    return g;
}

// ---- models -----------------------------------------------------------
// power law:  R0 * (E/E0)^(-DGAMMA)        par: [0]=R0 [1]=E0 [2]=DGAMMA
// broken PL:  k*(E/Eb)^Delta for E<=Eb, k above   par: [0]=k [1]=Eb [2]=Delta
TF1* make_pl(const char* name, double xlo, double xhi, double E0, bool fixSlope) {
    auto* f = new TF1(name, "[0]*TMath::Power(x/[1], -[2])", xlo, xhi);
    f->SetParameters(1., E0, 0.);
    f->FixParameter(1, E0);
    if (fixSlope) f->FixParameter(2, 0.);
    return f;
}

TF1* make_broken(const char* name, double xlo, double xhi) {
    return new TF1(name, "[0]*(x<=[1] ? TMath::Power(x/[1],[2]) : 1.)", xlo, xhi);
}

int do_fit(TGraphErrors* g, TF1* f, const char* opt = "RNSQ") {
    TFitResultPtr res = g->Fit(f, opt);
    int st = (int)res;
    if (st != 0) std::cerr << "Warning: fit of " << f->GetName() << " returned status " << st << "\n";
    return st;
}

// Broken PL: scan Eb on a log grid (fit k, Delta with Eb fixed), then release Eb.
// Needs a freshly constructed f (no stale limits).
void fit_broken(TGraphErrors* g, TF1* f, const std::vector<double>& x, double k0) {
    const int n = (int)x.size(), NSCAN = 40;
    const double lo = x[1], hi = x[n - 2];
    double bestChi2 = 1e300, bestEb = lo;
    for (int k = 0; k <= NSCAN; k++) {
        double Eb = lo * TMath::Power(hi / lo, double(k) / NSCAN);
        f->SetParameters(k0, Eb, 0.);
        f->FixParameter(1, Eb);
        g->Fit(f, "RNQ");
        if (f->GetChisquare() < bestChi2) { bestChi2 = f->GetChisquare(); bestEb = Eb; }
    }
    f->SetParameters(k0, bestEb, 0.);
    f->ReleaseParameter(1);
    f->SetParLimits(1, lo, hi);
    g->Fit(f, "RNQ");
}

// Toy MC: how often does a *constant* truth give a chi2 improvement from the
// broken PL at least as large as the data's? (accounts for the Eb scan / look-elsewhere)
double broken_toy_pvalue(const RatioData& d, const TF1* f0, double dchi2_data, int nToys,
                         double xlo, double xhi, double E0) {
    TRandom3 rng(12345);
    int nExceed = 0;
    for (int t = 0; t < nToys; t++) {
        RatioData td = d;
        for (size_t i = 0; i < d.n(); i++) td.r[i] = rng.Gaus(f0->Eval(d.x[i]), d.er[i]);
        std::unique_ptr<TGraphErrors> g(make_graph(td, false));
        std::unique_ptr<TF1> fc(make_pl("ft_c", xlo, xhi, E0, true));
        std::unique_ptr<TF1> fb(make_broken("ft_b", xlo, xhi));
        g->Fit(fc.get(), "RNQ");
        fit_broken(g.get(), fb.get(), td.x, fc->GetParameter(0));
        if (fc->GetChisquare() - fb->GetChisquare() >= dchi2_data) nExceed++;
    }
    return double(nExceed) / nToys;
}

TGraphErrors* make_pull(const RatioData& d, const TF1* f, double xshift,
                        int mstyle, int color) {
    auto* p = new TGraphErrors((int)d.n());
    for (size_t i = 0; i < d.n(); i++) {
        p->SetPoint(i, d.x[i] * xshift, (d.r[i] - f->Eval(d.x[i])) / d.er[i]);
        p->SetPointError(i, 0, 1);
    }
    p->SetMarkerStyle(mstyle); p->SetMarkerColor(color); p->SetLineColor(color);
    return p;
}

}  // namespace

void fit_ratio_ArS(const char* fAr = "/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/ARGON_ANALYSIS/unfold_results/unfold_results_argon_6bin_low1_up2_v2.root",                                         
                   const char* fS  = "/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/SULFUR_ANALYSIS/unfold_results/unfold_results_flux_spectrum_sulfur_120months_6bin_smooth_2e5up_1e5low_bkg_v2.root",
                   double Efit_min = 2*1e2, double Efit_max = 1e5,
                   bool draw = true, const char* outname = "ratio_ArS",
                   int nToys = 0)
{
    const bool save = (outname && outname[0] != '\0');
    const bool wasBatch = gROOT->IsBatch();
    if (!draw) gROOT->SetBatch(kTRUE);

    // ---------------- input ----------------
    std::unique_ptr<TFile> a(TFile::Open(fAr)), s(TFile::Open(fS));
    if (!a || a->IsZombie() || !s || s->IsZombie()) { std::cerr << "Error opening file\n"; return; }
    std::unique_ptr<TH1D> hAr(load_hist(a.get(), "flux")),      hS(load_hist(s.get(), "flux"));
    std::unique_ptr<TH1D> hArP(load_hist(a.get(), "flux_pow")), hSP(load_hist(s.get(), "flux_pow"));
    a->Close(); s->Close();
    if (!hAr || !hS) { std::cerr << "Histogram 'flux' not found\n"; return; }
    if (!same_binning(hAr.get(), hS.get())) return;

    RatioData d = build_ratio(hAr.get(), hS.get(), Efit_min, Efit_max);
    const int n = (int)d.n();
    if (n < 2) { std::cerr << "Too few points (" << n << ")\n"; return; }

    // Pivot energy = log-centre of the range: decorrelates R0 and DGAMMA
    const double E0  = TMath::Sqrt(d.x.front() * d.x.back());
    const double xlo = d.x.front() * 0.9, xhi = d.x.back() * 1.1;

    std::unique_ptr<TGraphErrors> g(make_graph(d, false));   // fit graph
    TGraphErrors* gd = make_graph(d, true);                  // drawing graph (kept alive for the canvas)

    // ---------------- fits ----------------
    TF1* f0 = make_pl("f_const", xlo, xhi, E0, true);
    do_fit(g.get(), f0);

    TF1* f1 = make_pl("f_pl", xlo, xhi, E0, false);
    f1->SetParameter(0, f0->GetParameter(0));
    do_fit(g.get(), f1);

    const bool doBreak = (n >= 5);
    TF1* f2 = make_broken("f_break", xlo, xhi);
    double dchi2_b = 0, pval_b = -1;
    if (doBreak) {
        fit_broken(g.get(), f2, d.x, f0->GetParameter(0));
        dchi2_b = f0->GetChisquare() - f2->GetChisquare();
        if (nToys > 0) pval_b = broken_toy_pvalue(d, f0, dchi2_b, nToys, xlo, xhi, E0);
    }

    const double dchi2 = f0->GetChisquare() - f1->GetChisquare();
    const double pval  = TMath::Prob(dchi2, 1);
    const double nsig  = TMath::Sqrt(TMath::Max(dchi2, 0.));    // exact for 1 dof; no underflow at tiny p

    auto report = [&](const char* tag, TF1* f, int npar) {
        double chi2 = f->GetChisquare();
        std::cout << "[" << tag << "]  chi2/ndf = " << chi2 << "/" << f->GetNDF()
                  << "  (GOF p = " << TMath::Prob(chi2, f->GetNDF()) << ")"
                  << "  AIC = " << chi2 + 2. * npar
                  << "  BIC = " << chi2 + npar * TMath::Log(n) << "\n";
    };
    std::cout << "N points = " << n << ",  E0 = " << E0 << " GeV\n";
    report("const", f0, 1);
    std::cout << "    R0 = " << f0->GetParameter(0) << " +- " << f0->GetParError(0) << "\n";
    report("power law", f1, 2);
    std::cout << "    R0 = " << f1->GetParameter(0) << " +- " << f1->GetParError(0)
              << "  DGAMMA = " << f1->GetParameter(2) << " +- " << f1->GetParError(2) << "\n";
    std::cout << "    Delta chi2 (const vs PL) = " << dchi2 << ",  p(1 dof) = " << pval
              << " (" << nsig << " sigma)\n";
    if (doBreak) {
        report("broken PL", f2, 3);
        std::cout << "    k = " << f2->GetParameter(0) << " +- " << f2->GetParError(0)
                  << "  Eb = " << f2->GetParameter(1) << " +- " << f2->GetParError(1)
                  << "  Delta = " << f2->GetParameter(2) << " +- " << f2->GetParError(2) << "\n"
                  << "    Delta chi2 (const vs break) = " << dchi2_b << "\n";
        if (pval_b >= 0)
            std::cout << "    toy-MC p-value (look-elsewhere included) = " << pval_b
                      << "  [" << nToys << " toys]\n";
        else
            std::cout << "    NOTE: chi2 is not chi2-distributed for a break fit (Eb is not defined under\n"
                         "          the null); set nToys > 0 for a proper p-value.\n";
    }

    if (!draw && !save) { gROOT->SetBatch(wasBatch); return; }

    // ---------------- graphics ----------------
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TGraphErrors* p0 = make_pull(d, f0, 0.98, 24, kBlue + 1);
    TGraphErrors* p1 = make_pull(d, f1, 1.02, 20, kRed + 1);
    TGraphErrors* p2 = doBreak ? make_pull(d, f2, 1.04, 21, kGreen + 2) : nullptr;

    gd->SetMarkerStyle(20); gd->SetMarkerSize(1.0); gd->SetLineWidth(2);
    f0->SetLineColor(kBlue + 1);  f0->SetLineStyle(2); f0->SetLineWidth(2);
    f1->SetLineColor(kRed + 1);   f1->SetLineStyle(1); f1->SetLineWidth(2);
    f2->SetLineColor(kGreen + 2); f2->SetLineStyle(1); f2->SetLineWidth(3);

    // Explicit pads: the bottom pad is taller so that its *plot area* matches the others
    // (its bottom margin holds the x-axis labels/title).
    TCanvas* c = new TCanvas("c_ratio", "Ratio Ar/S", 800, 1000);
    const double RM = 0.04;                       // right margin: keeps the right frame visible
    TPad* pd1 = new TPad("pd1", "", 0, 0.66, 1, 1.00);
    TPad* pd2 = new TPad("pd2", "", 0, 0.36, 1, 0.66);
    TPad* pd3 = new TPad("pd3", "", 0, 0.00, 1, 0.36);
    pd1->SetTopMargin(0.04);    pd1->SetBottomMargin(0);
    pd2->SetTopMargin(0);       pd2->SetBottomMargin(0);
    pd3->SetTopMargin(0);       pd3->SetBottomMargin(0.17);
    for (TPad* p : {pd1, pd2, pd3}) {
        p->SetLeftMargin(0.14); p->SetRightMargin(RM);
        p->SetLogx(); p->SetGrid(); p->Draw();
    }

    // --- pad 1: flux x E^2.7 ---
    pd1->cd();
    gPad->SetLogy();
    if (hArP && hSP) {
        TH1D* cAr = (TH1D*)hArP->Clone("cAr"); cAr->SetDirectory(nullptr);
        TH1D* cS  = (TH1D*)hSP->Clone("cS");   cS->SetDirectory(nullptr);
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
        TLegend* l1 = new TLegend(0.2, 0.4, 0.9, 0.6);
        l1->SetBorderSize(0);
        l1->AddEntry((TObject*)nullptr, "flux_pow not found in files", "");
        l1->Draw();
    }

    // --- pad 2: ratio + fits ---
    pd2->cd();
    double rmin = 1e300, rmax = -1e300;
    for (int i = 0; i < n; i++) {
        rmin = std::min(rmin, d.r[i] - d.er[i]); rmax = std::max(rmax, d.r[i] + d.er[i]);
    }
    const double pad = 0.25 * (rmax - rmin);
    const double ylo = std::max(0., rmin - pad), yhi = rmax + pad;
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
    l2->AddEntry(f0, Form("constant: R_{0}=%.3g #pm %.2g,  #chi^{2}/ndf=%.2f/%d",
                          f0->GetParameter(0), f0->GetParError(0), f0->GetChisquare(), f0->GetNDF()), "l");
    l2->AddEntry(f1, Form("power law: R_{0}=%.3g #pm %.2g, #Delta#gamma=%.3f #pm %.3f, #chi^{2}/ndf=%.2f/%d",
                          f1->GetParameter(0), f1->GetParError(0), f1->GetParameter(2),
                          f1->GetParError(2), f1->GetChisquare(), f1->GetNDF()), "l");
    if (doBreak)
        l2->AddEntry(f2, Form("broken PL: k=%.3g #pm %.2g, E_{b}=%.3g GeV, #Delta=%.3f #pm %.3f, #chi^{2}/ndf=%.2f/%d",
                              f2->GetParameter(0), f2->GetParError(0), f2->GetParameter(1),
                              f2->GetParameter(2), f2->GetParError(2), f2->GetChisquare(), f2->GetNDF()), "l");
    l2->AddEntry((TObject*)nullptr,
                 Form("E_{0}=%.3g GeV,  #Delta#chi^{2}=%.2f (p=%.2g, %.1f#sigma)", E0, dchi2, pval, nsig), "");
    l2->Draw();

    // --- pad 3: pulls ---
    pd3->cd();
    TH1F* fp = gPad->DrawFrame(xlo, -10., xhi, 10.);
    fp->GetXaxis()->SetTitle("Energy [GeV]");
    fp->GetXaxis()->SetTitleSize(0.07); fp->GetXaxis()->SetTitleOffset(1.2);
    fp->GetXaxis()->SetLabelSize(0.06);
    fp->GetYaxis()->SetTitle("(data - fit) / #sigma");
    fp->GetYaxis()->SetTitleSize(0.06); fp->GetYaxis()->SetTitleOffset(1.1);
    fp->GetYaxis()->SetLabelSize(0.05);
    TLine* l0 = new TLine(xlo, 0, xhi, 0); l0->SetLineStyle(2); l0->Draw();
    p0->Draw("P SAME"); p1->Draw("P SAME");
    if (p2) p2->Draw("P SAME");
    TLegend* l3 = new TLegend(0.16, 0.88, 0.96, 0.98);
    l3->SetNColumns(3); l3->SetTextSize(0.04);
    l3->SetBorderSize(0); l3->SetFillStyle(0);
    l3->AddEntry(p0, "constant", "p"); l3->AddEntry(p1, "power law", "p");
    if (p2) l3->AddEntry(p2, "broken power law", "p");
    l3->Draw();

    c->Update();
/*
    if (save) {
        c->SaveAs(Form("%s.pdf", outname));
        c->SaveAs(Form("%s.png", outname));
        std::cout << "Plots saved to " << outname << ".pdf / .png\n";
    }
*/
    gROOT->SetBatch(wasBatch);
}