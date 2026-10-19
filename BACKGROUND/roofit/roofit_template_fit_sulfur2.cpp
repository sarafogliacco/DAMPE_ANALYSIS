#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLine.h"
#include "TColor.h"
#include "TSystem.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooAbsReal.h"
#include "RooFormulaVar.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooArgList.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

using namespace RooFit;
using namespace std;

// ---------------------------------------------------------------------
// One template (one nuclear species)
// ---------------------------------------------------------------------
struct Elem {
	string name;      // short name used in output ("Si", "S", ...)
	string label;     // legend label
	string path;      // ROOT file with the MC histograms
	string color;     // hex colour
	bool   optional;  // if true and file is missing -> skip silently
	// filled at runtime
	TH1D        *h   = nullptr;
	RooDataHist *dh  = nullptr;
	RooHistPdf  *pdf = nullptr;
	RooRealVar  *n   = nullptr;
	RooAbsReal  *I   = nullptr;
};

// ---------------------------------------------------------------------
// Usage:
//   .x roofit_template_fit_sulfur.cpp(4, 2, false)
//   range_low_req / range_high_req are snapped to bin edges. Repeat the fit
//   with e.g. 13.5, 14.0, 14.5 to estimate the systematic uncertainty.
// ---------------------------------------------------------------------
void roofit_template_fit_sulfur(int iBin = 4, int iRebin = 2, bool writefile = false,
                                double range_low_req = 14.0, double range_high_req = 21.0){

	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPadLeftMargin(0.09);
	gStyle->SetPadRightMargin(0.19);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadBottomMargin(0.09);

	// ------------------------------------------------------------------
	// Energy bin
	// ------------------------------------------------------------------
	double E = 0.0;
	TString title;

	if(iBin == 1) { E = 121.1527961; title = "100 GeV < E_{BGO} < 178 GeV"  ; }
	if(iBin == 2) { E = 177.8277918; title = "178 GeV < E_{BGO} < 316 GeV"  ; }
	if(iBin == 3) { E = 261.015534 ; title = "316 GeV < E_{BGO} < 562 GeV"  ; }
	if(iBin == 4) { E = 383.1188748; title = "562 GeV < E_{BGO} < 1.0 TeV"  ; }
	if(iBin == 5) { E = 562.3413673; title = "1.0 TeV < E_{BGO} < 1.8 TeV"  ; }
	if(iBin == 6) { E = 825.4041434; title = "1.8 TeV < E_{BGO} < 3.2 TeV"  ; }
	if(iBin == 7) { E = 1467.79767 ; title = "3.2 TeV < E_{BGO} < 5.6 TeV"  ; }
	if(iBin == 8) { E = 3162.274615; title = "5.6 TeV < E_{BGO} < 10.0 TeV" ; }
	if(iBin == 9) { E = 6812.921547; title = "10.0 TeV < E_{BGO} < 31.6 TeV"; }
	if(iBin ==10) { E = 31622.7766 ; title = "31.6 TeV < E_{BGO} < 100.0 TeV";}

	if (E <= 0.0) { cerr << "Invalid iBin = " << iBin << endl; return; }

	// Sulfur purity window (09/2026)
	double L = log10(E);
	double S_high = 15.915774  + 0.6569678*L  - 0.167878*pow(L,2) + 0.0179829*pow(L,3);
	double S_low  = 15.1162445 + 0.53318365*L - 0.167878*pow(L,2) + 0.0179829*pow(L,3);

	cout << "High limit: " << S_high << endl;
	cout << "Low limit: "  << S_low  << endl;

	// Histogram name: h01..h09, h10
	string hname = (iBin < 10) ? Form("h0%d", iBin) : Form("h%d", iBin);

	// ------------------------------------------------------------------
	// Templates (order = order in the legend and in the output file)
	// To add a species (e.g. Al, P) just add a line here.
	// ------------------------------------------------------------------
	const string base  = "/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/";
	vector<Elem> el = {
		{"Si", "Silicon",   base+"BACKGROUND/out_root/MC_FTFP_SILICON_SATcorr_STKcut1200_PSDcutXY_CORRECTED_mcfit.root",                       "#6a4c93", false},
		// Phosphorus (Z=15): set the correct path. If the file does not exist it is skipped.
		{"P",  "Phosphorus",base+"secondaries/out_root/MC_FTFP_PHOSPHORUS_SATcorr_STKcut1200_PSDcutXY_CORRECTED.root",                        "#e6d800", true },
		{"S",  "Sulfur",    base+"SULFUR_ANALYSIS/out_root/MC_FTFP_SULFUR_SATcorr_STKcut1200_PSDcutXY_CORRECTED_1PeV_mcfit.root",             "#006627", false},
		{"Cl", "Chlorine",  base+"secondaries/out_root/MC_FTFP_CLORO_SATcorr_STKcut1200_PSDcutXY_CORRECTED.root",                              "#0bb4ff", false},
		{"Ar", "Argon",     base+"ARGON_ANALYSIS/out_root/MC_FTFP_ARGON_SATcorr_STKcut1200_PSDcutXY_CORRECTED_mcfit.root",                    "#50e991", false},
		{"K",  "Potassium", base+"secondaries/out_root/MC_FTFP_POTASSIUM_SATcorr_STKcut1200_PSDcutXY_CORRECTED.root",                          "#007585", false},
		{"Ca", "Calcium",   base+"CALCIUM_ANALYSIS/out_root/MC_FTFP_CALCIUM_SATcorr_STKcut1200_PSDcutXY_CORRECTED_100GeV-100TeV_mcfit.root", "#ff0000", false},
	};

	// ------------------------------------------------------------------
	// Flight data
	// ------------------------------------------------------------------
	TFile *f_data = TFile::Open((base+"BACKGROUND/out_root/Skim_120months_SATcorr_STKcut1200_PSDcutXY_10bin_v2.root").c_str());
	if (!f_data || f_data->IsZombie()) { cerr << "Cannot open data file" << endl; return; }
	TH1D *h_data = (TH1D*)f_data->Get(hname.c_str());
	if (!h_data) { cerr << "Data histogram not found: " << hname << endl; return; }
	h_data = (TH1D*)h_data->Clone("h_data_clone");
	h_data->SetDirectory(0);
	h_data->Rebin(iRebin);

	// ------------------------------------------------------------------
	// Fit range snapped to bin edges of the (rebinned) data histogram
	// ------------------------------------------------------------------
	TAxis *ax = h_data->GetXaxis();
	int bLow  = ax->FindBin(range_low_req + 1e-6);
	int bHigh = ax->FindBin(range_high_req - 1e-6);
	double range_low  = ax->GetBinLowEdge(bLow);
	double range_high = ax->GetBinUpEdge(bHigh);
	cout << "Fit range (snapped to bin edges): " << range_low << " - " << range_high << endl;

	// ------------------------------------------------------------------
	// Load MC histograms
	// ------------------------------------------------------------------
	printf("DATA nbins=%d [%.4f, %.4f] width=%.5f\n", h_data->GetNbinsX(),
	       ax->GetXmin(), ax->GetXmax(), h_data->GetBinWidth(1));

	// Resample a histogram onto the binning of a reference histogram.
	// Content is redistributed by bin overlap (uniform inside each source bin),
	// so the number of events is conserved.
	auto resample = [](TH1D *src, TH1D *ref, const char *name) -> TH1D* {
		TH1D *out = (TH1D*)ref->Clone(name);
		out->Reset();
		out->SetDirectory(0);
		TAxis *as = src->GetXaxis();
		for (int i = 1; i <= ref->GetNbinsX(); ++i) {
			double lo = ref->GetXaxis()->GetBinLowEdge(i);
			double hi = ref->GetXaxis()->GetBinUpEdge(i);
			double sum = 0;
			int j0 = max(1, as->FindBin(lo));
			int j1 = min(src->GetNbinsX(), as->FindBin(hi));
			for (int j = j0; j <= j1; ++j) {
				double slo = as->GetBinLowEdge(j), shi = as->GetBinUpEdge(j);
				double ov  = min(hi, shi) - max(lo, slo);
				if (ov > 0) sum += src->GetBinContent(j) * ov / (shi - slo);
			}
			out->SetBinContent(i, sum);
			out->SetBinError(i, sqrt(fabs(sum)));
		}
		return out;
	};

	vector<Elem> used;
	for (auto &e : el) {
		if (e.optional && gSystem->AccessPathName(e.path.c_str())) {
			cout << "[INFO] " << e.name << " template not found, skipped: " << e.path << endl;
			continue;
		}
		TFile *f = TFile::Open(e.path.c_str());
		if (!f || f->IsZombie()) { cerr << "Cannot open MC file for " << e.name << ": " << e.path << endl; return; }
		TH1D *h = (TH1D*)f->Get(hname.c_str());
		if (!h) { cerr << "Histogram " << hname << " not found for " << e.name << endl; f->ls(); return; }
		printf("%-3s original: nbins=%d [%.4f, %.4f] width=%.5f  -> resampled onto data binning\n",
		       e.name.c_str(), h->GetNbinsX(), h->GetXaxis()->GetXmin(),
		       h->GetXaxis()->GetXmax(), h->GetBinWidth(1));
		// no Rebin on the MC: resample the original histogram directly onto the data axis
		e.h = resample(h, h_data, Form("h_mc_%s", e.name.c_str()));
		printf("    integral in fit range = %.1f\n", e.h->Integral(bLow, bHigh));
		if (e.h->Integral(bLow, bHigh) <= 0) {
			cerr << "Template " << e.name << " is empty in the fit range." << endl; return;
		}
		// remove negative bins
		for (int i = 1; i <= e.h->GetNbinsX(); i++)
			if (e.h->GetBinContent(i) < 0) e.h->SetBinContent(i, 0);

		used.push_back(e);
	}
	el = used;

	// ------------------------------------------------------------------
	// RooFit objects
	// ------------------------------------------------------------------
	RooRealVar x("x", "PSD Charge ", range_low, range_high);
	RooDataHist data("data", "data", x, h_data);

	double N_data = h_data->Integral(bLow, bHigh);
	double seed_tot = 0;
	for (auto &e : el) seed_tot += e.h->Integral(bLow, bHigh);

	RooArgList pdfList, coefList;
	for (auto &e : el) {
		e.dh  = new RooDataHist(Form("mc%s", e.name.c_str()), Form("mc %s", e.name.c_str()), x, e.h);
		e.pdf = new RooHistPdf(Form("modelmc%s", e.name.c_str()), Form("modelmc%s", e.name.c_str()), x, *e.dh);
		double seed = e.h->Integral(bLow, bHigh) / seed_tot * N_data;
		e.n   = new RooRealVar(Form("mc%s_frac", e.name.c_str()), Form("Yield of %s", e.name.c_str()),
		                       seed, 0.0, 10.0 * N_data);
		pdfList.add(*e.pdf);
		coefList.add(*e.n);
	}

	RooAddPdf model("model", "model", pdfList, coefList);

	// ------------------------------------------------------------------
	// Ranges and fit
	// ------------------------------------------------------------------
	x.setRange("FULL",     range_low, range_high);
	x.setRange("Sulfur_2", S_low, S_high);

	RooFitResult *fitres = model.fitTo(data, Range("FULL"), Save(), Extended(true),
	                                   PrintLevel(0), Strategy(2));

	fitres->Print("v");
	cout << "status  = " << fitres->status()  << endl;
	cout << "covQual = " << fitres->covQual() << endl;
	cout << "EDM     = " << fitres->edm()     << endl;
	cout << "chi2/ndf (x frame, see below) will be printed after plotting" << endl;

	// ------------------------------------------------------------------
	// Purity / contamination in the window "Sulfur_2"
	// ------------------------------------------------------------------
	RooArgList denomArgs;
	string denomExpr;
	for (size_t i = 0; i < el.size(); ++i) {
		el[i].I = el[i].pdf->createIntegral(x, NormSet(x), Range("Sulfur_2"));
		denomExpr += Form("@%zu*@%zu", 2*i, 2*i+1);
		if (i < el.size()-1) denomExpr += " + ";
		denomArgs.add(*el[i].n);
		denomArgs.add(*el[i].I);
	}
	RooFormulaVar Ntot_S("Ntot_S", denomExpr.c_str(), denomArgs);

	Elem *eS = nullptr;
	for (auto &e : el) if (e.name == "S") eS = &e;

	cout << "\n=== Contaminations in S window [" << S_low << ", " << S_high << "] ===" << endl;
	vector<pair<double,double>> conts;   // (value, error) of every non-S species, in order
	for (auto &e : el) {
		if (e.name == "S") continue;
		RooFormulaVar cont(Form("cont_%s", e.name.c_str()), "(@0*@1)/@2",
		                   RooArgList(*e.n, *e.I, Ntot_S));
		double val = cont.getVal();
		double err = cont.getPropagatedError(*fitres);
		conts.push_back({val, err});
		cout << "Contamination of " << e.name << " in S = " << val << " ± " << err << endl;
	}

	RooFormulaVar purity_S("purity_S", "(@0*@1)/@2", RooArgList(*eS->n, *eS->I, Ntot_S));
	double purity     = purity_S.getVal();
	double purity_err = purity_S.getPropagatedError(*fitres);
	cout << "Purity of S          = " << purity << " ± " << purity_err << endl;
	cout << "Sum of contaminations = " << (1.0 - purity) << " ± " << purity_err << endl;

	// ------------------------------------------------------------------
	// Plot
	// ------------------------------------------------------------------
	TCanvas *c1 = new TCanvas("c1", "c1", 1300, 980);
	c1->SetLogy();
	RooPlot *dframe = x.frame(range_low, range_high, ax->GetNbins() > 0 ? bHigh - bLow + 1 : 100);
	dframe->SetTitle(" ");

	data.plotOn(dframe);
	model.plotOn(dframe, Name("modelmc_postfit"), LineColor(kRed), LineWidth(4));
	for (auto &e : el) {
		int col = TColor::GetColor(e.color.c_str());
		model.plotOn(dframe, Components(*e.pdf), Name(Form("%s_postfit", e.pdf->GetName())),
		             LineStyle(1), LineColor(col), FillColor(col), FillStyle(3003),
		             DrawOption("FL"), RooFit::VLines(), MoveToBack());
	}
	double chi2ndf = dframe->chiSquare("modelmc_postfit", "h_data", (int)(el.size()) - 1);
	cout << "chi2/ndf = " << chi2ndf << endl;
	dframe->Draw();

	TLine *lmax = new TLine(S_high, dframe->GetMinimum(), S_high, dframe->GetMaximum());
	lmax->SetLineColor(kGray+4); lmax->SetLineStyle(9); lmax->SetLineWidth(2); lmax->Draw();
	TLine *lmin = new TLine(S_low, dframe->GetMinimum(), S_low, dframe->GetMaximum());
	lmin->SetLineColor(kGray+4); lmin->SetLineStyle(9); lmin->SetLineWidth(2); lmin->Draw();

	TLegend *l = new TLegend(0.82, 0.33, 0.99, 0.9);
	for (auto &e : el)
		l->AddEntry(Form("%s_postfit", e.pdf->GetName()), e.label.c_str(), "f");
	l->AddEntry("modelmc_postfit", "MC total ", "l");
	l->AddEntry("data", "Flight data", "epl");
	l->SetLineColor(0);
	l->SetTextSize(0.035);
	l->Draw("same");

	TLatex latex; latex.SetTextSize(0.04); latex.SetNDC(true); latex.SetTextFont(42);
	latex.DrawLatex(0.5, 0.88, title);

	// ------------------------------------------------------------------
	// Output files
	// ------------------------------------------------------------------
	if (writefile) {
		TString outfileName = "contaminazioni_sulfur_10bins_ago26.txt";
		ofstream outFile(outfileName.Data(), ios::app);
		if (!outFile.is_open()) {
			cerr << "ERROR: unable to open " << outfileName << endl;
		} else {
			outFile << fixed << setprecision(6);
			outFile << iBin << " " << E << " ";
			for (auto &c : conts) outFile << c.first << " " << c.second << " ";
			outFile << purity << " " << purity_err << " " << range_low << endl;
			outFile.close();
			cout << "TXT saved: " << outfileName << endl;
			cout << "Columns: iBin E ";
			for (auto &e : el) if (e.name != "S") cout << e.name << " err ";
			cout << "purity err range_low" << endl;
		}

		TString plot_name = "PLOTS/template_sulfur_10bins_ago26_" + TString(hname.c_str())
		                    + Form("_from%.1f", range_low) + ".png";
		cout << "Trying to save: " << plot_name << endl;
		c1->SaveAs(plot_name.Data());
		cout << "PNG saved." << endl;
	}
}