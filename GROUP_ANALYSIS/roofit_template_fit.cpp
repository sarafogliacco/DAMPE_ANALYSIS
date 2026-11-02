// =====================================================================
// Template fit of the PSD charge spectrum and background study for the
// GROUP  S - Cl - Ar - K - Ca   (Z = 16 ... 20)
//
// Signal      = every species of the group, integrated in the group window
//               [cut_low (lower edge of S) , cut_high_grp (upper edge of Ca)]
// Background  = every other species (C ... P below the window, Sc above it)
// Purity      = N_group(window) / N_total(window)
//
// Usage (root -l):
//   .L roofit_template_fit_SCa_group.C
//   roofit_template_fit_SCa_group(iBin, iRebin, writefile, addApprox, addSc,
//                                 approxVariant, fixFarSpecies)
//
//   approxVariant : how the missing odd-Z templates are built (use it for the
//                   systematic uncertainty, quote the spread of the background)
//        0 = shift of one neighbour (your original choice)
//        1 = shift of the other neighbour
//        2 = average of the two shifted neighbours
// =====================================================================
#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLine.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooFormulaVar.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooArgList.h"
#include "RooArgSet.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

using namespace RooFit;
using namespace std;

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------
struct SpeciesInfo {
	string name;      // C, O, Ne ...
	string label;     // legend
	string file;      // MC file ("" = template built from neighbours)
	string colorHex;
	int    fillStyle;
	bool   approx;    // template built by shifting neighbours
	bool   inGroup;   // part of the S-Ca signal group
	bool   fixedZero;
	TH1D*  h;
	RooDataHist* dh;
	RooHistPdf*  pdf;
	RooRealVar*  yield;
	RooAbsReal*  integ;   // fraction of the normalised pdf inside the group window
	RooProduct*  nwin;    // yield * integ = expected events inside the window
};

static TH1D* loadHist(const string& fname, const string& hname, int rebin)
{
	TFile* f = TFile::Open(fname.c_str());
	if (!f || f->IsZombie()) { cerr << "ERROR: cannot open " << fname << endl; return nullptr; }
	TH1D* h = (TH1D*)f->Get(hname.c_str());
	if (!h) { cerr << "ERROR: histogram " << hname << " not found in " << fname << endl; return nullptr; }
	h->SetDirectory(0);          // survive closing the file
	f->Close();
	h->Rebin(rebin);
	return h;
}

// Shift a template by dx (in charge units) with linear interpolation
static TH1D* shiftHist(const TH1D* h, double dx, const char* name)
{
	TH1D* hs = (TH1D*)h->Clone(name);
	hs->SetDirectory(0);
	hs->Reset();
	const int nb = h->GetNbinsX();
	for (int i = 1; i <= nb; ++i) {
		double c = h->GetBinCenter(i) - dx;
		double v = 0.;
		if (c >= h->GetBinCenter(1) && c <= h->GetBinCenter(nb)) v = h->Interpolate(c);
		hs->SetBinContent(i, v > 0. ? v : 0.);
	}
	return hs;
}

// Shifted template rescaled to the integral of its source
static TH1D* shiftScaled(const TH1D* src, double dx, const char* name)
{
	TH1D* t = shiftHist(src, dx, name);
	if (t->Integral() > 0.) t->Scale(src->Integral() / t->Integral());
	return t;
}

// Resample a histogram on a new regular grid [xlo,xhi] with nb bins,
// redistributing content by overlap (assumes uniform density inside each source bin).
static TH1D* resample(const TH1D* h, const char* name, int nb, double xlo, double xhi)
{
	TH1D* r = new TH1D(name, h->GetTitle(), nb, xlo, xhi);
	r->SetDirectory(0);
	const TAxis* ax = h->GetXaxis();
	for (int j = 1; j <= h->GetNbinsX(); ++j) {
		const double c = h->GetBinContent(j);
		if (c == 0.) continue;
		const double a = ax->GetBinLowEdge(j), b = ax->GetBinUpEdge(j);
		for (int i = 1; i <= nb; ++i) {
			const double A = r->GetXaxis()->GetBinLowEdge(i), B = r->GetXaxis()->GetBinUpEdge(i);
			const double ov = std::min(b, B) - std::max(a, A);
			if (ov > 0.) r->AddBinContent(i, c * ov / (b - a));
		}
	}
	r->SetEntries(h->GetEntries());
	return r;
}

// ---------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------
void roofit_template_fit_SCa_group(int iBin = 4, int iRebin = 2, bool writefile = false,
                                   bool addApprox = true,     // approximate N, F, Al, P, Cl, K templates (no MC)
                                   bool addSc = true,         // approximate Sc template (Ca shifted by +1): background above the window
                                   int  approxVariant = 0,    // 0 / 1 / 2, see header
                                   bool fixFarSpecies = false)// true: C, O, Ne fixed to 0 (use with fit_low ~ 11)
{
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPadLeftMargin(0.09);
	gStyle->SetPadRightMargin(0.19);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadBottomMargin(0.09);

	if (iBin < 1 || iBin > 10) { cerr << "ERROR: iBin must be in 1..10" << endl; return; }
	if (approxVariant < 0 || approxVariant > 2) { cerr << "ERROR: approxVariant must be 0, 1 or 2" << endl; return; }

	// -----------------------------------------------------------------
	// Energy bin: 6 bins of 1/6 decade (100 GeV-1 TeV), 3 of 1/3 decade (1-10 TeV),
	// 1 of 1 decade (10-100 TeV). BGO_E = geometric centre, used for the window.
	// -----------------------------------------------------------------
	double E_min, E_max;
	if      (iBin <= 6) { E_min = pow(10., 2. + (iBin - 1) / 6.); E_max = pow(10., 2. + iBin / 6.); }
	else if (iBin <= 9) { E_min = pow(10., 3. + (iBin - 7) / 3.); E_max = pow(10., 3. + (iBin - 6) / 3.); }
	else                { E_min = 1.0e4;                          E_max = 1.0e5; }
	const double BGO_E = sqrt(E_min * E_max);

	auto fmtE = [](double E) -> TString { return (E < 999.5) ? Form("%.0f GeV", E) : Form("%.1f TeV", E / 1000.); };
	TString title = fmtE(E_min) + " < E_{BGO} < " + fmtE(E_max);

	// -----------------------------------------------------------------
	// Group window
	//  - lower edge : lower edge of the S window (your cut_low polynomial)
	//  - upper edge : upper edge of the Ca window.
	//    If you have a Ca parametrisation, set haveCaParam = true and fill ca_c[].
	//    Otherwise the S upper edge is moved by (peak_Ca - peak_S) of the MC templates
	//    (computed below, once the templates are loaded).
	// -----------------------------------------------------------------
	const double lE = log10(BGO_E);
	const double cut_high = 17.491863 + 3.08338485*lE - 0.976515*pow(lE,2) + 0.10604*pow(lE,3);
	const double cut_low    = 14.2341715 + 1.46564565*lE - 0.48318*pow(lE,2) + 0.0520877*pow(lE,3);

	const bool   haveCaParam = false;
	const double ca_c[4]     = {0., 0., 0., 0.};   // TODO: Ca upper-edge polynomial in log10(E): c0 + c1*lE + c2*lE^2 + c3*lE^3

	// Full histogram range and range actually used in the fit
	double       range_low  = 5.5;    
	double       range_high = 22.0;   
	const double fit_low    = 5.5;    // 5.5 = full fit; with fit_low > ~10.5 you can use fixFarSpecies = true
	const double fit_high   = 22.0;

	// -----------------------------------------------------------------
	// Species list (ordered by Z)
	// -----------------------------------------------------------------
	const string base = "/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/";
	const string bkg  = base + "BACKGROUND/out_root/";
	const string tail = "_SATcorr_STKcut1200_PSDcutXY_CORRECTED";

	const vector<string> grp = {"S", "Cl", "Ar", "K", "Ca"};
	auto isGroup = [&](const string& n) { return std::find(grp.begin(), grp.end(), n) != grp.end(); };

	vector<SpeciesInfo> sp;
	auto add = [&](const string& n, const string& l, const string& f, const string& col, int fs, bool approx) {
		SpeciesInfo s;
		s.name = n; s.label = l; s.file = f; s.colorHex = col; s.fillStyle = fs;
		s.approx = approx; s.inGroup = isGroup(n); s.fixedZero = false;
		s.h = nullptr; s.dh = nullptr; s.pdf = nullptr; s.yield = nullptr; s.integ = nullptr; s.nwin = nullptr;
		sp.push_back(s);
	};

	add("C",  "Carbon",    bkg + "MC_FTFP_CARBON"    + tail + "_mcfit.root", "#6a4c93", 3003, false);
	if (addApprox) add("N", "Nitrogen (approx)", "", "#b5651d", 3004, true);
	add("O",  "Oxygen",    bkg + "MC_FTFP_OXYGEN"    + tail + "_mcfit.root", "#006627", 3003, false);
	if (addApprox) add("F", "Fluorine (approx)", "", "#2ca02c", 3004, true);
	add("Ne", "Neon",      bkg + "MC_FTFP_NEON"      + tail + "_mcfit.root", "#0bb4ff", 3003, false);
	add("Na", "Sodium",    bkg + "MC_FTFP_SODIUM"    + tail + "_mcfit.root", "#50e991", 3002, false);
	add("Mg", "Magnesium", bkg + "MC_FTFP_MAGNESIUM" + tail + "_mcfit.root", "#e6d800", 3003, false);
	if (addApprox) add("Al", "Aluminium (approx)", "", "#8c564b", 3004, true);
	add("Si", "Silicon",   bkg + "MC_FTFP_SILICON"   + tail + "_mcfit.root", "#81008d", 3003, false);
	if (addApprox) add("P",  "Phosphorus (approx)", "", "#1f3a93", 3004, true);
	add("S",  "Sulfur",    base + "SULFUR_ANALYSIS/out_root/MC_FTFP_SULFUR" + tail + "_1PeV_mcfit.root", "#ff6b6b", 3001, false);
	if (addApprox) add("Cl", "Chlorine (approx)", "", "#00a6a6", 3004, true);
	add("Ar", "Argon",     base + "ARGON_ANALYSIS/out_root/MC_FTFP_ARGON" + tail + "_mcfit.root", "#a775f8", 3003, false);
	if (addApprox) add("K", "Potassium (approx)", "", "#c0c0c0", 3004, true);
	add("Ca", "Calcium",   base + "CALCIUM_ANALYSIS/out_root/MC_FTFP_CALCIUM" + tail + "_100GeV-100TeV_mcfit.root", "#ffa500", 3003, false);


	auto idx = [&](const string& n) -> int {
		for (size_t i = 0; i < sp.size(); ++i) if (sp[i].name == n) return (int)i;
		return -1;
	};

	if (idx("S") < 0 || idx("Ca") < 0) { cerr << "ERROR: S and Ca templates are required" << endl; return; }

	// -----------------------------------------------------------------
	// Load histograms
	// -----------------------------------------------------------------
	string hname = Form(iBin < 10 ? "h0%d" : "h%d", iBin);

	TH1D* h_data = loadHist(bkg + "Skim_120months_SATcorr_STKcut1200_PSDcutXY_all_12histo_v2.root", hname, iRebin);
	if (!h_data) return;

	for (auto& s : sp) {
		if (s.file.empty()) continue;
		s.h = loadHist(s.file, hname, iRebin);
		if (!s.h) return;
	}

	// ---- Binning diagnostics and common grid --------------------------
	cout << "Binning DATA: " << h_data->GetNbinsX() << " bins, ["
	     << h_data->GetXaxis()->GetXmin() << ", " << h_data->GetXaxis()->GetXmax()
	     << "], width " << h_data->GetBinWidth(1) << endl;
	for (auto& s : sp) {
		if (!s.h) continue;
		cout << "Binning " << s.name << ": " << s.h->GetNbinsX() << " bins, ["
		     << s.h->GetXaxis()->GetXmin() << ", " << s.h->GetXaxis()->GetXmax()
		     << "], width " << s.h->GetBinWidth(1) << endl;
	}

	range_high = std::min(range_high, h_data->GetXaxis()->GetXmax());
	const double wbin = h_data->GetBinWidth(1);
	// align range_low to a data bin edge so data counts stay integer (Poisson)
	range_low = h_data->GetXaxis()->GetBinLowEdge(h_data->FindBin(range_low + 1e-6));
	const int    nbG  = TMath::Nint((range_high - range_low) / wbin);
	const double xhiG = range_low + nbG * wbin;

	h_data = resample(h_data, "h_data_common", nbG, range_low, xhiG);
	for (auto& s : sp) {
		if (!s.h) continue;
		s.h = resample(s.h, ("h_mc_" + s.name + "_common").c_str(), nbG, range_low, xhiG);
	}
	cout << "Common grid: " << nbG << " bins in [" << range_low << ", " << xhiG << "]" << endl;

	// -----------------------------------------------------------------
	// Approximate templates for the species without MC
	//   dx = +1 : start from the neighbour with Z-1 ; dx = -1 : from the neighbour with Z+1
	// -----------------------------------------------------------------
	struct Approx { const char* miss; const char* srcA; double dxA; const char* srcB; double dxB; };
	const vector<Approx> apTab = {
		// missing  A (default)    B (alternative)
		{"N",  "O",  -1.0,   "C",  +1.0},
		{"F",  "Ne", -1.0,   "O",  +1.0},
		{"Al", "Mg", +1.0,   "Si", -1.0},
		{"P",  "Si", +1.0,   "S",  -1.0},
		{"Cl", "Ar", -1.0,   "S",  +1.0},
		{"K",  "Ar", +1.0,   "Ca", -1.0},
	};
	for (auto& a : apTab) {
		const int im = idx(a.miss);
		if (im < 0 || !sp[im].approx) continue;
		const int iA = idx(a.srcA);
		const int iB = a.srcB ? idx(a.srcB) : -1;
		if (iA < 0) continue;

		const string nm = string("h_mc_") + a.miss + "_approx";
		TH1D* t = nullptr;
		if (approxVariant == 1 && iB >= 0) {
			t = shiftScaled(sp[iB].h, a.dxB, nm.c_str());
		} else if (approxVariant == 2 && iB >= 0) {
			TH1D* ta = shiftScaled(sp[iA].h, a.dxA, (nm + "_a").c_str());
			TH1D* tb = shiftScaled(sp[iB].h, a.dxB, (nm + "_b").c_str());
			t = (TH1D*)ta->Clone(nm.c_str());
			t->SetDirectory(0);
			t->Add(ta, tb, 0.5, 0.5);
			delete ta; delete tb;
		} else {
			t = shiftScaled(sp[iA].h, a.dxA, nm.c_str());
		}
		sp[im].h = t;
	}

	// -----------------------------------------------------------------
	// Group window upper edge
	// -----------------------------------------------------------------
	const TH1D* hS  = sp[idx("S")].h;
	const TH1D* hCa = sp[idx("Ca")].h;
	const double pkS  = hS->GetBinCenter(hS->GetMaximumBin());
	const double pkCa = hCa->GetBinCenter(hCa->GetMaximumBin());

	//double cut_high_grp;
	//if (haveCaParam) cut_high_grp = ca_c[0] + ca_c[1]*lE + ca_c[2]*pow(lE,2) + ca_c[3]*pow(lE,3);
	//else             cut_high_grp = cut_high_S + (pkCa - pkS);   // fallback: shift S edge by the Ca-S peak distance

	//cout << "S window         : [" << cut_low << ", " << cut_high << "]" << endl;
	cout << "Peak S / Ca (MC) : " << pkS << " / " << pkCa << endl;
	cout << "GROUP window     : [" << cut_low << ", " << cut_high << "]"
	     << (haveCaParam ? "  (Ca parametrisation)" : "  (fallback: S edge + peak shift, replace with the Ca parametrisation)") << endl;
	if (cut_high >= xhiG) {
		cerr << "WARNING: group upper edge (" << cut_high << ") is beyond the histogram range (" << xhiG
		     << "): the window is truncated." << endl;
		//cut_high = xhiG - 1e-6;
	}

	// ---- Template diagnostics ------------------------------------------
	{
		const int w1 = h_data->FindBin(cut_low), w2 = h_data->FindBin(cut_high);
		cout << "DATA : integral=" << h_data->Integral()
		     << "  in window=" << h_data->Integral(w1, w2)
		     << "  peak@" << h_data->GetBinCenter(h_data->GetMaximumBin()) << endl;
		for (auto& s : sp) {
			cout << "TEMPL " << s.name << (s.approx ? " (approx)" : "") << ": integral=" << s.h->Integral()
			     << "  max=" << s.h->GetMaximum()
			     << "  peak@" << s.h->GetBinCenter(s.h->GetMaximumBin())
			     << "  frac in window=" << (s.h->Integral() > 0 ? s.h->Integral(w1, w2) / s.h->Integral() : 0.)
			     << endl;
		}
	}

	// Avoid exactly-empty template bins: floor RELATIVE to each template's maximum
	for (auto& s : sp) {
		const double fl = 1e-9 * s.h->GetMaximum();
		for (int i = 1; i <= s.h->GetNbinsX(); ++i)
			if (s.h->GetBinContent(i) <= fl) s.h->SetBinContent(i, fl);
	}

	// -----------------------------------------------------------------
	// RooFit objects
	// -----------------------------------------------------------------
	RooRealVar x("x", "PSD Charge ", range_low, xhiG);
	const double fit_hi_eff = std::min(fit_high, xhiG);
	x.setRange("FULL", range_low, xhiG);
	x.setRange("FITR", fit_low,   fit_hi_eff);
	x.setRange("grp",  cut_low,   cut_high);

	RooDataHist data("data", "data", x, h_data);

	// Seeds computed on the same range as the fit window
	const int b1 = h_data->FindBin(fit_low  + 1e-6);
	const int b2 = h_data->FindBin(fit_hi_eff - 1e-6);
	const double N_data = h_data->Integral(b1, b2);

	for (auto& s : sp)
		s.fixedZero = fixFarSpecies && (s.name == "C" || s.name == "O" || s.name == "Ne");

	// Seed from the DATA around each template peak (template integrals have very different scales)
	auto seedFor = [&](SpeciesInfo& s) -> double {
		if (s.fixedZero) return 0.;
		const double pk = s.h->GetBinCenter(s.h->GetMaximumBin());
		const int lo = h_data->FindBin(pk - 0.4), hi = h_data->FindBin(pk + 0.4);
		const double dwin = h_data->Integral(lo, hi);
		const double tw   = s.h->Integral(lo, hi);
		const double tt   = s.h->Integral(b1, b2);
		if (tw <= 0. || tt <= 0.) return 0.01 * N_data;
		return std::min(0.3 * dwin * tt / tw, 0.5 * N_data);
	};

	RooArgList pdfList, coefList;
	for (auto& s : sp) {
		s.dh  = new RooDataHist(("mc" + s.name).c_str(), ("mc " + s.name).c_str(), x, s.h);
		s.pdf = new RooHistPdf(("modelmc" + s.name).c_str(), ("modelmc" + s.name).c_str(), x, *s.dh);

		const double seed = seedFor(s);
		s.yield = new RooRealVar(("mc" + s.name + "_frac").c_str(), ("Yield of " + s.name).c_str(),
		                         seed, 0.0, 10.0 * N_data);
		if (s.fixedZero) { s.yield->setVal(0.); s.yield->setConstant(kTRUE); }

		pdfList.add(*s.pdf);
		coefList.add(*s.yield);
	}

	RooAddPdf model("model", "model", pdfList, coefList);

	// -----------------------------------------------------------------
	// Fit (extended, in FITR)
	// NB: the pdfs are normalised on the full x range, so yield * integral(window)
	//     is the number of events of that species inside the group window.
	// -----------------------------------------------------------------
	RooFitResult* fitres = model.fitTo(data, Range("FITR"), Extended(kTRUE), SumW2Error(false),
	                                   Strategy(2), Save(), PrintLevel(0));

	// -----------------------------------------------------------------
	// Purity / contaminations in the S-Ca group window
	// (RooProduct / RooAddition: no limit on the number of species,
	//  error propagation uses the full covariance matrix)
	// -----------------------------------------------------------------
	RooArgList allWin, sigWin;
	for (auto& s : sp) {
		s.integ = s.pdf->createIntegral(x, NormSet(x), Range("grp"));
		s.nwin  = new RooProduct(("Nwin_" + s.name).c_str(), ("Nwin " + s.name).c_str(),
		                         RooArgList(*s.yield, *s.integ));
		allWin.add(*s.nwin);
		if (s.inGroup) sigWin.add(*s.nwin);
	}
	RooAddition Ntot_sample("Ntot_sample", "Ntot in window", allWin);
	RooAddition Nsig_sample("Nsig_sample", "Nsig (S-Ca) in window", sigWin);

	RooFormulaVar purity_sample("purity_sample", "@0/@1", RooArgList(Nsig_sample, Ntot_sample));
	const double purity     = purity_sample.getVal();
	const double purity_err = purity_sample.getPropagatedError(*fitres);   // 1-purity has the same error

	// fraction of a species in the window sample, with error
	auto fracOf = [&](const SpeciesInfo& s, RooAbsReal& denom, double& val, double& err) {
		RooFormulaVar f(("frac_" + s.name).c_str(), "@0/@1", RooArgList(*s.nwin, denom));
		val = f.getVal();
		err = f.getPropagatedError(*fitres);
	};

	cout << endl << "=== Window composition: contaminations (species OUTSIDE S-Ca) ===" << endl;
	for (auto& s : sp) {
		if (s.inGroup) continue;
		double v, e; fracOf(s, Ntot_sample, v, e);
		cout << "  " << setw(3) << s.name << (s.approx ? " (approx)" : "         ") << " : " << v << " ± " << e << endl;
	}
	cout << "=== Window composition: group members (fraction of total) ===" << endl;
	for (auto& s : sp) {
		if (!s.inGroup) continue;
		double v, e; fracOf(s, Ntot_sample, v, e);
		cout << "  " << setw(3) << s.name << (s.approx ? " (approx)" : "         ") << " : " << v << " ± " << e << endl;
	}
	cout << "Purity of the S-Ca group      = " << purity << " ± " << purity_err << endl;
	cout << "Total contamination (1-purity)= " << (1.0 - purity) << " ± " << purity_err << endl;
	cout << "N total in window (fit)       = " << Ntot_sample.getVal() << endl;

	fitres->Print("v");
	cout << "status  = " << fitres->status()  << endl;
	cout << "covQual = " << fitres->covQual() << endl;
	cout << "EDM     = " << fitres->edm()     << endl;

	// -----------------------------------------------------------------
	// Plot
	// -----------------------------------------------------------------
	TCanvas* c1 = new TCanvas("c1", "c1", 1800, 980);
	c1->SetLogy();

	RooPlot* dframe = x.frame(range_low, xhiG, nbG);
	dframe->SetTitle(" ");

	data.plotOn(dframe, Name("data_pts"));
	model.plotOn(dframe, Name("modelmc_postfit"), LineColor(kRed), LineWidth(4));

	for (auto& s : sp) {
		if (s.fixedZero) continue;
		Int_t col = TColor::GetColor(s.colorHex.c_str());
		model.plotOn(dframe, Components(*s.pdf), Name(("modelmc" + s.name + "_postfit").c_str()),
		             LineStyle(1), LineColor(col), FillColor(col), FillStyle(s.fillStyle),
		             DrawOption("FL"), RooFit::VLines(), MoveToBack());
	}
	dframe->SetMinimum(0.5);
	dframe->Draw();

	TLine* lmin = new TLine(cut_low, dframe->GetMinimum(), cut_low, dframe->GetMaximum());
	lmin->SetLineColor(kGray+4); lmin->SetLineStyle(9); lmin->SetLineWidth(2); lmin->Draw();

	TLine* lmax = new TLine(cut_high, dframe->GetMinimum(), cut_high, dframe->GetMaximum());
	lmax->SetLineColor(kGray+4); lmax->SetLineStyle(9); lmax->SetLineWidth(2); lmax->Draw();

	TLine* lfit1 = new TLine(fit_low, dframe->GetMinimum(), fit_low, dframe->GetMaximum());
	lfit1->SetLineColor(kBlack); lfit1->SetLineStyle(2); lfit1->SetLineWidth(1); lfit1->Draw();

	TLegend* l = new TLegend(0.82, 0.20, 0.99, 0.93);
	for (auto& s : sp) {
		if (s.fixedZero) continue;
		l->AddEntry(("modelmc" + s.name + "_postfit").c_str(), s.label.c_str(), "f");
	}
	l->AddEntry("modelmc_postfit", "MC total ", "l");
	l->AddEntry("data_pts", "Flight data", "epl");
	l->SetLineColor(0);
	l->SetTextSize(0.027);
	l->Draw("same");

	TLatex latex; latex.SetTextSize(0.04); latex.SetNDC(true); latex.SetTextFont(42);
	latex.DrawLatex(0.45, 0.88, title);
	latex.SetTextSize(0.03);
	latex.DrawLatex(0.45, 0.83, Form("S-Ca window [%.2f, %.2f]   purity = %.4f #pm %.4f", cut_low, cut_high, purity, purity_err));

	// -----------------------------------------------------------------
	// Output files
	// -----------------------------------------------------------------
	if (writefile) {

		// CSV: one row per energy bin, appended; header written if the file is new/empty.
		// One file per approxVariant so that the systematic spread can be computed afterwards.
		TString outfileName = Form("txt/contaminazioni_SCa_group_var%d_21oct26.csv", approxVariant);
		const vector<string> csvOrder = {"C","N","O","F","Ne","Na","Mg","Al","Si","P"};   // species outside S-Ca

		bool needHeader = true;
		{
			std::ifstream test(outfileName.Data());
			if (test.is_open() && test.peek() != std::ifstream::traits_type::eof()) needHeader = false;
		}

		std::ofstream outFile(outfileName.Data(), std::ios::app);
		if (!outFile.is_open()) {
			cerr << "ERROR: unable to open " << outfileName << endl;
		} else {
			if (needHeader) {
				outFile << "#E_min,E_max,win_low,win_high";
				for (auto& n : csvOrder) outFile << "," << n << "," << n << "_err";
				outFile << ",Tot,Tot_err,purity,fit_status,covQual" << endl;
			}

			outFile << std::scientific << std::setprecision(6);
			outFile << E_min << "," << E_max << "," << cut_low << "," << cut_high;

			for (auto& n : csvOrder) {
				double val = 0., err = 0.;          // 0,0 = species not in the fit (NOT "no contamination")
				const int k = idx(n);
				if (k >= 0 && !sp[k].fixedZero) fracOf(sp[k], Ntot_sample, val, err);
				outFile << "," << val << "," << err;
			}
			outFile << "," << (1.0 - purity) << "," << purity_err << "," << purity
			        << "," << fitres->status() << "," << fitres->covQual() << endl;
			outFile.close();
			cout << "CSV saved: " << outfileName << endl;
		}

		TString plot_name1 = Form("PLOTS/template_SCa_group_var%d_21oct26_%s.png", approxVariant, hname.c_str());
		cout << "Trying to save: " << plot_name1 << endl;
		c1->SaveAs(plot_name1.Data());
		cout << "PNG saved." << endl;
	}
}