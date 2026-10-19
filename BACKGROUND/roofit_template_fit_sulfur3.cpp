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

using namespace RooFit;
using namespace std;

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------
struct SpeciesInfo {
	string name;      // C, O, Ne ...
	string label;     // legend
	string file;      // MC file ("" = template built by shifting neighbours)
	string colorHex;
	int    fillStyle;
	bool   fixedZero;
	TH1D*  h;
	RooDataHist* dh;
	RooHistPdf*  pdf;
	RooRealVar*  yield;
	RooAbsReal*  integ;
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

/*
// Approximate missing odd-Z template: average of (lower +1) and (upper -1), same integral scale
static TH1D* approxTemplate(const TH1D* hLow, const TH1D* hUp, const char* name)
{
	TH1D* a = shiftHist(hLow, +1.0, Form("%s_a", name));
	TH1D* b = shiftHist(hUp,  -1.0, Form("%s_b", name));
	if (a->Integral() > 0) a->Scale(hLow->Integral() / a->Integral());
	if (b->Integral() > 0) b->Scale(hUp->Integral()  / b->Integral());
	TH1D* r = (TH1D*)a->Clone(name);
	r->SetDirectory(0);
	r->Add(a, b, 0.5, 0.5);
	delete a; delete b;
	return r;
}
*/



// ---------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------
void roofit_template_fit_sulfur(int iBin = 4, int iRebin = 2, bool writefile = false,
                                bool addApproxAlP = true,   // approximate Al and P templates (missing MC)
                                bool fixFarSpecies = false) // true: C, O, Ne fixed to 0 (use together with fit_low ~ 11)
{
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPadLeftMargin(0.09);
	gStyle->SetPadRightMargin(0.19);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadBottomMargin(0.09);

	double BGO_E = 0.0;
	TString title;

	// 10 bins (6 - 3 - 1)
	// NB: CONTROLLARE che energie e titoli siano coerenti (vedi bin 7, 8, 9, 10):
	// S_low/S_high dipendono da BGO_E, i titoli sono solo estetici.
	if(iBin == 1) { BGO_E = 121.1527961; title = "100 GeV < E_{BGO} < 178 GeV"  ; }
	if(iBin == 2) { BGO_E = 177.8277918; title = "178 GeV < E_{BGO} < 316 GeV"  ; }
	if(iBin == 3) { BGO_E = 261.015534 ; title = "316 GeV < E_{BGO} < 562 GeV"  ; }
	if(iBin == 4) { BGO_E = 383.1188748; title = "562 GeV < E_{BGO} < 1.0 TeV"  ; }
	if(iBin == 5) { BGO_E = 562.3413673; title = "1.0 TeV < E_{BGO} < 1.8 TeV"  ; }
	if(iBin == 6) { BGO_E = 825.4041434; title = "1.8 TeV < E_{BGO} < 3.2 TeV"  ; }
	if(iBin == 7) { BGO_E = 1467.79767 ; title = "3.2 TeV < E_{BGO} < 5.6 TeV"  ; }
	if(iBin == 8) { BGO_E = 3162.274615; title = "5.6 TeV < E_{BGO} < 10.0 TeV" ; }
	if(iBin == 9) { BGO_E = 6812.921547; title = "10.0 TeV < E_{BGO} < 31.6 TeV"; }
	if(iBin ==10) { BGO_E = 31622.7766 ; title = "31.6 TeV < E_{BGO} < 100.0 TeV";}
	if (BGO_E == 0.0) { cerr << "ERROR: iBin must be in 1..10" << endl; return; }

	// Bin edges in GeV: 6 bins of 1/6 decade (100 GeV-1 TeV), 3 of 1/3 decade (1-10 TeV), 1 of 1 decade (10-100 TeV).
	// (BGO_E above are the geometric centres of these bins; the old titles did not match them for iBin >= 4.)
	double E_min, E_max;
	if      (iBin <= 6) { E_min = pow(10., 2. + (iBin - 1) / 6.); E_max = pow(10., 2. + iBin / 6.); }
	else if (iBin <= 9) { E_min = pow(10., 3. + (iBin - 7) / 3.); E_max = pow(10., 3. + (iBin - 6) / 3.); }
	else                { E_min = 1.0e4;                          E_max = 1.0e5; }
	auto fmtE = [](double E) -> TString { return (E < 999.5) ? Form("%.0f GeV", E) : Form("%.1f TeV", E / 1000.); };
	title = fmtE(E_min) + " < E_{BGO} < " + fmtE(E_max);

	// SULFUR 09/2026 : charge window used for the purity
	const double lE = log10(BGO_E);
	const double S_high =  15.5823475+ 1.0482565*lE  - 0.26595*pow(lE,2) +  0.0260414*pow(lE,3);
	//15.915774  + 0.6569678*lE   - 0.167878*pow(lE,2) + 0.0179829*pow(lE,3);
	const double S_low  = 14.8049915 + 0.8536717*lE -  0.26595*pow(lE,2) + 0.0260414*pow(lE,3);
		//15.1162445 + 0.53318365*lE  - 0.167878*pow(lE,2) + 0.0179829*pow(lE,3);

	cout << "High limit: " << S_high << endl;
	cout << "Low limit: "  << S_low  << endl;

	// Full histogram range and range actually used in the fit
	double       range_low  = 5.5;   // verra allineato a un bordo di bin dei dati
	const double range_high = 21.0;
	const double fit_low    = 5.5;    // inizio finestra di fit (5.5 = fit completo; con fit_low>~10.5 puoi usare fixFarSpecies=true)
	const double fit_high   = 21.0;

	// -----------------------------------------------------------------
	// Species list (ordered by Z)
	// -----------------------------------------------------------------
	const string base = "/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/";
	const string bkg  = base + "BACKGROUND/out_root/";
	const string tail = "_SATcorr_STKcut1200_PSDcutXY_CORRECTED";

	vector<SpeciesInfo> sp;
	auto add = [&](const string& n, const string& l, const string& f, const string& col, int fs) {
		SpeciesInfo s;
		s.name = n; s.label = l; s.file = f; s.colorHex = col; s.fillStyle = fs;
		s.fixedZero = false; s.h = nullptr; s.dh = nullptr; s.pdf = nullptr; s.yield = nullptr; s.integ = nullptr;
		sp.push_back(s);
	};

	add("C",  "Carbon",    bkg + "MC_FTFP_CARBON"    + tail + "_mcfit.root", "#6a4c93", 3003);
	if (addApproxAlP) add("N", "Nitrogen (approx)", "", "#b5651d", 3004);
	add("O",  "Oxygen",    bkg + "MC_FTFP_OXYGEN"    + tail + "_mcfit.root", "#006627", 3003);
	if (addApproxAlP) add("F", "Fluorine (approx)", "", "#2ca02c", 3004);
	add("Ne", "Neon",      bkg + "MC_FTFP_NEON"      + tail + "_mcfit.root", "#0bb4ff", 3003);
	add("Na", "Sodium",    bkg + "MC_FTFP_SODIUM"    + tail + "_mcfit.root", "#50e991", 3002);
	add("Mg", "Magnesium", bkg + "MC_FTFP_MAGNESIUM" + tail + "_mcfit.root", "#e6d800", 3003);
	if (addApproxAlP) add("Al", "Aluminium (approx)", "", "#8c564b", 3004);
	add("Si", "Silicon",   bkg + "MC_FTFP_SILICON"   + tail + "_mcfit.root", "#81008d", 3003);
	if (addApproxAlP) add("P",  "Phosphorus (approx)", "", "#1f3a93", 3004);
	add("S",  "Sulfur",    base + "SULFUR_ANALYSIS/out_root/MC_FTFP_SULFUR" + tail + "_1PeV_mcfit.root", "#ff6b6b", 3001);
	//if (addApproxAlP) add("Cl", "Chlorine (approx)", "", "#00a6a6", 3004);
	add("Cl", "Chlorine",  base + "secondaries/out_root/MC_FTFP_CLORO" + tail + ".root", "#00a6a6", 3003);
	add("Ar", "Argon",     base + "ARGON_ANALYSIS/out_root/MC_FTFP_ARGON" + tail + "_mcfit.root", "#a775f8", 3003);
	//if (addApproxAlP) add("K", "Potassium (approx)", "", "#c0c0c0", 3004);
	add("K",  "Potassium", base + "secondaries/out_root/MC_FTFP_POTASSIUM" + tail + ".root", "#c0c0c0", 3003);
	add("Ca", "Calcium",   base + "CALCIUM_ANALYSIS/out_root/MC_FTFP_CALCIUM" + tail + "_100GeV-100TeV_mcfit.root", "#ffa500", 3003);



	auto find = [&](const string& n) -> int {
		for (size_t i = 0; i < sp.size(); ++i) if (sp[i].name == n) return (int)i;
		return -1;
	};

	// -----------------------------------------------------------------
	// Load histograms
	// -----------------------------------------------------------------
	// Histogram naming: h01..h09, h10, h11, ...
	string hname = (iBin < 10) ? Form("h0%d", iBin) : Form("h%d", iBin);

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

	// Common grid: data bin width, from range_low to range_high
	const double wbin = h_data->GetBinWidth(1);
	// align range_low to a data bin edge so data counts stay integer (Poisson)
	range_low = h_data->GetXaxis()->GetBinLowEdge(h_data->FindBin(range_low + 1e-6));
	const int    nbG  = TMath::Nint((range_high - range_low) / wbin);
	const double xhiG = range_low + nbG * wbin;
	TH1D* h_data_old = h_data;
	h_data = resample(h_data_old, "h_data_common", nbG, range_low, xhiG);
	for (auto& s : sp) {
		if (!s.h) continue;
		TH1D* old = s.h;
		s.h = resample(old, ("h_mc_" + s.name + "_common").c_str(), nbG, range_low, xhiG);
	}
	cout << "Common grid: " << nbG << " bins in [" << range_low << ", " << xhiG << "]" << endl;

		if (addApproxAlP) {
		// dx = +1: parte dal vicino con Z-1 ; dx = -1: parte dal vicino con Z+1
		struct Approx { const char* miss; const char* src; double dx; };
		const vector<Approx> ap = {
			{"N",  "O",  -1.0},   // oppure {"N",  "O",  -1.0}
			{"F",  "Ne", -1.0},   // oppure {"F",  "O",  +1.0}
			{"Al", "Mg", +1.0},   // oppure {"Al", "Si", -1.0}
			{"P",  "Si", +1.0},   // oppure {"P",  "S",  -1.0}
			{"Cl", "Ar", -1.0},   // oppure {"Cl", "S",  +1.0}
			{"K",  "Ar", +1.0}    // oppure {"K",  "Ca", -1.0}
		};
		for (auto& a : ap) {
			const int im = find(a.miss), is = find(a.src);
			if (im < 0 || is < 0) continue;        // specie non presente: salta
			TH1D* t = shiftHist(sp[is].h, a.dx, Form("h_mc_%s_approx", a.miss));
			if (t->Integral() > 0) t->Scale(sp[is].h->Integral() / t->Integral());
			sp[im].h = t;
		}
	}


/*
	if (addApproxAlP) {
		sp[find("N")].h  = approxTemplate(sp[find("C")].h,  sp[find("O")].h,  "h_mc_N_approx");
		sp[find("F")].h  = approxTemplate(sp[find("O")].h,  sp[find("Ne")].h, "h_mc_F_approx");
		sp[find("Al")].h = approxTemplate(sp[find("Mg")].h, sp[find("Si")].h, "h_mc_Al_approx");
		sp[find("P")].h  = approxTemplate(sp[find("Si")].h, sp[find("S")].h,  "h_mc_P_approx");
		sp[find("Cl")].h = approxTemplate(sp[find("S")].h,  sp[find("Ar")].h, "h_mc_Cl_approx");
		sp[find("K")].h  = approxTemplate(sp[find("Ar")].h, sp[find("Ca")].h, "h_mc_K_approx");
	}
*/
	// ---- Template diagnostics (scale, peak position, fraction in S window) ----
	{
		const int w1 = h_data->FindBin(S_low), w2 = h_data->FindBin(S_high);
		cout << "DATA : integral=" << h_data->Integral()
		     << "  peak@" << h_data->GetBinCenter(h_data->GetMaximumBin()) << endl;
		for (auto& s : sp) {
			cout << "TEMPL " << s.name << ": integral=" << s.h->Integral()
			     << "  max=" << s.h->GetMaximum()
			     << "  peak@" << s.h->GetBinCenter(s.h->GetMaximumBin())
			     << "  frac in S window=" << (s.h->Integral() > 0 ? s.h->Integral(w1, w2) / s.h->Integral() : 0.)
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
	x.setRange("Sulfur_2", S_low, S_high);

	RooDataHist data("data", "data", x, h_data);

	// Seeds computed on the same range as the fit window
	const int b1 = h_data->FindBin(fit_low  + 1e-6);
	const int b2 = h_data->FindBin(fit_hi_eff - 1e-6);
	const double N_data = h_data->Integral(b1, b2);

	for (auto& s : sp) {
		s.fixedZero = fixFarSpecies && (s.name == "C" || s.name == "O" || s.name == "Ne");
	}

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

		double seed = seedFor(s);
		s.yield = new RooRealVar(("mc" + s.name + "_frac").c_str(), ("Yield of " + s.name).c_str(),
		                         seed, 0.0, 10.0 * N_data);
		if (s.fixedZero) { s.yield->setVal(0.); s.yield->setConstant(kTRUE); }

		pdfList.add(*s.pdf);
		coefList.add(*s.yield);
	}

	RooAddPdf model("model", "model", pdfList, coefList);

	// -----------------------------------------------------------------
	// Fit (extended, only in the window around S)
	// -----------------------------------------------------------------
	RooFitResult* fitres = model.fitTo(data, Range("FITR"), Extended(kTRUE), SumW2Error(false), Strategy(2), Save(), PrintLevel(0));

	// -----------------------------------------------------------------
	// Purity / contaminations in the sulfur window
	// -----------------------------------------------------------------
	for (auto& s : sp)
		s.integ = s.pdf->createIntegral(x, NormSet(x), Range("Sulfur_2"));

	RooArgList denomArgs;
	string denomExpr;
	for (size_t i = 0; i < sp.size(); ++i) {
		denomExpr += "@" + to_string(2*i) + "*@" + to_string(2*i + 1);
		if (i < sp.size() - 1) denomExpr += " + ";
		denomArgs.add(*sp[i].yield);
		denomArgs.add(*sp[i].integ);
	}
	RooFormulaVar Ntot_S("Ntot_S", denomExpr.c_str(), denomArgs);

	for (auto& s : sp) {
		if (s.name == "S") continue;
		RooFormulaVar cont("cont_" + TString(s.name.c_str()), "(@0*@1)/@2",
		                   RooArgList(*s.yield, *s.integ, Ntot_S));
		cout << "Contamination of " << s.name << " in S = "
		     << cont.getVal() << " ± " << cont.getPropagatedError(*fitres) << endl;
	}

	const int iS = find("S");
	RooFormulaVar purity_S("purity_S", "(@0*@1)/@2", RooArgList(*sp[iS].yield, *sp[iS].integ, Ntot_S));
	const double purity     = purity_S.getVal();
	const double purity_err = purity_S.getPropagatedError(*fitres);

	cout << "Purity S              = " << purity << " ± " << purity_err << endl;
	cout << "Sum of contaminations = " << (1.0 - purity) << " ± " << purity_err << endl;

	fitres->Print("v");
	cout << "status  = " << fitres->status()  << endl;
	cout << "covQual = " << fitres->covQual() << "   (atteso 3)" << endl;
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

	TLine *lmax = new TLine(S_high, dframe->GetMinimum(), S_high, dframe->GetMaximum());
	lmax->SetLineColor(kGray+4); lmax->SetLineStyle(9); lmax->SetLineWidth(2); lmax->Draw();

	TLine *lmin = new TLine(S_low, dframe->GetMinimum(), S_low, dframe->GetMaximum());
	lmin->SetLineColor(kGray+4); lmin->SetLineStyle(9); lmin->SetLineWidth(2); lmin->Draw();

	// Fit window
	TLine *lfit1 = new TLine(fit_low, dframe->GetMinimum(), fit_low, dframe->GetMaximum());
	lfit1->SetLineColor(kBlack); lfit1->SetLineStyle(2); lfit1->SetLineWidth(1); lfit1->Draw();

	TLegend *l = new TLegend(0.82, 0.25, 0.99, 0.92);
	for (auto& s : sp) {
		if (s.fixedZero) continue;
		l->AddEntry(("modelmc" + s.name + "_postfit").c_str(), s.label.c_str(), "f");
	}
	l->AddEntry("modelmc_postfit", "MC total ", "l");
	l->AddEntry("data_pts", "Flight data", "epl");
	l->SetLineColor(0);
	l->SetTextSize(0.03);
	l->Draw("same");

	TLatex latex; latex.SetTextSize(0.04); latex.SetNDC(true); latex.SetTextFont(42);
	latex.DrawLatex(0.5, 0.88, title);

	// -----------------------------------------------------------------
	// Output files
	// -----------------------------------------------------------------
	if (writefile) {

		// ---- CSV: one row per energy bin, appended; header written if the file is new/empty ----
		TString outfileName = "txt/contaminazioni_sulfur_13oct26_v2.csv";
		const vector<string> csvOrder = {"C","O","N","F","Ne","Na","Mg","Al","Si","P","Cl","Ar","K","Ca"};

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
				outFile << "#E_min,E_max";
				for (auto& n : csvOrder) outFile << "," << n << "," << n << "_err";
				outFile << ",Tot,Tot_err" << endl;
			}

			outFile << std::scientific << std::setprecision(6);
			outFile << E_min << "," << E_max;

			for (auto& n : csvOrder) {
				double val = 0., err = 0.;          // 0,0 if the species is not in the fit
				const int k = find(n);
				if (k >= 0 && !sp[k].fixedZero) {
					RooFormulaVar cont("cont_" + TString(n.c_str()), "(@0*@1)/@2",
					                   RooArgList(*sp[k].yield, *sp[k].integ, Ntot_S));
					val = cont.getVal();
					err = cont.getPropagatedError(*fitres);
				}
				outFile << "," << val << "," << err;
			}
			outFile << "," << (1.0 - purity) << "," << purity_err << endl;
			outFile.close();
			cout << "CSV saved: " << outfileName << endl;
		}

		TString plot_name1 = "PLOTS/SULFUR/template_sulfur_13oct26_v2_" + TString(hname.c_str()) + ".png";
		cout << "Trying to save: " << plot_name1 << endl;
		c1->SaveAs(plot_name1.Data());
		cout << "PNG saved." << endl;
	}
}