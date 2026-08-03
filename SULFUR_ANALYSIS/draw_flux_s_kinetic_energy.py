import sys
from ROOT import gStyle, TGraph, TGraphErrors, TGraphAsymmErrors, TLatex, TLegend, TLine, TBox, TCanvas, gPad, kBlack, kGray, kRed, kBlue, kGreen, kAzure, kOrange, kMagenta, TPad, TGaxis
import array as ary
import numpy as np
import math

def make_flux_graph_DAMPE(filename, color, marker, size, n_drop_last=0):
    Emean    = np.loadtxt(filename, skiprows=1, usecols=(0,), unpack=True)
    #Elow     = np.loadtxt(filename, skiprows=1, usecols=(1,), unpack=True)
    #Eup      = np.loadtxt(filename, skiprows=1, usecols=(2,), unpack=True)
    Flux_2   = np.loadtxt(filename, skiprows=1, usecols=(1,), unpack=True) 
    Flux_stat= np.loadtxt(filename, skiprows=1, usecols=(2,), unpack=True)
    #Flux_sysA= np.loadtxt(filename, skiprows=1, usecols=(5,), unpack=True)
    #Flux_sysH= np.loadtxt(filename, skiprows=1, usecols=(6,), unpack=True)

    if n_drop_last > 0:
        Emean     = Emean[:-n_drop_last]
        Flux_2    = Flux_2[:-n_drop_last]
        Flux_stat = Flux_stat[:-n_drop_last]

    Flux    = (Flux_2) * Emean**2.7 * 0.8
    Flux_low= (Flux_stat) * Emean**2.7 *0.8
    Flux_up = (Flux_stat) * Emean**2.7 *0.8
    #syst    = (Flux_sysA) * Emean**2.7
    #syst_had= (Flux_sysH) * Emean**2.7

    null = np.zeros(len(Emean))

    gr = TGraphAsymmErrors(len(Emean), Emean, Flux, null, null, Flux_low, Flux_up)
    gr.SetLineColor(color)
    gr.SetMarkerColor(color)
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(size)

    # NOTE: the original code also built gr_in / gr_out (systematic-band
    # graphs) using an undefined "Eup" variable (its source column was
    # commented out) and never actually drew them anywhere. That dead code
    # has been removed. If you want to re-enable the systematic band, load
    # Elow/Eup/Flux_sysA/Flux_sysH again and rebuild gr_in/gr_out from them.

    return gr

def make_flux_graph_AMS(filename, color, marker, size, Z, m_S, gamma=2.7):
    # File:
    # Quantity R-GV bin-range Flux-(GV m^2 s sr)^-1 stat-err syst-err

    R        = np.loadtxt(filename, skiprows=1, usecols=(1,), unpack=True)
    Flux_R   = np.loadtxt(filename, skiprows=1, usecols=(4,), unpack=True)
    stat_err = np.loadtxt(filename, skiprows=1, usecols=(5,), unpack=True)
    sys_err  = np.loadtxt(filename, skiprows=1, usecols=(6,), unpack=True)

    # Rigidity -> kinetic energy
    Emean = np.sqrt((Z * R)**2 + m_S**2) - m_S
    dR_dE = np.sqrt((Z * R)**2 + m_S**2) / (Z**2 * R)

    Flux = Flux_R * dR_dE * Emean**gamma

    Flux_low = np.sqrt(stat_err**2 + sys_err**2) * dR_dE * Emean**gamma
    Flux_up  = Flux_low.copy()

    null = np.zeros(len(Emean))

    gr = TGraphAsymmErrors(
        len(Emean),
        Emean,
        Flux,
        null, null,
        Flux_low,
        Flux_up
    )

    gr.SetLineColor(color)
    gr.SetMarkerColor(color)
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(size)

    return gr

def make_flux_graph_TRACER(filename, color, marker, size, A, gamma=2.7):

    Emean_n = np.loadtxt(filename, skiprows=1, usecols=(3,), delimiter=",", unpack=True)
    Flux_n  = np.loadtxt(filename, skiprows=1, usecols=(5,), delimiter=",", unpack=True)
    Err_up  = np.loadtxt(filename, skiprows=1, usecols=(6,), delimiter=",", unpack=True)
    Err_dw  = np.loadtxt(filename, skiprows=1, usecols=(7,), delimiter=",", unpack=True)

    Emean = Emean_n * A

    Flux = Flux_n * Emean_n**gamma * A**(gamma-1)
    Flux_low = Err_dw * Emean_n**gamma * A**(gamma-1)
    Flux_up  = Err_up * Emean_n**gamma * A**(gamma-1)

    null = np.zeros(len(Emean))

    gr = TGraphAsymmErrors(
        len(Emean),
        Emean,
        Flux,
        null, null,
        Flux_low,
        Flux_up
    )

    gr.SetLineColor(color)
    gr.SetMarkerColor(color)
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(size)

    return gr

def make_flux_graph_TRACER99(filename, color, marker, size, A, gamma=2.7):

    Emean_n = np.loadtxt(filename, skiprows=1, usecols=(1,), unpack=True)
    Flux_n  = np.loadtxt(filename, skiprows=1, usecols=(3,), unpack=True)
    Err     = np.loadtxt(filename, skiprows=1, usecols=(4,), unpack=True)

    Emean = Emean_n * A

    Flux = Flux_n * Emean_n**gamma * A**(gamma-1)
    Flux_err = Err * Emean_n**gamma * A**(gamma-1)

    null = np.zeros(len(Emean))

    gr = TGraphAsymmErrors(
        len(Emean),
        Emean,
        Flux,
        null, null,
        Flux_err,
        Flux_err
    )

    gr.SetLineColor(color)
    gr.SetMarkerColor(color)
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(size)

    return gr

def make_flux_graph_CRISIS(filename, color, marker, size, A, gamma=2.7):

    Emean_n = np.loadtxt(filename, skiprows=1, usecols=(1,), unpack=True)
    Flux_n  = np.loadtxt(filename, skiprows=1, usecols=(4,), unpack=True)
    Err     = np.loadtxt(filename, skiprows=1, usecols=(5,), unpack=True)

    Emean = Emean_n * A

    Flux = Flux_n * Emean_n**gamma * A**(gamma-1)
    Flux_err = Err * Emean_n**gamma * A**(gamma-1)

    null = np.zeros(len(Emean))

    gr = TGraphAsymmErrors(
        len(Emean),
        Emean,
        Flux,
        null, null,
        Flux_err,
        Flux_err
    )

    gr.SetLineColor(color)
    gr.SetMarkerColor(color)
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(size)

    return gr


def draw_panel(pad, graph, title, Ymin, Ymax, show_x=False):

    pad.cd()

    frame = pad.DrawFrame(1e2, Ymin, 1e6, Ymax)  # ora Ymin/Ymax variabili

    # Y axis: label solo sul primo pannello
    frame.GetYaxis().SetTitle("")
    frame.GetYaxis().SetLabelSize(0.045)
    frame.GetYaxis().SetTitleOffset(0.7)

    # X axis: label solo in basso
    frame.GetXaxis().SetTitle("")
    frame.GetXaxis().SetLabelSize(0.045 if show_x else 0)

    graph.Draw("P SAME")

    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.07)
    latex.DrawLatex(0.2, 0.85, title)



if __name__ == '__main__':

    #file_Si_DAMPE = './DATA_POINTS/flux_spectrum_silicon_21bins_nobkg_norm34_5iter_toy10000.dat'
    file_S_DAMPE = 'data/flux_spectrum_sulfur_120months_6-4bin_smooth.dat'  #flux_spectrum_sulfur_120months_6-4bin_nobkg_smooth.dat'
    gr_S_DAMPE = make_flux_graph_DAMPE(file_S_DAMPE, kRed+1, 20, 1.3, n_drop_last=4)

    file_S_AMS = '/mnt/c/Users/saraf/Desktop/dampe/sulfur/AMS/sulfur_ams.txt'
    gr_S_AMS = make_flux_graph_AMS(file_S_AMS, kGray+3, 22, 1.5, 16, 26.06)

    file_S_CRISIS = '/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/SULFUR_ANALYSIS/data/other_exp/S_CRISIS_1977.txt'
    gr_S_CRISIS = make_flux_graph_CRISIS(file_S_CRISIS, kBlue-8, 27, 1.8, 32.)

    file_S_TRACER99 = '/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/SULFUR_ANALYSIS/data/other_exp/S_TRACER_1999.txt'
    gr_S_TRACER99 = make_flux_graph_TRACER99(file_S_TRACER99, kRed-8, 29, 1.4, 28.)

    file_S_TRACER2008 = '/mnt/c/Users/saraf/Desktop/DAMPE_ANALYSIS/SULFUR_ANALYSIS/data/other_exp/S_TRACER_2008.txt'
    gr_S_TRACER2008 = make_flux_graph_TRACER(file_S_TRACER2008, kGray+1, 25, 1.4, 28.)

    cc = TCanvas("cc", "Flux", 900, 750)
    cc.SetLeftMargin(0.13)
    cc.SetRightMargin(0.04)
    cc.SetTopMargin(0.05)
    cc.SetBottomMargin(0.12)
    cc.SetTicks(1,1)
    cc.SetLogx()

    frame = cc.DrawFrame(1e1, 0, 1e6, 1500)

    frame.GetXaxis().SetTitle("Kinetic energy [GeV]")
    frame.GetYaxis().SetTitle("E^{2.7} Flux [m^{-2} s^{-1} sr^{-1} (GeV)^{1.7}]")

    frame.GetXaxis().SetLabelSize(0.035)
    frame.GetYaxis().SetLabelSize(0.035)

    frame.GetXaxis().SetTitleSize(0.035)
    frame.GetYaxis().SetTitleSize(0.035)

    frame.GetXaxis().SetTitleOffset(1.4)
    frame.GetYaxis().SetTitleOffset(1.8)

    frame.GetXaxis().CenterTitle()
    frame.GetYaxis().CenterTitle()

    #gr_S_TRACER99.Draw("P SAME")
    gr_S_TRACER2008.Draw("P SAME")
    gr_S_CRISIS.Draw("P SAME")
    gr_S_AMS.Draw("P SAME")
    gr_S_DAMPE.Draw("P SAME")

    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextFont(42)
    latex.DrawLatex(0.2, 0.86, "Sulfur")

    # ------------------- LEGEND

    leg = TLegend(0.19, 0.48, 0.42, 0.84)  # x1,y1,x2,y2 in NDC pad1
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.025)

    #leg.AddEntry(gr_S_TRACER99,"TRACER 99","pl")
    leg.AddEntry(gr_S_TRACER2008,"TRACER 2008","pl")
    leg.AddEntry(gr_S_CRISIS, "CRISIS [1977]", "PL")
    leg.AddEntry(gr_S_AMS,   "AMS-02 [2023]", "PL")
    leg.AddEntry(gr_S_DAMPE, "DAMPE (this work)", "PL")
    leg.Draw()

    Z = 16
    # oxygen softening
    #E_soft  = 15.4e3 * Z   # GeV
    #dE_soft = 4.8e3 * Z    # GeV
    # general softening
    E_soft  = 15.e3 * Z   # GeV
    dE_soft = 1.6e3 * Z    # GeV
    Ymin = 400
    Ymax = 3400

    #soft = TBox(E_soft - dE_soft, Ymin, E_soft + dE_soft, Ymax)
    #soft.SetFillColorAlpha(kAzure-9, 0.2)
    #soft.SetFillStyle(3002)
    #soft.SetLineColor(0)
    #soft.Draw("same")
    #line_15TV = TLine(E_soft , Ymin, E_soft , Ymax)
    #line_15TV.SetLineColor(12)
    #line_15TV.SetLineWidth(2)
    #line_15TV.SetLineStyle(9)
    #line_15TV.Draw("same")

    title = TLatex()
    title.SetTextAlign(13)
    title.SetTextSize(0.06)
    title.SetTextFont(42)
    title.DrawLatex(4000,1400,"Preliminary")


    cc.Update()

    #cc.SaveAs('PLOTS/flux_SULFUR_noDAMPE_v2.png')
    cc.SaveAs('PLOTS/flux_SULFUR_2026_Orb120Month_MLions_20bins_smooth_PLOT.eps')
    #cc.SaveAs('PLOTS/flux_SULFUR_2026_Orb120Month_MLions_20bins_5iterMax_smooth_PLOT_noHad_d.pdf')

    input("Press enter..")