#!/usr/bin/env python


from argparse import ArgumentParser
import os
import re

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal cmd-line options
r.gROOT.SetBatch(True)
r.gStyle.SetOptStat(False)
import sys
#sys.dont_write_bytecode = True

import math

r.TEventList.__init__._creates = False
r.TH1F.__init__._creates = False
r.TGraphErrors.__init__._creates = False
r.TGraphAsymmErrors.__init__._creates = False


import tools.plot_utils as pu
import tools.utils as utils
import tools.signal as signal
import tools.samples as sample
import tools.region as region
import tools.plot as plot


################################################################################
## Unused functions
#def getSystHists(plot, reg, b, nom_yield, nom_hist) :
#    for s in b.systList :
#        hist_name = ""
#        if "abs" in plot.variable and "DPB_vSS" not in plot.variable :
#            replace_var = plot.variable.replace("abs(","")
#            replace_var = replace_var.replace(")","")
#            hist_name = replace_var
#        elif plot.variable == "DPB_vSS - 0.9*abs(cosThetaB)" :
#            hist_name = "DPB_minus_COSB"
#        else : hist_name = plot.variable
#        h_up = pu.th1d("h_"+b.treename+"_"+hist_name+"_"+s.name+"_up", "", int(plot.nbins), plot.x_range_min, plot.x_range_max, plot.x_label, plot.y_label)
#        h_dn = pu.th1d("h_"+b.treename+"_"+hist_name+"_"+s.name+"_dn", "", int(plot.nbins), plot.x_range_min, plot.x_range_max, plot.x_label, plot.y_label)
#
#        h_up.SetMinimum(plot.y_range_min)
#        h_up.SetMaximum(plot.y_range_max)
#        h_up.GetXaxis().SetTitle(plot.x_label)
#        h_up.GetXaxis().SetTitleFont(42)
#        h_up.GetXaxis().SetLabelFont(42)
#        h_up.GetXaxis().SetLabelSize(0.035)
#        h_up.GetXaxis().SetTitleSize(0.048 * 0.85)
#        h_up.GetXaxis().SetLabelOffset(-999)
#        h_up.GetXaxis().SetTitleOffset(-999)
#
#        h_up.GetYaxis().SetTitle(plot.y_label)
#        h_up.GetYaxis().SetTitleFont(42)
#        h_up.GetYaxis().SetLabelFont(42)
#        h_up.GetYaxis().SetTitleOffset(1.4)
#        h_up.GetYaxis().SetLabelOffset(0.013)
#        h_up.GetYaxis().SetLabelSize(1.2 * 0.035)
#        h_up.GetYaxis().SetTitleSize(0.055 * 0.85)
#
#        h_dn.SetMinimum(plot.y_range_min)
#        h_dn.SetMaximum(plot.y_range_max)
#        h_dn.GetXaxis().SetTitle(plot.x_label)
#        h_dn.GetXaxis().SetTitleFont(42)
#        h_dn.GetXaxis().SetLabelFont(42)
#        h_dn.GetXaxis().SetLabelSize(0.035)
#        h_dn.GetXaxis().SetTitleSize(0.048 * 0.85)
#        h_dn.GetXaxis().SetLabelOffset(-999)
#        h_dn.GetXaxis().SetTitleOffset(-999)
#
#        h_dn.GetYaxis().SetTitle(plot.y_label)
#        h_dn.GetYaxis().SetTitleFont(42)
#        h_dn.GetYaxis().SetLabelFont(42)
#        h_dn.GetYaxis().SetTitleOffset(1.4)
#        h_dn.GetYaxis().SetLabelOffset(0.013)
#        h_dn.GetYaxis().SetLabelSize(1.2 * 0.035)
#        h_dn.GetYaxis().SetTitleSize(0.055 * 0.85)
#
#        #for hsys in [h_up, h_dn] :
#        #    yax = hsys.GetYaxis()
#        #    xax = hsys.GetXaxis()
#
#        #    yax.SetTitleSize(0.05)
#        #    yax.SetLabelSize(0.045)
#        #    yax.SetLabelOffset(0.008)
#        #    yax.SetTitleOffset(1.2)
#        #    yax.SetLabelFont(42)
#        #    yax.SetTitleFont(42)
#        #    yax.SetNdivisions(5)
#
#        if s.isWeightSys() :
#            name_up = s.up_name
#            name_up = "syst_" + name_up.replace('syst_', "")
#            name_dn = s.down_name
#            name_dn = "syst_" + name_dn.replace('syst_', "")
#            weight_up = ""
#            weight_dn = ""
#            if "PILEUPUP" in name_up :
#                weight_up = " eventweightNOPUPW * pupw_up "
#            else :
#                #print " +++ getSystHists isWeightSys UP not applying PRW +++ "
#                #weight_up = " eventweightNOPUPW * %s"%(str(name_dn))
#                weight_up = " eventweight * %s"%(str(name_up))
#
#                if "vv" in b.name and "SF" in reg.name :
#                    print "SYS add mu_VVSF to VV"
#                    weight_up = weight_up + " * 1.02"
#                elif "vv" in b.name and "SF" not in reg.name :
#                    print "SYS add mu_VVDF to VV"
#                    weight_up = weight_up + " * 1.02"
#                elif "ttbar" in b.name :
#                    print "SYS add mu_TTBAR to %s"%b.name
#                    weight_up = weight_up + " * 1.06"
#
#                #if "ttbar" in b.name :
#                #    weight_up = weight_up + " * 0.99"
#                #elif "vv" in b.name :
#                #    if "sf" in reg.name.lower() :
#                #        weight_up = weight_up + " * 1.23"
#                #    else :
#                #        weight_up = weight_up + " * 1.27"
#
#            if "PILEUPDOWN" in name_dn :
#                weight_dn = " eventweightNOPUPW * pupw_down "
#            else :
#                #print " +++ getSystHists isWeightSys DOWN not applying PRW +++ "
#                #weight_dn = " eventweightNOPUPW * %s"%(str(name_dn))
#                weight_dn = " eventweight * %s"%(str(name_dn))
#
#                if "vv" in b.name and "SF" in reg.name :
#                    print "SYS add mu_VVSF to VV"
#                    weight_dn = weight_dn + " * 1.02"
#                elif "vv" in b.name and "SF" not in reg.name :
#                    print "SYS add mu_VVDF to VV"
#                    weight_dn = weight_dn + " * 1.02"
#                elif "ttbar" in b.name :
#                    print "SYS add mu_TTBAR to %s"%b.name
#                    weight_dn = weight_dn + " * 1.06"
#
#                #if "ttbar" in b.name :
#                #    weight_dn = weight_dn + " * 0.99"
#                #elif "vv" in b.name :
#                #    if "sf" in reg.name.lower() :
#                #        weight_dn = weight_dn + " * 1.23"
#                #    else :
#                #        weight_dn = weight_dn + " * 1.27"
#
#            cut_up = "(" + reg.tcut + ") * %s * %s"%(weight_up, str(b.scale_factor))
#            cut_dn = "(" + reg.tcut + ") * %s * %s"%(weight_dn, str(b.scale_factor))
#            #cut_up = "(" + reg.tcut + ") * eventweight * " + str(name_up) + " * " + str(b.scale_factor)
#            #cut_dn = "(" + reg.tcut + ") * eventweight * " + str(name_dn) + " * " + str(b.scale_factor)
#
#            cut_up = r.TCut(cut_up)
#            cut_dn = r.TCut(cut_dn)
#            sel = r.TCut("1")
#
#            cmd_up = "%s>>%s"%(plot.variable, h_up.GetName())
#            cmd_dn = "%s>>%s"%(plot.variable, h_dn.GetName())
#
#            s.tree.Draw(cmd_up, cut_up * sel)
#            s.tree.Draw(cmd_dn, cut_dn * sel)
#
#            # add overflow to these guys' last bins
#            pu.add_overflow_to_lastbin(h_up)
#            pu.add_overflow_to_lastbin(h_dn)
#
#            print "    %s   (+%.2f, -%.2f)"%(s.name, h_up.Integral(0,-1)-nom_yield, nom_yield-h_dn.Integral(0,-1))
#
#            s.up_histo = h_up
#            s.down_histo = h_dn
#
#        elif s.isKinSys() :
#            #print " +++ getSystHists isKinSys not applying PRW +++ "
#            #cut = "(" + reg.tcut + ") * eventweightNOPUPW * " + str(b.scale_factor)
#            cut = "(" + reg.tcut + ") * eventweight * " + str(b.scale_factor)
#
#            if "vv" in b.name and "SF" in reg.name :
#                print "SYS add mu_VVSF to VV"
#                cut = cut + " * 1.02"
#            elif "vv" in b.name and "SF" not in reg.name :
#                print "SYS add mu_VVDF to VV"
#                cut = cut + " * 1.02"
#            elif "ttbar" in b.name :
#                cut = cut + " * 1.06"
#            #if "ttbar" in b.name :
#            #    cut = cut + " * 0.99"
#            #elif "vv" in b.name :
#            #    if "sf" in reg.name.lower() :
#            #        cut = cut + " * 1.23"
#            #    else :
#            #        cut = cut + " * 1.27"
#
#            cut = r.TCut(cut)
#            sel = r.TCut("1")
#            cmd_up = "%s>>%s"%(plot.variable, h_up.GetName())
#            cmd_dn = "%s>>%s"%(plot.variable, h_dn.GetName())
#
#            s.tree_up.Draw(cmd_up, cut * sel)
#            pu.add_overflow_to_lastbin(h_up)
#            s.up_histo = h_up
#
#            is_one_side = False
#            if "JER" in syst.name or ("MET_" in syst.name and "Reso" in syst.name)  : is_one_side = True
#            if "ResoPerp" in syst.name or "ResoPara" in syst.name : is_one_side = True
#            if not s.isOneSided() :
#            #if not is_one_side :
#                s.tree_down.Draw(cmd_dn, cut * sel)
#                pu.add_overflow_to_lastbin(h_dn)
#                s.down_histo = h_dn
#            else :
#                s.down_histo = nom_hist.Clone("%s_down_hist"%s.name)
#                h_dn = s.down_histo
#
#            if s.isOneSided() :
#                print "    %s (+%.2f, -%.2f)"%(s.name, h_up.Integral(0,-1)-nom_yield, nom_yield-h_dn.Integral(0,-1))
#            else :
#                print "    %s  (+%.2f, -%.2f)"%(s.name, h_up.Integral(0,-1)-nom_yield, nom_yield-h_dn.Integral(0,-1))
#
#def check_2d_consistency(plot, data, backgrounds) :
#    if plot.sample == "Data" and data.treename == "":
#        print 'check_2d_consistency ERROR    Requested sample is ("Data") and the data sample is empty. Exitting.'
#        sys.exit()
#    sample_in_backgrounds = False
#    for b in backgrounds :
#        if plot.sample == b.name and plot.sample != "Data" : sample_in_backgrounds = True
#    if not sample_in_backgrounds and plot.sample != "Data" :
#        print 'check_2d_consistency ERROR    Requested sample ("%s") is not in the backgrounds list:'%plot.sample
#        print backgrounds
#        print 'check_2d_consistency ERROR    Exitting.'
#        sys.exit()
#
#    print 'check_2d_consistency    Samples consistent with plot.'

################################################################################
def histos_for_legend(histos) :
    '''
    rearrange histogram list for legend

    param:
        histos : list(TH1F)
            histograms in plot ordered by total events (assumes len >= 4)

    returns:
        list(TH1F)
    '''

    if len(histos) == 7 :
        indices = [0, 4, 1, 5, 2, 6, 3]
    elif len(histos) == 6 :
        indices = [0, 3, 1, 4, 2, 5]
    elif len(histos) == 5 :
        indices = [0, 3, 1, 4, 2]
    elif len(histos) == 4 :
        indices = [0, 2, 1, 3]
    else:
        return histos

    return [histos[idx] for idx in indices]

################################################################################
def make_1dprofileRMS(plot, reg, data, backgrounds ) :
    print "make_1dprofileRMS    Plotting %s"%plot.name

    r.gStyle.SetOptStat(1100)

    # grab the plot's canvas and pretty up
    c = plot.canvas
    c.cd()
    c.SetGridx(1)
    c.SetGridy(1)
    c.SetFrameFillColor(0)
    c.SetFillColor(0)
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.14)
    c.SetBottomMargin(1.3*c.GetBottomMargin())

    name_on_plot = ""
    if plot.sample != "Data" :
        for b in backgrounds :
            if b.name != plot.sample : continue
            name_on_plot += b.displayname

            hist_name_x = ""
            hist_name_y = ""
            if "abs" in plot.xVariable :
                x_repl = plot.xVariable.replace("abs(","")
                x_repl = x_repl.replace(")","")
                hist_name_x = x_repl
            else : hist_name_x = plot.xVariable

            if "abs" in plot.yVariable :
                y_repl = plot.yVariable.replace("abs(","")
                y_repl = y_repl.replace(")","")
                hist_name_y = y_repl
            else : hist_name_y = plot.yVariable

            hx = pu.th1d("h_"+b.treename+"_"+hist_name_x, "", int(plot.n_binsX), plot.x_range_min, plot.x_range_max, plot.x_label, plot.y_label)
            hx.Sumw2
            cut = "(" + reg.tcut + ") * eventweight * " + str(b.scale_factor)
            cut = r.TCut(cut)
            sel = r.TCut("1")
            cmd = "%s>>+%s"%(plot.xVariable, hx.GetName())
            b.tree.Draw(cmd, cut * sel)

            gr = r.TGraphErrors()

            for i in range(hx.GetNbinsX()) :
                hy = pu.th1d("h_" + b.treename + "_" + hist_name_y, "", 100, -200, 200,  plot.y_label, "")
                cut_up = hx.GetBinLowEdge(i+1) + hx.GetBinWidth(i+1)
                cut_down = hx.GetBinLowEdge(i+1)
                cut = "(" + reg.tcut + " && ( %s >= %s && %s <= %s)"%(plot.xVariable, cut_down, plot.xVariable, cut_up)  + ") * eventweight * " + str(b.scale_factor)
                cut = r.TCut(cut)
                sel = r.TCut("1")
                cmd = "%s>>%s"%(plot.yVariable, hy.GetName())
                b.tree.Draw(cmd, cut * sel)
                rms = hy.GetRMS()
                rms_err = hy.GetRMSError()
                gr.SetPoint(i, hx.GetBinCenter(i+1), rms)
                gr.SetPointError(i, 0.5*hx.GetBinWidth(i+1), rms_err)
                hy.Delete()

            gr.SetMarkerStyle(8)
            gr.SetMarkerSize(1.15 * gr.GetMarkerSize())
            gr.SetMarkerColor(r.TColor.GetColor("#E67067"))
            gr.GetYaxis().SetRangeUser(plot.y_range_min, plot.y_range_max)

            gr.Draw("ap")
            r.gPad.Update()
            gr.Draw("ap")
            gr.GetYaxis().SetTitle(plot.y_label)
            gr.GetXaxis().SetTitle(plot.x_label)


    if plot.sample == "Data" :
        name_on_plot = "Data"
        hist_name_x = ""
        hist_name_y = ""
        if "abs" in plot.xVariable :
            x_repl = plot.xVariable.replace("abs(","")
            x_repl = x_repl.replace(")","")
            hist_name_x = x_repl
        else : hist_name_x = plot.xVariable

        if "abs" in plot.yVariable :
            y_repl = plot.yVariable.replace("abs(","")
            y_repl = y_repl.replace(")","")
            hist_name_y = y_repl
        else : hist_name_y = plot.yVariable

        hx = pu.th1d("h_data_"+hist_name_x, "", int(plot.n_binsX), plot.x_range_min, plot.x_range_max, plot.x_label, plot.y_label)
        hx.Sumw2
        cut = "(" + reg.tcut + ") * eventweight"
        cut = r.TCut(cut)
        sel = r.TCut("1")
        cmd = "%s>>+%s"%(plot.xVariable, hx.GetName())
        data.tree.Draw(cmd, cut * sel)

        g = r.TGraphErrors()

        for i in range(hx.GetNbinsX()) :
            hy = pu.th1d("h_data_" + hist_name_y, "", 100, -200, 200,  plot.y_label, "")
            cut_up = hx.GetBinLowEdge(i+1) + hx.GetBinWidth(i+1)
            cut_down = hx.GetBinLowEdge(i+1)
            cut = "(" + reg.tcut + " && ( %s >= %s && %s <= %s)"%(plot.xVariable, cut_down, plot.xVariable, cut_up)  + ") * eventweight"
            cut = r.TCut(cut)
            sel = r.TCut("1")
            cmd = "%s>>%s"%(plot.yVariable, hy.GetName())
            data.tree.Draw(cmd, cut * sel)
            rms = hy.GetRMS()
            rms_err = hy.GetRMSError()
            gr.SetPoint(i, hx.GetBinCenter(i+1), rms)
            gr.SetPointError(i, 0.5*hx.GetBinWidth(i+1), rms_err)
            hy.Delete()

        gr.SetMarkerStyle(8)
        gr.SetMarkerSize(1.15 * gr.GetMarkerSize())
        gr.SetMarkerColor(r.TColor.GetColor("#5E9AD6"))
        gr.GetYaxis().SetRangeUser(plot.y_range_min, plot.y_range_max)

        gr.Draw("ap")
        r.gPad.Update()
        gr.Draw("ap")
        gr.GetYaxis().SetTitle(plot.y_label)
        gr.GetXaxis().SetTitle(plot.x_label)

    pu.draw_text_on_top(text="%s : #bf{%s}"%(plot.name, name_on_plot))
    c.Update()
    r.gPad.RedrawAxis()

    # set output
    outname = plot.name + ".eps"
    out = g.plots_dir + g_outdir
    c.SaveAs(outname)
    utils.mv_file_to_dir(outname, out, True)
    fullname = out + "/" + outname
    print "%s saved to : %s"%(outname, os.path.abspath(fullname))

def make_1dprofile(plot, reg, data, backgrounds) :

    print "make_1dprofile    Plotting %s"%plot.name

    # get the stats box back to show the mean and err
    r.gStyle.SetOptStat(1100)

    # grab the plot's canvas and pretty up
    c = plot.canvas
    c.cd()
    c.SetGridx(1)
    c.SetGridy(1)
    c.SetFrameFillColor(0)
    c.SetFillColor(0)
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.14)
    c.SetBottomMargin(1.3*c.GetBottomMargin())

    # name_on_plot : text to put on top of the plot pad in addition to the region name so that
    # we know which sample is being plotted on the profile
    name_on_plot = ""

    if plot.sample != "Data" :
        for b in backgrounds :
            if b.name != plot.sample : continue
            name_on_plot += b.displayname
            h = r.TProfile("hprof_"+b.name+"_"+plot.xVariable+"_"+plot.yVariable,";%s;%s"%(plot.x_label,plot.y_label), int(plot.n_binsX), plot.x_range_min, plot.x_range_max, plot.y_range_min, plot.y_range_max)

            cut = "(" + reg.tcut + ") * eventweight * " + str(b.scale_factor)
            cut = r.TCut(cut)
            sel = r.TCut("1")
            cmd = "%s:%s>>+%s"%(plot.yVariable, plot.xVariable, h.GetName())
            b.tree.Draw(cmd, cut * sel, "prof")
            h.SetMarkerColor(r.TColor.GetColor("#E67067"))

    if plot.sample == "Data" :
        name_on_plot = "Data"
        h = r.TProfile("hprof_Data_"+plot.xVariable+"_"+plot.yVariable,";%s;%s"%(plot.x_label,plot.y_label), int(plot.n_binsX), plot.x_range_min, plot.x_range_max, plot.y_range_min, plot.y_range_max)

        cut = "(" + reg.tcut + ")"
        cut = r.TCut(cut)
        sel = r.TCut("1")
        cmd = "%s:%s>>+%s"%(plot.yVariable, plot.xVariable, h.GetName())
        data.tree.Draw(cmd, cut * sel, "prof")
        h.SetMarkerColor(r.TColor.GetColor("#5E9AD6"))

    #  r.TGaxis.SetMaxDigits(2)
    h.GetYaxis().SetTitleOffset(1.6 * h.GetYaxis().GetTitleOffset())
    h.GetXaxis().SetTitleOffset(1.2 * h.GetXaxis().GetTitleOffset())
    h.SetMarkerStyle(8)
    h.SetLineColor(r.kBlack)
    h.SetMarkerSize(1.15*h.GetMarkerSize())
    h.GetYaxis().SetRangeUser(plot.y_range_min,plot.y_range_max)

    h.Draw()
    r.gPad.Update()

    # move the stats box so it does not block the text on top
    st = h.FindObject("stats")
    st.SetY1NDC(0.93 * st.GetY1NDC())
    st.SetY2NDC(0.93 * st.GetY2NDC())

    # draw
    h.Draw()

    # draw descriptive text on top
  #  pu.draw_text_on_top(text="%s : #bf{%s}"%(plot.name, name_on_plot),pushup=1.035)
    pu.draw_text_on_top(text="%s : #bf{%s}"%(plot.name, name_on_plot))
    c.Update()

    r.gPad.RedrawAxis()

    # set output
    outname =  plot.name + ".eps"
    out = g.plots_dir + g_outdir
    c.SaveAs(outname)
    utils.mv_file_to_dir(outname, out, True)
    fullname = out + "/" + outname
    print "%s saved to : %s"%(outname, os.path.abspath(fullname))

def make_plotsRatio(plot, reg, data, backgrounds, plot_i, n_plots) :
    '''
    Make a ratio plot with background stack, signal, and data on top
    canvas and a ratio of data to SM background on the bottom

    param:
        plot : plot class
            1D plot defined in configuration file
        reg : region class
            region defined in configuration file
        data : background class
            data sample defined in configuration file
        backgrounds : list(background class)
            background samples defined in configuration file
        plot_i : int
            current plot number
        n_plots : int
            total number of plots

    returns (None)
    '''


    # add_scale_factors = False
    # if not add_scale_factors :
    #     print " *** NOT ADDING S2L SCALE FACTORS *** "

    print 20*"-","Plotting [%d/%d] %s"%(plot_i, n_plots, plot.name), 20*'-'

    ############################################################################
    # Intialize plot components

    # Canvases
    rcan = plot.ratioCanvas
    rcan.canvas.cd()
    rcan.upper_pad.cd()

    if plot.isLog() : rcan.upper_pad.SetLogy(True)
    rcan.upper_pad.Update()

    # Stack for MC
    hax = r.TH1F("axes", "", int(plot.nbins), plot.x_range_min, plot.x_range_max)
    hax.SetMinimum(plot.y_range_min)
    hax.SetMaximum(plot.y_range_max)
    hax.GetXaxis().SetTitle(plot.x_label)
    hax.GetXaxis().SetTitleFont(42)
    hax.GetXaxis().SetLabelFont(42)
    hax.GetXaxis().SetLabelSize(0.035)
    hax.GetXaxis().SetTitleSize(0.048 * 0.85)
    # ADD: Move x_axis out of way for ratio plot
    hax.GetXaxis().SetTitleOffset(-999)
    hax.GetXaxis().SetLabelOffset(-999)

    hax.GetYaxis().SetTitle(plot.y_label)
    hax.GetYaxis().SetTitleFont(42)
    hax.GetYaxis().SetLabelFont(42)
    hax.GetYaxis().SetTitleOffset(1.4)
    hax.GetYaxis().SetLabelOffset(0.013)
    hax.GetYaxis().SetLabelSize(1.2 * 0.035)
    hax.GetYaxis().SetTitleSize(0.055 * 0.85)
    hax.Draw()
    rcan.upper_pad.Update()

    stack = r.THStack("stack_"+plot.name, "")

    # Legend
    # ADD: change legend location
    leg = pu.default_legend(xl=0.55,yl=0.71,xh=0.93,yh=0.90)
    leg.SetNColumns(2)

    # Yield table
    yield_tbl = []

    ############################################################################
    # Get background MC stack

    # Initilize lists and defaults
    histos = []
    all_histos = []
    avoid_bkg = []
    #h_nom_fake = None

    n_total_sm_yield = n_total_sm_error = 0.
    hists_to_clear = []

    # Add MC backgrounds to stack
    for b in backgrounds :
        # Initilize histogram
        h_name = "h_"+reg.name+'_'+b.treename+"_"+plot.variable
        hists_to_clear.append(h_name)
        h = pu.th1d(h_name, "", int(plot.nbins), plot.x_range_min, plot.x_range_max, plot.x_label, plot.y_label)

        h.SetLineColor(b.color)
        h.GetXaxis().SetLabelOffset(-999)
        h.SetFillColor(b.color)
        h.SetFillStyle(1001)
        h.Sumw2


        # cut and make the sample weighted, applying the scale_factor
        weight_str = "eventweight"
        #weight_str = "1"
        # weight_str = ""
        # if add_scale_factors :
        #     if "fakes" in b.name :
        #         weight_str = "FakeWeight"
        #     elif "vv" in b.name and "SF" in reg.name :
        #         print "adding mu_VVSF"
        #         weight_str = "eventweight * 1.02"
        #     elif "vv" in b.name and "SF" not in reg.name :
        #         print "adding mu_VVDF"
        #         weight_str = "eventweight * 1.02"
        #     elif "ttbar" in b.name :
        #         print "adding mu_TTBAR"
        #         weight_str = "eventweight * 1.06"
        #     else :
        #         weight_str = "eventweight"
        #     #elif "vv" in b.name and "crv" in reg.name and "SF" not in reg.name :
        #     #    weight_str = "eventweight * 1.27"
        #     #elif "vv" in b.name and "crvSF" in reg.name :
        #     #    weight_str = "eventweight * 1.22"
        # else :
        #     #print "+++ not applying pileup reweighting to sample %s +++"%b.name
        #     #weight_str = "eventweightNOPUPW"
        #     weight_str = "eventweight"

        # Draw final histogram (i.e. selections and weights applied)
        cut = "(" + reg.tcut + ") * %s * "%weight_str + str(b.scale_factor)
        cut = r.TCut(cut)
        sel = r.TCut("1")
        draw_cmd = "%s>>+%s"%(plot.variable, h.GetName())
        b.tree.Draw(draw_cmd, cut * sel, "goff")

        # Yield +/- stat error
        stat_err = r.Double(0.0)
        integral = h.IntegralAndError(0,-1,stat_err)
        if not b.isSignal():
            n_total_sm_yield += float(integral)
            n_total_sm_error = (n_total_sm_error**2 + float(stat_err)**2)**0.5
        yield_tbl.append("%10s: %.2f +/- %.2f"%(b.name, integral, stat_err))

        # Add overflow
        if plot.add_overflow:
            pu.add_overflow_to_lastbin(h)

        # Record all and non-empty histograms
        all_histos.append(h)
        histos.append(h) if integral > 0 else avoid_bkg.append(b.name)
        rcan.upper_pad.Update()

    if not len(histos):
        print "make_plotsRatio ERROR :: All SM hists are empty. Skipping"
        return
    #yield_tbl.append("%10s: %.2f"%('Total SM', n_total_sm_yield))
    yield_tbl.append("%10s: %.2f +/- %.2f"%('Total SM', n_total_sm_yield, n_total_sm_error))

    # Order the hists by total events
    histos = sorted(histos, key=lambda h: h.Integral())
    for h in histos :
        #if "fake" in h.GetName() or "Fake" in h.GetName() : continue
        stack.Add(h)
        # add items to legend in order of stack
        # name_for_legend = ""
        # for b in backgrounds :
        #     if b.treename in h.GetName() :
        #         name_for_legend = b.displayname
        # leg.AddEntry(h, name_for_legend, "f")
    rcan.upper_pad.Update()

    # removes = {}
    # removal_list = []
    # try :
    #     removal_list = removes[reg.name]
    # except :
    #     print "S2L NOT REMOVING SAMPLES FROM LEGEND"
    #     pass

    #tmp_leg_histos = []
    #for h in all_histos :
    #    keep_histo = True
    #    for rem in removal_list :
    #        if rem in h.GetName() :
    #            keep_histo = False
    #            print "not adding %s to legend"%rem
    #    if keep_histo :
    #        tmp_leg_histos.append(h)
    # print "length of tmp_leg_histos = %d"%len(tmp_leg_histos)



    h_leg = sorted(all_histos, key=lambda h: h.Integral(), reverse=True)
    histos_for_leg = histos_for_legend(h_leg)

    # # we do this down below now so that data and SM can be on top
    # for h in histos_for_leg :
    #     # add items to legend in order of stack
    #     name_for_legend = ""
    #     for b in backgrounds :
    #         if b.treename in h.GetName() :
    #             name_for_legend = b.displayname
    #     leg.AddEntry(h, name_for_legend, "f")
    # rcan.upper_pad.Update()


    ############################################################################
    # Get data hist

    hd_name = "h_"+reg.name+'_data_'+plot.variable
    hd = pu.th1d(hd_name, "", int(plot.nbins), plot.x_range_min, plot.x_range_max, plot.x_label, plot.y_label)
    hd.Sumw2

    cut = "(" + reg.tcut + ")"
    cut = r.TCut(cut)
    sel = r.TCut("1")
    draw_cmd = "%s>>%s"%(plot.variable, hd.GetName())
    data.tree.Draw(draw_cmd, cut * sel, "goff")
    hd.GetXaxis().SetLabelOffset(-999)

    # print the yield +/- stat error
    stat_err = r.Double(0.0)
    integral = hd.IntegralAndError(0,-1,stat_err)
    #yield_tbl.append("%10s: %.2f +/- %.2f"%('Data',integral, stat_err))
    yield_tbl.append("%10s: %.2f"%('Data',integral))
    yield_tbl.append("%10s: %.2f"%("Data/MC", integral/n_total_sm_yield))
    #print "Data: %.2f +/- %.2f"%(integral, stat_err)

    # Add overflow
    if plot.add_overflow:
        pu.add_overflow_to_lastbin(hd)

    gdata = pu.convert_errors_to_poisson(hd)
    #gdata.SetLineWidth(2)
    #uglify
    gdata.SetLineWidth(1)
    gdata.SetMarkerStyle(20)
    gdata.SetMarkerSize(1.5)
    gdata.SetLineColor(1)
    leg.AddEntry(gdata, "Data", "p")
    rcan.upper_pad.Update()


    ############################################################################
    # systematics loop
    r.gStyle.SetHatchesSpacing(0.9)

    # dummy histo for legend
    #mcError = r.TH1F("mcError", "mcError", 2,0,2)
    #mcError.SetLineWidth(3)
    #mcError.SetFillStyle(3354)
    #mcError.SetFillColor(r.kBlack)
    #mcError.SetLineColor(r.TColor.GetColor("#FC0F1D"))
    #leg.AddEntry(mcError, "Total SM", "fl")
    # uglify for Moriond
    mcError = r.TH1F("mcError", "mcError", 2,0,2)
    #mcError.SetLineWidth(2)
    #uglify
    mcError.SetLineWidth(3)
    mcError.SetFillStyle(3345)
    mcError.SetFillColor(r.kBlue)
    mcError.SetLineColor(r.kBlack)
    #leg.AddEntry(mcError, "Total SM", "fl")
    #uglify for Moriond
    leg.AddEntry(mcError, "Standard Model", "fl")

    # now add backgrounds to legend
    for h in histos_for_leg :
        name_for_legend = ""
        for b in backgrounds :
            if b.treename in h.GetName()[7:] : #TESTING
                name_for_legend = b.displayname
        leg.AddEntry(h, name_for_legend, "f")

    # histogram for total stack
    #totalSM = stack.GetStack().Last().Clone("totalSM")
    #nominalAsymErrors = pu.th1_to_tgraph(totalSM)
    #nominalAsymErrors.SetMarkerSize(0)
    #nominalAsymErrors.SetLineWidth(0)
    #nominalAsymErrors.SetFillStyle(3354)
    #nominalAsymErrors.SetFillColor(r.kGray + 3)
   # leg.AddEntry(nominalAsymErrors, "Bkg. Uncert.", "f")
    # uglify for stop-2l Moriond
    totalSM = stack.GetStack().Last().Clone("totalSM")
    nominalAsymErrors = pu.th1_to_tgraph(totalSM)
    nominalAsymErrors.SetMarkerSize(0)
    nominalAsymErrors.SetLineWidth(0)
    nominalAsymErrors.SetFillStyle(3345)
    nominalAsymErrors.SetFillColor(r.kBlue)

    # if g_doSys :
    #     # totalSystHisto will hold each samples'
    #     # variation
    #     totalSysHisto = totalSM.Clone()
    #     totalSysHisto.Reset("ICE")
    #     transient = r.TGraphAsymmErrors()

    #     # add to the error band the contribution from the up-variations
    #     systematics_up = [s.up_name for s in backgrounds[0].systList]
    #     for up_sys in systematics_up :
    #         transient = r.TGraphAsymmErrors()
    #         transient.Clear()
    #         for b in backgrounds :
    #             if b.name in avoid_bkg : continue
    #             if b.isSignal() : continue
    #             if 'fakes' in b.name : continue
    #             for syst in b.systList :
    #                 if syst.up_name != up_sys : continue
    #                 #print "[%s] adding %s to up histo : %.2f (%.2f)"%(up_sys, b.name, syst.up_histo.Integral(), totalSM.Integral())
    #                 totalSysHisto.Add(syst.up_histo)
    #         print "NOT ADDDING FAKE TO SYS HISTOS"
    #         #if h_nom_fake :
    #         #    totalSysHisto.Add(h_nom_fake)
    #         transient = pu.th1_to_tgraph(totalSysHisto)
    #         print " > %s"%up_sys
    #         pu.add_to_band(transient, nominalAsymErrors)
    #         totalSysHisto.Reset()

    #     # add to the error band the contribution from the down-variations
    #     systematics_down = [s.down_name for s in backgrounds[0].systList]
    #     for dn_sys in systematics_down :
    #         transient = r.TGraphAsymmErrors()
    #         transient.Clear()
    #         for b in backgrounds :
    #             if b.name in avoid_bkg : continue
    #             if b.isSignal() : continue
    #             if 'fakes' in b.name : continue
    #             for syst in b.systList :
    #                 if syst.down_name != dn_sys : continue
    #                 #print "[%s] adding %s to down histo : %.2f"%(dn_sys, b.name, syst.down_histo.Integral())
    #                # if syst.isOneSided() : continue

    #                # if "JER" in syst.name : continue
    #                # if "ResoPerp" in syst.name or "ResoPara" in syst.name : continue
    #                 totalSysHisto.Add(syst.down_histo)
    #         print "NOT ADDDING FAKE TO SYS HISTOS"
    #         #if h_nom_fake :
    #         #    totalSysHisto.Add(h_nom_fake)
    #         transient = pu.th1_to_tgraph(totalSysHisto)
    #         print " > %s"%dn_sys
    #         pu.add_to_band(transient, nominalAsymErrors)
    #         totalSysHisto.Reset()

    ############################################################################
    # draw the MC stack and do cosmetics
    stack.SetMinimum(plot.y_range_min)

    # Determine y range for plot
    max_mult = 2.0 if any([b.isSignal() for b in backgrounds]) else 1.66

    maxy = max(hd.GetMaximum(), stack.GetMaximum())
    if not plot.isLog() :
        hax.SetMaximum(max_mult*maxy)
        hax.Draw()
        rcan.upper_pad.Update()
        stack.SetMaximum(max_mult*maxy)
    else :
        hax.SetMaximum(1e3*plot.y_range_max)
        hax.Draw()
        rcan.upper_pad.Update()
        stack.SetMaximum(1e3*plot.y_range_max)
    #    #stack.SetMaximum(plot.y_range_max)

    # Add stack to canvas
    stack.Draw("HIST SAME")
    rcan.upper_pad.Update()

    # symmetrize the errors
    for i in xrange(nominalAsymErrors.GetN()) :
        ehigh = nominalAsymErrors.GetErrorYhigh(i)
        elow  = nominalAsymErrors.GetErrorYlow(i)


        error_sym = r.Double(0.0)
        error_sym = (ehigh + elow) / 2.0

        if ehigh != error_sym:
            print "initial error (+%.2f,-%.2f), symmetrized = (+%.2f,-%.2f)"%(ehigh,elow, error_sym, error_sym)


        nominalAsymErrors.SetPointEYhigh(i,0.0)
        nominalAsymErrors.SetPointEYhigh(i, error_sym)
        nominalAsymErrors.SetPointEYlow(i,0.0)
        nominalAsymErrors.SetPointEYlow(i,error_sym)

    # draw the error band
    nominalAsymErrors.Draw("same && E2")

    # draw the total bkg line
    hist_sm = stack.GetStack().Last().Clone("hist_sm")
    #hist_sm.SetLineColor(r.TColor.GetColor("#FC0F1D"))
    #uglify for Stp2-l Moriond
    hist_sm.SetLineColor(r.kBlack)
    hist_sm.SetLineWidth(mcError.GetLineWidth())
    hist_sm.SetLineStyle(1)
    hist_sm.SetFillStyle(0)
    #uglify
    hist_sm.SetLineWidth(3)
    hist_sm.Draw("hist same")

    ############################################################################
    # Get signal hist
    leg_sig = pu.default_legend(xl=0.55, yl=0.6, xh=0.91, yh=0.71)
    leg_sig.SetNColumns(1)

    sig_histos = []
    for s in backgrounds :
        if not s.isSignal() : continue

        hist_name = ""
        if "abs" in plot.variable and "DPB" not in plot.variable :
            replace_var = plot.variable.replace("abs(","")
            replace_var = replace_var.replace(")","")
            hist_name = replace_var
        elif plot.variable == "DPB_vSS - 0.9*abs(cosThetaB)" :
            hist_name = "DPB_minus_COSB"
        else : hist_name = plot.variable

        h = pu.th1d("h_"+s.treename+"_"+hist_name, "", int(plot.nbins), plot.x_range_min, plot.x_range_max, plot.x_label, plot.y_label)
        h.SetLineWidth(2)
        h.SetLineStyle(2)
        h.SetLineColor(s.color)
        h.GetXaxis().SetLabelOffset(-999)
        h.SetFillStyle(0)
        h.Sumw2

        #cut = "(" + new_tcut + ") * eventweightNOPUPW * susy3BodyRightPol *" + str(s.scale_factor)
        cut = "(" + reg.tcut + ") * eventweight *" + str(s.scale_factor)
        cut = r.TCut(cut)
        sel = r.TCut("1")
        cmd = "%s>>+%s"%(plot.variable, h.GetName())
        s.tree.Draw(cmd, cut * sel, "goff")

        # print the yield +/- stat error
        stat_err = r.Double(0.0)
        integral = h.IntegralAndError(0,-1,stat_err)
        yield_tbl.append("%10s: %.2f +/- %.2f"%(s.name, integral, stat_err))
        #print "%s: %.2f +/- %.2f"%(s.name, integral, stat_err)

        # add overflow
        pu.add_overflow_to_lastbin(h)

        #leg.AddEntry(h, s.displayname, "l")
        leg_sig.AddEntry(h, s.displayname, "l")
        sig_histos.append(h)
        rcan.upper_pad.Update()

    ############################################################################
    #
    # Print yield table
    print "Yields for %s region"%reg.displayname
    print 30*'-'
    for yld in yield_tbl:
        if 'total' in yld.lower() or 'Data/MC' in yld:
            print 20*'-'
        print yld
    print 30*'-'

    #draw the signals
    for hsig in sig_histos :
        hsig.Draw("hist same")

    # draw the data graph
    gdata.Draw("option same pz 0")
    # draw the legend
    leg.Draw()
    leg_sig.Draw()
    r.gPad.RedrawAxis()

    # add some text/labels
    #uglify
    pu.draw_text(text="ATLAS",x=0.18,y=0.85,size=0.05,font=72)
    pu.draw_text(text="Internal",x=0.325,y=0.85,size=0.05,font=42)
    #pu.draw_text(text="Internal",x=0.325,y=0.85,size=0.06,font=42)
    pu.draw_text(text="#sqrt{s} = 13 TeV, 36.1 fb^{-1}", x=0.18, y=0.79, size=0.04)
    #pu.draw_text(text="L = 36 fb^{-1}, #sqrt{s} = 13 TeV",x=0.18,y=0.79, size=0.04)
    pu.draw_text(text="Higgs LFV", x=0.18, y=0.74, size=0.04)
    pu.draw_text(text=reg.displayname,      x=0.18,y=0.68, size=0.04)
    #pu.draw_text(text=reg.displayname,      x=0.18,y=0.74, size=0.04)

    r.gPad.SetTickx()
    r.gPad.SetTicky()

    rcan.upper_pad.Update()

    #### Lower Pad
    rcan.lower_pad.cd()

    rcan.lower_pad.SetTickx()
    rcan.lower_pad.SetTicky()
    rcan.lower_pad.Update()

    # get the total SM histo
    h_sm = stack.GetStack().Last().Clone("h_sm")

    # yaxis
    yax = h_sm.GetYaxis()
    yax.SetRangeUser(0,2)
    yax.SetTitle("Data / SM")
    yax.SetTitleSize(0.14 * 0.83)
    yax.SetLabelSize(0.13 * 0.81)
    yax.SetLabelOffset(0.98 * 0.013 * 1.08)
    yax.SetTitleOffset(0.45 * 1.2)
    yax.SetLabelFont(42)
    yax.SetTitleFont(42)
    yax.SetNdivisions(5)
    # xaxis
    xax = h_sm.GetXaxis()
    xax.SetTitleSize(1.1 * 0.14)
    xax.SetLabelSize(yax.GetLabelSize())
    xax.SetLabelOffset(1.15*0.02)
    xax.SetTitleOffset(0.85 * xax.GetTitleOffset())
    xax.SetLabelFont(42)
    xax.SetTitleFont(42)

    h_sm.SetTickLength(0.06)
    h_sm.Draw("AXIS")
    rcan.lower_pad.Update()

    # get the ratio-error band
    ratioBand = r.TGraphAsymmErrors(nominalAsymErrors)
    pu.buildRatioErrorBand(nominalAsymErrors, ratioBand)

    # draw lines
    #pu.draw_line(plot.x_range_min, 1.0, plot.x_range_max, 1.0,color=r.kRed,style=2,width=1)
    #uglify
    pu.draw_line(plot.x_range_min, 1.0, plot.x_range_max, 1.0,color=r.kBlack,style=2,width=1)
    pu.draw_line(plot.x_range_min, 0.5, plot.x_range_max, 0.5,style=3,width=1)
    pu.draw_line(plot.x_range_min, 1.5, plot.x_range_max, 1.5,style=3,width=1)

    # convert to tgraphs to get the ratio
    g_data = pu.convert_errors_to_poisson(hd)
    g_sm = pu.th1_to_tgraph(h_sm)
    g_ratio = pu.tgraphAsymmErrors_divide(g_data, g_sm)

    # For Data/MC only use the statistical error for data
    # since we explicity draw the MC error band
    nominalAsymErrorsNoSys = r.TGraphAsymmErrors(nominalAsymErrors)
    for i in xrange(nominalAsymErrorsNoSys.GetN()) :
        nominalAsymErrorsNoSys.SetPointError(i-1,0,0,0,0)
    ratio_raw = pu.tgraphAsymmErrors_divide(g_data, nominalAsymErrorsNoSys)
    #ratio_raw = pu.tgraphAsymmErrors_divide(nominalAsymErrorsNoSys,g_data)
    ratio = r.TGraphAsymmErrors()

    x1, y1 = r.Double(0.0), r.Double(0.0)
    index = 0
    for i in xrange(ratio_raw.GetN()) :
        ratio_raw.GetPoint(i, x1, y1)
        if y1 > 0. :
            ratio.SetPoint(index, x1, y1)
            ratio.SetPointError(index, ratio_raw.GetErrorXlow(i), ratio_raw.GetErrorXhigh(i), ratio_raw.GetErrorYlow(i), ratio_raw.GetErrorYhigh(i))
            index+=1
    #ratio.SetLineWidth(2)
    #uglify
    ratio.SetLineWidth(1)
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(1.5)
    ratio.SetLineColor(1)
    ratio.SetMarkerSize(1.5)
    rcan.lower_pad.Update()

    ratioBand.Draw("E2")
    #ratioBand.Draw("same && E2")
    ratio.Draw("option same pz 0")
    rcan.lower_pad.Update()

    rcan.canvas.Update()
    ############################################################################

    print "\n"
    outname = plot.name + ".eps"
    rcan.canvas.SaveAs(outname)
    out = g.plots_dir + g_outdir
    utils.mv_file_to_dir(outname, out, True)
    fullname = out + "/" + outname

    # outname = plot.name + ".root"
    # rcan.canvas.SaveAs(outname)
    # out = g.plots_dir + g_outdir
    # utils.mv_file_to_dir(outname, out, True)

    print "%s saved to : %s"%(outname, os.path.abspath(fullname))
    ############################################################################

def make_plotsComparison(plot, reg, data, backgrounds) :
    print "make_plotsComparison    Plotting %s"%plot.name

    c = plot.canvas
    c.cd()
    c.SetFrameFillColor(0)
    c.SetFillColor(0)
    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(1.3*c.GetBottomMargin())

    if plot.isLog() : c.SetLogy(True)
    c.Update()

    leg = None
    if plot.leg_is_left : leg = pu.default_legend(xl=0.2,yl=0.7,xh=0.47, yh=0.87)
    elif plot.leg_is_bottom_right : leg = pu.default_legend(xl=0.7, yl=0.17,xh=0.97,yh=0.41)
    elif plot.leg_is_bottom_left : leg = pu.default_legend(xl=0.2,yl=0.2,xh=0.47,yh=0.37)
    else : leg = pu.default_legend(xl=0.7,yl=0.65,xh=0.97,yh=0.87)

    histos = []
    maxy = []

    for b in backgrounds :
        hist_name = ""

        if "abs" in plot.variable and "DPB_vSS" not in plot.variable :
            replace_var = plot.variable.replace("abs(","")
            replace_var = replace_var.replace(")","")
            hist_name = replace_var
        elif plot.variable == "DPB_vSS - 0.9*abs(cosThetaB)" :
            hist_name = "DPB_minus_COSB"
        else : hist_name = plot.variable

        h = pu.th1d("h_"+b.treename+"_"+hist_name, "", int(plot.nbins), plot.x_range_min, plot.x_range_max, plot.x_label, plot.y_label)
        h.SetLineColor(b.color)
        h.SetLineWidth(2)
        h.SetLineStyle(b.line_style)
        h.SetFillColor(0)
        h.Sumw2

        # cut and make the sample weighted, applying any scale_factor
        weight_str = ""
        if b.isSignal() :
            #weight_str = "eventweightNOPUPW * susy3BodyRightPol"
            weight_str = "eventweight"
        else :
            #weight_str = "1"
            weight_str = "eventweight"
        cut = "(" + reg.tcut + ") * %s *"%weight_str + str(b.scale_factor)
        cut = r.TCut(cut)
        sel = r.TCut("1")
        cmd = "%s>>+%s"%(plot.variable, h.GetName())
        b.tree.Draw(cmd, cut * sel)

        # print the yield +/- stat error
        stat_err = r.Double(0.0)
        integral = h.IntegralAndError(0,-1, stat_err)
        print "%s: %.2f +/- %.2f"%(b.name, integral, stat_err)

        # normalize since we only care about shapes
        if integral != 0 : h.Scale(1/integral)

        # add overflow
        pu.add_overflow_to_lastbin(h)

        # setup the axes
        x = h.GetXaxis()
        x.SetTitle(plot.x_label)
        x.SetTitleSize(0.048 * 0.85)
        x.SetLabelSize(0.035)
        x.SetLabelOffset(1.15 * 0.02)
        x.SetTitleOffset(0.95 * x.GetTitleOffset())
        x.SetLabelFont(42)
        x.SetTitleFont(42)

        y = h.GetYaxis()
        y.SetTitle("Arb. units")
        y.SetTitleSize(0.055 * 0.85)
        y.SetLabelSize(1.2 * 0.035)
        y.SetLabelOffset(0.013)
        y.SetTitleOffset(1.4)
        y.SetLabelFont(42)
        y.SetTitleFont(42)

        leg.AddEntry(h, b.displayname, "l")
        histos.append(h)
        c.Update()
        maxy.append(h.GetMaximum())

    maxy_ = 1.25*max(maxy)

    is_first = True
    for hist in histos :
        if is_first :
            is_first = False
            hist.SetMaximum(maxy_)
            hist.Draw("hist")
        hist.Draw("hist same")
        hist.SetMaximum(maxy_)

    # legend
    leg.Draw()

    pu.draw_text(text="#it{ATLAS} Internal",x=0.18,y=0.83, size = 0.06)
    #pu.draw_text(text="SF Pre-selection (==0 b-jets)",x=0.18,y=0.78)
    #pu.draw_text_on_top(text=plot.name)

    outname = plot.name + ".eps"
    c.SaveAs(outname)
    out = g.plots_dir + g_outdir
    utils.mv_file_to_dir(outname, out, True)
    fullname = out + "/" + outname
    print "%s saved to : %s"%(outname, os.path.abspath(fullname))

def make_plotsStack(plot, reg, data, backgrounds, plot_i, n_plots, for_ratio=True):
    '''
    Make a stack plot of all backgrounds and data (if defined)

    param:
        plot : plot class
            1D plot defined in configuration file
        reg : region class
            region defined in configuration file
        data : background class
            data sample defined in configuration file
        backgrounds : list(background class)
            background samples defined in configuration file
        plot_i : int
            current plot number
        n_plots : int
            total number of plots

    returns canvas
    '''
    # Sanity Check
    assert not for_ratio or data, "ERROR :: Need data for ratio plot"

    print 20*"-","Plotting [%d/%d] %s"%(plot_i, n_plots, plot.name), 20*'-'

    ############################################################################
    # Intialize plot components

    # Canvas
    can = plot.canvas
    can.cd()
    can.SetFrameFillColor(0)
    can.SetFillColor(0)
    can.SetLeftMargin(0.14)
    can.SetRightMargin(0.05)
    can.SetBottomMargin(1.3*can.GetBottomMargin())

    if plot.isLog() : can.SetLogy(True)
    can.Update()


    hax = r.TH1F("axes", "", int(plot.nbins), plot.x_range_min, plot.x_range_max)
    hax.SetMinimum(plot.y_range_min)
    hax.SetMaximum(plot.y_range_max)
    xax = hax.GetXaxis()
    xax.SetTitle(plot.x_label)
    xax.SetTitleFont(42)
    xax.SetLabelFont(42)
    xax.SetLabelSize(0.035)
    xax.SetTitleSize(0.048 * 0.85)
    if for_ratio:
        hax.GetXaxis().SetTitleOffset(-999)
        hax.GetXaxis().SetLabelOffset(-999)
    else:
        xax.SetLabelOffset(1.15 * 0.02)
        xax.SetTitleOffset(1.5 * xax.GetTitleOffset())

    yax = hax.GetYaxis()
    yax.SetTitle(plot.y_label)
    yax.SetTitleFont(42)
    yax.SetLabelFont(42)
    yax.SetTitleOffset(1.4)
    yax.SetLabelOffset(0.013)
    yax.SetLabelSize(1.2 * 0.035)
    yax.SetTitleSize(0.055 * 0.85)
    if for_ratio:
        hax.Draw()
    else:
        hax.Draw("AXIS")
    can.Update()

    stack = r.THStack("stack_"+plot.name, "")

    # Legend
    leg = pu.default_legend(xl=0.55,yl=0.71,xh=0.93,yh=0.90)
    leg.SetNColumns(2)
    leg_sig = pu.default_legend(xl=0.55, yl=0.6, xh=0.91, yh=0.71)
    leg_sig.SetNColumns(1)

    # Yield table
    yield_tbl = []

    ############################################################################
    # Get background MC stack

    # Initilize lists and defaults
    histos = []
    all_histos = []
    sig_histos = []
    avoid_bkg = []

    #h_nom_fake = None

    n_total_sm_yield = n_total_sm_error = 0.
    hists_to_clear = []

    # Add MC backgrounds to stack
    for b in backgrounds :
        # Initilize histogram
        h_name_tmp = re.sub(r'[(){}[\]]+','',plot.variable)
        h_name = "h_"+reg.name+'_'+b.treename+"_"+h_name_tmp
        hists_to_clear.append(h_name)
        h = pu.th1d(h_name, "", int(plot.nbins),
                    plot.x_range_min, plot.x_range_max,
                    plot.x_label, plot.y_label)

        h.SetLineColor(b.color)
        h.GetXaxis().SetLabelOffset(-999)
        if b.isSignal():
            h.SetLineWidth(2)
            h.SetLineStyle(2)
            h.SetFillStyle(0)
        else:
            h.SetFillColor(b.color)
            h.SetFillStyle(1001)
        h.Sumw2

        weight_str = "eventweight"

        # Draw final histogram (i.e. selections and weights applied)
        cut = "(" + reg.tcut + ") * %s * "%weight_str + str(b.scale_factor)
        cut = r.TCut(cut)
        sel = r.TCut("1")
        draw_cmd = "%s>>+%s"%(plot.variable, h.GetName())
        b.tree.Draw(draw_cmd, cut * sel, "goff")


        # Yield +/- stat error
        stat_err = r.Double(0.0)
        integral = h.IntegralAndError(0,-1,stat_err)
        if not b.isSignal():
            n_total_sm_yield += float(integral)
            n_total_sm_error = (n_total_sm_error**2 + float(stat_err)**2)**0.5
        yield_tbl.append("%10s: %.2f +/- %.2f"%(b.name, integral, stat_err))

        # Add overflow
        if plot.add_overflow:
            pu.add_overflow_to_lastbin(h)
        if plot.add_underflow:
            pu.add_underflow_to_firstbin(h)

        # Record all and non-empty histograms
        if b.isSignal():
            leg_sig.AddEntry(h, b.displayname, "l")
            sig_histos.append(h)
        else:
            all_histos.append(h)
            histos.append(h) if integral > 0 else avoid_bkg.append(b.name)
        can.Update()

    if not len(histos):
        print "make_plotsStack ERROR :: All SM hists are empty. Skipping"
        return

    # Order the hists by total events
    histos = sorted(histos, key=lambda h: h.Integral())
    for h in histos :
        if "fake" in h.GetName() or "Fake" in h.GetName() : continue
        if plot.isNorm(): h.Scale(1.0/float(n_total_sm_yield))
        stack.Add(h)
    can.Update()


    yield_tbl.append("%10s: %.2f +/- %.2f"%('Total SM', n_total_sm_yield, n_total_sm_error))

    h_leg = sorted(all_histos, key=lambda h: h.Integral(), reverse=True)
    histos_for_leg = histos_for_legend(h_leg)

    ############################################################################
    # Get data hist
    if data:
        hd_name = "h_"+reg.name+'_data_'+plot.variable
        hd = pu.th1d(hd_name, "", int(plot.nbins),
                                  plot.x_range_min, plot.x_range_max,
                                  plot.x_label, plot.y_label)
        hd.Sumw2

        cut = "(" + reg.tcut + ")"
        cut = r.TCut(cut)
        sel = r.TCut("1")
        draw_cmd = "%s>>%s"%(plot.variable, hd.GetName())
        data.tree.Draw(draw_cmd, cut * sel, "goff")
        hd.GetXaxis().SetLabelOffset(-999)

        # print the yield +/- stat error
        stat_err = r.Double(0.0)
        integral = hd.IntegralAndError(0,-1,stat_err)
        yield_tbl.append("%10s: %.2f +/- %.2f"%('Data',integral, stat_err))
        yield_tbl.append("%10s: %.2f"%("Data/MC", integral/n_total_sm_yield))
        #print "Data: %.2f +/- %.2f"%(integral, stat_err)

        if plot.isNorm():
            hd.Scale(1.0/float(integral))

        # Add overflow
        if plot.add_overflow:
            pu.add_overflow_to_lastbin(h)
        if plot.add_underflow:
            pu.add_underflow_to_firstbin(h)

        gdata = pu.convert_errors_to_poisson(hd)
        #gdata.SetLineWidth(2)
        #uglify
        gdata.SetLineWidth(1)
        gdata.SetMarkerStyle(20)
        gdata.SetMarkerSize(1.5)
        gdata.SetLineColor(1)
        leg.AddEntry(gdata, "Data", "p")
        can.Update()

    ############################################################################
    # systematics loop
    r.gStyle.SetHatchesSpacing(0.9)

    mcError = r.TH1F("mcError", "mcError", 2,0,2)
    mcError.SetLineWidth(3)
    mcError.SetFillStyle(3345)
    mcError.SetFillColor(r.kBlue)
    mcError.SetLineColor(r.kBlack)
    leg.AddEntry(mcError, "Standard Model", "fl")

    # now add backgrounds to legend
    for h in histos_for_leg :
        name_for_legend = ""
        for b in backgrounds :
            if b.treename in h.GetName() :
                name_for_legend = b.displayname
        leg.AddEntry(h, name_for_legend, "f")

    totalSM = stack.GetStack().Last().Clone("totalSM")
    nominalAsymErrors = pu.th1_to_tgraph(totalSM)
    nominalAsymErrors.SetMarkerSize(0)
    nominalAsymErrors.SetLineWidth(0)
    nominalAsymErrors.SetFillStyle(3345)
    nominalAsymErrors.SetFillColor(r.kBlue)

    ############################################################################
    # draw the MC stack and do cosmetics
    stack.SetMinimum(plot.y_range_min)

    # Determine y range for plot
    max_mult = 2.0 if any([b.isSignal() for b in backgrounds]) else 1.66
    if data:
        maxy = max(hd.GetMaximum(), stack.GetMaximum())
    else:
        maxy = stack.GetMaximum()
    if not plot.isLog() :
        hax.SetMaximum(max_mult*maxy)
        hax.Draw()
        can.Update()
        stack.SetMaximum(max_mult*maxy)
    else :
        hax.SetMaximum(1e3*plot.y_range_max)
        hax.Draw()
        can.Update()
        stack.SetMaximum(1e3*plot.y_range_max)

    # Add stack to canvas
    stack.Draw("HIST SAME")
    can.Update()

    # symmetrize the errors
    for i in xrange(nominalAsymErrors.GetN()) :
        ehigh = nominalAsymErrors.GetErrorYhigh(i)
        elow  = nominalAsymErrors.GetErrorYlow(i)


        error_sym = r.Double(0.0)
        error_sym = (ehigh + elow) / 2.0

        if ehigh != error_sym:
            print "initial error (+%.2f,-%.2f), symmetrized = (+%.2f,-%.2f)"%(ehigh,elow, error_sym, error_sym)


        nominalAsymErrors.SetPointEYhigh(i,0.0)
        nominalAsymErrors.SetPointEYhigh(i, error_sym)
        nominalAsymErrors.SetPointEYlow(i,0.0)
        nominalAsymErrors.SetPointEYlow(i,error_sym)

    # draw the error band
    nominalAsymErrors.Draw("same && E2")

    # draw the total bkg line
    hist_sm = stack.GetStack().Last().Clone("hist_sm")
    hist_sm.SetLineColor(r.kBlack)
    hist_sm.SetLineWidth(mcError.GetLineWidth())
    hist_sm.SetLineStyle(1)
    hist_sm.SetFillStyle(0)
    hist_sm.SetLineWidth(3)
    hist_sm.Draw("hist same")

    #draw the signals
    for hsig in sig_histos :
        hsig.Draw("hist same")

    ############################################################################
    # Print yield table
    print "Yields for %s region"%reg.displayname
    print 30*'-'
    for yld in yield_tbl:
        if 'total' in yld.lower() or 'Data/MC' in yld:
            print 20*'-'
        print yld
    print 30*'-'

    ############################################################################
    # Finalize the plot

    # draw the data graph
    if data:
        gdata.Draw("option same pz 0")

    # draw the legend
    leg.Draw()
    leg_sig.Draw()
    r.gPad.RedrawAxis()

    # add some text/labels
    pu.draw_text(text="ATLAS",x=0.18,y=0.85,size=0.05,font=72)
    pu.draw_text(text="Internal",x=0.325,y=0.85,size=0.05,font=42)
    pu.draw_text(text="#sqrt{s} = 13 TeV, 36.1 fb^{-1}", x=0.18, y=0.79, size=0.04)
    pu.draw_text(text="Higgs LFV", x=0.18, y=0.74, size=0.04)
    pu.draw_text(text=reg.displayname,      x=0.18,y=0.68, size=0.04)

    r.gPad.SetTickx()
    r.gPad.SetTicky()

    can.Update()
    ############################################################################
    # Save histogram
    if for_ratio:
        return can, stack, gdata, nominalAsymErrors
    else:
        print "\n"
        outname = plot.name + ".eps"
        can.SaveAs(outname)
        out = g.plots_dir + g_outdir
        utils.mv_file_to_dir(outname, out, True)
        fullname = out + "/" + outname
        #print "%s saved to : %s"%(outname, os.path.abspath(fullname))
        #return True
################################################################################
def make_plots2D(plot, reg, data, backgrounds) :

    # check if we want to do profile vs. TH2F plots
    if plot.do_profile :
        make_1dprofile(plot, reg, data, backgrounds)
        return
    elif plot.do_profileRMS :
        make_1dprofileRMS(plot, reg, data, backgrounds)
        return
    print "make_plots2D    Plotting %s"%plot.name

    # set the palette/colors
    #pu.set_palette(name="redbluevector")
    pu.set_palette(name="")

    # get the canvas from the plot
    c = plot.canvas
    c.cd()
    # these should be in the default cavnas setting in utils/plot.py
    c.SetFrameFillColor(0)
    c.SetFillColor(0)
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.14)
    c.SetBottomMargin(1.3*c.GetBottomMargin())

    # name_on_plot : text to put on top of the plot pad in addition to the region name so that
    # we know which sample is being plotted on the TH2D
    name_on_plot = ""

    ### MC
    if plot.sample != "Data" :
        for b in backgrounds :
            if b.name != plot.sample : continue
            name_on_plot = b.displayname

            hist_name_x = ""
            hist_name_y = ""
            if "abs" in plot.xVariable :
                x_repl = plot.xVariable.replace("abs(","")
                x_repl = x_repl.replace(")","")
                hist_name_x = x_repl
            else : hist_name_x = plot.xVariable

            if "abs" in plot.yVariable :
                y_repl = plot.yVariable.replace("abs(","")
                y_repl = y_repl.replace(")","")
                hist_name_y = y_repl
            else : hist_name_y = plot.yVariable

            h = pu.th2f("h_"+b.name+"_"+hist_name_x+"_"+hist_name_y, "", int(plot.n_binsX), plot.x_range_min, plot.x_range_max, int(plot.n_binsY), plot.y_range_min, plot.y_range_max, plot.x_label, plot.y_label)

            # get the cut, and weight the sample (for TH2 the scale_factor is not so important if we normalize)
            cut = "(" + reg.tcut + ") * eventweight * " + str(b.scale_factor)
            cut = r.TCut(cut)
            sel = r.TCut("1")
            cmd = "%s:%s>>+%s"%(plot.yVariable,plot.xVariable,h.GetName())
            b.tree.Draw(cmd, cut * sel)

    ### DATA
    if plot.sample == "Data" :
        name_on_plot = "Data"
        h = pu.th2f("h_"+data.name+"_"+plot.xVariable+"_"+plot.yVariable, "", int(plot.n_binsX), plot.x_range_min, plot.x_range_max, int(plot.n_binsY), plot.y_range_min, plot.y_range_max, plot.x_label, plot.y_label)

        cut = "(" + reg.tcut + ")"
        cut = r.TCut(cut)
        sel = r.TCut("1")
        cmd = "%s:%s>>+%s"%(plot.yVariable,plot.xVariable,h.GetName())
        data.tree.Draw(cmd, cut * sel)

    # normalize the histgrams to unity
    if h.Integral()!=0 : h.Scale(1/h.Integral())

    h.Smooth()
    h.GetYaxis().SetTitleOffset(1.52 * h.GetYaxis().GetTitleOffset())
    h.GetXaxis().SetTitleOffset(1.2 * h.GetXaxis().GetTitleOffset())
    #r.TGaxis.SetMaxDigits(2)
    #h.GetZaxis().SetLabelSize(0.75 * h.GetLabelSize())

    gr = r.TGraph2D(1)
    gr.SetMarkerStyle(r.kFullSquare)
    gr.SetMarkerSize(2.0 * gr.GetMarkerSize())
    gr.SetHistogram(h)
    gr.Draw(plot.style)

    h.Draw(plot.style)

    # write descriptive text on top of the pad
    pu.draw_text_on_top(text="DF pre-selection + b-veto: #bf{%s}"%(name_on_plot))
    #pu.draw_text_on_top(text="%s : #bf{%s}"%(plot.name,name_on_plot))

    c.Update()

    # redraw the now-invisible axes ticks
    r.gPad.RedrawAxis()

    outname = plot.name+".eps"
    out = g.plots_dir + g_outdir
    c.SaveAs(outname)
    utils.mv_file_to_dir(outname, out, True)
    fullname = out + "/" + outname
    print "%s saved to : %s"%(outname, os.path.abspath(fullname))

def make_plots1D(plot, reg, data, backgrounds, plot_i, n_plots) :
    '''
    Determine 1D plot type and run appropriate 1D plot method

    param:
        plot : plot class
            1D plot defined in configuration file
        reg : region class
            region defined in configuration file
        data : background class
            data sample defined in configuration file
        backgrounds : list(background class)
            background samples defined in configuration file
        plot_i : int
            current plot number
        n_plots : int
            total number of plots

    returns (None)
    '''
    if plot.ratioCanvas :
        make_plotsRatio(plot, reg, data, backgrounds, plot_i, n_plots)
    elif plot.is_comparison :
        make_plotsComparison(plot, reg, data, backgrounds)
    else :
        make_plotsStack(plot, reg, data, backgrounds, plot_i, n_plots, for_ratio=False)

################################################################################
def make_plots(plots, regions, data, backgrounds) :
    '''
    Make sure that the plots are not asking for undefined region

    param:
        plots : list(plot class)
            plots defined in config file
        regions : list(region class)
            regions defined in config file
        data : background class
        backgrounds: list(background class)

    returns (None)
    '''
    for reg in regions:
        # first check that there are plots for the given region
        plots_with_region = [p for p in plots if p.region == reg.name]
        if not len(plots_with_region): continue

        print 20*'-', "Plots for %s region"%reg.displayname, 20*'-','\n'
        ########################################################################
        # Set event lists, if they already exist load it.
        # Otherwise make it and save

        print "Setting EventLists for %s"%reg.name
        cut = r.TCut(reg.tcut)
        sel = r.TCut("1")

        # EventLists for Backgrounds
        for b in backgrounds :
            save_cut  = r.TCut(sel*cut)
            list_name = "list_" + reg.name + "_" + b.treename
            save_name = g.event_list_dir + list_name + ".root"
            save_name = os.path.normpath(save_name)
            if g_dbg: print 'List name =',list_name

            # Reset event list
            b.tree.SetEventList(0)

            # check if the list already exists
            load_eventlist = True
            if os.path.isfile(save_name) :
                rfile = r.TFile.Open(save_name)
                print "%10s : EventList found at %s"%(b.name, save_name)
                if save_cut != rfile.Get("cut"):
                    print "EventList cuts have changed. Remaking EventList"
                    load_eventlist = False
            else:
                load_eventlist = False

            if load_eventlist:
                list = rfile.Get(list_name)
                b.tree.SetEventList(list)
                if g_dbg : list.Print()
            else:
                print "Creating TEventList for", b.name
                rfile = r.TFile(save_name,'recreate')
                draw_list = ">> " + list_name
                b.tree.Draw(draw_list, sel*cut)
                list = r.gROOT.FindObject(list_name).Clone()
                b.tree.SetEventList(list)
                list.Write(list_name)
                save_cut.Write("cut")
                rfile.Close()


        # systematics trees
        # if g_doSys :
        #     cut = reg.tcut
        #     cut = r.TCut(cut)
        #     sel = r.TCut("1")
        #     for b in backgrounds :
        #         if 'fakes' in b.name : continue
        #         if b.isSignal() : continue
        #         for s in b.systList :
        #             if not s.isKinSys() : continue

        #             # up variation
        #             list_name_up = "list_" + reg.name + "_" + b.treename + "_" + s.name + "_UP"
        #             save_name_up = g.event_list_dir + list_name_up + ".root"

        #             if os.path.isfile(save_name_up) :
        #                 rfile = r.TFile.Open(save_name_up)
        #                 list = rfile.Get(list_name_up)
        #                 print "%s : EventList found at %s"%(b.name, os.path.abspath(save_name_up))
        #                 if g_dbg : list.Print()
        #                 s.tree_up.SetEventList(list)
        #             else :
        #                 draw_list = ">> " + list_name_up
        #                 s.tree_up.Draw(draw_list, sel*cut)
        #                 list = r.gROOT.FindObject(list_name_up)
        #                 print list_name_up
        #                 s.tree_up.SetEventList(list)
        #                 list.SaveAs(save_name_up)

        #             # down variation

        #             list_name_dn = "list_" + reg.name + "_" + b.treename + "_" + s.name + "_DN"
        #             save_name_dn = g.event_list_dir + list_name_dn + ".root"

        #             if not s.isOneSided() :
        #                 if os.path.isfile(save_name_dn) :
        #                     rfile = r.TFile.Open(save_name_dn)
        #                     list = rfile.Get(list_name_dn)
        #                     print "%s : EventList found at %s"%(b.name, os.path.abspath(save_name_dn))
        #                     s.tree_down.SetEventList(list)
        #                 else :
        #                     draw_list = ">> " + list_name_dn
        #                     s.tree_down.Draw(draw_list, sel*cut)
        #                     list = r.gROOT.FindObject(list_name_dn)
        #                     s.tree_down.SetEventList(list)
        #                     list.SaveAs(save_name_dn)

        # EventLists for Data
        if data :
            data_list_name = "list_" + reg.name + "_" + data.treename
            data_save_name = g.event_list_dir + data_list_name + ".root"
            data_save_name = os.path.normpath(data_save_name)

            # Reset event list
            data.tree.SetEventList(0)

            if os.path.isfile(data_save_name) :
                rfile = r.TFile.Open(data_save_name)
                data_list = rfile.Get(data_list_name)
                print "%10s : EventList found at %s"%('Data', data_save_name)
                if g_dbg : data_list.Print()
                data.tree.SetEventList(data_list)
            else :
                draw_list = ">> " + data_list_name
                data.tree.Draw(draw_list, sel * cut)
                data_list = r.gROOT.FindObject(data_list_name)
                data.tree.SetEventList(data_list)
                data_list.SaveAs(data_save_name)

        ########################################################################
        # Make 1D and 2D plots
        n_total = len(plots_with_region)
        for n_current, p in enumerate(plots_with_region, 1):
            if not p.is2D : make_plots1D(p, reg, data, backgrounds, n_current, n_total)
            elif p.is2D : make_plots2D(p, reg, data, backgrounds)
        print 30*'-', '%s region completed'%reg.displayname, 30*'-'
    print 30*'-', 'PLOTS COMPLETED', 30*'-','\n'

def check_for_consistency(plots, regions) :
    '''
    Make sure that the plots are not asking for undefined region

    param:
        plots : list(plot class)
            plots defined in config file
        regions : list(region class)
            regions defined in config file

    returns (None)
    '''
    region_names = [r.name for r in regions]
    bad_regions = [p.region for p in plots if p.region not in region_names]
    if len(bad_regions) > 0 :
        print 'check_for_consistency ERROR    '\
        'You have configured a plot for a region that is not defined. '\
        'Here is the list of "bad regions":'
        print bad_regions
        print 'check_for_consistency ERROR    The regions that are defined in the configuration ("%s") are:'%g_plotConfig
        print region_names
        print "check_for_consistency ERROR    Exiting."
        sys.exit()
    else :
        print "check_for_consistency    Plots and regions consistent."

def get_plotConfig(conf) :
    """
    Modify path and check for file

    param:
        conf : str
            name of configuration file

    returns:
        str:
            path to configuration file
    """
    configuration_file = conf
    if not configuration_file.endswith(".py"):
        configuration_file += '.py'
    configuration_file = os.path.normpath(configuration_file)
    if os.path.isfile(configuration_file) :
        return configuration_file
    else :
        print 'get_plotConfig ERROR    '\
              'Input g_plotConfig ("%s") is not found in the directory/path.'\
              'Does it exist? Exitting.'%(configuration_file)
        sys.exit()

################################################################################
if __name__=="__main__" :
    global g_plotConfig, g_requestRegion, g_requestPlot, g_outdir, g_dbg

    parser = ArgumentParser()
    parser.add_argument("-c", "--plotConfig",
                            default="",
                            help='provide the name of the config file')
    parser.add_argument("-s", "--doSys", action="store_true",
                            default=False)
    parser.add_argument("-r", "--requestRegion",
                            default="",
                            help='request a region to plot -- will make all plots in the config that are in this region')
    parser.add_argument("-p", "--requestPlot",
                            default="",
                            help='request a specific plot -- provide the name of the plot')
    parser.add_argument("-o", "--outdir",
                            default="./",
                            help='provide the name of the output directory to save the plots. Will be under \'indir/plots/\'')
    parser.add_argument("-d", "--dbg", action="store_true",
                            default=False,
                            help='set debug/verbosity to True')

    args = parser.parse_args()
    g_plotConfig      = args.plotConfig
    g_doSys           = args.doSys
    g_requestRegion   = args.requestRegion
    g_requestPlot     = args.requestPlot
    g_outdir          = args.outdir
    g_dbg             = args.dbg

    if g_requestRegion != "" and g_requestPlot != "" :
        print 'ERROR    You have requested both a reagion ("%s") AND a plot ("%s").'%(g_requestRegion, g_requestPlot)
        print 'ERROR    You may only request one at a time. Exitting.'
        sys.exit()

    print " ++ ------------------------- ++ "
    print "      plotter                    "
    print "                                 "
    print " plot config      :  %s          "%g_plotConfig
    print " requested region :  %s          "%g_requestRegion
    print " requested plot   :  %s          "%g_requestPlot
    print " systematics      :  %s          "%g_doSys
    print " output directory :  %s          "%g_outdir
    print " debug            :  %s          "%g_dbg
    print "                                 "
    print " ++ ------------------------- ++ \n"

    # get the config file
    conf_file = get_plotConfig(g_plotConfig)
    print "Found the configuration file: %s"%conf_file

    # variables to be filled
    plots = []
    data = None
    backgrounds = []
    systematics = []
    regions = []
    execfile(conf_file)

    check_for_consistency(plots, regions)

    # print out the loaded backgrounds and plots
    if g_dbg :
        for p in plots :
            p.Print()
    print "+-----------------------+ "
    print "Loaded backgrounds:    "
    for b in backgrounds :
        print '\t',
        b.Print()
    print "Loaded data sample:    "
    if data :
        print '\t',
        data.Print()
    print "+-----------------------+ "
    # if g_doSys :
    #     print "+-----------------------+ "
    #     print "  Loaded systematics:     "
    #     for s in systematics :
    #         s.check()
    #         s.Print()
    #         for b in backgrounds :
    #             if b.isSignal() : continue
    #             b.addSys(s)
    #             if g_dbg :
    #                 for s in b.systList :
    #                     print s.tree
    #     print "+-----------------------+ "

    # make the plots
    if g_requestRegion != "" :
        requested_plots = []
        for p in plots :
            if p.region == g_requestRegion : requested_plots.append(p)
        make_plots(requested_plots, regions, data, backgrounds)
    elif g_requestPlot != "" :
        requested_plots = []
        for p in plots :
            if p.name == g_requestPlot : requested_plots.append(p)
        make_plots(requested_plots, regions, data, backgrounds)
    else :
        make_plots(plots, regions, data, backgrounds)
