#!/usr/bin/env python
import pdb
"""
================================================================================
As comprehensive a plotting script as there is this side of the Mississipi

Examples
    ./plotter.py -c config_file.py

Description:
    Creates various histogram plots (e.g. stack, ratio plot, overlay, etc.)
    composed of several samples (i.e. data, background, and signal). The script
    requires a detailed configuration file be provided in that defines
    backgrounds, data, systematics, regions, and plots.


Author:
    Alex Armstrong <alarmstr@cern.ch>
    ... with lots of ideas borrowed from Danny Antrim <dantrim@cern.ch>

License:
    Copyright: (C) <May 16th, 2018>; University of California, Irvine
================================================================================
"""
# General python
import sys, os, traceback, argparse
import time
import importlib
import re

# Root data analysis framework
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal cmd-line options
r.gROOT.SetBatch(True)
r.gStyle.SetOptStat(False)

r.TEventList.__init__._creates = False
r.TH1F.__init__._creates = False
r.TGraphErrors.__init__._creates = False
r.TGraphAsymmErrors.__init__._creates = False

# Local classes for plotting
import tools.plot_utils as pu
import tools.utils as utils
import tools.signal as signal
import tools.background as background
import tools.region as region
import tools.plot as plot
from global_variables import event_list_dir, plots_dir

################################################################################
def main ():
    """ Main Function """

    global args
    global yield_tbl
    # Slim plots to make if requested
    rReg = args.requestRegion
    rPlot = args.requestPlot
    plots = conf.plots
    if not rReg and not rPlot:
        pass # Use all plots
    elif rReg and not rPlot:
        plots = [p for p in plots if p.region == rReg]
    elif not rReg and rPlot:
        plots = [p for p in plots if p.name == rPlot]
    elif rReg and rPlot:
        plots = [p for p in plots if p.name == rPlot and p.region == rReg]

    make_plots(plots, conf.regions, conf.data, conf.backgrounds)

################################################################################
# PLOTTING FUNCTIONS
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

    '''
    for reg in regions:
        # first check that there are plots for the given region
        plots_with_region = [p for p in plots if p.region == reg.name]
        if not len(plots_with_region): continue

        print '\n', 20*'-', "Plots for %s region"%reg.displayname, 20*'-', '\n'

        ########################################################################
        # Set event lists, if they already exist load it.
        # Otherwise make it and save

        #TODO: Modularize out
        # - Loop over regions
        # - Loop over samples = backgrounds + data
        # - event_lists = get_event_list(sample, selection, list_dir)
        # - b.tree.SetEventLists(event_lists)

        print "Setting EventLists for %s"%reg.name
        cut = r.TCut(reg.tcut)
        sel = r.TCut("1")

        # EventLists for Backgrounds
        for b in backgrounds :
            save_cut  = r.TCut(sel*cut)
            list_name = "list_" + reg.name + "_" + b.treename
            save_name = event_list_dir + list_name + ".root"
            save_name = os.path.normpath(save_name)
            if args.verbose: print 'List name =',list_name

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
                if args.verbose : list.Print()
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

        # EventLists for Data
        if data :
            data_list_name = "list_" + reg.name + "_" + data.treename
            data_save_name = event_list_dir + data_list_name + ".root"
            data_save_name = os.path.normpath(data_save_name)

            # Reset event list
            data.tree.SetEventList(0)

            if os.path.isfile(data_save_name) :
                rfile = r.TFile.Open(data_save_name)
                data_list = rfile.Get(data_list_name)
                print "%10s : EventList found at %s"%('Data', data_save_name)
                if args.verbose : data_list.Print()
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
            if p.is2D :
                make_plots2D(p, reg, data, backgrounds)
            else:
                make_plots1D(p, reg, data, backgrounds, n_current, n_total)


    print 30*'-', 'PLOTS COMPLETED', 30*'-','\n'

################################################################################
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
        make_plotsStack(plot, reg, data, backgrounds, plot_i, n_plots)

def make_plots2D(plot, reg, data, backgrounds):
        print "TODO"

################################################################################
class PlotParts():
    def __init__(self):
        self.hist_axis = None
        self.mc_stack = None
        self.mc_hist = None
        self.signal_hists = []
        self.data_hist = None
        self.error_graph = None
        self.leg = None
    def draw_all(self, region_name):
        if not self.hist_axis:
            print "ERROR :: Hist axis not provided. Cannot draw "
            return
        self.hist_axis.Draw('AXIS')
        if self.mc_stack: self.mc_stack.Draw("hist && same")
        if self.error_graph: self.error_graph.Draw("E2 && same")
        if self.mc_hist: self.mc_hist.Draw('hist same')
        for hsig in self.signal_hists: hsig.Draw("hist same")
        if self.data_hist: self.data_hist.Draw("option same pz 0")
        if self.leg: self.leg.Draw('same')

        # add some text/labels
        pu.draw_text(text="ATLAS",x=0.18,y=0.85,size=0.05,font=72)
        pu.draw_text(text="Internal",x=0.325,y=0.85,size=0.05,font=42)
        pu.draw_text(text="#sqrt{s} = 13 TeV, 36.1 fb^{-1}", x=0.18, y=0.79, size=0.04)
        pu.draw_text(text="Higgs LFV", x=0.18, y=0.74, size=0.04)
        pu.draw_text(text=region_name,      x=0.18,y=0.68, size=0.04)

        # Reset axis
        r.gPad.RedrawAxis()
        r.gPad.SetTickx()
        r.gPad.SetTicky()

    def Print(self):
        print "self.hist_axis = ", self.hist_axis
        print "self.mc_stack = ", self.mc_stack
        print "self.error_graph = ", self.error_graph
        print "self.mc_hist = ", self.mc_hist
        print "self.signal_hists = ", self.signal_hists
        print "self.data_hist = ", self.data_hist
        print "self.leg = ", self.leg

def make_plotsStack(plot, reg, data, backgrounds, plot_i, n_plots):
    ''' '''
    print 20*"-","Plotting [%d/%d] %s"%(plot_i, n_plots, plot.name), 20*'-'

    # Prepare canvas
    can = plot.canvas
    can.cd()
    can.SetFrameFillColor(0)
    can.SetFillColor(0)
    can.SetLeftMargin(0.14)
    can.SetRightMargin(0.05)
    can.SetBottomMargin(1.3*can.GetBottomMargin())
    if plot.isLog() : can.SetLogy(True)

    # Initialize plot components
    # - allows them to persist when function scope ends
    # - thereby preserving them on the canvas
    stack_plot_parts = PlotParts()
    #stack_plot_parts = None
    #hist_axis = None
    #mc_stack = None
    #mc_hist = None
    #signal_hists = []
    #data_hist = None
    #error_graph = None
    #leg = None
    draw_stack(can, plot, reg, data, backgrounds, stack_plot_parts)
    #stack_plot_parts.draw_all()
    can.SaveAs("InProgress5.pdf")
    can.Update()

    can.SaveAs("InProgress6.pdf")
    save_plot(can, plot.name+".eps")

def make_plotsRatio(plot, reg, data, backgrounds, plot_i, n_plots) :
    ''' '''
    print 20*"-","Plotting [%d/%d] %s"%(plot_i, n_plots, plot.name), 20*'-'
    ############################################################################
    # Intialize plot components
    print "Got Nothing Yet"

################################################################################
def draw_stack(can, plot, reg, data, backgrounds, plot_parts):
    ''' '''
    #TODO: Remove dependence on passing canvas to functions
    # Get histograms
    #global hist_axis, mc_stack, mc_hist, signal_hists, data_hist, error_graph, leg
    make_stack_hists(can, plot, reg, data, backgrounds, plot_parts)

    #TODO: Seperate histogram making into seperate functions
    #mc_stack = make_hist1D(plot, reg, backgrounds, stack=true)
    #mc_total_hist = stack_to_total(mc_stack)
    #signals_hists = make_hist1D(plot, reg, signals)
    #data_hist = make_hist1D(plot, reg, data)
    #error_graph = make_error_graph(mc_hist)
    #
    ## Prepare legend
    #leg = make_legend(data=data_hist, stack=mc_stack, signals=signals_hists)
    #
    # Draw onto current canvas
    #can.cd()
    #hist_axis.Draw('AXIS')
    #mc_stack.Draw("hist && same")
    #error_graph.Draw("E2 && same")
    #mc_hist.Draw('hist same')
    #for hsig in signal_hists: hsig.Draw("hist same")
    #data_hist.Draw("option same pz 0")
    #leg.Draw('same')
    plot_parts.draw_all(reg.displayname)
    #can.Update()
    can.SaveAs("InProgress4.pdf")

def save_plot(can, outname):
    can.SaveAs("InProgress7.pdf")
    can.SaveAs(outname)
    out = plots_dir + args.outdir
    utils.mv_file_to_dir(outname, out, True)
    #fullname = out + "/" + outname
    #print "%s saved to : %s"%(outname, os.path.abspath(fullname))

################################################################################
def make_stack_hists(can, plot, reg, data, backgrounds, plot_parts, for_ratio=False):
    ''' '''

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
        if b.isSignal() : continue
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
        #TODO: return NONE and add flags downstream to handle
        sys.exit()

    # Order the hists by total events
    histos = sorted(histos, key=lambda h: h.Integral())
    for h in histos :
        if "fake" in h.GetName() or "Fake" in h.GetName() : continue
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

    can.Update()
    can.SaveAs("InProgress3.pdf")

    plot_parts.hist_axis = hax
    plot_parts.mc_stack = stack
    plot_parts.mc_hist = hist_sm
    plot_parts.signal_hists = sig_histos
    plot_parts.data_hist = gdata
    plot_parts.error_graph = nominalAsymErrors
    plot_parts.leg = leg
    #return (hax, stack, hist_sm, sig_histos, gdata, nominalAsymErrors, leg)
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
# SETUP FUNCTIONS
def check_args(args):
    """ Check the input arguments are as expected """
    configuration_file = os.path.normpath(args.plotConfig)
    if not os.path.exists(configuration_file):
        print "ERROR :: Cannot find config file:", configuration_file
        sys.exit()

def check_environment():
    """ Check if the shell environment is setup as expected """
    python_ver = sys.version_info[0] + 0.1*sys.version_info[1]
    if python_ver < 2.7:
        print "ERROR :: Running old version of python\n", sys.version
        sys.exit()

def print_inputs(args):
    """ Print the program inputs """
    full_path = os.path.abspath(__file__)
    prog_name = os.path.basename(full_path)
    prog_dir = os.path.dirname(full_path)

    print " ==================================================================\n"
    print " Program : %s "%prog_name
    print " Run from: %s "%prog_dir
    print ""
    print " Options:"
    print "     plot config      :  %s "%args.plotConfig
    print "     specific plot    :  %s "%args.requestPlot
    print "     specific region  :  %s "%args.requestRegion
    print "     output directory :  %s "%args.outdir
    print "     verbose          :  %s "%args.verbose
    print ""
    print "===================================================================\n"

    # print out the loaded backgrounds and plots
    print " ============================"
    if conf.backgrounds :
        print "Loaded backgrounds:    "
        for b in conf.backgrounds :
            print '\t',
            b.Print()
    if conf.data :
        print "Loaded data sample:    "
        print '\t',
        conf.data.Print()
    if args.verbose:
        print "Loaded plots:"
        for p in conf.plots :
            p.Print()
    print " ============================"

def check_for_consistency(plots, regions) :
    '''
    Make sure that the plots are not asking for undefined region

    param:
        plots : list(plot class)
            plots defined in config file
        regions : list(region class)
            regions defined in config file
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

def print_hello_world(param=''):
    """
    function synopsis
    args:
        param (type) - description [default: '']
    returns:
        (type) - description
    """
    print "Hello World! " + param
    return True if param else False

def get_plotConfig(configuration_file) :
    """
    Modify path and check for file

    param:
        configuration_file : str
            name of configuration file

    returns:
        str:
            path to configuration file
    """
    configuration_file = os.path.normpath(configuration_file)
    if os.path.isfile(configuration_file) :
        return configuration_file
    else :
        print 'get_plotConfig ERROR    '\
              'Input g_plotConfig ("%s") is not found in the directory/path.'\
              'Does it exist? Exitting.'%(configuration_file)
        sys.exit()

################################################################################
# Run main when not imported
if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = argparse.ArgumentParser(
                description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument("-c", "--plotConfig",
                                default="",
                                help='name of the config file')
        parser.add_argument("-r", "--requestRegion",
                                default="",
                                help='request a region to plot -- will make all plots in the config that are in this region')
        parser.add_argument("-p", "--requestPlot",
                                default="",
                                help='request a specific plot -- provide the name of the plot')
        parser.add_argument("-o", "--outdir",
                                default="./",
                                help='name of the output directory to save plots.')
        parser.add_argument("-v", "--verbose",
                                action="store_true",
                                help='set verbosity mode')
        args = parser.parse_args()

        if args.verbose:
            print '>'*40
            print 'Running {}...'.format(os.path.basename(__file__))
            print time.asctime()

        check_args(args)
        check_environment()

        # Import configuration file
        import_conf = args.plotConfig.replace(".py","")
        conf = importlib.import_module(import_conf)

        check_for_consistency(conf.plots, conf.regions)
        print_inputs(args)

        yield_tbl = []
        main()

        if args.verbose:
            print time.asctime()
            time = (time.time() - start_time)
            print 'TOTAL TIME: %fs'%time,
            print ''
            print '<'*40
    except KeyboardInterrupt, e: # Ctrl-C
        print 'Program ended by keyboard interruption'
        raise e
    except SystemExit, e: # sys.exit()
        print 'Program ended by system exit'
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        os._exit(1)

