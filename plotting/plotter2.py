#!/usr/bin/env python
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
from collections import OrderedDict

# Root data analysis framework
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal cmd-line options
r.gROOT.SetBatch(True)
r.gStyle.SetOptStat(False)

# Prevent root from taking ownership when creating classes
r.TEventList.__init__._creates = False
r.TH1F.__init__._creates = False
r.TGraphErrors.__init__._creates = False
r.TGraphAsymmErrors.__init__._creates = False
r.TLegend.__init__._creates = False

# Local classes for plotting
import tools.plot_utils as pu
import tools.utils as m_utils
import tools.sample as m_sample
import tools.region as m_region
import tools.plot as m_plot
from global_variables import event_list_dir, plots_dir

################################################################################
def main ():
    """ Main Function """

    global args
    # Slim plots to make if requested
    rReg = args.requestRegion
    rPlot = args.requestPlot
    if not rReg and not rPlot:
        pass # Use all plots
    elif rReg and not rPlot:
        PLOTS = [p for p in PLOTS if p.region == rReg]
    elif not rReg and rPlot:
        PLOTS = [p for p in PLOTS if p.name == rPlot]
    elif rReg and rPlot:
        PLOTS = [p for p in PLOTS if p.name == rPlot and p.region == rReg]

    make_plots()

################################################################################
# PLOTTING FUNCTIONS
def make_plots() :
    ''' '''
    for reg in REGIONS:
        # Get plots for this region
        plots_with_region = [p for p in PLOTS if p.region == reg.name]
        if not len(plots_with_region): continue

        print '\n', 20*'-', "Plots for %s region"%reg.displayname, 20*'-', '\n'

        ########################################################################
        print "Setting EventLists for %s"%reg.name
        cut = r.TCut(reg.tcut)
        for sample in SAMPLES :
            list_name = "list_" + reg.name + "_" + sample.name
            sample.set_event_list(cut, list_name, event_list_dir)

        ########################################################################
        # Make 1D and 2D plots
        n_total = len(plots_with_region)
        for n_current, plot in enumerate(plots_with_region, 1):
            print 20*"-","Plotting [%d/%d] %s"%(n_current, n_total, plot.name), 20*'-'
            if p.is2D :
                make_plots2D(plot, reg)
            else:
                make_plots1D(plot, reg)

    print 30*'-', 'PLOTS COMPLETED', 30*'-','\n'

################################################################################
def make_plots1D(plot, reg) :
    ''' '''
    if plot.ptype == m_plot.Types.ratio :
        make_plotsRatio(plot, reg)
    elif plot.ptype == m_plot.Types.comparison :
        make_plotsComparison(plot, reg)
    elif plot.ptype == m_plot.Types.stack :
        make_plotsStack(plot, reg)
    else:
        print "ERROR :: Unknown plot type,", plot.ptype.name

def make_plots2D(plot, reg):
        print "TODO"

################################################################################
def make_plotsStack(plot, reg):
    ''' '''
    #TODO: Modularize so this can be reused for ratio plot code
    # Get Canvas
    can = plot.pads.canvas
    can.cd()
    if plot.doLogY : can.SetLogy(True)

    # Track yields for printing
    # TODO: Create a Yield Table class
    yield_tbl = []

    # Get plotting primitives
    # Not all primitives are for drawing. Some are for preserving pointers
    legend = make_stack_legend()
    axis = make_stack_axis(plot)
    mc_stack, mc_total, signals, hists = add_stack_backgrounds(plot, legend, reg, yield_tbl)
    data, data_hist = add_stack_data(plot, legend, reg, yield_tbl, mc_total.Integral())
    error_leg, mc_errors = add_stack_mc_errors(plot, legend, hists, mc_stack)
    reformat_axis(plot, legend, mc_stack, data, axis, signals)

    # Checks - Move on to next plot in case of failure
    if not mc_stack: return

    # Print yield table
    print "Yields for %s region"%reg.displayname
    print 30*'-'
    for yld in yield_tbl:
        if 'total' in yld.lower() or 'Data/MC' in yld:
            print 20*'-'
        print yld
    print 30*'-'

    # Draw the histograms
    draw_stack(axis, mc_stack, mc_errors, mc_total, signals, data, legend, reg.displayname)

    # Reset axis
    can.RedrawAxis()
    can.SetTickx()
    can.SetTicky()
    can.Update()

    # Save the histogram
    save_plot(can, plot.name+".eps")

    # Clean up
    root_delete([axis, mc_stack, mc_errors, mc_total, error_leg, signals, data,
                 data_hist, legend, hists])

def make_plotsRatio(plot, reg) :
    ''' '''
    ############################################################################
    # Intialize plot components
    # Canvases
    rcan = plot.pads
    rcan.canvas.cd()
    rcan.upper_pad.cd()

    if plot.doLogY : rcan.upper_pad.SetLogy(True)
    rcan.upper_pad.Update()

################################################################################
#Stack Plot Functions
################################################################################
def make_stack_legend():
    leg = pu.default_legend(xl=0.55,yl=0.71,xh=0.93,yh=0.90)
    leg.SetNColumns(2)
    leg_sig = pu.default_legend(xl=0.55, yl=0.6, xh=0.91, yh=0.71)
    leg_sig.SetNColumns(1)
    return leg

def make_stack_axis(plot, for_ratio=False):
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

    return hax

def add_stack_backgrounds(plot, leg, reg, yield_tbl):
    stack = r.THStack("stack_"+plot.name, "")

    # Initilize lists and defaults
    mc_samples = [s for s in SAMPLES if s.isMC]
    histos = []
    all_histos = []
    sig_histos = []
    avoid_bkg = []

    n_total_sm_yield = n_total_sm_error = 0.
    hists_to_clear = []

    # Make MC sample hists
    for mc_sample in mc_samples :
        # Initilize histogram
        h_name_tmp = re.sub(r'[(){}[\]]+','',plot.variable)
        h_name = "h_"+reg.name+'_'+mc_sample.name+"_"+h_name_tmp
        hists_to_clear.append(h_name)
        h = pu.th1d(h_name, "", int(plot.nbins),
                    plot.x_range_min, plot.x_range_max,
                    plot.x_label, plot.y_label)
        h.leg_name = mc_sample.displayname #dynamic classes...ooo yeah!

        h.SetLineColor(mc_sample.color)
        h.GetXaxis().SetLabelOffset(-999)
        if mc_sample.isSignal:
            h.SetLineWidth(2)
            h.SetLineStyle(2)
            h.SetFillStyle(0)
        else:
            h.SetFillColor(mc_sample.color)
            h.SetFillStyle(1001)
        h.Sumw2

        # Draw final histogram (i.e. selections and weights applied)
        cut = "(%s) * %s * %s"%(reg.tcut, mc_sample.weight_str, str(mc_sample.scale_factor))
        cut = r.TCut(cut)
        sel = r.TCut("1")
        draw_cmd = "%s>>+%s"%(plot.variable, h.GetName())
        mc_sample.tree.Draw(draw_cmd, cut * sel, "goff")


        # Yield +/- stat error
        stat_err = r.Double(0.0)
        integral = h.IntegralAndError(0,-1,stat_err)
        if not mc_sample.isSignal:
            n_total_sm_yield += float(integral)
            n_total_sm_error = (n_total_sm_error**2 + float(stat_err)**2)**0.5
        yield_tbl.append("%10s: %.2f +/- %.2f"%(mc_sample.name, integral, stat_err))

        # Add overflow
        if plot.add_overflow:
            pu.add_overflow_to_lastbin(h)
        if plot.add_underflow:
            pu.add_underflow_to_firstbin(h)

        # Record all and non-empty histograms
        if mc_sample.isSignal:
            leg_sig.AddEntry(h, mc_sample.displayname, "l")
            sig_histos.append(h)
        else:
            all_histos.append(h)
            histos.append(h) if integral > 0 else avoid_bkg.append(mc_sample.name)

    if not len(histos):
        print "ERROR (make_stack_background) :: All SM hists are empty. Skipping"
        return None, None, sig_histos, all_histos

    # Order the hists by total events
    histos = sorted(histos, key=lambda h: h.Integral())
    for h in histos :
        if "fake" in h.GetName() or "Fake" in h.GetName() : continue
        stack.Add(h)

    yield_tbl.append("%10s: %.2f +/- %.2f"%('Total SM', n_total_sm_yield, n_total_sm_error))

    h_leg = sorted(all_histos, key=lambda h: h.Integral(), reverse=True)
    histos_for_leg = histos_for_legend(h_leg)

    # draw the total bkg line
    hist_sm = stack.GetStack().Last().Clone("hist_sm")
    hist_sm.SetLineColor(r.kBlack)
    hist_sm.SetLineWidth(3)
    hist_sm.SetLineStyle(1)
    hist_sm.SetFillStyle(0)
    hist_sm.SetLineWidth(3)

    return stack, hist_sm, sig_histos, histos_for_leg 

def add_stack_data(plot, leg, reg, yield_tbl, sm_yield):
    #TODO: Look for a way to combine with backgrounds
    data = [s for s in SAMPLES if not s.isMC]
    assert len(data) <= 1, "ERROR :: Multiple data samples"
    if len(data) == 0: return None, None
    data = data[0]

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
    yield_tbl.append("%10s: %.2f"%("Data/MC", integral/sm_yield))
    #print "Data: %.2f +/- %.2f"%(integral, stat_err)

    # Add overflow
    if plot.add_overflow:
        pu.add_overflow_to_lastbin(hd)
    if plot.add_underflow:
        pu.add_underflow_to_firstbin(hd)

    gdata = pu.convert_errors_to_poisson(hd)
    #gdata.SetLineWidth(2)
    #uglify
    gdata.SetLineWidth(1)
    gdata.SetMarkerStyle(20)
    gdata.SetMarkerSize(1.5)
    gdata.SetLineColor(1)
    leg.AddEntry(gdata, "Data", "p")

    return gdata, hd

def add_stack_mc_errors(plot, leg, hists, stack):
    r.gStyle.SetHatchesSpacing(0.9)

    mcError = r.TH1F("mcError", "mcError", 2,0,2)
    mcError.SetLineWidth(3)
    mcError.SetFillStyle(3345)
    mcError.SetFillColor(r.kBlue)
    mcError.SetLineColor(r.kBlack)
    leg.AddEntry(mcError, "Standard Model", "fl")

    # now add backgrounds to legend
    for h in hists :
        leg.AddEntry(h, h.leg_name, "f")

    totalSM = stack.GetStack().Last().Clone("totalSM")
    nominalAsymErrors = pu.th1_to_tgraph(totalSM)
    totalSM.Delete()
    nominalAsymErrors.SetMarkerSize(0)
    nominalAsymErrors.SetLineWidth(0)
    nominalAsymErrors.SetFillStyle(3345)
    nominalAsymErrors.SetFillColor(r.kBlue)

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

    return mcError, nominalAsymErrors

def reformat_axis(plot, leg, stack, hd, hax, signals):
    # draw the MC stack and do cosmetics
    stack.SetMinimum(plot.y_range_min)

    # Determine y range for plot
    max_mult = 2.0 if len(signals) else 1.66
    if hd:
        maxy = max(hd.GetMaximum(), stack.GetMaximum())
    else:
        maxy = stack.GetMaximum()
    if not plot.doLogY:
        hax.SetMaximum(max_mult*maxy)
        stack.SetMaximum(max_mult*maxy)
    else:
        hax.SetMaximum(1e3*plot.y_range_max)
        stack.SetMaximum(1e3*plot.y_range_max)

def draw_stack(axis, mc_stack, mc_errors, mc_total, signals, data, legend, reg_name) :

    axis.Draw()
    mc_stack.Draw("HIST SAME")
    mc_errors.Draw("E2 same")
    mc_total.Draw('hist same')
    for hsig in signals: hsig.Draw("hist same")
    if data: data.Draw("option same pz 0")
    legend.Draw()
    pu.draw_atlas_label('Internal','Higgs LFV', reg_name)

def save_plot(can, outname):
    save_path = os.path.join(plots_dir, args.outdir, outname)
    save_path = os.path.normpath(save_path)
    can.SaveAs(save_path)
    print "%s saved to : %s"%(outname, os.path.dirname(save_path))

def root_delete(root_objects):
    for ro in root_objects:
        if type(ro) == list:
            root_delete(ro)
        elif issubclass(type(ro), r.TObject):
            ro.Delete()
        else:
            print "ERROR :: Unknown primitive type", type(ro)

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

    # print out the loaded samples and plots
    print " ============================"
    if SAMPLES :
        print "Loaded samples:    "
        for sample in SAMPLES :
            print '\t',
            sample.Print()
    if args.verbose:
        print "Loaded plots:"
        for plot in PLOTS :
            plot.Print()
    print " ============================"

def check_for_consistency() :
    '''
    Make sure that the plots are not asking for undefined region

    param:
        plots : list(plot class)
            plots defined in config file
        regions : list(region class)
            regions defined in config file
    '''
    region_names = [r.name for r in REGIONS]
    bad_regions = [p.region for p in PLOTS if p.region not in region_names]
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

        SAMPLES = conf.backgrounds + [conf.data]
        REGIONS = conf.regions
        PLOTS = conf.plots

        check_for_consistency()
        print_inputs(args)

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

