#!/usr/bin/env python
import pdb
"""
================================================================================
Make fake factor files and plots
Examples
    make_fake_factor_plots.py -c config_file.conf

Author:
    Alex Armstrong <alarmstr@cern.ch>
Licence:
    Copyright: (C) <June 1st, 2018>; University of California, Irvine
================================================================================
"""

# General python
import sys, os, traceback, argparse
import time
from re import sub
from copy import copy
from collections import defaultdict
from importlib import import_module
import subprocess
from contextlib import contextmanager


# Root data analysis framework
import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal cmd-line options
r.gROOT.SetBatch(True)
r.gStyle.SetOptStat(False)

# Turn off root ownership when creating classes
r.TEventList.__init__._creates = False
r.TH1F.__init__._creates = False
r.TGraphErrors.__init__._creates = False
r.TGraphAsymmErrors.__init__._creates = False
r.TLegend.__init__._creates = False
r.THStack.__init__._creates = False

# Local classes for plotting
import tools.plot_utils as pu
from tools.plot import Types
from tools.YieldTable import UncFloat
from global_variables import event_list_dir, plots_dir
from plotter2 import make_stack_axis

@contextmanager
def open_root(f_name, f_mode):
    ofile = r.TFile(f_name, f_mode)
    try:
        yield ofile
    finally:
        ofile.Write()
        ofile.Close()

class KeyManager(object) :
    bkg_mc_str = 'truth_mc_3ID_lep'
    fake_mc_str = 'truth_mc_2ID_1aID_lep'

    def __init__(self):
        # TODO: move globals into class
        self._mc_stack = None
        self.mc_stack_hists = {}
        self.mc_hist = None
        self._mc_truth_bkg_stack = None
        self.mc_truth_bkg_stack_hists = {}
        self.mc_truth_bkg_hist = None
        self._mc_truth_fake_stack = None
        self.mc_truth_fake_stack_hists = {}
        self.mc_truth_fake_hist = None
        self.data_hist = None
        self.data_corr_hist = None
        self.fake_factor_keys = {}

    ############################################################################
    @property
    def data_fake_factor(self):
        return self.fake_factor_keys[self.data_hist]

    @property
    def data_corr_fake_factor(self):
        return self.fake_factor_keys[self.data_corr_hist]

    @property
    def mc_fake_factor(self):
        return self.fake_factor_keys[self.mc_truth_fake_hist]

    # TODO: Remove property after testing that there is no error
    @property
    def mc_stack(self):
        return self._mc_stack
    @mc_stack.setter
    def mc_stack(self, new_key_name):
        old_key_name = self._mc_stack
        assert not old_key_name or old_key_name == new_key_name, (
            "ERROR :: stack name generator not consistent", old_key_name, new_key_name)
        self._mc_stack = new_key_name

    @property
    def mc_truth_bkg_stack(self):
        return self._mc_truth_bkg_stack
    @mc_truth_bkg_stack.setter
    def mc_truth_bkg_stack(self, new_key_name):
        old_key_name = self._mc_truth_bkg_stack
        assert not old_key_name or old_key_name == new_key_name, (
            "ERROR :: stack name generator not consistent", old_key_name, new_key_name)
        self._mc_truth_bkg_stack = new_key_name

    @property
    def mc_truth_fake_stack(self):
        return self._mc_truth_fake_stack
    @mc_truth_fake_stack.setter
    def mc_truth_fake_stack(self, new_key_name):
        old_key_name = self._mc_truth_fake_stack
        assert not old_key_name or old_key_name == new_key_name, (
            "ERROR :: stack name generator not consistent", old_key_name, new_key_name)
        self._mc_truth_fake_stack = new_key_name

    ############################################################################
    def generate_hist_key(self, sample, region, cut):
        sample_name = remove_num_den(sample.name)
        if region.truth_fake_sel in cut:
            key = self.fake_mc_str + "_" + sample_name
            self.mc_truth_fake_stack_hists[sample_name] = key
        elif region.truth_bkg_sel in cut:
            key = self.bkg_mc_str + "_" + sample_name
            self.mc_truth_bkg_stack_hists[sample_name] = key
        elif sample.isMC:
            key = sample_name
            self.mc_stack_hists[sample_name] = key
        else:
            key = sample_name
            self.data_hist = key
        return key

    def generate_stack_key(self, stack_hist_key):
        if self.is_bkg_mc(stack_hist_key):
            stack_key = "stack_%s"%self.bkg_mc_str
            self.mc_truth_bkg_stack = stack_key
        elif self.is_fake_mc(stack_hist_key):
            stack_key = "stack_%s"%self.fake_mc_str
            self.mc_truth_fake_stack = stack_key
        else:
            stack_key = "stack_mc"
            self.mc_stack = stack_key
        return stack_key

    def generate_mc_hist_key(self, hist_key):
        if self.is_bkg_mc(hist_key):
            hist_key = "total_%s_hist"%self.bkg_mc_str
            self.mc_truth_bkg_hist = hist_key
        elif self.is_fake_mc(hist_key):
            hist_key = "total_%s_hist"%self.fake_mc_str
            self.mc_truth_fake_hist = hist_key
        else:
            hist_key = "total_mc_hist"
            self.mc_hist = hist_key
        return hist_key

    def generate_data_corr_key(self):
        data_corr_key = "%s_bkgd_subtracted"%self.data_hist
        self.data_corr_hist = data_corr_key
        return data_corr_key

    def generate_fake_factor_key(self, hist_key):
        ff_key = "fake_factor_" + hist_key
        self.fake_factor_keys[hist_key] = ff_key
        return ff_key

    ############################################################################
    def get_data_keys(self):
        keys = [self.data_hist, self.data_corr_hist]
        keys = [k for k in keys if k]
        return keys

    def get_stack_keys(self):
        keys = [self.mc_stack, self.mc_truth_bkg_stack, self.mc_truth_fake_stack]
        keys = [k for k in keys if k]
        return keys

    def get_total_mc_hist_keys(self):
        keys = [self.mc_hist, self.mc_truth_bkg_hist, self.mc_truth_fake_hist]
        keys = [k for k in keys if k]
        return keys

    def get_truth_bkg_keys(self):
        keys = self.mc_truth_bkg_stack_hists.values()
        keys += [self.mc_truth_bkg_stack, self.mc_truth_bkg_hist]
        keys = [k for k in keys if k]
        return keys

    def get_truth_fake_keys(self):
        keys = self.mc_truth_fake_stack_hists.values()
        keys += [self.mc_truth_fake_stack, self.mc_truth_fake_hist]
        keys = [k for k in keys if k]
        return keys

    def get_fake_factor_input_keys(self):
        keys = self.get_data_keys() + [self.mc_truth_fake_hist] 
        keys = [k for k in keys if k]
        return keys


    ############################################################################
    def is_data(self, hist_key):
        return hist_key in self.get_data_keys()

    def is_stack(self, hist_key):
        return hist_key in self.get_stack_keys()

    def is_fake_mc(self, hist_key):
        return hist_key in self.get_truth_fake_keys()

    def is_bkg_mc(self, hist_key):
        return hist_key in self.get_truth_bkg_keys()

    def is_fake_factor(self, hist_key):
        return hist_key in self.fake_factor_keys


################################################################################
def main ():
    """ Main Function """

    global args

    # Create hist container
    # Organization : dict["chnanel"]["sample"]["num/den"] = TH1D
    hists = defaultdict(lambda: defaultdict(lambda: defaultdict(r.TH1D)))
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
        # Make hists
        for plot in plots_with_region:
            print "Running on %s plot"%plot.variable
            add_ff_hist_primitives(plot, hists, reg)
    print "Identified channels:", hists.keys()
    print "Samples in each channel:",sorted(hists[hists.keys()[0]][conf.DEN_STR].keys(), key=len)
    print "\nFormatting histograms"
    format_and_combine_hists(hists)
    print "Making fake factor histograms"
    ff_hists = get_fake_factor_hists(hists)
    print "Making plots"
    save_and_write_hists(ff_hists, hists)
    print 30*'-', 'PLOTS COMPLETED', 30*'-','\n'
################################################################################
# Fake factor functions
def add_ff_hist_primitives(plot, hists, reg):
    channel_name = remove_num_den(reg.name)
    num_or_den = get_fake_channel(reg.name)
    samples = [s for s in SAMPLES if num_or_den in s.name]

    for sample in samples:
        # Get histogram(s) for each plot-sample pair
        # For MC samples, make additional hists for truth matched selections

        # Get cuts
        cuts = []
        if sample.isMC:
            base_cut = "((%s) * %s * %f)"%(reg.tcut, sample.weight_str, sample.scale_factor)
            weight = "%s * %f"%(sample.weight_str, sample.scale_factor)
            cuts.append("(%s) * %s"%(reg.tcut, weight))
            cuts.append("(%s && %s) * %s"%(reg.tcut, reg.truth_fake_sel, weight))
            cuts.append("(%s && %s) * %s"%(reg.tcut, reg.truth_bkg_sel, weight))
        else: # Data
            cuts.append("(" + reg.tcut + ")")

        # Apply each cut to the sample and fill hist
        for cut in cuts:
            hist_key = KEYS.generate_hist_key(sample, reg, cut)
            # histogram name must be unique relative to all hists made by script
            h_name = "h_"+reg.name+'_'+hist_key+"_"+plot.variable
            hist = build_hist(h_name, plot, sample, cut)
            print "Base histograms (Sample: %s): %s created"%(sample.name, h_name)
            hist.displayname = sample.displayname
            hist.plot = plot
            hist.color = sample.color
            hists[channel_name][num_or_den][hist_key] = hist

def format_and_combine_hists(hists):
    '''
    Make all the combined histograms.

    1) THStack plots for all sets of MC samples
    2) Total SM plots for all sets of MC samples
    3) Data with non-fake background subtracted

    args:
        hists (dict[dict[dict[TH1D]]] : histograms organized by channel,
        sample or combined samples, and then by numerator or denominator
        selections

    '''
    # First loop to make MC stack histograms
    for channel_name, ch_dict in hists.items():
        for num_or_den, sample_dict in ch_dict.items():
            for hist_key, hist in sample_dict.items():
                #TODO: Loop over keys from key manager instead?
                if KEYS.is_data(hist_key): continue
                stack_key = KEYS.generate_stack_key(hist_key)
                if stack_key not in sample_dict:
                    print "Stack histograms (Channel: %s [%s]): %s initialized"%(
                            channel_name, num_or_den, stack_key)
                    stack_name = stack_key+'_'+num_or_den
                    sample_dict[stack_key] = r.THStack(stack_name, "")
                    sample_dict[stack_key].plot = hist.plot

                sample_dict[stack_key].Add(hist)

        # Second loop to make MC total histograms
        for num_or_den, sample_dict in ch_dict.items():
            for hist_key, hist in sample_dict.items():
                #TODO: Use KeyManager to grab correct hist
                if not KEYS.is_stack(hist_key): continue
                mc_hist_key = KEYS.generate_mc_hist_key(hist_key)
                hist_name = mc_hist_key+'_'+num_or_den
                mc_total_hist = hist.GetStack().Last().Clone(hist_name)
                mc_total_hist.plot = hist.plot
                if KEYS.is_bkg_mc(hist_key):
                    mc_total_hist.displayname = KEYS.bkg_mc_str.replace("_"," ").title()
                elif KEYS.is_fake_mc(hist_key):
                    mc_total_hist.displayname = KEYS.fake_mc_str.replace("_"," ").title()
                else:
                    mc_total_hist.displayname = 'Total MC'
                mc_total_hist.SetLineColor(r.kBlack)
                mc_total_hist.SetLineWidth(3)
                mc_total_hist.SetLineStyle(1)
                mc_total_hist.SetFillStyle(0)
                mc_total_hist.SetLineWidth(3)
                mc_total_hist.is_total = True

                sample_dict[mc_hist_key] = mc_total_hist
                print "Total MC histogram  (Channel: %s [%s]): %s created"%(
                       channel_name, num_or_den, mc_hist_key)

            # Grab hists
            data_hist = sample_dict[KEYS.data_hist]
            mc_background_hist = sample_dict[KEYS.mc_truth_bkg_hist]
            data_corr_hist_key = KEYS.generate_data_corr_key()

            # Create and store background-subtracted data histogram
            data_corrected_name = "%s_%s"%(data_corr_hist_key, num_or_den)
            data_corrected_hist = data_hist.Clone(data_corrected_name)
            data_corrected_hist.Add(mc_background_hist, -1)
            data_corrected_hist.displayname = "Data (bkgd subtracted)"
            data_corrected_hist.plot = data_hist.plot
            sample_dict[data_corr_hist_key] = data_corrected_hist
            print "Data corrected histogram (Channel: %s [%s]): %s created"%(
                    channel_name, num_or_den, data_corr_hist_key)

def get_fake_factor_hists(hists):
    ff_hists = defaultdict(lambda: defaultdict(r.TH1D))
    for channel_name, ch_dict in hists.iteritems():
        for hist_key in KEYS.get_fake_factor_input_keys():
            fake_factor_key = KEYS.generate_fake_factor_key(hist_key)
            fake_factor_name = channel_name + "_" + fake_factor_key
            num_hist = ch_dict[conf.NUM_STR][hist_key]
            ff_hist = num_hist.Clone(fake_factor_name)
            ff_hist.Divide(ch_dict[conf.DEN_STR][hist_key])

            # Append some information   
            ff_hist.displayname = hist_key.replace("_"," ")
            ff_hist.plot = num_hist.plot
            ff_hist.plot = copy(num_hist.plot)

            # Format the hists
            ff_hist.plot.update(doLogY = False, doNorm = True)

            ff_hists[channel_name][fake_factor_key] = ff_hist
    return ff_hists

def save_and_write_hists(ff_hists_dict, hists):
# Writing fake factor hists to root file
#with open_root(args.ofile_name,"RECREATE") as ofile:
#    for channel_name, ff_hists in ff_hists_dict.iteritems():
#        ff_hists[KEYS.data_corr_fake_factor].Write()

# Saving plots of fake factor hists
    for channel_name, ff_hists in ff_hists_dict.iteritems():
        data_corr_ff_hist = ff_hists[KEYS.data_corr_fake_factor]
        data_ff_hist = ff_hists[KEYS.data_fake_factor]
        mc_ff_hist = ff_hists[KEYS.mc_fake_factor]
        plot_title = 'Fake Factor (Ch: %s)'%channel_name
        hists_to_plot = [mc_ff_hist, data_ff_hist, data_corr_ff_hist]
        plot = data_corr_ff_hist.plot
        save_hist(plot_title, plot, channel_name, hists_to_plot)

    # Save all other desired plots
    for channel_name, ch_dict in hists.iteritems():
        for num_or_den, sample_dict in ch_dict.iteritems():
            # Save MC Stacks
            mc_stack = sample_dict[KEYS.mc_stack]
            mc_hist = sample_dict[KEYS.mc_hist]
            plot_title = 'MC Backgrounds'
            plot_title += ' (%s)'%num_or_den
            hists_to_plot = [mc_stack, mc_hist]
            plot = mc_stack.plot
            save_hist(plot_title, plot, channel_name, hists_to_plot)

            mc_truth_bkg_stack = sample_dict[KEYS.mc_truth_bkg_stack]
            mc_truth_bkg_hist = sample_dict[KEYS.mc_truth_bkg_hist]
            plot_title = 'MC Backgrounds with 3 ID truth-matched prompt leptons'
            plot_title += ' (%s)'%num_or_den
            hists_to_plot = [mc_truth_bkg_stack, mc_truth_bkg_hist]
            plot = mc_truth_bkg_stack.plot
            save_hist(plot_title, plot, channel_name, hists_to_plot)

            mc_truth_fake_stack = sample_dict[KEYS.mc_truth_fake_stack]
            mc_truth_fake_hist = sample_dict[KEYS.mc_truth_fake_hist]
            plot_title = 'MC Backgrounds with anti-ID truth-matched fake lepton'
            plot_title += ' (%s)'%num_or_den
            hists_to_plot = [mc_truth_fake_stack, mc_truth_fake_hist]
            plot = mc_truth_fake_stack.plot
            save_hist(plot_title, plot, channel_name, hists_to_plot)

            # Overlay of data before and after correction with MC truth
            # background stack
            data_hist = sample_dict[KEYS.data_hist]
            data_corr_hist = sample_dict[KEYS.data_corr_hist]
            plot_title = 'Data before and after background substraction'
            plot_title += ' (%s)'%num_or_den
            hists_to_plot = [mc_truth_bkg_hist, data_hist, data_corr_hist]
            plot = data_hist.plot
            save_hist(plot_title, plot, channel_name, hists_to_plot)

            # Overlay of fake MC with data after MC truth background
            # subtractio
            plot_title = 'Resulting fake estimates in Data and MC'
            plot_title += ' (%s)'%num_or_den
            hists_to_plot = [mc_truth_fake_stack, data_corr_hist]
            plot = mc_truth_fake_stack.plot
            save_hist(plot_title, plot, channel_name, hists_to_plot)

            # Overlay of data with stack of total MC truth background and total
            # MC truth fake
            stack = r.THStack("bkg_and_fake_mc","")
            stack.Add(mc_truth_bkg_hist)
            stack.Add(mc_truth_fake_hist)
            plot_title = 'MC breakdown of fake and non-fake backgrounds'
            plot_title += ' (%s)'%num_or_den
            hists_to_plot = [stack, data_hist]
            plot = data_hist.plot
            save_hist(plot_title, plot, channel_name, hists_to_plot)

MARKER_COLORS = [r.kBlack, r.kBlue, r.kGreen+3, r.kRed]

def get_marker_color():
    for color in MARKER_COLORS:
        yield color

def save_hist(title, plot, reg_name, hist_list):
    can = plot.pads.canvas
    can.cd()
    can.SetTitle(title)
    if plot.doLogY : can.SetLogy(True)

    # Make Axis
    axis = make_stack_axis(plot)

    # Format Primitives
    # TODO: Sort THStack by integral

    # Format primitives and fill legend
    legend = pu.default_legend(xl=0.55,yl=0.71,xh=0.93,yh=0.90)
    legend.SetNColumns(2)
   
    color_counter = 0
    for hist in hist_list:
        if isinstance(hist, r.THStack):
            for stack_hist in hist.GetHists():
                if hasattr(stack_hist, 'color'):
                    color = stack_hist.color
                else:
                    color = MARKER_COLORS[color_counter]
                    color_counter = (color_counter + 1)%len(MARKER_COLORS)
                stack_hist.SetFillColor(color)
                stack_hist.SetLineColor(color)
                stack_hist.SetFillStyle(1001)
                legend.AddEntry(stack_hist, stack_hist.displayname, "f")
            hist.Modified()
        else:
            hist.SetLineWidth(1)
            hist.SetLineColor(r.kBlack)
            hist.SetFillStyle(0)
            hist.SetMarkerStyle(r.kFullCircle)
            hist.SetMarkerSize(1.5)
            hist.SetMarkerColor(MARKER_COLORS[color_counter])
            color_counter = (color_counter + 1)%len(MARKER_COLORS)
            legend.AddEntry(hist, hist.displayname, "p")

    # Draw primitives to canvas
    axis.Draw()
    for hist in hist_list:
        if isinstance(hist, r.THStack):
            hist.Draw("HIST SAME")
        elif hasattr(hist, 'is_total'):
            hist.Draw("HIST SAME")
        else:
            hist.Draw("pE1 same")

    legend.Draw()
    pu.draw_atlas_label('Internal','Higgs LFV', reg_name)

    # Finalize
    can.RedrawAxis()
    can.SetTickx()
    can.SetTicky()
    can.Update()

    # Save
    outname = title + ".eps"
    outname = outname.replace(" ","_")
    outname = sub(r'[:\-(){}[\]]+','', outname)
    save_path = os.path.join(plots_dir, args.dir_name, outname)
    save_path = os.path.normpath(save_path)
    can.SaveAs(save_path)
    axis.Delete()


################################################################################
# Fake factor functions
def build_hist(h_name, plot, sample, cut):
    hist = pu.th1d(h_name, "", int(plot.nbins),
                plot.xmin, plot.xmax,
                plot.ylabel, plot.ylabel)
    hist.Sumw2
    hist.SetLineColor(sample.color)
    cut = r.TCut(cut)
    draw_cmd = "%s>>+%s"%(plot.variable, hist.GetName())
    sample.tree.Draw(draw_cmd, cut, "goff")
    print "%s: Yield -> %.2f"%(sample.name, hist.Integral())
    if plot.add_overflow:
        pu.add_overflow_to_lastbin(hist)
    if plot.add_underflow:
        pu.add_underflow_to_firstbin(hist)
    return hist

def remove_num_den(name):
    name = name.replace(conf.NUM_STR,"")
    name = name.replace(conf.DEN_STR,"")
    name = name.replace("__","_")
    if name.endswith("_"): name = name[:-1]
    if name.startswith("_"): name = name[1:]
    return name


def get_fake_channel(reg_name):
    return conf.DEN_STR if conf.DEN_STR in reg_name else conf.NUM_STR

################################################################################
# Check functions
def check_environment():
    """ Check if the shell environment is setup as expected """
    assert os.environ['USER'], "USER variable not set"

    python_ver = sys.version_info[0] + 0.1*sys.version_info[1]
    assert python_ver >= 2.7, ("Running old version of python\n", sys.version)

def check_input_args():
    """
    Check that user inputs are as expected
    """
    if not os.path.isfile(args.config):
        print "ERROR :: configuration file not found: %s"%(args.config)
        sys.exit()

    if not os.path.exists(args.dir_name):
        print "ERROR :: output directory not found: %s"%(args.dir_name)
        sys.exit()

    of = os.path.join(args.dir_name, args.ofile_name)
    if os.path.exists(of):
        if not os.path.exists("%s.bu"%of):
            print "Renaming old output file %s -> %s.bu"%(of, of)
            mv_cmd = 'mv %s %s.bu'%(of, of)
            subprocess.call(mv_cmd, shell=True)
        else:
            print "WARNING :: Output file already exists: %s"%of
            print "\tConsider deleting it or its backup (%s.bu)"%of
            sys.exit()

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

def check_import_globals():

    #TODO: Allow plotting of multiple variables
    assert len(set([p.variable for p in conf.PLOTS])) == 1, (
        "ERROR :: Currently cannot handle plotting multiple variables")

    reg_names = [p.region for p in conf.PLOTS]
    print "TESTING :: plot regions =", reg_names
    assert all(conf.NUM_STR in n or conf.DEN_STR in n for n in reg_names),(
        "ERROR :: Non fake factor region defined")
    
    num_regions = [n for n in reg_names if conf.NUM_STR in n]
    den_regions = [n for n in reg_names if conf.DEN_STR in n]
    assert len(num_regions) == len(den_regions), (
        "ERROR :: Number of numerator and denominator regions do not match")
    
    num_regions = [n.replace(conf.NUM_STR, "") for n in num_regions]
    den_regions = [n.replace(conf.DEN_STR, "") for n in den_regions]
    assert sorted(num_regions) == sorted(den_regions), (
        "ERROR :: numerator and denominator regions do not matched")


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
    print "     plot config      :  %s "%args.config
    print "     output file      :  %s "%args.ofile_name
    print "     output directory :  %s "%args.dir_name
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

################################################################################
# Run main when not imported
if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = argparse.ArgumentParser(
                description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument("-c", "--config",
                            default="",
                            help='path to config file')
        parser.add_argument('-o', '--ofile_name',
                            default="FF_hists.root",
                            help="Output file name")
        parser.add_argument('-d', '--dir_name',
                            default="./",
                            help="Output directory")
        parser.add_argument('-v', '--verbose',
                            action='store_true', default=False,
                            help='verbose output')
        args = parser.parse_args()

        check_environment()
        check_input_args()

        import_conf = args.config.replace(".py","")
        conf = import_module(import_conf)
        check_import_globals()

        SAMPLES = conf.SAMPLES
        REGIONS = conf.REGIONS
        PLOTS = conf.PLOTS
        YIELD_TBL = conf.YIELD_TBL
        YIELD_TABLES = []
        KEYS = KeyManager()

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


