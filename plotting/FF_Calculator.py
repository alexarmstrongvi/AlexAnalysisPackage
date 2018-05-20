#!/usr/bin/env python
"""
================================================================================
TODO: Synopsis
Examples
    python <TODO:script name>
TODO: Description (for long scripts)

Author: 
    Alex Armstrong <alarmstr@cern.ch>
TODO: Licence:
    This script is in the public domain, free from copyrights or restrictions.
    Copyright: (C) <TODO:date>; University of California, Irvine 
================================================================================
"""

import sys, os, traceback, argparse
import time
import itertools
import array 
import ROOT as r
import global_variables as g
from collections import defaultdict

r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal cmd-line options
r.gROOT.SetBatch(True)
r.gStyle.SetOptStat(False)

import tools.plot_utils as pu
import tools.utils as utils
import tools.signal as signal
import tools.background as background
import tools.region as region
import tools.plot as plot

################################################################################
def main ():
    """ Main Function """
    
    global args
    global plots, regions, backgrounds, data
    ############################################################################
    # Prepare samples
    ############################################################################

    # Make histograms for FF calculations
    hists = get_FF_hists(data, backgrounds, regions, plots)

    # Substract contaminating MC
    hists = prepare_hists(hists)
    for ch_name, ch_hists in hists.iteritems():
        for h_name, hist in ch_hists.iteritems():
            print "INFO :: %s - %s - Integral = %.2f"%(ch_name, h_name, hist.Integral())
    
    # Divide to get FF
    FF_hists = make_FF_hists(hists)
    for ch_name, ch_hists in FF_hists.iteritems():
        for h_name, hist in ch_hists.iteritems():
            print "INFO :: %s - %s - Integral = %.2f"%(ch_name, h_name, hist.Integral())
            print "INFO :: Bin values : [",
            for xbin in range(hist.GetXaxis().GetNbins()):
                print "%.2f, "%hist.GetBinContent(xbin),
            print "]"

    # Save FF
    finalize_hists(hists, FF_hists)


################################################################################
# FUNCTIONS
def check_input_args():
    """ 
    Check that user inputs are as expected
    """
    if not os.path.isfile(args.config):
        print "ERROR :: configuration file not found: %s"%(args.config) 
        sys.exit()

def set_eventlists(data, backgrounds, reg):
    """
    Set/create the event lists for the input samples and their region cuts

    """
    cut = r.TCut(reg.tcut)
    sel = r.TCut("1")
    samples = backgrounds+[data]

    for sample in samples:
        list_name = "list_" + reg.name + "_" + sample.treename
        save_name = g.event_list_dir + list_name + ".root"
        save_name = os.path.normpath(save_name)
    
        # Reset event list
        sample.tree.SetEventList(0)
    
        # check if the list already exists
        if os.path.isfile(save_name) :
            rfile = r.TFile.Open(save_name)
            elist = rfile.Get(list_name)
            print "%10s : EventList found at %s"%(sample.name, save_name)
            sample.tree.SetEventList(elist)
        else :
            draw_list = ">> " + list_name
            sample.tree.Draw(draw_list, sel*cut)
            elist = r.gROOT.FindObject(list_name)
            sample.tree.SetEventList(elist)
            elist.SaveAs(save_name)

def get_FF_hists(data, backgrounds, regions, plots):
    hists = {}
    samples = backgrounds + [data]
    for reg in regions:
        print "INFO :: Making plots for %s region"%reg.displayname
        set_eventlists(data, backgrounds, reg)
        channel_name = get_hist_channel(reg.name)
        channel_hists = {}
        for sample, plot in itertools.product(samples, plots):
            if reg.name != plot.region: continue
            h_name = "h_"+reg.name+'_'+sample.treename+"_"+plot.variable
            print "INFO :: Makeing histogram:", h_name
            
            h = pu.th1d(h_name, "", int(plot.nbins), 
                        plot.x_range_min, plot.x_range_max, 
                        plot.x_label, plot.y_label)
            h.SetLineColor(sample.color)
            h.Sumw2
            if sample.name == "Data":
                cut = "(" + reg.tcut + ")"
            else:
                cut = "(" + reg.tcut + ") * eventweight * " + str(sample.scale_factor)
           
            cut = r.TCut(cut)
            sel = r.TCut("1")
            draw_cmd = "%s>>+%s"%(plot.variable, h.GetName())
            sample.tree.Draw(draw_cmd, cut * sel, "goff")
            if plot.add_overflow:
                pu.add_overflow_to_lastbin(h)
            channel_hists[h_name] = h.Clone()
            del h
        if channel_name in hists:
            hists[channel_name] = combine_dicts(hists[channel_name], channel_hists)
        else:
            hists[channel_name] = channel_hists
    print "INFO :: %d identified channels ->"%len(hists.keys()), hists.keys()
    return hists

def prepare_hists(hists):
    # Rebin histograms
    #newBinsEl = array.array('d',[0, 15, 20, 25, 35, 100])
    #newBinsMu = array.array('d',[0, 15, 20, 25, 35, 100])
    newBinsEl = array.array('d',[0, 5, 10, 15, 20, 23, 26, 29, 32, 35, 40, 45, 50, 100])
    newBinsMu = array.array('d',[0, 5, 10, 15, 20, 23, 26, 29, 32, 35, 40, 45, 50, 100])
    samples_to_subtract = ['WZ','ZZ','ttbar']
    
    rebin_hists = defaultdict(dict)
    bkgd_subtracted = defaultdict(dict)
    for ch_name, ch_hists in hists.iteritems():
        print "INFO :: Modifying hists in channel", ch_name
        non_fake_bkgd_hists = [[],[]]
        for h_name, hist in ch_hists.iteritems():
            # Rebin
            ch_flav = get_channel_flav(ch_name)
            if ch_flav == 'e':
                rebin_hists[ch_name][h_name] = hist.Rebin(len(newBinsEl)-1, h_name, newBinsEl)
            else:
                rebin_hists[ch_name][h_name] = hist.Rebin(len(newBinsMu)-1, h_name, newBinsMu)


            # Capture relevent hists
            if 'data' in h_name:
                h_tmp = rebin_hists[ch_name][h_name].Clone(hist.GetName()+"_minus_bkgd")
                bkgd_subtracted[ch_name]['%s_minus_bkgd'%h_name] = h_tmp 
                print "INFO :: Data yield for %s before background subtraction: %.2f"%(h_name, h_tmp.Integral())
            if any(x in h_name for x in samples_to_subtract):
                if hist_is_num(h_name):
                    non_fake_bkgd_hists[0].append(rebin_hists[ch_name][h_name])
                else:
                    non_fake_bkgd_hists[1].append(rebin_hists[ch_name][h_name])

        for h_name, hist in bkgd_subtracted[ch_name].iteritems():
            if 'minus_bkgd' in h_name:
                num_or_den = 0 if hist_is_num(h_name) else 1
                for bkgd in non_fake_bkgd_hists[num_or_den]:
                    hist.Add(bkgd, -1)
                print "INFO :: Data yield after background subtraction: ", hist.Integral()
    hists = combine_dicts(rebin_hists, bkgd_subtracted)
    return hists

def make_FF_hists(hists):
    hists_to_divide = ['data','zee', 'zmumu']
    FF_hists = {}
    for ch_name, ch_hists in hists.iteritems():
        FF_hists_ch = {}
        for h_name, hist in ch_hists.iteritems():
            if not any(x in h_name for x in hists_to_divide): continue
            if hist_is_num(h_name):  
                ff_name = h_name.replace("num","")
                FF_hists_ch[ff_name] = hist.Clone(ff_name)
        for h_name, hist in ch_hists.iteritems():
            if not any(x in h_name for x in hists_to_divide): continue
            if not hist_is_num(h_name):  
                ff_name = h_name.replace('den','')
                FF_hists_ch[ff_name].Divide(hist)
        FF_hists[ch_name] = FF_hists_ch
    return FF_hists

def finalize_hists(hists, ff_hists):
    can = r.TCanvas("Fake Factor")
    can.cd()
    can.SetFrameFillColor(0)
    can.SetFillColor(0)
    can.SetLeftMargin(0.14)
    can.SetRightMargin(0.05)
    can.SetBottomMargin(1.3*can.GetBottomMargin())
    can.Update()
    
    ofile = r.TFile(args.ofile_name,'RECREATE')
    all_hists = combine_dicts(hists, ff_hists)
    for ch_name, ch_hists in all_hists.iteritems():
        for h_name, hist in ch_hists.iteritems():
            hist.Draw("pE1")
            hist.Write()
            if not ('num' in h_name or 'den' in h_name):
                can.SaveAs('%s%s.eps'%(g.plots_dir,h_name))
    ofile.Write()
    ofile.Close()

def hist_is_num(hist_name):
    return 'num' in hist_name
def get_hist_channel(hist_name):
    return hist_name.split('_')[-1]
def get_channel_flav(ch_name):
    return ch_name[-1]
def combine_dicts(d1, d2):
    dic = d1.copy()
    if not d1 and d2: return d2
    elif d1 and not d2: return d1

    for x, y in zip(d1, d2):
        t1, t2 = type(d1[x]), type(d2[y])
        if t1 == t2 and t1 == dict:
            dic[x] = combine_dicts(d1[x], d2[y])
        elif t1 != dict and t2 != dict:
            if any(a == b for a, b in zip(d1, d2)):
                print "WARNING :: Overwriting dictionary entries from:"
                print d1
                print d2
            dic.update(d2)
    return dic



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
        parser.add_argument('-v', '--verbose', 
                            action='store_true', default=False, 
                            help='verbose output')
        args = parser.parse_args()
        
        check_input_args()
        plots = []
        data = None
        backgrounds = []
        regions = []
        execfile(args.config)
        region_names = [p.region for p in plots]
        regions = [reg for reg in regions if reg.name in region_names]

        assert (data and backgrounds and regions and plots), (
            "ERROR :: configuration file not loaded correctly"
        )
        
        if args.verbose: 
            print '>'*40
            print 'Running {}...'.format(os.path.basename(__file__))
            print time.asctime()
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

