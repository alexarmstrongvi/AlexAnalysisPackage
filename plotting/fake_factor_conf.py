"""
================================================================================
An example configuration file for plotter.py

Author:
    Alex Armstrong <alarmstr@cern.ch>
    ... with lots of ideas borrowed from Danny Antrim <dantrim@cern.ch>

License:
    Copyright: (C) <May 20th, 2018>; University of California, Irvine
================================================================================
"""

# General python
import sys, os
import re
from collections import namedtuple
from copy import deepcopy

# Root data analysis framework
import ROOT
# Prevent root from printing anything to the screen
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# Local classes for plotting
from PlotTools.plot import Plot1D, Plot2D, Types
from PlotTools.sample import Sample, MCsample, Data, Background, Signal
from PlotTools.region import Region
from PlotTools.YieldTable import YieldTable, UncFloat


################################################################################
# Useful globals
################################################################################
import global_variables as g
# Path to directories with flat ntuples
bkg_ntuple_dir    = g.mc_ntuples
signal_ntuple_dir = g.mc_ntuples
data_ntuple_dir   = g.data_ntuples
fake_ntuple_dir   = g.data_ntuples

# Luminosity options
# Description | 2015-16 (ipb) |  2015-16 (ifb) | 2015 (ipb) | 2016 (ipb) | dummy
# Value       |     36180     |      36.01     |    320     |    32971   |   1
lumi = 36180

# Strings for plotting
MCsample.weight_str = 'eventweight'
MCsample.scale_factor = lumi
Sample.input_file_treename = 'superNt'
Plot1D.auto_set_ylimits = False
Plot1D.doLogY = False
Plot1D.ymax = 6e3
NUM_STR = "num"
DEN_STR = "den"

################################################################################
# Samples
# - Create all the samples that will be put into plots (e.g. backgrounds,
#   signal, data)
################################################################################
SAMPLES = []

for num_den in [NUM_STR, DEN_STR]:
    # Data
    data = Data('data_%s'%num_den)
    SAMPLES.append(data)

    #######################################
    # Initialize all backgrounds
    # Higgs (Htt + HWW + HV + ttH)
    higgs = Background("higgs_%s"%num_den, "Higgs")
    higgs.color = ROOT.kRed - 7
    SAMPLES.append(higgs)

    # ttbar
    ttbar = Background("ttbar_%s"%num_den, "t#bar{t}")
    ttbar.color = ROOT.kOrange+2
    SAMPLES.append(ttbar)

    # singletop
    stop = Background("st_%s"%num_den, "Single top")
    stop.color = ROOT.kOrange+1
    SAMPLES.append(stop)

    # W+top
    wtop = Background("wt_%s"%num_den, "Wt")
    wtop.color = ROOT.kOrange+8
    SAMPLES.append(wtop)

    # ttbar + X
    ttbar_x = Background("ttbar_x_%s"%num_den, "t#bar{t} + X")
    ttbar_x.color = ROOT.kOrange+5
    SAMPLES.append(ttbar_x)

    # Top (ttbar + singletop + Wt)
    top = Background("top_%s"%num_den, "Top")
    top.color = ROOT.kYellow+1
    SAMPLES.append(top)
    
    # VV combined
    VV = Background("VV_%s"%num_den, "VV")
    VV.color = ROOT.kSpring-7
    VV.scale_factor *= 1.05
    SAMPLES.append(VV)

    # VVV combined
    VVV = Background("VVV_%s"%num_den, "VVV")
    VVV.color = ROOT.kSpring-7
    SAMPLES.append(VVV)

    # Zll
    zll = Background("zll_%s"%num_den, "Zll")
    zll.color = ROOT.kAzure-7
    SAMPLES.append(zll)

    # Zee
    zee = Background("zee_%s"%num_den, "Zee")
    zee.color = ROOT.kAzure-7
    SAMPLES.append(zee)

    # Zmumu
    zmumu = Background("zmumu_%s"%num_den, "Zmumu")
    zmumu.color = ROOT.kAzure-9
    SAMPLES.append(zmumu)

    # Ztt
    ztt = Background("ztt_%s"%num_den, "Z#tau#tau")
    ztt.color = ROOT.kAzure+8
    SAMPLES.append(ztt)

    # Wjets
    wjets = Background("wjets_%s"%num_den, "W+jets")
    wjets.color = ROOT.kOrange + 2
    SAMPLES.append(wjets)

    # W+gamma
    wgamma = Background("wgamma_%s"%num_den, "W+gamma")
    wgamma.color = ROOT.kOrange-1
    SAMPLES.append(wgamma)

    # Z+gamma
    zgamma = Background("zgamma_%s"%num_den, "Z+gamma")
    zgamma.color = ROOT.kOrange-3
    SAMPLES.append(zgamma)

    #######################################
    ## Build the TChain/TTree for each sample
    # To remove sample from plot, comment out the line setting its TChain
    # Samples with empty TChains get removed below
    data_dsids = g.groups['data15']+g.groups['data16']
    data_ntuple_dir2 = data_ntuple_dir[:-1] + '_' + num_den + "/"
    bkg_ntuple_dir2 = bkg_ntuple_dir[:-1] + '_' + num_den + "/"

    data.set_chain_from_dsid_list(data_dsids, data_ntuple_dir2)
    higgs.set_chain_from_dsid_list(g.groups['higgs'], bkg_ntuple_dir2)
    top.set_chain_from_dsid_list(g.groups['top'], bkg_ntuple_dir2)
    #ttbar.set_chain_from_dsid_list(g.groups['ttbar'], bkg_ntuple_dir2)
    #stop.set_chain_from_dsid_list(g.groups['singletop'], bkg_ntuple_dir2)
    #wtop.set_chain_from_dsid_list(g.groups['Wt'], bkg_ntuple_dir2)
    ttbar_x.set_chain_from_dsid_list(g.groups['ttbar_lep'], bkg_ntuple_dir2)
    VV.set_chain_from_dsid_list(g.groups['VV'], bkg_ntuple_dir2)
    VVV.set_chain_from_dsid_list(g.groups['VVV'], bkg_ntuple_dir2)
    zll.set_chain_from_dsid_list(g.groups['zll'], bkg_ntuple_dir2)
    #zee.set_chain_from_dsid_list(g.groups['zee'], bkg_ntuple_dir2)
    #zmumu.set_chain_from_dsid_list(g.groups['zmumu'], bkg_ntuple_dir2)
    ztt.set_chain_from_dsid_list(g.groups['ztt'], bkg_ntuple_dir2)
    wjets.set_chain_from_dsid_list(g.groups['wjets'], bkg_ntuple_dir2)
    wgamma.set_chain_from_dsid_list(g.groups['wgamma'], bkg_ntuple_dir2)
    zgamma.set_chain_from_dsid_list(g.groups['zgamma'], bkg_ntuple_dir2)

SAMPLES = [s for s in SAMPLES if s.is_setup()]
assert SAMPLES, "ERROR :: No samples are setup"
################################################################################
# Regions
# - Create all the regions that will be applied to samples
################################################################################

########################################
# Define important selections

# General
emu = 'dilep_flav == 0'
mue = 'dilep_flav == 1'
ee = 'dilep_flav == 2'
mumu = 'dilep_flav == 3'
DF_OS = 'dilep_flav <= 1 && l_q[0]*l_q[1]<0'
SF_OS = 'dilep_flav > 1 && l_q[0]*l_q[1]<0'
Z_emu = 'Z_dilep_flav == 0'
Z_mue = 'Z_dilep_flav == 1'
Z_ee = 'Z_dilep_flav == 2'
Z_mumu = 'Z_dilep_flav == 3'
Z_DF_OS = 'Z_dilep_flav <= 1 && l_q[0]*l_q[1]<0'
Z_SF_OS = 'Z_dilep_flav > 1 && l_q[0]*l_q[1]<0'

# Same flavor lepton triggers
ee15_trig   = 'HLT_2e12_lhloose_L12EM10VH'
mumu15_trig = 'HLT_mu18_mu8noL1'
ee16_trig   = 'HLT_2e17_lhvloose_nod0'
mumu16_trig = 'HLT_mu22_mu8noL1'
# Different flavor lepton triggs
emu15_trig  = 'HLT_e17_lhloose_mu14'
emu16_trig  = 'HLT_e17_lhloose_nod0_mu14'
# Single lepton triggers
e15_trig    = '(HLT_e24_lhmedium_L1EM20VH || HLT_e60_lhmedium || HLT_e120_lhloose)'
mu15_trig   = '(HLT_mu20_iloose_L1MU15 || HLT_mu40)'
e16_trig    = '(HLT_e26_lhtight_nod0_ivarloose || HLT_e60_lhmedium_nod0 || HLT_e140_lhloose_nod0)'
mu16_trig   = '(HLT_mu26_ivarmedium || HLT_mu50)'

# Triggers with added pT requirements
e15_trig_emu_pT   = '(%s && dilep_flav == 0 && l_pt[0] >= 25)'%e15_trig
mu15_trig_emu_pT  = '(%s && dilep_flav == 0 && l_pt[0] < 25 && l_pt[1] >= 21)'%mu15_trig
e16_trig_emu_pT   = '(%s && dilep_flav == 0 && l_pt[0] >= 27)'%e16_trig
mu16_trig_emu_pT  = '(%s && dilep_flav == 0 && l_pt[0] < 27 && l_pt[1] >= 28)'%mu16_trig
emu15_trig_emu_pT = '(%s && dilep_flav == 0 && 18 <= l_pt[0] && l_pt[0] < 25 && 15 <= l_pt[1] && l_pt[1] < 21)'%emu15_trig
emu16_trig_emu_pT = '(%s && dilep_flav == 0 && 18 <= l_pt[0] && l_pt[0] < 27 && 15 <= l_pt[1] && l_pt[1] < 28)'%emu16_trig
e15_trig_mue_pT   = '(%s && dilep_flav == 1 && l_pt[1] >= 25)'%e15_trig
mu15_trig_mue_pT  = '(%s && dilep_flav == 1 && l_pt[1] < 25 && l_pt[0] >= 21)'%mu15_trig
e16_trig_mue_pT   = '(%s && dilep_flav == 1 && l_pt[1] >= 27)'%e16_trig
mu16_trig_mue_pT  = '(%s && dilep_flav == 1 && l_pt[1] < 27 && l_pt[0] >= 28)'%mu16_trig
emu15_trig_mue_pT = '(%s && dilep_flav == 1 && 18 <= l_pt[1] && l_pt[1] < 25 && 15 <= l_pt[0] && l_pt[0] < 21)'%emu15_trig
emu16_trig_mue_pT = '(%s && dilep_flav == 1 && 18 <= l_pt[1] && l_pt[1] < 27 && 15 <= l_pt[0] && l_pt[0] < 28)'%emu16_trig
e15_trig_ee_pT   = '(%s && Z_dilep_flav == 2 && l_pt[0] >= 25)'%e15_trig
e16_trig_ee_pT   = '(%s && Z_dilep_flav == 2 && l_pt[0] >= 27)'%e16_trig
mu15_trig_mumu_pT  = '(%s && Z_dilep_flav == 3 && l_pt[0] >= 21)'%mu15_trig
mu16_trig_mumu_pT  = '(%s && Z_dilep_flav == 3 && l_pt[0] >= 27)'%mu16_trig

# Combined triggers
e15_trig_pT = '(%s || %s || %s)'%(e15_trig_emu_pT, e15_trig_mue_pT, e15_trig_ee_pT)
mu15_trig_pT = '(%s || %s || %s)'%(mu15_trig_emu_pT, mu15_trig_mue_pT,mu15_trig_mumu_pT)
e16_trig_pT = '(%s || %s || %s)'%(e16_trig_emu_pT, e16_trig_mue_pT, e16_trig_ee_pT)
mu16_trig_pT = '(%s || %s || %s)'%(mu16_trig_emu_pT, mu16_trig_mue_pT, mu16_trig_mumu_pT)
emu15_trig_pT = '(%s || %s)'%(emu15_trig_emu_pT, emu15_trig_mue_pT)
emu16_trig_pT = '(%s || %s)'%(emu16_trig_emu_pT, emu16_trig_mue_pT)

dilep15_trig      = '(treatAsYear==2015 && %s)'%emu15_trig
dilep16_trig     = '(treatAsYear==2016 && %s)'%emu16_trig
singlelep15_trig = '(treatAsYear==2015 && (%s || %s))'%(e15_trig,mu15_trig)
singlelep16_trig = '(treatAsYear==2016 && (%s || %s))'%(e16_trig,mu16_trig)
dilep_trig = '(%s || %s)'%(dilep15_trig, dilep16_trig)
singlelep_trig = '(%s || %s)'%(singlelep15_trig, singlelep16_trig)
lepton_trig = '(%s || %s)'%(singlelep_trig, dilep_trig)

dilep15_trig_pT      = '(treatAsYear==2015 && %s)'%emu15_trig_pT
dilep16_trig_pT     = '(treatAsYear==2016 && %s)'%emu16_trig_pT
singlelep15_trig_pT = '(treatAsYear==2015 && (%s || %s))'%(e15_trig_pT,mu15_trig_pT)
singlelep16_trig_pT = '(treatAsYear==2016 && (%s || %s))'%(e16_trig_pT,mu16_trig_pT)
dilep_trig_pT = '(%s || %s)'%(dilep15_trig_pT, dilep16_trig_pT)
singlelep_trig_pT = '(%s || %s)'%(singlelep15_trig_pT, singlelep16_trig_pT)
lepton_trig_pT = '(%s || %s)'%(singlelep_trig_pT, dilep_trig_pT)

########################################
# Create regions
REGIONS = []

zjets_FF_CR = '1'
zjets_FF_CR += ' && %s'%singlelep_trig_pT
zjets_FF_CR += ' && fabs(lep_d0sigBSCorr[0]) < 15 && fabs(lep_d0sigBSCorr[1]) < 15 && fabs(lep_d0sigBSCorr[2]) < 15'
zjets_FF_CR += ' && fabs(lep_z0SinTheta[0]) < 15 && fabs(lep_z0SinTheta[1]) < 15 && fabs(lep_z0SinTheta[2]) < 15'
zjets_FF_CR += ' && (80 < Z_MLL && Z_MLL < 100)'
zjets_FF_CR += ' && l_mT[2] < 40'
zjets_FF_CR += ' && (Z2_MLL < 80 || 100 < Z2_MLL)'
zjets_FF_CR += ' && MET < 60'
#zjets_FF_CR += ' && nBJets == 0'
#zjets_FF_CR += ' && nLJets <= 2'
#zjets_FF_CR += ' && fabs(l_eta[2]) < 1.45'
#zjets_FF_CR += ' && fabs(l_eta[2]) >= 1.45'
zjets_FF_truth_base = '(!isMC || (0 < l_truthClass[0] && l_truthClass[0] <= 2))' #Prompt Leading Lepton
zjets_FF_truth_base += ' && (!isMC || (0 < l_truthClass[1] && l_truthClass[1] <= 2))' #Prompt Subleading Lepton
zjets_FF_truth_num = zjets_FF_truth_base\
                   +' && (!isMC || (0 < l_truthClass[2] && l_truthClass[2] <= 2))' #Prompt Probe Lepton
zjets_FF_truth_den = zjets_FF_truth_base\
                   + ' && (!isMC || (l_truthClass[2] <= 0 || 2 < l_truthClass[2]))' #Fake Probe Lepton
num_den_dict = {'den' : 'nLepID == 2 && nLepAntiID >= 1',
                'num' : 'nLepID == 3'}
chan_dict = {'eee' : ['ee','e','Z_dilep_flav==2 && l_flav[2]==0'],
             'mme' : ['mumu','e','Z_dilep_flav==3 && l_flav[2]==0'],
             'eem' : ['ee','m','Z_dilep_flav==2 && l_flav[2]==1'],
             'mmm' : ['mumu','m','Z_dilep_flav==3 && l_flav[2]==1'],
             'm' : ['ll','m','l_flav[2]==1'],
             'e' : ['ll','e','l_flav[2]==0'],
        }
for num_den, num_den_sel in num_den_dict.iteritems():
    for chan, ops in chan_dict.iteritems():
        id_or_aid = 'ID' if num_den == 'num' else 'anti-ID'
        chan_name = '%s + %s %s'%(ops[0], id_or_aid, ops[1])

        name = 'zjets_FF_CR%s_%s'%(num_den, chan)
        displayname = 'Z+jets CR (%s)'%(chan_name)
        #displayname = 'Z+jets CR (%s, |#eta| < 1.4)'%(chan_name)
        #displayname = 'Z+jets CR (%s, |#eta| #ge 1.4)'%(chan_name)
        REGIONS.append(Region(name, displayname))
        REGIONS[-1].tcut = ' && '.join([ops[2], zjets_FF_CR])
        REGIONS[-1].truth_fake_sel = zjets_FF_truth_den
        REGIONS[-1].truth_bkg_sel = zjets_FF_truth_num





################################################################################
# Variables
# - Create all the variables that one might like to plot
################################################################################

# Improved plot defining setup
# These plots will be copied into new plots for each region being run
plot_defaults = {
    'l_pt[2]'            : Plot1D( bin_range=[0.0, 1000.0],  bin_width=0.5, ptype=Types.stack, doLogY=False, add_overflow = False, xunits='GeV', xlabel='Fake candidate lepton p_{T}'),
    'l_eta[2]'           : Plot1D( bin_range=[-3.0, 3.0, 0, 4000],   bin_width=0.01, ptype=Types.stack, doLogY= False, add_underflow = True, xlabel='Fake candidate lepton #eta'),
    'fabs(l_eta[2]):l_pt[2]'   : Plot2D( bin_range=[0, 1000,0, 3.0], xbin_width = 1, ybin_width = 0.01, ylabel='Fake probe lepton |#eta|', xunits='GeV', xlabel='Leading lepton p_{T}')
}
plot_defaults['l_pt[2]'].rebin_bins = [0,10,11,15,20,25,35,1000]
plot_defaults['l_eta[2]'].rebin_bins = [-3.0, -2.5, -1.45, 0, 1.45, 2.5, 3.0]
plot_defaults['fabs(l_eta[2]):l_pt[2]'].rebin_xbins = plot_defaults['l_pt[2]'].rebin_bins 
plot_defaults['fabs(l_eta[2]):l_pt[2]'].rebin_ybins = [0, 1.45, 2.5, 3.0] 

region_plots = {}
################################################################################
# Toggle options for execution
# - Options expected to change often between different plottings
################################################################################


#######################################
# Yield Table
YIELD_TBL = YieldTable()
# Add formulas to yield table
# Key values will be the labels on the table
# Formulas should use sample names and these symbols: +, -, *, /, (, ), [0-9]
#YIELD_TBL.formulas['Data/MC'] = "data/MC"
#YIELD_TBL.formulas['1 - (WZ+ZZ)/Data'] = "1 - (wz + zz)/data"
#YIELD_TBL.formulas['Zll/Data'] = "zll/data"

#######################################
# What regions to plot
region_ops = []
region_ops += ['zjets_FF_CRden_e', 'zjets_FF_CRnum_e']
region_ops += ['zjets_FF_CRden_m', 'zjets_FF_CRnum_m']
#region_ops += ['zjets_FF_CRden_eem', 'zjets_FF_CRnum_eem']
#region_ops += ['zjets_FF_CRden_mmm', 'zjets_FF_CRnum_mmm']
#region_ops += ['zjets_FF_CRden_eee', 'zjets_FF_CRnum eee']
#region_ops += ['zjets_FF_CRden_mme', 'zjets_FF_CRnum_mme']
#region_ops += ['zjets_FF_CRden_eem']
#######################################
# What variables to plot
vars_to_plot = []
#vars_to_plot += ['l_pt[2]']
#vars_to_plot += ['l_eta[2]']
vars_to_plot = ['fabs(l_eta[2]):l_pt[2]']

# Remove duplicate names
vars_to_plot = list(set(vars_to_plot))

################################################################################
# Create plots
################################################################################
PLOTS = []
for var in vars_to_plot:

    # Use defualt settings, if needed, when plotting a specific element of a
    # vector (i.e. lepton_flav[1] uses lepton_flav settings if lepton_flav[1]
    # is not defined separately)
    if '[' in var and var not in plot_defaults:
        key = var.split('[')[0]
    else:
        key = var

    # Create plot for each region
    for region in region_ops:

        # Grab the default plot unless a region specific one is defined
        if region in region_plots and key in region_plots[region]:
            p = deepcopy(region_plots[region][key])
        elif key in plot_defaults:
            p = deepcopy(plot_defaults[key])
        else:
            assert False, ("ERROR :: requested plot not defined:", var)

        if p.is2D:
            varx = var.split(':')[1]
            vary = var.split(':')[0]
            p.update(region, varx, vary)
        else:
            p.update(region, var)

        # Set plot type if not already set
        if p.ptype == Types.default:
            n_bkgds = len([s for s in SAMPLES if s.isMC and not s.isSignal])
            n_signal = len([s for s in SAMPLES if s.isMC and s.isSignal])
            n_data = len([s for s in SAMPLES if not s.isMC])

            if n_bkgds and n_data:
                p.ptype = Types.ratio
            elif n_bkgds or n_data:
                 p.ptype = Types.stack
            else:
                assert False, "ERROR :: No samples are setup"

        # Setup correct canvas for plotting
        if p.ptype == Types.ratio:
            p.setRatioPads(p.name)
        elif p.ptype == Types.stack:
            p.setStackPads(p.name)
        elif p.ptype == Types.two_dim:
            p.set2DPads(p.name)
        else:
            print "WARNING :: %s plots are not yet setup"%p.ptype.name
            continue

        PLOTS.append(p)

print "Configuration File End"
