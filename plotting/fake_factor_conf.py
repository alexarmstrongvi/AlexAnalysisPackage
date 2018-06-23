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
from tools.plot import Plot1D, Types
from tools.sample import Sample, MCsample, Data, Background, Signal
from tools.region import Region
import tools.systematic as systematic
from tools.YieldTable import YieldTable, UncFloat


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
Plot1D.ymax = 8e3
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

    # Top (ttbar + singletop + Wt)
    top = Background("top_%s"%num_den, "Top")
    top.color = ROOT.kYellow+1
    SAMPLES.append(top)
    
    # VV combined
    VV = Background("ggllvv_%s"%num_den, "ggllvv")
    VV.color = ROOT.kSpring-7
    SAMPLES.append(VV)

    # WW
    WW = Background("ww_%s"%num_den, "WW")
    WW.color = ROOT.kTeal-5
    SAMPLES.append(WW)

    # ZZ
    ZZ = Background("zz_%s"%num_den, "ZZ")
    ZZ.color = ROOT.kMagenta-5
    SAMPLES.append(ZZ)

    # WZ
    WZ = Background("wz_%s"%num_den, "WZ")
    WZ.color = ROOT.kSpring+9
    SAMPLES.append(WZ)

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

    # Higgs -> tau tau
    htt = Background("htt_%s"%num_den, "H#tau#tau")
    htt.color = ROOT.kRed - 7
    SAMPLES.append(htt)

    # Higgs -> W W
    hww = Background("hww_%s"%num_den, "HWW")
    hww.color = ROOT.kBlue+3
    SAMPLES.append(hww)

    #######################################
    ## Build the TChain/TTree for each sample
    # To remove sample from plot, comment out the line setting its TChain
    # Samples with empty TChains get removed below

    data_dsids = g.groups['data15']+g.groups['data16']
    data_ntuple_dir2 = data_ntuple_dir[:-1] + '_' + num_den + "/"
    bkg_ntuple_dir2 = bkg_ntuple_dir[:-1] + '_' + num_den + "/"

    data.set_chain_from_dsid_list(data_dsids, data_ntuple_dir2)
    top.set_chain_from_dsid_list(g.groups['top'], bkg_ntuple_dir2)
    #ttbar.set_chain_from_dsid_list(g.groups['ttbar'], bkg_ntuple_dir2)
    #stop.set_chain_from_dsid_list(g.groups['singletop'], bkg_ntuple_dir2)
    #wtop.set_chain_from_dsid_list(g.groups['Wt'], bkg_ntuple_dir2)
    VV.set_chain_from_dsid_list(g.groups['ggllvv'], bkg_ntuple_dir2)
    WW.set_chain_from_dsid_list(g.groups['ww'], bkg_ntuple_dir2)
    ZZ.set_chain_from_dsid_list(g.groups['zz'], bkg_ntuple_dir2)
    WZ.set_chain_from_dsid_list(g.groups['wz'], bkg_ntuple_dir2)
    zll.set_chain_from_dsid_list(g.groups['zll'], bkg_ntuple_dir2)
    #zee.set_chain_from_dsid_list(g.groups['zee'], bkg_ntuple_dir2)
    #zmumu.set_chain_from_dsid_list(g.groups['zmumu'], bkg_ntuple_dir2)
    ztt.set_chain_from_dsid_list(g.groups['ztt'], bkg_ntuple_dir2)
    wjets.set_chain_from_dsid_list(g.groups['wjets'], bkg_ntuple_dir2)
    #wgamma.set_chain_from_dsid_list(g.groups['w_gamma'], bkg_ntuple_dir2)
    #htt.set_chain_from_dsid_list(g.groups['htt'], bkg_ntuple_dir2)
    #hww.set_chain_from_dsid_list(g.groups['hww'], bkg_ntuple_dir2)

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
e15_trig_emu_pT   = '(%s && dilep_flav == 0 && Lep0Pt >= 25)'%e15_trig
mu15_trig_emu_pT  = '(%s && dilep_flav == 0 && Lep0Pt < 25 && Lep1Pt >= 21)'%mu15_trig
e16_trig_emu_pT   = '(%s && dilep_flav == 0 && Lep0Pt >= 27)'%e16_trig
mu16_trig_emu_pT  = '(%s && dilep_flav == 0 && Lep0Pt < 27 && Lep1Pt >= 28)'%mu16_trig
emu15_trig_emu_pT = '(%s && dilep_flav == 0 && 18 <= Lep0Pt && Lep0Pt < 25 && 15 <= Lep1Pt && Lep1Pt < 21)'%emu15_trig
emu16_trig_emu_pT = '(%s && dilep_flav == 0 && 18 <= Lep0Pt && Lep0Pt < 27 && 15 <= Lep1Pt && Lep1Pt < 28)'%emu16_trig
e15_trig_mue_pT   = '(%s && dilep_flav == 1 && Lep1Pt >= 25)'%e15_trig
mu15_trig_mue_pT  = '(%s && dilep_flav == 1 && Lep1Pt < 25 && Lep0Pt >= 21)'%mu15_trig
e16_trig_mue_pT   = '(%s && dilep_flav == 1 && Lep1Pt >= 27)'%e16_trig
mu16_trig_mue_pT  = '(%s && dilep_flav == 1 && Lep1Pt < 27 && Lep0Pt >= 28)'%mu16_trig
emu15_trig_mue_pT = '(%s && dilep_flav == 1 && 18 <= Lep1Pt && Lep1Pt < 25 && 15 <= Lep0Pt && Lep0Pt < 21)'%emu15_trig
emu16_trig_mue_pT = '(%s && dilep_flav == 1 && 18 <= Lep1Pt && Lep1Pt < 27 && 15 <= Lep0Pt && Lep0Pt < 28)'%emu16_trig
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

zjets_FF_CR_add = '1'
zjets_FF_CR_add += ' && %s'%singlelep_trig_pT
zjets_FF_CR_add += ' && (75 < Z_MLL && Z_MLL < 105)'
zjets_FF_CR_add += ' && nBJets == 0'
zjets_FF_CR_add += ' && l_mT[2] < 50'
zjets_FF_CR_add += ' && (Z2_MLL < 80 || 100 < Z2_MLL)'
zjets_FF_CR_add += ' && MET < 50'
zjets_FF_CR_add += ' && (l_flav[2]!=1 || nLepID!=2 || fabs(lep_d0sigBSCorr[2]) < 15)'
zjets_FF_CR_add += ' && (l_flav[2]!=1 || nLepID!=3 || fabs(lep_d0sigBSCorr[2]) < 3)'
zjets_FF_CR_add += ' && (l_flav[2]!=0 || fabs(lep_d0sigBSCorr[2]) < 5)'
zjets_FF_CR_add += ' && (fabs(lep_z0SinTheta[2]) < 0.5)'
#zjets_FF_CR_add += ' && (l_flav[2] != 1 || l_ID[2] > 1)'
zjets_FF_CR_add += ' && (l_flav[2] != 1 || l_ID[2] <= 1)'
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
        displayname = 'Z+jets FF CR (%s)'%(chan_name)
        REGIONS.append(Region(name, displayname))
        REGIONS[-1].tcut = ' && '.join([ops[2], zjets_FF_CR_add])
        REGIONS[-1].truth_fake_sel = zjets_FF_truth_den
        REGIONS[-1].truth_bkg_sel = zjets_FF_truth_num


# Wjet fake regions
wjets_FF_CRden_base = 'nLepID == 1 && nLepAntiID >= 1'
wjets_FF_CRnum_base = 'nLepID == 2'
wjets_FF_CRden_add = '1'
wjets_FF_CRden_add += ' && l_q[0] == aID_Lep0Q && aID_dilep_flav <= 1' # reject SFOS zjets
wjets_FF_CRden_add += ' && (aID_MLL < 70 || 110 < aID_MLL)' # reject zjets
wjets_FF_CRden_add += ' && nBJets==0' # reject top
FF_CRden_emu_ch  = 'aID_dilep_flav == 0'
FF_CRden_mue_ch  = 'aID_dilep_flav == 1'
FF_CRden_ee_ch   = 'aID_dilep_flav == 2'
FF_CRden_mumu_ch = 'aID_dilep_flav == 3'

REGIONS.append(Region('wjets_FF_CRden_emu', 'wjets FF CR (anti-ID el)'))
REGIONS[-1].tcut = wjets_FF_CRden_base + '&&' + wjets_FF_CRden_add + '&&' + FF_CRden_emu_ch

REGIONS.append(Region('wjets_FF_CRden_mue', 'wjets FF CR (anti-ID mu)'))
REGIONS[-1].tcut = wjets_FF_CRden_base + '&&' + wjets_FF_CRden_add + '&&' + FF_CRden_mue_ch

wjets_FF_CRnum_base = 'nLepID == 2' # require final state
wjets_FF_CRnum_add = '1'
wjets_FF_CRnum_add += ' && l_q[0] == l_q[1]' # reject zjets
wjets_FF_CRnum_add += ' && (MLL < 70 || 110 < MLL)' # reject zjets
wjets_FF_CRnum_add += ' && nBJets==0' # reject top
FF_CRnum_emu_ch  = 'dilep_flav == 0'
FF_CRnum_mue_ch  = 'dilep_flav == 1'
FF_CRnum_ee_ch   = 'dilep_flav == 2'
FF_CRnum_mumu_ch = 'dilep_flav == 3'

REGIONS.append(Region('wjets_FF_CRnum_emu', 'wjets FF CR (ID el)'))
REGIONS[-1].tcut = wjets_FF_CRnum_base + '&&' + wjets_FF_CRnum_add + '&&' + FF_CRnum_emu_ch

REGIONS.append(Region('wjets_FF_CRnum_mue', 'wjets FF CR (ID mu)'))
REGIONS[-1].tcut = wjets_FF_CRnum_base + '&&' + wjets_FF_CRnum_add + '&&' + FF_CRnum_emu_ch


################################################################################
# Variables
# - Create all the variables that one might like to plot
################################################################################

# Improved plot defining setup
# These plots will be copied into new plots for each region being run
plot_defaults = {
    'l_pt[2]'            : Plot1D( bin_range=[0.0, 100.0],  bin_width=5, ptype=Types.stack, doLogY=False, add_overflow = False, xunits='GeV', xlabel='Fake candidate lepton p_{T}'),
    'l_eta[2]'           : Plot1D( bin_range=[-3.0, 3.0],   nbins=20, ptype=Types.stack, doLogY= False, add_underflow = True, xlabel='Fake candidate lepton #eta'),
}
plot_defaults['l_pt[2]'].rebin_bins = [0,5,10,15,20,25,30,35,100]
plot_defaults['l_eta[2]'].rebin_bins = [-3.0, -2.5, -1.5, 1.5, 2.5, 3.0]

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
#region_ops += ['zjets_FF_CRden_e', 'zjets_FF_CRnum_e']
region_ops += ['zjets_FF_CRden_m', 'zjets_FF_CRnum_m']
#region_ops += ['zjets_FF_CRden_eem', 'zjets_FF_CRnum_eem']
#region_ops += ['zjets_FF_CRden_mmm', 'zjets_FF_CRnum_mmm']
#region_ops += ['zjets_FF_CRden_eee', 'zjets_FF_CRnum eee']
#region_ops += ['zjets_FF_CRden_mme', 'zjets_FF_CRnum_mme']
#region_ops += ['zjets_FF_CRden_eem']
#######################################
# What variables to plot
vars_to_plot = []
vars_to_plot += ['l_pt[2]']
#vars_to_plot += ['l_eta[2]']

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
        else:
            print "WARNING :: %s plots are not yet setup"%p.ptype.name
            continue

        PLOTS.append(p)

print "Configuration File End"
