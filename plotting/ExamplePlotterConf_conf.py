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

################################################################################
# Samples
# - Create all the samples that will be put into plots (e.g. backgrounds,
#   signal, data)
################################################################################
SAMPLES = []

# Data
data = Data()
SAMPLES.append(data)

################################################################################
# Fakes
## Fakes
fakes = Background("fakes", "Fakes")
fakes.scale_factor = 1 # Derived from data so no need for scaling
fakes.color = ROOT.kGray
SAMPLES.append(fakes)

################################################################################
# Signal
signal_branching_ratio = 0.01
signal_SF = 1
signal_label = "Higgs LFV" if signal_SF == 1 else "Higgs LFV (%dX)"%signal_SF

signal = Signal("higgs_lfv", signal_label)
signal.scale_factor = lumi * signal_branching_ratio * signal_SF
signal.color = ROOT.kGreen
SAMPLES.append(signal)

################################################################################
# Backgrounds

#######################################
# Initialize all backgrounds
# ttbar
ttbar = Background("ttbar", "t#bar{t}")
ttbar.color = ROOT.kOrange+2
SAMPLES.append(ttbar)

# singletop
stop = Background("st", "Single top")
stop.color = ROOT.kOrange+1
SAMPLES.append(stop)

# W+top
wtop = Background("wt", "Wt")
wtop.color = ROOT.kOrange+8
SAMPLES.append(wtop)

# VV combined
VV = Background("ggllvv", "ggllvv")
VV.color = ROOT.kSpring-7
SAMPLES.append(VV)

# WW
WW = Background("ww", "WW")
WW.color = ROOT.kSpring-6
SAMPLES.append(WW)

# ZZ
ZZ = Background("zz", "ZZ")
ZZ.color = ROOT.kSpring-4
SAMPLES.append(ZZ)

# WZ
WZ = Background("wz", "WZ")
WZ.color = ROOT.kSpring-5
SAMPLES.append(WZ)

# Zll
zll = Background("zll", "Zll")
zll.color = ROOT.kAzure-9
SAMPLES.append(zll)

# Zee
zee = Background("zee", "Zee")
zee.color = ROOT.kAzure-7
SAMPLES.append(zee)

# Zmumu
zmumu = Background("zmumu", "Zmumu")
zmumu.color = ROOT.kAzure-9
SAMPLES.append(zmumu)

# Ztt
ztt = Background("ztt", "Z#tau#tau")
ztt.color = ROOT.kAzure-5
SAMPLES.append(ztt)

# Wjets
wjets = Background("wjets", "W+jets")
wjets.color = ROOT.kOrange
SAMPLES.append(wjets)

# W+gamma
wgamma = Background("wgamma", "W+gamma")
wgamma.color = ROOT.kOrange-1
SAMPLES.append(wgamma)

# Higgs -> tau tau
htt = Background("htt", "H#tau#tau")
htt.color = ROOT.kRed
SAMPLES.append(htt)

# Higgs -> W W
hww = Background("hww", "HWW")
hww.color = ROOT.kBlue+3
SAMPLES.append(hww)


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

# Region building blocks
Baseline_Sel = ('Lep0Pt >= 45 && Lep1Pt >= 15 '
              + '&& (30 < MLL && MLL < 150) '
              + '&& nBJets==0 '
              + '&& ( !'+mue+' || el1pT_trackclus_ratio < 1.2) '
              + '&&' + DF_OS
              #+ '&& dilep_flav <= 1 '
              + '&&' + lepton_trig_pT)
VBF_stripped = "JetN_g30 >= 2 && j_pt[0] > 40 && Mjj > 400 && DEtaJJ > 3"


########################################
# Create regions
REGIONS = []

REGIONS.append(Region(name = "no_sel", displayname = "No Selections"))
REGIONS[-1].tcut = '1' #DF_OS

REGIONS.append(Region("trig_only", "Lepton Triggers"))
REGIONS[-1].tcut = lepton_trig_pT + "&&" + DF_OS

REGIONS.append(Region("topCR", "Top CR"))
REGIONS[-1].tcut = "nBJets >= 1 && MET > 40 &&" + DF_OS + " &&" + lepton_trig_pT

zll_cr_base = "75 < Z_MLL && Z_MLL < 105 && " + Z_SF_OS
zll_cr_add  = singlelep_trig_pT
#zll_cr_add += " && nLepID == 2"
#zll_cr_add += " && nLepAntiID >= 1"
zll_cr_add += " && nLepID == 3"
#zll_cr_add += " && Z_Lep2_flav==1" # electron
zll_cr_add += " && Z_Lep2_flav==0" # muon
#zll_cr_add += " && Z_dilep_flav==2" # el-el
#zll_cr_add += " && Z_dilep_flav==3" # mu-mu
REGIONS.append(Region("zCR_ee", "Z CR (Channel: El-El)"))
REGIONS[-1].tcut = zll_cr_base + " && " + zll_cr_add + " && " + Z_ee

REGIONS.append(Region("zCR_mumu", "Z CR (Channel: Mu-Mu)"))
REGIONS[-1].tcut = zll_cr_base + " && " + zll_cr_add + " && " + Z_mumu

ZTauTau_CR =  ('Lep0Pt >= 30 && Lep0Pt < 45 && Lep1Pt >= 15 '
              + '&& (30 < MLL && MLL < 150) '
              + '&& nBJets==0 '
              + '&& ( !'+mue+' || el1pT_trackclus_ratio < 1.2) '
              + '&&' + DF_OS)
REGIONS.append(Region("ztautauCR", "Ztautau CR"))
REGIONS[-1].tcut = ZTauTau_CR

# Z+Jets fake regions
zjets_FF_CRden_base = 'nLepID == 2 && nLepAntiID >= 1'
zjets_FF_CRnum_base = 'nLepID == 3'
zjets_FF_CR_add = singlelep_trig_pT
zjets_FF_CR_add += '&& (75 < Z_MLL && Z_MLL < 105)'
zjets_FF_CR_add += ' && nBJets == 0'
zjets_FF_CR_add += ' && Z_Lep2_mT < 50'
zjets_FF_CR_add += ' && (Z2_MLL < 80 || 100 < Z2_MLL)'
zjets_FF_CR_add += ' && MET < 50'
zjets_FF_CR_add += ' && (!isMC || (0 < l_truthClass[0] && l_truthClass[0] <= 2))' #Prompt Leading Lepton
zjets_FF_CR_add += ' && (!isMC || (0 < l_truthClass[1] && l_truthClass[1] <= 2))' #Prompt Subleading Lepton
#zjets_FF_CR_add += ' && (!isMC || (0 < l_truthClass[2] && l_truthClass[2] <= 2))' #Prompt Probe Lepton
zjets_FF_CR_add += ' && (!isMC || (l_truthClass[2] <= 0 || 2 < l_truthClass[2]))' #Fake Probe Lepton
num_den_dict = {'den' : 'nLepID == 2 && nLepAntiID >= 1',
                'num' : 'nLepID == 3'}
chan_dict = {'eee' : ['ee','e','Z_dilep_flav==2 && Z_Lep2_flav==1'],
             'mme' : ['mumu','e','Z_dilep_flav==3 && Z_Lep2_flav==1'],
             'eem' : ['ee','m','Z_dilep_flav==2 && Z_Lep2_flav==0'],
             'mmm' : ['mumu','m','Z_dilep_flav==3 && Z_Lep2_flav==0'],
             'm' : ['ll','m','Z_Lep2_flav==0'],
             'e' : ['ll','e','Z_Lep2_flav==1'],
        }
for num_den, num_den_sel in num_den_dict.iteritems():
    for chan, ops in chan_dict.iteritems():
        id_or_aid = 'ID' if num_den == 'num' else 'anti-ID'
        chan_name = '%s + %s %s'%(ops[0], id_or_aid, ops[1])

        name = 'zjets_FF_CR%s_%s'%(num_den, chan)
        displayname = 'Z+jets FF CR (%s)'%(chan_name)
        REGIONS.append(Region(name, displayname))
        REGIONS[-1].tcut = ' && '.join([num_den_sel, zjets_FF_CR_add, ops[2], singlelep_trig])

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

# Baseline regions
REGIONS.append(Region("baseline", "Baseline"))
REGIONS[-1].tcut = Baseline_Sel

REGIONS.append(Region("baseline_emu", "Baseline (Channel: El-Mu)"))
REGIONS[-1].tcut = Baseline_Sel + " && " + emu

REGIONS.append(Region("baseline_mue", "Baseline (Channel: Mu-El)"))
REGIONS[-1].tcut = Baseline_Sel + " && " + mue

REGIONS.append(Region("vbf", "VBF"))
REGIONS[-1].tcut = "(%s) && (%s)"%(Baseline_Sel, VBF_stripped)

REGIONS.append(Region("optimized", "Optimized"))
REGIONS[-1].tcut = "(%s) && !(%s) && DphiLep1MET < 1 && MtLep0 > 50 && MtLep1 < 40 && ((MET+Lep1Pt)/Lep1Pt) > 0.5"%(Baseline_Sel, VBF_stripped)

################################################################################
# Variables
# - Create all the variables that one might like to plot
################################################################################

# Improved plot defining setup
# These plots will be copied into new plots for each region being run
plot_defaults = {
    # Event level
    'eventweight'          : Plot1D( bin_range=[-0.001, 0.002], nbins=100, xlabel='Event weight'),
    'treatAsYear'          : Plot1D( bin_range=[2007.5, 2018.5], bin_width=1, xlabel='treatAsYear'),
    'isMC'                 : Plot1D( bin_range=[-1.5, 2.5],   bin_width=1, ptype=Types.stack, xlabel='is Monte Carlo'),
    # Multiplicity
    'n_preLeptons'         : Plot1D( bin_range=[-0.5, 10.5],  bin_width=1, xlabel='N_{pre-leptons}'),
    'n_baseLeptons'        : Plot1D( bin_range=[-0.5, 10.5],  bin_width=1, xlabel='N_{baseline leptons}'),
    'n_leptons'            : Plot1D( bin_range=[-0.5, 10.5],  bin_width=1, xlabel='N_{signal leptons}'),
    # Leptons
    'l_pt'                 : Plot1D( bin_range=[0.0, 500.0],  nbins=25, xunits='GeV', xlabel='Lepton p_{T}'),
    'aID_Lep0Pt'           : Plot1D( bin_range=[0.0, 200.0],  nbins=40, xunits='GeV', xlabel='p_{T}^{leading lep} anti-ID'),
    ## Electrons
    'el1_track_pt'         : Plot1D( bin_range=[0.0, 500.0],  nbins=25, xunits='GeV', xlabel='Leading electron track p_{T}'),
    'el1_clus_pt'          : Plot1D( bin_range=[0.0, 500.0],  nbins=25, xunits='GeV', xlabel='Leading electron cluster p_{T}'),
    'preEl_Et'             : Plot1D( bin_range=[0.0, 50.0],   nbins=50, add_overflow=False, xunits='GeV', xlabel='E_{T}^{preElectron}'),
    'preEl_pt'             : Plot1D( bin_range=[0.0, 50.0],   nbins=50, xunits='GeV', xlabel='p_{T}^{preElectron}'),
    'preEl_clusEtaBE'      : Plot1D( bin_range=[-3.0, 3.0],   nbins=60, xlabel='#eta_{preEl clusterBE}'),
    'preEl_eta'            : Plot1D( bin_range=[-3.0, 3.0],   nbins=60, xlabel='#eta_{preElectron}'),
    'baseEl_ID'            : Plot1D( bin_range=[-0.5, 8.5],   bin_width=1, xlabel='el_{Baseline}^{ID}(non-inclusive)'),
    'El_ID'                : Plot1D( bin_range=[-1.5, 5.5],   bin_width=1, doLogY=False, xlabel='electron ID(non-inclusive)'),
    'el_type'              : Plot1D( bin_range=[-1.5, 39.5],  nbins=41, xlabel='electron truth type'),
    'el_origin'            : Plot1D( bin_range=[-1.5, 46.5],  nbins=48, xlabel='electron truth origin'),
    'el_d0sigBSCorr'       : Plot1D( bin_range=[-6, 6],       nbins=60, add_underflow=True, xunits='mm', xlabel='Electron d_{0}/#sigma_{d_{0}} BSCorr'),
    'el_z0SinTheta'        : Plot1D( bin_range=[-0.6, 0.6],   nbins=60, add_underflow=True, xunits='mm', xlabel='Electron z_{0}sin(#theta)'),
    ## Muons
    'preMu_pt'             : Plot1D( bin_range=[0, 50.0],     nbins=50, add_overflow=False, xunits='GeV', xlabel='p_{T}^{preMuon}'),
    'preMu_ID'             : Plot1D( bin_range=[-1.5, 5.5],   bin_width=1, xlabel='mu_{pre}^{ID} (non-inclusive)'),
    'baseMu_pt'            : Plot1D( bin_range=[0.0, 200.0],  nbins=20, xunits='GeV', xlabel='p_{T}^{BaselineMuon}'),
    'baseMu_eta'           : Plot1D( bin_range=[-3.0, 3.0],   nbins=20, xlabel='#eta_{BaselineMuon}'),
    'baseMu_ID'            : Plot1D( bin_range=[-1.5, 5.5],   bin_width=1, xlabel='mu_{Baseline}^{ID}(non-inclusive)'),
    'Mu_ID'                : Plot1D( bin_range=[-1.5, 5.5],   bin_width=1, doLogY=False, xlabel='muon ID (non-inclusive)'),
    'mu_type'              : Plot1D( bin_range=[-1.5, 39.5],  bin_width=1, xlabel='muon truth type'),
    'mu_origin'            : Plot1D( bin_range=[-1.5, 46.5],  bin_width=1, xlabel='muon truth origin'),
    'mu_d0sigBSCorr'       : Plot1D( bin_range=[-6, 6],       nbins=60, add_underflow=True, xunits='mm', xlabel='Muon d_{0}/#sigma_{d_{0}} BSCorr'),
    'mu_z0SinTheta'        : Plot1D( bin_range=[-0.6, 0.6],   nbins=60, add_underflow=True, xunits='mm', xlabel='Muon z_{0}sin(#theta)'),
    ## Taus
    'preTau_q'             : Plot1D( bin_range=[-5.5, 5.5],   bin_width=1, xlabel='Tau charge'),
    'preTau_nTracks'       : Plot1D( bin_range=[-1.5, 8.5],   bin_width=1, xlabel='preTau nTracks'),
    'baseTau_pT'           : Plot1D( bin_range=[0.0, 200.0],  nbins=20, xunits='GeV', xlabel='p_{T}^{BaselineTau}'),
    'baseTau_eta'          : Plot1D( bin_range=[-3.0, 3.0],   nbins=20, xlabel='#eta_{BaselineMuon}'),
    'baseTau_nTracks'      : Plot1D( bin_range=[-1.5, 8.5],   bin_width=1, xlabel='Baseline Tau nTracks'),
    'baseTau_ID'           : Plot1D( bin_range=[-1.5, 5.5],   bin_width=1, xlabel='tau_{Baseline}^{ID}(non-inclusive)'),

    ## General lepton
    'l_eta'                : Plot1D( bin_range=[-3.0, 3.0],   nbins=20, xlabel='Lepton #eta'),
    'l_phi'                : Plot1D( bin_range=[0.0, 3.15],   nbins=30, xlabel='Lepton #phi'),
    'l_flav'               : Plot1D( bin_range=[-0.5, 4.5],   bin_width=1, xlabel='Lepton flavor (0: e, 1: m)'),
    'l_type'               : Plot1D( bin_range=[-1.5, 39.5],  bin_width=1, ptype=Types.stack, xlabel='Lepton type'),
    'l_origin'             : Plot1D( bin_range=[-1.5, 46.5],  bin_width=1, ptype=Types.stack, xlabel='Lepton origin'),
    'l_BkgMotherPdgId'     : Plot1D( bin_range=[-20.5, 20.5], bin_width=1, add_underflow=True, ptype=Types.stack, xlabel='Lepton Mother PdgID'),
    'l_truthClass'         : Plot1D( bin_range=[-1.5, 11.5],  bin_width=1, doNorm=True, doLogY=False, ptype=Types.stack, xlabel='Lepton truth classification'),
    'l_truthClass[0]'      : Plot1D( bin_range=[-1.5, 11.5],  bin_width=1, doNorm=True, doLogY=False, ptype=Types.stack, xlabel='Leading lepton truth classification'),
    'l_truthClass[1]'      : Plot1D( bin_range=[-1.5, 11.5],  bin_width=1, doNorm=True, doLogY=False, ptype=Types.stack, xlabel='Subleading lepton truth classification'),
    'l_truthClass[2]'      : Plot1D( bin_range=[-1.5, 11.5],  bin_width=1, doNorm=True, doLogY=False, ptype=Types.stack, xlabel='Fake candidate lepton Truth Classification'),
    'l_type[0]'            : Plot1D( bin_range=[-1.5, 39.5],  bin_width=1, ptype=Types.stack, xlabel='Leading Z Lepton type'),
    'l_origin[0]'          : Plot1D( bin_range=[-1.5, 46.5],  bin_width=1, ptype=Types.stack, xlabel='Leading Z Lepton origin'),
    'l_type[1]'            : Plot1D( bin_range=[-1.5, 39.5],  bin_width=1, ptype=Types.stack, xlabel='Subleading Z Lepton type'),
    'l_origin[1]'          : Plot1D( bin_range=[-1.5, 46.5],  bin_width=1, ptype=Types.stack, xlabel='Subleading Z Lepton origin'),
    'l_type[2]'            : Plot1D( bin_range=[-1.5, 39.5],  bin_width=1, ptype=Types.stack, xlabel='Fake Candidate Lepton type'),
    'l_origin[2]'          : Plot1D( bin_range=[-1.5, 46.5],  bin_width=1, ptype=Types.stack, xlabel='Fake Candidate Lepton origin'),
    'l_BkgMotherPdgId[2]'  : Plot1D( bin_range=[-20.5, 20.5], bin_width=1, doNorm=True, add_underflow=True, ptype=Types.stack, xlabel='Fake Candidate Lepton Mother PdgID'),
    'Lep_Iso'              : Plot1D( bin_range=[-1.5, 7.5],   bin_width=1, xlabel='Lepton Isolation (non-inclusive)'),
    'l_q'                  : Plot1D( bin_range=[-1.5, 1.5],   bin_width=1, xlabel='Lepton charge'),
    'LepLepSign'           : Plot1D( bin_range=[-1.5, 1.5],   bin_width=1, xlabel='Leptons sign product'),
    'l_pt[0]'              : Plot1D( bin_range=[0.0, 200.0],  nbins=40, xunits='GeV', xlabel='p_{T}^{leading lep}'),
    'l_pt[1]'              : Plot1D( bin_range=[0.0, 200.0],  nbins=40, xunits='GeV', xlabel='p_{T}^{subleading lep}'),
    'Lep0Pt'               : Plot1D( bin_range=[0.0, 200.0],  nbins=40, xunits='GeV', xlabel='p_{T}^{leading lep}'),
    'Lep1Pt'               : Plot1D( bin_range=[0.0, 200.0],  nbins=40, xunits='GeV', xlabel='p_{T}^{subleading lep}'),
    'Lep0Eta'              : Plot1D( bin_range=[-3.0, 3.0],   nbins=20, xlabel='#eta^{leading lep}'),
    'Lep1Eta'              : Plot1D( bin_range=[-3.0, 3.0],   nbins=20, xlabel='#eta^{subleading lep}'),
    'Lep0Phi'              : Plot1D( bin_range=[0.0, 3.15],   nbins=30, xlabel='#phi^{leading lep}'),
    'Lep1Phi'              : Plot1D( bin_range=[0.0, 3.15],   nbins=30, xlabel='#phi^{subleading lep}'),
    'MLep0'                : Plot1D( bin_range=[0.0, -1],     nbins=25, xunits='GeV', xlabel='M_{l0}'),
    'MLep1'                : Plot1D( bin_range=[0.0, -1],     nbins=25, xunits='GeV', xlabel='M_{l1}'),
    'DEtaLL'               : Plot1D( bin_range=[0.0, 6.0],    nbins=20, xlabel='#Delta#eta_{ll}'),
    'DphiLL'               : Plot1D( bin_range=[0.0, 3.15],   nbins=30, xlabel='#Delta#phi_{ll}'),
    'drll'                 : Plot1D( bin_range=[0.0, 6.0],    nbins=60, xlabel='#DeltaR_{ll}'),
    'dRy_sEl_bMu_noCalo'   : Plot1D( bin_range=[0.0, 6.0],    nbins=60, xlabel='#DeltaR_{sig E, base nonCalo Mu}'),
    'dRy_sEl_bMu_Calo'     : Plot1D( bin_range=[0.0, 6.0],    nbins=60, xlabel='#DeltaR_{sig E, base Calo Mu}'),
    'isCaloTagged'         : Plot1D( bin_range=[-1.5, 3.5],   bin_width=1, xlabel='Calo-Tagged Muon'),
    'dilep_flav'           : Plot1D( bin_range=[-0.5, 4.5],   bin_width=1, xlabel='Dilepton flavor'),
    'isEM'                 : Plot1D( bin_range=[-1.5, 3.5],   bin_width=1, xlabel='Dilepton flavor is el mu'),
    'isME'                 : Plot1D( bin_range=[-1.5, 3.5],   bin_width=1, xlabel='Dilepton flavor is mu el'),
    'MtLep0'               : Plot1D( bin_range=[0.0, 250.0],  nbins=15, xunits='GeV', xlabel='m_{T}(l_{0},MET)'),
    'MtLep1'               : Plot1D( bin_range=[0.0, 140.0],  nbins=20, xunits='GeV', xlabel='m_{T}(l_{1},MET)'),
    #'MLL'                 : Plot1D( bin_range=[50.0, 130],   nbins=80, xunits='GeV', xlabel='M_{ll}'),
    'MLL'                  : Plot1D( bin_range=[0.0, 300.0],  nbins=100, xunits='GeV', xlabel='M_{ll}'),
    'ptll'                 : Plot1D( bin_range=[0.0, 500.0],  nbins=50, xunits='GeV', xlabel='pT_{ll}'),
    # MET + leptons
    'MET'                  : Plot1D( bin_range=[0.0, 200.0],  bin_width=4, xunits='GeV', xlabel='E_{T}^{miss}'),
    'METPhi'               : Plot1D( bin_range=[0.0, 3.15],   nbins=30, xlabel='MET_{#phi}'),
    'MCollASym'            : Plot1D( bin_range=[0.0, 250.0],  nbins=25, xunits='GeV', xlabel='LFV Collinear Mass m_{coll}'),
    'dpt_ll'               : Plot1D( bin_range=[0.0, 150.0],  nbins=20, xunits='GeV', xlabel='#Deltap_{T}^{ll}'),
    'DphiLep0MET'          : Plot1D( bin_range=[-3.15, 3.15], nbins=63, add_underflow=True, xlabel='#Delta#phi(l_{0},MET)'),
    'DphiLep1MET'          : Plot1D( bin_range=[-3.15, 3.15], nbins=63, add_underflow=True, xlabel='#Delta#phi(l_{1},MET)'),
    'tau_pT'               : Plot1D( bin_range=[0.0, 200.0],  nbins=20, xunits='GeV', xlabel='p_{T}^{subleading lep + MET}'),
    'taulep1_pT_ratio'     : Plot1D( bin_range=[0.0, 3],      nbins=20, xlabel='p_{T}^{subleading lep + MET} / p_{T}^{leading lep}'),
    # Jets
    'preJet_pt'            : Plot1D( bin_range=[0.0, 100.0],  nbins=20, xunits='GeV', xlabel='p_{T}^{preJet}'),
    'preJet_eta'           : Plot1D( bin_range=[-5.0, 5.0],   nbins=100, xlabel='#eta_{preJet}'),
    'preJet_JVT'           : Plot1D( bin_range=[-0.2, 1.1],   nbins=39, xlabel='preJet JVT (|eta|<=2.4 & pT < 60)'),
    'baseJet_eta'          : Plot1D( bin_range=[-5.0, 5.0],   nbins=100, xlabel='#eta_{baseJet}'),
    'baseJet_mv2c10'       : Plot1D( bin_range=[-2, 2],       nbins=80, xlabel='mv2c10_{baseJet}'),
    'jetN'                 : Plot1D( bin_range=[-0.5, 7.5],   bin_width=1, xlabel='N_{base jet}'),
    'jet_N2p4Eta25Pt'      : Plot1D( bin_range=[-0.5, 7.5],   bin_width=1, xlabel='N_{jet} (p_{T}>25GeV, |#eta|<2.5)'),
    'signalJetN'           : Plot1D( bin_range=[-0.5, 7.5],   bin_width=1, xlabel='N_{sig jets}'),
    'jetN_g30'             : Plot1D( bin_range=[-0.5, 7.5],   bin_width=1, xlabel='N_{jet} (p_{T}>30GeV)'),
    'nLJets'               : Plot1D( bin_range=[-0.5, 7.5],   bin_width=1, xlabel='N_{light jet}'),
    'nBJets'               : Plot1D( bin_range=[-0.5, 7.5],   bin_width=1, xlabel='N_{Bjet}'),
    'btag'                 : Plot1D( bin_range=[-1.5, 3.5],   bin_width=1, xlabel='B-tagged jet'),
    'nForwardJets'         : Plot1D( bin_range=[-0.5, 7.5],   bin_width=1, xlabel='N_{F jet}'),
    'j_pt[0]'              : Plot1D( bin_range=[0.0, 500.0],  bin_width=20, xunits='GeV', xlabel='p_{T}^{leading jet}'),
    'j_pt[1]'              : Plot1D( bin_range=[0.0, 500.0],  nbins=20, xunits='GeV', xlabel='p_{T}^{subleading jet}'),
    'j_pt[2]'              : Plot1D( bin_range=[0.0, 500.0],  nbins=20, xunits='GeV', xlabel='p_{T}^{3rd leading jet}'),
    'j_pt[3]'              : Plot1D( bin_range=[0.0, 500.0],  nbins=20, xunits='GeV', xlabel='p_{T}^{4th leading jet}'),
    'j_pt'                 : Plot1D( bin_range=[0.0, 500.0],  nbins=20, xunits='GeV', xlabel='Jet p_{T}'),
    'mjj'                  : Plot1D( bin_range=[0.0, -1],     nbins=25, xunits='GeV', xlabel='Dijet mass'),
    'dEtaJJ'               : Plot1D( bin_range=[0.0, 6.0],    nbins=20, xlabel='#Delta#eta(j0,j1)'),
    'j_eta'                : Plot1D( bin_range=[-5.0, 5.0],   nbins=100, xlabel='Jet #eta'),
    'j_jvt'                : Plot1D( bin_range=[0.0, -1],     nbins=25, xlabel='Jet JVT'),
    'j_jvf'                : Plot1D( bin_range=[0.0, -1],     nbins=25, xlabel='Jet JVF'),
    'j_phi'                : Plot1D( bin_range=[0.0, 3.15],   nbins=30, xlabel='Jet #phi'),
    'j_flav'               : Plot1D( bin_range=[-0.5, 4.5],   nbins=5, xlabel='Jet flavor (0:NA,1:CL,2:CB,3:F)'),
    # Leptons
    'preEl_EcaloClus'      : Plot1D( bin_range=[ 0, 50],      nbins=50,  add_overflow=False, xunits='GeV', xlabel='E_{preEl CaloCluster}'),
    'preEl_etaCaloClus'    : Plot1D( bin_range=[-3.0, 3.0],   nbins=20, xlabel='#eta_{preEl CaloCluster}'),
    'baseEl_etconetopo20'  : Plot1D( bin_range=[ -7, 15],     nbins=50,  doLogY=False, xunits='GeV', xlabel='E_{T}^{baselineEl} conetopo20' ),
    'baseEl_ptvarcone20'   : Plot1D( bin_range=[ -1, 4],      nbins=25,  doLogY=False, xunits='GeV', xlabel='p_{T}^{baselineEl} ptvarcone20' ),
    'baseMu_etconetopo20'  : Plot1D( bin_range=[ -3, 3],      nbins=60,  xunits='GeV', xlabel='E_{T}^{baselineMu} conetopo20' ),
    'baseMu_ptvarcone30'   : Plot1D( bin_range=[ -3, 3],      nbins=60,  doLogY=False, xunits='GeV', xlabel='p_{T}^{baselineMu} ptvarcone30' ),
    'el1pT_trackclus_ratio': Plot1D( bin_range=[0, 3],        nbins=30, xlabel='el_{subleading pT}^{track} / el_{subleading pT}^{cluster}'),
    # Fakes
    'nLepID'               : Plot1D( bin_range=[-1.5, 7.5],   bin_width=1, xlabel='ID lepton multiplicity'),
    'nLepAntiID'           : Plot1D( bin_range=[-1.5, 7.5],   bin_width=1, xlabel='anti-ID lepton multiplicity'),
    'aID_Lep0Pt'           : Plot1D( bin_range=[0.0, 200.0],  nbins=40, xunits='GeV', xlabel='p_{T}^{leading lep} anti-ID'),
    'aID_Lep1Pt'           : Plot1D( bin_range=[0.0, 100.0],  nbins=20, xunits='GeV', xlabel='p_{T}^{subleading lep} anti-ID'),
    'aID_Lep0Eta'          : Plot1D( bin_range=[-3.0, 3.0],   nbins=20, xlabel='#eta^{leading lep} anti-ID'),
    'aID_Lep1Eta'          : Plot1D( bin_range=[-3.0, 3.0],   nbins=20, xlabel='#eta^{subleading lep} anti-ID'),
    'aID_MLL'              : Plot1D( bin_range=[0.0, 300.0],  nbins=100, xunits='GeV', xlabel='M_{ll}(ID, anti-ID)'),
    'aID_Lep0Flav'         : Plot1D( bin_range=[-1.5, 5.5],   bin_width=1, xlabel='Leading antiID flavor'),
    'aID_Lep1Flav'         : Plot1D( bin_range=[-1.5, 5.5],   bin_width=1, xlabel='Subleading antiID flavor'),
    'aID_Lep0Q'            : Plot1D( bin_range=[-1.5, 1.5],   bin_width=1, xlabel='Leading anti-ID charge'),
    'aID_Lep1Q'            : Plot1D( bin_range=[-1.5, 1.5],   bin_width=1, xlabel='Subleading anti-ID charge'),
    'Z_MLL'                : Plot1D( bin_range=[65, 115.0],   bin_width=1, xunits='GeV', xlabel='M_{ll} (Z pair)'),
    'Z2_MLL'               : Plot1D( bin_range=[0, 200.0],    bin_width=5, logy_max=1e6, xunits='GeV', xlabel='M_{ll} (2nd Z pair)'),
    'aID_dpt_ll'           : Plot1D( bin_range=[0.0, 150.0],  nbins=20, xunits='GeV', xlabel='#Deltap_{T}^{ll}(ID, anti-ID)'),
    'aID_drll'             : Plot1D( bin_range=[0.0, 6.0],    nbins=60, xlabel='#DeltaR_{ll}(ID, anti-ID)'),
    'dR_Zl'                : Plot1D( bin_range=[0.0, 6.0],    nbins=60, xlabel='#DeltaR(Z, lep)'),
    'Z_dilep_flav'         : Plot1D( bin_range=[-1.5, 5.5],   nbins=7, xlabel='Z dilepton flavor'),
    'Z2_dilep_flav'        : Plot1D( bin_range=[-1.5, 5.5],   nbins=7, xlabel='2nd Z dilepton flavor'),
    'Z_Lep2_pT'            : Plot1D( bin_range=[0.0, 200.0],  bin_width=5, xunits='GeV', xlabel='Fake candidate lepton p_{T}'),
    'Z_Lep2_eta'           : Plot1D( bin_range=[-3.0, 3.0],   nbins=20, xlabel='Fake candidate lepton #eta'),
    'Z_Lep2_flav'          : Plot1D( bin_range=[-1.5, 2.5],   bin_width=1, xlabel='Fake candidate flavor'),
    'Z_dilep_sign'         : Plot1D( bin_range=[-2.5, 2.5],   bin_width=1, xlabel='Z Dilepton Sign : OS(-1) SS(1)'),
    'Z2_dilep_sign'        : Plot1D( bin_range=[-2.5, 2.5],   bin_width=1, xlabel='2nd Z Dilepton Sign : OS(-1) SS(1)'),
    'Z_Lep2_dPhi_MET'      : Plot1D( bin_range=[-3.15, 3.15], nbins=63, add_underflow=True, xlabel='#Delta#phi(l_{3},MET)'),
    'Z_Lep2_mT'            : Plot1D( bin_range=[0.0, 200.0],  bin_width=4, xunits='GeV', xlabel='3rd leading lepton m_{T}'),
}

# Add any labels to the plots
l_truthClass_labels = [ 'Undefined',' ','prompt El','prompt Mu', 'prompt Pho', 'Bkg El from conv','prompt El from FSR','hadron','Mu as e','HF tau','HF B','HF C',' ']
plot_defaults['l_truthClass[0]'].bin_labels = l_truthClass_labels
plot_defaults['l_truthClass[1]'].bin_labels = l_truthClass_labels
plot_defaults['l_truthClass[2]'].bin_labels = l_truthClass_labels
plot_defaults['isMC'].bin_labels = [' ', 'Data', 'MC', ' ']
plot_defaults['Z_Lep2_flav'].bin_labels = [' ', 'Muon', 'Electron', ' ']


# To alter the plot properties for a specific region
# Deep copy the default plot into new plot dictionary
# Edit that copy as needed

region_plots = {}
#region_plots['zjets_FF_CRden_eem'] = {
#  'Z_MLL' : deepcopy(plot_defaults['Z_MLL'])
#}
#region_plots['zjets_FF_CRden_eem']['Z_MLL'].update(bin_range=[80, 100], bin_width=1)

### Define variables for plots


################################################################################
# Toggle options for execution
# - Options expected to change often between different plottings
################################################################################

#######################################
## Build the TChain/TTree for each sample
# To remove sample from plot, comment out the line setting its TChain
# Samples with empty TChains get removed below
data_dsids = g.groups['data15']+g.groups['data16']
data.set_chain_from_dsid_list(data_dsids, data_ntuple_dir, exclude_strs=['FFest'])
#fakes.set_chain_from_dsid_list(data_dsids, fake_ntuple_dir, search_strs=['FFest'])
ttbar.set_chain_from_dsid_list(g.groups['ttbar'], bkg_ntuple_dir)
stop.set_chain_from_dsid_list(g.groups['singletop'], bkg_ntuple_dir)
wtop.set_chain_from_dsid_list(g.groups['Wt'], bkg_ntuple_dir)
VV.set_chain_from_dsid_list(g.groups['ggllvv'], bkg_ntuple_dir)
WW.set_chain_from_dsid_list(g.groups['ww'], bkg_ntuple_dir)
ZZ.set_chain_from_dsid_list(g.groups['zz'], bkg_ntuple_dir)
WZ.set_chain_from_dsid_list(g.groups['wz'], bkg_ntuple_dir)
#zll.set_chain_from_dsid_list(g.groups['zll'], bkg_ntuple_dir)
zee.set_chain_from_dsid_list(g.groups['zee'], bkg_ntuple_dir)
zmumu.set_chain_from_dsid_list(g.groups['zmumu'], bkg_ntuple_dir)
ztt.set_chain_from_dsid_list(g.groups['ztt'], bkg_ntuple_dir)
wjets.set_chain_from_dsid_list(g.groups['wjets'], bkg_ntuple_dir)
#wgamma.set_chain_from_dsid_list(g.groups['w_gamma'], bkg_ntuple_dir)
htt.set_chain_from_dsid_list(g.groups['htt'], bkg_ntuple_dir)
hww.set_chain_from_dsid_list(g.groups['hww'], bkg_ntuple_dir)
#signal.set_chain_from_dsid_list(g.groups['higgs_lfv'], signal_ntuple_dir)
SAMPLES = [s for s in SAMPLES if s.is_setup()]
assert SAMPLES, "ERROR :: No samples are setup"

#######################################
# Yield Table
YIELD_TBL = YieldTable()
# Add formulas to yield table
# Key values will be the labels on the table
# Formulas should use sample names and these symbols: +, -, *, /, (, ), [0-9]
YIELD_TBL.formulas['(WZ+ZZ)/Data'] = "(wz + zz)/data"
YIELD_TBL.formulas['Zll/Data'] = "(zee + zmumu)/data"

#######################################
# What regions to plot
region_ops = []
#region_ops += ['zjets_FF_CRden_e', 'zjets_FF_CRnum_e']
#region_ops += ['zjets_FF_CRden_m', 'zjets_FF_CRnum_m']
#region_ops += ['zjets_FF_CRden_eem']
#region_ops += ['zjets_FF_CRden_eem', 'zjets_FF_CRden_mmm']
#region_ops += ['zjets_FF_CRden_eee', 'zjets_FF_CRden_mme']
#region_ops += ['zjets_FF_CRnum_eem', 'zjets_FF_CRnum_mmm']
#region_ops += ['zjets_FF_CRnum_eee', 'zjets_FF_CRnum_mme']
#region_ops += ['zCR_ee', 'zCR_mumu']

#######################################
# What variables to plot
vars_to_plot = []
#vars_to_plot += ['nLepID']
#vars_to_plot += ['nLepAntiID']
#vars_to_plot += ['Z_Lep2_flav']
#vars_to_plot += ['Z_MLL']
#vars_to_plot += ['Z_MLL','nBJets','Z_Lep2_mT','Z2_MLL','MET']
vars_to_plot += ['l_pt[0]','l_pt[1]','Z_Lep2_pT']
vars_to_plot += ['l_truthClass[0]','l_truthClass[1]','l_truthClass[2]']

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
            assert n_data <= 1, "ERROR :: More than one data sample setup"

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
