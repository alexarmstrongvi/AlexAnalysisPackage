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

# Root data analysis framework
import ROOT
# Prevent root from printing anything to the screen
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# Local classes for plotting
import tools.plot as plot
from tools.sample import Data, Background, Signal
import tools.region as region
import tools.systematic as systematic


################################################################################
# Useful globals
################################################################################
import global_variables as g
# Path to directories with flat ntuples
#TODO: rawdir -> mc_ntuples_dir, signal_ntuples_dir, data_ntuples_dir
rawdir        = g.mc_ntuples
signal_rawdir = g.mc_ntuples
data_rawdir   = g.data_ntuples

# Luminosity options
# Description | 2015-16 (ipb) |  2015-16 (ifb) | 2015 (ipb) | 2016 (ipb) | dummy
# Value       |     36180     |      36.01     |    320     |    32971   |   1
lumi_ = 36180

################################################################################
# Samples
# - Create all the samples that will be put into plots (e.g. backgrounds,
#   signal, data)
################################################################################
# Data
data = Data()
data.treename = "data"

################################################################################
# Fakes
## Fakes
fakes = Background("fakes", "Fakes")
fakes.color = ROOT.kGray
fakes.treename = "fakes"

################################################################################
# Signal
signal_branching_ratio = 0.01
signal_SF = 1
signal_label = "Higgs LFV" if signal_SF == 1 else "Higgs LFV (%dX)"%signal_SF

signal = Signal("higgs_lfv", signal_label)
signal.scale_factor = lumi_ * signal_branching_ratio * signal_SF
signal.color = ROOT.kGreen
signal.treename = "signal"

################################################################################
# Backgrounds

#######################################
# Initialize all backgrounds
# ttbar
ttbar = Background("ttbar", "t#bar{t}")
ttbar.scale_factor = lumi_
ttbar.color = ROOT.kOrange+2
ttbar.treename = "ttbar"

# singletop
stop = Background("st", "Single top")
stop.scale_factor = lumi_
stop.color = ROOT.kOrange+1
stop.treename = "ST"

# W+top
wtop = Background("wt", "Wt")
wtop.scale_factor = lumi_
wtop.color = ROOT.kOrange+8
wtop.treename = "wt"

# WW
WW = Background("ww", "WW")
WW.scale_factor = lumi_
WW.color = ROOT.kSpring-6
WW.treename = "WW"

# ZZ
ZZ = Background("zz", "ZZ")
ZZ.scale_factor = lumi_
ZZ.color = ROOT.kSpring-4
ZZ.treename = "ZZ"

# WZ
WZ = Background("wz", "WZ")
WZ.scale_factor = lumi_
WZ.color = ROOT.kSpring-5
WZ.treename = "WZ"

# Zll
zll = Background("zll", "Zll")
zll.scale_factor = lumi_
zll.color = ROOT.kAzure-9
zll.treename = "zll"

# Zee
zee = Background("zee", "Zee")
zee.scale_factor = lumi_
zee.color = ROOT.kAzure-7
zee.treename = "zee"

# Zmumu
zmumu = Background("zmumu", "Zmumu")
zmumu.scale_factor = lumi_
zmumu.color = ROOT.kAzure-9
zmumu.treename = "zmumu"

# Ztt
ztt = Background("ztt", "Z#tau#tau")
ztt.scale_factor = lumi_
ztt.color = ROOT.kAzure-5
ztt.treename = "ztt"

# Wjets
wjets = Background("wjets", "W+jets")
wjets.scale_factor = lumi_
wjets.color = ROOT.kOrange
wjets.treename = "wjets"

# W+gamma
wgamma = Background("wgamma", "W+gamma")
wgamma.scale_factor = lumi_
wgamma.color = ROOT.kOrange-1
wgamma.treename = "w_gamma"

# Higgs -> tau tau
htt = Background("htt", "H#tau#tau")
htt.scale_factor = lumi_
htt.color = ROOT.kRed
htt.treename = "htt"

# Higgs -> W W
hww = Background("hww", "HWW")
hww.scale_factor = lumi_
hww.color = ROOT.kBlue+3
hww.treename = "hww"


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

# Combined triggers
e15_trig_pT = '(%s || %s)'%(e15_trig_emu_pT, e15_trig_mue_pT)
mu15_trig_pT = '(%s || %s)'%(mu15_trig_emu_pT, mu15_trig_mue_pT)
e16_trig_pT = '(%s || %s)'%(e16_trig_emu_pT, e16_trig_mue_pT)
mu16_trig_pT = '(%s || %s)'%(mu16_trig_emu_pT, mu16_trig_mue_pT)
emu15_trig_pT = '(%s || %s)'%(emu15_trig_emu_pT, emu15_trig_mue_pT)
emu16_trig_pT = '(%s || %s)'%(emu16_trig_emu_pT, emu16_trig_mue_pT)

dilep15_trig      = '(treatAsYear==2015 && %s)'%emu15_trig
dilep16_trig     = '(treatAsYear==2016 && %s)'%emu16_trig
singlelep15_trig = '(treatAsYear==2015 && (%s || %s))'%(e15_trig,mu15_trig)
singlelep16_trig = '(treatAsYear==2016 && (%s || %s))'%(e16_trig,mu16_trig)
dilep_trig = '(%s || %s)'%(dilep15_trig, dilep16_trig)
singlelep_trig = '(%s || %s)'%(singlelep15_trig, singlelep16_trig)

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
regions = []

reg = region.Region()
reg.name = "no_sel"
reg.displayname = "No Selections"
reg.tcut = '1' #DF_OS
regions.append(reg)

reg = region.Region()
reg.name = "trig_only"
reg.displayname = "Lepton Triggers"
reg.tcut = lepton_trig_pT + "&&" + DF_OS
regions.append(reg)

reg = region.Region()
reg.name = "topCR"
reg.displayname = "Top CR"
reg.tcut = "nBJets >= 1 && MET > 40 &&" + DF_OS + " &&" + lepton_trig_pT
regions.append(reg)

reg = region.Region()
reg.name = "zCR_ee"
reg.displayname = "Z CR (Channel: El-El)"
reg.tcut = "75 < MLL && MLL < 105 && " + SF_OS + " && " + ee + " && " + lepton_trig_pT
regions.append(reg)

reg = region.Region()
reg.name = "zCR_mumu"
reg.displayname = "Z CR (Channel: Mu-Mu)"
reg.tcut = "75 < MLL && MLL < 105 && " + SF_OS + " && " + mumu + " && " + lepton_trig_pT
regions.append(reg)

ZTauTau_CR =  ('Lep0Pt >= 30 && Lep0Pt < 45 && Lep1Pt >= 15 '
              + '&& (30 < MLL && MLL < 150) '
              + '&& nBJets==0 '
              + '&& ( !'+mue+' || el1pT_trackclus_ratio < 1.2) '
              + '&&' + DF_OS)
reg = region.Region()
reg.name = "ztautauCR"
reg.displayname = "Ztautau CR"
reg.tcut = ZTauTau_CR
regions.append(reg)

# Z+Jets fake regions
zjets_FF_CRden_base = 'nLepID == 2 && nLepAntiID >= 1'
zjets_FF_CRnum_base = 'nLepID == 3'
zjets_FF_CR_add = '70 < Z_MLL && Z_MLL < 110'
zjets_FF_CR_add += ' && nBJets == 0'
zjets_FF_CR_add += ' && MET < 50'
zjets_FF_CR_add += ' && Z_Lep2_mT < 50'
zjets_FF_CR_add += ' && (Z2_MLL < 80 || 100 < Z2_MLL)'
#zjets_FF_CR_add += ' && !isMC || (l_truthClass[2] == 1 || l_truthClass[2] == 2)' # Only prompt leptons
#zjets_FF_CR_add += ' && !isMC || !(l_truthClass[2] == 1 || l_truthClass[2] == 2)'
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

        reg = region.Region()
        reg.name = 'zjets_FF_CR%s_%s'%(num_den, chan)
        reg.displayname = 'Z+jets FF CR (%s)'%(chan_name)
        reg.tcut = ' && '.join([num_den_sel, zjets_FF_CR_add, ops[2], singlelep_trig])
        regions.append(reg)

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

reg = region.Region()
reg.name = 'wjets_FF_CRden_emu'
reg.displayname = 'wjets FF CR (anti-ID el)'
reg.tcut = wjets_FF_CRden_base + '&&' + wjets_FF_CRden_add + '&&' + FF_CRden_emu_ch
regions.append(reg)

reg = region.Region()
reg.name = 'wjets_FF_CRden_mue'
reg.displayname = 'wjets FF CR (anti-ID mu)'
reg.tcut = wjets_FF_CRden_base + '&&' + wjets_FF_CRden_add + '&&' + FF_CRden_mue_ch
regions.append(reg)

wjets_FF_CRnum_base = 'nLepID == 2' # require final state
wjets_FF_CRnum_add = '1'
wjets_FF_CRnum_add += ' && l_q[0] == l_q[1]' # reject zjets
wjets_FF_CRnum_add += ' && (MLL < 70 || 110 < MLL)' # reject zjets
wjets_FF_CRnum_add += ' && nBJets==0' # reject top
FF_CRnum_emu_ch  = 'dilep_flav == 0'
FF_CRnum_mue_ch  = 'dilep_flav == 1'
FF_CRnum_ee_ch   = 'dilep_flav == 2'
FF_CRnum_mumu_ch = 'dilep_flav == 3'

reg = region.Region()
reg.name = 'wjets_FF_CRnum_emu'
reg.displayname = 'wjets FF CR (ID el)'
reg.tcut = wjets_FF_CRnum_base + '&&' + wjets_FF_CRnum_add + '&&' + FF_CRnum_emu_ch
regions.append(reg)

reg = region.Region()
reg.name = 'wjets_FF_CRnum_mue'
reg.displayname = 'wjets FF CR (ID mu)'
reg.tcut = wjets_FF_CRnum_base + '&&' + wjets_FF_CRnum_add + '&&' + FF_CRnum_emu_ch
regions.append(reg)

# Baseline regions
reg = region.Region()
reg.name = "baseline"
reg.displayname = "Baseline"
reg.tcut = Baseline_Sel
regions.append(reg)

reg = region.Region()
reg.name = "baseline_emu"
reg.displayname = "Baseline (Channel: El-Mu)"
reg.tcut = Baseline_Sel + " && " + emu
regions.append(reg)

reg = region.Region()
reg.name = "baseline_mue"
reg.displayname = "Baseline (Channel: Mu-El)"
reg.tcut = Baseline_Sel + " && " + mue
regions.append(reg)

reg = region.Region()
reg.name = "vbf"
reg.displayname = "VBF"
reg.tcut = "(%s) && (%s)"%(Baseline_Sel, VBF_stripped)
regions.append(reg)

reg = region.Region()
reg.name = "optimized"
reg.displayname = "Optimized"
reg.tcut = "(%s) && !(%s) && DphiLep1MET < 1 && MtLep0 > 50 && MtLep1 < 40 && ((MET+Lep1Pt)/Lep1Pt) > 0.5"%(Baseline_Sel, VBF_stripped)
regions.append(reg)

################################################################################
# Variables
# - Create all the variables that one might like to plot
################################################################################
### Define variables for plots
HistOp1D = namedtuple('HistOp1D', 'regions,   nBinsX, x0, x1, y0,   y1,   logY, norm   xUnits, xLabel, yLabel, add_overflow, add_underflow, ratioPlot')
HistOp1D.__new__.__defaults__=   ('baseline', 12,     -2,  10, None, None, False, False, '',     '',     'Events', True, False, True)

HistOpMap = {
    # Event level
    'RunNumber'       : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Event run number'),
    'event_number'    : HistOp1D(nBinsX=1000, x0=0,  x1=1000000000,    xUnits='',    xLabel='Event number'),
    'isMC'            : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5,   xUnits='',    xLabel='is Monte Carlo'),
    'eventweight'     : HistOp1D(nBinsX=100, x0=-0.001,  x1=0.002,    xUnits='',    xLabel='Event weight', logY=True),
    'dsid'            : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Sample DSID'),
    'treatAsYear'     : HistOp1D(nBinsX=11, x0=2007.5,x1=2018.5,xUnits='',   xLabel='treatAsYear'),
    # Multiplicity
    'n_preLeptons'    : HistOp1D(nBinsX=11, x0=-0.5, x1=10.5,  xUnits='',    xLabel='N_{pre-leptons}'),
    'n_baseLeptons'   : HistOp1D(nBinsX=11, x0=-0.5, x1=10.5,  xUnits='',    xLabel='N_{baseline leptons}'),
    'n_leptons'       : HistOp1D(nBinsX=11, x0=-0.5, x1=10.5,  xUnits='',    xLabel='N_{signal leptons}'),
    # Leptons
    'l_pt'            : HistOp1D(nBinsX=25, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='Lepton p_{T}'),
    'aID_Lep0Pt'          : HistOp1D(nBinsX=40, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{leading lep} anti-ID', logY=True),
    ## Electrons
    'el1_track_pt'    : HistOp1D(nBinsX=25, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='Leading electron track p_{T}'),
    'el1_clus_pt'     : HistOp1D(nBinsX=25, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='Leading electron cluster p_{T}'),
    'preEl_Et'        : HistOp1D(nBinsX=50, x0=0.0,  x1=50.0, xUnits='GeV', xLabel='E_{T}^{preElectron}', add_overflow=False),
    'preEl_pt'        : HistOp1D(nBinsX=50, x0=0.0,  x1=50.0, xUnits='GeV', xLabel='p_{T}^{preElectron}'),
    'preEl_clusEtaBE' : HistOp1D(nBinsX=60, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta_{preEl clusterBE}'),
    'preEl_eta'       : HistOp1D(nBinsX=60, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta_{preElectron}'),
    'baseEl_ID'       : HistOp1D(nBinsX=9,  x0=-0.5, x1=8.5,   xUnits='',    xLabel='el_{Baseline}^{ID}(non-inclusive)', logY=True),
    'El_ID'           : HistOp1D(nBinsX=7,  x0=-1.5, x1=5.5,   xUnits='',    xLabel='electron ID(non-inclusive)', logY=False),
    'el_type'         : HistOp1D(nBinsX=41,  x0=-1.5, x1=39.5,   xUnits='',    xLabel='electron truth type', logY=True),
    'el_origin'       : HistOp1D(nBinsX=48,  x0=-1.5, x1=46.5,   xUnits='',    xLabel='electron truth origin', logY=True),
    'el_d0sigBSCorr'  : HistOp1D(nBinsX=60,  x0=-6, x1=6,   xUnits='mm',    xLabel='Electron d_{0}/#sigma_{d_{0}} BSCorr', add_underflow=True),
    'el_z0SinTheta'  : HistOp1D(nBinsX=60,  x0=-0.6, x1=0.6,   xUnits='mm',    xLabel='Electron z_{0}sin(#theta)', add_underflow=True),
    ## Muons
    'preMu_pt'        : HistOp1D(nBinsX=50, x0=0,  x1=50.0, xUnits='GeV', xLabel='p_{T}^{preMuon}', add_overflow=False),
    'preMu_ID'      : HistOp1D(nBinsX=7,  x0=-1.5, x1=5.5,   xUnits='',    xLabel='mu_{pre}^{ID} (non-inclusive)'),
    'baseMu_pt'       : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{BaselineMuon}'),
    'baseMu_eta'      : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta_{BaselineMuon}'),
    'baseMu_ID'     : HistOp1D(nBinsX=7,  x0=-1.5, x1=5.5,   xUnits='',    xLabel='mu_{Baseline}^{ID}(non-inclusive)', logY=True),
    'Mu_ID'           : HistOp1D(nBinsX=7,  x0=-1.5, x1=5.5,   xUnits='',    xLabel='muon ID (non-inclusive)', logY=False),
    'mu_type'         : HistOp1D(nBinsX=41,  x0=-1.5, x1=39.5,   xUnits='',    xLabel='muon truth type', logY=True),
    'mu_origin'       : HistOp1D(nBinsX=48,  x0=-1.5, x1=46.5,   xUnits='',    xLabel='muon truth origin', logY=True),
    'mu_d0sigBSCorr'  : HistOp1D(nBinsX=60,  x0=-6, x1=6,   xUnits='mm',    xLabel='Muon d_{0}/#sigma_{d_{0}} BSCorr', add_underflow=True),
    'mu_z0SinTheta'  : HistOp1D(nBinsX=60,  x0=-0.6, x1=0.6,   xUnits='mm',    xLabel='Muon z_{0}sin(#theta)', add_underflow=True),
    ## Taus
    'preTau_q'        : HistOp1D(nBinsX=33,  x0=-5.5, x1=5.5,   xUnits='',    xLabel='Tau charge'),
    'preTau_nTracks'  : HistOp1D(nBinsX=10,  x0=-1.5, x1=8.5,   xUnits='',    xLabel='preTau nTracks'),
    'baseTau_pT'      : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{BaselineTau}'),
    'baseTau_eta'     : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta_{BaselineMuon}'),
    'baseTau_nTracks' : HistOp1D(nBinsX=10,  x0=-1.5, x1=8.5,   xUnits='',    xLabel='Baseline Tau nTracks'),
    'baseTau_ID'      : HistOp1D(nBinsX=7,  x0=-1.5, x1=5.5,   xUnits='',    xLabel='tau_{Baseline}^{ID}(non-inclusive)', logY=True),

    ## General lepton
    'l_eta'           : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='Lepton #eta'),
    'l_phi'           : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='Lepton #phi'),
    'l_flav'          : HistOp1D(nBinsX=5,  x0=-0.5, x1=4.5,   xUnits='',    xLabel='Lepton flavor (0: e, 1: m)'),
    'l_type'          : HistOp1D(nBinsX=41,  x0=-1.5, x1=39.5,   xUnits='',    xLabel='Lepton type', logY=True, ratioPlot=False),
    'l_origin'        : HistOp1D(nBinsX=48,  x0=-1.5, x1=46.5,    xUnits='',    xLabel='Lepton origin', logY=True, ratioPlot=False),
    'l_BkgMotherPdgId' : HistOp1D(nBinsX=41,  x0=-20.5, x1=20.5,    xUnits='',    xLabel='Lepton Mother PdgID', logY=True, add_underflow=True, ratioPlot=False),
    'l_truthClass'    : HistOp1D(nBinsX=12,  x0=-1.5, x1=10.5,    xUnits='',    xLabel='Lepton Truth Classification', logY=True, norm=True, ratioPlot=False),
    'l_type[0]'       : HistOp1D(nBinsX=41,  x0=-1.5, x1=39.5,   xUnits='',    xLabel='Leading Z Lepton type', logY=True, ratioPlot=False),
    'l_origin[0]'     : HistOp1D(nBinsX=48,  x0=-1.5, x1=46.5,    xUnits='',    xLabel='Leading Z Lepton origin', logY=True, ratioPlot=False),
    'l_type[1]'       : HistOp1D(nBinsX=41,  x0=-1.5, x1=39.5,   xUnits='',    xLabel='Subleading Z Lepton type', logY=True, ratioPlot=False),
    'l_origin[1]'     : HistOp1D(nBinsX=48,  x0=-1.5, x1=46.5,    xUnits='',    xLabel='Subleading Z Lepton origin', logY=True, ratioPlot=False),
    'l_type[2]'       : HistOp1D(nBinsX=41,  x0=-1.5, x1=39.5,   xUnits='',    xLabel='Fake Candidate Lepton type',  logY=True, ratioPlot=False),
    'l_origin[2]'     : HistOp1D(nBinsX=48,  x0=-1.5, x1=46.5,    xUnits='',    xLabel='Fake Candidate Lepton origin', logY=True, ratioPlot=False),
    'l_BkgMotherPdgId[2]': HistOp1D(nBinsX=41,  x0=-20.5, x1=20.5,    xUnits='',    xLabel='Fake Candidate Lepton Mother PdgID', logY=True, norm=True, add_underflow=True, ratioPlot=False),
    'Lep_Iso'         : HistOp1D(nBinsX=9,  x0=-1.5, x1=7.5,   xUnits='',    xLabel='Lepton Isolation (non-inclusive)'),
    'l_q'             : HistOp1D(nBinsX=3,  x0=-1.5, x1=1.5,   xUnits='',    xLabel='Lepton charge'),
    'LepLepSign'      : HistOp1D(nBinsX=3,  x0=-1.5, x1=1.5,   xUnits='',    xLabel='Leptons sign product'),
    'Lep0Pt'          : HistOp1D(nBinsX=40, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{leading lep}', logY=True),
    #'Lep1Pt'          : HistOp1D(nBinsX=40, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{subleading lep}', logY=True),
    #'Lep0Pt'          : HistOp1D(nBinsX=50, x0=0.0,  x1=50.0, xUnits='GeV', xLabel='p_{T}^{leading lep}', add_overflow=False),
    'Lep1Pt'          : HistOp1D(nBinsX=50, x0=0.0,  x1=50.0, xUnits='GeV', xLabel='p_{T}^{subleading lep}', add_overflow=False),
    'Lep0Eta'         : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta^{leading lep}'),
    'Lep1Eta'         : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta^{subleading lep}'),
    'Lep0Phi'         : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='#phi^{leading lep}'),
    'Lep1Phi'         : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='#phi^{subleading lep}'),
    'MLep0'           : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='GeV', xLabel='M_{l0}'),
    'MLep1'           : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='GeV', xLabel='M_{l1}'),
    'DEtaLL'          : HistOp1D(nBinsX=20, x0=0.0,  x1=6.0,   xUnits='',    xLabel='#Delta#eta_{ll}'),
    'DphiLL'          : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='#Delta#phi_{ll}'),
    'drll'            : HistOp1D(nBinsX=60, x0=0.0,  x1=6.0,   xUnits='',    xLabel='#DeltaR_{ll}'),
    'dRy_sEl_bMu_noCalo' : HistOp1D(nBinsX=60, x0=0.0,  x1=6.0,   xUnits='',    xLabel='#DeltaR_{sig E, base nonCalo Mu}'),
    'dRy_sEl_bMu_Calo'   : HistOp1D(nBinsX=60, x0=0.0,  x1=6.0,   xUnits='',    xLabel='#DeltaR_{sig E, base Calo Mu}'),
    'isCaloTagged'    : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5,   xUnits='',    xLabel='Calo-Tagged Muon', logY=True),
    'dilep_flav'      : HistOp1D(nBinsX=5,  x0=-0.5, x1=4.5,   xUnits='',    xLabel='Dilepton flavor'),
    'isEM'            : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5,   xUnits='',    xLabel='Dilepton flavor is el mu'),
    'isME'            : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5,   xUnits='',    xLabel='Dilepton flavor is mu el'),
    'MtLep0'          : HistOp1D(nBinsX=15, x0=0.0,  x1=250.0, xUnits='GeV', xLabel='m_{T}(l_{0},MET)'),
    'MtLep1'          : HistOp1D(nBinsX=20, x0=0.0,  x1=140.0, xUnits='GeV', xLabel='m_{T}(l_{1},MET)'),
    #'MLL'             : HistOp1D(nBinsX=80, x0=50.0, x1=130, xUnits='GeV', xLabel='M_{ll}', logY=True),
    'MLL'             : HistOp1D(nBinsX=100, x0=0.0, x1=300.0, xUnits='GeV', xLabel='M_{ll}'),
    'ptll'            : HistOp1D(nBinsX=50, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='pT_{ll}', logY=True),
    # MET + leptons
    'MET'             : HistOp1D(nBinsX=50, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='E_{T}^{miss}', logY=True),
    'METPhi'          : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='MET_{#phi}'),
    'MCollASym'       : HistOp1D(nBinsX=25, x0=0.0,  x1=250.0, xUnits='GeV', xLabel='LFV Collinear Mass m_{coll}'),
    'dpt_ll'          : HistOp1D(nBinsX=20, x0=0.0,  x1=150.0, xUnits='GeV', xLabel='#Deltap_{T}^{ll}'),
    'DphiLep0MET'     : HistOp1D(nBinsX=63, x0=-3.15,  x1=3.15,  xUnits='',    xLabel='#Delta#phi(l_{0},MET)', add_underflow=True, logY=True),
    'DphiLep1MET'     : HistOp1D(nBinsX=63, x0=-3.15,  x1=3.15,  xUnits='',    xLabel='#Delta#phi(l_{1},MET)', add_underflow=True, logY=True),
    'tau_pT'          : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{subleading lep + MET}'),
    'taulep1_pT_ratio': HistOp1D(nBinsX=20, x0=0.0,  x1=3,     xUnits='',    xLabel='p_{T}^{subleading lep + MET} / p_{T}^{leading lep}'),
    # Jets
    'preJet_pt'       : HistOp1D(nBinsX=20, x0=0.0,  x1=100.0, xUnits='GeV', xLabel='p_{T}^{preJet}'),
    'preJet_eta'      : HistOp1D(nBinsX=100, x0=-5.0, x1=5.0,   xUnits='',    xLabel='#eta_{preJet}'),
    'preJet_JVT'      : HistOp1D(nBinsX=39, x0=-0.2,  x1=1.1,     xUnits='',    xLabel='preJet JVT (|eta|<=2.4 & pT < 60)', logY=True),
    'baseJet_eta'     : HistOp1D(nBinsX=100, x0=-5.0, x1=5.0,   xUnits='',    xLabel='#eta_{baseJet}'),
    'baseJet_mv2c10'  : HistOp1D(nBinsX=80, x0=-2,   x1=2,    xUnits='',    xLabel='mv2c10_{baseJet}', logY=True),
    'jetN'            : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{base jet}'),
    'jet_N2p4Eta25Pt' : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{jet} (p_{T}>25GeV, |#eta|<2.5)'),
    'signalJetN'      : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{sig jets}'),
    'jetN_g30'        : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{jet} (p_{T}>30GeV)'),
    'nLJets'          : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{CL jet}', logY=True),
    'nBJets'          : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{CB jet}', logY=True),
    'btag'            : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5,   xUnits='',    xLabel='B-tagged jet'),
    'nForwardJets'    : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{F jet}'),
    'j_pt[0]'         : HistOp1D(nBinsX=20, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='p_{T}^{leading jet}', logY=True),
    'j_pt[1]'         : HistOp1D(nBinsX=20, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='p_{T}^{subleading jet}', logY=True),
    'j_pt[2]'         : HistOp1D(nBinsX=20, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='p_{T}^{3rd leading jet}', logY=True),
    'j_pt[3]'         : HistOp1D(nBinsX=20, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='p_{T}^{4th leading jet}', logY=True),
    'j_pt'            : HistOp1D(nBinsX=20, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='Jet p_{T}'),
    'mjj'             : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='GeV', xLabel='Dijet mass'),
    'dEtaJJ'          : HistOp1D(nBinsX=20, x0=0.0,  x1=6.0,   xUnits='',    xLabel='#Delta#eta(j0,j1)'),
    'j_eta'           : HistOp1D(nBinsX=100, x0=-5.0, x1=5.0,   xUnits='',    xLabel='Jet #eta'),
    'j_jvt'           : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Jet JVT'),
    'j_jvf'           : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Jet JVF'),
    'j_phi'           : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='Jet #phi'),
    'j_flav'          : HistOp1D(nBinsX=5,  x0=-0.5, x1=4.5,   xUnits='',    xLabel='Jet flavor (0:NA,1:CL,2:CB,3:F)'),
    # Leptons
    'preEl_EcaloClus'           : HistOp1D(nBinsX=50, x0= 0,   x1=50, xUnits='GeV', xLabel='E_{preEl CaloCluster}', add_overflow=False),
    'preEl_etaCaloClus'         : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0, xUnits='',    xLabel='#eta_{preEl CaloCluster}'),
    'baseEl_etconetopo20'       : HistOp1D(nBinsX=50, x0= -7,   x1=15, xUnits='GeV', xLabel='E_{T}^{baselineEl} conetopo20' , logY=False),
    'baseEl_ptvarcone20'        : HistOp1D(nBinsX=25, x0= -1,   x1=4, xUnits='GeV', xLabel='p_{T}^{baselineEl} ptvarcone20' , logY=False),
    'baseMu_etconetopo20'       : HistOp1D(nBinsX=60, x0= -3,   x1=3, xUnits='GeV', xLabel='E_{T}^{baselineMu} conetopo20' ),
    'baseMu_ptvarcone30'        : HistOp1D(nBinsX=60, x0= -3,   x1=3, xUnits='GeV', xLabel='p_{T}^{baselineMu} ptvarcone30' , logY=False),
    'el1pT_trackclus_ratio'     : HistOp1D(nBinsX=30, x0=0,    x1=3,   xUnits='', xLabel='el_{subleading pT}^{track} / el_{subleading pT}^{cluster}', logY=True),

    # Fakes
    'nLepID' : HistOp1D(nBinsX=9, x0=-1.5,    x1=7.5,   xUnits='', xLabel='ID lepton multiplicity', logY=True),
    'nLepAntiID' : HistOp1D(nBinsX=9, x0=-1.5,    x1=7.5,   xUnits='', xLabel='anti-ID lepton multiplicity', logY=True),
    'aID_Lep0Pt'          : HistOp1D(nBinsX=40, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{leading lep} anti-ID', logY=True),
    'aID_Lep1Pt'          : HistOp1D(nBinsX=20, x0=0.0,  x1=100.0, xUnits='GeV', xLabel='p_{T}^{subleading lep} anti-ID', logY=True),
    'aID_Lep0Eta'         : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta^{leading lep} anti-ID', logY=True),
    'aID_Lep1Eta'         : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta^{subleading lep} anti-ID', logY=True),
    'aID_MLL'             : HistOp1D(nBinsX=100, x0=0.0, x1=300.0, xUnits='GeV', xLabel='M_{ll}(ID, anti-ID)'),
    'aID_Lep0Flav'        : HistOp1D(nBinsX=7,  x0=-1.5, x1=5.5,   xUnits='',    xLabel='Leading antiID flavor', logY=True),
    'aID_Lep1Flav'        : HistOp1D(nBinsX=7,  x0=-1.5, x1=5.5,   xUnits='',    xLabel='Subleading antiID flavor', logY=True),
    'aID_Lep0Q'           : HistOp1D(nBinsX=3,  x0=-1.5, x1=1.5,   xUnits='',    xLabel='Leading anti-ID charge'),
    'aID_Lep1Q'           : HistOp1D(nBinsX=3,  x0=-1.5, x1=1.5,   xUnits='',    xLabel='Subleading anti-ID charge'),
    'Z_MLL'               : HistOp1D(nBinsX=50, x0=65, x1=115.0, xUnits='GeV', xLabel='M_{ll} (best lep pair)'),
    'Z2_MLL'               : HistOp1D(nBinsX=50, x0=0, x1=200.0, xUnits='GeV', xLabel='M_{ll} (2nd best lep pair)'),
    'aID_dpt_ll'          : HistOp1D(nBinsX=20, x0=0.0,  x1=150.0, xUnits='GeV', xLabel='#Deltap_{T}^{ll}(ID, anti-ID)'),
    'aID_drll'            : HistOp1D(nBinsX=60, x0=0.0,  x1=6.0,   xUnits='',    xLabel='#DeltaR_{ll}(ID, anti-ID)'),
    'dR_Zl'               : HistOp1D(nBinsX=60, x0=0.0,  x1=6.0,   xUnits='',    xLabel='#DeltaR(Z, lep)'),
    'Z_dilep_flav'        : HistOp1D(nBinsX=7,  x0=-1.5, x1=5.5,   xUnits='',    xLabel='Z dilepton flavor', logY=True),
    'Z2_dilep_flav'       : HistOp1D(nBinsX=7,  x0=-1.5, x1=5.5,   xUnits='',    xLabel='2nd Z dilepton flavor', logY=True),
    'Z_Lep2_pT'         : HistOp1D(nBinsX=40, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='3rd leading lepton p_{T}', logY=True),
    #'Z_Lep2_pT'         : HistOp1D(nBinsX=50, x0=0.0,  x1=50.0, xUnits='GeV', xLabel='3rd leading lepton p_{T}', logY=False),
    'Z_Lep2_eta'        : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='3rd leading lepton #eta', logY=True),
    'Z_dilep_sign'       :HistOp1D(nBinsX=5, x0=-2.5, x1=2.5,   xUnits='',    xLabel='Z Dilepton Sign : OS(-1) SS(1)', logY=True),
    'Z2_dilep_sign'       :HistOp1D(nBinsX=5, x0=-2.5, x1=2.5,   xUnits='',    xLabel='2nd Z Dilepton Sign : OS(-1) SS(1)', logY=True),
    'Z_Lep2_dPhi_MET'     : HistOp1D(nBinsX=63, x0=-3.15,  x1=3.15,  xUnits='',    xLabel='#Delta#phi(l_{3},MET)', add_underflow=True, logY=True),
    'Z_Lep2_mT'           : HistOp1D(nBinsX=50, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='3rd leading lepton m_{T}', logY=False),
}

################################################################################
# Toggle options for execution
# - Options expected to change often between different plottings
################################################################################

#######################################
## Build the TChain/TTree for each sample
# To remove sample from plot, comment out its formation here and where the
# background gets appended to the samples list
data_dsids = g.groups['data15']+g.groups['data16']
data.set_chain_from_dsid_list(data_dsids, data_rawdir, exclude_strs='FFest')
#fakes.set_chain_from_dsid_list(data_dsids, data_rawdir, search_strs='FFest')
ttbar.set_chain_from_dsid_list(g.groups['ttbar'], rawdir)
stop.set_chain_from_dsid_list(g.groups['singletop'], rawdir)
wtop.set_chain_from_dsid_list(g.groups['Wt'], rawdir)
WW.set_chain_from_dsid_list(g.groups['ww'], rawdir)
ZZ.set_chain_from_dsid_list(g.groups['zz'], rawdir)
WZ.set_chain_from_dsid_list(g.groups['wz'], rawdir)
#zll.set_chain_from_dsid_list(g.groups['zll'], rawdir)
zee.set_chain_from_dsid_list(g.groups['zee'], rawdir)
zmumu.set_chain_from_dsid_list(g.groups['zmumu'], rawdir)
ztt.set_chain_from_dsid_list(g.groups['ztt'], rawdir)
wjets.set_chain_from_dsid_list(g.groups['wjets'], rawdir)
#wgamma.set_chain_from_dsid_list(g.groups['w_gamma'], rawdir)
htt.set_chain_from_dsid_list(g.groups['htt'], rawdir)
hww.set_chain_from_dsid_list(g.groups['hww'], rawdir)
#signal.set_chain_from_dsid_list(g.groups['higgs_lfv'], signal_rawdir)

## Add backgrounds to be included in plots
backgrounds = []
backgrounds.append(ttbar)
backgrounds.append(stop)
backgrounds.append(wtop)
backgrounds.append(WW)
backgrounds.append(ZZ)
backgrounds.append(WZ)
#backgrounds.append(zll)
backgrounds.append(zee)
backgrounds.append(zmumu)
backgrounds.append(ztt)
backgrounds.append(wjets)
#backgrounds.append(wgamma)
backgrounds.append(htt)
backgrounds.append(hww)
#backgrounds.append(fakes)
#backgrounds.append(signal)

# What regions to plot
region_ops = []
region_ops += ['zjets_FF_CRden_e', 'zjets_FF_CRden_m']

# What variables to plot
vars_to_plot = []
vars_to_plot += ['Lep0Pt']

# Remove duplicate names
vars_to_plot = list(set(vars_to_plot))

################################################################################
# Create plots
################################################################################
plots = []
for var in vars_to_plot:

    # Use vector settings when plotting a specific element of the vector
    # i.e. lepton_flav[1] uses lepton_flav settings if lepton_flav[1]
    # is not defined
    if '[' in var and var not in HistOpMap:
        var_trim = var.split('[')[0]
        ops = HistOpMap[var_trim]
    else:
        ops = HistOpMap[var]


    # Set contingent defaults
    y0_def = 0.1 if ops.logY and not ops.norm else 1e-4 if ops.norm else 0
    y1_def = 1e7 if ops.logY and not ops.norm else 1 if ops.norm else 5000

    # Determine plot properties
    name_ = re.sub(r'[(){}[\]]+','',var)
    bin_width = (ops.x1 - ops.x0) / float(ops.nBinsX)
    width_label = str(round(bin_width, 2))
    if not ops.xUnits and width_label == '1.0':
        yLabel = ops.yLabel
    elif width_label == '1.0':
        yLabel = "%s / %s"%(ops.yLabel, ops.xUnits)
    else:
        yLabel = "%s / %s %s"%(ops.yLabel,width_label,ops.xUnits)
    y0 = ops.y0 if ops.y0 else y0_def
    y1 = ops.y1 if ops.y1 else y1_def

    # Create plot for each indicated region
    #for region in ops.regions.split(','):
    for region in region_ops:
        region = region.strip()
        p = plot.Plot1D()
        p.initialize(region, var, name="%s_%s"%(region, name_))
        p.doLogY = ops.logY
        p.doNorm = ops.norm
        p.add_overflow = ops.add_overflow
        p.labels(x=ops.xLabel, y=yLabel)
        p.xax(bin_width, ops.x0, ops.x1)
        p.yax(y0, y1)
        assert data or backgrounds, "No data or backgrounds defined"
        #if data and backgrounds and ops.ratioPlot:
        #    p.setRatioCanvas(p.name)
        #else:
        p.setDefaultCanvas(p.name)

        plots.append(p)

print "Configuration File End"
