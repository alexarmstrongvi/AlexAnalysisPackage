import sys, os
import re
import ROOT
from collections import namedtuple
import tools.plot as plot
import tools.background as background
import tools.region as region
import tools.systematic as systematic
import global_variables as g

#######################################
# samples
########################################

tag = g.tag
analysis_dir = g.analysis_dir
analysis_run_dir = g.analysis_run_dir

# Path to directories with flat ntuples
rawdir        = g.analysis_run_dir + "ntuples/test_bjetOR/"
signal_rawdir = g.mc_ntuples
data_rawdir   = g.data_ntuples

# Path to directory with condor filelists
filelist_dir = g.input_files

# Luminosity options [2015-16 (ipb), 2015-16 (ifb), 2015 (ipb), 2016 (ipb)]
lumi_val = 0
lumi_ = [36180, 36.01, 3209, 32971, 1][lumi_val] 

backgrounds = []
######## MC
## ttbar
ttbar = background.Background("ttbar", "t#bar{t}")
ttbar.set_debug()
ttbar.scale_factor = lumi_
ttbar.set_fillStyle(0)
ttbar.setLineStyle(1)
ttbar.set_color(ROOT.kOrange+2)
ttbar.set_treename("ttbar")
#ttbar.set_chain_from_file_list(['miniNtuple_410009_STT.root'], rawdir)
#ttbar.set_chain_from_file_list(['miniNtuple_410009_SFF.root'], rawdir)
#ttbar.set_chain_from_file_list(['miniNtuple_410009_SFT.root'], rawdir)
#ttbar.set_chain_from_file_list(['miniNtuple_410009_STF.root'], rawdir)
#ttbar.set_chain_from_file_list(['miniNtuple_410009_BTT.root'], rawdir)
#ttbar.set_chain_from_file_list(['miniNtuple_410009_Base.root'], rawdir)
#ttbar.set_chain_from_file_list(['miniNtuple_410009_1.root'], rawdir)
#ttbar.set_chain_from_file_list(['miniNtuple_410009_2.root'], rawdir)
#ttbar.set_chain_from_file_list(['miniNtuple_410009_3.root'], rawdir)
#ttbar.set_chain_from_file_list(['miniNtuple_410009_4.root'], rawdir)
#ttbar.set_chain_from_file_list(['miniNtuple_410009_5.root'], rawdir)
ttbar.set_chain_from_file_list(['miniNtuple_410009_test.root'], rawdir)
#ttbar.CheckForDuplicates()
backgrounds.append(ttbar)

#### DATA
#data = background.Data()
#data.set_color(ROOT.kBlack)
#data.set_treename("data")
##data.set_chain_from_list_CONDOR(filelist_dir + "data/", data_rawdir)
#data.set_chain_from_list_CONDOR2([filelist_dir + "data15/", filelist_dir + "data16/"] , [data_rawdir,data_rawdir])
##data.CheckForDuplicates()


##########################################
# systematics
##########################################
systematics = []

##########################################
# regions
##########################################
regions = []

### Cuts ###
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



# Region building blocks
Baseline_Sel = ('Lep0Pt >= 45 && Lep1Pt >= 15 '
              + '&& (30 < MLL && MLL < 150) '
              + '&& nBJets==0 '
              + '&& ( !'+mue+' || el1pT_trackclus_ratio < 1.2) '
              + '&&' + DF_OS
              + '&&' + singlelep_trig_pT) 
VBF_stripped = "JetN_g30 >= 2 && j_pt[0] > 40 && Mjj > 400 && DEtaJJ > 3"

# Define regions

reg = region.Region()
reg.name = "no_sel"
reg.displayname = "No Selections"
reg.tcut = DF_OS
regions.append(reg)

reg = region.Region()
reg.name = "trig_only"
reg.displayname = "Single Lepton Triggers"
reg.tcut = singlelep_trig_pT + "&&" + DF_OS
regions.append(reg)

#reg = region.Region()
#reg.name = "topCR"
#reg.displayname = "Top CR"
#reg.tcut = "nBJets >= 1 && MET > 40 &&" + DF_OS + " &&" + singlelep_trig
#regions.append(reg)
#
#reg = region.Region()
#reg.name = "zCR_ee"
#reg.displayname = "Z CR (Channel: El-El)"
#reg.tcut = "75 < MLL && MLL < 105 && " + SF_OS + " && " + ee + " && " + singlelep_trig
#regions.append(reg)
#
#reg = region.Region()
#reg.name = "zCR_mumu"
#reg.displayname = "Z CR (Channel: Mu-Mu)"
#reg.tcut = "75 < MLL && MLL < 105 && " + SF_OS + " && " + mumu + " && " + singlelep_trig
#regions.append(reg)

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


##########################################
# plots
##########################################
plots = []

vars = {}

### Define variables for plots
HistOp1D = namedtuple('HistOp1D', 'regions,   nBinsX, x0, x1, y0,   y1,   logY,   xUnits, xLabel, yLabel, add_overflow')
HistOp1D.__new__.__defaults__=   ('baseline', 12,     -2,  10, None, None, False,  '',     '',     'Events', True)

HistOpMap = {
    'MLL' : HistOp1D(nBinsX=18, x0=20, x1=160, xUnits='GeV', xLabel='M_{ll}', regions='baseline, baseline_emu, baseline_mue', add_overflow=False),
    'treatAsYear' : HistOp1D(nBinsX=10, x0=2010.5, x1=2020.5, xUnits='', xLabel='Year', regions='baseline, baseline_emu, baseline_mue'),
}

vars_to_plot = []
## Variable to get Yields
#vars_to_plot += ['treatAsYear']

## Quick Plots
#vars_to_plot += ['nCentralBJets', 'nCentralLJets', 'j_pt[0]', 'j_pt[1]', 'j_pt[2]', 'j_pt[3]']
vars_to_plot += ['MLL']

## Z CR
#vars_to_plot += ['Lep0Pt', 'Lep1Pt', 'MLL']

## Top CR
#vars_to_plot += ['nCentralBJets', 'Lep0Pt', 'Lep1Pt', 'j_pt[0]', 'j_pt[1]', 'MET']

## Trigger Variables
#vars_to_plot += ["Lep0Pt", "Lep0Phi", "Lep0Eta", "Lep1Pt", "Lep1Phi", "Lep1Eta"]
#vars_to_plot += ["El_ID", "Mu_ID", "El_etconetopo20", "El_ptvarcone20", "Mu_ptvarcone30", "Mu_etconetopo20"]

## Baseline Selection
#vars_to_plot += ['Lep0Pt', 'Lep1Pt', 'MLL', 'nCentralBJets', 'el1pT_trackclus_ratio']

## Multiplicity Plots
#vars_to_plot += ['nCentralLJets', 'nCentralBJets', 'nForwardJets', 'n_leptons']

## Basic Kinematics
#vars_to_plot += ['Lep0Pt', 'Lep1Pt', 'MET', 'j_pt', 'Lep0Eta', 'Lep1Eta', 'j_eta',] 

## Object Definitions Plots
#vars_to_plot += ['n_preLeptons',
#                'preEl_EcaloClus', 'preEl_etaCaloClus', 'preEl_Et', 'preEl_pt',
#                'preEl_clusEtaBE', 'preEl_eta', 'baseEl_etconetopo20',
#                'baseEl_ptvarcone20', 'baseEl_ID', 'el_type', 'el_origin',
#                'el1_track_pt', 'el1_clus_pt', 
#                'preMu_pt', 'preMu_ID', 'baseMu_pt', 'baseMu_eta',
#                'baseMu_etconetopo20', 'baseMu_ptvarcone30', 'baseMu_ID',
#                'mu_type', 'mu_origin', 'preTau_q', 'preTau_nTracks',
#                'baseTau_pT', 'baseTau_eta', 'baseTau_nTracks', 'baseTau_ID',
#                'tau_pT', 'taulep1_pT_ratio', 'preJet_pt', 'preJet_eta',
#                'preJet_JVT', 'baseJet_eta', 'baseJet_mv2c10']

# Remove duplicate names
vars_to_plot = list(set(vars_to_plot))

for var in vars_to_plot:
    ops = HistOpMap[var]

    # Set contingent defaults
    y0_def = 0.1 if ops.logY else 0
    y1_def = 1e7 if ops.logY else 5000

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
    for region in ops.regions.split(','):
        region = region.strip()
        print "Plotting for region", region
        p = plot.Plot1D()
        p.initialize(region, var, name="%s_%s"%(region, name_))
        p.doLogY = ops.logY
        p.add_overflow = ops.add_overflow
        p.labels(x=ops.xLabel, y=yLabel)
        p.xax(bin_width, ops.x0, ops.x1)
        p.yax(y0, y1)
        #p.setRatioCanvas(p.name)
        p.setDefaultCanvas(p.name)
        plots.append(p)


