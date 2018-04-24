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
rawdir        = g.mc_ntuples
signal_rawdir = g.mc_ntuples
data_rawdir   = g.data_ntuples

# Path to directory with condor filelists
filelist_dir = g.input_files

# Luminosity options [2015-16 (ipb), 2015-16 (ifb), 2015 (ipb), 2016 (ipb)]
lumi_val = 0
lumi_ = [36180, 36.01, 3209, 32971, 1][lumi_val] 

backgrounds = []
####### MC
# ttbar
ttbar = background.Background("ttbar", "t#bar{t}")
ttbar.set_debug()
ttbar.scale_factor = lumi_
ttbar.set_fillStyle(0)
ttbar.setLineStyle(1)
ttbar.set_color(ROOT.kOrange+2)
ttbar.set_treename("ttbar")
ttbar.set_chain_from_dsid_list(g.groups['ttbar'], rawdir)

# singletop
stop = background.Background("st", "Single top")
stop.set_debug()
stop.scale_factor = lumi_
stop.set_fillStyle(0)
stop.setLineStyle(1)
stop.set_color(ROOT.kOrange+1)
stop.set_treename("ST")
stop.set_chain_from_dsid_list(g.groups['singletop'], rawdir)

# W+top
wtop = background.Background("wt", "Wt")
wtop.set_debug()
wtop.scale_factor = lumi_
wtop.set_fillStyle(0)
wtop.setLineStyle(1)
wtop.set_color(ROOT.kOrange+8)
wtop.set_treename("wt")
wtop.set_chain_from_dsid_list(g.groups['Wt'], rawdir)

# WW
WW = background.Background("ww", "WW")
WW.set_debug()
WW.scale_factor = lumi_ 
WW.set_fillStyle(0)
WW.setLineStyle(1)
WW.set_color(ROOT.kSpring-6)
WW.set_treename("WW")
WW.set_chain_from_dsid_list(g.groups['ww'], rawdir)

# ZZ
ZZ = background.Background("zz", "ZZ")
ZZ.set_debug()
ZZ.scale_factor = lumi_ 
ZZ.set_fillStyle(0)
ZZ.setLineStyle(1)
ZZ.set_color(ROOT.kSpring-4)
ZZ.set_treename("ZZ")
ZZ.set_chain_from_dsid_list(g.groups['zz'], rawdir)

# WZ
WZ = background.Background("wz", "WZ")
WZ.set_debug()
WZ.scale_factor = lumi_ 
WZ.set_fillStyle(0)
WZ.setLineStyle(1)
WZ.set_color(ROOT.kSpring-5)
WZ.set_treename("WZ")
WZ.set_chain_from_dsid_list(g.groups['wz'], rawdir)

# Zll
#zll = background.Background("zll", "Zll")
#zll.set_debug()
#zll.scale_factor = lumi_
#zll.set_fillStyle(0)
#zll.setLineStyle(1)
#zll.set_color(ROOT.kAzure-9)
#zll.set_treename("zll")
#zll.set_chain_from_dsid_list(g.groups['zll'], rawdir)

# Zee
zee = background.Background("zee", "Zee")
zee.set_debug()
zee.scale_factor = lumi_
zee.set_fillStyle(0)
zee.setLineStyle(1)
zee.set_color(ROOT.kAzure-7)
zee.set_treename("zee")
zee.set_chain_from_dsid_list(g.groups['zee'], rawdir)

# Zmumu
zmumu = background.Background("zmumu", "Zmumu")
zmumu.set_debug()
zmumu.scale_factor = lumi_
zmumu.set_fillStyle(0)
zmumu.setLineStyle(1)
zmumu.set_color(ROOT.kAzure-9)
zmumu.set_treename("zmumu")
zmumu.set_chain_from_dsid_list(g.groups['zmumu'], rawdir)

# Ztt
ztt = background.Background("ztt", "Z#tau#tau")
ztt.set_debug()
ztt.scale_factor = lumi_
ztt.set_fillStyle(0)
ztt.setLineStyle(1)
ztt.set_color(ROOT.kAzure-5)
ztt.set_treename("ztt")
ztt.set_chain_from_dsid_list(g.groups['ztt'], rawdir)

# Wjets
wjets = background.Background("wjets", "W+jets")
wjets.set_debug()
wjets.scale_factor = lumi_
wjets.set_fillStyle(0)
wjets.setLineStyle(1)
wjets.set_color(ROOT.kOrange)
wjets.set_treename("wjets")
wjets.set_chain_from_dsid_list(g.groups['wjets'], rawdir)

# W+gamma
wgamma = background.Background("wgamma", "W+gamma")
wgamma.set_debug()
wgamma.scale_factor = lumi_
wgamma.set_fillStyle(0)
wgamma.setLineStyle(1)
wgamma.set_color(ROOT.kOrange-1)
wgamma.set_treename("w_gamma")
wgamma.set_chain_from_dsid_list(g.groups['w_gamma'], rawdir)

## Higgs -> tau tau
htt = background.Background("htt", "H#tau#tau")
#higgs = background.Background("higgs", "Higgs")
htt.scale_factor = lumi_
htt.set_fillStyle(0)
htt.setLineStyle(1)
htt.set_color(ROOT.kRed)
htt.set_treename("htt")
htt.set_chain_from_dsid_list(g.groups['htt'], rawdir)

## Higgs -> W W
hww = background.Background("hww", "HWW")
#higgs = background.Background("higgs", "Higgs")
hww.scale_factor = lumi_
hww.set_fillStyle(0)
hww.setLineStyle(1)
hww.set_color(ROOT.kBlue+3)
hww.set_treename("hww")
hww.set_chain_from_dsid_list(g.groups['hww'], rawdir)

signal_branching_ratio = 0.01
signal_SF = 1
signal_label = "Higgs LFV" if signal_SF == 1 else "Higgs LFV (%dX)"%signal_SF
signal = background.Background("higgs_lfv", signal_label)
signal.setSignal()
signal.scale_factor = lumi_ * signal_branching_ratio * signal_SF 
signal.set_fillStyle(0)
signal.set_color(ROOT.kGreen)
signal.set_treename("signal")
signal.set_chain_from_dsid_list(g.groups['higgs_lfv'], signal_rawdir)


#### DATA
data = background.Data()
data.set_color(ROOT.kBlack)
data.set_treename("data")
#data.set_chain_from_list_CONDOR(filelist_dir + "data/", data_rawdir)
data.set_chain_from_list_CONDOR2([filelist_dir + "data15/", filelist_dir + "data16/"] , [data_rawdir,data_rawdir])


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

# Define regions

reg = region.Region()
reg.name = "no_sel"
reg.displayname = "No Selections"
reg.tcut = DF_OS
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
zjets_FF_CR_add += ' && Z_Lep2_mT < 50'
zjets_FF_CR_add += ' && (Z2_MLL < 80 || 100 < Z2_MLL)'
FF_CR_eee_ch = 'Z_dilep_flav==2 && Z_Lep2_flav==1'
FF_CR_mme_ch = 'Z_dilep_flav==2 && Z_Lep2_flav==1'
FF_CR_eem_ch = 'Z_dilep_flav==3 && Z_Lep2_flav==0'
FF_CR_mmm_ch = 'Z_dilep_flav==3 && Z_Lep2_flav==0'
num_den_dict = {'den' : 'nLepID == 2 && nLepAntiID >= 1',
                'num' : 'nLepID == 3'}
chan_dict = {'eee' : ['ee','e','Z_dilep_flav==2 && Z_Lep2_flav==1'],
             'mme' : ['mumu','e','Z_dilep_flav==3 && Z_Lep2_flav==1'],
             'eem' : ['ee','m','Z_dilep_flav==2 && Z_Lep2_flav==0'],
             'mmm' : ['mumu','m','Z_dilep_flav==3 && Z_Lep2_flav==0']
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
wjets_FF_CRden_base = 'nLepID == 1 && nLepAntiID >= 1' # require final state
wjets_FF_CRnum_base = 'nLepID == 2' # require final state
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


##########################################
# plots
##########################################
plots = []

vars = {}

### Define variables for plots
HistOp1D = namedtuple('HistOp1D', 'regions,   nBinsX, x0, x1, y0,   y1,   logY,   xUnits, xLabel, yLabel, add_overflow, add_underflow')
HistOp1D.__new__.__defaults__=   ('baseline', 12,     -2,  10, None, None, False,  '',     '',     'Events', True, False)

HistOpMap = {
    # Event level
    'RunNumber'       : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Event run number',                  regions='baseline'),
    'event_number'    : HistOp1D(nBinsX=1000, x0=0,  x1=1000000000,    xUnits='',    xLabel='Event number',              regions='baseline_emu, baseline_mue'),
    'isMC'            : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5,   xUnits='',    xLabel='is Monte Carlo',                    regions='baseline'),
    'eventweight'     : HistOp1D(nBinsX=100, x0=-0.001,  x1=0.002,    xUnits='',    xLabel='Event weight',                      regions='no_sel, baseline', logY=True),
    'dsid'            : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Sample DSID',                       regions='baseline'),
    'treatAsYear'     : HistOp1D(nBinsX=11, x0=2007.5,x1=2018.5,xUnits='',   xLabel='treatAsYear',                       regions='trig_only, topCR, ztautauCR, baseline, baseline_mue, baseline_emu'),
    # Multiplicity
    'n_preLeptons'    : HistOp1D(nBinsX=11, x0=-0.5, x1=10.5,  xUnits='',    xLabel='N_{pre-leptons}',                   regions='no_sel'),
    'n_baseLeptons'   : HistOp1D(nBinsX=11, x0=-0.5, x1=10.5,  xUnits='',    xLabel='N_{baseline leptons}',              regions='baseline'),
    'n_leptons'       : HistOp1D(nBinsX=11, x0=-0.5, x1=10.5,  xUnits='',    xLabel='N_{signal leptons}',                regions='baseline, baseline_emu, baseline_mue'),
    # Leptons
    'l_pt'            : HistOp1D(nBinsX=25, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='Lepton p_{T}',                      regions='baseline, baseline_emu, baseline_mue'),
    'aID_Lep0Pt'          : HistOp1D(nBinsX=40, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{leading lep} anti-ID', logY=True),
    ## Electrons
    'el1_track_pt'    : HistOp1D(nBinsX=25, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='Leading electron track p_{T}',      regions='baseline'),
    'el1_clus_pt'     : HistOp1D(nBinsX=25, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='Leading electron cluster p_{T}',    regions='baseline'),
    'preEl_Et'        : HistOp1D(nBinsX=50, x0=0.0,  x1=50.0, xUnits='GeV', xLabel='E_{T}^{preElectron}',               regions='no_sel', add_overflow=False),
    'preEl_pt'        : HistOp1D(nBinsX=50, x0=0.0,  x1=50.0, xUnits='GeV', xLabel='p_{T}^{preElectron}',               regions='no_sel'),
    'preEl_clusEtaBE' : HistOp1D(nBinsX=60, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta_{preEl clusterBE}',            regions='no_sel'),
    'preEl_eta'       : HistOp1D(nBinsX=60, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta_{preElectron}',                regions='no_sel'),
    'baseEl_ID'       : HistOp1D(nBinsX=9,  x0=-0.5, x1=8.5,   xUnits='',    xLabel='el_{Baseline}^{ID}(non-inclusive)', regions='no_sel', logY=True),
    'El_ID'           : HistOp1D(nBinsX=7,  x0=-1.5, x1=5.5,   xUnits='',    xLabel='electron ID(non-inclusive)', regions='baseline_emu, baseline_mue', logY=False),
    'el_type'         : HistOp1D(nBinsX=41,  x0=-1.5, x1=39.5,   xUnits='',    xLabel='electron truth type',               regions='baseline, baseline_emu, baseline_mue', logY=True),
    'el_origin'       : HistOp1D(nBinsX=48,  x0=-1.5, x1=46.5,   xUnits='',    xLabel='electron truth origin',             regions='baseline, baseline_emu, baseline_mue', logY=True),
    'el_d0sigBSCorr'  : HistOp1D(nBinsX=60,  x0=-6, x1=6,   xUnits='mm',    xLabel='Electron d_{0}/#sigma_{d_{0}} BSCorr',     regions='baseline', add_underflow=True),
    'el_z0SinTheta'  : HistOp1D(nBinsX=60,  x0=-0.6, x1=0.6,   xUnits='mm',    xLabel='Electron z_{0}sin(#theta)',     regions='baseline', add_underflow=True),
    ## Muons
    'preMu_pt'        : HistOp1D(nBinsX=50, x0=0,  x1=50.0, xUnits='GeV', xLabel='p_{T}^{preMuon}',                   regions='no_sel', add_overflow=False),
    'preMu_ID'      : HistOp1D(nBinsX=7,  x0=-1.5, x1=5.5,   xUnits='',    xLabel='mu_{pre}^{ID} (non-inclusive)',     regions='no_sel'),
    'baseMu_pt'       : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{BaselineMuon}',              regions='no_sel'),
    'baseMu_eta'      : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta_{BaselineMuon}',               regions='no_sel'),
    'baseMu_ID'     : HistOp1D(nBinsX=7,  x0=-1.5, x1=5.5,   xUnits='',    xLabel='mu_{Baseline}^{ID}(non-inclusive)', regions='no_sel', logY=True),
    'Mu_ID'           : HistOp1D(nBinsX=7,  x0=-1.5, x1=5.5,   xUnits='',    xLabel='muon ID (non-inclusive)', regions='baseline_mue, baseline_emu', logY=False),
    'mu_type'         : HistOp1D(nBinsX=41,  x0=-1.5, x1=39.5,   xUnits='',    xLabel='muon truth type',                   regions='baseline, baseline_emu, baseline_mue', logY=True),
    'mu_origin'       : HistOp1D(nBinsX=48,  x0=-1.5, x1=46.5,   xUnits='',    xLabel='muon truth origin',                 regions='baseline, baseline_emu, baseline_mue', logY=True),
    'mu_d0sigBSCorr'  : HistOp1D(nBinsX=60,  x0=-6, x1=6,   xUnits='mm',    xLabel='Muon d_{0}/#sigma_{d_{0}} BSCorr',     regions='baseline', add_underflow=True),
    'mu_z0SinTheta'  : HistOp1D(nBinsX=60,  x0=-0.6, x1=0.6,   xUnits='mm',    xLabel='Muon z_{0}sin(#theta)',     regions='baseline', add_underflow=True),
    ## Taus
    'preTau_q'        : HistOp1D(nBinsX=33,  x0=-5.5, x1=5.5,   xUnits='',    xLabel='Tau charge',                        regions='no_sel'),
    'preTau_nTracks'  : HistOp1D(nBinsX=10,  x0=-1.5, x1=8.5,   xUnits='',    xLabel='preTau nTracks',                    regions='no_sel'),
    'baseTau_pT'      : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{BaselineTau}',               regions='no_sel'),
    'baseTau_eta'     : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta_{BaselineMuon}',               regions='no_sel'),
    'baseTau_nTracks' : HistOp1D(nBinsX=10,  x0=-1.5, x1=8.5,   xUnits='',    xLabel='Baseline Tau nTracks',              regions='no_sel'),
    'baseTau_ID'      : HistOp1D(nBinsX=7,  x0=-1.5, x1=5.5,   xUnits='',    xLabel='tau_{Baseline}^{ID}(non-inclusive)',regions='no_sel', logY=True),

    ## General lepton
    'l_eta'           : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='Lepton #eta',                       regions='baseline'),
    'l_phi'           : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='Lepton #phi',                       regions='baseline'),
    'l_flav'          : HistOp1D(nBinsX=5,  x0=-0.5, x1=4.5,   xUnits='',    xLabel='Lepton flavor (0: e, 1: m)',        regions='baseline'),
    'l_type'          : HistOp1D(nBinsX=41,  x0=-1.5, x1=39.5,   xUnits='',    xLabel='Lepton type',                       regions='baseline, baseline_emu, baseline_mue', logY=True),
    'l_origin'        : HistOp1D(nBinsX=48,  x0=-1.5, x1=46.5,    xUnits='',    xLabel='Lepton origin',                     regions='baseline, baseline_emu, baseline_mue', logY=True),
    'l_type[0]'       : HistOp1D(nBinsX=41,  x0=-1.5, x1=39.5,   xUnits='',    xLabel='Lepton type',                       regions='baseline, baseline_emu, baseline_mue', logY=True),
    'l_origin[0]'     : HistOp1D(nBinsX=48,  x0=-1.5, x1=46.5,    xUnits='',    xLabel='Lepton origin',                     regions='baseline, baseline_emu, baseline_mue', logY=True),
    'l_type[1]'       : HistOp1D(nBinsX=41,  x0=-1.5, x1=39.5,   xUnits='',    xLabel='Lepton type',                       regions='baseline, baseline_emu, baseline_mue', logY=True),
    'l_origin[1]'     : HistOp1D(nBinsX=48,  x0=-1.5, x1=46.5,    xUnits='',    xLabel='Lepton origin',                     regions='baseline, baseline_emu, baseline_mue', logY=True),
    'Lep_Iso'         : HistOp1D(nBinsX=9,  x0=-1.5, x1=7.5,   xUnits='',    xLabel='Lepton Isolation (non-inclusive)', regions='baseline'),
    'l_q'             : HistOp1D(nBinsX=3,  x0=-1.5, x1=1.5,   xUnits='',    xLabel='Lepton charge',                     regions='baseline'),
    'LepLepSign'      : HistOp1D(nBinsX=3,  x0=-1.5, x1=1.5,   xUnits='',    xLabel='Leptons sign product',              regions='baseline'),
    #'Lep0Pt'          : HistOp1D(nBinsX=40, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{leading lep}',               regions='no_sel, trig_only, ztautauCR, zCR_mumu, zCR_ee, topCR, baseline, baseline_emu, baseline_mue', logY=True),
    #'Lep1Pt'          : HistOp1D(nBinsX=40, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{subleading lep}',            regions='no_sel, trig_only, ztautauCR, zCR_mumu, zCR_ee, topCR, baseline, baseline_emu, baseline_mue', logY=True),
    'Lep0Pt'          : HistOp1D(nBinsX=50, x0=0.0,  x1=50.0, xUnits='GeV', xLabel='p_{T}^{leading lep}', add_overflow=False),
    'Lep1Pt'          : HistOp1D(nBinsX=50, x0=0.0,  x1=50.0, xUnits='GeV', xLabel='p_{T}^{subleading lep}', add_overflow=False),
    'Lep0Eta'         : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta^{leading lep}',                regions='no_sel, trig_only, baseline'),
    'Lep1Eta'         : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta^{subleading lep}',             regions='no_sel, trig_only, baseline'),
    'Lep0Phi'         : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='#phi^{leading lep}',                regions='no_sel, trig_only'),
    'Lep1Phi'         : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='#phi^{subleading lep}',             regions='no_sel, trig_only'),
    'MLep0'           : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='GeV', xLabel='M_{l0}',                            regions='baseline'),
    'MLep1'           : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='GeV', xLabel='M_{l1}',                            regions='baseline'),
    'DEtaLL'          : HistOp1D(nBinsX=20, x0=0.0,  x1=6.0,   xUnits='',    xLabel='#Delta#eta_{ll}',                   regions='baseline'),
    'DphiLL'          : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='#Delta#phi_{ll}',                   regions='baseline'),
    'drll'            : HistOp1D(nBinsX=60, x0=0.0,  x1=6.0,   xUnits='',    xLabel='#DeltaR_{ll}',                      regions='baseline_emu, baseline_mue'),
    'dRy_sEl_bMu_noCalo' : HistOp1D(nBinsX=60, x0=0.0,  x1=6.0,   xUnits='',    xLabel='#DeltaR_{sig E, base nonCalo Mu}', regions='baseline_emu, baseline_mue'),
    'dRy_sEl_bMu_Calo'   : HistOp1D(nBinsX=60, x0=0.0,  x1=6.0,   xUnits='',    xLabel='#DeltaR_{sig E, base Calo Mu}', regions='baseline_emu, baseline_mue'),
    'isCaloTagged'    : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5,   xUnits='',    xLabel='Calo-Tagged Muon',                  regions='no_sel', logY=True),
    'dilep_flav'      : HistOp1D(nBinsX=5,  x0=-0.5, x1=4.5,   xUnits='',    xLabel='Dilepton flavor',                   regions='baseline'),
    'isEM'            : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5,   xUnits='',    xLabel='Dilepton flavor is el mu',          regions='baseline'),
    'isME'            : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5,   xUnits='',    xLabel='Dilepton flavor is mu el',          regions='baseline'),
    'MtLep0'          : HistOp1D(nBinsX=15, x0=0.0,  x1=250.0, xUnits='GeV', xLabel='m_{T}(l_{0},MET)',                  regions='baseline'),
    'MtLep1'          : HistOp1D(nBinsX=20, x0=0.0,  x1=140.0, xUnits='GeV', xLabel='m_{T}(l_{1},MET)',                  regions='baseline'),
    #'MLL'             : HistOp1D(nBinsX=80, x0=50.0, x1=130, xUnits='GeV', xLabel='M_{ll}', logY=True),
    'MLL'             : HistOp1D(nBinsX=100, x0=0.0, x1=300.0, xUnits='GeV', xLabel='M_{ll}',                            regions='trig_only, ztautauCR, zCR_mumu, zCR_ee, baseline, baseline_emu, baseline_mue'),
    'ptll'            : HistOp1D(nBinsX=50, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='pT_{ll}',                           regions='baseline', logY=True),
    # MET + leptons
    'MET'             : HistOp1D(nBinsX=50, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='E_{T}^{miss}',                      regions='baseline, trig_only, zCR_mumu, zCR_ee, topCR, ztautauCR', logY=True),
    'METPhi'          : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='MET_{#phi}',                        regions='baseline'),
    'MCollASym'       : HistOp1D(nBinsX=25, x0=0.0,  x1=250.0, xUnits='GeV', xLabel='LFV Collinear Mass m_{coll}',       regions='baseline'),
    'dpt_ll'          : HistOp1D(nBinsX=20, x0=0.0,  x1=150.0, xUnits='GeV', xLabel='#Deltap_{T}^{ll}',                  regions='baseline'),
    'DphiLep0MET'     : HistOp1D(nBinsX=63, x0=-3.15,  x1=3.15,  xUnits='',    xLabel='#Delta#phi(l_{0},MET)',             regions='baseline', add_underflow=True, logY=True),
    'DphiLep1MET'     : HistOp1D(nBinsX=63, x0=-3.15,  x1=3.15,  xUnits='',    xLabel='#Delta#phi(l_{1},MET)',             regions='baseline', add_underflow=True, logY=True),
    'tau_pT'          : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{subleading lep + MET}',      regions='baseline'),
    'taulep1_pT_ratio': HistOp1D(nBinsX=20, x0=0.0,  x1=3,     xUnits='',    xLabel='p_{T}^{subleading lep + MET} / p_{T}^{leading lep}',  regions='baseline'),
    # Jets
    'preJet_pt'       : HistOp1D(nBinsX=20, x0=0.0,  x1=100.0, xUnits='GeV', xLabel='p_{T}^{preJet}',                    regions='no_sel'),
    'preJet_eta'      : HistOp1D(nBinsX=100, x0=-5.0, x1=5.0,   xUnits='',    xLabel='#eta_{preJet}',                     regions='no_sel'),
    'preJet_JVT'      : HistOp1D(nBinsX=39, x0=-0.2,  x1=1.1,     xUnits='',    xLabel='preJet JVT (|eta|<=2.4 & pT < 60)', regions='no_sel', logY=True),
    'baseJet_eta'     : HistOp1D(nBinsX=100, x0=-5.0, x1=5.0,   xUnits='',    xLabel='#eta_{baseJet}',                    regions='no_sel'),
    'baseJet_mv2c10'  : HistOp1D(nBinsX=80, x0=-2,   x1=2,    xUnits='',    xLabel='mv2c10_{baseJet}',                  regions='no_sel', logY=True),
    'jetN'            : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{base jet}',                      regions='baseline'),
    'jet_N2p4Eta25Pt' : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{jet} (p_{T}>25GeV, |#eta|<2.5)', regions='baseline'),
    'signalJetN'      : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{sig jets}',                      regions='baseline'),
    'jetN_g30'        : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{jet} (p_{T}>30GeV)',             regions='baseline'),
    'nLJets'          : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{CL jet}',                        regions='trig_only, baseline', logY=True),
    'nBJets'          : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{CB jet}',                        regions='trig_only, topCR, baseline, baseline_emu, baseline_mue', logY=True),
    'btag'            : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5,   xUnits='',    xLabel='B-tagged jet',                      regions='baseline'),
    'nForwardJets'    : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{F jet}',                         regions='baseline, baseline_emu, baseline_mue'),
    'j_pt[0]'         : HistOp1D(nBinsX=20, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='p_{T}^{leading jet}',               regions='baseline, topCR', logY=True),
    'j_pt[1]'         : HistOp1D(nBinsX=20, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='p_{T}^{subleading jet}',            regions='baseline, topCR', logY=True),
    'j_pt[2]'         : HistOp1D(nBinsX=20, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='p_{T}^{3rd leading jet}',           regions='baseline', logY=True),
    'j_pt[3]'         : HistOp1D(nBinsX=20, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='p_{T}^{4th leading jet}',           regions='baseline', logY=True),
    'j_pt'            : HistOp1D(nBinsX=20, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='Jet p_{T}',                         regions='baseline, baseline_emu, baseline_mue'),
    'mjj'             : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='GeV', xLabel='Dijet mass',                        regions='baseline'),
    'dEtaJJ'          : HistOp1D(nBinsX=20, x0=0.0,  x1=6.0,   xUnits='',    xLabel='#Delta#eta(j0,j1)',                 regions='baseline'),
    'j_eta'           : HistOp1D(nBinsX=100, x0=-5.0, x1=5.0,   xUnits='',    xLabel='Jet #eta',                          regions='baseline, baseline_emu, baseline_mue'),
    'j_jvt'           : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Jet JVT',                           regions='baseline'),
    'j_jvf'           : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Jet JVF',                           regions='baseline'),
    'j_phi'           : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='Jet #phi',                          regions='baseline'),
    'j_flav'          : HistOp1D(nBinsX=5,  x0=-0.5, x1=4.5,   xUnits='',    xLabel='Jet flavor (0:NA,1:CL,2:CB,3:F)',   regions='baseline'),
    # Leptons
    'preEl_EcaloClus'           : HistOp1D(nBinsX=50, x0= 0,   x1=50, xUnits='GeV', xLabel='E_{preEl CaloCluster}',          regions='no_sel', add_overflow=False),
    'preEl_etaCaloClus'         : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0, xUnits='',    xLabel='#eta_{preEl CaloCluster}',       regions='no_sel'),
    'baseEl_etconetopo20'       : HistOp1D(nBinsX=50, x0= -7,   x1=15, xUnits='GeV', xLabel='E_{T}^{baselineEl} conetopo20' ,   regions='no_sel, baseline_emu, baseline_mue', logY=False),
    'baseEl_ptvarcone20'        : HistOp1D(nBinsX=25, x0= -1,   x1=4, xUnits='GeV', xLabel='p_{T}^{baselineEl} ptvarcone20' ,  regions='no_sel, baseline_emu, baseline_mue', logY=False),
    'baseMu_etconetopo20'       : HistOp1D(nBinsX=60, x0= -3,   x1=3, xUnits='GeV', xLabel='E_{T}^{baselineMu} conetopo20' ,   regions='no_sel, baseline_emu, baseline_mue'),
    'baseMu_ptvarcone30'        : HistOp1D(nBinsX=60, x0= -3,   x1=3, xUnits='GeV', xLabel='p_{T}^{baselineMu} ptvarcone30' ,  regions='no_sel, baseline_emu, baseline_mue', logY=False),    
    'el1pT_trackclus_ratio'     : HistOp1D(nBinsX=30, x0=0,    x1=3,   xUnits='', xLabel='el_{subleading pT}^{track} / el_{subleading pT}^{cluster}', regions='trig_only, baseline_mue', logY=True),

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
    #'Z_Lep2_pT'         : HistOp1D(nBinsX=40, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='3rd leading lepton p_{T}', logY=True),
    'Z_Lep2_pT'         : HistOp1D(nBinsX=50, x0=0.0,  x1=50.0, xUnits='GeV', xLabel='3rd leading lepton p_{T}', logY=False),
    'Z_Lep2_eta'        : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='3rd leading lepton #eta', logY=True),
    'Z_dilep_sign'       :HistOp1D(nBinsX=5, x0=-2.5, x1=2.5,   xUnits='',    xLabel='Z Dilepton Sign : OS(-1) SS(1)', logY=True),
    'Z2_dilep_sign'       :HistOp1D(nBinsX=5, x0=-2.5, x1=2.5,   xUnits='',    xLabel='2nd Z Dilepton Sign : OS(-1) SS(1)', logY=True),
    'Z_Lep2_dPhi_MET'     : HistOp1D(nBinsX=63, x0=-3.15,  x1=3.15,  xUnits='',    xLabel='#Delta#phi(l_{3},MET)', add_underflow=True, logY=True),
    'Z_Lep2_mT'           : HistOp1D(nBinsX=50, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='3rd leading lepton m_{T}', logY=False),

}

## Backgrounds to include
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
backgrounds.append(wgamma)
backgrounds.append(htt)
backgrounds.append(hww)
#backgrounds.append(signal)
#data = None

region_ops = []
#region_ops += ['wjets_FF_CRden_emu', 'wjets_FF_CRden_mue']
#region_ops += ['wjets_FF_CRnum_emu', 'wjets_FF_CRnum_mue']
#region_ops += ['zjets_FF_CRden_eee', 'zjets_FF_CRden_mme']
#region_ops += ['zjets_FF_CRden_eem', 'zjets_FF_CRden_mmm']
#region_ops += ['zjets_FF_CRnum_eee', 'zjets_FF_CRnum_mme']
#region_ops += ['zjets_FF_CRnum_eem', 'zjets_FF_CRnum_mmm']
region_ops += ['baseline_mue', 'baseline_emu']

vars_to_plot = []
## Variable to get Yields
vars_to_plot += ['treatAsYear']

## Quick Plots
#vars_to_plot += ['Lep0Pt', 'Lep1Pt']
#vars_to_plot += ['Z_Lep2_pT']

## Z CR
#vars_to_plot += ['Lep0Pt', 'Lep1Pt', 'MLL']

## Top CR
#vars_to_plot += ['nBJets', 'Lep0Pt', 'Lep1Pt', 'j_pt[0]', 'j_pt[1]', 'MET']

## Baseline Selection
#vars_to_plot += ['Lep0Pt', 'Lep1Pt', 'MLL', 'MET']

## Trigger Variables
#vars_to_plot += ["Lep0Pt", "Lep0Phi", "Lep0Eta", "Lep1Pt", "Lep1Phi", "Lep1Eta"]

## Fake Rejection
#vars_to_plot += ["El_ID", "Mu_ID"] 
#vars_to_plot += ["Lep_Iso"]
#vars_to_plot += ["baseEl_etconetopo20", "baseEl_ptvarcone20", "baseMu_ptvarcone30", "baseMu_etconetopo20"]
#vars_to_plot += ["mu_z0SinTheta", "mu_d0sigBSCorr", "el_z0SinTheta", "el_d0sigBSCorr"]
#vars_to_plot += ["isCaloTagged", "dRy_sEl_bMu_noCalo", "dRy_sEl_bMu_Calo", 'drll']

## Baseline Selection
#vars_to_plot += ['Lep0Pt', 'Lep1Pt', 'MLL', 'nBJets', 'el1pT_trackclus_ratio', 'MET']

## Multiplicity Plots
#vars_to_plot += ['nLJets', 'nBJets', 'nForwardJets', 'n_leptons']

## Basic Kinematics
#vars_to_plot += ['Lep0Pt', 'Lep1Pt', 'MET', 'j_pt', 'Lep0Eta', 'Lep1Eta', 'j_eta',] 

## Fake plots
#vars_to_plot += ['nLepID', 'nLepAntiID', 'nBJets']
    # Wjets
#vars_to_plot += ['Lep0Pt', 'aID_Lep0Pt', 'aID_MLL', 'aID_drll', 'aID_Lep0Eta', 'MET']
#vars_to_plot += ['Lep0Pt', 'Lep1Pt', 'MLL', 'drll', 'Lep1Eta', 'MET', 'dilep_flav', 'LepLepSign']
    # Z+Jets
#vars_to_plot += ['Z_MLL','nBJets','Z_Lep2_mT']
#vars_to_plot += ['Z_MLL', 'Z_Lep2_pT', 'dR_Zl', 'Z_Lep2_eta', 'MET', 'DphiLep0MET', 'Z_Lep2_dPhi_MET', 'Z2_dilep_sign', 'Z2_MLL']

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
    #for region in ops.regions.split(','):
    for region in region_ops:
        region = region.strip()
        print "Plotting for region", region
        p = plot.Plot1D()
        p.initialize(region, var, name="%s_%s"%(region, name_))
        p.doLogY = ops.logY
        p.add_overflow = ops.add_overflow
        p.labels(x=ops.xLabel, y=yLabel)
        p.xax(bin_width, ops.x0, ops.x1)
        p.yax(y0, y1)
        assert data or backgrounds, "No data or backgrounds defined"
        if data and backgrounds:
            p.setRatioCanvas(p.name)
        else:
            p.setDefaultCanvas(p.name)
            
        plots.append(p)


