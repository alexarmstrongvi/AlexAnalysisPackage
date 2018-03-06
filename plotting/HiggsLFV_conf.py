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

lumi_ = [36180, 36.01, 3209, 32971] #[2015-16 (ipb), 2015-16 (ifb), 2015 (ipb), 2016 (ipb)]
lumi_val = 0

backgrounds = []
######## MC
## ttbar
ttbar = background.Background("ttbar", "t#bar{t}")
ttbar.set_debug()
ttbar.scale_factor = lumi_[lumi_val]#
ttbar.set_fillStyle(0)
ttbar.setLineStyle(1)
ttbar.set_color(ROOT.kOrange+2)
ttbar.set_treename("ttbar")
ttbar.set_chain_from_list_CONDOR(filelist_dir+ "ttbar/", rawdir)
backgrounds.append(ttbar)

# singletop
stop = background.Background("st", "Single top")
stop.set_debug()
stop.scale_factor = lumi_[lumi_val]
stop.set_fillStyle(0)
stop.setLineStyle(1)
stop.set_color(ROOT.kOrange+1)
stop.set_treename("ST")
stop.set_chain_from_list_CONDOR(filelist_dir+ "singletop/", rawdir)
backgrounds.append(stop)

# WW
WW = background.Background("ww", "WW")
WW.set_debug()
WW.scale_factor = lumi_[lumi_val] #* 1.06
WW.set_fillStyle(0)
WW.setLineStyle(1)
WW.set_color(ROOT.kSpring-6)
WW.set_treename("WW")
WW.set_chain_from_list_CONDOR(filelist_dir+ "ww/", rawdir)
backgrounds.append(WW)

# ZZ
ZZ = background.Background("zz", "ZZ")
ZZ.set_debug()
ZZ.scale_factor = lumi_[lumi_val] #* 1.06
ZZ.set_fillStyle(0)
ZZ.setLineStyle(1)
ZZ.set_color(ROOT.kSpring-4)
ZZ.set_treename("ZZ")
ZZ.set_chain_from_list_CONDOR(filelist_dir+ "zz/", rawdir)
backgrounds.append(ZZ)

# WZ
WZ = background.Background("wz", "WZ")
WZ.set_debug()
WZ.scale_factor = lumi_[lumi_val] #* 1.06
WZ.set_fillStyle(0)
WZ.setLineStyle(1)
WZ.set_color(ROOT.kSpring-5)
WZ.set_treename("WZ")
WZ.set_chain_from_list_CONDOR(filelist_dir+ "wz/", rawdir)
backgrounds.append(WZ)

## diboson
#diboson = background.Background("vv", "VV")
#diboson.set_debug()
#diboson.scale_factor = lumi_[lumi_val] #* 1.06
#diboson.set_fillStyle(0)
#diboson.setLineStyle(1)
#diboson.set_color(ROOT.kSpring-6)
#diboson.set_treename("diboson")
#diboson.set_chain_from_list_CONDOR(filelist_dir+ "diboson/", rawdir)
##diboson.set_chain_from_list_CONDOR(filelist_dir+ "diboson_new_check/", diboson_rawdir_SF)
#backgrounds.append(diboson)

# Zll
zll = background.Background("zll", "Zll")
zll.set_debug()
zll.scale_factor = lumi_[lumi_val]
zll.set_fillStyle(0)
zll.setLineStyle(1)
zll.set_color(ROOT.kAzure-9)
zll.set_treename("zll")
zll.set_chain_from_list_CONDOR(filelist_dir+ "zll/", rawdir)
backgrounds.append(zll)

# Ztt
ztt = background.Background("ztt", "Z#tau#tau")
ztt.set_debug()
ztt.scale_factor = lumi_[lumi_val]
ztt.set_fillStyle(0)
ztt.setLineStyle(1)
ztt.set_color(ROOT.kAzure-5)
ztt.set_treename("ztt")
ztt.set_chain_from_list_CONDOR(filelist_dir+ "ztt/", rawdir)
backgrounds.append(ztt)

# Wjets
wjets = background.Background("wjets", "W+jets")
wjets.set_debug()
wjets.scale_factor = lumi_[lumi_val]
wjets.set_fillStyle(0)
wjets.setLineStyle(1)
wjets.set_color(ROOT.kOrange)
wjets.set_treename("wjets")
wjets.set_chain_from_list_CONDOR(filelist_dir+ "wjets/", rawdir)
backgrounds.append(wjets)

## Higgs -> tau tau
htt = background.Background("htt", "H#tau#tau")
#higgs = background.Background("higgs", "Higgs")
htt.scale_factor = lumi_[lumi_val]
htt.set_fillStyle(0)
htt.setLineStyle(1)
htt.set_color(ROOT.kRed)
htt.set_treename("htt")
htt.set_chain_from_list_CONDOR(filelist_dir+ "htt/", rawdir)
backgrounds.append(htt)

## Higgs -> W W
hww = background.Background("hww", "HWW")
#higgs = background.Background("higgs", "Higgs")
hww.scale_factor = lumi_[lumi_val]
hww.set_fillStyle(0)
hww.setLineStyle(1)
hww.set_color(ROOT.kBlue+3)
hww.set_treename("hww")
hww.set_chain_from_list_CONDOR(filelist_dir+ "hww/", rawdir)
backgrounds.append(hww)

#fakes = background.Background("fakes", "FNP")
#fakes.scale_factor = 1.0 # Feb 7 2017 - fakes are full dataset, yo
##fakes.scale_factor = 2.95 # scale from 12.2/fb to 36/fb
#fakes.set_treename("superNt")
##fakes.set_file(fake_rawdir + "physics_Main_276262_303560_FakesInclusive.root")
##fakes.set_file(fake_rawdir + "fakes_3body.root")
#fakes.set_file(fake_rawdir + "fakes_3body_mar10_361ifb.root")
#fakes.set_merged_tree("superNt")
#fakes.set_color(94) # stop-2l uglify for Moriond
##fakes.set_color(r.kOrange+7)
#fakes.set_fillStyle(0)
#fakes.setLineStyle(1)
#backgrounds.append(fakes)

# Samples currently need to be remade
#s0 = background.Background("higgs_lfv", "Higgs LFV")
#s0.setSignal()
#s0.scale_factor = lumi_[lumi_val] * 10
#s0.set_fillStyle(0)
#s0.set_color(ROOT.kGreen)
#s0.set_treename("s0")
#s0.set_chain_from_list_CONDOR(filelist_dir + "higgs_lfv", signal_rawdir)
#backgrounds.append(s0)


#### DATA
data = background.Data()
data.set_color(r.kBlack)
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
lep_eta_cut = 'fabs(l_eta[0]) <= 2.4 && fabs(l_eta[1]) <= 2.4'

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
mu15_trig   = '(HLT_mu20_iloose_L1MU15 || HLT_mu50)'
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

# Combined triggers
e15_trig_pT = '(%s || %s)'%(e15_trig_emu_pT, e15_trig_mue_pT)
mu15_trig_pT = '(%s || %s)'%(mu15_trig_emu_pT, mu15_trig_mue_pT)
e16_trig_pT = '(%s || %s)'%(e16_trig_emu_pT, e16_trig_mue_pT)
mu16_trig_pT = '(%s || %s)'%(mu16_trig_emu_pT, mu16_trig_mue_pT)
emu15_trig_pT = '(%s || %s)'%(emu15_trig_emu_pT, emu15_trig_mue_pT)
emu16_trig_pT = '(%s || %s)'%(emu16_trig_emu_pT, emu16_trig_mue_pT)

dilep15_trig      = '(treatAsYear==2015 && %s)'%emu15_trig_pT
dilep16_trig     = '(treatAsYear==2016 && %s)'%emu16_trig_pT
singlelep15_trig = '(treatAsYear==2015 && (%s || %s))'%(e15_trig_pT,mu15_trig_pT)
singlelep16_trig = '(treatAsYear==2016 && (%s || %s))'%(e16_trig_pT,mu16_trig_pT)

dilep_trig = '(%s || %s)'%(dilep15_trig, dilep16_trig)
singlelep_trig = '(%s || %s)'%(singlelep15_trig, singlelep16_trig)


# Region building blocks
Baseline_Sel = 'l_pt[0] >= 45 && l_pt[1] >= 15 '\
              + '&& 30 < MLL && MLL < 150 '\
              + '&& nCentralBJets==0 '\
              + '&& (dilep_flav != 0 || (el1pT_trackclus_ratio < 1.2)) '\
              + '&&' + singlelep_trig
VBF_stripped = "JetN_g30 >= 2 && j_pt[0] > 40 && Mjj > 400 && DEtaJJ > 3"

# Define regions

reg = region.Region()
reg.name = "no_sel"
reg.displayname = "No Selections"
reg.tcut = "1"
regions.append(reg)
reg = region.Region()

reg.name = "trig_only"
reg.displayname = "Single Lepton Triggers"
reg.tcut = singlelep_trig
regions.append(reg)

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
reg.tcut = "(%s) && !(%s) && DphiLep1MET < 1 && MtLep0 > 50 && MtLep1 < 40 && ((MET+l_pt[1])/l_pt[1]) > 0.5"%(Baseline_Sel, VBF_stripped)
regions.append(reg)


##########################################
# plots
##########################################
plots = []

vars = {}

### Define variables for plots
# Standard kinematics
HistOp1D = namedtuple('HistOp1D', 'regions,   nBinsX, x0, x1, y0,   y1,   logY,   xUnits, xLabel, yLabel')
HistOp1D.__new__.__defaults__=   ('baseline', 12,     -2,  10, None, None, False,  '',     '',     'Events')

HistOpMap = {
    # Event level
    'RunNumber'       : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Event run number',                  regions='baseline'),
    'event_number'    : HistOp1D(nBinsX=1000, x0=0,  x1=1000000000,    xUnits='',    xLabel='Event number',              regions='baseline_emu, baseline_mue'),
    'isMC'            : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5,   xUnits='',    xLabel='is Monte Carlo',                    regions='baseline'),
    'eventweight'     : HistOp1D(nBinsX=100, x0=-0.001,  x1=0.002,    xUnits='',    xLabel='Event weight',                      regions='no_sel, baseline', logY=True),
    'dsid'            : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Sample DSID',                       regions='baseline'),
    'treatAsYear'     : HistOp1D(nBinsX=11, x0=2007.5,x1=2018.5,xUnits='',   xLabel='treatAsYear',                       regions='no_sel'),
    # Multiplicity
    'n_baseLeptons'   : HistOp1D(nBinsX=11, x0=-0.5, x1=10.5,  xUnits='',    xLabel='N_{baseline leptons}',              regions='baseline'),
    'n_leptons'       : HistOp1D(nBinsX=11, x0=-0.5, x1=10.5,  xUnits='',    xLabel='N_{signal leptons}',                regions='baseline'),
    # Leptons
    'l_pt'            : HistOp1D(nBinsX=25, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='Lepton p_{T}',                      regions='baseline'),
    ## Electrons
    'el1_track_pt'    : HistOp1D(nBinsX=25, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='Leading electron track p_{T}',      regions='baseline'),
    'el1_clus_pt'     : HistOp1D(nBinsX=25, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='Leading electron cluster p_{T}',    regions='baseline'),
    'preEl_Et'        : HistOp1D(nBinsX=20, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='E_{T}^{preElectron}',               regions='trig_only'),
    'preEl_pt'        : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{preElectron}',               regions='trig_only'),
    'preEl_clusEtaBE' : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta_{preEl clusterBE}',            regions='trig_only'),
    'preEl_eta'       : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta_{preElectron}',                regions='trig_only'),
    'baseEl_ID'       : HistOp1D(nBinsX=9,  x0=-0.5, x1=8.5,   xUnits='',    xLabel='el_{Baseline}^{ID}(non-inclusive)', regions='trig_only'),
    'el_type'         : HistOp1D(nBinsX=9,  x0=-0.5, x1=8.5,   xUnits='',    xLabel='electron truth type',               regions='trig_only'),
    'el_origin'       : HistOp1D(nBinsX=9,  x0=-0.5, x1=8.5,   xUnits='',    xLabel='electron truth origin',             regions='trig_only'),
    ## Muons
    'preMu_pt'        : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{preMuon}',                   regions='trig_only'),
    'preMuon_ID'      : HistOp1D(nBinsX=9,  x0=-0.5, x1=8.5,   xUnits='',    xLabel='mu_{pre}^{ID} (non-inclusive)',     regions='trig_only'),
    'baseMu_pt'       : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{BaselineMuon}',              regions='trig_only'),
    'baseMu_eta'      : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta_{BaselineMuon}',               regions='trig_only'),
    'baseMuon_ID'     : HistOp1D(nBinsX=9,  x0=-0.5, x1=8.5,   xUnits='',    xLabel='mu_{Baseline}^{ID}(non-inclusive)', regions='trig_only'),
    'mu_type'         : HistOp1D(nBinsX=9,  x0=-0.5, x1=8.5,   xUnits='',    xLabel='muon truth type',                   regions='trig_only'),
    'mu_origin'       : HistOp1D(nBinsX=9,  x0=-0.5, x1=8.5,   xUnits='',    xLabel='muon truth origin',                 regions='trig_only'),
    ## Taus
    'preTau_q'        : HistOp1D(nBinsX=3,  x0=-1.5, x1=1.5,   xUnits='',    xLabel='Tau charge',                        regions='trig_only'),
    'preTau_nTracks'  : HistOp1D(nBinsX=3,  x0=-1.5, x1=1.5,   xUnits='',    xLabel='preTau nTracks',                    regions='trig_only'),
    'baseTau_pT'      : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{BaselineTau}',               regions='trig_only'),
    'baseTau_eta'     : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta_{BaselineMuon}',               regions='trig_only'),
    'baseTau_nTracks' : HistOp1D(nBinsX=3,  x0=-1.5, x1=1.5,   xUnits='',    xLabel='Baseline Tau nTracks',              regions='trig_only'),
    'baseTau_ID'      : HistOp1D(nBinsX=9,  x0=-0.5, x1=8.5,   xUnits='',    xLabel='tau_{Baseline}^{ID}(non-inclusive)',regions='trig_only'),

    ## General lepton
    'l_eta'           : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='Lepton #eta',                       regions='baseline'),
    'l_phi'           : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='Lepton #phi',                       regions='baseline'),
    'l_flav'          : HistOp1D(nBinsX=5,  x0=-0.5, x1=4.5,   xUnits='',    xLabel='Lepton flavor (0: e, 1: m)',        regions='baseline'),
    'l_type'          : HistOp1D(nBinsX=5,  x0=-0.5, x1=4.5,   xUnits='',    xLabel='Lepton type',                       regions='baseline'),
    'l_origin'        : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Lepton origin',                     regions='baseline'),
    'l_q'             : HistOp1D(nBinsX=3,  x0=-1.5, x1=1.5,   xUnits='',    xLabel='Lepton charge',                     regions='baseline'),
    'LepLepSign'      : HistOp1D(nBinsX=3,  x0=-1.5, x1=1.5,   xUnits='',    xLabel='Leptons sign product',              regions='baseline'),
    'Lep0Pt'          : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{leading lep}',               regions='trig_only, baseline'),
    'Lep1Pt'          : HistOp1D(nBinsX=20, x0=0.0,  x1=100.0, xUnits='GeV', xLabel='p_{T}^{subleading lep}',            regions='trig_only, baseline'),
    'Lep0Eta'         : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta^{leading lep}',                regions='trig_only, baseline'),
    'Lep1Eta'         : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta^{subleading lep}',             regions='trig_only, baseline'),
    'Lep0Phi'         : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='#phi^{leading lep}',                regions='baseline'),
    'Lep1Phi'         : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='#phi^{subleading lep}',             regions='baseline'),
    'MLep0'           : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='GeV', xLabel='M_{l0}',                            regions='baseline'),
    'MLep1'           : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='GeV', xLabel='M_{l1}',                            regions='baseline'),
    'DEtaLL'          : HistOp1D(nBinsX=20, x0=0.0,  x1=6.0,   xUnits='',    xLabel='#Delta#eta_{ll}',                   regions='baseline'),
    'DphiLL'          : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='#Delta#phi_{ll}',                   regions='baseline'),
    'drll'            : HistOp1D(nBinsX=20, x0=0.0,  x1=6.0,   xUnits='',    xLabel='#DeltaR_{ll}',                      regions='baseline'),
    'dilep_flav'      : HistOp1D(nBinsX=5,  x0=-0.5, x1=4.5,   xUnits='',    xLabel='Dilepton flavor',                   regions='baseline'),
    'isEM'            : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5,   xUnits='',    xLabel='Dilepton flavor is el mu',          regions='baseline'),
    'isME'            : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5,   xUnits='',    xLabel='Dilepton flavor is mu el',          regions='baseline'),
    'MtLep0'          : HistOp1D(nBinsX=15, x0=0.0,  x1=250.0, xUnits='GeV', xLabel='m_{T}(l_{0},MET)',                  regions='baseline'),
    'MtLep1'          : HistOp1D(nBinsX=20, x0=0.0,  x1=140.0, xUnits='GeV', xLabel='m_{T}(l_{1},MET)',                  regions='baseline'),
    'MLL'             : HistOp1D(nBinsX=28, x0=20.0, x1=300.0, xUnits='GeV', xLabel='M_{ll}',                            regions='no_sel, trig_only, baseline, baseline_mue, baseline_emu'),
    'ptll'            : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='pT_{ll}',                           regions='baseline'),
    # MET + leptons
    'MET'             : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='E_{T}^{miss}',                      regions='baseline'),
    'METPhi'          : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='MET_{#phi}',                        regions='baseline'),
    'MCollASym'       : HistOp1D(nBinsX=25, x0=0.0,  x1=250.0, xUnits='GeV', xLabel='LFV Collinear Mass m_{coll}',       regions='baseline'),
    'dpt_ll'          : HistOp1D(nBinsX=20, x0=0.0,  x1=150.0, xUnits='GeV', xLabel='#Deltap_{T}^{ll}',                  regions='baseline'),
    'DphiLep0MET'     : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='#Delta#phi(l_{0},MET)',             regions='baseline'),
    'DphiLep1MET'     : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='#Delta#phi(l_{1},MET)',             regions='baseline'),
    'tau_pT'          : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{subleading lep + MET}',      regions='baseline'),
    'taulep1_pT_ratio': HistOp1D(nBinsX=20, x0=0.0,  x1=3,     xUnits='',    xLabel='p_{T}^{subleading lep + MET} / p_{T}^{leading lep}',  regions='baseline'),
    # Jets
    'preJet_pt'       : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{preJet}',                    regions='trig_only'),
    'preJet_eta'      : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta_{preJet}',                     regions='trig_only'),
    'preJet_JVT'      : HistOp1D(nBinsX=20, x0=0.0,  x1=3,     xUnits='',    xLabel='preJet JVT (|eta|<=2.4 & pT < 60)', regions='trig_only'),
    'baseJet_eta'     : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta_{baseJet}',                    regions='trig_only'),
    'baseJet_mv2c10'  : HistOp1D(nBinsX=20, x0=-1,   x1=10,    xUnits='',    xLabel='mv2c10_{baseJet}',                  regions='trig_only'),
    'JetN'            : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{base jet}',                      regions='baseline'),
    'Jet_N2p4Eta25Pt' : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{jet} (p_{T}>25GeV, |#eta|<2.5)', regions='baseline'),
    'SignalJetN'      : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{sig jets}',                      regions='baseline'),
    'JetN_g30'        : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{jet} (p_{T}>30GeV)',             regions='baseline'),
    'nCentralLJets'   : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{CL jet}',                        regions='baseline'),
    'nCentralBJets'   : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{CB jet}',                        regions='trig_only, baseline'),
    'Btag'            : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5,   xUnits='',    xLabel='B-tagged jet',                      regions='baseline'),
    'nForwardJets'    : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{F jet}',                         regions='baseline'),
    'j_pt[0]'         : HistOp1D(nBinsX=20, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='p_{T}^{leading jet}',               regions='baseline'),
    'j_pt'            : HistOp1D(nBinsX=20, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='Jet p_{T}',                         regions='baseline'),
    'Mjj'             : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='GeV', xLabel='Dijet mass',                        regions='baseline'),
    'DEtaJJ'          : HistOp1D(nBinsX=20, x0=0.0,  x1=6.0,   xUnits='',    xLabel='#Delta#eta(j0,j1)',                 regions='baseline'),
    'j_eta'           : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='Jet #eta',                          regions='baseline'),
    'j_jvt'           : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Jet JVT',                           regions='baseline'),
    'j_jvf'           : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Jet JVF',                           regions='baseline'),
    'j_phi'           : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='Jet #phi',                          regions='baseline'),
    'j_flav'          : HistOp1D(nBinsX=5,  x0=-0.5, x1=4.5,   xUnits='',    xLabel='Jet flavor (0:NA,1:CL,2:CB,3:F)',   regions='baseline'),
    # Leptons
    'preEl_EcaloClus'           : HistOp1D(nBinsX=20, x0= 0,   x1=200, xUnits='GeV', xLabel='E_{preEl CaloCluster}',          regions='trig_only'),
    'preEl_etaCaloClus'         : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0, xUnits='',    xLabel='#eta_{preEl CaloCluster}',       regions='trig_only'),
    'baseEl_etconetopo20'       : HistOp1D(nBinsX=20, x0= -5,   x1=100, xUnits='GeV', xLabel='E_T^{baselineEl} conetopo20' ,   regions='trig_only'),
    'baseEl_ptvarcone20'        : HistOp1D(nBinsX=20, x0= 0,   x1=200, xUnits='GeV', xLabel='p_T^{baselineEl} ptvarcone20' ,  regions='trig_only'),
    'el1pT_trackclus_ratio'     : HistOp1D(nBinsX=25, x0=0,    x1=3,   xUnits='', xLabel='el_{subleading pT}^{track} / el_{subleading pT}^{cluster}', regions='trig_only, baseline_mue'),
    'baseMu_etconetopo20'       : HistOp1D(nBinsX=20, x0= 0,   x1=200, xUnits='GeV', xLabel='E_T^{baselineMu} conetopo20' ,   regions='trig_only'),
    'baseMu_ptvarcone30'        : HistOp1D(nBinsX=20, x0= 0,   x1=200, xUnits='GeV', xLabel='p_T^{baselineMu} ptvarcone30' ,  regions='trig_only'),

    'HLT_mu18_mu8noL1'                 : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='mu18 mu8noL1 trigger',                 regions='baseline'),
    'HLT_2e12_lhloose_L12EM10VH'       : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='2e12 lhloose L12EM10VH trigger',       regions='baseline'),
    'HLT_e17_lhloose_mu14'             : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='e17 lhloose mu14 trigger',             regions='baseline'),
    'HLT_mu20_mu8noL1'                 : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='mu20 mu8noL1 trigger',                 regions='baseline'),
    'HLT_2e15_lhvloose_L12EM13VH'      : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='2e15 lhvloose L12EM13VH trigger',      regions='baseline'),
    'HLT_2e17_lhvloose_nod0'           : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='2e17 lhvloose nod0 trigger',           regions='baseline'),
    'HLT_mu22_mu8noL1'                 : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='mu22 mu8noL1 trigger',                 regions='baseline'),
    'HLT_e17_lhloose_nod0_mu14'        : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='e17 lhloose nod0 mu14 trigger',        regions='baseline'),
    'HLT_e60_lhmedium'                 : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='e60 lhmedium trigger',                 regions='baseline'),
    'HLT_e24_lhmedium_L1EM20VH'        : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='e24 lhmedium L1EM20VH trigger',        regions='baseline'),
    'HLT_e24_lhmedium_iloose_L1EM18VH' : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='e24 lhmedium iloose L1EM18VH trigger', regions='baseline'),
    'HLT_mu20_iloose_L1MU15'           : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='mu20 iloose L1MU15 trigger',           regions='baseline'),
    'HLT_mu24_imedium'                 : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='mu24 imedium trigger',                 regions='baseline'),
    'HLT_mu26_imedium'                 : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='mu26 imedium trigger',                 regions='baseline'),
    'HLT_e24_lhtight_nod0_ivarloose'   : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='e24 lhtight nod0 ivarloose trigger',   regions='baseline'),
    'HLT_e26_lhtight_nod0_ivarloose'   : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='e26 lhtight nod0 ivarloose trigger',   regions='baseline'),
    'HLT_e60_lhmedium_nod0'            : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='e60 lhmedium nod0 trigger',            regions='baseline'),
    'HLT_e120_lhloose'                 : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='e120 lhloose trigger',                 regions='baseline'),
    'HLT_e140_lhloose_nod0'            : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='e140 lhloose nod0 trigger',            regions='baseline'),
    'HLT_mu24_iloose'                  : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='mu24 iloose trigger',                  regions='baseline'),
    'HLT_mu24_iloose_L1MU15'           : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='mu24 iloose L1MU15 trigger',           regions='baseline'),
    'HLT_mu24_ivarloose'               : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='mu24 ivarloose trigger',               regions='baseline'),
    'HLT_mu24_ivarloose_L1MU15'        : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='mu24 ivarloose L1MU15 trigger',        regions='baseline'),
    'HLT_mu24_ivarmedium'              : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='mu24 ivarmedium trigger',              regions='baseline'),
    'HLT_mu26_ivarmedium'              : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='mu26 ivarmedium trigger',              regions='baseline'),
    'HLT_mu50'                         : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5, xUnits='', xLabel='mu50 trigger',                         regions='baseline')
}

#vars_to_plot = ['eventweight', 'Lep0Pt', 'Lep0Eta', 'nCentralLJets']
#vars_to_plot = ['Lep0Pt', 'Lep1Pt', 'MET', 'MLL', 'nCentralLJets', 'nCentralBJets', 'j_pt[0]']
#vars_to_plot = ['treatAsYear', 'MLL']
vars_to_plot = ['eventweight']
#vars_to_plot = ['preEl_EcaloClus', 'preEl_etaCaloClus', 'preEl_Et', 'preEl_pt',
#                'preEl_clusEtaBE', 'preEl_eta', 'baseEl_etconetopo20',
#                'baseEl_ptvarcone20', 'baseEl_ID', 'el_type', 'el_origin',
#                'el1_track_pt', 'el1_clus_pt', 'el1pT_trackclus_ratio',
#                'preMu_pt', 'preMuon_ID', 'baseMu_pt', 'baseMu_eta',
#                'baseMu_etconetopo20', 'baseMu_ptvarcone30', 'baseMuon_ID',
#                'mu_type', 'mu_origin', 'preTau_q', 'preTau_nTracks',
#                'baseTau_pT', 'baseTau_eta', 'baseTau_nTracks', 'baseTau_ID',
#                'tau_pT', 'taulep1_pT_ratio', 'preJet_pt', 'preJet_eta',
#                'preJet_JVT', 'baseJet_eta', 'baseJet_mv2c10']

for var in vars_to_plot:
    ops = HistOpMap[var]

    # Set contingent defaults
    y0_def = 0.1 if ops.logY else 0
    y1_def = 1e7 if ops.logY else 5000

    # Determine plot properties
    name_ = re.sub(r'[(){}[\]]+','',var)
    bin_width = (ops.x1 - ops.x0) / float(ops.nBinsX)
    width_label = str(round(bin_width, 2))
    if not ops.xUnits:
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
        p.labels(x=ops.xLabel, y=yLabel)
        p.xax(bin_width, ops.x0, ops.x1)
        p.yax(y0, y1)
        p.setRatioCanvas(p.name)
        #p.setDefaultCanvas(p.name)
        plots.append(p)


