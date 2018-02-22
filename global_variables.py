#!/bin/bash/env python
from collections import OrderedDict, namedtuple, defaultdict
import ROOT


tag = "n0235"
analysis_dir = "/data/uclhc/uci/user/armstro1/SusyNt/analysis_%s/"%tag
analysis_run_dir = "/data/uclhc/uci/user/armstro1/SusyNt/analysis_%s_run/"%tag
#input_files = analysis_dir+'inputs_LFV/'
input_files = analysis_dir+'inputs_LFV_check/'
output_dir = analysis_run_dir+'outputs/'
dsid_dir = analysis_run_dir+'dsid_filelist/'
LOCAL_DSID_DIR = '/data/uclhc/uci/user/dantrim/susynt_productions/n0235/'
LOCAL_DSID_SUBDIRS = {'data':'data','mc':'mc'}

S2L_trigger = '(((pass_HLT_2e12_lhloose_L12EM10VH||pass_HLT_e17_lhloose_mu14||pass_HLT_mu18_mu8noL1)&&treatAsYear==2015)||((pass_HLT_2e17_lhvloose_nod0||pass_HLT_e17_lhloose_nod0_mu14||pass_HLT_mu22_mu8noL1)&&treatAsYear==2016))'
S2L_ptCuts  = 'Lep0Pt>25.&&Lep1Pt>20.&&MLL>40.'
S2L_isOS    = '(l_q[0]*l_q[1])<0'
S2L_jetVeto = 'nCentralBJets==0 && nForwardJets==0 && nCentralLJets==0'
DF_OS = 'dilep_flav <= 1 && l_q[0]*l_q[1]<0'
SF_OS = 'dilep_flav > 1 && l_q[0]*l_q[1]<0'
lep_eta_cut = 'fabs(l_eta[0]) <= 2.4 && fabs(l_eta[1]) <= 2.4'

# Same flavor lepton triggers
ee15_trig   = 'pass_HLT_2e12_lhloose_L12EM10VH'
mumu15_trig = 'pass_HLT_mu18_mu8noL1'
ee16_trig   = 'pass_HLT_2e17_lhvloose_nod0'
mumu16_trig = 'pass_HLT_mu22_mu8noL1'
# Different flavor lepton triggs
emu15_trig  = 'pass_HLT_e17_lhloose_mu14'
emu16_trig  = 'pass_HLT_e17_lhloose_nod0_mu14'
# Single lepton triggers
e15_trig    = '(pass_HLT_e24_lhmedium_L1EM20VH || pass_HLT_e60_lhmedium || pass_HLT_e120_lhloose)'
mu15_trig   = '(pass_HLT_mu20_iloose_L1MU15 || pass_HLT_mu50)'
e16_trig    = '(pass_HLT_e26_lhtight_nod0_ivarloose || pass_HLT_e60_lhmedium_nod0 || pass_HLT_e140_lhloose_nod0)'
mu16_trig   = '(pass_HLT_mu26_ivarmedium || pass_HLT_mu50)'
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

dilep15_trig     = '(treatAsYear==2015 && %s)'%emu15_trig_pT
dilep16_trig     = '(treatAsYear==2016 && %s)'%emu16_trig_pT
singlelep15_trig = '(treatAsYear==2015 && (%s || %s))'%(e15_trig_pT,mu15_trig_pT)
singlelep16_trig = '(treatAsYear==2016 && (%s || %s))'%(e16_trig_pT,mu16_trig_pT)

Sel = {
    ''             : '1',
    'all'       : '1',
    'dilep_trig': '(%s || %s)'%(dilep15_trig, dilep16_trig),
    'singlelep_trig': '(%s || %s)'%(singlelep15_trig,singlelep16_trig),
    'Baseline'  : 'Lep0Pt >= 45 && Lep1Pt >= 15 '
                   + '&& 30 < MLL && MLL < 150 '
                   + '&& nCentralBJets==0 '
                   + '&& (dilep_flav != 0 || (el0_track_pt/el0_clus_pt) < 1.2) ',
    'emu'       : 'dilep_flav == 0',
    'mue'       : 'dilep_flav == 1',
    'ee'        : 'dilep_flav == 2',
    'mumu'      : 'dilep_flav == 3',
    'tt'        : 'nLepID == 2', #'dilepID = 1',
    'lt'        : 'nLepID == 1 && nLepAntiID == 1', #'dilepID = 2',
    'tl'        : 'nLepID == 1 && nLepAntiID == 1', #'dilepID = 3',
    'll'        : 'nLepAntiID == 2', #'dilepID = 4',
    }
Sel['DF'] = "((%s) || (%s))"%(Sel['mue'], Sel['emu'])
Sel['SF'] = "((%s) || (%s))"%(Sel['mumu'], Sel['ee'])
Sel['VBF_raw'] = "JetN_g30 >= 2 && j_pt[0] > 40 && Mjj > 400 && DEtaJJ > 3"
Sel['VBF'] = "(%s) && (%s)"%(Sel['Baseline'],Sel['VBF_raw'])
Sel['Optimized'] = "(%s) && !(%s) && DphiLep1MET < 1 && MtLep0 > 50 && MtLep1 < 40 && ((MET+Lep1Pt)/Lep1Pt) > 0.5"%(Sel['Baseline'],Sel['VBF_raw'])
Sel['Fake_Base'] = "(nLepID+nLepAntiID)==2 && " + Sel['SF']
Sel['Fake'] = Sel['Fake_Base']+" && fabs(MLL-91.188)>10 && nCentralBJets==0"
trigger_sel_name = 'singlelep_trig'
trigger_selection = Sel[trigger_sel_name].replace("pass_","")

# Define the luminosity you want
luminosity = "36180" # ipb
#luminosity = "3209" # ipb for data15
#luminosity = "32971" # ipb for data16

# Signal Branching Ratio
BR = '0.01'

histMaxY = {}

HistOp1D = namedtuple('HistOp1D', 'nBinsX, x0, x1, xUnits, xLabel, yLabel')
HistOp1D.__new__.__defaults__= (25,0,-1,'','','Events')

HistOpMap = {
    'RunNumber'       : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Event run number',                  yLabel='Events'),
    'event_number'    : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Event number',                      yLabel='Events'),
    'isMC'            : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5,   xUnits='',    xLabel='is Monte Carlo',                    yLabel='Events'),
    'event_weight'    : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Event weight',                      yLabel='Events'),
    'dsid'            : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Sample DSID',                       yLabel='Events'),
    'treatAsYear'     : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='treatAsYear',                       yLabel='Events'),
    'nBaseLeptons'    : HistOp1D(nBinsX=11, x0=-0.5, x1=10.5,  xUnits='',    xLabel='N_{baseline leptons}',              yLabel='Events'),
    'nSignalLeptons'  : HistOp1D(nBinsX=11, x0=-0.5, x1=10.5,  xUnits='',    xLabel='N_{signal leptons}',                yLabel='Events'),
    'nBaseTaus'       : HistOp1D(nBinsX=11, x0=-0.5, x1=10.5,  xUnits='',    xLabel='N_{baseline taus}',                 yLabel='Events'),
    'nSignalTaus'     : HistOp1D(nBinsX=11, x0=-0.5, x1=10.5,  xUnits='',    xLabel='N_{signal taus}',                   yLabel='Events'),
    'l_pt'            : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='GeV', xLabel='Lepton p_{T}',                      yLabel='Events'),
    'el0_track_pt'    : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='GeV', xLabel='Leading electron track p_{T}',      yLabel='Events'),
    'el0_clus_pt'     : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='GeV', xLabel='Leading electron cluster p_{T}',    yLabel='Events'),
    'l_eta'           : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='Lepton #eta',                       yLabel='Events'),
    'l_phi'           : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='Lepton #phi',                       yLabel='Events'),
    'l_flav'          : HistOp1D(nBinsX=5,  x0=-0.5, x1=4.5,   xUnits='',    xLabel='Lepton flavor (0: e, 1: m)',        yLabel='Events'),
    'l_type'          : HistOp1D(nBinsX=5,  x0=-0.5, x1=4.5,   xUnits='',    xLabel='Lepton type',                       yLabel='Events'),
    'l_origin'        : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Lepton origin',                     yLabel='Events'),
    'l_q'             : HistOp1D(nBinsX=3,  x0=-1.5, x1=1.5,   xUnits='',    xLabel='Lepton charge',                     yLabel='Events'),
    'LepLepSign'      : HistOp1D(nBinsX=3,  x0=-1.5, x1=1.5,   xUnits='',    xLabel='Leptons sign product',              yLabel='Events'),
    'MET'             : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='E_{T}^{miss}',                      yLabel='Events'),
    'met_Et'          : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='E_{T}^{miss}',                      yLabel='Events'),
    'met_phi'         : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='MET_{#phi}',                        yLabel='Events'),
    'METPhi'          : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='MET_{#phi}',                        yLabel='Events'),
    'Lep0Pt'          : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='p_{T}^{leading lep}',               yLabel='Events'),
    'Lep1Pt'          : HistOp1D(nBinsX=20, x0=0.0,  x1=100.0, xUnits='GeV', xLabel='p_{T}^{subleading lep}',            yLabel='Events'),
    'Lep0Eta'         : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta^{leading lep}',                yLabel='Events'),
    'Lep1Eta'         : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='#eta^{subleading lep}',             yLabel='Events'),
    'Lep0Phi'         : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='#phi^{leading lep}',                yLabel='Events'),
    'Lep1Phi'         : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='#phi^{subleading lep}',             yLabel='Events'),
    'MLep0'           : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='GeV', xLabel='M_{l0}',                            yLabel='Events'),
    'MLep1'           : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='GeV', xLabel='M_{l1}',                            yLabel='Events'),
    'DEtaLL'          : HistOp1D(nBinsX=20, x0=0.0,  x1=6.0,   xUnits='',    xLabel='#Delta#eta_{ll}',                   yLabel='Events'),
    'DphiLL'          : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='#Delta#phi_{ll}',                   yLabel='Events'),
    'drll'            : HistOp1D(nBinsX=20, x0=0.0,  x1=6.0,   xUnits='',    xLabel='#DeltaR_{ll}',                      yLabel='Events'),
    'dilep_flav'      : HistOp1D(nBinsX=5,  x0=-0.5, x1=4.5,   xUnits='',    xLabel='Dilepton flavor',                   yLabel='Events'),
    'isEM'            : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5,   xUnits='',    xLabel='Dilepton flavor is el mu',          yLabel='Events'),
    'isME'            : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5,   xUnits='',    xLabel='Dilepton flavor is mu el',          yLabel='Events'),
    'MCollASym'       : HistOp1D(nBinsX=25, x0=0.0,  x1=250.0, xUnits='GeV', xLabel='LFV Collinear Mass m_{coll}',       yLabel='Events'),
    'MtLep0'          : HistOp1D(nBinsX=15, x0=0.0,  x1=250.0, xUnits='GeV', xLabel='m_{T}(l_{0},MET)',                  yLabel='Events'),
    'MtLep1'          : HistOp1D(nBinsX=20, x0=0.0,  x1=140.0, xUnits='GeV', xLabel='m_{T}(l_{1},MET)',                  yLabel='Events'),
    'MLL'             : HistOp1D(nBinsX=25, x0=0.0,  x1=400.0, xUnits='GeV', xLabel='M_{ll}',                            yLabel='Events'),
    'ptll'            : HistOp1D(nBinsX=20, x0=0.0,  x1=200.0, xUnits='GeV', xLabel='pT_{ll}',                           yLabel='Events'),
    'dpt_ll'          : HistOp1D(nBinsX=20, x0=0.0,  x1=150.0, xUnits='GeV', xLabel='#Deltap_{T}^{ll}',                  yLabel='Events'),
    'DphiLep0MET'     : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='#Delta#phi(l_{0},MET)',             yLabel='Events'),
    'DphiLep1MET'     : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='#Delta#phi(l_{1},MET)',             yLabel='Events'),
    'JetN'            : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{base jet}',                      yLabel='Events'),
    'Jet_N2p4Eta25Pt' : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{jet} (p_{T}>25GeV, |#eta|<2.5)', yLabel='Events'),
    'SignalJetN'      : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{sig jets}',                      yLabel='Events'),
    'JetN_g30'        : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{jet} (p_{T}>30GeV)',             yLabel='Events'),
    'nCentralLJets'   : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{CL jet}',                        yLabel='Events'),
    'nCentralBJets'   : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{CB jet}',                        yLabel='Events'),
    'Btag'            : HistOp1D(nBinsX=5,  x0=-1.5, x1=3.5,   xUnits='',    xLabel='B-tagged jet',                      yLabel='Events'),
    'nForwardJets'    : HistOp1D(nBinsX=8,  x0=-0.5, x1=7.5,   xUnits='',    xLabel='N_{F jet}',                         yLabel='Events'),
    'j_pt[0]'         : HistOp1D(nBinsX=20, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='p_{T}^{leading jet}',               yLabel='Events'),
    'j_pt'            : HistOp1D(nBinsX=20, x0=0.0,  x1=500.0, xUnits='GeV', xLabel='Jet p_{T}',                         yLabel='Events'),
    'Mjj'             : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='GeV', xLabel='Dijet mass',                        yLabel='Events'),
    'DEtaJJ'          : HistOp1D(nBinsX=20, x0=0.0,  x1=6.0,   xUnits='',    xLabel='#Delta#eta(j0,j1)',                 yLabel='Events'),
    'j_eta'           : HistOp1D(nBinsX=20, x0=-3.0, x1=3.0,   xUnits='',    xLabel='Jet #eta',                          yLabel='Events'),
    'j_jvt'           : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Jet JVT',                           yLabel='Events'),
    'j_jvf'           : HistOp1D(nBinsX=25, x0=0.0,  x1=-1,    xUnits='',    xLabel='Jet JVF',                           yLabel='Events'),
    'j_phi'           : HistOp1D(nBinsX=30, x0=0.0,  x1=3.15,  xUnits='',    xLabel='Jet #phi',                          yLabel='Events'),
    'j_flav'          : HistOp1D(nBinsX=5,  x0=-0.5, x1=4.5,   xUnits='',    xLabel='Jet flavor (0:NA,1:CL,2:CB,3:F)',   yLabel='Events'),
    'pass_HLT_mu18_mu8noL1'                 : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='mu18 mu8noL1 trigger bit',                 yLabel='Events'),
    'pass_HLT_2e12_lhloose_L12EM10VH'       : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='2e12 lhloose L12EM10VH trigger bit',       yLabel='Events'),
    'pass_HLT_e17_lhloose_mu14'             : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='e17 lhloose mu14 trigger bit',             yLabel='Events'),
    'pass_HLT_mu20_mu8noL1'                 : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='mu20 mu8noL1 trigger bit',                 yLabel='Events'),
    'pass_HLT_2e15_lhvloose_L12EM13VH'      : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='2e15 lhvloose L12EM13VH trigger bit',      yLabel='Events'),
    'pass_HLT_2e17_lhvloose_nod0'           : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='2e17 lhvloose nod0 trigger bit',           yLabel='Events'),
    'pass_HLT_mu22_mu8noL1'                 : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='mu22 mu8noL1 trigger bit',                 yLabel='Events'),
    'pass_HLT_e17_lhloose_nod0_mu14'        : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='e17 lhloose nod0 mu14 trigger bit',        yLabel='Events'),
    'pass_HLT_e60_lhmedium'                 : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='e60 lhmedium trigger bit',                 yLabel='Events'),
    'pass_HLT_e24_lhmedium_L1EM20VH'        : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='e24 lhmedium L1EM20VH trigger bit',        yLabel='Events'),
    'pass_HLT_e24_lhmedium_iloose_L1EM18VH' : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='e24 lhmedium iloose L1EM18VH trigger bit', yLabel='Events'),
    'pass_HLT_mu20_iloose_L1MU15'           : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='mu20 iloose L1MU15 trigger bit',           yLabel='Events'),
    'pass_HLT_mu24_imedium'                 : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='mu24 imedium trigger bit',                 yLabel='Events'),
    'pass_HLT_mu26_imedium'                 : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='mu26 imedium trigger bit',                 yLabel='Events'),
    'pass_HLT_e24_lhtight_nod0_ivarloose'   : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='e24 lhtight nod0 ivarloose trigger bit',   yLabel='Events'),
    'pass_HLT_e26_lhtight_nod0_ivarloose'   : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='e26 lhtight nod0 ivarloose trigger bit',   yLabel='Events'),
    'pass_HLT_e60_lhmedium_nod0'            : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='e60 lhmedium nod0 trigger bit',            yLabel='Events'),
    'pass_HLT_e120_lhloose'                 : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='e120 lhloose trigger bit',                 yLabel='Events'),
    'pass_HLT_e140_lhloose_nod0'            : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='e140 lhloose nod0 trigger bit',            yLabel='Events'),
    'pass_HLT_mu24_iloose'                  : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='mu24 iloose trigger bit',                  yLabel='Events'),
    'pass_HLT_mu24_iloose_L1MU15'           : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='mu24 iloose L1MU15 trigger bit',           yLabel='Events'),
    'pass_HLT_mu24_ivarloose'               : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='mu24 ivarloose trigger bit',               yLabel='Events'),
    'pass_HLT_mu24_ivarloose_L1MU15'        : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='mu24 ivarloose L1MU15 trigger bit',        yLabel='Events'),
    'pass_HLT_mu24_ivarmedium'              : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='mu24 ivarmedium trigger bit',              yLabel='Events'),
    'pass_HLT_mu26_ivarmedium'              : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='mu26 ivarmedium trigger bit',              yLabel='Events'),
    'pass_HLT_mu50'                         : HistOp1D(nBinsX=5, x0=-1.5, x1=3.5, xUnits='', xLabel='mu50 trigger bit',                         yLabel='Events')
}

# Samples recommended by Duc Bao Ta (xiedebao@googlemail.com)
# subset of full sample list used in Htt analysis (ATL-COM-PHYS-2017-446)
# additions where indicated
groups = {}
groups['higgs_lfv'] = [
    344084, 344085, 344086, 344087, 344088, 344089, 344090, 344091,
    345077, 345078, 
    345124, 345125, 
    345213, 345214, 345215, 345216,
    345218, 345219
]
groups['htt'] = [
    341122, 
    341155,
    342170, 342171, 342172,
    345073,
    345120,
    #343323,
    #343334,
    #343942,
    #343953,
]
groups['hww'] = [
    341079, 341080,
]
groups['ttbar'] = [
    410009, 
]
groups['singletop'] = [
    410011, 410012, 
    410015, 410016, 
    410025, 410026,
]
groups['diboson'] = [
    361600, 361601, 
    361603, 361604, 
    361607, 
    361610,
]
groups['zll'] = [
    364100, 364101, 364102, 364103, 364104, 364105, 364106, 364107, 364108,
    364109, 364110, 364111, 364112, 364113,
    364114, 364115, 364116, 364117, 364118, 364119, 364120, 364121, 364122,
    364123, 364124, 364125, 364126, 364127,
    364198, 364199, 364200, 364201, 364202, 364203,
    364204, 364205, 364206, 364207, 364208, 364209,
    #304018, 304019, 304020, replaced with 344441-3
    344441, 344442, 344443,
    345099, 345100, 345101, 345102,
]
groups['ztt'] = [
    344772, 344776, 344780,
    #304021,
    #364128, 364129, 364130, 364131, 364132, 364133, 364134, 364135, 364136,
    364137, 364138, 364139, 364140, 364141, 364210, 364211, 364212, 364213,
    364214, 364215,
] 
groups['wjets'] = [
    # Temporary fake sample
    361100, 361101, 361102, 361103, 361104, 361105,
]
groups['data15'] = [
    276262, 276329, 276336, 276416, 276511, 276689, 276778, 276790, 276952, 
    276954, 278880, 278912, 278968, 279169, 279259, 279279, 279284, 279345, 
    279515, 279598, 279685, 279813, 279867, 279928, 279932, 279984, 280231, 
    280273, 280319, 280368, 280423, 280464, 280500, 280520, 280614, 280673, 
    280753, 280853, 280862, 280950, 280977, 281070, 281074, 281075, 281317, 
    281385, 281411, 282625, 282631, 282712, 282784, 282992, 283074, 283155, 
    283270, 283429, 283608, 283780, 284006, 284154, 284213, 284285, 284420, 
    284427, 284484,
]
groups['data16'] = [
    297730, 298595, 298609, 298633, 298687, 298690, 298771, 298773, 298862, 
    298967, 299055, 299144, 299147, 299184, 299243, 299584, 300279, 300345, 
    300415, 300418, 300487, 300540, 300571, 300600, 300655, 300687, 300784, 
    300800, 300863, 300908, 301912, 301918, 301932, 301973, 302053, 302137, 
    302265, 302269, 302300, 302347, 302380, 302391, 302393, 302737, 302831, 
    302872, 302919, 302925, 302956, 303007, 303079, 303201, 303208, 303264, 
    303266, 303291, 303304, 303338, 303421, 303499, 303560, 303638, 303832, 
    303846, 303892, 303943, 304006, 304008, 304128, 304178, 304198, 304211, 
    304243, 304308, 304337, 304409, 304431, 304494, 305380, 305543, 305571, 
    305618, 305671, 305674, 305723, 305727, 305735, 305777, 305811, 305920, 
    306269, 306278, 306310, 306384, 306419, 306442, 306448, 306451, 307126, 
    307195, 307259, 307306, 307354, 307358, 307394, 307454, 307514, 307539, 
    307569, 307601, 307619, 307656, 307710, 307716, 307732, 307861, 307935, 
    308047, 308084, 309375, 309390, 309440, 309516, 309640, 309674, 309759, 
    310015, 310247, 310249, 310341, 310370, 310405, 310468, 310473, 310634, 
    310691, 310738, 310809, 310863, 310872, 310969, 311071, 311170, 311244, 
    311287, 311321, 311365, 311402, 311473, 311481,
]
def get_mc_groups():
    mc_groups = {}
    for key, lst in groups.iteritems():
        if 'data' in key.lower():
            continue
        mc_groups[key] = lst
    return mc_groups
def get_bkgd_groups():
    bkgd_groups = {}
    for key, lst in groups.iteritems():
        if 'data' in key.lower() or 'signal' in key.lower():
            continue
        bkgd_groups[key] = lst
    return bkgd_groups
def get_data_groups():
    data_groups = {}
    for key, lst in groups.iteritems():
        if 'data' not in key.lower():
            continue
        data_groups[key] = lst
    return data_groups

# Set ATLAS style
def setATLASStyle(path="/data/uclhc/uci/user/armstro1/ATLASStyle/"):
    ROOT.gROOT.SetMacroPath(path)
    ROOT.gROOT.LoadMacro("AtlasStyle.C")
    ROOT.gROOT.LoadMacro("AtlasLabels.C")
    ROOT.gROOT.LoadMacro("AtlasUtils.C")
    ROOT.SetAtlasStyle()

# Dummy ratio histogram
def dummifyHistogram(histo):
    ratio_histo = histo.Clone();
    ratio_histo.Reset();
    ratio_histo.SetMarkerSize(1.2);
    ratio_histo.SetMarkerStyle(20);
    ratio_histo.SetLineWidth(2);
    ratio_histo.GetYaxis().SetTitle("#frac{Data}{SM}");
    ratio_histo.GetXaxis().SetLabelSize(0.1);
    ratio_histo.GetXaxis().SetLabelOffset(0.02);
    ratio_histo.GetXaxis().SetTitleSize(0.12);
    ratio_histo.GetXaxis().SetTitleOffset(1.);
    ratio_histo.GetYaxis().SetRangeUser(0.001,2);
    ratio_histo.GetYaxis().SetLabelSize(0.1);
    ratio_histo.GetYaxis().SetTitleSize(0.12);
    ratio_histo.GetYaxis().SetTitleOffset(0.5);
    ratio_histo.GetYaxis().SetNdivisions(5);
    return ratio_histo
