////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) <2018> by Alex Armstrong
///
/// @file makeFlatNtuple.h
/// @author Alex Armstrong <alarmstr@cern.ch>
/// @date <May 27th, 2018>
///
////////////////////////////////////////////////////////////////////////////////

#pragma once

// general cpp
#include <iostream>
#include <string>
#include <cmath>
#include <limits.h>

// analysis include(s)
#include "Superflow/Superflow.h"
#include "Superflow/Superlink.h"

#include "SusyNtuple/ChainHelper.h"
#include "SusyNtuple/string_utils.h"
#include "SusyNtuple/KinematicTools.h"
#include "SusyNtuple/MCTruthClassifierDefs.h"
//#include "SusyNtuple/SusyNt.h"
//#include "SusyNtuple/SusyDefs.h"
//#include "SusyNtuple/SusyNtObject.h"
//#include "SusyNtuple/SusyNtTools.h"

#include "AlexAnalysisPackage/ApplyFakeFactor.h"

/// @brief Available run selections
enum Sel {FAKE_NUM, FAKE_DEN, BASELINE, BASE_DEN};

using std::string;

TChain* setup_chain();
ApplyFakeFactor* initialize_fake_factor_tool();
Superflow* get_cutflow(TChain* chain, Sel sel_type);

//
Superflow* initialize_superflow(TChain *chain);
void add_xtau_cuts(Superflow* superflow, Sel sel_type);
void add_cleaing_cuts(Superflow* superflow, Sel sel_type);
void add_fake_var_defs(Superflow* superflow, Sel sel_type);
void add_lepton_cuts(Superflow* superflow, Sel sel_type);
void add_final_cuts(Superflow* superflow, Sel sel_type);
void add_shortcut_variables(Superflow* superflow, Sel sel_type);
void add_event_variables(Superflow* superflow, Sel sel_type);
void add_trigger_variables(Superflow* superflow, Sel sel_type);
void add_lepton_variables(Superflow* superflow, Sel sel_type);
void add_met_variables(Superflow* superflow, Sel sel_type);
void add_jet_variables(Superflow* superflow, Sel sel_type);
void add_fake_var_defs(Superflow* superflow, Sel sel_type);
void add_shortcut_variables_reset(Superflow* superflow, Sel sel_type);

//

// Important globals
TChain* m_chain = new TChain("susyNt");
ApplyFakeFactor* m_applyFakeFactorTool = new ApplyFakeFactor("FF_tool");

bool m_do_fakes = false;

// TODO: Final check on globals
// Search by removing m_ and see if anything turns up in search
// Search for m_ that doesn't belong to global

// All globals must be initialized here and reset in the source file
int m_cutflags = 0;  ///< Cutflags used for cleaning cuts

LeptonVector m_preLeptons;
LeptonVector m_baseLeptons;
LeptonVector m_signalLeptons;
ElectronVector m_preElectrons;
ElectronVector m_baseElectrons;
ElectronVector m_signalElectrons;
MuonVector m_preMuons;
MuonVector m_baseMuons;
MuonVector m_signalMuons;
TauVector m_preTaus;
TauVector m_baseTaus;
TauVector m_signalTaus;
Susy::Met  m_met;

uint m_lepID_n = 0,
uint m_lepAntiID_n = 0;
int m_antiID_idx0 = -1;  // TODO: Change to a Susy::Lepton m_antiID_lep0
int m_antiID_idx1 = -1;  // TODO: Change to a Susy::Lepton m_antiID_lep1
TLorentzVector m_LepFake0;
TLorentzVector m_LepFake1;

// Find leptons paired with Z
// Indices 0 and 1 are closest Z pair
// Index 2 is leading anti-ID or 3rd leading ID lepton
// Indices 3 and 4 are second closest Z pair
LeptonVector m_Zlep(5, nullptr); //TODO: change name to m_Zordered_leptons
LeptonVector m_selectLeptons;

TLorentzVector m_lepton0;
TLorentzVector m_lepton1;
TLorentzVector m_dileptonP4;
Susy::Electron* m_el0 = nullptr;
Susy::Electron* m_el1 = nullptr;
Susy::Muon* m_mu0 = nullptr;
Susy::Muon* m_mu1 = nullptr;

JetVector m_preJets;
JetVector m_baseJets;
JetVector m_signalJets;
JetVector m_lightJets;
JetVector m_BJets;
JetVector m_forwardJets;
TLorentzVector m_JetP4, m_Jet1, m_Jet0;

////////////////////////////////////////////////////////////////////////////////
// Configuration settings
////////////////////////////////////////////////////////////////////////////////

SuperflowRunMode m_run_mode = SuperflowRunMode::nominal; ///< The mode in which
            ///< Superflow will process events (e.g. compute weights or not)

////////////////////////////////////////////////////////////////////////////////
// Tool for applying the fake factor

// set path to where the fake factor histograms are located
string m_path_to_FF_file = "/data/uclhc/uci/user/armstro1/SusyNt/"
                            "SusyNt_n0235_LFV_analysis/analysis/"
                            "AlexAnalysisPackage/plotting/FF_hists/";
string fake_file;
string m_el_FF_hist, m_mu_FF_hist;
// Toggle which FF hists to use
if (false) {
    // My own fake factors
    fake_file = "FF_hists.root";
    m_el_FF_hist = "h_zjets_FF_CR_e_data_Z_Lep2_pT_minus_bkgd";
    m_mu_FF_hist = "h_zjets_FF_CR_m_data_Z_Lep2_pT_minus_bkgd";
} else if (false) {
    // Fake Factors from another analysis
    fake_file = "FakeFactors_April12_2017_2.4.29_IncludingJVTSFs.root";
    m_el_FF_hist = "FakeFactor_el";
    m_mu_FF_hist = "FakeFactor_mu";
} else {
    // Fake factors hists that always return 1
    fake_file = "dummy_FF_hists.root";
    m_el_FF_hist = "h_zjets_FF_CR_e_data_Z_Lep2_pT_minus_bkgd";
    m_mu_FF_hist = "h_zjets_FF_CR_m_data_Z_Lep2_pT_minus_bkgd";
}
m_fake_path = m_path_to_FF_file + m_fake_file;




