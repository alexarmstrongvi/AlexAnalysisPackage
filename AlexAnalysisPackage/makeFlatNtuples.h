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
using sflow::Superflow;
using sflow::SuperflowRunMode;
using sflow::Superlink;

void setup_chain(TChain* chain, string iname);
void initialize_fake_factor_tool(ApplyFakeFactor* applyFakeFactorTool);
Superflow* get_cutflow(TChain* chain, Sel sel_type);
string determine_suffix(string user_suffix, Sel sel_type, bool apply_ff);
Superflow* initialize_superflow(TChain *chain, string name_suffix);

void add_shortcut_variables(Superflow* superflow, Sel sel_type);
void add_pre_cuts(Superflow* superflow, Sel sel_type);
void add_cleaing_cuts(Superflow* superflow);
void add_baseline_lepton_cuts(Superflow* superflow);
void add_baseline_den_lepton_cuts(Superflow* superflow);
void add_fake_num_lepton_cuts(Superflow* superflow);
void add_fake_den_lepton_cuts(Superflow* superflow);
void add_final_cuts(Superflow* superflow, Sel sel_type);
void add_event_variables(Superflow* superflow);
void add_trigger_variables(Superflow* superflow);
void add_prelepton_variables(Superflow* superflow);
void add_baselepton_variables(Superflow* superflow);
void add_signallepton_variables(Superflow* superflow);
void add_xtau_lepton_variables(Superflow* superflow);
void add_tau_variables(Superflow* superflow);
void add_other_lepton_variables(Superflow* superflow);
void add_met_variables(Superflow* superflow);
void add_jet_variables(Superflow* superflow);
void add_fake_variables(Superflow* superflow);
void add_shortcut_variables_reset(Superflow* superflow);


void add_fake_shortcut_variables(Superlink* sl);
bool is_ID_lepton(Superlink* sl, Susy::Lepton* lepton);
bool is_antiID_lepton(Susy::Lepton* lepton);
bool is_1lep_trig_matched(Superlink* sl, string trig_name, LeptonVector leptons);
void add_SFOS_lepton_cut(Superflow* superflow);
void add_DFOS_lepton_cut(Superflow* superflow);
int get_lepton_truth_class(Susy::Lepton* lepton);

// Important globals
TChain* m_chain = new TChain("susyNt");
ApplyFakeFactor* m_applyFakeFactorTool = new ApplyFakeFactor("FF_tool");

bool m_denominator_selection = false;

// All globals must be initialized here and reset in the source file
int m_cutflags = 0;  ///< Cutflags used for cleaning cuts

TLorentzVector m_MET;

uint m_lepID_n = 0;
uint m_lepAntiID_n = 0;
Susy::Lepton* m_antiID_lep0 = nullptr;
Susy::Lepton* m_antiID_lep1 = nullptr;
TLorentzVector m_antiID_lep0_TLV;
TLorentzVector m_antiID_lep1_TLV;

// Find leptons paired with Z
// Indices 0 and 1 are closest Z pair
// Index 2 is leading anti-ID or 3rd leading ID lepton
// Indices 3 and 4 are second closest Z pair
LeptonVector m_Zlep(5, nullptr); //TODO: change name to m_Zordered_leptons
LeptonVector m_selectLeptons;
LeptonVector m_triggerLeptons;

TLorentzVector m_lepton0;
TLorentzVector m_lepton1;
TLorentzVector m_dileptonP4;
Susy::Electron* m_el0 = nullptr;
Susy::Electron* m_el1 = nullptr;
Susy::Muon* m_mu0 = nullptr;
Susy::Muon* m_mu1 = nullptr;

JetVector m_lightJets;
JetVector m_BJets;
JetVector m_forwardJets;
TLorentzVector m_Dijet_TLV, m_Jet1_TLV, m_Jet0_TLV;

////////////////////////////////////////////////////////////////////////////////
// Configuration settings
////////////////////////////////////////////////////////////////////////////////

SuperflowRunMode m_run_mode = SuperflowRunMode::nominal;  ///< The mode in which
            ///< Superflow will process events (e.g. compute weights or not)

////////////////////////////////////////////////////////////////////////////////
// Tool for applying the fake factor

// set path to where the fake factor histograms are located
string m_path_to_FF_file = "/data/uclhc/uci/user/armstro1/SusyNt/SusyNt_n0235_LFV_analysis/analysis/AlexAnalysisPackage/plotting/FF_hists/";

//TODO: Get a better toggle setup than just commenting out lines
//// My own fake factors
//fake_file = "FF_hists.root";
//m_el_FF_hist = "h_zjets_FF_CR_e_data_Z_Lep2_pT_minus_bkgd";
//m_mu_FF_hist = "h_zjets_FF_CR_m_data_Z_Lep2_pT_minus_bkgd";
//// Fake Factors from another analysis
//fake_file = "FakeFactors_April12_2017_2.4.29_IncludingJVTSFs.root";
//m_el_FF_hist = "FakeFactor_el";
//m_mu_FF_hist = "FakeFactor_mu";
// Fake factors hists that always return 1
string fake_file = "dummy_FF_hists.root";
string m_el_FF_hist = "h_zjets_FF_CR_e_data_Z_Lep2_pT_minus_bkgd";
string m_mu_FF_hist = "h_zjets_FF_CR_m_data_Z_Lep2_pT_minus_bkgd";
//}
string m_fake_path = m_path_to_FF_file + fake_file;




