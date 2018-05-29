////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) <2018> by Alex Armstrong with much code borrowed from
///     Aleattin Mete (Alaettin.Serhan.Mete@cern.ch)
///
/// @file makeFlatNtuple.cxx
/// @author Alex Armstrong <alarmstr@cern.ch>
/// @date <May 27th, 2018>
/// @brief Make flat ntuples from SusyNts
///
////////////////////////////////////////////////////////////////////////////////
// TODO:
// - See the effect of adding variables on ntuple production time
// - Fix dilep trigger matching
// - Improve all the variable nameing schemes

#include "AlexAnalysisPackage/makeFlatNtuples.h"

using std::string;
using std::cout;
using std::cerr;
using std::vector;

// Various superflow classes
using namespace sflow;

////////////////////////////////////////////////////////////////////////////////
// Declarations
////////////////////////////////////////////////////////////////////////////////
// Useful macros
// TODO: Run over a vector of leptons
#define ADD_1LEP_TRIGGER_VAR(trig_name, leptons) { \
    *superflow << NewVar(#trig_name" trigger bit"); { \
        *superflow << HFTname(#trig_name); \
        *superflow << [=](Superlink* sl, var_bool*) -> bool { \
            return is_1lep_trig_matched(sl, #trig_name, leptons); \
        }; \
        *superflow << SaveVar(); \
    } \
}

//#define ADD_1LEP_TRIGGER_VAR(trig_name, lep) { \
//    *superflow << NewVar(#trig_name" trigger bit"); { \
//        *superflow << HFTname(#trig_name); \
//        *superflow << [=](Superlink* sl, var_bool*) -> bool { \
//            if(!lep) return false;\
//            bool trig_fired = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, #trig_name); \
//            if (!trig_fired) return false;\
//            bool trig_matched = sl->tools->triggerTool().lepton_trigger_match(lep, #trig_name);\
//            return trig_matched; }; \
//        *superflow << SaveVar(); \
//    } \
//}
// Trig Matching for dilep triggers is buggy
// so currently not trigger matching
#define ADD_2LEP_TRIGGER_VAR(trig_name, lep0, lep1) { \
  *superflow << NewVar(#trig_name" trigger bit"); { \
      *superflow << HFTname(#trig_name); \
      *superflow << [=](Superlink* sl, var_bool*) -> bool { \
          if(!lep0 || ! lep1) return false;\
          bool trig_fired = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, #trig_name); \
          return trig_fired;\
          if (!trig_fired) return false; \
          bool trig_matched = sl->tools->triggerTool().dilepton_trigger_match(sl->nt->evt(), lep0, lep1, #trig_name);\
          return trig_matched; }; \
      *superflow << SaveVar(); \
  } \
}
#define ADD_MULTIPLICITY_VAR(container) { \
  *superflow << NewVar("number of "#container); { \
    *superflow << HFTname("n_"#container); \
    *superflow << [=](Superlink* sl, var_int*) -> int { \
      return sl->container->size(); \
    }; \
    *superflow << SaveVar(); \
  } \
}

/// @brief Command line argument parser
struct Args {
    ////////////////////////////////////////////////////////////////////////////
    // Initialize arguments
    ////////////////////////////////////////////////////////////////////////////
    string PROG_NAME;  // always required

    // Required arguments
    string input_name;  ///< Input file name

    // Optional arguments with defaults
    unsigned int n_events  = -1; ///< number of events processed
    unsigned int n_skipped = 0; ///< number of initila entries to skip
    string name_suffix = ""; ///< suffix appended to output name
    bool baseline_sel = false;  ///< Apply baseline denominator selection
    bool baseline_den = false;  ///< Apply baseline denominator selection
    bool fake_num = false;  ///< Apply fake numerator selection
    bool fake_den = false;  ///< Apply fake denominator selection
    bool apply_ff = false; ///< Apply fake factor event weight to denominator events

    ////////////////////////////////////////////////////////////////////////////
    // Print information
    ////////////////////////////////////////////////////////////////////////////
    void print_usage() const {
        printf("===========================================================\n");
        printf(" %s\n", PROG_NAME.c_str());
        printf(" Makeing analysis flat ntuples from SusyNts\n");
        printf("===========================================================\n");
        printf("Required Parameters:\n");
        printf("\t-i, --input     Input file name\n");
        printf("\nOptional Parameters:\n");
        printf("\t-n, --nevents   number of events to process\n");
        printf("\t-k, --nskipped  number of initial entries to skip\n");
        printf("\t-s, --suffix    suffix appended to output name\n");
        printf("\t--baseline_sel  use baseline event selection [default selection]\n");
        printf("\t--baseline_den  use baseline denominator event selection \n");
        printf("\t--fake_num      use Z+jets numerator event selection \n");
        printf("\t--fake_den      use Z+jets denominator event selection \n");
        printf("\t--apply_ff      apply fake factor to denominator events\n");
        printf("\t-h, --help      print this help message\n");
        printf("===========================================================\n");
    }

    void print() const {
        printf("===========================================================\n");
        printf(" %s Configuration\n", PROG_NAME.c_str());
        printf("===========================================================\n");
        printf("\tInput file name : %s\n", input_name.c_str());
        printf("\tnevents         : %i\n", n_events);
        printf("\tnskipped        : %i\n", n_skipped);
        printf("\tsuffix          : %s\n", name_suffix.c_str());
        printf("\tbaseline_sel    : %i\n", baseline_sel);
        printf("\tbaseline_den    : %i\n", baseline_den);
        printf("\tfake_num        : %i\n", fake_num);
        printf("\tfake_den        : %i\n", fake_den);
        printf("\tapply_ff        : %i\n", apply_ff);
        printf("===========================================================\n");
    }

    ////////////////////////////////////////////////////////////////////////////
    // Parser
    ////////////////////////////////////////////////////////////////////////////
    bool parse(int argc, char* argv[]) {
        PROG_NAME = argv[0];

        // Parse arguments
        for (int i = 0; i< argc; ++i) {
            // Grab arguments
            string arg = argv[i];
            string arg_value = argc > i+1 ? argv[i+1] : "";
            // Skip if arg set to arg value and not arg name
            if (arg.at(0) != '-') continue;

            // Check for required arguments
            if (arg == "-i" || arg == "--input") {
                input_name = arg_value;
            } else if (arg == "-n" || arg == "--nevents") {
                n_events = atoi(arg_value.c_str());
            } else if (arg == "-k" || arg == "--nskipped") {
                n_skipped = atoi(arg_value.c_str());
            } else if (arg == "-s" || arg == "--suffix") {
                name_suffix = true;
            } else if (arg == "--baseline_sel") {
                baseline_sel = true;
            } else if (arg == "--baseline_den") {
                baseline_den = true;
            } else if (arg == "--fake_num") {
                fake_num = true;
            } else if (arg == "--fake_den") {
                fake_den = true;
            } else if (arg == "--apply_ff") {
                apply_ff = true;
            } else if (arg == "-h" || arg == "--help") {
                print_usage();
                return false;
            } else {
                cerr << "ERROR :: Unrecognized input argument: "
                     << arg << " -> " << arg_value << '\n';
                print_usage();
            }
        }

        // Check arguments
        if (input_name.size() == 0) {
            cerr << "ERROR :: No input source given\n";
            return false;
        } else if (apply_ff && !(baseline_den || fake_den)) {
            cerr << "ERROR :: trying to apply fake factor without using any "
                 << "denominator selection options\n";
            return false;
        }
        return true;
    }
} args;

////////////////////////////////////////////////////////////////////////////////
/// @brief Main function
///
/// Run with help option (-h, --help) to see available parameters
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
    cout << "\n====== RUNNING " << argv[0] << " ====== \n";
    if (args.parse(argc, argv)) {
        // Finished parsing arguments
        args.print();
    } else {
        // Failed to parse arguments or help requested
        return 1;
    }
    m_denominator_selection = args.fake_den || args.baseline_den;

    ////////////////////////////////////////////////////////////////////////////
    // Main implementation
    ////////////////////////////////////////////////////////////////////////////
    // Build list of cutflows to run
    setup_chain(m_chain, args.input_name);
    printf("makeFlatNtuple :: Total events available : %lli\n",m_chain->GetEntries());

    if (args.apply_ff) initialize_fake_factor_tool(m_applyFakeFactorTool);

    // Run the chosen cutflows
    if (args.baseline_sel) {
        cout << "\n\n Creating baseline cutflow \n\n";
        Superflow* sf = get_cutflow(m_chain, BASELINE);
        m_chain->Process(sf, args.input_name.c_str(), args.n_events, args.n_skipped);
        delete sf; sf = 0;
    }
    if (args.baseline_den) {
        cout << "\n\n Creating baseline denominator cutflow \n\n";
        Superflow* sf = get_cutflow(m_chain, BASE_DEN);
        m_chain->Process(sf, args.input_name.c_str(), args.n_events, args.n_skipped);
        delete sf; sf = 0;
    }
    if (args.fake_den) {
        cout << "\n\n Creating Z+jets denominator cutflow \n\n";
        Superflow* sf = get_cutflow(m_chain, FAKE_DEN);
        m_chain->Process(sf, args.input_name.c_str(), args.n_events, args.n_skipped);
        delete sf; sf = 0;
    }
    if (args.fake_num) {
        cout << "\n\n Creating Z+jets numerator cutflow \n\n";
        Superflow* sf = get_cutflow(m_chain, FAKE_NUM);
        m_chain->Process(sf, args.input_name.c_str(), args.n_events, args.n_skipped);
        delete sf; sf = 0;
    }
    
    //TODO: ttbar, zll, W+jets, Diboson, Ztautau Regions

    delete m_chain;
    cout << "\n====== SUCCESSFULLY RAN " << argv[0] << " ====== \n";
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
////////////////////////////////////////////////////////////////////////////////
void setup_chain(TChain* chain, string iname) {
    chain->SetDirectory(0);  // Remove ROOT ownership

    bool inputIsFile = Susy::utils::endswith(iname, ".root");
    bool inputIsList = Susy::utils::endswith(iname, ".txt");
    bool inputIsDir  = Susy::utils::endswith(iname, "/");

    if (inputIsFile) {
        ChainHelper::addFile(chain, iname);
    } else if (inputIsList) {
      // If a list of ROOT files
        ChainHelper::addFileList(chain, iname);
    } else if (inputIsDir) {
        ChainHelper::addFileDir(chain, iname);
    } else {
        printf("ERROR (initialize_chain) :: Unrecognized input %s", iname.c_str());
    }
}

void initialize_fake_factor_tool(ApplyFakeFactor* applyFakeFactorTool) {
    applyFakeFactorTool->set_savedFakeFactorFileName(m_fake_path);
    applyFakeFactorTool->set_saved_elFakeFactorName(m_el_FF_hist);
    applyFakeFactorTool->set_saved_muFakeFactorName(m_mu_FF_hist);
    applyFakeFactorTool->initialize().ignore();
}

Superflow* get_cutflow(TChain* chain, Sel sel_type) {
    args.name_suffix = determine_suffix(args.name_suffix, sel_type, args.apply_ff);

    Superflow* superflow = initialize_superflow(chain, args.name_suffix);

    // How to add cutflow entry:
    // Create lambda function that returns bool value of cut.
    // Pass that with "<<" into the function CutName("Cut name").
    // Pass that with "<<" into the dereferenced cutflow object.

    // *IMPORTANT* The order that cuts are added is very important as that is
    // order in which cuts are applied and global superflow variables filled

    ////////////////////////////////////////////////////////////////////////////
    // Define helpful variables for use inside superflow variable definitions
    add_shortcut_variables(superflow, sel_type);

    ////////////////////////////////////////////////////////////////////////////
    // Add cuts
    add_pre_cuts(superflow, sel_type);
    add_cleaing_cuts(superflow);
    if (sel_type == BASELINE) {
        add_baseline_lepton_cuts(superflow);
    } else if (sel_type == BASE_DEN) {
        add_baseline_den_lepton_cuts(superflow);
    } else if (sel_type == FAKE_NUM) {
        add_fake_num_lepton_cuts(superflow);
    } else if (sel_type == FAKE_DEN) {
        add_fake_den_lepton_cuts(superflow);
    }

    add_final_cuts(superflow, sel_type);

    ////////////////////////////////////////////////////////////////////////////
    // Add superflow variables
    add_event_variables(superflow);
    add_trigger_variables(superflow);
    add_met_variables(superflow);
    add_prelepton_variables(superflow);
    add_baselepton_variables(superflow);
    add_signallepton_variables(superflow);
    add_xtau_lepton_variables(superflow);
    add_tau_variables(superflow);
    add_other_lepton_variables(superflow);
    add_jet_variables(superflow);
    add_fake_variables(superflow);

    ////////////////////////////////////////////////////////////////////////////
    // Reset helpful variables
    add_shortcut_variables_reset(superflow);

    return superflow;
}

// TODO: Switch to using references instead of pointers
string determine_suffix(string user_suffix, Sel sel_type, bool apply_ff) {
    string full_suffix = user_suffix;

    if (user_suffix != "") full_suffix += "_";

    if (apply_ff) full_suffix += "fakes_";

    if (sel_type == FAKE_NUM) full_suffix += "zjets_num";
    else if (sel_type == FAKE_DEN) full_suffix += "zjets_den";
    else if (sel_type == BASELINE) full_suffix += "baseline";
    else if (sel_type == BASE_DEN) full_suffix += "baseline_den";

    return full_suffix;
}
Superflow* initialize_superflow(TChain *chain, string name_suffix) {
    // Move run_mode to globals

    Superflow* superflow = new Superflow();        // initialize the cutflow
    superflow->setAnaName("SuperflowAna");         // arbitrary
    superflow->setAnaType(AnalysisType::Ana_HLFV); // analysis type, passed to SusyNt
    superflow->setLumi(1.0);                       // set the MC normalized to X pb-1
    superflow->setSampleName(args.input_name);     // sample name, check to make sure it's set OK
    superflow->setRunMode(m_run_mode);             // make configurable via run_mode
    superflow->setCountWeights(true);              // print the weighted cutflows
    if(name_suffix != "") superflow->setFileSuffix(name_suffix);
    superflow->setChain(chain);
    superflow->nttools().initTriggerTool(ChainHelper::firstFile(args.input_name, 0.));

    return superflow;
}

////////////////////////////////////////////////////////////////////////////////
// Add lambda expression to superflow
void add_shortcut_variables(Superflow* superflow, Sel sel_type) {
    *superflow << CutName("read in") << [=](Superlink* sl) -> bool {
        m_cutflags = sl->nt->evt()->cutFlags[NtSys::NOM];

        ////////////////////////////////////////////////////////////////////////
        // Fake shortcuts
        add_fake_shortcut_variables(sl);

        // Pick the right set and ordingering of leptons for the given selection
        // type. Affects trigger matching, dilepton requirements, etc...
        //
        // The first two leptons in the container are assumed to be those that
        // should be used in determining dilepton properties (e.g. MLL, dphi)
        // and requirements (e.g. SFOS, DFOS). The third is an additional,
        // usually probe lepton.
        // TODO: Add grabbing of lepton to its own function for each selection type
        if (sel_type == FAKE_DEN || sel_type == FAKE_NUM) {
            m_selectLeptons = m_Zlep;
            m_triggerLeptons.push_back(m_Zlep.at(0));
            m_triggerLeptons.push_back(m_Zlep.at(1));
        } else if (sel_type == BASELINE ) {
            m_selectLeptons = *sl->leptons;
            if (sl->leptons->size() > 0) m_triggerLeptons.push_back(sl->leptons->at(0));
        } else if (sel_type == BASE_DEN){
            // TODO: add function to return vector of leptons for this selection
            cout << "ERROR :: Not configured for this region yet\n";
        } else {
            cout << "ERROR :: Region configuration not yet defined\n";
        }
            
            
            // TODO: Add case of baseline denominator

        int n_selectleps = 0;
        for (Susy::Lepton* lep : m_selectLeptons) {
            if (!lep) continue;
            n_selectleps++;

            // For dilepton trigger matching
            if (lep->isEle()) {
                if (!m_el0) m_el0 = dynamic_cast<Susy::Electron*>(lep);
                else if (!m_el1) m_el1 = dynamic_cast<Susy::Electron*>(lep);
            }
            else if (lep->isMu()) {
                if (!m_mu0) m_mu0 = dynamic_cast<Susy::Muon*>(lep);
                else if (!m_mu1) m_mu1 = dynamic_cast<Susy::Muon*>(lep);
            }
        }

        if (m_selectLeptons.size() > 0 and m_selectLeptons.at(0)) {
            m_lepton0 = *m_selectLeptons[0];
        }
        if (m_selectLeptons.size() > 1 and m_selectLeptons.at(1)) {
            m_lepton1 = *m_selectLeptons[1];
        }
        m_dileptonP4 = m_lepton0 + m_lepton1;

        ////////////////////////////////////////////////////////////////////////
        // Jet shortcuts
        if (sl->baseJets->size() > 0) {
            m_Jet0_TLV = *sl->baseJets->at(0);
            if (sl->baseJets->size() > 1) {
               m_Jet1_TLV = *sl->baseJets->at(1);
               m_Dijet_TLV = m_Jet0_TLV + m_Jet1_TLV;
            }
        }
        // TODO: replace auto with explicit class
        for (auto& jet : *sl->baseJets) {
              if (sl->tools->jetSelector().isLight(jet))  {
                   m_lightJets.push_back(jet);
              } else if (sl->tools->jetSelector().isB(jet)) {
                  m_BJets.push_back(jet);
              } else if (sl->tools->jetSelector().isForward(jet))  {
                  m_forwardJets.push_back(jet);
              }
        }
        std::sort(m_lightJets.begin()  , m_lightJets.end()  , comparePt);
        std::sort(m_BJets.begin()      , m_BJets.end()      , comparePt);
        std::sort(m_forwardJets.begin(), m_forwardJets.end(), comparePt);

        return true; // All events pass this cut
    };
}
// TODO: Replace [=] with only the needed variables to see if it improves performance
void add_pre_cuts(Superflow* superflow, Sel sel_type) {
  // xTauFW Cut
  if (sel_type == BASELINE || sel_type == BASE_DEN) {
    *superflow << CutName("xTau: 2+ Loose Leptons") << [](Superlink* sl) -> bool {
        uint nLooseLeptons = 0;
        for (const auto* mu : *sl->preMuons) {if (mu->loose) nLooseLeptons++;}
        for (const auto* ele : *sl->preElectrons) {if (ele->looseLLH) nLooseLeptons++;}
        return nLooseLeptons >= 2;
    };
  }
}
void add_cleaing_cuts(Superflow* superflow) {
  *superflow << CutName("Pass GRL") << [=](Superlink* sl) -> bool {
      return (sl->tools->passGRL(m_cutflags));
  };
  *superflow << CutName("LAr error") << [=](Superlink* sl) -> bool {
      return (sl->tools->passLarErr(m_cutflags));
  };
  *superflow << CutName("Tile error") << [=](Superlink* sl) -> bool {
      return (sl->tools->passTileErr(m_cutflags));
  };
  *superflow << CutName("TTC veto") << [=](Superlink* sl) -> bool {
      return (sl->tools->passTTC(m_cutflags));
  };
  *superflow << CutName("SCT err") << [=](Superlink* sl) -> bool {
      return (sl->tools->passSCTErr(m_cutflags));
  };
}
void add_baseline_lepton_cuts(Superflow* superflow) {
    *superflow << CutName("nBaselineLep = nSignalLep") << [](Superlink* sl) -> bool {
        return (sl->leptons->size() == sl->baseLeptons->size());
    };
    *superflow << CutName("2+ leptons") << [](Superlink* sl) -> bool {
        return (sl->leptons->size() >= 2);
    };
    add_DFOS_lepton_cut(superflow);
}
void add_baseline_den_lepton_cuts(Superflow* superflow) {
    *superflow << CutName("1-ID Lepton and 1 Anti-ID Lepton") << [=](Superlink* /*sl*/) -> bool {
        return (m_lepID_n == 1 && m_lepAntiID_n == 1);
    };
    add_DFOS_lepton_cut(superflow);
}
void add_fake_num_lepton_cuts(Superflow* superflow) {
    *superflow << CutName("3-ID Leptons") << [=](Superlink* /*sl*/) -> bool {
        return (m_lepID_n == 3);
    };
    add_SFOS_lepton_cut(superflow);
}
void add_fake_den_lepton_cuts(Superflow* superflow) {
    *superflow << CutName("2-ID Leptons and 1+ Anti-ID Lepton") << [=](Superlink* /*sl*/) -> bool {
        return (m_lepID_n == 2 && m_lepAntiID_n >= 1);
    };
    add_SFOS_lepton_cut(superflow);
}
void add_final_cuts(Superflow* superflow, Sel sel_type) {
    *superflow << CutName("pass good vertex") << [=](Superlink* sl) -> bool {
        return (sl->tools->passGoodVtx(m_cutflags));
    };
    *superflow << CutName("jet cleaning") << [](Superlink* sl) -> bool {
        return (sl->tools->passJetCleaning(sl->baseJets));
    };
    if (sel_type == BASELINE) {
      *superflow << CutName("2 leptons") << [](Superlink* sl) -> bool {
          return (sl->leptons->size() == 2);
      };
    }
    *superflow << CutName("Tau veto") << [](Superlink* sl) -> bool {
        return (sl->taus->size() == 0);
    };
}
void add_event_variables(Superflow* superflow) {
  *superflow << NewVar("Event run number"); {
    *superflow << HFTname("RunNumber");
    *superflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->run; };
    *superflow << SaveVar();
  }
  *superflow << NewVar("Event number"); {
    *superflow << HFTname("event_number");
    *superflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->eventNumber; };
    *superflow << SaveVar();
  }
  *superflow << NewVar("is Monte Carlo"); {
    *superflow << HFTname("isMC");
    *superflow << [](Superlink* sl, var_bool*) -> bool { return sl->nt->evt()->isMC ? true : false; };
    *superflow << SaveVar();
  }
  *superflow << NewVar("event weight"); {
    *superflow << HFTname("eventweight");
    *superflow << [=](Superlink* sl, var_double*) -> double {
        // TODO: add fake factor tool to superflow
        float fakeFactor = 1;
        if (args.apply_ff && m_denominator_selection) {
            bool probeIsEl = m_antiID_lep0->isEle();
            LepEnum::LepType typeOfLep = probeIsEl ? LepEnum::Electron : LepEnum::Muon;
            fakeFactor = m_applyFakeFactorTool->apply(m_antiID_lep0_TLV.Pt(), typeOfLep);
        }
        return sl->weights->product() * sl->nt->evt()->wPileup * fakeFactor;
    };
    *superflow << SaveVar();
  }
  *superflow << NewVar("sample DSID"); {
    *superflow << HFTname("dsid");
    *superflow << [](Superlink* sl, var_int*) -> int {return sl->nt->evt()->mcChannel;};
    *superflow << SaveVar();
  }
  *superflow << NewVar("treatAsYear"); {
      // 15/16 Year ID
      *superflow << HFTname("treatAsYear");
      *superflow << [](Superlink* sl, var_double*) -> int { return sl->nt->evt()->treatAsYear; };
      *superflow << SaveVar();
  }
}
void add_trigger_variables(Superflow* superflow) {
    ////////////////////////////////////////////////////////////////////////////
    // Trigger Variables
    // ADD_*_TRIGGER_VAR preprocessor defined

    ////////////////////////////////////////////////////////////////////////////
    // 2015
    ADD_2LEP_TRIGGER_VAR(HLT_e17_lhloose_mu14, m_el0, m_mu0)
    ADD_2LEP_TRIGGER_VAR(HLT_e24_lhmedium_L1EM20VHI_mu8noL1, m_el0, m_mu0)
    ADD_2LEP_TRIGGER_VAR(HLT_e7_lhmedium_mu24, m_mu0, m_el0)
    ADD_2LEP_TRIGGER_VAR(HLT_2e12_lhloose_L12EM10VH, m_el0, m_el1)
    // TODO: Add to SusyNts HLT_2mu10)

    // Single Electron Triggers
    ADD_1LEP_TRIGGER_VAR(HLT_e24_lhmedium_L1EM20VH, m_triggerLeptons)
    ADD_1LEP_TRIGGER_VAR(HLT_e60_lhmedium, m_triggerLeptons)
    ADD_1LEP_TRIGGER_VAR(HLT_e120_lhloose, m_triggerLeptons)

    // Single Muon Triggers
    ADD_1LEP_TRIGGER_VAR(HLT_mu20_iloose_L1MU15, m_triggerLeptons)
    ADD_1LEP_TRIGGER_VAR(HLT_mu40, m_triggerLeptons)
    

    ////////////////////////////////////////////////////////////////////////////
    // 2016

    // Dilepton Triggers
    ADD_2LEP_TRIGGER_VAR(HLT_e17_lhloose_nod0_mu14, m_el0, m_mu0)
    // TODO: Add to SusyNts HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1, m_el0, m_mu0)
    ADD_2LEP_TRIGGER_VAR(HLT_e7_lhmedium_nod0_mu24, m_mu0, m_el0)
    ADD_2LEP_TRIGGER_VAR(HLT_2e17_lhvloose_nod0, m_el0, m_el1)
    // TODO: Add to SusyNts HLT_2mu14)

    // Single Electron Triggers
    ADD_1LEP_TRIGGER_VAR(HLT_e26_lhtight_nod0_ivarloose, m_triggerLeptons)
    ADD_1LEP_TRIGGER_VAR(HLT_e60_lhmedium_nod0, m_triggerLeptons)
    ADD_1LEP_TRIGGER_VAR(HLT_e140_lhloose_nod0, m_triggerLeptons)

    // Single Muon Triggers
    ADD_1LEP_TRIGGER_VAR(HLT_mu26_ivarmedium, m_triggerLeptons)
    ADD_1LEP_TRIGGER_VAR(HLT_mu50, m_triggerLeptons)
}
void add_met_variables(Superflow* superflow) {
  //////////////////////////////////////////////////////////////////////////////
  // MET

  // Fill MET variable inside Et var
  *superflow << NewVar("transverse missing energy (Et)"); {
    *superflow << HFTname("MET");
    *superflow << [=](Superlink* sl, var_double*) -> double {
      m_MET.SetPxPyPzE(sl->met->Et*cos(sl->met->phi),
                     sl->met->Et*sin(sl->met->phi),
                     0.,
                     sl->met->Et);
      return m_MET.Pt();
    };
    *superflow << SaveVar();
  }
  *superflow << NewVar("transverse missing energy (Phi)"); {
    *superflow << HFTname("METPhi");
    *superflow << [=](Superlink* /*sl*/, var_double*) -> double {
      return m_MET.Phi();
    };
    *superflow << SaveVar();
  }
}
void add_prelepton_variables(Superflow* superflow) {
    ADD_MULTIPLICITY_VAR(preLeptons)
    ADD_MULTIPLICITY_VAR(preElectrons)
    ADD_MULTIPLICITY_VAR(preMuons)
    //////////////////////////////////////////////////////////////////////////////
    // PreElectrons
    *superflow << NewVar("preElectron E_caloClus"); {
      *superflow << HFTname("preEl_EcaloClus");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : *sl->preElectrons) {
          out.push_back(el->clusE);
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("preElectron eta_caloClus"); {
      *superflow << HFTname("preEl_etaCaloClus");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : *sl->preElectrons) {
          out.push_back(el->clusEta);
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("preElectron Et"); {
      *superflow << HFTname("preEl_Et");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : *sl->preElectrons) {
          float E_caloClus = el->clusE;
          float eta_caloClus = el->clusEta;
          out.push_back(E_caloClus / std::cosh(eta_caloClus));
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("preElectron pt"); {
      *superflow << HFTname("preEl_pt");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : *sl->preElectrons) {out.push_back(el->Pt()); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("preElectron clusterEtaBE"); {
      *superflow << HFTname("preEl_clusEtaBE");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : *sl->preElectrons) {out.push_back(el->clusEtaBE);}
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("preElectron eta"); {
      *superflow << HFTname("preEl_eta");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : *sl->preElectrons) {out.push_back(el->Eta());}
        return out;
      };
      *superflow << SaveVar();
    }
    //////////////////////////////////////////////////////////////////////////////
    // PreMuons

    *superflow << NewVar("preMuon pt"); {
      *superflow << HFTname("preMu_pt");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& mu : *sl->preMuons) {out.push_back(mu->Pt()); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("isCaloTagged"); {
      *superflow << HFTname("isCaloTagged");
      *superflow << [](Superlink* sl, var_int_array*) -> vector<int> {
        vector<int> out;
        for (auto& mu : *sl->preMuons) { out.push_back(mu->isCaloTagged);}
        return out;
      };
      *superflow << SaveVar();
    }

    *superflow << NewVar("PreMuon ID (non-inclusive)"); {
      *superflow << HFTname("preMu_ID");
      *superflow << [](Superlink* sl, var_int_array*) -> vector<int> {
        vector<int> out;
        for (auto& mu : *sl->preMuons) {
          if (mu->tight) out.push_back(0);
          else if (mu->medium) out.push_back(1);
          else if (mu->loose) out.push_back(2);
          else if (mu->veryLoose) out.push_back(3);
          else out.push_back(4);
        }
        return out;
      };
      *superflow << SaveVar();
    }
}
void add_baselepton_variables(Superflow* superflow) {
    ADD_MULTIPLICITY_VAR(baseLeptons)
    ADD_MULTIPLICITY_VAR(baseElectrons)
    ADD_MULTIPLICITY_VAR(baseMuons)

    //////////////////////////////////////////////////////////////////////////////
    // Baseline Electrons
    *superflow << NewVar("Baseline Electron etconetopo20"); {
      *superflow << HFTname("baseEl_etconetopo20");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : *sl->baseElectrons) {out.push_back(el->etconetopo20); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Baseline Electron ptvarcone20"); {
      *superflow << HFTname("baseEl_ptvarcone20");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : *sl->baseElectrons) {out.push_back(el->ptvarcone20); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Baseline Electron ID (non-inclusive)"); {
      *superflow << HFTname("baseEl_ID");
      *superflow << [](Superlink* sl, var_int_array*) -> vector<int> {
        vector<int> out;
        for (auto& el : *sl->baseElectrons) {
          if (el->tightLLH) out.push_back(0);
          else if (el->mediumLLH) out.push_back(1);
          else if (el->looseLLHBLayer) out.push_back(3);
          else if (el->looseLLH) out.push_back(2);
          else if (el->veryLooseLLH) out.push_back(4);
          else out.push_back(5);
        }
        return out;
      };
      *superflow << SaveVar();
    }
    //////////////////////////////////////////////////////////////////////////////
    // Baseline Muons
    *superflow << NewVar("baseline Muon pt"); {
      *superflow << HFTname("baseMu_pt");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& mu : *sl->baseMuons) {out.push_back(mu->Pt()); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Baseline Muon eta"); {
      *superflow << HFTname("baseMu_eta");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& mu : *sl->baseMuons) {out.push_back(mu->Eta()); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Baseline Muon etconetopo20"); {
      *superflow << HFTname("baseMu_etconetopo20");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& mu : *sl->baseMuons) {out.push_back(mu->etconetopo20); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Baseline Muon ptvarcone30"); {
      *superflow << HFTname("baseMu_ptvarcone30");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& mu : *sl->baseMuons) {out.push_back(mu->ptvarcone30); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Baseline Muon ID (non-inclusive)"); {
      *superflow << HFTname("baseMu_ID");
      *superflow << [](Superlink* sl, var_int_array*) -> vector<int> {
        vector<int> out;
        for (auto& mu : *sl->baseMuons) {
          if (mu->tight) out.push_back(0);
          else if (mu->medium) out.push_back(1);
          else if (mu->loose) out.push_back(2);
          else if (mu->veryLoose) out.push_back(3);
          else out.push_back(4);
        }
        return out;
      };
      *superflow << SaveVar();
    }
}
void add_signallepton_variables(Superflow* superflow) {
    ADD_MULTIPLICITY_VAR(leptons)
    ADD_MULTIPLICITY_VAR(electrons)
    ADD_MULTIPLICITY_VAR(muons)
    //////////////////////////////////////////////////////////////////////////////
    // Signal Electrons
    *superflow << NewVar("Electron type"); {
      *superflow << HFTname("el_type");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : *sl->electrons) {out.push_back(el->mcType); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Electron origin"); {
      *superflow << HFTname("el_origin");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : *sl->electrons) {out.push_back(el->mcOrigin); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Electron ID (non-inclusive)"); {
      *superflow << HFTname("El_ID");
      *superflow << [](Superlink* sl, var_int_array*) -> vector<int> {
        vector<int> out;
        for (auto& el : *sl->electrons) {
          if (el->tightLLH) out.push_back(0);
          else if (el->mediumLLH) out.push_back(1);
          else if (el->looseLLHBLayer) out.push_back(3);
          else if (el->looseLLH) out.push_back(2);
          else if (el->veryLooseLLH) out.push_back(4);
          else out.push_back(5);
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("dRy(sigEl, preMu)"); {
      *superflow << HFTname("dRy_sEl_pMu");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : *sl->electrons) {
          for (auto& mu : *sl->preMuons) {
            out.push_back(el->DeltaRy(*mu));
          }
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("dRy(sigEl, preMu_shrdTrk)"); {
      *superflow << HFTname("dRy_sEl_pMu_shrdtrk");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : *sl->electrons) {
          for (auto& mu : *sl->preMuons) {
            if(!(el->sharedMuTrk.size()>0)) continue;
            if(mu->idx > ((int)el->sharedMuTrk.size()-1)) continue;
            if(el->sharedMuTrk[mu->idx]!=1) continue;
            out.push_back(el->DeltaRy(*mu));
          }
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("dRy(sigEl, baseMu not Calo)"); {
      *superflow << HFTname("dRy_sEl_bMu_noCalo");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : *sl->electrons) {
          for (auto& mu : *sl->preMuons) {
            if (!sl->tools->muonSelector().isBaseline(mu)) continue;
            if (sl->tools->muonSelector().isSignal(mu)) continue;
            if (mu->isCaloTagged) continue;
            out.push_back(el->DeltaRy(*mu));
          }
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("dRy(sigEl, baseMu CaloTag)"); {
      *superflow << HFTname("dRy_sEl_bMu_Calo");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : *sl->electrons) {
          for (auto& mu : *sl->preMuons) {
            if (!sl->tools->muonSelector().isBaseline(mu)) continue;
            if (sl->tools->muonSelector().isSignal(mu)) continue;
            if (!mu->isCaloTagged) continue;
            out.push_back(el->DeltaRy(*mu));
          }
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Electron d0sigBSCorr"); {
      *superflow << HFTname("el_d0sigBSCorr");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : *sl->electrons) {out.push_back(el->d0sigBSCorr); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Electron z0SinTheta"); {
      *superflow << HFTname("el_z0SinTheta");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : *sl->electrons) {out.push_back(el->z0SinTheta()); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("subleading electron track pt"); {
      *superflow << HFTname("el1_track_pt");
      *superflow << [](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() > 1 && sl->leptons->at(1)->isEle()) {
            const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(sl->leptons->at(1));
            return ele ? ele->trackPt : -DBL_MAX;
        } else {
            return -DBL_MAX;
        }
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("subleading electron clus pt"); {
      *superflow << HFTname("el1_clus_pt");
      *superflow << [](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() > 1 && sl->leptons->at(1)->isEle()) {
            const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(sl->leptons->at(1));
            return ele ? ele->clusE : -DBL_MAX;
        } else {
            return -DBL_MAX;
        }
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Subleading Electron pT track-cluster ratio"); {
      *superflow << HFTname("el1pT_trackclus_ratio");
      *superflow << [](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() > 1 && sl->leptons->at(1)->isEle()) {
            const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(sl->leptons->at(1));
            return ele ? ele->trackPt / ele->clusE : -DBL_MAX;
        } else {
            return -DBL_MAX;
        }
      };
      *superflow << SaveVar();
    }


    //////////////////////////////////////////////////////////////////////////////
    // Signal Muons
    *superflow << NewVar("Muon ID (non-inclusive)"); {
      *superflow << HFTname("Mu_ID");
      *superflow << [](Superlink* sl, var_int_array*) -> vector<int> {
        vector<int> out;
        for (auto& mu : *sl->muons) {
          if (!mu->veryLoose) out.push_back(4);
          else if (mu->veryLoose && !mu->loose) out.push_back(3);
          else if (mu->loose && !mu->medium) out.push_back(2);
          else if (mu->medium && !mu->tight) out.push_back(1);
          else if (mu->tight) out.push_back(0);
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Muon type"); {
      *superflow << HFTname("mu_type");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& mu : *sl->muons) {out.push_back(mu->mcType); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Muon origin"); {
      *superflow << HFTname("mu_origin");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& mu : *sl->muons) {out.push_back(mu->mcOrigin); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Muon d0sigBSCorr"); {
      *superflow << HFTname("mu_d0sigBSCorr");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& mu : *sl->muons) {out.push_back(mu->d0sigBSCorr); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Muon z0SinTheta"); {
      *superflow << HFTname("mu_z0SinTheta");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& mu : *sl->muons) {out.push_back(mu->z0SinTheta()); }
        return out;
      };
      *superflow << SaveVar();
    }
    //////////////////////////////////////////////////////////////////////////////
    // Signal Leptons
    *superflow << NewVar("Lepton Iso (non-inclusive)"); {
      *superflow << HFTname("Lep_Iso");
      *superflow << [=](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for (auto& lep :  m_selectLeptons) {
          if (!lep) continue;
          bool flag = false;
          out.push_back(-1);  // for tracking all entries and normalizing bins
          if (lep->isoGradient) {flag=true, out.push_back(0);}
          if (lep->isoGradientLoose) {flag=true, out.push_back(1);}
          if (lep->isoLoose) {flag=true, out.push_back(2);}
          if (lep->isoLooseTrackOnly) {flag=true, out.push_back(3);}
          if (lep->isoFixedCutTightTrackOnly) {flag=true; out.push_back(4);}
          if (!flag) out.push_back(5);
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("lepton pt"); {
      *superflow << HFTname("l_pt");
      *superflow << [=](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->Pt());
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("lepton eta"); {
      *superflow << HFTname("l_eta");
      *superflow << [=](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->Eta());
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("lepton phi"); {
      *superflow << HFTname("l_phi");
      *superflow << [=](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->Phi());
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("lepton flavor (0: e, 1: m)"); {
      *superflow << HFTname("l_flav");
      *superflow << [=](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->isEle() ? 0 : 1);
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("lepton type"); {
      *superflow << HFTname("l_type");
      *superflow << [=](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->mcType);
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("lepton origin"); {
      *superflow << HFTname("l_origin");
      *superflow << [=](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->mcOrigin);
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("lepton charge"); {
      *superflow << HFTname("l_q");
      *superflow << [=](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->q);
        }
        return out;
        };
      *superflow << SaveVar();
    }
    *superflow << NewVar("lepton BkgMotherPdgId"); {
      *superflow << HFTname("l_BkgMotherPdgId");
      *superflow << [=](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->mcBkgMotherPdgId);
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("lepton BkgTruthOrigin"); {
      *superflow << HFTname("l_BkgTruthOrigin");
      *superflow << [=](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->mcBkgTruthOrigin);
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("lepton matched2TruthLepton"); {
      *superflow << HFTname("l_matched2TruthLepton");
      *superflow << [=](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->matched2TruthLepton);
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("lepton classification"); {
        *superflow << HFTname("l_truthClass");
        *superflow << [=](Superlink* /*sl*/, var_int_array*) -> vector<int> {
            vector<int> out;
            for(auto& lepton : m_selectLeptons) {
                out.push_back(get_lepton_truth_class(lepton));
            }
            return out;
        };
        *superflow << SaveVar();
    }
}
void add_xtau_lepton_variables(Superflow* superflow) {
    // xTauFW variable
    *superflow << NewVar("leptons sign product"); {
      *superflow << HFTname("LepLepSign");
      *superflow << [=](Superlink* sl, var_int*) -> int {
        if (sl->leptons->size() < 2) return -INT_MAX;
        return m_selectLeptons.at(0)->q*m_selectLeptons.at(1)->q;
      };
      *superflow << SaveVar();
    }
    // Lepton variables
    // xTauFW variable
    *superflow << NewVar("Leading lepton pT"); {
      *superflow << HFTname("Lep0Pt");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() < 1) return -DBL_MAX;
        return m_lepton0.Pt(); };
      *superflow << SaveVar();
    }
    // xTauFW variable
    *superflow << NewVar("Subleading lepton pT"); {
      *superflow << HFTname("Lep1Pt");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() < 2) return -DBL_MAX;
        return m_lepton1.Pt(); };
      *superflow << SaveVar();
    }
    // xTauFW variable
    *superflow << NewVar("Leading lepton eta"); {
      *superflow << HFTname("Lep0Eta");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() < 1) return -DBL_MAX;
        return m_lepton0.Eta(); };
      *superflow << SaveVar();
    }
    // xTauFW variable
    *superflow << NewVar("Subleading lepton eta"); {
      *superflow << HFTname("Lep1Eta");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() < 2) return -DBL_MAX;
        return m_lepton1.Eta(); };
      *superflow << SaveVar();
    }
    // xTauFW variable
    *superflow << NewVar("Leading lepton phi"); {
      *superflow << HFTname("Lep0Phi");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() < 1) return -DBL_MAX;
        return m_lepton0.Phi(); };
      *superflow << SaveVar();
    }
    // xTauFW variable
    *superflow << NewVar("Subleading lepton phi"); {
      *superflow << HFTname("Lep1Phi");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() < 2) return -DBL_MAX;
        return m_lepton1.Phi(); };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Leading lepton invariant mass"); {
      *superflow << HFTname("MLep0");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() < 1) return -DBL_MAX;
        return m_lepton0.M(); };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Subleading lepton invariant mass"); {
      *superflow << HFTname("MLep1");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() < 2) return -DBL_MAX;
        return m_lepton1.M(); };
      *superflow << SaveVar();
    }
}
void add_tau_variables(Superflow* superflow) {
    ADD_MULTIPLICITY_VAR(preTaus)
    ADD_MULTIPLICITY_VAR(baseTaus)
    ADD_MULTIPLICITY_VAR(taus)
    //////////////////////////////////////////////////////////////////////////////
    // PreTaus
    *superflow << NewVar("preTau charge"); {
      *superflow << HFTname("preTau_q");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& tau : *sl->taus) {out.push_back(tau->q); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("preTau nTracks"); {
      *superflow << HFTname("preTau_nTracks");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& tau : *sl->taus) {out.push_back(tau->nTrack); }
        return out;
      };
      *superflow << SaveVar();
    }

    //////////////////////////////////////////////////////////////////////////////
    // Baseline Taus
    *superflow << NewVar("Baseline Tau pT"); {
      *superflow << HFTname("baseTau_pT");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& tau : *sl->taus) {out.push_back(tau->Pt()); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Baseline Tau eta"); {
      *superflow << HFTname("baseTau_eta");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& tau : *sl->taus) {out.push_back(tau->Eta()); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Baseline Tau nTracks"); {
      *superflow << HFTname("baseTau_nTracks");
      *superflow << [](Superlink* sl, var_int_array*) -> vector<int> {
        vector<int> out;
        for (auto& tau : *sl->taus) {out.push_back(tau->nTrack); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Baseline Tau ID (non-inclusive)"); {
      *superflow << HFTname("baseTau_ID");
      *superflow << [](Superlink* sl, var_int_array*) -> vector<int> {
        vector<int> out;
        for (auto& tau : *sl->taus) {
          if (!tau->loose) out.push_back(3);
          else if (tau->loose && !tau->medium) out.push_back(2);
          else if (tau->medium && !tau->tight) out.push_back(1);
          else if (tau->tight) out.push_back(0);
        }
        return out;
      };
      *superflow << SaveVar();
    }
}
void add_other_lepton_variables(Superflow* superflow) {
    //////////////////////////////////////////////////////////////////////////////
    // Two-lepton properties
    *superflow << NewVar("Delta eta between leptons"); {
      *superflow << HFTname("DEtaLL");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() < 2) return -DBL_MAX;
        return fabs(m_lepton0.Eta() - m_lepton1.Eta()); };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Delta phi between leptons"); {
      *superflow << HFTname("DphiLL");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() < 2) return -DBL_MAX;
        return fabs(m_lepton0.DeltaPhi(m_lepton1));};
      *superflow << SaveVar();
    }
    *superflow << NewVar("delta R of di-lepton system"); {
      *superflow << HFTname("drll");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() < 2) return -DBL_MAX;
        return m_lepton0.DeltaR(m_lepton1); };
      *superflow << SaveVar();
    }
    *superflow << NewVar("dilepton flavor"); {
      *superflow << HFTname("dilep_flav");
      *superflow << [=](Superlink* sl, var_int*) -> int {
        if(sl->leptons->size()<2) return -INT_MAX;
        if(sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isMu()){return 0;}       // e mu  case
        else if(sl->leptons->at(0)->isMu() && sl->leptons->at(1)->isEle()){return 1;}  // mu e  case
        else if(sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isEle()){return 2;} // e e   case
        else if(sl->leptons->at(0)->isMu() && sl->leptons->at(1)->isMu()){return 3;}   // mu mu case
        return 4;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("dilepton flavor is e mu"); {
      *superflow << HFTname("isEM");
      *superflow << [=](Superlink* sl, var_int*) -> int {
          if(sl->leptons->size()<2) return -INT_MAX;
        return sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isMu();  // e mu  case
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("dilepton flavor is mu e"); {
      *superflow << HFTname("isME");
      *superflow << [=](Superlink* sl, var_int*) -> int {
          if(sl->leptons->size()<2) return -INT_MAX;
          return sl->leptons->at(0)->isMu() && sl->leptons->at(1)->isEle();  // mu e  case
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("collinear mass, M_coll"); {
      *superflow << HFTname("MCollASym");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() < 2) return -DBL_MAX;
        double deta = fabs(m_lepton0.Eta()-m_lepton1.Eta());
          double dphi = m_lepton0.DeltaPhi(m_lepton1);
          return sqrt(2*m_lepton0.Pt()*(m_lepton1.Pt()+sl->met->Et)*(cosh(deta)-cos(dphi)));
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("mass of di-lepton system, M_ll"); {
      *superflow << HFTname("MLL");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() < 2) return -DBL_MAX;
        return m_dileptonP4.M(); };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Pt of di-lepton system, Pt_ll"); {
      *superflow << HFTname("ptll");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() < 2) return -DBL_MAX;
        return m_dileptonP4.Pt(); };
      *superflow << SaveVar();
    }
    *superflow << NewVar("diff_pt between leading and sub-leading lepton"); {
      *superflow << HFTname("dpt_ll");
      *superflow << [=](Superlink* sl, var_double*) -> double {
          if (sl->leptons->size() < 2) return -DBL_MAX;
          return m_lepton0.Pt() - m_lepton1.Pt();
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("dphi between leading lepton and MET"); {
      *superflow << HFTname("DphiLep0MET");
      *superflow << [=](Superlink* sl, var_double*) -> double {
          if (sl->leptons->size() < 1) return -DBL_MAX;
          return m_lepton0.DeltaPhi(m_MET);
      };
      *superflow << SaveVar();
    }

    //////////////////////////////////////////////////////////////////////////////
    // Optimized Region Cuts
    *superflow << NewVar("dphi between sub-leading lepton and MET"); {
      *superflow << HFTname("DphiLep1MET");
      *superflow << [=](Superlink* sl, var_double*) -> double {
          if (sl->leptons->size() < 2) return -DBL_MAX;
          return m_lepton1.DeltaPhi(m_MET);
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("transverse mass of leading lepton"); {
      *superflow << HFTname("MtLep0");
      *superflow << [=](Superlink* sl, var_double*) -> double {
          if (sl->leptons->size() < 1) return -DBL_MAX;
          return sqrt( 2*m_lepton0.Pt()*m_MET.Pt()*(1-cos (m_lepton0.DeltaPhi(m_MET)) ) );
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("transverse mass of subleading lepton"); {
      *superflow << HFTname("MtLep1");
      *superflow << [=](Superlink* sl, var_double*) -> double {
          if (sl->leptons->size() < 2) return -DBL_MAX;
          return sqrt( 2*m_lepton1.Pt()*m_MET.Pt()*(1-cos (m_lepton1.DeltaPhi(m_MET)) ) );
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("approximate tau pT"); {
      *superflow << HFTname("tau_pT");
      *superflow << [=](Superlink* sl, var_double*) -> double {
          if (sl->leptons->size() < 2) return -DBL_MAX;
          return (m_MET + m_lepton1).Pt();
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("ratio of tau pT to subleading lepton pT"); {
      *superflow << HFTname("taulep1_pT_ratio");
      *superflow << [=](Superlink* sl, var_double*) -> double {
          if (sl->leptons->size() < 2) return -DBL_MAX;
          TLorentzVector tau_estimate = m_MET + m_lepton1;
          return tau_estimate.Pt() / m_lepton0.Pt();
      };
      *superflow << SaveVar();
    }
}
void add_jet_variables(Superflow* superflow) {
    ////////////////////////////////////////////////////////////////////////////
    // Multiplicity variables
    ADD_MULTIPLICITY_VAR(preJets)
    ADD_MULTIPLICITY_VAR(baseJets)
    ADD_MULTIPLICITY_VAR(jets)
    //xTauFW variable
    *superflow << NewVar("number of baseline jets"); {
      *superflow << HFTname("JetN");
      *superflow << [](Superlink* sl, var_int*) -> int { return sl->baseJets->size(); };
      *superflow << SaveVar();
    }
    //xTauFW variable
    *superflow << NewVar("number of baseline jets (2p4Eta25Pt)"); {
      *superflow << HFTname("Jet_N2p4Eta25Pt");
      *superflow << [](Superlink* sl, var_int*) -> int {
          uint nPassedJets = 0;
          for (auto& jet : *sl->baseJets) {
              if (fabs(jet->Eta()) < 2.4 && jet->Pt() > 25){
                  nPassedJets += 1;
              }
          }
          return nPassedJets;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("number of signal jets"); {
      *superflow << HFTname("SignalJetN");
      *superflow << [](Superlink* sl, var_int*) -> int { return sl->jets->size(); };
      *superflow << SaveVar();
    }
    *superflow << NewVar("number of light jets"); {
      *superflow << HFTname("nLJets");
      *superflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfLJets(*sl->jets)/*(**sl->baseJets)*/; };
      *superflow << SaveVar();
    }
    *superflow << NewVar("number of b jets"); {
      *superflow << HFTname("nBJets");
      *superflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfBJets(*sl->jets)/*(**sl->baseJets)*/; };
      *superflow << SaveVar();
    }
    //xTauFW variable
    *superflow << NewVar("b-tagged jet"); {
      *superflow << HFTname("Btag");
      *superflow << [](Superlink* sl, var_bool*) -> bool { return sl->tools->numberOfBJets(*sl->jets) > 0;};
      *superflow << SaveVar();
    }
    *superflow << NewVar("number of forward jets"); {
      *superflow << HFTname("nForwardJets");
      *superflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfFJets(*sl->jets)/*(**sl->baseJets)*/; };
      *superflow << SaveVar();
    }
    //////////////////////////////////////////////////////////////////////////////
    // Pre Jets
    *superflow << NewVar("preJet pt"); {
      *superflow << HFTname("preJet_pt");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& jet : *sl->preJets) {out.push_back(jet->Pt()); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("preJet eta"); {
      *superflow << HFTname("preJet_eta");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& jet : *sl->preJets) {out.push_back(jet->Eta()); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("preJet JVT if eta<=2.4 and pT<60"); {
      *superflow << HFTname("preJet_JVT");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& jet : *sl->preJets) {
          if (jet->Pt() < 60 && fabs(jet->Eta()) <= 2.4) out.push_back(jet->jvt);
          else
            out.push_back(-DBL_MAX);
        }
        return out;
      };
      *superflow << SaveVar();
    }

    //////////////////////////////////////////////////////////////////////////////
    // Baseline Jets
    *superflow << NewVar("Baseline Jet eta"); {
      *superflow << HFTname("baseJet_eta");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& jet : *sl->baseJets) {out.push_back(jet->Eta()); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Baseline Jet mv2c10"); {
      *superflow << HFTname("baseJet_mv2c10");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& jet : *sl->baseJets) {out.push_back(jet->mv2c10); }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("jet eta"); {
      *superflow << HFTname("j_eta");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for(auto& jet : *sl->baseJets) {
          out.push_back(jet->Eta());
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("jet JVT"); {
      *superflow << HFTname("j_jvt");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for(auto& jet : *sl->baseJets) {
          out.push_back(jet->jvt);
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("jet JVF"); {
      *superflow << HFTname("j_jvf");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for(auto& jet : *sl->baseJets) {
          out.push_back(jet->jvf);
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("jet phi"); {
      *superflow << HFTname("j_phi");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for(auto& jet : *sl->baseJets) {
          out.push_back(jet->Phi());
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("jet flavor (0: NA, 1: L, 2: B, 3: F)"); {
      *superflow << HFTname("j_flav");
      *superflow << [](Superlink* sl, var_int_array*) -> vector<int> {
        vector<int> out; int flav = 0;
        for(auto& jet : *sl->baseJets) {
          if(sl->tools->jetSelector().isLight(jet))  { flav = 1; }
          else if(sl->tools->jetSelector().isB(jet)) { flav = 2; }
          else if(sl->tools->jetSelector().isForward(jet))  { flav = 3; }
          out.push_back(flav);
          flav=0;
        }
        return out;
      };
      *superflow << SaveVar();
    }

    //////////////////////////////////////////////////////////////////////////////
    // VBF Region Cuts
    *superflow << NewVar("number of jets pT > 30"); {
      *superflow << HFTname("JetN_g30");
      *superflow << [](Superlink* sl, var_int*) -> int {
          uint nJets = 0;
          for (auto& jet : *sl->baseJets) {
              if (jet->Pt() > 30) nJets++;
          }
          return nJets;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("jet pt"); {
      *superflow << HFTname("j_pt");
      *superflow << [](Superlink* sl, var_float_array*) -> vector<double> {
        vector<double> out;
        for(auto& jet : *sl->baseJets) {
          out.push_back(jet->Pt());
        }
        return out;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("dijet mass"); {
      *superflow << HFTname("Mjj");
      *superflow << [](Superlink* sl, var_double*) -> double {
          if (sl->baseJets->size() < 2) return -DBL_MAX;
          return m_Dijet_TLV.M();
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("DeltaEta between two leading jets"); {
      *superflow << HFTname("DEtaJJ");
      *superflow << [](Superlink* sl, var_double*) -> double {
          if (sl->baseJets->size() < 2) return -DBL_MAX;
          return fabs(m_Jet0_TLV.Eta() - m_Jet1_TLV.Eta());
      };
      *superflow << SaveVar();
    }
}
void add_fake_variables(Superflow* superflow) {
  //////////////////////////////////////////////////////////////////////////////
  // Fake variables
  //////////////////////////////////////////////////////////////////////////////
    *superflow << NewVar("Number of ID leptons"); {
      *superflow << HFTname("nLepID");
      *superflow << [=](Superlink* /*sl*/, var_int*) -> int {
          return m_lepID_n;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Number of Anti-ID leptons"); {
      *superflow << HFTname("nLepAntiID");
      *superflow << [=](Superlink* /*sl*/, var_int*) -> int {
          return m_lepAntiID_n;
      };
      *superflow << SaveVar();
    }
    //////////////////////////////////////////////////////////////////////////////
    // 1 antiID lepton case
    *superflow << NewVar("Leading anti-ID lepton charge"); {
        *superflow << HFTname("aID_Lep0Q");
        *superflow << [=](Superlink* /*sl*/, var_int*) -> int {
            // TODO: move lepton requirements outside of larger function and
            // remove them from the individual lambda expressions
            if (m_lepAntiID_n < 1) return -INT_MAX;
            return m_antiID_lep0->q;
        };
        *superflow << SaveVar();
    }
    *superflow << NewVar("Leading aID-lepton pT"); {
      *superflow << HFTname("aID_Lep0Pt");
      *superflow << [=](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 1) return -DBL_MAX;
        return m_antiID_lep0_TLV.Pt(); };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Leading aID-lepton eta"); {
      *superflow << HFTname("aID_Lep0Eta");
      *superflow << [=](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 1) return -DBL_MAX;
        return m_antiID_lep0_TLV.Eta(); };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Leading aID-lepton invariant mass"); {
      *superflow << HFTname("aID_MLep0");
      *superflow << [=](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 1) return -DBL_MAX;
        return m_antiID_lep0_TLV.M(); };
      *superflow << SaveVar();
    }
    *superflow << NewVar("transverse mass of leading aID-lepton"); {
      *superflow << HFTname("aID_MtLep0");
      *superflow << [=](Superlink* /*sl*/, var_double*) -> double {
          if (m_lepAntiID_n < 1) return -DBL_MAX;
          return sqrt( 2*m_antiID_lep0_TLV.Pt()*m_MET.Pt()*(1-cos (m_antiID_lep0_TLV.DeltaPhi(m_MET)) ) );
      };
      *superflow << SaveVar();
    }
    // Z + jets variables
    *superflow << NewVar("Z leptons' flavor"); {
      *superflow << HFTname("Z_dilep_flav");
      *superflow << [=](Superlink* sl, var_int*) -> int {
        if(sl->leptons->size() < 2) return -INT_MAX;
        if(m_Zlep.at(0)->isEle() && m_Zlep.at(1)->isMu()){return 0;}       // e mu  case
        else if(m_Zlep.at(0)->isMu() && m_Zlep.at(1)->isEle()){return 1;}  // mu e  case
        else if(m_Zlep.at(0)->isEle() && m_Zlep.at(1)->isEle()){return 2;} // e e   case
        else if(m_Zlep.at(0)->isMu() && m_Zlep.at(1)->isMu()){return 3;}   // mu mu case
        return 4;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Z dilepton sign"); {
      *superflow << HFTname("Z_dilep_sign");
      *superflow << [=](Superlink* sl, var_int*) -> int {
        if(sl->leptons->size() < 2) return -INT_MAX;
        return m_Zlep.at(0)->q * m_Zlep.at(1)->q;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("MLL of closest Z pair"); {
      *superflow << HFTname("Z_MLL");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() < 2) return -DBL_MAX;
        return (*m_Zlep.at(0)+*m_Zlep.at(1)).M();
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("2nd Z leptons' flavor"); {
      *superflow << HFTname("Z2_dilep_flav");
      *superflow << [=](Superlink* sl, var_int*) -> int {
        if(sl->leptons->size() < 2) return -INT_MAX;
        else if (!m_Zlep.at(3) or !m_Zlep.at(4)) return -2;
        if(m_Zlep.at(3)->isEle() && m_Zlep.at(4)->isMu()){return 0;}       // e mu  case
        else if(m_Zlep.at(3)->isMu() && m_Zlep.at(4)->isEle()){return 1;}  // mu e  case
        else if(m_Zlep.at(3)->isEle() && m_Zlep.at(4)->isEle()){return 2;} // e e   case
        else if(m_Zlep.at(3)->isMu() && m_Zlep.at(4)->isMu()){return 3;}   // mu mu case
        return 4;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("2nd Z dilepton sign"); {
      *superflow << HFTname("Z2_dilep_sign");
      *superflow << [=](Superlink* sl, var_int*) -> int {
        if(sl->leptons->size() < 2) return -INT_MAX;
        else if (!m_Zlep.at(3) or !m_Zlep.at(4)) return -2;
        return m_Zlep.at(3)->q * m_Zlep.at(4)->q;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("MLL of 2nd-closest Z pair"); {
      *superflow << HFTname("Z2_MLL");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (sl->leptons->size() < 2) return -DBL_MAX;
        else if (!m_Zlep.at(3) or !m_Zlep.at(4)) return -5.0;
        return (*m_Zlep.at(3)+*m_Zlep.at(4)).M();
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Non-Z lepton pT"); {
      *superflow << HFTname("Z_Lep2_pT");
      *superflow << [=](Superlink* /*sl*/, var_double*) -> double {
        if (m_Zlep.at(2)) return m_Zlep.at(2)->Pt();
        return -DBL_MAX;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Non-Z lepton eta"); {
      *superflow << HFTname("Z_Lep2_eta");
      *superflow << [=](Superlink* /*sl*/, var_double*) -> double {
        if (m_Zlep.at(2)) return m_Zlep.at(2)->Eta();
        return -DBL_MAX;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Non-Z lepton flavor"); {
      *superflow << HFTname("Z_Lep2_flav");
      *superflow << [=](Superlink* /*sl*/, var_int*) -> int {
        if (m_Zlep.at(2)) return m_Zlep.at(2)->isEle();
        return -INT_MAX;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Non-Z lepton charge"); {
      *superflow << HFTname("Z_Lep2_q");
      *superflow << [=](Superlink* /*sl*/, var_int*) -> int {
        if (m_Zlep.at(2)) return m_Zlep.at(2)->q;
        return -INT_MAX;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("DeltaPhi of Non-Z lep and MET"); {
      *superflow << HFTname("Z_Lep2_dPhi_MET");
      *superflow << [=](Superlink* /*sl*/, var_double*) -> double {
        if (m_Zlep.at(2)) return m_Zlep.at(2)->DeltaPhi(m_MET);
        return -DBL_MAX;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Non-Z lepton mT"); {
      *superflow << HFTname("Z_Lep2_mT");
      *superflow << [=](Superlink* /*sl*/, var_double*) -> double {
        if (m_Zlep.at(2)) {
            return sqrt( 2*m_Zlep.at(2)->Pt()*m_MET.Pt()*(1-cos (m_Zlep.at(2)->DeltaPhi(m_MET)) ) );
        }
        return -DBL_MAX;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("delta R of Z and aID-Lepton"); {
      *superflow << HFTname("dR_Zl");
      *superflow << [=](Superlink* /*sl*/, var_double*) -> double {
        if (m_Zlep.at(2)) return m_Zlep.at(2)->DeltaR(*m_Zlep.at(0)+*m_Zlep.at(1));
        return -DBL_MAX;
      };
      *superflow << SaveVar();
    }
    // W + jets variables
    *superflow << NewVar("aID-ID dilepton flavor"); {
      *superflow << HFTname("aID_dilep_flav");
      *superflow << [=](Superlink* sl, var_int*) -> int {
        if (m_lepAntiID_n < 1 || sl->leptons->size() < 1) return -INT_MAX;
        auto* lep0 = sl->leptons->at(0);
        if(m_antiID_lep0->isEle() && lep0->isMu()){return 0;}       // e mu  case
        else if(m_antiID_lep0->isMu() && lep0->isEle()){return 1;}  // mu e  case
        else if(m_antiID_lep0->isEle() && lep0->isEle()){return 2;} // e e   case
        else if(m_antiID_lep0->isMu() && lep0->isMu()){return 3;}   // mu mu case
        return 4;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("aID-ID sign product"); {
      *superflow << HFTname("aID_LepLepSign");
      *superflow << [=](Superlink* sl, var_int*) -> int {
        if (m_lepAntiID_n < 1 || sl->leptons->size() < 1) return -INT_MAX;
        return m_antiID_lep0->q * sl->leptons->at(0)->q;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("delta R of ID-lepton and aID-lepton"); {
      *superflow << HFTname("aID_drll");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (m_lepAntiID_n < 1 || sl->leptons->size() < 1) return -DBL_MAX;
        return m_antiID_lep0_TLV.DeltaR(m_lepton0); };
      *superflow << SaveVar();
    }
    *superflow << NewVar("ID-aID collinear mass, M_coll"); {
      *superflow << HFTname("aID_MCollASym");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (m_lepAntiID_n < 1 || sl->leptons->size() < 1) return -DBL_MAX;
        double deta = fabs(m_antiID_lep0_TLV.Eta()-m_lepton0.Eta());
        double dphi = m_antiID_lep0_TLV.DeltaPhi(m_lepton0);
        return sqrt(2*m_antiID_lep0_TLV.Pt()*(m_lepton0.Pt()+sl->met->Et)*(cosh(deta)-cos(dphi)));
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("mass of ID-aID-dilepton system, M_ll"); {
      *superflow << HFTname("aID_MLL");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (m_lepAntiID_n < 1 || sl->leptons->size() < 1) return -DBL_MAX;
        return (m_antiID_lep0_TLV+m_lepton0).M();
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Pt of ID-aID-dilepton system, Pt_ll"); {
      *superflow << HFTname("aID_ptll");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (m_lepAntiID_n < 1 || sl->leptons->size() < 1) return -DBL_MAX;
        return (m_antiID_lep0_TLV+m_lepton0).Pt(); };
      *superflow << SaveVar();
    }
    *superflow << NewVar("diff_pt between ID-aID-leptons"); {
      *superflow << HFTname("aID_dpt_ll");
      *superflow << [=](Superlink* sl, var_double*) -> double {
          if (m_lepAntiID_n < 1 || sl->leptons->size() < 1) return -DBL_MAX;
          return m_antiID_lep0_TLV.Pt() - m_lepton0.Pt();
      };
      *superflow << SaveVar();
    }

    //////////////////////////////////////////////////////////////////////////////
    // 2 antiID lepton case (QCD)
    *superflow << NewVar("Subleading anti-ID lepton flavor"); {
      *superflow << HFTname("aID_Lep1Flav");
      *superflow << [=](Superlink* /*sl*/, var_int*) -> int {
          if (m_lepAntiID_n < 2) return -INT_MAX;
          return  m_antiID_lep1->isEle() ? 0 : 1;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Subleading anti-ID lepton charge"); {
      *superflow << HFTname("aID_Lep1Q");
      *superflow << [=](Superlink* /*sl*/, var_int*) -> int {
        if (m_lepAntiID_n < 2) return -INT_MAX;
        return m_antiID_lep1->q;
        };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Subleading aID-lepton pT"); {
      *superflow << HFTname("aID2_Lep1Pt");
      *superflow << [=](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 2) return -DBL_MAX;
        return m_antiID_lep1_TLV.Pt(); };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Subleading aID-lepton eta"); {
      *superflow << HFTname("aID2_Lep1Eta");
      *superflow << [=](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 2) return -DBL_MAX;
        return m_antiID_lep1_TLV.Eta(); };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Subleading aID-lepton invariant mass"); {
      *superflow << HFTname("aID2_MLep1");
      *superflow << [=](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 2) return -DBL_MAX;
        return m_antiID_lep1_TLV.M(); };
      *superflow << SaveVar();
    }
    *superflow << NewVar("delta R of aID-leptons"); {
      *superflow << HFTname("aID2_drll");
      *superflow << [=](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 2) return -DBL_MAX;
        return m_antiID_lep0_TLV.DeltaR(m_antiID_lep1_TLV); };
      *superflow << SaveVar();
    }
    *superflow << NewVar("aID-dilepton flavor"); {
      *superflow << HFTname("aID2_dilep_flav");
      *superflow << [=](Superlink* /*sl*/, var_int*) -> int {
        if(m_lepAntiID_n < 2) return -INT_MAX;
        auto* lep0 = m_antiID_lep0;
        auto* lep1 = m_antiID_lep1;
        if(lep0->isEle() && lep1->isMu()){return 0;}       // e mu  case
        else if(lep0->isMu() && lep1->isEle()){return 1;}  // mu e  case
        else if(lep0->isEle() && lep1->isEle()){return 2;} // e e   case
        else if(lep0->isMu() && lep1->isMu()){return 3;}   // mu mu case
        return 4;
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("aID collinear mass, M_coll"); {
      *superflow << HFTname("aID2_MCollASym");
      *superflow << [=](Superlink* sl, var_double*) -> double {
        if (m_lepAntiID_n < 2) return -DBL_MAX;
        double deta = fabs(m_antiID_lep0_TLV.Eta()-m_antiID_lep1_TLV.Eta());
        double dphi = m_antiID_lep0_TLV.DeltaPhi(m_antiID_lep1_TLV);
        return sqrt(2*m_antiID_lep0_TLV.Pt()*(m_antiID_lep1_TLV.Pt()+sl->met->Et)*(cosh(deta)-cos(dphi)));
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("mass of aID-dilepton system, M_ll"); {
      *superflow << HFTname("aID2_MLL");
      *superflow << [=](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 2) return -DBL_MAX;
        return (m_antiID_lep0_TLV+m_antiID_lep1_TLV).M(); };
      *superflow << SaveVar();
    }
    *superflow << NewVar("Pt of aID-dilepton system, Pt_ll"); {
      *superflow << HFTname("aID2_ptll");
      *superflow << [=](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 2) return -DBL_MAX;
        return (m_antiID_lep0_TLV+m_antiID_lep1_TLV).Pt(); };
      *superflow << SaveVar();
    }
    *superflow << NewVar("diff_pt between aID-leptons"); {
      *superflow << HFTname("aID2_dpt_ll");
      *superflow << [=](Superlink* /*sl*/, var_double*) -> double {
          if (m_lepAntiID_n < 2) return -DBL_MAX;
          return m_antiID_lep0_TLV.Pt() - m_antiID_lep1_TLV.Pt();
      };
      *superflow << SaveVar();
    }
    *superflow << NewVar("transverse mass of subleading aID-lepton"); {
      *superflow << HFTname("aID2_MtLep1");
      *superflow << [=](Superlink* /*sl*/, var_double*) -> double {
          if (m_lepAntiID_n < 2) return -DBL_MAX;
          return sqrt( 2*m_antiID_lep1_TLV.Pt()*m_MET.Pt()*(1-cos (m_antiID_lep1_TLV.DeltaPhi(m_MET)) ) );
      };
      *superflow << SaveVar();
    }
}
void add_shortcut_variables_reset(Superflow* superflow) {
  *superflow << [=](Superlink* /*sl*/, var_void*) {
    m_cutflags = 0;
    m_MET = {};
    m_lepID_n = m_lepAntiID_n = 0;
    m_antiID_lep0_TLV = m_antiID_lep1_TLV = {};
    std::fill(m_Zlep.begin(), m_Zlep.end(), nullptr);
    m_selectLeptons.clear(); m_triggerLeptons.clear();
    m_dileptonP4 = m_lepton0 = m_lepton1 = {};
    m_el0 = m_el1 = 0;
    m_mu0 = m_mu1 = 0;
    m_lightJets.clear(); m_BJets.clear(); m_forwardJets.clear();
    m_Dijet_TLV = m_Jet1_TLV = m_Jet0_TLV = {};
  };
}

////////////////////////////////////////////////////////////////////////////////
// Useful functions
void add_fake_shortcut_variables(Superlink* sl) {

    vector<int> lepID_vec, lepAntiID_vec;
    // ID Leptons - same as signal leptons
    for (Susy::Lepton* lepton : *sl->baseLeptons) {
        bool lepID_bool = is_ID_lepton(sl, lepton);
        lepID_vec.push_back(lepID_bool);
        if (lepID_bool) m_lepID_n++;
    }
    // Anti-ID Leptons
    for (Susy::Lepton* lepton : *sl->baseLeptons) {
        bool lepAntiID_bool = is_antiID_lepton(lepton);
        lepAntiID_vec.push_back(lepAntiID_bool);
      if (lepAntiID_bool) m_lepAntiID_n++;
    }

    std::vector<int>::iterator begin = lepAntiID_vec.begin();
    std::vector<int>::iterator end = lepAntiID_vec.end();
    int antiID_idx0 = -1;
    int antiID_idx1 = -1;
    if (m_lepAntiID_n >= 1) {
      antiID_idx0 = std::find(begin, end, true) - begin;
      m_antiID_lep0 = sl->baseLeptons->at(antiID_idx0);
      m_antiID_lep0_TLV = *m_antiID_lep0;
    }
    if (m_lepAntiID_n >= 2) {
      antiID_idx1 = std::find(begin + antiID_idx0 + 1, end, true) - begin;
      m_antiID_lep1 = sl->baseLeptons->at(antiID_idx1);
      m_antiID_lep1_TLV = *m_antiID_lep1;
    }

    // Find ID lepton combination with MLL closest to Z-peak
    float Z_diff = FLT_MAX;
    int zlep_idx2 = -1;
    for (uint ii = 0; ii < sl->leptons->size(); ++ii) {
        Susy::Lepton *lep_ii = sl->leptons->at(ii);
        for (uint jj = ii+1; jj < sl->leptons->size(); ++jj) {
            Susy::Lepton *lep_jj = sl->leptons->at(jj);
            float Z_diff_cf = fabs((*lep_ii+*lep_jj).M() - 91.2);
            if (Z_diff_cf < Z_diff) {
                Z_diff = Z_diff_cf;
                m_Zlep.at(0) = lep_ii;
                m_Zlep.at(1) = lep_jj;
                zlep_idx2 = ii > 0 ? 0 : jj > 1 ? 1 : 2;
            }
        }
    }
    if (sl->leptons->size() == 2 && m_lepAntiID_n >= 1) {
        m_Zlep.at(2) = m_antiID_lep0;
    } else if (sl->leptons->size() >= 3) {
        m_Zlep.at(2) = sl->leptons->at(zlep_idx2);
    }
    // Find base lepton combination with MLL closest to Z-peak
    // Exclude ID leptons already matched with Z
    Z_diff = FLT_MAX;
    for (uint ii = 0; ii < sl->baseLeptons->size(); ++ii) {
        Susy::Lepton *lep_ii = sl->baseLeptons->at(ii);
        if (lep_ii == m_Zlep.at(0) || lep_ii == m_Zlep.at(1)) continue;
        for (uint jj = ii+1; jj < sl->baseLeptons->size(); ++jj) {
            Susy::Lepton *lep_jj = sl->baseLeptons->at(jj);
            if (lep_jj == m_Zlep.at(0) || lep_jj == m_Zlep.at(1)) continue;
            float Z_diff_cf = fabs((*lep_ii+*lep_jj).M() - 91.2);
            if (Z_diff_cf < Z_diff) {
                Z_diff = Z_diff_cf;
                m_Zlep.at(3) = lep_ii;
                m_Zlep.at(4) = lep_jj;
            }
        }
    }
}
bool is_ID_lepton(Superlink* sl, Susy::Lepton* lepton) {
    for (Susy::Lepton* sig_lepton : *sl->leptons) {
        if (lepton == sig_lepton) {
            return true;
        }
    }
    return false;
}

bool is_antiID_lepton(Susy::Lepton* lepton) {
    bool pt_pass = 0, eta_pass = 0, iso_pass = 0, id_pass = 0;
    bool passID_cuts = 0, passAntiID_cuts = 0;
    if (lepton->isEle()) {
        const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(lepton);
        float absEtaBE = fabs(ele->clusEtaBE);
        pt_pass  = ele->pt > 15;
        eta_pass = absEtaBE < 1.37 || 1.52 < absEtaBE;
        iso_pass = ele->isoGradient;
        id_pass  = ele->mediumLLH;
        passID_cuts = iso_pass && id_pass;
        passAntiID_cuts = ele->veryLooseLLH;
    } else if (lepton->isMu()) {
        const Susy::Muon* mu = dynamic_cast<const Susy::Muon*>(lepton);
        pt_pass  = mu->pt > 10;
        eta_pass = fabs(mu->eta) < 2.47;
        iso_pass = mu->isoGradient;
        id_pass  = mu->medium;
        //TODO: Study effect of removing id_pass so muons only fail iso requirements
        passID_cuts = iso_pass && id_pass;
        passAntiID_cuts = 1;
    }
    return pt_pass && eta_pass && passAntiID_cuts && !passID_cuts;
}

bool is_1lep_trig_matched(Superlink* sl, string trig_name, LeptonVector leptons) {
    for (Susy::Lepton* lep : leptons) {
        if(!lep) continue;
        bool trig_fired = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, trig_name);
        if (!trig_fired) continue;
        bool trig_matched = sl->tools->triggerTool().lepton_trigger_match(lep, trig_name);
        if (trig_matched) return true;
    }
    return false;
}

void add_SFOS_lepton_cut(Superflow* superflow) {
    *superflow << CutName("SFOS leptons") << [=](Superlink* /*sl*/) -> bool {
        bool SF = m_selectLeptons.at(0)->isEle() == m_selectLeptons.at(1)->isEle();
        bool OS = m_selectLeptons.at(0)->q != m_selectLeptons.at(1)->q;
        return SF && OS;
    };
}

void add_DFOS_lepton_cut(Superflow* superflow) {
    *superflow << CutName("DFOS leptons") << [=](Superlink* /*sl*/) -> bool {
        bool DF = m_selectLeptons.at(0)->isEle() != m_selectLeptons.at(1)->isEle();
        bool OS = m_selectLeptons.at(0)->q != m_selectLeptons.at(1)->q;
        return DF && OS;
    };
}

int get_lepton_truth_class(Susy::Lepton* lepton) {
    if (lepton==nullptr) return -INT_MAX;

    // Get Truth information
    int T = lepton->mcType;
    int O = lepton->mcOrigin;
    int MO = lepton->mcBkgTruthOrigin;
    int MT = 0; // Not stored in SusyNt::Lepton
    int M_ID = lepton->mcBkgMotherPdgId;

    using namespace MCTruthPartClassifier;

    bool mother_is_el = fabs(M_ID) == 11;
    //bool mother_is_piZero = fabs(M_ID) == 111;
    bool bkgEl_from_phoConv = T==BkgElectron && O==PhotonConv;
    bool bkgEl_from_EMproc = T==BkgElectron && O==ElMagProc;
    bool fromSMBoson = O==WBoson || O==ZBoson || O==Higgs || O==DiBoson;
    bool MfromSMBoson = MO==WBoson || MO==ZBoson || MO==Higgs || MO==DiBoson;
    //bool noChargeFlip = M_ID*lepton->q < 0;
    //bool chargeFlip = M_ID*lepton->q > 0;

    // Defs from https://indico.cern.ch/event/725960/contributions/2987219/attachments/1641430/2621432/TruthDef_April242018.pdf
    bool promptEl1 = T==IsoElectron; //&& noChargeFlip;
    bool promptEl2 = (bkgEl_from_phoConv && mother_is_el); //&& noChargeFlip;
    bool promptEl3 = bkgEl_from_EMproc && MT==IsoElectron && (MO==top || MfromSMBoson);
    bool promptEl = promptEl1 || promptEl2 || promptEl3;

    bool promptEl_from_FSR1 = bkgEl_from_phoConv && MO==FSRPhot;
    bool promptEl_from_FSR2 = T==NonIsoPhoton && O==FSRPhot;
    bool promptEl_from_FSR = promptEl_from_FSR1 || promptEl_from_FSR2;

    //bool promptChargeFlipEl1 = T==IsoElectron && chargeFlip;
    //bool promptChargeFlipEl2 = (bkgEl_from_phoConv && mother_is_el) && chargeFlip;
    //bool promptChargeFlipEl = promptChargeFlipEl1 || promptChargeFlipEl2;

    bool promptMuon = T==IsoMuon && (O==top || fromSMBoson || O==HiggsMSSM || O==MCTruthPartClassifier::SUSY);

    bool promptPho1 = T==IsoPhoton && O==PromptPhot;
    bool promptPho2 = bkgEl_from_phoConv && MT==IsoPhoton && MO==PromptPhot;
    bool promptPho3 = bkgEl_from_EMproc  && MT==IsoPhoton && MO==PromptPhot;
    bool promptPho4 = bkgEl_from_phoConv && MT==BkgPhoton && MO==UndrPhot;
    bool promptPho5 = T==BkgPhoton && O==UndrPhot;
    bool promptPho = promptPho1 || promptPho2 || promptPho3 || promptPho4 || promptPho5;

    bool hadDecay1 = T==BkgElectron && (O==DalitzDec || O==ElMagProc || O==LightMeson || O==StrangeMeson);
    bool hadDecay2 = bkgEl_from_phoConv && MT==BkgPhoton && (MO==PiZero || MO==LightMeson || MO==StrangeMeson);
    bool hadDecay3 = bkgEl_from_EMproc && ((MT==BkgElectron && MO==StrangeMeson) || (MT==BkgPhoton && MO==PiZero));
    bool hadDecay4 = T==BkgPhoton && (O==LightMeson || O==PiZero);
    bool hadDecay5 = T==BkgMuon && (O==LightMeson || O==StrangeMeson || O==PionDecay || O==KaonDecay);
    bool hadDecay6 = T==Hadron;
    bool hadDecay = hadDecay1 || hadDecay2 || hadDecay3 || hadDecay4 || hadDecay5 || hadDecay6;

    bool Mu_as_e1 = (T==NonIsoElectron || T==NonIsoPhoton) && O==Mu;
    bool Mu_as_e2 = bkgEl_from_EMproc && MT==NonIsoElectron && MO==Mu;
    bool Mu_as_e3 = bkgEl_from_phoConv && MT==NonIsoPhoton && MO==Mu;
    bool Mu_as_e = Mu_as_e1 || Mu_as_e2 || Mu_as_e3;

    bool HF_tau1 =  (T==NonIsoElectron || T==NonIsoPhoton) && O==TauLep;
    bool HF_tau2 =  bkgEl_from_phoConv && MT==NonIsoPhoton && MO==TauLep;
    bool HF_tau3 =  T==NonIsoMuon && O==TauLep;
    bool HF_tau =  HF_tau1 || HF_tau2 || HF_tau3;

    bool HF_B1 = T==NonIsoElectron && (O==BottomMeson || O==BBbarMeson || O==BottomBaryon);
    bool HF_B2 = T==BkgPhoton && O==BottomMeson;
    bool HF_B3 = bkgEl_from_phoConv && MT==BkgPhoton && MO==BottomMeson;
    bool HF_B4 = (T==IsoMuon || T==NonIsoMuon) && (O==BottomMeson || O==BBbarMeson || O==BottomBaryon);
    bool HF_B = HF_B1 || HF_B2 || HF_B3 || HF_B4;

    bool HF_C1 = T==NonIsoElectron && (O==CharmedMeson || O==CharmedBaryon || O==CCbarMeson);
    bool HF_C2 = T==BkgElectron && O==CCbarMeson;
    bool HF_C3 = T==BkgPhoton && (O==CharmedMeson || O==CCbarMeson);
    bool HF_C4 = bkgEl_from_phoConv && MT==BkgPhoton && (MO==CharmedMeson || MO==CCbarMeson);
    bool HF_C5 = T==NonIsoMuon && (O==CharmedMeson || O==CharmedBaryon || O==CCbarMeson);
    bool HF_C6 = (T==IsoMuon || T==BkgMuon) && (O==CCbarMeson || MO==CCbarMeson);
    bool HF_C =  HF_C1 || HF_C2 || HF_C3 || HF_C4 || HF_C5 || HF_C6;


    if      (promptEl)           return 1;
    else if (promptMuon)         return 2;
    else if (promptPho)          return 3;
    else if (bkgEl_from_phoConv) return 4;
    else if (promptEl_from_FSR)  return 5;
    else if (hadDecay)           return 6;
    else if (Mu_as_e)            return 7;
    else if (HF_tau)             return 8;
    else if (HF_B)               return 9;
    else if (HF_C)               return 10;
    // else if (promptChargeFlipEl) return 2;
    else if (T && O && M_ID) {
        cout << "Unexpected Truth Class: "
             << "T = " << T << ", "
             << "O = " << O << ", "
             << "MT = " << MT << ", "
             << "MO = " << MO << ", "
             << "M_ID = " << M_ID << endl;
    }
    return -1;
}



