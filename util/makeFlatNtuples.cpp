////////////////////////////////////////////////////////////////////////////////
/// Copyright (c) <2018> by Alex Armstrong with much code borrowed from
///     Aleattin Mete (Alaettin.Serhan.Mete@cern.ch)
///
/// @file makeFlatNtuple.cpp
/// @author Alex Armstrong <alarmstr@cern.ch>
/// @date <May 27th, 2018>
/// @brief Make flat ntuples from SusyNts
///
////////////////////////////////////////////////////////////////////////////////
// TODO:
// - See the effect of adding variables on ntuple production time
// - Fix dilep trigger matching

#include "AlexAnalysisPackage/makeFlatNtuple.h"

using std::string, std::cout, std::cerr, std::vector;

////////////////////////////////////////////////////////////////////////////////
// Declarations
////////////////////////////////////////////////////////////////////////////////
// Useful macros
// TODO: Run over a vector of leptons
#define ADD_1LEP_TRIGGER_VAR(trig_name, lep) { \
  *cutflow << NewVar(#trig_name" trigger bit"); { \
      *cutflow << HFTname(#trig_name); \
      *cutflow << [&](Superlink* sl, var_bool*) -> bool { \
          if(!lep) return false;\
          bool trig_fired = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, #trig_name); \
          if (!trig_fired) return false;\
          bool trig_matched = sl->tools->triggerTool().lepton_trigger_match(lep, #trig_name);\
          return trig_matched; }; \
      *cutflow << SaveVar(); \
  } \
}
// Trig Matching for dilep triggers is buggy
// so currently not trigger matching
#define ADD_2LEP_TRIGGER_VAR(trig_name, lep0, lep1) { \
  *cutflow << NewVar(#trig_name" trigger bit"); { \
      *cutflow << HFTname(#trig_name); \
      *cutflow << [&](Superlink* sl, var_bool*) -> bool { \
          if(!lep0 || ! lep1) return false;\
          bool trig_fired = sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, #trig_name); \
          return trig_fired;\
          if (!trig_fired) return false; \
          bool trig_matched = sl->tools->triggerTool().dilepton_trigger_match(sl->nt->evt(), lep0, lep1, #trig_name);\
          return trig_matched; }; \
      *cutflow << SaveVar(); \
  } \
}
#define ADD_MULTIPLICITY_VAR(container) { \
  *cutflow << NewVar("number of "#container); { \
    *cutflow << HFTname("n_"#container); \
    *cutflow << [&](Superlink* sl, var_int*) -> int { \
      return sl->container->size(); \
    }; \
    *cutflow << SaveVar(); \
  } \
}

/// @brief Command line argument parser
struct Args {
    ////////////////////////////////////////////////////////////////////////////
    // Initialize arguments
    ////////////////////////////////////////////////////////////////////////////
    string PROG_NAME;  // always required

    // Required arguments
    string ifile_name;  ///< Input file name

    // Optional arguments with defaults
    unsigned int n_events  = -1; ///< number of events processed
    string name_suffix = ""; ///< suffix appended to output name
    bool baseline_sel = true;  ///< Apply baseline denominator selection
    bool baseline_den = false;  ///< Apply baseline denominator selection
    bool fake_num = false;  ///< Apply fake numerator selection
    bool fake_den = false;  ///< Apply fake denominator selection
    bool apply_ff = false; ///< TODO: change to apply_ff
    // TODO:
    // fake_num, fake_den, -fake_est, -base_sel
    // Produces ntuples for up to 4 options
    // no option -> base_sel, fake option -> base_sel = false unless explicitely requested

    ////////////////////////////////////////////////////////////////////////////
    // Print information
    // TODO: remove single char option for fake stuff
    ////////////////////////////////////////////////////////////////////////////
    void print_usage() const {
        printf("===========================================================\n");
        printf(" %s\n", PROG_NAME.c_str());
        printf(" TODO : Brief description of executable\n");
        printf("===========================================================\n");
        printf("Required Parameters:\n");
        printf("\t-i, --input     Input file name\n");
        printf("\nOptional Parameters:\n");
        printf("\t-n, --nevents   number of events to process\n");
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
        printf("\tInput file name : %s\n", ifile_name.c_str());
        printf("\tnevents         : %i\n", n_events);
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
                ifile_name = arg_value;
            } else if (arg == "-n" || arg == "--nevents") {
                n_events = arg_value;
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
        bool denominator_selection = baseline_den || fake_den;
        bool numerator_selection = baseline_sel || fake_num;

        if (ifile_name.size() == 0) {
            cerr << "ERROR :: No input source given\n";
            return false;
        } else if (apply_ff && !denominator_selection) {
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

    cout << "TESTING :: Stop here to check" << std::endl;
    return 0;
    ////////////////////////////////////////////////////////////////////////////
    // Main implementation
    ////////////////////////////////////////////////////////////////////////////
    // Build list of cutflows to run
    setup_chain(m_chain, args.input_name);

    if (args.apply_ff) initialize_fake_factor_tool(m_applyFakeFactorTool);


    vector<Superflow*> superflows;
    if (base_sel) {
        superflows.push_back(get_cutflow(m_chain, BASE_SEL))
    } else if (args.fake_num) {
        superflows.push_back(get_cutflow(m_chain, FAKE_NUM))
    } else if (args.fake_den) {
        superflows.push_back(get_cutflow(m_chain, FAKE_DEN))
    } else if (args.base_den) {
        superflows.push_back(get_cutflow(m_chain, BASE_DEN))
    } else {
        cout << "ERROR :: No configuration for that region yet\n";
        //TODO: ttbar, zll, W+jets, Diboson, Ztautau
    }

    // Run the chosen cutflows
    for (uint ii = 0; ii < superflows.size(); i++) {
        printf("\n\n [%d/%d] Running cutflow \n\n",ii+1, superflows.size())
        Superflow *sf = superflows.at(ii)
        chain->Process(sf, input_name, n_events, n_skip_events);
        delete sf;
        delete chain;
    }
    superflows.clear()
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
        printf("ERROR (initialize_chain) :: Unrecognized input %s", iname);
        return 0;
    }
}

ApplyFakeFactor* initialize_fake_factor_tool(ApplyFakeFactor* applyFakeFactorTool) {
    //////////////////////////////////////////////////////////////////////////////
    // Initialize relevant fake variables. Set inside cutflow
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
    add_pre_cuts(superflow);
    add_cleaing_cuts(superflow);
    if (sel_type = BASELINE) {
        add_baseline_lepton_cuts(superflow);
    } else if (self_type = BASE_DEN) {
        add_baseline_den_lepton_cuts(superflow);
    } else if (sel_type = FAKE_NUM) {
        add_fake_num_lepton_cuts(superflow);
    } else if (sel_type = FAKE_DEN) {
        add_fake_den_lepton_cuts(superflow);

    add_final_cuts(superflow, sel_type);

    ////////////////////////////////////////////////////////////////////////////
    // Add superflow variables
    add_event_variables(superflow);
    add_trigger_variables(superflow);
    add_lepton_variables(superflow);

    add_met_variables(superflow);
    add_jet_variables(superflow);
    if (m_do_fakes) {
        add_fake_variables(superflow, sel_type);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Reset helpful variables
    add_shortcut_variables_reset(superflow);

    return superflow;
}

////////////////////////////////////////////////////////////////////////////////
// TODO: Switch to using references instead of pointers
string determine_suffix(string user_suffix, Sel sel_type, bool apply_ff) {
    string full_suffix = user_suffix;

    if (user_suffix != "") full_suffix += "_";

    if (apply_ff) full_suffix += 'fakes_';

    if (sel_type = FAKE_NUM) full_suffix += "zjets_num";
    else if (sel_type = FAKE_DEN) full_suffix += "zjets_den";
    else if (sel_type = BASE_SEL) full_suffix +=
    else if (self_type = BASE_DEN) full_suffix += "baseline_den";

    return full_suffix
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
void add_shortcut_variables(Superflow* superflow, Sel sel_type) {
    *cutflow << CutName("(No Cuts) Define shortcut variables")
             << [&](Superlink* sl, var_void*) -> bool {
        m_cutflags = sl->nt->evt()->cutFlags[NtSys::NOM];

        ////////////////////////////////////////////////////////////////////////
        // Lepton shortcuts
        m_preLeptons = *sl->preLeptons;
        m_baseLeptons =  *sl->baseLeptons;
        m_signalLeptons = *sl->leptons;

        m_preElectrons = *sl->preElectrons;
        m_baseElectrons =  *sl->baseElectrons;
        m_signalElectrons = *sl->electrons;

        m_preMuons = *sl->preMuons;
        m_baseMuons =  *sl->baseMuons;
        m_signalMuons = *sl->muons;

        m_preTaus = *sl->preTaus;
        m_baseTaus =  *sl->baseTaus;
        m_signalTaus = *sl->taus;

        ////////////////////////////////////////////////////////////////////////
        // MET shortcuts
        m_met = *sl->met;

        ////////////////////////////////////////////////////////////////////////
        // Fake shortcuts
        add_fake_shortcut_variables(superflow, sel_type);

        if (m_signalLeptons.size() >= 1) m_lepton0 = *m_signalLeptons.at(0);
        if (m_signalLeptons.size() >= 2) m_lepton1 = *m_signalLeptons.at(1);
        m_dileptonP4 = m_lepton0 + m_lepton1;
        m_el0 = m_signalElectrons.size() >= 1 ? m_signalElectrons.at(0) : nullptr;
        m_el1 = m_signalElectrons.size() >= 2 ? m_signalElectrons.at(1) : nullptr;
        m_mu0 = m_signalMuons.size() >= 1 ? m_signalMuons.at(0) : nullptr;
        m_mu1 = m_signalMuons.size() >= 2 ? m_signalMuons.at(1) : nullptr;

        // Pick the right set and ordingering of leptons for the given selection
        // type. Affects trigger matching, dilepton requirements, etc...
        if (sel_type == FAKE_DEN || sel_type == FAKE_NUM) {
            m_selectLeptons = m_Zlep;
        } else {
            m_selectLeptons = m_signalLeptons;
        }

        ////////////////////////////////////////////////////////////////////////
        // Jet shortcuts
        m_preJets = *sl->preJets;
        m_baseJets = *sl->baseJets;
        m_signalJets = *sl->jets;

        if (m_baseJets.size() > 0) {
            m_Jet0 = *m_baseJets.at(0);
            if (m_baseJets.size() > 1) {
               m_Jet1 = *m_baseJets.at(1);
               m_JetP4 = m_Jet0 + m_Jet1;
            }
        }
        // TODO: replace auto with explicit class
        for (auto& jet : m_baseJets) {
              if (sl->tools->m_jetSelector->isLight(jet))  {
                   m_lightJets.push_back(jet);
              } else if (sl->tools->m_jetSelector->isB(jet)) {
                  m_BJets.push_back(jet);
              } else if (sl->tools->m_jetSelector->isForward(jet))  {
                  m_forwardJets.push_back(jet);
              }
        }
        std::sort(m_lightJets.begin()  , m_lightJets.end()  , comparePt);
        std::sort(m_BJets.begin()      , m_BJets.end()      , comparePt);
        std::sort(m_forwardJets.begin(), m_forwardJets.end(), comparePt);

        return true; // All events pass this cut
    }
}

void add_pre_cuts(Superflow* superflow) {
  // All events before cuts
  *cutflow << CutName("read in") << [](Superlink* /*sl*/) -> bool { return true; };

  // xTauFW Cut
  if (!do_fakes_zjets) {
    *cutflow << CutName("xTau: 2+ Loose Leptons") << [&](Superlink* sl) -> bool {
        uint nLooseLeptons = 0;
        for (const auto* mu : m_preMuons) {if (mu->loose) nLooseLeptons++;}
        for (const auto* ele : m_preElectrons) {if (ele->looseLLH) nLooseLeptons++;}
        return nLooseLeptons >= 2;
    };
  }
}
void add_cleaing_cuts(Superflow* superflow) {
  *cutflow << CutName("Pass GRL") << [&](Superlink* sl) -> bool {
      return (sl->tools->passGRL(m_cutflags));
  };
  *cutflow << CutName("LAr error") << [&](Superlink* sl) -> bool {
      return (sl->tools->passLarErr(m_cutflags));
  };
  *cutflow << CutName("Tile error") << [&](Superlink* sl) -> bool {
      return (sl->tools->passTileErr(m_cutflags));
  };
  *cutflow << CutName("TTC veto") << [&](Superlink* sl) -> bool {
      return (sl->tools->passTTC(m_cutflags));
  };
  *cutflow << CutName("SCT err") << [&](Superlink* sl) -> bool {
      return (sl->tools->passSCTErr(m_cutflags));
  };
}

void add_baseline_lepton_cuts(Superflow* superflow) {
    *cutflow << CutName("nBaselineLep = nSignalLep") << [](Superlink* sl) -> bool {
        return (m_signalLeptons.size() == m_baseLeptons->size());
    };
    *cutflow << CutName("2+ leptons") << [](Superlink* sl) -> bool {
        return (m_signalLeptons.size() >= 2);
    };
    add_DFOS_lepton_cut(superflow);
}
void add_baseline_den_lepton_cuts(Superflow* superflow) {
    *cutflow << CutName("1-ID Lepton and 1 Anti-ID Lepton") << [&](Superlink* /*sl*/) -> bool {
        return (m_lepID_n == 1 && m_lepAntiID_n == 1);
    };
    add_DFOS_lepton_cut(superflow);
}
void add_fake_num_lepton_cuts(Superflow* superflow) {
    *cutflow << CutName("3-ID Leptons") << [&](Superlink* /*sl*/) -> bool {
        return (m_lepID_n == 3);
    };
    add_SFOS_lepton_cut(superflow);
}
void add_fake_den_lepton_cuts(Superflow* superflow) {
    *cutflow << CutName("2-ID Leptons and 1+ Anti-ID Lepton") << [&](Superlink* /*sl*/) -> bool {
        return (m_lepID_n == 2 && m_lepAntiID_n >= 1);
    };
    add_SFOS_lepton_cut(superflow);
}

void add_final_cuts(Superflow* superflow, Sel sel_type) {
    *cutflow << CutName("pass good vertex") << [&](Superlink* sl) -> bool {
        return (sl->tools->passGoodVtx(cutflags));
    };
    *cutflow << CutName("jet cleaning") << [&](Superlink* sl) -> bool {
        return (sl->tools->passJetCleaning(&m_baseJets));
    };
    if (sel_type == BASELINE) {
      *cutflow << CutName("2 leptons") << [](Superlink* sl) -> bool {
          return (m_signalLeptons.size() == 2);
      };
    }
    *cutflow << CutName("Tau veto") << [](Superlink* sl) -> bool {
        return (m_signalTaus.size() == 0);
    };
}

void add_event_variables(Superflow* superflow, Sel sel_type) {
  *cutflow << NewVar("Event run number"); {
    *cutflow << HFTname("RunNumber");
    *cutflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->run; };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Event number"); {
    *cutflow << HFTname("event_number");
    *cutflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->eventNumber; };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("is Monte Carlo"); {
    *cutflow << HFTname("isMC");
    *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->nt->evt()->isMC ? true : false; };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("event weight"); {
    *cutflow << HFTname("eventweight");
    *cutflow << [&](Superlink* sl, var_double*) -> double {
        // TODO: add fake factor tool to superflow
        float fakeFactor = 1;
        if (args.apply_ff && m_denominator_selection) {
            bool probeIsEl = m_baseLeptons.at(m_antiID_idx0)->isEle();
            LepEnum::LepType typeOfLep = probeIsEl ? LepEnum::Electron : LepEnum::Muon;
            fakeFactor = m_applyFakeFactorTool->apply(m_LepFake0.Pt(), typeOfLep);
        }
        return sl->weights->product() * sl->nt->evt()->wPileup * fakeFactor;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("sample DSID"); {
    *cutflow << HFTname("dsid");
    *cutflow << [](Superlink* sl, var_int*) -> int {return sl->nt->evt()->mcChannel;};
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("treatAsYear"); {
      // 15/16 Year ID
      *cutflow << HFTname("treatAsYear");
      *cutflow << [](Superlink* sl, var_double*) -> int { return sl->nt->evt()->treatAsYear; };
      *cutflow << SaveVar();
  }
}
void add_trigger_variables(Superflow* superflow, Sel sel_type) {
  //////////////////////////////////////////////////////////////////////////////
  // Trigger Variables
  // ADD_*_TRIGGER_VAR preprocessor defined

  //////////////////////////////////////////////////////////////////////////////
  // 2015

  // Dilepton Triggers
  ADD_2LEP_TRIGGER_VAR(HLT_e17_lhloose_mu14, m_el0, m_mu0)
  ADD_2LEP_TRIGGER_VAR(HLT_e24_lhmedium_L1EM20VHI_mu8noL1, m_el0, m_mu0)
  ADD_2LEP_TRIGGER_VAR(HLT_e7_lhmedium_mu24, m_mu0, m_el0)
  ADD_2LEP_TRIGGER_VAR(HLT_2e12_lhloose_L12EM10VH, m_el0, m_el1)
  // TODO: Add to SusyNts HLT_2mu10)

  // Single Electron Triggers
  // TODO: Trigger match a provided vector of electrons instead of just one
  ADD_1LEP_TRIGGER_VAR(HLT_e24_lhmedium_L1EM20VH, m_el0)
  ADD_1LEP_TRIGGER_VAR(HLT_e60_lhmedium, m_el0)
  ADD_1LEP_TRIGGER_VAR(HLT_e120_lhloose, m_el0)

  // Single Muon Triggers
  ADD_1LEP_TRIGGER_VAR(HLT_mu20_iloose_L1MU15, m_mu0)
  ADD_1LEP_TRIGGER_VAR(HLT_mu40, m_mu0)

  //////////////////////////////////////////////////////////////////////////////
  // 2016

  // Dilepton Triggers
  ADD_2LEP_TRIGGER_VAR(HLT_e17_lhloose_nod0_mu14, m_el0, m_mu0)
  // TODO: Add to SusyNts HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1, m_el0, m_mu0)
  ADD_2LEP_TRIGGER_VAR(HLT_e7_lhmedium_nod0_mu24, m_mu0, m_el0)
  ADD_2LEP_TRIGGER_VAR(HLT_2e17_lhvloose_nod0, m_el0, m_el1)
  // TODO: Add to SusyNts HLT_2mu14)

  // Single Electron Triggers
  ADD_1LEP_TRIGGER_VAR(HLT_e26_lhtight_nod0_ivarloose, m_el0)
  ADD_1LEP_TRIGGER_VAR(HLT_e60_lhmedium_nod0, m_el0)
  ADD_1LEP_TRIGGER_VAR(HLT_e140_lhloose_nod0, m_el0)

  // Single Muon Triggers
  ADD_1LEP_TRIGGER_VAR(HLT_mu26_ivarmedium, m_mu0)
  ADD_1LEP_TRIGGER_VAR(HLT_mu50, m_mu0)
}

void add_lepton_variables(Superflow* superflow, Sel sel_type) {
    //////////////////////////////////////////////////////////////////////////////
    // Lepton informtion
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // Get multiplicity variables
    // ADD_MULTIPLICITY_VAR preprocessor defined
    ADD_MULTIPLICITY_VAR(preLeptons)
    ADD_MULTIPLICITY_VAR(baseLeptons)
    ADD_MULTIPLICITY_VAR(leptons)
    ADD_MULTIPLICITY_VAR(preElectrons)
    ADD_MULTIPLICITY_VAR(baseElectrons)
    ADD_MULTIPLICITY_VAR(electrons)
    ADD_MULTIPLICITY_VAR(preMuons)
    ADD_MULTIPLICITY_VAR(baseMuons)
    ADD_MULTIPLICITY_VAR(muons)
    ADD_MULTIPLICITY_VAR(preTaus)
    ADD_MULTIPLICITY_VAR(baseTaus)
    ADD_MULTIPLICITY_VAR(taus)

    //////////////////////////////////////////////////////////////////////////////
    // PreElectrons
    *cutflow << NewVar("preElectron E_caloClus"); {
      *cutflow << HFTname("preEl_EcaloClus");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : m_preElectrons) {
          out.push_back(el->clusE);
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("preElectron eta_caloClus"); {
      *cutflow << HFTname("preEl_etaCaloClus");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : m_preElectrons) {
          out.push_back(el->clusEta);
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("preElectron Et"); {
      *cutflow << HFTname("preEl_Et");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : m_preElectrons) {
          float E_caloClus = el->clusE;
          float eta_caloClus = el->clusEta;
          out.push_back(E_caloClus / std::cosh(eta_caloClus));
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("preElectron pt"); {
      *cutflow << HFTname("preEl_pt");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : m_preElectrons) {out.push_back(el->Pt()); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("preElectron clusterEtaBE"); {
      *cutflow << HFTname("preEl_clusEtaBE");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : m_preElectrons) {out.push_back(el->clusEtaBE);}
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("preElectron eta"); {
      *cutflow << HFTname("preEl_eta");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : m_preElectrons) {out.push_back(el->Eta());}
        return out;
      };
      *cutflow << SaveVar();
    }

    //////////////////////////////////////////////////////////////////////////////
    // Baseline Electrons
    *cutflow << NewVar("Baseline Electron etconetopo20"); {
      *cutflow << HFTname("baseEl_etconetopo20");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : m_baseElectrons) {out.push_back(el->etconetopo20); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Baseline Electron ptvarcone20"); {
      *cutflow << HFTname("baseEl_ptvarcone20");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : m_baseElectrons) {out.push_back(el->ptvarcone20); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Baseline Electron ID (non-inclusive)"); {
      *cutflow << HFTname("baseEl_ID");
      *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for (auto& el : m_baseElectrons) {
          if (el->tightLLH) out.push_back(0);
          else if (el->mediumLLH) out.push_back(1);
          else if (el->looseLLHBLayer) out.push_back(3);
          else if (el->looseLLH) out.push_back(2);
          else if (el->veryLooseLLH) out.push_back(4);
          else out.push_back(5);
        }
        return out;
      };
      *cutflow << SaveVar();
    }

    //////////////////////////////////////////////////////////////////////////////
    // Signal Electrons
    *cutflow << NewVar("Electron type"); {
      *cutflow << HFTname("el_type");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : m_signalElectrons) {out.push_back(el->mcType); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Electron origin"); {
      *cutflow << HFTname("el_origin");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : m_signalElectrons) {out.push_back(el->mcOrigin); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Electron ID (non-inclusive)"); {
      *cutflow << HFTname("El_ID");
      *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for (auto& el : m_signalElectrons) {
          if (el->tightLLH) out.push_back(0);
          else if (el->mediumLLH) out.push_back(1);
          else if (el->looseLLHBLayer) out.push_back(3);
          else if (el->looseLLH) out.push_back(2);
          else if (el->veryLooseLLH) out.push_back(4);
          else out.push_back(5);
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("dRy(sigEl, preMu)"); {
      *cutflow << HFTname("dRy_sEl_pMu");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : m_signalElectrons) {
          for (auto& mu : m_preMuons) {
            out.push_back(el->DeltaRy(*mu));
          }
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("dRy(sigEl, preMu_shrdTrk)"); {
      *cutflow << HFTname("dRy_sEl_pMu_shrdtrk");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : m_signalElectrons) {
          for (auto& mu : m_preMuons) {
            if(!(el->sharedMuTrk.size()>0)) continue;
            if(mu->idx > ((int)el->sharedMuTrk.size()-1)) continue;
            if(el->sharedMuTrk[mu->idx]!=1) continue;
            out.push_back(el->DeltaRy(*mu));
          }
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("dRy(sigEl, baseMu not Calo)"); {
      *cutflow << HFTname("dRy_sEl_bMu_noCalo");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : m_signalElectrons) {
          for (auto& mu : m_preMuons) {
            if (!cutflow->nttools().m_muonSelector->isBaseline(mu)) continue;
            if (cutflow->nttools().isSignal(mu)) continue;
            if (mu->isCaloTagged) continue;
            out.push_back(el->DeltaRy(*mu));
          }
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("dRy(sigEl, baseMu CaloTag)"); {
      *cutflow << HFTname("dRy_sEl_bMu_Calo");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : m_signalElectrons) {
          for (auto& mu : m_preMuons) {
            if (!cutflow->nttools().m_muonSelector->isBaseline(mu)) continue;
            if (cutflow->nttools().isSignal(mu)) continue;
            if (!mu->isCaloTagged) continue;
            out.push_back(el->DeltaRy(*mu));
          }
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Electron d0sigBSCorr"); {
      *cutflow << HFTname("el_d0sigBSCorr");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : m_signalElectrons) {out.push_back(el->d0sigBSCorr); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Electron z0SinTheta"); {
      *cutflow << HFTname("el_z0SinTheta");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& el : m_signalElectrons) {out.push_back(el->z0SinTheta()); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("subleading electron track pt"); {
      *cutflow << HFTname("el1_track_pt");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() > 1 && m_signalLeptons.at(1)->isEle()) {
            const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(m_signalLeptons.at(1));
            return ele ? ele->trackPt : -DBL_MAX;
        } else {
            return -DBL_MAX;
        }
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("subleading electron clus pt"); {
      *cutflow << HFTname("el1_clus_pt");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() > 1 && m_signalLeptons.at(1)->isEle()) {
            const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(m_signalLeptons.at(1));
            return ele ? ele->clusE : -DBL_MAX;
        } else {
            return -DBL_MAX;
        }
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Subleading Electron pT track-cluster ratio"); {
      *cutflow << HFTname("el1pT_trackclus_ratio");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() > 1 && m_signalLeptons.at(1)->isEle()) {
            const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(m_signalLeptons.at(1));
            return ele ? ele->trackPt / ele->clusE : -DBL_MAX;
        } else {
            return -DBL_MAX;
        }
      };
      *cutflow << SaveVar();
    }
    //////////////////////////////////////////////////////////////////////////////
    // PreMuons

    *cutflow << NewVar("preMuon pt"); {
      *cutflow << HFTname("preMu_pt");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& mu : m_preMuons) {out.push_back(mu->Pt()); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("isCaloTagged"); {
      *cutflow << HFTname("isCaloTagged");
      *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for (auto& mu : m_preMuons) { out.push_back(mu->isCaloTagged);}
        return out;
      };
      *cutflow << SaveVar();
    }

    *cutflow << NewVar("PreMuon ID (non-inclusive)"); {
      *cutflow << HFTname("preMu_ID");
      *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for (auto& mu : m_preMuons) {
          if (mu->tight) out.push_back(0);
          else if (mu->medium) out.push_back(1);
          else if (mu->loose) out.push_back(2);
          else if (mu->veryLoose) out.push_back(3);
          else out.push_back(4);
        }
        return out;
      };
      *cutflow << SaveVar();
    }

    //////////////////////////////////////////////////////////////////////////////
    // Baseline Muons
    *cutflow << NewVar("baseline Muon pt"); {
      *cutflow << HFTname("baseMu_pt");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& mu : m_baseMuons) {out.push_back(mu->Pt()); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Baseline Muon eta"); {
      *cutflow << HFTname("baseMu_eta");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& mu : m_baseMuons) {out.push_back(mu->Eta()); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Baseline Muon etconetopo20"); {
      *cutflow << HFTname("baseMu_etconetopo20");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& mu : m_baseMuons) {out.push_back(mu->etconetopo20); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Baseline Muon ptvarcone30"); {
      *cutflow << HFTname("baseMu_ptvarcone30");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& mu : m_baseMuons) {out.push_back(mu->ptvarcone30); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Baseline Muon ID (non-inclusive)"); {
      *cutflow << HFTname("baseMu_ID");
      *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for (auto& mu : m_baseMuons) {
          if (mu->tight) out.push_back(0);
          else if (mu->medium) out.push_back(1);
          else if (mu->loose) out.push_back(2);
          else if (mu->veryLoose) out.push_back(3);
          else out.push_back(4);
        }
        return out;
      };
      *cutflow << SaveVar();
    }

    *cutflow << NewVar("Muon ID (non-inclusive)"); {
      *cutflow << HFTname("Mu_ID");
      *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for (auto& mu : m_signalMuons) {
          if (!mu->veryLoose) out.push_back(4);
          else if (mu->veryLoose && !mu->loose) out.push_back(3);
          else if (mu->loose && !mu->medium) out.push_back(2);
          else if (mu->medium && !mu->tight) out.push_back(1);
          else if (mu->tight) out.push_back(0);
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    //////////////////////////////////////////////////////////////////////////////
    // Signal Muons
    *cutflow << NewVar("Muon type"); {
      *cutflow << HFTname("mu_type");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& mu : m_signalMuons) {out.push_back(mu->mcType); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Muon origin"); {
      *cutflow << HFTname("mu_origin");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& mu : m_signalMuons) {out.push_back(mu->mcOrigin); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Muon d0sigBSCorr"); {
      *cutflow << HFTname("mu_d0sigBSCorr");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& mu : m_signalMuons) {out.push_back(mu->d0sigBSCorr); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Muon z0SinTheta"); {
      *cutflow << HFTname("mu_z0SinTheta");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& mu : m_signalMuons) {out.push_back(mu->z0SinTheta()); }
        return out;
      };
      *cutflow << SaveVar();
    }

    //////////////////////////////////////////////////////////////////////////////
    // PreTaus
    *cutflow << NewVar("preTau charge"); {
      *cutflow << HFTname("preTau_q");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& tau : m_preTaus) {out.push_back(tau->q); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("preTau nTracks"); {
      *cutflow << HFTname("preTau_nTracks");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& tau : m_preTaus) {out.push_back(tau->nTrack); }
        return out;
      };
      *cutflow << SaveVar();
    }

    //////////////////////////////////////////////////////////////////////////////
    // Baseline Taus
    *cutflow << NewVar("Baseline Tau pT"); {
      *cutflow << HFTname("baseTau_pT");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& tau : m_preTaus) {out.push_back(tau->Pt()); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Baseline Tau eta"); {
      *cutflow << HFTname("baseTau_eta");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& tau : m_preTaus) {out.push_back(tau->Eta()); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Baseline Tau nTracks"); {
      *cutflow << HFTname("baseTau_nTracks");
      *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for (auto& tau : m_baseTaus) {out.push_back(tau->nTrack); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Baseline Tau ID (non-inclusive)"); {
      *cutflow << HFTname("baseTau_ID");
      *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for (auto& tau : m_baseTaus) {
          if (!tau->loose) out.push_back(3);
          else if (tau->loose && !tau->medium) out.push_back(2);
          else if (tau->medium && !tau->tight) out.push_back(1);
          else if (tau->tight) out.push_back(0);
        }
        return out;
      };
      *cutflow << SaveVar();
    }

    //////////////////////////////////////////////////////////////////////////////
    // Signal Leptons
    *cutflow << NewVar("Lepton Iso (non-inclusive)"); {
      *cutflow << HFTname("Lep_Iso");
      *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for (auto& lep : m_signalLeptons) {
          bool flag = false;
          out.push_back(-1); // for tracking all entries and normalizing bins
          if (lep->isoGradient){flag=true, out.push_back(0);}
          if (lep->isoGradientLoose){flag=true, out.push_back(1);}
          if (lep->isoLoose){flag=true, out.push_back(2);}
          if (lep->isoLooseTrackOnly){flag=true, out.push_back(3);}
          if (lep->isoFixedCutTightTrackOnly){flag=true; out.push_back(4);}
          if (!flag) out.push_back(5);
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton pt"); {
      *cutflow << HFTname("l_pt");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->Pt());
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton eta"); {
      *cutflow << HFTname("l_eta");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->Eta());
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton phi"); {
      *cutflow << HFTname("l_phi");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->Phi());
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton flavor (0: e, 1: m)"); {
      *cutflow << HFTname("l_flav");
      *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->isEle() ? 0 : 1);
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton type"); {
      *cutflow << HFTname("l_type");
      *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->mcType);
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton origin"); {
      *cutflow << HFTname("l_origin");
      *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->mcOrigin);
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton charge"); {
      *cutflow << HFTname("l_q");
      *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->q);
        }
        return out;
        };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton BkgMotherPdgId"); {
      *cutflow << HFTname("l_BkgMotherPdgId");
      *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->mcBkgMotherPdgId);
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton BkgTruthOrigin"); {
      *cutflow << HFTname("l_BkgTruthOrigin");
      *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->mcBkgTruthOrigin);
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton matched2TruthLepton"); {
      *cutflow << HFTname("l_matched2TruthLepton");
      *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for(auto& lepton : m_selectLeptons) {
          if (lepton) out.push_back(lepton->matched2TruthLepton);
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("lepton classification"); {
      *cutflow << HFTname("l_truthClass");
      *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for(auto& lepton : m_selectLeptons) {
          if (!lepton) continue;
          // Get Truth information
          int T = lepton->mcType;
          int O = lepton->mcOrigin;
          int MO = lepton->mcBkgTruthOrigin;
          int MT = 0; // Not stored in SusyNt::Lepton
          int M_ID = lepton->mcBkgMotherPdgId;

          uint lep_class = 0;
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
          //TODO: Remove photon conversion prompt electrons
          bool promptEl1 = T==IsoElectron; //&& noChargeFlip;
          bool promptEl2 = (bkgEl_from_phoConv && mother_is_el); //&& noChargeFlip;
          bool promptEl3 = bkgEl_from_EMproc && MT==IsoElectron && (MO==top || MfromSMBoson);
          bool promptEl4 = bkgEl_from_phoConv && MO==FSRPhot;
          bool promptEl5 = T==NonIsoPhoton && O==FSRPhot;
          bool promptEl = promptEl1 || promptEl2 || promptEl3 || promptEl4 || promptEl5;

          //bool promptChargeFlipEl1 = T==IsoElectron && chargeFlip;
          //bool promptChargeFlipEl2 = (bkgEl_from_phoConv && mother_is_el) && chargeFlip;
          //bool promptChargeFlipEl = promptChargeFlipEl1 || promptChargeFlipEl2;

          bool promptMuon = T==IsoMuon && (
              O==top || fromSMBoson || O==HiggsMSSM || O==MCTruthPartClassifier::SUSY);

          bool promptPho1 = T==IsoPhoton && O==PromptPhot;
          bool promptPho2 = bkgEl_from_phoConv && MT==IsoPhoton && MO==PromptPhot;
          bool promptPho3 = bkgEl_from_EMproc  && MT==IsoPhoton && MO==PromptPhot;
          bool promptPho4 = bkgEl_from_phoConv && MT==BkgPhoton && MO==UndrPhot;
          bool promptPho5 = T==BkgPhoton && O==UndrPhot;
          bool promptPho = promptPho1 || promptPho2 || promptPho3 || promptPho4 || promptPho5;

          bool hadDecay1 = T==BkgElectron && (
              O==DalitzDec || O==ElMagProc || O==LightMeson || O==StrangeMeson);
          bool hadDecay2 = bkgEl_from_phoConv && MT==BkgPhoton && (
              MO==PiZero || MO==LightMeson || MO==StrangeMeson);
          bool hadDecay3 = bkgEl_from_EMproc &&
              ((MT==BkgElectron && MO==StrangeMeson) || (MT==BkgPhoton && MO==PiZero));
          bool hadDecay4 = T==BkgPhoton && (O==LightMeson || O==PiZero);
          bool hadDecay5 = T==BkgMuon && (
              O==LightMeson || O==StrangeMeson || O==PionDecay || O==KaonDecay);
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

          if (promptEl) lep_class = 1;
         // else if (promptChargeFlipEl) lep_class = 2;
          else if (promptMuon) lep_class = 2;
          else if (promptPho) lep_class = 3;
          else if (hadDecay) lep_class = 4;
          // Stand-in while Mother type is not available
          else if (Mu_as_e) lep_class = 5;
          else if (HF_tau) lep_class = 6;
          else if (HF_B) lep_class = 7;
          else if (HF_C) lep_class = 8;
          else if (bkgEl_from_phoConv) lep_class = -1;
          else if (T && O && M_ID) {
              cout << "Unexpected Truth Class: "
                   << "T = " << T << ", "
                   << "O = " << O << ", "
                   << "MT = " << MT << ", "
                   << "MO = " << MO << ", "
                   << "M_ID = " << M_ID << endl;
          }
          out.push_back(lep_class);
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    // xTauFW variable
    *cutflow << NewVar("leptons sign product"); {
      *cutflow << HFTname("LepLepSign");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if (m_signalLeptons.size() < 2) return -INT_MAX;
        return m_signalLeptons.at(0)->q*m_signalLeptons.at(1)->q;
      };
      *cutflow << SaveVar();
    }

    //////////////////////////////////////////////////////////////////////////////
    // Two-lepton properties

    // Lepton variables
    // xTauFW variable
    *cutflow << NewVar("Leading lepton pT"); {
      *cutflow << HFTname("Lep0Pt");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() < 1) return -DBL_MAX;
        return m_lepton0.Pt(); };
      *cutflow << SaveVar();
    }
    // xTauFW variable
    *cutflow << NewVar("Subleading lepton pT"); {
      *cutflow << HFTname("Lep1Pt");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() < 2) return -DBL_MAX;
        return m_lepton1.Pt(); };
      *cutflow << SaveVar();
    }
    // xTauFW variable
    *cutflow << NewVar("Leading lepton eta"); {
      *cutflow << HFTname("Lep0Eta");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() < 1) return -DBL_MAX;
        return m_lepton0.Eta(); };
      *cutflow << SaveVar();
    }
    // xTauFW variable
    *cutflow << NewVar("Subleading lepton eta"); {
      *cutflow << HFTname("Lep1Eta");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() < 2) return -DBL_MAX;
        return m_lepton1.Eta(); };
      *cutflow << SaveVar();
    }
    // xTauFW variable
    *cutflow << NewVar("Leading lepton phi"); {
      *cutflow << HFTname("Lep0Phi");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() < 1) return -DBL_MAX;
        return m_lepton0.Phi(); };
      *cutflow << SaveVar();
    }
    // xTauFW variable
    *cutflow << NewVar("Subleading lepton phi"); {
      *cutflow << HFTname("Lep1Phi");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() < 2) return -DBL_MAX;
        return m_lepton1.Phi(); };
      *cutflow << SaveVar();
    }

    *cutflow << NewVar("Leading lepton invariant mass"); {
      *cutflow << HFTname("MLep0");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() < 1) return -DBL_MAX;
        return m_lepton0.M(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Subleading lepton invariant mass"); {
      *cutflow << HFTname("MLep1");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() < 2) return -DBL_MAX;
        return m_lepton1.M(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Delta eta between leptons"); {
      *cutflow << HFTname("DEtaLL");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() < 2) return -DBL_MAX;
        return fabs(m_lepton0.Eta() - m_lepton1.Eta()); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Delta phi between leptons"); {
      *cutflow << HFTname("DphiLL");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() < 2) return -DBL_MAX;
        return fabs(m_lepton0.DeltaPhi(m_lepton1));};
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta R of di-lepton system"); {
      *cutflow << HFTname("drll");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() < 2) return -DBL_MAX;
        return m_lepton0.DeltaR(m_lepton1); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("dilepton flavor"); {
      *cutflow << HFTname("dilep_flav");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if(m_signalLeptons.size()<2) return -INT_MAX;
        if(m_signalLeptons.at(0)->isEle() && m_signalLeptons.at(1)->isMu()){return 0;}       // e mu  case
        else if(m_signalLeptons.at(0)->isMu() && m_signalLeptons.at(1)->isEle()){return 1;}  // mu e  case
        else if(m_signalLeptons.at(0)->isEle() && m_signalLeptons.at(1)->isEle()){return 2;} // e e   case
        else if(m_signalLeptons.at(0)->isMu() && m_signalLeptons.at(1)->isMu()){return 3;}   // mu mu case
        return 4;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("dilepton flavor is e mu"); {
      *cutflow << HFTname("isEM");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
          if(m_signalLeptons.size()<2) return -INT_MAX;
        return m_signalLeptons.at(0)->isEle() && m_signalLeptons.at(1)->isMu();  // e mu  case
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("dilepton flavor is mu e"); {
      *cutflow << HFTname("isME");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
          if(m_signalLeptons.size()<2) return -INT_MAX;
          return m_signalLeptons.at(0)->isMu() && m_signalLeptons.at(1)->isEle();  // mu e  case
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("collinear mass, M_coll"); {
      *cutflow << HFTname("MCollASym");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() < 2) return -DBL_MAX;
        double deta = fabs(m_lepton0.Eta()-m_lepton1.Eta());
          double dphi = m_lepton0.DeltaPhi(m_lepton1);
          return sqrt(2*m_lepton0.Pt()*(m_lepton1.Pt()+m_met.Et)*(cosh(deta)-cos(dphi)));
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("mass of di-lepton system, M_ll"); {
      *cutflow << HFTname("MLL");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() < 2) return -DBL_MAX;
        return m_dileptonP4.M(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Pt of di-lepton system, Pt_ll"); {
      *cutflow << HFTname("ptll");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() < 2) return -DBL_MAX;
        return m_dileptonP4.Pt(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("diff_pt between leading and sub-leading lepton"); {
      *cutflow << HFTname("dpt_ll");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
          if (m_signalLeptons.size() < 2) return -DBL_MAX;
          return m_lepton0.Pt() - m_lepton1.Pt();
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("dphi between leading lepton and MET"); {
      *cutflow << HFTname("DphiLep0MET");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
          if (m_signalLeptons.size() < 1) return -DBL_MAX;
          return m_lepton0.DeltaPhi(MET);
      };
      *cutflow << SaveVar();
    }

    //////////////////////////////////////////////////////////////////////////////
    // Optimized Region Cuts
    *cutflow << NewVar("dphi between sub-leading lepton and MET"); {
      *cutflow << HFTname("DphiLep1MET");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
          if (m_signalLeptons.size() < 2) return -DBL_MAX;
          return m_lepton1.DeltaPhi(MET);
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("transverse mass of leading lepton"); {
      *cutflow << HFTname("MtLep0");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
          if (m_signalLeptons.size() < 1) return -DBL_MAX;
          return sqrt( 2*m_lepton0.Pt()*MET.Pt()*(1-cos (m_lepton0.DeltaPhi(MET)) ) );
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("transverse mass of subleading lepton"); {
      *cutflow << HFTname("MtLep1");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
          if (m_signalLeptons.size() < 2) return -DBL_MAX;
          return sqrt( 2*m_lepton1.Pt()*MET.Pt()*(1-cos (m_lepton1.DeltaPhi(MET)) ) );
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("approximate tau pT"); {
      *cutflow << HFTname("tau_pT");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
          if (m_signalLeptons.size() < 2) return -DBL_MAX;
          return (MET + m_lepton1).Pt();
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("ratio of tau pT to subleading lepton pT"); {
      *cutflow << HFTname("taulep1_pT_ratio");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
          if (m_signalLeptons.size() < 2) return -DBL_MAX;
          TLorentzVector tau_estimate = MET + m_lepton1;
          return tau_estimate.Pt() / m_lepton0.Pt();
      };
      *cutflow << SaveVar();
    }
}
void add_met_variables(Superflow* superflow, Sel sel_type) {
  //////////////////////////////////////////////////////////////////////////////
  // MET

  // Fill MET variable inside Et var
  TLorentzVector MET;
  *cutflow << NewVar("transverse missing energy (Et)"); {
    *cutflow << HFTname("MET");
    *cutflow << [&](Superlink* sl, var_double*) -> double {
      MET.SetPxPyPzE(m_met.Et*cos(m_met.phi),
                     m_met.Et*sin(m_met.phi),
                     0.,
                     m_met.Et);
      return MET.Pt();
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("transverse missing energy (Phi)"); {
    *cutflow << HFTname("METPhi");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
      return MET.Phi();
    };
    *cutflow << SaveVar();
  }
}
void add_jet_variables(Superflow* superflow, Sel sel_type) {
    ////////////////////////////////////////////////////////////////////////////
    // Multiplicity variables
    ADD_MULTIPLICITY_VAR(preJets)
    ADD_MULTIPLICITY_VAR(baseJets)
    ADD_MULTIPLICITY_VAR(jets)
    //xTauFW variable
    *cutflow << NewVar("number of baseline jets"); {
      *cutflow << HFTname("JetN");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return m_baseJets.size(); };
      *cutflow << SaveVar();
    }
    //xTauFW variable
    *cutflow << NewVar("number of baseline jets (2p4Eta25Pt)"); {
      *cutflow << HFTname("Jet_N2p4Eta25Pt");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
          uint nPassedJets = 0;
          for (auto& jet : m_baseJets) {
              if (fabs(jet->Eta()) < 2.4 && jet->Pt() > 25){
                  nPassedJets += 1;
              }
          }
          return nPassedJets;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of signal jets"); {
      *cutflow << HFTname("SignalJetN");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return m_signalJets.size(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of light jets"); {
      *cutflow << HFTname("nLJets");
      *cutflow << [&](Superlink* sl, var_int*) -> int { return sl->tools->numberOfLJets(m_signalJets)/*(*m_baseJets)*/; };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of b jets"); {
      *cutflow << HFTname("nBJets");
      *cutflow << [&](Superlink* sl, var_int*) -> int { return sl->tools->numberOfBJets(m_signalJets)/*(*m_baseJets)*/; };
      *cutflow << SaveVar();
    }
    //xTauFW variable
    *cutflow << NewVar("b-tagged jet"); {
      *cutflow << HFTname("Btag");
      *cutflow << [&](Superlink* sl, var_bool*) -> bool { return sl->tools->numberOfBJets(m_signalJets) > 0;};
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("number of forward jets"); {
      *cutflow << HFTname("nForwardJets");
      *cutflow << [&](Superlink* sl, var_int*) -> int { return sl->tools->numberOfFJets(m_signalJets)/*(*m_baseJets)*/; };
      *cutflow << SaveVar();
    }
    //////////////////////////////////////////////////////////////////////////////
    // Pre Jets
    *cutflow << NewVar("preJet pt"); {
      *cutflow << HFTname("preJet_pt");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& jet : m_preJets) {out.push_back(jet->Pt()); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("preJet eta"); {
      *cutflow << HFTname("preJet_eta");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& jet : m_preJets) {out.push_back(jet->Eta()); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("preJet JVT if eta<=2.4 and pT<60"); {
      *cutflow << HFTname("preJet_JVT");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& jet : m_preJets) {
          if (jet->Pt() < 60 && fabs(jet->Eta()) <= 2.4) out.push_back(jet->jvt);
          else
            out.push_back(-DBL_MAX);
        }
        return out;
      };
      *cutflow << SaveVar();
    }

    //////////////////////////////////////////////////////////////////////////////
    // Baseline Jets
    *cutflow << NewVar("Baseline Jet eta"); {
      *cutflow << HFTname("baseJet_eta");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& jet : m_baseJets) {out.push_back(jet->Eta()); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Baseline Jet mv2c10"); {
      *cutflow << HFTname("baseJet_mv2c10");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for (auto& jet : m_baseJets) {out.push_back(jet->mv2c10); }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("jet eta"); {
      *cutflow << HFTname("j_eta");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for(auto& jet : m_baseJets) {
          out.push_back(jet->Eta());
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("jet JVT"); {
      *cutflow << HFTname("j_jvt");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for(auto& jet : m_baseJets) {
          out.push_back(jet->jvt);
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("jet JVF"); {
      *cutflow << HFTname("j_jvf");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for(auto& jet : m_baseJets) {
          out.push_back(jet->jvf);
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("jet phi"); {
      *cutflow << HFTname("j_phi");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for(auto& jet : m_baseJets) {
          out.push_back(jet->Phi());
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("jet flavor (0: NA, 1: L, 2: B, 3: F)"); {
      *cutflow << HFTname("j_flav");
      *cutflow << [&](Superlink* sl, var_int_array*) -> vector<int> {
        vector<int> out; int flav = 0;
        for(auto& jet : m_baseJets) {
          if(sl->tools->m_jetSelector->isLight(jet))  { flav = 1; }
          else if(sl->tools->m_jetSelector->isB(jet)) { flav = 2; }
          else if(sl->tools->m_jetSelector->isForward(jet))  { flav = 3; }
          out.push_back(flav);
          flav=0;
        }
        return out;
      };
      *cutflow << SaveVar();
    }

    //////////////////////////////////////////////////////////////////////////////
    // VBF Region Cuts
    *cutflow << NewVar("number of jets pT > 30"); {
      *cutflow << HFTname("JetN_g30");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
          uint nJets = 0;
          for (auto& jet : m_baseJets) {
              if (jet->Pt() > 30) nJets++;
          }
          return nJets;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("jet pt"); {
      *cutflow << HFTname("j_pt");
      *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
        vector<double> out;
        for(auto& jet : m_baseJets) {
          out.push_back(jet->Pt());
        }
        return out;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("dijet mass"); {
      *cutflow << HFTname("Mjj");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
          if (m_baseJets.size() < 2) return -DBL_MAX;
          return m_JetP4.M();
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("DeltaEta between two leading jets"); {
      *cutflow << HFTname("DEtaJJ");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
          if (m_baseJets.size() < 2) return -DBL_MAX;
          return fabs(m_Jet0.Eta() - m_Jet1.Eta());
      };
      *cutflow << SaveVar();
    }
}
void add_fake_variables(Superflow* superflow, Sel sel_type) {
  //////////////////////////////////////////////////////////////////////////////
  // Fake variables
  //////////////////////////////////////////////////////////////////////////////
    *cutflow << NewVar("Number of ID leptons"); {
      *cutflow << HFTname("nLepID");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
          return m_lepID_n;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Number of Anti-ID leptons"); {
      *cutflow << HFTname("nLepAntiID");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
          return m_lepAntiID_n;
      };
      *cutflow << SaveVar();
    }
    //////////////////////////////////////////////////////////////////////////////
    // 1 antiID lepton case
    *cutflow << NewVar("Leading anti-ID lepton charge"); {
      *cutflow << HFTname("aID_Lep0Q");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if (m_lepAntiID_n < 1) return -INT_MAX;
        // TODO: Add antiID_lep0 lepton object
        auto* fakelep0 = m_baseLeptons.at(m_antiID_idx0);
        return fakelep0->q;
        };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Leading aID-lepton pT"); {
      *cutflow << HFTname("aID_Lep0Pt");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 1) return -DBL_MAX;
        return m_LepFake0.Pt(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Leading aID-lepton eta"); {
      *cutflow << HFTname("aID_Lep0Eta");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 1) return -DBL_MAX;
        return m_LepFake0.Eta(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Leading aID-lepton invariant mass"); {
      *cutflow << HFTname("aID_MLep0");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 1) return -DBL_MAX;
        return m_LepFake0.M(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("transverse mass of leading aID-lepton"); {
      *cutflow << HFTname("aID_MtLep0");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
          if (m_lepAntiID_n < 1) return -DBL_MAX;
          return sqrt( 2*m_LepFake0.Pt()*MET.Pt()*(1-cos (m_LepFake0.DeltaPhi(MET)) ) );
      };
      *cutflow << SaveVar();
    }
    // Z + jets variables
    *cutflow << NewVar("Z leptons' flavor"); {
      *cutflow << HFTname("Z_dilep_flav");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if(m_signalLeptons.size() < 2) return -INT_MAX;
        if(m_Zlep.at(0)->isEle() && m_Zlep.at(1)->isMu()){return 0;}       // e mu  case
        else if(m_Zlep.at(0)->isMu() && m_Zlep.at(1)->isEle()){return 1;}  // mu e  case
        else if(m_Zlep.at(0)->isEle() && m_Zlep.at(1)->isEle()){return 2;} // e e   case
        else if(m_Zlep.at(0)->isMu() && m_Zlep.at(1)->isMu()){return 3;}   // mu mu case
        return 4;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Z dilepton sign"); {
      *cutflow << HFTname("Z_dilep_sign");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if(m_signalLeptons.size() < 2) return -INT_MAX;
        return m_Zlep.at(0)->q * m_Zlep.at(1)->q;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("MLL of closest Z pair"); {
      *cutflow << HFTname("Z_MLL");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() < 2) return -DBL_MAX;
        return (*m_Zlep.at(0)+*m_Zlep.at(1)).M();
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("2nd Z leptons' flavor"); {
      *cutflow << HFTname("Z2_dilep_flav");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if(m_signalLeptons.size() < 2) return -INT_MAX;
        else if (!m_Zlep.at(3) or !m_Zlep.at(4)) return -2;
        if(m_Zlep.at(3)->isEle() && m_Zlep.at(4)->isMu()){return 0;}       // e mu  case
        else if(m_Zlep.at(3)->isMu() && m_Zlep.at(4)->isEle()){return 1;}  // mu e  case
        else if(m_Zlep.at(3)->isEle() && m_Zlep.at(4)->isEle()){return 2;} // e e   case
        else if(m_Zlep.at(3)->isMu() && m_Zlep.at(4)->isMu()){return 3;}   // mu mu case
        return 4;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("2nd Z dilepton sign"); {
      *cutflow << HFTname("Z2_dilep_sign");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if(m_signalLeptons.size() < 2) return -INT_MAX;
        else if (!m_Zlep.at(3) or !m_Zlep.at(4)) return -2;
        return m_Zlep.at(3)->q * m_Zlep.at(4)->q;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("MLL of 2nd-closest Z pair"); {
      *cutflow << HFTname("Z2_MLL");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_signalLeptons.size() < 2) return -DBL_MAX;
        else if (!m_Zlep.at(3) or !m_Zlep.at(4)) return -5.0;
        return (*m_Zlep.at(3)+*m_Zlep.at(4)).M();
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Non-Z lepton pT"); {
      *cutflow << HFTname("Z_Lep2_pT");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_Zlep.at(2)) return m_Zlep.at(2)->Pt();
        return -DBL_MAX;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Non-Z lepton eta"); {
      *cutflow << HFTname("Z_Lep2_eta");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_Zlep.at(2)) return m_Zlep.at(2)->Eta();
        return -DBL_MAX;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Non-Z lepton flavor"); {
      *cutflow << HFTname("Z_Lep2_flav");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if (m_Zlep.at(2)) return m_Zlep.at(2)->isEle();
        return -INT_MAX;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Non-Z lepton charge"); {
      *cutflow << HFTname("Z_Lep2_q");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if (m_Zlep.at(2)) return m_Zlep.at(2)->q;
        return -INT_MAX;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("DeltaPhi of Non-Z lep and MET"); {
      *cutflow << HFTname("Z_Lep2_dPhi_MET");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_Zlep.at(2)) return m_Zlep.at(2)->DeltaPhi(MET);
        return -DBL_MAX;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Non-Z lepton mT"); {
      *cutflow << HFTname("Z_Lep2_mT");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_Zlep.at(2)) {
            return sqrt( 2*m_Zlep.at(2)->Pt()*MET.Pt()*(1-cos (m_Zlep.at(2)->DeltaPhi(MET)) ) );
        }
        return -DBL_MAX;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta R of Z and aID-Lepton"); {
      *cutflow << HFTname("dR_Zl");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_Zlep.at(2)) return m_Zlep.at(2)->DeltaR(*m_Zlep.at(0)+*m_Zlep.at(1));
        return -DBL_MAX;
      };
      *cutflow << SaveVar();
    }
    // W + jets variables
    *cutflow << NewVar("aID-ID dilepton flavor"); {
      *cutflow << HFTname("aID_dilep_flav");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if (m_lepAntiID_n < 1 || m_signalLeptons.size() < 1) return -INT_MAX;
        auto* fakelep0 = m_baseLeptons.at(m_antiID_idx0);
        auto* lep0 = m_signalLeptons.at(0);
        if(fakelep0->isEle() && lep0->isMu()){return 0;}       // e mu  case
        else if(fakelep0->isMu() && lep0->isEle()){return 1;}  // mu e  case
        else if(fakelep0->isEle() && lep0->isEle()){return 2;} // e e   case
        else if(fakelep0->isMu() && lep0->isMu()){return 3;}   // mu mu case
        return 4;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("aID-ID sign product"); {
      *cutflow << HFTname("aID_LepLepSign");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if (m_lepAntiID_n < 1 || m_signalLeptons.size() < 1) return -INT_MAX;
        auto* fakelep0 = m_baseLeptons.at(m_antiID_idx0);
        auto* lep0 = m_signalLeptons.at(0);
        return fakelep0->q * lep0->q;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta R of ID-lepton and aID-lepton"); {
      *cutflow << HFTname("aID_drll");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 1 || m_signalLeptons.size() < 1) return -DBL_MAX;
        return m_LepFake0.DeltaR(m_lepton0); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("ID-aID collinear mass, M_coll"); {
      *cutflow << HFTname("aID_MCollASym");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 1 || m_signalLeptons.size() < 1) return -DBL_MAX;
        double deta = fabs(m_LepFake0.Eta()-m_lepton0.Eta());
        double dphi = m_LepFake0.DeltaPhi(m_lepton0);
        return sqrt(2*m_LepFake0.Pt()*(m_lepton0.Pt()+m_met.Et)*(cosh(deta)-cos(dphi)));
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("mass of ID-aID-dilepton system, M_ll"); {
      *cutflow << HFTname("aID_MLL");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 1 || m_signalLeptons.size() < 1) return -DBL_MAX;
        return (m_LepFake0+m_lepton0).M();
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Pt of ID-aID-dilepton system, Pt_ll"); {
      *cutflow << HFTname("aID_ptll");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 1 || m_signalLeptons.size() < 1) return -DBL_MAX;
        return (m_LepFake0+m_lepton0).Pt(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("diff_pt between ID-aID-leptons"); {
      *cutflow << HFTname("aID_dpt_ll");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
          if (m_lepAntiID_n < 1 || m_signalLeptons.size() < 1) return -DBL_MAX;
          return m_LepFake0.Pt() - m_lepton0.Pt();
      };
      *cutflow << SaveVar();
    }

    //////////////////////////////////////////////////////////////////////////////
    // 2 antiID lepton case (QCD)
    *cutflow << NewVar("Subleading anti-ID lepton flavor"); {
      *cutflow << HFTname("aID_Lep1Flav");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
          if (m_lepAntiID_n < 2) return -INT_MAX;
          auto* fakelep1 = m_baseLeptons.at(m_antiID_idx1);
          return  fakelep1->isEle() ? 0 : 1;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Subleading anti-ID lepton charge"); {
      *cutflow << HFTname("aID_Lep1Q");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if (m_lepAntiID_n < 2) return -INT_MAX;
        auto* fakelep1 = m_baseLeptons.at(m_antiID_idx1);
        return fakelep1->q;
        };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Subleading aID-lepton pT"); {
      *cutflow << HFTname("aID2_Lep1Pt");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 2) return -DBL_MAX;
        return m_LepFake1.Pt(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Subleading aID-lepton eta"); {
      *cutflow << HFTname("aID2_Lep1Eta");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 2) return -DBL_MAX;
        return m_LepFake1.Eta(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Subleading aID-lepton invariant mass"); {
      *cutflow << HFTname("aID2_MLep1");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 2) return -DBL_MAX;
        return m_LepFake1.M(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta R of aID-leptons"); {
      *cutflow << HFTname("aID2_drll");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 2) return -DBL_MAX;
        return m_LepFake0.DeltaR(m_LepFake1); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("aID-dilepton flavor"); {
      *cutflow << HFTname("aID2_dilep_flav");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if(m_lepAntiID_n < 2) return -INT_MAX;
        auto* fakelep0 = m_baseLeptons.at(m_antiID_idx0);
        auto* fakelep1 = m_baseLeptons.at(m_antiID_idx1);
        if(fakelep0->isEle() && fakelep1->isMu()){return 0;}       // e mu  case
        else if(fakelep0->isMu() && fakelep1->isEle()){return 1;}  // mu e  case
        else if(fakelep0->isEle() && fakelep1->isEle()){return 2;} // e e   case
        else if(fakelep0->isMu() && fakelep1->isMu()){return 3;}   // mu mu case
        return 4;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("aID collinear mass, M_coll"); {
      *cutflow << HFTname("aID2_MCollASym");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 2) return -DBL_MAX;
        double deta = fabs(m_LepFake0.Eta()-m_LepFake1.Eta());
        double dphi = m_LepFake0.DeltaPhi(m_LepFake1);
        return sqrt(2*m_LepFake0.Pt()*(m_LepFake1.Pt()+m_met.Et)*(cosh(deta)-cos(dphi)));
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("mass of aID-dilepton system, M_ll"); {
      *cutflow << HFTname("aID2_MLL");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 2) return -DBL_MAX;
        return (m_LepFake0+m_LepFake1).M(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Pt of aID-dilepton system, Pt_ll"); {
      *cutflow << HFTname("aID2_ptll");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (m_lepAntiID_n < 2) return -DBL_MAX;
        return (m_LepFake0+m_LepFake1).Pt(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("diff_pt between aID-leptons"); {
      *cutflow << HFTname("aID2_dpt_ll");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
          if (m_lepAntiID_n < 2) return -DBL_MAX;
          return m_LepFake0.Pt() - m_LepFake1.Pt();
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("transverse mass of subleading aID-lepton"); {
      *cutflow << HFTname("aID2_MtLep1");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
          if (m_lepAntiID_n < 2) return -DBL_MAX;
          return sqrt( 2*m_LepFake1.Pt()*MET.Pt()*(1-cos (m_LepFake1.DeltaPhi(MET)) ) );
      };
      *cutflow << SaveVar();
    }
}
void add_shortcut_variables_reset(Superflow* superflow, Sel sel_type) {
  *cutflow << [&](Superlink* /*sl*/, var_void*) {
    m_preLeptons.clear(); m_baseLeptons.clear(); m_signalLeptons.clear(), m_selectLeptons.clear();
    m_preElectrons.clear(); m_baseElectrons.clear(); m_signalElectrons.clear();
    m_preMuons.clear(); m_baseMuons.clear(); m_signalMuons.clear();
    m_dileptonP4 = m_lepton0 = m_lepton1 = {};
    m_preTaus.clear(); m_baseTaus.clear(); m_signalTaus.clear();
    m_preJets.clear(); m_baseJets.clear(); m_signalJets.clear();
    m_lightJets.clear(); m_BJets.clear(); m_forwardJets.clear();
    m_JetP4 = m_Jet1 = m_Jet0 = {};
    m_lepID_n = m_lepAntiID_n = 0;
    m_antiID_idx0 = m_antiID_idx1 = -1;
    m_LepFake0 = m_LepFake1 = {};
    std::fill(m_Zlep.begin(), m_Zlep.end(), nullptr);
  };
}

////////////////////////////////////////////////////////////////////////////////
void add_fake_shortcut_variables(Superflow* superflow, Sel sel_type) {

    // ID Leptons - same as signal leptons
    for (Susy::Lepton* lepton : m_baseLeptons) {
      if (is_ID_lepton(lepton)) m_lepID_n++;
    }
    // Anti-ID Leptons
    for (Susy::Lepton* lepton : m_baseLeptons) {
      if (is_antiID_lepton(lepton)) m_lepAntiID_n++;
    }

    std::vector<int>::iterator begin = lepAntiID_vec.begin();
    std::vector<int>::iterator end = lepAntiID_vec.end();
    if (m_lepAntiID_n >= 1) {
      m_antiID_idx0 = std::find(begin, end, true) - begin;
      m_LepFake0 = m_baseLeptons.at(m_antiID_idx0);
    }
    if (m_lepAntiID_n >= 2) {
      m_antiID_idx1 = std::find(begin + m_antiID_idx0 + 1, end, true) - begin;
      m_LepFake1 = m_baseLeptons.at(m_antiID_idx1);
    }

    // Find ID lepton combination with MLL closest to Z-peak
    float Z_diff = FLT_MAX;
    int zlep_idx2 = -1;
    for (uint ii = 0; ii < m_signalLeptons.size(); ++ii) {.at(1] = lep_jj;
                zlep_idx2 = ii > 0 ? 0 : jj > 1 ? 1 : 2;
            }
        })    }

    if (m_signalLeptons.size() == 2 && m_lepAntiID_n >= 1) {.at(2] = m_baseLeptons.at(m_antiID_idx0);)    } else if (m_signalLeptons.size() >= 3) {.at(2] = m_signalLeptons.at(zlep_idx2);)    }
    // Find base lepton combination with MLL closest to Zpeak
    // Exclude ID leptons already matched with Z
    Z_diff = FLT_MAX;
    for (uint ii = 0; ii < m_baseLeptons.size(); ++ii) {.at(4] = lep_jj;
            }
        })    }
}

bool is_ID_lepton(Susy::Lepton* leptn) {ton* sig_lepton : m_signalLeptons) {
        if (lepton == sig_lepton) {
            return true;
        }
    }
    return false;
}

bool is_antiID_lepton(Susy::Leptn* lepton) {sy::Electron*>(lepton);
        float absEtaBE = fabs(ele->clusEtaBE);
        pt_pass  = ele->pt > 15;
        eta_pass = absEtaBE < 1.37 || 1.52 < absEtaBE;
        iso_pass = ele->isoGradient;
        id_pass  = ele->mediumLLH;
        passedID_cuts = iso_pass && id_pass;
        passedAntiID_cuts = ele->veryLooseLLH;
    } else if (lepton->isMu()) {
        const Susy::Muon* mu = dynamic_cast<const Susy::Muon*>(lepton);
        pt_pass  = mu->pt > 10;
        eta_pass = fabs(mu->eta) < 2.47;
        iso_pass = mu->isoGradient;
        id_pass  = mu->medium;
        //TODO: Study effect of removing id_pass so muons only fail iso requirements
        passedID_cuts = iso_pass && id_pass;
        passedAntiID_cuts = 1;
    }
    bool lepAntiID_bool = pt_pass && eta_pass
                        && passedAntiID_cuts && !passedID_cuts;

    return lepAntiID_bool;
}

void add_SFOS_lepton_cut(Superflow* superflow) {
    *cutflow << CutName("SFOS leptons") << [&](Superlink /*sl*/) -> bool {
        bool SF = m_selectLeptons.at(0)->isEle() == m_selectLeptons.at(1)->isEle();
        bool OS = m_selectLeptons.at(0)->q != m_selectLeptons.at(1)->q;
        return SF && OS;
    };
}

void add_DFOS_lepton_cut(Superflow* superflow) {
    *cutflow << CutName("DFOS leptons") << [&](Superlink* sl) -> bool {
        bool DF = m_selectLeptons.at(0)->isEle() != m_selectLeptons.at(1)->isEle();
        bool OS = m_selectLeptons.at(0)->q != m_selectLeptons.at(1)->q;
        return DF && OS;
    };
}




