//////////////////////////////////////////////////////////////////////////////////////////
// \project:    ATLAS Experiment at CERN's LHC
// \package:    MyAnalysis
// \file:       $Id$
// \author:     Alaettin.Serhan.Mete@cern.ch
// \history:    N/A
//
// Copyright (C) 2016 University of California, Irvine
//////////////////////////////////////////////////////////////////////////////////////////

#include <cmath>

// analysis include(s)
#include "Superflow/Superflow.h"
#include "Superflow/Superlink.h"
#include "SusyNtuple/ChainHelper.h"
#include "SusyNtuple/string_utils.h"
#include "SusyNtuple/KinematicTools.h"
//#include "SusyNtuple/SusyNt.h"
//#include "SusyNtuple/SusyDefs.h"
//#include "SusyNtuple/SusyNtObject.h"
//#include "SusyNtuple/SusyNtTools.h"

#define ADD_TRIGGER_VAR(trig_name) { \
  *cutflow << NewVar(#trig_name" trigger bit"); { \
      *cutflow << HFTname(#trig_name); \
      *cutflow << [](Superlink* sl, var_bool*) -> bool { \
          return sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, #trig_name); }; \
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
using namespace sflow;
//using namespace Susy;

///////////////////////////////////////////////////////////////////////
// Usage
///////////////////////////////////////////////////////////////////////
void usage(std::string progName)
{
  printf("=================================================================\n");
  printf("%s [options]\n",progName.c_str());
  printf("=================================================================\n");
  printf("Options:\n");
  printf("-h        Print this help\n");
  printf("-n        Number of events to be processed (default: -1)\n");
  printf("-f        Input file as *.root, list of *.root in a *.txt,\n");
  printf("          or a DIR/ containing *.root (default: none)\n");
  printf("=================================================================\n");
}

///////////////////////////////////////////////////////////////////////
// Print event information
///////////////////////////////////////////////////////////////////////
void printEventInformation(Superlink* sl) {
  // event
  printf("==============================================\n");
  sl->nt->evt()->print();
  // pre-obj
  printf("==============================================\n");
  printf("Pre-objects: \n");
  for(auto &obj : *(sl->preElectrons))  { obj->print(); }
  printf("\n");
  for(auto &obj : *(sl->preMuons))      { obj->print(); }
  printf("\n");
  for(auto &obj : *(sl->preJets))       { obj->print(); }
  // base-obj
  printf("==============================================\n");
  printf("Base-objects: \n");
  for(auto &obj : *(sl->baseElectrons)) { obj->print(); }
  printf("\n");
  for(auto &obj : *(sl->baseMuons))     { obj->print(); }
  printf("\n");
  for(auto &obj : *(sl->baseJets))      { obj->print(); }
  // signal-obj
  printf("==============================================\n");
  printf("Signal-objects: \n");
  for(auto &obj : *(sl->electrons))     { obj->print(); }
  printf("\n");
  for(auto &obj : *(sl->muons))         { obj->print(); }
  printf("\n");
  for(auto &obj : *(sl->jets))          { obj->print(); }
  printf("==============================================\n");
  return;
}

///////////////////////////////////////////////////////////////////////
// Main function
///////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  //////////////////////////////////////////////////////////////////////////////
  // Read user inputs - NOT super safe so be careful :)
  unsigned int n_events  = -1;
  unsigned int n_skip_events  = 0;
  char *input_file  = nullptr;
  char *name_suffix = nullptr;
  SuperflowRunMode run_mode = SuperflowRunMode::nominal; // SuperflowRunMode::all_syst; //SuperflowRunMode::nominal;
  int c;

  opterr = 0;
  while ((c = getopt (argc, argv, "f:s:n:h")) != -1)
    switch (c)
      {
      case 'f':
        input_file = optarg;
        break;
      case 's':
        name_suffix = optarg;
        break;
      case 'n':
        n_events = atoi(optarg);
        break;
      case 'h':
        usage("makeMiniNtuple");
        return 1;
      case '?':
        if (optopt == 'f')
          fprintf (stderr, "makeMiniNtuple\t Option -%c requires an argument.\n", optopt);
        else if (optopt == 'n')
          fprintf (stderr, "makeMiniNtuple\t Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "makeMiniNtuple\t Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "makeMiniNtuple\t Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();
      }

  // Catch problems or cast
  for (int index = optind; index < argc; index++)
    printf ("makeMiniNtuple\t Non-option argument %s\n", argv[index]);
  if (input_file==nullptr) {
    printf("makeMiniNtuple\t An input file must be provided with option -f (a list, a DIR or single file)\n");
    return 0;
  }

  // Print information
  printf("makeMiniNtuple\t =================================================================\n");
  printf("makeMiniNtuple\t Running MyAnalysis/makeMiniNtuple\n");
  printf("makeMiniNtuple\t =================================================================\n");
  printf("makeMiniNtuple\t   Flags:\n");
  printf("makeMiniNtuple\t     Input file (-f)         : %s\n",input_file);
  printf("makeMiniNtuple\t     Number of events (-n)   : %i\n",n_events );
  printf("makeMiniNtuple\t =================================================================\n");
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  // Setup Superflow
  TChain* chain = new TChain("susyNt");
  chain->SetDirectory(0);

  bool inputIsFile = Susy::utils::endswith(input_file, ".root");
  bool inputIsList = Susy::utils::endswith(input_file, ".txt");
  bool inputIsDir  = Susy::utils::endswith(input_file, "/");

  if(inputIsFile) {
    ChainHelper::addFile(chain, input_file);
  } else if (inputIsList) {
    // If a list of ROOT files
    ChainHelper::addFileList(chain, input_file);
  } else if (inputIsDir) {
    ChainHelper::addFileDir(chain, input_file);
  } else {
    printf("makeMiniNtuple\t Cannot understand input %s",input_file);
    return 0;
  }

  Superflow* cutflow = new Superflow(); // initialize the cutflow
  cutflow->setAnaName("SuperflowAna");                // arbitrary
  cutflow->setAnaType(AnalysisType::Ana_HLFV);        // analysis type, passed to SusyNt
  cutflow->setLumi(1.0);                              // set the MC normalized to X pb-1
  cutflow->setSampleName(input_file);                 // sample name, check to make sure it's set OK
  cutflow->setRunMode(run_mode);                      // make configurable via run_mode
  cutflow->setCountWeights(true);                     // print the weighted cutflows
  if(name_suffix != nullptr) cutflow->setFileSuffix(name_suffix);
  cutflow->setChain(chain);
  cutflow->nttools().initTriggerTool(ChainHelper::firstFile(input_file,0.));

  printf("makeMiniNtuple\t Total events available : %lli\n",chain->GetEntries());

  //////////////////////////////////////////////////////////////////////////////
  //
  // SUPERFLOW METHODS BEGIN HERE
  //
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  //
  // Add Cuts to cutflow
  //
  //////////////////////////////////////////////////////////////////////////////

  // All events before cuts
  *cutflow << CutName("read in") << [](Superlink* /*sl*/) -> bool { return true; };

  //  Cleaning Cuts
  // How to add cutflow entry:
  // Create lambda function that returns bool value of cut.
  // Pass that with "<<" into the function CutName("Cut name").
  // Pass that with "<<" into the dereferenced cutflow object.

  // xTauFW Cut
  *cutflow << CutName("xTau: 2 Loose Leptons") << [&](Superlink* sl) -> bool {
      uint nLooseLeptons = 0;
      for (const auto* mu : *sl->preMuons) {if (mu->loose) nLooseLeptons++;}
      for (const auto* ele : *sl->preElectrons) {if (ele->looseLLH) nLooseLeptons++;}
      return nLooseLeptons == 2;
  };

  int cutflags = 0;
  *cutflow << CutName("Pass GRL") << [&](Superlink* sl) -> bool {
      cutflags = sl->nt->evt()->cutFlags[NtSys::NOM];
      return (sl->tools->passGRL(cutflags));
  };
  *cutflow << CutName("LAr error") << [&](Superlink* sl) -> bool {
      return (sl->tools->passLarErr(cutflags));
  };
  *cutflow << CutName("Tile error") << [&](Superlink* sl) -> bool {
      return (sl->tools->passTileErr(cutflags));
  };
  *cutflow << CutName("TTC veto") << [&](Superlink* sl) -> bool {
      return (sl->tools->passTTC(cutflags));
  };
  *cutflow << CutName("SCT err") << [&](Superlink* sl) -> bool {
      return (sl->tools->passSCTErr(cutflags));
  };
  *cutflow << CutName("nBaselineLep = nSignalLep") << [](Superlink* sl) -> bool {
      return (sl->leptons->size() == sl->baseLeptons->size());
  };
  *cutflow << CutName("2+ leptons") << [](Superlink* sl) -> bool {
      return (sl->leptons->size() >= 2);
  };
  *cutflow << CutName("DFOS leptons") << [&](Superlink* sl) -> bool {
      bool DF = sl->leptons->at(0)->isEle() != sl->leptons->at(1)->isEle();
      bool OS = sl->leptons->at(0)->q != sl->leptons->at(1)->q;
      return DF && OS;
  };
  *cutflow << CutName("pass good vertex") << [&](Superlink* sl) -> bool {
      return (sl->tools->passGoodVtx(cutflags));
  };
  *cutflow << CutName("jet cleaning") << [&](Superlink* sl) -> bool {
      return (sl->tools->passJetCleaning(sl->baseJets));
  };
  *cutflow << CutName("bad muon veto") << [&](Superlink* sl) -> bool {
      return (sl->tools->passBadMuon(sl->preMuons));
  };
  *cutflow << CutName("2 leptons") << [](Superlink* sl) -> bool {
      return (sl->leptons->size() == 2);
  };
  *cutflow << CutName("Tau veto") << [](Superlink* sl) -> bool {
      return (sl->taus->size() == 0);
  };

  //////////////////////////////////////////////////////////////////////////////
  //
  // Variables stored in output nTuple
  //
  // - Accesing objects through superlink and assigning them to internal variables
  //
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  // Event level details
  //////////////////////////////////////////////////////////////////////////////
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
    *cutflow << [](Superlink* sl, var_double*) -> double {
        return sl->weights->product() * sl->nt->evt()->wPileup;
        //return sl->weights->product();
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

  //////////////////////////////////////////////////////////////////////////////
  // Trigger Variables
  // ADD_TRIGGER_VAR preprocessor defined
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  // 2015

  // Dilepton Triggers
  ADD_TRIGGER_VAR(HLT_e17_lhloose_mu14)
  ADD_TRIGGER_VAR(HLT_e24_lhmedium_L1EM20VHI_mu8noL1)
  ADD_TRIGGER_VAR(HLT_e7_lhmedium_mu24)

  // Single Electron Triggers
  ADD_TRIGGER_VAR(HLT_e24_lhmedium_L1EM20VH)
  ADD_TRIGGER_VAR(HLT_e60_lhmedium)
  ADD_TRIGGER_VAR(HLT_e120_lhloose)

  // Single Muon Triggers
  ADD_TRIGGER_VAR(HLT_mu20_iloose_L1MU15)
  ADD_TRIGGER_VAR(HLT_mu40)

  //////////////////////////////////////////////////////////////////////////////
  // 2016

  // Dilepton Triggers
  ADD_TRIGGER_VAR(HLT_e17_lhloose_nod0_mu14)
  ADD_TRIGGER_VAR(HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1)
  ADD_TRIGGER_VAR(HLT_e7_lhmedium_nod0_mu24)

  // Single Electron Triggers
  ADD_TRIGGER_VAR(HLT_e26_lhtight_nod0_ivarloose)
  ADD_TRIGGER_VAR(HLT_e60_lhmedium_nod0)
  ADD_TRIGGER_VAR(HLT_e140_lhloose_nod0)

  // Single Muon Triggers
  ADD_TRIGGER_VAR(HLT_mu26_ivarmedium)
  ADD_TRIGGER_VAR(HLT_mu50)


  //////////////////////////////////////////////////////////////////////////////
  // Lepton informtion
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  // Define some space saving variables
  LeptonVector preLeptons, baseLeptons, signalLeptons;
  ElectronVector preElectrons, baseElectrons, signalElectrons;
  MuonVector preMuons, baseMuons, signalMuons;
  TauVector preTaus, baseTaus, signalTaus;
  Susy::Met met;
  *cutflow << [&](Superlink* sl, var_void*) {
    preLeptons = *sl->preLeptons;
    baseLeptons =  *sl->baseLeptons;
    signalLeptons = *sl->leptons;

    preElectrons = *sl->preElectrons;
    baseElectrons =  *sl->baseElectrons;
    signalElectrons = *sl->electrons;

    preMuons = *sl->preMuons;
    baseMuons =  *sl->baseMuons;
    signalMuons = *sl->muons;

    preTaus = *sl->preTaus;
    baseTaus =  *sl->baseTaus;
    signalTaus = *sl->taus;

    met = *sl->met;
  };

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
      for (auto& el : preElectrons) {
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
      for (auto& el : preElectrons) {
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
      for (auto& el : preElectrons) {
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
      for (auto& el : preElectrons) {out.push_back(el->Pt()); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("preElectron clusterEtaBE"); {
    *cutflow << HFTname("preEl_clusEtaBE");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for (auto& el : preElectrons) {out.push_back(el->clusEtaBE);}
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("preElectron eta"); {
    *cutflow << HFTname("preEl_eta");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for (auto& el : preElectrons) {out.push_back(el->clusEtaBE);}
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
      for (auto& el : baseElectrons) {out.push_back(el->etconetopo20); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Baseline Electron ptvarcone20"); {
    *cutflow << HFTname("baseEl_ptvarcone20");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for (auto& el : baseElectrons) {out.push_back(el->ptvarcone20); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Baseline Electron ID (non-inclusive)"); {
    *cutflow << HFTname("baseEl_ID");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for (auto& el : baseElectrons) {
        if (!el->veryLooseLLH) out.push_back(0);
        else if (el->veryLooseLLH && !el->looseLLH) out.push_back(1);
        else if (el->looseLLH && !el->looseLLHBLayer) out.push_back(2);
        else if (el->looseLLHBLayer && !el->mediumLLH) out.push_back(3);
        else if (el->mediumLLH && !el->tightLLH) out.push_back(4);
        else if (el->tightLLH) out.push_back(5);
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
      for (auto& el : signalElectrons) {out.push_back(el->mcType); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Electron origin"); {
    *cutflow << HFTname("el_origin");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for (auto& el : signalElectrons) {out.push_back(el->mcOrigin); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("subleading electron track pt"); {
    *cutflow << HFTname("el1_track_pt");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
      if (signalLeptons.size() > 1 && signalLeptons.at(1)->isEle()) {
          const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(signalLeptons.at(1));
          return ele ? ele->trackPt : -2;
      } else {
          return -1;
      }
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("subleading electron clus pt"); {
    *cutflow << HFTname("el1_clus_pt");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
      if (signalLeptons.size() > 1 && signalLeptons.at(1)->isEle()) {
          const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(signalLeptons.at(1));
          return ele ? ele->clusE : -2;
      } else {
          return -1;
      }
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Subleading Electron pT track-cluster ratio"); {
    *cutflow << HFTname("el1pT_trackclus_ratio");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
      if (signalLeptons.size() > 1 && signalLeptons.at(1)->isEle()) {
          const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(signalLeptons.at(1));
          return ele ? ele->trackPt / ele->clusE : -2;
      } else {
          return -1;
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
      for (auto& mu : preMuons) {out.push_back(mu->Pt()); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("PreMuon ID (non-inclusive)"); {
    *cutflow << HFTname("preMuon_ID");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for (auto& mu : preMuons) {
        if (!mu->veryLoose) out.push_back(0);
        else if (mu->veryLoose && !mu->loose) out.push_back(1);
        else if (mu->loose && !mu->medium) out.push_back(2);
        else if (mu->medium && !mu->tight) out.push_back(3);
        else if (mu->tight) out.push_back(4);
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
      for (auto& mu : baseMuons) {out.push_back(mu->Pt()); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Baseline Muon eta"); {
    *cutflow << HFTname("baseMu_eta");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for (auto& mu : baseMuons) {out.push_back(mu->Eta()); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Baseline Muon etconetopo20"); {
    *cutflow << HFTname("baseMu_etconetopo20");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for (auto& mu : baseMuons) {out.push_back(mu->etconetopo20); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Baseline Muon ptvarcone30"); {
    *cutflow << HFTname("baseMu_ptvarcone30");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for (auto& mu : baseMuons) {out.push_back(mu->ptvarcone30); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Baseline Muon ID (non-inclusive)"); {
    *cutflow << HFTname("baseMuon_ID");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for (auto& mu : baseMuons) {
        if (!mu->veryLoose) out.push_back(0);
        else if (mu->veryLoose && !mu->loose) out.push_back(1);
        else if (mu->loose && !mu->medium) out.push_back(2);
        else if (mu->medium && !mu->tight) out.push_back(3);
        else if (mu->tight) out.push_back(4);
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
      for (auto& mu : signalMuons) {out.push_back(mu->mcType); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Muon origin"); {
    *cutflow << HFTname("mu_origin");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for (auto& mu : signalMuons) {out.push_back(mu->mcOrigin); }
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
      for (auto& tau : preTaus) {out.push_back(tau->q); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("preTau nTracks"); {
    *cutflow << HFTname("preTau_nTracks");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for (auto& tau : preTaus) {out.push_back(tau->nTrack); }
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
      for (auto& tau : preTaus) {out.push_back(tau->Pt()); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Baseline Tau eta"); {
    *cutflow << HFTname("baseTau_eta");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for (auto& tau : preTaus) {out.push_back(tau->Eta()); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Baseline Tau nTracks"); {
    *cutflow << HFTname("baseTau_nTracks");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for (auto& tau : baseTaus) {out.push_back(tau->nTrack); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Baseline Tau ID (non-inclusive)"); {
    *cutflow << HFTname("baseTau_ID");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for (auto& tau : baseTaus) {
        if (!tau->loose) out.push_back(0);
        else if (tau->loose && !tau->medium) out.push_back(1);
        else if (tau->medium && !tau->tight) out.push_back(2);
        else if (tau->tight) out.push_back(3);
      }
      return out;
    };
    *cutflow << SaveVar();
  }

  //////////////////////////////////////////////////////////////////////////////
  // Signal Leptons
  *cutflow << NewVar("lepton pt"); {
    *cutflow << HFTname("l_pt");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for (auto& lepton : signalLeptons) {
        out.push_back(lepton->Pt());
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("lepton eta"); {
    *cutflow << HFTname("l_eta");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
        out.push_back(lepton->Eta());
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("lepton phi"); {
    *cutflow << HFTname("l_phi");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
        out.push_back(lepton->Phi());
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("lepton flavor (0: e, 1: m)"); {
    *cutflow << HFTname("l_flav");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : signalLeptons) {
        out.push_back(lepton->isEle() ? 0 : 1);
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("lepton type"); {
    *cutflow << HFTname("l_type");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
        out.push_back(lepton->mcType);
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("lepton origin"); {
    *cutflow << HFTname("l_origin");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
        out.push_back(lepton->mcOrigin);
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("lepton charge"); {
    *cutflow << HFTname("l_q");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
          out.push_back(lepton->q);
      }
      return out;
      };
    *cutflow << SaveVar();
  }
  //xTauFW variable
  *cutflow << NewVar("leptons sign product"); {
    *cutflow << HFTname("LepLepSign");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
      return signalLeptons.at(0)->q*signalLeptons.at(1)->q;
    };
    *cutflow << SaveVar();
  }

  //////////////////////////////////////////////////////////////////////////////
  // MET

  // Fill MET variable inside Et var
  TLorentzVector MET;
  *cutflow << NewVar("transverse missing energy (Et)"); {
    *cutflow << HFTname("MET");
    *cutflow << [&](Superlink* sl, var_float*) -> double {
      MET.SetPxPyPzE(sl->met->Et*cos(sl->met->phi),
                     sl->met->Et*sin(sl->met->phi),
                     0.,
                     sl->met->Et);
      return MET.Pt();
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("transverse missing energy (Phi)"); {
    *cutflow << HFTname("METPhi");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return MET.Phi(); };
    *cutflow << SaveVar();
  }

  //////////////////////////////////////////////////////////////////////////////
  // Two-lepton properties

  // Gather useful information
  TLorentzVector lepton0 ;
  TLorentzVector lepton1 ;
  TLorentzVector dileptonP4 ;
  *cutflow << [&](Superlink* /*sl*/, var_void*) {
    lepton0 = *signalLeptons.at(0);
    lepton1 = *signalLeptons.at(1);
    dileptonP4 = lepton0 + lepton1;
  };

  // Lepton variables
  //xTauFW variable
  *cutflow << NewVar("Leading lepton pT"); {
    *cutflow << HFTname("Lep0Pt");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return lepton0.Pt(); };
    *cutflow << SaveVar();
  }
  //xTauFW variable
  *cutflow << NewVar("Subleading lepton pT"); {
    *cutflow << HFTname("Lep1Pt");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return lepton1.Pt(); };
    *cutflow << SaveVar();
  }
  //xTauFW variable
  *cutflow << NewVar("Leading lepton eta"); {
    *cutflow << HFTname("Lep0Eta");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return lepton0.Eta(); };
    *cutflow << SaveVar();
  }
  //xTauFW variable
  *cutflow << NewVar("Subleading lepton eta"); {
    *cutflow << HFTname("Lep1Eta");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return lepton1.Eta(); };
    *cutflow << SaveVar();
  }
  //xTauFW variable
  *cutflow << NewVar("Leading lepton phi"); {
    *cutflow << HFTname("Lep0Phi");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return lepton0.Phi(); };
    *cutflow << SaveVar();
  }
  //xTauFW variable
  *cutflow << NewVar("Subleading lepton phi"); {
    *cutflow << HFTname("Lep1Phi");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return lepton1.Phi(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("Leading lepton invariant mass"); {
    *cutflow << HFTname("MLep0");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return lepton0.M(); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Subleading lepton invariant mass"); {
    *cutflow << HFTname("MLep1");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return lepton1.M(); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Delta eta between leptons"); {
    *cutflow << HFTname("DEtaLL");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
       return fabs(lepton0.Eta() - lepton1.Eta()); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Delta phi between leptons"); {
    *cutflow << HFTname("DphiLL");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
       return fabs(lepton0.DeltaPhi(lepton1));};
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("delta R of di-lepton system"); {
    *cutflow << HFTname("drll");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return lepton0.DeltaR(lepton1); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("dilepton flavor"); {
    *cutflow << HFTname("dilep_flav");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
	    if(signalLeptons.at(0)->isEle() && signalLeptons.at(1)->isMu()){return 0;}       // e mu  case
        else if(signalLeptons.at(0)->isMu() && signalLeptons.at(1)->isEle()){return 1;}  // mu e  case
        else if(signalLeptons.at(0)->isEle() && signalLeptons.at(1)->isEle()){return 2;} // e e   case
        else if(signalLeptons.at(0)->isMu() && signalLeptons.at(1)->isMu()){return 3;}   // mu mu case
	    else{return 4;}
        };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("dilepton flavor is e mu"); {
    *cutflow << HFTname("isEM");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
	    return signalLeptons.at(0)->isEle() && signalLeptons.at(1)->isMu();  // e mu  case
        };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("dilepton flavor is mu e"); {
    *cutflow << HFTname("isME");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        return signalLeptons.at(0)->isMu() && signalLeptons.at(1)->isEle();  // mu e  case
        };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("collinear mass, M_coll"); {
    *cutflow << HFTname("MCollASym");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
	  double deta = fabs(lepton0.Eta()-lepton1.Eta());
	  double dphi = lepton0.DeltaPhi(lepton1);
	  return sqrt(2*lepton0.Pt()*(lepton1.Pt()+met.Et)*(cosh(deta)-cos(dphi)));
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("mass of di-lepton system, M_ll"); {
    *cutflow << HFTname("MLL");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return dileptonP4.M(); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Pt of di-lepton system, Pt_ll"); {
    *cutflow << HFTname("ptll");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return dileptonP4.Pt(); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("diff_pt between leading and sub-leading lepton"); {
    *cutflow << HFTname("dpt_ll");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
        return lepton0.Pt() - lepton1.Pt();
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("dphi between leading lepton and MET"); {
    *cutflow << HFTname("DphiLep0MET");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
        return lepton0.DeltaPhi(MET);
    };
    *cutflow << SaveVar();
  }

  //////////////////////////////////////////////////////////////////////////////
  // Optimized Region Cuts
  *cutflow << NewVar("dphi between sub-leading lepton and MET"); {
    *cutflow << HFTname("DphiLep1MET");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
        return lepton1.DeltaPhi(MET);
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("transverse mass of leading lepton"); {
    *cutflow << HFTname("MtLep0");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
        return sqrt( 2*lepton0.Pt()*MET.Pt()*(1-cos (lepton0.DeltaPhi(MET)) ) );
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("transverse mass of subleading lepton"); {
    *cutflow << HFTname("MtLep1");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
        return sqrt( 2*lepton1.Pt()*MET.Pt()*(1-cos (lepton1.DeltaPhi(MET)) ) );
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("approximate tau pT"); {
    *cutflow << HFTname("tau_pT");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
        return (MET + lepton1).Pt();
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("ratio of tau pT to subleading lepton pT"); {
    *cutflow << HFTname("taulep1_pT_ratio");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
        TLorentzVector tau_estimate = MET + lepton1;
        return tau_estimate.Pt() / lepton0.Pt();
    };
    *cutflow << SaveVar();
  }


  //////////////////////////////////////////////////////////////////////////////
  // Jet Variables
  //////////////////////////////////////////////////////////////////////////////

  // Useful variables
  JetVector preJets, baseJets, signalJets;
  JetVector centralLightJets, centralBJets, forwardJets;
  TLorentzVector JetP4, Jet1, Jet0;

  *cutflow << [&](Superlink* sl, var_void*) {
    preJets = *sl->preJets;
    baseJets = *sl->baseJets;
    signalJets = *sl->jets;

    if (baseJets.size() > 0) {
        Jet0 = *baseJets.at(0);
        if (baseJets.size() > 1) {
           Jet1 = *baseJets.at(1);
           JetP4 = Jet0 + Jet1;
        }
    }
    for (auto& jet : baseJets) {
          if (sl->tools->m_jetSelector->isCentralLight(jet))  {
               centralLightJets.push_back(jet);
          } else if (sl->tools->m_jetSelector->isCentralB(jet)) {
              centralBJets.push_back(jet);
          } else if (sl->tools->m_jetSelector->isForward(jet))  {
              forwardJets.push_back(jet);
          }
    }
    std::sort(centralLightJets.begin()  , centralLightJets.end()  , comparePt);
    std::sort(centralBJets.begin()      , centralBJets.end()      , comparePt);
    std::sort(forwardJets.begin()       , forwardJets.end()       , comparePt);
  };

  //////////////////////////////////////////////////////////////////////////////
  // Multiplicity variables
  ADD_MULTIPLICITY_VAR(preJets)
  ADD_MULTIPLICITY_VAR(baseJets)
  ADD_MULTIPLICITY_VAR(jets)
  //xTauFW variable
  *cutflow << NewVar("number of baseline jets"); {
    *cutflow << HFTname("JetN");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return baseJets.size(); };
    *cutflow << SaveVar();
  }
  //xTauFW variable
  *cutflow << NewVar("number of baseline jets (2p4Eta25Pt)"); {
    *cutflow << HFTname("Jet_N2p4Eta25Pt");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        uint nPassedJets = 0;
        for (auto& jet : baseJets) {
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
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return signalJets.size(); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("number of central light jets"); {
    *cutflow << HFTname("nCentralLJets");
    *cutflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfCLJets(*sl->baseJets)/*(*baseJets)*/; };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("number of central b jets"); {
    *cutflow << HFTname("nCentralBJets");
    *cutflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfCBJets(*sl->baseJets)/*(*baseJets)*/; };
    *cutflow << SaveVar();
  }
  //xTauFW variable
  *cutflow << NewVar("b-tagged jet"); {
    *cutflow << HFTname("Btag");
    *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->tools->numberOfCBJets(*sl->baseJets) > 0;};
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("number of forward jets"); {
    *cutflow << HFTname("nForwardJets");
    *cutflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfFJets(*sl->baseJets)/*(*baseJets)*/; };
    *cutflow << SaveVar();
  }
  //////////////////////////////////////////////////////////////////////////////
  // Pre Jets
  *cutflow << NewVar("preJet pt"); {
    *cutflow << HFTname("preJet_pt");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for (auto& jet : preJets) {out.push_back(jet->Pt()); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("preJet eta"); {
    *cutflow << HFTname("preJet_eta");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for (auto& jet : preJets) {out.push_back(jet->Eta()); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("preJet JVT if eta<=2.4 and pT<60"); {
    *cutflow << HFTname("preJet_JVT");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for (auto& jet : preJets) {
        if (jet->Pt() < 60 && fabs(jet->Eta()) <= 2.4) out.push_back(jet->jvt);
        else
          out.push_back(0);
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
      for (auto& jet : baseJets) {out.push_back(jet->Eta()); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Baseline Jet mv2c10"); {
    *cutflow << HFTname("baseJet_mv2c10");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for (auto& jet : baseJets) {out.push_back(jet->mv2c10); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("jet eta"); {
    *cutflow << HFTname("j_eta");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& jet : baseJets) {
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
      for(auto& jet : baseJets) {
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
      for(auto& jet : baseJets) {
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
      for(auto& jet : baseJets) {
        out.push_back(jet->Phi());
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("jet flavor (0: NA, 1: CL, 2: CB, 3: F)"); {
    *cutflow << HFTname("j_flav");
    *cutflow << [&](Superlink* sl, var_int_array*) -> vector<int> {
      vector<int> out; int flav = 0;
      for(auto& jet : baseJets) {
        if(sl->tools->m_jetSelector->isCentralLight(jet))  { flav = 1; }
        else if(sl->tools->m_jetSelector->isCentralB(jet)) { flav = 2; }
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
        for (auto& jet : baseJets) {
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
      for(auto& jet : baseJets) {
        out.push_back(jet->Pt());
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("dijet mass"); {
    *cutflow << HFTname("Mjj");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
        if (baseJets.size() > 1) {
            return JetP4.M();
        } else {
            return -1.0;
        }
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("DeltaEta between two leading jets"); {
    *cutflow << HFTname("DEtaJJ");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
        if (baseJets.size() > 1) {
            return fabs(Jet0.Eta() - Jet1.Eta());
        } else {
            return -1.0;
        }
    };
    *cutflow << SaveVar();
  }


  ///////////////////////////////////////////////////////////////////////
  // Clear internal variables
  ///////////////////////////////////////////////////////////////////////
  *cutflow << [&](Superlink* /*sl*/, var_void*) {
    preLeptons.clear(); baseLeptons.clear(); signalLeptons.clear();
    preElectrons.clear(); baseElectrons.clear(); signalElectrons.clear();
    preMuons.clear(); baseMuons.clear(); signalMuons.clear();
    preTaus.clear(); baseTaus.clear(); signalTaus.clear();
    preJets.clear(); baseJets.clear(); signalJets.clear();
    centralLightJets.clear(); centralBJets.clear(); forwardJets.clear();
  };

  ///////////////////////////////////////////////////////////////////////
  //
  // Systematics
  //
  ///////////////////////////////////////////////////////////////////////

  // Here comes the systematics...

  ///////////////////////////////////////////////////////////////////////
  //
  // SUPERFLOW METHODS END HERE
  //
  ///////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  // Initialize the cutflow and start the event loop.
  chain->Process(cutflow, input_file, n_events, n_skip_events);

  delete cutflow;
  delete chain;

  ///////////////////////////////////////////////////////////////////////
  // Print information
  printf("makeMiniNtuple\t =================================================================\n");
  printf("makeMiniNtuple\t All done!\n");
  printf("makeMiniNtuple\t =================================================================\n");
  return 0;
}
