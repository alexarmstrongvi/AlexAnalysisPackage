//////////////////////////////////////////////////////////////////////////////////////////
// \project:    ATLAS Experiment at CERN's LHC
// \package:    MyAnalysis
// \file:       $Id$
// \author:     Alaettin.Serhan.Mete@cern.ch
// \history:    N/A
//
// Copyright (C) 2016 University of California, Irvine
//////////////////////////////////////////////////////////////////////////////////////////

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

bool passTrigger(Superlink* sl, std::string trig_name) {
  return sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, trig_name);
}

///////////////////////////////////////////////////////////////////////
// Main function
///////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  ///////////////////////////////////////////////////////////////////////
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
  ///////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
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

  ///////////////////////////////////////////////////////////////////////
  // Superflow methods begin here
  ///////////////////////////////////////////////////////////////////////

  *cutflow << CutName("read in") << [](Superlink* /*sl*/) -> bool { return true; };

  //  Cleaning Cuts
  // How to add cutflow entry:
  // Create lambda function that returns bool value of cut.
  // Pass that with "<<" into the function CutName("Cut name").
  // Pass that with "<<" into the dereferenced cutflow object.

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
  *cutflow << CutName("pass good vertex") << [&](Superlink* sl) -> bool {
      return (sl->tools->passGoodVtx(cutflags));
  };
  *cutflow << CutName("jet cleaning") << [&](Superlink* sl) -> bool {
      return (sl->tools->passJetCleaning(sl->baseJets));
  };
  *cutflow << CutName("bad muon veto") << [&](Superlink* sl) -> bool {
      return (sl->tools->passBadMuon(sl->preMuons));
  };


  //  Output Ntuple Setup
  //      > Ntuple variables

  // Storing a bunch of numbers...
  *cutflow << NewVar("Event run number"); {
    *cutflow << HFTname("RunNumber");
    *cutflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->run; };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("Event number"); {
    *cutflow << HFTname("EventNumber");
    *cutflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->eventNumber; };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("is Monte Carlo"); {
    *cutflow << HFTname("isMC");
    *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->nt->evt()->isMC ? true : false; };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("event weight"); {
    *cutflow << HFTname("event_weight");
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

  // Dipleton tiggers
  *cutflow << NewVar("HLT_mu18_mu8noL1 trigger bit"); {
      *cutflow << HFTname("HLT_mu18_mu8noL1");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
        return passTrigger(sl, "HLT_mu18_mu8noL1");
      };
      *cutflow << SaveVar();
  }

  *cutflow << NewVar("HLT_2e12_lhloose_L12EM10VH trigger bit"); {
      *cutflow << HFTname("HLT_2e12_lhloose_L12EM10VH");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
        return passTrigger(sl, "HLT_2e12_lhloose_L12EM10VH");
      };
      *cutflow << SaveVar();
  }

  *cutflow << NewVar("HLT_e17_lhloose_mu14 trigger bit"); {
      *cutflow << HFTname("HLT_e17_lhloose_mu14");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
        return passTrigger(sl, "HLT_e17_lhloose_mu14");
      };
      *cutflow << SaveVar();
  }

  *cutflow << NewVar("HLT_mu20_mu8noL1 trigger bit"); {
      *cutflow << HFTname("HLT_mu20_mu8noL1");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
        return passTrigger(sl, "HLT_mu20_mu8noL1");
      };
      *cutflow << SaveVar();
  }

  *cutflow << NewVar("HLT_2e15_lhvloose_L12EM13VH trigger bit"); {
      *cutflow << HFTname("HLT_2e15_lhvloose_L12EM13VH");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
        return passTrigger(sl, "HLT_2e15_lhvloose_L12EM13VH");
      };
      *cutflow << SaveVar();
  }

  *cutflow << NewVar("HLT_2e17_lhvloose_nod0 trigger bit"); {
      *cutflow << HFTname("HLT_2e17_lhvloose_nod0");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
        return passTrigger(sl, "HLT_2e17_lhvloose_nod0");
      };
      *cutflow << SaveVar();
  }

  *cutflow << NewVar("HLT_mu22_mu8noL1 trigger bit"); {
      *cutflow << HFTname("HLT_mu22_mu8noL1");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
        return passTrigger(sl, "HLT_mu22_mu8noL1");
      };
      *cutflow << SaveVar();
  }

  *cutflow << NewVar("HLT_e17_lhloose_nod0_mu14 trigger bit"); {
      *cutflow << HFTname("HLT_e17_lhloose_nod0_mu14");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
        return passTrigger(sl, "HLT_e17_lhloose_nod0_mu14");
      };
      *cutflow << SaveVar();
  }

   // Single Lepton Triggers
  *cutflow << NewVar("HLT_e60_lhmedium trigger bit"); {
      *cutflow << HFTname("HLT_e60_lhmedium");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
          return passTrigger(sl, "HLT_e60_lhmedium");
        };
      *cutflow << SaveVar();
  }
  *cutflow << NewVar("HLT_e24_lhmedium_L1EM20VH trigger bit"); {
      *cutflow << HFTname("HLT_e24_lhmedium_L1EM20VH");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
          return passTrigger(sl, "HLT_e24_lhmedium_L1EM20VH");
        };
      *cutflow << SaveVar();
  }
  *cutflow << NewVar("HLT_e24_lhmedium_iloose_L1EM18VH trigger bit"); {
      *cutflow << HFTname("HLT_e24_lhmedium_iloose_L1EM18VH");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
          return passTrigger(sl, "HLT_e24_lhmedium_iloose_L1EM18VH");
        };
      *cutflow << SaveVar();
  }
  *cutflow << NewVar("HLT_mu20_iloose_L1MU15 trigger bit"); {
      *cutflow << HFTname("HLT_mu20_iloose_L1MU15");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
          return passTrigger(sl, "HLT_mu20_iloose_L1MU15");
        };
      *cutflow << SaveVar();
  }
  *cutflow << NewVar("HLT_mu24_imedium trigger bit"); {
      *cutflow << HFTname("HLT_mu24_imedium");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
          return passTrigger(sl, "HLT_mu24_imedium");
        };
      *cutflow << SaveVar();
  }
  *cutflow << NewVar("HLT_mu26_imedium trigger bit"); {
      *cutflow << HFTname("HLT_mu26_imedium");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
          return passTrigger(sl, "HLT_mu26_imedium");
        };
      *cutflow << SaveVar();
  }
  *cutflow << NewVar("HLT_e24_lhtight_nod0_ivarloose trigger bit"); {
      *cutflow << HFTname("HLT_e24_lhtight_nod0_ivarloose");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
          return passTrigger(sl, "HLT_e24_lhtight_nod0_ivarloose");
        };
      *cutflow << SaveVar();
  }
  *cutflow << NewVar("HLT_e26_lhtight_nod0_ivarloose trigger bit"); {
      *cutflow << HFTname("HLT_e26_lhtight_nod0_ivarloose");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
          return passTrigger(sl, "HLT_e26_lhtight_nod0_ivarloose");
        };
      *cutflow << SaveVar();
  }
  *cutflow << NewVar("HLT_e60_lhmedium_nod0 trigger bit"); {
      *cutflow << HFTname("HLT_e60_lhmedium_nod0");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
          return passTrigger(sl, "HLT_e60_lhmedium_nod0");
        };
      *cutflow << SaveVar();
  }
  *cutflow << NewVar("HLT_e120_lhloose trigger bit"); {
      *cutflow << HFTname("HLT_e120_lhloose");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
          return passTrigger(sl, "HLT_e120_lhloose");
        };
      *cutflow << SaveVar();
  }
  *cutflow << NewVar("HLT_e140_lhloose_nod0 trigger bit"); {
      *cutflow << HFTname("HLT_e140_lhloose_nod0");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
          return passTrigger(sl, "HLT_e140_lhloose_nod0");
        };
      *cutflow << SaveVar();
  }
  *cutflow << NewVar("HLT_mu24_iloose trigger bit"); {
      *cutflow << HFTname("HLT_mu24_iloose");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
          return passTrigger(sl, "HLT_mu24_iloose");
        };
      *cutflow << SaveVar();
  }
  *cutflow << NewVar("HLT_mu24_iloose_L1MU15 trigger bit"); {
      *cutflow << HFTname("HLT_mu24_iloose_L1MU15");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
          return passTrigger(sl, "HLT_mu24_iloose_L1MU15");
        };
      *cutflow << SaveVar();
  }
  //*cutflow << NewVar("HLT_mu24_ivarloose trigger bit"); {
  //    *cutflow << HFTname("HLT_mu24_ivarloose");
  //    *cutflow << [](Superlink* sl, var_bool*) -> bool {
  //        return passTrigger(sl, "HLT_mu24_ivarloose");
  //    };
  //    *cutflow << SaveVar();
  //}
  //*cutflow << NewVar("HLT_mu24_ivarloose_L1MU15 trigger bit"); {
  //    *cutflow << HFTname("HLT_mu24_ivarloose_L1MU15");
  //    *cutflow << [](Superlink* sl, var_bool*) -> bool {
  //        return passTrigger(sl, "HLT_mu24_ivarloose_L1MU15");
  //    };
  //    *cutflow << SaveVar();
  //}
  //*cutflow << NewVar("HLT_mu24_ivarmedium trigger bit"); {
  //    *cutflow << HFTname("HLT_mu24_ivarmedium");
  //    *cutflow << [](Superlink* sl, var_bool*) -> bool {
  //        return passTrigger(sl, "HLT_mu24_ivarmedium");
  //    };
  //    *cutflow << SaveVar();
  //}
  *cutflow << NewVar("HLT_mu26_ivarmedium trigger bit"); {
      *cutflow << HFTname("HLT_mu26_ivarmedium");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
          return passTrigger(sl, "HLT_mu26_ivarmedium");
        };
      *cutflow << SaveVar();
  }
  *cutflow << NewVar("HLT_mu40"); {
      *cutflow << HFTname("HLT_mu40");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
          return passTrigger(sl, "HLT_mu40");
        };
      *cutflow << SaveVar();
  }
  *cutflow << NewVar("HLT_mu50"); {
      *cutflow << HFTname("HLT_mu50");
      *cutflow << [](Superlink* sl, var_bool*) -> bool {
          return passTrigger(sl, "HLT_mu50");
        };
      *cutflow << SaveVar();
  }

  // 15/16 Year ID
  *cutflow << NewVar("treatAsYear"); {
      *cutflow << HFTname("treatAsYear");
      *cutflow << [](Superlink* sl, var_double*) -> int {
        return sl->nt->evt()->treatAsYear;
      };
      *cutflow << SaveVar();
  }


  // Accesing objects through superlink and assigning them to internal variables
  LeptonVector preLeptons, baseLeptons, signalLeptons;
  ElectronVector signalElectrons;
  Susy::Met met;
  *cutflow << [&](Superlink* sl, var_void*) {
    preLeptons = *sl->preLeptons;
    baseLeptons   = *sl->baseLeptons;
    signalLeptons = *sl->leptons;
    signalElectrons = *sl->electrons;
    met = *sl->met;
  };

  *cutflow << NewVar("number of pre leptons"); {
    *cutflow << HFTname("nPreLeptons");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return preLeptons.size(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of baseline leptons"); {
    *cutflow << HFTname("nBaselineLeptons");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return baseLeptons.size(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of signal leptons"); {
    *cutflow << HFTname("nSignalLeptons");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return signalLeptons.size(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of baseline taus"); {
    *cutflow << HFTname("nBaseTaus");
    *cutflow << [&](Superlink* sl, var_int*) -> int { return sl->baseTaus->size(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of signal taus"); {
    *cutflow << HFTname("nSignalTaus");
    *cutflow << [&](Superlink* sl, var_int*) -> int { return sl->taus->size(); };
    *cutflow << SaveVar();
  }

  // Storing a bunch of vectors
  uint lepID_n = 0, lepAntiID_n = 0;
  vector<int> lepID_vec, lepAntiID_vec; 
  // ID Leptons
  *cutflow << [&](Superlink* /*sl*/, var_void*) {
    for(auto& lepton : preLeptons) {
      bool pt_pass = 0, eta_pass = 0, iso_pass = 0, id_pass = 0, OR_pass = 0;
      if (lepton->isEle()) {
        const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(lepton);
        float absEtaBE = fabs(ele->clusEtaBE);
        pt_pass  = ele->pt > 15;
        eta_pass = absEtaBE < 1.37 || 1.52 < absEtaBE;
        OR_pass  = ele->isBaseline; // TODO: Improve
        iso_pass = ele->isoGradient;
        id_pass  = ele->mediumLLH;
      } else if (lepton->isMu()) {
        const Susy::Muon* mu = dynamic_cast<const Susy::Muon*>(lepton);
        pt_pass  = mu->pt > 10;
        eta_pass = 1;
        OR_pass  = mu->isBaseline; // TODO: Improve
        iso_pass = mu->isoGradient;
        id_pass  = mu->medium;
      }
      bool lepID_bool = pt_pass && eta_pass && iso_pass && id_pass && OR_pass;
      lepID_vec.push_back(lepID_bool);
      if(lepID_bool) lepID_n++;
    }
  };
  // Anti-ID Leptons
  *cutflow << [&](Superlink* /*sl*/, var_void*) {
    for(auto& lepton : preLeptons) {
      bool pt_pass = 0, eta_pass = 0, iso_pass = 0, id_pass = 0, OR_pass = 0;
      bool passedID_cuts = 0, passedAntiID_cuts = 0;
      if (lepton->isEle()) {
        const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(lepton);
        float absEtaBE = fabs(ele->clusEtaBE);
        pt_pass  = ele->pt > 15;
        eta_pass = absEtaBE < 1.37 || 1.52 < absEtaBE;
        OR_pass  = ele->isBaseline; // TODO: Improve  
        iso_pass = ele->isoGradient;
        id_pass  = ele->mediumLLH;
        passedID_cuts = iso_pass && id_pass; 
        passedAntiID_cuts = ele->veryLooseLLH;
      } else if (lepton->isMu()) {
        const Susy::Muon* mu = dynamic_cast<const Susy::Muon*>(lepton);
        pt_pass  = mu->pt > 10;
        eta_pass = 1;
        OR_pass  = 1; // TODO: Improve
        iso_pass = mu->isoGradient;
        id_pass  = mu->medium;
        passedID_cuts = iso_pass && id_pass; 
        passedAntiID_cuts = 1;
      }

      bool lepAntiID_bool = pt_pass && eta_pass && OR_pass
                        && passedAntiID_cuts && !passedID_cuts;
      lepAntiID_vec.push_back(lepAntiID_bool);
      if(lepAntiID_bool) lepAntiID_n++;
    }
  };

  *cutflow << NewVar("ID Leptons"); {
    *cutflow << HFTname("lepID");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {return lepID_vec;};
    *cutflow << SaveVar();
  };
  *cutflow << NewVar("Anti ID Leptons"); {
    *cutflow << HFTname("lepAntiID");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {return lepAntiID_vec;};
    *cutflow << SaveVar();
  };
  *cutflow << NewVar("Lepton ID status (ID(1), antiID(-1), Other(0))"); {
    *cutflow << HFTname("lepID_status");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
        vector<int> out;
        for (uint i = 0; i < preLeptons.size(); i++) {    
            // ID -> 1, AntiID -> -1, Other -> 0
            if (lepID_vec.at(i)) out.push_back(1);
            else if (lepAntiID_vec.at(i)) out.push_back(-1);
            else out.push_back(0);
        }
        return out;
    };
    *cutflow << SaveVar();
  };
  *cutflow << NewVar("Number of ID leptons"); {
    *cutflow << HFTname("nLepID");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return lepID_n; };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Number of Anti-ID leptons"); {
    *cutflow << HFTname("nLepAntiID");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return lepAntiID_n; };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("lepton flavor (1: e, 2: m)"); {
    *cutflow << HFTname("lepFlavor");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : preLeptons) {
        out.push_back(lepton->isEle() ? 1 : 2);
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("lepton charge"); {
    *cutflow << HFTname("lepCharge");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : preLeptons) {
        out.push_back(lepton->q);
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("lepton passes overlap removal"); {
    *cutflow << HFTname("lepPassOR");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : preLeptons) {
        out.push_back(lepton->isBaseline); //TODO: Improve
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("lepTruthDetailed"); {
    *cutflow << HFTname("lepTruthDetailed");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : preLeptons) {
        out.push_back(lepton->matched2TruthLepton);
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("lepMatchesTrigger"); {
    *cutflow << HFTname("lepMatchesTrigger");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : preLeptons) {
        (void)lepton;
        out.push_back(1); // TODO: Add
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Gradient"); {
    *cutflow << HFTname("Gradient");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : preLeptons) {
        out.push_back(lepton->isoGradient);
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("GradientLoose"); {
    *cutflow << HFTname("GradientLoose");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : preLeptons) {
        out.push_back(lepton->isoGradientLoose);
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("LooseTrackOnly"); {
    *cutflow << HFTname("LooseTrackOnly");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : preLeptons) {
        out.push_back(lepton->isoLooseTrackOnly);
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("lepEleEtaBE"); {
    *cutflow << HFTname("lepEleEtaBE");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : preLeptons) {
        if (lepton->isEle()) {
          const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(lepton);
          out.push_back(ele->clusEtaBE);
        } else {
          out.push_back(0);
        }
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("TightLH"); {
    *cutflow << HFTname("TightLH");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : preLeptons) {
        if (lepton->isEle()) {
          const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(lepton);
          out.push_back(ele->tightLLH);
        } else {
          out.push_back(0);
        }
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("MediumLH"); {
    *cutflow << HFTname("MediumLH");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : preLeptons) {
        if (lepton->isEle()) {
          const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(lepton);
          out.push_back(ele->mediumLLH);
        } else {
          out.push_back(0);
        }
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("LooseLH"); {
    *cutflow << HFTname("LooseLH");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : preLeptons) {
        if (lepton->isEle()) {
          const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(lepton);
          out.push_back(ele->looseLLH);
        } else {
          out.push_back(0);
        }
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("LooseAndBLayerLH."); {
    *cutflow << HFTname("LooseAndBLayerLH.");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : preLeptons) {
        if (lepton->isEle()) {
          const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(lepton);
          out.push_back(ele->looseLLHBLayer);
        } else {
          out.push_back(0);
        }
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("VeryLooseLH"); {
    *cutflow << HFTname("VeryLooseLH");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : preLeptons) {
        if (lepton->isEle()) {
          const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(lepton);
          out.push_back(ele->veryLooseLLH);
        } else {
          out.push_back(0);
        }
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Tight"); {
    *cutflow << HFTname("Tight");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : preLeptons) {
        if (lepton->isMu()) {
          const Susy::Muon* mu = dynamic_cast<const Susy::Muon*>(lepton);
          out.push_back(mu->tight);
        } else {
          out.push_back(0);
        }
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Medium"); {
    *cutflow << HFTname("Medium");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : preLeptons) {
        if (lepton->isMu()) {
          const Susy::Muon* mu = dynamic_cast<const Susy::Muon*>(lepton);
          out.push_back(mu->medium);
        } else {
          out.push_back(0);
        }
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Loose"); {
    *cutflow << HFTname("Loose");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : preLeptons) {
        if (lepton->isMu()) {
          const Susy::Muon* mu = dynamic_cast<const Susy::Muon*>(lepton);
          out.push_back(mu->loose);
        } else {
          out.push_back(0);
        }
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("VeryLoose"); {
    *cutflow << HFTname("VeryLoose");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : preLeptons) {
        if (lepton->isMu()) {
          const Susy::Muon* mu = dynamic_cast<const Susy::Muon*>(lepton);
          out.push_back(mu->veryLoose);
        } else {
          out.push_back(0);
        }
      }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("lepton pt"); {
    *cutflow << HFTname("lepPt");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : preLeptons) {
        out.push_back(lepton->Pt());
      }
      return out;
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("lepton eta"); {
    *cutflow << HFTname("lepEta");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : preLeptons) {
        out.push_back(lepton->Eta());
      }
      return out;
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("lepton phi"); {
    *cutflow << HFTname("lepPhi");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : preLeptons) {
        out.push_back(lepton->Phi());
      }
      return out;
    };
    *cutflow << SaveVar();
  }


  *cutflow << NewVar("lepton type"); {
    *cutflow << HFTname("lepType");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for(auto& lepton : preLeptons) {
        out.push_back(lepton->mcType);
      }
      return out;
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("lepton origin"); {
    *cutflow << HFTname("lepOrigin");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : preLeptons) {
        out.push_back(lepton->mcOrigin);
      }
      return out;
    };
    *cutflow << SaveVar();
  }

  // MET
  TLorentzVector MET;
  *cutflow << NewVar("transverse missing energy (Et)"); {
    *cutflow << HFTname("met_Et");
    *cutflow << [&](Superlink* sl, var_float*) -> float {
      MET.SetPxPyPzE(sl->met->Et*cos(sl->met->phi),
                     sl->met->Et*sin(sl->met->phi),
                     0.,
                     sl->met->Et);
      return MET.Pt();
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("transverse missing energy (Phi)"); {
    *cutflow << HFTname("met_phi");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> float { return MET.Phi(); };
    *cutflow << SaveVar();
  }
  // Interesting Leptons
  TLorentzVector lepton0, lepton1, dileptonP4;
  int lepton0_idx = -1, lepton1_idx = -1; 
  *cutflow << [&](Superlink* /*sl*/, var_void*) {
    float lep_pt0 = 0, lep_pt1 = 0;
    for(uint i = 0; i < preLeptons.size(); i++) {
       if (!(lepID_vec.at(i) || lepAntiID_vec.at(i))) continue;
       auto& lepton = preLeptons.at(i);
       double pt = lepton->pt;
       if (pt < lep_pt1) continue;

       // Build TLVector
       double eta = lepton->eta;
       double phi = lepton->phi;
       // Determine lepton mass in GeV
       double mass = lepton->isEle() ? 0.000511 : 0.105658;
       TLorentzVector lep; lep.SetPtEtaPhiM(pt, eta, phi, mass);
       // Update
       if (lep.Pt() > lep_pt0) {
           lepton0_idx = i;
           lep_pt0 = lep.Pt();
           lepton0 = lep;
       } else if (lep.Pt() > lep_pt1) {
           lepton1_idx = i;
           lep_pt1 = lep.Pt();
           lepton1 = lep;
       }
       dileptonP4 = lepton0 + lepton1;
    } 
  };

  *cutflow << NewVar("Leading lepton pT"); {
    *cutflow << HFTname("Lep0Pt");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return lepton0.Pt(); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Subleading lepton pT"); {
    *cutflow << HFTname("Lep1Pt");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return lepton1.Pt(); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Leading lepton eta"); {
    *cutflow << HFTname("Lep0Eta");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return lepton0.Eta(); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Subleading lepton eta"); {
    *cutflow << HFTname("Lep1Eta");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return lepton1.Eta(); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Leading lepton phi"); {
    *cutflow << HFTname("Lep0Phi");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return lepton0.Phi(); };
    *cutflow << SaveVar();
  }
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
        if (lepton0_idx < 0 || lepton1_idx < 0) {return -DBL_MAX;}
       return fabs(lepton0.Eta() - lepton1.Eta()); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Delta phi between leptons"); {
    *cutflow << HFTname("DphiLL");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
        if (lepton0_idx < 0 || lepton1_idx < 0) {return -DBL_MAX;}
       return fabs(lepton0.DeltaPhi(lepton1));};
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("delta R of di-lepton system"); {
    *cutflow << HFTname("drll");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { 
        if (lepton0_idx < 0 || lepton1_idx < 0) {return -DBL_MAX;}
        return lepton0.DeltaR(lepton1); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("collinear mass, M_coll"); {
    *cutflow << HFTname("MCollASym");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
        if (lepton0_idx < 0 || lepton1_idx < 0) {return -DBL_MAX;}
        double deta = fabs(lepton0.Eta()-lepton1.Eta());
        double dphi = lepton0.DeltaPhi(lepton1);
        return sqrt(2*lepton0.Pt()*(lepton1.Pt()+met.Et)*(cosh(deta)-cos(dphi)));
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("transverse mass of leading lepton"); {
    *cutflow << HFTname("MtLep0");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
        if (lepton0_idx < 0 || lepton1_idx < 0) {return -DBL_MAX;}
        return sqrt( 2*lepton0.Pt()*MET.Pt()*(1-cos (lepton0.DeltaPhi(MET)) ) );
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("transverse mass of subleading lepton"); {
    *cutflow << HFTname("MtLep1");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
        if (lepton0_idx < 0 || lepton1_idx < 0) {return -DBL_MAX;}
        return sqrt( 2*lepton1.Pt()*MET.Pt()*(1-cos (lepton1.DeltaPhi(MET)) ) );
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("mass of di-lepton system, M_ll"); {
    *cutflow << HFTname("MLL");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { 
        if (lepton0_idx < 0 || lepton1_idx < 0) {return -DBL_MAX;}
        return dileptonP4.M(); 
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("Pt of di-lepton system, Pt_ll"); {
    *cutflow << HFTname("ptll");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { 
        if (lepton0_idx < 0 || lepton1_idx < 0) {return -DBL_MAX;}
        return dileptonP4.Pt(); 
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("diff_pt between leading and sub-leading lepton"); {
    *cutflow << HFTname("dpt_ll");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
        if (lepton0_idx < 0 || lepton1_idx < 0) {return -DBL_MAX;}
        return lepton0.Pt() - lepton1.Pt();
    };
    *cutflow << SaveVar();
 }

  *cutflow << NewVar("dphi between leading lepton and MET"); {
    *cutflow << HFTname("DphiLep0MET");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
        if (lepton0_idx < 0) {return -DBL_MAX;}
        return lepton0.DeltaPhi(MET);
    };
    *cutflow << SaveVar();
 }

  *cutflow << NewVar("dphi between sub-leading lepton and MET"); {
    *cutflow << HFTname("DphiLep1MET");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
        if (lepton1_idx < 0) {return -DBL_MAX;}
        return lepton1.DeltaPhi(MET);
    };
    *cutflow << SaveVar();
 }
  *cutflow << NewVar("dilepton flavor"); {
    *cutflow << HFTname("dilep_flav");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if (lepton0_idx < 0 || lepton1_idx < 0) {return -1;}

        if(preLeptons.at(lepton0_idx)->isEle() && preLeptons.at(lepton1_idx)->isMu()){return 0;}       // e mu  case
        else if(preLeptons.at(lepton0_idx)->isMu() && preLeptons.at(lepton1_idx)->isEle()){return 1;}  // mu e  case
        else if(preLeptons.at(lepton0_idx)->isEle() && preLeptons.at(lepton1_idx)->isEle()){return 2;} // e e   case
        else if(preLeptons.at(lepton0_idx)->isMu() && preLeptons.at(lepton1_idx)->isMu()){return 3;}   // mu mu case
	    else{return 4;}
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("dilepton flavor is e mu"); {
    *cutflow << HFTname("isEM");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if (lepton0_idx < 0 || lepton1_idx < 0) return 0;
	    return preLeptons.at(lepton0_idx)->isEle() && preLeptons.at(lepton1_idx)->isMu();  // e mu  case
        };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("dilepton flavor is mu e"); {
    *cutflow << HFTname("isME");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if (lepton0_idx < 0 || lepton1_idx < 0) return 0;
        return preLeptons.at(lepton0_idx)->isMu() && preLeptons.at(lepton1_idx)->isEle();  // mu e  case
        };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Dilepton Charge Sign -1 (OS), 1 (SS), 0 (N/A) "); {
    *cutflow << HFTname("qq");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if (lepton0_idx < 0 || lepton1_idx < 0) return 0;
	    return preLeptons.at(lepton0_idx)->q*preLeptons.at(lepton1_idx)->q;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("ID channel: 0-N/A, 1-tt, 2-lt, 3-tl, 4-ll "); {
    *cutflow << HFTname("dilepID");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if (lepton0_idx < 0 || lepton1_idx < 0) return 0;
        else if (lepID_vec.at(lepton0_idx) && lepID_vec.at(lepton1_idx)) return 1;
        else if (lepAntiID_vec.at(lepton0_idx) && lepID_vec.at(lepton1_idx)) return 2;
        else if (lepID_vec.at(lepton0_idx) && lepAntiID_vec.at(lepton1_idx)) return 3;
        else if (lepAntiID_vec.at(lepton0_idx) && lepAntiID_vec.at(lepton1_idx)) return 4;
        else return 0;
    };
    *cutflow << SaveVar();
  }

  // Jet Variables
  JetVector baseJets, signalJets, lightJets, BJets, forwardJets;
  TLorentzVector JetP4, Jet1, Jet0;

  *cutflow << [&](Superlink* sl, var_void*) {
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
          if (sl->tools->m_jetSelector->isLight(jet))  {
               lightJets.push_back(jet);
          } else if (sl->tools->m_jetSelector->isB(jet)) {
              BJets.push_back(jet);
          } else if (sl->tools->m_jetSelector->isForward(jet))  {
              forwardJets.push_back(jet);
          }
    }
    std::sort(lightJets.begin()  , lightJets.end()  , comparePt);
    std::sort(BJets.begin()      , BJets.end()      , comparePt);
    std::sort(forwardJets.begin()       , forwardJets.end()       , comparePt);
  };

  *cutflow << NewVar("number of baseline jets"); {
    *cutflow << HFTname("JetN");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return baseJets.size(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of select baseline jets"); {
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

  *cutflow << NewVar("number of light jets"); {
    *cutflow << HFTname("nLJets");
    *cutflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfLJets(*sl->baseJets)/*(*baseJets)*/; };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of B jets"); {
    *cutflow << HFTname("nBJets");
    *cutflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfBJets(*sl->baseJets)/*(*baseJets)*/; };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("b-tagged jet"); {
    *cutflow << HFTname("Btag");
    *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->tools->numberOfBJets(*sl->baseJets) > 0;};
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of forward jets"); {
    *cutflow << HFTname("nForwardJets");
    *cutflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfFJets(*sl->baseJets)/*(*baseJets)*/; };
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

  // *cutflow << NewVar("jet JVF"); {
  //   *cutflow << HFTname("j_jvf");
  //   *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
  //     vector<double> out;
  //     for(auto& jet : baseJets) {
  //       out.push_back(jet->jvf);
  //     }
  //     return out;
  //   };
  //   *cutflow << SaveVar();
  // }

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

  *cutflow << NewVar("jet flavor (0: NA, 1: L, 2: B, 3: F)"); {
    *cutflow << HFTname("j_flav");
    *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
      vector<double> out; int flav = 0;
      for(auto& jet : baseJets) {
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


  // Clear internal variables
  *cutflow << [&](Superlink* /*sl*/, var_void*) {
    preLeptons.clear(); baseLeptons.clear();signalLeptons.clear();
    baseJets.clear();signalJets.clear();
    lightJets.clear();BJets.clear();forwardJets.clear();
    lepID_n = 0; lepAntiID_n = 0;
    lepton0_idx = -1, lepton1_idx = -1; 
    lepID_vec.clear(); lepAntiID_vec.clear(); 

  };

  ///////////////////////////////////////////////////////////////////////
  // Systematics
  ///////////////////////////////////////////////////////////////////////

  // Here comes the systematics...

  ///////////////////////////////////////////////////////////////////////
  // Superflow methods end here
  ///////////////////////////////////////////////////////////////////////

  // Initialize the cutflow and start the event loop.
  chain->Process(cutflow, input_file, n_events, n_skip_events);

  delete cutflow;
  delete chain;

  // Print information
  printf("makeMiniNtuple\t =================================================================\n");
  printf("makeMiniNtuple\t All done!\n");
  printf("makeMiniNtuple\t =================================================================\n");
  return 0;
}
