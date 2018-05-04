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
  printf("-s        Append suffix to output filename\n");
  printf("-i        Input file as *.root, list of *.root in a *.txt,\n");
  printf("          or a DIR/ containing *.root (default: none)\n");
  printf("-f        make fake samples [zjets, wjets]\n");
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
  char *do_fakes = nullptr;
  SuperflowRunMode run_mode = SuperflowRunMode::nominal; // SuperflowRunMode::all_syst; //SuperflowRunMode::nominal;
  int c;

  opterr = 0;
  while ((c = getopt (argc, argv, "i:s:n:f:h")) != -1)
    switch (c)
      {
      case 'i':
        input_file = optarg;
        break;
      case 's':
        name_suffix = optarg;
        break;
      case 'n':
        n_events = atoi(optarg);
        break;
      case 'f':
        do_fakes = optarg;        
        break;
      case 'h':
        usage("makeMiniNtuple");
        return 1;
      case '?':
        if (optopt == 'i')
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
    printf("makeMiniNtuple\t An input file must be provided with option -i (a list, a DIR or single file)\n");
    return 0;
  }
  std::string fake_op = do_fakes ? do_fakes : " - ";
  bool do_fakes_zjets = false, do_fakes_wjets = false;
  if (fake_op == " - ");
  else if (fake_op == "zjets") do_fakes_zjets = true;
  else if (fake_op == "wjets") do_fakes_wjets = true;
  else {
    cout << "makeMiniNtuple\t Unrecognized fake options: " << fake_op << '\n';
    return 0;
  }

  // Print information
  printf("makeMiniNtuple\t =================================================================\n");
  printf("makeMiniNtuple\t Running MyAnalysis/makeMiniNtuple\n");
  printf("makeMiniNtuple\t =================================================================\n");
  printf("makeMiniNtuple\t   Flags:\n");
  printf("makeMiniNtuple\t     Input file (-i)         : %s\n",input_file);
  printf("makeMiniNtuple\t     Number of events (-n)   : %i\n",n_events );
  printf("makeMiniNtuple\t     Adding fakes (-f)       : %s\n",do_fakes );
  printf("makeMiniNtuple\t     Appending suffix (-s)   : %s\n",name_suffix );
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
  // Initialize relevant fake variables. Set inside cutflow 
  vector<int> lepID_vec, lepAntiID_vec;
  uint lepID_n = 0, lepAntiID_n = 0;
  int antiID_idx0 = -1, antiID_idx1 = -1;
  // Find leptons paired with Z
  // Indices 0 and 1 are closest Z pair
  // Index 2 is leading anti-ID or 3rd leading ID lepton
  // Indices 3 and 4 are second closest Z pair
  LeptonVector Zlep(5, nullptr);
  TLorentzVector LepFake0, LepFake1;


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
  if (!do_fakes) {
    *cutflow << CutName("xTau: 2+ Loose Leptons") << [&](Superlink* sl) -> bool {
        uint nLooseLeptons = 0;
        for (const auto* mu : *sl->preMuons) {if (mu->loose) nLooseLeptons++;}
        for (const auto* ele : *sl->preElectrons) {if (ele->looseLLH) nLooseLeptons++;}
        return nLooseLeptons >= 2;
    };
  }

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
  if (do_fakes) {
    //////////////////////////////////////////////////////////////////////////////
    // Set relevant fake variables. This cut does not reject any events
    *cutflow << CutName("NoCut (Define Fake Vars)") << [&](Superlink* sl) -> bool {
      // ID Leptons
      for(auto& lepton : *sl->baseLeptons) {
        bool lepID_bool = false;
        for (auto& sig_lepton : *sl->leptons) {
            if (sig_lepton == lepton) {
                lepID_bool = true;
                break;
            }
        }
        lepID_vec.push_back(lepID_bool);
        if(lepID_bool) lepID_n++;
      }
      // Anti-ID Leptons
      for(auto& lepton : *sl->baseLeptons) {
        bool pt_pass = 0, eta_pass = 0, iso_pass = 0, id_pass = 0;
        bool passedID_cuts = 0, passedAntiID_cuts = 0;
        if (lepton->isEle()) {
          const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(lepton);
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
          passedID_cuts = iso_pass && id_pass;
          passedAntiID_cuts = 1;
        }
        bool lepAntiID_bool = pt_pass && eta_pass
                              && passedAntiID_cuts && !passedID_cuts;
        lepAntiID_vec.push_back(lepAntiID_bool);
        if(lepAntiID_bool) lepAntiID_n++;
      }
      std::vector<int>::iterator begin = lepAntiID_vec.begin();
      std::vector<int>::iterator end = lepAntiID_vec.end();
      if (lepAntiID_n >= 1) {
        antiID_idx0 = std::find(begin, end, true) - begin;
        LepFake0 = *sl->baseLeptons->at(antiID_idx0);
      }
      if (lepAntiID_n >= 2) {
        antiID_idx1 = std::find(begin + antiID_idx0 + 1, end, true) - begin;
        LepFake1 = *sl->baseLeptons->at(antiID_idx1);
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
                  Zlep[0] = lep_ii;
                  Zlep[1] = lep_jj;
                  zlep_idx2 = ii > 0 ? 0 : jj > 1 ? 1 : 2;
              }
          }
      }
      if (sl->leptons->size() == 2 && lepAntiID_n >= 1) { 
        Zlep[2] = sl->baseLeptons->at(antiID_idx0);
      } else if (sl->leptons->size() >= 3) {
        Zlep[2] = sl->leptons->at(zlep_idx2);
      }
      // Find base lepton combination with MLL closest to Z-peak
      // Exclude ID leptons already matched with Z
      Z_diff = FLT_MAX;
      for (uint ii = 0; ii < sl->baseLeptons->size(); ++ii) {
          Susy::Lepton *lep_ii = sl->baseLeptons->at(ii); 
          if (lep_ii == Zlep[0] || lep_ii == Zlep[1]) continue;
          for (uint jj = ii+1; jj < sl->baseLeptons->size(); ++jj) {
              Susy::Lepton *lep_jj = sl->baseLeptons->at(jj); 
              if (lep_jj == Zlep[0] || lep_jj == Zlep[1]) continue;
              float Z_diff_cf = fabs((*lep_ii+*lep_jj).M() - 91.2);
              if (Z_diff_cf < Z_diff) {
                  Z_diff = Z_diff_cf;
                  Zlep[3] = lep_ii;
                  Zlep[4] = lep_jj;
              }
          }
      }

      // All events pass this "cut".
      return true;
    };
  } else {
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
  }
  if (do_fakes_zjets) {
    *cutflow << CutName("2-ID Leptons + 1 Lepton") << [&](Superlink* /*sl*/) -> bool {
        return (lepID_n >= 2 && lepID_n+lepAntiID_n >= 3);
    };
    *cutflow << CutName("SFOS leptons") << [&](Superlink *sl) -> bool {
        (void)sl;
        bool SF = Zlep[0]->isEle() == Zlep[1]->isEle();
        bool OS = Zlep[0]->q != Zlep[1]->q;
        return SF && OS;
    };
  } else if (do_fakes_wjets) {
    *cutflow << CutName("1-ID Lepton + 1 Lepton") << [&](Superlink* /*sl*/) -> bool {
        return (lepID_n >= 1 && lepID_n+lepAntiID_n >= 2);
    };
    //*cutflow << CutName("DFOS leptons") << [&](Superlink* sl) -> bool {
    //    bool DF = sl->leptons->at(0)->isEle() != sl->baseLeptons->at(antiID_idx0)->isEle();
    //    bool OS = sl->leptons->at(0)->q != sl->baseLeptons->at(antiID_idx0)->q;
    //    return DF;
    //};
  }
  *cutflow << CutName("pass good vertex") << [&](Superlink* sl) -> bool {
      return (sl->tools->passGoodVtx(cutflags));
  };
  *cutflow << CutName("jet cleaning") << [&](Superlink* sl) -> bool {
      return (sl->tools->passJetCleaning(sl->baseJets));
  };
  *cutflow << CutName("bad muon veto") << [&](Superlink* sl) -> bool {
      //return (sl->tools->passBadMuon(sl->preMuons));
      // Temporary change to be in agreement with xTau
      return true;
  };
  if (!do_fakes) {
    *cutflow << CutName("2 leptons") << [](Superlink* sl) -> bool {
        return (sl->leptons->size() == 2);
    };
  }
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
  // Define some space saving variables
  LeptonVector preLeptons, baseLeptons, signalLeptons, selectLeptons;
  ElectronVector preElectrons, baseElectrons, signalElectrons;
  MuonVector preMuons, baseMuons, signalMuons;
  TauVector preTaus, baseTaus, signalTaus;
  Susy::Met met;
  Susy::Lepton *el0, *el1, *mu0, *mu1;
  *cutflow << [&](Superlink* sl, var_void*) {
    preLeptons = *sl->preLeptons;
    baseLeptons =  *sl->baseLeptons;
    signalLeptons = *sl->leptons;
    
    if (do_fakes_zjets) selectLeptons = Zlep;
    else selectLeptons = signalLeptons;

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

    el0 = signalElectrons.size() >= 1 ? signalElectrons.at(0) : nullptr;
    el1 = signalElectrons.size() >= 2 ? signalElectrons.at(1) : nullptr;
    mu0 = signalMuons.size() >= 1 ? signalMuons.at(0) : nullptr;
    mu1 = signalMuons.size() >= 2 ? signalMuons.at(1) : nullptr;
  };

  //////////////////////////////////////////////////////////////////////////////
  // Trigger Variables
  // ADD_*_TRIGGER_VAR preprocessor defined
  //////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////////////
  // 2015

  // Dilepton Triggers
  ADD_2LEP_TRIGGER_VAR(HLT_e17_lhloose_mu14, el0, mu0)
  ADD_2LEP_TRIGGER_VAR(HLT_e24_lhmedium_L1EM20VHI_mu8noL1, el0, mu0)
  ADD_2LEP_TRIGGER_VAR(HLT_e7_lhmedium_mu24, mu0, el0)
  ADD_2LEP_TRIGGER_VAR(HLT_2e12_lhloose_L12EM10VH, el0, el1)
  // TODO: Add to SusyNts HLT_2mu10)

  // Single Electron Triggers
  ADD_1LEP_TRIGGER_VAR(HLT_e24_lhmedium_L1EM20VH, el0)
  ADD_1LEP_TRIGGER_VAR(HLT_e60_lhmedium, el0)
  ADD_1LEP_TRIGGER_VAR(HLT_e120_lhloose, el0)

  // Single Muon Triggers
  ADD_1LEP_TRIGGER_VAR(HLT_mu20_iloose_L1MU15, mu0)
  ADD_1LEP_TRIGGER_VAR(HLT_mu40, mu0)

  //////////////////////////////////////////////////////////////////////////////
  // 2016

  // Dilepton Triggers
  ADD_2LEP_TRIGGER_VAR(HLT_e17_lhloose_nod0_mu14, el0, mu0)
  // TODO: Add to SusyNts HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1, el0, mu0)
  ADD_2LEP_TRIGGER_VAR(HLT_e7_lhmedium_nod0_mu24, mu0, el0)
  ADD_2LEP_TRIGGER_VAR(HLT_2e17_lhvloose_nod0, el0, el1)
  // TODO: Add to SusyNts HLT_2mu14)

  // Single Electron Triggers
  ADD_1LEP_TRIGGER_VAR(HLT_e26_lhtight_nod0_ivarloose, el0)
  ADD_1LEP_TRIGGER_VAR(HLT_e60_lhmedium_nod0, el0)
  ADD_1LEP_TRIGGER_VAR(HLT_e140_lhloose_nod0, el0)

  // Single Muon Triggers
  ADD_1LEP_TRIGGER_VAR(HLT_mu26_ivarmedium, mu0)
  ADD_1LEP_TRIGGER_VAR(HLT_mu50, mu0)


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
      for (auto& el : preElectrons) {out.push_back(el->Eta());}
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
  *cutflow << NewVar("Electron ID (non-inclusive)"); {
    *cutflow << HFTname("El_ID");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for (auto& el : signalElectrons) {
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
      for (auto& el : signalElectrons) {
        for (auto& mu : preMuons) {
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
      for (auto& el : signalElectrons) {
        for (auto& mu : preMuons) {
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
      for (auto& el : signalElectrons) {
        for (auto& mu : preMuons) {
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
      for (auto& el : signalElectrons) {
        for (auto& mu : preMuons) {
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
      for (auto& el : signalElectrons) {out.push_back(el->d0sigBSCorr); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Electron z0SinTheta"); {
    *cutflow << HFTname("el_z0SinTheta");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for (auto& el : signalElectrons) {out.push_back(el->z0SinTheta()); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("subleading electron track pt"); {
    *cutflow << HFTname("el1_track_pt");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
      if (signalLeptons.size() > 1 && signalLeptons.at(1)->isEle()) {
          const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(signalLeptons.at(1));
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
      if (signalLeptons.size() > 1 && signalLeptons.at(1)->isEle()) {
          const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(signalLeptons.at(1));
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
      if (signalLeptons.size() > 1 && signalLeptons.at(1)->isEle()) {
          const Susy::Electron* ele = dynamic_cast<const Susy::Electron*>(signalLeptons.at(1));
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
      for (auto& mu : preMuons) {out.push_back(mu->Pt()); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("isCaloTagged"); {
    *cutflow << HFTname("isCaloTagged");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for (auto& mu : preMuons) { out.push_back(mu->isCaloTagged);}
      return out;
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("PreMuon ID (non-inclusive)"); {
    *cutflow << HFTname("preMu_ID");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for (auto& mu : preMuons) {
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
    *cutflow << HFTname("baseMu_ID");
    *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
      vector<int> out;
      for (auto& mu : baseMuons) {
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
      for (auto& mu : signalMuons) {
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
  *cutflow << NewVar("Muon d0sigBSCorr"); {
    *cutflow << HFTname("mu_d0sigBSCorr");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for (auto& mu : signalMuons) {out.push_back(mu->d0sigBSCorr); }
      return out;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Muon z0SinTheta"); {
    *cutflow << HFTname("mu_z0SinTheta");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for (auto& mu : signalMuons) {out.push_back(mu->z0SinTheta()); }
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
      for (auto& lep : signalLeptons) {
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
      for(auto& lepton : selectLeptons) {
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
      for(auto& lepton : selectLeptons) {
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
      for(auto& lepton : selectLeptons) {
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
      for(auto& lepton : selectLeptons) {
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
      for(auto& lepton : selectLeptons) {
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
      for(auto& lepton : selectLeptons) {
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
      for(auto& lepton : selectLeptons) {
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
      for(auto& lepton : selectLeptons) {
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
      for(auto& lepton : selectLeptons) {
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
      for(auto& lepton : selectLeptons) {
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
      for(auto& lepton : selectLeptons) {
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
        bool mother_is_piZero = fabs(M_ID) == 111;
        bool bkgEl_from_phoConv = T==BkgElectron && O==PhotonConv;
        //bool noChargeFlip = M_ID*lepton->q < 0;
        //bool chargeFlip = M_ID*lepton->q > 0;
        
        bool promptEl1 = T==IsoElectron; //&& noChargeFlip;
        bool promptEl2 = (bkgEl_from_phoConv && mother_is_el); //&& noChargeFlip;
        bool promptEl3 = bkgEl_from_phoConv && MO==FSRPhot;
        bool promptEl4 = T==NonIsoPhoton && O==FSRPhot;
        bool promptEl = promptEl1 || promptEl2 || promptEl3 || promptEl4;

        //bool promptChargeFlipEl1 = T==IsoElectron && chargeFlip;
        //bool promptChargeFlipEl2 = (bkgEl_from_phoConv && mother_is_el) && chargeFlip;
        //bool promptChargeFlipEl = promptChargeFlipEl1 || promptChargeFlipEl2;
        
        bool promptMuon = T==IsoMuon && (
            O==top || O==WBoson || O==ZBoson || O==Higgs || O==HiggsMSSM || 
            O==MCTruthPartClassifier::SUSY || O==DiBoson);

        bool promptPho1 = T==IsoPhoton && O==PromptPhot;
        bool promptPho2 = bkgEl_from_phoConv && MT==IsoPhoton && MO==PromptPhot;
        bool promptPho3 = bkgEl_from_phoConv && MT==BkgPhoton && MO==UndrPhot;
        bool promptPho = promptPho1 || promptPho2 || promptPho3;

        bool hadDecay1 = T==BkgElectron && (
            O==DalitzDec || O==ElMagProc || O==LightMeson || O==StrangeMeson);
        bool hadDecay2 = bkgEl_from_phoConv && MT==BkgPhoton && (
            MO==PiZero || MO==LightMeson || MO==StrangeMeson);
        bool hadDecay3 = T==BkgPhoton && (O==LightMeson || O==PiZero);
        bool hadDecay4 = T==BkgMuon && (
            O==LightMeson || O==StrangeMeson || O==PionDecay || O==KaonDecay);
        bool hadDecay5 = T==Hadron;
        bool hadDecay = hadDecay1 || hadDecay2 || hadDecay3 || hadDecay4 || hadDecay5;

        bool HF_tau_mu1 =  (T==NonIsoElectron || T==NonIsoPhoton) && O==TauLep;
        bool HF_tau_mu2 =  bkgEl_from_phoConv && MT==NonIsoPhoton && MO==TauLep;
        bool HF_tau_mu3 =  T==NonIsoMuon && O==TauLep;
        bool HF_tau_mu4 =  (T==NonIsoElectron || T==NonIsoPhoton) && O==Mu;
        bool HF_tau_mu5 =  bkgEl_from_phoConv && MT==NonIsoPhoton && MO==Mu;
        bool HF_tau_mu =  HF_tau_mu1 || HF_tau_mu2 || HF_tau_mu3 || HF_tau_mu4 || HF_tau_mu5;

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
        else if (bkgEl_from_phoConv && mother_is_piZero) lep_class = 5;
        else if (HF_tau_mu) lep_class = 6;
        else if (HF_B) lep_class = 7;
        else if (HF_C) lep_class = 8;
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
  //xTauFW variable
  *cutflow << NewVar("leptons sign product"); {
    *cutflow << HFTname("LepLepSign");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
      if (signalLeptons.size() < 2) return -INT_MAX;
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
    *cutflow << [&](Superlink* sl, var_double*) -> double {
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
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
      return MET.Phi();
    };
    *cutflow << SaveVar();
  }

  //////////////////////////////////////////////////////////////////////////////
  // Two-lepton properties

  // Gather useful information
  TLorentzVector lepton0 ;
  TLorentzVector lepton1 ;
  TLorentzVector dileptonP4 ;
  *cutflow << [&](Superlink* /*sl*/, var_void*) {
    if (signalLeptons.size() >= 1) lepton0 = *signalLeptons.at(0);
    if (signalLeptons.size() >= 2) lepton1 = *signalLeptons.at(1);
    dileptonP4 = lepton0 + lepton1;
  };

  // Lepton variables
  //xTauFW variable
  *cutflow << NewVar("Leading lepton pT"); {
    *cutflow << HFTname("Lep0Pt");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
      if (signalLeptons.size() < 1) return -DBL_MAX;
      return lepton0.Pt(); };
    *cutflow << SaveVar();
  }
  //xTauFW variable
  *cutflow << NewVar("Subleading lepton pT"); {
    *cutflow << HFTname("Lep1Pt");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
      if (signalLeptons.size() < 2) return -DBL_MAX;
      return lepton1.Pt(); };
    *cutflow << SaveVar();
  }
  //xTauFW variable
  *cutflow << NewVar("Leading lepton eta"); {
    *cutflow << HFTname("Lep0Eta");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
      if (signalLeptons.size() < 1) return -DBL_MAX;
      return lepton0.Eta(); };
    *cutflow << SaveVar();
  }
  //xTauFW variable
  *cutflow << NewVar("Subleading lepton eta"); {
    *cutflow << HFTname("Lep1Eta");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
      if (signalLeptons.size() < 2) return -DBL_MAX;
      return lepton1.Eta(); };
    *cutflow << SaveVar();
  }
  //xTauFW variable
  *cutflow << NewVar("Leading lepton phi"); {
    *cutflow << HFTname("Lep0Phi");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
      if (signalLeptons.size() < 1) return -DBL_MAX;
      return lepton0.Phi(); };
    *cutflow << SaveVar();
  }
  //xTauFW variable
  *cutflow << NewVar("Subleading lepton phi"); {
    *cutflow << HFTname("Lep1Phi");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
      if (signalLeptons.size() < 2) return -DBL_MAX;
      return lepton1.Phi(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("Leading lepton invariant mass"); {
    *cutflow << HFTname("MLep0");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
      if (signalLeptons.size() < 1) return -DBL_MAX;
      return lepton0.M(); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Subleading lepton invariant mass"); {
    *cutflow << HFTname("MLep1");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
      if (signalLeptons.size() < 2) return -DBL_MAX;
      return lepton1.M(); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Delta eta between leptons"); {
    *cutflow << HFTname("DEtaLL");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
      if (signalLeptons.size() < 2) return -DBL_MAX;
      return fabs(lepton0.Eta() - lepton1.Eta()); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Delta phi between leptons"); {
    *cutflow << HFTname("DphiLL");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
      if (signalLeptons.size() < 2) return -DBL_MAX;
      return fabs(lepton0.DeltaPhi(lepton1));};
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("delta R of di-lepton system"); {
    *cutflow << HFTname("drll");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
      if (signalLeptons.size() < 2) return -DBL_MAX;
      return lepton0.DeltaR(lepton1); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("dilepton flavor"); {
    *cutflow << HFTname("dilep_flav");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
      if(signalLeptons.size()<2) return -INT_MAX;
      if(signalLeptons.at(0)->isEle() && signalLeptons.at(1)->isMu()){return 0;}       // e mu  case
      else if(signalLeptons.at(0)->isMu() && signalLeptons.at(1)->isEle()){return 1;}  // mu e  case
      else if(signalLeptons.at(0)->isEle() && signalLeptons.at(1)->isEle()){return 2;} // e e   case
      else if(signalLeptons.at(0)->isMu() && signalLeptons.at(1)->isMu()){return 3;}   // mu mu case
      return 4;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("dilepton flavor is e mu"); {
    *cutflow << HFTname("isEM");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
	    if(signalLeptons.size()<2) return -INT_MAX;
      return signalLeptons.at(0)->isEle() && signalLeptons.at(1)->isMu();  // e mu  case
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("dilepton flavor is mu e"); {
    *cutflow << HFTname("isME");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if(signalLeptons.size()<2) return -INT_MAX;
        return signalLeptons.at(0)->isMu() && signalLeptons.at(1)->isEle();  // mu e  case
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("collinear mass, M_coll"); {
    *cutflow << HFTname("MCollASym");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
      if (signalLeptons.size() < 2) return -DBL_MAX;
      double deta = fabs(lepton0.Eta()-lepton1.Eta());
	    double dphi = lepton0.DeltaPhi(lepton1);
	    return sqrt(2*lepton0.Pt()*(lepton1.Pt()+met.Et)*(cosh(deta)-cos(dphi)));
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("mass of di-lepton system, M_ll"); {
    *cutflow << HFTname("MLL");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
      if (signalLeptons.size() < 2) return -DBL_MAX;
      return dileptonP4.M(); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Pt of di-lepton system, Pt_ll"); {
    *cutflow << HFTname("ptll");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
      if (signalLeptons.size() < 2) return -DBL_MAX;
      return dileptonP4.Pt(); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("diff_pt between leading and sub-leading lepton"); {
    *cutflow << HFTname("dpt_ll");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (signalLeptons.size() < 2) return -DBL_MAX;
        return lepton0.Pt() - lepton1.Pt();
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("dphi between leading lepton and MET"); {
    *cutflow << HFTname("DphiLep0MET");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (signalLeptons.size() < 1) return -DBL_MAX;
        return lepton0.DeltaPhi(MET);
    };
    *cutflow << SaveVar();
  }

  //////////////////////////////////////////////////////////////////////////////
  // Optimized Region Cuts
  *cutflow << NewVar("dphi between sub-leading lepton and MET"); {
    *cutflow << HFTname("DphiLep1MET");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (signalLeptons.size() < 2) return -DBL_MAX;
        return lepton1.DeltaPhi(MET);
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("transverse mass of leading lepton"); {
    *cutflow << HFTname("MtLep0");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (signalLeptons.size() < 1) return -DBL_MAX;
        return sqrt( 2*lepton0.Pt()*MET.Pt()*(1-cos (lepton0.DeltaPhi(MET)) ) );
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("transverse mass of subleading lepton"); {
    *cutflow << HFTname("MtLep1");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (signalLeptons.size() < 2) return -DBL_MAX;
        return sqrt( 2*lepton1.Pt()*MET.Pt()*(1-cos (lepton1.DeltaPhi(MET)) ) );
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("approximate tau pT"); {
    *cutflow << HFTname("tau_pT");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (signalLeptons.size() < 2) return -DBL_MAX;
        return (MET + lepton1).Pt();
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("ratio of tau pT to subleading lepton pT"); {
    *cutflow << HFTname("taulep1_pT_ratio");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (signalLeptons.size() < 2) return -DBL_MAX;
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
  JetVector lightJets, BJets, forwardJets;
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
  *cutflow << NewVar("number of light jets"); {
    *cutflow << HFTname("nLJets");
    *cutflow << [&](Superlink* sl, var_int*) -> int { return sl->tools->numberOfLJets(signalJets)/*(*baseJets)*/; };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("number of b jets"); {
    *cutflow << HFTname("nBJets");
    *cutflow << [&](Superlink* sl, var_int*) -> int { return sl->tools->numberOfBJets(signalJets)/*(*baseJets)*/; };
    *cutflow << SaveVar();
  }
  //xTauFW variable
  *cutflow << NewVar("b-tagged jet"); {
    *cutflow << HFTname("Btag");
    *cutflow << [&](Superlink* sl, var_bool*) -> bool { return sl->tools->numberOfBJets(signalJets) > 0;};
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("number of forward jets"); {
    *cutflow << HFTname("nForwardJets");
    *cutflow << [&](Superlink* sl, var_int*) -> int { return sl->tools->numberOfFJets(signalJets)/*(*baseJets)*/; };
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
  *cutflow << NewVar("jet flavor (0: NA, 1: L, 2: B, 3: F)"); {
    *cutflow << HFTname("j_flav");
    *cutflow << [&](Superlink* sl, var_int_array*) -> vector<int> {
      vector<int> out; int flav = 0;
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
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (baseJets.size() < 2) return -DBL_MAX;
        return JetP4.M();
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("DeltaEta between two leading jets"); {
    *cutflow << HFTname("DEtaJJ");
    *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (baseJets.size() < 2) return -DBL_MAX;
        return fabs(Jet0.Eta() - Jet1.Eta());
    };
    *cutflow << SaveVar();
  }

  //////////////////////////////////////////////////////////////////////////////
  // Fake variables
  //////////////////////////////////////////////////////////////////////////////
  if (do_fakes) {
    *cutflow << NewVar("Lepton ID status (ID(1), antiID(-1))"); {
      *cutflow << HFTname("lepID_status");
      *cutflow << [&](Superlink* /*sl*/, var_int_array*) -> vector<int> {
          vector<int> out;
          for (uint i = 0; i < baseLeptons.size(); i++) {
              // ID -> 1, AntiID -> -1
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
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { 
          return lepID_n; 
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Number of Anti-ID leptons"); {
      *cutflow << HFTname("nLepAntiID");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { 
          return lepAntiID_n; 
      };
      *cutflow << SaveVar();
    }
    //////////////////////////////////////////////////////////////////////////////
    // 1 antiID lepton case
    *cutflow << NewVar("Leading anti-ID lepton charge"); {
      *cutflow << HFTname("aID_Lep0Q");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if (lepAntiID_n < 1) return -INT_MAX; 
        auto* fakelep0 = baseLeptons.at(antiID_idx0);
        return fakelep0->q;
        };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Leading aID-lepton pT"); {
      *cutflow << HFTname("aID_Lep0Pt");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (lepAntiID_n < 1) return -DBL_MAX;
        return LepFake0.Pt(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Leading aID-lepton eta"); {
      *cutflow << HFTname("aID_Lep0Eta");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (lepAntiID_n < 1) return -DBL_MAX;
        return LepFake0.Eta(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Leading aID-lepton invariant mass"); {
      *cutflow << HFTname("aID_MLep0");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (lepAntiID_n < 1) return -DBL_MAX;
        return LepFake0.M(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("transverse mass of leading aID-lepton"); {
      *cutflow << HFTname("aID_MtLep0");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
          if (lepAntiID_n < 1) return -DBL_MAX;
          return sqrt( 2*LepFake0.Pt()*MET.Pt()*(1-cos (LepFake0.DeltaPhi(MET)) ) );
      };
      *cutflow << SaveVar();
    }
    // Z + jets variables
    *cutflow << NewVar("Z leptons' flavor"); {
      *cutflow << HFTname("Z_dilep_flav");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if(signalLeptons.size() < 2) return -INT_MAX;
        if(Zlep[0]->isEle() && Zlep[1]->isMu()){return 0;}       // e mu  case
        else if(Zlep[0]->isMu() && Zlep[1]->isEle()){return 1;}  // mu e  case
        else if(Zlep[0]->isEle() && Zlep[1]->isEle()){return 2;} // e e   case
        else if(Zlep[0]->isMu() && Zlep[1]->isMu()){return 3;}   // mu mu case
        return 4;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Z dilepton sign"); {
      *cutflow << HFTname("Z_dilep_sign");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if(signalLeptons.size() < 2) return -INT_MAX;
        return Zlep[0]->q * Zlep[1]->q;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("MLL of closest Z pair"); {
      *cutflow << HFTname("Z_MLL");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (signalLeptons.size() < 2) return -DBL_MAX;
        return (*Zlep[0]+*Zlep[1]).M(); 
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("2nd Z leptons' flavor"); {
      *cutflow << HFTname("Z2_dilep_flav");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if(signalLeptons.size() < 2) return -INT_MAX;
        else if (!Zlep[3] or !Zlep[4]) return -2;
        if(Zlep[3]->isEle() && Zlep[4]->isMu()){return 0;}       // e mu  case
        else if(Zlep[3]->isMu() && Zlep[4]->isEle()){return 1;}  // mu e  case
        else if(Zlep[3]->isEle() && Zlep[4]->isEle()){return 2;} // e e   case
        else if(Zlep[3]->isMu() && Zlep[4]->isMu()){return 3;}   // mu mu case
        return 4;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("2nd Z dilepton sign"); {
      *cutflow << HFTname("Z2_dilep_sign");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if(signalLeptons.size() < 2) return -INT_MAX;
        else if (!Zlep[3] or !Zlep[4]) return -2;
        return Zlep[3]->q * Zlep[4]->q;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("MLL of 2nd-closest Z pair"); {
      *cutflow << HFTname("Z2_MLL");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (signalLeptons.size() < 2) return -DBL_MAX;
        else if (!Zlep[3] or !Zlep[4]) return -5.0;
        return (*Zlep[3]+*Zlep[4]).M(); 
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Non-Z lepton pT"); {
      *cutflow << HFTname("Z_Lep2_pT");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (Zlep[2]) return Zlep[2]->Pt();
        return -DBL_MAX;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Non-Z lepton eta"); {
      *cutflow << HFTname("Z_Lep2_eta");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (Zlep[2]) return Zlep[2]->Eta();
        return -DBL_MAX;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Non-Z lepton flavor"); {
      *cutflow << HFTname("Z_Lep2_flav");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if (Zlep[2]) return Zlep[2]->isEle();
        return -INT_MAX;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Non-Z lepton charge"); {
      *cutflow << HFTname("Z_Lep2_q");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if (Zlep[2]) return Zlep[2]->q;
        return -INT_MAX;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("DeltaPhi of Non-Z lep and MET"); {
      *cutflow << HFTname("Z_Lep2_dPhi_MET");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (Zlep[2]) return Zlep[2]->DeltaPhi(MET);
        return -DBL_MAX;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Non-Z lepton mT"); {
      *cutflow << HFTname("Z_Lep2_mT");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (Zlep[2]) {
            return sqrt( 2*Zlep[2]->Pt()*MET.Pt()*(1-cos (Zlep[2]->DeltaPhi(MET)) ) );
        }
        return -DBL_MAX;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta R of Z and aID-Lepton"); {
      *cutflow << HFTname("dR_Zl");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (Zlep[2]) return Zlep[2]->DeltaR(*Zlep[0]+*Zlep[1]);
        return -DBL_MAX;
      }; 
      *cutflow << SaveVar();
    }
    // W + jets variables
    *cutflow << NewVar("aID-ID dilepton flavor"); {
      *cutflow << HFTname("aID_dilep_flav");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if (lepAntiID_n < 1 || signalLeptons.size() < 1) return -INT_MAX; 
        auto* fakelep0 = baseLeptons.at(antiID_idx0);
        auto* lep0 = signalLeptons.at(0);
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
        if (lepAntiID_n < 1 || signalLeptons.size() < 1) return -INT_MAX; 
        auto* fakelep0 = baseLeptons.at(antiID_idx0);
        auto* lep0 = signalLeptons.at(0);
        return fakelep0->q * lep0->q; 
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta R of ID-lepton and aID-lepton"); {
      *cutflow << HFTname("aID_drll");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (lepAntiID_n < 1 || signalLeptons.size() < 1) return -DBL_MAX;
        return LepFake0.DeltaR(lepton0); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("ID-aID collinear mass, M_coll"); {
      *cutflow << HFTname("aID_MCollASym");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (lepAntiID_n < 1 || signalLeptons.size() < 1) return -DBL_MAX;
        double deta = fabs(LepFake0.Eta()-lepton0.Eta());
        double dphi = LepFake0.DeltaPhi(lepton0);
        return sqrt(2*LepFake0.Pt()*(lepton0.Pt()+met.Et)*(cosh(deta)-cos(dphi)));
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("mass of ID-aID-dilepton system, M_ll"); {
      *cutflow << HFTname("aID_MLL");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (lepAntiID_n < 1 || signalLeptons.size() < 1) return -DBL_MAX;
        return (LepFake0+lepton0).M(); 
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Pt of ID-aID-dilepton system, Pt_ll"); {
      *cutflow << HFTname("aID_ptll");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (lepAntiID_n < 1 || signalLeptons.size() < 1) return -DBL_MAX;
        return (LepFake0+lepton0).Pt(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("diff_pt between ID-aID-leptons"); {
      *cutflow << HFTname("aID_dpt_ll");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
          if (lepAntiID_n < 1 || signalLeptons.size() < 1) return -DBL_MAX;
          return LepFake0.Pt() - lepton0.Pt();
      };
      *cutflow << SaveVar();
    }

    //////////////////////////////////////////////////////////////////////////////
    // 2 antiID lepton case (QCD)
    *cutflow << NewVar("Subleading anti-ID lepton flavor"); {
      *cutflow << HFTname("aID_Lep1Flav");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { 
          if (lepAntiID_n < 2) return -INT_MAX;
          auto* fakelep1 = baseLeptons.at(antiID_idx1);
          return  fakelep1->isEle() ? 0 : 1;
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Subleading anti-ID lepton charge"); {
      *cutflow << HFTname("aID_Lep1Q");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if (lepAntiID_n < 2) return -INT_MAX; 
        auto* fakelep1 = baseLeptons.at(antiID_idx1);
        return fakelep1->q;
        };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Subleading aID-lepton pT"); {
      *cutflow << HFTname("aID2_Lep1Pt");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (lepAntiID_n < 2) return -DBL_MAX;
        return LepFake1.Pt(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Subleading aID-lepton eta"); {
      *cutflow << HFTname("aID2_Lep1Eta");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (lepAntiID_n < 2) return -DBL_MAX;
        return LepFake1.Eta(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Subleading aID-lepton invariant mass"); {
      *cutflow << HFTname("aID2_MLep1");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (lepAntiID_n < 2) return -DBL_MAX;
        return LepFake1.M(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("delta R of aID-leptons"); {
      *cutflow << HFTname("aID2_drll");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (lepAntiID_n < 2) return -DBL_MAX;
        return LepFake0.DeltaR(LepFake1); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("aID-dilepton flavor"); {
      *cutflow << HFTname("aID2_dilep_flav");
      *cutflow << [&](Superlink* /*sl*/, var_int*) -> int {
        if(lepAntiID_n < 2) return -INT_MAX;
        auto* fakelep0 = baseLeptons.at(antiID_idx0);
        auto* fakelep1 = baseLeptons.at(antiID_idx1);
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
        if (lepAntiID_n < 2) return -DBL_MAX;
        double deta = fabs(LepFake0.Eta()-LepFake1.Eta());
        double dphi = LepFake0.DeltaPhi(LepFake1);
        return sqrt(2*LepFake0.Pt()*(LepFake1.Pt()+met.Et)*(cosh(deta)-cos(dphi)));
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("mass of aID-dilepton system, M_ll"); {
      *cutflow << HFTname("aID2_MLL");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (lepAntiID_n < 2) return -DBL_MAX;
        return (LepFake0+LepFake1).M(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("Pt of aID-dilepton system, Pt_ll"); {
      *cutflow << HFTname("aID2_ptll");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
        if (lepAntiID_n < 2) return -DBL_MAX;
        return (LepFake0+LepFake1).Pt(); };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("diff_pt between aID-leptons"); {
      *cutflow << HFTname("aID2_dpt_ll");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
          if (lepAntiID_n < 2) return -DBL_MAX;
          return LepFake0.Pt() - LepFake1.Pt();
      };
      *cutflow << SaveVar();
    }
    *cutflow << NewVar("transverse mass of subleading aID-lepton"); {
      *cutflow << HFTname("aID2_MtLep1");
      *cutflow << [&](Superlink* /*sl*/, var_double*) -> double {
          if (lepAntiID_n < 2) return -DBL_MAX;
          return sqrt( 2*LepFake1.Pt()*MET.Pt()*(1-cos (LepFake1.DeltaPhi(MET)) ) );
      };
      *cutflow << SaveVar();
    }
  }
  ///////////////////////////////////////////////////////////////////////
  // Clear internal variables
  ///////////////////////////////////////////////////////////////////////
  // NOTE!: Superflow assumes this expression is the last to be added. 
  *cutflow << [&](Superlink* /*sl*/, var_void*) {
    preLeptons.clear(); baseLeptons.clear(); signalLeptons.clear(), selectLeptons.clear();
    preElectrons.clear(); baseElectrons.clear(); signalElectrons.clear();
    preMuons.clear(); baseMuons.clear(); signalMuons.clear();
    preTaus.clear(); baseTaus.clear(); signalTaus.clear();
    preJets.clear(); baseJets.clear(); signalJets.clear();
    lightJets.clear(); BJets.clear(); forwardJets.clear();
    if (do_fakes) {
        lepID_vec.clear(); lepAntiID_vec.clear();
        lepID_n = lepAntiID_n = 0;
        antiID_idx0 = antiID_idx1 = -1;
        LepFake0 = LepFake1 = {};
        std::fill(Zlep.begin(), Zlep.end(), nullptr);
    }
  };

  ///////////////////////////////////////////////////////////////////////
  // END OF CUTFLOW EXPRESSIONS 
  ///////////////////////////////////////////////////////////////////////

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
