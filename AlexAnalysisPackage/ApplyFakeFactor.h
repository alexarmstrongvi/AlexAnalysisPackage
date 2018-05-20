/**
 * @file LeptonFakeFactorTools/ApplyFakeFactor.h
 * @authors Joseph Reichert <Joseph.Reichert@cern.ch>
 * @date August 2016
 * @brief Tool to read fake factor histograms and return 
 * the binned fake factor value
**/

#ifndef ApplyFakeFactor_H
#define ApplyFakeFactor_H

#include <AsgTools/MessageCheck.h>
#include <AsgTools/AsgMessaging.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <fstream>
#include <iostream>
#include <string> // for string
#include <vector> // for vector
#include <cmath>

namespace LepEnum{
  enum LepType{
    Electron,
    Muon,
    Tau,
    NonLepton
  };
}

class ApplyFakeFactor : public asg::AsgMessaging
{

public: 
  /// Standard constructor
  ApplyFakeFactor(const std::string& name = "ApplyFakeFactor");
  
  /// Standard destructor
  ~ApplyFakeFactor();

  // Main methods
public:
  // Initialize this class
  StatusCode initialize();

  // Finalize this class; everything that should be done after the 
  // event loop should go here
  StatusCode finalize();

  // The main apply method -- assumes pt in GeV
  // FIXME May want to add methods for pt in GeV vs. MeV
  // May also want an xAOD method
  double apply(double pt, LepEnum::LepType typeOfLep);
  double apply(double pt, double eta, LepEnum::LepType typeOfLep);

  double getStatErr(double pt, LepEnum::LepType typeOfLep);
  double getStatErr(double pt, double eta, LepEnum::LepType typeOfLep);

  double getSystErr(double pt, LepEnum::LepType typeOfLep);
  double getSystErr(double pt, double eta, LepEnum::LepType typeOfLep);

  void set_savedFakeFactorFileName(TString s){m_savedFakeFactorFileName = s;};

  void set_saved_elFakeFactorName(TString s){m_saved_elFakeFactorName = s;};
  void set_saved_muFakeFactorName(TString s){m_saved_muFakeFactorName = s;};
  void set_saved_tauFakeFactorName(TString s){m_saved_tauFakeFactorName = s;};

  void set_saved_elFakeFactorName_2D(TString s){m_saved_elFakeFactorName_2D = s;};
  void set_saved_muFakeFactorName_2D(TString s){m_saved_muFakeFactorName_2D = s;};
  void set_saved_tauFakeFactorName_2D(TString s){m_saved_tauFakeFactorName_2D = s;};

private:

  int getBin(TH1F* hist, double pt);
  int getBinX(TH2F* hist, double pt);
  int getBinY(TH2F* hist, double eta);

  TH1F* get1DFakeFactorHist(LepEnum::LepType typeOfLep, bool getSystErrHist=false);
  TH2F* get2DFakeFactorHist(LepEnum::LepType typeOfLep, bool getSystErrHist=false);
  
  TH1F* Load1DHistogram(TFile* fakeFactorFile, TString name);
  TH2F* Load2DHistogram(TFile* fakeFactorFile, TString name);

  TString m_savedFakeFactorFileName = "";

  TString m_saved_elFakeFactorName = "";
  TString m_saved_muFakeFactorName = "";
  TString m_saved_tauFakeFactorName = "";

  TString m_saved_elFakeFactorName_2D = "";
  TString m_saved_muFakeFactorName_2D = "";
  TString m_saved_tauFakeFactorName_2D = "";

  TH1F* m_saved_elFakeFactor = 0;
  TH1F* m_saved_muFakeFactor = 0;
  TH1F* m_saved_tauFakeFactor = 0;

  TH2F* m_saved_elFakeFactor_2D = 0;
  TH2F* m_saved_muFakeFactor_2D = 0;
  TH2F* m_saved_tauFakeFactor_2D = 0;

  TH1F* m_saved_elFakeFactorSyst = 0;
  TH1F* m_saved_muFakeFactorSyst = 0;
  TH1F* m_saved_tauFakeFactorSyst = 0;

  TH2F* m_saved_elFakeFactorSyst_2D = 0;
  TH2F* m_saved_muFakeFactorSyst_2D = 0;
  TH2F* m_saved_tauFakeFactorSyst_2D = 0;
};

//----------------------------------------------------------------------------------------
#endif
