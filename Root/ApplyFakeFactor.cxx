#include "AlexAnalysisPackage/ApplyFakeFactor.h"

//=============================================================================
// Constructor
//=============================================================================
ApplyFakeFactor::ApplyFakeFactor(const std::string& name) :
  asg::AsgMessaging(name)
{
  return;
}

//=============================================================================
// Destructor
//=============================================================================
ApplyFakeFactor::~ApplyFakeFactor()
{
  return;
}


StatusCode ApplyFakeFactor::initialize()
{
  StatusCode sc = StatusCode::SUCCESS;

  // Open saved fake factor file, and read in the histograms
  TFile* fakeFactorFile = TFile::Open(m_savedFakeFactorFileName, "READ");

  if(!fakeFactorFile){
    ATH_MSG_FATAL("Invalid saved fake factor file! Exiting. Expected name: " << m_savedFakeFactorFileName);
    abort();
  }

  m_saved_elFakeFactor  = Load1DHistogram(fakeFactorFile, m_saved_elFakeFactorName);
  m_saved_muFakeFactor  = Load1DHistogram(fakeFactorFile, m_saved_muFakeFactorName);
  m_saved_tauFakeFactor = Load1DHistogram(fakeFactorFile, m_saved_tauFakeFactorName);

  m_saved_elFakeFactor_2D  = Load2DHistogram(fakeFactorFile, m_saved_elFakeFactorName_2D);
  m_saved_muFakeFactor_2D  = Load2DHistogram(fakeFactorFile, m_saved_muFakeFactorName_2D);
  m_saved_tauFakeFactor_2D = Load2DHistogram(fakeFactorFile, m_saved_tauFakeFactorName_2D);

  m_saved_elFakeFactorSyst  = Load1DHistogram(fakeFactorFile, m_saved_elFakeFactorName+"_Syst");
  m_saved_muFakeFactorSyst  = Load1DHistogram(fakeFactorFile, m_saved_muFakeFactorName+"_Syst");
  m_saved_tauFakeFactorSyst = Load1DHistogram(fakeFactorFile, m_saved_tauFakeFactorName+"_Syst");

  m_saved_elFakeFactorSyst_2D  = Load2DHistogram(fakeFactorFile, m_saved_elFakeFactorName_2D+"_Syst");
  m_saved_muFakeFactorSyst_2D  = Load2DHistogram(fakeFactorFile, m_saved_muFakeFactorName_2D+"_Syst");
  m_saved_tauFakeFactorSyst_2D = Load2DHistogram(fakeFactorFile, m_saved_tauFakeFactorName_2D+"_Syst");

  fakeFactorFile->Close();

  return sc;
}

StatusCode ApplyFakeFactor::finalize()
{
  StatusCode sc = StatusCode::SUCCESS;

  std::vector<TH1F*> hists = {
    m_saved_elFakeFactor,
    m_saved_muFakeFactor,
    m_saved_tauFakeFactor
  };

  std::vector<TH2F*> hists2D = {
    m_saved_elFakeFactor_2D,
    m_saved_muFakeFactor_2D,
    m_saved_tauFakeFactor_2D
  };

  for(auto hist : hists){
    if(hist) delete hist;
  }

  for(auto hist : hists2D){
    if(hist) delete hist;
  }

  return sc;
}

double ApplyFakeFactor::apply(double pt, LepEnum::LepType typeOfLep){

  TH1F* hist = get1DFakeFactorHist(typeOfLep);
  if(!hist) return 0;
  int bin = getBin(hist,pt);
  return hist->GetBinContent(bin);
}

double ApplyFakeFactor::apply(double pt, double eta, LepEnum::LepType typeOfLep){

  TH2F* hist = get2DFakeFactorHist(typeOfLep);
  if(!hist) return 0;
  int xbin = getBinX(hist, pt);
  int ybin = getBinY(hist, eta);
  return hist->GetBinContent(xbin, ybin);
}

double ApplyFakeFactor::getStatErr(double pt, LepEnum::LepType typeOfLep){

  TH1F* hist = get1DFakeFactorHist(typeOfLep);
  if(!hist) return 0;
  int bin = getBin(hist,pt);
  return hist->GetBinError(bin);
}

double ApplyFakeFactor::getStatErr(double pt, double eta, LepEnum::LepType typeOfLep){

  TH2F* hist = get2DFakeFactorHist(typeOfLep);
  if(!hist) return 0;
  int xbin = getBinX(hist, pt);
  int ybin = getBinY(hist, eta);
  return hist->GetBinError(xbin, ybin);
}

double ApplyFakeFactor::getSystErr(double pt, LepEnum::LepType typeOfLep){

  TH1F* hist = get1DFakeFactorHist(typeOfLep, true);
  if(!hist) return 0;
  int bin = getBin(hist,pt);
  return hist->GetBinContent(bin);
}

double ApplyFakeFactor::getSystErr(double pt, double eta, LepEnum::LepType typeOfLep){

  TH2F* hist = get2DFakeFactorHist(typeOfLep, true);
  if(!hist) return 0;
  int xbin = getBinX(hist, pt);
  int ybin = getBinY(hist, eta);
  return hist->GetBinContent(xbin, ybin);
}

int ApplyFakeFactor::getBin(TH1F* hist, double pt){

  // don't read values from over/underflow bins
  int bin = hist->FindBin(pt);
  if(bin == 0) bin += 1;
  if(bin == hist->GetNbinsX()+1) bin -= 1;

  return bin;
}
int ApplyFakeFactor::getBinX(TH2F* hist, double pt){

  // don't read values from over/underflow bins
  int xbin = hist->GetXaxis()->FindBin(pt);
  if(xbin == 0) xbin += 1;
  if(xbin == hist->GetNbinsX()+1) xbin -= 1;

  return xbin;
}
int ApplyFakeFactor::getBinY(TH2F* hist, double eta){

  eta = fabs(eta); // seems safe to assume things are binned in abs(eta)

  // don't read values from over/underflow bins
  int ybin = hist->GetYaxis()->FindBin(eta);
  if(ybin == 0) ybin += 1;
  if(ybin == hist->GetNbinsY()+1) ybin -= 1;

  return ybin;
}

TH1F* ApplyFakeFactor::get1DFakeFactorHist(LepEnum::LepType typeOfLep, bool getSystErrHist){

  TH1F* hist = 0;

  if(getSystErrHist){
    if(typeOfLep == LepEnum::Electron)  hist = m_saved_elFakeFactorSyst;
    else if(typeOfLep == LepEnum::Muon) hist = m_saved_muFakeFactorSyst;
    else if(typeOfLep == LepEnum::Tau)  hist = m_saved_tauFakeFactorSyst;
    else{
      ATH_MSG_ERROR("Lepton type " << typeOfLep << " not recognized!");
    }

    if(!hist){
      ATH_MSG_ERROR("Saved fake factor syst histogram for lepton type " << typeOfLep << " cannot be found!");
    }
  }
  else{
    if(typeOfLep == LepEnum::Electron)  hist = m_saved_elFakeFactor;
    else if(typeOfLep == LepEnum::Muon) hist = m_saved_muFakeFactor;
    else if(typeOfLep == LepEnum::Tau)  hist = m_saved_tauFakeFactor;
    else{
      ATH_MSG_ERROR("Lepton type " << typeOfLep << " not recognized!");
    }

    if(!hist){
      ATH_MSG_ERROR("Saved fake factor histogram for lepton type " << typeOfLep << " cannot be found!");
    }
  }

  return hist;

}
TH2F* ApplyFakeFactor::get2DFakeFactorHist(LepEnum::LepType typeOfLep, bool getSystErrHist){

  TH2F* hist = 0;

  if(getSystErrHist){
    if(typeOfLep == LepEnum::Electron)  hist = m_saved_elFakeFactorSyst_2D;
    else if(typeOfLep == LepEnum::Muon) hist = m_saved_muFakeFactorSyst_2D;
    else if(typeOfLep == LepEnum::Tau)  hist = m_saved_tauFakeFactorSyst_2D;
    else{
      ATH_MSG_ERROR("Lepton type " << typeOfLep << " not recognized!");
    }

    if(!hist){
      ATH_MSG_ERROR("Saved fake factor syst histogram for lepton type " << typeOfLep << " cannot be found!");
    }
  }
  else{
    if(typeOfLep == LepEnum::Electron)  hist = m_saved_elFakeFactor_2D;
    else if(typeOfLep == LepEnum::Muon) hist = m_saved_muFakeFactor_2D;
    else if(typeOfLep == LepEnum::Tau)  hist = m_saved_tauFakeFactor_2D;
    else{
      ATH_MSG_ERROR("Lepton type " << typeOfLep << " not recognized!");
    }

    if(!hist){
      ATH_MSG_ERROR("Saved fake factor histogram for lepton type " << typeOfLep << " cannot be found!");
    }
  }

  return hist;
}

TH1F* ApplyFakeFactor::Load1DHistogram(TFile* fakeFactorFile, TString name){

  // if the hist name was not specified, don't try to load it into memory
  if(name == "" || name == "_Syst"){
    return 0;
  }

  TH1F* hist = (TH1F*)fakeFactorFile->Get(name);
  if(!hist){
    ATH_MSG_WARNING("Unable to open " << name << " histogram in " << m_savedFakeFactorFileName << ". Skipping.");
  }
  else{
    hist->SetDirectory(0);
  }

  return hist;
}

TH2F* ApplyFakeFactor::Load2DHistogram(TFile* fakeFactorFile, TString name){

  // if the hist name was not specified, don't try to load it into memory
  if(name == "" || name == "_Syst"){
    return 0;
  }

  TH2F* hist = (TH2F*)fakeFactorFile->Get(name);
  if(!hist){
    ATH_MSG_WARNING("Unable to open " << name << " histogram in " << m_savedFakeFactorFileName << ". Skipping.");
  }
  else{
    hist->SetDirectory(0);
  }

  return hist;
}
