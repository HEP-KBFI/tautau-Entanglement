
#include "DataFormats/FWLite/interface/InputSource.h"                                           // fwlite::InputSource
#include "DataFormats/FWLite/interface/OutputFiles.h"                                           // fwlite::OutputFiles
#include "FWCore/ParameterSet/interface/ParameterSet.h"                                         // edm::ParameterSet
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"                             // edm::readPSetsFrom()
#include "FWCore/PluginManager/interface/PluginManager.h"                                       // edmplugin::PluginManager::configure()
#include "FWCore/PluginManager/interface/standard.h"                                            // edmplugin::standard::config()
#include "PhysicsTools/FWLite/interface/TFileService.h"                                         // fwlite::TFileService

#include "TauAnalysis/Entanglement/interface/cmsException.h"                                    // cmsException

#include <TBenchmark.h>                                                                         // TBenchmark
#include <TError.h>                                                                             // gErrorAbortLevel, kError
#include <TTree.h>                                                                              // TTree
#include <TMatrixD.h>                                                                           // TMatrixD
#include <TObject.h>                                                                            // TObject
#include <TMinuit.h>                                                                            // TMinuit
#include <TString.h>                                                                            // Form()

#include <assert.h>                                                                             // assert()
#include <cstdlib>                                                                              // EXIT_SUCCESS, EXIT_FAILURE
#include <fstream>                                                                              // std::ofstream
#include <iostream>                                                                             // std::cout
#include <string>                                                                               // std::string
#include <vector>                                                                               // std::vector<>
#include <cmath>                                                                                // std::fabs()

class EntanglementData
{
 public:
  EntanglementData(float hPlus_r, float hPlus_n, float hPlus_k, float hMinus_r, float hMinus_n, float hMinus_k, float evtWeight)
    : hPlus_r_(hPlus_r)
    , hPlus_n_(hPlus_n)
    , hPlus_k_(hPlus_k)
    , hMinus_r_(hMinus_r)
    , hMinus_n_(hMinus_n)
    , hMinus_k_(hMinus_k)
    , evtWeight_(evtWeight)
  {}
  ~EntanglementData()
  {}

  inline
  float
  get_hPlus_r() const
  {
    return hPlus_r_;
  }
  inline
  float
  get_hPlus_n() const
  {
    return hPlus_n_;
  }
  inline
  float
  get_hPlus_k() const
  {
    return hPlus_k_;
  }
  inline
  float
  get_hMinus_r() const
  {
    return hMinus_r_;
  }
  inline
  float
  get_hMinus_n() const
  {
    return hMinus_n_;
  }
  inline
  float
  get_hMinus_k() const
  {
    return hMinus_k_;
  }

  inline
  float
  get_evtWeight() const
  {
    return evtWeight_;
  }

 private:
  float hPlus_r_;
  float hPlus_n_;
  float hPlus_k_;
  float hMinus_r_;
  float hMinus_n_;
  float hMinus_k_;

  float evtWeight_;
};

class EntanglementDataset : public TObject
{
 public:
  EntanglementDataset()
  {}
  ~EntanglementDataset()
  {}

  void
  push_back(float hPlus_r, float hPlus_n, float hPlus_k, float hMinus_r, float hMinus_n, float hMinus_k, float evtWeight = 1.)
  {
    data_.push_back(EntanglementData(hPlus_r, hPlus_n, hPlus_k, hMinus_r, hMinus_n, hMinus_k, evtWeight));
  }

  inline
  size_t
  size() const
  {
    return data_.size();
  }

  inline
  const EntanglementData&
  operator[](size_t idx) const
  {
    return data_[idx];
  }

  inline
  const EntanglementData&
  at(size_t idx) const
  {
    return data_[idx];
  }

 private:
  std::vector<EntanglementData> data_;
};

void mlfit_fcn(int& npar, double* gin, double& f, double* par, int iflag)
{
  assert(npar == 15);
  double Bp_r = par[0];
  double Bp_n = par[1];
  double Bp_k = par[2];

  double Bm_r = par[3];
  double Bm_n = par[4];
  double Bm_k = par[5];

  double C_rr = par[6];
  double C_rn = par[7];
  double C_rk = par[8];
  double C_nr = par[9];
  double C_nn = par[10];
  double C_nk = par[11];
  double C_kr = par[12];
  double C_kn = par[13];
  double C_kk = par[14];

  const EntanglementDataset* mlfitData = dynamic_cast<EntanglementDataset*>(gMinuit->GetObjectFit());
  assert(mlfitData);

  double logL = 0.;

  size_t numEntries = mlfitData->size();
  for ( size_t idxEntry = 0; idxEntry < numEntries; ++idxEntry )
  {
    const EntanglementData& entry = mlfitData->at(idxEntry);

    // CV: Equation for probability p of tau pair to be in given spin state
    //     taken from Eq. (2.6) of the paper Comput.Phys.Commun. 64 (1990) 275.
    //     The signs for the terms Bp_X and C_XX, in which the tau+ enters, are flipped, 
    //     because helicity frame defined with respect to tau- (and not tau+) direction in EntanglementNtupleProducer.cc code.
    double p = 1. 
      + (Bp_r*entry.get_hPlus_r() + Bp_n*entry.get_hPlus_n()  + Bp_k*entry.get_hPlus_k())
      - (Bm_r*entry.get_hMinus_r() + Bm_n*entry.get_hMinus_n() + Bm_k*entry.get_hMinus_k())
      - (C_rr*entry.get_hPlus_r()*entry.get_hMinus_r() + C_rn*entry.get_hPlus_r()*entry.get_hMinus_n() + C_rk*entry.get_hPlus_r()*entry.get_hMinus_k())
      - (C_nr*entry.get_hPlus_n()*entry.get_hMinus_r() + C_nn*entry.get_hPlus_n()*entry.get_hMinus_n() + C_nk*entry.get_hPlus_n()*entry.get_hMinus_k())
      - (C_kr*entry.get_hPlus_k()*entry.get_hMinus_r() + C_kn*entry.get_hPlus_k()*entry.get_hMinus_n() + C_kk*entry.get_hPlus_k()*entry.get_hMinus_k());
    const double epsilon = 1.e-12;
    if ( p < epsilon ) p = epsilon;

    logL -= 2.*entry.get_evtWeight()*log(p);
    //logL -= 2.*log(p);
  }

  f = logL;
}

int main(int argc, char* argv[])
{
//--- throw an exception in case ROOT encounters an error
  gErrorAbortLevel = kError;

//--- stop ROOT from keeping track of all histograms
  TH1::AddDirectory(false);

//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]\n";
    return EXIT_FAILURE;
  }

  std::cout << "<analyzeEntanglementNtuple>:\n";

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("analyzeEntanglementNtuple");

//--- read python configuration parameters
  std::cout << "Reading config file " << argv[1] << "\n";
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") )
    throw cmsException("analyzeEntanglementNtuple", __LINE__) << "No ParameterSet 'process' found in config file !!";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameterSet("process");

  edm::ParameterSet cfg_input = cfg.getParameterSet("fwliteInput");
  int maxEvents_beforeCuts = cfg_input.getParameter<int>("maxEvents_beforeCuts");
  int maxEvents_afterCuts = cfg_input.getParameter<int>("maxEvents_afterCuts");
  std::cout << " maxEvents: beforeCuts = " << maxEvents_beforeCuts << ", afterCuts = " << maxEvents_afterCuts << "\n";

  edm::ParameterSet cfg_analyze = cfg.getParameterSet("analyzeEntanglementNtuple");
  std::string treeName = cfg_analyze.getParameter<std::string>("treeName");
  std::cout << " treeName = " << treeName << "\n";
  std::string mode = cfg_analyze.getParameter<std::string>("mode");
  std::cout << " mode = " << mode << "\n";
  float minVisTauPt = cfg_analyze.getParameter<double>("minVisTauPt");
  std::cout << " minVisTauPt = " << minVisTauPt << "\n";
  float maxAbsVisTauEta = cfg_analyze.getParameter<double>("maxAbsVisTauEta");
  std::cout << " maxAbsVisTauEta = " << maxAbsVisTauEta << "\n";
  int maxNumChargedKaons = cfg_analyze.getParameter<int>("maxNumChargedKaons");
  std::cout << " maxNumChargedKaons = " << maxNumChargedKaons << "\n";
  int maxNumNeutralKaons = cfg_analyze.getParameter<int>("maxNumNeutralKaons");
  std::cout << " maxNumNeutralKaons = " << maxNumNeutralKaons << "\n";
  int maxNumPhotons = cfg_analyze.getParameter<int>("maxNumPhotons");
  std::cout << " maxNumPhotons = " << maxNumPhotons << "\n";
  float maxSumPhotonEn = cfg_analyze.getParameter<double>("maxSumPhotonEn");
  std::cout << " maxSumPhotonEn = " << maxSumPhotonEn << "\n";
  std::string branchName_evtWeight = cfg_analyze.getParameter<std::string>("branchName_evtWeight");
  std::cout << " branchName_evtWeight = " << branchName_evtWeight << "\n";
  //bool isDEBUG = cfg_analyze.getParameter<bool>("isDEBUG");

  fwlite::InputSource inputFiles(cfg);
  unsigned reportEvery = inputFiles.reportAfter();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().c_str());

  std::vector<std::string> inputFileNames = inputFiles.files();
  size_t numInputFiles = inputFileNames.size();
  std::cout << "Loaded " << numInputFiles << " file(s).\n";

  TDirectory* dir = fs.getBareDirectory();
  dir->cd();

  TMatrixD C(3, 3);

  EntanglementDataset mlfitData; 

  int analyzedEntries = 0;
  double analyzedEntries_weighted = 0.;
  int selectedEntries = 0;
  double selectedEntries_weighted = 0.;
  int processedInputFiles = 0;
  bool STOP = false;
  for ( size_t idxInputFile = 0; idxInputFile < numInputFiles && !STOP; ++idxInputFile )
  {
    const std::string & inputFileName = inputFileNames.at(idxInputFile);
    std::cout << "Opening #" << idxInputFile << " file " << inputFileName << '\n';
    TFile* inputFile = TFile::Open(inputFileName.c_str());
    if ( !inputFile )
      throw cmsException("analyzeEntanglementNtuple", __LINE__) 
        << "The file " << inputFileName << " failed to open !!";
   
    TTree* inputTree = dynamic_cast<TTree*>(inputFile->Get(treeName.c_str()));
    if ( !inputTree )
      throw cmsException("analyzeEntanglementNtuple", __LINE__) 
        << "The file " << inputFileName << " does not contain a TTree named '" << treeName << "' !!";
    std::cout << "The file " << inputFileName << " contains " << inputTree->GetEntries() << " entries\n";

    ++processedInputFiles;

    Float_t hPlus_r, hPlus_n, hPlus_k;
    inputTree->SetBranchAddress(Form("%s_hPlus_r", mode.c_str()), &hPlus_r);
    inputTree->SetBranchAddress(Form("%s_hPlus_n", mode.c_str()), &hPlus_n);
    inputTree->SetBranchAddress(Form("%s_hPlus_k", mode.c_str()), &hPlus_k);
    Float_t visTauPlus_pt, visTauPlus_eta;
    inputTree->SetBranchAddress(Form("%s_visPlus_pt", mode.c_str()), &visTauPlus_pt);
    inputTree->SetBranchAddress(Form("%s_visPlus_eta", mode.c_str()), &visTauPlus_eta);
    Int_t tauPlus_nChargedKaons, tauPlus_nNeutralKaons, tauPlus_nPhotons;
    Float_t tauPlus_sumPhotonEn;
    inputTree->SetBranchAddress("gen_tauPlus_nChargedKaons", &tauPlus_nChargedKaons);
    inputTree->SetBranchAddress("gen_tauPlus_nNeutralKaons", &tauPlus_nNeutralKaons);
    inputTree->SetBranchAddress("gen_tauPlus_nPhotons", &tauPlus_nPhotons);
    inputTree->SetBranchAddress("gen_tauPlus_sumPhotonEn", &tauPlus_sumPhotonEn);

    Float_t hMinus_r, hMinus_n, hMinus_k;
    inputTree->SetBranchAddress(Form("%s_hMinus_r", mode.c_str()), &hMinus_r);
    inputTree->SetBranchAddress(Form("%s_hMinus_n", mode.c_str()), &hMinus_n);
    inputTree->SetBranchAddress(Form("%s_hMinus_k", mode.c_str()), &hMinus_k);
    Float_t visTauMinus_pt, visTauMinus_eta;
    inputTree->SetBranchAddress(Form("%s_visMinus_pt", mode.c_str()), &visTauMinus_pt);
    inputTree->SetBranchAddress(Form("%s_visMinus_eta", mode.c_str()), &visTauMinus_eta);
    Int_t tauMinus_nChargedKaons, tauMinus_nNeutralKaons, tauMinus_nPhotons;
    Float_t tauMinus_sumPhotonEn;
    inputTree->SetBranchAddress("gen_tauMinus_nChargedKaons", &tauMinus_nChargedKaons);
    inputTree->SetBranchAddress("gen_tauMinus_nNeutralKaons", &tauMinus_nNeutralKaons);
    inputTree->SetBranchAddress("gen_tauMinus_nPhotons", &tauMinus_nPhotons);
    inputTree->SetBranchAddress("gen_tauMinus_sumPhotonEn", &tauMinus_sumPhotonEn);

    Float_t evtWeight = 1.;
    if ( branchName_evtWeight != "" )
    {
      inputTree->SetBranchAddress(branchName_evtWeight.c_str(), &evtWeight);
    }

    int numEntries = inputTree->GetEntries();
    for ( int idxEntry = 0; idxEntry < numEntries && !STOP; ++idxEntry )
    {
      inputTree->GetEntry(idxEntry);

      ++analyzedEntries;
      analyzedEntries_weighted += evtWeight;
      if ( (analyzedEntries % reportEvery) == 0 )
      {
        std::cout << "processing Entry " << analyzedEntries << "\n";
      }

      if ( !(visTauPlus_pt  > minVisTauPt && std::fabs(visTauPlus_eta)  < maxAbsVisTauEta) ) continue;
      if ( maxNumChargedKaons != -1  && tauPlus_nChargedKaons  > maxNumChargedKaons ) continue;
      if ( maxNumNeutralKaons != -1  && tauPlus_nNeutralKaons  > maxNumNeutralKaons ) continue;
      if ( maxNumPhotons      != -1  && tauPlus_nPhotons       > maxNumPhotons      ) continue;
      if ( maxSumPhotonEn     >=  0. && tauPlus_sumPhotonEn    > maxSumPhotonEn     ) continue;
      if ( !(visTauMinus_pt > minVisTauPt && std::fabs(visTauMinus_eta) < maxAbsVisTauEta) ) continue;
      if ( maxNumChargedKaons != -1  && tauMinus_nChargedKaons > maxNumChargedKaons ) continue;
      if ( maxNumNeutralKaons != -1  && tauMinus_nNeutralKaons > maxNumNeutralKaons ) continue;
      if ( maxNumPhotons      != -1  && tauMinus_nPhotons      > maxNumPhotons      ) continue;
      if ( maxSumPhotonEn     >=  0. && tauMinus_sumPhotonEn   > maxSumPhotonEn     ) continue;

      // CV: compute matrix C according to Eq. (25)
      //     in the paper arXiv:2211.10513
      double c = -9.*evtWeight;
      C[0][0] += c*hPlus_r*hMinus_r;
      C[0][1] += c*hPlus_r*hMinus_n;
      C[0][2] += c*hPlus_r*hMinus_k;
      C[1][0] += c*hPlus_n*hMinus_r;
      C[1][1] += c*hPlus_n*hMinus_n;
      C[1][2] += c*hPlus_n*hMinus_k;
      C[2][0] += c*hPlus_k*hMinus_r;
      C[2][1] += c*hPlus_k*hMinus_n;
      C[2][2] += c*hPlus_k*hMinus_k;
 
      mlfitData.push_back(hPlus_r, hPlus_n, hPlus_k, hMinus_r, hMinus_n, hMinus_k, evtWeight);

      ++selectedEntries;
      selectedEntries_weighted += evtWeight;

      if ( maxEvents_beforeCuts != -1 && analyzedEntries >= maxEvents_beforeCuts ) STOP = true;
      if ( maxEvents_afterCuts  != -1 && selectedEntries >= maxEvents_afterCuts  ) STOP = true;
    }

    delete inputTree;
    delete inputFile;
  }

  C *= (1./selectedEntries_weighted);

  std::cout << "Processing Summary:\n";
  std::cout << " processedInputFiles = " << processedInputFiles << " (out of " << numInputFiles << ")\n";
  std::cout << " analyzedEntries = " << analyzedEntries << " (weighted = " << analyzedEntries_weighted << ")\n";
  std::cout << " selectedEntries = " << selectedEntries << " (weighted = " << selectedEntries_weighted << ")\n";

  std::cout << "Matrix C (measured using Eq. (25) in arXiv:2211.10513):\n";
  C.Print();

  std::cout << "Standard Model expectation (given by Eq. (69) in arXiv:2208:11723):\n";
  TMatrixD C_exp(3, 3);
  C_exp[0][0] = +1.;
  C_exp[1][1] = +1.;
  C_exp[2][2] = -1.;
  C_exp.Print();

  // define parameters for maximum-likelihood (ML) fit
  std::map<size_t, std::string> parNames;
  parNames[0]  = "Bp_r";
  parNames[1]  = "Bp_n";
  parNames[2]  = "Bp_k";
  parNames[3]  = "Bm_r";
  parNames[4]  = "Bm_n";
  parNames[5]  = "Bm_k";
  parNames[6]  = "C_rr";
  parNames[7]  = "C_rn";
  parNames[8]  = "C_rk";
  parNames[9]  = "C_nr";
  parNames[10] = "C_nn";
  parNames[11] = "C_nk";
  parNames[12] = "C_kr";
  parNames[13] = "C_kn";
  parNames[14] = "C_kk";

  // initialize Minuit
  TMinuit mlfit(15);
  for ( size_t idxPar = 0; idxPar < 15; ++idxPar )
  {
    double par0 = 0.;
    if ( idxPar >= 6 && idxPar <= 14 )
    {
      size_t idxRow = (idxPar - 6) / 3;
      size_t idxCol = (idxPar - 6) % 3;
      par0 = C[idxRow][idxCol];
    }
    mlfit.DefineParameter(idxPar, parNames[idxPar].c_str(), par0, 0.1, -2., +2.);
  }

  // set function pointer
  mlfit.SetFCN(mlfit_fcn);

  // set minimization strategy: 1 = standard, 2 = try to improve minimum (slower)
  double arglist[10];
  arglist[0] = 2;
  int errorFlag = 0;
  mlfit.mnexcm("SET STR", arglist, 1, errorFlag);

  // set output level: -1 = no output, 1 = standard output
  mlfit.SetPrintLevel(-1);

  // disable warnings
  mlfit.mnexcm("SET NOW", arglist, 1, errorFlag);

  // set polarimetric vectors of tau+ and tau- in processed dataset
  gMinuit = &mlfit;
  gMinuit->SetObjectFit(&mlfitData);

  mlfit.SetMaxIterations(10000);
  mlfit.Migrad();

  std::cout << "Fit Results:\n";
  for ( size_t idxPar = 0; idxPar < 15; ++idxPar )
  {
    double parValue, parError;
    mlfit.GetParameter(idxPar, parValue, parError);
    std::cout << parNames[idxPar] << " = " << parValue << " +/- " << parError << "\n";
  }

  mlfit.mnhess();
  std::cout << "Hesse matrix:\n";
  mlfit.mnmatu(2);

  clock.Show("analyzeEntanglementNtuple");

  return EXIT_SUCCESS;
}
