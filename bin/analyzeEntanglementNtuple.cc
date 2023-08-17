
#include "DataFormats/FWLite/interface/InputSource.h"                             // fwlite::InputSource
#include "DataFormats/FWLite/interface/OutputFiles.h"                             // fwlite::OutputFiles
#include "FWCore/ParameterSet/interface/ParameterSet.h"                           // edm::ParameterSet
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"               // edm::readPSetsFrom()
#include "FWCore/PluginManager/interface/PluginManager.h"                         // edmplugin::PluginManager::configure()
#include "FWCore/PluginManager/interface/standard.h"                              // edmplugin::standard::config()
#include "PhysicsTools/FWLite/interface/TFileService.h"                           // fwlite::TFileService

#include "TauAnalysis/Entanglement/interface/cmsException.h"                      // cmsException
#include "TauAnalysis/Entanglement/interface/comp_EigenVectors_and_EigenValues.h" // comp_EigenVectors_and_EigenValues()
#include "TauAnalysis/Entanglement/interface/format_vT.h"                         // format_vint(), vdouble, vint
#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h"                 // Matrix3x3, Vector3
#include "TauAnalysis/Entanglement/interface/passesStatusSelection.h"             // passesStatusSelection()
#include "TauAnalysis/Entanglement/interface/printEigenVectors_and_EigenValues.h" // printEigenVectors_and_EigenValues()

#include "Math/Functions.h"                                                       // ROOT::Math::Transpose() 
#include <TBenchmark.h>                                                           // TBenchmark
#include <TError.h>                                                               // gErrorAbortLevel, kError
#include <TMath.h>                                                                // TMath::Nint()
#include <TObject.h>                                                              // TObject
#include <TRandom3.h>                                                             // TRandom3
#include <TString.h>                                                              // Form()
#include <TTree.h>                                                                // TTree

#include <algorithm>                                                              // std::sort()
#include <assert.h>                                                               // assert()
#include <cmath>                                                                  // std::fabs()
#include <cstdlib>                                                                // EXIT_SUCCESS, EXIT_FAILURE
#include <fstream>                                                                // std::ofstream
#include <iostream>                                                               // std::cout
#include <string>                                                                 // std::string
#include <utility>                                                                // std::make_pair(), std::pair<>
#include <vector>                                                                 // std::vector<>

const size_t npar = 15;

typedef std::vector<double> vdouble;
typedef std::map<int, vdouble> map_vdouble;
typedef std::map<int, map_vdouble> map2_vdouble;

class EntanglementData
{
 public:
  EntanglementData(float hPlus_r, float hPlus_n, float hPlus_k, 
                   float hMinus_r, float hMinus_n, float hMinus_k,
                   float evtWeight)
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
  EntanglementDataset(const std::vector<double>& par_gen)
    : par_gen_(par_gen)
  {}
  EntanglementDataset(const EntanglementDataset& dataset, int maxEvents_afterCuts = -1)
    : par_gen_(dataset.par_gen_)
  {
    // CV: If maxEvents_afterCuts == 0, an empty dataset will be created.
    //     This is the intended behaviour and is used for building bootstrap samples. 
    if ( maxEvents_afterCuts != 0 )
    {
      if ( maxEvents_afterCuts > 0 && maxEvents_afterCuts < (int)dataset.data_.size() )
      {
        size_t numEntries = dataset.data_.size();
        for ( size_t idxEntry = 0; idxEntry < numEntries; ++idxEntry )
        {
          const EntanglementData& entry = dataset.data_.at(idxEntry);
          data_.push_back(entry);
        }
      }
      else
      {
        data_ = dataset.data_;
      }
    }
  }
  ~EntanglementDataset()
  {}

  void
  push_back(const EntanglementData& entry)
  {
    data_.push_back(entry);
  }

  inline
  size_t
  size() const
  {
    return data_.size();
  }

  inline
  const EntanglementData&
  at(size_t idx) const
  {
    return data_[idx];
  }

  inline
  double
  get_par_gen(size_t idx) const
  {
    return par_gen_[idx];
  }

 private:
   std::vector<EntanglementData> data_;

   std::vector<double> par_gen_;
};

EntanglementDataset
build_bootstrap_sample(const EntanglementDataset& dataset, TRandom& rnd, int maxEvents_afterCuts = -1)
{
  EntanglementDataset bootstrap_sample(dataset, 0);
  size_t sampleSize = ( maxEvents_afterCuts > 0 ) ? maxEvents_afterCuts : dataset.size();
  for ( size_t idxSample = 0; idxSample < sampleSize; ++idxSample )
  {
    size_t idxEntry = rnd.Integer(dataset.size());
    const EntanglementData& entry = dataset.at(idxEntry);
    bootstrap_sample.push_back(entry);
  }
  return bootstrap_sample;
}

double
comp_Rchsh(const math::Matrix3x3& C, int verbosity = -1)
{
  if ( verbosity >= 1 )
  {
    std::cout << "<comp_Rchsh>:\n";
  }
  math::Matrix3x3 CT = ROOT::Math::Transpose(C);
  math::Matrix3x3 CT_times_C = CT*C;
  std::vector<std::pair<TVectorD, double>> EigenVectors_and_EigenValues = comp_EigenVectors_and_EigenValues(CT_times_C);
  if ( verbosity >= 1 )
  {
    printEigenVectors_and_EigenValues(EigenVectors_and_EigenValues);
  }
  assert(EigenVectors_and_EigenValues.size() == 3);
  double Rchsh = EigenVectors_and_EigenValues[0].second + EigenVectors_and_EigenValues[1].second;
  if ( verbosity >= 1 )
  {
    std::cout << "Rchsh = " << Rchsh << "\n";
  }
  return Rchsh;
}

class Measurement
{
 public:
  Measurement(const math::Vector3& Bp, const math::Vector3& BpErr,
              const math::Vector3& Bm, const math::Vector3& BmErr,
              const math::Matrix3x3& C, const math::Matrix3x3& CErr,
              int verbosity = -1)
    : Bp_(Bp)
    , BpErr_(BpErr)
    , Bm_(Bm)
    , BmErr_(BmErr)
    , C_(C)
    , CErr_(CErr)
    , Rchsh_(0.)
  {
    Rchsh_ = comp_Rchsh(C, verbosity);
  }
  ~Measurement()
  {}

  void
  set_BpErr(const math::Vector3& BpErr)
  {
    BpErr_ = BpErr;
  }
  void
  set_BmErr(const math::Vector3& BmErr)
  {
    BmErr_ = BmErr;
  }
  void
  set_CErr(const math::Matrix3x3& CErr)
  {
    CErr_ = CErr;
  }
  void
  set_RchshErr(double Rchsh)
  {
    Rchsh_ = Rchsh;
  }

  const math::Vector3&
  get_Bp() const
  {
    return Bp_;
  }
  const math::Vector3&
  get_BpErr() const
  {
    return BpErr_;
  }
  const math::Vector3&
  get_Bm() const
  {
    return Bm_;
  }
  const math::Vector3&
  get_BmErr() const
  {
    return BmErr_;
  }

  const math::Matrix3x3&
  get_C() const
  {
    return C_;
  }
  const math::Matrix3x3&
  get_CErr() const
  {
    return CErr_;
  }

  double
  get_Rchsh() const
  {
    return Rchsh_;
  }
  double
  get_RchshErr() const
  {
    return RchshErr_;
  }

 private:
  math::Vector3 Bp_;
  math::Vector3 BpErr_;
  math::Vector3 Bm_;
  math::Vector3 BmErr_;
  math::Matrix3x3 C_;
  math::Matrix3x3 CErr_;
  double Rchsh_;
  double RchshErr_;
};

Measurement
comp_Bp_Bm_C(const EntanglementDataset& dataset, 
             int verbosity = -1)
{
  math::Vector3 Bp;
  math::Vector3 BpErr;
  math::Vector3 Bm;
  math::Vector3 BmErr;
  math::Matrix3x3 C;
  math::Matrix3x3 CErr;

  double evtWeight_sum = 0.;

  size_t numEntries = dataset.size();
  for ( size_t idxEntry = 0; idxEntry < numEntries; ++idxEntry )
  {
    const EntanglementData& entry = dataset.at(idxEntry);

    double hPlus_r = entry.get_hPlus_r();
    double hPlus_n = entry.get_hPlus_n();
    double hPlus_k = entry.get_hPlus_k();
    
    double hMinus_r = entry.get_hMinus_r();
    double hMinus_n = entry.get_hMinus_n();
    double hMinus_k = entry.get_hMinus_k();

    double evtWeight = entry.get_evtWeight();

    Bp(0) += evtWeight*hPlus_r;
    Bp(1) += evtWeight*hPlus_n;
    Bp(2) += evtWeight*hPlus_k;

    Bm(0) += evtWeight*hMinus_r;
    Bm(1) += evtWeight*hMinus_n;
    Bm(2) += evtWeight*hMinus_k;
    
    // CV: compute matrix C according to Eq. (25)
    //     in the paper arXiv:2211.10513
    double c = -9.*evtWeight;
    C(0,0) += c*hPlus_r*hMinus_r;
    C(0,1) += c*hPlus_r*hMinus_n;
    C(0,2) += c*hPlus_r*hMinus_k;
    C(1,0) += c*hPlus_n*hMinus_r;
    C(1,1) += c*hPlus_n*hMinus_n;
    C(1,2) += c*hPlus_n*hMinus_k;
    C(2,0) += c*hPlus_k*hMinus_r;
    C(2,1) += c*hPlus_k*hMinus_n;
    C(2,2) += c*hPlus_k*hMinus_k;

    evtWeight_sum += evtWeight;
  }

  Bp *= (1./evtWeight_sum);
  Bm *= (1./evtWeight_sum);

  C *= (1./evtWeight_sum);

  Measurement measurement(Bp, BpErr, Bm, BmErr, C, CErr, verbosity);
  return measurement;
}

std::pair<double, double>
comp_median_and_Err(const std::vector<double>& measuredValues)
{
  std::vector<double> tmp = measuredValues;
  // CV: sort measured values into ascending order
  std::sort(tmp.begin(), tmp.end());
  size_t numMeasurements = measuredValues.size();
  int idxMedian = TMath::Nint(0.5*numMeasurements);
  double median = tmp[idxMedian];
  int idxPlus1Sigma = TMath::Nint(0.84*numMeasurements);
  int idxMinus1Sigma = TMath::Nint(0.16*numMeasurements);
  double Err = tmp[idxPlus1Sigma] - tmp[idxMinus1Sigma];
  assert(Err >= 0.);
  return std::make_pair(median, Err);
}

std::pair<math::Vector3, math::Vector3>
comp_median_and_Err(const std::vector<math::Vector3>& measuredVectors)
{
  map_vdouble tmp;
  size_t numMeasurements = measuredVectors.size();
  for ( size_t idxMeasurement = 0; idxMeasurement < numMeasurements; ++idxMeasurement )
  {
    const math::Vector3& measuredVector = measuredVectors[idxMeasurement];
    for ( int idxElement = 0; idxElement < 3; ++idxElement )
    {
      tmp[idxElement].push_back(measuredVector(idxElement));
    }
  }
  math::Vector3 median;
  math::Vector3 Err;
  for ( size_t idxElement = 0; idxElement < 3; ++idxElement )
  {
    std::pair<double,double> median_and_Err = comp_median_and_Err(tmp[idxElement]);
    median(idxElement) = median_and_Err.first;
    Err(idxElement) = median_and_Err.second;
  }
  return std::make_pair(median, Err);
}

std::pair<math::Matrix3x3, math::Matrix3x3>
comp_median_and_Err(const std::vector<math::Matrix3x3>& measuredMatrices)
{
  map2_vdouble tmp;
  size_t numMeasurements = measuredMatrices.size();
  for ( size_t idxMeasurement = 0; idxMeasurement < numMeasurements; ++idxMeasurement )
  {
    const math::Matrix3x3& measuredMatrix = measuredMatrices[idxMeasurement];
    for ( int idxRow = 0; idxRow < 3; ++idxRow )
    {
      for ( int idxColumn = 0; idxColumn < 3; ++idxColumn )
      {
        tmp[idxRow][idxColumn].push_back(measuredMatrix(idxRow,idxColumn));
      }
    }
  }
  math::Matrix3x3 median;
  math::Matrix3x3 Err;
  for ( int idxRow = 0; idxRow < 3; ++idxRow )
  {
    for ( int idxColumn = 0; idxColumn < 3; ++idxColumn )
    {
      std::pair<double,double> median_and_Err = comp_median_and_Err(tmp[idxRow][idxColumn]);
      median(idxRow,idxColumn) = median_and_Err.first;
      Err(idxRow,idxColumn) = median_and_Err.second;
    }
  }
  return std::make_pair(median, Err);
}

void
comp_median_and_Err(const std::vector<Measurement>& measurements,
                    math::Vector3& Bp_median, math::Vector3& BpErr,
                    math::Vector3& Bm_median, math::Vector3& BmErr,
                    math::Matrix3x3& C_median, math::Matrix3x3& CErr,
                    double& Rchsh_median, double& RchshErr)
{
  std::vector<math::Vector3> measuredBp;
  std::vector<math::Vector3> measuredBm;
  std::vector<math::Matrix3x3> measuredC;
  std::vector<double> measuredRchsh;
  for ( const Measurement& measurement : measurements )
  {
    measuredBp.push_back(measurement.get_Bp());
    measuredBm.push_back(measurement.get_Bm());
    measuredC.push_back(measurement.get_C());
    measuredRchsh.push_back(measurement.get_Rchsh());
  }
  std::pair<math::Vector3, math::Vector3> Bp_median_and_Err = comp_median_and_Err(measuredBp);
  Bp_median = Bp_median_and_Err.first;
  BpErr = Bp_median_and_Err.second;
  std::pair<math::Vector3, math::Vector3> Bm_median_and_Err = comp_median_and_Err(measuredBm);
  Bm_median = Bm_median_and_Err.first;
  BmErr = Bm_median_and_Err.second;
  std::pair<math::Matrix3x3, math::Matrix3x3> C_median_and_Err = comp_median_and_Err(measuredC);
  C_median = C_median_and_Err.first;
  CErr = C_median_and_Err.second;
  std::pair<double, double> Rchsh_median_and_Err = comp_median_and_Err(measuredRchsh);
  Rchsh_median = Rchsh_median_and_Err.first;
  RchshErr = Rchsh_median_and_Err.second;
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
  float minTauTIP = cfg_analyze.getParameter<double>("minTauTIP");
  std::cout << " minTauTIP = " << minTauTIP << "\n";
  int maxNumChargedKaons = cfg_analyze.getParameter<int>("maxNumChargedKaons");
  std::cout << " maxNumChargedKaons = " << maxNumChargedKaons << "\n";
  int maxNumNeutralKaons = cfg_analyze.getParameter<int>("maxNumNeutralKaons");
  std::cout << " maxNumNeutralKaons = " << maxNumNeutralKaons << "\n";
  int maxNumPhotons = cfg_analyze.getParameter<int>("maxNumPhotons");
  std::cout << " maxNumPhotons = " << maxNumPhotons << "\n";
  float maxSumPhotonEn = cfg_analyze.getParameter<double>("maxSumPhotonEn");
  std::cout << " maxSumPhotonEn = " << maxSumPhotonEn << "\n";
  float maxChi2 = cfg_analyze.getParameter<double>("maxChi2");
  std::cout << " maxChi2 = " << maxChi2 << "\n";
  vint statusSelection = cfg_analyze.getParameter<vint>("statusSelection");
  std::cout << " statusSelection = " << format_vint(statusSelection) << "\n";
  std::string branchName_evtWeight = cfg_analyze.getParameter<std::string>("branchName_evtWeight");
  std::cout << " branchName_evtWeight = " << branchName_evtWeight << "\n";
  
  std::vector<double> par_gen = cfg_analyze.getParameter<vdouble>("par_gen");
  if ( par_gen.size() != npar )
    throw cmsException("analyzeEntanglementNtuple", __LINE__) << "Invalid Configuration parameter 'par_gen' !!";

  unsigned numBootstrapSamples = cfg_analyze.getParameter<unsigned>("numBootstrapSamples");
  TRandom3 rnd;

  bool scanLikelihood = cfg_analyze.getParameter<bool>("scanLikelihood");
  std::cout << " scanLikelihood = " << scanLikelihood << "\n";
  
  int verbosity = cfg_analyze.getUntrackedParameter<int>("verbosity");

  fwlite::InputSource inputFiles(cfg);
  unsigned reportEvery = inputFiles.reportAfter();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().c_str());

  std::vector<std::string> inputFileNames = inputFiles.files();
  size_t numInputFiles = inputFileNames.size();
  std::cout << "Loaded " << numInputFiles << " file(s).\n";
  
  EntanglementDataset dataset(par_gen); 

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
    Float_t visPlus_pt, visPlus_eta;
    inputTree->SetBranchAddress(Form("%s_visPlus_pt", mode.c_str()), &visPlus_pt);
    inputTree->SetBranchAddress(Form("%s_visPlus_eta", mode.c_str()), &visPlus_eta);
    Float_t tauPlus_tip;
    inputTree->SetBranchAddress(Form("%s_tauPlus_tip", mode.c_str()), &tauPlus_tip);
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
    Float_t visMinus_pt, visMinus_eta;
    inputTree->SetBranchAddress(Form("%s_visMinus_pt", mode.c_str()), &visMinus_pt);
    inputTree->SetBranchAddress(Form("%s_visMinus_eta", mode.c_str()), &visMinus_eta);
    Float_t tauMinus_tip;
    inputTree->SetBranchAddress(Form("%s_tauMinus_tip", mode.c_str()), &tauMinus_tip);
    Int_t tauMinus_nChargedKaons, tauMinus_nNeutralKaons, tauMinus_nPhotons;
    Float_t tauMinus_sumPhotonEn;
    inputTree->SetBranchAddress("gen_tauMinus_nChargedKaons", &tauMinus_nChargedKaons);
    inputTree->SetBranchAddress("gen_tauMinus_nNeutralKaons", &tauMinus_nNeutralKaons);
    inputTree->SetBranchAddress("gen_tauMinus_nPhotons", &tauMinus_nPhotons);
    inputTree->SetBranchAddress("gen_tauMinus_sumPhotonEn", &tauMinus_sumPhotonEn);

    Float_t kinFit_chi2;
    inputTree->SetBranchAddress("kinFit_chi2", &kinFit_chi2);
    Int_t kinFit_status;
    inputTree->SetBranchAddress("kinFit_status", &kinFit_status);

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

      if ( !(visPlus_pt  > minVisTauPt && std::fabs(visPlus_eta)  < maxAbsVisTauEta) ) continue;
      if ( !(tauPlus_tip > minTauTIP) ) continue;
      if ( maxNumChargedKaons     != -1  && tauPlus_nChargedKaons  > maxNumChargedKaons            ) continue;
      if ( maxNumNeutralKaons     != -1  && tauPlus_nNeutralKaons  > maxNumNeutralKaons            ) continue;
      if ( maxNumPhotons          != -1  && tauPlus_nPhotons       > maxNumPhotons                 ) continue;
      if ( maxSumPhotonEn         >=  0. && tauPlus_sumPhotonEn    > maxSumPhotonEn                ) continue;
      if ( !(visMinus_pt > minVisTauPt && std::fabs(visMinus_eta) < maxAbsVisTauEta) ) continue;
      if ( !(tauMinus_tip > minTauTIP) ) continue;
      if ( maxNumChargedKaons     != -1  && tauMinus_nChargedKaons > maxNumChargedKaons            ) continue;
      if ( maxNumNeutralKaons     != -1  && tauMinus_nNeutralKaons > maxNumNeutralKaons            ) continue;
      if ( maxNumPhotons          != -1  && tauMinus_nPhotons      > maxNumPhotons                 ) continue;
      if ( maxSumPhotonEn         >=  0. && tauMinus_sumPhotonEn   > maxSumPhotonEn                ) continue;
      if ( maxChi2                != -1  && kinFit_chi2            > maxChi2                       ) continue;
      if ( statusSelection.size() >   0  && !passesStatusSelection(kinFit_status, statusSelection) ) continue;

      dataset.push_back(EntanglementData(hPlus_r, hPlus_n, hPlus_k, hMinus_r, hMinus_n, hMinus_k, evtWeight));

      ++selectedEntries;
      selectedEntries_weighted += evtWeight;

      if ( maxEvents_beforeCuts != -1 && analyzedEntries >= maxEvents_beforeCuts ) STOP = true;
    }

    delete inputTree;
    delete inputFile;
  }

  std::cout << "Processing Summary:\n";
  std::cout << " processedInputFiles = " << processedInputFiles << " (out of " << numInputFiles << ")\n";
  std::cout << " analyzedEntries = " << analyzedEntries << " (weighted = " << analyzedEntries_weighted << ")\n";
  std::cout << " selectedEntries = " << selectedEntries << " (weighted = " << selectedEntries_weighted << ")\n";

  EntanglementDataset nominal_sample(dataset, maxEvents_afterCuts);

  Measurement nominal_measurement = comp_Bp_Bm_C(nominal_sample, verbosity);
 
  // CV: estimate uncertainties on Bp and Bm vectors, on tau spin correlation matrix C,
  //     and on Entanglement observables with bootstrap samples
  std::vector<Measurement> bootstrap_measurements;
  for ( size_t idxBootstrapSample = 0; idxBootstrapSample < numBootstrapSamples; ++idxBootstrapSample )
  {
    EntanglementDataset bootstrap_sample = build_bootstrap_sample(dataset, rnd, maxEvents_afterCuts);
    Measurement measurement = comp_Bp_Bm_C(bootstrap_sample, -1);
    bootstrap_measurements.push_back(measurement);
  }

  math::Vector3 Bp_median, BpErr, Bm_median, BmErr;
  math::Matrix3x3 C_median, CErr;
  double Rchsh_median, RchshErr;
  comp_median_and_Err(bootstrap_measurements, 
    Bp_median, BpErr,
    Bm_median, BmErr,
    C_median, CErr,
    Rchsh_median, RchshErr);
  std::cout << "Matrix C (measured using Eq. (25) of arXiv:2211.10513):\n";
  std::cout << nominal_measurement.get_C() << "\n";
  std::cout << "+/-\n";
  std::cout << CErr << "\n";
  std::cout << "Rchsh = " << nominal_measurement.get_Rchsh() << " +/- " << RchshErr << "\n";

  if ( verbosity >= 1 )
  {
    std::cout << "Standard Model expectation (given by Eq. (69) of arXiv:2208:11723):\n";
    TMatrixD C_exp(3,3);
    C_exp[0][0] = +1.;
    C_exp[1][1] = +1.;
    C_exp[2][2] = -1.;
    C_exp.Print();
  }

  clock.Show("analyzeEntanglementNtuple");

  return EXIT_SUCCESS;
}
