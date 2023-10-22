
#include "DataFormats/FWLite/interface/InputSource.h"                    // fwlite::InputSource
#include "DataFormats/FWLite/interface/OutputFiles.h"                    // fwlite::OutputFiles
#include "FWCore/ParameterSet/interface/ParameterSet.h"                  // edm::ParameterSet
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"      // edm::readPSetsFrom()
#include "PhysicsTools/FWLite/interface/TFileService.h"                  // fwlite::TFileService

#include "TauAnalysis/Entanglement/interface/bookHistogram1d.h"          // bookHistogram1d()
#include "TauAnalysis/Entanglement/interface/bookHistogram2d.h"          // bookHistogram2d()
#include "TauAnalysis/Entanglement/interface/cmsException.h"             // cmsException
#include "TauAnalysis/Entanglement/interface/comp_concurrence.h"         // comp_concurrence()
#include "TauAnalysis/Entanglement/interface/comp_Ek.h"                  // comp_Ek()
#include "TauAnalysis/Entanglement/interface/comp_Rchsh.h"               // comp_Rchsh()
#include "TauAnalysis/Entanglement/interface/comp_steerability.h"        // comp_steerability()
#include "TauAnalysis/Entanglement/interface/Dataset.h"                  // spin::Dataset
#include "TauAnalysis/Entanglement/interface/format_vT.h"                // format_vint(), vdouble, vint
#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h"        // math::Matrix3x3, math::Vector3
#include "TauAnalysis/Entanglement/interface/Measurement.h"              // spin::Measurement
#include "TauAnalysis/Entanglement/interface/passesStatusSelection.h"    // passesStatusSelection()
#include "TauAnalysis/Entanglement/interface/SpinAlgo_by_mlfit.h"        // SpinAlgo_by_mlfit
#include "TauAnalysis/Entanglement/interface/SpinAnalyzer.h"             // spin::SpinAnalyzer
#include "TauAnalysis/Entanglement/interface/dumpJSON.h"                 // dumpJSON()
#include "TauAnalysis/Entanglement/interface/BranchAddressInitializer.h" // BranchAddressInitializer

#include <TAxis.h>                                                       // TAxis
#include <TBenchmark.h>                                                  // TBenchmark
#include <TError.h>                                                      // gErrorAbortLevel, kError
#include <TH1.h>                                                         // TH1
#include <TH2.h>                                                         // TH2
#include <TString.h>                                                     // Form(), TString
#include <TTree.h>                                                       // TTree

#include <algorithm>                                                     // std::sort()
#include <assert.h>                                                      // assert()
#include <cmath>                                                         // std::fabs()
#include <cstdlib>                                                       // EXIT_SUCCESS, EXIT_FAILURE, std::system()
#include <filesystem>                                                    // std::filesystem::current_path()
#include <iostream>                                                      // std::cout
#include <string>                                                        // std::string
#include <vector>                                                        // std::vector<>
#include <fstream>                                                       // std::ofstream

const size_t npar = 15;

void
addToOutputFile(fwlite::TFileService& fs, const TH1* histogram)
{
  const TAxis* xAxis = histogram->GetXaxis();
  int numBinsX = xAxis->GetNbins();
  double xMin = xAxis->GetXmin();
  double xMax = xAxis->GetXmax();
  TH1* output = bookHistogram1d(fs, histogram->GetName(), numBinsX, xMin, xMax);
  for ( int idxBinX = 1; idxBinX <= numBinsX; ++idxBinX )
  {
    output->SetBinContent(idxBinX, histogram->GetBinContent(idxBinX));
    output->SetBinError(idxBinX, histogram->GetBinError(idxBinX));
  }
}

void
addToOutputFile(fwlite::TFileService& fs, const TH2* histogram)
{
  const TAxis* xAxis = histogram->GetXaxis();
  int numBinsX = xAxis->GetNbins();
  double xMin = xAxis->GetXmin();
  double xMax = xAxis->GetXmax();
  const TAxis* yAxis = histogram->GetYaxis();
  int numBinsY = yAxis->GetNbins();
  double yMin = yAxis->GetXmin();
  double yMax = yAxis->GetXmax();
  TH2* output = bookHistogram2d(fs, histogram->GetName(), numBinsX, xMin, xMax, numBinsY, yMin, yMax);
  for ( int idxBinX = 1; idxBinX <= numBinsX; ++idxBinX )
  {
    for ( int idxBinY = 1; idxBinY <= numBinsY; ++idxBinY )
    {
      output->SetBinContent(idxBinX, idxBinY, histogram->GetBinContent(idxBinX, idxBinY));
      output->SetBinError(idxBinX, idxBinY, histogram->GetBinError(idxBinX, idxBinY));
    }
  }
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
  float minVisTauZ = cfg_analyze.getParameter<double>("minVisTauZ");
  std::cout << " minVisTauZ = " << minVisTauZ << "\n";
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
  
  bool apply_evtWeight = cfg_analyze.getParameter<bool>("apply_evtWeight");
  const double absCosTheta_cut = cfg_analyze.getParameter<double>("absCosTheta_cut");

  std::vector<double> par_gen = cfg_analyze.getParameter<vdouble>("par_gen");
  if ( par_gen.size() != npar )
    throw cmsException("analyzeEntanglementNtuple", __LINE__) 
      << "Invalid Configuration parameter 'par_gen' !!\n";

  cfg_analyze.addParameter<int>("maxEvents_afterCuts", maxEvents_afterCuts);
  spin::SpinAnalyzer spinAnalyzer(cfg_analyze);
  std::string spinAnalyzer_algo = cfg_analyze.getParameter<std::string>("spinAnalyzer");
  
  int verbosity = cfg_analyze.getUntrackedParameter<int>("verbosity");

  fwlite::InputSource inputFiles(cfg);
  unsigned reportEvery = inputFiles.reportAfter();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().c_str());

  std::vector<std::string> inputFileNames = inputFiles.files();
  size_t numInputFiles = inputFileNames.size();
  std::cout << "Loaded " << numInputFiles << " file(s).\n";
  
  spin::Dataset dataset_passed;
  spin::Dataset dataset_failed;

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
    BranchAddressInitializer bai(inputTree);

    Float_t hPlus_n, hPlus_r, hPlus_k;
    bai.setBranchAddress(hPlus_n, Form("%s_hPlus_n", mode.c_str()));
    bai.setBranchAddress(hPlus_r, Form("%s_hPlus_r", mode.c_str()));
    bai.setBranchAddress(hPlus_k, Form("%s_hPlus_k", mode.c_str()));
    Float_t visPlus_pt, visPlus_eta;
    bai.setBranchAddress(visPlus_pt, Form("%s_visPlus_pt", mode.c_str()));
    bai.setBranchAddress(visPlus_eta, Form("%s_visPlus_eta", mode.c_str()));
    Float_t tauPlus_tip;
    bai.setBranchAddress(tauPlus_tip, Form("%s_tauPlus_tip", mode.c_str()));
    Int_t tauPlus_nChargedKaons, tauPlus_nNeutralKaons, tauPlus_nPhotons;
    Float_t tauPlus_sumPhotonEn;
    bai.setBranchAddress(tauPlus_nChargedKaons, "gen_tauPlus_nChargedKaons");
    bai.setBranchAddress(tauPlus_nNeutralKaons, "gen_tauPlus_nNeutralKaons");
    bai.setBranchAddress(tauPlus_nPhotons, "gen_tauPlus_nPhotons");
    bai.setBranchAddress(tauPlus_sumPhotonEn, "gen_tauPlus_sumPhotonEn");
    
    Float_t hMinus_n, hMinus_r, hMinus_k;
    bai.setBranchAddress(hMinus_n, Form("%s_hMinus_n", mode.c_str()));
    bai.setBranchAddress(hMinus_r, Form("%s_hMinus_r", mode.c_str()));
    bai.setBranchAddress(hMinus_k, Form("%s_hMinus_k", mode.c_str()));
    Float_t visMinus_pt, visMinus_eta;
    bai.setBranchAddress(visMinus_pt, Form("%s_visMinus_pt", mode.c_str()));
    bai.setBranchAddress(visMinus_eta, Form("%s_visMinus_eta", mode.c_str()));
    Float_t tauMinus_tip;
    bai.setBranchAddress(tauMinus_tip, Form("%s_tauMinus_tip", mode.c_str()));
    Int_t tauMinus_nChargedKaons, tauMinus_nNeutralKaons, tauMinus_nPhotons;
    Float_t tauMinus_sumPhotonEn;
    bai.setBranchAddress(tauMinus_nChargedKaons, "gen_tauMinus_nChargedKaons");
    bai.setBranchAddress(tauMinus_nNeutralKaons, "gen_tauMinus_nNeutralKaons");
    bai.setBranchAddress(tauMinus_nPhotons, "gen_tauMinus_nPhotons");
    bai.setBranchAddress(tauMinus_sumPhotonEn, "gen_tauMinus_sumPhotonEn");

    Float_t kinFit_chi2;
    bai.setBranchAddress(kinFit_chi2, "kinFit_chi2");
    Int_t kinFit_status;
    bai.setBranchAddress(kinFit_status, "kinFit_status");

    Float_t zPlus, zMinus, cosThetaStar;
    bai.setBranchAddress(zPlus, Form("%s_zPlus", mode.c_str()));
    bai.setBranchAddress(zMinus, Form("%s_zMinus", mode.c_str()));
    bai.setBranchAddress(cosThetaStar, Form("%s_cosThetaStar", mode.c_str()));

    Float_t evtWeight = 1.;
    if ( branchName_evtWeight != "" && apply_evtWeight )
    {
      bai.setBranchAddress(evtWeight, branchName_evtWeight);
    }

    inputTree->SetBranchStatus("*", 0);
    for(const std::string & boundBranchName: bai.getBoundBranchNames())
    {
      inputTree->SetBranchStatus(boundBranchName.c_str(), 1);
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

      if ( !(zPlus > minVisTauZ) ) continue;
      if ( !(tauPlus_tip > minTauTIP) ) continue;
      if ( maxNumChargedKaons       != -1  && tauPlus_nChargedKaons  > maxNumChargedKaons            ) continue;
      if ( maxNumNeutralKaons       != -1  && tauPlus_nNeutralKaons  > maxNumNeutralKaons            ) continue;
      if ( maxNumPhotons            != -1  && tauPlus_nPhotons       > maxNumPhotons                 ) continue;
      if ( maxSumPhotonEn           >=  0. && tauPlus_sumPhotonEn    > maxSumPhotonEn                ) continue;
      if ( !(zMinus > minVisTauZ) ) continue;
      if ( !(tauMinus_tip > minTauTIP) ) continue;
      if ( maxNumChargedKaons       != -1  && tauMinus_nChargedKaons > maxNumChargedKaons            ) continue;
      if ( maxNumNeutralKaons       != -1  && tauMinus_nNeutralKaons > maxNumNeutralKaons            ) continue;
      if ( maxNumPhotons            != -1  && tauMinus_nPhotons      > maxNumPhotons                 ) continue;
      if ( maxSumPhotonEn           >=  0. && tauMinus_sumPhotonEn   > maxSumPhotonEn                ) continue;
      if ( mode == "startPos" || mode == "kinFit" )
      {
        if ( maxChi2                != -1  && kinFit_chi2            > maxChi2                       ) continue;
        if ( statusSelection.size() >   0  && !passesStatusSelection(kinFit_status, statusSelection) ) continue;
      }
      if ( absCosTheta_cut > 0. && std::fabs(cosThetaStar) > absCosTheta_cut )                         continue;

      spin::Data entry(hPlus_n, hPlus_r, hPlus_k, hMinus_n, hMinus_r, hMinus_k, evtWeight);

      if ( visPlus_pt  > minVisTauPt && std::fabs(visPlus_eta)  < maxAbsVisTauEta &&
           visMinus_pt > minVisTauPt && std::fabs(visMinus_eta) < maxAbsVisTauEta )
      {
        dataset_passed.push_back(entry);

        ++selectedEntries;
        selectedEntries_weighted += evtWeight;
      }
      else
      {
        dataset_failed.push_back(entry);
      }

      if ( maxEvents_beforeCuts != -1 && analyzedEntries >= maxEvents_beforeCuts ) STOP = true;
    }

    delete inputTree;
    delete inputFile;
  }

  clock.Start("spinAnalyzer");
  std::cout << "Processing Summary:\n";
  std::cout << " processedInputFiles = " << processedInputFiles << " (out of " << numInputFiles << ")\n";
  std::cout << " analyzedEntries = " << analyzedEntries << " (weighted = " << analyzedEntries_weighted << ")\n";
  std::cout << " selectedEntries = " << selectedEntries << " (weighted = " << selectedEntries_weighted << ")\n";
  
  if ( selectedEntries == 0 )
  {
    if ( analyzedEntries == 0 ) std::cerr << "WARNING: Tree '" << treeName << "' contains no events !!\n";
    else                        std::cerr << "WARNING: No events pass selection !!\n";
    std::cerr << "Please check that the tau decay mode you are analyzing was enabled in the Ntuple production.\n";
    return EXIT_SUCCESS;
  }

  spin::SpinAlgo_by_mlfit::set_dataset_norm_passed(&dataset_passed);
  spin::SpinAlgo_by_mlfit::set_dataset_norm_failed(&dataset_failed);

  spin::Measurement measurement = spinAnalyzer(dataset_passed);
  std::cout << "Matrix C:\n";
  std::cout << measurement.get_C() << "\n";
  std::cout << "+/-\n";
  std::cout << measurement.get_CErr() << "\n";
  std::cout << "Ek = " << measurement.get_Ek() << " +/- " << measurement.get_EkErr() << "\n";
  std::cout << "concurrence = " << measurement.get_concurrence() << " +/- " << measurement.get_concurrenceErr() << "\n";
  std::cout << "steerability = " << measurement.get_steerability() << " +/- " << measurement.get_steerabilityErr() << "\n";
  std::cout << "Rchsh = " << measurement.get_Rchsh() << " +/- " << measurement.get_RchshErr() << "\n";

  if ( verbosity >= 1 )
  {
    math::Vector3 Bp_exp;
    math::Vector3 Bm_exp;
    for ( size_t idxElement = 0; idxElement < 3; ++idxElement )
    {
      Bp_exp(idxElement) = 0.;
      Bm_exp(idxElement) = 0.;
    }
    math::Matrix3x3 C_exp;
    for ( int idxRow = 0; idxRow < 3; ++idxRow )
    {
      for ( int idxColumn = 0; idxColumn < 3; ++idxColumn )
      {
        C_exp(idxRow,idxColumn) = par_gen[6 + 3*idxRow + idxColumn];
      }
    }
    std::cout << "Standard Model expectation:\n";
    std::cout << "Matrix C:\n";
    std::cout << C_exp << "\n";
    double Ek_exp = comp_Ek(C_exp);
    std::cout << "Ek = " << Ek_exp << "\n";
    double concurrence_exp = comp_concurrence(Bp_exp, Bm_exp, C_exp, verbosity);
    std::cout << "concurrence = " << concurrence_exp << "\n";
    double Rchsh_exp = comp_Rchsh(C_exp);
    std::cout << "Rchsh = " << Rchsh_exp << "\n";
    double steerability_exp = comp_steerability(C_exp, 360, 360);
    std::cout << "steerability = " << steerability_exp << "\n";
  }
  clock.Show("spinAnalyzer");

  // Save the results to a JSON file
  const std::string jsonOutoutFileName = cfg_analyze.getParameter<std::string>("jsonOutoutFileName");
  std::ofstream jsonOutputFile(jsonOutoutFileName);
  if ( ! jsonOutputFile )
  {
    throw cmsException(argv[0], __LINE__) << "Could not create file: " << jsonOutoutFileName;
  }
  jsonOutputFile << dumpJSON(measurement, absCosTheta_cut) << '\n';
  jsonOutputFile.close();

  clock.Show("analyzeEntanglementNtuple");

  return EXIT_SUCCESS;
}
