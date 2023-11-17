
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

bool
contains(const vstring & list,
         const std::string & element)
{
  return std::find(list.begin(), list.cend(), element) != list.end();
}

bool
isSelected_decayMode(int tau_decayMode, const vint& selected_decayModes)
{
  bool retVal = false;
  for ( int selected_decayMode : selected_decayModes )
  {
    if ( tau_decayMode == selected_decayMode )
    {
      retVal = true;
      break;
    }
  }
  return retVal;
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
  vint tauPlus_decayModes = cfg_analyze.getParameter<vint>("tauPlus_decayModes");
  std::cout << " tauPlus_decayModes = " << format_vint(tauPlus_decayModes) << "\n";
  vint tauMinus_decayModes = cfg_analyze.getParameter<vint>("tauMinus_decayModes");
  std::cout << " tauMinus_decayModes = " << format_vint(tauMinus_decayModes) << "\n";
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
  bool apply_acceptanceCuts = cfg_analyze.getParameter<bool>("apply_acceptanceCuts");
  std::cout << " apply_acceptanceCuts = " << apply_acceptanceCuts << "\n";
  bool apply_startPos_isCorrectSignCut = cfg_analyze.getParameter<bool>("apply_startPos_isCorrectSignCut");
  if ( apply_startPos_isCorrectSignCut )
  {
    std::cerr << "WARNING: apply_startPos_isCorrectSignCut enabled in Configuration file - this setting is to be used for debugging purposes only !!\n";
  }
  std::string branchName_evtWeight = cfg_analyze.getParameter<std::string>("branchName_evtWeight");
  std::cout << " branchName_evtWeight = " << branchName_evtWeight << "\n";
  
  bool apply_evtWeight = cfg_analyze.getParameter<bool>("apply_evtWeight");
  const double absCosTheta_cut = cfg_analyze.getParameter<double>("absCosTheta_cut");
  vstring binned_measurements = cfg_analyze.getParameter<vstring>("binned_measurements");

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

  spin::BinnedDataset1d * binnedDataset_zPlus        = contains(binned_measurements, "zPlus")        ? new spin::BinnedDataset1d("zPlus",        20,  0.,  1.) : nullptr;
  spin::BinnedDataset1d * binnedDataset_zMinus       = contains(binned_measurements, "zMinus")       ? new spin::BinnedDataset1d("zMinus",       20,  0.,  1.) : nullptr;
  spin::BinnedDataset1d * binnedDataset_cosThetaStar = contains(binned_measurements, "cosThetaStar") ? new spin::BinnedDataset1d("cosThetaStar", 20, -1., +1.) : nullptr;

  spin::BinnedDataset2d * binnedDataset_zPlus_vs_cosThetaStar   = contains(binned_measurements, "zPlus_vs_cosThetaStar")   ? new spin::BinnedDataset2d("zPlus_vs_cosThetaStar",   10, -1., +1., 10, 0., 1.) : nullptr;
  spin::BinnedDataset2d * binnedDataset_zMinus_vs_cosThetaStar  = contains(binned_measurements, "zMinus_vs_cosThetaStar")  ? new spin::BinnedDataset2d("zMinus_vs_cosThetaStar",  10, -1., +1., 10, 0., 1.) : nullptr;
  spin::BinnedDataset2d * binnedDataset_zPlus_vs_zMinus         = contains(binned_measurements, "zPlus_vs_zMinus")         ? new spin::BinnedDataset2d("zPlus_vs_zMinus",         10,  0.,  1., 10, 0., 1.) : nullptr;
  spin::BinnedDataset2d * binnedDataset_visPlusPt_vs_visMinusPt = contains(binned_measurements, "visPlusPt_vs_visMinusPt") ? new spin::BinnedDataset2d("visPlusPt_vs_visMinusPt", 12,  0.,  6., 12, 0., 6.) : nullptr;

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
    Int_t tauPlus_decayMode;
    bai.setBranchAddress(tauPlus_decayMode, "gen_tauPlus_decayMode");
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
    Int_t tauMinus_decayMode;
    bai.setBranchAddress(tauMinus_decayMode, "gen_tauMinus_decayMode");
    Int_t tauMinus_nChargedKaons, tauMinus_nNeutralKaons, tauMinus_nPhotons;
    Float_t tauMinus_sumPhotonEn;
    bai.setBranchAddress(tauMinus_nChargedKaons, "gen_tauMinus_nChargedKaons");
    bai.setBranchAddress(tauMinus_nNeutralKaons, "gen_tauMinus_nNeutralKaons");
    bai.setBranchAddress(tauMinus_nPhotons, "gen_tauMinus_nPhotons");
    bai.setBranchAddress(tauMinus_sumPhotonEn, "gen_tauMinus_sumPhotonEn");

    Int_t startPos_isCorrectSign;
    bai.setBranchAddress(startPos_isCorrectSign, "startPos_isCorrectSign");

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

    Bool_t passesAcceptanceCuts = true;
    bai.setBranchAddress(passesAcceptanceCuts, "passesAcceptanceCuts");

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
      if ( tauPlus_decayModes.size() > 0  && !isSelected_decayMode(tauPlus_decayMode, tauPlus_decayModes)   ) continue;
      if ( maxNumChargedKaons       != -1  && tauPlus_nChargedKaons  > maxNumChargedKaons            ) continue;
      if ( maxNumNeutralKaons       != -1  && tauPlus_nNeutralKaons  > maxNumNeutralKaons            ) continue;
      if ( maxNumPhotons            != -1  && tauPlus_nPhotons       > maxNumPhotons                 ) continue;
      if ( maxSumPhotonEn           >=  0. && tauPlus_sumPhotonEn    > maxSumPhotonEn                ) continue;
      if ( !(zMinus > minVisTauZ) ) continue;
      if ( !(tauMinus_tip > minTauTIP) ) continue;
      if ( tauMinus_decayModes.size() > 0 && !isSelected_decayMode(tauMinus_decayMode, tauMinus_decayModes) ) continue;
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
      if ( apply_acceptanceCuts && !passesAcceptanceCuts )                                             continue;
      if ( apply_startPos_isCorrectSignCut && startPos_isCorrectSign != +1 )                           continue;

      spin::Data entry(hPlus_n, hPlus_r, hPlus_k, hMinus_n, hMinus_r, hMinus_k, evtWeight);

      if ( visPlus_pt  > minVisTauPt && std::fabs(visPlus_eta)  < maxAbsVisTauEta &&
           visMinus_pt > minVisTauPt && std::fabs(visMinus_eta) < maxAbsVisTauEta )
      {
        dataset_passed.push_back(entry);

        if(binnedDataset_zPlus)        { binnedDataset_zPlus->push_back(zPlus, entry); }
        if(binnedDataset_zMinus)       { binnedDataset_zMinus->push_back(zMinus, entry); }
        if(binnedDataset_cosThetaStar) { binnedDataset_cosThetaStar->push_back(cosThetaStar, entry); }

        if(binnedDataset_zPlus_vs_cosThetaStar)   { binnedDataset_zPlus_vs_cosThetaStar->push_back(cosThetaStar, zPlus, entry); }
        if(binnedDataset_zMinus_vs_cosThetaStar)  { binnedDataset_zMinus_vs_cosThetaStar->push_back(cosThetaStar, zMinus, entry); }
        if(binnedDataset_zPlus_vs_zMinus)         { binnedDataset_zPlus_vs_zMinus->push_back(zMinus, zPlus,  entry); }
        if(binnedDataset_visPlusPt_vs_visMinusPt) { binnedDataset_visPlusPt_vs_visMinusPt->push_back(visMinus_pt, visPlus_pt, entry); }

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

  // CV: compute binned measurements only if spinAnalyzer is set to 'by_summation' mode,
  //     as running the binned measurements in 'by_mlfit' mode is too time-consuming
  if(binnedDataset_zPlus)
  {
    std::cout << "Processing binned measurement as function of zPlus...\n";
    spin::BinnedMeasurement1d binnedMeasurement_zPlus = spinAnalyzer(*binnedDataset_zPlus);
    TH1* histogram_Rchsh_vs_zPlus = binnedMeasurement_zPlus.get_histogram("Rchsh");
    addToOutputFile(fs, histogram_Rchsh_vs_zPlus);
    std::cout << " Done.\n";
  }
  if(binnedDataset_zMinus)
  {
    std::cout << "Processing binned measurement as function of zMinus...\n";
    spin::BinnedMeasurement1d binnedMeasurement_zMinus = spinAnalyzer(*binnedDataset_zMinus);
    TH1* histogram_Rchsh_vs_zMinus = binnedMeasurement_zMinus.get_histogram("Rchsh");
    addToOutputFile(fs, histogram_Rchsh_vs_zMinus);
    std::cout << " Done.\n";
  }
  if(binnedDataset_cosThetaStar)
  {
    std::cout << "Processing binned measurement as function of cosThetaStar...\n";
    spin::BinnedMeasurement1d binnedMeasurement_cosThetaStar = spinAnalyzer(*binnedDataset_cosThetaStar);
    TH1* histogram_Rchsh_vs_cosThetaStar = binnedMeasurement_cosThetaStar.get_histogram("Rchsh");
    addToOutputFile(fs, histogram_Rchsh_vs_cosThetaStar);
    std::cout << " Done.\n";
  }

  if(binnedDataset_zPlus_vs_cosThetaStar)
  {
    std::cout << "Processing binned measurement as function of zPlus and cosThetaStar...\n";
    spin::BinnedMeasurement2d binnedMeasurement_zPlus_vs_cosThetaStar = spinAnalyzer(*binnedDataset_zPlus_vs_cosThetaStar);
    TH2* histogram_Rchsh_vs_zPlus_and_cosThetaStar = binnedMeasurement_zPlus_vs_cosThetaStar.get_histogram("Rchsh");
    addToOutputFile(fs, histogram_Rchsh_vs_zPlus_and_cosThetaStar);
    std::cout << " Done.\n";
  }
  if(binnedDataset_zMinus_vs_cosThetaStar)
  {
    std::cout << "Processing binned measurement as function of zMinus and cosThetaStar...\n";
    spin::BinnedMeasurement2d binnedMeasurement_zMinus_vs_cosThetaStar = spinAnalyzer(*binnedDataset_zMinus_vs_cosThetaStar);
    TH2* histogram_Rchsh_vs_zMinus_vs_cosThetaStar = binnedMeasurement_zMinus_vs_cosThetaStar.get_histogram("Rchsh");
    addToOutputFile(fs, histogram_Rchsh_vs_zMinus_vs_cosThetaStar);
    std::cout << " Done.\n";
  }
  if(binnedDataset_zPlus_vs_zMinus)
  {
    std::cout << "Processing binned measurement as function of zPlus and zMinus...\n";
    spin::BinnedMeasurement2d binnedMeasurement_zPlus_vs_zMinus = spinAnalyzer(*binnedDataset_zPlus_vs_zMinus);
    TH2* histogram_Rchsh_vs_zPlus_vs_zMinus = binnedMeasurement_zPlus_vs_zMinus.get_histogram("Rchsh");
    addToOutputFile(fs, histogram_Rchsh_vs_zPlus_vs_zMinus);
    std::cout << " Done.\n";
  }
  if(binnedDataset_visPlusPt_vs_visMinusPt)
  {
    std::cout << "Processing binned measurement as function of visPlusPt and visMinusPt...\n";
    spin::BinnedMeasurement2d binnedMeasurement_visPlusPt_vs_visMinusPt = spinAnalyzer(*binnedDataset_visPlusPt_vs_visMinusPt);
    TH2* histogram_Rchsh_vs_visPlusPt_vs_visMinusPt = binnedMeasurement_visPlusPt_vs_visMinusPt.get_histogram("Rchsh");
    addToOutputFile(fs, histogram_Rchsh_vs_visPlusPt_vs_visMinusPt);
    std::cout << " Done.\n";
  }

  clock.Show("analyzeEntanglementNtuple");

  return EXIT_SUCCESS;
}
