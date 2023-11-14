
#include "DataFormats/FWLite/interface/InputSource.h"                             // fwlite::InputSource
#include "DataFormats/FWLite/interface/OutputFiles.h"                             // fwlite::OutputFiles
#include "FWCore/ParameterSet/interface/ParameterSet.h"                           // edm::ParameterSet
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"               // edm::readPSetsFrom()
#include "FWCore/PluginManager/interface/PluginManager.h"                         // edmplugin::PluginManager::configure()
#include "FWCore/PluginManager/interface/standard.h"                              // edmplugin::standard::config()
#include "PhysicsTools/FWLite/interface/TFileService.h"                           // fwlite::TFileService

#include "TauAnalysis/Entanglement/interface/cmsException.h"                      // cmsException
#include "TauAnalysis/Entanglement/interface/comp_BandC.h"                        // comp_C()
#include "TauAnalysis/Entanglement/interface/comp_EigenVectors_and_EigenValues.h" // comp_EigenVectors_and_EigenValues()
#include "TauAnalysis/Entanglement/interface/comp_Rchsh.h"                        // comp_Rchsh()
#include "TauAnalysis/Entanglement/interface/convert_to_TMatrixD.h"               // convert_to_TMatrixD()
#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h"                 // math::Matrix3x3
#include "TauAnalysis/Entanglement/interface/printEigenVectors_and_EigenValues.h" // printEigenVectors_and_EigenValues()
#include "TauAnalysis/Entanglement/interface/BranchAddressInitializer.h"          // BranchAddressInitializer

#include <TBenchmark.h>                                                           // TBenchmark
#include <TError.h>                                                               // gErrorAbortLevel, kError
#include <TString.h>                                                              // Form()
#include <TTree.h>                                                                // TTree

#include <cmath>                                                                  // std::fabs()
#include <cstdlib>                                                                // EXIT_SUCCESS, EXIT_FAILURE
#include <iostream>                                                               // std::cout
#include <string>                                                                 // std::string
#include <vector>                                                                 // std::vector<>

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

  std::cout << "<dumpEntanglementNtuple>:\n";

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("dumpEntanglementNtuple");

//--- read python configuration parameters
  std::cout << "Reading config file " << argv[1] << "\n";
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") )
    throw cmsException("dumpEntanglementNtuple", __LINE__) << "No ParameterSet 'process' found in config file !!";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameterSet("process");

  edm::ParameterSet cfg_input = cfg.getParameterSet("fwliteInput");
  int maxEvents_beforeCuts = cfg_input.getParameter<int>("maxEvents_beforeCuts");
  int maxEvents_afterCuts = cfg_input.getParameter<int>("maxEvents_afterCuts");
  std::cout << " maxEvents: beforeCuts = " << maxEvents_beforeCuts << ", afterCuts = " << maxEvents_afterCuts << "\n";

  edm::ParameterSet cfg_dump = cfg.getParameterSet("dumpEntanglementNtuple");
  std::string treeName = cfg_dump.getParameter<std::string>("treeName");
  std::cout << " treeName = " << treeName << "\n";
  std::string mode = cfg_dump.getParameter<std::string>("mode");
  std::cout << " mode = " << mode << "\n";
  float minVisTauPt = cfg_dump.getParameter<double>("minVisTauPt");
  std::cout << " minVisTauPt = " << minVisTauPt << "\n";
  float maxAbsVisTauEta = cfg_dump.getParameter<double>("maxAbsVisTauEta");
  std::cout << " maxAbsVisTauEta = " << maxAbsVisTauEta << "\n";
  float minVisTauZ = cfg_dump.getParameter<double>("minVisTauZ");
  std::cout << " minVisTauZ = " << minVisTauZ << "\n";
  float minTauTIP = cfg_dump.getParameter<double>("minTauTIP");
  std::cout << " minTauTIP = " << minTauTIP << "\n";
  int maxNumChargedKaons = cfg_dump.getParameter<int>("maxNumChargedKaons");
  std::cout << " maxNumChargedKaons = " << maxNumChargedKaons << "\n";
  int maxNumNeutralKaons = cfg_dump.getParameter<int>("maxNumNeutralKaons");
  std::cout << " maxNumNeutralKaons = " << maxNumNeutralKaons << "\n";
  int maxNumPhotons = cfg_dump.getParameter<int>("maxNumPhotons");
  std::cout << " maxNumPhotons = " << maxNumPhotons << "\n";
  float maxSumPhotonEn = cfg_dump.getParameter<double>("maxSumPhotonEn");
  std::cout << " maxSumPhotonEn = " << maxSumPhotonEn << "\n";
  std::string branchName_evtWeight = cfg_dump.getParameter<std::string>("branchName_evtWeight");
  std::cout << " branchName_evtWeight = " << branchName_evtWeight << "\n";

  bool apply_evtWeight = cfg_dump.getParameter<bool>("apply_evtWeight");

  //int verbosity = cfg_dump.getUntrackedParameter<int>("verbosity");

  fwlite::InputSource inputFiles(cfg);
  unsigned reportEvery = inputFiles.reportAfter();

  //fwlite::OutputFiles outputFile(cfg);
  //fwlite::TFileService fs = fwlite::TFileService(outputFile.file().c_str());

  std::vector<std::string> inputFileNames = inputFiles.files();
  size_t numInputFiles = inputFileNames.size();
  std::cout << "Loaded " << numInputFiles << " file(s).\n";

  int analyzedEntries = 0;
  double analyzedEntries_weighted = 0.;
  int selectedEntries = 0;
  double selectedEntries_weighted = 0.;
  math::Matrix3x3 C_sum;
  int processedInputFiles = 0;
  bool isFirst = true;
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

    UInt_t run, lumi;
    ULong64_t event;
    bai.setBranchAddress(run, "run");
    bai.setBranchAddress(lumi, "lumi");
    bai.setBranchAddress(event, "event");

    ULong64_t entry;
    bai.setBranchAddress(entry, "entry");

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

    Float_t zPlus, zMinus, cosThetaStar;
    bai.setBranchAddress(zPlus, Form("%s_zPlus", mode.c_str()));
    bai.setBranchAddress(zMinus, Form("%s_zMinus", mode.c_str()));
    bai.setBranchAddress(cosThetaStar, Form("%s_cosThetaStar", mode.c_str()));

    Float_t evtWeight = 1.;
    if ( branchName_evtWeight != "" && apply_evtWeight )
    {
      bai.setBranchAddress(evtWeight, branchName_evtWeight.c_str());
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

      if ( !(visPlus_pt  > minVisTauPt && std::fabs(visPlus_eta)  < maxAbsVisTauEta) ) continue;
      if ( !(zPlus > minVisTauZ) ) continue;
      if ( !(tauPlus_tip > minTauTIP) ) continue;
      if ( maxNumChargedKaons       != -1  && tauPlus_nChargedKaons  > maxNumChargedKaons            ) continue;
      if ( maxNumNeutralKaons       != -1  && tauPlus_nNeutralKaons  > maxNumNeutralKaons            ) continue;
      if ( maxNumPhotons            != -1  && tauPlus_nPhotons       > maxNumPhotons                 ) continue;
      if ( maxSumPhotonEn           >=  0. && tauPlus_sumPhotonEn    > maxSumPhotonEn                ) continue;
      if ( !(visMinus_pt > minVisTauPt && std::fabs(visMinus_eta) < maxAbsVisTauEta) ) continue;
      if ( !(zPlus > minVisTauZ) ) continue;
      if ( !(tauMinus_tip > minTauTIP) ) continue;
      if ( maxNumChargedKaons       != -1  && tauMinus_nChargedKaons > maxNumChargedKaons            ) continue;
      if ( maxNumNeutralKaons       != -1  && tauMinus_nNeutralKaons > maxNumNeutralKaons            ) continue;
      if ( maxNumPhotons            != -1  && tauMinus_nPhotons      > maxNumPhotons                 ) continue;
      if ( maxSumPhotonEn           >=  0. && tauMinus_sumPhotonEn   > maxSumPhotonEn                ) continue;

      // CV: compute spin correlation matrix C according to Eq. (25)
      //     in the paper arXiv:2211.10513
      math::Matrix3x3 C = comp_C(hPlus_n, hPlus_r, hPlus_k, hMinus_n, hMinus_r, hMinus_k);

      C_sum += evtWeight*C;

      if ( isFirst )
      {
        std::cout << "--------------------------------------------------------------------------------\n";
        isFirst = false;
      }
      std::cout << "event #" << entry << " (" << run << ":" << lumi << ":" << event << "):\n";
      std::cout << "C:\n";
      std::cout << C << "\n";
      std::cout << "cos(theta*) = " << cosThetaStar << "\n";
      std::cout << "--------------------------------------------------------------------------------\n";

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

  C_sum *= 1./selectedEntries_weighted;

  std::cout << "<C>:\n";
  std::cout << C_sum << "\n";

  math::Matrix3x3 CT_sum = ROOT::Math::Transpose(C_sum);
  std::vector<std::pair<TVectorD, double>> EigenVectors_and_EigenValues = comp_EigenVectors_and_EigenValues(convert_to_TMatrixD(CT_sum*C_sum));
  printEigenVectors_and_EigenValues(EigenVectors_and_EigenValues);

  double Rchsh = comp_Rchsh(C_sum);
  std::cout << "Rchsh = " << Rchsh << "\n";

  clock.Show("dumpEntanglementNtuple");

  return EXIT_SUCCESS;
}
