
#include "DataFormats/FWLite/interface/InputSource.h"                 // fwlite::InputSource
#include "DataFormats/FWLite/interface/OutputFiles.h"                 // fwlite::OutputFiles
#include "FWCore/ParameterSet/interface/ParameterSet.h"               // edm::ParameterSet
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"   // edm::readPSetsFrom()
#include "FWCore/PluginManager/interface/PluginManager.h"             // edmplugin::PluginManager::configure()
#include "FWCore/PluginManager/interface/standard.h"                  // edmplugin::standard::config()
#include "PhysicsTools/FWLite/interface/TFileService.h"               // fwlite::TFileService

#include "TauAnalysis/Entanglement/interface/bookHistogram1d.h"       // bookHistogram1d()
#include "TauAnalysis/Entanglement/interface/bookHistogram2d.h"       // bookHistogram2d()
#include "TauAnalysis/Entanglement/interface/cmsException.h"          // cmsException
#include "TauAnalysis/Entanglement/interface/format_vT.h"             // format_vint(), vdouble, vint
#include "TauAnalysis/Entanglement/interface/passesStatusSelection.h" // passesStatusSelection()
#include "TauAnalysis/Entanglement/interface/showHistogram1d.h"       // showHistogram1d()
#include "TauAnalysis/Entanglement/interface/showHistogram2d.h"       // showHistogram2d()

#include <TBenchmark.h>                                               // TBenchmark
#include <TError.h>                                                   // gErrorAbortLevel, kError
#include <TH1.h>                                                      // TH1
#include <TH2.h>                                                      // TH2
#include <TString.h>                                                  // Form()
#include <TTree.h>                                                    // TTree

#include <assert.h>                                                   // assert()
#include <cmath>                                                      // std::fabs()
#include <cstdlib>                                                    // EXIT_SUCCESS, EXIT_FAILURE
#include <fstream>                                                    // std::ofstream
#include <iostream>                                                   // std::cout
#include <string>                                                     // std::string
#include <vector>                                                     // std::vector<>

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

  std::cout << "<makeControlPlots>:\n";

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("makeControlPlots");

//--- read python configuration parameters
  std::cout << "Reading config file " << argv[1] << "\n";
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") )
    throw cmsException("makeControlPlots", __LINE__) << "No ParameterSet 'process' found in config file !!";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameterSet("process");

  edm::ParameterSet cfg_ctrlPlots = cfg.getParameterSet("makeControlPlots");
  std::string treeName = cfg_ctrlPlots.getParameter<std::string>("treeName");
  std::cout << " treeName = " << treeName << "\n";
  std::string mode = cfg_ctrlPlots.getParameter<std::string>("mode");
  std::cout << " mode = " << mode << "\n";
  float minVisTauPt = cfg_ctrlPlots.getParameter<double>("minVisTauPt");
  std::cout << " minVisTauPt = " << minVisTauPt << "\n";
  float maxAbsVisTauEta = cfg_ctrlPlots.getParameter<double>("maxAbsVisTauEta");
  std::cout << " maxAbsVisTauEta = " << maxAbsVisTauEta << "\n";
  float minTauTIP = cfg_ctrlPlots.getParameter<double>("minTauTIP");
  std::cout << " minTauTIP = " << minTauTIP << "\n";
  int maxNumChargedKaons = cfg_ctrlPlots.getParameter<int>("maxNumChargedKaons");
  std::cout << " maxNumChargedKaons = " << maxNumChargedKaons << "\n";
  int maxNumNeutralKaons = cfg_ctrlPlots.getParameter<int>("maxNumNeutralKaons");
  std::cout << " maxNumNeutralKaons = " << maxNumNeutralKaons << "\n";
  int maxNumPhotons = cfg_ctrlPlots.getParameter<int>("maxNumPhotons");
  std::cout << " maxNumPhotons = " << maxNumPhotons << "\n";
  float maxSumPhotonEn = cfg_ctrlPlots.getParameter<double>("maxSumPhotonEn");
  std::cout << " maxSumPhotonEn = " << maxSumPhotonEn << "\n";
  float maxChi2 = cfg_ctrlPlots.getParameter<double>("maxChi2");
  std::cout << " maxChi2 = " << maxChi2 << "\n";
  vint statusSelection = cfg_ctrlPlots.getParameter<vint>("statusSelection");
  std::cout << " statusSelection = " << format_vint(statusSelection) << "\n";
  std::string branchName_evtWeight = cfg_ctrlPlots.getParameter<std::string>("branchName_evtWeight");
  std::cout << " branchName_evtWeight = " << branchName_evtWeight << "\n";

  bool apply_evtWeight = cfg_ctrlPlots.getParameter<bool>("apply_evtWeight");

  fwlite::InputSource inputFiles(cfg);
  int maxEvents = inputFiles.maxEvents();
  std::cout << " maxEvents = " << maxEvents << "\n";
  unsigned reportEvery = inputFiles.reportAfter();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().c_str());

  std::vector<std::string> inputFileNames = inputFiles.files();
  size_t numInputFiles = inputFileNames.size();
  std::cout << "Loaded " << numInputFiles << " file(s).\n";

  TH1* histogram_tauPlusPt     = bookHistogram1d(fs, "tauPlusPt",    40,  0., 200.);
  TH1* histogram_tauPlusEta    = bookHistogram1d(fs, "tauPlusEta",   50, -5.,  +5.);
  TH1* histogram_tauPlusTIP    = bookHistogram1d(fs, "tauPlusTIP",   50,  0.,   0.025);
  TH1* histogram_visPlusPt     = bookHistogram1d(fs, "visPlusPt",    40,  0., 200.);
  TH1* histogram_visPlusEta    = bookHistogram1d(fs, "visPlusEta",   50, -5.,  +5.);
  TH1* histogram_tauMinusPt    = bookHistogram1d(fs, "tauMinusPt",   40,  0., 200.);
  TH1* histogram_tauMinusEta   = bookHistogram1d(fs, "tauMinusEta",  50, -5.,  +5.);
  TH1* histogram_tauMinusTIP   = bookHistogram1d(fs, "tauMinusTIP",  50,  0.,   0.025);
  TH1* histogram_visMinusPt    = bookHistogram1d(fs, "visMinusPt",   40,  0., 200.);
  TH1* histogram_visMinusEta   = bookHistogram1d(fs, "visMinusEta",  50, -5.,  +5.);
  TH1* histogram_mTauTau       = bookHistogram1d(fs, "mTauTau",      40,  0., 200.);
  TH1* histogram_mVis          = bookHistogram1d(fs, "mVis",         40,  0., 200.);
  TH1* histogram_cosThetaStar  = bookHistogram1d(fs, "cosThetaStar", 40, -1.,  +1.);
  TH1* histogram_chi2          = bookHistogram1d(fs, "chi2",         50,  0.,  50.);

  TH1* histogram_Bp_r          = bookHistogram1d(fs, "Bp_r",         40, -1.,  +1.);
  TH1* histogram_Bp_n          = bookHistogram1d(fs, "Bp_n",         40, -1.,  +1.);
  TH1* histogram_Bp_k          = bookHistogram1d(fs, "Bp_k",         40, -1.,  +1.);
  TH1* histogram_Bm_r          = bookHistogram1d(fs, "Bm_r",         40, -1.,  +1.);
  TH1* histogram_Bm_n          = bookHistogram1d(fs, "Bm_n",         40, -1.,  +1.);
  TH1* histogram_Bm_k          = bookHistogram1d(fs, "Bm_k",         40, -1.,  +1.);
  TH1* histogram_C_rr          = bookHistogram1d(fs, "C_rr",         72, -9.,  +9.);
  TH1* histogram_C_nn          = bookHistogram1d(fs, "C_nn",         72, -9.,  +9.);
  TH1* histogram_C_kk          = bookHistogram1d(fs, "C_kk",         72, -9.,  +9.);

  TH2* histogram_zPlus_vs_zMinus = bookHistogram2d(fs, "zPlus_vs_zMinus", 20, 0., 1., 20, 0., 1.);

  int analyzedEntries = 0;
  double analyzedEntries_weighted = 0.;
  int selectedEntries = 0;
  double selectedEntries_weighted = 0.;
  int processedInputFiles = 0;
  bool STOP = false;
  double avEvtWeight = 0.;
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

    Float_t tauPlus_pt, tauPlus_eta, tauPlus_tip;
    inputTree->SetBranchAddress(Form("%s_tauPlus_pt", mode.c_str()), &tauPlus_pt);
    inputTree->SetBranchAddress(Form("%s_tauPlus_eta", mode.c_str()), &tauPlus_eta);
    inputTree->SetBranchAddress(Form("%s_tauPlus_tip", mode.c_str()), &tauPlus_tip);
    Float_t visPlus_pt, visPlus_eta;
    inputTree->SetBranchAddress(Form("%s_visPlus_pt", mode.c_str()), &visPlus_pt);
    inputTree->SetBranchAddress(Form("%s_visPlus_eta", mode.c_str()), &visPlus_eta);
    Int_t tauPlus_nChargedKaons, tauPlus_nNeutralKaons, tauPlus_nPhotons;
    Float_t tauPlus_sumPhotonEn;
    inputTree->SetBranchAddress("gen_tauPlus_nChargedKaons", &tauPlus_nChargedKaons);
    inputTree->SetBranchAddress("gen_tauPlus_nNeutralKaons", &tauPlus_nNeutralKaons);
    inputTree->SetBranchAddress("gen_tauPlus_nPhotons", &tauPlus_nPhotons);
    inputTree->SetBranchAddress("gen_tauPlus_sumPhotonEn", &tauPlus_sumPhotonEn);
    Float_t hPlus_r, hPlus_n, hPlus_k;
    inputTree->SetBranchAddress(Form("%s_hPlus_r", mode.c_str()), &hPlus_r);
    inputTree->SetBranchAddress(Form("%s_hPlus_n", mode.c_str()), &hPlus_n);
    inputTree->SetBranchAddress(Form("%s_hPlus_k", mode.c_str()), &hPlus_k);

    Float_t tauMinus_pt, tauMinus_eta, tauMinus_tip;
    inputTree->SetBranchAddress(Form("%s_tauMinus_pt", mode.c_str()), &tauMinus_pt);
    inputTree->SetBranchAddress(Form("%s_tauMinus_eta", mode.c_str()), &tauMinus_eta);
    inputTree->SetBranchAddress(Form("%s_tauMinus_tip", mode.c_str()), &tauMinus_tip);
    Float_t visMinus_pt, visMinus_eta;
    inputTree->SetBranchAddress(Form("%s_visMinus_pt", mode.c_str()), &visMinus_pt);
    inputTree->SetBranchAddress(Form("%s_visMinus_eta", mode.c_str()), &visMinus_eta);
    Int_t tauMinus_nChargedKaons, tauMinus_nNeutralKaons, tauMinus_nPhotons;
    Float_t tauMinus_sumPhotonEn;
    inputTree->SetBranchAddress("gen_tauMinus_nChargedKaons", &tauMinus_nChargedKaons);
    inputTree->SetBranchAddress("gen_tauMinus_nNeutralKaons", &tauMinus_nNeutralKaons);
    inputTree->SetBranchAddress("gen_tauMinus_nPhotons", &tauMinus_nPhotons);
    inputTree->SetBranchAddress("gen_tauMinus_sumPhotonEn", &tauMinus_sumPhotonEn);
    Float_t hMinus_r, hMinus_n, hMinus_k;
    inputTree->SetBranchAddress(Form("%s_hMinus_r", mode.c_str()), &hMinus_r);
    inputTree->SetBranchAddress(Form("%s_hMinus_n", mode.c_str()), &hMinus_n);
    inputTree->SetBranchAddress(Form("%s_hMinus_k", mode.c_str()), &hMinus_k);

    Float_t mTauTau, mVis, cosThetaStar;
    inputTree->SetBranchAddress(Form("%s_mTauTau", mode.c_str()), &mTauTau);
    inputTree->SetBranchAddress(Form("%s_mVis", mode.c_str()), &mVis);
    inputTree->SetBranchAddress(Form("%s_cosThetaStar", mode.c_str()), &cosThetaStar);

    Float_t zPlus, zMinus;
    inputTree->SetBranchAddress(Form("%s_zPlus", mode.c_str()), &zPlus);
    inputTree->SetBranchAddress(Form("%s_zMinus", mode.c_str()), &zMinus);

    Float_t kinFit_chi2;
    inputTree->SetBranchAddress("kinFit_chi2", &kinFit_chi2);
    Int_t kinFit_status;
    inputTree->SetBranchAddress("kinFit_status", &kinFit_status);

    Float_t evtWeight = 1.;
    if ( branchName_evtWeight != "" && apply_evtWeight )
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

      histogram_tauPlusPt->Fill(tauPlus_pt, evtWeight);
      histogram_tauPlusEta->Fill(tauPlus_eta, evtWeight);
      histogram_tauPlusTIP->Fill(tauPlus_tip, evtWeight);
      histogram_visPlusPt->Fill(visPlus_pt, evtWeight);
      histogram_visPlusEta->Fill(visPlus_eta, evtWeight);

      histogram_tauMinusPt->Fill(tauMinus_pt, evtWeight);
      histogram_tauMinusEta->Fill(tauMinus_eta, evtWeight);
      histogram_tauMinusTIP->Fill(tauMinus_tip, evtWeight);
      histogram_visMinusPt->Fill(visMinus_pt, evtWeight);
      histogram_visMinusEta->Fill(visMinus_eta, evtWeight);

      histogram_mTauTau->Fill(mTauTau, evtWeight);
      histogram_mVis->Fill(mVis, evtWeight);
      histogram_cosThetaStar->Fill(cosThetaStar, evtWeight);

      histogram_chi2->Fill(kinFit_chi2, evtWeight);

      histogram_Bp_r->Fill(hPlus_r, evtWeight);
      histogram_Bp_n->Fill(hPlus_n, evtWeight);
      histogram_Bp_k->Fill(hPlus_k, evtWeight);

      histogram_Bm_r->Fill(hMinus_r, evtWeight);
      histogram_Bm_n->Fill(hMinus_n, evtWeight);
      histogram_Bm_k->Fill(hMinus_k, evtWeight);

      double c = -9.;
      histogram_C_rr->Fill(c*hPlus_r*hMinus_r, evtWeight);
      histogram_C_nn->Fill(c*hPlus_n*hMinus_n, evtWeight);
      histogram_C_kk->Fill(c*hPlus_k*hMinus_k, evtWeight);

      histogram_zPlus_vs_zMinus->Fill(zMinus, zPlus, evtWeight);

      ++selectedEntries;
      selectedEntries_weighted += evtWeight;

      if ( maxEvents != -1 && analyzedEntries >= maxEvents ) STOP = true;

      avEvtWeight += std::fabs(evtWeight);
    }

    delete inputTree;
    delete inputFile;
  }

  std::cout << "Processing Summary:\n";
  std::cout << " processedInputFiles = " << processedInputFiles << " (out of " << numInputFiles << ")\n";
  std::cout << " analyzedEntries = " << analyzedEntries << " (weighted = " << analyzedEntries_weighted << ")\n";
  std::cout << " selectedEntries = " << selectedEntries << " (weighted = " << selectedEntries_weighted << ")\n";

  avEvtWeight /= selectedEntries;

  showHistogram1d(histogram_tauPlusPt,       "#tau^{+} p_{T} [GeV]",     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_tauPlusEta,      "#tau^{+} #eta",            1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());
  showHistogram1d(histogram_tauPlusTIP,      "#tau^{+} d_{IP} [cm]",     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());
  showHistogram1d(histogram_visPlusPt,       "#tau^{+}_{h} p_{T} [GeV]", 1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_visPlusEta,      "#tau^{+}_{h} #eta",        1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());

  showHistogram1d(histogram_tauMinusPt,      "#tau^{-} p_{T} [GeV]",     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_tauMinusEta,     "#tau^{-} #eta",            1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());
  showHistogram1d(histogram_tauMinusTIP,     "#tau^{-} d_{IP} [cm]",     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());
  showHistogram1d(histogram_visMinusPt,      "#tau^{-}_{h} p_{T} [GeV]", 1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_visMinusEta,     "#tau^{-}_{h} #eta",        1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());

  showHistogram1d(histogram_mTauTau,         "m_{#tau#tau} [GeV]",       1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_mVis,            "m_{vis} [GeV]",            1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_cosThetaStar,    "cos(#theta^{*})",          1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());

  showHistogram1d(histogram_Bp_r,            "Bp_r",                     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_Bp_n,            "Bp_n",                     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_Bp_k,            "Bp_k",                     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());

  showHistogram1d(histogram_Bm_r,            "Bm_r",                     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_Bm_n,            "Bm_n",                     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_Bm_k,            "Bm_k",                     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());

  showHistogram1d(histogram_C_rr,            "C_rr",                     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_C_nn,            "C_nn",                     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_C_kk,            "C_kk",                     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());

  showHistogram2d(histogram_zPlus_vs_zMinus, "z^{-}",                    1.2,                    "z^{+}",  1.3, avEvtWeight, false, outputFile.file());

  clock.Show("makeControlPlots");

  return EXIT_SUCCESS;
}
