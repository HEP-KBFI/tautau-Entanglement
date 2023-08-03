
#include "DataFormats/Candidate/interface/Candidate.h"                    // reco::Candidate::LorentzVector
#include "DataFormats/FWLite/interface/InputSource.h"                     // fwlite::InputSource
#include "DataFormats/FWLite/interface/OutputFiles.h"                     // fwlite::OutputFiles
#include "FWCore/ParameterSet/interface/ParameterSet.h"                   // edm::ParameterSet
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"       // edm::readPSetsFrom()
#include "FWCore/PluginManager/interface/PluginManager.h"                 // edmplugin::PluginManager::configure()
#include "FWCore/PluginManager/interface/standard.h"                      // edmplugin::standard::config()
#include "PhysicsTools/FWLite/interface/TFileService.h"                   // fwlite::TFileService

#include "TauAnalysis/Entanglement/interface/bookHistogram1d.h"           // bookHistogram1d()
#include "TauAnalysis/Entanglement/interface/cmsException.h"              // cmsException
#include "TauAnalysis/Entanglement/interface/get_localCoordinateSystem.h" // get_localCoordinateSystem()
#include "TauAnalysis/Entanglement/interface/showHistogram1d.h"           // showHistogram1d()
#include "TauAnalysis/Entanglement/interface/showHistogram1d.h"           // showHistogram1d()
#include "TauAnalysis/Entanglement/interface/square.h"                    // square()

#include <TBenchmark.h>                                                   // TBenchmark
#include <TError.h>                                                       // gErrorAbortLevel, kError
#include <TH1.h>                                                          // TH1D
#include <TString.h>                                                      // Form()
#include <TTree.h>                                                        // TTree

#include <assert.h>                                                       // assert()
#include <cmath>                                                          // std::acos(), std::cos(), std::fabs(), std::sin(), std::sinh(), std::sqrt()
#include <cstdlib>                                                        // EXIT_SUCCESS, EXIT_FAILURE
#include <fstream>                                                        // std::ofstream
#include <iostream>                                                       // std::cout
#include <string>                                                         // std::string
#include <vector>                                                         // std::vector<>

reco::Candidate::LorentzVector
build_p4(double pt, double eta, double phi, double mass)
{
  // CV: formulas for (px,py,pz) as function of (pt,eta,phi) taken from 
  //       https://en.wikipedia.org/wiki/Pseudorapidity
  double px = pt*std::cos(phi);
  double py = pt*std::sin(phi);
  double pz = pt*std::sinh(eta);
  double energy = std::sqrt(square(px) + square(py) + square(pz) + square(mass));
  reco::Candidate::LorentzVector p4(px, py, pz, energy);
  return p4;
}

reco::Candidate::Point
build_point(double x, double y, double z)
{
  reco::Candidate::Point point(x, y, z);
  return point;
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

  std::cout << "<makeResolutionPlots>:\n";

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("makeResolutionPlots");

//--- read python configuration parameters
  std::cout << "Reading config file " << argv[1] << "\n";
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") )
    throw cmsException("makeResolutionPlots", __LINE__) << "No ParameterSet 'process' found in config file !!";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameterSet("process");

  edm::ParameterSet cfg_resPlots = cfg.getParameterSet("makeResolutionPlots");
  std::string treeName = cfg_resPlots.getParameter<std::string>("treeName");
  std::cout << " treeName = " << treeName << "\n";
  std::string mode = cfg_resPlots.getParameter<std::string>("mode");
  std::cout << " mode = " << mode << "\n";
  float minVisTauPt = cfg_resPlots.getParameter<double>("minVisTauPt");
  std::cout << " minVisTauPt = " << minVisTauPt << "\n";
  float maxAbsVisTauEta = cfg_resPlots.getParameter<double>("maxAbsVisTauEta");
  std::cout << " maxAbsVisTauEta = " << maxAbsVisTauEta << "\n";
  int maxNumChargedKaons = cfg_resPlots.getParameter<int>("maxNumChargedKaons");
  std::cout << " maxNumChargedKaons = " << maxNumChargedKaons << "\n";
  int maxNumNeutralKaons = cfg_resPlots.getParameter<int>("maxNumNeutralKaons");
  std::cout << " maxNumNeutralKaons = " << maxNumNeutralKaons << "\n";
  int maxNumPhotons = cfg_resPlots.getParameter<int>("maxNumPhotons");
  std::cout << " maxNumPhotons = " << maxNumPhotons << "\n";
  float maxSumPhotonEn = cfg_resPlots.getParameter<double>("maxSumPhotonEn");
  std::cout << " maxSumPhotonEn = " << maxSumPhotonEn << "\n";
  std::string branchName_evtWeight = cfg_resPlots.getParameter<std::string>("branchName_evtWeight");
  std::cout << " branchName_evtWeight = " << branchName_evtWeight << "\n";
  //bool isDEBUG = cfg_analyze.getParameter<bool>("isDEBUG");

  fwlite::InputSource inputFiles(cfg);
  int maxEvents = inputFiles.maxEvents();
  std::cout << " maxEvents = " << maxEvents << "\n";
  unsigned reportEvery = inputFiles.reportAfter();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().c_str());

  std::vector<std::string> inputFileNames = inputFiles.files();
  size_t numInputFiles = inputFileNames.size();
  std::cout << "Loaded " << numInputFiles << " file(s).\n";

  TH1* histogram_res_pv_x         = bookHistogram1d(fs, "res_pv_x",         100,  -0.01, +0.01);
  TH1* histogram_res_pv_y         = bookHistogram1d(fs, "res_pv_y",         100,  -0.01, +0.01);
  TH1* histogram_res_pv_z         = bookHistogram1d(fs, "res_pv_z",         100,  -0.01, +0.01);

  TH1* histogram_res_tauPlus_pt   = bookHistogram1d(fs, "res_tauPlus_pt",   100, -50.,  +50.);
  TH1* histogram_res_tauPlus_eta  = bookHistogram1d(fs, "res_tauPlus_eta",  100,  -0.5,  +0.5);
  TH1* histogram_res_tauPlus_phi  = bookHistogram1d(fs, "res_tauPlus_phi",  100,  -0.5,  +0.5);
  TH1* histogram_res_svPlus_r     = bookHistogram1d(fs, "res_svPlus_r",     100,  -0.01, +0.01);
  TH1* histogram_res_svPlus_n     = bookHistogram1d(fs, "res_svPlus_n",     100,  -0.01, +0.01);
  TH1* histogram_res_svPlus_k     = bookHistogram1d(fs, "res_svPlus_k",     100,  -0.25, +0.25);
  TH1* histogram_res_visPlus_pt   = bookHistogram1d(fs, "res_visPlus_pt",   100, -50.,  +50.);
  TH1* histogram_res_visPlus_eta  = bookHistogram1d(fs, "res_visPlus_eta",  100,  -0.5,  +0.5);
  TH1* histogram_res_visPlus_phi  = bookHistogram1d(fs, "res_visPlus_phi",  100,  -0.5,  +0.5);
  TH1* histogram_res_nuPlus_pt    = bookHistogram1d(fs, "res_nuPlus_pt",    100, -50. , +50.);
  TH1* histogram_res_nuPlus_eta   = bookHistogram1d(fs, "res_nuPlus_eta",   100,  -0.5,  +0.5);
  TH1* histogram_res_nuPlus_phi   = bookHistogram1d(fs, "res_nuPlus_phi",   100,  -0.5,  +0.5);
  TH1* histogram_res_nuPlus_px    = bookHistogram1d(fs, "res_nuPlus_px",    100, -50.,  +50.);
  TH1* histogram_res_nuPlus_py    = bookHistogram1d(fs, "res_nuPlus_py",    100, -50.,  +50.);
  TH1* histogram_res_nuPlus_pz    = bookHistogram1d(fs, "res_nuPlus_pz",    100, -50.,  +50.);
  TH1* histogram_res_hPlus_r      = bookHistogram1d(fs, "res_hPlus_r",      100,  -0.5,  +0.5);
  TH1* histogram_res_hPlus_n      = bookHistogram1d(fs, "res_hPlus_n",      100,  -0.5,  +0.5);
  TH1* histogram_res_hPlus_k      = bookHistogram1d(fs, "res_hPlus_k",      100,  -0.5,  +0.5);
  TH1* histogram_res_hPlus_angle  = bookHistogram1d(fs, "res_hPlus_angle",  100,   0.,    0.5);

  TH1* histogram_res_tauMinus_pt  = bookHistogram1d(fs, "res_tauMinus_pt",  100, -50.,  +50.);
  TH1* histogram_res_tauMinus_eta = bookHistogram1d(fs, "res_tauMinus_eta", 100,  -0.5,  +0.5);
  TH1* histogram_res_tauMinus_phi = bookHistogram1d(fs, "res_tauMinus_phi", 100,  -0.5,  +0.5);
  TH1* histogram_res_svMinus_r    = bookHistogram1d(fs, "res_svMinus_r",    100,  -0.01, +0.01);
  TH1* histogram_res_svMinus_n    = bookHistogram1d(fs, "res_svMinus_n",    100,  -0.01, +0.01);
  TH1* histogram_res_svMinus_k    = bookHistogram1d(fs, "res_svMinus_k",    100,  -0.25, +0.25);
  TH1* histogram_res_nuMinus_pt   = bookHistogram1d(fs, "res_nuMinus_Pt",   100, -50.,  +50.);
  TH1* histogram_res_nuMinus_eta  = bookHistogram1d(fs, "res_nuMinus_Eta",  100,  -0.5,  +0.5);
  TH1* histogram_res_nuMinus_phi  = bookHistogram1d(fs, "res_nuMinus_Phi",  100,  -0.5,  +0.5);
  TH1* histogram_res_nuMinus_px   = bookHistogram1d(fs, "res_nuMinus_Px",   100, -50.,  +50.);
  TH1* histogram_res_nuMinus_py   = bookHistogram1d(fs, "res_nuMinus_Py",   100, -50.,  +50.);
  TH1* histogram_res_nuMinus_pz   = bookHistogram1d(fs, "res_nuMinus_Pz",   100, -50.,  +50.);
  TH1* histogram_res_visMinus_pt  = bookHistogram1d(fs, "res_visMinus_pt",  100, -50.,  +50.);
  TH1* histogram_res_visMinus_eta = bookHistogram1d(fs, "res_visMinus_eta", 100,  -0.5,  +0.5);
  TH1* histogram_res_visMinus_phi = bookHistogram1d(fs, "res_visMinus_phi", 100,  -0.5,  +0.5);
  TH1* histogram_res_hMinus_r     = bookHistogram1d(fs, "res_hMinus_r",     100,  -0.5,  +0.5);
  TH1* histogram_res_hMinus_n     = bookHistogram1d(fs, "res_hMinus_n",     100,  -0.5,  +0.5);
  TH1* histogram_res_hMinus_k     = bookHistogram1d(fs, "res_hMinus_k",     100,  -0.5,  +0.5);
  TH1* histogram_res_hMinus_angle = bookHistogram1d(fs, "res_hMinus_angle", 100,   0.,    0.5);

  TH1* histogram_res_higgs_pt     = bookHistogram1d(fs, "res_higgs_pt",     100, -50.,  +50.);
  TH1* histogram_res_higgs_eta    = bookHistogram1d(fs, "res_higgs_eta",    100,  -0.5,  +0.5);
  TH1* histogram_res_higgs_phi    = bookHistogram1d(fs, "res_higgs_phi",    100,  -0.5,  +0.5);
  TH1* histogram_res_higgs_px     = bookHistogram1d(fs, "res_higgs_px",     100, -50.,  +50.);
  TH1* histogram_res_higgs_py     = bookHistogram1d(fs, "res_higgs_py",     100, -50.,  +50.);
  TH1* histogram_res_higgs_pz     = bookHistogram1d(fs, "res_higgs_pz",     100, -50.,  +50.);

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

    Float_t gen_pv_x, gen_pv_y, gen_pv_z;
    inputTree->SetBranchAddress("gen_pv_x", &gen_pv_x);
    inputTree->SetBranchAddress("gen_pv_y", &gen_pv_y);
    inputTree->SetBranchAddress("gen_pv_z", &gen_pv_z);

    Float_t gen_tauPlus_pt, gen_tauPlus_eta, gen_tauPlus_phi, gen_tauPlus_mass;
    inputTree->SetBranchAddress("gen_tauPlus_pt", &gen_tauPlus_pt);
    inputTree->SetBranchAddress("gen_tauPlus_eta", &gen_tauPlus_eta);
    inputTree->SetBranchAddress("gen_tauPlus_phi", &gen_tauPlus_phi);
    inputTree->SetBranchAddress("gen_tauPlus_mass", &gen_tauPlus_mass);
    Int_t gen_tauPlus_nChargedKaons, gen_tauPlus_nNeutralKaons, gen_tauPlus_nPhotons;
    Float_t gen_tauPlus_sumPhotonEn;
    inputTree->SetBranchAddress("gen_tauPlus_nChargedKaons", &gen_tauPlus_nChargedKaons);
    inputTree->SetBranchAddress("gen_tauPlus_nNeutralKaons", &gen_tauPlus_nNeutralKaons);
    inputTree->SetBranchAddress("gen_tauPlus_nPhotons", &gen_tauPlus_nPhotons);
    inputTree->SetBranchAddress("gen_tauPlus_sumPhotonEn", &gen_tauPlus_sumPhotonEn);
    Float_t gen_svPlus_x, gen_svPlus_y, gen_svPlus_z;
    inputTree->SetBranchAddress("gen_svPlus_x", &gen_svPlus_x);
    inputTree->SetBranchAddress("gen_svPlus_y", &gen_svPlus_y);
    inputTree->SetBranchAddress("gen_svPlus_z", &gen_svPlus_z);
    Float_t gen_visPlus_pt, gen_visPlus_eta, gen_visPlus_phi, gen_visPlus_mass;
    inputTree->SetBranchAddress("gen_visPlus_pt", &gen_visPlus_pt);
    inputTree->SetBranchAddress("gen_visPlus_eta", &gen_visPlus_eta);
    inputTree->SetBranchAddress("gen_visPlus_phi", &gen_visPlus_phi);
    inputTree->SetBranchAddress("gen_visPlus_mass", &gen_visPlus_mass);
    Float_t gen_nuPlus_pt, gen_nuPlus_eta, gen_nuPlus_phi, gen_nuPlus_mass;
    inputTree->SetBranchAddress("gen_nuPlus_pt", &gen_nuPlus_pt);
    inputTree->SetBranchAddress("gen_nuPlus_eta", &gen_nuPlus_eta);
    inputTree->SetBranchAddress("gen_nuPlus_phi", &gen_nuPlus_phi);
    inputTree->SetBranchAddress("gen_nuPlus_mass", &gen_nuPlus_mass);
    Float_t gen_hPlus_r, gen_hPlus_n, gen_hPlus_k;
    inputTree->SetBranchAddress("gen_hPlus_r", &gen_hPlus_r);
    inputTree->SetBranchAddress("gen_hPlus_n", &gen_hPlus_n);
    inputTree->SetBranchAddress("gen_hPlus_k", &gen_hPlus_k);

    Float_t gen_tauMinus_pt, gen_tauMinus_eta, gen_tauMinus_phi, gen_tauMinus_mass;
    inputTree->SetBranchAddress("gen_tauMinus_pt", &gen_tauMinus_pt);
    inputTree->SetBranchAddress("gen_tauMinus_eta", &gen_tauMinus_eta);
    inputTree->SetBranchAddress("gen_tauMinus_phi", &gen_tauMinus_phi);
    inputTree->SetBranchAddress("gen_tauMinus_mass", &gen_tauMinus_mass);
    Int_t gen_tauMinus_nChargedKaons, gen_tauMinus_nNeutralKaons, gen_tauMinus_nPhotons;
    Float_t gen_tauMinus_sumPhotonEn;
    inputTree->SetBranchAddress("gen_tauMinus_nChargedKaons", &gen_tauMinus_nChargedKaons);
    inputTree->SetBranchAddress("gen_tauMinus_nNeutralKaons", &gen_tauMinus_nNeutralKaons);
    inputTree->SetBranchAddress("gen_tauMinus_nPhotons", &gen_tauMinus_nPhotons);
    inputTree->SetBranchAddress("gen_tauMinus_sumPhotonEn", &gen_tauMinus_sumPhotonEn);
    Float_t gen_svMinus_x, gen_svMinus_y, gen_svMinus_z;
    inputTree->SetBranchAddress("gen_svMinus_x", &gen_svMinus_x);
    inputTree->SetBranchAddress("gen_svMinus_y", &gen_svMinus_y);
    inputTree->SetBranchAddress("gen_svMinus_z", &gen_svMinus_z);
    Float_t gen_visMinus_pt, gen_visMinus_eta, gen_visMinus_phi, gen_visMinus_mass;
    inputTree->SetBranchAddress("gen_visMinus_pt", &gen_visMinus_pt);
    inputTree->SetBranchAddress("gen_visMinus_eta", &gen_visMinus_eta);
    inputTree->SetBranchAddress("gen_visMinus_phi", &gen_visMinus_phi);
    inputTree->SetBranchAddress("gen_visMinus_mass", &gen_visMinus_mass);
    Float_t gen_nuMinus_pt, gen_nuMinus_eta, gen_nuMinus_phi, gen_nuMinus_mass;
    inputTree->SetBranchAddress("gen_nuMinus_pt", &gen_nuMinus_pt);
    inputTree->SetBranchAddress("gen_nuMinus_eta", &gen_nuMinus_eta);
    inputTree->SetBranchAddress("gen_nuMinus_phi", &gen_nuMinus_phi);
    inputTree->SetBranchAddress("gen_nuMinus_mass", &gen_nuMinus_mass);
    Float_t gen_hMinus_r, gen_hMinus_n, gen_hMinus_k;
    inputTree->SetBranchAddress("gen_hMinus_r", &gen_hMinus_r);
    inputTree->SetBranchAddress("gen_hMinus_n", &gen_hMinus_n);
    inputTree->SetBranchAddress("gen_hMinus_k", &gen_hMinus_k);

    Float_t rec_pv_x, rec_pv_y, rec_pv_z;
    inputTree->SetBranchAddress(Form("%s_pv_x", mode.c_str()), &rec_pv_x);
    inputTree->SetBranchAddress(Form("%s_pv_y", mode.c_str()), &rec_pv_y);
    inputTree->SetBranchAddress(Form("%s_pv_z", mode.c_str()), &rec_pv_z);

    Float_t rec_tauPlus_pt, rec_tauPlus_eta, rec_tauPlus_phi, rec_tauPlus_mass;
    inputTree->SetBranchAddress(Form("%s_tauPlus_pt", mode.c_str()), &rec_tauPlus_pt);
    inputTree->SetBranchAddress(Form("%s_tauPlus_eta", mode.c_str()), &rec_tauPlus_eta);
    inputTree->SetBranchAddress(Form("%s_tauPlus_phi", mode.c_str()), &rec_tauPlus_phi);
    inputTree->SetBranchAddress(Form("%s_tauPlus_mass", mode.c_str()), &rec_tauPlus_mass);
    Float_t rec_svPlus_x, rec_svPlus_y, rec_svPlus_z;
    inputTree->SetBranchAddress(Form("%s_svPlus_x", mode.c_str()), &rec_svPlus_x);
    inputTree->SetBranchAddress(Form("%s_svPlus_y", mode.c_str()), &rec_svPlus_y);
    inputTree->SetBranchAddress(Form("%s_svPlus_z", mode.c_str()), &rec_svPlus_z);
    Float_t rec_visPlus_pt, rec_visPlus_eta, rec_visPlus_phi, rec_visPlus_mass;
    inputTree->SetBranchAddress(Form("%s_visPlus_pt", mode.c_str()), &rec_visPlus_pt);
    inputTree->SetBranchAddress(Form("%s_visPlus_eta", mode.c_str()), &rec_visPlus_eta);
    inputTree->SetBranchAddress(Form("%s_visPlus_phi", mode.c_str()), &rec_visPlus_phi);
    inputTree->SetBranchAddress(Form("%s_visPlus_mass", mode.c_str()), &rec_visPlus_mass);
    Float_t rec_nuPlus_pt, rec_nuPlus_eta, rec_nuPlus_phi, rec_nuPlus_mass;
    inputTree->SetBranchAddress(Form("%s_nuPlus_pt", mode.c_str()), &rec_nuPlus_pt);
    inputTree->SetBranchAddress(Form("%s_nuPlus_eta", mode.c_str()), &rec_nuPlus_eta);
    inputTree->SetBranchAddress(Form("%s_nuPlus_phi", mode.c_str()), &rec_nuPlus_phi);
    inputTree->SetBranchAddress(Form("%s_nuPlus_mass", mode.c_str()), &rec_nuPlus_mass);
    Float_t rec_hPlus_r, rec_hPlus_n, rec_hPlus_k;
    inputTree->SetBranchAddress(Form("%s_hPlus_r", mode.c_str()), &rec_hPlus_r);
    inputTree->SetBranchAddress(Form("%s_hPlus_n", mode.c_str()), &rec_hPlus_n);
    inputTree->SetBranchAddress(Form("%s_hPlus_k", mode.c_str()), &rec_hPlus_k);

    Float_t rec_tauMinus_pt, rec_tauMinus_eta, rec_tauMinus_phi, rec_tauMinus_mass;
    inputTree->SetBranchAddress(Form("%s_tauMinus_pt", mode.c_str()), &rec_tauMinus_pt);
    inputTree->SetBranchAddress(Form("%s_tauMinus_eta", mode.c_str()), &rec_tauMinus_eta);
    inputTree->SetBranchAddress(Form("%s_tauMinus_phi", mode.c_str()), &rec_tauMinus_phi);
    inputTree->SetBranchAddress(Form("%s_tauMinus_mass", mode.c_str()), &rec_tauMinus_mass);
    Float_t rec_svMinus_x, rec_svMinus_y, rec_svMinus_z;
    inputTree->SetBranchAddress(Form("%s_svMinus_x", mode.c_str()), &rec_svMinus_x);
    inputTree->SetBranchAddress(Form("%s_svMinus_y", mode.c_str()), &rec_svMinus_y);
    inputTree->SetBranchAddress(Form("%s_svMinus_z", mode.c_str()), &rec_svMinus_z);
    Float_t rec_visMinus_pt, rec_visMinus_eta, rec_visMinus_phi, rec_visMinus_mass;
    inputTree->SetBranchAddress(Form("%s_visMinus_pt", mode.c_str()), &rec_visMinus_pt);
    inputTree->SetBranchAddress(Form("%s_visMinus_eta", mode.c_str()), &rec_visMinus_eta);
    inputTree->SetBranchAddress(Form("%s_visMinus_phi", mode.c_str()), &rec_visMinus_phi);
    inputTree->SetBranchAddress(Form("%s_visMinus_mass", mode.c_str()), &rec_visMinus_mass);
    Float_t rec_nuMinus_pt, rec_nuMinus_eta, rec_nuMinus_phi, rec_nuMinus_mass;
    inputTree->SetBranchAddress(Form("%s_nuMinus_pt", mode.c_str()), &rec_nuMinus_pt);
    inputTree->SetBranchAddress(Form("%s_nuMinus_eta", mode.c_str()), &rec_nuMinus_eta);
    inputTree->SetBranchAddress(Form("%s_nuMinus_phi", mode.c_str()), &rec_nuMinus_phi);
    inputTree->SetBranchAddress(Form("%s_nuMinus_mass", mode.c_str()), &rec_nuMinus_mass);
    Float_t rec_hMinus_r, rec_hMinus_n, rec_hMinus_k;
    inputTree->SetBranchAddress(Form("%s_hMinus_r", mode.c_str()), &rec_hMinus_r);
    inputTree->SetBranchAddress(Form("%s_hMinus_n", mode.c_str()), &rec_hMinus_n);
    inputTree->SetBranchAddress(Form("%s_hMinus_k", mode.c_str()), &rec_hMinus_k);

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

      if ( !(gen_visPlus_pt  > minVisTauPt && std::fabs(gen_visPlus_eta)  < maxAbsVisTauEta) ) continue;
      if ( maxNumChargedKaons != -1  && gen_tauPlus_nChargedKaons  > maxNumChargedKaons ) continue;
      if ( maxNumNeutralKaons != -1  && gen_tauPlus_nNeutralKaons  > maxNumNeutralKaons ) continue;
      if ( maxNumPhotons      != -1  && gen_tauPlus_nPhotons       > maxNumPhotons      ) continue;
      if ( maxSumPhotonEn     >=  0. && gen_tauPlus_sumPhotonEn    > maxSumPhotonEn     ) continue;
      if ( !(gen_visMinus_pt > minVisTauPt && std::fabs(gen_visMinus_eta) < maxAbsVisTauEta) ) continue;
      if ( maxNumChargedKaons != -1  && gen_tauMinus_nChargedKaons > maxNumChargedKaons ) continue;
      if ( maxNumNeutralKaons != -1  && gen_tauMinus_nNeutralKaons > maxNumNeutralKaons ) continue;
      if ( maxNumPhotons      != -1  && gen_tauMinus_nPhotons      > maxNumPhotons      ) continue;
      if ( maxSumPhotonEn     >=  0. && gen_tauMinus_sumPhotonEn   > maxSumPhotonEn     ) continue;

      reco::Candidate::LorentzVector gen_tauPlus_p4  = build_p4(gen_tauPlus_pt,  gen_tauPlus_eta,  gen_tauPlus_phi,  gen_tauPlus_mass);
      reco::Candidate::Point         gen_svPlus      = build_point(gen_svPlus_x, gen_svPlus_y, gen_svPlus_z);
      reco::Candidate::Vector tauPlus_r, tauPlus_n, tauPlus_k;
      get_localCoordinateSystem(gen_tauPlus_p4, nullptr, nullptr, kBeam, tauPlus_r, tauPlus_n, tauPlus_k);
      double gen_svPlus_r = gen_svPlus.x()*tauPlus_r.x() + gen_svPlus.y()*tauPlus_r.y() + gen_svPlus.z()*tauPlus_r.z();
      double gen_svPlus_n = gen_svPlus.x()*tauPlus_n.x() + gen_svPlus.y()*tauPlus_n.y() + gen_svPlus.z()*tauPlus_n.z();
      double gen_svPlus_k = gen_svPlus.x()*tauPlus_k.x() + gen_svPlus.y()*tauPlus_k.y() + gen_svPlus.z()*tauPlus_k.z();
      reco::Candidate::LorentzVector gen_nuPlus_p4   = build_p4(gen_nuPlus_pt,   gen_nuPlus_eta,   gen_nuPlus_phi,   gen_nuPlus_mass);
      reco::Candidate::LorentzVector gen_tauMinus_p4 = build_p4(gen_tauMinus_pt, gen_tauMinus_eta, gen_tauMinus_phi, gen_tauMinus_mass);
      reco::Candidate::Point         gen_svMinus     = build_point(gen_svMinus_x, gen_svMinus_y, gen_svMinus_z); 
      reco::Candidate::Vector tauMinus_r, tauMinus_n, tauMinus_k;
      get_localCoordinateSystem(gen_tauMinus_p4, nullptr, nullptr, kBeam, tauMinus_r, tauMinus_n, tauMinus_k);
      double gen_svMinus_r = gen_svMinus.x()*tauMinus_r.x() + gen_svMinus.y()*tauMinus_r.y() + gen_svMinus.z()*tauMinus_r.z();
      double gen_svMinus_n = gen_svMinus.x()*tauMinus_n.x() + gen_svMinus.y()*tauMinus_n.y() + gen_svMinus.z()*tauMinus_n.z();
      double gen_svMinus_k = gen_svMinus.x()*tauMinus_k.x() + gen_svMinus.y()*tauMinus_k.y() + gen_svMinus.z()*tauMinus_k.z();
      reco::Candidate::LorentzVector gen_nuMinus_p4  = build_p4(gen_nuMinus_pt,  gen_nuMinus_eta,  gen_nuMinus_phi,  gen_nuMinus_mass);
      reco::Candidate::LorentzVector gen_higgs_p4    = gen_tauPlus_p4 + gen_tauMinus_p4;

      reco::Candidate::LorentzVector rec_tauPlus_p4  = build_p4(rec_tauPlus_pt,  rec_tauPlus_eta,  rec_tauPlus_phi,  rec_tauPlus_mass);
      reco::Candidate::Point         rec_svPlus      = build_point(rec_svPlus_x, rec_svPlus_y, rec_svPlus_z);
      double rec_svPlus_r = rec_svPlus.x()*tauPlus_r.x() + rec_svPlus.y()*tauPlus_r.y() + rec_svPlus.z()*tauPlus_r.z();
      double rec_svPlus_n = rec_svPlus.x()*tauPlus_n.x() + rec_svPlus.y()*tauPlus_n.y() + rec_svPlus.z()*tauPlus_n.z();
      double rec_svPlus_k = rec_svPlus.x()*tauPlus_k.x() + rec_svPlus.y()*tauPlus_k.y() + rec_svPlus.z()*tauPlus_k.z();
      reco::Candidate::LorentzVector rec_nuPlus_p4   = build_p4(rec_nuPlus_pt,   rec_nuPlus_eta,   rec_nuPlus_phi,   rec_nuPlus_mass);
      double hPlus_angle = std::acos(gen_hPlus_r*rec_hPlus_r + gen_hPlus_n*rec_hPlus_n + gen_hPlus_k*rec_hPlus_k); 
      reco::Candidate::LorentzVector rec_tauMinus_p4 = build_p4(rec_tauMinus_pt, rec_tauMinus_eta, rec_tauMinus_phi, rec_tauMinus_mass);
      reco::Candidate::Point         rec_svMinus     = build_point(rec_svMinus_x, rec_svMinus_y, rec_svMinus_z);
      double rec_svMinus_r = rec_svMinus.x()*tauMinus_r.x() + rec_svMinus.y()*tauMinus_r.y() + rec_svMinus.z()*tauMinus_r.z();
      double rec_svMinus_n = rec_svMinus.x()*tauMinus_n.x() + rec_svMinus.y()*tauMinus_n.y() + rec_svMinus.z()*tauMinus_n.z();
      double rec_svMinus_k = rec_svMinus.x()*tauMinus_k.x() + rec_svMinus.y()*tauMinus_k.y() + rec_svMinus.z()*tauMinus_k.z();
      reco::Candidate::LorentzVector rec_nuMinus_p4  = build_p4(rec_nuMinus_pt,  rec_nuMinus_eta,  rec_nuMinus_phi,  rec_nuMinus_mass);
      double hMinus_angle = std::acos(gen_hMinus_r*rec_hMinus_r + gen_hMinus_n*rec_hMinus_n + gen_hMinus_k*rec_hMinus_k); 
      reco::Candidate::LorentzVector rec_higgs_p4    = rec_tauPlus_p4 + rec_tauMinus_p4;

      histogram_res_pv_x->Fill(rec_pv_x - gen_pv_x, evtWeight);
      histogram_res_pv_y->Fill(rec_pv_y - gen_pv_y, evtWeight);
      histogram_res_pv_z->Fill(rec_pv_z - gen_pv_z, evtWeight);

      histogram_res_tauPlus_pt->Fill(rec_tauPlus_pt - gen_tauPlus_pt, evtWeight);
      histogram_res_tauPlus_eta->Fill(rec_tauPlus_eta - gen_tauPlus_eta, evtWeight);
      histogram_res_tauPlus_phi->Fill(rec_tauPlus_phi - gen_tauPlus_phi, evtWeight);
      histogram_res_svPlus_r->Fill(rec_svPlus_r - gen_svPlus_r, evtWeight);
      histogram_res_svPlus_n->Fill(rec_svPlus_n - gen_svPlus_n, evtWeight);
      histogram_res_svPlus_k->Fill(rec_svPlus_k - gen_svPlus_k, evtWeight);
      histogram_res_visPlus_pt->Fill(rec_visPlus_pt - gen_visPlus_pt, evtWeight);
      histogram_res_visPlus_eta->Fill(rec_visPlus_eta - gen_visPlus_eta, evtWeight);
      histogram_res_visPlus_phi->Fill(rec_visPlus_phi - gen_visPlus_phi, evtWeight);
      histogram_res_nuPlus_pt->Fill(rec_nuPlus_pt - gen_nuPlus_pt, evtWeight);
      histogram_res_nuPlus_eta->Fill(rec_nuPlus_eta - gen_nuPlus_eta, evtWeight);
      histogram_res_nuPlus_phi->Fill(rec_nuPlus_phi - gen_nuPlus_phi, evtWeight);
      histogram_res_nuPlus_px->Fill(rec_nuPlus_p4.px() - gen_nuPlus_p4.px(), evtWeight);
      histogram_res_nuPlus_py->Fill(rec_nuPlus_p4.py() - gen_nuPlus_p4.py(), evtWeight);
      histogram_res_nuPlus_pz->Fill(rec_nuPlus_p4.pz() - gen_nuPlus_p4.pz(), evtWeight);
      histogram_res_hPlus_r->Fill(rec_hPlus_r - gen_hPlus_r, evtWeight);
      histogram_res_hPlus_n->Fill(rec_hPlus_n - gen_hPlus_n, evtWeight);
      histogram_res_hPlus_k->Fill(rec_hPlus_k - gen_hPlus_k, evtWeight);
      histogram_res_hPlus_angle->Fill(hPlus_angle, evtWeight);

      histogram_res_tauMinus_pt->Fill(rec_tauMinus_pt - gen_tauMinus_pt, evtWeight);
      histogram_res_tauMinus_eta->Fill(rec_tauMinus_eta - gen_tauMinus_eta, evtWeight);
      histogram_res_tauMinus_phi->Fill(rec_tauMinus_phi - gen_tauMinus_phi, evtWeight);
      histogram_res_svMinus_r->Fill(rec_svMinus_r - gen_svMinus_r, evtWeight);
      histogram_res_svMinus_n->Fill(rec_svMinus_n - gen_svMinus_n, evtWeight);
      histogram_res_svMinus_k->Fill(rec_svMinus_k - gen_svMinus_k, evtWeight);
      histogram_res_visMinus_pt->Fill(rec_visMinus_pt - gen_visMinus_pt, evtWeight);
      histogram_res_visMinus_eta->Fill(rec_visMinus_eta - gen_visMinus_eta, evtWeight);
      histogram_res_visMinus_phi->Fill(rec_visMinus_phi - gen_visMinus_phi, evtWeight);
      histogram_res_nuMinus_pt->Fill(rec_nuMinus_pt - gen_nuMinus_pt, evtWeight);
      histogram_res_nuMinus_eta->Fill(rec_nuMinus_eta - gen_nuMinus_eta, evtWeight);
      histogram_res_nuMinus_phi->Fill(rec_nuMinus_phi - gen_nuMinus_phi, evtWeight);
      histogram_res_nuMinus_px->Fill(rec_nuMinus_p4.px() - gen_nuMinus_p4.px(), evtWeight);
      histogram_res_nuMinus_py->Fill(rec_nuMinus_p4.py() - gen_nuMinus_p4.py(), evtWeight);
      histogram_res_nuMinus_pz->Fill(rec_nuMinus_p4.pz() - gen_nuMinus_p4.pz(), evtWeight);
      histogram_res_hMinus_r->Fill(rec_hMinus_r - gen_hMinus_r, evtWeight);
      histogram_res_hMinus_n->Fill(rec_hMinus_n - gen_hMinus_n, evtWeight);
      histogram_res_hMinus_k->Fill(rec_hMinus_k - gen_hMinus_k, evtWeight);
      histogram_res_hMinus_angle->Fill(hMinus_angle, evtWeight);

      histogram_res_higgs_pt->Fill(rec_higgs_p4.pt() - gen_higgs_p4.pt(), evtWeight);
      histogram_res_higgs_eta->Fill(rec_higgs_p4.eta() - gen_higgs_p4.eta(), evtWeight);
      histogram_res_higgs_phi->Fill(rec_higgs_p4.phi() - gen_higgs_p4.phi(), evtWeight);
      histogram_res_higgs_px->Fill(rec_higgs_p4.px() - gen_higgs_p4.px(), evtWeight);
      histogram_res_higgs_py->Fill(rec_higgs_p4.py() - gen_higgs_p4.py(), evtWeight);
      histogram_res_higgs_pz->Fill(rec_higgs_p4.pz() - gen_higgs_p4.pz(), evtWeight);

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

  showHistogram1d(histogram_res_pv_x,         "PV x^{rec} - x^{gen} [cm]",                        1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_pv_y,         "PV y^{rec} - y^{gen} [cm]",                        1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_pv_z,         "PV z^{rec} - z^{gen} [cm]",                        1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());

  showHistogram1d(histogram_res_tauPlus_pt,   "#tau^{+} p_{T}^{rec} - p_{T}^{gen} [GeV]",         1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_tauPlus_eta,  "#tau^{+} #eta^{rec} - #eta^{gen}",                 1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());
  showHistogram1d(histogram_res_tauPlus_phi,  "#tau^{+} #phi^{rec} - #phi^{gen}",                 1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());
  showHistogram1d(histogram_res_svPlus_r,     "SV(#tau^{+}) r^{rec} - r^{gen} [cm]",              1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_svPlus_n,     "SV(#tau^{+}) n^{rec} - n^{gen} [cm]",              1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_svPlus_k,     "SV(#tau^{+}) k^{rec} - k^{gen} [cm]",              1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_visPlus_pt,   "#tau^{+}_{h} p_{T}^{rec} - p_{T}^{gen} [GeV]",     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_visPlus_eta,  "#tau^{+}_{h} #eta^{rec} - #eta^{gen}",             1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_visPlus_phi,  "#tau^{+}_{h} #phi^{rec} - #eta^{gen}",             1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_nuPlus_pt,    "#bar{#nu}_{#tau} p_{T}^{rec} - p_{T}^{gen} [GeV]", 1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_nuPlus_eta,   "#bar{#nu}_{#tau} #eta^{rec} - #eta^{gen}",         1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_nuPlus_phi,   "#bar{#nu}_{#tau} #phi^{rec} - #phi^{gen}",         1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_nuPlus_px,    "#bar{#nu}_{#tau} p_{x}^{rec} - p_{x}^{gen} [GeV]", 1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_nuPlus_py,    "#bar{#nu}_{#tau} p_{y}^{rec} - p_{y}^{gen} [GeV]", 1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_nuPlus_pz,    "#bar{#nu}_{#tau} p_{z}^{rec} - p_{z}^{gen} [GeV]", 1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_hPlus_r,      "h^{+} r^{rec} - r^{gen}",                          1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());
  showHistogram1d(histogram_res_hPlus_n,      "h^{+} n^{rec} - n^{gen}",                          1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());
  showHistogram1d(histogram_res_hPlus_k,      "h^{+} k^{rec} - k^{gen}",                          1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());
  showHistogram1d(histogram_res_hPlus_angle,  "angle(h^{+,rec},h^{+,gen})",                       1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());

  showHistogram1d(histogram_res_tauMinus_pt,  "#tau^{-} p_{T}^{rec} - p_{T}^{gen} [GeV]",         1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_tauMinus_eta, "#tau^{-} #eta^{rec} - #eta^{gen}",                 1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());
  showHistogram1d(histogram_res_tauMinus_phi, "#tau^{-} #phi^{rec} - #phi^{gen}",                 1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());
  showHistogram1d(histogram_res_svMinus_r,    "SV(#tau^{-}) r^{rec} - r^{gen} [cm]",              1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_svMinus_n,    "SV(#tau^{-}) n^{rec} - n^{gen} [cm]",              1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_svMinus_k,    "SV(#tau^{-}) k^{rec} - k^{gen} [cm]",              1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_visMinus_pt,  "#tau^{-}_{h} p_{T}^{rec} - p_{T}^{gen} [GeV]",     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_visMinus_eta, "#tau^{-}_{h} #eta^{rec} - #eta^{gen}",             1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_visMinus_phi, "#tau^{-}_{h} #phi^{rec} - #eta^{gen}",             1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_nuMinus_pt,   "#bar{#nu}_{#tau} p_{T}^{rec} - p_{T}^{gen} [GeV]", 1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_nuMinus_eta,  "#bar{#nu}_{#tau} #eta^{rec} - #eta^{gen}",         1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_nuMinus_phi,  "#bar{#nu}_{#tau} #phi^{rec} - #phi^{gen}",         1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_nuMinus_px,   "#bar{#nu}_{#tau} p_{x}^{rec} - p_{x}^{gen} [GeV]", 1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_nuMinus_py,   "#bar{#nu}_{#tau} p_{y}^{rec} - p_{y}^{gen} [GeV]", 1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_nuMinus_pz,   "#bar{#nu}_{#tau} p_{z}^{rec} - p_{z}^{gen} [GeV]", 1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_hMinus_r,     "h^{-} r^{rec} - r^{gen}",                          1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());
  showHistogram1d(histogram_res_hMinus_n,     "h^{-} n^{rec} - n^{gen}",                          1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());
  showHistogram1d(histogram_res_hMinus_k,     "h^{-} k^{rec} - k^{gen}",                          1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());
  showHistogram1d(histogram_res_hMinus_angle, "angle(h^{-,rec},h^{-,gen})",                       1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());

  showHistogram1d(histogram_res_higgs_pt,     "H p_{T}^{rec} - p_{T}^{gen} [GeV]",                1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_higgs_eta,    "H #eta^{rec} - #eta^{gen}",                        1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_higgs_phi,    "H #phi^{rec} - #phi^{gen}",                        1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_higgs_px,     "H p_{x}^{rec} - p_{x}^{gen} [GeV]",                1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_higgs_py,     "H p_{y}^{rec} - p_{y}^{gen} [GeV]",                1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_res_higgs_pz,     "H p_{z}^{rec} - p_{z}^{gen} [GeV]",                1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());

  clock.Show("makeResolutionPlots");

  return EXIT_SUCCESS;
}
