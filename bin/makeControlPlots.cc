
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
#include <TCanvas.h>                                                                            // TCanvas
#include <TH1.h>                                                                                // TH1D
#include <TH2.h>                                                                                // TH2D
#include <TAxis.h>                                                                              // TAxis

#include <assert.h>                                                                             // assert()
#include <cstdlib>                                                                              // EXIT_SUCCESS, EXIT_FAILURE
#include <fstream>                                                                              // std::ofstream
#include <iostream>                                                                             // std::cout
#include <string>                                                                               // std::string
#include <vector>                                                                               // std::vector
#include <cmath>                                                                                // std::fabs

TH1*
bookHistogram1d(fwlite::TFileService& fs, const std::string& name, int numBinsX, double xMin, double xMax)
{
  TH1* histogram = fs.make<TH1D>(name.c_str(), name.c_str(), numBinsX, xMin, xMax);
  histogram->Sumw2();
  return histogram;
}

void showHistogram1d(TH1* histogram, 
                     const std::string& xAxisTitle, double xAxisOffset, 
                     bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
                     const std::string& outputFileName)
{
  double integral = histogram->Integral();
  if ( integral > 0. ) histogram->Scale(1./integral);

  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2); 
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);
  canvas->SetLogy(useLogScale);
  
  histogram->SetTitle("");
  histogram->SetStats(false);
  histogram->SetMinimum(yMin);
  histogram->SetMaximum(yMax);

  TAxis* xAxis = histogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);  

  TAxis* yAxis = histogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  histogram->SetLineColor(1);
  histogram->SetLineWidth(1);
  histogram->SetLineStyle(1);
  histogram->SetMarkerColor(1);
  histogram->SetMarkerSize(2);
  histogram->SetMarkerStyle(8);
  histogram->Draw("E1P");

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  //canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  
  delete canvas;  
}

TH2*
bookHistogram2d(fwlite::TFileService& fs, const std::string& name, int numBinsX, double xMin, double xMax, int numBinsY, double yMin, double yMax)
{
  TH2* histogram = fs.make<TH2D>(name.c_str(), name.c_str(), numBinsX, xMin, xMax, numBinsY, yMin, yMax);
  histogram->Sumw2();
  return histogram;
}

void showHistogram2d(TH2* histogram, 
                     const std::string& xAxisTitle, double xAxisOffset, 
                     const std::string& yAxisTitle, double yAxisOffset,
                     const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2); 
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);
  
  histogram->SetTitle("");
  histogram->SetStats(false);

  TAxis* xAxis = histogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);  

  TAxis* yAxis = histogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  histogram->Draw("BOX");

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  //canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  
  delete canvas;  
}

int main(int argc, char* argv[])
{
//--- throw an exception in case ROOT encounters an error
  gErrorAbortLevel = kError;

//--- stop ROOT from keeping track of all histograms
  TH1::AddDirectory(false);

//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "<makeControlPlots>:" << std::endl;

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("makeControlPlots");

//--- read python configuration parameters
  std::cout << "Reading config file " << argv[1] << std::endl;
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") )
    throw cmsException("makeControlPlots", __LINE__) << "No ParameterSet 'process' found in config file !!";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameterSet("process");

  edm::ParameterSet cfg_ctrlPlots = cfg.getParameterSet("makeControlPlots");
  std::string treeName = cfg_ctrlPlots.getParameter<std::string>("treeName");
  std::cout << " treeName = " << treeName << std::endl;
  float minVisTauPt = cfg_ctrlPlots.getParameter<double>("minVisTauPt");
  std::cout << " minVisTauPt = " << minVisTauPt << std::endl;
  float maxAbsVisTauEta = cfg_ctrlPlots.getParameter<double>("maxAbsVisTauEta");
  std::cout << " maxAbsVisTauEta = " << maxAbsVisTauEta << std::endl;
  std::string branchName_evtWeight = cfg_ctrlPlots.getParameter<std::string>("branchName_evtWeight");
  std::cout << " branchName_evtWeight = " << branchName_evtWeight << std::endl;
  //bool isDEBUG = cfg_analyze.getParameter<bool>("isDEBUG");

  fwlite::InputSource inputFiles(cfg);
  int maxEvents = inputFiles.maxEvents();
  std::cout << " maxEvents = " << maxEvents << std::endl;
  unsigned reportEvery = inputFiles.reportAfter();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

  std::vector<std::string> inputFileNames = inputFiles.files();
  size_t numInputFiles = inputFileNames.size();
  std::cout << "Loaded " << numInputFiles << " file(s).\n";

  TH1* histogram_tauPlusPt       = bookHistogram1d(fs, "tauPlusPt",       40,  0., 200.);
  TH1* histogram_tauPlusEta      = bookHistogram1d(fs, "tauPlusEta",      50, -5.,  +5.);
  TH1* histogram_visTauPlusPt    = bookHistogram1d(fs, "visTauPlusPt",    40,  0., 200.);
  TH1* histogram_visTauPlusEta   = bookHistogram1d(fs, "visTauPlusEta",   50, -5.,  +5.);
  TH1* histogram_tauMinusPt      = bookHistogram1d(fs, "tauMinusPt",      40,  0., 200.);
  TH1* histogram_tauMinusEta     = bookHistogram1d(fs, "tauMinusEta",     50, -5.,  +5.);
  TH1* histogram_visTauMinusPt   = bookHistogram1d(fs, "visTauMinusPt",   40,  0., 200.);
  TH1* histogram_visTauMinusEta  = bookHistogram1d(fs, "visTauMinusEta",  50, -5.,  +5.);
  TH1* histogram_mTauTau         = bookHistogram1d(fs, "mTauTau",         40,  0., 200.);
  TH1* histogram_mVis            = bookHistogram1d(fs, "mVis",            40,  0., 200.);
  TH1* histogram_cosTheta        = bookHistogram1d(fs, "cosTheta",        40, -1.,  +1.);

  TH2* histogram_zPlus_vs_zMinus = bookHistogram2d(fs, "zPlus_vs_zMinus", 40, 0., 1., 40, 0., 1.);

  int analyzedEntries = 0;
  float analyzedEntries_weighted = 0.;
  int selectedEntries = 0;
  float selectedEntries_weighted = 0.;
  int processedInputFiles = 0;
  bool STOP = false;
  for ( size_t idxInputFile = 0; idxInputFile < numInputFiles && !STOP; ++idxInputFile )
  {
    const std::string & inputFileName = inputFileNames.at(idxInputFile);
    std::cout << "Opening #" << idxInputFile << " file " << inputFileName << '\n';
    TFile* inputFile = TFile::Open(inputFileName.data());
    if ( !inputFile )
      throw cmsException("analyzeEntanglementNtuple", __LINE__) 
        << "The file " << inputFileName << " failed to open !!";
   
    TTree* inputTree = dynamic_cast<TTree*>(inputFile->Get(treeName.data()));
    if ( !inputTree )
      throw cmsException("analyzeEntanglementNtuple", __LINE__) 
        << "The file " << inputFileName << " does not contain a TTree named '" << treeName << "' !!";
    std::cout << "The file " << inputFileName << " contains " << inputTree->GetEntries() << " entries\n";

    ++processedInputFiles;

    Float_t tauPlus_pt, tauPlus_eta;
    inputTree->SetBranchAddress("tauPlus_pt", &tauPlus_pt);
    inputTree->SetBranchAddress("tauPlus_eta", &tauPlus_eta);
    Float_t visTauPlus_pt, visTauPlus_eta;
    inputTree->SetBranchAddress("visTauPlus_pt", &visTauPlus_pt);
    inputTree->SetBranchAddress("visTauPlus_eta", &visTauPlus_eta);

    Float_t tauMinus_pt, tauMinus_eta;
    inputTree->SetBranchAddress("tauMinus_pt", &tauMinus_pt);
    inputTree->SetBranchAddress("tauMinus_eta", &tauMinus_eta);
    Float_t visTauMinus_pt, visTauMinus_eta;
    inputTree->SetBranchAddress("visTauMinus_pt", &visTauMinus_pt);
    inputTree->SetBranchAddress("visTauMinus_eta", &visTauMinus_eta);

    Float_t mTauTau, mVis, cosTheta;
    inputTree->SetBranchAddress("mTauTau", &mTauTau);
    inputTree->SetBranchAddress("mVis", &mVis);
    inputTree->SetBranchAddress("cosTheta", &cosTheta);

    Float_t zPlus, zMinus;
    inputTree->SetBranchAddress("zPlus", &zPlus);
    inputTree->SetBranchAddress("zMinus", &zMinus);

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
        std::cout << "processing Entry " << analyzedEntries << std::endl;
      }

      if ( !(visTauPlus_pt  > minVisTauPt && std::fabs(visTauPlus_eta)  < maxAbsVisTauEta) ) continue;
      if ( !(visTauMinus_pt > minVisTauPt && std::fabs(visTauMinus_eta) < maxAbsVisTauEta) ) continue;

      histogram_tauPlusPt->Fill(tauPlus_pt, evtWeight);
      histogram_tauPlusEta->Fill(tauPlus_eta, evtWeight);
      histogram_visTauPlusPt->Fill(visTauPlus_pt, evtWeight);
      histogram_visTauPlusEta->Fill(visTauPlus_eta, evtWeight);

      histogram_tauMinusPt->Fill(tauMinus_pt, evtWeight);
      histogram_tauMinusEta->Fill(tauMinus_eta, evtWeight);
      histogram_visTauMinusPt->Fill(visTauMinus_pt, evtWeight);
      histogram_visTauMinusEta->Fill(visTauMinus_eta, evtWeight);

      histogram_mTauTau->Fill(mTauTau, evtWeight);
      histogram_mVis->Fill(mVis, evtWeight);
      histogram_cosTheta->Fill(cosTheta, evtWeight);

      histogram_zPlus_vs_zMinus->Fill(zMinus, zPlus, evtWeight);

      ++selectedEntries;
      selectedEntries_weighted += evtWeight;

      if ( maxEvents != -1 && analyzedEntries >= maxEvents ) STOP = true;
    }

    delete inputTree;
    delete inputFile;
  }

  std::cout << "Processing Summary:\n";
  std::cout << " processedInputFiles = " << processedInputFiles << " (out of " << numInputFiles << ")\n";
  std::cout << " analyzedEntries = " << analyzedEntries << " (weighted = " << analyzedEntries_weighted << ")\n";
  std::cout << " selectedEntries = " << selectedEntries << " (weighted = " << selectedEntries_weighted << ")\n";

  showHistogram1d(histogram_tauPlusPt,       "#tau^{+} p_{T} [GeV]",     1.2, true, 1.e-3, 1.e0, "Events", 1.3, "makeControlPlots_tauPlusPt.png");
  showHistogram1d(histogram_tauPlusEta,      "#tau^{+} #eta",            1.2, true, 1.e-3, 1.e0, "Events", 1.3, "makeControlPlots_tauPlusEta.png");
  showHistogram1d(histogram_visTauPlusPt,    "#tau^{+}_{h} p_{T} [GeV]", 1.2, true, 1.e-3, 1.e0, "Events", 1.3, "makeControlPlots_visTauPlusPt.png");
  showHistogram1d(histogram_visTauPlusEta,   "#tau^{+}_{h} #eta",        1.2, true, 1.e-3, 1.e0, "Events", 1.3, "makeControlPlots_visTauPlusEta.png");

  showHistogram1d(histogram_tauMinusPt,      "#tau^{-} p_{T} [GeV]",     1.2, true, 1.e-3, 1.e0, "Events", 1.3, "makeControlPlots_tauMinusPt.png");
  showHistogram1d(histogram_tauMinusEta,     "#tau^{-} #eta",            1.2, true, 1.e-3, 1.e0, "Events", 1.3, "makeControlPlots_tauMinusEta.png");
  showHistogram1d(histogram_visTauMinusPt,   "#tau^{-}_{h} p_{T} [GeV]", 1.2, true, 1.e-3, 1.e0, "Events", 1.3, "makeControlPlots_visTauMinusPt.png");
  showHistogram1d(histogram_visTauMinusEta,  "#tau^{-}_{h} #eta",        1.2, true, 1.e-3, 1.e0, "Events", 1.3, "makeControlPlots_visTauMinusEta.png");

  showHistogram1d(histogram_mTauTau,         "m_{#tau#tau} [GeV]",       1.2, true, 1.e-3, 1.e0, "Events", 1.3, "makeControlPlots_mTauTau.png");
  showHistogram1d(histogram_mVis,            "m_{vis} [GeV]",            1.2, true, 1.e-3, 1.e0, "Events", 1.3, "makeControlPlots_mVis.png");
  showHistogram1d(histogram_cosTheta,        "cos(#theta)",              1.2, true, 1.e-3, 1.e0, "Events", 1.3, "makeControlPlots_cosTheta.png");

  showHistogram2d(histogram_zPlus_vs_zMinus, "z^{-}",                    1.2,                    "z^{+}",  1.2, "makeControlPlots_zPlus_vs_zMinus.png");

  clock.Show("makeControlPlots");

  return EXIT_SUCCESS;
}
