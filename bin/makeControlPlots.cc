
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
                     double avEvtWeight,
                     bool showStatsBox,
                     const std::string& outputFileName)
{
  //double integral = histogram->Integral();
  //if ( integral > 0. ) histogram->Scale(1./integral);
  if ( avEvtWeight > 0. )
  {
    histogram->Scale(1./avEvtWeight);
  }

  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2); 
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);
  canvas->SetLogy(useLogScale);
  
  histogram->SetTitle("");
  histogram->SetStats(showStatsBox);
  histogram->SetMinimum(yMin*histogram->Integral());
  histogram->SetMaximum(yMax*histogram->Integral());

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
  histogram->SetMarkerSize(1);
  histogram->SetMarkerStyle(8);
  histogram->Draw("E1P");

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  outputFileName_plot.append("_");
  outputFileName_plot.append(histogram->GetName());
  canvas->Print(std::string(outputFileName_plot).append(".png").c_str());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").c_str());
  
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
                     double avEvtWeight,
                     bool showStatsBox,
                     const std::string& outputFileName)
{
  //double integral = histogram->Integral();
  //if ( integral > 0. ) histogram->Scale(1./integral);
  if ( avEvtWeight > 0. )
  {
    histogram->Scale(1./avEvtWeight);
  }

  TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 800);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2); 
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);
  
  histogram->SetTitle("");
  histogram->SetStats(showStatsBox);

  TAxis* xAxis = histogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.c_str());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);  

  TAxis* yAxis = histogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.c_str());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  histogram->Draw("BOX");

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  outputFileName_plot.append("_");
  outputFileName_plot.append(histogram->GetName());
  canvas->Print(std::string(outputFileName_plot).append(".png").c_str());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").c_str());
  
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
  float minVisTauPt = cfg_ctrlPlots.getParameter<double>("minVisTauPt");
  std::cout << " minVisTauPt = " << minVisTauPt << "\n";
  float maxAbsVisTauEta = cfg_ctrlPlots.getParameter<double>("maxAbsVisTauEta");
  std::cout << " maxAbsVisTauEta = " << maxAbsVisTauEta << "\n";
  int maxNumChargedKaons = cfg_ctrlPlots.getParameter<int>("maxNumChargedKaons");
  std::cout << " maxNumChargedKaons = " << maxNumChargedKaons << "\n";
  int maxNumNeutralKaons = cfg_ctrlPlots.getParameter<int>("maxNumNeutralKaons");
  std::cout << " maxNumNeutralKaons = " << maxNumNeutralKaons << "\n";
  int maxNumPhotons = cfg_ctrlPlots.getParameter<int>("maxNumPhotons");
  std::cout << " maxNumPhotons = " << maxNumPhotons << "\n";
  float maxSumPhotonEn = cfg_ctrlPlots.getParameter<double>("maxSumPhotonEn");
  std::cout << " maxSumPhotonEn = " << maxSumPhotonEn << "\n";
  std::string branchName_evtWeight = cfg_ctrlPlots.getParameter<std::string>("branchName_evtWeight");
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

  TH1* histogram_C_rr            = bookHistogram1d(fs, "C_rr",            72, -9.,  +9.);
  TH1* histogram_C_nn            = bookHistogram1d(fs, "C_nn",            72, -9.,  +9.);
  TH1* histogram_C_kk            = bookHistogram1d(fs, "C_kk",            72, -9.,  +9.);

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

    Float_t tauPlus_pt, tauPlus_eta;
    inputTree->SetBranchAddress("tauPlus_pt", &tauPlus_pt);
    inputTree->SetBranchAddress("tauPlus_eta", &tauPlus_eta);
    Float_t visTauPlus_pt, visTauPlus_eta;
    inputTree->SetBranchAddress("visTauPlus_pt", &visTauPlus_pt);
    inputTree->SetBranchAddress("visTauPlus_eta", &visTauPlus_eta);
    Int_t tauPlus_nChargedKaons, tauPlus_nNeutralKaons, tauPlus_nPhotons;
    Float_t tauPlus_sumPhotonEn;
    inputTree->SetBranchAddress("tauPlus_nChargedKaons", &tauPlus_nChargedKaons);
    inputTree->SetBranchAddress("tauPlus_nNeutralKaons", &tauPlus_nNeutralKaons);
    inputTree->SetBranchAddress("tauPlus_nPhotons", &tauPlus_nPhotons);
    inputTree->SetBranchAddress("tauPlus_sumPhotonEn", &tauPlus_sumPhotonEn);
    Float_t hPlus_r, hPlus_n, hPlus_k;
    inputTree->SetBranchAddress("hPlus_r", &hPlus_r);
    inputTree->SetBranchAddress("hPlus_n", &hPlus_n);
    inputTree->SetBranchAddress("hPlus_k", &hPlus_k);

    Float_t tauMinus_pt, tauMinus_eta;
    inputTree->SetBranchAddress("tauMinus_pt", &tauMinus_pt);
    inputTree->SetBranchAddress("tauMinus_eta", &tauMinus_eta);
    Float_t visTauMinus_pt, visTauMinus_eta;
    inputTree->SetBranchAddress("visTauMinus_pt", &visTauMinus_pt);
    inputTree->SetBranchAddress("visTauMinus_eta", &visTauMinus_eta);
    Int_t tauMinus_nChargedKaons, tauMinus_nNeutralKaons, tauMinus_nPhotons;
    Float_t tauMinus_sumPhotonEn;
    inputTree->SetBranchAddress("tauMinus_nChargedKaons", &tauMinus_nChargedKaons);
    inputTree->SetBranchAddress("tauMinus_nNeutralKaons", &tauMinus_nNeutralKaons);
    inputTree->SetBranchAddress("tauMinus_nPhotons", &tauMinus_nPhotons);
    inputTree->SetBranchAddress("tauMinus_sumPhotonEn", &tauMinus_sumPhotonEn);
    Float_t hMinus_r, hMinus_n, hMinus_k;
    inputTree->SetBranchAddress("hMinus_r", &hMinus_r);
    inputTree->SetBranchAddress("hMinus_n", &hMinus_n);
    inputTree->SetBranchAddress("hMinus_k", &hMinus_k);

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
  showHistogram1d(histogram_visTauPlusPt,    "#tau^{+}_{h} p_{T} [GeV]", 1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_visTauPlusEta,   "#tau^{+}_{h} #eta",        1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());

  showHistogram1d(histogram_tauMinusPt,      "#tau^{-} p_{T} [GeV]",     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_tauMinusEta,     "#tau^{-} #eta",            1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());
  showHistogram1d(histogram_visTauMinusPt,   "#tau^{-}_{h} p_{T} [GeV]", 1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_visTauMinusEta,  "#tau^{-}_{h} #eta",        1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());

  showHistogram1d(histogram_mTauTau,         "m_{#tau#tau} [GeV]",       1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_mVis,            "m_{vis} [GeV]",            1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_cosTheta,        "cos(#theta)",              1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, false, outputFile.file());

  showHistogram1d(histogram_C_rr,            "C_rr",                     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_C_nn,            "C_nn",                     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());
  showHistogram1d(histogram_C_kk,            "C_kk",                     1.2, true, 1.e-3, 1.e0, "Events", 1.3, avEvtWeight, true,  outputFile.file());

  showHistogram2d(histogram_zPlus_vs_zMinus, "z^{-}",                    1.2,                    "z^{+}",  1.3, avEvtWeight, false, outputFile.file());

  clock.Show("makeControlPlots");

  return EXIT_SUCCESS;
}
