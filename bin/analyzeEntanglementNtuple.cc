
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
#include <TVectorD.h>                                                                           // TVectorD
#include <TObject.h>                                                                            // TObject
#include <Minuit2/Minuit2Minimizer.h>                                                           // ROOT::Minuit2::Minuit2Minimizer
#include <Math/Functor.h>                                                                       // ROOT::Math::Functor
#include <TCanvas.h>                                                                            // TCanvas
#include <TGraph.h>                                                                             // TGraph
#include <TH1.h>                                                                                // TH1F
#include <TLegend.h>                                                                            // TLegend
#include <TString.h>                                                                            // Form()

#include <assert.h>                                                                             // assert()
#include <cstdlib>                                                                              // EXIT_SUCCESS, EXIT_FAILURE
#include <fstream>                                                                              // std::ofstream
#include <iostream>                                                                             // std::cout
#include <string>                                                                               // std::string
#include <vector>                                                                               // std::vector<>
#include <cmath>                                                                                // std::fabs()

const size_t npar = 15;

typedef std::vector<double> vdouble;

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
  ~EntanglementDataset()
  {}

  void
  push_back(float hPlus_r, float hPlus_n, float hPlus_k, 
            float hMinus_r, float hMinus_n, float hMinus_k,
            float evtWeight = 1.)
  {
    EntanglementData entry = EntanglementData(
      hPlus_r, hPlus_n, hPlus_k, 
      hMinus_r, hMinus_n, hMinus_k, 
      evtWeight);
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

  double
  get_par_gen(size_t idx) const
  {
    return par_gen_[idx];
  }

 private:
   std::vector<EntanglementData> data_;

   std::vector<double> par_gen_;
};

static EntanglementDataset* gEntanglementDataset = nullptr;

double
get_p(const double* par, const EntanglementData& entry)
{
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

  return p;
}

double
mlfit_fcn(const double* par)
{
  const EntanglementDataset* mlfitData = gEntanglementDataset;
  assert(mlfitData);

  size_t numEntries = mlfitData->size();

  double par_gen[npar];
  for ( size_t idxPar = 0; idxPar < npar; ++idxPar )
  {
    par_gen[idxPar] = mlfitData->get_par_gen(idxPar);
  }

  double norm = 0.;
  for ( size_t idxEntry = 0; idxEntry < numEntries; ++idxEntry )
  {
    const EntanglementData& entry = mlfitData->at(idxEntry);

    norm += entry.get_evtWeight()*get_p(par, entry)/get_p(par_gen, entry);
    //norm += get_p(par, entry)/get_p(par_gen, entry);
  }

  double logL = 0.;
  for ( size_t idxEntry = 0; idxEntry < numEntries; ++idxEntry )
  {
    const EntanglementData& entry = mlfitData->at(idxEntry);

    double p = get_p(par, entry);

    logL -= 2.*entry.get_evtWeight()*log(p/norm);
    //logL -= 2.*log(p/norm);
  }

  return logL;
}

TGraph*
scan_mlfit_fixed(const std::vector<double>& parValues,
                 size_t parToScan, const std::string& parName, size_t numPoints, double xMin, double xMax,
                 double yMin)
{
  std::cout << "<scan_mlfit_fixed>:\n";

  assert(xMax > xMin && numPoints >= 2);
  double dx = (xMax - xMin)/(numPoints - 1);

  double par[npar];
  for ( size_t idxPar = 0; idxPar < npar; ++idxPar )
  {
    par[idxPar] = parValues[idxPar];
  }

  TGraph* graph = new TGraph(numPoints);
  for ( size_t idxPoint = 0; idxPoint < numPoints; ++idxPoint )
  {
    double x = xMin + idxPoint*dx;
    par[parToScan] = x;
    double y = mlfit_fcn(par);
    std::cout << " point #" << idxPoint << ": x = " << x << ", y = " << y - yMin << "\n";
    graph->SetPoint(idxPoint, x, y - yMin);
  }
  return graph;
}

TGraph*
scan_mlfit_profiled(ROOT::Math::Minimizer* mlfit,
                    const std::vector<double>& parValues,
                    size_t parToScan, const std::string& parName, size_t numPoints, double xMin, double xMax,
                    double yMin)
{
  std::cout << "<scan_mlfit_profiled>:\n";

  assert(xMax > xMin && numPoints >= 2);
  double dx = (xMax - xMin)/(numPoints - 1);

  for ( size_t idxPar = 0; idxPar < npar; ++idxPar )
  {
    mlfit->SetVariableValue(idxPar, parValues[idxPar]);
  }

  mlfit->FixVariable(parToScan);

  TGraph* graph = new TGraph(numPoints);
  for ( size_t idxPoint = 0; idxPoint < numPoints; ++idxPoint )
  {
    double x = xMin + idxPoint*dx;
    mlfit->SetFixedVariable(parToScan, parName.c_str(), x);
    mlfit->Minimize();
    double y = mlfit->MinValue();
    std::cout << " point #" << idxPoint << ": x = " << x << ", y = " << y - yMin << "\n";
    graph->SetPoint(idxPoint, x, y - yMin);
  }
  return graph;
}

void showGraphs(double canvasSizeX, double canvasSizeY,
                TGraph* graph1, const std::string& legendEntry1,
                TGraph* graph2, const std::string& legendEntry2,
                double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, 
                double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
                bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
                const std::string& outputFileName)
{
  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetTopMargin(0.05);
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);
  canvas->SetRightMargin(0.05);
  canvas->SetLogy(useLogScale);
  
  canvas->SetGridx(1);
  canvas->SetGridy(1);

  if ( !graph1 ) 
  {
    std::cerr << "<showGraphs>: graph1 = NULL --> skipping !!" << std::endl;
    return;
  }

  TH1* dummyHistogram = new TH1F("dummyHistogram", "dummyHistogram", 10, xMin, xMax);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);

  TAxis* xAxis = dummyHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleSize(0.045);
  xAxis->SetTitleOffset(xAxisOffset);  
  if ( xMax > xMin ) 
  {
    std::cout << "limiting x-axis range to " << xMin << ".." << xMax << std::endl;
    xAxis->SetRangeUser(xMin, xMax);
  }

  TAxis* yAxis = dummyHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleSize(0.045);
  yAxis->SetTitleOffset(yAxisOffset);

  dummyHistogram->Draw();
  //dummyHistogram->Draw("axis");

  graph1->SetMarkerColor(1);
  graph1->SetMarkerSize(1);
  graph1->SetMarkerStyle(8);
  graph1->SetLineColor(1);
  graph1->SetLineWidth(1);
  graph1->SetLineStyle(1);
  graph1->Draw("PL");

  if ( graph2 ) 
  {
    graph2->SetMarkerColor(2);
    graph2->SetMarkerSize(1);
    graph2->SetMarkerStyle(8);
    graph2->SetLineColor(2);
    graph2->SetLineWidth(1);
    graph2->SetLineStyle(1);
    graph2->Draw("PL");
  }

  TLegend* legend = 0;
  if ( legendEntry1 != "" ) 
  {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, "", "brNDC"); 
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextSize(legendTextSize);
    legend->AddEntry(graph1, legendEntry1.data(), "L");
    if ( graph2 ) legend->AddEntry(graph2, legendEntry2.data(), "L");
    legend->Draw();
  }

  dummyHistogram->Draw("axissame");

  canvas->RedrawAxis();

  canvas->Update();
  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  //if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  canvas->Print(std::string(outputFileName_plot).append(".png").data());
  //canvas->Print(std::string(outputFileName_plot).append(".pdf").data());

  delete dummyHistogram;
  delete legend;
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
  
  std::vector<double> par_gen = cfg_analyze.getParameter<vdouble>("par_gen");
  if ( par_gen.size() != npar )
    throw cmsException("analyzeEntanglementNtuple", __LINE__) << "Invalid Configuration parameter 'par_gen' !!";

  bool scanLikelihood = cfg_analyze.getParameter<bool>("scanLikelihood");
  std::cout << " scanLikelihood = " << scanLikelihood << "\n";
  
  //bool isDEBUG = cfg_analyze.getParameter<bool>("isDEBUG");

  fwlite::InputSource inputFiles(cfg);
  unsigned reportEvery = inputFiles.reportAfter();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().c_str());

  std::vector<std::string> inputFileNames = inputFiles.files();
  size_t numInputFiles = inputFileNames.size();
  std::cout << "Loaded " << numInputFiles << " file(s).\n";

  TVectorD Bp(3);
  TVectorD BpErr(3);
  TVectorD Bm(3);
  TVectorD BmErr(3);
  TMatrixD C(3,3);
  TMatrixD CErr(3,3);
  
  TTree* fitResult = fs.make<TTree>("fitResult", "fitResult");
  fitResult->Branch("Bp",    "TVectorD", &Bp);
  fitResult->Branch("BpErr", "TVectorD", &BpErr);
  fitResult->Branch("Bm",    "TVectorD", &Bm);
  fitResult->Branch("BmErr", "TVectorD", &BmErr);
  fitResult->Branch("C",     "TMatrixD", &C);
  fitResult->Branch("CErr",  "TMatrixD", &CErr);

  EntanglementDataset mlfitData(par_gen); 

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
  TMatrixD C_exp(3,3);
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
  ROOT::Math::Minimizer* mlfit = new ROOT::Minuit2::Minuit2Minimizer();
  for ( size_t idxPar = 0; idxPar < npar; ++idxPar )
  {
    double par0 = 0.;
    if ( idxPar >= 6 && idxPar <= 14 )
    {
      size_t idxRow = (idxPar - 6) / 3;
      size_t idxCol = (idxPar - 6) % 3;
      par0 = C[idxRow][idxCol];
      //par0 = C_exp[idxRow][idxCol];
    }
    mlfit->SetLimitedVariable(idxPar, parNames[idxPar].c_str(), par0, 0.1, -2., +2.);
  }

  mlfit->SetMaxFunctionCalls(10000);
  mlfit->SetTolerance(1.e-3);
  mlfit->SetPrintLevel(-1);

  // set function pointer
  ROOT::Math::Functor f(&mlfit_fcn, npar);
  mlfit->SetFunction(f);

  gEntanglementDataset = &mlfitData;

  mlfit->Minimize();
  mlfit->Hesse();

  std::cout << "Fit Results:\n";
  mlfit->PrintResults();

  std::vector<double> parValues(npar);
  std::vector<double> parErrors(npar);
  for ( size_t idxPar = 0; idxPar < npar; ++idxPar )
  {
    const double* X = mlfit->X();
    parValues[idxPar] = X[idxPar];
    parErrors[idxPar] = sqrt(mlfit->CovMatrix(idxPar, idxPar));
    std::cout << parNames[idxPar] << " = " << parValues[idxPar] << " +/- " << parErrors[idxPar] << "\n";  
  }

  for ( size_t idxPar = 0; idxPar < npar; ++idxPar )
  {
    if ( idxPar <= 2 )
    {
      size_t idx = idxPar;
      Bp[idx] = parValues[idxPar];
      BpErr[idx] = parErrors[idxPar];
    }
    else if ( idxPar <= 5 )
    {
      size_t idx = idxPar - 3;
      Bm[idx] = parValues[idxPar];
      BmErr[idx] = parErrors[idxPar];
    }
    else if ( idxPar <= 14 )
    {	
      size_t idxRow = (idxPar - 6) / 3;
      size_t idxCol = (idxPar - 6) % 3;
      C[idxRow][idxCol] = parValues[idxPar];
      CErr[idxRow][idxCol] = parErrors[idxPar];
    } else assert(0);
  }
  fitResult->Fill();

  if ( scanLikelihood )
  {
    double mlfitMin = mlfit->MinValue();
    std::vector<size_t> parsToScan = { 6, 10, 14 };
    for ( size_t i = 0; i < parsToScan.size(); ++i )
    {
      size_t parToScan = parsToScan[i];
      std::cout << "Scanning likelihood as function of parameter " << parNames[parToScan] << "...\n";
      double xMin = parValues[parToScan] - std::max(0.01, 5.*parErrors[parToScan]);
      double xMax = parValues[parToScan] + std::max(0.01, 5.*parErrors[parToScan]);
      TGraph* mlfitScan_fixed = scan_mlfit_fixed(parValues, parToScan, parNames[parToScan], 100, xMin, xMax, mlfitMin);
      TGraph* mlfitScan_profiled = scan_mlfit_profiled(mlfit, parValues, parToScan, parNames[parToScan], 100, xMin, xMax, mlfitMin);
      std::string outputFileName = (const char*)TString(outputFile.file().c_str()).ReplaceAll(".root", Form("_mlfitScan_%s.png", parNames[parToScan].c_str()));
      showGraphs(1150, 950,
                 mlfitScan_fixed,    "Fixed",
                 mlfitScan_profiled, "Profiled",
                 0.040, 0.48, 0.78, 0.22, 0.14, 
	         xMin, xMax, parNames[parToScan], 1.2,
                 false, 0., 2.5e+1, "-2 log(L)", 1.4,
                 outputFileName);
      std::cout << " Done.\n";
    }
  }

  delete mlfit;

  clock.Show("analyzeEntanglementNtuple");

  return EXIT_SUCCESS;
}
