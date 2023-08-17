
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
#include <Math/Functor.h>                                                         // ROOT::Math::Functor
#include <Minuit2/Minuit2Minimizer.h>                                             // ROOT::Minuit2::Minuit2Minimizer
#include <TBenchmark.h>                                                           // TBenchmark
#include <TCanvas.h>                                                              // TCanvas
#include <TError.h>                                                               // gErrorAbortLevel, kError
#include <TGraph.h>                                                               // TGraph
#include <TH1.h>                                                                  // TH1F
#include <TLegend.h>                                                              // TLegend
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

static const EntanglementDataset* gEntanglementDataset = nullptr;

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
comp_Bp_Bm_C_by_summation(const EntanglementDataset& dataset, 
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

std::map<size_t, std::string>
get_parNames()
{
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
  return parNames;
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

Measurement
comp_Bp_Bm_C_by_mlfit(const EntanglementDataset& dataset, ROOT::Math::Minimizer* mlfit, 
                      bool scanLikelihood = false, const std::string& outputFileName = "", 
                      int verbosity = -1)
{
  // initialize fit parameters
  Measurement startpos = comp_Bp_Bm_C_by_summation(dataset);
  const math::Matrix3x3& startpos_C = startpos.get_C();
  std::map<size_t, std::string> parNames = get_parNames();
  size_t npar = parNames.size();
  for ( size_t idxPar = 0; idxPar < npar; ++idxPar )
  {
    double par0 = 0.;
    if ( idxPar >= 6 && idxPar <= 14 )
    {
      size_t idxRow = (idxPar - 6) / 3;
      size_t idxCol = (idxPar - 6) % 3;
      par0 = startpos_C(idxRow,idxCol);
    }
    mlfit->SetLimitedVariable(idxPar, parNames[idxPar].c_str(), par0, 0.1, -2., +2.);
  }

  gEntanglementDataset = &dataset;

  mlfit->Minimize();
  mlfit->Hesse();

  if ( verbosity >= 1 )
  {
    std::cout << "Fit Results:\n";
    mlfit->PrintResults();
  }

  std::vector<double> parValues(npar);
  std::vector<double> parErrors(npar);
  for ( size_t idxPar = 0; idxPar < npar; ++idxPar )
  {
    const double* X = mlfit->X();
    parValues[idxPar] = X[idxPar];
    parErrors[idxPar] = sqrt(mlfit->CovMatrix(idxPar, idxPar));
    if ( verbosity >= 1 )
    {
      std::cout << parNames[idxPar] << " = " << parValues[idxPar] << " +/- " << parErrors[idxPar] << "\n";
    }
  }

  if ( scanLikelihood )
  {
    double mlfitMin = mlfit->MinValue();
    std::map<size_t, std::string> parNames = get_parNames();
    std::vector<size_t> parsToScan = { 6, 10, 14 };
    for ( size_t i = 0; i < parsToScan.size(); ++i )
    {
      size_t parToScan = parsToScan[i];      
      std::cout << "Scanning likelihood as function of parameter " << parNames[parToScan] << "...\n";
      double xMin = parValues[parToScan] - std::max(0.01, 5.*parErrors[parToScan]);
      double xMax = parValues[parToScan] + std::max(0.01, 5.*parErrors[parToScan]);
      TGraph* mlfitScan_fixed = scan_mlfit_fixed(parValues, parToScan, parNames[parToScan], 100, xMin, xMax, mlfitMin);
      TGraph* mlfitScan_profiled = scan_mlfit_profiled(mlfit, parValues, parToScan, parNames[parToScan], 100, xMin, xMax, mlfitMin);
      std::string outputFileName_plot = (const char*)TString(outputFileName.c_str()).ReplaceAll(".root", Form("_mlfitScan_%s.png", parNames[parToScan].c_str()));
      showGraphs(1150, 950,
                 mlfitScan_fixed,    "Fixed",
                 mlfitScan_profiled, "Profiled",
                 0.040, 0.48, 0.78, 0.22, 0.14, 
	         xMin, xMax, parNames[parToScan], 1.2,
                 false, 0., 2.5e+1, "-2 log(L)", 1.4,
                 outputFileName_plot);
      std::cout << " Done.\n";
    }
  }

  math::Vector3 Bp;
  math::Vector3 BpErr;
  math::Vector3 Bm;
  math::Vector3 BmErr;
  math::Matrix3x3 C;
  math::Matrix3x3 CErr;
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

  Measurement nominal_measurement_by_summation = comp_Bp_Bm_C_by_summation(nominal_sample, verbosity);
 
  // initialize Minuit
  ROOT::Math::Minimizer* mlfit = new ROOT::Minuit2::Minuit2Minimizer();
  mlfit->SetMaxFunctionCalls(10000);
  mlfit->SetTolerance(1.e-3);
  mlfit->SetPrintLevel(-1);

  // set function pointer
  ROOT::Math::Functor f(&mlfit_fcn, npar);
  mlfit->SetFunction(f);

  Measurement nominal_measurement_by_mlfit = comp_Bp_Bm_C_by_mlfit(nominal_sample, mlfit, scanLikelihood, outputFile.file(), verbosity);

  // CV: estimate uncertainties on Bp and Bm vectors, on tau spin correlation matrix C,
  //     and on Entanglement observables with bootstrap samples
  std::vector<Measurement> bootstrap_measurements_by_summation;
  std::vector<Measurement> bootstrap_measurements_by_mlfit;
  for ( size_t idxBootstrapSample = 0; idxBootstrapSample < numBootstrapSamples; ++idxBootstrapSample )
  {
    EntanglementDataset bootstrap_sample = build_bootstrap_sample(dataset, rnd, maxEvents_afterCuts);
  
    Measurement measurement_by_summation = comp_Bp_Bm_C_by_summation(bootstrap_sample, -1);
    bootstrap_measurements_by_summation.push_back(measurement_by_summation);

    Measurement measurement_by_mlfit = comp_Bp_Bm_C_by_mlfit(bootstrap_sample, mlfit, false, "", -1);
    bootstrap_measurements_by_mlfit.push_back(measurement_by_mlfit);
  }

  math::Vector3 Bp_median_by_summation, BpErr_by_summation, Bm_median_by_summation, BmErr_by_summation;
  math::Matrix3x3 C_median_by_summation, CErr_by_summation;
  double Rchsh_median_by_summation, RchshErr_by_summation;
  comp_median_and_Err(bootstrap_measurements_by_summation, 
    Bp_median_by_summation, BpErr_by_summation,
    Bm_median_by_summation, BmErr_by_summation,
    C_median_by_summation, CErr_by_summation,
    Rchsh_median_by_summation, RchshErr_by_summation);
  std::cout << "Matrix C (measured using Eq. (25) of arXiv:2211.10513):\n";
  std::cout << nominal_measurement_by_summation.get_C() << "\n";
  std::cout << "+/-\n";
  std::cout << CErr_by_summation << "\n";
  std::cout << "Rchsh = " << nominal_measurement_by_summation.get_Rchsh() << " +/- " << RchshErr_by_summation << "\n";

  math::Vector3 Bp_median_by_mlfit, BpErr_by_mlfit, Bm_median_by_mlfit, BmErr_by_mlfit;
  math::Matrix3x3 C_median_by_mlfit, CErr_by_mlfit;
  double Rchsh_median_by_mlfit, RchshErr_by_mlfit;
  comp_median_and_Err(bootstrap_measurements_by_mlfit, 
    Bp_median_by_mlfit, BpErr_by_mlfit,
    Bm_median_by_mlfit, BmErr_by_mlfit,
    C_median_by_mlfit, CErr_by_mlfit,
    Rchsh_median_by_mlfit, RchshErr_by_mlfit);
  std::cout << "Matrix C (measured using maximum-likelihood fit):\n";
  std::cout << nominal_measurement_by_mlfit.get_C() << "\n";
  std::cout << "+/-\n";
  std::cout << CErr_by_mlfit << "\n";
  if ( verbosity >= 1 )
  {
    std::cout << "Uncertainty estimated by maximum-likelihood fit (for comparison):\n";
    std::cout << nominal_measurement_by_mlfit.get_CErr() << "\n";
  }
  std::cout << "Rchsh = " << nominal_measurement_by_mlfit.get_Rchsh() << " +/- " << RchshErr_by_mlfit << "\n";

  if ( verbosity >= 1 )
  {
    std::cout << "Standard Model expectation (given by Eq. (69) of arXiv:2208:11723):\n";
    TMatrixD C_exp(3,3);
    C_exp[0][0] = +1.;
    C_exp[1][1] = +1.;
    C_exp[2][2] = -1.;
    C_exp.Print();
  }

  delete mlfit;

  clock.Show("analyzeEntanglementNtuple");

  return EXIT_SUCCESS;
}
