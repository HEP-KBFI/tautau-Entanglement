#include "TauAnalysis/Entanglement/interface/SpinAlgo_by_mlfit.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h" // cmsException
#include "TauAnalysis/Entanglement/interface/comp_Rchsh.h"   // comp_Rchsh()

#include <Math/Functor.h>                                    // ROOT::Math::Functor
#include <TCanvas.h>                                         // TCanvas
#include <TGraph.h>                                          // TGraph
#include <TH1.h>                                             // TH1F
#include <TLegend.h>                                         // TLegend
#include <TString.h>                                         // Form()

#include <assert.h>                                          // assert()
#include <map>                                               // std::map<>

const size_t npar = 15;

static const spin::Dataset* gDataset = nullptr;
static const std::vector<double>* gpar_gen = nullptr;

using namespace spin;

SpinAlgo_by_mlfit::SpinAlgo_by_mlfit(const edm::ParameterSet& cfg)
  : SpinAlgoBase(cfg)
  , mlfit_(nullptr)
  , par_gen_(cfg.getParameter<std::vector<double>>("par_gen"))
  , scanLikelihood_(cfg.getParameter<bool>("mlfit_scan_likelihood"))
  , outputFileName_(cfg.getParameter<std::string>("mlfit_outputFileName"))
  , algo_by_summation_(new SpinAlgo_by_summation(cfg))
{
  // read elements of Bp and Bm vector and of tau spin correlation matrix C,
  // that were used to produce the analyzed Monte Carlo sample.
  // The values are used for the normalization of the likelihood function,
  // to avoid biases arising from the event selection.
  if ( par_gen_.size() != npar )
    throw cmsException("SpinAlgo_by_mlfit", __LINE__) 
      << "Invalid Configuration parameter 'par_gen' !!\n";

  // initialize Minuit
  mlfit_ = new ROOT::Minuit2::Minuit2Minimizer();
 
  mlfit_->SetMaxFunctionCalls(10000);
  mlfit_->SetTolerance(1.e-3);
  mlfit_->SetPrintLevel(-1);
}

SpinAlgo_by_mlfit::~SpinAlgo_by_mlfit()
{
  delete algo_by_summation_;
  delete mlfit_;
}

namespace
{
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

  double
  get_p(const double* par, const spin::Data& entry)
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
    assert(gDataset);
    const spin::Dataset* mlfitData = gDataset;

    size_t numEntries = mlfitData->size();

    assert(gpar_gen);
    double par_gen[npar];
    for ( size_t idxPar = 0; idxPar < npar; ++idxPar )
    {
      par_gen[idxPar] = gpar_gen->at(idxPar);
    }

    double norm = 0.;
    for ( size_t idxEntry = 0; idxEntry < numEntries; ++idxEntry )
    {
      const spin::Data& entry = mlfitData->at(idxEntry);

      norm += entry.get_evtWeight()*get_p(par, entry)/get_p(par_gen, entry);
      //norm += get_p(par, entry)/get_p(par_gen, entry);
    }

    double logL = 0.;
    for ( size_t idxEntry = 0; idxEntry < numEntries; ++idxEntry )
    {
      const spin::Data& entry = mlfitData->at(idxEntry);

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
}

spin::Measurement
SpinAlgo_by_mlfit::operator()(const spin::Dataset& dataset)
{
  // initialize fit parameters
  spin::Measurement startpos = (*algo_by_summation_)(dataset);
  const math::Matrix3x3& startpos_C = startpos.get_C();
  std::map<size_t, std::string> parNames = get_parNames();
  assert(parNames.size() == npar);
  for ( size_t idxPar = 0; idxPar < npar; ++idxPar )
  {
    double par0 = 0.;
    if ( idxPar >= 6 && idxPar <= 14 )
    {
      size_t idxRow = (idxPar - 6) / 3;
      size_t idxCol = (idxPar - 6) % 3;
      par0 = startpos_C(idxRow,idxCol);
    }
    mlfit_->SetLimitedVariable(idxPar, parNames[idxPar].c_str(), par0, 0.1, -2., +2.);
  }

  // set function pointer
  ROOT::Math::Functor f(&mlfit_fcn, npar);
  mlfit_->SetFunction(f);

  gDataset = &dataset;
  gpar_gen = &par_gen_;

  mlfit_->Minimize();
  mlfit_->Hesse();

  if ( verbosity_ >= 1 )
  {
    std::cout << "Fit Results:\n";
    mlfit_->PrintResults();
  }

  std::vector<double> parValues(npar);
  std::vector<double> parErrors(npar);
  for ( size_t idxPar = 0; idxPar < npar; ++idxPar )
  {
    const double* X = mlfit_->X();
    parValues[idxPar] = X[idxPar];
    parErrors[idxPar] = sqrt(mlfit_->CovMatrix(idxPar, idxPar));
  }

  if ( scanLikelihood_ )
  {
    double mlfitMin = mlfit_->MinValue();
    std::map<size_t, std::string> parNames = get_parNames();
    std::vector<size_t> parsToScan = { 6, 10, 14 };
    for ( size_t i = 0; i < parsToScan.size(); ++i )
    {
      size_t parToScan = parsToScan[i];      
      std::cout << "Scanning likelihood as function of parameter " << parNames[parToScan] << "...\n";
      double xMin = parValues[parToScan] - std::max(0.01, 5.*parErrors[parToScan]);
      double xMax = parValues[parToScan] + std::max(0.01, 5.*parErrors[parToScan]);
      TGraph* mlfitScan_fixed = scan_mlfit_fixed(parValues, parToScan, parNames[parToScan], 100, xMin, xMax, mlfitMin);
      TGraph* mlfitScan_profiled = scan_mlfit_profiled(mlfit_, parValues, parToScan, parNames[parToScan], 100, xMin, xMax, mlfitMin);
      std::string mlfitScan_outputFileName = (const char*)TString(outputFileName_.c_str()).ReplaceAll(".root", Form("_mlfitScan_%s.png", parNames[parToScan].c_str()));
      showGraphs(1150, 950,
                 mlfitScan_fixed,    "Fixed",
                 mlfitScan_profiled, "Profiled",
                 0.040, 0.48, 0.78, 0.22, 0.14, 
	         xMin, xMax, parNames[parToScan], 1.2,
                 false, 0., 2.5e+1, "-2 log(L)", 1.4,
                 mlfitScan_outputFileName);
      std::cout << " Done.\n";
    }
  }

  math::Vector3 Bp, Bm;
  math::Matrix3x3 C;
  for ( size_t idxPar = 0; idxPar < npar; ++idxPar )
  {
    if ( idxPar <= 2 )
    {
      size_t idx = idxPar;
      Bp[idx] = parValues[idxPar];
    }
    else if ( idxPar <= 5 )
    {
      size_t idx = idxPar - 3;
      Bm[idx] = parValues[idxPar];
    }
    else if ( idxPar <= 14 )
    {	
      size_t idxRow = (idxPar - 6) / 3;
      size_t idxCol = (idxPar - 6) % 3;
      C[idxRow][idxCol] = parValues[idxPar];
    } else assert(0);
  }

  double Rchsh = comp_Rchsh(C, verbosity_);

  spin::Measurement measurement(Bp, Bm, C, Rchsh);
  return measurement;
}
