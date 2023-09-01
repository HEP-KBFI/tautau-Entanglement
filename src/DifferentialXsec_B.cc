#include "TauAnalysis/Entanglement/interface/DifferentialXsec_B.h"

#include "RooArgSet.h"     // RooArgSet
#include "RooGenericPdf.h" // RooGenericPdf
#include "RooRealVar.h"    // RooRealVar
#include "RooDataHist.h"   // RooDataHist

#include <TString.h>       // Form()

using namespace spin;

DifferentialXsec_B::DifferentialXsec_B(const std::string& label, int numBins)
  : label_(label)
  , histogram_(nullptr)
{
  std::string histogramName = Form("histogram1d_%s", label.c_str());
  histogram_ = new TH1D(histogramName.c_str(), histogramName.c_str(), numBins, -1., +1.);
  histogram_->Sumw2();
}

DifferentialXsec_B::~DifferentialXsec_B()
{
  delete histogram_;
}

void
DifferentialXsec_B::fill(double cosTheta, double evtWeight)
{
  histogram_->Fill(cosTheta, evtWeight);
}

const std::string&
DifferentialXsec_B::get_label() const
{
  return label_;
}

const TH1*
DifferentialXsec_B::get_histogram() const
{
  return histogram_;
}

namespace spin
{

double
fit_DifferentialXsec_B(const DifferentialXsec_B& Xsec, int verbosity)
{
  if ( verbosity >= 2 )
  {
    std::cout << "<fit_DifferentialXsec_B>: fitting " << Xsec.get_label() << "\n";
  }

  RooRealVar cosTheta("cosTheta", "cosTheta", -1., +1.);
    
  RooDataHist data("data", "data", RooArgSet(cosTheta, cosTheta), Xsec.get_histogram());

  // CV: need to restrict range of fit variable B_ii to the interval [-1,+1],
  //     to avoid that PDF becomes less than zero
  //    (which is "unphysical" behaviour for a differential cross section and causes warnings from RooFit)
  RooRealVar B_i("B_i", "B_i", 0., -1., +1.);

  RooGenericPdf pdf("pdf", "pdf", "0.50*(1. + B_i*cosTheta)", RooArgSet(B_i, cosTheta));    
  pdf.fitTo(data, RooFit::PrintLevel(-1));
  if ( verbosity >= 2 )
  {
    std::cout << " " << Xsec.get_label() << " = " << B_i.getVal() << " +/- " << B_i.getError() << "\n";
  }

  return B_i.getVal();
}

}
