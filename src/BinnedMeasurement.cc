#include "TauAnalysis/Entanglement/interface/BinnedMeasurement.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h" // cmsException

#include <TString.h>                                         // Form()

#include <utility>                                           // std::pair<>

using namespace spin;

BinnedMeasurement::BinnedMeasurement(const std::string& name)
  : name_(name)
{}

BinnedMeasurement::~BinnedMeasurement()
{}

double
BinnedMeasurement::get_value(const Measurement& entry, const std::string& observable)
{
  if      ( observable == "C_rr"  ) return entry.get_C()(0,0);
  else if ( observable == "C_nn"  ) return entry.get_C()(1,1);
  else if ( observable == "C_kk"  ) return entry.get_C()(2,2);
  else if ( observable == "Rchsh" ) return entry.get_Rchsh();
  else throw cmsException("BinnedMeasurement::get_value", __LINE__)
    << "Invalid function argument 'observable' = " << observable << " !!\n";
}

double
BinnedMeasurement::get_Err(const Measurement& entry, const std::string& observable)
{
  if      ( observable == "C_rr"  ) return entry.get_CErr()(0,0);
  else if ( observable == "C_nn"  ) return entry.get_CErr()(1,1);
  else if ( observable == "C_kk"  ) return entry.get_CErr()(2,2);
  else if ( observable == "Rchsh" ) return entry.get_RchshErr();
  else throw cmsException("BinnedMeasurement::get_Err", __LINE__)
    << "Invalid function argument 'observable' = " << observable << " !!\n";
}

BinnedMeasurement1d::BinnedMeasurement1d(const std::string& name, int numBinsX, double xMin, double xMax)
  : BinnedMeasurement(name)
  , numBinsX_(numBinsX)
  , xMin_(xMin)
  , xMax_(xMax)
  , binning_(nullptr)
{
  binning_ = new TH1D(name.c_str(), name.c_str(), numBinsX, xMin, xMax);
  binning_->Sumw2();
}

BinnedMeasurement1d::~BinnedMeasurement1d()
{
  delete binning_;
}

TH1*
BinnedMeasurement1d::get_histogram(const std::string& observable)
{
  std::string histogramName = Form("%s_%s", name_.c_str(), observable.c_str());
  TH1* histogram = (TH1*)binning_->Clone(histogramName.c_str());
  histogram->Reset();
  for ( const std::pair<int, Measurement>& entry : entries_ )
  {
    histogram->SetBinContent(entry.first, get_value(entry.second, observable));
    histogram->SetBinError(entry.first, get_Err(entry.second, observable));
  }
  return histogram;
}

BinnedMeasurement2d::BinnedMeasurement2d(const std::string& name, int numBinsX, double xMin, double xMax, int numBinsY, double yMin, double yMax)
  : BinnedMeasurement(name)
  , numBinsX_(numBinsX)
  , xMin_(xMin)
  , xMax_(xMax)
  , numBinsY_(numBinsY)
  , yMin_(yMin)
  , yMax_(yMax)
  , binning_(nullptr)
{
  binning_ = new TH2D(name.c_str(), name.c_str(), numBinsX, xMin, xMax, numBinsY, yMin, yMax);
  binning_->Sumw2();
}

BinnedMeasurement2d::~BinnedMeasurement2d()
{
  delete binning_;
}

TH2*
BinnedMeasurement2d::get_histogram(const std::string& observable)
{
  std::string histogramName = Form("%s_%s", name_.c_str(), observable.c_str());
  TH2* histogram = (TH2*)binning_->Clone(histogramName.c_str());
  histogram->Reset();
  for ( const std::pair<int, Measurement>& entry : entries_ )
  {
    histogram->SetBinContent(entry.first, get_value(entry.second, observable));
    histogram->SetBinError(entry.first, get_Err(entry.second, observable));
  }
  return histogram;
}
