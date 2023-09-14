#include "TauAnalysis/Entanglement/interface/SpinAnalyzer.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h"                   // cmsException
#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h"              // Matrix3x3, Vector3
#include "TauAnalysis/Entanglement/interface/SpinAlgo_by_asymmetry.h"          // SpinAlgo_by_asymmetry
#include "TauAnalysis/Entanglement/interface/SpinAlgo_by_differentialXsec1d.h" // SpinAlgo_by_differentialXsec1d
#include "TauAnalysis/Entanglement/interface/SpinAlgo_by_differentialXsec2d.h" // SpinAlgo_by_differentialXsec2d
#include "TauAnalysis/Entanglement/interface/SpinAlgo_by_mlfit.h"              // SpinAlgo_by_mlfit
#include "TauAnalysis/Entanglement/interface/SpinAlgo_by_summation.h"          // SpinAlgo_by_summation

#include <TMath.h>                                                             // TMath::Nint()

#include <algorithm>                                                           // std::sort()
#include <map>                                                                 // std::map<>
#include <string>                                                              // std::string
#include <utility>                                                             // std::make_pair(), std::pair<>

using namespace spin;

SpinAnalyzer::SpinAnalyzer(const edm::ParameterSet& cfg)
  : spinAnalyzer_(cfg.getParameter<std::string>("spinAnalyzer"))
  , algo_(nullptr)
  , numBootstrapSamples_(cfg.getParameter<unsigned>("numBootstrapSamples"))
  , maxEvents_afterCuts_(cfg.getParameter<int>("maxEvents_afterCuts"))
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
{
  if      ( spinAnalyzer_ == "by_asymmetry"          ) algo_ = new SpinAlgo_by_asymmetry(cfg);
  else if ( spinAnalyzer_ == "by_differentialXsec1d" ) algo_ = new SpinAlgo_by_differentialXsec1d(cfg);
  else if ( spinAnalyzer_ == "by_differentialXsec2d" ) algo_ = new SpinAlgo_by_differentialXsec2d(cfg);
  else if ( spinAnalyzer_ == "by_mlfit"              ) algo_ = new SpinAlgo_by_mlfit(cfg);
  else if ( spinAnalyzer_ == "by_summation"          ) algo_ = new SpinAlgo_by_summation(cfg);
  else throw cmsException("SpinAnalyzer", __LINE__)
    << "Invalid Configuration parameter 'spinAnalyzer' = " << spinAnalyzer_ << " !!\n";
}

SpinAnalyzer::~SpinAnalyzer()
{
  delete algo_;
}

void
SpinAnalyzer::set_verbosity(int verbosity)
{
  verbosity_ = verbosity;
  algo_->set_verbosity(verbosity);
}

namespace
{
  spin::Dataset
  build_bootstrap_sample(const spin::Dataset& dataset, TRandom& rnd, int maxEvents_afterCuts = -1)
  {
    spin::Dataset bootstrap_sample(dataset, 0);
    size_t sampleSize = ( maxEvents_afterCuts > 0 ) ? maxEvents_afterCuts : dataset.size();
    for ( size_t idxSample = 0; idxSample < sampleSize; ++idxSample )
    {
      size_t idxEntry = rnd.Integer(dataset.size());
      const spin::Data& entry = dataset.at(idxEntry);
      bootstrap_sample.push_back(entry);
    }
    return bootstrap_sample;
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
    std::map<int, std::vector<double>> tmp;
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
    std::map<int, std::map<int, std::vector<double>>> tmp;
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
  comp_median_and_Err(const std::vector<spin::Measurement>& measurements,
                      math::Vector3& Bp_median, math::Vector3& BpErr,
                      math::Vector3& Bm_median, math::Vector3& BmErr,
                      math::Matrix3x3& C_median, math::Matrix3x3& CErr,
                      double& concurrence_median, double& concurrenceErr,
                      double& Ek_median, double& EkErr,
                      double& Rchsh_median, double& RchshErr,
                      double& steerability_median, double& steerabilityErr)
  {
    std::vector<math::Vector3> measuredBp;
    std::vector<math::Vector3> measuredBm;
    std::vector<math::Matrix3x3> measuredC;
    std::vector<double> measured_concurrence;
    std::vector<double> measuredEk;
    std::vector<double> measuredRchsh;
    std::vector<double> measured_steerability;
    for ( const spin::Measurement& measurement : measurements )
    {
      measuredBp.push_back(measurement.get_Bp());
      measuredBm.push_back(measurement.get_Bm());
      measuredC.push_back(measurement.get_C());
      measured_concurrence.push_back(measurement.get_concurrence());
      measuredEk.push_back(measurement.get_Ek());
      measuredRchsh.push_back(measurement.get_Rchsh());
      measured_steerability.push_back(measurement.get_steerability());
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
    std::pair<double, double> concurrence_median_and_Err = comp_median_and_Err(measured_concurrence);
    concurrence_median = concurrence_median_and_Err.first;
    concurrenceErr = concurrence_median_and_Err.second;
    std::pair<double, double> Ek_median_and_Err = comp_median_and_Err(measuredEk);
    Ek_median = Ek_median_and_Err.first;
    EkErr = Ek_median_and_Err.second;
    std::pair<double, double> Rchsh_median_and_Err = comp_median_and_Err(measuredRchsh);
    Rchsh_median = Rchsh_median_and_Err.first;
    RchshErr = Rchsh_median_and_Err.second;
    std::pair<double, double> steerability_median_and_Err = comp_median_and_Err(measured_steerability);
    steerability_median = steerability_median_and_Err.first;
    steerabilityErr = steerability_median_and_Err.second;
  }
}

spin::Measurement
SpinAnalyzer::build_measurement(const spin::Dataset& dataset, int maxEvents_afterCuts, int verbosity) const
{
  spin::Dataset nominal_sample(dataset, maxEvents_afterCuts);
  algo_->set_verbosity(verbosity);
  spin::Measurement nominal_measurement = (*algo_)(nominal_sample);

  // CV: estimate uncertainties on Bp and Bm vectors, on tau spin correlation matrix C,
  //     and on Entanglement observables with bootstrap samples
  if ( verbosity >= 1 )
  {
    std::cout << "Generating bootstrap samples...\n";
  }
  std::vector<spin::Measurement> bootstrap_measurements;
  for ( size_t idxBootstrapSample = 0; idxBootstrapSample < numBootstrapSamples_; ++idxBootstrapSample )
  {
    if ( verbosity >= 1 )
    {
      if ( idxBootstrapSample > 0 && (idxBootstrapSample % 100) == 0 )
      {
        std::cout << " Processing " << idxBootstrapSample << "th sample\n";
      }
    }
    spin::Dataset bootstrap_sample = build_bootstrap_sample(dataset, rnd_, maxEvents_afterCuts);
    algo_->set_verbosity(-1);
    spin::Measurement bootstrap_measurement = (*algo_)(bootstrap_sample);
    bootstrap_measurements.push_back(bootstrap_measurement);
  }
  if ( verbosity >= 1 )
  {
    std::cout << " Done.\n";
  }
  math::Vector3 Bp_median, BpErr, Bm_median, BmErr;
  math::Matrix3x3 C_median, CErr;
  double concurrence_median, concurrenceErr, Ek_median, EkErr, Rchsh_median, RchshErr, steerability_median, steerabilityErr;
  comp_median_and_Err(bootstrap_measurements, 
    Bp_median, BpErr,
    Bm_median, BmErr,
    C_median, CErr,
    concurrence_median, concurrenceErr,
    Ek_median, EkErr,
    Rchsh_median, RchshErr,
    steerability_median, steerabilityErr);
  nominal_measurement.set_BpErr(BpErr);
  nominal_measurement.set_BmErr(BmErr);
  nominal_measurement.set_CErr(CErr);
  nominal_measurement.set_concurrenceErr(concurrenceErr);
  nominal_measurement.set_EkErr(EkErr);
  nominal_measurement.set_RchshErr(RchshErr);
  nominal_measurement.set_steerabilityErr(steerabilityErr);
  
  return nominal_measurement;
}

spin::Measurement
SpinAnalyzer::operator()(const spin::Dataset& dataset) const
{
  return build_measurement(dataset, maxEvents_afterCuts_, verbosity_);
}

spin::BinnedMeasurement1d
SpinAnalyzer::operator()(const spin::BinnedDataset1d& dataset) const
{
  if ( verbosity_ >= 2 )
  {
    std::cout << "<SpinAnalyzer::operator()>:\n";
  }
  spin::BinnedMeasurement1d measurement(dataset.name_, dataset.numBinsX_, dataset.xMin_, dataset.xMax_);
  for ( int idxBinX = 1; idxBinX <= dataset.binning_->GetNbinsX(); ++idxBinX )
  {
    measurement.binning_->SetBinContent(idxBinX, dataset.binning_->GetBinContent(idxBinX));
    measurement.binning_->SetBinError(idxBinX, dataset.binning_->GetBinError(idxBinX));
  }
  int numEntries_unallocated = dataset.numEntries_;
  double prob_unallocated = 1.;
  for ( const std::pair<int, spin::Dataset>& entry : dataset.entries_ )
  {    
    double prob = entry.second.numEntries_/(double)dataset.numEntries_;
    int maxEvents = 0;
    if ( prob < prob_unallocated )
    {
      maxEvents = rnd_.Binomial(numEntries_unallocated, prob/prob_unallocated);
      numEntries_unallocated -= maxEvents;
      prob_unallocated -= prob;
    }
    else
    {
      maxEvents = numEntries_unallocated;
      numEntries_unallocated = 0;
      prob_unallocated = 0.;
    }
    if ( maxEvents > 0 ) measurement.entries_[entry.first] = build_measurement(entry.second, maxEvents, -1);
    else                 measurement.entries_[entry.first] = Measurement();
    if ( verbosity_ >= 2 )
    {
      int idxBinX, idxBinY, idxBinZ;
      dataset.binning_->GetBinXYZ(entry.first, idxBinX, idxBinY, idxBinZ);
      double x = dataset.binning_->GetXaxis()->GetBinCenter(idxBinX);
      std::cout << "x = " << x << " (#entries = " << dataset.entries_.find(entry.first)->second.numEntries_ << "):"
                << " Rchsh = " << measurement.entries_[entry.first].get_Rchsh() 
                << " +/- " << measurement.entries_[entry.first].get_RchshErr() << "\n";
      std::cout << "C:\n";
      std::cout << measurement.entries_[entry.first].get_C() << "\n";
    }
  }
  return measurement;
}

spin::BinnedMeasurement2d
SpinAnalyzer::operator()(const spin::BinnedDataset2d& dataset) const
{
  if ( verbosity_ >= 2 )
  {
    std::cout << "<SpinAnalyzer::operator()>:\n";
  }
  spin::BinnedMeasurement2d measurement(dataset.name_, dataset.numBinsX_, dataset.xMin_, dataset.xMax_, dataset.numBinsY_, dataset.yMin_, dataset.yMax_);
  for ( int idxBinX = 1; idxBinX <= dataset.binning_->GetNbinsX(); ++idxBinX )
  {
    for ( int idxBinY = 1; idxBinY <= dataset.binning_->GetNbinsY(); ++idxBinY )
    {
      measurement.binning_->SetBinContent(idxBinX, idxBinY, dataset.binning_->GetBinContent(idxBinX, idxBinY));
      measurement.binning_->SetBinError(idxBinX, idxBinY, dataset.binning_->GetBinError(idxBinX, idxBinY));
    }
  }
  int numEntries_unallocated = dataset.numEntries_;
  double prob_unallocated = 1.;
  for ( const std::pair<int, spin::Dataset>& entry : dataset.entries_ )
  {
    double prob = entry.second.numEntries_/(double)dataset.numEntries_;
    int maxEvents = 0;
    if ( prob < prob_unallocated )
    {
      maxEvents = rnd_.Binomial(numEntries_unallocated, prob/prob_unallocated);
      numEntries_unallocated -= maxEvents;
      prob_unallocated -= prob;
    }
    else
    {
      maxEvents = numEntries_unallocated;
      numEntries_unallocated = 0;
      prob_unallocated = 0.;
    }
    if ( maxEvents > 0 ) measurement.entries_[entry.first] = build_measurement(entry.second, maxEvents, -1);
    else                 measurement.entries_[entry.first] = Measurement();
    if ( verbosity_ >= 2 )
    {
      int idxBinX, idxBinY, idxBinZ;
      dataset.binning_->GetBinXYZ(entry.first, idxBinX, idxBinY, idxBinZ);
      double x = dataset.binning_->GetXaxis()->GetBinCenter(idxBinX);
      double y = dataset.binning_->GetYaxis()->GetBinCenter(idxBinY);
      std::cout << "x = " << x << ", y = " << y << " (#entries = " << dataset.entries_.find(entry.first)->second.numEntries_ << "):" 
                << " Rchsh = " << measurement.entries_[entry.first].get_Rchsh() 
                << " +/- " << measurement.entries_[entry.first].get_RchshErr() << "\n";
      std::cout << "C:\n";
      std::cout << measurement.entries_[entry.first].get_C() << "\n";
    }
  }
  return measurement;
}
