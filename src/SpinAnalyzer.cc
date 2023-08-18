#include "TauAnalysis/Entanglement/interface/SpinAnalyzer.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h"          // cmsException
#include "TauAnalysis/Entanglement/interface/Matrix_and_Vector.h"     // Matrix3x3, Vector3
#include "TauAnalysis/Entanglement/interface/SpinAlgo_by_mlfit.h"     // SpinAlgo_by_mlfit
#include "TauAnalysis/Entanglement/interface/SpinAlgo_by_summation.h" // SpinAlgo_by_summation

#include <TMath.h>                                                    // TMath::Nint()

#include <algorithm>                                                  // std::sort()
#include <string>                                                     // std::string
#include <utility>                                                    // std::make_pair(), std::pair<>

SpinAnalyzer::SpinAnalyzer(const edm::ParameterSet& cfg)
  : spinAnalyzer_(cfg.getParameter<std::string>("spinAnalyzer"))
  , algo_(nullptr)
  , numBootstrapSamples_(cfg.getParameter<unsigned>("numBootstrapSamples"))
  , maxEvents_afterCuts_(cfg.getParameter<int>("maxEvents_afterCuts"))
  , verbosity_(cfg.getUntrackedParameter<int>("verbosity"))
{
  if      ( spinAnalyzer_ == "by_mlfit"     ) algo_ = new SpinAlgo_by_mlfit(cfg);
  else if ( spinAnalyzer_ == "by_summation" ) algo_ = new SpinAlgo_by_summation(cfg);
  else throw cmsException("SpinAnalyzer", __LINE__)
    << "Invalid Configuration parameter 'spinAnalyzer' = " << spinAnalyzer_ << " !!\n";
}

SpinAnalyzer::~SpinAnalyzer()
{
  delete algo_;
}

namespace
{
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
                      double& Rchsh_median, double& RchshErr)
  {
    std::vector<math::Vector3> measuredBp;
    std::vector<math::Vector3> measuredBm;
    std::vector<math::Matrix3x3> measuredC;
    std::vector<double> measuredRchsh;
    for ( const spin::Measurement& measurement : measurements )
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
}

spin::Measurement
SpinAnalyzer::operator()(const EntanglementDataset& dataset)
{
  EntanglementDataset nominal_sample(dataset, maxEvents_afterCuts_);
  algo_->set_verbosity(verbosity_);
  spin::Measurement nominal_measurement = (*algo_)(nominal_sample);

  // CV: estimate uncertainties on Bp and Bm vectors, on tau spin correlation matrix C,
  //     and on Entanglement observables with bootstrap samples
  std::cout << "Generating bootstrap samples...\n";
  std::vector<spin::Measurement> bootstrap_measurements;
  for ( size_t idxBootstrapSample = 0; idxBootstrapSample < numBootstrapSamples_; ++idxBootstrapSample )
  {
    if ( idxBootstrapSample > 0 && (idxBootstrapSample % 100) == 0 )
    {
      std::cout << " Processing " << idxBootstrapSample << "th sample.\n";
    }
    EntanglementDataset bootstrap_sample = build_bootstrap_sample(dataset, rnd_, maxEvents_afterCuts_);
    algo_->set_verbosity(-1);
    spin::Measurement bootstrap_measurement = (*algo_)(bootstrap_sample);
    bootstrap_measurements.push_back(bootstrap_measurement);
  }
  std::cout << " Done.\n";
  math::Vector3 Bp_median, BpErr, Bm_median, BmErr;
  math::Matrix3x3 C_median, CErr;
  double Rchsh_median, RchshErr;
  comp_median_and_Err(bootstrap_measurements, 
    Bp_median, BpErr,
    Bm_median, BmErr,
    C_median, CErr,
    Rchsh_median, RchshErr);
  nominal_measurement.set_BpErr(BpErr);
  nominal_measurement.set_BmErr(BmErr);
  nominal_measurement.set_CErr(CErr);
  nominal_measurement.set_RchshErr(RchshErr);
  
  return nominal_measurement;
}
