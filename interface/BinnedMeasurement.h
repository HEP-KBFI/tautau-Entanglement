#ifndef TauAnalysis_Entanglement_BinnedMeasurement_h
#define TauAnalysis_Entanglement_BinnedMeasurement_h

#include "TauAnalysis/Entanglement/interface/Measurement.h" // spin::Measurement

#include <TH1.h>                                            // TH1D
#include <TH2.h>                                            // TH2D

#include <map>                                              // std::map<>

namespace spin
{

class BinnedMeasurement
{
 public:
  BinnedMeasurement(const std::string& name);
  ~BinnedMeasurement();

  friend class SpinAnalyzer;

 protected:
  double
  get_value(const Measurement& entry, const std::string& observable);
  double
  get_Err(const Measurement& entry, const std::string& observable);

  std::string name_;

  std::map<int, Measurement> entries_;
};

class BinnedMeasurement1d : public BinnedMeasurement
{
 public:
  BinnedMeasurement1d(const std::string& name, int numBinsX, double xMin, double xMax);
  ~BinnedMeasurement1d();

  TH1*
  get_histogram(const std::string& observable);

  friend class SpinAnalyzer;

 protected:
  int numBinsX_;
  double xMin_;
  double xMax_;
  TH1* binning_;
};

class BinnedMeasurement2d : public BinnedMeasurement
{
 public:
  BinnedMeasurement2d(const std::string& name, int numBinsX, double xMin, double xMax, int numBinsY, double yMin, double yMax);
  ~BinnedMeasurement2d();

  TH2*
  get_histogram(const std::string& observable);

  friend class SpinAnalyzer;

 protected:
  int numBinsX_;
  double xMin_;
  double xMax_;
  int numBinsY_;
  double yMin_;
  double yMax_;
  TH2* binning_;
};

}

#endif // TauAnalysis_Entanglement_BinnedMeasurement_h
