#ifndef TauAnalysis_Entanglement_BinnedDataset_h
#define TauAnalysis_Entanglement_BinnedDataset_h

#include "TauAnalysis/Entanglement/interface/Data.h"    // spin::Data
#include "TauAnalysis/Entanglement/interface/Dataset.h" // spin::Dataset

#include <TH1.h>                                        // TH1D
#include <TH2.h>                                        // TH2D

#include <map>                                          // std::map<>
#include <vector>                                       // std::vector<>

namespace spin
{

class BinnedDataset
{
 public:
  BinnedDataset();
  BinnedDataset(const std::string& name);
  ~BinnedDataset();

  friend class SpinAnalyzer;

 protected:
  void
  push_back(int idxBin, const Data& entry);

  std::string name_;

  std::map<int, Dataset> entries_;
  size_t numEntries_;
};

class BinnedDataset1d : public BinnedDataset
{
 public:
  BinnedDataset1d(const std::string& name, int numBinsX, double xMin, double xMax);
  ~BinnedDataset1d();

  void
  push_back(double x, const Data& entry);

  friend class SpinAnalyzer;

 private:
  int numBinsX_;
  double xMin_;
  double xMax_;
  TH1* binning_;
};

class BinnedDataset2d : public BinnedDataset
{
 public:
  BinnedDataset2d(const std::string& name, int numBinsX, double xMin, double xMax, int numBinsY, double yMin, double yMax);
  ~BinnedDataset2d();

  void
  push_back(double x, double y, const Data& entry);

  friend class SpinAnalyzer;

 private:
  int numBinsX_;
  double xMin_;
  double xMax_;
  int numBinsY_;
  double yMin_;
  double yMax_;
  TH2* binning_;
};

}

#endif // TauAnalysis_Entanglement_BinnedDataset_h
