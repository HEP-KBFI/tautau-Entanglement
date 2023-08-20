#include "TauAnalysis/Entanglement/interface/BinnedDataset.h"

using namespace spin;

BinnedDataset::BinnedDataset()
{}

BinnedDataset::BinnedDataset(const std::string& name)
  : name_(name)
  , numEntries_(0)
{}
  
BinnedDataset::~BinnedDataset()
{}

void
BinnedDataset::push_back(int idxBin, const Data& entry)
{
  entries_[idxBin].push_back(entry);
  ++numEntries_;
}

BinnedDataset1d::BinnedDataset1d(const std::string& name, int numBinsX, double xMin, double xMax)
  : BinnedDataset(name)
  , numBinsX_(numBinsX)
  , xMin_(xMin)
  , xMax_(xMax)
  , binning_(nullptr)
{
  binning_ = new TH1D(name.c_str(), name.c_str(), numBinsX, xMin, xMax);
  binning_->Sumw2();
}

BinnedDataset1d::~BinnedDataset1d()
{
  delete binning_;
}

void
BinnedDataset1d::push_back(double x, const Data& entry)
{
  binning_->Fill(x, entry.get_evtWeight());
  int idxBin = binning_->FindBin(x);
  BinnedDataset::push_back(idxBin, entry);
}

BinnedDataset2d::BinnedDataset2d(const std::string& name, int numBinsX, double xMin, double xMax, int numBinsY, double yMin, double yMax)
  : BinnedDataset(name)
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

BinnedDataset2d::~BinnedDataset2d()
{
  delete binning_;
}

void
BinnedDataset2d::push_back(double x, double y, const Data& entry)
{
  binning_->Fill(x, y, entry.get_evtWeight());
  int idxBin = binning_->FindBin(x, y);
  BinnedDataset::push_back(idxBin, entry);
}
