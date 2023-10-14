#include "TauAnalysis/Entanglement/interface/Dataset.h"

using namespace spin;

Dataset::Dataset()
  : TObject()
  , numEntries_(0)
{}

Dataset::Dataset(const Dataset& dataset, int maxEvents_afterCuts)
  : TObject(dataset)
{
  // CV: If maxEvents_afterCuts == 0, an empty dataset will be created.
  //     This is the intended behaviour and is used for building bootstrap samples. 
  if ( maxEvents_afterCuts != 0 )
  {
    if ( maxEvents_afterCuts > 0 && maxEvents_afterCuts < (int)dataset.data_.size() )
    {
      // Karl: is this intended behavior that we're not copying the first maxEvents_afterCuts elements?
      std::copy(dataset.data_.begin(), dataset.data_.end(), data_.begin());
    }
    else
    {
      data_ = dataset.data_;
    }
    numEntries_ = data_.size();
  }
}

Dataset::~Dataset()
{}

void
Dataset::push_back(const Data& entry)
{
  data_.push_back(entry);
  ++numEntries_;
}
