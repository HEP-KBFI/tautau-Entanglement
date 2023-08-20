#include "TauAnalysis/Entanglement/interface/Dataset.h"

using namespace spin;

Dataset::Dataset()
  : numEntries_(0)
{}

Dataset::Dataset(const Dataset& dataset, int maxEvents_afterCuts)
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
        const Data& entry = dataset.data_.at(idxEntry);
        data_.push_back(entry);
      }
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
