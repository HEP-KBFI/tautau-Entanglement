#include "TauAnalysis/Entanglement/interface/EntanglementDataset.h"

EntanglementDataset::EntanglementDataset()
{}

EntanglementDataset::EntanglementDataset(const EntanglementDataset& dataset, int maxEvents_afterCuts)
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
        const EntanglementData& entry = dataset.data_.at(idxEntry);
        data_.push_back(entry);
      }
    }
    else
    {
      data_ = dataset.data_;
    }
  }
}

EntanglementDataset::~EntanglementDataset()
{}

void
EntanglementDataset::push_back(const EntanglementData& entry)
{
  data_.push_back(entry);
}
