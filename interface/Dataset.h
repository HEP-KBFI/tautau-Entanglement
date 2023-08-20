#ifndef TauAnalysis_Entanglement_Dataset_h
#define TauAnalysis_Entanglement_Dataset_h

#include "TauAnalysis/Entanglement/interface/Data.h" // spin::Data

#include <TObject.h>                                 // TObject

#include <vector>                                    // std::vector<>

namespace spin
{

class Dataset : public TObject
{
 public:
  Dataset();
  Dataset(const Dataset& dataset, int maxEvents_afterCuts = -1);
  ~Dataset();
  
  void
  push_back(const Data& entry);

  inline
  size_t
  size() const
  {
    return data_.size();
  }

  inline
  const Data&
  at(size_t idx) const
  {
    return data_[idx];
  }

  friend class SpinAnalyzer;

 private:
  std::vector<Data> data_;
  size_t numEntries_;
};

}

#endif // TauAnalysis_Entanglement_Dataset_h
