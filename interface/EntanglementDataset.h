#ifndef TauAnalysis_Entanglement_EntanglementDataset_h
#define TauAnalysis_Entanglement_EntanglementDataset_h

#include "TauAnalysis/Entanglement/interface/EntanglementData.h" // EntanglementData

#include <TObject.h>                                             // TObject

#include <vector>                                                // std::vector<>

class EntanglementDataset : public TObject
{
 public:
  EntanglementDataset();
  EntanglementDataset(const EntanglementDataset& dataset, int maxEvents_afterCuts = -1);
  ~EntanglementDataset();
  
  void
  push_back(const EntanglementData& entry);

  inline
  size_t
  size() const
  {
    return data_.size();
  }

  inline
  const EntanglementData&
  at(size_t idx) const
  {
    return data_[idx];
  }

 private:
  std::vector<EntanglementData> data_;
};

#endif // TauAnalysis_Entanglement_EntanglementDataset_h
