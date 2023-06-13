#ifndef TauAnalysis_Entanglement_Resolutions_h
#define TauAnalysis_Entanglement_Resolutions_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"        // edm::ParameterSet

class Resolutions
{
 public:
  Resolutions(const edm::ParameterSet& cfg);
  ~Resolutions();

  double
  get_recoilResolutionPx() const;

  double
  get_recoilResolutionPy() const;
  
  double
  get_recoilResolutionPz() const;

  double
  get_recoilResolutionE() const;

  double
  get_pvResolutionXY() const;
  
  double
  get_pvResolutionZ() const;
  
  double
  get_svResolutionParl() const;
  
  double
  get_svResolutionPerp() const;
  
  double
  get_tipResolutionPerp() const;

 private:
  double recoilResolutionPx_; // [GeV]
  double recoilResolutionPy_; // [GeV]
  double recoilResolutionPz_; // [GeV]
  double recoilResolutionE_;  // [GeV]
  double pvResolutionXY_;     // [cm]
  double pvResolutionZ_;      // [cm]
  double svResolutionParl_;   // [cm]
  double svResolutionPerp_;   // [cm]
  double tipResolutionPerp_;  // [cm]
};

#endif // TauAnalysis_Entanglement_Resolutions_h
