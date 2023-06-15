#ifndef TauAnalysis_Entanglement_get_localCoordinateSystem_h
#define TauAnalysis_Entanglement_get_localCoordinateSystem_h

#include "DataFormats/Candidate/interface/Candidate.h"         // Candidate::LorentzVector, Candidate::Point, Candidate::Vector

#include "TauAnalysis/Entanglement/interface/KinematicEvent.h" // KinematicEvent

#include <Math/Boost.h>                                        // Boost

enum { kBeam, kHiggs };

void
get_localCoordinateSystem(const reco::Candidate::LorentzVector& p4,
                          const reco::Candidate::LorentzVector* recoilP4, const ROOT::Math::Boost* boost_ttrf,
                          int hAxis,
                          reco::Candidate::Vector& r, reco::Candidate::Vector& n, reco::Candidate::Vector& k,
                          int verbosity = 0);

#endif // TauAnalysis_Entanglement_get_localCoordinateSystem_h
