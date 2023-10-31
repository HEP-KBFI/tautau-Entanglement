#ifndef TauAnalysis_Entanglement_comp_PCA_line2line_h
#define TauAnalysis_Entanglement_comp_PCA_line2line_h

#include "DataFormats/Candidate/interface/Candidate.h" // reco::Candidate::LorentzVector, reco::Candidate::Point

std::pair<reco::Candidate::Point, reco::Candidate::Point>
comp_PCA_line2line(const reco::Candidate::Point& P1,
                   const reco::Candidate::Vector& V1,
                   const reco::Candidate::Point& p2,
                   const reco::Candidate::Vector& V2,
                   int verbosity = -1);

reco::Candidate::Point
comp_PCA_line2line(const reco::Candidate::Point& pv,
                   const reco::Candidate::LorentzVector& tauP4,
                   const reco::Candidate::Point& sv,
                   const reco::Candidate::LorentzVector& visTauP4,
                   int verbosity = -1);

#endif // TauAnalysis_Entanglement_comp_PCA_line2line_h

