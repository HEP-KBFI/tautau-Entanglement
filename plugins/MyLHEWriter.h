#ifndef TauAnalysis_Entanglement_MyLHEWriter_h
#define TauAnalysis_Entanglement_MyLHEWriter_h

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Candidate/interface/Candidate.h"                    // reco::Candidate
#include "DataFormats/Candidate/interface/CandidateFwd.h"                 // reco::CandidateView
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h" // LHERunInfoProduct
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"   // LHEEventProduct

#include <fstream>                                                        // std::ofstream
#include <string>                                                         // std::string

class MyLHEWriter : public edm::one::EDAnalyzer<edm::one::WatchRuns> 
{
 public:
  explicit MyLHEWriter(const edm::ParameterSet& cfg);
  ~MyLHEWriter() override;
 
 protected:
  void beginRun(const edm::Run& run, const edm::EventSetup& es) override;
  void endRun(const edm::Run& run, const edm::EventSetup& es) override;
  void analyze(const edm::Event& event, const edm::EventSetup& es) override;

 private:
  std::ofstream file_;
  std::ofstream file1_;

  edm::InputTag srcLHERunInfo_;
  edm::EDGetTokenT<LHERunInfoProduct> tokenLHERunInfo_;

  edm::InputTag srcLHEEvent_;
  edm::EDGetTokenT<LHEEventProduct> tokenLHEEvent_;

  edm::InputTag srcGenParticles_;
  edm::EDGetTokenT<reco::CandidateView> tokenGenParticles_;

  std::string outputFileName_;
  std::string outputFileName1_;
};

#endif // TauAnalysis_Entanglement_MyLHEWriter_h


