#include "TauAnalysis/Entanglement/plugins/MyLHEWriter.h"

#include <TString.h>                                      // Form(), ReplaceAll(), TString

#include <algorithm>                                      // std::find()
#include <iomanip>                                        // std::setw()
#include <memory>                                         // std::copy()
#include <vector>                                         // std::vector<>

MyLHEWriter::MyLHEWriter(const edm::ParameterSet& cfg)
{
  srcLHERunInfo_ = cfg.getParameter<edm::InputTag>("srcLHERunInfo");
  tokenLHERunInfo_ = consumes<LHERunInfoProduct, edm::InRun>(srcLHERunInfo_);

  srcLHEEvent_ = cfg.getParameter<edm::InputTag>("srcLHEEvent");
  tokenLHEEvent_ = consumes<LHEEventProduct>(srcLHEEvent_);

  srcGenParticles_ = cfg.getParameter<edm::InputTag>("srcGenParticles");
  tokenGenParticles_ = consumes<reco::CandidateView>(srcGenParticles_);

  outputFileName_ = cfg.getParameter<std::string>("outputFileName");
  outputFileName1_ = TString(outputFileName_.c_str()).ReplaceAll(".", "1.");
}

MyLHEWriter::~MyLHEWriter()
{}

void
MyLHEWriter::beginRun(const edm::Run& run, const edm::EventSetup& es)
{
  file_.open(outputFileName_, std::fstream::out | std::fstream::trunc);
  file1_.open(outputFileName1_, std::fstream::out | std::fstream::trunc);
}

void
MyLHEWriter::endRun(const edm::Run& run, const edm::EventSetup& es)
{
  edm::Handle<LHERunInfoProduct> lheRunInfo = run.getHandle(tokenLHERunInfo_);

  std::copy(lheRunInfo->begin(), lheRunInfo->end(), std::ostream_iterator<std::string>(file_));

  file1_ << LHERunInfoProduct::endOfFile();
  file1_.close();

  file_.close();

  system(Form("cat %s >> %s", outputFileName1_.c_str(), outputFileName_.c_str()));
  system(Form("rm -rf %s", outputFileName1_.c_str()));
}

namespace
{
  int
  getIdx(const std::vector<const reco::Candidate*>& candidatePtrs, const reco::Candidate* candidate)
  {
    auto found = std::find(candidatePtrs.begin(), candidatePtrs.end(), candidate);
    if ( found != candidatePtrs.end() ) return found - candidatePtrs.begin();
    else                                return -1;
  }

  std::string
  format_float(float x, int precision)
  {
    std::string format = std::string("%1.") + Form("%i", precision) + "E";
    return Form(format.c_str(), x);
  }
}

void
MyLHEWriter::analyze(const edm::Event& event, const edm::EventSetup& es)
{
  edm::Handle<LHEEventProduct> lheEvent = event.getHandle(tokenLHEEvent_);

  edm::Handle<reco::CandidateView> genParticles = event.getHandle(tokenGenParticles_);

  file1_ << "<event>\n";

  const lhef::HEPEUP& hepeup = lheEvent->hepeup();
  file1_ << std::setw(4) << genParticles->size() << " " 
         << std::setw(6) << hepeup.IDPRUP << " " 
         << format_float(hepeup.XWGTUP, 7) << " " 
         << format_float(hepeup.SCALUP, 8) << " " 
         << format_float(hepeup.AQEDUP, 8) << " " << format_float(hepeup.AQCDUP, 8) << "\n";

  std::vector<const reco::Candidate*> genParticlePtrs;
  for ( reco::CandidateView::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); ++genParticle )
  {
    genParticlePtrs.push_back(&*genParticle);
  }

  for ( reco::CandidateView::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); ++genParticle )
  {
    // indices to mother particles
    int numMothers = genParticle->numberOfMothers();
    int idxMother1 = getIdx(genParticlePtrs, genParticle->mother(0));
    int idxMother2 = getIdx(genParticlePtrs, genParticle->mother(numMothers - 1));

    // indices to daughter particles
    int numDaughters = genParticle->numberOfDaughters();
    int idxDaughter1 = getIdx(genParticlePtrs, genParticle->daughter(0));
    int idxDaughter2 = getIdx(genParticlePtrs, genParticle->daughter(numDaughters - 1));

    // CV: increase all indices to mother and daughter particles by +1
    //     in order to match the index convention used by Luca's LHE file reader
    file1_ << std::setw(8) << genParticle->pdgId() << " "
           << std::setw(2) << genParticle->status() << " "
           << std::setw(4) << 1 + idxMother1 << " " << std::setw(4) << 1 + idxMother2 << " "
           << std::setw(4) << 1 + idxDaughter1 << " " << std::setw(4) << 1 + idxDaughter2 << " "
           << format_float(genParticle->px(), 10) << " " << format_float(genParticle->py(), 10) << " " << format_float(genParticle->pz(), 10) << " "
           << format_float(genParticle->energy(), 10) << " " << format_float(genParticle->mass(), 10) << " "
           << format_float(0., 3) << " " << 0 << "\n";           
  }

  file1_ << "</event>\n";
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(MyLHEWriter);
