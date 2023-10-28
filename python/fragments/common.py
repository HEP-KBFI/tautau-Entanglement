from TauAnalysis.Entanglement.fragments.tools import COM_ENERGY_BELLE2, read_pythia_cfg
import FWCore.ParameterSet.Config as cms

externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
  args = cms.vstring(),
  nEvents = cms.untracked.uint32(5000),
  numberOfParameters = cms.uint32(1),
  outputFile = cms.string('cmsgrid_final.lhe'),
  scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh')
)
generator = cms.EDFilter(
  "Pythia8HadronizerFilter",
  maxEventsToPrint = cms.untracked.int32(1),
  pythiaPylistVerbosity = cms.untracked.int32(1),
  filterEfficiency = cms.untracked.double(1.0),
  pythiaHepMCVerbosity = cms.untracked.bool(True),
  comEnergy = cms.double(COM_ENERGY_BELLE2),
  PythiaParameters = cms.PSet(
    processParameters = cms.vstring(read_pythia_cfg()),
    parameterSets = cms.vstring('processParameters')
  )
)
