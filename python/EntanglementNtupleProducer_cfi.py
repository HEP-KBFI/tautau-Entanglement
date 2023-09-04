import FWCore.ParameterSet.Config as cms

from TauAnalysis.Entanglement.smearing_cfi import smearing
ntupleProducer = cms.EDAnalyzer("EntanglementNtupleProducer",
    src = cms.InputTag(''),
    collider = cms.string(""),
    hAxis = cms.string(""),
    resolutions = cms.PSet(),
    smearing = smearing.clone(
        rndSeed = cms.uint64(1)
    ),
    #applySmearing = cms.bool(False),
    applySmearing = cms.bool(True),
    startPosFinder = cms.PSet(
        algos = cms.vint32(1),
        applyHiggsMassConstraint = cms.bool(True),
        applyRecoilEnergy_and_PzConstraint = cms.bool(True),
        skip = cms.bool(False)
    ),
    kinematicFit = cms.PSet(
        applyLifetimeConstraint = cms.bool(False),
        skip = cms.bool(False)
    ),
    srcEvtWeights = cms.VInputTag('genWeight'),
    verbosity = cms.untracked.int32(-1),
    #verbosity = cms.untracked.int32(3),
    cartesian = cms.untracked.bool(True)
)
