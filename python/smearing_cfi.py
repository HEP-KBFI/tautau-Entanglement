import FWCore.ParameterSet.Config as cms

smearing = cms.PSet(
    applySmearing_recoilPx = cms.bool(False),
    applySmearing_recoilPy = cms.bool(False),
    applySmearing_recoilPz = cms.bool(False),
    applySmearing_recoilE  = cms.bool(False),
    applySmearing_pvXY     = cms.bool(False),
    applySmearing_pvZ      = cms.bool(False),
    applySmearing_svPerp   = cms.bool(False),
    applySmearing_svParl   = cms.bool(False),
    applySmearing_tipPerp  = cms.bool(False)
)
