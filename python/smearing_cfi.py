import FWCore.ParameterSet.Config as cms

smearing = cms.PSet(
    recoilSmearPx = cms.bool(False),
    recoilSmearPy = cms.bool(False),
    recoilSmearPz = cms.bool(False),
    recoilSmearE  = cms.bool(False),
    pvSmearXY     = cms.bool(False),
    pvSmearZ      = cms.bool(False),
    svSmearPerp   = cms.bool(False),
    svSmearParl   = cms.bool(False),
    tipSmearPerp  = cms.bool(False)
)
