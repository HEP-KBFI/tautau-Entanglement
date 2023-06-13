import FWCore.ParameterSet.Config as cms

resolutions = cms.PSet(
    recoilResolutionPx = cms.double(5.),     # [GeV]
    recoilResolutionPy = cms.double(5.),     # [GeV]
    recoilResolutionPz = cms.double(5.),     # [GeV]
    recoilResolutionE  = cms.double(5.),     # [GeV]
    pvResolutionXY     = cms.double(0.0015), # [cm]
    pvResolutionZ      = cms.double(0.0030), # [cm]
    svResolutionParl   = cms.double(0.1000), # [cm]
    svResolutionPerp   = cms.double(0.0020), # [cm]
    tipResolutionPerp  = cms.double(0.0020)  # [cm]
)
