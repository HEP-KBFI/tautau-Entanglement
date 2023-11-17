import FWCore.ParameterSet.Config as cms

# CV: acceptance cuts for LHC
#     taken from config files for hadronic tau reconstruction in CMS
acceptanceCuts_LHC = cms.PSet(
  chargedHadrons = cms.PSet(
      minEta = cms.double(-2.5),
      maxEta = cms.double(+2.5),
      minPt  = cms.double( 1.)
  ),
  photons = cms.PSet(
      minEta = cms.double(-2.5),
      maxEta = cms.double(+2.5),
      minPt  = cms.double( 0.5)
  )
)

# CV: acceptance cuts for charged hadrons and for photons at Belle II
#     taken from the paper arXiv:2010.15361
acceptanceCuts_SuperKEKB = cms.PSet(
  chargedHadrons = cms.PSet(
      minTheta = cms.double(0.2967), #  17 degrees
      maxTheta = cms.double(2.6180), # 150 degrees
      minPt  = cms.double(0.1)       # [GeV]
  ),
  photons = cms.PSet(
      minTheta = cms.double(0.2967), #  17 degrees
      maxTheta = cms.double(2.6180), # 150 degrees
      minEnergy  = cms.double(0.1)   # [GeV]
  )
)
