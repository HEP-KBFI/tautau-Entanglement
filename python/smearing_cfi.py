import FWCore.ParameterSet.Config as cms

smearing = cms.PSet(
    applySmearing_recoil_px     = cms.bool(False),
    applySmearing_recoil_py     = cms.bool(False),
    applySmearing_recoil_pz     = cms.bool(False),
    applySmearing_recoil_energy = cms.bool(False),

    applySmearing_pv_xy         = cms.bool(False),
    applySmearing_pv_z          = cms.bool(False),

    applySmearing_track_pt      = cms.bool(False),
    applySmearing_track_theta   = cms.bool(False),
    applySmearing_track_phi     = cms.bool(False),

    applySmearing_ecal_energy   = cms.bool(False),
    applySmearing_ecal_theta    = cms.bool(False),
    applySmearing_ecal_phi      = cms.bool(False),

    applySmearing_sv_perp       = cms.bool(False),
    applySmearing_sv_parl       = cms.bool(False),

    applySmearing_tip_perp      = cms.bool(False),

    rndSeed                     = cms.uint64(0)
)
