import FWCore.ParameterSet.Config as cms

resolutions = cms.PSet(
    recoilResolution_px     = cms.double(1.),     # [GeV]
    recoilResolution_py     = cms.double(1.),     # [GeV]
    recoilResolution_pz     = cms.double(5.),     # [GeV]
    recoilResolution_energy = cms.double(5.1),     # [GeV]

    pvResolution_xy         = cms.double(0.0015), # [cm]
    pvResolution_z          = cms.double(0.0030), # [cm]
    
    trackResolution_pt      = cms.double(2.e-5),  # resolution on 1/pT in units of GeV^-1
    trackResolution_theta   = cms.double(3.e-4),  # [rad]
    trackResolution_phi     = cms.double(3.e-4),  # [rad]

    ecalResolution_energy_a = cms.double(0.166),  # coefficient a in resolution function sigma_E = a*sqrt(E) + b*E, where E is in units of GeV
    ecalResolution_energy_b = cms.double(0.011),  # coefficient b in resolution function sigma_E = a*sqrt(E) + b*E, where E is in units of GeV
    ecalResolution_theta    = cms.double(1.7e-3), # [rad]
    ecalResolution_phi      = cms.double(1.7e-3), # [rad]

    svResolution_parl       = cms.double(0.1000), # [cm]
    svResolution_perp       = cms.double(0.0020), # [cm]

    tipResolution_perp      = cms.double(0.0020)  # [cm]
)
