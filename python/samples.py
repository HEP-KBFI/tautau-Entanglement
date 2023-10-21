# CV: define Standard Model expectation for matrix C, given by Eq. (69) of arXiv:2208:11723,
#     for comparison with measured matrix elements
PAR_GEN_LHC = [ 0., 0., 0., 0., 0., 0., +1., 0., 0., 0., +1., 0., 0., 0., -1. ]

# Karl: see https://github.com/HEP-KBFI/tautau-Entanglement/issues/2#issuecomment-1731618553
PAR_GEN_SUPERKEKB = [ 0., 0., 0., 0., 0., 0., -0.41989, 0., 0., 0., 0.526703, 0., 0., 0., 0.893186 ]

MAX_SUM_PHOTON_EN_LHC = 5.
MAX_SUM_PHOTON_EN_SUPERKEKB = 0.2

# CV: define samples for tau spin analysis @ LHC
samples_LHC = {
  'ggH_htt_pythia8' : {
    'inputFilePath' : '/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v3/100000/',
    'numJobs' : 10,
    'apply_evtWeight' : True
  },
  'ggH_htt_tauola' : {
   'inputFilePath' : '/store/mc/RunIISummer16MiniAODv2/GluGluToHToTauTau_M125_13TeV_amcatnloFXFX_pythia8_tauola/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/',
   'numJobs' : 10,
   'apply_evtWeight' : True
  },
  'dy_lo_pythia8' : {
   'inputFilePath' : '/store/mc/RunIISummer20UL18MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/270000/',
   'numJobs' : 10,
   'apply_evtWeight' : False
  },
  'dy_nlo_pythia8' : {
   'inputFilePath' : '/store/mc/RunIISummer20UL18MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/230000/',
   'numJobs' : 10,
   'apply_evtWeight' : False
  },
}

# CV: define samples for tau spin analysis @ SuperKEKB (Belle)
samples_SuperKEKB = {
  'dy_lo_pythia8' : {
    'inputFilePath' : '/local/karl/belle_eeToTauTau/aod/unwgt_pythia_internal4/',
    'numJobs' : 100,
    'apply_evtWeight' : False,
  },
  'dy_lo_tauola' : {
    'inputFilePath' : '/local/karl/belle_eeToTauTau/aod/unwgt_tauola_sanc_fine_v2_10M/',
    'numJobs' : 100,
    'apply_evtWeight' : False,
  },
  'dy_lo_taudecay' : {
    'inputFilePath' : '/local/karl/belle_eeToTauTau/aod/unwgt_taudecay/',
    'numJobs' : 100,
    'apply_evtWeight' : False,
  },
  'dy_lo_kkmc_orig' : {
    'inputFilePath' : '/local/karl/belle_eeToTauTau/kkmc/orig_100j_100Kevs/',
    'numJobs' : 100,
    'apply_evtWeight' : False,
  },
  'dy_lo_kkmc_bbb' : {
    'inputFilePath' : '/local/karl/belle_eeToTauTau/kkmc/bbb_100j_100Kevs/',
    'numJobs' : 100,
    'apply_evtWeight' : False,
  },
  'dy_lo_kkmc_pythia' : {
    'inputFilePath' : '/local/karl/belle_eeToTauTau/kkmc/pythia_100j_100Kevs/',
    'numJobs' : 100,
    'apply_evtWeight' : False,
  },
  'dy_lo_pythia8_ext' : {
    'inputFilePath' : '/local/karl/belle_eeToTauTau/aod/unwgt_pythia_extended/',
    'numJobs' : 2000,
    'apply_evtWeight' : False,
  },
}
