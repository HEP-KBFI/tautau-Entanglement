
# CV: define samples for tau spin analysis @ LHC
samples_LHC = {
  'ggH_htt_pythia8' : {
    'inputFilePath' : '/store/mc/RunIISummer20UL18MiniAODv2/GluGluHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v3/100000/',
    'numJobs' : 10,
    'process' : "ggH_htt_pythia8",
    'par_gen' : [ 0., 0., 0., 0., 0., 0., +1., 0., 0., 0., +1., 0., 0., 0., -1. ],
    'apply_evtWeight' : True
  },
  ##'ggH_htt_tauola' : {
  ##  'inputFilePath' : '/store/mc/RunIISummer16MiniAODv2/GluGluToHToTauTau_M125_13TeV_amcatnloFXFX_pythia8_tauola/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/',
  ##  'numJobs' : 10,
  ##  'process' : "ggH_htt_tauola",
  ##  'par_gen' : [ 0., 0., 0., 0., 0., 0., +1., 0., 0., 0., +1., 0., 0., 0., -1. ],
  ##  'apply_evtWeight' : True
  ##},
  #'dy_lo_pythia8' : {
  #  'inputFilePath' : '/store/mc/RunIISummer20UL18MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/270000/',
  #  'numJobs' : 10,
  #  'process' : "dy_lo_pythia8",
  #  'par_gen' : [ 0., 0., 0., 0., 0., 0., +1., 0., 0., 0., +1., 0., 0., 0., -1. ],
  #  'apply_evtWeight' : False
  #},
  #'dy_nlo_pythia8' : {
  #  'inputFilePath' : '/store/mc/RunIISummer20UL18MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/230000/',
  #  'numJobs' : 10,
  #  'process' : "dy_nlo_pythia8",
  #  'par_gen' : [ 0., 0., 0., 0., 0., 0., +1., 0., 0., 0., +1., 0., 0., 0., -1. ],
  #  'apply_evtWeight' : False
  #},
}

# CV: define samples for tau spin analysis @ SuperKEKB (Belle)
samples_SuperKEKB = {
  'dy_lo_pythia8' : {
    'inputFilePath' : '/local/karl/belle_eeToTauTau/aod/unwgt_pythia_internal4/',
    'numJobs' : 100,
    'process' : "dy_lo_pythia8",
    'par_gen' : [ 0., 0., 0., 0., 0., 0., -0.249686, 0.00508408, 0.00391348, 0.00653464, 0.138683, -0.276996, -0.00186421, 0.243658, 0.885224 ],
    'apply_evtWeight' : False,
  },
  'dy_lo_tauola' : {
    'inputFilePath' : '/local/karl/belle_eeToTauTau/aod/unwgt_tauola_sanc_fine/',
    'numJobs' : 100,
    'process' : "dy_lo_tauola",
    'par_gen' : [ 0., 0., 0., 0., 0., 0., -0.249686, 0.00508408, 0.00391348, 0.00653464, 0.138683, -0.276996, -0.00186421, 0.243658, 0.885224 ],
    'apply_evtWeight' : False,
  },
  'dy_lo_taudecay' : {
    'inputFilePath' : '/local/karl/belle_eeToTauTau/aod/unwgt_taudecay/',
    'numJobs' : 100,
    'process' : "dy_lo_taudecay",
    'par_gen' : [ 0., 0., 0., 0., 0., 0., -0.249686, 0.00508408, 0.00391348, 0.00653464, 0.138683, -0.276996, -0.00186421, 0.243658, 0.885224 ],
    'apply_evtWeight' : False,
  },
}
