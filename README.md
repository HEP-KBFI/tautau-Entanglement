# tautau-Entanglement
Code for studies of quantum entanglement in Z->tautau and H->tautau events


# To checkout the code:
source /cvmfs/cms.cern.ch/cmsset_default.sh

mkdir Entanglement

cd Entanglement

export SCRAM_ARCH=slc7_amd64_gcc700

cmsrel CMSSW_10_6_20

cd CMSSW_10_6_20/src/

cmsenv

git clone https://github.com/HEP-KBFI/tautau-Entanglement TauAnalysis/Entanglement

cd TauAnalysis/Entanglement

scram b -j 8

# To produce Ntuples:
source /cvmfs/cms.cern.ch/cmsset_default.sh

cd Entanglement/CMSSW_10_6_20/src/TauAnalysis/Entanglement/test

./runJobs.py


# To analyze the Ntuples and determine the spin correlation matrix C:
source /cvmfs/cms.cern.ch/cmsset_default.sh

cd Entanglement/CMSSW_10_6_20/src/TauAnalysis/Entanglement/test

analyzeEntanglementNtuple analyzeEntanglementNtuple_cfg.py
