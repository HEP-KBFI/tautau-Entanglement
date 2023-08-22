# tautau-Entanglement
Code for studies of quantum entanglement in Z->tautau and H->tautau events


# To checkout the code:
source /cvmfs/cms.cern.ch/cmsset_default.sh

mkdir Entanglement

cd Entanglement

export SCRAM_ARCH=slc7_amd64_gcc10

cmsrel CMSSW_12_4_8

cd CMSSW_12_4_8/src/

cmsenv

git clone https://github.com/HEP-KBFI/tautau-Entanglement TauAnalysis/Entanglement

git remote set-url origin git+ssh://git@github.com/HEP-KBFI/tautau-Entanglement

cd TauAnalysis/Entanglement

scram b -j 8

# To produce Ntuples:
source /cvmfs/cms.cern.ch/cmsset_default.sh

cd Entanglement/CMSSW_12_4_8/src/TauAnalysis/Entanglement/test

./produceNtuples.py

Execute make command as instructed in terminal window

# To run the analysis code on these Ntuples:
source /cvmfs/cms.cern.ch/cmsset_default.sh

cd Entanglement/CMSSW_12_4_8/src/TauAnalysis/Entanglement/test

./analyzeNtuples.py

Execute make command as instructed in terminal window

