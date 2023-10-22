# tautau-Entanglement

Code for studies of quantum entanglement in Z->tautau and H->tautau events

## Setup

```bash
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
```

## Ntuple production

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd Entanglement/CMSSW_12_4_8/src/TauAnalysis/Entanglement
./test/produceNtuples.py -v 2023Oct06_wSmearing -s dy_lo_pythia8 -j local # or -j cluster
```

Execute make command as instructed in terminal window

## Run the analysis code on the Ntuples

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd Entanglement/CMSSW_12_4_8/src/TauAnalysis/Entanglement
# run inclusive measurement + scan in |cos(theta)| on all decay modes (except had_had) on the cluster
./test/analyzeNtuples.py -v 2023Oct22        -V 2023Oct06_wSmearing -s dy_lo_pythia8_ext -j cluster \
  -M inclusive scan
# run inclusive measurement + scan in |cos(theta)| on fully hadronic decay modes (had_had) locally
./test/analyzeNtuples.py -v 2023Oct22_hadHad -V 2023Oct06_wSmearing -s dy_lo_pythia8_ext -j local   \
  -M inclusive scan -d had_had
```

## Plot the results

```bash
plot_entanglement_vars.py -i ~/Entanglement/analysis/SuperKEKB/2023Oct22 \
  -o ~/Entanglement/analysis/SuperKEKB/2023Oct22/json_plots -E png
```

## Latex tables

```bash
# note: multiple arguments
get_latex_tables.py -i ~/Entanglement/analysis/SuperKEKB/2023Oct22 \
                       ~/Entanglement/analysis/SuperKEKB/2023Oct22_hadHad
# to produce tables for a particular cut in |cos(theta)|
get_latex_tables.py -i ~/Entanglement/analysis/SuperKEKB/2023Oct22 \
                       ~/Entanglement/analysis/SuperKEKB/2023Oct22_hadHad \
                    -C absCosTheta0p55 # or "opt" if you have run the analysis with -M opt
```