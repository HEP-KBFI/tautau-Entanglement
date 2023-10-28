# MC production

## Generating gridpacks

Execute the following outside of CMSSW:

```bash
git clone -b ee2tt https://github.com/HEP-KBFI/genproductions.git
cd genproductions/bin/MadGraph5_aMCatNLO

# to generate e+ e- > ta+ ta- without decaying the taus
./gridpack_generation.sh ee2tt cards/production/ee2tt

# to generate e+ e- > ta+ ta- and decay the taus with TauDecay
./gridpack_generation.sh ee2tt_TauDecay cards/production/ee2tt_TauDecay
```

Copy the resulting `*.tar.xz` files to some common directory and edit the `GRIDPACKS` dictionary in `tools.py` in the current directory.

### TauDecay commands

The [`TauDecay`](https://arxiv.org/abs/1212.6247) commands are generated with the following snippet:

```python
dms = {
  'pi'   : 'pi+ vt~',
  'rho'  : 'pi+ vt~ pi0',
  'a1_1' : 'pi+ vt~ pi0 pi0',
  'a1_2' : 'pi+ vt~ pi+ pi-', # requires EFT=1 at the very end
}

cmds = []
for dm_tauPlus in dms:
  cmd_tauPlus = dms[dm_tauPlus]
  for dm_tauMinus in dms:
    cmd_tauMinus = dms[dm_tauMinus].replace('~', '').translate(str.maketrans('+-', '-+'))
    pfix = 'add process' if cmds else 'generate'
    sfix = 'EFT=1' if dm_tauPlus == 'a1_2' or dm_tauMinus == 'a1_2' else ''
    cmd = f'{pfix} e+ e- > ta+ ta-, (ta+ > {cmd_tauPlus}), (ta- > {cmd_tauMinus}) {sfix}'
    cmds.append(cmd.lstrip())

assert(len(cmds) == len(dms)**2)
print('\n'.join(cmds))
```

## Running MC sample production

### MadGraph + Pythia or Tauola

Example commands:

```bash
# Pythia, 200M events
generate_mc_production_jobs.sh 2000 100000 ee2tt_fragment_pythia.py \
  /local/$USER/belle_eeToTauTau/aod/unwgt_pythia_extended aod run

# Tauola, 10M events
generate_mc_production_jobs.sh 100 100000 ee2tt_fragment_tauola.py \
  /local/$USER/belle_eeToTauTau/aod/unwgt_tauola aod run

# TauDecay, 10M events
generate_mc_production_jobs.sh 100 100000 ee2tt_fragment_TauDecay.py \
  /local/$USER/belle_eeToTauTau/aod/unwgt_taudecay aod run
```

Note that if you run these commands in some other directory, you'd have to modify the path to the Pythia fragment (i.e., the 3rd argument) accordingly.
Each of those commands creates a script called `submit.sh` (or whatever you would specify as the 7th argument), which you run in order to submit the jobs to the cluster.
If you just want to generate config files but not run the jobs, then replace `run` with `test` in the above examples.

NB! If you're running Tauola and would want to increase the number of bins in cosTheta to reduce errors caused by the interpolation, do the following (assuming that you have CMSSW already set up):

```bash
cd $HOME # or some other place outside of $CMSSW_BASE

# Download and unpack Tauola++
wget http://tauolapp.web.cern.ch/tauolapp/resources/TAUOLA.1.1.8/TAUOLA.1.1.8-LHC.tar.gz
tar xzf TAUOLA.1.1.8-LHC.tar.gz
cd TAUOLA

# Set up the environment
cd $CMSSW_BASE
export HEPMC_ROOT=$(scram tool info hepmc | grep "^HEPMC_BASE=" | tr '=' ' ' | awk '{print $2}')
export LHAPDF_ROOT=$(scram tool info lhapdf | grep "^LHAPDF_BASE=" | tr '=' ' ' | awk '{print $2}')
export PYTHIA8_ROOT=$(scram tool info pythia8 | grep "^PYTHIA8_BASE=" | tr '=' ' ' | awk '{print $2}')
export BOOST_ROOT=$(scram tool info boost | grep "^BOOST_BASE=" | tr '=' ' ' | awk '{print $2}')
cd -

# Build
CPPFLAGS="-I${BOOST_ROOT}/include" ./configure \
  --prefix=$PWD/prefix \
  --without-hepmc3 \
  --with-hepmc=$HEPMC_ROOT \
  --with-pythia8=$PYTHIA8_ROOT \
  --with-lhapdf=$LHAPDF_ROOT
make -j8

# Install
make install
cp -v $PWD/prefix/lib/* $CMSSW_BASE/lib/$SCRAM_ARCH

# Configure CMSSW accordingly
cd $CMSSW_BASE/src
scram setup tauolapp
scram build ExternalLinks ToolUpdated
scram build echo_tauolapp_USED_BY -j8
scram build echo_tauolapp_USES
```

If you want to use finer binning for the SANC tables, then
1. Recompile Tauola++ by changing `NCOS` to `200` in `src/tauolaCInterfaces/Tauola.h`:
```cpp
static const int NS1=100,NS2=100,NS3=100,NCOS=21; 
```
2. Uncomment this line in `run_mc_production_job.sh`:
```bash
#TABLE="/local/karl/gridpacks/table11-11.txt.long" # works only if you recompile Tauola++
```

#### Generating the SANC tables

Assuming that you have already downloaded Tauola++ (but necessarily compiled or installed) as instructed in the previous section, you can produce your own SANC tables as follows:

```bash
# While you're still in the TAUOLA directory
cd SANC

# Apply changes to the makefiles
patch Makefile $CMSSW_BASE/src/TauAnalysis/Entanglement/data/sanc.patch
patch LoopTools-2.1/makefile $CMSSW_BASE/src/TauAnalysis/Entanglement/data/looptools.patch

# Build
make
make tables

# Create the SANC tables
./SANCinterface.exe
```

### KKMC production

Create a script that submits KKMC production jobs to the cluster:

```bash
# Create a script called submit.sh that submits 10 KKMC production jobs to SLURM that each generate 10000 events.
# The output files will be stored in the directory that's specified as the 4th argument in the following.
create_kkmc_jobs.py default 10 10000 /local/$USER/belle_eeToTauTau/kkmc/default_10jobs_10Kevents run
```

This would generate a file called `submit.sh` (or whatever you provide as the 5th argument to `create_kkmc_jobs.py`), which you run in order to submit the jobs.
You can replace `default` with `orig` or `bbb` (read [this](https://github.com/HEP-KBFI/tautau-Entanglement/issues/2#issuecomment-1746796984) to learn more about their differences).

The jobs can fail if they are unable to fetch payloads via cURL, which does happen quite often.
If you just want to run the production interactively, use `test` instead of `run` as the 5th argument (cf the above example).
