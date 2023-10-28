#!/bin/bash

# Example usage:
# run_mc_production_job.sh 1 1000 python/fragments/ee2tt_fragment_pythia.py /local/$USER/results aod run
# This would run a MC production job for 1000 events with seed 1 using ee2tt_fragment_pythia.py fragment.
# The resulting AOD file ends up in /local/$USER/results.
# The fubak argument instructs whether to run the job ("run") or just generate the cfg file in $PWD ("test")

if [ "$#" -lt 6 ]; then
  echo "Too few arguments";
  exit 1;
fi

SEED=$1
NEVENTS=$2
FRAGMENT=$3
DSTDIR=$4
TIER=$5
MODE=$6

TABLE="/local/karl/gridpacks/table11-11.txt.orig"
#TABLE="/local/karl/gridpacks/table11-11.txt.long" # works only if you recompile Tauola++

if [ ! -f $TABLE ]; then
  echo "No such file: $TABLE"
  exit 1;
fi

if [[ ! -d "$DSTDIR" ]]; then
  echo "No such directory: $DSTDIR";
  exit 1;
fi

JOB_NAME="job_${SEED}_$(basename $DSTDIR)";
echo "Host name is: `hostname`"

# Make a temporary directory, go there
if [ $HOSTNAME != manivald ]; then
  OUTDIR=/scratch/local/$USER/mc_production/$JOB_NAME;
else
  OUTDIR=/scratch/persistent/$USER/mc_production/$JOB_NAME;
fi

if [ "$MODE" = "run" ]; then
  rm -rfv $OUTDIR;
  mkdir -pv $OUTDIR;
  cd $OUTDIR;
  ln -s $TABLE table11-11.txt;
elif [ "$MODE" = "test" ]; then
  :
else
  echo "Invalid mode: $MODE";
  exit 1;
fi

FIRST_EVENT=$(( $NEVENTS * ($SEED - 1) + 1 ));

CUSTOM_CMDS="process.RandomNumberGeneratorService.externalLHEProducer.initialSeed=$SEED;";
CUSTOM_CMDS+="process.RandomNumberGeneratorService.generator.initialSeed=$SEED;";
CUSTOM_CMDS+="process.source.firstRun=cms.untracked.uint32(1);";
CUSTOM_CMDS+="process.source.firstLuminosityBlock=cms.untracked.uint32($SEED);";
CUSTOM_CMDS+="process.source.firstEvent=cms.untracked.uint32($FIRST_EVENT);";
CUSTOM_CMDS+="process.source.numberEventsInLuminosityBlock=cms.untracked.uint32($NEVENTS);";
CUSTOM_CMDS+="process.MessageLogger.cerr.FwkReport.reportEvery=1;";

STEP="LHE,GEN";
if [ "$TIER" = "aod" ]; then
  EVENTCONTENT="AODSIM";
  DATATIER="AODSIM";
elif [ "$TIER" = "nanoaod" ]; then
  EVENTCONTENT="NANOAODGEN";
  DATATIER="NANOAOD";
  STEP="$STEP,NANOGEN";
else
  echo "Invalid tier: $TIER";
  exit 1;
fi

OUTPUT="${TIER}sim_${SEED}.root"
CFG="cfg_${SEED}.py"

cmsDriver.py $FRAGMENT                                            \
  --fileout file:$OUTPUT --mc --eventcontent $EVENTCONTENT        \
  --datatier $DATATIER --conditions auto:mc --step LHE,GEN        \
  --no_exec --python_filename=$CFG --number=$NEVENTS --nThreads=1 \
  --customise_commands "$CUSTOM_CMDS";

if [ "$MODE" = "run" ]; then
  /usr/bin/time --verbose cmsRun $CFG
  cp -v $OUTPUT $DSTDIR

  cd -
  rm -rf $OUTDIR
fi
