#!/bin/bash

# Example usage:
# generate_mc_production_jobs.sh 100 1000 python/fragments/ee2tt_fragment_pythia.py /local/$USER/results aod run
# Generates 100 AOD production jobs that each produce 1000 events from the pythia fragment.
# The script generates a file called submit.sh (which you can rename by adding a 7th argument to the above).
# There are three modes, specified by the 6th argument:
# - run -- runs the MC production jobs on the cluster
# - test -- generates cfg files for local testing

if [ "$#" -lt 6 ]; then
  echo "Too few arguments";
  exit 1;
fi

NOF_JOBS=$1
NOF_EVENTS_PER_JOB=$2
FRAGMENT=$(realpath $3)
DSTDIR=$4
TIER=$5
MODE=$6

SUBMISSION_SCRIPT="submit.sh";
if [ ! -z $7 ]; then
  echo "Changing submission script name from $SUBMISSION_SCRIPT to $7";
  SUBMISSION_SCRIPT=$7;
fi

if [[ ! $SUBMISSION_SCRIPT == *.sh ]]; then
  echo "Invalid name for the submission script: $SUBMISSION_SCRIPT";
  exit 1;
fi

if [ -f $SUBMISSION_SCRIPT ]; then
  echo "Submission script already exists: $SUBMISSION_SCRIPT";
  echo "Either remove the file or provide a new name as 7th argument.";
  exit 1;
fi

if [[ ! $DSTDIR == /local/* ]]; then
  echo "Invalid name for the destination directory: $DSTDIR";
  exit 1;
fi

if ! [[ "$NOF_JOBS" =~ ^[0-9]+$ && $NOF_JOBS -gt 0 ]]; then
  echo "Invalid number of jobs: $NOF_JOBS";
  exit 1;
fi

if ! [[ "$NOF_EVENTS_PER_JOB" =~ ^[0-9]+$ && $NOF_EVENTS_PER_JOB -gt 0 ]]; then
  echo "Invalid number of events per job: $NOF_EVENTS_PER_JOB";
  exit 1;
fi

if [ ! -f $FRAGMENT ]; then
  echo "No such fragment found: $3";
  exit 1;
fi

CMSSW_BASE_SRC="$CMSSW_BASE/src";
FRAGMENT_BASE="${FRAGMENT#$CMSSW_BASE_SRC/}"

mkdir -pv $DSTDIR;
if [ "$MODE" = "run" ]; then
  LOGDIR="${DSTDIR/local/home}/logs";
  mkdir -pv $LOGDIR;
elif [ "$MODE" = "test" ]; then
  :
else
  echo "Invalid mode: $MODE";
  exit 1;
fi

echo -e "#!/bin/bash\n" > $SUBMISSION_SCRIPT;
for i in `seq 1 $NOF_JOBS`; do
  if [ "$MODE" = "run" ]; then
    SUBMISSION_COMMAND="sbatch --partition=main --output=$LOGDIR/out_${i}.log "; # note the trailing space
  elif [ "$MODE" = "test" ]; then
    SUBMISSION_COMMAND="";
  fi
  echo "${SUBMISSION_COMMAND}run_mc_production_job.sh $i $NOF_EVENTS_PER_JOB $FRAGMENT_BASE $DSTDIR $TIER $MODE" >> $SUBMISSION_SCRIPT;
done
chmod +x $SUBMISSION_SCRIPT;

echo "Run: $(realpath $SUBMISSION_SCRIPT)"
