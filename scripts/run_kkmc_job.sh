#!/bin/bash

# cmd line arguments
MODE=$1
SEED=$2
NOF_EVENTS=$3
DSTDIR=$4

if [[ ! -d "$DSTDIR" ]]; then
  echo "No such directory: $DSTDIR";
  exit 1;
fi

# Make a temporary directory, go there
JOB_NAME="job_kkmc_${SEED}_$(basename $DSTDIR)";
echo "Host name is: `hostname`"

if [ $HOSTNAME != manivald ]; then
  OUTDIR=/scratch/local/$USER/mc_production/$JOB_NAME;
else
  OUTDIR=/scratch/persistent/$USER/mc_production/$JOB_NAME;
fi
rm -rfv $OUTDIR
mkdir -pv $OUTDIR
cd $_

# Run MC production
OUTFN=${MODE}_${SEED}

singularity exec -B /local/ -B /scratch/local \
  docker://ktht/basf2-centos7:v01-12-01-1df2691 \
  /home/karl/sandbox/kkmc/run.sh $OUTFN $SEED $NOF_EVENTS $MODE

# Tarball the HepMC output file
tar czf ${OUTFN}.tgz ${OUTFN}.hepmc

# Convert the HepMC file into a ROOT file
hepmc2gen.py inputFiles=${OUTFN}.hepmc outputFile=${OUTFN}.root

# Copy the files
for ext in root txt tgz; do
  cp -v ${OUTFN}.${ext} $DSTDIR
done

# Cleanup
cd -
rm -rfv $OUTDIR
